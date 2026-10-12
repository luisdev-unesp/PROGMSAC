#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <omp.h>
#include "PROGMSAC.h"

/* =========================================================================
   LEITURA E GERENCIAMENTO DE SEQUÊNCIAS
   ========================================================================= */

char** ler_sequencias(int quantidade, int tam_max) {
    char **seqs = malloc(sizeof(char*) * quantidade);
    if (!seqs) return NULL;

    for (int i = 0; i < quantidade; i++) {
        seqs[i] = malloc(sizeof(char) * (tam_max + 1));
        if (!seqs[i]) {
            for (int k = 0; k < i; k++) free(seqs[k]);
            free(seqs);
            return NULL;
        }
        if (scanf("%1000010s", seqs[i]) != 1) {
            printf("\nErro ao ler a sequencia %d\n", i + 1);
            for (int k = 0; k <= i; k++) free(seqs[k]);
            free(seqs);
            return NULL;
        }
    }
    return seqs;
}

int carregar_matriz(matriz_sub* m) {
    char nome_arquivo[100];
    printf("\nde o nome do arquivo:");
    
    if (scanf("%99s", nome_arquivo) != 1) {
        printf("\nErro ao ler o nome do arquivo.\n");
        return 0;
    }

    FILE *arquivo = fopen(nome_arquivo, "r");
    if (arquivo == NULL) {
        printf("\nerro ao abrir o arquivo");
        return 0;
    }

    for (int i = 0; i < 256; i++) {
        m->hash[i] = -1;
    }

    char buffer[1024];

    while (fgets(buffer, sizeof(buffer), arquivo)) {
        char *ptr = buffer;
        while (*ptr == ' ' || *ptr == '\t') ptr++;

        if (*ptr != '#' && *ptr != '\n' && *ptr != '\r' && *ptr != '\0') {
            break;
        }
    }

    int qtd = 0;
    char *token = strtok(buffer, " \t\r\n");
    while (token != NULL) {
        unsigned char c = (unsigned char)token[0];
        m->hash[c] = qtd++;
        token = strtok(NULL, " \t\r\n");
    }

    m->tam = qtd;

    if (m->tam == 0) {
        printf("Cabecalho invalido ou nao encontrado.\n");
        fclose(arquivo);
        return 0;
    }

    m->valores = malloc(m->tam * sizeof(int*));
    for (int i = 0; i < m->tam; i++) {
        m->valores[i] = malloc(m->tam * sizeof(int));
    }

    for (int i = 0; i < m->tam; i++) {
        char rotulo_linha[16];
        if (fscanf(arquivo, "%15s", rotulo_linha) != 1) {
            printf("\nErro ao ler rotulo da linha %d\n", i);
            fclose(arquivo);
            return 0;
        }

        for (int j = 0; j < m->tam; j++) {
            if (fscanf(arquivo, "%d", &m->valores[i][j]) != 1) {
                printf("\nErro ao ler valor da matriz em [%d][%d]\n", i, j);
                fclose(arquivo);
                return 0;
            }
        }
    }

    fclose(arquivo);
    return 1;   
}

int indice_hash(char amino, matriz_sub* m) {
    return m->hash[(unsigned char)amino];
}

int obter_score(char a, char b, int flag_proteina, const matriz_sub* m) {
    char ca = toupper((unsigned char)a);
    char cb = toupper((unsigned char)b);

    if (flag_proteina) {
        int indice_1 = m->hash[(unsigned char)ca];
        int indice_2 = m->hash[(unsigned char)cb];
        if (indice_1 != -1 && indice_2 != -1) {
            return m->valores[indice_1][indice_2]; 
        }
        return mismatch;
    }
    return (ca == cb) ? match : mismatch;
}

/* =========================================================================
   MATRIZ DE DISTÂNCIA NORMALIZADA
   ========================================================================= */

double** calcular_matriz_distancia(char **seqs, int qtd, int flag_proteina, const matriz_sub *m) {
    double **matriz_dist = malloc(sizeof(double*) * qtd);
    for (int i = 0; i < qtd; i++) {
        matriz_dist[i] = malloc(sizeof(double) * qtd);
    }

    double *max_scores = malloc(sizeof(double) * qtd);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < qtd; i++) {
        max_scores[i] = (double)calcular_score_gotoh_1d(seqs[i], seqs[i], flag_proteina, m);
    }

    int total_pares = (qtd * (qtd - 1)) / 2;
    int pares_concluidos = 0;

    printf("Calculando matriz de distância (%d pares em paralelo)...\n", total_pares);
    fflush(stdout);

    #pragma omp parallel for schedule(dynamic) shared(pares_concluidos)
    for (int i = 0; i < qtd; i++) {
        matriz_dist[i][i] = 0.0;
        for (int j = i + 1; j < qtd; j++) {
            int score = calcular_score_gotoh_1d(seqs[i], seqs[j], flag_proteina, m);

            double max_s = (max_scores[i] > max_scores[j]) ? max_scores[i] : max_scores[j];
            double dist = max_s - (double)score;
            if (dist < 0.0) dist = 0.0;

            matriz_dist[i][j] = dist;
            matriz_dist[j][i] = dist;

            #pragma omp atomic
            pares_concluidos++;

            if (pares_concluidos % 25 == 0 || pares_concluidos == total_pares) {
                #pragma omp critical
                {
                    printf("Progresso: %d / %d pares (%.1f%%)\n",
                           pares_concluidos, total_pares,
                           ((double)pares_concluidos / total_pares) * 100.0);
                    fflush(stdout);
                }
            }
        }
    }

    free(max_scores);
    return matriz_dist;
}

int calcular_score_gotoh_1d(const char *seq1, const char *seq2, int flag_proteina, const matriz_sub *m) {
    int tam_1 = strlen(seq1);
    int tam_2 = strlen(seq2);
    int INF = 1e9;

    int *M_prev = malloc((tam_2 + 1) * sizeof(int));
    int *D_prev = malloc((tam_2 + 1) * sizeof(int));
    int *I_prev = malloc((tam_2 + 1) * sizeof(int));

    int *M_curr = malloc((tam_2 + 1) * sizeof(int));
    int *D_curr = malloc((tam_2 + 1) * sizeof(int));
    int *I_curr = malloc((tam_2 + 1) * sizeof(int));

    M_prev[0] = 0;
    D_prev[0] = -INF;
    I_prev[0] = -INF;

    for (int j = 1; j <= tam_2; j++) {
        M_prev[j] = gapopen + (j - 1) * gapextend;
        D_prev[j] = -INF;
        I_prev[j] = gapopen + (j - 1) * gapextend;
    }

    for (int i = 1; i <= tam_1; i++) {
        M_curr[0] = gapopen + (i - 1) * gapextend;
        D_curr[0] = gapopen + (i - 1) * gapextend;
        I_curr[0] = -INF;

        for (int j = 1; j <= tam_2; j++) {
            int custo_troca = obter_score(seq1[i - 1], seq2[j - 1], flag_proteina, m);

            int d1 = M_prev[j] + gapopen + gapextend;
            int d2 = D_prev[j] + gapextend;
            D_curr[j] = (d1 > d2) ? d1 : d2;

            int i1 = M_curr[j - 1] + gapopen + gapextend;
            int i2 = I_curr[j - 1] + gapextend;
            I_curr[j] = (i1 > i2) ? i1 : i2;

            int m_diag = M_prev[j - 1] + custo_troca;
            int maior = m_diag;
            if (D_curr[j] > maior) maior = D_curr[j];
            if (I_curr[j] > maior) maior = I_curr[j];
            M_curr[j] = maior;
        }

        int *tmp;
        tmp = M_prev; M_prev = M_curr; M_curr = tmp;
        tmp = D_prev; D_prev = D_curr; D_curr = tmp;
        tmp = I_prev; I_prev = I_curr; I_curr = tmp;
    }

    int score_final = M_prev[tam_2];

    free(M_prev); free(D_prev); free(I_prev);
    free(M_curr); free(D_curr); free(I_curr);

    return score_final;
}

/* =========================================================================
   NEIGHBOR-JOINING (ÁRVORE FILOGENÉTICA)
   ========================================================================= */

double* neighbor_step1(double **matriz, int n, int *ativo) {
    double *vetor = malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) {
        vetor[i] = 0.0;
        if (!ativo[i]) continue;
        for (int j = 0; j < n; j++) {
            if (i != j && ativo[j]) { 
                vetor[i] += matriz[i][j];
            }
        } 
        vetor[i] /= (n > 2) ? (n - 2) : 1;
    }
    return vetor;
}

void neighbor_step2(double *vetor, int n, double **matriz, int *ativo, int *idx_matriz_u, int *idx_matriz_v) {
    double menor_elem = 1e9; 
    *idx_matriz_u = -1;
    *idx_matriz_v = -1;

    for (int i = 0; i < n; i++) {
        if (!ativo[i]) continue;
        for (int j = i + 1; j < n; j++) {
            if (!ativo[j]) continue;
            double q = (matriz[i][j] - vetor[i] - vetor[j]);
            if (q < menor_elem) {
                menor_elem = q;
                *idx_matriz_u = i;
                *idx_matriz_v = j;
            }
        }
    }
}

void neighbor_step3(No **nos, double **matriz, double *vetor, int i, int j, int novo_id) {
    double dist_i = 0.5 * (matriz[i][j] + (vetor[i] - vetor[j]));
    double dist_j = 0.5 * (matriz[i][j] + (vetor[j] - vetor[i]));

    if (dist_i < 0) dist_i = 0.0;
    if (dist_j < 0) dist_j = 0.0;

    nos[novo_id] = criar_no_interno(novo_id, nos[i], nos[j], dist_i, dist_j);
}

No* arvore_neighbor(double **dist_inicial, int n) {
    int maximo_nos = n * 2;

    double **matriz_dist_mod = malloc(sizeof(double*) * maximo_nos);
    for (int i = 0; i < maximo_nos; i++) {
        matriz_dist_mod[i] = calloc(maximo_nos, sizeof(double));
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matriz_dist_mod[i][j] = dist_inicial[i][j];
        }
    }

    No **nos = malloc(sizeof(No*) * maximo_nos);
    int *ativo = calloc(maximo_nos, sizeof(int));
    for (int i = 0; i < n; i++) {
        nos[i] = criar_folha(i, i);
        ativo[i] = 1;
    }

    int r = n;
    int proximo_id = n;
    while (r > 2) {
        double *vet_aux = neighbor_step1(matriz_dist_mod, proximo_id, ativo);

        int u, v;
        neighbor_step2(vet_aux, proximo_id, matriz_dist_mod, ativo, &u, &v);

        neighbor_step3(nos, matriz_dist_mod, vet_aux, u, v, proximo_id);

        double dist_uv = matriz_dist_mod[u][v];
        for (int k = 0; k < proximo_id; k++) {
            if (ativo[k] && k != u && k != v) {
                double dist_novo = 0.5 * (matriz_dist_mod[u][k] + matriz_dist_mod[v][k] - dist_uv);
                if (dist_novo < 0) dist_novo = 0.0;
                matriz_dist_mod[proximo_id][k] = dist_novo;
                matriz_dist_mod[k][proximo_id] = dist_novo;
            }
        }

        ativo[u] = 0;
        ativo[v] = 0;
        ativo[proximo_id] = 1;

        proximo_id++;
        r--;
        free(vet_aux);
    }

    int u = -1, v = -1;
    for (int i = 0; i < proximo_id; i++) {
        if (ativo[i]) {
            if (u == -1) u = i;
            else v = i;
        }
    }

    double dist_final = matriz_dist_mod[u][v] / 2.0;
    No *raiz = criar_no_interno(proximo_id, nos[u], nos[v], dist_final, dist_final);

    for (int i = 0; i < maximo_nos; i++) {
        free(matriz_dist_mod[i]);
    }
    free(matriz_dist_mod);
    free(nos);
    free(ativo);

    return raiz;
}

No* criar_folha(int id, int seq_idx) {
    No *n = malloc(sizeof(No));
    n->id = id;
    n->seq_idx = seq_idx;
    n->dist_pai = 0.0;
    n->esq = NULL;
    n->dir = NULL;
    return n;
}

No* criar_no_interno(int id, No *filho_esq, No *filho_dir, double dist_esq, double dist_dir) {
    No *n = malloc(sizeof(No));
    n->id = id;
    n->seq_idx = -1;
    n->dist_pai = 0.0;
    n->esq = filho_esq;
    n->dir = filho_dir;
    filho_esq->dist_pai = dist_esq;
    filho_dir->dist_pai = dist_dir;
    return n;
}

/* =========================================================================
   CÁLCULO OTIMIZADO DE SCORE DE COLUNA
   ========================================================================= */

double calculo_score_coluna(Perfil *p1, Perfil *p2, int col1, int col2, const matriz_sub *m, int flag_proteina) {
    double score_total = 0.0;

    for (int i = 0; i < p1->num_seqs; i++) {
        char c1 = p1->sequencias[i][col1];
        for (int j = 0; j < p2->num_seqs; j++) {
            char c2 = p2->sequencias[j][col2];
            
            if (c1 == '-' && c2 == '-') {
                score_total += 0.0; 
            } else if (c1 == '-' || c2 == '-') {
                score_total += gapextend; 
            } else {
                score_total += obter_score(c1, c2, flag_proteina, m);
            }
        }
    }
    return score_total;
}

/* =========================================================================
   ALINHAMENTO E FUSÃO DE PERFIS (OTIMIZADO)
   ========================================================================= */

Perfil* fundir_perfil(Perfil *perfil1, Perfil *perfil2, const matriz_sub *m, int flag_proteina) {
    int l1 = perfil1->comprimento;
    int l2 = perfil2->comprimento;
    int n1 = perfil1->num_seqs;
    int n2 = perfil2->num_seqs;

    int *counts1 = calloc((size_t)l1 * 256, sizeof(int));
    int *counts2 = calloc((size_t)l2 * 256, sizeof(int));

    unsigned char chars1[256], chars2[256];
    int num_c1 = 0, num_c2 = 0;
    unsigned char pres1[256] = {0}, pres2[256] = {0};

    for (int s = 0; s < n1; s++) {
        for (int col = 0; col < l1; col++) {
            unsigned char c = (unsigned char)perfil1->sequencias[s][col];
            counts1[col * 256 + c]++;
            if (!pres1[c]) { pres1[c] = 1; chars1[num_c1++] = c; }
        }
    }
    for (int s = 0; s < n2; s++) {
        for (int col = 0; col < l2; col++) {
            unsigned char c = (unsigned char)perfil2->sequencias[s][col];
            counts2[col * 256 + c]++;
            if (!pres2[c]) { pres2[c] = 1; chars2[num_c2++] = c; }
        }
    }

    size_t rows = l1 + 1;
    size_t cols = l2 + 1;

    double *M_data = malloc(rows * cols * sizeof(double));
    double *D_data = malloc(rows * cols * sizeof(double));
    double *I_data = malloc(rows * cols * sizeof(double));

    double **M = malloc(rows * sizeof(double*));
    double **D = malloc(rows * sizeof(double*));
    double **I = malloc(rows * sizeof(double*));

    for (size_t i = 0; i < rows; i++) {
        M[i] = M_data + i * cols;
        D[i] = D_data + i * cols;
        I[i] = I_data + i * cols;
    }

    M[0][0] = 0;
    D[0][0] = -1e9;
    I[0][0] = -1e9;

    for (int i = 1; i <= l1; i++) { 
        double val = (gapopen + (i - 1) * gapextend) * n1 * n2;
        M[i][0] = val; D[i][0] = val; I[i][0] = -1e9;
    }
    for (int j = 1; j <= l2; j++) { 
        double val = (gapopen + (j - 1) * gapextend) * n1 * n2;
        M[0][j] = val; D[0][j] = -1e9; I[0][j] = val;
    }

    for (int i = 1; i <= l1; i++) {
        int *c1_ptr = &counts1[(i - 1) * 256];
        for (int j = 1; j <= l2; j++) {
            int *c2_ptr = &counts2[(j - 1) * 256];
            double score_col = 0.0;

            for (int k1 = 0; k1 < num_c1; k1++) {
                unsigned char char1 = chars1[k1];
                int cnt1 = c1_ptr[char1];
                if (cnt1 == 0) continue;

                for (int k2 = 0; k2 < num_c2; k2++) {
                    unsigned char char2 = chars2[k2];
                    int cnt2 = c2_ptr[char2];
                    if (cnt2 == 0) continue;

                    double s_pair;
                    if (char1 == '-' && char2 == '-') s_pair = 0.0;
                    else if (char1 == '-' || char2 == '-') s_pair = gapextend;
                    else s_pair = obter_score(char1, char2, flag_proteina, m);

                    score_col += (double)(cnt1 * cnt2) * s_pair;
                }
            }

            double d1 = M[i - 1][j] + (gapopen + gapextend) * n1 * n2;
            double d2 = D[i - 1][j] + gapextend * n1 * n2;
            D[i][j] = (d1 > d2) ? d1 : d2;

            double i1 = M[i][j - 1] + (gapopen + gapextend) * n1 * n2;
            double i2 = I[i][j - 1] + gapextend * n1 * n2;
            I[i][j] = (i1 > i2) ? i1 : i2;

            double m_diag = M[i - 1][j - 1] + score_col;
            double maior = m_diag;
            if (D[i][j] > maior) maior = D[i][j];
            if (I[i][j] > maior) maior = I[i][j];
            M[i][j] = maior;
        }
    }

    int i = l1, j = l2;
    int capacidade_max = l1 + l2 + 10;
    int total_seqs = n1 + n2;
    
    char **temp_seqs = malloc(sizeof(char*) * total_seqs);
    for (int s = 0; s < total_seqs; s++) {
        temp_seqs[s] = malloc(capacidade_max);
        temp_seqs[s][0] = '\0';
    }

    int idx_str = 0;
    while (i > 0 || j > 0) {
        double score_col = 0.0;
        if (i > 0 && j > 0) {
            int *c1_ptr = &counts1[(i - 1) * 256];
            int *c2_ptr = &counts2[(j - 1) * 256];
            for (int k1 = 0; k1 < num_c1; k1++) {
                unsigned char char1 = chars1[k1];
                int cnt1 = c1_ptr[char1];
                if (cnt1 == 0) continue;
                for (int k2 = 0; k2 < num_c2; k2++) {
                    unsigned char char2 = chars2[k2];
                    int cnt2 = c2_ptr[char2];
                    if (cnt2 == 0) continue;
                    double s_pair;
                    if (char1 == '-' && char2 == '-') s_pair = 0.0;
                    else if (char1 == '-' || char2 == '-') s_pair = gapextend;
                    else s_pair = obter_score(char1, char2, flag_proteina, m);
                    score_col += (double)(cnt1 * cnt2) * s_pair;
                }
            }
        }

        if (i > 0 && j > 0 && M[i][j] == M[i - 1][j - 1] + score_col) {
            for (int s = 0; s < n1; s++){ 
                temp_seqs[s][idx_str] = perfil1->sequencias[s][i - 1];
            }
            for (int s = 0; s < n2; s++){ 
                temp_seqs[n1 + s][idx_str] = perfil2->sequencias[s][j - 1];
            }
            i--; 
            j--;
        } else if (i > 0 && M[i][j] == D[i][j]) {
            for (int s = 0; s < n1; s++){ 
                temp_seqs[s][idx_str] = perfil1->sequencias[s][i - 1];
            }
            for (int s = 0; s < n2; s++){ 
                temp_seqs[n1 + s][idx_str] = '-'; 
            }
            i--;
        } else {
            for (int s = 0; s < n1; s++){ 
                temp_seqs[s][idx_str] = '-'; 
            }
            for (int s = 0; s < n2; s++){ 
                temp_seqs[n1 + s][idx_str] = perfil2->sequencias[s][j - 1];
            }
            j--;
        }
        idx_str++;
    }

    Perfil *novo_perfil = malloc(sizeof(Perfil));
    novo_perfil->num_seqs = total_seqs;
    novo_perfil->comprimento = idx_str;
    novo_perfil->sequencias = malloc(sizeof(char*) * total_seqs);

    for (int s = 0; s < total_seqs; s++) {
        novo_perfil->sequencias[s] = malloc(idx_str + 1);
        for (int k = 0; k < idx_str; k++) {
            novo_perfil->sequencias[s][k] = temp_seqs[s][idx_str - 1 - k];
        }
        novo_perfil->sequencias[s][idx_str] = '\0';
        free(temp_seqs[s]);
    }

    free(temp_seqs);
    free(counts1);
    free(counts2);
    free(M_data); free(D_data); free(I_data);
    free(M); free(D); free(I);

    return novo_perfil;
}

Perfil* alinhar_arvore_recur(No *no, char **seq_main, const matriz_sub *m, int flag_proteina) {
    if (no->esq == NULL && no->dir == NULL) {
        Perfil *p = malloc(sizeof(Perfil));
        p->num_seqs = 1;
        p->comprimento = strlen(seq_main[no->seq_idx]);
        p->sequencias = malloc(sizeof(char*));
        size_t len = strlen(seq_main[no->seq_idx]);
        p->sequencias[0] = malloc(len + 1);
        strcpy(p->sequencias[0], seq_main[no->seq_idx]);
        return p;
    }

    Perfil *perfil_esq = alinhar_arvore_recur(no->esq, seq_main, m, flag_proteina);
    Perfil *perfil_dir = alinhar_arvore_recur(no->dir, seq_main, m, flag_proteina);

    printf("Fundindo perfis: (%d seqs) + (%d seqs)...\n", perfil_esq->num_seqs, perfil_dir->num_seqs);
    fflush(stdout);

    Perfil *perfil_fundido = fundir_perfil(perfil_esq, perfil_dir, m, flag_proteina);

    liberar_perfil(perfil_esq);
    liberar_perfil(perfil_dir);

    return perfil_fundido;
}

/* =========================================================================
   EXPORTAÇÃO DE RESULTADOS E AUXILIARES
   ========================================================================= */

void salvar_fasta_alinhado(Perfil *perfil, const char *nome_arquivo) {
    FILE *f = fopen(nome_arquivo, "w");
    if (!f) {
        printf("\nErro ao criar o arquivo FASTA.\n");
        return;
    }
    for (int s = 0; s < perfil->num_seqs; s++) {
        fprintf(f, ">sequencia_%d\n", s + 1);
        fprintf(f, "%s\n", perfil->sequencias[s]);
    }
    fclose(f);
}

void salvar_clustal(Perfil *perfil, const char *nome_arquivo, int tamanho_linha) {
    FILE *f = fopen(nome_arquivo, "w");
    if (!f) {
        printf("\nErro ao criar o arquivo Clustal.\n");
        return;
    }
    int comprimento_total = perfil->comprimento;
    int offset = 0;

    fprintf(f, "CLUSTAL multiple sequence alignment\n\n");

    while (offset < comprimento_total) {
        fprintf(f, "\n");
        for (int s = 0; s < perfil->num_seqs; s++) {
            fprintf(f, "seq_%-2d  ", s + 1);
            for (int k = 0; k < tamanho_linha && (offset + k) < comprimento_total; k++) {
                fputc(perfil->sequencias[s][offset + k], f);
            }
            fprintf(f, "\n");
        }
        fprintf(f, "        ");
        for (int k = 0; k < tamanho_linha && (offset + k) < comprimento_total; k += 10) {
            fprintf(f, ".         ");
        }
        fprintf(f, " %d\n", offset + (tamanho_linha < (comprimento_total - offset) ? tamanho_linha : (comprimento_total - offset)));
        
        offset += tamanho_linha;
    }
    fclose(f);
}

/* =========================================================================
   LIBERAÇÃO DE MEMÓRIA
   ========================================================================= */

void liberar_matriz_distancia(double **matriz_dist, int qtd) {
    if (!matriz_dist) return;
    for (int i = 0; i < qtd; i++) {
        free(matriz_dist[i]);
    }
    free(matriz_dist);
}

void liberar_sequencias(char **seqs, int quantidade) {
    if (!seqs) return;
    for (int i = 0; i < quantidade; i++) {
        free(seqs[i]);
    }
    free(seqs);
}

void liberar_matriz(matriz_sub *m) {
    if (m->valores != NULL) {
        for (int i = 0; i < m->tam; i++) {
            free(m->valores[i]);
        }
        free(m->valores);
        m->valores = NULL;
    }
}

void liberar_perfil(Perfil *perfil) {
    if (!perfil) return;
    for (int i = 0; i < perfil->num_seqs; i++) {
        free(perfil->sequencias[i]);
    }
    free(perfil->sequencias);
    free(perfil);
}

void liberar_arvore(No *n) {
    if (n == NULL) return;
    liberar_arvore(n->esq); 
    liberar_arvore(n->dir); 
    free(n);                
}