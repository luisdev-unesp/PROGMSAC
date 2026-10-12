#ifndef PROGMSAC_H
#define PROGMSAC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <omp.h>

#define match 2
#define mismatch -3
#define gapopen -10
#define gapextend -2

/* =========================================================================
   ESTRUTURAS (STRUCTS)
   ========================================================================= */

typedef struct {
    int tam;
    int** valores;
    int hash[256];
} matriz_sub;

typedef struct No {
    int id;               
    int seq_idx;          
    double dist_pai;      
    struct No *esq;       
    struct No *dir;       
} No;

typedef struct {
    int num_seqs;      
    int comprimento;   
    char **sequencias; 
} Perfil;

/* =========================================================================
   DECLARAÇÕES DAS FUNÇÕES (PROTÓTIPOS)
   ========================================================================= */

// Leitura e gerenciamento de sequências
char** ler_sequencias(int quantidade, int tam_max);

// Operações com Matriz de Substituição e Score
int carregar_matriz(matriz_sub* m);
int indice_hash(char amino, matriz_sub* m);
int obter_score(char a, char b, int flag_proteina, const matriz_sub* m);

// Alinhamento Par a Par e Matriz de Distância
double** calcular_matriz_distancia(char **seqs, int qtd, int flag_proteina, const matriz_sub *m);
int calcular_score_gotoh_1d(const char *seq1, const char *seq2, int flag_proteina, const matriz_sub *m);

// Algoritmo Neighbor-Joining (Construção de Árvore)
double* neighbor_step1(double **matriz, int n, int *ativo);
void neighbor_step2(double *vetor, int n, double **matriz, int *ativo, int *idx_matriz_u, int *idx_matriz_v);
void neighbor_step3(No **nos, double **matriz, double *vetor, int i, int j, int novo_id);
No* arvore_neighbor(double **dist_inicial, int n);

// Estrutura de Nós da Árvore
No* criar_folha(int id, int seq_idx);
No* criar_no_interno(int id, No *filho_esq, No *filho_dir, double dist_esq, double dist_dir);

// Alinhamento Múltiplo de Sequências (Perfis)
Perfil* fundir_perfil(Perfil *perfil1, Perfil *perfil2, const matriz_sub *m, int flag_proteina);
Perfil* alinhar_arvore_recur(No *no, char **seq_main, const matriz_sub *m, int flag_proteina);
double calculo_score_coluna(Perfil *p1, Perfil *p2, int col1, int col2, const matriz_sub *m, int flag_proteina);

// Exportação de Resultados
void salvar_fasta_alinhado(Perfil *perfil, const char *nome_arquivo);
void salvar_clustal(Perfil *perfil, const char *nome_arquivo, int tamanho_linha);

// Liberação de Memória
void liberar_matriz_distancia(double **matriz_dist, int qtd);
void liberar_sequencias(char **seqs, int quantidade);
void liberar_matriz(matriz_sub *m);
void liberar_arvore(No *n);
void liberar_perfil(Perfil *perfil);

#endif // PROGMSAC_H