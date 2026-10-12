#include <stdio.h>
#include <stdlib.h>
#include "PROGMSAC.h"

int main() {
    printf("Vai ser feito o alinhamento de Proteinas (p) ou DNA/RNA (d): ");
    char escolha;
    if (scanf(" %c", &escolha) != 1) return 1;

    int flag_proteina = (escolha == 'p');
    matriz_sub m = {0, NULL, {0}};

    if (flag_proteina && !carregar_matriz(&m)) {
        return 1;
    }

    int quantidade;
    printf("Quantas sequencias vão ser verificadas: ");
    if (scanf("%d", &quantidade) != 1 || quantidade < 1) return 1;

    if (quantidade < 3) {
        printf("\nSão necessárias pelo menos 3 sequencias para construir a arvore.\n");
        if (flag_proteina) liberar_matriz(&m);
        return 0;
    }

    // 1. Leitura das sequências
    char **seq_main = ler_sequencias(quantidade, 1000010);
    if (!seq_main) {
        if (flag_proteina) liberar_matriz(&m);
        return 1;
    }

    // 2. Cálculo da Matriz de Distância
    double **matriz_dist = calcular_matriz_distancia(seq_main, quantidade, flag_proteina, &m);

    // 3. Construção e Alinhamento Múltiplo via Árvore
    No *raiz_filogenetica = arvore_neighbor(matriz_dist, quantidade);
    printf("\nArvore filogenetica calculada com sucesso!\n");

    Perfil *perfil_final = alinhar_arvore_recur(raiz_filogenetica, seq_main, &m, flag_proteina);

    // 4. Exportação dos Resultados
    salvar_fasta_alinhado(perfil_final, "alinhamento.fasta");
    salvar_clustal(perfil_final, "alinhamento.aln", 60);

    printf("\nSucesso! Os seguintes arquivos foram criados:\n");
    printf("1. alinhamento.fasta\n");
    printf("2. alinhamento.aln (Clustal)\n");

    // 5. Liberação de Memória
    liberar_perfil(perfil_final);
    liberar_arvore(raiz_filogenetica);
    liberar_matriz_distancia(matriz_dist, quantidade);
    liberar_sequencias(seq_main, quantidade);

    if (flag_proteina) {
        liberar_matriz(&m);
    }

    return 0;
}