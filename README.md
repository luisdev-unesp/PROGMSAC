
# PROGMSAC

Alinhamento múltiplo de sequências (MSA) **progressivo**, escrito em C, para **DNA/RNA** e **proteínas**.
Usa alinhamento par a par com gap afim (Gotoh), árvore guia por Neighbor-Joining e fusão progressiva de perfis, com a etapa mais pesada paralelizada com OpenMP.

## Como funciona

1. **Matriz de distância** — todos os pares de sequências são alinhados com o algoritmo de Gotoh (gap afim, memória linear). A distância de cada par é o maior score de auto-alinhamento entre as duas sequências menos o score do par (nunca negativa). Essa etapa roda em paralelo (OpenMP).
2. **Árvore guia** — a matriz de distância alimenta o algoritmo **Neighbor-Joining**.
3. **Alinhamento progressivo** — a árvore é percorrida das folhas até a raiz, fundindo perfis (perfil × perfil) com Gotoh. O score de uma coluna é a soma dos scores par a par entre os resíduos das duas colunas (gap × gap = 0, gap × resíduo = penalidade de extensão).
4. **Saída** — o alinhamento final é gravado em FASTA e em Clustal.

### Pontuação

| Tipo | Substituição |
|------|--------------|
| DNA/RNA | match `+2`, mismatch `-3` |
| Proteína | matriz **BLOSUM62** (`BLOSUM62.txt`) |

Penalidades de gap: abertura `-10`, extensão `-2`. Todos os valores ficam em `src/PROGMSAC.h` (`match`, `mismatch`, `gapopen`, `gapextend`).

## Estrutura do repositório

```
.
├── src/
│   ├── main.c          # entrada/saída e fluxo principal
│   ├── PROGMSAC.c      # Gotoh, Neighbor-Joining, fusão de perfis, exportação
│   └── PROGMSAC.h      # structs, protótipos e parâmetros de pontuação
├── scripts/
│   ├── baixar_dna.py   # gera entradas de teste de DNA a partir do NCBI
│   └── baixar_prot.py  # gera entradas de teste de proteínas a partir do NCBI
├── BLOSUM62.txt        # matriz de substituição (formato NCBI)
└── README.md
```

## Requisitos

- `gcc` com suporte a OpenMP (Linux, WSL, ou macOS com `brew install gcc`)
- Python 3 (apenas biblioteca padrão) e acesso à internet, somente para os scripts de download
- Opcional: [T-Coffee](https://tcoffee.org/) para visualizar o alinhamento em HTML

## Compilação

```bash
gcc -O3 -march=native -fopenmp src/main.c src/PROGMSAC.c -o progmsa -lm
```

## Formato de entrada

O programa lê da entrada padrão, então pode ser usado de forma interativa ou com redirecionamento (`< arquivo`).
As sequências são **cruas** (sem cabeçalho FASTA), uma por linha, sem espaços. São necessárias **pelo menos 3** sequências, com até 1.000.010 caracteres cada.

**DNA/RNA**

```
d
3
ACGTACGTAC
ACGTTCGTAC
ACGAACGTAC
```

**Proteínas** (a segunda linha é o caminho da matriz BLOSUM62)

```
p
BLOSUM62.txt
3
MKVLAAGIVG
MKVLSAGIVG
MKILAAGLVG
```

Qualquer caractere diferente de `p` na primeira linha é tratado como DNA/RNA.

## Gerando dados de teste

Os scripts baixam sequências do NCBI (E-utilities) e já gravam os arquivos no formato de entrada acima. Rode a partir da raiz do repositório:

```bash
python3 scripts/baixar_dna.py
python3 scripts/baixar_prot.py
```

| Cenário | DNA (arquivo / qtd × tamanho) | Proteína (arquivo / qtd × tamanho) |
|---------|-------------------------------|------------------------------------|
| 1. Poucas e pequenas | `teste1_poucas_pequenas.txt` — 5 × 200 pb | `teste1_prot_poucas_pequenas.txt` — 5 × 100 aa |
| 2. Poucas e grandes | `teste2_poucas_ultralongas.txt` — 5 × 1.000.000 pb | `teste2_prot_poucas_grandes.txt` — 5 × 1500 aa |
| 3. Muitas e pequenas | `teste3_varias_pequenas.txt` — 200 × 200 pb | `teste3_prot_muitas_pequenas.txt` — 150 × 100 aa |
| 4. Muitas e grandes | `teste4_varias_grandes.txt` — 50 × 20.000 pb | `teste4_prot_muitas_grandes.txt` — 40 × 1500 aa |

Se a busca retornar menos sequências do que o pedido, o script repete a primeira; sequências curtas são completadas com `A`.

## Executando

```bash
./progmsa < teste1_poucas_pequenas.txt
```

Arquivos gerados no diretório atual:

- `alinhamento.fasta` — alinhamento em FASTA
- `alinhamento.aln` — alinhamento em formato Clustal (60 colunas por bloco)

Para visualizar com o T-Coffee:

```bash
t_coffee -infile alinhamento.fasta -output html
```

## Limitações conhecidas

- **Memória na fusão de perfis:** cada fusão aloca três matrizes de `double` de tamanho `(L1+1) × (L2+1)`, ou seja, memória O(L1·L2). Para duas sequências de 20.000 posições isso já passa de 9 GB, e o cenário de ~1 milhão de pb (`teste2_poucas_ultralongas.txt`) não cabe em memória.
- **Tempo:** o alinhamento par a par é O(L²) por par de sequências, e a matriz de distância tem n(n−1)/2 pares.
- **Ordem na saída:** as sequências no FASTA/Clustal aparecem na ordem das folhas da árvore guia, e não na ordem de entrada; os rótulos `sequencia_N` / `seq_N` seguem essa ordem de saída.
- A primeira linha da entrada define o tipo (`p` ou outro caractere para DNA/RNA); não há validação dos caracteres das sequências.
- **Mudança em sequẽncias de RNA:** no alinhamento de sequẽncias de proteínas é preciso ajustar os valores de match, mismatch, gap e gapextends.
