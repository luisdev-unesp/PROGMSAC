import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
import time

def descarregar_cenario_proteina(nome_cenario, ficheiro_saida, query_term, slen_range, max_seqs, tamanho_alvo, matriz_arquivo="BLOSUM62.txt"):
    """
    Baixa sequências de proteínas do NCBI E-utilities (db=protein) e formata para o arquivo de entrada do programa em C.
    """
    print(f"\n==========================================")
    print(f"Processando Cenário (Proteínas): {nome_cenario}")
    print(f"Arquivo de saída: {ficheiro_saida}")
    print(f"Buscando {max_seqs} proteínas com tamanho ~{tamanho_alvo} aa...")

    # Busca no banco de dados 'protein' usando a faixa de tamanho em aminoácidos [SLEN]
    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        f"db=protein&term={urllib.parse.quote(query_term)}+AND+{slen_range}[SLEN]"
        f"&retmax={max_seqs * 4}&usehistory=y"
    )

    try:
        resposta_busca = urllib.request.urlopen(search_url).read()
        xml_root = ET.fromstring(resposta_busca)

        webenv_element = xml_root.find('WebEnv')
        query_key_element = xml_root.find('QueryKey')

        if webenv_element is None or query_key_element is None:
            raise Exception("Erro ao recuperar WebEnv ou QueryKey do NCBI.")

        webenv = webenv_element.text
        query_key = query_key_element.text

        fetch_url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
            f"db=protein&query_key={query_key}&WebEnv={webenv}&rettype=fasta&retmode=text"
        )
        
        print("Baixando dados FASTA de proteínas do NCBI...")
        dados_fasta = urllib.request.urlopen(fetch_url).read().decode('utf-8', errors='ignore')

        sequencias = []
        seq_atual = []

        # Aminoácidos padrão aceitos pela BLOSUM62
        aminoacidos_validos = set("ARNDCQEGHILKMFPSTWYV")

        for linha in dados_fasta.split('\n'):
            linha = linha.strip()
            if linha.startswith('>'):
                if seq_atual:
                    seq_completa = ''.join(seq_atual).upper()
                    # Filtra caracteres ambíguos/desconhecidos (ex: 'X', 'U', 'Z', 'B', '*')
                    if all(c in aminoacidos_validos for c in seq_completa) and len(seq_completa) >= int(tamanho_alvo * 0.5):
                        sequencias.append(seq_completa)
                    seq_atual = []
            elif linha:
                seq_atual.append(linha)
                
        if seq_atual:
            seq_completa = ''.join(seq_atual).upper()
            if all(c in aminoacidos_validos for c in seq_completa) and len(seq_completa) >= int(tamanho_alvo * 0.5):
                sequencias.append(seq_completa)

        sequencias = sequencias[:max_seqs]

        if len(sequencias) == 0:
            print(f"Aviso: Nenhuma sequência válida encontrada para {nome_cenario}.")
            return

        # Preenche com replicação caso o número retornado seja menor que o solicitado
        while len(sequencias) < max_seqs:
            sequencias.append(sequencias[0])

        # Escreve o arquivo no formato de entrada aceito pelo programa C (main.c)
        with open(ficheiro_saida, "w") as f:
            f.write("p\n")               # Opção 'p' para proteína
            f.write(f"{matriz_arquivo}\n")# Nome do arquivo da matriz BLOSUM62
            f.write(f"{len(sequencias)}\n") # Quantidade total de sequências
            
            for s in sequencias:
                # Ajusta a sequência ao tamanho alvo (preenche com 'A' - Alanina - se for menor)
                s_formatada = s[:tamanho_alvo].ljust(tamanho_alvo, 'A')
                f.write(s_formatada + "\n")

        print(f"Sucesso! Arquivo '{ficheiro_saida}' gerado com {len(sequencias)} proteínas de ~{tamanho_alvo} aa.")

    except Exception as e:
        print(f"Erro ao processar {nome_cenario}: {e}")

# ==============================================================================
# CONFIGURAÇÃO DOS 4 CENÁRIOS DE TESTE COM PROTEÍNAS
# ==============================================================================

cenarios = [
    {
        "nome": "1. Poucas Sequências Pequenas (Poucas Qtd, Pequeno Tam)",
        "ficheiro": "teste1_prot_poucas_pequenas.txt",
        "query": "bacteria[filter]",
        "slen": "50:150",
        "qtd": 5,          # Poucas sequências
        "tamanho": 100     # Sequência pequena (100 aminoácidos)
    },
    {
        "nome": "2. Poucas Sequências Grandes (Poucas Qtd, Grande Tam)",
        "ficheiro": "teste2_prot_poucas_grandes.txt",
        "query": "bacteria[filter]",
        "slen": "1000:3000",
        "qtd": 5,          # Poucas sequências
        "tamanho": 1500    # Sequência grande (1500 aminoácidos)
    },
    {
        "nome": "3. Muitas Sequências Pequenas (Muitas Qtd, Pequeno Tam)",
        "ficheiro": "teste3_prot_muitas_pequenas.txt",
        "query": "bacteria[filter]",
        "slen": "50:150",
        "qtd": 150,        # Muitas sequências
        "tamanho": 100     # Sequência pequena
    },
    {
        "nome": "4. Muitas Sequências Grandes (Muitas Qtd, Grande Tam)",
        "ficheiro": "teste4_prot_muitas_grandes.txt",
        "query": "bacteria[filter]",
        "slen": "1000:3000",
        "qtd": 40,         # Muitas sequências
        "tamanho": 1500    # Sequência grande
    }
]

# Execução em lote dos cenários
for c in cenarios:
    descarregar_cenario_proteina(
        nome_cenario=c["nome"],
        ficheiro_saida=c["ficheiro"],
        query_term=c["query"],
        slen_range=c["slen"],
        max_seqs=c["qtd"],
        tamanho_alvo=c["tamanho"],
        matriz_arquivo="BLOSUM62.txt"
    )
    time.sleep(2)
