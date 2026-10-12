import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
import time

def descarregar_cenario(nome_cenario, ficheiro_saida, query_term, slen_range, max_seqs, tamanho_alvo):
    """
    Baixa sequências do NCBI E-utilities e formata para o arquivo de teste.
    """
    print(f"\n==========================================")
    print(f"Processando Cenário: {nome_cenario}")
    print(f"Arquivo de saída: {ficheiro_saida}")
    print(f"Buscando {max_seqs} sequências com tamanho ~{tamanho_alvo} pb...")

    # Aumentamos retmax para ter margem de filtragem (sem Ns)
    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        f"db=nuccore&term={urllib.parse.quote(query_term)}+AND+{slen_range}[SLEN]"
        f"&retmax={max_seqs * 3}&usehistory=y"
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
            f"db=nuccore&query_key={query_key}&WebEnv={webenv}&rettype=fasta&retmode=text"
        )
        
        print("Baixando dados FASTA do NCBI...")
        dados_fasta = urllib.request.urlopen(fetch_url).read().decode('utf-8', errors='ignore')

        sequencias = []
        seq_atual = []

        for linha in dados_fasta.split('\n'):
            linha = linha.strip()
            if linha.startswith('>'):
                if seq_atual:
                    seq_completa = ''.join(seq_atual).upper()
                    # Filtra ruídos ('N') e garante comprimento mínimo razoável
                    if 'N' not in seq_completa and len(seq_completa) >= int(tamanho_alvo * 0.5):
                        sequencias.append(seq_completa)
                    seq_atual = []
            elif linha:
                seq_atual.append(linha)
                
        if seq_atual:
            seq_completa = ''.join(seq_atual).upper()
            if 'N' not in seq_completa and len(seq_completa) >= int(tamanho_alvo * 0.5):
                sequencias.append(seq_completa)

        sequencias = sequencias[:max_seqs]

        if len(sequencias) == 0:
            print(f"Aviso: Nenhuma sequência válida encontrada para {nome_cenario}.")
            return

        # Preenche com replicação se a busca retornar menos que o solicitado
        while len(sequencias) < max_seqs:
            sequencias.append(sequencias[0])

        # Escreve no formato esperado pelo código C
        with open(ficheiro_saida, "w") as f:
            f.write("d\n")  # Tipo de dado (DNA)
            f.write(f"{len(sequencias)}\n")  # Número total de sequências
            
            for s in sequencias:
                s_formatada = s[:tamanho_alvo].ljust(tamanho_alvo, 'A')
                f.write(s_formatada + "\n")

        print(f"Sucesso! Arquivo '{ficheiro_saida}' gerado com {len(sequencias)} sequências de {tamanho_alvo} pb.")

    except Exception as e:
        print(f"Erro ao processar {nome_cenario}: {e}")

# ==============================================================================
# CONFIGURAÇÃO DOS 4 CENÁRIOS DE TESTE
# ==============================================================================

cenarios = [
    {
        "nome": "1. Poucas Sequências Pequenas",
        "ficheiro": "teste1_poucas_pequenas.txt",
        "query": "bacteria[filter]",
        "slen": "100:300",
        "qtd": 5,
        "tamanho": 200
    },
    {
        "nome": "2. Poucas Sequências Ultralongas (~1 Milhão pb)",
        "ficheiro": "teste2_poucas_ultralongas.txt",
        "query": "\"complete genome\"[Title] AND bacteria[filter]",
        "slen": "900000:1200000",
        "qtd": 5,  # Alterado para 5 sequências
        "tamanho": 1000000
    },
    {
        "nome": "3. Várias Sequências Pequenas",
        "ficheiro": "teste3_varias_pequenas.txt",
        "query": "bacteria[filter]",
        "slen": "100:300",
        "qtd": 200,
        "tamanho": 200
    },
    {
        "nome": "4. Várias Sequências Grandes",
        "ficheiro": "teste4_varias_grandes.txt",
        "query": "\"complete sequence\"[Title] AND bacteria[filter]",
        "slen": "15000:30000",
        "qtd": 50,
        "tamanho": 20000
    }
]

# Execução em lote com pausa para não exceder limites de requisições da API
for c in cenarios:
    descarregar_cenario(
        nome_cenario=c["nome"],
        ficheiro_saida=c["ficheiro"],
        query_term=c["query"],
        slen_range=c["slen"],
        max_seqs=c["qtd"],
        tamanho_alvo=c["tamanho"]
    )
    time.sleep(2)
