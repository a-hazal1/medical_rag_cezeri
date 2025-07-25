from Bio import Entrez
import json

Entrez.email = "hazallarikann@gmail.com"

def get_pubmed_abstracts(query="diabcanceretes", max_results=50):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    ids = Entrez.read(handle)["IdList"]

    handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="medline", retmode="xml")
    records = Entrez.read(handle)

    results = []
    for rec in records["PubmedArticle"]:
        try:
            title = rec["MedlineCitation"]["Article"]["ArticleTitle"]
            abstract = rec["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0]
            results.append({"title": title, "abstract": abstract})
        except:
            continue

    with open("data/pubmed_docs.json", "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)

get_pubmed_abstracts("heart disease", 100)