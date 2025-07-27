from Bio import Entrez
import json
import xml.etree.ElementTree as ET

Entrez.email = "hazallarikann@gmail.com"

def get_pubmed_abstracts(query="diabetes", max_results=50):
    print(f" '{query}' için PubMed araması yapılıyor...")

    # 1. E-Search
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    search_results = Entrez.read(handle)
    ids = search_results["IdList"]

    if not ids:
        print(" Hiç sonuç bulunamadı!")
        return

    print(f" {len(ids)} ID bulundu, detaylar çekiliyor...")

    # 2. E-Fetch (XML olarak al, sonra elle parse edeceğiz)
    handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="medline", retmode="xml")
    xml_data = handle.read()
    root = ET.fromstring(xml_data)

    # 3. XML'den başlık ve özet çıkar
    ns = {'ns': 'http://www.ncbi.nlm.nih.gov/pubmed'}
    results = []
    for article in root.findall(".//PubmedArticle"):
        try:
            title = article.findtext(".//ArticleTitle")
            abstract = " ".join([elem.text for elem in article.findall(".//AbstractText") if elem.text]) or None

            results.append({
                "title": title,
                "abstract": abstract
            })
        except Exception as e:
            print(f"⚠️ Hata: {e}")

    # 4. JSON'a kaydet
    with open("data/pubmed_docs.json", "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)

    print(f"{len(results)} makale başarıyla kaydedildi → data/pubmed_docs.json")

if __name__ == "__main__":
    get_pubmed_abstracts("diabetes", 100)
