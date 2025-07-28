from datasets import load_dataset
import json
import os

def get_dataset(max_samples=10000):
    print(f" HuggingFace üzerinden {max_samples} ornek indiriliyo")

    dataset = load_dataset("JYumeko/processed_pubmed_scientific_papers", split="train") #veri setini yükle
    dataset = dataset.select(range(min(max_samples, len(dataset))))#alt küme oluşturma 
    results = []
    for sample in dataset:
        abstract = sample.get("abstract")
        summary = sample.get("summary")
        if abstract and summary:
            results.append({
                "summary": summary.strip(),    
                "abstract": abstract.strip()
            })

    os.makedirs("data", exist_ok=True)
    with open("data/dataset.json", "w", encoding="utf-8") as f: #dataset bu adrestedir -> https://drive.google.com/drive/u/0/folders/1Yp20WxNSLKWAO4WRPA3_HOThAy2m19SZ
        json.dump(results, f, ensure_ascii=False, indent=2)

    print(f"{len(results)} örnek başariyla kaydedildi")

if __name__ == "__main__":
    get_dataset(10000)
