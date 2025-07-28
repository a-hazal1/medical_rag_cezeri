from sentence_transformers import SentenceTransformer
import faiss
import numpy as np
import json
import os

#embedding modeli
MODEL = "sentence-transformers/all-MiniLM-L6-v2"
model = SentenceTransformer(MODEL)

with open("data/dataset.json", "r", encoding="utf-8") as f:
    docs = json.load(f)

#embed edilecek içerikleri hazırla
corpus = [doc["summary"] for doc in docs if "summary" in doc and doc["abstract"].strip()]
titles = [doc["abstract"] for doc in docs if "summary" in doc and doc["abstract"].strip()]

#embedding hesapla
embeddings = model.encode(corpus, show_progress_bar=True)

# FAISS index oluştur
dimension = embeddings.shape[1]
index = faiss.IndexFlatL2(dimension)
index.add(np.array(embeddings).astype("float32"))

#kayıtet
os.makedirs("vector_index", exist_ok=True)
faiss.write_index(index, "vector_index/faiss_index.index")
with open("vector_index/docs.json", "w", encoding="utf-8") as f:
    json.dump([
        {"abstract": t, "summary": a} for t, a in zip(titles, corpus)
    ], f, ensure_ascii=False, indent=2)

print("FAISS vektör veritabanı ve belge JSON dosyası oluşturuldu")
