from fastapi import FastAPI
from pydantic import BaseModel
from sentence_transformers import SentenceTransformer
from transformers import AutoTokenizer, AutoModelForCausalLM
import torch
import faiss
import json
import numpy as np

app = FastAPI()

# FAISS index'i ve belgeleri yükle
index = faiss.read_index("vector_index/faiss_index.index")
with open("vector_index/docs.json", "r", encoding="utf-8") as f:
    documents = json.load(f)

# Embedding modeli
embedder = SentenceTransformer("sentence-transformers/all-MiniLM-L6-v2")

# BioGPT modeli
BIOGPT_MODEL = "microsoft/BioGPT"
tokenizer = AutoTokenizer.from_pretrained(BIOGPT_MODEL)
model = AutoModelForCausalLM.from_pretrained(BIOGPT_MODEL)

class Query(BaseModel):
    question: str

@app.post("/query")
def answer_question(query: Query):
    # Sorguyu embed et
    vec = embedder.encode([query.question])
    vec = np.array(vec).astype("float32")

    # En yakın 3 belgeyi bul
    distances, indices = index.search(vec, 3)

    # Abstract'ları birleştir
    context = "\n".join([documents[i].get("abstract", "") for i in indices[0]])

    # BioGPT'ye prompt oluştur
    prompt = f"Context:\n{context}\n\nQuestion: {query.question}\nAnswer:"
    input_ids = tokenizer(prompt, return_tensors="pt").input_ids
    output = model.generate(input_ids, max_new_tokens=150)
    answer = tokenizer.decode(output[0], skip_special_tokens=True)

    return {
        "answer": answer, # Generated answer from BioGPT
        "sources": [documents[i] for i in indices[0]] # Include the sources of the answer
    }
