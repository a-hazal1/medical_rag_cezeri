from fastapi import FastAPI #fastapi kütüphanesi
from pydantic import BaseModel #sorgu modeli için pydantic
from sentence_transformers import SentenceTransformer #embedding için sentence-transformers
import faiss #faiss kütüphanesi
import json 
import numpy as np 
import requests 
import time 
from fastapi.middleware.cors import CORSMiddleware 
from fastapi.staticfiles import StaticFiles #frontend dosyaları için

app = FastAPI()

# CORS ayarları
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# FAISS index ve belgeler (sadece arama için, çıktıdan kaldırıyoruz)
index = faiss.read_index("vector_index/faiss_index.index")
with open("vector_index/docs.json", "r", encoding="utf-8") as f:
    documents = json.load(f)

#embedding modeli
embedder = SentenceTransformer("sentence-transformers/all-MiniLM-L6-v2")

#ollama ayarları
OLLAMA_URL = "http://host.docker.internal:11434/api/generate"
OLLAMA_MODEL = "mistral"

class Query(BaseModel):#sorgu modeli
    question: str

#sorgu uzunluk sınıflandırması
def get_query_length_label(query: str) -> str:
    wc = len(query.split()) # kelime sayısı
    if wc <= 5:
        return "Short"
    elif wc <= 12:
        return "Medium"
    else:
        return "Long"

@app.post("/query")
def answer_question(query: Query): #sorgu işleme fonksiyonu
    start_time = time.time()

    #Embed + FAISS
    t0 = time.time()
    vec = embedder.encode([query.question]).astype("float32")
    distances, indices = index.search(vec, 3)
    retrieval_time = round((time.time() - t0) * 1000, 2)

    #context oluştur 
    context = "\n".join(documents[i].get("summary", "") for i in indices[0])

    #prompt
    prompt = f"""
    Context:
    {context}

    Question: {query.question}
    Answer:"""

    #llm çağrısı
    t1 = time.time()
    resp = requests.post(OLLAMA_URL, json={
        "model": OLLAMA_MODEL,
        "prompt": prompt,
        "stream": False,
    })
    t2 = time.time()
    generation_time = round((t2 - t1) * 1000, 2)
    total_time = round((time.time() - start_time) * 1000, 2)

    #json cevabı
    return {
        "query":             query.question,
        "query_length":      len(query.question.split()),
        "query_length_label": get_query_length_label(query.question),
        "retrieval_time_ms": retrieval_time,
        "generation_time_ms": generation_time,
        "total_time_ms":     total_time,
        "answer":            resp.json().get("response", "__no_response__").strip()
    }


app.mount("/", StaticFiles(directory="static", html=True), name="static") #frontend dosyaları
