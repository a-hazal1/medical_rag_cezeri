Medical RAG: Bilimsel Soru-Cevap Sistemi
Bu proje, tıbbi sorulara bilimsel cevaplar üreten bir Retrieval-Augmented Generation (RAG) sistemidir. HuggingFace üzerinden indirilen PubMed özet verisini işler, vektör veritabanı oluşturur ve bir büyük dil modeli (LLM) ile doğru, açıklayıcı yanıtlar üretir.

Sistem, veri hazırlama işlemleri Python ortamında, API ve değerlendirme işlemleri ise Docker içinde gerçekleştirilir.

Özellikler
1. HuggingFace PubMed özet verisiyle çalışma
2. SentenceTransformer ile embedding işlemi
3. FAISS ile hızlı vektör arama
4. Mistral (Ollama) üzerinden LLM kullanımı
5. FastAPI 
6. BLEU, ROUGE, METEOR, BERTScore ile değerlendirme
7. Docker ile kolay servis kurulumu

Proje Yapısı
.
├── Dockerfile                   # API sunucusu için yapılandırma
├── docker-compose.yml           # Uygulama orkestrasyonu
├── data_loader.py               # PubMed verisi indirme (sadece Python)
├── ingest.py                    # FAISS dizini oluşturma (sadece Python)
├── rag_api.py                   # RAG sistem REST API
├── results/
    └── time_results.py              # Sorgu sürelerini ölçer
    └── query_results.py             # Model çıktısını kaydeder
    └── scores.py                    # Otomatik değerlendirme metrikleri
    └── query.csv                    # Test sorguları
    └── groundtruth.csv              # Doğru cevaplar
├── requirements.txt             # Python bağımlılıkları
├── static/
    └── index.html               #frontend bağlantısı
├── data/                        # İndirilen veri dosyası
    └── vector_index/                # FAISS index ve belge verisi

Kurulum Adımları
Gerekli Yazılımlar
1. Python ≥ 3.8
2. Docker
3. Ollama
    ollama run mistral

Veri Seti ve FAISS Dizini (Python Ortamında)
Aşağıdaki iki adım sadece ilk çalıştırmada yapılmalıdır:

1. Veri Setini İndir
    python data_loader.py
Bu işlem HuggingFace üzerinden PubMed özet verisini indirir ve data/dataset.json dosyasına kaydeder.

2. Vektör Dizini Oluştur
    python ingest.py
Bu komut:
all-MiniLM-L6-v2 modeliyle özetleri gömlemeye çevirir.
FAISS kullanarak vector_index/faiss_index.index dosyasını oluşturur.
İlgili belgeleri vector_index/docs.json içinde saklar.

Docker ile API Servisini Başlat
API'yi çalıştırmak için:
docker-compose up --build
localhost:8000 üzerinde FastAPI çalışır.
Ollama üzerinden mistral LLM modelini kullanır.

🔎 API Kullanımı
POST /query
Örnek İstek:
{
  "question": "What causes high blood pressure?"
}
Örnek Yanıt:
{
  "query": "What causes high blood pressure?",
  "query_length": 5,
  "query_length_label": "Short",
  "retrieval_time_ms": 43.8,
  "generation_time_ms": 691.1,
  "total_time_ms": 734.9,
  "answer": "High blood pressure can be caused by genetics, diet, lack of exercise, and stress."
}

Değerlendirme
1. Zaman Ölçümü
    python time_results.py
2. Model Yanıtlarını Üret
    python query_results.py
3. Otomatik Metriklerle Skorla
    python scores.py

Kullanılan Teknolojiler
  Amaç	            Teknoloji / Model
Embedding	    sentence-transformers/all-MiniLM-L6-v2
Vektör Arama	FAISS
Yanıt Üretimi	Mistral (Ollama)
API	            FastAPI
Değerlendirme	evaluate (BLEU, ROUGE, METEOR, BERTScore)

Uyarı
Bu sistem klinik veya medikal tanı amaçlı kullanılamaz. Sadece akademik ve araştırma kullanımına yöneliktir.

Katkı ve İletişim
Katkı sağlamak veya öneri göndermek isterseniz; geliştirici: [hazallarikann@gmail.com]

