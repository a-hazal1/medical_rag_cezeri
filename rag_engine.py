from llama_index.core import SimpleDirectoryReader, Settings, VectorStoreIndex, PromptTemplate
from llama_index.embeddings.huggingface import HuggingFaceEmbedding
from llama_index.llms.ollama import Ollama

# Veri yükleme
loader = SimpleDirectoryReader(
    input_dir="database",
    required_exts=[".pdf"],
    recursive=True
)
docs = loader.load_data()

# Gömme işlemi
embed_model = HuggingFaceEmbedding(
    model_name="sentence-transformers/all-MiniLM-L6-v2",# HuggingFace model adı
    trust_remote_code=True
)
Settings.embed_model = embed_model

# Index oluşturuluyor
index = VectorStoreIndex.from_documents(docs)

# LLM kurulumu
llm = Ollama(model="llama3", request_timeout=120.0)
Settings.llm = llm

# Query Engine
query_engine = index.as_query_engine(streaming=True, similarity_top_k=4)

# Prompt Template
qa_prompt_tmpl_str = (
    "You are a knowledgeable and reliable medical assistant. "
    "Based on the medical context provided below, please give an accurate and clear answer.\n"
    "---------------------\n"
    "{context_str}\n"
    "---------------------\n"
    "If the information is insufficient to answer, say 'I don't know!'.\n"
    "Query: {query_str}\n"
    "Answer: "
)

qa_prompt_tmpl = PromptTemplate(template=qa_prompt_tmpl_str)
query_engine.update_prompts({"response_synthesizer:text_qa_template": qa_prompt_tmpl})


# Ana sorgulama fonksiyonu
def ask_medical_question(query: str) -> str:
    response = query_engine.query(query)
    return str(response)
