from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from rag_engine import ask_medical_question

app = FastAPI()

class Question(BaseModel):
    query: str

@app.post("/ask")
async def ask(question: Question):
    try:
        answer = ask_medical_question(question.query)
        return {"answer": answer}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
