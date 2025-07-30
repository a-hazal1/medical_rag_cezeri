import csv
import requests

input_csv = "query.csv"
output_csv = "query_results.csv"

#csv'den sorguları oku
queries = []
with open(input_csv, newline='', encoding='utf-8') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        query = row["query"]
        queries.append(query)


with open(output_csv, mode="w", newline='', encoding="utf-8") as file:
    writer = csv.writer(file)
    writer.writerow(["query", "answer"])  

    for query in queries:
        print(f"Sending: {query}")

        try:
            response = requests.post(
                "http://localhost:8000/query",
                json={"question": query},
                timeout=300
            )
            data = response.json()

            
            answer = data.get("answer", "Cevap bulunamadı")#modelin cevsbı

        except Exception as e:
            print(f"Error for query: {query} -> {e}")
            answer = f"Hata: {e}"

        #csv'ye yaz
        writer.writerow([query, answer])
