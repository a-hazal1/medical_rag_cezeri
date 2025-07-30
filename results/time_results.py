import csv
import time
import requests

input_csv = "query.csv"
output_csv = "timing_results.csv"

#sorgu uzunluk sınıflandırması
def get_query_length_label(query: str) -> str:
    wc = len(query.split()) # kelime sayısı
    if wc <= 5:
        return "Short"
    elif wc <= 12:
        return "Medium"
    else:
        return "Long"

# CSV'den sorguları oku
queries = []
with open(input_csv, newline='', encoding='utf-8') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        query = row["query"]
        q_length_label = get_query_length_label(query)
        queries.append((query, q_length_label))

# Sonuçları yaz
with open(output_csv, mode="w", newline='', encoding="utf-8") as file:
    writer = csv.writer(file)
    writer.writerow(["Query", "Query_Length", "Retrieval_Time_MS", "Generation_Time_MS", "Total_Time_MS"])

    for query, q_length in queries:
        print(f"Sending: {query}")
        start_time = time.time()

        try:
            response = requests.post(
                "http://localhost:8000/query",
                json={"question": query},
                timeout=300
            )
            total_time = (time.time() - start_time) * 1000  # toplam süre (ms)
            data = response.json()

            retrieval_time = data.get("retrieval_time_ms", -1)
            generation_time = data.get("generation_time_ms", -1)

        except Exception as e:
            print(f"Error for query: {query} -> {e}")
            retrieval_time = generation_time = total_time = -1

        writer.writerow([query, q_length, int(retrieval_time), int(generation_time), int(total_time)])
