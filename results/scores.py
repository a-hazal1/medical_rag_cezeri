import pandas as pd
import matplotlib.pyplot as plt
from evaluate import load
import logging

logging.getLogger("transformers").setLevel(logging.ERROR)
logging.getLogger("transformers.modeling_utils").setLevel(logging.ERROR)

#değerlendirme metriklerini yükle
bleu = load("bleu")
rouge = load("rouge")
meteor = load("meteor")
bertscore = load("bertscore")

def normalize(text):
    return str(text).strip().lower()

pred_df = pd.read_csv("query_results.csv", on_bad_lines='skip')
ref_df = pd.read_csv("groundtruth.csv", on_bad_lines='skip')

pred_dict = {normalize(q): a for q, a in zip(pred_df["query"], pred_df["answer"])}
ref_dict  = {normalize(q): a for q, a in zip(ref_df["query"], ref_df["answer"])}

#ortak sorguları bul
common_queries = list(set(pred_dict.keys()) & set(ref_dict.keys()))
print(f"Toplam eşleşen soru: {len(common_queries)}")

if not common_queries:
    raise ValueError("Hiç ortak soru yok. CSV dosyalarını kontrol et.")

#referanslar ve tahminler listesi
predictions = [pred_dict[q] for q in common_queries]
references  = [[ref_dict[q]] for q in common_queries]

#metrikleri hesapla
bleu_score    = bleu.compute(predictions=predictions, references=references)
rouge_score   = rouge.compute(predictions=predictions, references=references)
meteor_score  = meteor.compute(predictions=predictions, references=references)
bertscore_val = bertscore.compute(predictions=predictions, references=[r[0] for r in references], lang="en")

metrics = {
    "BLEU":        bleu_score["bleu"],
    "ROUGE1":      rouge_score["rouge1"],
    "ROUGE2":      rouge_score["rouge2"],
    "ROUGEL":      rouge_score["rougeL"],
    "METEOR":      meteor_score["meteor"],
    "BERTScore_F1": sum(bertscore_val["f1"]) / len(bertscore_val["f1"])
}

# Skorları CSV'ye yaz
df_scores = pd.DataFrame({
    "Metric": list(metrics.keys()),
    "Score": list(metrics.values())
})
df_scores.to_csv("metrics.csv", index=False)
print("Skorlar metrics.csv dosyasına kaydedildi.")

#skorları görselleştir
plt.figure(figsize=(10, 6))
bars = plt.bar(df_scores["Metric"], df_scores["Score"], color="lightseagreen", edgecolor="black")
plt.ylim(0, 1.05)
plt.title("Değerlendirme Metrikleri")
plt.ylabel("Skor (0–1)")
plt.xlabel("Metrik")

for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2.0, yval + 0.01, f"{yval:.3f}", ha="center", va="bottom")

plt.tight_layout()
plt.savefig("metrics_plot.png", dpi=300)
plt.show()
print("Grafik metrics_plot.png olarak kaydedildi.")
