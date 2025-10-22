import pandas as pd
import glob, os

quant_dirs = glob.glob("/Users/zainabbutul/maize_rnaseq_project/QUANT/*")
count_tables = []

for qd in quant_dirs:
    sample_name = os.path.basename(qd)
    df = pd.read_csv(f"{qd}/abundance.tsv", sep="\t")
    df = df[["target_id", "est_counts"]].rename(columns={"est_counts": sample_name})
    count_tables.append(df.set_index("target_id"))

count_matrix = pd.concat(count_tables, axis=1)
count_matrix.to_csv("/Users/zainabbutul/maize_rnaseq_project/count_matrix.csv")

print("✅ Count matrix saved at /Users/zainabbutul/maize_rnaseq_project/count_matrix.csv")