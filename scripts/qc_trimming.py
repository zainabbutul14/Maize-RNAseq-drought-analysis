import pyfastx
import glob
import os
import pandas as pd
import subprocess

# --------------------------
# Paths
# --------------------------
fastq_folder = "/Users/zainabbutul/maize_rnaseq_project/ALL_FASTQ"
trimmed_folder = "/Users/zainabbutul/maize_rnaseq_project/TRIMMED"
os.makedirs(trimmed_folder, exist_ok=True)

# Get all FASTQ files
fastq_files = glob.glob(f"{fastq_folder}/*.fastq")

# --------------------------
# Step 1: QC
# --------------------------
qc_data = []
print("🔹 Running QC...\n")
for fq_path in fastq_files:
    fq = pyfastx.Fastq(fq_path)
    total_reads = len(fq)
    read_lengths = [len(read.seq) for read in fq]
    avg_length = sum(read_lengths) / len(read_lengths)

    qc_data.append({
        "sample": os.path.basename(fq_path),
        "total_reads": total_reads,
        "avg_read_length": avg_length
    })

    print(f"📂 {fq_path}")
    print(f"   Total Reads: {total_reads}")
    print(f"   Avg Read Length: {avg_length:.1f}\n")

# Save QC stats to CSV
qc_df = pd.DataFrame(qc_data)
qc_df.to_csv("/Users/zainabbutul/maize_rnaseq_project/QC_stats.csv", index=False)
print("✅ QC stats saved to QC_stats.csv\n")

# --------------------------
# Step 2: Trimming
# --------------------------
print("🔹 Trimming reads...\n")
for fq_path in fastq_files:
    sample_name = os.path.basename(fq_path).replace(".fastq", "_trimmed.fastq")
    output_path = os.path.join(trimmed_folder, sample_name)

    print(f"✂️ Trimming {fq_path} -> {output_path}")

    # Run cutadapt as a subprocess
    cmd = [
        "cutadapt",
        "-q", "20",  # trim low-quality bases
        "-m", "30",  # discard reads < 30 bases
        "-o", output_path,  # output file
        fq_path  # input FASTQ
    ]

    subprocess.run(cmd, check=True)

print("\n✅ Trimming completed. Trimmed files are in TRIMMED folder.")