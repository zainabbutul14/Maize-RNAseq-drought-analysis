import os
import pandas as pd

# Path to your GTF file
gtf_file = "/Users/zainabbutul/Downloads/ncbi_dataset 2/ncbi_dataset/data/GCF_902167145.1/genomic.gtf"

# Read GTF manually (tab-separated, skip comments)
gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=[
    "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"
])

# Filter for transcript features
transcripts = gtf[gtf["feature"] == "transcript"].copy()

# Extract transcript_id and gene_id from the attribute column
transcripts["transcript_id"] = transcripts["attribute"].str.extract('transcript_id "([^"]+)"')
transcripts["gene_id"] = transcripts["attribute"].str.extract('gene_id "([^"]+)"')

# Keep only transcript_id and gene_id columns
tx2gene = transcripts[["transcript_id", "gene_id"]]

# Save as CSV
tx2gene.to_csv("/Users/zainabbutul/maize_rnaseq_project/reference/transcript2gene.csv", index=False)

print("✅ transcript2gene.csv created successfully!")



import shutil
import gzip

# Paths
downloaded_file = "/Users/zainabbutul/Downloads/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.cdna.fa.gz"
reference_folder = "/Users/zainabbutul/maize_rnaseq_project/reference"
final_fasta = os.path.join(reference_folder, "transcripts.fa")


print(f"✅ FASTA ready at: {final_fasta}")
