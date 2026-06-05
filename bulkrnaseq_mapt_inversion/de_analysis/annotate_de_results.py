import csv
from pathlib import Path

# ===============================================================
# CONFIG
# ===============================================================
# Gene-id -> (symbol, type) mapping produced by make_gene_id_mapping.py.
mapping_file_path = "gene_id_to_symbol_and_type_mapping.csv"

# DESeq2 result CSVs to annotate. For each file `X.csv` this writes
# `X_with_name_and_type.csv` next to it (the form the plotting scripts read).
# Add one entry per differential-expression result you want to annotate.
input_files = [
    "path/to/all_genes_result_cov_rin_sex_ancestry_flipped_eur_vs_aa.csv",
]
# ===============================================================

gene_mapping = {}
with open(mapping_file_path, 'r') as infile:
    reader = csv.reader(infile)
    next(reader)
    for rows in reader:
        gene_mapping[rows[0]] = (rows[1], rows[2])

def process_file(input_path):
    input_path = Path(input_path)
    output_path = input_path.parent / f"{input_path.stem}_with_name_and_type{input_path.suffix}"
    
    with open(input_path, 'r') as infile, open(output_path, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        
        header = next(reader)
        new_header = ['gene_id', 'Gene_Symbol', 'Gene_Type'] + header[:-1]
        writer.writerow(new_header)
        
        for row in reader:
            gene_id = row[-1]
            gene_symbol, gene_type = gene_mapping.get(gene_id, ('Not found', 'Not found'))
            new_row = [gene_id, gene_symbol, gene_type] + row[:-1]
            writer.writerow(new_row)

# Process each file
for file_path in input_files:
    print(f"Processing: {Path(file_path).name}")
    process_file(file_path)
    print(f"Completed: {Path(file_path).name}")