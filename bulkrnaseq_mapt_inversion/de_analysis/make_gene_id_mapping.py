import csv

# Check how many lines we are reading and parsing
line_count = 0
gene_count = 0
gene_mapping = {}  # To store the mapping

# CONFIG — GENCODE annotation GTF (v38 used here) and the mapping CSV to write.
gtf_file_path = "path/to/gencode.v38.annotation.gtf"

# Reading the GTF file
with open(gtf_file_path, "r") as gtf:
    for line in gtf:
        line_count += 1
        if line.startswith("#"):  # Skip header lines
            continue
        parts = line.strip().split("\t")
        if parts[2] == "gene":  # Focus on gene entries
            attributes_section = parts[8]
            attributes = dict(attribute.strip().split(' ') for attribute in attributes_section.split('; ') if attribute)
            attributes = {k: v.strip('"') for k, v in attributes.items()}
            
            if 'gene_id' in attributes and 'gene_name' in attributes and 'gene_type' in attributes:
                gene_id = attributes['gene_id']
                gene_name = attributes['gene_name']
                gene_type = attributes['gene_type']
                gene_mapping[gene_id] = (gene_name, gene_type)  # Store both name and type
                gene_count += 1

print(f"Processed {line_count} lines and found {gene_count} genes.")

# Add this at the beginning of your script
output_file_path = "gene_id_to_symbol_and_type_mapping.csv"

# Write the header and each mapping
with open(output_file_path, 'w') as out_file:
    out_file.write("ENSEMBL_ID,Gene_Symbol,Gene_Type\n")  # Add Gene_Type to the header
    for gene_id, (gene_name, gene_type) in gene_mapping.items():
        out_file.write(f"{gene_id},{gene_name},{gene_type}\n")
