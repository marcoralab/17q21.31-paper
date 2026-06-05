import pandas as pd
import os
import glob
from collections import defaultdict
import re

# ===============================================================
# CONFIG — edit these paths for your environment.
# After downloading the GEO data, point each path below at the
# corresponding files/directories on your machine.
# ===============================================================

# Directories of PacBio iso-seq BED files (one BED file per sample).
BED_DIRS = [
    "/path/to/PB000732_20250227_pool1_analyzed_isoforms",
    "/path/to/PB000732_20250227_pool2_analyzed_isoforms",
    "/path/to/Kat_ISOseq_ALLDATA_dec2020/BEDfiles",
    "/path/to/PB000673/20240613_delivery",
]

# Directories containing per-sample classification/annotation .txt files
# (paired with the BED files above).
ANNOTATION_DIRS = [
    "/path/to/PB000732_20250227_pool1_analyzed_isoforms",
    "/path/to/PB000732_20250227_pool2_analyzed_isoforms",
    "/path/to/PB000673/20240613_delivery",
    "/path/to/Kat_ISOseq_ALLDATA_dec2020/Annotations",
]

# Sample metadata CSV. Required columns:
#   File_Name, HARMONIZED_NAMES, Ancestry, haplotype, Mutation, Cell_Type,
#   Brain region_Caudate, Brain region_MFG
METADATA_FILE = "/path/to/All_Meta_ISOSEQ_all_combined.csv"

# Reference transcript-to-exon mapping CSVs (H1 and H2 alts).
H1_TRANSCRIPT_EXONS_FILE = "/path/to/H1_transcirpt_name_with_exons.csv"
H2_TRANSCRIPT_EXONS_FILE = "/path/to/H2_transcirpt_name_with_exons.csv"

# Haplotype-to-isoform reference CSVs.
ISOFORM_MAPPING_FILE_H1 = "/path/to/unique_filtered_h1h1_haplotypes_isoforms.csv"
ISOFORM_MAPPING_FILE_H2 = "/path/to/unique_filtered_h2_haplotypes_isoforms.csv"

# Some samples come from an earlier PacBio call set whose BED and
# classification files use a different schema (header at the bottom,
# PB-IDs that need a PB->G remap). The processing code detects these
# "older-format" files by checking whether the file path contains these
# directory fragments. Set them to your older-format directories, or to a
# non-matching sentinel string if you have no older-format data.
OLDER_FORMAT_BED_DIR        = "/path/to/Kat_ISOseq_ALLDATA_dec2020/BEDfiles"
OLDER_FORMAT_ANNOTATION_DIR = "/path/to/Kat_ISOseq_ALLDATA_dec2020/Annotations"

# Path fragments used to detect the current PacBio iso-seq directories
# (pool1 / pool2). Used only by some logging branches.
PB_DIRECTORIES = [
    "/path/to/PB000732_20250227_pool1_analyzed_isoforms",
    "/path/to/PB000732_20250227_pool2_analyzed_isoforms",
]

# All processed outputs land here.
OUTPUT_DIR = "results/exact_exon_analysis"

# Per-sample cached PB-ID -> G-number CSVs (for the older-format dataset).
PB_TO_G_MAPPING_DIR = os.path.join(OUTPUT_DIR, "pb_to_g_mappings")
# ===============================================================

# BED file directories (one BED file per sample)
bed_dirs = BED_DIRS

# Annotation directories for per-isoform frequency data
annotation_dirs = ANNOTATION_DIRS

# Sample metadata CSV
metadata_file = METADATA_FILE

# After loading the metadata file
metadata = pd.read_csv(metadata_file)
# This analysis covers wild-type (WT) samples only.
metadata = metadata[metadata['Mutation'] == 'WT'].reset_index(drop=True)
print(f"Loaded metadata with {len(metadata)} samples")

# Define the output directory
output_dir = OUTPUT_DIR
os.makedirs(output_dir, exist_ok=True)

# Define the chromosome identifiers to check in the BED files
chromosome_identifiers = [
    'chr17_GL000258v2_alt',  # H2 chromosome
    'GL000258.2',            # H2 alternate name
    'chr17_KI270908v1_alt',  # H1 chromosome 
    'KI270908.1',            # H1 alternate name
    'chr17',                 # Generic chromosome 17
    '17'                     # Another possible chromosome 17 identifier
]

# Define the exon position to exon number mapping - H1 chromosome
h1_exon_mapping = {
    "664515-664598": "1",
    "673976-674063": "2",
    "676502-676589": "3",
    "680492-680558": "4",
    "684984-686049": "4aL",
    "685296-686049": "4a", 
    "689158-689214": "5",
    "691997-692194": "6",
    "693578-693705": "7",
    "698517-698783": "9",
    "712428-712521": "10",
    "716362-716443": "11",
    "719587-719710": "11i",
    "720731-720844": "12",
    "726067-726242": "13"
}

# Define the exon position to exon number mapping - H2 chromosome
h2_exon_mapping = {
    "825869-825952": "1",
    "816505-816592": "2",
    "813979-814066": "3", 
    "810008-810074": "4",
    "804487-805552": "4aL",
    "804487-805240": "4a",
    "801323-801379": "5",
    "798313-798511": "6",
    "796802-796929": "7",
    "791725-791991": "9",
    "778232-778325": "10",
    "774310-774392": "11",
    "769883-769996": "12",
    "764488-764663": "13"
}

# Isoform mapping file paths
isoform_mapping_file_h1 = ISOFORM_MAPPING_FILE_H1
isoform_mapping_file_h2 = ISOFORM_MAPPING_FILE_H2

# Path to the folder with PB to G mappings
pb_to_g_mapping_dir = PB_TO_G_MAPPING_DIR

# Load the metadata file
metadata = pd.read_csv(metadata_file)
# This analysis covers wild-type (WT) samples only.
metadata = metadata[metadata['Mutation'] == 'WT'].reset_index(drop=True)
print(f"Loaded metadata with {len(metadata)} samples")

# Add this function after loading the metadata file
def determine_cell_type(row):
    if pd.notna(row['Cell_Type']) and row['Cell_Type'] != 'Brain':
        return row['Cell_Type']
    elif pd.notna(row['Brain region_Caudate']) and row['Brain region_Caudate'] == 'x':
        return 'Brain_Caudate'
    elif pd.notna(row['Brain region_MFG']) and row['Brain region_MFG'] == 'x':
        return 'Brain_MFG'
    else:
        return 'Unknown'

# Apply this function to metadata and create a new mapping dictionary
metadata['Derived_Cell_Type'] = metadata.apply(determine_cell_type, axis=1)

# Add to your existing dictionary mapping code
sample_to_derived_cell_type = {}
for _, row in metadata.iterrows():
    if 'File_Name' in row and pd.notna(row['File_Name']):
        sample_name = row['File_Name']
        derived_cell_type = row['Derived_Cell_Type']
        
        # Add to mapping dictionaries
        sample_to_derived_cell_type[sample_name] = derived_cell_type
        
        # Also add to partial/numeric mappings as you do with ancestry and haplotype
        name_parts = sample_name.split('_')
        for i in range(len(name_parts)):
            partial_name = '_'.join(name_parts[i:])
            if partial_name and partial_name != sample_name:
                sample_to_derived_cell_type[partial_name] = derived_cell_type
                
        # Map numeric portions
        numbers = re.findall(r'\d+', sample_name)
        for num in numbers:
            if len(num) >= 3:  # Only use reasonably long numbers
                sample_to_derived_cell_type[num] = derived_cell_type

# Create a dictionary to map sample names to ancestry, haplotype, and cell type
sample_to_ancestry = {}
sample_to_haplotype = {}
sample_to_Cell_Type = {}

for _, row in metadata.iterrows():
    if 'File_Name' in row and pd.notna(row['File_Name']):
        sample_name = row['File_Name']
        ancestry = row.get('Ancestry', 'Unknown')
        haplotype = row.get('haplotype', 'Unknown')
        Cell_Type = row.get('Cell_Type', 'Unknown')

        # Store the original sample name mapping
        sample_to_ancestry[sample_name] = ancestry
        sample_to_haplotype[sample_name] = haplotype
        sample_to_Cell_Type[sample_name] = Cell_Type
        
        # For sample names like "Mayo_5_cau", also map "95_Mayo_5_cau" format
        # by creating mappings for each part of the name
        name_parts = sample_name.split('_')
        for i in range(len(name_parts)):
            # Create mappings for partial matches (Mayo, Mayo_5, etc.)
            partial_name = '_'.join(name_parts[i:])
            if partial_name and partial_name != sample_name:
                sample_to_ancestry[partial_name] = ancestry
                sample_to_haplotype[partial_name] = haplotype
                sample_to_Cell_Type[partial_name] = Cell_Type
                
        # Also map numeric portions
        import re
        numbers = re.findall(r'\d+', sample_name)
        for num in numbers:
            if len(num) >= 3:  # Only use reasonably long numbers
                sample_to_ancestry[num] = ancestry
                sample_to_haplotype[num] = haplotype
                sample_to_Cell_Type[num] = Cell_Type
                
print(f"Created mapping for {len(sample_to_ancestry)} samples with ancestry and haplotype information")

# Function to get PB to G mapping for a sample from cache files
def get_pb_to_g_mapping_for_sample(sample_name):
    """Get the PB to G mapping for a specific sample from the cached CSV files."""
    # Create a cache file name pattern to look for
    if not isinstance(sample_name, str):
        sample_name = str(sample_name) 
    
    # If we have a direct match for the sample name
    cache_file = os.path.join(pb_to_g_mapping_dir, f"{sample_name}_pb_to_g_mapping.csv")
    if os.path.exists(cache_file):
        try:
            mapping_df = pd.read_csv(cache_file)
            return dict(zip(mapping_df['PB_Number'], mapping_df['G_Number']))
        except Exception as e:
            print(f"Error loading cached mapping for {sample_name}: {e}")
    
    # If not, check if we have a related mapping file
    if os.path.exists(pb_to_g_mapping_dir):
        try:
            # Get list of all mapping files
            all_files = [f for f in os.listdir(pb_to_g_mapping_dir) if f.endswith('_pb_to_g_mapping.csv')]
            
            # Look for a file that contains a part of the sample name
            sample_parts = sample_name.split('_')
            for part in sample_parts:
                if len(part) >= 3:  # Only use parts that are long enough to be meaningful
                    for file in all_files:
                        if part in file:
                            try:
                                mapping_df = pd.read_csv(os.path.join(pb_to_g_mapping_dir, file))
                                print(f"Using mapping from {file} for sample {sample_name}")
                                return dict(zip(mapping_df['PB_Number'], mapping_df['G_Number']))
                            except Exception as e:
                                print(f"Error loading cached mapping from {file}: {e}")
        except Exception as e:
            print(f"Error searching for mapping files: {e}")
    
    # If we reach here, no mapping was found
    return {}

# Function to extract exon structure from a BED file row
def get_exon_structure(block_sizes, block_starts, start_pos):
    if block_sizes and block_starts:
        sizes = [int(s) for s in block_sizes.split(',') if s]
        starts = [int(s) for s in block_starts.split(',') if s]
        if len(sizes) == len(starts):
            exons = []
            for size, start in zip(sizes, starts):
                exon_start = start_pos + start
                exon_end = exon_start + size
                exons.append((exon_start, exon_end))
            return exons
    return []

# Function to get the annotation file path for a sample
def get_annotation_file_path(sample_name, verbose=True):
    # Ensure sample_name is a string
    if sample_name is None:
        sample_name = "unknown"
    elif not isinstance(sample_name, str):
        sample_name = str(sample_name)
    
    # Try direct matching for annotation files
    classification_patterns = [
        f"{sample_name}.classification.final.txt",
        f"{sample_name}.clustered.classification.final.txt"
    ]
    
    # First, search for classification files specifically
    for annot_dir in annotation_dirs:
        for pattern in classification_patterns:
            # Skip wildcard patterns at this stage
            if '*' in pattern:
                continue
                
            annot_file = os.path.join(annot_dir, pattern)
            if os.path.exists(annot_file):
                # Determine the file type safely
                try:
                    with open(annot_file, 'r') as f:
                        lines = f.readlines()
                        
                        # Check if the file is empty
                        if not lines:
                            if verbose:
                                print(f"Empty file: {annot_file}")
                            continue
                        
                        # Check if header might be at the bottom
                        last_line = lines[-1].strip()
                        header_at_bottom = ('isoform' in last_line.lower() and 'FL' in last_line.lower()) or \
                                          ('transcript' in last_line.lower() and 'supporting_reads' in last_line.lower())
                        
                        if header_at_bottom:
                            if verbose:
                                print(f"Found classification file with header at bottom: {annot_file}")
                            return annot_file, "bottom_header"
                        
                        # Check if this is a older-format file
                        is_older_format = OLDER_FORMAT_ANNOTATION_DIR in annot_file
                        if is_older_format:
                            # For older-format files, check if header contains 'FL'
                            file_type = 'older-format'
                            if verbose:
                                print(f"Found older-format classification file: {annot_file} (type: {file_type})")
                        else:
                            # For non-older-format files, use the first line header
                            header_line = lines[0].strip()
                            num_fields = len(header_line.split('\t'))
                            file_type = 'first' if num_fields < 20 else 'second'
                            if verbose:
                                print(f"Found classification file: {annot_file} (type: {file_type}, cols: {num_fields})")
                    return annot_file, file_type
                except Exception as e:
                    if verbose:
                        print(f"Error reading {annot_file}: {e}")
                    continue
    
    # If more relaxed matching needed, look for files containing sample_name
    for annot_dir in annotation_dirs:
        if os.path.exists(annot_dir):
            try:
                all_files = [f for f in os.listdir(annot_dir) if sample_name in f and f.endswith('.txt')]
                for file in all_files:
                    annot_file = os.path.join(annot_dir, file)
                    try:
                        with open(annot_file, 'r') as f:
                            lines = f.readlines()
                            if not lines:
                                continue
                            
                            # Check header location
                            last_line = lines[-1].strip()
                            header_at_bottom = ('isoform' in last_line.lower() and 'FL' in last_line.lower()) or \
                                              ('transcript' in last_line.lower() and 'supporting_reads' in last_line.lower())
                            
                            if header_at_bottom:
                                return annot_file, "bottom_header"
                            
                            is_older_format = OLDER_FORMAT_ANNOTATION_DIR in annot_file
                            if is_older_format:
                                return annot_file, 'older-format'
                            else:
                                header_line = lines[0].strip()
                                num_fields = len(header_line.split('\t'))
                                file_type = 'first' if num_fields < 20 else 'second'
                                return annot_file, file_type
                    except:
                        continue
            except:
                continue
    
    return None, None

# Function to process frequency data from annotation files
def process_annotation_file(file_path, file_type, pb_to_g_mapping=None):
    """
    Processes the annotation file, calculates frequencies, and returns a dictionary mapping
    transcript names to their frequency data.
    """
    if file_path is None:
        return None
        
    # Extract sample name from file path for debugging
    sample_name = os.path.basename(file_path).split('.')[0]
    
    # Check if this is a older-format file
    is_older_format = OLDER_FORMAT_ANNOTATION_DIR in file_path
    if is_older_format:
        print(f"Processing older-format annotation file: {file_path}")
    
    # Check if this is from one of the PB ID directories
    is_pb_directory = any(path in file_path for path in PB_DIRECTORIES)
    
    if is_pb_directory:
        print(f"Processing file from PB directory: {file_path}")
    
    try:
        # Read the entire file into memory to check for header location
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            print(f"Warning: Empty file {file_path}")
            return None
            
        # Check if the last line contains a header (common in these files)
        last_line = lines[-1].strip().lower()
        header_at_bottom = ('isoform' in last_line and 'fl' in last_line) or \
                          ('transcript' in last_line and 'supporting_reads' in last_line)
        
        if header_at_bottom:
            # Use the last line as header and all other lines as data
            header = lines[-1].strip().split('\t')
            data_lines = [line.strip().split('\t') for line in lines[:-1]]
            
            # Create DataFrame manually
            df = pd.DataFrame(data_lines, columns=header)
            print(f"Created DataFrame with header from bottom of file. Columns: {df.columns.tolist()}")
            
            # Determine correct column names regardless of case
            transcript_col = None
            length_col = None
            read_count_col = None
            
            # Find the transcript column
            for col in df.columns:
                if 'isoform'.lower() == str(col).lower() or 'transcript'.lower() in str(col).lower():
                    transcript_col = col
                    break
            if transcript_col is None:
                transcript_col = df.columns[0]  # First column is usually the transcript
                
            # Find the length column
            for col in df.columns:
                if 'length'.lower() == str(col).lower():
                    length_col = col
                    break
            if length_col is None and len(df.columns) > 3:
                length_col = df.columns[3]  # Length is often in the 4th column
                
            # Find the read count column - check for both 'FL' and 'fl' 
            for col in df.columns:
                if 'fl'.lower() == str(col).lower() or 'supporting_reads'.lower() in str(col).lower():
                    read_count_col = col
                    break
            if read_count_col is None:
                read_count_col = df.columns[-1]  # Last column is usually the read count
                
            print(f"Using columns: transcript={transcript_col}, length={length_col}, read_count={read_count_col}")
            
        elif is_older_format:
            # For older-format files, use pandas with header=0
            df = pd.read_csv(file_path, sep='\t', header=0)
            print(f"Created DataFrame from older-format file. Columns: {df.columns.tolist()}")
            
            # Enhanced column detection for older-format files
            transcript_col = None
            for possible_col in ['isoform', 'Isoform', 'transcript', 'Transcript', 'Transcript_Name']:
                if possible_col in df.columns:
                    transcript_col = possible_col
                    break
            
            length_col = None
            for possible_col in ['length', 'Length']:
                if possible_col in df.columns:
                    length_col = possible_col
                    break
            
            read_count_col = None
            for possible_col in ['FL', 'fl', 'Supporting_Reads', 'supporting_reads']:
                if possible_col in df.columns:
                    read_count_col = possible_col
                    break
            
            # If any columns weren't found, use positional fallbacks
            if transcript_col is None and len(df.columns) > 0:
                transcript_col = df.columns[0]
                print(f"Using positional fallback for transcript column: {transcript_col}")
            if length_col is None and len(df.columns) > 3:
                length_col = df.columns[3]
                print(f"Using positional fallback for length column: {length_col}")
            if read_count_col is None and len(df.columns) > 0:
                read_count_col = df.columns[-1]
                print(f"Using positional fallback for read count column: {read_count_col}")
                
            print(f"Using columns for older-format file: transcript={transcript_col}, length={length_col}, read_count={read_count_col}")
            
        else:
            # For other files, try standard parsing
            try:
                # First try with header=0 in case the file actually has a header
                df = pd.read_csv(file_path, sep='\t', header=0)
                first_col = df.columns[0]
                # Check if what was read as a header actually looks like data
                if first_col.isdigit() or first_col.startswith('PB') or first_col.startswith('CP') or first_col.startswith('G'):
                    # This looks like the first row is data, not header
                    # Read again with no header
                    df = pd.read_csv(file_path, sep='\t', header=None)
            except:
                # If that fails, try without a header
                df = pd.read_csv(file_path, sep='\t', header=None)
            
            # Assign column names based on file type
            if file_type == 'first' and len(df.columns) >= 9:
                column_names = ['Transcript_Name', 'Chromosome', 'Strand', 'Length', 'Exons', 
                               'Category', 'Gene', 'Gene_Type', 'Supporting_Reads']
                # Extend with placeholder names if needed
                if len(df.columns) > len(column_names):
                    column_names.extend([f'Extra_{i}' for i in range(len(column_names), len(df.columns))])
                
                df.columns = column_names[:len(df.columns)]
                print(f"Assigned first-type column names. Columns: {df.columns.tolist()}")
                
                transcript_col = 'Transcript_Name'
                length_col = 'Length'
                read_count_col = 'Supporting_Reads'
                
            elif file_type == 'second' and len(df.columns) >= 9:
                # For second type with FL column
                # Create basic column names
                column_names = ['Transcript_Name', 'Chromosome', 'Strand', 'Length', 'Exons']
                # Find the FL column position - usually last or at a specific index
                if len(df.columns) >= 23:  # If we have enough columns for standard second type
                    column_names.extend([f'Col_{i}' for i in range(len(column_names), 22)])
                    column_names.append('FL')
                    if len(df.columns) > 23:
                        column_names.extend([f'Extra_{i}' for i in range(23, len(df.columns))])
                else:
                    # For shorter files, add placeholders and assume last column is FL
                    column_names.extend([f'Col_{i}' for i in range(len(column_names), len(df.columns)-1)])
                    column_names.append('FL')
                
                df.columns = column_names
                print(f"Assigned second-type column names. Columns: {df.columns.tolist()}")
                
                transcript_col = 'Transcript_Name'
                length_col = 'Length'
                read_count_col = 'FL'
            else:
                print(f"Warning: Unknown file structure for {file_path} with {len(df.columns)} columns")
                if len(df.columns) >= 3:
                    # Attempt a basic mapping for unknown formats
                    print(f"Attempting basic column mapping for unknown format")
                    df.columns = ["Transcript_Name"] + [f"Col_{i}" for i in range(1, len(df.columns))]
                    transcript_col = 'Transcript_Name'
                    
                    # Try to identify length and read count columns by checking data types
                    numeric_cols = []
                    for col in df.columns[1:]:  # Skip the transcript column
                        try:
                            # Check if column has mostly numeric values
                            numeric_vals = pd.to_numeric(df[col], errors='coerce')
                            if numeric_vals.notna().mean() > 0.5:  # If more than 50% are numbers
                                numeric_cols.append(col)
                        except:
                            continue
                    
                    if len(numeric_cols) >= 2:
                        # Assume last numeric column is read count, second-to-last is length
                        read_count_col = numeric_cols[-1]
                        length_col = numeric_cols[-2]
                        print(f"Guessed columns: length={length_col}, read_count={read_count_col}")
                    else:
                        # Can't reliably guess, return empty
                        print(f"Unable to identify numeric columns for length and read count")
                        return {}
                else:
                    print(f"File has too few columns, cannot process")
                    return {}
        
        # For older-format files, if we still have issues with column identification
        if is_older_format and (transcript_col not in df.columns or length_col not in df.columns or read_count_col not in df.columns):
            # Try to identify numeric columns
            potential_numeric_cols = []
            for col in df.columns:
                try:
                    numeric_values = pd.to_numeric(df[col], errors='coerce')
                    non_null_pct = numeric_values.notna().mean()
                    if non_null_pct > 0.5:  # Column is mostly numeric
                        potential_numeric_cols.append((col, non_null_pct))
                except:
                    continue
            
            # Sort by percentage of numeric values
            potential_numeric_cols.sort(key=lambda x: x[1], reverse=True)
            
            # If we found potential numeric columns, use them
            if len(potential_numeric_cols) >= 2:
                if read_count_col not in df.columns:
                    read_count_col = potential_numeric_cols[0][0]
                    print(f"Auto-detected read count column: {read_count_col}")
                if length_col not in df.columns:
                    length_col = potential_numeric_cols[1][0]
                    print(f"Auto-detected length column: {length_col}")
        
        # Make sure the columns exist
        if transcript_col not in df.columns:
            print(f"Warning: Transcript column '{transcript_col}' not found in {file_path}")
            print(f"Available columns: {df.columns.tolist()}")
            return {}
            
        if length_col not in df.columns:
            print(f"Warning: Length column '{length_col}' not found in {file_path}")
            print(f"Available columns: {df.columns.tolist()}")
            return {}
            
        if read_count_col not in df.columns:
            print(f"Warning: Read count column '{read_count_col}' not found in {file_path}")
            print(f"Available columns: {df.columns.tolist()}")
            return {}
        
        # Convert columns to numeric
        df[length_col] = pd.to_numeric(df[length_col], errors='coerce')
        df[read_count_col] = pd.to_numeric(df[read_count_col], errors='coerce')
        
        # Calculate frequencies
        df['Raw_Frequency'] = df[read_count_col]
        df['Normalized_Frequency'] = df[read_count_col] / (df[length_col] + 1e-10)
        
        # Create a mapping of transcript to frequency data
        transcript_frequencies = {}
        for _, row in df.iterrows():
            # Get transcript name and strip quotes if present
            transcript_name = str(row[transcript_col]).strip('"\'')
            if not transcript_name:
                continue

            # Create mappings for both the full name and the base name (without version)
            transcript_frequencies[transcript_name] = {
                'Raw_Frequency': row['Raw_Frequency'],
                'Normalized_Frequency': row['Normalized_Frequency'],
                'Read_Count': row[read_count_col],
                'Length': row[length_col]
            }
            
            # Also add a mapping for the base name (without version) if different
            if '.' in transcript_name:
                base_name = transcript_name.split('.')[0]
                if base_name != transcript_name:
                    transcript_frequencies[base_name] = transcript_frequencies[transcript_name]
        
        print(f"Created frequency mapping for {len(transcript_frequencies)} transcripts")
        
        # Count by prefix type
        prefix_counts = {
            'PB': sum(1 for t in transcript_frequencies if str(t).startswith('PB')),
            'CP': sum(1 for t in transcript_frequencies if str(t).startswith('CP')),
            'G': sum(1 for t in transcript_frequencies if str(t).startswith('G')),
            'other': sum(1 for t in transcript_frequencies if not any(str(t).startswith(p) for p in ['PB', 'CP', 'G']))
        }
        print(f"Transcript prefix distribution: {prefix_counts}")
        
        # Apply PB to G mapping if available and this is a older-format file
        if is_older_format and pb_to_g_mapping:
            # Create a new mapping with G numbers when available
            transcript_frequencies_with_g = {}
            mapped_count = 0
            
            for transcript_name, freq_info in transcript_frequencies.items():
                # Add the original mapping
                transcript_frequencies_with_g[transcript_name] = freq_info.copy()
                
                # If this transcript has a G number mapping, add that too
                if transcript_name in pb_to_g_mapping:
                    g_number = pb_to_g_mapping[transcript_name]
                    # Add G number as an additional key
                    transcript_frequencies_with_g[g_number] = freq_info.copy()
                    mapped_count += 1
                
                # Also try with base name if original has a version number
                elif '.' in transcript_name:
                    base_name = transcript_name.split('.')[0]
                    if base_name in pb_to_g_mapping:
                        g_number = pb_to_g_mapping[base_name]
                        transcript_frequencies_with_g[g_number] = freq_info.copy()
                        mapped_count += 1
            
            print(f"Added {mapped_count} G number mappings for older-format transcripts")
            
            # Replace the original dictionary with the expanded one
            transcript_frequencies = transcript_frequencies_with_g
        
        return transcript_frequencies
        
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        import traceback
        traceback.print_exc()
        return {}

def map_exons_to_names(exons, chromosome, tolerance_percentage=0.1):
    """
    Maps exon coordinates to their symbolic names (1, 2, 3, 4, etc.)
    with a tolerance of ±tolerance_percentage of exon length.
    
    Args:
        exons: List of tuples [(start1, end1), (start2, end2), ...]
        chromosome: Chromosome identifier
        tolerance_percentage: Percentage of exon length to allow for variation (default: 0.1)
        
    Returns:
        String of comma-separated exon names in original order
    """
    # Choose the appropriate mapping based on chromosome
    if any(h2_id in chromosome for h2_id in ["GL000258v2_alt", "GL000258.2"]):
        mapping = h2_exon_mapping
    else:
        # Default to H1 mapping for other cases
        mapping = h1_exon_mapping
    
    result = []
    
    # Process each exon in the original order
    for start, end in exons:
        found_match = False
        exon_length = end - start
        
        # Check each mapping entry
        for key, value in mapping.items():
            ref_start, ref_end = map(int, key.split('-'))
            ref_length = ref_end - ref_start
            
            # Calculate tolerance based on average length
            avg_length = (exon_length + ref_length) / 2
            tolerance = int(avg_length * tolerance_percentage)
            
            # Check if within tolerance
            if abs(start - ref_start) <= tolerance and abs(end - ref_end) <= tolerance:
                result.append(value)
                found_match = True
                break
        
        # If no match found, add the original coordinates in brackets
        if not found_match:
            result.append(f"[{start}-{end}]")
    
    # Return comma-separated string
    return ", ".join(result) if result else ""

def determine_isoform_type(exons_present):
    """
    Determines the isoform type based on the presence of specific exons:
    - 2N/1N/0N based on exons 2 and 3
    - 4a and 4aL suffixes if present
    - 4R/3R suffix based on exons 9, 10, 11, and 12
    
    Args:
        exons_present: String of comma-separated exon names
    
    Returns:
        String representing the isoform type
    """
    if not exons_present:
        return "Unknown"
    
    # Split into individual exons
    exons = [e.strip() for e in exons_present.split(',')]
    
    # Check for exons 2 and 3
    has_exon2 = "2" in exons
    has_exon3 = "3" in exons
    
    # Determine N status
    if has_exon2 and has_exon3:
        n_status = "2N"
    elif has_exon2 or has_exon3:
        n_status = "1N"
    else:
        n_status = "0N"
    
    # Check for 4a and 4aL
    has_4a = "4a" in exons
    has_4aL = "4aL" in exons
    
    # Check for R-region exons (9, 10, 11, 12)
    has_exon9 = "9" in exons
    has_exon10 = "10" in exons
    has_exon11 = "11" in exons
    has_exon12 = "12" in exons
    
    # Build the full isoform type
    isoform_type = n_status
    
    # Determine R status - need to have at least exons 9, 11, and 12 to count as R
 #   if has_exon9 and has_exon11 and has_exon12:
 #       # 4R if all four, 3R if missing exon 10
 #       if has_exon10:
 #           isoform_type += "4R"
 #       else:
 #           isoform_type += "3R"
 
    # Determine R status based ONLY on presence/absence of exon 10
    if has_exon10:
        isoform_type += "4R"
    else:
        isoform_type += "3R"
    
    # Add 4a/4aL suffixes if present
    if has_4a:
        isoform_type += "4a"
    
    if has_4aL:
        isoform_type += "4aL"
    
    return isoform_type


# Function to process a single BED file
def process_bed_file(file_path):
    import re
    
    transcripts = []
    filename = os.path.basename(file_path)
    # Extract sample name as the first token before a dot
    sample_name = filename.split('.')[0]
    
    print(f"\nProcessing BED file for sample: {sample_name}")

    # Check if this is a older-format bed file
    is_older_format_bed = OLDER_FORMAT_BED_DIR in file_path
    pb_to_g_mapping = {}
    
    # Get ancestry and haplotype information for this sample
    ancestry = sample_to_ancestry.get(sample_name, "Unknown")
    haplotype = sample_to_haplotype.get(sample_name, "Unknown")
    Cell_Type = sample_to_Cell_Type.get(sample_name, "Unknown")
    
    Derived_Cell_Type = sample_to_derived_cell_type.get(sample_name, "Unknown")

    ancestry_haplotype = f"{ancestry}_{haplotype}" if ancestry != "Unknown" and haplotype != "Unknown" else "Unknown"
    # Also try alternate sample name formats if not found
    if ancestry == "Unknown":
        # Try removing numeric prefix if present (e.g., "95_Mayo_5_cau" → "Mayo_5_cau")
        if '_' in sample_name:
            # Check if it starts with a number followed by underscore
            prefix_match = re.match(r'^\d+_(.+)$', sample_name)
            # For example, in the numeric prefix matching section:
            if prefix_match:
                # Extract the part after the number_
                without_prefix = prefix_match.group(1)
                ancestry = sample_to_ancestry.get(without_prefix, "Unknown")
                haplotype = sample_to_haplotype.get(without_prefix, "Unknown")
                Cell_Type = sample_to_Cell_Type.get(without_prefix, "Unknown")
                Derived_Cell_Type = sample_to_derived_cell_type.get(without_prefix, "Unknown")
                # Update the combined field
                ancestry_haplotype = f"{ancestry}_{haplotype}" if ancestry != "Unknown" and haplotype != "Unknown" else "Unknown"
                if ancestry != "Unknown":
                    print(f"Found ancestry info after removing numeric prefix: {without_prefix} for sample {sample_name}")

        # If still unknown, try with first part after underscore
        if ancestry == "Unknown" and '_' in sample_name:
            alt_sample_name = sample_name.split('_')[0]
            ancestry = sample_to_ancestry.get(alt_sample_name, "Unknown")
            haplotype = sample_to_haplotype.get(alt_sample_name, "Unknown")

    # Also try matching by just the numeric part
    if ancestry == "Unknown":
        import re
        # Extract numbers from the sample name
        number_match = re.findall(r'\d+', sample_name)
        if number_match:
            for num in number_match:
                # Only try with reasonably long number sequences
                if len(num) >= 3:  
                    alt_sample_name = num
                    ancestry = sample_to_ancestry.get(alt_sample_name, "Unknown")
                    haplotype = sample_to_haplotype.get(alt_sample_name, "Unknown")
                    if ancestry != "Unknown":
                        print(f"Found ancestry info using numeric pattern {num} for sample {sample_name}")
                        break
        
    if is_older_format_bed:
        print(f"Processing older-format BED file: {file_path}")
        # Get PB to G mapping for older-format BED files 
        pb_to_g_mapping = get_pb_to_g_mapping_for_sample(sample_name)
        if pb_to_g_mapping:
            print(f"Found PB to G mapping with {len(pb_to_g_mapping)} entries for older-format BED file")
        else:
            print(f"WARNING: No PB to G mapping found for older-format BED file {file_path}")
            # Try alternate sample name formats
            if '_' in sample_name:
                alt_sample_name = sample_name.split('_')[0]
                print(f"Trying alternative sample name: {alt_sample_name}")
                pb_to_g_mapping = get_pb_to_g_mapping_for_sample(alt_sample_name)
                if pb_to_g_mapping:
                    print(f"Found PB to G mapping with {len(pb_to_g_mapping)} entries using alternate name")
    
    # Check for a corresponding renaming key file (e.g., 181_HBCC_1314_DLPFC.renaming_key.txt)
    key_file = os.path.join(os.path.dirname(file_path), sample_name + ".renaming_key.txt")
    renaming_mapping = {}
    cp_to_pb_mapping = {}  # The reverse mapping we need (CP to PB)

    if os.path.exists(key_file):
        try:
            print(f"Found renaming key file: {key_file}")
            with open(key_file, 'r') as kf:
                for line in kf:
                    line = line.strip()
                    if line:
                        # Expecting each line to be: PB_name  CP_name
                        items = line.split()
                        if len(items) >= 2:
                            pb_name, cp_name = items[0], items[1]
                            # Store the mapping in both directions
                            renaming_mapping[cp_name] = pb_name  # CP → PB (what we primarily need)
                            # Also add mapping for CP in other formats (with and without dot)
                            cp_base = cp_name.split('.')[0] if '.' in cp_name else cp_name
                            renaming_mapping[cp_base] = pb_name
            
            # Debug output
            print(f"Loaded {len(renaming_mapping)} CP to PB ID mappings from {key_file}")
            sample_keys = list(renaming_mapping.keys())[:5]
            if sample_keys:
                print("Sample mappings:")
                for key in sample_keys:
                    print(f"  {key} → {renaming_mapping[key]}")
        except Exception as e:
            print(f"Error reading renaming key file {key_file}: {e}")
    else:
        print(f"No renaming key file found for {sample_name}")

    # Get frequency data for this sample - try multiple sample name variants
    frequency_data = {}
    frequency_sources = []

    # Try with original sample name
    annot_file, file_type = get_annotation_file_path(sample_name)

    # Check if this is a older-format file
    is_older_format = False
    if annot_file:
        is_older_format = OLDER_FORMAT_ANNOTATION_DIR in annot_file

    # Get PB to G mapping if this is a older-format file (reuse the mapping from above if exists)
    if is_older_format and not pb_to_g_mapping:
        pb_to_g_mapping = get_pb_to_g_mapping_for_sample(sample_name)
            
        if pb_to_g_mapping:
            print(f"Found PB to G mapping with {len(pb_to_g_mapping)} entries for older-format sample {sample_name}")
        else:
            print(f"No PB to G mapping found for older-format sample {sample_name} - may affect mapping")
            
    if annot_file:
        frequency_data = process_annotation_file(annot_file, file_type, pb_to_g_mapping) or {}
        if frequency_data:
            frequency_sources.append(f"original name ({sample_name})")
            print(f"Found frequency data for {len(frequency_data)} transcripts using sample name: {sample_name}")
    
    # If no data found and sample name contains underscores, try with first part
    if not frequency_data and '_' in sample_name:
        alt_sample_name = sample_name.split('_')[0]
        annot_file, file_type = get_annotation_file_path(alt_sample_name)
        if annot_file:
            # Check if this is a older-format file
            is_older_format_alt = OLDER_FORMAT_ANNOTATION_DIR in annot_file
            if is_older_format_alt and not pb_to_g_mapping:
                pb_to_g_mapping = get_pb_to_g_mapping_for_sample(alt_sample_name)
                
            frequency_data = process_annotation_file(annot_file, file_type, pb_to_g_mapping) or {}
            if frequency_data:
                frequency_sources.append(f"alternate name ({alt_sample_name})")
                print(f"Found frequency data for {len(frequency_data)} transcripts using alternate sample name: {alt_sample_name}")
    
    # If still no data found and filename contains numbers, try with just numbers
    if not frequency_data:
        import re
        numbers = re.findall(r'\d+', sample_name)
        if numbers:
            for num in numbers:
                if len(num) >= 3:  # Only try with reasonably long number sequences
                    annot_file, file_type = get_annotation_file_path(num)
                    if annot_file:
                        # Check if this is a older-format file
                        is_older_format_num = OLDER_FORMAT_ANNOTATION_DIR in annot_file
                        if is_older_format_num and not pb_to_g_mapping:
                            pb_to_g_mapping = get_pb_to_g_mapping_for_sample(num)
                            
                        frequency_data = process_annotation_file(annot_file, file_type, pb_to_g_mapping) or {}
                        if frequency_data:
                            frequency_sources.append(f"number ({num})")
                            print(f"Found frequency data for {len(frequency_data)} transcripts using number: {num}")
                            break
    
    # Log frequency data keys for debugging
    if frequency_data:
        print(f"Frequency data source: {', '.join(frequency_sources)}")
        print(f"First 5 frequency data keys: {list(frequency_data.keys())[:5]}")
        
        # Count prefix types in frequency data
        prefix_counts = {
            'PB': sum(1 for t in frequency_data if str(t).startswith('PB')),
            'CP': sum(1 for t in frequency_data if str(t).startswith('CP')),
            'G': sum(1 for t in frequency_data if str(t).startswith('G')),
            'other': sum(1 for t in frequency_data if not any(str(t).startswith(p) for p in ['PB', 'CP', 'G']))
        }
        print(f"Frequency data prefix distribution: {prefix_counts}")
    else:
        print(f"No frequency data found for sample {sample_name}")
    
    try:
        # Read the BED file
        bed_data = pd.read_csv(file_path, sep='\t', header=None)
        if len(bed_data.columns) < 12:
            print(f"Warning: {file_path} has fewer than 12 columns. Skipping...")
            return []
        # Assign column names
        bed_data.columns = ["chrom", "start", "end", "name", "score", "strand", "thick_start", 
                              "thick_end", "rgb", "block_count", "block_sizes", "block_starts"] + \
                              [f"extra{i}" for i in range(len(bed_data.columns) - 12)]
        
        # Extract transcripts from chromosomes of interest
        matched_chroms = []
        for chrom in chromosome_identifiers:
            region_data = bed_data[bed_data['chrom'] == chrom].copy()
            if not region_data.empty:
                matched_chroms.append(f"{chrom} ({len(region_data)})")
                # Inside the loop that processes each transcript
                for _, row in region_data.iterrows():
                    # Get the transcript name from the BED file
                    original_name = row['name']
                    
                    # Initialize variables
                    transcript_name = None
                    pb_id_for_freq = None
                    cp_id = None
                    
                    # Check if this is a CP ID from one of the PB directories
                    if not is_older_format_bed and original_name.startswith('CP.'):
                        cp_id = original_name  # Store the original CP ID
                        
                        # First use CP ID for frequency lookup
                        pb_id_for_freq = cp_id
                        
                        # Check if we have a mapping to convert CP to PB
                        mapped_pb = None
                        if cp_id in renaming_mapping:
                            mapped_pb = renaming_mapping[cp_id]
                            print(f"  Mapped CP to PB: {cp_id} → {mapped_pb}")
                            
                            # For CP IDs, we'll use the PB ID as the transcript name
                            transcript_name = mapped_pb
                        else:
                            # If no mapping found, use the CP ID
                            transcript_name = cp_id
                            
                    # For older-format files, extract as before
                    elif is_older_format_bed and '|' in original_name:
                        parts = original_name.split('|')
                        if len(parts) >= 2:
                            transcript_name = parts[0]  # Use G identifier for older-format files
                            pb_id_for_freq = parts[1]   # Use PB ID for frequency lookup
                        else:
                            transcript_name = parts[0]
                            pb_id_for_freq = parts[0]
                    else:
                        # For other files, use standard handling
                        transcript_name = original_name.split('|')[0] if '|' in original_name else original_name
                        pb_id_for_freq = transcript_name
                    
                    # CRITICAL: If we're using PB ID as transcript name but CP ID for lookup,
                    # we need to add both to possible_names to ensure we find frequency data
                    possible_names = []
                    
                    # Add original IDs first
                    if cp_id:
                        possible_names.append(cp_id)  # Add CP ID for frequency lookup
                    
                    possible_names.append(pb_id_for_freq)  # Add PB ID (or whatever we're using for freq) 
                    
                    # For non-older-format files, we need to try both with the original CP ID and the mapped PB ID
                    if not is_older_format_bed and cp_id and transcript_name != cp_id:
                        # Make sure both are in possible_names
                        if transcript_name not in possible_names:
                            possible_names.append(transcript_name)
                    
                    # Add variant without version number if it has one
                    for name in list(possible_names):  # Create a copy of the list to avoid modification during iteration
                        if '.' in name:
                            base_name = name.split('.')[0]
                            if base_name not in possible_names:
                                possible_names.append(base_name)
                    
                    # Special handling for older-format BED files - add G number variants
                    g_number = None
                    if is_older_format_bed and (pb_id_for_freq.startswith('PB') or pb_id_for_freq.startswith('pb')):
                        # Look for the transcript in the mapping
                        for pb_variant in [pb_id_for_freq, pb_id_for_freq.upper(), pb_id_for_freq.lower()]:
                            if pb_variant in pb_to_g_mapping:
                                g_number = pb_to_g_mapping[pb_variant]
                                possible_names.append(g_number)
                                print(f"Added G number for older-format BED transcript: {pb_id_for_freq} → {g_number}")
                                break
                        
                        # Also try with base name if no match found and has version
                        if g_number is None and '.' in pb_id_for_freq:
                            base_name = pb_id_for_freq.split('.')[0]
                            if base_name in pb_to_g_mapping:
                                g_number = pb_to_g_mapping[base_name]
                                possible_names.append(g_number)
                                print(f"Added G number for base name: {base_name} → {g_number}")
                    
                    # Look up frequency data using all possible name variants
                    freq_info = None
                    matched_name = None
                    for name_var in possible_names:
                        if name_var in frequency_data:
                            freq_info = frequency_data[name_var]
                            matched_name = name_var
                            print(f"  Found frequency match: {name_var} in frequency data")
                            break
                    
                    # For CP transcripts that still have no freq match, try prefix matching
                    if not freq_info and cp_id and '.' in cp_id:
                        cp_parts = cp_id.split('.')
                        for freq_key in frequency_data.keys():
                            freq_parts = freq_key.split('.')
                            if len(cp_parts) >= 3 and len(freq_parts) >= 3 and freq_key.startswith('CP.'):
                                if cp_parts[0] == freq_parts[0] and cp_parts[1] == freq_parts[1]:
                                    # Match CP prefix and subtype (CP.H1, etc)
                                    freq_info = frequency_data[freq_key]
                                    matched_name = freq_key
                                    print(f"  Found frequency match by CP pattern: {cp_id} → {freq_key}")
                                    break
                    
                    # Extract exon structure
                    exons = get_exon_structure(row['block_sizes'], row['block_starts'], row['start'])
                    
                    # Default values for frequency data
                    raw_freq = None
                    norm_freq = None
                    read_count = None
                    length = None
                    
                    # If we found frequency data, extract the values
                    if freq_info:
                        raw_freq = freq_info.get('Raw_Frequency')
                        norm_freq = freq_info.get('Normalized_Frequency')
                        read_count = freq_info.get('Read_Count')
                        length = freq_info.get('Length')
                    
                    if exons:
                        
                        # Map exons to their symbolic names
                        exons_present = map_exons_to_names(exons, chrom, tolerance_percentage=0.1)

                        transcripts.append({
                            'sample': sample_name,
                            'source_file': filename,
                            'chromosome': chrom,
                            'transcript': transcript_name,  # This will be PB ID for CP files if mapped
                            'original_transcript': original_name,
                            'cp_id': cp_id,  # Store the CP ID explicitly
                            'pb_id': pb_id_for_freq if pb_id_for_freq != transcript_name else None,
                            'matched_freq_name': matched_name,
                            'g_number': g_number,  # Add G number if we found one
                            'start': row['start'],
                            'end': row['end'],
                            'strand': row['strand'],
                            'exons': exons,
                            'exon_count': len(exons),
                            'exons_present': exons_present,
                            'Raw_Frequency': raw_freq,
                            'Normalized_Frequency': norm_freq,
                            'Read_Count': read_count,
                            'Length': length,
                            'is_older_format': is_older_format_bed,  # Flag older-format transcripts
                            'ancestry': ancestry,         # Add ancestry
                            'haplotype': haplotype,
                            'Cell_Type': Cell_Type,
                            'Derived_Cell_Type': Derived_Cell_Type,
                            'ancestry_haplotype': ancestry_haplotype,
                        })
        
        print(f"Matched chromosomes: {', '.join(matched_chroms)}")
        print(f"Found {len(transcripts)} transcripts with exon structure")
        
        # Count how many transcripts have frequency data
        freq_count = sum(1 for t in transcripts if t.get('Raw_Frequency') is not None)
        if len(transcripts) > 0:
            print(f"Found frequency data for {freq_count}/{len(transcripts)} transcripts ({freq_count/len(transcripts)*100:.1f}%)")
            
            # For older-format files, count how many have G numbers
            if is_older_format_bed:
                g_count = sum(1 for t in transcripts if t.get('g_number') is not None)
                print(f"Found G number mappings for {g_count}/{len(transcripts)} older-format transcripts ({g_count/len(transcripts)*100:.1f}%)")
                
                # If older-format file has 0 frequency matches, try direct matching as a last resort
                if freq_count == 0 and frequency_data:
                    print(f"Attempting direct older-format transcript matching for {sample_name}")
                    matches_found = 0
                    
                    for i, transcript_data in enumerate(transcripts):
                        current_transcript = transcript_data['transcript']
                        
                        # Try each frequency key - checking for matching PB.X.Y pattern
                        for freq_key in frequency_data.keys():
                            # Compare only the relevant parts (PB.1.X)
                            current_parts = current_transcript.split('.')
                            freq_parts = freq_key.split('.')
                            
                            if (len(current_parts) >= 3 and len(freq_parts) >= 3 and
                                current_parts[0] == freq_parts[0] and 
                                current_parts[1] == freq_parts[1]):
                                
                                # Match found! Use this frequency data
                                freq_info = frequency_data[freq_key]
                                transcripts[i]['Raw_Frequency'] = freq_info.get('Raw_Frequency')
                                transcripts[i]['Normalized_Frequency'] = freq_info.get('Normalized_Frequency')
                                transcripts[i]['Read_Count'] = freq_info.get('Read_Count')
                                transcripts[i]['Length'] = freq_info.get('Length')
                                transcripts[i]['matched_freq_name'] = freq_key
                                matches_found += 1
                                break
                    
                    # Report results after direct matching
                    if matches_found > 0:
                        print(f"  Direct matching found {matches_found}/{len(transcripts)} matches!")
                        # Recalculate frequency counts
                        freq_count = sum(1 for t in transcripts if t.get('Raw_Frequency') is not None)
                        print(f"  Updated: Found frequency data for {freq_count}/{len(transcripts)} transcripts ({freq_count/len(transcripts)*100:.1f}%)")
        else:
            print("No transcripts with exon structure found")
        
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        import traceback
        traceback.print_exc()
    
    return transcripts

# Main processing function
def main():
    print("Starting transcript exon analysis...")
    
    all_transcripts = []
    all_files = []
    
    # Gather all BED files from the directories
    for bed_dir in bed_dirs:
        if os.path.exists(bed_dir):
            bed_files = glob.glob(os.path.join(bed_dir, "*.bed")) + \
                        glob.glob(os.path.join(bed_dir, "*.bed12"))
            print(f"Found {len(bed_files)} BED files in {bed_dir}")
            all_files.extend(bed_files)
    
    print(f"Total BED files found: {len(all_files)}")
    
    # Counter for frequency data
    transcript_freq_found = 0
    total_transcripts = 0
        
    # Process each BED file
    for i, file_path in enumerate(all_files):
        if i % 20 == 0 or i == len(all_files) - 1:
            print(f"Processing file {i+1}/{len(all_files)}: {os.path.basename(file_path)}")
        transcripts = process_bed_file(file_path)
        
        # Count frequency data
        for transcript in transcripts:
            total_transcripts += 1
            if transcript.get('Raw_Frequency') is not None:
                transcript_freq_found += 1
                
        all_transcripts.extend(transcripts)
    
    print(f"\nProcessed {len(all_files)} BED files and found {len(all_transcripts)} transcripts")
    if total_transcripts > 0:
        print(f"Found frequency data for {transcript_freq_found} out of {total_transcripts} transcripts ({transcript_freq_found/total_transcripts*100:.2f}% coverage)")
    else:
        print("No transcripts found in the BED files.")
    
    # Convert to DataFrame for easier analysis
    # Convert to DataFrame for easier analysis
    df = pd.DataFrame(all_transcripts)
    if df.empty:
        print("No transcripts found. Exiting.")
        return

    # Check ancestry and haplotype distribution
    if 'ancestry' in df.columns:
        ancestry_counts = df['ancestry'].value_counts()
        print("\n--- Ancestry Distribution ---")
        print(ancestry_counts)
        
        if 'haplotype' in df.columns:
            haplotype_counts = df['haplotype'].value_counts()
            print("\n--- Haplotype Distribution ---")
            print(haplotype_counts)
            
            # Count ancestry-haplotype combinations
            ancestry_haplotype_counts = df.groupby(['ancestry', 'haplotype']).size()
            print("\n--- Ancestry-Haplotype Combinations ---")
            print(ancestry_haplotype_counts)
        
    # Create an intermediary output file to debug frequency extraction
    raw_transcript_data = os.path.join(output_dir, "all_transcripts_raw_detailed.csv")
    # Convert exons to string for CSV export
    df_detailed = df.copy()
    if 'exons' in df_detailed.columns:
        df_detailed['exons'] = df_detailed['exons'].apply(lambda x: '; '.join([f"{start}-{end}" for start, end in x]) if isinstance(x, list) else x)
    df_detailed.to_csv(raw_transcript_data, index=False)
    print(f"Saved detailed raw transcript data to {raw_transcript_data} for debugging")
        
    # Ensure frequency columns have proper numeric data
    frequency_columns = ['Raw_Frequency', 'Normalized_Frequency', 'Read_Count', 'Length']
    for col in frequency_columns:
        if col in df.columns:
            # Convert to numeric, handling errors
            df[col] = pd.to_numeric(df[col], errors='coerce')
            
            # Count non-null values
            non_null_count = df[col].notna().sum()
            print(f"Column '{col}' has {non_null_count} non-null values out of {len(df)} ({non_null_count/len(df)*100:.2f}%)")

    # --- Sample-level statistics ---
    print("\n--- Sample-Level Statistics ---")
    
    # Count transcripts per sample
    sample_counts = df.groupby('sample').size()
    print(f"Found transcripts in {len(sample_counts)} unique samples")
    print(f"Average transcripts per sample: {sample_counts.mean():.1f}")
    print(f"Samples with most transcripts: {sample_counts.nlargest(3).to_dict()}")
    print(f"Samples with fewest transcripts: {sample_counts.nsmallest(3).to_dict()}")
    
    # Count frequency data by sample
    if 'Raw_Frequency' in df.columns:
        freq_by_sample = df.groupby('sample')['Raw_Frequency'].apply(lambda x: x.notna().sum())
        samples_with_freq = sum(freq_by_sample > 0)
        print(f"Samples with any frequency data: {samples_with_freq}/{len(sample_counts)} ({samples_with_freq/len(sample_counts)*100:.1f}%)")
        
        # Show samples with best frequency coverage
        coverage_by_sample = df.groupby('sample').apply(lambda x: x['Raw_Frequency'].notna().sum() / len(x) * 100)
        print(f"Samples with highest frequency coverage: {coverage_by_sample.nlargest(5).to_dict()}")
        
        # Show samples with worst frequency coverage (but some data)
        non_zero_coverage = coverage_by_sample[coverage_by_sample > 0]
        if len(non_zero_coverage) > 0:
            print(f"Samples with lowest non-zero frequency coverage: {non_zero_coverage.nsmallest(5).to_dict()}")

    # --- Build grouping dictionaries based on identical exon structure ---
    print("\n--- Building Exon Structure Groups ---")
    transcript_to_exon_structure = {}
    transcript_to_samples = defaultdict(set)
    transcript_to_chromosomes = defaultdict(set)
    
    for _, row in df.iterrows():
        transcript_name = row['transcript']
        sample = row['sample']
        chromosome = row['chromosome']
        
        # Skip if exons is not a valid list
        if not isinstance(row['exons'], list):
            continue
            
        # Use frozenset of exons for hashability
        exon_structure = frozenset((start, end) for start, end in row['exons'])
        if transcript_name not in transcript_to_exon_structure:
            transcript_to_exon_structure[transcript_name] = exon_structure
        transcript_to_samples[transcript_name].add(sample)
        transcript_to_chromosomes[transcript_name].add(chromosome)
    
    # Map exon structures to transcript names
    exon_structure_to_transcripts = defaultdict(list)
    for transcript, exon_structure in transcript_to_exon_structure.items():
        exon_structure_to_transcripts[exon_structure].append(transcript)
    
    print(f"Found {len(transcript_to_exon_structure)} unique transcripts")
    print(f"Found {len(exon_structure_to_transcripts)} unique exon structures")
    
    # Count transcripts by prefix
    prefix_counts = {
        'PB': sum(1 for t in transcript_to_exon_structure if str(t).startswith('PB')),
        'CP': sum(1 for t in transcript_to_exon_structure if str(t).startswith('CP')),
        'G': sum(1 for t in transcript_to_exon_structure if str(t).startswith('G')),
        'other': sum(1 for t in transcript_to_exon_structure if not any(str(t).startswith(p) for p in ['PB', 'CP', 'G']))
    }
    print(f"Transcript prefix distribution: {prefix_counts}")
    
    # --- Add CP_transcript and G_transcript columns ---
    print("\n--- Adding Transcript Cross-References ---")
    def get_identical_cols(transcript):
        exon_structure = transcript_to_exon_structure.get(transcript)
        if not exon_structure:
            return "None", "None"
        # Get all transcripts sharing this exon structure (excluding the current one)
        identical_transcripts = [t for t in exon_structure_to_transcripts[exon_structure] if t != transcript]
        cp = [t for t in identical_transcripts if str(t).startswith("CP")]
        g = [t for t in identical_transcripts if str(t).startswith("G")]
        cp_str = ", ".join(cp) if cp else "None"
        g_str = ", ".join(g) if g else "None"
        return cp_str, g_str
    
    cp_transcripts = []
    g_transcripts = []
    for t in df['transcript']:
        cp_str, g_str = get_identical_cols(t)
        cp_transcripts.append(cp_str)
        g_transcripts.append(g_str)
    df['CP_transcript'] = cp_transcripts
    df['G_transcript'] = g_transcripts

    # --- Add a harmonized_transcript_name column ---
    print("\n--- Creating Harmonized Transcript Names ---")

    # First, add CP_transcript and G_transcript columns as before
    print("\n--- Adding Transcript Cross-References ---")
    def get_identical_cols(transcript):
        exon_structure = transcript_to_exon_structure.get(transcript)
        if not exon_structure:
            return "None", "None"
        # Get all transcripts sharing this exon structure (excluding the current one)
        identical_transcripts = [t for t in exon_structure_to_transcripts[exon_structure] if t != transcript]
        cp = [t for t in identical_transcripts if str(t).startswith("CP")]
        g = [t for t in identical_transcripts if str(t).startswith("G")]
        cp_str = ", ".join(cp) if cp else "None"
        g_str = ", ".join(g) if g else "None"
        return cp_str, g_str

    cp_transcripts = []
    g_transcripts = []
    for t in df['transcript']:
        cp_str, g_str = get_identical_cols(t)
        cp_transcripts.append(cp_str)
        g_transcripts.append(g_str)
    df['CP_transcript'] = cp_transcripts
    df['G_transcript'] = g_transcripts

    # Next, create the Arbitrary_Isoform column
    print("\n--- Creating Arbitrary Isoform Names ---")
    # Function to generate alphabetical labels (A, B, C, ..., Z, AA, AB, ...)
    def generate_labels():
        # Single letters first
        for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
            yield c
        
        # Then double letters
        for c1 in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
            for c2 in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                yield c1 + c2
        
        # Triple letters if needed
        for c1 in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
            for c2 in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                for c3 in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                    yield c1 + c2 + c3

    # Create a dictionary to store exon structure to arbitrary name mapping
    exon_structure_to_arbitrary_name = {}
    # Create a label generator
    label_gen = generate_labels()

    # Add the new column
    df['Arbitrary_Isoform'] = ""


    def are_similar_exons(exons1, exons2, tolerance_percentage=0.1):
        """
        Check if two exon structures are similar within a percentage-based tolerance.
        Each exon's start and end can vary by ±tolerance_percentage of the exon length.
        """
        if len(exons1) != len(exons2):
            return False
        
        # Sort both exon lists by start position
        sorted_exons1 = sorted(exons1, key=lambda x: x[0])
        sorted_exons2 = sorted(exons2, key=lambda x: x[0])
        
        # Compare each exon pair
        for (start1, end1), (start2, end2) in zip(sorted_exons1, sorted_exons2):
            # Calculate exon lengths
            length1 = end1 - start1
            length2 = end2 - start2
            
            # Use average length to calculate tolerance
            avg_length = (length1 + length2) / 2
            tolerance = int(avg_length * tolerance_percentage)
            
            # Check if start positions are within tolerance
            if abs(start1 - start2) > tolerance:
                return False
            # Check if end positions are within tolerance
            if abs(end1 - end2) > tolerance:
                return False
        
        return True

    # Store exon structures for each arbitrary name
    arbitrary_name_exon_structures = {}

    # Process all rows for flexible exon structure comparison
    for idx, row in df.iterrows():
        # Skip if exons is not a valid list
        if not isinstance(row['exons'], list) or len(row['exons']) == 0:
            continue
        
        # Get current exon structure
        current_exons = [(start, end) for start, end in row['exons']]
        
        # Check if this structure is similar to any existing one
        found_match = False
        for name, existing_structures in arbitrary_name_exon_structures.items():
            for existing_exons in existing_structures:
                if are_similar_exons(current_exons, existing_exons, tolerance_percentage=0.1):
                    # Found a similar structure, use its arbitrary name
                    df.at[idx, 'Arbitrary_Isoform'] = name
                    # Add this structure to the list for this name (for future comparisons)
                    arbitrary_name_exon_structures[name].append(current_exons)
                    found_match = True
                    break
            if found_match:
                break
        
        # If no similar structure found, create a new arbitrary name
        if not found_match:
            new_name = next(label_gen)
            df.at[idx, 'Arbitrary_Isoform'] = new_name
            # Initialize the list of structures for this new name
            arbitrary_name_exon_structures[new_name] = [current_exons]

    # Count unique arbitrary isoforms
    arbitrary_isoform_count = df['Arbitrary_Isoform'].nunique()
    print(f"Created {arbitrary_isoform_count} unique arbitrary isoform names")

    # Now create the harmonized_transcript_name with Arbitrary_Isoform as fallback
    print("\n--- Creating Harmonized Transcript Names ---")
    def get_harmonized(row):
        # First check if CP_transcript has a value
        if row['CP_transcript'] != "None":
            # Use first CP transcript if multiple are listed
            return row['CP_transcript'].split(", ")[0]
        # Next check if G_transcript has a value
        elif row['G_transcript'] != "None":
            # Use first G transcript if multiple are listed
            return row['G_transcript'].split(", ")[0]
        # Otherwise use the Arbitrary_Isoform
        else:
            return row['Arbitrary_Isoform']

    df['harmonized_transcript_name'] = df.apply(get_harmonized, axis=1)
    
    
    # --- Map to Isoforms based on exact exon structure ---
    print("\n--- Mapping Transcripts to Isoforms Based on Exact Exon Structure ---")

    # Define paths to the transcript exon mapping files
    h1_transcript_exons_file = H1_TRANSCRIPT_EXONS_FILE
    h2_transcript_exons_file = H2_TRANSCRIPT_EXONS_FILE

    # Create a dictionary to map exon structures to isoform names
    exon_structure_to_isoform = {}

    # Function to process the exon structure mapping files
    def load_exon_to_isoform_mapping(file_path):
        try:
            mapping_df = pd.read_csv(file_path)
            count = 0
            
            for _, row in mapping_df.iterrows():
                # Get the exon structure from the file
                if 'exons' in row and 'new_transcript_name' in row:
                    try:
                        # Parse the exon string into a list of tuples
                        exon_str = row['exons']
                        if isinstance(exon_str, str):
                            # Extract the tuples from the string format [(x, y), (a, b), ...]
                            exon_tuples = []
                            import re
                            matches = re.findall(r'\((\d+),\s*(\d+)\)', exon_str)
                            for match in matches:
                                start, end = int(match[0]), int(match[1])
                                exon_tuples.append((start, end))
                            
                            # Sort the exons by start position for consistent comparison
                            exon_tuples.sort(key=lambda x: x[0])
                            
                            # Create a frozenset for hashable lookup
                            exon_structure = frozenset(exon_tuples)
                            
                            # Map this structure to the isoform name
                            isoform_name = row['new_transcript_name']
                            exon_structure_to_isoform[exon_structure] = isoform_name
                            count += 1
                    except Exception as e:
                        print(f"Error parsing exon structure: {e}")
                        continue
            
            print(f"Loaded {count} exon structure to isoform mappings from {file_path}")
            return count
        except Exception as e:
            print(f"Error loading exon structure mapping file {file_path}: {e}")
            return 0

    # Load both mapping files
    h1_count = load_exon_to_isoform_mapping(h1_transcript_exons_file)
    h2_count = load_exon_to_isoform_mapping(h2_transcript_exons_file)
    print(f"Total exon structure to isoform mappings: {len(exon_structure_to_isoform)}")


    def map_to_isoform_by_exons(row, tolerance_percentage=0.1):
        # Skip if exons is not a valid list
        if not isinstance(row['exons'], list) or len(row['exons']) == 0:
            return ""
        
        # Get the sorted exon structure for this transcript
        sorted_exons = sorted([(start, end) for start, end in row['exons']], key=lambda x: x[0])
        
        # Check each known isoform exon structure
        for known_structure, isoform_name in exon_structure_to_isoform.items():
            # Convert frozenset to sorted list for ordered comparison
            known_exons = sorted(list(known_structure), key=lambda x: x[0])
            
            # Skip if different number of exons
            if len(sorted_exons) != len(known_exons):
                continue
            
            # Check if structures match within tolerance
            is_match = True
            for (start1, end1), (start2, end2) in zip(sorted_exons, known_exons):
                # Calculate exon lengths and tolerance
                length1 = end1 - start1
                length2 = end2 - start2
                avg_length = (length1 + length2) / 2
                tolerance = int(avg_length * tolerance_percentage)
                
                # Check if start and end positions are within tolerance
                if abs(start1 - start2) > tolerance or abs(end1 - end2) > tolerance:
                    is_match = False
                    break
            
            if is_match:
                return isoform_name
                
        # No match found within tolerance
        return ""

    # Apply the mapping function to each row
    #df['Isoform'] = df.apply(map_to_isoform_by_exons, axis=1)
    # Apply the mapping function to each row
   # df['Isoform'] = df.apply(lambda row: map_to_isoform_by_exons(row, tolerance=20), axis=1)
    df['Isoform'] = df.apply(lambda row: map_to_isoform_by_exons(row, tolerance_percentage=0.1), axis=1)
    # Check isoform mapping success
    mapped_count = sum(1 for iso in df['Isoform'] if iso)  # Count non-empty strings
    if len(df) > 0:
        print(f"Successfully mapped {mapped_count}/{len(df)} transcripts to isoforms based on exact exon structure ({mapped_count/len(df)*100:.2f}%)")
    
    # --- Merge frequency data with exon structure ---
    print("\n--- Merging Frequency Data with Exon Structure ---")
    # Calculate average frequency by harmonized transcript name
    if 'Raw_Frequency' in df.columns and df['Raw_Frequency'].notna().sum() > 0:
        avg_freq = df.groupby('harmonized_transcript_name')['Raw_Frequency'].mean().reset_index()
        avg_freq.columns = ['harmonized_transcript_name', 'Average_Raw_Frequency']
        df = pd.merge(df, avg_freq, on='harmonized_transcript_name', how='left')
        
        # Fill in Average_Raw_Frequency where Raw_Frequency is missing
        df['Final_Frequency'] = df['Raw_Frequency'].fillna(df['Average_Raw_Frequency'])
        
        # Count how many transcripts have frequency data now
        final_freq_count = df['Final_Frequency'].notna().sum()
        print(f"After merging: {final_freq_count}/{len(df)} transcripts have frequency data ({final_freq_count/len(df)*100:.2f}%)")
    
    # --- Create harmonized_isoform column (combining Isoform and Arbitrary_Isoform) ---
    print("\n--- Creating Harmonized Isoform Column ---")

    def get_harmonized_isoform(row):
        # First check if Isoform has a value
        if row['Isoform'] and not pd.isna(row['Isoform']):
            return row['Isoform']
        # Otherwise use the Arbitrary_Isoform
        else:
            return row['Arbitrary_Isoform']

    df['harmonized_isoform'] = df.apply(get_harmonized_isoform, axis=1)

    # Check how many rows use each source
    isoform_count = sum(1 for iso in df['Isoform'] if iso and not pd.isna(iso))
    arbitrary_count = len(df) - isoform_count
    print(f"harmonized_isoform source breakdown:")
    print(f"  - From Isoform: {isoform_count}/{len(df)} ({isoform_count/len(df)*100:.2f}%)")
    print(f"  - From Arbitrary_Isoform: {arbitrary_count}/{len(df)} ({arbitrary_count/len(df)*100:.2f}%)")
    
    # --- Add Isoform_type column based on presence of specific exons ---
    print("\n--- Adding Isoform Type Classification ---")

    # Apply the function to determine isoform type
    df['Isoform_type'] = df['exons_present'].apply(determine_isoform_type)

    # Print summary of isoform types
    isoform_type_counts = df['Isoform_type'].value_counts()
    print(f"Isoform type distribution:")
    for isoform_type, count in isoform_type_counts.items():
        print(f"  {isoform_type}: {count} transcripts ({count/len(df)*100:.2f}%)")
    
    
    # --- Add read count normalization columns ---
    print("\n--- Adding Read Count Normalization Columns ---")

    # Create a dictionary to store total read count per sample
    sample_to_total_read_count = {}

    # Calculate total read count for each sample
    for sample in df['sample'].unique():
        sample_reads = df[df['sample'] == sample]['Read_Count'].sum()
        sample_to_total_read_count[sample] = sample_reads

    # Add the total read count column
    df['total_read_count_sample'] = df['sample'].map(sample_to_total_read_count)

    # Calculate normalized read count
    df['read_count_over_total_sample_read_count'] = df.apply(
        lambda row: row['Read_Count'] / row['total_read_count_sample'] 
        if pd.notna(row['Read_Count']) and row['total_read_count_sample'] > 0 
        else None, axis=1
    )

    # Print summary of the new columns
    print(f"Added total read count per sample and normalized read count columns")
    print(f"Non-null normalized read counts: {df['read_count_over_total_sample_read_count'].notna().sum()}/{len(df)}")
    
    # --- Save the processed transcript data ---
    # Convert the exons list to a string for CSV export
    df_for_csv = df.copy()
    if 'exons' in df_for_csv.columns:
        df_for_csv['exons'] = df_for_csv['exons'].apply(lambda x: '; '.join([f"{start}-{end}" for start, end in x]) if isinstance(x, list) else x)
    
    # Save to CSV
    raw_output = os.path.join(output_dir, "all_transcripts_processed.csv")
    df_for_csv.to_csv(raw_output, index=False)
    print(f"Processed transcript data saved to {raw_output}")
    
    # --- Create isoform-specific summary ---
    if 'Isoform' in df.columns and df['Isoform'].notna().sum() > 0:
        print("\n--- Creating Isoform Summary ---")
        # Keep only rows with valid isoforms
        df_isoforms = df[df['Isoform'].notna() & (df['Isoform'] != "")]
        
        if not df_isoforms.empty:
            # Calculate summary statistics by isoform
            isoform_summary = df_isoforms.groupby('Isoform').agg({
                'sample': 'nunique',
                'Final_Frequency': ['mean', 'median', 'std', 'count'],
                'exon_count': ['mean', 'median', 'min', 'max']
            }).reset_index()
            
            # Flatten the multi-index columns
            isoform_summary.columns = ['_'.join(col).strip('_') for col in isoform_summary.columns.values]
            
            # Save isoform summary
            isoform_output = os.path.join(output_dir, "isoform_summary.csv")
            isoform_summary.to_csv(isoform_output, index=False)
            print(f"Isoform summary saved to {isoform_output}")
            
            # Create a visualization of isoform frequencies
            try:
                print("Creating isoform frequency visualization...")
                import matplotlib.pyplot as plt
                
                # Get top 15 isoforms by frequency
                top_isoforms = isoform_summary.nlargest(15, 'Final_Frequency_mean')
                
                plt.figure(figsize=(12, 8))
                plt.bar(top_isoforms['Isoform'], top_isoforms['Final_Frequency_mean'])
                plt.title('Average Frequency of Top 15 Isoforms')
                plt.xlabel('Isoform')
                plt.ylabel('Average Frequency')
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()
                
                # Save the plot
                plot_output = os.path.join(output_dir, "top_isoforms_frequency.png")
                plt.savefig(plot_output, dpi=300)
                plt.close()
                print(f"Isoform frequency visualization saved to {plot_output}")
            except Exception as e:
                print(f"Error creating visualization: {e}")

    # Validate using the same leniency criteria
    print("\n--- Validating Arbitrary Isoform Names with ±10 base leniency ---")

    validation_issues = []

    for idx, row in df.iterrows():
        arbitrary_name = row['Arbitrary_Isoform']
        
        # Skip if no arbitrary name or exons not a list
        if not arbitrary_name or not isinstance(row['exons'], list):
            continue
        
        # Get current exon structure
        current_exons = [(start, end) for start, end in row['exons']]
        
        # Check if this matches any of the stored structures for this name
        if arbitrary_name in arbitrary_name_exon_structures:
            # Check against all stored structures for this name
#            is_valid = any(are_similar_exons(current_exons, stored_exons, tolerance=20) 
#                        for stored_exons in arbitrary_name_exon_structures[arbitrary_name])
            
            is_valid = any(are_similar_exons(current_exons, stored_exons, tolerance_percentage=0.1) 
                            for stored_exons in arbitrary_name_exon_structures[arbitrary_name])
            
            if not is_valid:
                # Found a mismatch
                validation_issues.append({
                    'arbitrary_name': arbitrary_name,
                    'row_index': idx,
                    'transcript': row['transcript'],
                    'sample': row['sample'],
                    'exons': current_exons,
                    'stored_exons_count': len(arbitrary_name_exon_structures[arbitrary_name])
                })

    # Print validation results
    if validation_issues:
        print(f"WARNING: Found {len(validation_issues)} instances where rows with the same Arbitrary_Isoform have different exon structures:")
        for issue in validation_issues[:10]:  # Print first 10 issues to avoid excessive output
            print(f"  Arbitrary Isoform: {issue['arbitrary_name']}")
            print(f"  Transcript: {issue['transcript']} (Sample: {issue['sample']})")
            print(f"  Expected exons: {issue['expected_exons']}")
            print(f"  Actual exons: {issue['actual_exons']}")
            print("  " + "-" * 50)
        if len(validation_issues) > 10:
            print(f"  ... and {len(validation_issues) - 10} more issues.")
        
        # Save validation issues to file for further investigation
        validation_output = os.path.join(output_dir, "arbitrary_isoform_validation_issues.csv")
        validation_df = pd.DataFrame(validation_issues)
        
        # Convert frozensets to strings for CSV export
        validation_df['expected_exons'] = validation_df['expected_exons'].apply(lambda x: '; '.join([f"{start}-{end}" for start, end in x]))
        validation_df['actual_exons'] = validation_df['actual_exons'].apply(lambda x: '; '.join([f"{start}-{end}" for start, end in x]))
        
        validation_df.to_csv(validation_output, index=False)
        print(f"Saved validation issues to {validation_output}")
    else:
        print("Validation successful! All rows with the same Arbitrary_Isoform have identical exon structures.")
    
    # --- Check Isoform Consistency ---
    print("\n--- Checking Isoform to Exon Structure Consistency ---")

    # Create a dictionary to store unique exon structures for each isoform
    isoform_to_exon_structures = defaultdict(set)

    # Track isoforms and their exon structures
    for _, row in df.iterrows():
        isoform = row.get('Isoform', '')
        
        # Skip empty isoforms
        if not isoform:
            continue
        
        # Skip if exons is not a valid list
        if not isinstance(row['exons'], list) or len(row['exons']) == 0:
            continue
        
        # Create a consistent representation of the exon structure
        sorted_exons = tuple(sorted([(start, end) for start, end in row['exons']], key=lambda x: x[0]))
        # Add this exon structure to the set for this isoform
        isoform_to_exon_structures[isoform].add(sorted_exons)

    # Report results
    print(f"Found {len(isoform_to_exon_structures)} isoforms with exon structure information")
    inconsistent_isoforms = []

    for isoform, exon_structures in isoform_to_exon_structures.items():
        structure_count = len(exon_structures)
        if structure_count > 1:
            inconsistent_isoforms.append((isoform, structure_count))
            
    # Print results
    if inconsistent_isoforms:
        print(f"WARNING: Found {len(inconsistent_isoforms)} isoforms with multiple exon structures:")
        for isoform, count in sorted(inconsistent_isoforms, key=lambda x: x[1], reverse=True):
            print(f"  Isoform '{isoform}' has {count} different exon structures")
        
        # Save detailed report to file
        inconsistency_output = os.path.join(output_dir, "isoform_exon_inconsistencies.csv")
        
        # Prepare data for CSV
        inconsistency_data = []
        for isoform, structures in isoform_to_exon_structures.items():
            if len(structures) > 1:
                for i, structure in enumerate(structures):
                    # Convert the structure to a readable string
                    structure_str = "; ".join([f"{start}-{end}" for start, end in structure])
                    inconsistency_data.append({
                        'Isoform': isoform,
                        'Structure_Number': i+1,
                        'Total_Structures': len(structures),
                        'Exon_Structure': structure_str
                    })
        
        # Create and save DataFrame
        if inconsistency_data:
            pd.DataFrame(inconsistency_data).to_csv(inconsistency_output, index=False)
            print(f"Saved detailed inconsistency report to {inconsistency_output}")
    else:
        print("All isoforms have consistent exon structures (1 structure per isoform)")

    # Create a summary table
    summary_table = []
    for isoform, structures in isoform_to_exon_structures.items():
        structure_count = len(structures)
        summary_table.append({
            'Isoform': isoform,
            'Unique_Exon_Structures': structure_count,
            'Is_Consistent': 'Yes' if structure_count == 1 else 'No'
        })

    # Save summary table
    summary_output = os.path.join(output_dir, "isoform_exon_structure_summary.csv")
    pd.DataFrame(summary_table).to_csv(summary_output, index=False)
    print(f"Saved isoform exon structure summary to {summary_output}")
    
        
    # --- Check harmonized_transcript_name Consistency ---
    print("\n--- Checking harmonized_transcript_name to Exon Structure Consistency ---")

    # Create a dictionary to store unique exon structures for each harmonized transcript name
    harmonized_to_exon_structures = defaultdict(set)

    # Track harmonized names and their exon structures
    for _, row in df.iterrows():
        harmonized_name = row.get('harmonized_transcript_name', '')
        
        # Skip empty harmonized names
        if not harmonized_name:
            continue
        
        # Skip if exons is not a valid list
        if not isinstance(row['exons'], list) or len(row['exons']) == 0:
            continue
        
        # Create a consistent representation of the exon structure
        sorted_exons = tuple(sorted([(start, end) for start, end in row['exons']], key=lambda x: x[0]))
        # Add this exon structure to the set for this harmonized name
        harmonized_to_exon_structures[harmonized_name].add(sorted_exons)

    # Create a summary of transcripts with multiple exon structures
    harmonized_structure_counts = {}
    for name, structures in harmonized_to_exon_structures.items():
        structure_count = len(structures)
        harmonized_structure_counts[name] = structure_count

    # Report results
    inconsistent_harmonized = {name: count for name, count in harmonized_structure_counts.items() if count > 1}
    print(f"Found {len(inconsistent_harmonized)} harmonized transcript names with multiple exon structures:")

    # Print the top cases with the most different structures
    if inconsistent_harmonized:
        sorted_inconsistent = sorted(inconsistent_harmonized.items(), key=lambda x: x[1], reverse=True)
        for name, count in sorted_inconsistent[:10]:  # Show top 10
            print(f"  {name}: {count} different exon structures")
        
        # Save to CSV for further analysis
        inconsistent_output = os.path.join(output_dir, "harmonized_transcript_exon_inconsistencies.csv")
        inconsistent_df = pd.DataFrame([{"harmonized_transcript_name": name, "exon_structure_count": count} 
                                    for name, count in sorted_inconsistent])
        inconsistent_df.to_csv(inconsistent_output, index=False)
        print(f"Saved harmonized transcript name inconsistencies to {inconsistent_output}")
    else:
        print("All harmonized transcript names have consistent exon structures (1 structure per name)")


        
            
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()