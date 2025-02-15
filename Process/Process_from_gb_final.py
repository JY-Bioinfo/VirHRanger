import os
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def create_directory_structure(base_dir, base_name):
    """
    Creates the required directory structure:
    base_dir/base_name/
        ├── genbank/
        ├── genome/
        └── protein/
    """
    main_dir = os.path.join(base_dir, base_name)
    sub_dirs = ['genbank', 'genome', 'protein']
    
    for sub_dir in sub_dirs:
        path = os.path.join(main_dir, sub_dir)
        os.makedirs(path, exist_ok=True)
    
    return main_dir

def copy_genbank_file(source_file, dest_dir):
    """
    Copies the .gbff file to the genbank directory.
    """
    shutil.copy(source_file, dest_dir)

def extract_genome_fasta(records, output_path, header_suffix=""):
    """
    Extracts the genome sequences from multiple records and writes them to a single FASTA file.
    Each genome is a separate record in the FASTA.
    """
    genome_records = []
    for record in records:
        genome_record = SeqRecord(
            record.seq,
            id=f"{record.id}{header_suffix}",
            description=""
        )
        genome_records.append(genome_record)
    
    if genome_records:
        SeqIO.write(genome_records, output_path, "fasta")
        print(f"Genome FASTA written to {output_path}")
    else:
        print(f"No genome sequences found.")

def extract_protein_fasta(records, output_path, header_suffix=""):
    """
    Extracts all protein translations from multiple records and writes them to a single FASTA file.
    Each protein is a separate record in the FASTA.
    """
    protein_records = []
    for record in records:
        for feature in record.features:
            if feature.type == "CDS":
                protein_id = feature.qualifiers.get('protein_id', [None])[0]
                translation = feature.qualifiers.get('translation', [None])[0]
                
                if protein_id and translation:
                    protein_record = SeqRecord(
                        Seq(translation),
                        id=f"{protein_id}{header_suffix}",
                        description=""
                    )
                    protein_records.append(protein_record)
                else:
                    print(f"Warning: Missing protein_id or translation in {record.id}")
    
    if protein_records:
        SeqIO.write(protein_records, output_path, "fasta")
        print(f"Protein FASTA written to {output_path}")
    else:
        print(f"No protein translations found.")

def process_gbff_file(gbff_path, base_dir):
    """
    Processes a single .gbff file:
    - Creates directory structure
    - Copies the .gbff file
    - Extracts genome and protein FASTA files for all records
    """
    base_name = os.path.basename(gbff_path).replace('.gbff', '')
    main_dir = create_directory_structure(base_dir, base_name)
    
    # Copy the .gbff file to the genbank directory
    genbank_dir = os.path.join(main_dir, 'genbank')
    copy_genbank_file(gbff_path, genbank_dir)
    
    # Parse the .gbff file
    try:
        records = list(SeqIO.parse(gbff_path, "genbank"))
        if not records:
            print(f"No records found in {gbff_path}.")
            return
    except Exception as e:
        print(f"Error parsing {gbff_path}: {e}")
        return
    
    # Extract genome FASTA
    genome_fasta_path = os.path.join(main_dir, 'genome', f"{base_name}_genome.fasta")
    extract_genome_fasta(records, genome_fasta_path)
    
    # Extract protein FASTA
    protein_fasta_path = os.path.join(main_dir, 'protein', f"{base_name}_protein.fasta")
    extract_protein_fasta(records, protein_fasta_path)

def split_genbank(input_path, output_dir):
    """
    Splits a GenBank (.gbff) file into individual GenBank records.

    Parameters:
    - input_path: Path to the input .gbff file.
    - output_dir: Directory where individual .gbff files will be saved.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    try:
        records = SeqIO.parse(input_path, "genbank")
        count = 0
        for record in records:
            # Use record accession or another unique identifier for the filename
            filename = f"{record.id}.gbff"
            output_path = os.path.join(output_dir, filename)
            SeqIO.write(record, output_path, "genbank")
            count += 1
        print(f"  Successfully split GenBank file into {count} records.")
        return count > 0
    except Exception as e:
        print(f"  Error splitting GenBank file: {e}")
        return False

def split_fasta(input_path, output_dir, file_extension="fasta"):
    """
    Splits a FASTA file into individual FASTA records.

    Parameters:
    - input_path: Path to the input FASTA file.
    - output_dir: Directory where individual FASTA files will be saved.
    - file_extension: Extension for the output files (default is 'fasta').
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    try:
        records = SeqIO.parse(input_path, "fasta")
        count = 0
        for record in records:
            # Use record ID for the filename
            # Replace any characters that are invalid in filenames
            #safe_id = "".join([c if c.isalnum() or c in ('-', '_') else "_" for c in record.id])
            safe_id = record.id
            filename = f"{safe_id}.{file_extension}"
            output_path = os.path.join(output_dir, filename)
            SeqIO.write(record, output_path, "fasta")
            count += 1
        print(f"  Successfully split FASTA file into {count} records.")
        return count > 0
    except Exception as e:
        print(f"  Error splitting FASTA file: {e}")
        return False

def process_subfolder(parent_folder, subfolder_name):
    """
    Processes a specific subfolder (genbank, genome, or protein) within a parent folder.

    Parameters:
    - parent_folder: Path to the parent folder containing the subfolder.
    - subfolder_name: Name of the subfolder to process (genbank, genome, or protein).
    """
    subfolder_path = os.path.join(parent_folder, subfolder_name)
    if not os.path.isdir(subfolder_path):
        print(f"  Subfolder '{subfolder_name}' does not exist in '{parent_folder}'. Skipping.")
        return

    # Identify files based on subfolder type
    if subfolder_name == "genbank":
        file_extensions = ('.gbff', '.gbk')
        split_function = split_genbank
        output_extension = "gbff"
    elif subfolder_name == "genome":
        file_extensions = ('.fasta', '.fa', '.fna')
        split_function = split_fasta
        output_extension = "fasta"
    elif subfolder_name == "protein":
        file_extensions = ('.fasta', '.fa', '.faa')
        split_function = split_fasta
        output_extension = "fasta"
    else:
        print(f"  Unknown subfolder '{subfolder_name}'. Skipping.")
        return

    # Find relevant files in the subfolder
    files = [f for f in os.listdir(subfolder_path) if f.lower().endswith(file_extensions)]
    
    if not files:
        print(f"  No relevant files found in '{subfolder_path}'. Skipping.")
        return

    for file in files:
        input_file_path = os.path.join(subfolder_path, file)
        output_dir = subfolder_path
        print(f"  Processing file: {input_file_path}")
        success = split_function(input_file_path, output_dir) if subfolder_name != "genbank" else split_function(input_file_path, output_dir)
        if success:
            try:
                os.remove(input_file_path)
                print(f"    Removed original file: {input_file_path}")
            except Exception as e:
                print(f"    Error removing original file '{input_file_path}': {e}")

def process_all_folders(parent_directory):
    """
    Iterates over all subdirectories in the parent directory and processes each.

    Parameters:
    - parent_directory: Path to the parent directory containing the 100 folders.
    """
    if not os.path.isdir(parent_directory):
        # print(f"Error: The provided path '{parent_directory}' is not a directory or does not exist.")
        return

    subfolders = [os.path.join(parent_directory, d) for d in os.listdir(parent_directory) if os.path.isdir(os.path.join(parent_directory, d))]
    
    if not subfolders:
        # print(f"No subdirectories found in '{parent_directory}'.")
        return

    total_genbank = total_genome = total_protein = 0
    for folder in subfolders:
        folder_name = os.path.basename(folder)
        # print(f"Processing folder: {folder_name}")
        
        # Process 'genbank' subfolder
        # print(f"  Processing 'genbank' subfolder...")
        genbank_initial = total_genbank
        process_subfolder(folder, "genbank")
        # No return value captured; could be enhanced to aggregate counts

        # Process 'genome' subfolder
        # print(f"  Processing 'genome' subfolder...")
        genome_initial = total_genome
        process_subfolder(folder, "genome")
        
        # Process 'protein' subfolder
        # print(f"  Processing 'protein' subfolder...")
        protein_initial = total_protein
        process_subfolder(folder, "protein")
        
        # print(f"Completed processing folder: {folder_name}\n")

    # print("All folders have been processed.")
    
    
    
    
def main():
    import argparse

    parser = argparse.ArgumentParser(description="Process GenBank .gbff files to extract genome and protein FASTA sequences.")
    parser.add_argument('-i', '--input_dir', required=True, help="Path to the directory containing .gbff files.")
    parser.add_argument('-o', '--output_dir', required=True, help="Path to the directory where output folders will be created.")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir

    # Validate input directory
    if not os.path.isdir(input_dir):
        print(f"Error: Input directory '{input_dir}' does not exist.")
        return

    # Create output directory if it doesn't exist
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Get list of .gbff files
    gbff_files = [f for f in os.listdir(input_dir) if f.endswith('.gbff')]

    if not gbff_files:
        print(f"No .gbff files found in '{input_dir}'.")
        return

    print(f"Found {len(gbff_files)} .gbff file(s). Starting processing...")

    for gbff_file in gbff_files:
        gbff_path = os.path.join(input_dir, gbff_file)
        print(f"\nProcessing file: {gbff_file}")
        process_gbff_file(gbff_path, output_dir)
    
    process_all_folders(output_dir)
        
    print("\nAll files have been processed successfully.")

if __name__ == "__main__":
    main()