import pandas as pd
from Bio import SeqIO, SearchIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def upload_fasta(file_path):
    """Load the fasta file."""
    sequences = list(SeqIO.parse(file_path, "fasta"))
    return sequences

def run_blastp(sequence, sequence_id):
    """Run BLASTP for the given sequence."""
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence.seq)
    blast_record = NCBIXML.read(result_handle)
    return blast_record

def process_blast_results(blast_record, cutoff=80.0):
    """Process the BLASTP results and filter hits based on sequence identity."""
    hits = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            identity = (hsp.identities / hsp.align_length) * 100
            if identity < cutoff:
                hits.append({
                    'query': blast_record.query,
                    'hit_id': alignment.hit_id,
                    'hit_def': alignment.hit_def,
                    'length': alignment.length,
                    'e_value': hsp.expect,
                    'identity': identity,
                    'alignment_length': hsp.align_length,
                    'query_start': hsp.query_start,
                    'query_end': hsp.query_end,
                    'hit_start': hsp.sbjct_start,
                    'hit_end': hsp.sbjct_end,
                    'sequence': hsp.sbjct
                })
    return pd.DataFrame(hits).sort_values(by='identity', ascending=True).head(5)

def save_results(results, output_csv, output_fasta):
    """Save the results to CSV and FASTA files."""
    results.to_csv(output_csv, index=False)
    
    fasta_records = []
    for index, row in results.iterrows():
        record = SeqRecord(
            Seq(row['sequence']),
            id=row['hit_id'],
            description=row['hit_def']
        )
        fasta_records.append(record)
    
    SeqIO.write(fasta_records, output_fasta, "fasta")

# Define file paths
input_fasta = r"C:\Users\Henry\Downloads\pat_5_abg.fasta"
output_csv = r"C:\Users\Henry\Downloads\Pat_files5_blasp.csv"
output_fasta = r"C:\Users\Henry\Downloads\Pat_files5_blasp.FASTA"

# Load sequences
sequences = upload_fasta(input_fasta)

# Define cutoff
identity_cutoff = 80.0

# Process each sequence
all_results = pd.DataFrame()
for seq in sequences:
    print(f"Running BLASTP for sequence: {seq.id}")
    blast_record = run_blastp(seq, seq.id)
    results = process_blast_results(blast_record, identity_cutoff)
    all_results = pd.concat([all_results, results], ignore_index=True)

# Save results
save_results(all_results, output_csv, output_fasta)
print("Results saved to CSV and FASTA files.")
