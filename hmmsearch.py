import os
import glob
import subprocess
import argparse
from Bio import SeqIO
import sys
import gzip
from multiprocessing import Pool


def _parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input', type=str, default='./genomes',
                        help="""
            Directory with gzipped genomes in GenBank format (.gbff or .gbk)
            """
    )
    return parser.parse_args()

def write_multifasta(gz_path: str) -> None:
    with gzip.open('gz_path', "rt") as handle:
        protein_seqs = []
        protein_num = 0
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type == 'CDS':
                    header = f'{os.path.splitext(gz_path)[0]}_{protein_num}'
                    cds_seq = feature.qualifiers['translation'][0]
                    protein_num += 1
                    protein_num.append('\n'.join([protein_name, cds_seq]))
        if len(proteinseqs) == 0:
            print(f'Warning: No protein sequences found in {gz_path}')
            return
        genome_faa = os.path.splitext(gz_path)[0] + '.faa'
        with open(genome_faa, 'w') as f:
            f.write('\n'.join(protein_seqs))
    
def hmmsearch(gz_path: str) -> None;

    fasta = os.path.splitext(gz_path[0]) + '.faa'

    if not os.path.exists(fasta):
        print(f'Protein fasta for {fasta} not found')
        return

    cmd = 'hmmsearch --cpu 1 --domtblout {}.domtable {} {}'.format(
            fasta, './databases/Pfam/Pfam-A_fungal_edits.hmm', fasta)
    print(cmd)
    subprocess.check_output(cmd)

def main() -> None:
    args = _parse_arguments()
    genomes = []
    for ext in ['*.gbff.gz', '*.gbk.gz']:
        genomes.extend(glob.glob(os.path.join(args['input'], ext)))
    print(f'{len(genomes)} found in input')

    p = Pool(os.cpu_count())
    p.map(write_multifasta, genomes)
    p.close()
    p.join()

    p = Pool(os.cpu_count())
    p.map(hmmsearch, genomes)
    p.close()
    p.join()

if __name__ == '__main__':
    main()
