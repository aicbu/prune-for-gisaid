#!/usr/bin/env python
from Bio import SeqIO
import pandas as pd
import argparse
import sys
import os

n_count=0
new_seqs=[]
new_ids=[]

parser = argparse.ArgumentParser(description='Count Ambigious bases and trim according to a cutoff')
parser.add_argument('fasta', type=str, help='input fasta')
parser.add_argument('amb', type=int, default='30', help='N% cutoff - default 30%')
args=parser.parse_args()

N_cutoff=args.amb
#Input fasta
input_fasta= os.path.join(args.fasta)
if not os.path.exists(args.fasta):
    sys.stderr.write(f"Error: can't find the fasta file\n")
    sys.exit(-1)
    

print("Sequences only have less than "+ str(100-N_cutoff) +"% coverage: \n")
for i in SeqIO.parse(open(input_fasta), 'fasta'):
    seq=i.seq
    n_count=round(seq.count('N')/29903*100,1)
    if n_count>=N_cutoff:
        print(i.id,"-",n_count)
    else:
        new_seqs.append(i)
        new_ids.append(i.id)
SeqIO.write(new_seqs,'output.fasta','fasta')
ids=pd.DataFrame(new_ids)
ids.to_csv('output.csv',index=False, header=False)
    
print("\nThese sequences will be removed and a new fasta and metadata files will be written..")
