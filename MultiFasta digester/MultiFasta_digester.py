from Bio import SeqIO
from pyteomics.parser import cleave, expasy_rules
from Bio.SeqUtils import molecular_weight
from collections import Counter
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    '-n',
    '--file_name',
    help=
    'Write the name of the fasta file (including the.fasta suffix)',
    type=str,
    action='append',
    nargs='*',
    required=True,
)

parser.add_argument(
    '-m',
    '--missed_cleavage',
    help=
    'How many missed cleavages are allowed',
    type=int,
    action='append',
    nargs='*',
    required=False,
)

args = parser.parse_args()  #store the arguments

seqs = dict()

for seq_record in SeqIO.parse(args.file_name[0][0], "fasta"):
    seqs[seq_record.id] = str(seq_record.seq)

digests = dict()

if args.missed_cleavage:
    for i in seqs:
        digests[i] = list(cleave(seqs[i],expasy_rules['trypsin'], args.missed_cleavage[0][0]))
else:
    for i in seqs:
        digests[i] = list(cleave(seqs[i],expasy_rules['trypsin']))

digests_list = list()

for i in digests:
    for j in digests[i]:
        if len(j) >= 7:
            digests_list.append(j)

c = Counter(digests_list)

for i in c.most_common(5):
    print(i[0] +' ---- '+str(i[1]) + ' ---- '+ str(molecular_weight(i[0],'protein')) )