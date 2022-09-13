from Bio import SeqIO
from pyteomics.parser import cleave, expasy_rules
from Bio.SeqUtils import molecular_weight
from collections import Counter
import argparse
import regex as re

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

parser.add_argument(
    '-c',
    '--compare_taxa',
    help=
    'The required taxa to be compared',
    type=str,
    action='append',
    nargs='*',
    required=False,
)


args = parser.parse_args()  #store the arguments
taxa = [i.split(',') for i in args.compare_taxa[0]]

seqs = dict()
for seq_record in SeqIO.parse(args.file_name[0][0], "fasta"):
    seqs[seq_record.description] = str(seq_record.seq)

def digester(sequence_dict):

    digests = dict()

    if args.missed_cleavage:
        for i in sequence_dict:
            digests[i] = list(cleave(sequence_dict[i],expasy_rules['trypsin'], args.missed_cleavage[0][0]))
    else:
        for i in sequence_dict:
            digests[i] = list(cleave(sequence_dict[i],expasy_rules['trypsin']))

    digests_list = list()

    for i in digests:
        for j in digests[i]:
            if (len(j) <= 7) or ('X' in j) or ('-' in j):
                continue
            else:
                digests_list.append(j)
    return digests_list


group_1 = dict()
group_2 = dict()
excl_group_1 = dict()
excl_group_2 = dict()

for key in seqs:
    if taxa[0][0] in key: 
        group_1[key] = seqs[key]

for key in seqs:
    if taxa[0][1] in key: 
        group_2[key] = seqs[key]

for key in seqs:
    if not taxa[0][0] in key: 
        excl_group_1[key] = seqs[key]

for key in seqs:
    if not taxa[0][1] in key:
        excl_group_2[key] = seqs[key]

group1_digests = digester(group_1)
group2_digests = digester(group_2) 
excl_group_1_digests = digester(excl_group_1)
excl_group_2_digests = digester(excl_group_2)

c_1 = Counter(group1_digests)
c_1_cent = [(i, c_1[i] / len(group_1) * 100.0) for i,count in c_1.most_common()]
c1top = [x[0] for x in c_1_cent if x[1]>50]
c1top_occ = [round(x[1],3) for x in c_1_cent if x[1]>50]

c_2 = Counter(group2_digests)
c_2_cent = [(i, c_2[i] / len(group_2) * 100.0) for i,count in c_2.most_common()]
c2top = [x[0] for x in c_2_cent if x[1]>50]
c2top_occ = [round(x[1],3) for x in c_2_cent if x[1]>50]

if len(taxa[0][0]) < 2:

    print('\n'+taxa[0][0]+ '\n'+'Common in most '+ taxa[0][0]+ '\n')
    for i in range(len(c1top)):
        temp = round(molecular_weight(c1top[i],'protein'),3)
        print(c1top[i] + f' (MW:{temp}) ' +f'----{c1top_occ[i]}')

print('\n'+ 'Present only in ' + taxa[0][0]+ '\n')

only_c1 = set(c1top) - set(excl_group_1_digests) #present only in mammalia top 50 but not in any other species. 

if len(only_c1) != 0:
    for i in only_c1:
        try:
            temp = round(molecular_weight(i,'protein'),3)
            print(i + f' (MW:{temp})')
        except ValueError:
            print('ValueError' + i)
            continue
else: 
    print('No unique marker found')

print('\n'+ f'Present in {taxa[0][0]} but not in {taxa[0][1]}' '\n')

c1notc2 = set(c1top) - set(group2_digests) #present in mammalia top 50 but does not match any of the aves digests. 

for i in c1notc2:
    try:
        temp = round(molecular_weight(i,'protein'),3)
        print(i + f' (MW:{temp})')
    except ValueError:
        print('ValueError' + i)
        continue

if len(taxa[0][1]) < 2:
    print('\n'+taxa[0][1]+ '\n'+'Common in most '+ taxa[0][1]+ '\n')
    for i in range(len(c2top)):
        temp = round(molecular_weight(c2top[i],'protein'),3)
        print(c2top[i] + f' (MW:{temp}) ' + f'----{c2top_occ[i]}')


print('\n'+ f'Present in {taxa[0][1]} but not in {taxa[0][0]}' '\n')

c2notc1 = set(c2top) - set(group1_digests) #present in aves top 50 but does not match any of the mammalia digests. 

for i in c2notc1:
    temp = round(molecular_weight(i,'protein'),3)
    print(i + f' (MW:{temp})')

print('\n'+'Present only in '+ taxa[0][1]+ '\n')

only_c2 = set(c2top) - set(excl_group_2_digests)  #present only in aves top 50 but not in any other species. 

if len(only_c2) != 0:
    for i in only_c2:
        temp = round(molecular_weight(i,'protein'),3)
        print(i + f' (MW:{temp})')
else:
    print('No unique marker found')

print('\n')

with open('test.txt', 'w') as f:

    print(taxa[0][0]+ '\n', file=f)

    for key,value in group_1.items():
        for j in  c1notc2:
            if not j in value:
                temp = round(molecular_weight(j,'protein'),3)
                print( j + f' (MW:{temp}) ' + '----------- not found in: ' + re.findall('\[.*\]', key)[0][1:-1], file=f)

    print('\n'+taxa[0][1]+ '\n', file=f)      

    for key,value in group_2.items():
        for j in  list(c2notc1):
            if not j in value:
                temp = round(molecular_weight(j,'protein'),3)
                print( j + f' (MW:{temp}) ' + '----------- not found in: ' + re.findall('\[.*\]', key)[0][1:-1], file=f)













