def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

modules = ['Bio','regex','pandas','json','textwrap','argparse']
for module in modules:
	install(module)
	print('Installed module')

from Bio import SeqIO
import regex as re
import pandas as pd
import json
import textwrap
import argparse

parser = argparse.ArgumentParser() # initiate the parser

parser.add_argument(  
    '-n',
    '--file_name',
    help=
    'Write the name of the file (including the suffix (.txt or .fasta))',
    type=str,
    action='append',
    nargs='*',
    required=True,
)

parser.add_argument(  
    '-o',
    '--output_file_name',
    help=
    'Write the name of the output file (including the suffix (.txt or .fasta))',
    type=str,
    action='append',
    nargs='*',
    required=True,
)

args = parser.parse_args()  #store the arguments

df_id2name = pd.read_csv('longnames.txt', sep = '|', encoding = 'unicode_escape')

with open ('hierarchy.txt', 'r') as f:
    hierarchy = f.readlines()

hierarchy = [i.split('|')[0] for i in hierarchy]
hierarchy = [h.split('-') for h in hierarchy]

syn = pd.read_csv('synonym_links.txt', sep = '|', encoding = 'unicode_escape')

with open (args.file_name[0][0], 'r') as f:
    seqs = f.readlines()

lineage = dict()
all_species = dict()
all_ID = list()

for line in seqs:
    if line[0] == '>':
        try: 
            species = re.findall('\[.*\]', line)[0][1:-1].replace('_',' ')
            ID_index = df_id2name.index[df_id2name['Taxon'] == str(species)].tolist()[0]
            ID = df_id2name['ID'].iloc[ID_index]
            all_species[species] = ID
            all_ID.append(ID)
            for hr in hierarchy:
                if str(ID) in hr:
                    lineage[species] = hr
                    break
           
        except IndexError:
            continue
            
lineage_key_ID = [all_species[x] for x in list(lineage.keys())]

other_lineage = dict()

for i in set(all_ID) - set(lineage_key_ID):
    new_ID = syn['new'].iloc[syn.index[syn['old'] == int(i)].tolist()[0]] 
    species = df_id2name['Taxon'].iloc[df_id2name.index[df_id2name['ID'] == int(new_ID)].tolist()[0]]
    for hr in hierarchy:
        if str(new_ID) in hr:
            other_lineage[species] = hr

all_entries = {**lineage, **other_lineage}
for species in all_entries:
    all_entries[species]=all_entries[species][7:]

dummy = dict()

for species in all_entries:
    temp = [df_id2name['Taxon'].iloc[df_id2name.index[df_id2name['ID'] == int(x)].tolist()[0]] for x in all_entries[species]]
    if len(temp[-1].split()) == 3:
        temp = temp[:-1]
    dummy[species] = temp

temp = ''

for species in dummy:
    for taxa in dummy[species]:
        temp = temp+taxa+'-'
    dummy[species] = temp
    temp = ''

db_prob = set()

with open (args.output_file_name[0][0], 'w') as f:
    for line in seqs:
        try:
            if line[0] == '>':
                species = re.findall('\[.*\]', line)[0][1:-1].replace('_',' ')
                line = line[0]+dummy[species]+line[1:]
            print(line, file=f)
        except KeyError:
            print(line, file=f)
            db_prob.add(species)
            continue

print('The following species need manual intervention:')
for i in db_prob:
	print(i)

