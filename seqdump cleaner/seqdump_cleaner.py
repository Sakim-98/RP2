def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

modules = ['Bio','regex','pandas','json','textwrap','argparse']
for module in modules:
	install(module)
	print(f'Installed {}')

from Bio import SeqIO
import textwrap
import regex as re
import os
import argparse

parser = argparse.ArgumentParser() # initiate the parser

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
    '-o',
    '--output_file_name',
    help=
    'Write the name of the fasta file (including the.fasta suffix)',
    type=str,
    action='append',
    nargs='*',
    required=True,
)

parser.add_argument(
    '-r',
    '--regex',
    help=
    'Write the name of the fasta file (including the.fasta suffix)',
    type=str,
    action='append',
    nargs='*',
    required=False,
)

parser.add_argument(
    '-u',
    '--unwanted_sequences',
    help=
    'Write the terms/words in the FASTA ID that you would like to be excluded',
    type=str,
    action='append',
    nargs='*',
    required=False,
)

args = parser.parse_args()  #store the arguments

# with open (args.file_name[0][0], 'r') as f: 
    
# 	file = f.readlines() # read the content of the file

pre_seqs = dict() # Entire sequences 

for seq_record in SeqIO.parse(args.file_name[0][0], "fasta"):
    pre_seqs[seq_record.description] = seq_record.seq # store each fasta sequence in dictionary where they keys are the FASTA identifiers

if args.unwanted_sequences:
    x = args.unwanted_sequences[0][0].split(',')

    for i in list(pre_seqs.keys()):

        if len(x) == 1:

            if x[0] in i or pre_seqs[i].count('X') >= 1: # remove unwanted sequences
                del pre_seqs[i]

        elif len(x) == 2:

            if x[0] in i or x[1] in i or pre_seqs[i].count('X') >= 1: # remove unwanted sequences
                del pre_seqs[i]

        elif len(x) == 3:

            if x[0] in i or x[1] in i or x[2] in i or pre_seqs[i].count('X') >= 1: # remove unwanted sequences
                del pre_seqs[i]
else:

    for i in list(pre_seqs.keys()):
        if 'partial' in i or pre_seqs[i].count('X') > 5: # remove unwanted sequences
            del pre_seqs[i]

seqs = dict() # only portion of the sequence involved in chain development

if args.regex:
    for i in pre_seqs:
        temp = re.findall(str(args.regex[0][0]), str(pre_seqs[i]))
        if not len(temp) == 0:
            seqs[i] = temp[0] # store sequence according to the rules
else:
    seqs = pre_seqs

seen = dict()
        
for key in seqs.copy():
    species = re.findall('\[.*\]', key)[0][1:-1]
    value = seqs[key]
    try:
        if seen[species] == value:
            del seqs[key]
    except KeyError:
        seen[species] = value 

with open (args.output_file_name[0][0] , 'w') as f: # write the entire thing to a file
    for i in seqs:
        f.write('>'+i + '\n' + '\n')
        wrapper = textwrap.TextWrapper(width=80)
        string = wrapper.fill(text=str(seqs[i]))
        f.write (string+'\n' + '\n')
