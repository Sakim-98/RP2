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

with open (args.file_name[0][0], 'r') as f: 
    
	file = f.readlines() # read the content of the file

with open ('placeholder.fasta', 'w') as f:
    
    for line in file: # replace the spaces in the ID with '_' for later manipulation in SeqIO operations
        if line[0] == '>':
            temp = '_'.join(line.split())
            f.write(temp+'\n'+'\n')
        else:
            f.write(line+'\n')

pre_seqs = dict() # Entire sequences 

for seq_record in SeqIO.parse("placeholder.fasta", "fasta"):
    pre_seqs[seq_record.id] = seq_record.seq # store each fasta sequence in dictionary where they keys are the FASTA identifiers

if args.unwanted_sequences:
    x = args.unwanted_sequences[0][0].split(',')

    for i in list(pre_seqs.keys()):

        if len(x) == 1:

            if x[0] in i or pre_seqs[i].count('X') > 5: # remove unwanted sequences
                del pre_seqs[i]

        elif len(x) == 2:

            if x[0] in i or x[1] in i or pre_seqs[i].count('X') > 5: # remove unwanted sequences
                del pre_seqs[i]

        elif len(x) == 3:

            if x[0] in i or x[1] in i or x[2] in i or pre_seqs[i].count('X') > 5: # remove unwanted sequences
                del pre_seqs[i]
else:

    for i in list(pre_seqs.keys()):
        if 'partial' in i or 'LOW_QUALITY' in i or pre_seqs[i].count('X') > 5: # remove unwanted sequences
            del pre_seqs[i]

seqs = dict() # only portion of the sequence involved in chain development

if args.regex:
    for i in pre_seqs:
        temp = re.findall(str(args.regex[0][0]), str(pre_seqs[i]))
        if not len(temp) == 0:
            seqs[i] = temp[0] # store sequence according to the rules
else:
    seqs = pre_seqs

species_list = list()

for i in list(seqs.keys()):
    temp = re.findall('\[.*\]', i)
    species_list.append(temp[0][1:-1]) # store the name of the species in a list

species_duplist = list()

for i in species_list:
    if species_list.count(i) > 1:
        species_duplist.append(i) # count if species name is present multiple times

temp_dict = dict()
temp_dict_sansdup = dict()

for i in seqs: 
    if re.findall('\[.*\]', i)[0][1:-1] in species_duplist:
        temp_dict[i] = seqs[i] # check if species is present multiple times
    else:
        temp_dict_sansdup[i] = seqs[i] # if not, add to the no duplicate dict
        
    for key,value in temp_dict.items(): # for multiply present species, get rid of one randomly if sequences already present
        if value not in temp_dict_sansdup.values():
            temp_dict_sansdup[key] = value

with open (args.output_file_name[0][0] , 'w') as f: # write the entire thing to a file
    for i in seqs:
        f.write('>'+i + '\n' + '\n')
        wrapper = textwrap.TextWrapper(width=80)
        string = wrapper.fill(text=str(seqs[i]))
        f.write (string+'\n' + '\n')

os.remove('placeholder.fasta') # delete the placeholder file