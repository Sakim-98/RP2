import requests
import pandas as pd
import os 
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    '-n',
    '--file_name',
    help=
    'Write the name of the csv file (including the.csv suffix)',
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

args = parser.parse_args()  #store the arguments

IDs = pd.read_csv(args.file_name[0][0]) #read the csv file downloaded from Ensembl (the search results)

with open('placeholder.fasta', 'w') as f:

    for i in range(len(IDs['id_with_url'])):

        r = requests.get(f"https://rest.ensembl.org/sequence/id/{IDs['id_with_url'][i]}?type=protein", headers={ "Content-Type" : "text/x-fasta"})

        if r.text[0] == '>':

            print('<'+IDs['species'][i], '\n', file=f)
            
            print(r.text, file=f)



with open("placeholder.fasta", "r") as f:
    
    lines = f.readlines()
    
with open(args.output_file_name[0][0], "w") as f:
    i = 0
    while i < len(lines):
        if lines[i][0] == '<':
            temp = lines[i+2].rstrip()+'_'+lines[i][1:].rstrip()
            print(temp, '\n', file=f)
            i+=2
        else:
            print(lines[i], file=f)
        i+=1

os.remove('placeholder.fasta')
