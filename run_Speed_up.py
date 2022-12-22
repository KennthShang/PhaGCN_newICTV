import numpy as np
import pandas as pd
import os
import Bio
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import pandas as pd
import subprocess
import argparse
import re

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--contigs', type=str, default = 'contigs.fa')
parser.add_argument('--len', type=int, default=8000)
parser.add_argument('--threads', type=int, default=8)
args = parser.parse_args()

threads = args.threads

if not os.path.exists("input"):
    _ = os.makedirs("input")
else:
    print("folder {0} exist... cleaning dictionary".format("input"))
    if os.listdir("input"):
        try:
            _ = subprocess.check_call("rm -rf {0}".format("input"), shell=True)
            _ = os.makedirs("input")
            print("Dictionary cleaned")
        except:
            print("Cannot clean your folder... permission denied")
            exit(1)

if not os.path.exists("pred"):
    _ = os.makedirs("pred")
else:
    print("folder {0} exist... cleaning dictionary".format("pred"))
    if os.listdir("pred"):
        try:
            _ = subprocess.check_call("rm -rf {0}".format("pred"), shell=True)
            _ = os.makedirs("pred")
            print("Dictionary cleaned")
        except:
            print("Cannot clean your folder... permission denied")
            exit(1)



if not os.path.exists("Split_files"):
    _ = os.makedirs("Split_files")
else:
    print("folder {0} exist... cleaning dictionary".format("Split_files"))
    if os.listdir("Split_files"):
        try:
            _ = subprocess.check_call("rm -rf {0}".format("Split_files"), shell=True)
            _ = os.makedirs("Split_files")
        except:
            print("Cannot clean your folder... permission denied")
            exit(1)


try:
    if os.path.exists('database/database.self-diamond.tab.abc'):
        print(f'Using preformatted DIAMOND database ...')
    else:
        make_diamond_cmd = 'diamond makedb --threads {threads} --in database/Caudovirales_protein.fasta -d database/database.dmnd'
        print("Creating Diamond database...")
        _ = subprocess.check_call(make_diamond_cmd, shell=True)
        diamond_cmd = 'diamond blastp --threads {threads} --sensitive -d database/database.dmnd -q database/Caudovirales_protein.fasta -o database/database.self-diamond.tab'
        print("Running Diamond...")
        _ = subprocess.check_call(diamond_cmd, shell=True)
        diamond_out_fp = "database/database.self-diamond.tab"
        database_abc_fp = "database/database.self-diamond.tab.abc"
        _ = subprocess.check_call("awk '$1!=$2 {{print $1,$2,$11}}' {0} > {1}".format(diamond_out_fp, database_abc_fp), shell=True)
except:
    print("create database failed")
    exit(1)


#####################################################################
##########################    Start Program  ########################
#####################################################################

#def special_match(strg, search=re.compile(r'[^ACGT]').search):
#    return not bool(search(strg))

##  Filter unknown family 
query_file = f"{args.contigs}"
db_virus_prefix = f"database/unknown_db/db"
output_file = f"unknown_out.tab"
virus_call = NcbiblastnCommandline(query=query_file,db=db_virus_prefix,out=output_file,outfmt="6 qseqid sseqid evalue pident length qlen", evalue=1e-10,
                                 task='megablast',perc_identity=95,num_threads=threads)
virus_call()


check_unknown = {}
with open(output_file) as file_out:
    for line in file_out.readlines():
        parse = line.replace("\n", "").split("\t")
        virus = parse[0]
        ident  = float(parse[-3])
        length = float(parse[-2])
        qlen   = float(parse[-1])
        if length/qlen > 0.95 and ident > 0.95:
            check_unknown[virus] = 1

try:
    with open('phage_do_have_family.csv', 'w') as file:
        file.write('Accession,Family,Score\n')
        for name in check_unknown:
            file.write(f'{name},no_family_avaliable,1\n')
except:
    pass
        
rec = []
for record in SeqIO.parse(f'{args.contigs}', 'fasta'):
    try:
        if check_unknown[record.id]:
            continue
    except:
        rec.append(record)

SeqIO.write(rec, f'filtered_contigs.fa', 'fasta')


cnt = 0
file_id = 0
records = []
for record in SeqIO.parse('filtered_contigs.fa', 'fasta'):
    if cnt !=0 and cnt%1000 == 0:
        SeqIO.write(records, "Split_files/contig_"+str(file_id)+".fasta","fasta") 
        records = []
        file_id+=1
        cnt = 0
    seq = str(record.seq)
    seq = seq.upper()
    if 1:#special_match(seq):
        if len(record.seq) > args.len:
            records.append(record)
            cnt+=1

SeqIO.write(records, "Split_files/contig_"+str(file_id)+".fasta","fasta")
file_id+=1

for i in range(file_id):
    cmd = "mv Split_files/contig_"+str(i)+".fasta input/"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("Moving file Error for file {0}".format("contig_"+str(i)))
        continue

    cmd = "python run_CNN.py"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("Pre-trained CNN Error for file {0}".format("contig_"+str(i)))
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        continue
        

    cmd = f"python run_KnowledgeGraph.py --threads {args.threads}"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("Knowledge Graph Error for file {0}".format("contig_"+str(i)))
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        continue
    
    cmd = "python run_GCN.py"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("GCN Error for file {0}".format("contig_"+str(i)))
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        continue

    # Clean files
    cmd = "rm input/*"
    out = subprocess.check_call(cmd, shell=True)

    name_list = pd.read_csv("name_list.csv")
    prediction = pd.read_csv("Cyber_data/prediction.csv")
    prediction = prediction.rename(columns={'contig_names':'idx'})
    contig_to_pred = pd.merge(name_list, prediction, on='idx')
    contig_to_pred.to_csv("pred/contig_"+str(i)+".csv", index = None)

    cmd = "rm name_list.csv"
    out = subprocess.check_call(cmd, shell=True)

df_list = []    
for file in sorted(list(os.listdir('pred'))):
    if 'contig' in file:
        df_list.append(pd.read_csv(f'pred/{file}'))

df = pd.concat(df_list)
df.to_csv('final_prediction.csv', index=False)
#cmd = "cat pred/* > final_prediction.csv"
#out = subprocess.check_call(cmd, shell=True)

