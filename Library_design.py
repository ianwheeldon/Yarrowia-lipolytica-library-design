#!/usr/bin/env python2.7
import gffutils
import itertools, csv
from subprocess import Popen, PIPE
from subprocess import check_output
import subprocess, io
import pandas as pd
from StringIO import StringIO
from scipy.stats import pearsonr
import string
import math

def isheader(line):
	return line[0] == '>'
def aspairs(f):
	seq_id = ''
	sequence = ''
	for header,group in itertools.groupby(f, isheader):
		if header:
			line = next(group)
			seq_id = line[1:].split()[0]
		else:
			sequence = ''.join(line.strip() for line in group)
			yield seq_id, sequence

gff_file = "GCF_001761485.1_ASM176148v1_genomic.gff"
fasta_file = "CLIB89/CLIB89.fasta"

#Reading the annotation file enabling us to find the start and end position of the CDSs
with open(gff_file, "rt") as file:
	data = file.read()
gffutils.create_db(data, dbfn='database.db', force=True, from_string=True,merge_strategy="merge", sort_attribute_values=True)
db = gffutils.FeatureDB('database.db', keep_order=True)
gene_ID = []
for i in list(db.features_of_type("gene")):
	gene_ID.append(i['ID'][0])

#Creating a dataframe composed of the gene IDs and the start and end location of each gene:
key = []
start = []
end = []
ID = []
chr = []
strand = []
for a in gene_ID:
	type = ['CDS', 'tRNA']
	gene = db[a]
	for j in type:
		for k in (db.children(gene, featuretype=j , order_by='start')):
			chr.append(k.seqid)
			start.append(k.start)
			end.append(k.end)
			key.append(k['ID'][0])
			ID.append(a.split('-')[1])
			strand.append(k.strand)

db_dict = {'Chromosome_ID': chr, 'Gene_ID': ID, 'CDS_ID': key, 'Start': start, 'End': end, 'Strand': strand}
db_df = pd.DataFrame(db_dict)

chrs = {}
with open(fasta_file, "rt") as file:
        seqs = aspairs(file)
        for seq in seqs:
                chr_id = seq[0]
                chr_str = seq[1]
                chr = []
                chr.append(chr_str)
                chrs[str(chr_id)] = chr

PAM = ['AGG', 'TGG', 'CGG', 'GGG']
#Running CHOPCHOP v3 to find all the sgRNAs in the first 300 bp of each CDS or tRNA gene:
chopchop_df = pd.DataFrame()
best_library = pd.DataFrame()

for index, row in db_df.iterrows():
	print row['Gene_ID']
	chromo = row['Chromosome_ID']
	if (row['End'] - row['Start'] < 500):
		start = row['Start']
		end = row['End']

	if (row['End'] - row['Start'] >= 500) and (row['Strand'] == '+'):
		start = row['Start']
		end = row['Start'] + 0.4 * (row['End'] - row['Start'])

	if (row['End'] - row['Start'] >= 500) and (row['Strand'] == '-'):
		start = row['End'] - 0.4 * (row['End'] - row['Start'])
		end = row['End']

	if start != end:
		ALL = check_output(["./chopchop.py", "-T", "1", "-M", "NGG", "--maxMismatches", "3", "-g", "20", "-G",  "CLIB89" , "-o", "Results" , "-Target", "%s:%d-%d"%(chromo, start, end), "--scoringMethod", "ALL"], universal_newlines=True)
		if len(ALL) != 0:
			ALL = ALL.decode("utf-8")
			data = io.StringIO(ALL)
			df = pd.read_csv(data, sep='\t')
			df['Target sequence'] = df['Target sequence'].apply(lambda x : x[0:20])
			df['Seed_sequence'] = df['Target sequence'].apply(lambda x : x[8:20])
			df['CDS_ID'] = [row['CDS_ID']] * df.shape[0]
			df['Gene_ID'] = [row['Gene_ID']] * df.shape[0]
			df['naive_score'] = df['XU_2015'] + df['DOENCH_2014'] + df['DOENCH_2016']/100 + df['MORENO_MATEOS_2015'] + df['ZHANG_2019']/100 + df['ALKAN_2018']/100

			#Calculating Seed_MM0 and the Quality score for each sgRNA using Bowtie:
			for a, b in df.iterrows():
				Seed_MM0 = 0
				for i in PAM:
					data = check_output(["bowtie/bowtie", "-a", "-p", "4", "-v", "0", "CLIB89/CLIB89", "-c", "%s"%(b['Seed_sequence']+i)], universal_newlines=True)
					if len(data) != 0:
						data = data.decode("utf-8")
						data = io.StringIO(data)
						dataframe = pd.read_csv(data, sep='\t', header=None, usecols=[0, 1, 2, 3, 4, 5, 6, 7], engine='python')
						Seed_MM0 += dataframe.shape[0]
				df.at[a, "Seed_MM0"] = int(Seed_MM0) - 1
			#Defining the Quality score for each sgRNA:
			for a, b in df.iterrows():
				if b['Seed_MM0'] == 0 and b['MM0'] == 0 and b['MM1'] == 0 and b['MM2'] == 0 and b['MM3'] == 0 and b['Self-complementarity'] == 0:
					df.at[a, "Quality_score"] = 1
				elif b['Seed_MM0'] == 0 and b['MM0'] == 0 and b['MM1'] == 0 and b['MM2'] == 0 and b['MM3'] == 1 and b['Self-complementarity'] == 0:
					df.at[a, "Quality_score"] = 2
				elif b['Seed_MM0'] == 0 and b['MM0'] == 0 and b['MM1'] == 0 and b['MM2'] == 0 and b['MM3'] == 2 and b['Self-complementarity'] == 0:
					df.at[a, "Quality_score"] = 3
				elif b['Seed_MM0'] == 0 and b['MM0'] == 0 and b['MM1'] == 0 and b['MM2'] == 0 and b['MM3'] == 3 and b['Self-complementarity'] == 0:
					df.at[a, "Quality_score"] = 4
				elif b['Seed_MM0'] == 0 and b['MM0'] == 0 and b['MM1'] == 0 and b['MM2'] == 1 and b['MM3'] == 0 and b['Self-complementarity'] == 0:
					df.at[a, "Quality_score"] = 5
				elif b['Seed_MM0'] == 0 and b['MM0'] == 0 and b['MM1'] == 0 and b['MM2'] == 1 and b['MM3'] == 1 and b['Self-complementarity'] == 0:
					df.at[a, "Quality_score"] = 6
				elif b['Seed_MM0'] == 0 and b['MM0'] == 0 and b['MM1'] == 0 and b['MM2'] == 1 and b['MM3'] == 2 and b['Self-complementarity'] == 0:
					df.at[a, "Quality_score"] = 7
				elif b['Seed_MM0'] == 0 and b['MM0'] == 0 and b['MM1'] == 0 and b['MM2'] == 1 and b['MM3'] == 3 and b['Self-complementarity'] == 0:
					df.at[a, "Quality_score"] = 8
				elif b['Seed_MM0'] == 0 and b['MM0'] == 0 and b['MM1'] == 0 and b['MM2'] == 2 and b['MM3'] == 0 and b['Self-complementarity'] == 0:
					df.at[a, "Quality_score"] = 9
				elif b['Seed_MM0'] == 0 and b['MM0'] == 0 and b['MM1'] == 0 and b['MM2'] == 2 and b['MM3'] == 1 and b['Self-complementarity'] == 0:
					df.at[a, "Quality_score"] = 10
				elif b['Seed_MM0'] == 0 and b['MM0'] == 0 and b['MM1'] == 0 and b['MM2'] == 2 and b['MM3'] == 2 and b['Self-complementarity'] == 0:
					df.at[a, "Quality_score"] = 11
				elif b['Seed_MM0'] == 0 and b['MM0'] == 0 and b['MM1'] == 0 and b['MM2'] == 2 and b['MM3'] == 3 and b['Self-complementarity'] == 0:
					df.at[a, "Quality_score"] = 12
				elif b['Seed_MM0'] == 1 and b['MM0'] == 1:
					df.at[a, "Quality_score"] = 13
				else:
					df.at[a, "Quality_score"] = 14

			#Ranking sgRNAs designed for each CDS based on their Quality_score and naive_score
			df.sort_values(by=['Quality_score', 'naive_score'], ascending=[True, False], inplace=True)
			chopchop_df = pd.concat([chopchop_df, df])
			chopchop_df.to_csv('CHOPCHOP_Total.csv')
			#Choosing the first best six sgRNAs for the final library
			best_library = pd.concat([best_library, df.loc[(df['Quality_score'] != 14)].head(6)])
			best_library.to_csv('BEST_LIBRARY.csv')

best_library.reset_index(drop=True, inplace=True)

#Adding a unique ID for each sgRNA in the best_library:
counter = 1
for index, row in best_library.iterrows():
	best_library.at[index, 'Unique_ID'] = 'sgRNA_' + str(counter)
	counter += 1

#Creating the non-targeting controls which correspond to 1% of the library size:
non_targeting_dict = {}
while len(non_targeting_dict.keys()) <= int(math.ceil(0.01 * best_library.shape[0])):
	ALL = check_output(["./control_guides.py", "--type", "Cas9", "--PAM", "NGG", "--how_many", "%d"%(1), "CLIB89/CLIB89.fasta"], universal_newlines=True)
	ALL = ALL.decode("utf-8")
	data = io.StringIO(ALL)
	df = pd.read_csv(data, sep='\t')
	for index, row in df.iterrows():
		bowtie_list = []
		for i in PAM:
                        data = check_output(["bowtie/bowtie", "-a", "-p", "4", "-v", "3", "CLIB89/CLIB89", "-c", "%s"%(row['Target sequence'][0:20] + i)], universal_newlines=True)
                        data1 = check_output(["bowtie/bowtie", "-a", "-p", "4", "-v", "0", "CLIB89/CLIB89", "-c", "%s"%(row['Target sequence'][8:20] + i)], universal_newlines=True)
                        data2 = check_output(["bowtie/bowtie", "-a", "-p", "4", "-v", "0", "CLIB89/CLIB89", "-c", "%s"%(row['Target sequence'][9:20] + i)], universal_newlines=True)
                        bowtie_list.extend([len(data), len(data1), len(data2)])
                if all(x == 0 for x in bowtie_list):
                        non_targeting_dict[row['Target sequence'][0:20]] = []
                        non_targeting_dict[row['Target sequence'][0:20]].extend([row['GC content (%)'], 'sgRNA_' + str(counter), 'Non-targeting control', 'Non-targeting control'])
                        counter += 1
non_targeting_df = pd.DataFrame.from_dict(non_targeting_dict, orient='index', columns=['GC content (%)', 'Unique_ID', 'Gene_ID', 'CDS_ID'])
non_targeting_df = non_targeting_df.reset_index()
non_targeting_df.rename(columns={'index': 'Target sequence'}, inplace=True)

non_targeting_df.to_csv('non_targeting.csv')
best_library.to_csv('BEST_LIBRARY.csv')

