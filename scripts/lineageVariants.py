'''
Part of pipeline to find all common variants in SARS-CoV-2 lineages. Takes a VCF containing SNP variants a list of
variants from minimap2 and paftools.js call and parses them and combines them into a single VCF filters the deletions
present at a lower frequency than the cutoff, and turns the VCF into a BED for both coding sequence variants and
nucleotides variants for use in the UCSC Genome Browser. Variants in coding regions will be named as the amino acid 
mutation name [Ref AA][AA coordinate in peptide][Alt AA] unless they are a indels.

USAGE: python parse.py -v [Variant List] -r [Reference Sequence]
'''

import pandas as pd
import argparse
from Bio import SeqIO
import pdb


class CommandLine() :
    '''
    argparse command line.
    '''

    def __init__(self, inOpts=None) :
        
        self.parser = argparse.ArgumentParser(description = 'Converts variant list to VCF and BED', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s -v [Vairant list] -f [Frequency Cutoff]' 
                                             )
        self.parser.add_argument('-va', '--vars', action = 'store', help='Path to file containing variant list')
        self.parser.add_argument('-vc', '--vcf', action = 'store', help='Path to file containing SNP VCF')
        self.parser.add_argument('-r', '--ref', action = 'store', help='Path to file containing SARS-CoV-2 reference sequence')
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)


def highFrequencyDeletions(data):
	'''
	Using list of vars from minimap2 -cs | paftools.js call (position ranges given are 0-based and half-open) creates
	dictionary containing only indels that have a frequency above the cutoff.
	'''
	cutoff = 0.85

	# First create dictionary that contains counts of all indels. Dictionary entries have form:
	# {([start], [end]) : {([ref], [alt]) : [count]}}
	variants = {}
	samples = len(set(data['query'].tolist()))
	for i in range(1, len(data)):
		start,end,ref,alt,depth  = data.iloc[i][['start', 'end', 'ref', 'alt', 'depth']]
		if '-' in [ref, alt]:
			if (start, end) not in variants.keys():
				variants[(start,end)] = {}
				variants[(start, end)][(ref, alt)] = 1
			else:
				if (ref, alt) not in variants[(start,end)].keys():
					variants[(start, end)][(ref, alt)] = 1
				else:
					variants[(start, end)][(ref, alt)] += 1

	# Go through deletion dictionary and add indels with a frequency above the cutoff to the high_freq dictionary.
	# high_freq dictionary has the form: {([start], [end]) : {([ref], [alt]) : ([count], [frequency])}}
	high_freq = {}
	for coords in variants.keys():
		for alleles in variants[coords].keys():
			if variants[coords][alleles]/samples > cutoff:
				high_freq[coords] = {}
				high_freq[coords][alleles] = variants[coords][alleles]
				high_freq[coords][alleles] = (high_freq[coords][alleles], round(high_freq[coords][alleles]/samples, 5))
	return high_freq


def createVCF(high_freq, chro, vars, vcf, samples, ref_seq):
	'''
	Creates a VCF file from the high frequency variants.
	'''

	# First reads in SNPs from previous part of pipeline
	# Adds AF (allele frequency) and CS (codon start) fields to INFO tag
	# IDs given according to the location in the genome and their effect (syn = synonomous)
	snpVCF = pd.read_csv(vcf, sep='\t', skiprows=3)
	syn = []
	for i in range(len(snpVCF)):
		code, codon_start, peptide = genome2Peptide(snpVCF.iloc[i]['POS']-1, snpVCF.iloc[i]['ALT'], snpVCF.iloc[i]['REF'], ref_seq)
		if code[0] == code[-1]:
			snpVCF.loc[i, 'ID'] = f'{peptide}:syn'
		else:
			snpVCF.loc[i, 'ID'] = f'{peptide}:{code}'
		info = snpVCF.iloc[i]['INFO']
		ac = int(info.split('=')[1].split(';')[0])
		an = int(info.split('=')[-1])
		snpVCF.loc[i, 'INFO'] = f'{info};AF={round(ac/an, 5)}'
		snpVCF.loc[i, 'INFO'] += f';CS={codon_start}'

	# Then goes through high_freq dictionary of high frequency indels and creates a pandas dictionary with same fields
	# as a VCF
	chrom = []
	pos = []
	idd = []
	ref = []
	alt = []
	qual = []
	filt = []
	info = []
	for coords in high_freq.keys():
		for alleles in high_freq[coords].keys():
			code, codon_start, peptide = genome2Peptide(coords[0], alleles[1], alleles[0], ref_seq)
			chrom.append(chro)
			idd.append(f'{peptide}:{code}')
			if alleles[1] == '-':
				pos.append(coords[0])
				new_alt, new_ref = vcfDeletionFormat(coords[0], alleles[0], ref_seq)
				ref.append(new_ref)
				alt.append(new_alt)
			elif alleles[0] == '-':
				pos.append(coords[0])
				new_alt, new_ref = vcfInsertionFormat(coords[0], alleles[1], ref_seq)
				ref.append(new_ref)
				alt.append(new_alt)
			else:
				pos.append(coords[0]+1)
				ref.append(alleles[0].upper())
				alt.append(alleles[1].upper())
			qual.append('.')
			filt.append('.')
			info.append(f'AC={high_freq[coords][alleles][0]};AN={samples};AF={high_freq[coords][alleles][1]};CS={codon_start}')
	del_data = pd.DataFrame({'#CHROM':chrom, 'POS':pos, 'ID':idd, 'REF':ref, 'ALT':alt, 'QUAL':qual, 'FILTER':filt, 'INFO':info})

	# Combine SNP data and indel data into single dataframe, write it to file, then add header lines.
	all_data = pd.concat([snpVCF, del_data])
	all_data = all_data.sort_values(['POS']).reset_index()
	all_data['POS'] = all_data['POS'].astype(int)
	all_data[all_data.columns[1:]].to_csv('.'.join(vars.split('.')[:-1])+'_nohead.vcf', sep='\t', index=None)
	with open(vcf, 'r') as file:
		header = file.readlines()[:3]
	with open('.'.join(vars.split('.')[:-1])+'_nohead.vcf', 'r') as file:
		lines = file.readlines()
	with open('.'.join(vars.split('.')[:-2])+'.vcf', 'w') as file:
		file.writelines(header)
		file.writelines(lines)



def createBED(high_freq, vars, vcf, samples, ref_seq):
	'''
	Create BED files from the high frequency variants. Creates two files: one containing amino acid variants (only 
	showing variants that change amino acid sequence) and one containing all nucleotide variants. Both files are 
	created separately using the VCF file containing all variants.
	'''

	# VCF data of all high frequency variants
	vcf_data = pd.read_csv('.'.join(vars.split('.')[:-2])+'.vcf', sep='\t', skiprows=3)

	# Create amino acid variants BED
	# Variants in the 5' or 3' UTR, intergenic regions or synonomous variants are not included
	# Variants are named according to their amino acid change and their codon position in the peptide
	# [Ref AA][AA coordinate in peptide][Alt AA]
	chrom = []
	start = []
	end = []
	name = []
	for i in range(len(vcf_data)):
		chro, pos, idd, ref, alt, info = vcf_data[['#CHROM','POS','ID', 'REF', 'ALT','INFO']].iloc[i]
		codon_start = int(info.split('=')[-1])
		if 'UTR' not in idd and 'intergenic' not in idd and 'syn' not in idd:
			chrom.append(chro)
			if len(alt) == 1 and len(ref) == 1:
				start.append(codon_start)
				end.append(codon_start+3)
			else:
				if len(ref) > 1:
					start.append(pos)
					end.append(pos+len(ref)-1)
				elif len(alt) > 1:
					start.append(pos)
					end.append(pos+1)
			name.append(idd)
	aa_bed_data = pd.DataFrame({'chrom':chrom, 'start':start, 'end':end, 'name':name})
	aa_bed_data.to_csv('.'.join(vars.split('.')[:-2])+'_aa.bed', sep='\t', index=None, header=False)

	# Create nucleotide variants BED
	# Includes all variants
	# Variants are named according to their reference and alternate nucleotides and their 1-based position in the
	# genome [Ref nuclotide][1-based genome coordinate][Alt nucleotide]
	chrom = []
	start = []
	end = []
	name = []
	for i in range(len(vcf_data)):
		chro, pos, idd, ref, alt, info = vcf_data[['#CHROM','POS','ID', 'REF', 'ALT','INFO']].iloc[i]
		chrom.append(chro)
		if len(alt) == 1 and len(ref) == 1:
			start.append(pos-1)
			end.append(pos)
			name.append(f'{ref}{pos}{alt}')
		else:
			if len(ref) > 1:
				start.append(pos)
				end.append(pos+len(ref)-1)
				name.append(f'del_{pos+1}')
			elif len(alt) > 1:
				start.append(pos)
				end.append(pos+1)
				name.append(f'ins_{pos+1}')
	nuc_bed_data = pd.DataFrame({'chrom':chrom, 'start':start, 'end':end, 'name':name})
	nuc_bed_data.to_csv('.'.join(vars.split('.')[:-2])+'_nuc.bed', sep='\t', index=None, header=False)


def genome2Peptide(genome_coord, alt, ref, ref_seq):
	'''
	Given a variant coordinate in the genome (0-based) plus the alternate allele, generates the amino acid mutation
	code [Ref AA][AA 1-based coordinate][Alt AA] and the genome coordinate (0-based) of the codon. Returns set of three 
	values: the mutation code ([Ref nuclotide][1-based genome coordinate][Alt nucleotide] for noncoding variants and 
	[Ref AA][AA coordinate in peptide][Alt AA] for coding variants or del_[1-based genome coordinate] or ins_[1-based genome coordinate]
	for indels), the 0-based genome coordinate of the codon start site for the variant (0 for noncoding or indel variants),
	and the peptide (or region) of the genome the variant is in. Accounts for ribosomal slippage in translation.
	'''

	# Dictionary containing 1-based start/end coordinates of the peptides in the SARS-CoV-2 genome
    # From https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2
	genes = {"5'UTR": (1,265),
		     'nsp1': (266,805),
			 'nsp2': (806,2719),
			 'nsp3': (2720,8554),
			 'nsp4': (8555,10054),
			 '3CL': (10055,10972),
			 'nsp6': (10973,11842),
			 'nsp7': (11843,12091),
			 'nsp8': (12092,12685),
			 'nsp9': (12686,13024),
			 'nsp10': (13025,13441),
			 'nsp11': (13442,13480),
			 'RdRp': (13442,16236),
			 'Hel': (16237,18039),
			 'ExoN': (18040,19620),
			 'nsp15': (19621,20658),
			 'MT': (20659,21552),
			 'S': (21563,25384),
			 'ORF3a': (25393,26220),
			 'E': (26245,26472),
			 'M': (26523,27191),
			 'ORF6': (27202,27387),
			 'ORF7a': (27394,27759),
			 'ORF7b': (27756,27887),
			 'ORF8': (27894,28259),
			 'N': (28274,29533),
			 'ORF10': (29558,29674),
			 "3'UTR": (29675, 29903)}
	RdRp_slip = 13467 # 1-based genome coordinate of RdRp ribosomal slippage site
	for record in SeqIO.parse(ref_seq, 'fasta'):
		seq = record.seq
	for peptide in genes.keys():
		if genome_coord >= genes[peptide][0] and genome_coord < genes[peptide][1]:
			if peptide in ["5'UTR", "3'UTR"]:
				return f'{seq[genome_coord+1]}{genome_coord}{alt}', 0, peptide
			elif alt == '-':
				return f'del_{genome_coord+1}', 0, peptide
			elif ref == '-':
				return f'ins_{genome_coord+1}', 0, peptide
			if peptide == 'RdRp':
				if genome_coord < RdRp_slip:
					codon = (genome_coord-(genes[peptide][0]-1))//3+1
					codon_start = (codon-1)*3 + (genes[peptide][0]-1)
					codon_seq = seq[codon_start:codon_start+3]
					ref = codon_seq.translate()
					diff = genome_coord-codon_start
					if diff == 2:
						altAA = (codon_seq[:diff]+alt.upper()).translate()
					else:
						altAA = (codon_seq[:diff]+alt.upper()+codon_seq[diff+1:]).translate()
				else:
					codon = (genome_coord-RdRp_slip)//3+1+9
					codon_start = (codon-1-9)*3 + (RdRp_slip)
					codon_seq = seq[codon_start:codon_start+3]
					ref = codon_seq.translate()
					diff = genome_coord-codon_start
					if diff == 2:
						altAA = (codon_seq[:diff]+alt.upper()).translate()
					else:
						altAA = (codon_seq[:diff]+alt.upper()+codon_seq[diff+1:]).translate()
					return f'{ref}{codon}{altAA}', codon_start, peptide

			codon = (genome_coord-(genes[peptide][0]-1))//3+1
			codon_start = (codon-1)*3 + (genes[peptide][0]-1)
			codon_seq = seq[codon_start:codon_start+3]
			ref = codon_seq.translate()
			diff = genome_coord-codon_start
			if diff == 2:
				altAA = (codon_seq[:diff]+alt.upper()).translate()
			else:
				altAA = (codon_seq[:diff]+alt.upper()+codon_seq[diff+1:]).translate()
			return f'{ref}{codon}{altAA}', codon_start, peptide
	if ref == '-':
		return f'ins_{genome_coord+1}', 0, 'intergenic'
	elif alt == '-':
		return f'del_{genome_coord+1}', 0, 'intergenic'	
	else:
		return f'{seq[genome_coord]}{genome_coord+1}{alt}', 0, 'intergenic'


def vcfDeletionFormat(pos, ref, ref_seq):
	'''
	Converts deletions to standard VCF format.
	'''
	for record in SeqIO.parse(ref_seq, 'fasta'):
		seq = record.seq
	new_alt = str(seq[pos-1])
	new_ref = new_alt+ref
	return new_alt, new_ref.upper()

def vcfInsertionFormat(pos, alt, ref_seq):
	'''
	Converts insertions to standard VCF format.
	'''
	for record in SeqIO.parse(ref_seq, 'fasta'):
		seq = record.seq
	new_ref = str(seq[pos-1])
	new_alt = new_ref+alt
	return new_alt.upper(), new_ref



def main():
	'''
	Reads variants list from minimap2, creates high_freq dictionary, then creates VCF and BED files.
	'''
	commands = CommandLine()
	data = pd.read_csv(commands.args.vars, sep='\t', names=[0, 'chr', 'start', 'end', 'depth', 'qual', 'ref', 'alt', 'query', 'query_start', 'query_end', 'strand'])
	high_freq = highFrequencyDeletions(data)
	createVCF(high_freq, data.iloc[0]['chr'], commands.args.vars, commands.args.vcf, len(set(data['query'].tolist())), commands.args.ref)
	createBED(high_freq, commands.args.vars, commands.args.vcf, len(set(data['query'].tolist())), commands.args.ref)


if __name__ == "__main__":
    main()
