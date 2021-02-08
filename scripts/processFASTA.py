'''
Replaces any ' ' in sample names with '_'
'''

from Bio import SeqIO
import sys

records = []
for record in SeqIO.parse(sys.argv[1], 'fasta'):
	record.id = record.description.replace(' ', '_')
	record.description = ''
	records.append(record)

SeqIO.write(records, f'{sys.argv[1]}.proc', 'fasta')