#Python file to convert text file with pvalues,fold-change, and sequence data into netphorest predictions

import pandas as pd
import os
import numpy as np

class SeqConvert:

	def __init__(self,filename,threshold=1.0):
		self.filename = filename
		self.threshold = threshold

	def get_trimmed_netphorest_frame(self):
		all_series = list()
		netphorest_frame, sequence_frame = self.get_netphorest_frame()

		print("\tTrimming NetPhorest DataFrame")
		
		for row in netphorest_frame.iterrows():
			peptide_full_sequence = row[0]
			peptide_info = sequence_frame.ix[peptide_full_sequence]
			sites = peptide_info['sites']
			curr_series = row[1]
			curr_site = curr_series['Position']
			if(curr_site in sites):
				all_series.append(curr_series)
		netphorest_frame = pd.DataFrame(all_series)
		netphorest_frame['pval'] = sequence_frame['pval']
		netphorest_frame['fc'] = sequence_frame['fc']
		return netphorest_frame

	#Returns a data frame containing all the relevant Netphorest Phosphorylation site predictions
	def get_netphorest_frame(self):
		sequence_frame = self.get_sequence_frame()
		print("\tGenerating NetPhorest DataFrame")
		self.write_sequence_file(".sequences.fasta",sequence_frame)
		os.system("cat .sequences.fasta | ./netphorest > .NetPhorest_Predictions.txt")
		netphorest_frame = pd.DataFrame.from_csv(".NetPhorest_Predictions.txt",sep="\t",index_col=0)	
		cols = netphorest_frame.columns.values
		renamed_columns = {cols[i - 1]: cols[i] for i in range(1,7)}
		renamed_columns['Classifier'] = 'Prediction'
		netphorest_frame = netphorest_frame.rename(columns=renamed_columns)
		netphorest_frame.index.name = 'Peptide #'
		return netphorest_frame,sequence_frame

	#Reads sequences, p-vaylues, fold-change values, and site information from text file.
	#Format of text file that this method expects will be documented later
	def get_sequence_frame(self):
		print("\tGenerating Sequence DataFrame")
		df = pd.DataFrame.from_csv(self.filename,sep=" ",index_col=-1)
		sequences,sites = self.get_sequences_and_sites(df)
		df['stripped_sequence'] = sequences
		df['sites'] = sites
		df = df.groupby(df.index,sort=False).first()
		return df

	def get_sequences_and_sites(self,sequence_frame):
		all_sequences = list()
		sites = list()
		for row in sequence_frame.iterrows():
			sequence = row[0]
			build_seq = str()
			build_sites = list()
			loc = 0
			in_brace = False
			for character in sequence:
				if(not in_brace and not character == '['):
					build_seq += character
					loc += 1
				elif(character == ']'):
					in_brace = False
				elif(character == '['):
					build_sites.append(loc)
					loc += 1
					in_brace = True
			all_sequences.append(build_seq)
			sites.append(build_sites)
		return all_sequences,sites


	def write_sequence_file(self, outfile, sequence_frame):
		sequence_writer = open(outfile,'wb')
		for row in sequence_frame.iterrows():
			sequence_with_site = row[0]
			stripped_sequence = row[1]['stripped_sequence']
			sequence_writer.write(">%s\n%s\n" % (sequence_with_site,stripped_sequence))
		sequence_writer.close()
		

def test():
	sq = SeqConvert("Testing_Data/test.txt")
	print(sq.get_trimmed_netphorest_frame())

test()