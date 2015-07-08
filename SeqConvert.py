#Python file to convert text file with pvalues,fold-change, and sequence data into netphorest predictions

import pandas as pd
import os
import numpy as np

class SeqConvert:

	def __init__(self,filename):
		self.filename = filename

	def get_trimmed_netphorest_frame(self):
		all_series = list()
		netphorest_frame, sequence_frame = self.get_netphorest_frame()
		for row in netphorest_frame.iterrows():
			peptide_num = row[0]
			peptide_info = sequence_frame.ix[int(peptide_num)]
			sites = peptide_info['sites']
			curr_series = row[1]
			curr_site = curr_series['Position']
			if(curr_site in sites):
				all_series.append(curr_series)
		netphorest_frame = pd.DataFrame(all_series)
		sequence_frame.index = [str(x) for x in range(0,len(sequence_frame.index))]
		netphorest_frame['pval'] = sequence_frame[['pval']]
		netphorest_frame['fc'] = sequence_frame[['fc']]
		return netphorest_frame

	#Returns a data frame containing all the relevant Netphorest Phosphorylation site predictions
	def get_netphorest_frame(self):
		sequence_frame = self.get_sequence_frame()
		self.write_sequence_file("sequences.fasta",sequence_frame)
		os.system("cat sequences.fasta | ./netphorest > NetPhorest_Predictions.txt")
		netphorest_frame = pd.DataFrame.from_csv("NetPhorest_Predictions.txt",sep="\t",index_col=0)
		cols = netphorest_frame.columns.values
		renamed_columns = {cols[i - 1]: cols[i] for i in range(1,7)}
		renamed_columns['Classifier'] = 'Prediction'
		netphorest_frame = netphorest_frame.rename(columns=renamed_columns)
		netphorest_frame.index.name = 'Peptide #'
		return netphorest_frame,sequence_frame

	#Reads sequences, p-vaylues, fold-change values, and site information from text file.
	#Format of text file that this method expects will be documented later
	def get_sequence_frame(self):
		df = pd.DataFrame.from_csv(self.filename,sep=" ",index_col=-1)
		df.drop_duplicates(inplace=True)
		indices = list(df.index)
		sequences, sites = self.strip_sites(indices)
		renamed_index = {indices[i]: sequence for i,sequence in enumerate(sequences)}
		df = df.rename(index=renamed_index)
		df['sites'] = sites
		return df

	#Removes Mass-Spec charge information from sequence strings
	#Returns a tuple containing:
	#	1) A list containing all sequence strings without the phosphorylation site data
	#	2) A list containing the phosphorylation sites for each of the sequnces
	def strip_sites(self,sequences_with_sites):
		all_sequences = list()
		sites = list()
		for curr_seq in sequences_with_sites:
			build_seq = str()
			build_sites = list()
			loc = 0
			in_brace = False
			for character in curr_seq:
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
		with open(outfile, 'wb') as sequences_file:
			for index, sequence in enumerate(list(sequence_frame.index)):
				sequences_file.write(">%s\n%s\n" % (str(index),sequence))