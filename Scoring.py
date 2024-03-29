#!/usr/bin/python

import SeqConvert as seq
import Cleanup as clean
import sys
import operator
import math
import os
import fnmatch
import argparse
import pandas as pd
import random
import time

#./Scoring.py Validation
#./Scoring.py Validation/459_Valid_HuangPKA/ValData33noise.txt Results_Folder

def main():

	parser = argparse.ArgumentParser(description='A phosphoproteome-wide kinase inference algorithm.')
	parser.add_argument('data_location', metavar='d', type=str,help='The filename or directory name where the data is stored.')
	parser.add_argument('--permutations',dest='permutations',type=int,help='The number of permutations to run.')
	parser.add_argument('--threshold',type=float,help='What p-value to threshold at.')
	parser.add_argument('--debug',dest='debug',action='store_true',help='Enabling option gives descriptive error message when program crashes.')

	parser.set_defaults(threshold=1.0,debug=False,permutations=0)
	args = parser.parse_args()

	data_location = args.data_location
	num_iterations = args.permutations
	debug_mode = args.debug
	threshold = args.threshold

	#Handle file case
	if(os.path.isfile(data_location)):
		generate_results(data_location,num_iterations,args.threshold,debug_mode)
	#Handle directory case
	elif(os.path.isdir(data_location)):
		for root, dirnames, filenames in os.walk(data_location):
			for filename in fnmatch.filter(filenames, '*.txt'):
				data_location = os.path.join(root,filename)
				generate_results(data_location,num_iterations,args.threshold,debug_mode)
	else:
		print("SysPhos could not find file or directory")

def generate_results(data_location,num_iterations,threshold,debug_mode):
	if(debug_mode):
		write_results(data_location,num_iterations,threshold)
	else:
		try:
			write_results(data_location,num_iterations,threshold)
		except:
			print("SysPhos encountered an error when processing %s" % data_location)

def write_results(data_location,num_iterations,threshold):
	print("Now working on %s" % data_location)
	savedir = "Results_%s/" % data_location[:data_location.rfind(".")]
	kinase_dir = savedir + "Kinase_Scores/"
	peptide_dir = savedir + "Peptide_Scores/"
	permutation_dir = savedir + "Permutation_Scores"

	clean.clean_existing(savedir)
	os.makedirs(savedir)
	os.makedirs(kinase_dir)
	os.makedirs(peptide_dir)
	os.makedirs(permutation_dir)

	seq_conv = seq.SeqConvert(data_location,threshold=threshold)
	netphorest_frame = seq_conv.get_trimmed_netphorest_frame()

	write_kinase_scores(kinase_dir,netphorest_frame)
	write_peptide_scores(peptide_dir,netphorest_frame)
	if(num_iterations > 0):
		write_permutation_scores(permutation_dir,netphorest_frame,num_iterations)
	clean.clean_all(savedir)
	write_readme(savedir,num_iterations,threshold)

def write_readme(savedir,num_iterations,threshold):
	print("\tWriting README.txt")
	writer = open("%s/%s" % (savedir, "README.txt"),"wb")
	date = time.strftime("%m/%d/%y")
	t = time.strftime("%I:%M:%S")
	writer.write("Results generated on %s at %s\nP-value threshold used: %s\nNumber of permutations performed: %s" % (date,t,threshold,num_iterations))
	writer.close()

def write_kinase_scores(savedir,netphorest_frame):
	print("\tCalculating Kinase Scores")
	for schema in list_all_schemas():
		kinase_outfile = "%s/Kinase_scores_%s.txt" % (savedir,schema[0])
		kinase_writer = open(kinase_outfile,'wb')

		kinase_scores = compute_kinase_scores(netphorest_frame,schema[1])
		sorted_scores = sorted(kinase_scores.items(),key=operator.itemgetter(1))
		for score in sorted_scores:
			kinase_prediction = str(score[0])
			kinase_score = str(score[1])
			kinase_writer.write("%s\t%s\n" % (kinase_prediction,kinase_score))
		kinase_writer.close()

def write_peptide_scores(savedir,netphorest_frame):
	print("\tCalculating Peptide Scores")
	for schema in list_all_schemas():
		peptide_outfile = "%s/Peptide_scores_%s.txt" % (savedir,schema[0])
		peptide_writer = open(peptide_outfile,'wb')
		peptide_power_scores,peptide_sig_scores,peptide_fc_scores,peptide_prediction_conf_scores = compute_peptide_scores(netphorest_frame,schema[1])
		sorted_scores = sorted(peptide_power_scores.items(),key=lambda item: sum(float(peptide) for peptide in item[1].values()))

		peptide_writer.write("Kinase\tPeptide\tScore\tSignificance\tFold-Change\tConfidence\n")
		for entry in sorted_scores:
			kinase_prediction = str(entry[0])
			kinase_score = str(sum(float(x) for x in entry[1].values()))
			peptide_writer.write("%s (Score: %s)\n" % (kinase_prediction,kinase_score))

			peptide_entries = sorted(entry[1].items(),key=lambda item: item[1])
			for peptide_entry in peptide_entries:
				peptide_sequence = peptide_entry[0]
				peptide_score = peptide_entry[1]
				significance = peptide_sig_scores[peptide_sequence]
				fold_change = peptide_fc_scores[peptide_sequence]
				confidence = peptide_prediction_conf_scores[(peptide_sequence,kinase_prediction)]

				peptide_writer.write("\t%s\t%s\t%s\t%s\t%s\n" % (peptide_sequence,peptide_score,significance,fold_change,confidence))
		
		peptide_writer.close()

def write_permutation_scores(savedir,netphorest_frame,num_iterations):
	print("\tWorking on Permutations")
	for schema in list_all_schemas():
		kinase_data = []
		for iteration in xrange(0,num_iterations):
			print("\t\tSchema: %s, Iteration: %s" % (schema[0],str(iteration)))
			kinase_scores = compute_kinase_scores(randomize_frame(netphorest_frame),schema[1])
			kinase_data.append(pd.Series([float(score) for score in kinase_scores.values()], index=kinase_scores.keys(),name="Permutation# %s" % str(iteration)))
		kinase_dataframe = pd.concat(kinase_data,axis=1,keys=[s.name for s in kinase_data])
		kinase_dataframe.to_csv("%s/Permutated_%s.txt" % (savedir,schema[0]), sep="\t")

def randomize_frame(dataframe):
		pvals = list(dataframe['pval'])
		fc = list(dataframe['fc'])

		combined = zip(pvals, fc)
		random.shuffle(combined)
		pvals[:], fc[:] = zip(*combined)

		df = dataframe.copy()
		df['pval'] = pvals
		df['fc'] = fc
		return df

def compute_kinase_scores(netphorest_frame,schema):
	kinase_scores = {}
	
	for row in netphorest_frame.iterrows():
		curr_series = row[1]
		curr_score = schema(curr_series['pval'], curr_series['fc'], curr_series['Posterior'])
		prediction = curr_series['Prediction']
		kinase_scores[prediction] = (kinase_scores[prediction] if prediction in kinase_scores else 0) + curr_score

	return kinase_scores

def compute_peptide_scores(netphorest_frame,schema):
	#maps kinase_name -> (peptide_name -> score)
	peptide_power_scores = {}

	#maps peptide_name -> p-value
	peptide_sig_scores = {}

	#maps peptide_name -> fc value
	peptide_fc_scores = {}

	#maps a tuple (peptide_name,kinase_name) -> confidence
	peptide_prediction_conf_scores = {}


	for row in netphorest_frame.iterrows():
		peptide_sequence = row[0]
		curr_series = row[1]

		prediction = curr_series['Prediction']
		p_value = curr_series['pval']
		fc = curr_series['fc']
		confidence = curr_series['Posterior']

		curr_score = schema(curr_series['pval'], curr_series['fc'], curr_series['Posterior'])

		peptide_power_dict = (peptide_power_scores[prediction] if prediction in peptide_power_scores else {})
		peptide_power_dict[peptide_sequence] = (peptide_power_dict[peptide_sequence] if peptide_sequence in peptide_power_dict else 0) + curr_score
		peptide_power_scores[prediction] = peptide_power_dict

		peptide_sig_scores[peptide_sequence] = p_value
		peptide_fc_scores[peptide_sequence] = fc
		peptide_prediction_conf_scores[(peptide_sequence,prediction)] = confidence

	return peptide_power_scores,peptide_sig_scores,peptide_fc_scores,peptide_prediction_conf_scores


def list_all_schemas():
	return [("Sig_Conf", lambda s,fc,c : math.log(s,.05) * c),
			("Fold_Conf", lambda s,fc,c : math.fabs(fc) * c),
			("Sig_Fold_Conf", lambda s,fc,c : math.log(s,.05) * math.fabs(fc) * c),
			("Fold_Conf_Preserve_Sign", lambda s,fc,c : fc * c),
			("Sig_Fold_Conf_Preserve_Sign", lambda s,fc,c : math.log(s,.05) * fc * c)]

if __name__ == '__main__':
	main()