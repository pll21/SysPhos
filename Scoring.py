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

#./Scoring.py Validation
#./Scoring.py Validation/459_Valid_HuangPKA/ValData33noise.txt Results_Folder

def main():

	parser = argparse.ArgumentParser(description='A phosphoproteome-wide kinase inference algorithm.')
	parser.add_argument('data_location', metavar='d', type=str,help='The filename or directory name where the data is stored.')
	parser.add_argument('--permutations',type=int,help='The number of permutations to run.')

	args = parser.parse_args()

	data_location = args.data_location
	num_iterations = args.permutations if args.permutations else 0

	#Handle file case
	if(os.path.isfile(data_location)):
		generate_all_results(data_location,num_iterations)
	#Handle directory case
	elif(os.path.isdir(data_location)):
		for root, dirnames, filenames in os.walk(data_location):
			for filename in fnmatch.filter(filenames, '*.txt'):
				generate_all_results(os.path.join(root,filename),num_iterations)
	else:
		print("SysPhos could not find file or directory")

def generate_all_results(data_location,num_iterations):
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

	seq_conv = seq.SeqConvert(data_location)
	netphorest_frame = seq_conv.get_trimmed_netphorest_frame()

	write_kinase_scores(kinase_dir,netphorest_frame)
	write_peptide_scores(peptide_dir,netphorest_frame)
	write_permutation_scores(permutation_dir,netphorest_frame,num_iterations)
	clean.clean_all(savedir)


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
		peptide_scores = compute_peptide_scores(netphorest_frame,schema[1])
		sorted_scores = sorted(peptide_scores.items(),key=lambda item: sum(float(peptide) for peptide in item[1].values()))
		for entry in sorted_scores:
			kinase_prediction = str(entry[0])
			kinase_score = str(sum(float(x) for x in entry[1].values()))
			peptide_writer.write("Peptide Scores for %s (Score: %s)\n" % (kinase_prediction,kinase_score))
			for peptide_entry in entry[1].items():
				peptide_writer.write("Peptide: %s\t Score: %s\n" % (peptide_entry[0],peptide_entry[1]))
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
		peptide_num = int(row[0])
		curr_series = row[1]
		curr_score = schema(curr_series['pval'], curr_series['fc'], curr_series['Posterior'])
		prediction = curr_series['Prediction']
		kinase_scores[prediction] = (kinase_scores[prediction] if prediction in kinase_scores else 0) + curr_score

	return kinase_scores

def compute_peptide_scores(netphorest_frame,schema):
	peptide_power_scores = {}

	for row in netphorest_frame.iterrows():
		peptide_num = int(row[0])
		curr_series = row[1]

		curr_score = schema(curr_series['pval'], curr_series['fc'], curr_series['Posterior'])
		prediction = curr_series['Prediction']
		peptide_power_dict = (peptide_power_scores[prediction] if prediction in peptide_power_scores else {})
		peptide_power_dict[peptide_num] = (peptide_power_dict[peptide_num] if peptide_num in peptide_power_dict else 0) + curr_score
		peptide_power_scores[prediction] = peptide_power_dict

	return peptide_power_scores


def list_all_schemas():
	return [("Sig_Conf", lambda s,fc,c : c / s),
			("Fold_Conf", lambda s,fc,c : math.fabs(fc) * c),
			("Sig_Fold_Conf", lambda s,fc,c : math.fabs(fc) * c / s),
			("Fold_Conf_Preserve_Sign", lambda s,fc,c : fc * c),
			("Sig_Fold_Conf_Preserve_Sign", lambda s,fc,c : fc * c / s),
			
			("Sig_Conf_Threshold", lambda s,fc,c : (c / s if s < .05 else 0.0)),
			("Sig_Fold_Conf_Threshold", lambda s,fc,c : (math.fabs(fc) * c / s if s < .05 else 0.0)),
			("Sig_Fold_Conf_Preserve_Sign_Threshold", lambda s,fc,c : (fc * c / s if s < .05 else 0.0))]

if __name__ == '__main__':
	main()