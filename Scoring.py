#!/usr/bin/python

import SeqConvert as seq
import Cleanup as clean
import sys
import operator
import math
import os
import fnmatch

#./Scoring.py Validation
#./Scoring.py Validation/459_Valid_HuangPKA/ValData33noise.txt Results_Folder

def main():

	data_location = sys.argv[1]
	num_iterations = (int(sys.argv[2]) if len(sys.argv) > 2 else 0)

	#Handle file case
	if(os.path.isfile(data_location)):
		print("Now working on %s" % data_location)
		savedir = "Results_%s" % data_location
		savedir = savedir[:savedir.rfind(".txt")]
		write_kinase_and_peptide_scores(data_location, savedir)
	#Handle directory case
	elif(os.path.isdir(data_location)):
		for root, dirnames, filenames in os.walk(data_location):
			for filename in fnmatch.filter(filenames, '*.txt'):
				filepath = os.path.join(root, filename)
				savedir = "Results_%s" % filepath[:filepath.rfind(".")]
				print("Now working on %s" % filepath)
				write_kinase_and_peptide_scores(filepath, savedir)
				permutation_test(filepath,savedir,num_iterations)
	else:
		print("SysPhos could not find file or directory")


def write_kinase_and_peptide_scores(infile,outdir):
	if(not os.path.exists(outdir)): os.makedirs(outdir)
	seq_conv = seq.SeqConvert(infile,False)
	netphorest_frame = seq_conv.get_trimmed_netphorest_frame()
	
	print("\tCalculating Scores")
	#print("\t\tOutdirectory: %s" % outdir)
	for schema in list_all_schemas():
		kinase_outfile = "%s/kinase_scores_%s.txt" % (outdir,schema[0])
		kinase_writer = open(kinase_outfile,'wb')

		peptide_outfile = "%s/peptide_scores_%s.txt" % (outdir,schema[0])
		peptide_writer = open(peptide_outfile,'wb')

		#print("\t\tKinase Outfile: %s\tPeptide Outfile: %s\t" % (kinase_outfile,peptide_outfile))

		kinase_scores,peptide_power_scores = compute_kinase_and_peptide_scores(netphorest_frame,schema[1])
		sorted_scores = sorted(kinase_scores.items(),key=operator.itemgetter(1))
		
		for kinase, peptide_power_dict in peptide_power_scores.iteritems():
			peptide_power_scores[kinase] = sorted(peptide_power_dict.items(),key=operator.itemgetter(1))
		
		for score in sorted_scores:
			kinase_prediction = str(score[0])
			kinase_score = str(score[1])

			kinase_writer.write("%s\t%s\n" % (kinase_prediction,kinase_score))
			
			peptide_writer.write("Peptide Powers for %s (Score: %s)\n" % (kinase_prediction,kinase_score))
			for power_entry in peptide_power_scores[kinase_prediction]:
				peptide_name = str(power_entry[0])
				peptide_score = str(power_entry[1])
				peptide_writer.write("\tPeptide: %s\tPower: %s\n" % (peptide_name,peptide_score))
				
	clean.clean_all(outdir)

def write_kinase_scores(infile,outdir):
	if(not os.path.exists(outdir)): os.makedirs(outdir)
	seq_conv = seq.SeqConvert(infile,True)
	netphorest_frame = seq_conv.get_trimmed_netphorest_frame()
	for schema in list_all_schemas():
		kinase_outfile = "%s/kinase_scores_%s.txt" % (outdir,schema[0])
		kinase_writer = open(kinase_outfile,'wb')

		kinase_scores = compute_kinase_scores(netphorest_frame,schema[1])
		sorted_scores = sorted(kinase_scores.items(),key=operator.itemgetter(1))

		for score in sorted_scores:
			kinase_prediction = str(score[0])
			kinase_score = str(score[1])

			kinase_writer.write("%s\t%s\n" % (kinase_prediction,kinase_score))


def permutation_test(infile,outdir,num_iterations):
	for x in xrange(0,num_iterations):
		print("Working on permutation #%s" % str(x))
		#write_kinase_scores(infile, outdir + str(x))
		write_kinase_scores(infile,"%s/Permutation%s" % (outdir,str(x)))

#Returns a dictionary containing the scores for each peptide
#Parameters: infile - file with list of peptides with significance and fold-change values
#Parameters: outdir - directory to save results to
#Parameters: scoring_scemas - list of tuples, where each tuple contains
#	1) Schema Name
#	2) Function that takes confidence, significance, and fold-change as parameters and gives score as result
def compute_kinase_and_peptide_scores(netphorest_frame,schema):
	kinase_scores = {}
	peptide_power_scores = {}

	for row in netphorest_frame.iterrows():
		peptide_num = int(row[0])
		curr_series = row[1]

		curr_score = schema(curr_series['pval'], curr_series['fc'], curr_series['Posterior'])
		prediction = curr_series['Prediction']
		
		kinase_scores[prediction] = (kinase_scores[prediction] if prediction in kinase_scores else 0) + curr_score
		
		peptide_power_dict = (peptide_power_scores[prediction] if prediction in peptide_power_scores else {})
		peptide_power_dict[peptide_num] = (peptide_power_dict[peptide_num] if peptide_num in peptide_power_dict else 0) + curr_score
		peptide_power_scores[prediction] = peptide_power_dict

	return kinase_scores, peptide_power_scores

def compute_kinase_scores(netphorest_frame,schema):
	kinase_scores = {}
	
	for row in netphorest_frame.iterrows():
		peptide_num = int(row[0])
		curr_series = row[1]
		curr_score = schema(curr_series['pval'], curr_series['fc'], curr_series['Posterior'])
		prediction = curr_series['Prediction']
		kinase_scores[prediction] = (kinase_scores[prediction] if prediction in kinase_scores else 0) + curr_score

	return kinase_scores

def list_all_schemas():
	return [("Sig_Conf", lambda s,fc,c : c / s),
			("Fold_Conf", lambda s,fc,c : math.fabs(fc) * c),
			("Sig_Fold_Conf", lambda s,fc,c : math.fabs(fc) * c / s),
			("Fold_Conf_Preserve_Sign", lambda s,fc,c : fc * c),
			("Sig_Fold_Conf_Preserve_Sign", lambda s,fc,c : fc * c / s)]

if __name__ == '__main__':
	main()