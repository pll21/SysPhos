#!/usr/bin/python

import SeqConvert as seq
import Cleanup as clean
import sys
import operator
import math
import os
import fnmatch
from multiprocessing import Lock, Process, Queue, current_process
import argparse

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
	clean.clean_existing(savedir)
	write_kinase_and_peptide_scores(data_location,savedir)
	permutation_test(data_location,savedir,num_iterations)
	


def write_kinase_and_peptide_scores(infile,outdir):
	if(not os.path.exists(outdir)): os.makedirs(outdir)
	seq_conv = seq.SeqConvert(infile,False)
	netphorest_frame = seq_conv.get_trimmed_netphorest_frame()
	
	print("\tCalculating Scores")
	for schema in list_all_schemas():
		kinase_outfile = "%s/kinase_scores_%s.txt" % (outdir,schema[0])
		kinase_writer = open(kinase_outfile,'wb')

		peptide_outfile = "%s/peptide_scores_%s.txt" % (outdir,schema[0])
		peptide_writer = open(peptide_outfile,'wb')

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
	print("\tCalculating Scores")
	for schema in list_all_schemas():
		kinase_outfile = "%s/kinase_scores_%s.txt" % (outdir,schema[0])
		kinase_writer = open(kinase_outfile,'wb')

		kinase_scores = compute_kinase_scores(netphorest_frame,schema[1])
		sorted_scores = sorted(kinase_scores.items(),key=operator.itemgetter(1))

		for score in sorted_scores:
			kinase_prediction = str(score[0])
			kinase_score = str(score[1])

			kinase_writer.write("%s\t%s\n" % (kinase_prediction,kinase_score))

	clean.clean_all(outdir)


def permutation_test(infile,outdir,num_iterations):
	print("Now performing permutation testing")
	if(not os.path.exists(outdir)): os.makedirs(outdir)
	seq_conv = seq.SeqConvert(infile,True)
	netphorest_frame = seq_conv.get_trimmed_netphorest_frame()

	#workers = 5
	#work_queue = Queue()
	#done_queue = Queue()
	#processes = []

	#for perm_number in xrange(0,num_iterations):
	#	work_queue.put(perm_number)

	#for w in xrange(workers):
	#	p = Process(target=worker, args=(infile,outdir,work_queue, done_queue))
	#	p.start()
	#	processes.append(p)
	#	work_queue.put('STOP')

	#for p in processes:
	#	p.join()


	#done_queue.put('STOP')

	for x in xrange(0,num_iterations):
		print("Working on permutation #%s" % str(x))

		write_kinase_scores(infile,"%s/Permutation%s" % (outdir,str(x)))
	if(num_iterations > 0):
		merge_permutations(outdir)

def worker(infile,outdir,work_queue,done_queue):
	try:
		for perm_num in iter(work_queue.get,'STOP'):
			write_kinase_scores(infile,"%s/Permutation%s" % (outdir,str(perm_num)))
			done_queue.put("Success: %s" % str(perm_num))
	except Exception, e:
		done_queue.put("Failure: %s" % str(perm_num))
	return True



def merge_permutations(outdir):
	#Get list of directories in outdir that start with Permutation
	permutation_directories = [x[0] for x in os.walk(outdir) if x[0][x[0].rfind("/") + 1:].startswith("Permutation")]
	#Get list of text file names in the first Permutation directory
	text_file_names = [y for y in [x[2] for x in os.walk(outdir)][1] if y.endswith(".txt")]
	#Create dictionary that maps text file name to a dictionary mapping kinase names to scores
	text_file_dict = {file_name : {} for file_name in text_file_names}
	for perm_dir in permutation_directories:
		for text_file_name in text_file_names:
			score_dict = text_file_dict[text_file_name]
			reader = open("%s/%s" % (perm_dir,text_file_name),"rb")
			for line in reader.readlines():
				kinase_name,kinase_score = tuple(line.strip().split("\t"))
				if(kinase_name in score_dict):
					score_dict[kinase_name].append(kinase_score)
				else:
					score_dict[kinase_name] = [kinase_score]
			text_file_dict[text_file_name] = score_dict
			reader.close()
	#For each entry in text_file_dict
	results_dir = "%s/Permutation_Results" % outdir
	if(not os.path.exists(results_dir)): os.makedirs(results_dir)
	for filename, score_dict in text_file_dict.iteritems():
		writer = open("%s/Permutation_Results/%s" % (outdir,filename),'wb')
		for kinase_name,kinase_scores in score_dict.iteritems():
			writer.write("%s%s\n" % (kinase_name,"\t" + "\t".join(kinase_scores)))
	clean.clean_directories(permutation_directories)

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