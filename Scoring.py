#!/usr/bin/python

import SeqConvert as seq
import Cleanup as clean
import sys
import operator
import math
import os
import fnmatch

#./Scoring.py Validation

def main():
	#abs_path = os.getcwd()
	#val_dirs = [name for name in os.listdir(sys.argv[1]) if os.path.isdir(os.path.join(sys.argv[1], name))]
	#for val_dir in val_dirs:
	#	files = [x for x in os.listdir(abs_path + "/" + sys.argv[1] + "/" + val_dir) if x.endswith(".txt")]
	#	data_dir = "Validation/%s" % val_dir
	#	for f in files:
	#		data_file = f
	#		inputdir = "%s/%s/%s" % (abs_path,data_dir,data_file)
	#		print("Now working on: %s" % inputdir)
	#		savedir = "%s/Results/%s/%s" % (abs_path,data_dir[data_dir.find("/") + 1:], data_file[:data_file.find(".txt")])
	#		compute_all_scores(inputdir,savedir)
	#		clean.clean_all(savedir)

	data_location = sys.argv[1]
	#Handle file case
	if(os.path.isfile(data_location)):
		print("Now working on %s" % data_location)
		savedir = "Results_%s" % data_location
		compute_score(data_location, savedir)
	#Handle directory case
	elif(os.path.isdir(data_location)):
		for root, dirnames, filenames in os.walk(data_location):
			for filename in fnmatch.filter(filenames, '*.txt'):
				filepath = os.path.join(root, filename)
				savedir = "Results_%s" % filepath[:filepath.rfind("/")]
				print("Now working on %s" % filepath)
				compute_score(filepath, savedir)
	else:
		print("SysPhos could not find file or directory")
	#First we need to get the filepaths of every .txt file in the directory





#Compute the peptide scores for all files in input. Cleanup by saving all files to their correct directories and deleting intermediate files.
#Parameters - input can be either a filename or directory name. 
#	1) If it is a filename we run get_scores() on just that file
#	2) If it is a directory we run get_scores() on every file in that directory
#Parameters - output is the directory to which we are saving our results
#def compute_all_scores(infile,outdir):
	#seq_conv = seq.SeqConvert(infile)
	#netphorest_frame = seq_conv.get_trimmed_netphorest_frame()
	#for schema in list_all_schemas():
	#	scores = get_scores(netphorest_frame,schema[1])
	#	if(not os.path.exists(outdir)): os.makedirs(outdir)
	#	outfile = outdir + "/" + schema[0] + ".txt"
	#	writer = open(outfile,'wb')
	#	sorted_scores = sorted(scores.items(),key=operator.itemgetter(1))
	#	for score in sorted_scores:
	#		writer.write(str(score[0]) + "\t" + str(score[1]) + "\n")

def compute_score(infile,outdir):
	if(not os.path.exists(outdir)): os.makedirs(outdir)

	seq_conv = seq.SeqConvert(infile)
	netphorest_frame = seq_conv.get_trimmed_netphorest_frame()
	for schema in list_all_schemas():
		scores = get_scores(netphorest_frame,schema[1])
		outfile = "%s/%s.txt" % (outdir,schema[0])
		writer = open(outfile,'wb')
		sorted_scores = sorted(scores.items(),key=operator.itemgetter(1))
		for score in sorted_scores:
			writer.write(str(score[0]) + "\t" + str(score[1]) + "\n")
	clean.clean_all(outdir)


#Returns a dictionary containing the scores for each peptide
#Parameters: infile - file with list of peptides with significance and fold-change values
#Parameters: outdir - directory to save results to
#Parameters: scoring_scemas - list of tuples, where each tuple contains
#	1) Schema Name
#	2) Function that takes confidence, significance, and fold-change as parameters and gives score as result
def get_scores(netphorest_frame,schema):
	scores = {}
	for row in netphorest_frame.iterrows():
		curr_series = row[1]
		prediction = curr_series['Prediction']
		scores[prediction] = (scores[prediction] if prediction in scores else 0) + schema(curr_series['pval'], curr_series['fc'], curr_series['Posterior'])
	return scores

def list_all_schemas():
	return [("Sig_Conf", lambda s,fc,c : math.log(s,.05) * c),
			("Fold_Conf", lambda s,fc,c : math.fabs(fc) * c),
			("Sig_Fold_Conf", lambda s,fc,c : math.log(s,.05) * math.fabs(fc) * c),
			("Fold_Conf_Preserve_Sign", lambda s,fc,c : fc * c),
			("Sig_Fold_Conf_Preserve_Sign", lambda s,fc,c : math.log(s,.05) * fc * c)]

if __name__ == '__main__':
	main()