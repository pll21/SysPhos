#This module simply moves all results into a given directory and deletes any files that were used in the process for input/output with NetPhorest
import os
import shutil

def clean_all(outdir):
	os.rename(".sequences.fasta", "%s/sequences.fasta" % outdir)
	os.remove(".NetPhorest_Predictions.txt")

def clean_directories(directories):
	for directory in directories:
		shutil.rmtree(directory)