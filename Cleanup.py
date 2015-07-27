#This module simply moves all results into a given directory and deletes any files that were used in the process for input/output with NetPhorest
import os
import shutil

def clean_all(outdir):
	os.remove(".sequences.fasta")
	os.remove(".NetPhorest_Predictions.txt")

def clean_directories(directories):
	for directory in directories:
		shutil.rmtree(directory)

def clean_existing(data_location):
	if(os.path.isdir(data_location)):
		shutil.rmtree(data_location)