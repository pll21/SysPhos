#This module simply moves all results into a given directory and deletes any files that were used in the process for input/output with NetPhorest
import os

def clean_all(savedir):
	os.remove("sequences.fasta")
	os.remove("NetPhorest_Predictions.txt")