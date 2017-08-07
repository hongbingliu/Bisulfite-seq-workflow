#!/usr/bin/env python2

import sys, getopt
import re
import os
import time
import pandas as pd
import numpy as np


def get_option():
	opts, args = getopt.getopt(sys.argv[1:], "hi:o:f:s:d:")
	input_file = "" #genes info file ( ps : split by '\t' )
	output_file = "" #output name
	ref = "" #gtf
	bismark = "" #bismark file
	divide = 1
	h = ""
	for op, value in opts:
		if op == "-i":
			input_file = value
		elif op == "-o":
			output_file = value
		elif op == "-f":
			ref = value
		elif op == "-s":
			bismark = value
		elif op == "-d":
			divide = value
		elif op == "-h":
			h = 'useages:\nBisulfite-seq data processing workflow.\n-i : genes info (split by tab)\n-o : output name\n-f : gene annotation file (only support gtf file)\n-s : bismark file\n-d : divide genes into <int> (default : 1)'
	return input_file,output_file,ref,bismark,divide,h

def sort_sam(bismark):
	time_start = time.time()
	print ("----------------- first step : sort bismark file -------------------")
	wk_dir = bismark + "-sort/"
	if os.path.exists(wk_dir):
		print ("----------------- bismark sorted file exist .. -------------------")
		print ("time: " + str (time.time()-time_start))
		return 1

	os.makedirs(wk_dir + "+")
	os.makedirs(wk_dir + "-")
	cmd = "cut -f2,3,4 " + bismark + "|sort -nk 3 | awk -F '" + "\\" + "t' '{if(NF>1){print $3 > " + '"' + wk_dir + '"' + "$1" + '"/"' + "$2" + '".txt"' + "}}' "
	print (cmd)
	
	os.system(cmd)

	print ("----------------- sort bismark file okay .. -------------------")
	print ("time: " + str (time.time()-time_start))
					
def genes_info_processing(input_file,ref,divide):
	time_start = time.time()
	print ("----------------- second step : gene info processing -------------------")
	
	wk_dir = input_file + "-sort"
	if os.path.exists(wk_dir):
		print ("----------------- gene info file exist .. -------------------")
		print ("time: " + str (time.time()-time_start))
		return 1

	os.makedirs(wk_dir)
	if divide:
		cmd = "sort -nrk 2 " + input_file + " |awk '{if($2>=1){print $0}}' > " + wk_dir + "/gene-sort.txt"
		print (cmd)
		os.system(cmd)
		f = open(wk_dir + "/gene-sort.txt")
		total = len(f.readlines())
		f.close()
		pair = total // int(divide)
		
		print ("divide into " + str(pair) + " genes each pair..")

		#processing gff3/gtf file

		cmd = ("awk '{if($3==" + '"gene"' + "){print $0}}' /copy1/a-reference/1-maize/0-all-genome/genes.gtf|awk -F '[\\t" + '"' + "]' '{print $1,$4,$5,$7,$10}' > " + wk_dir + "/gtf-sort.txt")
		print (cmd)
		os.system(cmd)

		#divide genes list..
		
		print ("divide genes list..")
		count = 1
		di = 1
		with open(wk_dir + "/gene-pair-" + str(di), 'a') as f:
			with open(wk_dir + "/gene-sort.txt") as fh:
				for i in fh:
					if count == pair:
						if di == int(divide):
							out = os.popen("grep '" + i.split("\t")[0] + "' " + wk_dir + "/gtf-sort.txt").read()
							f.write(out)
						else:
							out = os.popen("grep '" + i.split("\t")[0] + "' " + wk_dir + "/gtf-sort.txt").read()
							f.write(out)
							f.close()
							count = 0
							di += 1
							f = open(wk_dir + "/gene-pair-" + str(di), 'a')
					else:
						out = os.popen("grep '" + i.split("\t")[0] + "' " + wk_dir + "/gtf-sort.txt").read()
						f.write(out)
					count += 1

		print ("divide genes list okay..")

		#divide into different chrs..

		print ("divide into different chrs..")
		for i in range(1, int(divide) + 1):
			os.makedirs(wk_dir + "/pair-" + str(i) + ".+")
			os.makedirs(wk_dir + "/pair-" + str(i) + ".-")
			
			cmd = "sort -nk 2 " + wk_dir + "/gene-pair-" + str(i) + "|awk '{print $2" + '"\\t"' + '$3 > "' + wk_dir + "/pair-" + str(i) + '."$4"/' + '"$1".txt"' + "}'"
			print (cmd)
			os.system(cmd)

		print ("divide into different chrs okay..")
		print ("----------------- gene info processing okay .. -------------------")
		print ("time: " + str (time.time()-time_start))

def mapping_reads(input_file,bismark,output_file,divide):
	time_start = time.time()
	print ("----------------- third step : mapping reads -------------------")
	po_dir = input_file + "-sort/"
	read_dir = bismark + "-sort/"
	out_dir = output_file + "/"
	os.makedirs(output_file)
	
	for n in range(int(divide)):
		os.makedirs(out_dir + "pair-" + str(n + 1) + "/TSS")
		os.makedirs(out_dir + "pair-" + str(n + 1) + "/TTS")
		os.makedirs(out_dir + "pair-" + str(n + 1) + "/body")

	reads_file_list_po = os.popen("ls " + read_dir + "+/").read()[:-1].split("\n")
	reads_file_list_ne = os.popen("ls " + read_dir + "-/").read()[:-1].split("\n")

	for i in reads_file_list_po:
		diff_list_map(i,divide,po_dir,read_dir + "+/",out_dir)

	for i in reads_file_list_ne:
		diff_list_map(i,divide,po_dir,read_dir + "-/",out_dir)

	print ("okay..")
	print ("----------------- mapping reads okay .. -------------------")
	print ("time: " + str (time.time()-time_start))

def calculate_frequency(input_file,bismark,output_file,divide):
	time_start = time.time()
	print ("----------------- fourth step : calculating frequency in each bin -------------------")
	po_dir = input_file + "-sort/"
	read_dir = bismark + "-sort/"
	out_dir = output_file + "/"
	os.makedirs(out_dir + "TSS")
	os.makedirs(out_dir + "TTS")
	os.makedirs(out_dir + "body")
	
	# combine mapping results

	for n in range(int(divide)):
		body = out_dir + "pair-" + str(n + 1) + "/body/"
		TSS = out_dir + "pair-" + str(n + 1) + "/TSS/"
		TTS = out_dir + "pair-" + str(n + 1) + "/TTS/"
		
		# body

		da_unmethy = pd.read_table(body + "-.txt",header = None)
		da_unmethy = da_unmethy[0]
		counts = da_unmethy.value_counts()
		unmethy = {}
		for each in range(len(counts.index)):
			unmethy[str(int(counts.index[each]))] = str((counts.values[each]))
		#print (unmethy)

		da_methy = pd.read_table(body + "+.txt",header = None)
		da_methy = da_methy[0]
		counts = da_methy.value_counts()
		methy = {}
		for each in range(len(counts.index)):
			methy[str(int(counts.index[each]))] = str((counts.values[each]))
		#print (methy)

		f = open(out_dir + "body/pair-" + str(n + 1) + ".txt", 'w')
		for i in range(50):
			i = str(i)
			if i not in unmethy:
				unmethy[i] = 0
			if i not in methy:
				methy[i] = 0
			f.write((i + "\t" + str(int(methy[i])/(int(methy[i]) + int(unmethy[i])) * 100) + "\t" + "expression level " + str(n + 1) + "\n"))
		f.close()

		# TSS
			
		da_unmethy = pd.read_table(TSS + "-.txt",header = None)
		da_unmethy = da_unmethy[0]
		counts = da_unmethy.value_counts()
		unmethy = {}
		for each in range(len(counts.index)):
			unmethy[str(int(counts.index[each]))] = str((counts.values[each]))
	
		da_methy = pd.read_table(TSS + "+.txt",header = None)
		da_methy = da_methy[0]
		counts = da_methy.value_counts()
		methy = {}
		for each in range(len(counts.index)):
			methy[str(int(counts.index[each]))] = str((counts.values[each]))
		
		f = open(out_dir + "TSS/pair-" + str(n + 1) + ".txt", 'w')
		for i in range(-50,0):
			i = str(i)
			if i not in unmethy: 
				unmethy[i] = 0
			if i not in methy: 
				methy[i] = 0
			f.write((i + "\t" + str(int(methy[i])/(int(methy[i]) + int(unmethy[i])) * 100) + "\t" + "expression level " + str(n + 1) + "\n"))
		f.close()		

		# TTS
			
		da_unmethy = pd.read_table(TTS + "-.txt",header = None)
		da_unmethy = da_unmethy[0]
		counts = da_unmethy.value_counts()
		unmethy = {}
		for each in range(len(counts.index)):
			unmethy[str(int(counts.index[each]))] = str((counts.values[each]))
	
		da_methy = pd.read_table(TTS + "+.txt",header = None)
		da_methy = da_methy[0]
		counts = da_methy.value_counts()
		methy = {}
		for each in range(len(counts.index)):
			methy[str(int(counts.index[each]))] = str((counts.values[each]))
		
		f = open(out_dir + "TTS/pair-" + str(n + 1) + ".txt", 'w')
		for i in range(50,100):
			i = str(i)
			if i not in unmethy: 
				unmethy[i] = 0
			if i not in methy: 
				methy[i] = 0
			f.write((i + "\t" + str(int(methy[i])/(int(methy[i]) + int(unmethy[i])) * 100) + "\t" + "expression level " + str(n + 1) + "\n"))
		f.close()		

	cmd = "cat " + out_dir + "*/*.txt > " + out_dir	+ "all.result.txt"
	print (cmd)
	os.system(cmd)
	os.system("sed -i '1 iposition\tmethylation-level\tgroup' " + out_dir + "all.result.txt")

	print ("----------------- calculating frequency in each bin okay .. -------------------")

def diff_list_map(i,divide,po_dir,read_dir,out_dir):
		ismeth = read_dir[-2]
		print ("mapping chr " + str(i) + " ...")
		for n in range(int(divide)):

			print ("\tmapping to genes pair " + str(n + 1) + " ...")
			positive = po_dir + "pair-" + str(n + 1) + ".+/"
			negative = po_dir + "pair-" + str(n + 1) + ".-/"
			if i in os.popen("ls " + positive).read()[:-1].split("\n"):
				with open(read_dir + i) as f:#open reads file
					gi_list = open(positive + i, 'r').read()[:-1].split('\n')
					print (gi_list[:50])
					for each in f:
						read_po = int(each[:-1])
						if read_po < ( int(gi_list[0].split('\t')[0]) - 2000 ):
							continue
						elif read_po < ( int(gi_list[0].split('\t')[1]) + 2000 ):
							calculate_bin_positive(read_po, int(gi_list[0].split('\t')[0]), int(gi_list[0].split('\t')[1]), out_dir + "pair-" + str(n + 1), ismeth)
							for u in gi_list[1:]:
								if read_po < ( int(u.split('\t')[0]) - 2000 ):
									break
								elif read_po < ( int(u.split('\t')[1]) + 2000 ):
									calculate_bin_positive(read_po, int(u.split('\t')[0]), int(u.split('\t')[1]), out_dir + "pair-" + str(n + 1), ismeth)

						else:
							try:
								while read_po >= ( int(gi_list[0].split('\t')[1]) + 2000 ):
									gi_list.pop(0)
							except:
								break

							for u in gi_list:
								if read_po < ( int(u.split('\t')[0]) - 2000 ):
									break
								elif read_po < ( int(u.split('\t')[1]) + 2000 ):
									calculate_bin_positive(read_po, int(u.split('\t')[0]), int(u.split('\t')[1]), out_dir + "pair-" + str(n + 1), ismeth)
				print ("\t......positive is okay..")
				
			if i in os.popen("ls " + negative).read()[:-1].split("\n"):
				with open(read_dir + i) as f:#open reads file
					gi_list = open(negative + i, 'r').read()[:-1].split('\n')
					print (gi_list[:50])
					for each in f:
						read_po = int(each[:-1])
						if read_po <= ( int(gi_list[0].split('\t')[0]) - 2000 ):
							continue
						elif read_po <= ( int(gi_list[0].split('\t')[1]) + 2000 ):
							calculate_bin_negative(read_po, int(gi_list[0].split('\t')[0]), int(gi_list[0].split('\t')[1]), out_dir + "pair-" + str(n + 1), ismeth)
							for u in gi_list[1:]:
								if read_po <= ( int(u.split('\t')[0]) - 2000 ):
									break
								elif read_po <= ( int(u.split('\t')[1]) + 2000 ):
									calculate_bin_negative(read_po, int(u.split('\t')[0]), int(u.split('\t')[1]), out_dir + "pair-" + str(n + 1), ismeth)

						else:
							try:
								while read_po > ( int(gi_list[0].split('\t')[1]) + 2000 ):
									gi_list.pop(0)
							except:
								break

							for u in gi_list:
								if read_po <= ( int(u.split('\t')[0]) - 2000 ):
									break
								elif read_po <= ( int(u.split('\t')[1]) + 2000 ):
									calculate_bin_negative(read_po, int(u.split('\t')[0]), int(u.split('\t')[1]), out_dir + "pair-" + str(n + 1), ismeth)
				print ("\t......negative is okay..")

def calculate_bin_positive(a, b, c, d, e):

	#   calculate reads in which bin (positive)  #
	#   a : reads position                       #
	#   b : TSS position                         #
	#   c : TTS position                         #
	#   d : output dir                           #
	#   e : is it methylation ?                  #
	

	if a < b :
		f = open(d + "/TSS/" + e + ".txt", 'a')
		f.write(str (( a - b ) // 40 ) + "\n")
	elif a < c :
		f = open(d + "/body/" + e + ".txt", 'a')
		f.write(str ((( a - b ) / ( c - b )) // 0.02 ) + "\n")
	else:
		f = open(d + "/TTS/" + e + ".txt", 'a')
		f.write(str (( a - c ) // 40 + 50 ) + "\n")

def calculate_bin_negative(a, b, c, d, e):

	#   calculate reads in which bin (negative)  #
	#   a : reads position                       #
	#   b : TTS position                         #
	#   c : TSS position                         #
	#   d : output dir                           #
	#   e : is it methylation ?                  #
	
	if a <= b :
		f = open(d + "/TTS/" + e + ".txt", 'a')
		f.write(str (( b - a ) // 40 + 50 ) + "\n")
	elif a <= c :
		f = open(d + "/body/" + e + ".txt", 'a')
		f.write(str ((( c - a ) / ( c - b )) // 0.02 ) + "\n")
	else:
		f = open(d + "/TSS/" + e + ".txt", 'a')
		f.write(str (( c - a ) // 40 ) + "\n")

def main(input_file,output_file,ref,bismark,divide):
	sort_sam(bismark)
	genes_info_processing(input_file,ref,divide)
	mapping_reads(input_file,bismark,output_file,divide)
	calculate_frequency(input_file,bismark,output_file,divide)

if __name__ == "__main__":
	time_start = time.time()

	input_file,output_file,ref,bismark,divide,h = get_option()
	if str(h) == "":
		main(input_file,output_file,ref,bismark,divide)
		print ("time: " + str (time.time()-time_start))
	else:
		print (h)
