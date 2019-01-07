#!/usr/bin/env python

import os
import pybedtools
import statistics
import argparse
import sys
import shutil
from scipy import stats
import scipy.stats
from collections import Counter
import numpy as np
import pandas as pd
import csv
import itertools
import pprint
##################################################-------FUNCTIONS--------#########################################################
# This module will parse file of multiple outputs of scrmshaw to individual unique files for each training set and statisitcal method used. This will return three methods dictionaries with training sets as their keys
def parse_output(outfile,numparse): 
	# creating three separate dictionaries based on methods used  
	d_hexmcd={}
	d_imm={}
	d_pac={}
	x={''}
	global d_hexmcd_val
	global d_imm_val
	global d_pac_val
	d_hexmcd_val=['']
	d_imm_val=['']
	d_pac_val=['']
	
	with open (outfile) as file:
		rows=(line2.split('\t') for line2 in file)

		for row in rows:
		#based on the 14th column(names of different data sets) and 15th column (statistical method used) of scrmshaw_joined_output file giving values to each of the three method`s dictionaries
			if (row[15]=='hexmcd') and (int(row[16]) <= int(numparse)):
				#print(numparse)
				if row[14] not in d_hexmcd:
					myRow = [] # create a new list to use
					myRow.append(row) # add my new row to our new list
					d_hexmcd[row[14]] = myRow  #create a new entry if it isn't in the dictionary already

				else:
					d_hexmcd.get(row[14]).append(row)
					#count_hexmcd=count_hexmcd+1
			elif (row[15]=='imm')and (int(row[16]) <= int(numparse)):
				if row[14] not in d_imm:
					myRow = []
					myRow.append(row)
					d_imm[row[14]] = myRow
				else:
					d_imm.get(row[14]).append(row)
			elif (row[15]=='pac') and (int(row[16]) <= int(numparse) ):
				if row[14] not in d_pac:
					myRow = []
					myRow.append(row)
					d_pac[row[14]] = myRow
				else:
					d_pac.get(row[14]).append(row)
					#count_pac=count_pac+1

	#calculating number of keys(datasets) each dictionary ends up having		
		for key in d_hexmcd.keys():
			d_hexmcd_val.append(key)
		for key in d_imm.keys():
			d_imm_val.append(key)
		for key in d_pac.keys():
			d_pac_val.append(key)
			
	#creating separate files for each method w.r.t datasets, using the above newly created three dictionaries and moving them to tmp (temporary folder).
	
	#These individual unique files(based on methods and their data sets) will be used as input to perform the functions for getting the values to fill up the outfile file
	for key in d_hexmcd.keys():
		noOflines=len(d_hexmcd[key])
		
		with open(os.path.join(subdirectory,'hexmcd_'+key+'_fullLength.bed'),'w') as outfile:
			for i in d_hexmcd[key]:
				line=str(i)
				line=line.replace('\'','')
				line=line.strip('[')
				line=line.strip(']')
				line=line.strip(':\\n')
				outfile.write(" ".join(line.split()).replace(', ', '\t'))
				outfile.write("\n")
	
	for key in d_imm.keys():
		with open(os.path.join(subdirectory,'imm_'+key+'_fullLength.bed'),'w') as outfile:
			for i in d_imm[key]:
				line=str(i)
				line=line.replace('\'','')
				line=line.strip('[')
				line=line.strip(']')
				line=line.strip(':\\n')
				outfile.write(" ".join(line.split()).replace(', ', '\t'))
				outfile.write("\n")

	for key in d_pac.keys():
		with open(os.path.join(subdirectory,'pac_'+key+'_fullLength.bed'),'w') as outfile:
			for i in d_pac[key]:
				line=str(i)
				line=line.replace('\'','')
				line=line.strip('[')
				line=line.strip(']')
				line=line.strip(':\\n')
				outfile.write(" ".join(line.split()).replace(', ', '\t'))
				outfile.write("\n")
				
	return(d_imm_val,d_hexmcd_val,d_pac_val)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------


#This module will extract user specified number of predictions(or top good hits calculated from top good predictions) from the given training set's predictions

def parse_output_individual(indoutfile,parse,dataset,method): 
	# creating three separate dictionaries based on methods used  

	count=0
	#print("inside function")
		
	#print(parse)
	
	if nameOfmethod == 'SCRMshaw' or nameOfmethod =='SCRMSHAW' or nameOfmethod =='scrmshaw':
	#
		with open(subdirectory+"/"+method+"_"+dataset+"_fullLength.bed") as file2, open(subdirectory+"/"+method+"_"+dataset+"tmp.bed",'w') as tmp:		
			for line2 in file2:
			
				#print("opened file")
				row=line2.split('\t')
				lastCol=len(row)-1
				if int(row[lastCol]) <= parse:
					#print(row[16])
					tmp.write(line2)
					count+=1
				
	else:
		with open(indoutfile) as file2, open(subdirectory+"/"+method+"_"+dataset+"tmp.bed",'w') as tmp:		
			for line2 in file2:
		
				#print("opened file")
				row=line2.split('\t')
				lastCol=len(row)-1
				if int(row[lastCol]) <= parse:
					#print(row[16])
					tmp.write(line2)
					count+=1
				
		
	outfile=open(os.path.join(subdirectory,method+'_'+dataset+'.bed'),'w')
	#print(outfile)
	with open(subdirectory+"/"+method+"_"+dataset+"tmp.bed") as file6:	
		for line6 in file6:
			outfile.write(line6)
	
# 			
	pathThis=subdirectory+"/"+method+"_"+dataset+".bed"
	#return number of lines extracted
	#print(count)
	return count, pathThis
	
#---------------------------------------------------------------------------------------------------------------------------


#This module will find out the score cutoff and amplitude cutoff based on their cutoffs 
def topGoodPrediction(unsorted_indScrOutputfile):
	row=[]
	valuesScore=[]
	valuesAmp=[]
	
	with open(unsorted_indScrOutputfile) as sfile:
		for line in sfile:
			row=line.split('\t')
			valuesScore.append(float(row[4]))
			valuesAmp.append(float(row[3]))
	#print(values)
	# pull out the list from pandas frame
	valuesScore=list(valuesScore)
	valuesScore=sorted(valuesScore,key=float,reverse=True)
	#valuesAmp=sorted(valuesAmp,key=float,reverse=True)
	
	#for scores cutoff
	
	#get coordinates of all the points
	nPointsScore = len(valuesScore)
	allCoordScore = np.vstack((range(nPointsScore), valuesScore)).T
	#np.array([range(nPoints), values])

	# get the first point
	firstPointScore = allCoordScore[0]
	# get vector between first and last point - this is the line
	lineVecScore = allCoordScore[-1] - allCoordScore[0]
	lineVecNormScore = lineVecScore / np.sqrt(np.sum(lineVecScore**2))

	# find the distance from each point to the line:
	# vector between all points and first point
	vecFromFirstScore = allCoordScore - firstPointScore

	# To calculate the distance to the line, we split vecFromFirst into two 
	# components, one that is parallel to the line and one that is perpendicular 
	# Then, we take the norm of the part that is perpendicular to the line and 
	# get the distance.
	# We find the vector parallel to the line by projecting vecFromFirst onto 
	# the line. The perpendicular vector is vecFromFirst - vecFromFirstParallel
	# We project vecFromFirst by taking the scalar product of the vector with 
	# the unit vector that points in the direction of the line (this gives us 
	# the length of the projection of vecFromFirst onto the line). If we 
	# multiply the scalar product by the unit vector, we have vecFromFirstParallel
	scalarProductScore = np.sum(vecFromFirstScore * np.matlib.repmat(lineVecNormScore, nPointsScore, 1), axis=1)
	vecFromFirstParallelScore = np.outer(scalarProductScore, lineVecNormScore)
	vecToLineScore = vecFromFirstScore - vecFromFirstParallelScore
	# distance to line is the norm of vecToLine
	distToLineScore = np.sqrt(np.sum(vecToLineScore ** 2, axis=1))
	# knee/elbow is the point with max distance value
	idxOfBestPointScore = np.argmax(distToLineScore)
	
	#print("value of score at score elbow")
	#print(valuesScore[idxOfBestPointScore])
	scoreAtElbow=valuesScore[idxOfBestPointScore]
	
	#for amplitude cutoff
	
	valuesAmp=list(valuesAmp)
	#print(valuesAmp)
	#get coordinates of all the points
	nPointsAmp = len(valuesAmp)
	allCoordAmp = np.vstack((range(nPointsAmp), valuesAmp)).T
	#np.array([range(nPoints), values])

	# get the first point
	firstPointAmp = allCoordAmp[0]
	# get vector between first and last point - this is the line
	lineVecAmp = allCoordAmp[-1] - allCoordAmp[0]
	lineVecNormAmp = lineVecAmp / np.sqrt(np.sum(lineVecAmp**2))

	# find the distance from each point to the line:
	# vector between all points and first point
	vecFromFirstAmp = allCoordAmp - firstPointAmp
	scalarProductAmp = np.sum(vecFromFirstAmp * np.matlib.repmat(lineVecNormAmp, nPointsAmp, 1), axis=1)
	vecFromFirstParallelAmp = np.outer(scalarProductAmp, lineVecNormAmp)
	vecToLineAmp = vecFromFirstAmp - vecFromFirstParallelAmp
	# distance to line is the norm of vecToLine
	distToLineAmp = np.sqrt(np.sum(vecToLineAmp ** 2, axis=1))

	# knee/elbow is the point with max distance value
	idxOfBestPointAmp = np.argmax(distToLineAmp)
	#print("value at 843")
	#print(valuesAmp[idxOfBestPointAmp])
	ampAtElbow=valuesAmp[idxOfBestPointAmp]
	#print(idxOfBestPointScore,idxOfBestPointAmp)
	
	return idxOfBestPointScore,scoreAtElbow,ampAtElbow
	
	
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#This module will discard the peaks below the score cutoff determined by score curve
def extract_below_cutoff_vals(indoutfile,parse,dataset,method,numScore1,numAmp1):
	count=0
	with open(subdirectory+"/"+method+"_"+dataset+"tmp.bed") as fileA, open(subdirectory+"/"+method+"_"+dataset+"tmp2.bed",'w') as tmp2:
		for lineA in fileA:
			rowA=lineA.split('\t')
			#print(line2)
			#if the peak's score value is equal or above the score cutoff only then it will write to file
			if float(rowA[4]) >= numScore1:
				#print(row[16])
				tmp2.write(lineA)
				count+=1
	
	outfile=open(os.path.join(subdirectory,method+'_'+dataset+'.bed'),'w')
	with open(subdirectory+"/"+method+"_"+dataset+"tmp2.bed") as fileB:	
		for lineB in fileB:
			outfile.write(lineB)
		
	pathThis2=subdirectory+"/"+method+"_"+dataset+".bed"
	print("now")
	print(count)
	return count, pathThis2
#----------------------------------------------------------------------------------------------------------------------------


#################################################################--------------Main function----------------------####################################################################
def main():
	global subdirectory
	global d
	global d2
	global countd,countd_random
	global tsetBedOrNot
	global scrmshawOrNot
	global patternRecovery
	global t,x
	global nameOfmethod
	global continuous
	subdirectory='tmp' #Temporary directory
	nameOfmethod= 'SCRMshaw'
	t=1 #count of sets done
	global totalNumberOfCrmsKnownToCauseExpression
	totalNumberOfCrmsKnownToCauseExpression=0	
	#command line parsing files from users
	parser=argparse.ArgumentParser() 
	parser.add_argument('-po','--peaksJoinedOutputFile',help='Scrmshaw Output file')
	parser.add_argument('-listTset','--tsetListFormat',help='File of training set CRMs in the list format')
	
	args = parser.parse_args()
	
	peaksJoinedOutputFile=args.peaksJoinedOutputFile
	tsetListFormat=args.tsetListFormat
	
		#creating temporary directory to move stuff there
	my_path=os.getcwd()
	if not os.path.isdir(my_path+'/'+subdirectory):
		os.makedirs(my_path+'/'+subdirectory)
		
	

	#creating dictionary from the file that contain names of all crms against each training set; use to assess training set sensitivity; in the format name of sets being keys and their respective crms as key's values
	path_to_known_crms=os.path.abspath(tsetListFormat) 
	with open(path_to_known_crms) as fin:
		rows = (line.split('\t') for line in fin )
		d = {row[0].strip(':\:'):row[1:] for row in rows }			


	methods=['imm','hexmcd','pac']
	#Creating list to iterate through three methods
	methods_val=[None,None,None]

	#Parsing the output file into separate files for each training set and each method via creating three dictionaries for each method: keys being the names of data sets associated with that method 
	peaksJoinedOutputFile=os.path.abspath(peaksJoinedOutputFile) 
	methods_val[0],methods_val[1],methods_val[2]=parse_output(peaksJoinedOutputFile,35000)

	#This loop is used to iterate through all three method's dictionary:
	for num in range(len(methods)):
		print("Now method:"+methods[num])
		print("sets in this method "+str(methods_val[num]))
		#This loop is used to iterate through all the training sets in the given methods dictionary
		for x in methods_val[num]:
			if x in d.keys():
				#Getting tmp path of methods dictionary file
				for root, dirs, files in os.walk(os.getcwd()):
					for name in files:
						if name==methods[num]+'_'+x+'_fullLength.bed':
							unsorted_indScrmOutputfile=os.path.abspath(os.path.join(root,name))
	
				tab1=methods[num]	
				numParse,numScore,numAmp = topGoodPrediction(unsorted_indScrmOutputfile)
				print(x)
				print("Top pred:" + str(numParse))
				print('amplitude at cutoff'+str(numAmp))
				print('score at cutoff'+str(numScore))
				#extract number of lines of peaks given by the user
				noOflinesExtracted,unsorted_indScrmOutputfile=parse_output_individual(unsorted_indScrmOutputfile,numParse,x,tab1)
				
				
				noOflinesExtracted,unsorted_indScrmOutputfile=extract_below_cutoff_vals(unsorted_indScrmOutputfile,numParse,x,tab1,numScore,numAmp)
				print(noOflinesExtracted)
	
	
	
	
main()	
