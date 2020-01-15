#!/usr/bin/env python
import sys
import argparse
from random import randint

def parse_args() :
	parser = argparse.ArgumentParser(description='Takes an fstat file as input to generate the allele frequencies by populations, and input files in Outflank, Bayescan,  fdist2 and pcadapt hopefully')
	parser.add_argument('-i', '--input', required=True, help='The basename of input files. We expect to find a tped and a tfam files with this basename.', dest='input_basename')
	parser.add_argument('-o', '--output', required=True, help='The basename of the output files.', dest='output_basename')
	parser.add_argument('-s', '--samplesize', default='all', help='Number of individuals sampled per population. If all individuals are to be analysed type "all".', dest='samplesize')
	parser.add_argument('-k', '--keep', default='all', help='how many populations are to be kept?', dest='keep')
	parser.add_argument('-a', '--add', default='no', help='get the quantitaive trait', dest='add')
	args = parser.parse_args()
	return args.input_basename, args.output_basename,args.samplesize,args.keep,args.add

def sentences():# print info
	print
	print 'Takes a fstat file, and generates the allele frequencies per population, input files for Outflank, Bayescan, fdist2 and pcadapt'
	print 'Please indicates the input file name (after -i), the output file name (after -o) and the number of individuals you want to be sampled per population (after -s)'
	print 'Does not handle missing data'
	print



def fst(matrix012,output_basename,sample,nbreloci,nbrepop,nbreIndPop):
	filename1=output_basename+'fstat.dat'
	f = open(filename1,"w")
	print >> f, nbrepop,nbreloci,"2 1"
	nom="trait-1_locus-"	
	for i in range(0,nbreloci): 
		print >> f, nom+str(i+1)
	h=0
	for x in range(0,len(nbreIndPop)):
		for y in range(0,nbreIndPop[x]):
			print >> f, x+1,
			for u in range(0,nbreloci):
				if matrix012[h][u]==0:
					print >> f, 11,
				if matrix012[h][u]==1:
					print >> f, 12,
				if matrix012[h][u]==2:
					print >> f, 22,
			print >> f
			h=h+1	
	print "pegas input file is",filename1
	f.close

def info(nbreloci,nbrepop,matrix012,nbreIndPop,sample,keep): #print info
	print "In the studied sample, there are",nbreloci,"loci."
	if sample=='all':
		print "There are respectively",
		for i in range(0,nbrepop):
			print nbreIndPop[i],
		print "individuals in each population",
	else:
		print "You chose a sample of", int(sample),"individuals per population",
	if keep=='all':
		print "and you chose to keep all",nbrepop,"populations."
	else:
		print "and you chose to keep",keep,"populations."
	print

def selectInd(matrixtot,nbreloci,nbreIndPop,nbre,keep): #nbre=sample #if the users ask for a sample of the individuals
	matrixtotsampled=[]
	start=0
	tot=0
	nbreIndPop2=[]
	if keep!="all":
		for i in range(0,int(keep)):
			nbreIndPop2.append(nbreIndPop[i])
	else:
		nbreIndPop2=nbreIndPop
	choix=[]
	matrixtotsampled=[]
	for i in range(0,len(nbreIndPop2)):
		for u in range(0,nbre):	
			notok=0
			while notok==0:
				leChiffre=randint(start,start+nbreIndPop2[i]-1)
				if leChiffre not in choix:
     					choix.append(leChiffre)
					notok=1
			matrixtotsampled.append(matrixtot[choix[tot]])
			tot=tot+1

		start=start+nbreIndPop2[i]
	alleleEffectif=[]
	allelefrequence=[]
	nbreInd=0
	for i in range(0,len(nbreIndPop2)):
		alleleEffectif.append([])
		allelefrequence.append([])
		for k in range(0,nbreloci):
			alleleEffectif[i].append(0)
			allelefrequence[i].append(0)
		for u in range(0,nbre):
			for w in range(0,nbreloci):
				alleleEffectif[i][w]=alleleEffectif[i][w]+matrixtotsampled[nbreInd][w]
			nbreInd=nbreInd+1
	for a in range(0,len(nbreIndPop2)):
		for n in range(0,nbreloci) :
				allelefrequence[a][n]=(float(alleleEffectif[a][n])*100)/(nbre*2)
	nbreIndPop=len(nbreIndPop2)*[int(nbre)]
		
	return allelefrequence,nbreloci,matrixtotsampled,nbreIndPop,alleleEffectif



def addQuanti(nbrepop,nbreIndPop,allelefrequence,nbreloci,matrix012,alleleEffectif,input_basename,sample, keep,add): #add loci under selection to SNP matrix
	if add !='no':
		for i in range(0,len(input_basename)-3):
			if input_basename[i:i+4]=='ntrl':
				input_basenameq=input_basename[0:i]+'quanti'+input_basename[i+4:]
				print "\nQuantative trait loci are read in ",input_basenameq
				allelefrequenceq,nbrelociq,nbrepop,matrix012q,nbreIndPop,alleleEffectifq=readEntireFile(input_basenameq, keep)
				matrixtot=[]
				tot=0
				alleleEffectiftot=[]
				allelefrequencetot=[]
				for i in range(0,nbrepop):
					tot=tot+nbreIndPop[i]
					allelefrequencetot.append([])
					alleleEffectiftot.append([])
					for t in range (0,nbreloci):
						alleleEffectiftot[i].append(alleleEffectif[i][t])
						allelefrequencetot[i].append(allelefrequence[i][t])
					for j in range (0,nbrelociq):
						alleleEffectiftot[i].append(alleleEffectifq[i][j])
						allelefrequencetot[i].append(allelefrequenceq[i][j])
				for h in range(0,tot):
					matrixtot.append([])
					for a in range(0,nbreloci):
						matrixtot[h].append(matrix012[h][a])
					for b in range(0,nbrelociq):
						matrixtot[h].append(matrix012q[h][b])
				nbreloci=nbreloci+nbrelociq
	else:
		print "Quantitative trait(s) is/are not taken into account.\n"	
		matrixtot=matrix012
		alleleEffectiftot=alleleEffectif
		allelefrequencetot=allelefrequence
	return matrixtot,alleleEffectiftot,allelefrequencetot,nbreloci

def getAllelefrequencies(allelefrequence,nbreloci,nbrepop,filename): # calculate allele frequencies
	filename=filename+'freq.txt'
	print "Allele frequencies per population are in",filename
	f = open(filename,"w")
	for i in range(0,nbreloci):
		for h in range(0,nbrepop):
			print >> f, allelefrequence[h][i],
		print >> f 
	f.close

def fdistInput(filename,sample,nbreloci,nbrepop,nbreIndPop,alleleEffectif):# writes input file for Fdist
	#print alleleEffectif
	filename=filename+'fdist.txt'
	f = open(filename,"w")
	print >> f, 0
	print >> f, nbrepop
	print >> f, nbreloci
	print >> f
	for t in range(0,nbreloci):
		print >> f,2
		for y in range(0,nbrepop):
			print >> f, alleleEffectif[y][t], nbreIndPop[y]*2-alleleEffectif[y][t]
		print >> f
	print "The fdist input file is",filename
	f.close

def bayescanInput(filename,sample,nbreloci,nbrepop,nbreIndPop,alleleEffectif): # writes input file for bayescan
	filename=filename+'bayescan.txt'
	f = open(filename,"w")
	print >> f, '[loci]={0}\n'.format(nbreloci)
	print >> f, '[populations]={0}\n'.format(nbrepop)
	for y in range(0,nbrepop):
		print >> f, '[pop]={0}\n'.format(y+1),
		for t in range(0,nbreloci):
			print >> f, t+1, nbreIndPop[y]*2, 2, alleleEffectif[y][t], nbreIndPop[y]*2-alleleEffectif[y][t]
		print >> f
	print "The bayescan input file is",filename
	f.close



def outflankInput(matrix012,filename,sample,nbreloci,nbrepop,nbreIndPop) : # writes input file for OUTflank
	filename1=filename+'Outlocus.txt'
	f = open(filename1,"w")
	for i in range(0,nbreloci):
		print >> f, "locus"+str(i+1),
	f.close
	filename2=filename+'OutSNPmat.txt'
	f = open(filename2,"w")
	for h in range(0,len(matrix012)):
			for u in range(0,nbreloci):
				print >> f, matrix012[h][u],
			print >> f	
	f.close
	filename3=filename+'Outpop.txt'
	f = open(filename3,"w")
	for h in range(0,len(nbreIndPop)):
			for u in range(0,nbreIndPop[h]):
				print >> f, "pop"+str(h+1),	
	print "Outflank input files are",filename1,filename2,filename3
	f.close

def pcadaptInput(matrix012,filename,nbreloci) :# writes input file for PCADAPT
	filename=filename+'pcadapt.lfmm'
	f = open(filename,"w")
	for h in range(0,len(matrix012)):
			for u in range(0,nbreloci):
				print >> f, matrix012[h][u],
			print >> f	
	f.close
	print "The pcadapt input file is",filename



def readEntireFile(filename, keep): #read the quantinemo output file entirely
	f = open(filename,"r")
	#print(f.read())
	#f = safe_open()
	premiereLigne = 1
	nbrepop=0
	nbreloci=0
	indexInd=[]
	nbreInd=0
	matrix=[]
	matrix012=[]
	numerodepop=0
	popactuelle=0
	for idx, line in enumerate(f) :
		line = line.split()
		if premiereLigne == 1: # read the first line
			nbrepop=int(line[0]) # get the number of populations
			if keep=='all': #to keep only the populations we want
				toKeep=nbrepop
			else:
				toKeep=int(keep)
				nbrepop=int(keep)
			#print toKeep
			nbreloci=int(line[1]) # get the number of loci
			nbreIndPop=nbrepop*[0]

			alleleEffectif=[] #how many copies of one of the alleles for each locus per population
			for h in range(0,nbrepop):
				alleleEffectif.append([])
				for k in range(0,nbreloci):
					alleleEffectif[h].append(0)
			allelefrequence=[] #allele frequency in each population
			for h in range(0,nbrepop):
				allelefrequence.append([])
				for k in range(0,nbreloci):
					allelefrequence[h].append(0)
			premiereLigne =0
		elif len(line) != 1 and int(line[0])<(toKeep+1): # skip unwanted lines and unwanted populations
			
			indexInd.append(int(line[0])) #to know in which population we are and to count the number of individuals per population 
			if int(line[0]) != popactuelle:
				numerodepop=numerodepop+1
				popactuelle=int(line[0])
			nbreIndPop[numerodepop-1]=nbreIndPop[numerodepop-1]+1

			matrix.append([]) #same genotypes as in the source file
			matrix012.append([]) # genotypes coded as : 0 1 2
			for i in range(1,nbreloci+1) :
				matrix[nbreInd].append(int(line[i]))
				if int(line[i]) == 11:
					matrix012[nbreInd].append(0)
				elif int(line[i]) == 22:
					matrix012[nbreInd].append(2)
					alleleEffectif[numerodepop-1][i-1]=alleleEffectif[numerodepop-1][i-1]+2
				else :
					matrix012[nbreInd].append(1)
					alleleEffectif[numerodepop-1][i-1]=alleleEffectif[numerodepop-1][i-1]+1
			nbreInd=nbreInd+1
	for
 a in range(0,nbrepop):
				for n in range(0,nbreloci) :
					allelefrequence[a][n]=(float(alleleEffectif[a][n])*100)/(nbreIndPop[a]*2)
	f.close
	return allelefrequence,nbreloci,nbrepop, matrix012,nbreIndPop,alleleEffectif

def main():
	print
	ok=1
	if ok==1 :
		sentences()# print info
		input_basename, output_basename,sample, keep, add=parse_args() #get info from user
		allelefrequence,nbreloci,nbrepop,matrix012,nbreIndPop,alleleEffectif=readEntireFile(input_basename, keep)#read the quantinemo output file entirely 
		matrix012,alleleEffectif,allelefrequence,nbreloci=addQuanti(nbrepop,nbreIndPop,allelefrequence,nbreloci,matrix012,alleleEffectif,input_basename,sample, keep,add)#add loci under selection to SNP matrix
		if sample != "all": #if the users ask for a sample of the individuals, select ind 
			allelefrequence,nbreloci,matrix012,nbreIndPop,alleleEffectif=selectInd(matrix012,nbreloci,nbreIndPop,int(sample),keep) 
		info(nbreloci,nbrepop,matrix012,nbreIndPop,sample,keep) #print info 
		
		getAllelefrequencies(allelefrequence,nbreloci,nbrepop,output_basename) # calculate allele frequencies
		outflankInput(matrix012,output_basename,sample,nbreloci,nbrepop,nbreIndPop) # writes input file for OUTflank
		bayescanInput(output_basename,sample,nbreloci,nbrepop,nbreIndPop,alleleEffectif)# writes input file for Bayescan
		fdistInput(output_basename,sample,nbreloci,nbrepop,nbreIndPop,alleleEffectif) # writes input file for fdist
		pcadaptInput(matrix012,output_basename,nbreloci) # writes input file for PCAdapt
		fst(matrix012,output_basename,sample,nbreloci,nbrepop,nbreIndPop)
		if keep!='all': #when there is an outgroup
			keep='all'
			allelefrequence,nbreloci,nbrepop,matrix012,nbreIndPop,alleleEffectif=readEntireFile(input_basename, keep)
			matrix012,alleleEffectif,allelefrequence,nbreloci=addQuanti(nbrepop,nbreIndPop,allelefrequence,nbreloci,matrix012,alleleEffectif,input_basename,sample, keep,add)
			if sample != "all": 
				allelefrequence,nbreloci,matrix012,nbreIndPop,alleleEffectif=selectInd(matrix012,nbreloci,nbreIndPop,int(sample),keep)
			output_basename2=output_basename+"all"
			print "For all populations:",
			getAllelefrequencies(allelefrequence,nbreloci,nbrepop,output_basename2)
		
		


if __name__ == "__main__":

	main()

