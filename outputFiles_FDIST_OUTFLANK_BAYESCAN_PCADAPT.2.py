#!/usr/bin/env python
import sys
import argparse
def parse_args() :
	parser = argparse.ArgumentParser(description='By default for reading the result after outflank, pcadapt and fdist2 , for flk : type -fyes and for bayescan : -byes \n Read in : 1Tfdistresult.txt to 100Tfdistresult.txt, in 1FLKseuil975.txt or 995 & FLKtests.txt, in bayescan_fst.txt, in resultatPCadapt.txt, in resultatOutflank.txt')
	parser.add_argument('-q', '--quanti', default='1', help='get the number of quantitaive loci', dest='quanti')
	parser.add_argument('-b', '--bayescan', default='yes', help='was bayescan used ?', dest='bayescan')
	parser.add_argument('-o', '--outflank', default='yes', help='was outflank used ?', dest='outflank')
	parser.add_argument('-p', '--pcadapt', default='yes', help='was pcadapt used ?', dest='pcadapt')
	parser.add_argument('-d', '--fdist', default='yes', help='was fdist used ?', dest='fdist')
	args = parser.parse_args()
	return args.quanti, args.flk, args.bayescan, args.outflank, args.pcadapt, args.fdist





def findTrueBayescan(quantiLoci,n): #read Bayescan results
	inputfile="files/"+str(n)+"bayescanOutliers.txt"
	f = open(inputfile,"r")
	true=[]
	ligneOK=0
	for idx, line in enumerate(f) :
		line = line.split()
		if len(line)!=0:
			if line[0]=='$outliers':
				ligneOK=1
			elif ligneOK==1:
				for i in range(1,len(line)):
					true.append(line[i])
				ligneOK=0
	f.close()
	return true

def findTrueFdist(quantiLoci,n): #read FDist results
	inputfile="files/"+str(n)+"Tfdistresult.txt"
	f = open(inputfile,"r")
	true=[]
	for idx, line in enumerate(f) :
		line = line.split()
		if float(line[3])>=0.95:
			true.append(idx+1)
	f.close()
	return true

def findTruePCadapt10(quantiLoci,n): #read pcadapt results
	inputfile="files/"+str(n)+"resultatPCadapt10.txt"
	f = open(inputfile,"r")
	true=[]
	for idx, line in enumerate(f) :
		line = line.split()
		if len(line)!=0:
			if line[1]!="\"outliers\"":
				true.append(line[1])	

	f.close()

	return true
def findTruePCadapt5(quantiLoci,n): #read pcadapt results
	inputfile="files/"+str(n)+"resultatPCadapt5.txt"
	f = open(inputfile,"r")
	true=[]
	for idx, line in enumerate(f) :
		line = line.split()
		if len(line)!=0:
			if line[1]!="\"outliers\"":
				true.append(line[1])	

	f.close()

	return true

def findTrueOUTFLANK(quantiLoci,n): #read Outflank results
	inputfile="files/"+str(n)+"resultatOutflank.txt"
	f = open(inputfile,"r")
	flag=0
	colonne=666
	nbreTrue=0
	true=[]
	for idx, line in enumerate(f) :
		line = line.split()
		if len(line)!=0:
			if flag==0 :
				for i in range(len(line)) :
					if line[i]=='OutlierFlag':
						flag=1
						colonne=i+1
			else :
				if line[colonne]=="TRUE":
					nbreTrue=nbreTrue+1
					true.append(line[0])
	f.close()
	return true

def findTrueOUTFLANKL(quantiLoci,n): #OUTFLANK
	inputfile="files/"+str(n)+"resultatOutflankL.txt"
	f = open(inputfile,"r")
	flag=0
	colonne=666
	nbreTrue=0
	true=[]
	for idx, line in enumerate(f) :
		line = line.split()
		if len(line)!=0:
			if flag==0 :
				for i in range(len(line)) :
					if line[i]=='OutlierFlag':
						flag=1
						colonne=i+1
			else :
				if line[colonne]=="TRUE":
					nbreTrue=nbreTrue+1
					true.append(line[0])
	f.close()
	return true



def location(quanti_loci): # get the index of loci under selection
	quantiIndex=[]
	numero=float(100/(int(quanti_loci)+1))
	for i in range(0,int(quanti_loci)):
		quantiIndex.append(1000+i+1)
	return quantiIndex

def comparaison(quantiIndex,true):
	wellDetected=0
	notDetected=0
	falsePostive=0
	if len(true)!=0 :
		for i in range(0,len(quantiIndex)):
			ok=0
			for j in range(0,len(true)):
				if int(quantiIndex[i])==int(true[j]):
					wellDetected=wellDetected+1
					ok=1
			if ok==0:
				notDetected=notDetected+1
		falsePostive=len(true)-wellDetected
	else:
		notDetected=len(quantiIndex)

	return 	wellDetected,notDetected,falsePostive 

def tot(h,quantiLoci,nom,findtrue):
	wellDetectedM=0
	notDetectedM=0
	falsePostiveM=0
	for i in range(1,101):
		true=findtrue[h](quantiLoci,str(i))
		quantiIndex=location(quantiLoci)
		wellDetected,notDetected,falsePostive=comparaison(quantiIndex,true)
	
		
		wellDetectedM=wellDetectedM+wellDetected
		notDetectedM=notDetectedM+notDetected
		falsePostiveM=falsePostiveM+falsePostive
	if (falsePostiveM+wellDetectedM>0):
		tauxNom=float(falsePostiveM)/float(falsePostiveM+wellDetectedM)
	else:
		tauxNom='nd'
	return float(wellDetectedM)/100,float(notDetectedM)/100,float(falsePostiveM)/100, tauxNom,float(wellDetectedM)/float(quantiLoci*100),float(falsePostiveM)/(100*1000)

def writeResults(f,nom,wellDetectedM,notDetectedM,falsePostiveM,tauxNom,puissance,realFP):
	print >> f,nom, wellDetectedM,notDetectedM,falsePostiveM,tauxNom,puissance,realFP


def main():
	quantiLoci,flk,bayescan,outflank, pcadapt, fdist=parse_args() # get info from users
	quantiLoci=int(quantiLoci)
	softwares=[]
	findtrue=[]
	if (outflank=="yes"):
		softwares.append("Outflankl")
		findtrue.append(findTrueOUTFLANKL)#read results from outflank
	if (fdist=="yes"):
		softwares.append("fdist")
		findtrue.append(findTrueFdist)#read results from fdist
	if (pcadapt=="yes"):
		softwares.append("PCadapt5")
		findtrue.append(findTruePCadapt5)#read results from pcadapt
		softwares.append("PCadapt10")
		findtrue.append(findTruePCadapt10)
	if (bayescan=="yes"):
		softwares.append("Bayescan")
		findtrue.append(findTrueBayescan)#read results from bayescan
	print softwares," -> resultatTOT.txt"
	f = open("resultatTOT.txt","w")
	for i in range(0,len(softwares)):
		wellDetectedM,notDetectedM,falsePostiveM,tauxNom,puissance,realFP=tot(i,quantiLoci,softwares[i],findtrue)#calculate fasle postive true positive ...
		writeResults(f,softwares[i],wellDetectedM,notDetectedM,falsePostiveM,tauxNom,puissance,realFP) #write results
	f.close()
	for i in range(0,len(softwares)): #writes results software per software
		file1=softwares[i]+"res.txt"
		f = open(file1,"w")
		wellDetectedM,notDetectedM,falsePostiveM,tauxNom,puissance,realFP=tot(i,quantiLoci,softwares[i],findtrue)
		print >> f,wellDetectedM,notDetectedM,falsePostiveM,tauxNom,puissance,realFP
		f.close()

if __name__ == "__main__":
	main()

