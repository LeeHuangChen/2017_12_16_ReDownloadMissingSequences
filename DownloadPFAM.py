import Configurations as conf
import util
from cPickle import dump, load
from Bio import Entrez
import os
import sys
import time

#using the input of PFam Mappings, identify the set of proteins that need to be downloaded from NCBI
def identifyProteinSequences(pfamMappingContent):
	proteins=[]
	for line in pfamMappingContent.split("\n"):
		if len(line)>0 and ("PdbResNumStart" not in line): #if the line is not empty and not the header line:
			proteinID=line.split("\t")[0]
			proteins.append(proteinID)


	#clear duplicates
	proteins=list(set(proteins))

	return proteins

#download the sequence from pdb using bio.entrez
def fetchFASTASeqFromPDB(proteinQuery):
	#sets the email associated with this download
	Entrez.email = conf.EntrezEmail
	# #sets the maximum number of search results for the entrez search
	# retmax=3000

	#download the entry information using the pdb id
	handle = Entrez.esearch(db="protein",term=proteinQuery, retmax=1)
	record = Entrez.read(handle)
	handle.close()
	# #testing
	# print "Size:"
	# print sys.getsizeof(record), "bytes"
	# print sys.getsizeof(record)/conf.mbConst, "mb\n"
	
	

	uids=record["IdList"]
	# print uid
	handle = Entrez.efetch(db="protein", id=uids, rettype="fasta", retmode="text")
	result=handle.read().strip()
	handle.close()
	# #testing
	# print "Size:"
	# print sys.getsizeof(result), "bytes"
	# print sys.getsizeof(result)/conf.mbConst, "mb\n"
	return result

#takes care of adding protID to an array of queries such that all queries don't get larger then maxLen character count
def addToProtQueries(protQ, protID, maxLen):
	if len(protQ)==0:
		protQ.append(protID)
		return
	else:
		#i.e. we need to create another entry
		if len(protQ[len(protQ)-1])>maxLen:
			protQ.append(protID)
		else:
			protQ[len(protQ)-1]+= (" OR "+protID)

	
def downloadInOneGo():
	#create the output folder
	util.generateDirectories(conf.outputFolder)

	#input folder
	PFAMFolder=conf.PFAMFolder

	for infile in os.listdir(PFAMFolder):

		#the directory of the output infile is:
		outputDir=os.path.join(conf.outputFolder,infile.replace(conf.PFamExt,conf.outputExt))

		#create the output file
		open (outputDir,"w") 
		
		#identfy the proteins we need to download
		with open(os.path.join(PFAMFolder,infile),"r") as f:
			proteins=identifyProteinSequences(f.read())
		proteinDir=os.path.join(conf.outputFolder,infile.replace(conf.PFamExt,"_proteins.cPickle"))
		with open(proteinDir,"wb") as f:
			dump(proteins, f)
		print "Number of proteins",len(proteins)
		# #convert the proteins into a query
		# proteinQueries=[]
		# maxLen=250
		# for i, protID in enumerate(proteins):
		# 	addToProtQueries(proteinQueries, protID, maxLen)
		# 	# if i>0 and i<len(proteins)-1:
		# 	# 	proteinQuery+=" OR "

		for i, proteinQuery in enumerate(proteins):
			
			print "(",i+1,"/",len(proteins),")", int(i*100/float(len(proteins))),"%"

			results=fetchFASTASeqFromPDB(proteinQuery)
			with open (outputDir,"ab") as f:
				f.write(results)
				f.write("\n\n")
			#print "waiting..."
			time.sleep(.5)

def downloadSteps():
	#input folder
	PFAMFolder=conf.PFAMFolder

	for infile in os.listdir(PFAMFolder):
		proteinDir=os.path.join(conf.outputFolder,infile.replace(conf.PFamExt,"_proteins.cPickle"))
		progDir=os.path.join(conf.outputFolder,infile.replace(conf.PFamExt,"_progress.txt"))
		outputDir=os.path.join(conf.outputFolder,infile.replace(conf.PFamExt,conf.outputExt))
		with open(proteinDir,"rb") as f:
			proteinLeft=load(f)
		startIndex=0
		if os.path.isfile(progDir):
			with open(progDir,"r") as f:
				startIndex=int(f.read().strip())+1
				print "starting at index: ", startIndex
		for i in range(startIndex,len(proteinLeft)):
			protein=proteinLeft[i]
			print "(",i+1,"/",len(proteinLeft),")", int(i*100/float(len(proteinLeft))),"%"

			results=fetchFASTASeqFromPDB(protein)
			with open (outputDir,"ab") as f:
				f.write(">"+protein+":"+results.replace(">",""))
				f.write("\n\n")
			with open(progDir, "w") as f:
				f.write(str(i))
			#print "waiting..."
			#time.sleep(.12)


def main():
	downloadSteps()

		
			





if __name__=="__main__":
	main()