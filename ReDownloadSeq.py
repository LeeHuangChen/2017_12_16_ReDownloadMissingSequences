import util
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import Configurations as conf
import os
import time

#download the sequence from pdb using bio.entrez
def fetchFASTASeqFromPDB(proteinQuery, retmax):
	#sets the email associated with this download
	Entrez.email = conf.EntrezEmail
	# #sets the maximum number of search results for the entrez search
	# retmax=3000

	#download the entry information using the pdb id
	handle = Entrez.esearch(db="protein",term=proteinQuery, retmax=retmax)
	record = Entrez.read(handle)
	handle.close()
	
	uids=record["IdList"]
	# print uid
	handle = Entrez.efetch(db="protein", id=uids, rettype="fasta", retmode="text", retmax=retmax)
	result=handle.read().strip()
	handle.close()

	return result

def DownloadNewSeq(protID,retmax):
	result=fetchFASTASeqFromPDB(protID,retmax)
	# descs=[]
	# seqs=[]
	# isGoods=[]
	
	
	for record in result.split(">"):
		if len(record.strip())>0:
			desc=record.split("\n")[0]
			seq=record.replace(desc,"").strip()
			seqmerge=seq.replace("\n","").strip()
			#if(seqmerge!=len(seqmerge)*"X" and (protID in desc)):
			if(seqmerge!=len(seqmerge)*"X"):
				return protID+":"+record.strip()+"\n\n"

	if (len(result.split(">"))-1)==retmax:
		return "Retmax"
	else:
		return "Error"

def test():
	#DownloadNewSeq("2NUC")
	#DownloadNewSeq("3KGY")
	print DownloadNewSeq("4P00",10)
	#print 0
def reDownloadSeq():
	util.generateDirectories(conf.outputFolder)
	#input folder
	PFAMFolder=conf.PFAMFolder

	for infile in os.listdir(PFAMFolder):

		#the directory of the output infile is:
		outputDir=os.path.join(conf.outputFolder,infile.replace(conf.PFamExt,conf.outputExt))
		errorDir=os.path.join(conf.outputFolder,infile.replace(conf.PFamExt,"_error.txt"))
		inputDir=os.path.join(conf.PFAMFolder,infile)
		progDir=os.path.join(conf.outputFolder,infile.replace(conf.PFamExt,"_progress.txt"))
		progress=0
		if os.path.isfile(progDir):
			with open(progDir,"r") as f:
				progress=int(f.read().strip)
		#create the output file
		#open(outputDir,"w")
		open(errorDir,"w") 
		records=[]
		for record in SeqIO.parse(open(inputDir,"rU"),"fasta"):
			records.append(record)
		
		for i in range(progress,len(records)):
			record=records[i]
			seq = record.seq
			sid = record.id
			desc = record.description
			protID = desc.split(":")[0]
			seqmerge=str(seq).replace("\n","").strip()
			if seqmerge==len(seqmerge)*"X":#if sequence is all X
				#print seq
				print "downloading",i,"/",len(records), int(i*100/float(len(records))),"%"
				retmax=10
				strOut="Retmax"
				while (strOut=="Retmax" and retmax<1000):
					strOut=DownloadNewSeq(protID,retmax)
					time.sleep(.3)
					retmax=retmax*2
				
				if strOut=="Error" or retmax>=1000:
					print "Cannot find Seq for:",protID,"in",retmax,"downloads"
					with open(errorDir,"a") as f:
						f.write(protID+","+str(retmax)+"\n")

				else:
					with open(outputDir,"a") as f:
						f.write(strOut)
			else:
				with open(outputDir,"a") as f:
					f.write(">"+str(desc)+"\n"+str(seq)+"\n\n")
				#SeqIO.write(SeqRecord(seq,sid,record.name,desc),outputDir,"fasta")
				




def main():
	reDownloadSeq()
	#test()

if __name__ == '__main__':
	main()