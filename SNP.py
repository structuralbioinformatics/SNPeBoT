import subprocess,sys,os,re


class SNP:

	def __init__(self, SNPLine):
		self._SNPLine = SNPLine
		self._referenceSeq = ""
		self._refSite = ""
		self._altSite = ""
		self._alternateSeq = self.seperate_sequences()


	def get_SNPLine(self):
		""" return the SNP input line (tab seperated format) """
		return self._SNPLine

	def seperate_sequences(self):
		""" from the SNPLine set the reference and alternate sequences """
		refseq = self.get_SNPLine().split("\t")[1]
		altseq = self.get_SNPLine().split("\t")[2].strip()
		self._referenceSeq = refseq
		self._refSite = refseq[int(len(refseq)/2 - .5)-7: int(len(refseq)/2 - .5)+8]
		self._altSite = altseq[int(len(altseq)/2 - .5)-7: int(len(altseq)/2 - .5)+8]
		return altseq

	def get_referenceSeq(self):
		""" return the reference sequence """
		return self._referenceSeq

	def get_alternateSeq(self):
		""" return the alternate sequence """
		return self._alternateSeq

	def get_refSite(self):
		""" get the sequence 7 bp up and downstream of the reference allele """
		return self._refSite

	def get_altSite(self):
		""" get the sequence 7 bp up and downstream of the alternate allele """
		return self._altSite


	def set_motifBindSite(self,refmotif,altmotif):
		refmeme = "Resources/pwms/" + refmotif.strip() + ".meme"
		altmeme = "Resources/pwms/" + altmotif.strip() + ".meme"
		callstring = "grep -o " + "\"w= [^ ]*\" " + refmeme
		refMotLen = subprocess.check_output([callstring], shell = True).decode('utf-8').strip()
		refMotLen = int(refMotLen.split(" ")[1])
		ref2match = self.get_referenceSeq()[int(len(self.get_referenceSeq())/2 - .5)-(refMotLen-1): int(len(self.get_referenceSeq())/2 - .5)+(refMotLen)]
		if refmotif != altmotif:
			callstring = "grep -o " + "\"w= [^ ]*\" " + altmeme
			altMotLen = subprocess.check_output([callstring], shell = True).decode('utf-8').strip()
			altMotLen = int(altMotLen.split(" ")[1])
		else:
			altMotLen = refMotLen
		alt2match = self.get_alternateSeq()[int(len(self.get_alternateSeq())/2 - .5)-(altMotLen-1): int(len(self.get_alternateSeq())/2 - .5)+(altMotLen)]
		tempwriter = open("RefScoreWindow.fa","w")
		reffasta =  ">" + "ref" + "\n" + ref2match
		tempwriter.write(reffasta)
		tempwriter.close()
		tempwriter = open("AltScoreWindow.fa","w")
		altfasta =  ">" + "ref" + "\n" + alt2match
		tempwriter.write(altfasta)
		tempwriter.close()

