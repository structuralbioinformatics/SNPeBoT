import fimo,subprocess,sys,os,re,argparse,pickle,math
import pandas as pd
import numpy as np
from SNP import SNP
from keras.models import load_model
from sklearn import metrics


##############
#FILES NEEDED#
##############
# The row file for the Escore table ---> /home/pgohl/Work/pgohl/mutations/8mer.txt
# The column file for the Escore table ---> /home/pgohl/Work/pgohl/mutations/EscoreHeader.txt
# The Escore file ---> /sbi/users/interchange/boliva/ModCRE/pbm/CisBP_2019/Escores.txt
# The Cisbp meme files ---> Resources/pwms/*.meme

##############
# all found in Resources folder. Size: 2.9G
##############

###
### Bins for sequence similarities
###

# Dict_0-10.pickle
# Dict_10-20.pickle
# Dict_20-30.pickle
# Dict_30-40.pickle
# Dict_40-50.pickle
# Dict_50-60.pickle
# Dict_60-70.pickle 
# Dict_70-80.pickle
# Dict_80-90.pickle
# Dict_90-100.pickle

###
###
###



####
#### File import
####
with open('motifs.pickle', 'rb') as handle:
	Mdict = pickle.load(handle)





file3 = open("Resources/8mer.txt")
Merlist = file3.read()
Merlist = Merlist.split("\n")

fieldinfo =  open("Resources/EscoreHeader.txt","r")
scores = fieldinfo.read()
Escores = scores.split("\t")

print("loading Escore table")
df = pd.read_csv("Resources/Escores.txt", sep='\t', encoding='utf-8')
print("done")

# load the model from the local file
model = load_model("Resources/CNN_SNP.keras")
####
####
####



#######
#######
#######
"""
I am restructuring the prediction as a class to handle everything on a standalone basis 
"""
#######
#######
#######



class Prediction:
	motifs = Mdict
	"""
	### Input
	TF -> Name of the transcription factor with which to make the prediction
	lineFile -> location of the file storing the SNP input 
	###
	"""
	def __init__(self, TF, lineFile, Fimop_value):
		self._TF = TF
		self._pvalue = Fimop_value
		self._UniqMotifs = Prediction.motifs[self.get_TF().upper()]
		self._inputSNPlocation = lineFile
		self._sequences = self.set_sequences(self.get_InputLocation())
		self._Reference_Motifs = []
		self._Alternate_Motifs = []
		self._ConcordanceList = self.set_ConcordantData()
		self._Predictions = [] 
	"""
	GETTERS
	"""	
	def get_TF(self):
		"""
		Get the name of the transcription factor
		"""
		return self._TF
            
	def get_Fimo_pvalue(self):
		"""
		Get the pvalue to be used for the fimo filter
		"""
		return self._pvalue
        
	def get_InputLocation(self):
		"""
		Get the location of the input SNP file 
		"""
		return self._inputSNPlocation

	def get_sequences(self):
		"""
		Get the list of SNP objects from the file
		"""
		return self._sequences
	
	def get_UniqMotifs(self):
		"""
		Get the list of motifs that are associated to the given transcription factor 
		"""
		return self._UniqMotifs
            
	def get_ReferenceMotifs(self):
		"""
		Get the list of best matching motifs for the reference sequence
		"""
		return self._Reference_Motifs

	def get_AlternateMotifs(self):
		"""
		Get the list of best matching motifs for the alternate sequence
		"""
		return self._Alternate_Motifs

	def get_PredictionList(self):
		"""
		Get the list of Predictions 
		"""
		return self._Predictions

	def get_ConcordanceList(self):
		"""
		Get the list of FimoBased predictions
		"""
		return self._ConcordanceList

	"""
	FUNCTIONS
	"""
	
	def complementary(self, strand):
		"""
		Given a DNA sequence return the complementary strand 
		"""
		compliment = ""
		for base in strand[::-1]:
			if base == "A":
				compliment = compliment + "T"
			elif base == "T":
				compliment = compliment + "A"
			elif base == "G":
				compliment = compliment + "C"
			elif base == "C":
				compliment = compliment + "G"
		return compliment
	

	def set_sequences(self,LineFile):
		""" Get the individual lines """
		sequences = []
		with open(LineFile, "r") as fd:
			for line in fd:
				sequences.append(SNP(line))
		return sequences

	def find_motifs(self):
		for seqs in self.get_sequences():
			reffasta = ">" + self.get_TF() + "\n" + seqs.get_refSite()
			altfasta = ">" + self.get_TF() + "\n" + seqs.get_altSite()
			tempwriter = open("refFasta.fa","w")
			tempwriter.write(reffasta)
			tempwriter.close()
			tempwriter = open("altFasta.fa","w")
			tempwriter.write(altfasta)
			tempwriter.close()
			WildOptions = {}
			MutOptions = {}
			for motif in self.get_UniqMotifs():
				MotifMeasured = motif.strip() in scores 
				if not MotifMeasured:
					continue
				meme = "Resources/pwms/" + motif.strip() + ".meme"
				low = 1
				print(meme)
				Wildx = fimo.get_fimo_obj(meme,"refFasta.fa")
				for hit in Wildx.get_hits():
					if (hit.get_p_value() < low):
						low = hit.get_p_value()
						WildOptions[motif] = low
				low = 1
				Mutx = fimo.get_fimo_obj(meme,"altFasta.fa")
				for hit in Mutx.get_hits():
					if (hit.get_p_value() < low):
						low = hit.get_p_value()
						MutOptions[motif] = low

			try:
				self._Reference_Motifs.append(min(WildOptions, key=WildOptions.get))
				self._Alternate_Motifs.append(min(MutOptions, key=MutOptions.get))
			except:
				print("There are no motifs")
				exit(0)


	def set_ConcordantData(self):
		self.find_motifs()
		classificationList = []
		for x in range(0,len(self.get_sequences())):
			sequence = self.get_sequences()[x]
			refmotif = self.get_ReferenceMotifs()[x]
			altmotif = self.get_AlternateMotifs()[x]
			sequence.set_motifBindSite(refmotif, altmotif)
			refmeme =  "Resources/pwms/" + refmotif.strip() + ".meme"
			altmeme = "Resources/pwms/" + altmotif.strip() + ".meme"
			Wildx = fimo.get_fimo_obj(refmeme, "RefScoreWindow.fa", fimo_pvalue_threshold= self.get_Fimo_pvalue())
			Mutx = fimo.get_fimo_obj(altmeme, "AltScoreWindow.fa", fimo_pvalue_threshold= self.get_Fimo_pvalue())
			if len(Wildx.get_hits()) == 0 and len(Mutx.get_hits()) == 0:
				classification = "noHIT"
			elif len(Wildx.get_hits()) > len(Mutx.get_hits()):
				classification = "loss"
			elif len(Wildx.get_hits()) < len(Mutx.get_hits()):
				classification = "gain"
			else:
				classification = "no change"
			classificationList.append(classification)
		return classificationList



	def set_Escore(self):
		concords = self.get_ConcordanceList()
		messages = np.array(["gain", "loss", "no change"])
		for iteration in range(0,len(self.get_sequences())):
			scorelist = list()
			altscorelist = list()
			refmatch = list(filter(lambda x: self.get_ReferenceMotifs()[iteration] in x, Escores))
			altmatch = list(filter(lambda x: self.get_AlternateMotifs()[iteration] in x, Escores))
			refindexlist = []
			altindexlist = []
			for i in range(0,8):
				sequence = self.get_sequences()[iteration].get_refSite()[i:(i+8)]
				altsequence = self.get_sequences()[iteration].get_altSite()[i:(i+8)]
				for string in refmatch:
					refindexlist.append(Escores.index(string))
				for string in altmatch:
					altindexlist.append(Escores.index(string))
				index = "NA"
				if sequence.upper() in Merlist:
					index = next(i for i,v in enumerate(Merlist) if v.upper() == sequence.upper())
				else:
					sequence = self.complementary(sequence.upper())
					if sequence.upper() in Merlist:
						index = next(i for i,v in enumerate(Merlist) if v.upper() == sequence.upper())
				if index != "NA":
					scorelist.append(df[refmatch[0]][index-1])
				if altsequence.upper() in Merlist:
					altindex = next(i for i,v in enumerate(Merlist) if v.upper() == altsequence.upper())
				else:
					altsequence = self.complementary(altsequence.upper())
					if altsequence.upper() in Merlist:
						altindex = next(i for i,v in enumerate(Merlist) if v.upper() == altsequence.upper())
				if altindex != "NA":
					altscorelist.append(df[altmatch[0]][altindex-1])
				else:
					altscorelist.append("NA")
			combinedLists = [scorelist, altscorelist]
			refwindow = list()
			altwindow = list()
			for i in range(1,16):
				if(i > 8 ):
					window = 8 - (i - 8)
					altwindow.append(sum(altscorelist[-window:])/len(altscorelist[-window:]))
					refwindow.append(sum(scorelist[-window:])/len(scorelist[-window:]))
				else:
					window = i
					altwindow.append(sum(altscorelist[0:window])/len(altscorelist[0:window]))
					refwindow.append(sum(scorelist[0:window])/len(scorelist[0:window]))
			secondcombinedLists = refwindow
			secondcombinedLists.extend(altwindow)
			#print(refwindow, altwindow)
			secondcombinedLists = np.reshape(np.array(secondcombinedLists), (6,-1) )
			data = np.expand_dims(secondcombinedLists,axis=0)
			print(data)
			prediction = model.predict(data)
			prediction = np.argmax(prediction,axis=-1)
			prediction = np.where(prediction == 0, messages[0], np.where(prediction == 1, messages[1], messages[2]))
			if concords[iteration] == 'noHIT':
				pass
			#elif concords[iteration] == prediction[0]:
			#	self._Predictions.append((prediction[0],iteration))
			else:
				self._Predictions.append((prediction[0],iteration))
			#print(str(combinedLists).replace("[", "").replace("]",""),str(secondcombinedLists).replace("[", "").replace("]","") )


	def get_trueLabel(self, LineFile):
		fd = open(LineFile,"r")
		Truelabel = []
		for line in fd:
			line = line.split("\t")
			labeldict = {True : "gain", False : "loss"}
			if line[8] == "True":
				Truelabel.append(labeldict[int(line[6]) < int(line[7])])
			else:
				Truelabel.append("no change")
		return Truelabel

	def MeasurePerformance(self):
		labels = self.get_trueLabel(self.get_InputLocation())
		actual = []
		predicted = []
		for index in range(0, len(self.get_PredictionList())):
			predicted.append(self.get_PredictionList()[index][0])
			actual.append(labels[self.get_PredictionList()[index][1]])
		confusion_matrix = metrics.confusion_matrix(actual, predicted)
		return confusion_matrix

	def SavePerformance(self,OutputFile):
		cnf_matrix = self.MeasurePerformance()

		FP = cnf_matrix.sum(axis=0) - np.diag(cnf_matrix)  
		FN = cnf_matrix.sum(axis=1) - np.diag(cnf_matrix)
		TP = np.diag(cnf_matrix)
		TN = cnf_matrix.sum() - (FP + FN + TP)
	
		FP = FP.astype(float)
		FN = FN.astype(float)
		TP = TP.astype(float)
		TN = TN.astype(float)
		TPR = TP/(TP+FN)
		FPR = FP/(FP+TN)
		ACC = (TP+TN)/(TP+FP+FN+TN)
		MCC = 0
		for x in range(0,3):
			MCC += ((TP*TN - FP*FN))[x] / math.sqrt(((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))[x])

		MCC = MCC / 3 

		addition = self.get_TF() + ":" + "\t" + "FPR " + str(np.mean(FPR)) + "\t" + "TPR " + str(np.mean(TPR)) + "\t" + str(np.mean(ACC)) + "\t" + "MCC " + str(MCC) + "\n" 
		print(addition)
		with open(OutputFile,"a") as jd:
			jd.write(addition)







