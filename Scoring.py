import fimo,subprocess,sys,os,re,argparse
import pandas as pd

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

# get the program input 
parser = argparse.ArgumentParser(description='From an input of sequences in Proper format and motifs get the CNN input')
parser.add_argument('-t','--tf' , default=None, dest='transcription_factor', action='store', help='the transcription factor for which to retrieve values.')
parser.add_argument('-s','--seq' , default=None, dest='Sequence_File', action='store', help='the location of the mutated sequences.')
parser.add_argument('-m','--motifs' , default=None, dest='motifs', action='store', help='Motifs for fimo corresponding to the given transcription factor.(string of motifs seperated by newline character.)')
options = parser.parse_args()

# Saving the input to variables
TF = options.transcription_factor
lineFile = options.Sequence_File
uniqueMotifs = options.motifs

# The max amount of SNPs to handle per job
maxInput = 60

file3 = open("Resources/8mer.txt")
Merlist = file3.read()
Merlist = Merlist.split("\n")

fieldinfo =  open("Resources/EscoreHeader.txt","r")
scores = fieldinfo.read()
Escores = scores.split("\t")

print("loading Escore table")
df = pd.read_csv("Resources/Escores.txt", sep='\t', encoding='utf-8')
print("done")

numberofgains = 0
numberoflosses = 0 

with open("DataOutput.csv","w") as fd:
	fd.write("") 

with open("LineOutput.csv","a") as fd:
	fd.write("") 


def complementary(strand):
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

with open(lineFile, "r") as fd:
	number = 0
	for line in fd:
		if number >= maxInput:
			number += 1
			break
			#continue
		number += 1
		callstring = "sed -n " + str(number) + "p " + lineFile + " | cut -d'" + "\t" + "' -f 12 "
		wildSequence = subprocess.check_output([callstring], shell = True).decode('utf-8').strip()
		#callstring = "sed -n " + str(number) + "p " + lineFile + " | cut -d'" + "\t" + "' -f 7 "
		#wildbind = subprocess.check_output([callstring], shell = True).decode('utf-8').strip()
		wd = open("WildTemp.fa", "w")
		addition = ">" + TF + "\n" + wildSequence

		callstring = "sed -n " + str(number) + "p " + lineFile + " | cut -d'" + "\t" + "' -f 13 "
		mutSequence = subprocess.check_output([callstring], shell = True).decode('utf-8').strip()
		#callstring = "sed -n " + str(number) + "p " + lineFile + " | cut -d'" + "\t" + "' -f 8 "
		#mutbind = subprocess.check_output([callstring], shell = True).decode('utf-8').strip()

		md = open("MutTemp.fa", "w")
		wild2match = wildSequence[int(len(wildSequence)/2 - .5)-7: int(len(wildSequence)/2 - .5)+8]
		mut2match = mutSequence[int(len(mutSequence)/2 - .5)-7: int(len(mutSequence)/2 - .5)+8]
		wildfasta = ">" + TF + "\n" + wild2match
		tempwriter = open("wildFasta.fa","w")
		tempwriter.write(wildfasta)
		tempwriter.close()
		mutfasta = ">" + TF + "\n" + mut2match
		tempwriter = open("mutFasta.fa","w")
		tempwriter.write(mutfasta)
		tempwriter.close()
		WildOptions = {}
		MutOptions = {}

		for motif in uniqueMotifs.split("\n"):
			meme = "Resources/pwms/" + motif.strip() + ".meme"
			Wildx = fimo.get_fimo_obj(meme,"wildFasta.fa")
			low = 1
			for hit in Wildx.get_hits():
				if (hit.get_p_value() < low):
					low = hit.get_p_value()
					WildOptions[motif] = low
			low = 1
			Mutx = fimo.get_fimo_obj(meme,"mutFasta.fa")
			for hit in Mutx.get_hits():
				if (hit.get_p_value() < low):
					low = hit.get_p_value()
					MutOptions[motif] = low
		WildMotif = min(WildOptions, key=WildOptions.get)
		WildLine = uniqueMotifs.split("\n").index(WildMotif) + 1 
		MutantMotif = min(MutOptions, key=MutOptions.get)
		MutantLine = uniqueMotifs.split("\n").index(MutantMotif) + 1 
		""" Adding info for posthoc test on concordance """
		meme = "Resources/pwms/" + WildMotif.strip() + ".meme"
		callstring = "grep -o " + "\"w= [^ ]*\" " + meme
		refMotLen = subprocess.check_output([callstring], shell = True).decode('utf-8').strip()
		refMotLen = int(refMotLen.split(" ")[1])
		wild2match = wildSequence[int(len(wildSequence)/2 - .5)-(refMotLen-1): int(len(wildSequence)/2 - .5)+(refMotLen)]
		wildfasta = ">" + "ref" + "\n" + wild2match
		tempwriter = open("wildFasta.fa","w")
		tempwriter.write(wildfasta)
		tempwriter.close()
		altmeme = "Resources/pwms/" + MutantMotif.strip() + ".meme"
		callstring = "grep -o " + "\"w= [^ ]*\" " + altmeme
		altMotLen = subprocess.check_output([callstring], shell = True).decode('utf-8').strip()
		altMotLen = int(altMotLen.split(" ")[1])
		mut2match = mutSequence[int(len(mutSequence)/2 - .5)-(altMotLen-1): int(len(mutSequence)/2 - .5)+(altMotLen)]
		tempwriter = open("mutFasta.fa","w")
		mutfasta = ">" + "alt" + "\n" + mut2match
		tempwriter.write(mutfasta)
		tempwriter.close()
		Wildx = fimo.get_fimo_obj(meme,"wildFasta.fa", fimo_pvalue_threshold=0.05)
		Mutx = fimo.get_fimo_obj(altmeme,"mutFasta.fa", fimo_pvalue_threshold=0.05)
		if len(Wildx.get_hits()) == 0 and len(Mutx.get_hits()) == 0:
			classification = "noHIT" 
		elif len(Wildx.get_hits()) > len(Mutx.get_hits()):
			classification = "loss"
		elif len(Wildx.get_hits()) < len(Mutx.get_hits()):
			classification = "gain" 
		else:
			classification = "no change" 
   

		
		""" Here I have the info on what Motif is selected, WildMotif and MutantMotif, sequence info is in mutSequence and wildSequence
		In order to get non window scores I have to go through all octomers in each sequence and find the corresponding score for this motif."""
		scorelist = list()
		altscorelist = list()
		refmatch = list(filter(lambda x: WildMotif in x, Escores))
		altmatch = list(filter(lambda x: MutantMotif in x, Escores))
		refindexlist = []
		altindexlist = []
		for string in refmatch:
			refindexlist.append(Escores.index(string))
		for string in altmatch:
			altindexlist.append(Escores.index(string))
			
		for i in range(1,9):
			index = "NA"
			altindex = "NA"
			middle = i + 25
			start = middle - 8
			end = start + 8
			sequence = wildSequence[start:end]
			altsequence = mutSequence[start:end]
			""" getting all of the octomer scores for the reference sequence"""
                	
			if sequence.upper() in Merlist:
                		index = next(i for i,v in enumerate(Merlist) if v.upper() == sequence.upper())
			else:
				sequence = complementary(sequence.upper())
				if sequence.upper() in Merlist:
					index = next(i for i,v in enumerate(Merlist) if v.upper() == sequence.upper())
			if index != "NA":
				
				#callstring = "sed -n " + str(index+1) +"p /sbi/users/interchange/boliva/ModCRE/pbm/CisBP_2019/Escores.txt | cut -d'" + "\t"+ "' -f " + str(refindexlist[0]+1)
				try:
					#scorelist.append(float(subprocess.check_output([callstring], shell = True).decode('utf-8').strip()))
					scorelist.append(df[refmatch[0]][index-1])
				except:
					print(index )
					print(sequence.upper())
					exit(0)
			else:
				scorelist.append("NA")
			""" getting all of the octomer scores for the alternate sequence """
			if altsequence.upper() in Merlist:
				altindex = next(i for i,v in enumerate(Merlist) if v.upper() == altsequence.upper())
			else:
				altsequence = complementary(altsequence.upper())
				if altsequence.upper() in Merlist:
					altindex = next(i for i,v in enumerate(Merlist) if v.upper() == altsequence.upper())
			if altindex != "NA":
				#callstring = "sed -n " + str(altindex+1) +"p /sbi/users/interchange/boliva/ModCRE/pbm/CisBP_2019/Escores.txt | cut -d'" + "\t"+ "' -f " + str(altindexlist[0]+1)
				try:
					#altscorelist.append(float(subprocess.check_output([callstring], shell = True).decode('utf-8').strip()))
					altscorelist.append(df[altmatch[0]][altindex-1])
				except:
					print(altindex)
					print(altsequence.upper())
					exit(0)
			else:
				altscorelist.append("NA")
                		
		combinedLists = [scorelist, altscorelist]
		linetowrite = str(combinedLists).replace("[", "").replace("]","")
                
		""" Finished with the addition, I will just need to add the line to be written along with the averaged scores, change indicator and class """
                
                
                
                
                
                
		""" Here I am adding a feature saying if there was a change in motif"""
		if(WildMotif == MutantMotif):
			switchinMotif = 0
		else:
			switchinMotif = 1
		
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


		result = open("DataOutput.csv","a")
		line = line.split('\t')
		#addition_2 = line[0] + "\t" + str(int(line[1])-1) + "\t" + str(line[1]) + "\t" + line[0] + ":" + str(line[1]) + ":" + line[2] + ":" + line[3] + "\t" + "0" + "\t" + ".\n"
		#result2 = open("LineOutput.csv","a")
		addition =  str(refwindow).strip("[").strip("]") + ", " + str(altwindow).strip("[").strip("]") + ", " + str(switchinMotif) + ", " + linetowrite + ", " + classification + ", " + WildMotif + ", " + MutantMotif + "\n"
		
		result.write(addition)
		#result2.write(addition_2)
		result.close()
		#result2.close()
	print("Done")
	


