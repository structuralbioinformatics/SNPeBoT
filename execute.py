import argparse, LowerTest

parser = argparse.ArgumentParser(description="execute script")
parser.add_argument('-t','--tf' , default=None, dest='transcription_factor', action='store', help='input tf')
parser.add_argument('-p','--pval' , default=None, dest='p_val', action='store', help='input the fimo filter p-value to be used')
options = parser.parse_args()


TF = options.transcription_factor.strip()
Pval = options.p_val

linefile = "TF/" + TF + "/SNP.tsv"

x = LowerTest.Prediction(TF,linefile,Pval)


x.set_Escore()

OutputFileLocation = "Results/" + TF + "_Output.txt"

wd = open(OutputFileLocation,"w")
wd.write("SNP Line\tprediction\n")
wd.close()
wd = open(OutputFileLocation,"a")
for pred in x.get_PredictionList():
	pred = str(pred[1]) + "\t" + pred[0] + "\n"
	wd.write(pred)

#x.SavePerformance("TestUserOutput.tsv")









