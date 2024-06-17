import pickle, requests
import argparse

parser = argparse.ArgumentParser(description='From an input file with one rsid or vcf mutation per line retrieve the SNPeBoT input format')
parser.add_argument('-r','--rsids' , default=None, dest='rsids_input', action='store', help='a file containing input rsids')
parser.add_argument('-v','--vcf' , default=None, dest='vcf_input', action='store', help='a file containing input vcfs')
parser.add_argument('-t','--tf' , default=None, dest='TF_input', action='store', help='The TF for the input mutations')
parser.add_argument('-o','--out' , default="SNPList.tsv", dest='output', action='store', help='The location where the sequences are to be saved')
options = parser.parse_args()

Final_list = []
if options.rsids_input:

        import genopyc as gp
        
        with open(options.rsids_input) as fd:
                for line in fd:
                        Final_list.append(line.strip())
        print("getting coordinates")
        coordinates = gp.get_variants_position(Final_list)

        for rsid in coordinates:
                try:
                        url = 'https://api.genome.ucsc.edu/getData/sequence'
                        chrom = rsid[1]
                        snppos = rsid[2]
                        start = int(snppos) - 26
                        end = int(snppos) + 25
        
                        params = {"genome":"hg38","chrom":chrom,"start":start,"end":end}
                        try:
                                response = requests.get(url,params=params)
                        except:
                                print("chrom: " + chrom + "start: " + str(start) + " end: " + str(end))
                                exit(0)
                        output = response.json()
        
                        Sequence = output["dna"]
                        ref = gp.get_ancestral_allele(rsid[0])[rsid[0]]
                        alt = gp.get_variants_info(rsid[0])[rsid[0]]['mappings'][0]['allele_string'].split("/")[1].strip()
                        AltSequence = Sequence[:25] + alt.strip("\"") + Sequence[26:]
                        Sequence = Sequence[:25] + ref.strip("\"") + Sequence[26:]
                        wt_SNP = open(options.output,"a")
                
                        addition = options.TF_input + "\t" + Sequence + "\t" +  AltSequence + "\n"
                        print(options.TF_input + "\t" + Sequence + "\t" +  AltSequence)
                        wt_SNP.write(addition)
                except:
                        print("cant get sequence for: " , rsid[0])



if options.vcf_input:
        with open(options.vcf_input) as fd:
                for line in fd:
                        nexts = []
                        nexts.append("vcf")
                        nexts.append(line.strip().split("\t")[0])
                        nexts.append(line.strip().split("\t")[1])
                        nexts.append(line.strip().split("\t")[3])
                        nexts.append(line.strip().split("\t")[4])
                        Final_list.append(nexts)
        print("getting coordinates")
        coordinates = Final_list

        for vcf in coordinates:
                try:
                        url = 'https://api.genome.ucsc.edu/getData/sequence'
                        chrom = vcf[1]
                        snppos = vcf[2]
                        start = int(snppos) - 26
                        end = int(snppos) + 25
        
                        params = {"genome":"hg38","chrom":chrom,"start":start,"end":end}
                        try:
                                response = requests.get(url,params=params)
                        except:
                                print("chrom: " + chrom + "start: " + str(start) + " end: " + str(end))
                                exit(0)
                        output = response.json()
        
                        Sequence = output["dna"]
                        ref = vcf[3]
                        alt = vcf[4]
                        AltSequence = Sequence[:25] + alt.strip("\"") + Sequence[26:]
                        Sequence = Sequence[:25] + ref.strip("\"") + Sequence[26:]
                        wt_SNP = open(options.output,"a")
                
                        addition = options.TF_input + "\t" + Sequence + "\t" +  AltSequence + "\n"
                        print(options.TF_input + "\t" + Sequence + "\t" +  AltSequence)
                        wt_SNP.write(addition)
                except:
                        print("cant get sequence for: " , vcf[0])
        





