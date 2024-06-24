# SNPeBoT
Predictor of Transcription Factor Allele Specific Binding  

## Instructions for SNPeBoT CommandLine version

### Dependencies
SBILib  
fimo  
dssp  
genopyc 2.1.5  
tensorflow 2.13.0  
matplotlib 3.8.2  
jsonpickle 3.0.2  

After Downloading SNPeBoT the config.ini file must be edited  
open the config.ini  
under "dssp_path = ./dssp-4.4.2006/" and "meme_path = /soft/system/software/MEME/5.1.1-GCCcore-10.2.0-Python-3.8.6/bin/" change the paths to those corresponding to the executables on your machine  

pwms within the Resources folder must be decompressed  
In the Standalone/Resources folder execute the following command:  
```console
tar -xvzf pwms.tar.gz
```
Finally the Escores.txt file must be downloaded from CISBP and placed in the Resources folder  
To download the Escores file follow these steps (from the Resource folder):
```console
wget cisbp.ccbr.utoronto.ca/data/2.00/DataFiles/Bulk_downloads/EntireDataset/Escores.txt.zip
unzip Escores.txt.zip  
```
### Converting vcf or rsids

We have included a python script to enble the running of mutations in the format of rsids or vcf.
First all mutations must be stored in a file, then this file along with the relevant TF and an output file name should be passed to the command CreateInput.py like so:

for rsids:
```
python CreateInput.py -r {path to input rsid file} -t {TF name} -o {output file}
```

for vcf:
```
python CreateInput.py -v {path to input rsid file} -t {TF name} -o {output file}
```

### Running 

In order to run the program exeute the following command (from the Standalone folder) in the terminal  


```console
bash SNPeBoT.sh {path to input SNP list} {path to input TF list} {pval threshold}
```
The SNP input list is a tab seperated file with the following columns  
1) TF name that should be applied to this lines SNP  
2) The reference Sequence (51 bases with the central position being occupied be the mutation)  
3) The alternate Sequence (51 bases with the central position being occupied be the mutation)  

The TF input list is a txt file with one TF per row  


The output will be stored in the Results folder with the file name Output.txt prepended by the TF test  

### Additional PBM Experiments

It is possible to add PBM data not included in the CISBP Database. First the user would have to add E-score values to the Escores file downloaded in the set up step. To do this ensure that the tested DNA binding domain Ids are stored as column headers and that corresponding octamers align with the rownames as ordered in the Resources/8mer.txt file.
Next the Ids that were added as column headers will have to be added to a stored Dictionary (motifs.pickle). This dictionary stores as keys various transcription factors, with their values being a list of Ids from the E-score table relevant to that TF. Add the IDs to the list stored under the target TF key and save the dictionary. If the target TF doesn't exist as a key, create a new key value pair. Finally, store the PWM of the DNA binding domain in the Resources/pwms folder. Ensure that the pwm file is named with the column header id and bears the suffix ".meme". 

