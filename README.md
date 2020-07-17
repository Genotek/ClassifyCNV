## **ClassifyCNV**

ClassifyCNV is a tool that implements the 2019 ACMG guidelines to evaluate the pathogenicity of
germline duplications and deletions. 


### **Requirements**

ClassifyCNV runs on UNIX, Linux and MacOS. All of the necessary databases are included in the
repository.

Python 3.6 or above and BEDTools version 2.27.1 or above must be installed.

### **Input**

ClassifyCNV accepts a BED file as input. The file must include the following columns in
this order:
- chromosome
- CNV start position
- CNV end position
- CNV type (DEL or DUP)

CNVs on alternative contigs are not evaluated.
Both hg19 and hg38 coordinates are supported. 

### **Running the program**

Before running ClassifyCNV it is recommended that the ClinGen files are updated by executing
update_clingen.sh

Command to run ClassifyCNV:
```
python3 ClassifyCNV.py --infile YourCNVFile.bed --GenomeBuild {hg19,hg38}
```

Optional parameters:

```--cores```: number of threads to use; default is 1

```--precise```: should be used only if the exact CNV breakpoints are known; if this flag is on, the script will evaluate the effect of intragenic CNVs

```--outdir```: specify the name of a run-specific directory where the results will be saved to; the directory will be created inside the ClassifyCNV_results folder

### **Examples**

Sample datasets and a sample scoresheet are included in the Examples folder. 

### **Results**

Results are saved to the run-specific folder inside the ClassifyCNV_results folder.
The run-specific folder is named Result_dd_Mon_yyyy-hh-mm-ss.
The filled out table with the scores and the final classification 
of each CNV is named Scoresheet.txt. The column names correspond
to the section numbers of the ACMG rubrics available here:

http://cnvcalc.clinicalgenome.org/cnvcalc/cnv-loss

http://cnvcalc.clinicalgenome.org/cnvcalc/cnv-gain

The last two columns of the file include a list of dosage-sensitive genes contained within the CNV and a list of all protein coding genes in the CNV.

### **Citation**

If you use the software, please cite:


### **License**

The software provided herein is free for academic instruction and research use only. Commercial licenses are available to legal entities, including companies and organizations (both for-profit and non-profit), requiring the software for general commercial use. To obtain a commercial license please, contact us via e-mail: info@genotek.ru.

### **Disclaimer**

This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. The code should not be modified and/or redistributed without the permission of the authors.

### **Authors**

Tatiana Gurbich (Genotek Ltd)

Valery Ilinsky (Genotek Ltd)
