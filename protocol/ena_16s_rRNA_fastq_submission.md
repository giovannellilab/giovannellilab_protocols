![insilico](https://img.shields.io/badge/TYPE-in_silico-d3d7cf)
![author](https://img.shields.io/badge/AUTHOR-Matteo_Selci_and_Angelina_Cordone-ad7fa8)
![created](https://img.shields.io/badge/created-19/05/2021-lightgray)
[![giovannellilab](https://img.shields.io/badge/BY-Giovannelli_Lab-blue)](http://dgiovannelli.github.io)
 
# ENA submission of 16S rRNA fastq

>Protocol objective. 
- Submitting raw fastq files (Fw and Rw) on The European Nucleotide Archive (ENA)  

>**Duration: 1 h**

>**Adapted from:** https://ena-docs.readthedocs.io/en/latest/submit/reads.html

#### Equipment
- Computer
- Operating System: Linux (Ubuntu) 
- Stable WiFi connection

## Procedure
1. Register to ENA: google https://www.ebi.ac.uk/ena/browser/home, click on **SUBMIT** , click on "**SUBMIT TO ENA INTERACTIVELY**", click on "**REGISTRER**"
2. Once you are logged in, click on "**NEW SUBMISSION**"
3. Open a **Terminal Window** on your Linux System and type:
- 3.1 cd path/to/the/directory/with/your/fastq to submit
- 3.2  sudo apt-get install lftp
- 3.3  lftp webin2.ebi.ac.uk -u Webin-**YOUR ID NUMBER**
- 3.4 Enter your ENA Log-in password when promped
- 3.5 ls command to check the content of your drop box
- 3.6 mput **filename** command to upload files
- 3.6.1 To apply the same command on all fast.gz files type:  mput *fastq.gz
- 3.7 ls command to check of the uploaded files
- 3.8 bye command to exit

4. openssl md5 *R1_001.fastq.gz > md5_r1.csv  && openssl md5 *R2_001.fastq.gz > md5_r2.csv command to create two csv files for R1_fastq and R2_fastq, respectively. You will use these file in the step xxx
5.  Close the **Terminal** and come back on the **ENA SUBMISSION PAGE** 
6. Select **Submit sequence reads and experiments** option and click on **Next**
7. Click on **Create a new study**
8.  Compile all the lines as requested and click on **Next**.
9. Click on **Check list** and choose the appropriate format; i.e. for marine water samples choose: **GSC MIxS water** and click on **Next**
10. Click on **non-samples terms** and  select the option **pcr primers**
11. Click on **Download Template Spreadsheet**
### Template example  

#checklist_accession | xxxx | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp;
-- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | --
#unique_name_prefix | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp;
sample_alias | tax_id | scientific_name | common_name | sample_title | sample_description | project name | pcr primers | sequencing method | investigation type | collection date | geographic location (country and/or sea) | geographic location (latitude) | geographic location (longitude) | water environmental package | geographic location (depth) | environment (biome) | environment (feature) | environment (material)
#template | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | water | &nbsp; | &nbsp; | &nbsp; | &nbsp;
#units | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | DD | DD | &nbsp; | m | &nbsp; | &nbsp; | &nbsp;
1-xxxx_S118 | 358574 | 16S_rRNA | ST12 | ST12 | water | Ant17 | V3-V4 (515FB-926R) | Illumina miSeq | bacteria_archaea | 2017-01 | Ross Sea | -75.07186 | 163.70473 | water | 35 | water | water | water

**IMPORTANT:** The tax_id **358574**  is referred to the **uncultured microorganism [species]** indicated by **ENA** . Choose the more representative **tax_id** for your project.

12. Once the template is completed, save it as .tsv and click on **New Submission** and upload it typing on **Submit Completed Spreadsheet**
13. If errors are not present, click on **Next**
14. Click on **Two Fastq files (Paired)** and on **Next**
15. Click on **Download Template Spreadsheet**
### Template example

sample_alias | instrument_model | library_name | library_source | library_selection | library_strategy | design_description | library_construction_protocol | insert_size | forward_file_name | forward_file_md5 | reverse_file_name | reverse_file_md5
-- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | --
1-xxxx_S118 | Illumina MiSeq | &nbsp; | METAGENOMIC | PCR | AMPLICON | &nbsp; | &nbsp; | xxxx | 1-xxxx_S118_L001_R1_001.fastq.gz | calc md5 value | 1-xxxx_S118_L001_R2_001.fastq.gz | calc md5 value

**IMPORTANT** The md5 value is present within the file **md5_r1.csv** and **md5_r2.csv** (see above, point **4**).

16. Once the template is completed, save it as .tsv and click on **Upload Completed Spreadsheet** and then **Submit**.
17.  Click on **Runs** to check that everything is ok (i.e. click on **Show files** to check if both R1 and R2 fastq.gz are present)
