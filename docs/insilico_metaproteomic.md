![insilico](https://img.shields.io/badge/TYPE-in_silico-d3d7cf)
![author](https://img.shields.io/badge/AUTHOR-Donato_Giovannelli-ad7fa8)
![created](https://img.shields.io/badge/created-d28052020-lightgray)
[![giovannellilab](https://img.shields.io/badge/BY-Giovannelli_Lab-blue)](http://dgiovannelli.github.io)

 <img src="https://dgiovannelli.github.io//images/logopic/giovannellilab.png" width="100 px">

# Extended protocol name

**OBJ**: Basic metaproteomic processing pipeline using GalaxyP

>***Duration: variable, depending on dataset size***

The full set of tutorials to learn to work on Galaxy is available on the [GalaxyP website](https://galaxyproject.org/learn/).

#### Files and Equipment
- Access to internet
- The LC-MS/MS raw data, preferably in the .mgf format
- The corresponding assembled metagenome contigs in fasta nucleotide format

## Procedure
1. Access the [proteomics.usegalaxy.eu/](https://proteomics.usegalaxy.eu/) server
2. Log in with your account (or register if this is the first time accessing GalaxyP Europe)
3. Upload the MS data and the corresponding metagenomic contigs. For large files or large number of files use the FTP upload. You will need an FTP client like [FileZilla](https://filezilla-project.org/). The video tutorial is found in the [GalaxyP tutorial page](https://galaxyproject.org/tutorials/upload/)
4. Build a Dataset list for the MGF files you want to annotate together
5. From here you can move to section **Match peptide sequences** or to the section **Format your metagenome contigs** in case your fasta is not already composed of peptides.

### Format your metagenome contigs
A basic GalaxyP tutorial coviring how to get or create a peptide database can be found on the [GalaxyP tutotial website](https://galaxyproject.github.io/training-material/topics/proteomics/tutorials/database-handling/tutorial.html).
1. Upload your fasta nucleotide (if not already done in the previous step) containing your assembled contigs
2. Select the MetaGeneAnnotator (mga) tool from the Metagenomic Analysis menu. MetaGeneAnnotator will predict Open Reading Frames (ORF) from prokaryote and phage **//** Run Glimmer3 to predict the Open Reading Frames (ORF) on the contigs **//** Run *CAT contigs* to o predict the Open Reading Frames (ORF) on the contigs
4. Digest the predicted proteins using the *Digestor* tool to obtain the predicted peptides
5. Use the Protein Database Downloader tool is used to download the FASTA database from UniProt and cRAP (common Repository of Adventitious Proteins) database. To be able to distinguish contaminants from proteins of interest, you should add a tag to each contaminant protein. In order to do it run in order: 1. Run *FASTA-to-Tabular* tool on your crap database; 2. Run *Add Column* tool on the new output. In the field Add this value enter “CONTAMINANT” and execute. 4. Run *Tabular-to-FASTA* tool; Use column 1 and column 3 as Title columns and column 2 as sequence column; 4. Rename the *Tabular-to-FASTA* tool output to “Tagged cRAP database”
6. Depending on the search algorithm in use, you might need to merge all FASTA entries (i.e. proteins of interest and contaminants) in a single database. Run *FASTA Merge Files and Filter Unique Sequences* tool on the main database and the tagged cRAP database
7. Perform protein identification using *SearchGUI* to use various search engines and prepare results for input to Peptide Shake
8. Peptide Shaker Perform protein identification using various search engines based on results from SearchGUI

#### NOTES
The MetaGeneAnnotator is be used by [sixgill](https://github.com/dhmay/sixgill) to create the six frame translation to construct a database of 'metapeptides', short protein fragments for database search of LC-MS/MS metaproteomics data. Sixgill runs on the fastq file directly.

## Expected results (quantitative information, graphics, images)
Describe the expected results and if possible include a visual example (i.e. a picture, table or graph)

## Common problems, troubleshooting and solutions
- List all common problems and solutions

## Calculations
A detailed explanation of any formula used to convert or calculate the final data from the raw input should be included here. Please refer to any additional software needed or _in silico_ protocol. Please be sure to upload a .ods file with the example calculations along with this protocol if necessary. Use the same filename adding _calculation at the end. For example Wetbench_SedimentDNA_calculations.ods. You can render LaTeX mathematical expressions using [KaTeX](https://khan.github.io/KaTeX/) directly in StackEdit like:
$$
BAR=\frac{B-A}{B+A}
$$
A math cheatsheet to write formulas can be found in [KaTex documentation page](https://katex.org/docs/supported.html).

> Written with [[StackEdit](https://stackedit.io/)
