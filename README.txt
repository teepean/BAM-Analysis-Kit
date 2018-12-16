BAM Analysis Kit
----------------
Website: www.y-str.org
Project Page: www.y-str.org/2014/04/bam-analysis-kit.html

FEEDBACK / SUGGESTIONS
----------------------
Developer: Felix Immanuel <i@fi.id.au>
Feel free to contact me for feedback and suggestions.

SOURCE CODE
-----------
All sources are placed in 'src' folder. 
If any code is not found, it may be uploaded to GitHub at http://www.github.com/fiidau

DOCUMENTATION:
-------------
Execute 'BAM Analysis Kit', select the .BAM file and click 'Analysis'. 
After a few minutes to several hrs (or even days depending on your BAM
file and your computer processing power), your files will be ready
inside a sub-folder called 'out'

If you wish to execute from command-line, just execute console_bam.bat
The user interface 'BAM Analysis Kit.exe' simply modifies the bamkit.config
file and calls the console_bam.bat with input file as parameter. The 
batch file console_bam.bat reads the bamkit.config file and processes
accordingly.

CHANGE LOGS:
-----------
Change Log :1.8
* Bug Fix - filtered-x-chromosome-o37-results.csv is empty due to typo fixed.

Change Log :1.7
* Bug Fix - Unable to load BAM from folders with spaces fixed.

Change Log :1.6
* Human Genome Upgraded to GRCh37.75
* SAMtools upgraded to 1.2
* lobSTR upgraded to 3.0.3
* SNPs upgraded to dbSNP 142
* GATK upgraded to 3.4
* Picard upgraded to 1.134
* Filters SNPs tested by DNA testing companies
* Provides mtDNA mutations in both RSRS and rCRS
* Identifies mtDNA and ISOGG Y haplogroup
* Reports Population Admixture using dv3, globe13 and eurogenes36.
* UI/Speed improvements
* Minor code change and bug fixes
* lobSTR used to calculate only Y-STR but CODIS processing removed.
* SNPedia report
* Download size and operational disk space significantly reduced.

Change Log :1.5
* Bugfix - supports Y-STR/CODIS for all BAMs.

Change Log :1.4
* Upgraded lobSTR to v3.0.2, GATK to v3.2.2
* Output includes more accurate Y-STR values.
* Includes CODIS output.
* Separated Y-STR and CODIS as optional.
* Displays used software version for convenience and advanced use.
* Adds Read Group tags for BAM files without them.

Change Log :1.3
* Updated with lobSTR version 2.0.8 beta.

Change Log :1.2
* Positon offset fixed.

Change Log :1.1
* Some Y-DNA files don't get created if Y-Chromosome alone is selected - bug fixed.

Change Log :1.0
* Supports all BAM files with build 37 positions.
* mtDNA RSRS Values
* Y-STR values
* Y-SNPs values
* Autosomal Values
* X-DNA Values
* Telomere Length
* Processing Selected Chromosome