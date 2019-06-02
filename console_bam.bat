@echo off
title BAM Analysis Kit 2.09

REM     The MIT License (MIT)
REM     Copyright © 2013-2015 Felix Immanuel
REM     http://www.y-str.org
REM     
REM     Permission is hereby granted, free of charge, to any person obtaining a copy
REM     of this software and associated documentation files (the Softwareù), to deal
REM     in the Software without restriction, including without limitation the rights
REM     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
REM     copies of the Software, and to permit persons to whom the Software is furnished
REM     to do so, subject to the following conditions: The above copyright notice and
REM     this permission notice shall be included in all copies or substantial portions
REM     of the Software. THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY
REM     KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
REM     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
REM     EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES
REM     OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
REM     ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
REM     OTHER DEALINGS IN THE SOFTWARE.

echo.
echo *** BAM Analysis Kit 2.09 ***
echo.
echo Project Website: http://www.y-str.org/2014/04/bam-analysis-kit.html
echo Tools Used: SAMTools, picard, bamtools, lobSTR, telseq, Cygwin, GATK, Java, Yleaf v2, haplogrep
echo Script Developer: Felix Immanuel ^<i@fi.id.au^>
echo.

if [%1]==[] goto NOPARAM

REM - start reporting versions..
bin\cygwin\bin\bash.exe -c "echo -n 'lobSTR Version: ';/bin/lobSTR --version"
bin\cygwin\bin\bash.exe -c "/bin/bedtools --version"
bin\bow\samtools.exe --version
bin\cygwin\bin\bash.exe -c "echo -n 'Cygwin Version: ';/bin/uname -r"
bin\cygwin\bin\bash.exe -c "/bin/telseq --version|/bin/head -1"
bin\jre\bin\java -version 2^>^&1 ^| findstr /i version
bin\cygwin\bin\bash.exe -c "echo -n 'Picard Version: ';./bin/jre/bin/java.exe -jar bin/picard/picard.jar AddOrReplaceReadGroups --version"
REM - end reporting versions..

REM - saving old processing if it wasnt saved and accidentally started new processing
IF EXIST out.old\genome_complete.txt (
 bin\cygwin\bin\bash.exe -c "/bin/rm -fr out.old"
)

IF EXIST out\genome_complete.txt (
 bin\cygwin\bin\bash.exe -c "/bin/mv out out.old"
 mkdir out
)


SETLOCAL ENABLEDELAYEDEXPANSION
SET "CHR_1="
SET "CHR_2="
SET "CHR_3="
SET "CHR_4="
SET "CHR_5="
SET "CHR_6="
SET "CHR_7="
SET "CHR_8="
SET "CHR_9="
SET "CHR_10="
SET "CHR_11="
SET "CHR_12="
SET "CHR_13="
SET "CHR_14="
SET "CHR_15="
SET "CHR_16="
SET "CHR_17="
SET "CHR_18="
SET "CHR_19="
SET "CHR_20="
SET "CHR_21="
SET "CHR_22="
SET "CHR_X="
SET "CHR_Y="
SET "CHR_M="
SET "YSTR="
SET "TELOMERE="
SET "BAMKIT_JVM="
SET "BAMKIT_THREADS="
SET "ADMIXTURE="
SET "DEL_VCF="
SET "SNPEDIA="

for /f "tokens=1,2 delims== " %%X in (bamkit.config) do (
	IF NOT "%%X" == "#" (

	  IF NOT "%%Y" == "" (
	  set "%%X=%%Y"
	  )		  
	)
)

echo.
echo Configuration from bamkit.config
echo.
echo Process chr 1 set to %CHR_1%
echo Process chr 2 set to %CHR_2%
echo Process chr 3 set to %CHR_3%
echo Process chr 4 set to %CHR_4%
echo Process chr 5 set to %CHR_5%
echo Process chr 6 set to %CHR_6%
echo Process chr 7 set to %CHR_7%
echo Process chr 8 set to %CHR_8%
echo Process chr 9 set to %CHR_9%
echo Process chr 10 set to %CHR_10%
echo Process chr 11 set to %CHR_11%
echo Process chr 12 set to %CHR_12%
echo Process chr 13 set to %CHR_13%
echo Process chr 14 set to %CHR_14%
echo Process chr 15 set to %CHR_15%
echo Process chr 16 set to %CHR_16%
echo Process chr 17 set to %CHR_17%
echo Process chr 18 set to %CHR_18%
echo Process chr 19 set to %CHR_19%
echo Process chr 20 set to %CHR_20%
echo Process chr 21 set to %CHR_21%
echo Process chr 22 set to %CHR_22%
echo Process chr X set to %CHR_X%
echo Process chr Y set to %CHR_Y%
echo Process chr M set to %CHR_M%
echo Process Y-STR set to %YSTR%
echo Process telomere set to %TELOMERE%
echo JVM value is %BAMKIT_JVM%
echo No of parallel threads is %BAMKIT_THREADS%
echo Population admixture set to %ADMIXTURE%
echo Delete VCF after processing set to %DEL_VCF%
echo SNPedia Report set to %SNPEDIA%

echo.
echo.
echo Input BAM : %1
echo.

echo Pre-execution Cleanup ...

IF EXIST chrY_1.tab DEL /Q /F chrY_1.tab 
IF EXIST chrY_1.tmp DEL /Q /F chrY_1.tmp 
IF EXIST chrM_2.tmp DEL /Q /F chrM_2.tmp 
IF EXIST mtdna_max_pos DEL /Q /F mtdna_max_pos 
IF EXIST dbsnp_chrM.tab DEL /Q /F dbsnp_chrM.tab
IF EXIST snps.sorted DEL /Q /F snps.sorted
IF EXIST ref.sorted DEL /Q /F ref.sorted
IF EXIST snps.tmp DEL /Q /F  snps.tmp
IF EXIST chrY.tmp DEL /Q /F  chrY.tmp
IF EXIST lobSTR_CODIS.out DEL /Q /F  lobSTR_CODIS.out
IF EXIST lobSTR_Y-STR.out DEL /Q /F  lobSTR_Y-STR.out
IF EXIST inchr.bam DEL /Q /F  inchr.bam
IF EXIST bam_wh_tmp.bam DEL /Q /F  bam_wh_tmp.bam
IF EXIST ref.fa DEL /Q /F  ref.fa
IF EXIST ref.fa.fai DEL /Q /F  ref.fa.fai
IF EXIST ref.dict DEL /Q /F  ref.dict
IF EXIST bam_sorted.bam DEL /Q /F bam_sorted.bam
IF EXIST bam_sorted.bam.bai DEL /Q /F bam_sorted.bam.bai
IF EXIST bam_sorted_realigned.bam DEL /Q /F bam_sorted_realigned.bam
IF EXIST header DEL /Q /F  header
IF EXIST header01 DEL /Q /F  header01
IF EXIST header02 DEL /Q /F  header02
IF EXIST inchr.sam DEL /Q /F  inchr.sam
IF EXIST tmp.sam DEL /Q /F  tmp.sam
IF EXIST chr%%A.bam DEL /Q /F  chr%%A.bam
IF EXIST reads.bam DEL /Q /F  reads.bam
IF EXIST bam.intervals DEL /Q /F  bam.intervals
IF EXIST bam_sorted_realigned.bam.bai DEL /Q /F  bam_sorted_realigned.bam.bai
IF EXIST bam_sorted_realigned.bai DEL /Q /F  bam_sorted_realigned.bai
IF EXIST bam_wh.bam DEL /Q /F  bam_wh.bam
IF EXIST bam_out.vcf DEL /Q /F  bam_out.vcf
IF EXIST bam_out.vcf.idx DEL /Q /F  bam_out.vcf.idx
IF EXIST chr DEL /Q /F  chr
IF EXIST bam_chr%%A.vcf  DEL /Q /F  bam_chr%%A.vcf 
IF EXIST bam_chr%%A.vcf.gz  DEL /Q /F  bam_chr%%A.vcf.gz
IF EXIST bam_chr%%A.vcf.gz.tbi DEL /Q /F  bam_chr%%A.vcf.gz.tbi
IF EXIST bam_chr%%A.vcf.idx DEL /Q /F  bam_chr%%A.vcf.idx
IF EXIST chr%%A.tab DEL /Q /F  chr%%A.tab
IF EXIST chrM.seq DEL /Q /F  chrM.seq
IF EXIST chrM.tmp.tab DEL /Q /F  chrM.tmp.tab
IF EXIST chrM.tmp DEL /Q /F  chrM.tmp
IF EXIST genotype.txt DEL /Q /F genotype.txt
IF EXIST snps.txt DEL /Q /F snps.txt
IF EXIST bam_complete_sorted.bam DEL /Q /F bam_complete_sorted.bam
IF EXIST bam_complete_sorted.bam.bai DEL /Q /F bam_complete_sorted.bam.bai
IF EXIST ystr.filters DEL /Q /F ystr.filters
IF EXIST bed.a DEL /Q /F bed.a
IF EXIST bam_strs.aligned.stats DEL /Q /F bam_strs.aligned.stats
IF EXIST bam_strs_sorted.bam DEL /Q /F bam_strs_sorted.bam
IF EXIST bam_strs_sorted.bam.bai DEL /Q /F bam_strs_sorted.bam.bai
IF EXIST bam_ystrs.allelotype.stats DEL /Q /F bam_ystrs.allelotype.stats
IF EXIST bam_ystrs.vcf DEL /Q /F bam_ystrs.vcf
IF EXIST bam_strs.aligned.bam DEL /Q /F bam_strs.aligned.bam
IF EXIST bam_out_variants.vcf DEL /Q /F bam_out_variants.vcf
echo.

echo Sorting ...
REM - using windows version to get rid of escape seq for unix
bin\bow\samtools.exe sort -@ %BAMKIT_THREADS% %1 -o bam_complete_sorted.bam

echo.
echo Indexing the sorted BAM file ...
bin\bow\samtools.exe index -@ %BAMKIT_THREADS% bam_complete_sorted.bam

bin\cygwin\bin\bash.exe -c "/bin/echo -e '# rsid\tchr\tpos\tgenotype' > out/genome_full_snps.txt"
bin\cygwin\bin\bash.exe -c "/bin/echo -e '# chr\tpos\tgenotype\trsid' > out/genome_complete.txt"
bin\cygwin\bin\bash.exe -c "/bin/echo -e 'RSID,CHROMOSOME,POSITION,RESULT' > out/filtered-autosomal-o37-results.csv"
bin\cygwin\bin\bash.exe -c "/bin/echo -e 'RSID,CHROMOSOME,POSITION,RESULT' > out/filtered-x-chromosome-o37-results.csv"

REM --- Execute genome tools
FOR %%A IN (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M) DO (

IF "!CHR_%%A!" == "yes" (

echo.
echo ******* Processing Chromosome %%A *******
echo.

echo Splitting and preparing Chr %%A ...
REM it can be chr22 or 22 -- must be a way to detect.
bin\cygwin\bin\bash.exe -c "bin/bow/samtools.exe view -@ %BAMKIT_THREADS% -H bam_complete_sorted.bam|/bin/cut -f2|/bin/grep SN|/bin/grep %%A|/bin/cut -d':' -f2|/bin/head -1 > chr"
for /F "tokens=1" %%C in (chr) do (
bin\cygwin\bin\bash.exe -c "bin/bow/samtools.exe view -@ %BAMKIT_THREADS% -b bam_complete_sorted.bam %%C > chr%%A.bam"
)

copy ref\chr%%A.fa.gz ref.fa.gz > NUL

bin\cygwin\bin\bash.exe -c "/bin/gzip -d ref.fa.gz"

echo Checking and fixing ...
IF EXIST inchr.bam DEL /Q /F  inchr.bam
copy chr%%A.bam inchr.bam > NUL

echo.
echo Indexing Human Reference Genome ...
bin\bow\samtools.exe faidx ref.fa
echo.
echo Creating Sequence Dictionary for Human Reference Genome ...
bin\cygwin\bin\bash.exe -c "./bin/jre/bin/java.exe %BAMKIT_JVM% -jar bin/picard/picard.jar CreateSequenceDictionary R=ref.fa O=ref.dict"

echo.
echo Preparing BAM ...
bin\cygwin\bin\bash.exe -c "bin/bow/samtools.exe view -@ %BAMKIT_THREADS% inchr.bam | /bin/sed 's/\t/\tchr/2' > tmp.sam"
bin\cygwin\bin\bash.exe -c "/bin/cat tmp.sam | /bin/sed 's/\tchrchr/\tchr/' > inchr.sam"

if "%%A" == "M" (
  bin\cygwin\bin\bash.exe -c "/bin/cat inchr.sam | /bin/sed 's/\tchrMT/\tchrM/' > inchr_tmp.sam"
  DEL /Q /F inchr.sam
  copy inchr_tmp.sam inchr.sam > NUL
  DEL /Q /F inchr_tmp.sam
)

bin\cygwin\bin\bash.exe -c "bin/bow/samtools.exe view -@ %BAMKIT_THREADS% -bT ref.fa inchr.sam > reads.bam"
bin\cygwin\bin\bash.exe -c "bin/bow/samtools.exe view -@ %BAMKIT_THREADS% -H inchr.bam|/bin/grep -v SN > header01"
bin\cygwin\bin\bash.exe -c "bin/bow/samtools.exe view -@ %BAMKIT_THREADS% -H reads.bam > header02"
copy header01+header02 header /Y /B > NUL
bin\cygwin\bin\bash.exe -c "bin/bow/samtools.exe reheader header reads.bam > bam_wh_tmp.bam"

echo Adding or Replace Read Group Header ...
bin\jre\bin\java.exe %BAMKIT_JVM% -jar bin\picard\picard.jar AddOrReplaceReadGroups INPUT=bam_wh_tmp.bam OUTPUT=bam_wh.bam SORT_ORDER=coordinate RGID=rgid RGLB=rglib RGPL=illumina RGPU=rgpu RGSM=sample VALIDATION_STRINGENCY=SILENT

echo Sorting ...
bin\bow\samtools.exe sort -@ %BAMKIT_THREADS% bam_wh.bam -o bam_sorted.bam

echo.
echo Indexing the sorted BAM file ...
bin\bow\samtools.exe index -@ %BAMKIT_THREADS% bam_sorted.bam


	REM --------- YSTR  
  IF "%%A" == "Y" (
   IF "%YSTR%" == "yes" (
	echo Extracting Y STR Values ...
	echo lobSTR alignment ...
	bin\cygwin\bin\bash.exe -c "/bin/lobSTR --index-prefix ref/lobSTR/lobSTR_  -f bam_sorted.bam --rg-sample ggtsample --rg-lib ggtlibrary --fft-window-size 24 --fft-window-step 12 -o bam_strs -v --bam --noweb -p %BAMKIT_THREADS%"

	echo Sorting and Indexing ...
	bin\bow\samtools.exe sort -@ %BAMKIT_THREADS% bam_strs.aligned.bam -o bam_strs_sorted.bam
	bin\bow\samtools.exe index -@ %BAMKIT_THREADS% bam_strs_sorted.bam

	echo lobSTR allelotyper vcf...
	bin\cygwin\bin\bash.exe -c "export PATH=$PATH:/usr/local/bin:/usr/bin:/usr/lib/lapack; /bin/allelotype --command classify --index-prefix ref/lobSTR/lobSTR_ --out bam_ystrs --bam bam_strs_sorted.bam --noise_model /usr/local/share/lobSTR/models/illumina_v3.pcrfree  --min-border 5 --min-bp-before-indel 7 --maximal-end-match 15 --min-read-end-match 5 --strinfo ref/strinfo.tab -v --noweb"

	echo Generating Y-STR values ...
	bin\cygwin\bin\bash.exe -c "/bin/bedtools intersect -a bam_ystrs.vcf -b ref/lobSTR_ystr_hg19.bed -wa -wb | /bin/cut -f 1,2,10,14- | /bin/sed 's/:/\t/g' | /bin/cut -f 1,2,4,7,11- | /bin/sed 's/\//\t/' | /bin/cut -f 4 --complement | /bin/awk '{print $0 """\\t""" $6+$4/$5}' > lobSTR_Y-STR.out"
	bin\cygwin\bin\bash.exe -c "/bin/cat lobSTR_Y-STR.out|/bin/cut -f7,8|/bin/sed 's/\t/ = /g' >> out/Y-STR_Markers.txt"

	IF EXIST bam_strs.aligned.stats DEL /Q /F bam_strs.aligned.stats
	IF EXIST bam_strs_sorted.bam DEL /Q /F bam_strs_sorted.bam
	IF EXIST bam_strs_sorted.bam.bai DEL /Q /F bam_strs_sorted.bam.bai
	IF EXIST bam_ystrs.allelotype.stats DEL /Q /F bam_ystrs.allelotype.stats
	IF EXIST bam_ystrs.vcf DEL /Q /F bam_ystrs.vcf
	IF EXIST bam_strs.aligned.bam DEL /Q /F bam_strs.aligned.bam	
   ) 
  )
	REM --------- YSTR


echo.
echo Realignment of the sorted and indexed BAM file ...
  IF "%%A" == "M" (
	bin\cygwin\bin\bash.exe -c "./bin/jre/bin/java.exe %BAMKIT_JVM% -jar bin/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator   -window 3  -minReads 1 -R ref.fa -I bam_sorted.bam  -o bam.intervals"
	bin\cygwin\bin\bash.exe -c "./bin/jre/bin/java.exe %BAMKIT_JVM% -jar bin/gatk/GenomeAnalysisTK.jar -T IndelRealigner  -maxPosMove 10 -R ref.fa -I bam_sorted.bam -targetIntervals bam.intervals -o bam_sorted_realigned.bam"
  ) ELSE (
	bin\cygwin\bin\bash.exe -c "./bin/jre/bin/java.exe %BAMKIT_JVM% -jar bin/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ref.fa -I bam_sorted.bam  -o bam.intervals"
	bin\cygwin\bin\bash.exe -c "./bin/jre/bin/java.exe %BAMKIT_JVM% -jar bin/gatk/GenomeAnalysisTK.jar -T IndelRealigner -R ref.fa -I bam_sorted.bam -targetIntervals bam.intervals -o bam_sorted_realigned.bam"
)

echo.
echo Indexing the realigned BAM file ...
bin\bow\samtools.exe index -@ %BAMKIT_THREADS% bam_sorted_realigned.bam

echo.
echo Invoke the variant caller ...

 set "HAPLOID="
  IF "%%A" == "M" set HAPLOID=1  
  IF "%%A" == "Y" set HAPLOID=1  
  
  IF defined HAPLOID (
  
  IF "%%A" == "M" (
  
	bin\cygwin\bin\bash.exe -c "./bin/jre/bin/java.exe %BAMKIT_JVM% -jar bin/gatk/GenomeAnalysisTK.jar -l INFO -R ref.fa -T HaplotypeCaller -I bam_sorted_realigned.bam -nct %BAMKIT_THREADS% -o bam_chr%%A.vcf"
	bin\bow\bgzip.exe -@ %BAMKIT_THREADS% bam_chr%%A.vcf
	bin\bow\tabix.exe bam_chr%%A.vcf.gz	

	echo Extracting mtDNA markers ..
	bin\cygwin\bin\bash.exe -c "bin/bow/bcftools.exe query -f '%%CHROM\t%%POS\t[%%IUPACGT]\n' bam_chr%%A.vcf |/bin/sed 's/chr//g' | /bin/sed 's/\///g' > chr%%A.tab"	
	bin\cygwin\bin\bash.exe -c "bin/bow/bcftools.exe query -f '%%CHROM\t%%POS\t[%%IUPACGT]\n' bam_chr%%A.vcf |/bin/sed 's/chr//g' | /bin/sed 's/\///g'|/bin/cut -f2,3|/bin/awk '{if($2 ~ /^[ATGC]$/) printf $1\"\"$2\" \"}'|/bin/sed 's/\s$//g'|/bin/sed 's/\s/, /g' > out/rCRS_mtDNA.txt"

	REM - realign with RSRS to get accurate mtDNA markers in RSRS

	DEL /Q /F bam_chr%%A.*
	IF EXIST ref.fa DEL /Q /F ref.fa
	IF EXIST ref.dict DEL /Q /F ref.dict
	IF EXIST ref.fa.fai DEL /Q /F ref.fa.fai
	COPY /Y ref\chrM.RSRS.fasta ref.fa >NUL
	bin\bow\samtools.exe faidx ref.fa
	bin\cygwin\bin\bash.exe -c "./bin/jre/bin/java.exe %BAMKIT_JVM% -jar bin/picard/picard.jar CreateSequenceDictionary R=ref.fa O=ref.dict"
	IF EXIST bam.intervals DEL /Q /F bam.intervals
	bin\cygwin\bin\bash.exe -c "./bin/jre/bin/java.exe %BAMKIT_JVM% -jar bin/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator   -window 3  -minReads 1 -R ref.fa -I bam_sorted.bam  -o bam.intervals"
	IF EXIST bam_sorted_realigned.bam DEL /Q /F bam_sorted_realigned.bam
	bin\cygwin\bin\bash.exe -c "./bin/jre/bin/java.exe %BAMKIT_JVM% -jar bin/gatk/GenomeAnalysisTK.jar -T IndelRealigner  -maxPosMove 10 -R ref.fa -I bam_sorted.bam -targetIntervals bam.intervals -o bam_sorted_realigned.bam"
	bin\cygwin\bin\bash.exe -c "./bin/jre/bin/java.exe %BAMKIT_JVM% -jar bin/gatk/GenomeAnalysisTK.jar -l INFO -R ref.fa -T HaplotypeCaller -I bam_sorted_realigned.bam -nct %BAMKIT_THREADS% -o bam_chr%%A.vcf"
	bin\jre\bin\java.exe -jar bin\haplogrep\haplogrep-2.1.20.jar --in bam_chr%%A.vcf --format vcf --out out\mtDNA-haplogroup-haplogrep.txt
	bin\bow\bgzip.exe -@ %BAMKIT_THREADS% bam_chr%%A.vcf
	bin\bow\tabix.exe bam_chr%%A.vcf.gz	
	bin\cygwin\bin\bash.exe -c "bin/bow//bcftools.exe query -f '%%CHROM\t%%POS\t[%%IUPACGT]\n' bam_chr%%A.vcf |/bin/sed 's/chr//g' | /bin/sed 's/\///g' > chr%%A.tab"	
	bin\cygwin\bin\bash.exe -c "bin/bow/bcftools.exe query -f '%%REF\t%%POS\t[%%IUPACGT]\n' bam_chr%%A.vcf |/bin/sed 's/chr//g' | /bin/sed 's/\///g'|/bin/awk '{if($3 ~ /^[ATGC]$/) printf $1\"\"$2\"\"$3\" \"}'|/bin/sed 's/\s$//g'|/bin/sed 's/\s/, /g' > out/RSRS_mtDNA.txt"

  )
  
  IF "%%A" == "Y" (
	bin\cygwin\bin\bash.exe -c "./bin/jre/bin/java.exe %BAMKIT_JVM% -jar bin/gatk/GenomeAnalysisTK.jar -l INFO -R ref.fa -T UnifiedGenotyper -glm SNP -I bam_sorted_realigned.bam -rf BadCigar -nct %BAMKIT_THREADS% -o bam_chr%%A.vcf --output_mode EMIT_ALL_CONFIDENT_SITES"
	bin\bow\bgzip.exe -@ %BAMKIT_THREADS% bam_chr%%A.vcf
	bin\bow\tabix.exe bam_chr%%A.vcf.gz
	
	setlocal
	set PATH=bin\bow;bin\python373
	bin\python373\python bin\yleaf\Yleafw.py -bam bam_sorted_realigned.bam -ref hg19 -out out -q 10 -b 40 -t 8 -r 1
	move out\out.hg out\y-haplogroup.txt
	move out\bam_sorted_realigned\bam_sorted_realigned.chr out\y-haplogroup.chr
	move out\bam_sorted_realigned\bam_sorted_realigned.fmf out\y-haplogroup.fmf
	move out\bam_sorted_realigned\bam_sorted_realigned.log out\y-haplogroup.log
	move out\bam_sorted_realigned\bam_sorted_realigned.out out\y-haplogroup.out
	rmdir /S /Q out\bam_sorted_realigned
	endlocal

	echo Extracting Y-SNP markers ..
	bin\cygwin\bin\bash.exe -c "bin/bow/bcftools.exe query -f '%%CHROM\t%%POS\t[%%IUPACGT]\n' bam_chr%%A.vcf |/bin/sed 's/chr//g' > chr%%A.tab"
	
	bin\cygwin\bin\bash.exe -c "bin/bow/bcftools.exe query -f '%%POS\t%%REF\t[%%IUPACGT]\n' bam_chrY.vcf |/bin/sed 's/chr//g' | /bin/sort -k 1 > chrY.tmp"
	bin\bamkit\extract_ysnps.exe ref\ysnp_hg19.ref chrY.tmp
	bin\cygwin\bin\bash.exe -c "/bin/echo -e 'Position\tReference\tGenotype' > out/Variants_Y.txt"	
	
	
	
	bin\cygwin\bin\bash.exe -c "bin/bow/bcftools.exe query -f '%%POS\t%%REF\t[%%IUPACGT]\n' bam_chr%%A.vcf |/bin/sed 's/chr//g' > chrY_1.tab"
	
	bin\cygwin\bin\bash.exe -c "/bin/awk 'NR==FNR{a[$2]=$3;next} {if(a[$1] == \"\") print $1\"\t\"$2\"\t\"$3;}' <(/bin/gzip -dc ref/dbsnp_chrY.tab.gz) chrY_1.tab|/bin/grep -P -v 'A\tA'|/bin/grep -P -v 'T\tT'|/bin/grep -P -v 'G\tG'|/bin/grep -P -v 'C\tC' >> out/Variants_Y.txt"
	
	bin\cygwin\bin\bash.exe -c "/bin/cat out/Y_SNPs.txt | /bin/sed 's/\s//g'|/bin/sed 's/,/\n/g'|/bin/grep '+'|/bin/sed 's/\+/ = /g' > ystr.filters"
	bin\cygwin\bin\bash.exe -c "/bin/echo -e '# Mutation\tISOGG-Y-Haplogroup' > out/ISOGG_Y_Haplogroup.txt"
	bin\cygwin\bin\bash.exe -c "/bin/grep -f ystr.filters ref/isogg_ytree.ref >> out/ISOGG_Y_Haplogroup.txt"
  ) 
  
  ) ELSE (
	bin\cygwin\bin\bash.exe -c "./bin/jre/bin/java.exe %BAMKIT_JVM% -jar bin/gatk/GenomeAnalysisTK.jar -l INFO -R ref.fa -T UnifiedGenotyper -glm SNP -I bam_sorted_realigned.bam -rf BadCigar -nct %BAMKIT_THREADS% -o bam_chr%%A.vcf --output_mode EMIT_ALL_CONFIDENT_SITES"
	bin\bow\bgzip.exe -@ %BAMKIT_THREADS% bam_chr%%A.vcf
	bin\bow\tabix.exe bam_chr%%A.vcf.gz
    bin\cygwin\bin\bash.exe -c "bin/bow/bcftools.exe query -f '%%CHROM\t%%POS\t[%%TGT]\n' bam_chr%%A.vcf |/bin/sed 's/chr//g' > chr%%A.tab"
  )

bin\cygwin\bin\bash.exe -c "/bin/awk 'NR==FNR{a[$1,$2]=$3;next} ($1,$2) in a{ print a[$1,$2],$1,$2,$3}' <(/bin/gzip -dc ref/dbsnp_chr%%A.tab.gz) chr%%A.tab|/bin/sed 's/\\s/\t/g' > snps.tmp"
bin\cygwin\bin\bash.exe -c "/bin/cat snps.tmp|/bin/sed 's/\///g' >> out/genome_full_snps.txt"
bin\cygwin\bin\bash.exe -c "/bin/awk 'NR==FNR{a[$1,$2]=$3;next} { print $1,$2,$3,a[$1,$2]}' <(/bin/gzip -dc ref/dbsnp_chr%%A.tab.gz) chr%%A.tab|/bin/sed 's/\\s/\t/g' >> out/genome_complete.txt"

bin\cygwin\bin\bash.exe -c "/bin/cat snps.tmp | /bin/sort -k 1 > snps.sorted"

IF "%%A" == "X" (
bin\cygwin\bin\bash.exe -c "/bin/cat ref/snps_filtered/chr%%A | /bin/sort -k 1 > ref.sorted"
bin\cygwin\bin\bash.exe -c "/bin/join -t $'\t' snps.sorted ref.sorted |/bin/sed 's/\///g'|/bin/sort -n -k 3|/bin/awk '{print \"\\\"\"$1\"\\\"\"\",\"\"\\\"\"$2\"\\\"\"\",\"\"\\\"\"$3\"\\\"\"\",\"\"\\\"\"$4\"\\\"\"}' >> out/filtered-x-chromosome-o37-results.csv"
) ELSE (
IF NOT "%%A" == "Y" (
IF NOT "%%A" == "M" (
bin\cygwin\bin\bash.exe -c "/bin/cat ref/snps_filtered/chr%%A | /bin/sort -k 1 > ref.sorted"
bin\cygwin\bin\bash.exe -c "/bin/join -t $'\t' snps.sorted ref.sorted |/bin/sed 's/\///g'|/bin/sort -n -k 3|/bin/awk '{print \"\\\"\"$1\"\\\"\"\",\"\"\\\"\"$2\"\\\"\"\",\"\"\\\"\"$3\"\\\"\"\",\"\"\\\"\"$4\"\\\"\"}' >> out/filtered-autosomal-o37-results.csv"
)
)
)


REM -- final cleanup

IF EXIST chrY_1.tab DEL /Q /F chrY_1.tab 
IF EXIST chrY_1.tmp DEL /Q /F chrY_1.tmp 
IF EXIST chrM_2.tmp DEL /Q /F chrM_2.tmp 
IF EXIST mtdna_max_pos DEL /Q /F mtdna_max_pos 
IF EXIST dbsnp_chrM.tab DEL /Q /F dbsnp_chrM.tab
IF EXIST snps.sorted DEL /Q /F snps.sorted
IF EXIST ref.sorted DEL /Q /F ref.sorted
IF EXIST snps.tmp DEL /Q /F  snps.tmp
IF EXIST chrY.tmp DEL /Q /F  chrY.tmp
IF EXIST lobSTR_CODIS.out DEL /Q /F  lobSTR_CODIS.out
IF EXIST lobSTR_Y-STR.out DEL /Q /F  lobSTR_Y-STR.out
IF EXIST inchr.bam DEL /Q /F  inchr.bam
IF EXIST bam_wh_tmp.bam DEL /Q /F  bam_wh_tmp.bam
IF EXIST ref.fa DEL /Q /F  ref.fa
IF EXIST ref.fa.fai DEL /Q /F  ref.fa.fai
IF EXIST ref.dict DEL /Q /F  ref.dict
IF EXIST bam_sorted.bam DEL /Q /F bam_sorted.bam
IF EXIST bam_sorted.bam.bai DEL /Q /F bam_sorted.bam.bai
IF EXIST bam_sorted_realigned.bam DEL /Q /F bam_sorted_realigned.bam
IF EXIST header DEL /Q /F  header
IF EXIST header01 DEL /Q /F  header01
IF EXIST header02 DEL /Q /F  header02
IF EXIST inchr.sam DEL /Q /F  inchr.sam
IF EXIST tmp.sam DEL /Q /F  tmp.sam
IF EXIST chr%%A.bam DEL /Q /F  chr%%A.bam
IF EXIST reads.bam DEL /Q /F  reads.bam
IF EXIST bam.intervals DEL /Q /F  bam.intervals
IF EXIST bam_sorted_realigned.bam.bai DEL /Q /F  bam_sorted_realigned.bam.bai
IF EXIST bam_sorted_realigned.bai DEL /Q /F  bam_sorted_realigned.bai
IF EXIST bam_wh.bam DEL /Q /F  bam_wh.bam
IF EXIST bam_out.vcf DEL /Q /F  bam_out.vcf
IF EXIST bam_out.vcf.idx DEL /Q /F  bam_out.vcf.idx
IF EXIST chr DEL /Q /F  chr
IF EXIST bam_chr%%A.vcf  DEL /Q /F  bam_chr%%A.vcf 
IF "%DEL_VCF%" == "yes" (
	IF EXIST bam_chr%%A.vcf.gz  DEL /Q /F  bam_chr%%A.vcf.gz
) ELSE (
	IF EXIST bam_chr%%A.vcf.gz  MOVE bam_chr%%A.vcf.gz out >NUL
)
IF EXIST bam_chr%%A.vcf.gz.tbi DEL /Q /F  bam_chr%%A.vcf.gz.tbi
IF EXIST bam_chr%%A.vcf.idx DEL /Q /F  bam_chr%%A.vcf.idx
IF EXIST chr%%A.tab DEL /Q /F  chr%%A.tab
IF EXIST chrM.seq DEL /Q /F  chrM.seq
IF EXIST chrM.tmp.tab DEL /Q /F  chrM.tmp.tab
IF EXIST chrM.tmp DEL /Q /F  chrM.tmp
)
)

IF "%TELOMERE%" == "yes" (
	echo.
	echo Processing Telomere ...
	echo.	
	bin\cygwin\bin\bash.exe -c "/bin/telseq -H -o telseq.out bam_complete_sorted.bam"
	echo Telomere Length: > out\telomere.txt
	bin\cygwin\bin\bash.exe -c "/bin/cat telseq.out|/bin/cut -f7 >> out/telomere.txt"
	bin\cygwin\bin\bash.exe -c "/bin/conv --u2d out/telomere.txt"
	move telseq.out out > NUL
	type out\telomere.txt
)



REM SNPedia Report
IF "%SNPEDIA%" == "yes" (
echo.
echo Processing SNPedia Report ...
bin\bamkit\snpedia_report.exe ref\snpedia.txt.gz out\genome_full_snps.txt out\SNPedia_Report.html
echo.
)


IF "%ADMIXTURE%" == "yes" (

	echo.
	echo Processing Admixtures ...
	echo.

bin\cygwin\bin\bash.exe -c "/bin/cat out/genome_full_snps.txt |/bin/grep -v '#'|/bin/sed 's/\///g'| /bin/sort -k 1 > snps.txt"
bin\cygwin\bin\bash.exe -c "/bin/join -t $'\t' ref/admixture/snps.alleles snps.txt > genotype.txt"
echo.
echo Calculating Population Admixture - dv3 [K=12]
echo.
copy /Y ref\admixture\dv3.* . > NUL
bin\diydodecad\DIYDodecadWin.exe dv3.par
move genomewide.txt out\admix_result_dv3.txt > NUL
DEL /Q /F dv3.* > NUL

echo.
echo Calculating Population Admixture - eurogenes [K=36]
echo.
copy /Y ref\admixture\eurogenes.* . > NUL
bin\diydodecad\DIYDodecadWin.exe eurogenes.par
move genomewide.txt out\admix_result_eurogenes36.txt > NUL
DEL /Q /F eurogenes.* > NUL

echo.
echo Calculating Population Admixture - globe13 [K=13]
echo.
copy /Y ref\admixture\globe13.* . > NUL
bin\diydodecad\DIYDodecadWin.exe globe13.par
move genomewide.txt out\admix_result_globe13.txt > NUL
DEL /Q /F globe13.* > NUL


bin\bamkit\admixure_report.exe out\admix_result_dv3.txt dv3 "Dodecad [dv3]" out\Admixture_Dodecad_Report.html
bin\bamkit\admixure_report.exe out\admix_result_globe13.txt dv3 "Globe13" out\Admixture_Globe13_Report.html
bin\bamkit\admixure_report.exe out\admix_result_eurogenes36.txt eurogenes "Eurogenes [K=36]" out\Admixture_Eurogenes_Report.html

IF EXIST out\Admixture_Dodecad_Report.html DEL /Q /F out\admix_result_dv3.txt
IF EXIST out\Admixture_Globe13_Report.html DEL /Q /F out\admix_result_globe13.txt
IF EXIST out\Admixture_Eurogenes_Report.html DEL /Q /F out\admix_result_eurogenes36.txt

)

IF EXIST genotype.txt DEL /Q /F genotype.txt
IF EXIST snps.txt DEL /Q /F snps.txt
IF EXIST bam_complete_sorted.bam DEL /Q /F bam_complete_sorted.bam
IF EXIST bam_complete_sorted.bam.bai DEL /Q /F bam_complete_sorted.bam.bai
IF EXIST ystr.filters DEL /Q /F ystr.filters
IF EXIST bed.a DEL /Q /F bed.a
IF EXIST bam_strs.aligned.stats DEL /Q /F bam_strs.aligned.stats
IF EXIST bam_strs_sorted.bam DEL /Q /F bam_strs_sorted.bam
IF EXIST bam_strs_sorted.bam.bai DEL /Q /F bam_strs_sorted.bam.bai
IF EXIST bam_ystrs.allelotype.stats DEL /Q /F bam_ystrs.allelotype.stats
IF EXIST bam_ystrs.vcf DEL /Q /F bam_ystrs.vcf
IF EXIST bam_strs.aligned.bam DEL /Q /F bam_strs.aligned.bam
IF EXIST bam_out_variants.vcf DEL /Q /F bam_out_variants.vcf

explorer out
echo.
echo All Tasks Completed. Please find results in out subfolder.
echo Also check the logs/info in this window for errors (if any).
goto END
:NOPARAM
echo.
echo  Syntax:
echo     console_bam ^<bam-file^>
echo.
:END
pause