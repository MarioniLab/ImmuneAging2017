#
# Recipes for the result files of the scRNAseq analysis of CD4 cells of mice
#


#### Applications used
mapper = /path/to/gsnap/bin/bin/
samtools = /path/to/samtools-1.2
htseq_count = /path/to/HTSeq-0.6.1p1/scripts

#### Other variables
memlimit = 32000

#### Filenames of data sources 
ERCC = ERCC92
cur_folder = SS80
cur_folder_in = SS80
cur_genome = Mus_musculus

#### Data folders
# Reference and annotation
B6_reference = ./data/genome/Mus_musculus.GRCm38.dna.primary_assembly.fa
CAST_reference = ./data/genome/CAST_EiJ.chromosomes.unplaced.gt2k.fa
annotation_B6 = ./data/annotation/Mus_musculus.GRCm38.80.gtf
annotation_CAST = ./data/annotation/CAST_EiJ.consensus.protein_coding.gtf

ERCC_reference = ./data/ERCC/${ERCC}.fa
ERCC_annotation = ./data/ERCC/${ERCC}.gtf

# Merged files
merge = ./data/merge

# Input
cur_input = /path/to/fastq_files/${cur_genome}/${cur_folder_in}

# Output
cur_results = /path/to/bam_files/${cur_genome}/${cur_folder}

# Results
results = /path/to/results/${cur_genome}/${cur_folder}

#### Index
gsnap_index = /path/to/gsnap/database

#### Mapping
mapped_reads = $(addprefix ${cur_results}/,$(patsubst %1.fq.gz,%.bam,$(notdir ${fastq_files}))) 
fastq_files = $(shell cat all-files-fastq.txt)

#### Counting
read_count = $(addprefix ${results}/,$(patsubst %.bam,%.txt,$(notdir ${bam_files})))
bam_files = $(shell cat all-files-bam.txt)

# Rules to download data files for ERCC 
.PHONY: download_ERCC
download_ERCC:
	curl 'http://tools.lifetechnologies.com/content/sfs/manuals/ERCC92.zip' -o ./data/ERCC/${ERCC}.zip
	unzip ./data/ERCC/${ERCC}.zip -d ./data/ERCC
	rm ./data/ERCC/${ERCC}.zip


##### Targets with prerequisites
# Rules to concatenate genome and ERCC files

${merge}/Mus_musculus_ERCC.gtf: ${annotation_B6} ${ERCC_annotation}
	cat $< ${ERCC_annotation} > $@
${merge}/Mus_musculus_castaneus_ERCC.gtf: ${annotation_CAST} ${ERCC_annotation}
	cat $< ${ERCC_annotation} > $@

${merge}/Mus_musculus_ERCC.fa: ${B6_reference} ${ERCC_reference}
	cat $< ${ERCC_reference} > $@
${merge}/Mus_musculus_castaneus_ERCC.fa: ${CAST_reference} ${ERCC_reference}
	cat $< ${ERCC_reference} > $@


#### Indexing prior to mapping

.PHONY: B6-gmap-build
B6-gmap-build: ${merge}/Mus_musculus_ERCC.fa
	${bsub} -M 32000 -R 'rusage[mem=32000]' "${mapper}/gmap_build -d Mus_musculus_ERCC -D ${gsnap_index} -k 15 $<"

.PHONY: CAST-gmap-build
CAST-gmap-build: ${merge}/Mus_musculus_castaneus_ERCC.fa
	${bsub} -M 32000 -R 'rusage[mem=32000]' "${mapper}/gmap_build -d Mus_musculus_castaneus_ERCC -D ${gsnap_index} -k 15 $<"

#### Mapping
# Generate File for files to map
all-files-fastq.txt: 
	ls ${cur_input}/*1.fq.gz > $@

all-files-bam.txt:
	ls ${cur_results}/*.bam > $@

.PHONY: mapped-reads
mapped-reads: ${mapped_reads}

# For the mapping adjust the reference!
${cur_results}/%.bam: ${cur_input}/%1.fq.gz ${cur_input}/%2.fq.gz  
	bsub -M 32000 -n 10 -R 'select[gpfs]' -R 'rusage[mem=32000]' "${mapper}/gsnap -A sam -B 5 -t 10 -Q -D ${gsnap_index} -s Mus_musculus_ERCC -d Mus_musculus_ERCC --gunzip $< --gunzip $(word 2, $^) | ${samtools}/samtools view -bS -o $@ -"


#### Read counting
.PHONY: read-count
read-count: ${read_count}

${results}/%.txt: ${cur_results}/%.bam
	bsub -M 8000 -R 'rusage[mem=8000]' -R 'select[gpfs]' "${samtools}/samtools view $< | ${htseq_count}/htseq-count -s no - ${merge}/Mus_musculus_ERCC.gtf > $@"


