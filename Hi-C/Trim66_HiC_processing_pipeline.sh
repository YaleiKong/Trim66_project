#!/bin/bash
#SBATCH --job-name=sbatch
#SBATCH --partition=cpu
#SBATCH -n 1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

InputFolder=raw_data
outputFolder=output

for input in $InputFolder/*/*_R1.fq.gz
do
	sample_name=${input##*/}
    sample_name=${sample_name%%_R1.fq.gz}
    echo $sample_name

    fileR1=$input
    fileR2=${input/R1/R2}
    echo $fileR1
    echo $fileR2

    mkdir $outputFolder/$sample_name

    sbatch pairs_to_cooler.slurm $sample_name $outputFolder $fileR1 $fileR2
done


#!/bin/bash
#SBATCH --job-name=trim_align_pair
#SBATCH --partition=cpu
#SBATCH -n 35
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
thread=35
##########################
#### fastqc
##########################
module load fastqc
fastqc $3 $4 -o $2/$1

##########################
#### cutadapt
##########################
source activate cutadapt

cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o $2/$1/${1}_trimmed.R1.fastq.gz -p $2/$1/${1}_trimmed.R2.fastq.gz \
    -j 0 \
    $3 $4 &> $2/$1/"${1}.cutadapt.log"
# -a ADAPTER, --adapter ADAPTER Sequence of an adapter ligated to the 3' end (paired data: of the first read). 
# -A ADAPTER            3' adapter to be removed from second read in a pair.
# -o FILE, --output FILE Write trimmed reads to FILE. FASTQ or FASTA format is chosen depending on input. Summary report is sent to standard output. Use '{name}' for demultiplexing (see docs). Default: write to standard output
# -p FILE, --paired-output FILE Write second read in a pair to FILE.
# -j CORES, --cores CORES  Number of CPU cores to use. Use 0 to auto-detect. Default: 1
conda deactivate

##########################
#### bwa
#########################
module load bwa
module load samtools

index_prefix=/lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/GRCm38_primary_assembly_genome_related/bwa_index_Gencode_GRCm38_primary_assembly_genome/GRCm38.primary_assembly.genome.Gencode.fa

bwa mem -SP5M $index_prefix -t $thread $2/$1/${1}_trimmed.R1.fastq.gz $2/$1/${1}_trimmed.R2.fastq.gz > $2/$1/${1}.sam 2> $2/$1/${1}.log
# The -SP option is used to ensure the results are equivalent to that obtained by running bwa mem on each mate separately, while retaining the right formatting for paired-end reads. This option skips a step in bwa mem that forces alignment of a poorly aligned read given an alignment of its mate with the assumption that the two mates are part of a single genomic segment.
# The -5 option is used to report the 5' portion of chimeric alignments as the primary alignment. 
    # In Hi-C experiments, when a mate has chimeric alignments, typically, the 5' portion is the position of interest, while the 3' portion represents the same fragment as the mate. 
    # For chimeric alignments, bwa mem reports two alignments: one of them is annotated as primary and soft-clipped, retaining the full-length of the original sequence. 
    # The other end is annotated as hard-clipped and marked as either 'supplementary' or 'secondary'. The -5 option forces the 5'end to be always annotated as primary.
# The -M option is used to annotate the secondary/supplementary clipped reads as secondary rather than supplementary, for compatibility with some public software tools such as picard MarkDuplicates.
# The -t option is used for multi-threading and should not affect the result.

samtools view -bhS $2/$1/${1}.sam > $2/$1/${1}.bam
# -S: In older versions of SAMtools, this flag indicated the input is in SAM format. 
    # However, in more recent versions of SAMtools (1.0 and later), the -S option has been deprecated because the tool automatically detects the input format (SAM or BAM).
# -b: This option tells samtools view to output the result in BAM format. BAM is a binary format that is more compressed than SAM, which is plain text.
# -h, --with-header Include the header in the output.

samtools flagstat $2/$1/${1}.bam

samtools view -q 30 $2/$1/${1}.bam -o $2/$1/${1}.q30.bam
# -q 30: “Filter on read mapping quality (phred scale)": >=30 

samtools sort $2/$1/${1}.q30.bam -o $2/$1/${1}.q30.sort.bam
samtools index $2/$1/${1}.q30.sort.bam

rm $2/$1/${1}_trimmed.R1.fastq.gz $2/$1/${1}_trimmed.R2.fastq.gz $2/$1/${1}.sam $2/$1/${1}.bam 
#########################
### picard 
#########################
source activate picard
#### remove duplicates
picard MarkDuplicates \
      I=$2/$1/${1}.q30.sort.bam  \
      O=$2/$1/${1}.q30.sort.remove_duplicates.bam \
      REMOVE_DUPLICATES=true \
      ASSUME_SORT_ORDER=coordinate \
      M=$2/$1/"${1}_remove_dup_metrics.txt"
# INPUT=String
# I=String                      One or more input SAM or BAM files to analyze. Must be coordinate sorted.  Default value: null. This option may be specified 0 or more times. 
# OUTPUT=File
# O=File                        The output file to write marked records to  Required. 
# REMOVE_DUPLICATES=Boolean     If true do not write duplicates to the output file instead of writing them with appropriate flags set.  Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false} 
# ASSUME_SORT_ORDER=SortOrder
# ASO=SortOrder                 If not null, assume that the input file has this order even if the header says otherwise.  Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate, unknown}  Cannot be used in conjuction with option(s) ASSUME_SORTED (AS)
# METRICS_FILE=File
# M=File                        File to write duplication metrics to  Required. 

conda deactivate

#######################
## pairtools
#######################
chromosize=/lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/GRCm38_primary_assembly_genome_related/GRCm38_primary_assembly_genome_fasta_Gencode/GRCm38.primary_assembly.genome.Gencode.chrom.sizes
source activate pairtools
pairtools parse -c $chromosize --drop-sam --min-mapq 30 -o $2/$1/${1}.parsed.mapq30.pairs.gz $2/$1/${1}.q30.bam --nproc-in $thread --nproc-out $thread
# -c, --chroms-path <chroms_path>   Required Chromosome order used to flip interchromosomal mates: path to a chromosomes file (e.g. UCSC chrom.sizes or similar) whose first column lists scaffold names.
# -o, --output <output> If the path ends with .gz or .lz4, the output is bgzip-/lz4-compressed.
# --drop-sam    If specified, do not add sams to the output
# --min-mapq <min_mapq> The minimal MAPQ score to consider a read as uniquely mapped. Default: 1
# --nproc-in <nproc_in> Number of processes used by the auto-guessed input decompressing command. Default: 3.
# --nproc-out <nproc_out>   Number of processes used by the auto-guessed output compressing command. Default: 8.

pairtools sort -o $2/$1/${1}.parsed.mapq30.sorted.pairs.gz $2/$1/${1}.parsed.mapq30.pairs.gz --nproc $thread --memory 140G
# -o, --output <output>
# --nproc <nproc>   Number of processes to split the sorting work between. Default: 8
# --memory <memory> The amount of memory used by default. Default: 2G.
# --compress-program <compress_program> A binary to compress temporary sorted chunks. Must decompress input when the flag -d is provided. Suggested alternatives: gzip, lzop, lz4c, snzip. If “auto”, then use lz4c if available, and gzip otherwise. Default: auto.

# pairtools markasdup -o $2/$1/${1}_markasdup.pairs.gz $2/$1/${1}.parsed.mapq30.sorted.pairs.gz --nproc-in $thread --nproc-out $thread
pairtools dedup -o $2/$1/${1}.parsed.mapq30.sorted.deduped.pairs.gz $2/$1/${1}.parsed.mapq30.sorted.pairs.gz --output-stats $2/$1/${1}.parsed.mapq30.sorted.deduped.pairs.stats.txt -p $thread

conda deactivate 

#########################
#### pairix
#########################
source activate pairix
gunzip "$2/$1/${1}.parsed.mapq30.sorted.deduped.pairs.gz"
bgzip "$2/$1/${1}.parsed.mapq30.sorted.deduped.pairs"

pairix $2/$1/${1}.parsed.mapq30.sorted.deduped.pairs.gz
conda deactivate

#########################
#### cooler
#########################
chromosize=/lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/GRCm38_primary_assembly_genome_related/GRCm38_primary_assembly_genome_fasta_Gencode/GRCm38.primary_assembly.genome.Gencode.chrom.sizes
fa_file=/lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/GRCm38_primary_assembly_genome_related/GRCm38_primary_assembly_genome_fasta_Gencode/GRCm38.primary_assembly.genome.Gencode.fa

source activate cooler

cooler cload pairix $chromosize:5000 "$2/$1/${1}.parsed.mapq30.sorted.deduped.pairs.gz" "$2/$1/${1}.5000.cool"

cooler zoomify -r 10000,20000,25000,30000,40000,50000,100000,500000,1000000 --balance -o $2/$1/${1}.5000.zoomify.mcool $2/$1/${1}.5000.cool


# export bin information in the specific resolution
cooler dump --header -t bins $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/15000 | cut -f1-3 > $2/$1/${1}.bins.15000.tsv
cooler dump --header -t bins $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/50000 | cut -f1-3 > $2/$1/${1}.bins.50000.tsv
cooler dump --header -t bins $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/100000 | cut -f1-3 > $2/$1/${1}.bins.100000.tsv
cooler dump --header -t bins $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/500000 | cut -f1-3 > $2/$1/${1}.bins.500000.tsv
cooler dump --header -t bins $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/1000000 | cut -f1-3 > $2/$1/${1}.bins.1000000.tsv
# cooler dum: Dump a cooler’s data to a text stream.
# -H, --header  Print the header of column names as the first row. [default: False]
# -t, --table <table>   Which table to dump. [default: pixels]

conda deactivate

########################
### cooltools
########################
source activate cooltools

# Calculate expected-cis Hi-C signal
cooltools expected-cis $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/15000 -o $2/$1/${1}.expected_cis.15000.tsv -p $thread
cooltools expected-cis $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/50000 -o $2/$1/${1}.expected_cis.50000.tsv -p $thread
cooltools expected-cis $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/100000 -o $2/$1/${1}.expected_cis.100000.tsv -p $thread
cooltools expected-cis $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/500000 -o $2/$1/${1}.expected_cis.500000.tsv -p $thread
cooltools expected-cis $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/1000000 -o $2/$1/${1}.expected_cis.1000000.tsv -p $thread

# Calculate expected-cis Hi-C signal
# --smooth 
cooltools expected-cis $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/15000 -o $2/$1/${1}.expected_cis.smoothed.15000.tsv -p $thread --smooth
cooltools expected-cis $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/50000 -o $2/$1/${1}.expected_cis.smoothed.50000.tsv  -p $thread --smooth
cooltools expected-cis $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/100000 -o $2/$1/${1}.expected_cis.smoothed.100000.tsv  -p $thread --smooth

# Calculate expected-cis Hi-C signal
# --aggregate-smoothed
cooltools expected-cis $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/15000 -o $2/$1/${1}.expected_cis.aggregate_smoothed.15000.tsv -p $thread --smooth --aggregate-smoothed
cooltools expected-cis $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/50000 -o $2/$1/${1}.expected_cis.aggregate_smoothed.50000.tsv  -p $thread --smooth --aggregate-smoothed
cooltools expected-cis $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/100000 -o $2/$1/${1}.expected_cis.aggregate_smoothed.100000.tsv  -p $thread --smooth --aggregate-smoothed

# Calculate expected-trans Hi-C signal for trans regions of chromosomal interaction map
cooltools expected-trans $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/15000 -o $2/$1/${1}.expected_trans.15000.tsv -p $thread
cooltools expected-trans $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/50000 -o $2/$1/${1}.expected_trans.50000.tsv -p $thread
cooltools expected-trans $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/100000 -o $2/$1/${1}.expected_tran.100000.tsv  -p $thread

#### compartment (100000)
# Since the orientation of eigenvectors is determined up to a sign, the convention for Hi-C data anaylsis is to orient eigenvectors to be positively correlated with a binned profile of GC content as a ‘phasing track’.
# In humans and mice, GC content is useful for phasing because it typically has a strong correlation at the 100kb-1Mb bin level with the eigenvector.

# the phasing track
cooltools genome gc $2/$1/${1}.bins.100000.tsv $fa_file > $2/$1/${1}.bins.gc.100000.tsv
cooltools genome gc $2/$1/${1}.bins.500000.tsv $fa_file > $2/$1/${1}.bins.gc.500000.tsv
cooltools genome gc $2/$1/${1}.bins.1000000.tsv $fa_file > $2/$1/${1}.bins.gc.1000000.tsv

# Perform eigen value decomposition on a cooler matrix to calculate compartment signal by finding the eigenvector that correlates best with the phasing track.
cooltools eigs-cis --phasing-track $2/$1/${1}.bins.gc.100000.tsv -o $2/$1/${1}.100000 $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/100000
cooltools eigs-cis --phasing-track $2/$1/${1}.bins.gc.500000.tsv -o $2/$1/${1}.500000 $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/500000
cooltools eigs-cis --phasing-track $2/$1/${1}.bins.gc.1000000.tsv -o $2/$1/${1}.1000000 $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/1000000

#### saddle plots (100000)
# Calculate saddle statistics and generate saddle plots for an arbitrary signal track on the genomic bins of a contact matrix.
cooltools saddle $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/100000 $2/$1/${1}.100000.cis.vecs.tsv::E1 $2/$1/${1}.expected_cis.100000.tsv \
-o $2/$1/${1}.saddle.100000 --fig pdf --qrange 0.025 0.975 -n 38 --strength 

cooltools saddle $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/500000 $2/$1/${1}.500000.cis.vecs.tsv::E1 $2/$1/${1}.expected_cis.500000.tsv \
-o $2/$1/${1}.saddle.500000 --fig pdf --qrange 0.025 0.975 -n 38 --strength 

cooltools saddle $2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/1000000 $2/$1/${1}.1000000.cis.vecs.tsv::E1 $2/$1/${1}.expected_cis.1000000.tsv \
-o $2/$1/${1}.saddle.1000000 --fig pdf --qrange 0.025 0.975 -n 38 --strength 

#### insulation (resolution 50000)
cooltools insulation --threshold Li \
-o $2/$1/${1}.insulation.50000.window.150000.tsv \
$2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/50000 150000 -p $thread

cooltools insulation --threshold Li \
-o $2/$1/${1}.insulation.50000.window.250000.tsv \
$2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/50000 250000 -p $thread

cooltools insulation --threshold Li \
-o $2/$1/${1}.insulation.50000.window.500000.tsv \
$2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/50000 500000 -p $thread

cooltools insulation --threshold Li \
-o $2/$1/${1}.insulation.50000.window.1000000.tsv \
$2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/50000 1000000 -p $thread

cooltools insulation --threshold Li \
-o $2/$1/${1}.insulation.25000.window.75000.tsv \
$2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/25000 75000 -p $thread 

cooltools insulation --threshold Li \
-o $2/$1/${1}.insulation.25000.window.125000.tsv \
$2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/25000 125000 -p $thread 

cooltools insulation --threshold Li \
-o $2/$1/${1}.insulation.25000.window.250000.tsv \
$2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/25000 250000 -p $thread 

cooltools insulation --threshold Li \
-o $2/$1/${1}.insulation.25000.window.500000.tsv \
$2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/25000 500000 -p $thread 

#### dot (15000)
cooltools dots -p $thread \
-o $2/$1/${1}.dots.15000.tsv \
$2/$1/${1}.valid_pairs.5000.zoomify.mcool::resolutions/15000 $2/$1/${1}.expected_cis.15000.tsv

conda deactivate


