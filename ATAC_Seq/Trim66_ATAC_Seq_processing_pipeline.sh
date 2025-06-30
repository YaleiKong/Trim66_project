##########################################################################################################################
sbatch_fastp_bowtie2_samtools_picard_genrich.slurm
##########################################################################################################################
#!/bin/bash
#SBATCH --job-name=sbatch 
#SBATCH --partition=cpu
#SBATCH -n 1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

inputFolder="raw_data"
outputFolder="processed_fastp_bowtie2"

for input in $inputFolder/*R1.fq.gz
do
    cell=${input##*/}
    cell=${cell%%_R1.fq.gz}
    echo $cell

    fileR1=$input
    fileR2=${input/_R1/_R2}
    echo $fileR1
    echo $fileR2   

    mkdir $outputFolder/$cell 
    sbatch ./fastp_bowtie2_samtools_picard_genrich.slurm $cell $outputFolder $fileR1 $fileR2
done

##########################################################################################################################
fastp_bowtie2_samtools_picard_genrich.slurm
##########################################################################################################################
#!/bin/bash
#SBATCH --job-name=trim_align_genrich
#SBATCH --partition=cpu 
#SBATCH -n 5 
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3

#########
# fastp
#########
source activate fastp

fastp --in1 $3 --in2 $4 --out1 $2/$1/"${1}_R1.trimmed.fq.gz" --out2 $2/$1/"${1}_R2.trimmed.fq.gz" \
 -g -q 20 -u 40 -l 30 -g -h $2/$1/"${1}.fastp.html" &> $2/$1/"${1}.fastp.log"

# removes the polyG tails that arise from lack of signal in NextSeq/NovaSeq technologies. 
# -g: It is enabled for Nextseq/Novaseq data by default, and you can specify -g to enable it for any data, or specify -G to disable it.  
# -q: Quality threshold per base required. Default: 15, which means that a Phred quality score of at least 15 is required
# -u: Percent of bases allowed to be below the quality threshold to keep the read (0~100). Default 40 means 40% bases can fail the quality threshold. If more bases fail, the read is removed.
# -l 50: this specifies that if a read is shorter than 50 basepairs after all filters, it should be removed.

conda deactivate

#########
# QC
#########
source activate fastp
fastqc $3 $4
fastqc  $2/$1/"${1}_R1.trimmed.fq.gz"  $2/$1/"${1}_R2.trimmed.fq.gz" 
conda deactivate

#########
# bowtie2
#########
source activate bowtie2
bowtie2 --end-to-end --very-sensitive \
 -x /lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/GRCm38_primary_assembly_genome_related/bowtie2_index_Gencode_GRCm38_primary_assembly_genome/genome_GRCm38_primary_assembly_Gencode \
-1 $2/$1/"${1}_R1.trimmed.fq.gz"  -2 $2/$1/"${1}_R2.trimmed.fq.gz" \
-X 2000 -p 10 -S $2/$1/"${1}.fastp.bowtie2.sam" &> $2/$1/"${1}.bowtie2.log"
conda deactivate

# --end-to-end       entire read must align; no clipping (on)
#    OR
# --local            local alignment; ends might be soft clipped (off); trim some read characters to maximize the alignment score

# Presets:                 Same as:
#   For --end-to-end:
#    --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50
#    --fast                 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50
#    --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)
#    --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
# 
#   For --local:
#    --very-fast-local      -D 5 -R 1 -N 0 -L 25 -i S,1,2.00
#    --fast-local           -D 10 -R 2 -N 0 -L 22 -i S,1,1.75
#    --sensitive-local      -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)
#    --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50

# Paired-end:
      # -I/--minins <int>  minimum fragment length (0)
      # -X/--maxins <int>  maximum fragment length (500)
      # --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)
      # --no-mixed         suppress unpaired alignments for paired reads
      # --no-discordant    suppress discordant alignments for paired reads
      # --dovetail         concordant when mates extend past each other
      # --no-contain       not concordant when one mate alignment contains other
      # --no-overlap       not concordant when mates overlap at all
# Performance:
      # -p/--threads <int> number of alignment threads to launch (1)

##############################
# samtools (exclude chrM)
##############################
module load samtools
# single core

# find chr([1-9]|1[0-9]|X|Y) and =; add @ SQ lines present in the header
awk 'BEGIN {OFS = FS = "\t"} /^@/ || (/chr([1-9]|1[0-9]|X|Y)/ && /=/)' $2/$1/"${1}.fastp.bowtie2.sam" > $2/$1/"${1}.fastp.bowtie2.excluede_chrM.sam"
# BEGIN {OFS = FS = "\t"}: This initializes awk by setting the Output Field Separator (OFS) and the Field Separator (FS) to a tab character. This is done because SAM files are tab-delimited.
# /^@/: This matches any line starting with "@", which in SAM files are header lines. These lines are included in the output.
# /chr([1-9]|1[0-9]|X|Y)/: This part uses a regular expression to match lines where the chromosome field is 'chr1' to 'chr19', 'chrX', or 'chrY'. 
# && /=/: This checks for lines containing the '=' character, which in the context of a SAM file typically represents a match to the reference chromosome. 

samtools view -bS -F 4 -q 20 $2/$1/"${1}.fastp.bowtie2.excluede_chrM.sam"  -o $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.bam"
# -S: In older versions of SAMtools, this flag indicated the input is in SAM format. However, in more recent versions of SAMtools (1.0 and later), the -S option has been deprecated because the tool automatically detects the input format (SAM or BAM).
# -b: This option tells samtools view to output the result in BAM format. BAM is a binary format that is more compressed than SAM, which is plain text.
# -F 4: exclude segment unmapped
# -q 20: “Filter on read mapping quality (phred scale)": >=20

samtools sort $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.bam" -o $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.bam"

samtools index $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.bam"

################################
#### picard to remove duplicates
################################
source activate picard

#### remove duplicates
picard MarkDuplicates \
      I=$2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.bam"  \
      O=$2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam" \
      REMOVE_DUPLICATES=true \
      ASSUME_SORT_ORDER=coordinate \
      M=$2/$1/"${1}_remove_dup_metrics.txt"

#### check Insert Sizes
# picard CollectInsertSizeMetrics \
#       I=$2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam" \
#       O=$2/$1/"${1}_insert_size_metrics.txt" \
#       H=$2/$1/"${1}_insert_size_histogram.pdf" \
#       M=0.5
# Note: If processing a small file, set the minimum percentage option (M) to 0.5, otherwise an error may occur.

conda deactivate


################################
#### genrich 
################################
module load samtools 

samtools sort -n $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.bam" -o $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sortedByQueryname.bam"
# sorted by queryname (samtools sort -n) prepared for genrich


source activate genrich

#### genrich to remove chrM and duplicates
Genrich  -t $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sortedByQueryname.bam" \
 -o $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sortedByQueryname.genrich.narrowpeak" \
 -j  \
 -y  -r  -e chrM  \
 -p 0.01 -q 0.05 \
 -v \
 -f $2/$1/"${1}.fastp.bowtie2.excluede_chrM.mapped.Q20.sortedByQueryname..genrich.log" \
 &> $2/$1/"${1}.genrich.log"

# -j    ATAC-seq mode (must be specified)
# -y    Keep unpaired alignments (def. false)
# -r    Remove PCR duplicates
# -e <arg>  Chromosomes (reference sequences) to exclude. Can be a comma-separated list, e.g. -e chrM,chrY.
# -v    Option to print status updates/counts to stderr
# -f  log
# -d 100  expand cut sites to the given length (def. 100bp)
# -q <float>    Maximum q-value (FDR-adjusted p-value) for peak calling (default 0.05). An unadjusted p-value threshold can be used instead with -p <float>.
    # If a -q threshold is not specified, q-values are not calculated (reported as -1). 
# -a <float>    Minimum area under the curve (total significance) for a peak (default 20.0). Increasing this value results in fewer but higher confidence peaks.
# -v    verbose mode

conda deactivate

##########################################################################################################################
sbatch_bam_to_bw.slurm
##########################################################################################################################
#!/bin/bash
#SBATCH --job-name=sbatch 
#SBATCH --partition=cpu
#SBATCH -n 1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

inputFolder="raw_data"
outputFolder="processed_fastp_bowtie2"

for input in $outputFolder/*/*.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam
do
    cell=${input##*/}
    cell=${cell%%.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam}
    echo $cell

    sbatch ./bam_to_bw.slurm $outputFolder $input $cell 
done

##########################################################################################################################
bam_to_bw.slurm
##########################################################################################################################
#!/bin/bash
#SBATCH --job-name=bamCoverage 
#SBATCH --partition=cpu
#SBATCH -n 10
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load samtools 
samtools index $2 

module load miniconda3
source activate deeptools

bamCoverage -b $2 \
--normalizeUsing BPM \
--effectiveGenomeSize 2652783500 \
-o "$1/$3/${3}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bw" \
--numberOfProcessors 10 \
--binSize 10 --smoothLength 30 \
--ignoreForNormalization chrX \
--extendReads \
--centerReads

# https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
# --normalizeUsing BPM bins per million mapped reads (BPM)

conda deactivate

#########################################################################################################################
# calculate the TSS 
#########################################################################################################################
# https://github.com/jiananlin/TSS_enrichment_score_calculation
sbatch_TSS_enrichment_score_calculation_1.slurm
sbatch_TSS_enrichment_score_calculation_2.slurm
sbatch_TSS_enrichment_score_calculation_3.slurm
sbatch_TSS_enrichment_score_calculation_4.slurm





















#########################################################################################################################
calculate_FRiP.slurm
#########################################################################################################################
#!/bin/bash
#SBATCH --job-name=calculate_FRiP
#SBATCH --partition=cpu 
#SBATCH -n 20
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3

module load samtools
source activate bedtools

for sample in Trim66*
do
echo "${sample}_total_mapped_reads" >> total_reads_and_peak_reads.txt

samtools view -c -F 260 "${sample}/${sample}_fastp_bowtie2.mapped.Q20.excluede_chrM.sorted.remove_duplicates.bam" >> total_reads_and_peak_reads.txt
# to count the number of reads in each BAM file that are not secondary alignments (flag 256) and not unmapped (flag 4), then appends this count to the file.

cut -f1,2,3 "${sample}/${sample}_fastp_bowtie2.mapped.Q20.sortedByQueryname.genrich.sorted.narrowPeak" > "${sample}_peaks.bed"

echo "${sample}_peaks_reads" >> total_reads_and_peak_reads.txt

bedtools coverage -sorted -a "${sample}_peaks.bed" -b "${sample}/${sample}_fastp_bowtie2.mapped.Q20.excluede_chrM.sorted.remove_duplicates.bam" | awk '{sum+=$4} END {print sum}' >> total_reads_and_peak_reads.txt  
# It calculates the read coverage for each peak region using bedtools coverage, sums the coverage using awk, and then appends the total coverage to the file.

rm "${sample}_peaks.bed"
done

##########################################################################################################################
calculate_reads_coverage_of_enhancers.slurm
##########################################################################################################################
#!/bin/bash
#SBATCH --job-name=calculate_enhancer_reads_coverage
#SBATCH --partition=cpu 
#SBATCH -n 10
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load samtools


for input in 4*/*.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam
do
    cell=${input##*/}
    cell=${cell%%.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam}
    
    bedtools coverage -a /lustre/home/acct-medlqian/medlqian-loop3/database/enhancer/OR_63_enhancers_and_1_J_element_and_2_TAAR_enhancers.bed \
    -b $input -hist > "${cell}.enhancer.coverage.txt"

    # output:
    # chromosome, start position, end position, feature name, coverage count, feature length, total length of the region, and the coverage fraction.

    # -sorted (unsuitable)
    # -sorted -k1,1 -k2,2n bed_file:
    # The -k option specifies a sort key.
    # -k1,1 means "sort primarily by the first field". (lexicographical)
    # -k2,2n means "then sort numerically by the second field".
    # The n in 2n specifies a numerical sort, as opposed to a lexicographical sort.
    # samtools sort numerically by default for chromosomes and then lexicographically for other contigs.

    ### -counts, -d, -mean, and -hist are all mutually exclusive options.
    # -hist: Report a histogram of coverage for each feature in A as well as a summary histogram for _all_ features in A.
        # Output (tab delimited) after each feature in A:
        # 1) depth
        # 2) # bases at depth
        # 3) size of A
        # 4) % of A at depth
    # -d: Report the depth at each position in each A feature. Positions reported are one based. Each position and depth follow the complete A feature.
    # -counts: Only report the count of overlaps, don’t compute fraction, etc. Restricted by -f and -r.

    ### fraction
    # -f    Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp).
    # -F    Minimum overlap required as a fraction of B. Default is 1E-9 (i.e., 1bp).

    ### reciprocal
    # -r    Require that the fraction of overlap be reciprocal for A and B.
        # In other words, if -f is 0.90 and -r is used, this requires that B overlap at least 90% of A and that A also overlaps at least 90% of B.
    # -e    Require that the minimum fraction be satisfied for A _OR_ B. 
        # In other words, if -e is used with -f 0.90 and -F 0.10 this requires that either 90% of A is covered OR 10% of B is covered. Without -e, both fractions would have to be satisfied.

    ### strand
    # -s    Force “strandedness”. That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
    # -S    Require different strandedness. That is, only report hits in B that overlap A on the _opposite_ strand. By default, overlaps are reported without respect to strand.

    grep -v '^all' "${cell}.enhancer.coverage.txt" > "${cell}.enhancer.coverage.modified.txt"

    # the sum of fraction grouped by enhancer
    awk '{fraction[$4] += ($5 * $6) / ($3 - $2)} END {for (feature in fraction) print feature, fraction[feature]}'  "${cell}.enhancer.coverage.modified.txt" > "${cell}.enhancer.coverage.modified.sum.txt"
    # $4, $5, $7, $2, and $3 refer to the 4th (feature name), 5th (coverage count), 7th (feature length), 2nd (start position), and 3rd (end position) columns of your file, respectively.
    # fraction[$4] += ($5 * $7) / ($3 - $2) calculates the fraction for each line and accumulates it in an array fraction indexed by the feature name ($4).

    ### normalize by the number of all mapped reads
    input_library_size=$(samtools view -c -F 4 $input) 
    echo $input_library_size
    scale_factor=$(echo "scale=6; 1000000/$input_library_size" | bc)
    echo $scale_factor
    awk -v scale_factor="$scale_factor" 'BEGIN {OFS="\t"} {print $0, $2 * scale_factor}' "${cell}.enhancer.coverage.modified.sum.txt" > "${cell}.enhancer.coverage.modified.sum.RPM.txt"
done


# find the longest length of enhancer
# enhancer_length_baseline=$(awk 'BEGIN {max_length = 0} {len = $3 - $2; if (len > max_length) max_length = len} END {print max_length}' /lustre/home/acct-medlqian/medlqian-loop3/database/enhancer/OR_63_enhancers_and_1_J_element_and_2_TAAR_enhancers.bed)
# 671

### normalize by the length of feature

##########################################################################################################################
plot_reads_coverage_of_enhancers.R
##########################################################################################################################

##########################################################################################################################
bedtools_intersect_with_enhancer.slurm
##########################################################################################################################
#!/bin/bash
#SBATCH --job-name=bedtools_intersect_enhancer
#SBATCH --partition=cpu 
#SBATCH -n 20
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate bedtools
enhancer=/lustre/home/acct-medlqian/medlqian-loop3/database/enhancer/OR_63_enhancers_and_1_J_element_and_2_TAAR_enhancers.bed

bedtools intersect -wa -wb -a 431_Het/431_Het.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak -b $enhancer \
> 431_Het/431_Het.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak.enhancer_intersect.bed

bedtools intersect -wa -wb -a 435_Homo/435_Homo.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak -b $enhancer \
> 435_Homo/435_Homo.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak.enhancer_intersect.bed

bedtools intersect -wa -wb -a 455_Het/455_Het.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak -b $enhancer \
> 455_Het/455_Het.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak.enhancer_intersect.bed

bedtools intersect -wa -wb -a 464_Homo/464_Homo.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak -b $enhancer \
> 464_Homo/464_Homo.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak.enhancer_intersect.bed

conda deactivate

wc -l 431_Het/431_Het.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak.enhancer_intersect.bed
# 57
wc -l 435_Homo/435_Homo.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak.enhancer_intersect.bed
# 59
wc -l 455_Het/455_Het.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak.enhancer_intersect.bed
# 56 
wc -l 464_Homo/464_Homo.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak.enhancer_intersect.bed
# 59

##########################################################################################################################
sbatch_peak_heatmap_profile.slurm
##########################################################################################################################
#!/bin/bash
#SBATCH --job-name=sbatch 
#SBATCH --partition=cpu
#SBATCH -n 1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

outputFolder="processed_fastp_bowtie2"

for input in $outputFolder/*/*.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bw
do
    cell=${input##*/}
    cell=${cell%%.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bw}
    echo $cell

    sbatch ./peak_heatmap_profile.slurm $outputFolder $input $cell 
done

##########################################################################################################################
peak_heatmap_profile.slurm
##########################################################################################################################
#!/bin/bash
#SBATCH --job-name=peak_heatmap
#SBATCH --partition=cpu
#SBATCH -n 20
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
source activate deeptools

computeMatrix reference-point -S $2 \
-R "$1/$3/$3.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak" \
--referencePoint center -b 1000 -a 1000 -bs 10 \
-out "$1/$3/$3.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak.matrix.gz"

plotHeatmap -m "$1/$3/$3.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak.matrix.gz" \
-out "$1/$3/$3.fastp.bowtie2.mapped.Q20.sortedByQueryname.genrich.narrowpeak.matrix.gz.pdf" \
--refPointLabel "peak center" --samplesLabel  $3 --colorMap RdYlBu_r 
# --heatmapWidth 7 --yMax 25 -max 38

conda deactivate


##########################################################################################################################
merge.slurm
##########################################################################################################################
#!/bin/bash
#SBATCH --job-name=samtools_merge
#SBATCH --partition=cpu 
#SBATCH -n 1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load samtools 
samtools merge /lustre/home/acct-medlqian/medlqian-loop3/data/ATAC_Seq/20240712_Trim66_KO_2groups_ATAC_Seq_20240801_6G/processed_fastp_bowtie2/431_Het/431_Het.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam \
/lustre/home/acct-medlqian/medlqian-loop3/data/ATAC_Seq/20240712_Trim66_KO_2groups_ATAC_Seq_20240801_6G/processed_fastp_bowtie2/455_Het/455_Het.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam \
-o /lustre/home/acct-medlqian/medlqian-loop3/data/ATAC_Seq/20240712_Trim66_KO_2groups_ATAC_Seq_20240801_6G/merged_data/het/het.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam

samtools merge /lustre/home/acct-medlqian/medlqian-loop3/data/ATAC_Seq/20240712_Trim66_KO_2groups_ATAC_Seq_20240801_6G/processed_fastp_bowtie2/435_Homo/435_Homo.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam \
/lustre/home/acct-medlqian/medlqian-loop3/data/ATAC_Seq/20240712_Trim66_KO_2groups_ATAC_Seq_20240801_6G/processed_fastp_bowtie2/464_Homo/464_Homo.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam \
-o /lustre/home/acct-medlqian/medlqian-loop3/data/ATAC_Seq/20240712_Trim66_KO_2groups_ATAC_Seq_20240801_6G/merged_data/homo/homo.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam

##########################################################################################################################
sbatch_bam_to_bw.slurm
##########################################################################################################################
#!/bin/bash
#SBATCH --job-name=sbatch 
#SBATCH --partition=cpu
#SBATCH -n 1
#SBATCH --output=%j.out
#SBATCH --error=%j.err


outputFolder="merged_data"

for input in $outputFolder/*/*.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam
do
    cell=${input##*/}
    cell=${cell%%.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam}
    echo $cell

    sbatch ./bam_to_bw.slurm $outputFolder $input $cell 
done

##########################################################################################################################
bam_to_bw.slurm
##########################################################################################################################
#!/bin/bash
#SBATCH --job-name=bamCoverage 
#SBATCH --partition=cpu
#SBATCH -n 10
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load samtools 
samtools index $2 

module load miniconda3
source activate deeptools

# https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
bamCoverage -b $2 \
--normalizeUsing BPM \
--effectiveGenomeSize 2652783500 \
-o "$1/$3/${3}.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bw" \
--numberOfProcessors 10 \
--binSize 10 --smoothLength 30 \
--ignoreForNormalization chrX \
--extendReads \
--centerReads

conda deactivate

##########################################################################################################################
calculate_enhancers_rpm.slurm
##########################################################################################################################
#!/bin/bash
#SBATCH --job-name=calculate_enhancer_reads_coverage
#SBATCH --partition=cpu 
#SBATCH -n 30
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load miniconda3
module load samtools
source activate bedtools

for input in merged_data/*/*.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam
do
    cell=${input##*/}
    cell=${cell%%.fastp.bowtie2.excluede_chrM.mapped.Q20.sorted.remove_duplicates.bam}
    echo $cell 

    bedtools coverage -a /lustre/home/acct-medlqian/medlqian-loop3/database/enhancer/OR_63_enhancers_and_1_J_element_and_2_TAAR_enhancers.bed \
    -b $input -hist > "merged_data/${cell}/${cell}.enhancer.coverage.txt"

    grep -v '^all' "merged_data/${cell}/${cell}.enhancer.coverage.txt" > "merged_data/${cell}/${cell}.enhancer.coverage.modified.txt"

    # the sum of fraction grouped by enhancer
    awk '{fraction[$4] += ($5 * $6) / ($3 - $2)} END {for (feature in fraction) print feature, fraction[feature]}'  "merged_data/${cell}/${cell}.enhancer.coverage.modified.txt" > "merged_data/${cell}/${cell}.enhancer.coverage.modified.sum.txt"
    # $4, $5, $7, $2, and $3 refer to the 4th (feature name), 5th (coverage count), 7th (feature length), 2nd (start position), and 3rd (end position) columns of your file, respectively.
    # fraction[$4] += ($5 * $7) / ($3 - $2) calculates the fraction for each line and accumulates it in an array fraction indexed by the feature name ($4).

    ### normalize by the number of all mapped reads
    input_library_size=$(samtools view -c -F 4 $input) 
    echo $input_library_size
    scale_factor=$(echo "scale=6; 1000000/$input_library_size" | bc)
    echo $scale_factor
    awk -v scale_factor="$scale_factor" 'BEGIN {OFS="\t"} {print $0, $2 * scale_factor}' "merged_data/${cell}/${cell}.enhancer.coverage.modified.sum.txt" > "merged_data/${cell}/${cell}.enhancer.coverage.modified.sum.RPM.txt"
done


# find the longest length of enhancer
# enhancer_length_baseline=$(awk 'BEGIN {max_length = 0} {len = $3 - $2; if (len > max_length) max_length = len} END {print max_length}' /lustre/home/acct-medlqian/medlqian-loop3/database/enhancer/OR_63_enhancers_and_1_J_element_and_2_TAAR_enhancers.bed)
# 671

### normalize by the length of feature



