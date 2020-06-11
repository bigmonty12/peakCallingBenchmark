#!/usr/bin/env nextflow

/*
 * Defines parameters of ref genome and reads
 */
params.reads = "$baseDir/Fastq/*_{1,2}.fastq.gz"
params.genome_index = "/mnt/raid1/genomes/hg38_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index"
params.outdir = 'results'
params.goldPositive = "$baseDir/HumanRuns/GM12878/GoldStandards/hg38_ActiveRegions.bed"
params.goldNegative = "$baseDir/HumanRuns/GM12878/GoldStandards/hg38_Heterochrom_nomodifications.bed"
params.genomeInfo = "$baseDir/genome.info"

/*
 * Create `read_pairs_ch` channel that emits tuples containing three
 * elements: the pair ID, the first read-pair file and second read-pair file
 */
Channel
    .fromFilePairs( params.reads)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { ch_raw_reads_fastqc;
           ch_raw_reads_trimgalore }
            

/*
 * Step 1: FastQC
 */
process FastQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: { filename ->
                    filename.endsWith(".zip") ? "zips/$filename" : "$filename"
                }

    input:
    set val(name), file(reads) from ch_raw_reads_fastqc

    output:
    file "*.{zip,html}" into ch_fastq_reports_mqc

    script:
    """
    fastqc -q -t 8 ${name}_1.fastq.gz
    fastqc -q -t 8 ${name}_2.fastq.gz
    """
}

/*
 * Step 2: TrimGalore! 
 */
process trimGalore {
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: {filename ->
                    if (filename.endsWith(".html")) "fastqc/$filename"
                    else if (filename.endsWith(".zip")) "fastqc/zips/$filename"
                    else if (filename.endsWith("trimming_report.txt")) "logs/$filename"
                }

    input:
    set val(name), file(reads) from ch_raw_reads_trimgalore

    output:
    set val(name), file("*.fq.gz") into ch_trimmed_reads
    file "*.txt" into ch_trimgalore_results_mqc
    file "*.{zip,html}" into ch_trimgalore_fastqc_reports_mqc

    script:
    """
    trim_galore --paired --fastqc --gzip -j 4 ${name}_1.fastq.gz ${name}_2.fastq.gz
    """
}

/*
 * Step 3.1: Align with Bowtie2
 */
ch_bowtie_index = Channel
        .fromPath( "${params.genome_index}*", checkIfExists:true )
        .ifEmpty { exit 1, "Bowtie2 index directory not found: " }
process bowtie2 {
    input:
    // path(index) from params.genome_index
    set val(name), file(reads) from ch_trimmed_reads

    output:
    set val(name), file("*.bam") into ch_bowtie_bams

    script:
    """
    bowtie2 \\
        -p 16 \\
        -x $params.genome_index \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --very-sensitive \\
        -X 2000 \\
    | samtools view \\
        -@ 8 \\
        -b \\
        -f 2 \\
        -h \\
        -q 30 \\
        -O BAM \\
        -o ${name}.bam \\
        -
    """
}

/*
 * Step 3.2: Sort BAM by coordinate
 */
process SortBAM {
    publishDir "${params.outdir}/bams/sorted", mode: 'copy',
        saveAs: {filename ->
                    if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
                    else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
                    else if (filename.endsWith(".stats")) "samtools_stats/$filename"
                    else filename
                }
    input:
    set val(name), file(bam) from ch_bowtie_bams

    output:
    set val(name), file("*.sorted.{bam,bam.bai}") into ch_sort_bam
    set val(name), file("*.{flagstat,idxstats,stats}") into ch_sort_bam_mqc

    script:
    """
    samtools sort \\
        -@ 8 \\
        -m 4G \\
        -o ${name}.sorted.bam \\
        $bam
    samtools index ${name}.sorted.bam
    samtools flagstat ${name}.sorted.bam > ${name}.sorted.bam.flagstat
    samtools idxstats ${name}.sorted.bam > ${name}.sorted.bam.idxstats
    samtools stats ${name}.sorted.bam > ${name}.sorted.bam.stats
    """
}

/*
 * Step 3.3: Remove duplicates from BAM
 */
process RemoveDuplicates {
    publishDir "${params.outdir}/bams/dedup", mode: 'copy',
        saveAs: {filename -> filename
                }

    input:
    set val(name), file(bam_files) from ch_sort_bam

    output:
    set val(name), file("*.dedup.{bam,bam.bai}") into ch_dedup_bam

    script:
    """
    picard MarkDuplicates \\
        INPUT=${bam_files[0]} \\
        OUTPUT=${name}.sorted.dedup.bam \\
        ASSUME_SORTED=true \\
        REMOVE_DUPLICATES=true \\
        METRICS_FILE=${name}.MarkDuplicates.metrix.txt \\
        VALIDATION_STRINGENCY=LENIENT \\
        TMP_DIR=tmp
        
    samtools index ${name}.sorted.dedup.bam
    """
}
// Duplicate ch_dedup_bam into each peak caller channel
ch_dedup_bam .into{ch_dedup_bam_bam; ch_dedup_bam_bed; ch_dedup_bam_HMMR; ch_dedup_bam_HMMR_peaks; ch_dedup_bam_bampe; ch_dedup_bam_sort}

/*
 * Step 3.4: Sort BAM by read names
 */
process Genrich {
    publishDir "${params.outdir}/peaks/Genrich", mode: 'copy',
        saveAs: {filename -> 
                    if (filename.endsWith(".narrowPeak")) "narrowPeak/$filename"
                    }
                

    input:
    set val(name), file(bam_files) from ch_dedup_bam_sort

    output:
    set val(name), file("*.narrowPeak") into ch_narrowpeaks_G, ch_np_G
    file(bam_files) into ch_bam_G


    script:
    """
    samtools sort \\
        -n \\
        -@ 8 \\
        -m 4G \\
        -o ${name}.sorted_names.dedup.bam \\
        ${bam_files[0]}


    Genrich \\
        -t ${name}.sorted_names.dedup.bam \\
        -j \\
        -y \\
        -o ${name}.Genrich.narrowPeak
    """
}

ch_np_G .merge(ch_bam_G)
        .set{ ch_pbc_G }

/*
 * Step 4: Call peaks
 */
process CallHMMRATAC {
    publishDir "${params.outdir}/peaks/HMMRATAC", mode: 'copy',
        saveAs: {filename ->
                    if (filename.endsWith(".log")) "logs/$filename"
                    else if (filename.endsWith(".narrowPeak")) "narrowPeak/$filename"
                    else if (filename.endsWith(".bedgraph")) "bedgraphs/$filename"
                    else if (filename.endsWith(".gappedPeak")) "gappedPeak/$filename"
                    else if (filename.endsWith("_summits.bed")) "summits/$filename"
                    else if (filename.endsWith(".model")) "model/$filename"
                }

    maxForks 2

    input: 
    set val(name), file(bam_files) from ch_dedup_bam_HMMR

    output:
    set val(name), file("*.log") into ch_HMMRATAC_mqc
    // set val(name), file( "*summits.bed") into ch_HMMRATAC_summits
    file "*.{bedgraph,gappedPeak,model}" into ch_HMMRATAC_extra
    set val(name), file("*.narrowPeak") into ch_narrowpeaks_H, ch_np_H
    file(bam_files) into ch_bam_H
    file("genome.info") into ch_genomeinfo

    shell:
    '''
    samtools view -H !{bam_files[0]} | \
    perl -ne 'if(/^@SQ.*?SN:(\\w+)\\s+LN:(\\d+)/){print $1,"\\t",$2,"\\n"}' > genome.info

    HMMRATAC \
        -Xmx32G \
        -b !{bam_files[0]} \
        -i !{bam_files[1]} \
        -g genome.info \
        -p True \
        -o !{name}.HMMRATAC

    awk -v s=50 'BEGIN{OFS="\t"} {print $1, $2-s, $3+s, $4, $5}' !{name}.HMMRATAC_summits.bed > !{name}.HMMRATAC.narrowPeak
    echo !{bam_files[0]} !{bam_files[1]}
    '''
}

ch_np_H .merge(ch_bam_H)
        .set{ ch_pbc_H }


process CallPeaksBAMPE {
    publishDir "${params.outdir}/peaks/macs2", mode: 'copy',
        saveAs: {filename -> filename
                }

    input:
    set val(name), file(bam_files) from ch_dedup_bam_bampe

    output:
    set val(name), file("*.narrowPeak") into ch_narrowpeaks_bampe, ch_np_BAMPE
    file(bam_files) into ch_bam_BAMPE

    script:
    """
    macs2 callpeak \\
        -t ${bam_files[0]} \\
        --name ${name}.BAMPE \\
        -f BAMPE \\
        -g hs \\
        --keep-dup all \\
        --nomodel

    echo $bam_files
    """
}

ch_np_BAMPE .merge(ch_bam_BAMPE)
            .set{ ch_pbc_BAMPE }

process CallPeaksBAM {
    publishDir "${params.outdir}/peaks/macs2", mode: 'copy',
        saveAs: {filename -> filename
                }

    input:
    set val(name), file(bam_files) from ch_dedup_bam_bam

    output:
    set val(name), file("*.narrowPeak") into ch_narrowpeaks_bam, ch_np_BAM
    file(bam_files) into ch_bam_BAM

    script:
    """
    macs2 callpeak \\
    -t ${bam_files[0]} \\
    --name ${name}.BAM \\
    -f BAM \\
    -g hs \\
    --keep-dup all \\
    --nomodel \\
    --shift -100 \\
    --extsize 200

    echo $bam_files
    """
}

ch_np_BAM .merge(ch_bam_BAM)
          .set{ ch_pbc_BAM }

process CallPeaksBED {
    publishDir "${params.outdir}/peaks/macs2", mode: 'copy',
        saveAs: {filename -> filename
                }

    input:
    set val(name), file(bam_files) from ch_dedup_bam_bed

    output:
    set val(name), file("*.narrowPeak") into ch_narrowpeaks_bed, ch_np_BED
    file(bam_files) into ch_bam_BED

    script:
    """
    macs2 randsample \\
    -i ${bam_files[0]} \\
    -f BAMPE \\
    -p 100 \\
    -o ${name}.bed

    macs2 callpeak \\
    -t ${name}.bed \\
    --name ${name}.BED \\
    -f BED \\
    -g hs \\
    --keep-dup all \\
    --nomodel \\
    --shift -100 \\
    --extsize 200

    echo $bam_files
    """
}

ch_np_BED .merge(ch_bam_BED)
          .set{ ch_pbc_BED }


// Put all narrowPeaks into one channel
ch_narrowpeaks_bam .mix(ch_narrowpeaks_bed, ch_narrowpeaks_bampe, ch_narrowpeaks_G, ch_narrowpeaks_H)
                   .into { ch_positive_peaks; ch_negative_peaks}

process FindGoldStandardPositiveOverlaps {
    publishDir "${params.outdir}/GoldStandardOverlaps", mode: 'copy',
        saveAs: {filename -> filename
                }

    cache false

    input: 
    set val(name), file(peaks) from ch_positive_peaks

    output:
    file("*positive_overlap.txt") into ch_gold_positive_overlaps

    shell:
    '''
    sort -k1,1 -k2,2n !{peaks} | \\
    bedtools merge -i - > merged-peaks.bed
    bedtools intersect \\
        -a !{params.goldPositive} \\
        -b merged-peaks.bed \\
        -wo | \\
    awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$(NF) }END{print SUM}' > !{peaks}.positive_overlap.txt
    rm merged-peaks.bed
    '''
}

process FindGoldStandardNegativeOverlaps {
    publishDir "${params.outdir}/GoldStandardOverlaps", mode: 'copy',
        saveAs: {filename -> filename
                }

    cache false
    
    input: 
    set val(name), file(peaks) from ch_negative_peaks

    output:
    file("*negative_overlap.txt") into ch_gold_negative_overlaps

    shell:
    '''
    sort -k1,1 -k2,2n !{peaks} | \\
    bedtools merge -i - > merged-peaks.bed 
    bedtools intersect \\
        -a !{params.goldNegative} \\
        -b merged-peaks.bed \\
        -wo | \\
    awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$(NF) }END{print SUM}' > !{peaks}.negative_overlap.txt
    rm merged-peaks.bed
    '''
}


process PositiveOverlapsToFile {
    publishDir "${params.outdir}", mode: 'copy',
        saveAs: {filename -> filename}

    cache false

    input:
    file overlapList from ch_gold_positive_overlaps.toSortedList()

    output:
    file '*.txt' into ch_goldPositiveOverlaps

    shell:
    '''
    echo -e "PeakFile\tBasePairsOverlapped" > overlaps_goldPositive_sorted.txt
    for file in !{overlapList}
    do
        name=$(echo $file | sed -e "s,GoldStandardOverlaps/,," -e "s/.narrowPeak........._overlap.txt//")
        overlap=$(cat $file)
        echo -e "$name\t$overlap" >> overlaps_goldPositive.txt
    done
    sort -k1 overlaps_goldPositive.txt >> overlaps_goldPositive_sorted.txt
    rm overlaps_goldPositive.txt
    '''
}

process NegativeOverlapsToFile {
    publishDir "${params.outdir}", mode: 'copy',
        saveAs: {filename -> filename}

    cache false

    input:
    file overlapList from ch_gold_negative_overlaps.toSortedList()

    output:
    file '*.txt' into ch_goldNegativeOverlaps

    shell:
    '''
    echo -e "PeakFile\tBasePairsOverlapped" > overlaps_goldNegative_sorted.txt
    for file in !{overlapList}
    do
        name=$(echo $file | sed -e "s/.narrowPeak........._overlap.txt//")
        overlap=$(cat $file)
        echo -e "$name\t$overlap" >> overlaps_goldNegative.txt
    done
    sort -k1 overlaps_goldNegative.txt >> overlaps_goldNegative_sorted.txt
    rm overlaps_goldNegative.txt
    '''
}

process GoldCalculations {
    publishDir "${params.outdir}", mode: 'copy',
        saveAs: {filename -> filename}

    cache false

    input:
    file(negatives) from ch_goldNegativeOverlaps
    file(positives) from ch_goldPositiveOverlaps

    output:
    file '*.txt' into ch_calculations

    shell:
    '''
    echo -e "Name\tTP\tFP\tPrecision\tRecall\tFPR\tF1" > measures.txt
    rp=$(cat !{params.goldPositive} | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
    rn=$(cat !{params.goldNegative} | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
    paste overlaps_goldPositive_sorted.txt overlaps_goldNegative_sorted.txt | awk -v rp=$rp -v rn=$rn 'NR >= 2{print $1, $2,$4,$2/($2+$4), $2/rp, $4/rn, 2*((($2*$2)/(rp*($2+$4)))/(($2/($2+$4))+($2/rp)))}' OFS='\t' >> measures.txt
    '''
}

// Put all narrowPeak & bam combos into one channel
ch_pbc_G .mix( ch_pbc_H, ch_pbc_BAM, ch_pbc_BAMPE, ch_pbc_BED)
        .set{ ch_pbc_all }

process PerBaseCoverage {
    publishDir "${params.outdir}/perbasecoverage", mode: 'copy',
        saveAs: {filename -> filename}

    input:
    file(fileList) from ch_pbc_all

    output:
    file '*.txt' into ch_coverageStatistics

    shell:
    '''
    pre=$( echo !{fileList[1]} | cut -d . -f -2)
    bedtools sort -g !{params.genomeInfo} -i !{fileList[1]} > sorted.!{fileList[1]}
    bedtools coverage -a sorted.!{fileList[1]} -b !{fileList[2]} -d -sorted -g !{params.genomeInfo} > ${pre}.per-base-coverage.txt
    rm sorted.!{fileList[1]}
    '''
}

process CoverageStatistics {
    publishDir "${params.outdir}/", mode: 'copy',
        saveAs: {filename -> filename}

    input:
    file(fileList) from ch_coverageStatistics.toSortedList()

    output:
    file '*.txt' into ch_final

    shell:
    '''
    echo -e "Name\t0-2\t3-6\t>6" > coverageStatistics_sorted.txt
    for file in !{fileList}
    do
        echo $file
        pre=$(echo $file | cut -d . -f -2)
        bases=$(cat $file | wc -l)
        bad=$(awk '{if ($(NF) <= 2) { print $(NF) } }' $file | wc -l)
        bad_pct=$(echo $bad $bases | awk '{ print $1/$2 }')
        okay=$(awk '{if ($(NF) > 2 && $(NF) <= 6) { print $(NF) } }' $file | wc -l)
        okay_pct=$(echo $okay $bases | awk '{ print $1/$2 }')
        good=$(awk '{if ($(NF) >= 7) { print $(NF) } }' $file | wc -l)
        good_pct=$(echo $good $bases | awk '{ print $1/$2 }')
        echo -e "${pre}\t${bad}(${bad_pct})\t${okay}(${okay_pct})\t${good}(${good_pct})" >> coverageStatistics.txt
    done
    sort -k1,1 coverageStatistics.txt >> coverageStatistics_sorted.txt
    rm coverageStatistics.txt
   '''
}
