#!/usr/bin/env nextflow

/*
 * Defines parameters of ref genome and reads
 */
params.reads = "$baseDir/Fastq/*_{1,2}.fastq.gz"
params.genome_index = "/mnt/raid1/genomes/hg38_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index"
params.blacklist = "/mnt/raid1/genomes/hg38_genome/hg38-blacklist.bed"
params.outdir = 'results'
params.goldPositive = "$baseDir/HumanRuns/GM12878/GoldStandards/hg38_ActiveRegions.bed"
params.goldNegative = "$baseDir/HumanRuns/GM12878/GoldStandards/hg38_Heterochrom_nomodifications.bed"
params.genomeInfo = "$baseDir/genome.info"

/*
 * Create `read_pairs_ch` channel that emits tuples containing three
 * elements: the pair ID, the first read-pair file and second read-pair file
 */

ch_input = file(params.input)

process CheckDesign {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file design from ch_input

    output:
    file "*.csv" into ch_design_reads_csv

    script:
    """
    check_design.py $design design_reads.csv
    """
}

ch_design_reads_csv
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ] }
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
    [ ! -f ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
    [ ! -f ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
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
    [ ! -f ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
    [ ! -f ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
    trim_galore --paired --fastqc --gzip ${name}_1.fastq.gz ${name}_2.fastq.gz
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
ch_dedup_bam .into { ch_dedup_bam; ch_dedup_bam_sort}
/*
 * Step.3.4
 */

ch_dedup_bam
    .map { it -> [ it[0].split('_')[0..-3].join('_'), it[1] ] }
    .groupTuple(by: [0])
    .map { it ->  [ it[0], it[1].flatten() ] }
    .set { ch_dedup_bam }

process MergeReplicates {
    publishDir "${params.outdir}/bams/dedup", mode: 'copy',
        saveAs: {filename -> filename}

    input:
    set val(name), file(bams) from ch_dedup_bam

    output: 
    set val(name), file("*-merged.sorted.dedup.{bam,bam.bai}") into ch_merged_bams

    script:
    bam_files = bams.findAll { it.toString().endsWith('.bam') }.sort()
    """
    samtools merge ${name}-merged.sorted.dedup.bam ${bam_files.join(" ")}
    samtools index ${name}-merged.sorted.dedup.bam
    """
}

ch_dedup_bam_sort .mix(ch_merged_bams)
                  .into{ ch_dedup_bam_bam; ch_dedup_bam_bed; ch_dedup_bam_HMMR; ch_dedup_bam_HMMR_peaks; ch_dedup_bam_bampe; ch_dedup_bam_sort }
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
    file(bam_files) into ch_bam_GR
    file("*.sorted-names.dedup.bam") into ch_sorted_bam_GR


    shell:
    '''
    samtools sort \\
        -n \\
        -@ 8 \\
        -o !{name}.sorted-names.dedup.bam \\
        !{bam_files[0]}

    qVal=( 0001 001 01 05 1 25 50 )
    for i in "${qVal[@]}"
    do
        Genrich \\
        -t !{name}.sorted-names.dedup.bam \\
        -j \\
        -y \\
        -p 0.$i \\
        -E !{params.blacklist} \\
        -o !{name}.Genrich-p${i}.narrowPeak
    done 
    '''
}

ch_np_G .merge(ch_bam_G)
        .set{ ch_pbc_G }



ch_sorted_bam_GR
    .filter( ~/.*_.*/ )
    .map { it -> tuple(it.getName().split('_')[0], it) }
    .groupTuple()
    .set { ch_sorted_bam_GR }


process GenrichReplicates {
    publishDir "${params.outdir}/peaks/Genrich", mode: 'copy',
        saveAs: {filename -> 
                    if (filename.endsWith(".narrowPeak")) "narrowPeak/$filename"
                    }
    input:
    set val(name), file(bams) from ch_sorted_bam_GR
    

    output: 
    set val(name), file("*.narrowPeak") into ch_narrowpeaks_GR, ch_np_GR

    shell:
    bam_files = bams.findAll { it.toString().endsWith('.bam') }.sort()

    '''
    qVal=( 0001 001 01 05 1 25 50 )
    for i in "${qVal[@]}"
    do
        Genrich \\
            -t !{bam_files.join(",")} \\
            -j \\
            -y \\
            -E !{params.blacklist} \\
            -p 0.${i} \\
            -o !{name}.Genrich.Replicates-p${i}.narrowPeak
    done
    '''
}

ch_bam_GR
    .map { it -> it[0] }
    .filter( ~/^((?!_).)*$/ )
    .map { it -> tuple(it.getName().split('-')[0], it) }
    .set {ch_bam_GR}

ch_np_GR
    .join(ch_bam_GR)
    .set {ch_pbc_GR}
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
    file "*.{bedgraph,model}" into ch_HMMRATAC_extra
    set val(name), file("*.{gappedPeak,narrowPeak}") into ch_narrowpeaks_H, ch_np_H
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
        -e !{params.blacklist} \
        -p True \
        --threshold 0 \
        -o !{name}.HMMRATAC

    cutoffs=(0 05 10 15 20 25 30 35 40 45)
    for val in "${cutoffs[@]}"
    do
        awk -v val=$val '{ if ($13 >= val) { print } }' !{name}.HMMRATAC_peaks.gappedPeak > !{name}.HMMRATAC_peaks-p$val.gappedPeak
        awk -v val=$val '{ if ($5 >= val) { print } }' !{name}.HMMRATAC_summits.bed > !{name}.HMMRATAC_summits-p$val.bed
        awk -v s=50 'BEGIN{OFS="\t"} {print $1, $2-s, $3+s, $4, $5}' !{name}.HMMRATAC_summits-p$val.bed > !{name}.HMMRATAC-p$val.narrowPeak
    done
    rm !{name}.HMMRATAC_peaks.gappedPeak !{name}.HMMRATAC_summits.bed

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

    shell:
    '''
    qVal=( 0001 001 01 05 1 25 50 )
    for i in "${qVal[@]}"
    do
        macs2 callpeak \\
            -t !{bam_files[0]} \\
            --name !{name}.BAMPE-p${i} \\
            -f BAMPE \\
            -g hs \\
            --keep-dup all \\
            --nomodel \\
            -p 0.${i}
    done
    '''
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

    shell:
    '''
    qVal=( 0001 001 01 05 1 25 5 )
    for i in "${qVal[@]}"
    do
        macs2 callpeak \\
            -t !{bam_files[0]} \\
            --name !{name}.BAM-p${i} \\
            -f BAM \\
            -g hs \\
            --keep-dup all \\
            --nomodel \\
            -p 0.${i} \\
            --shift -100 \\
            --extsize 200
    done
    '''
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

    shell:
    '''
    macs2 randsample \\
    -i !{bam_files[0]} \\
    -f BAMPE \\
    -p 100 \\
    -o !{name}.bed

    qVal=( 0001 001 01 05 1 25 5 )
    for i in "${qVal[@]}"
    do
        macs2 callpeak \\
            -t !{name}.bed \\
            --name !{name}.BED-p${i} \\
            -f BED \\
            -g hs \\
            --keep-dup all \\
            --nomodel \\
            -p 0.${i} \\
            --shift -100 \\
            --extsize 200
    done

    '''
}

ch_np_BED .merge(ch_bam_BED)
          .set{ ch_pbc_BED }


// Put all narrowPeaks into one channel
ch_narrowpeaks_bam .mix(ch_narrowpeaks_bed, ch_narrowpeaks_bampe, ch_narrowpeaks_G, ch_narrowpeaks_H, ch_narrowpeaks_GR)
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
    for file in !{peaks}
    do
        awk '{if ($2>=0) print $0}' ${file} > no-negatives.${file}
        sort -k1,1 -k2,2n no-negatives.${file} | \\
        bedtools merge -i - > merged-peaks.bed
        bedtools intersect \\
            -a !{params.goldPositive} \\
            -b merged-peaks.bed \\
            -wo | \\
        awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$(NF) }END{print SUM}' > ${file}.positive_overlap.txt
        rm merged-peaks.bed no-negatives.${file}
    done
    '''
}
ch_gold_positive_overlaps .flatten()
                          .set {ch_gold_positive_overlaps}
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
    for file in !{peaks}
    do
        awk '{if ($2>=0) print $0}' ${file} > no-negatives.${file}
        sort -k1,1 -k2,2n no-negatives.${file} | \\
        bedtools merge -i - > merged-peaks.bed
        bedtools intersect \\
            -a !{params.goldNegative} \\
            -b merged-peaks.bed \\
            -wo | \\
        awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$(NF) }END{print SUM}' > ${file}.negative_overlap.txt
        rm merged-peaks.bed no-negatives.${file}
    done
    '''
}

ch_gold_negative_overlaps .flatten()
                          .set {ch_gold_negative_overlaps}
process PositiveOverlapsToFile {
    publishDir "${params.outdir}/GoldStandardOverlaps", mode: 'copy',
        saveAs: {filename -> filename}

    cache false

    input:
    file(overlapList) from ch_gold_positive_overlaps.toSortedList()

    output:
    file '*.txt' into ch_goldPositiveOverlaps

    shell:
    '''
    echo -e "PeakFile\tBasePairsOverlapped" >> overlaps_goldPositive_sorted.txt
    for file in !{overlapList}
    do
        name=$(echo $file | sed -e "s,GoldStandardOverlaps/,," -e "s/.narrowPeak........._overlap.txt/.summits/" -e "s/.gappedPeak........._overlap.txt/.gapped/")
        overlap=$(cat $file)
        echo -e "$name\t$overlap" >> overlaps_goldPositive.txt
    done
    sort -k1 overlaps_goldPositive.txt >> overlaps_goldPositive_sorted.txt
    rm overlaps_goldPositive.txt
    '''
}

process NegativeOverlapsToFile {
    publishDir "${params.outdir}/GoldStandardOverlaps", mode: 'copy',
        saveAs: {filename -> filename}

    cache false

    input:
    file(overlapList) from ch_gold_negative_overlaps.toSortedList()

    output:
    file '*.txt' into ch_goldNegativeOverlaps

    shell:
    '''
    echo -e "PeakFile\tBasePairsOverlapped" > overlaps_goldNegative_sorted.txt
    for file in !{overlapList}
    do
        name=$(echo $file | sed -e "s,GoldStandardOverlaps/,," -e "s/.narrowPeak........._overlap.txt/.summits/" -e "s/.gappedPeak........._overlap.txt/.gapped/") 
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
    echo -e "Peak.Caller\tTP\tFP\tPrecision\tRecall\tFPR\tF1" > measures.txt
    rp=$(cat !{params.goldPositive} | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
    rn=$(cat !{params.goldNegative} | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
    paste !{positives} !{negatives} | awk -v rp=$rp -v rn=$rn 'NR >= 2{print $1, $2,$4,$2/($2+$4), $2/rp, $4/rn, 2*((($2*$2)/(rp*($2+$4)))/(($2/($2+$4))+($2/rp)))}' OFS='\t' >> measures.txt
    rm !{positives} !{negatives}
    '''
}
ch_pbc_G. into{ ch_pbc_G; ch_print }
// Put all narrowPeak & bam combos into one channel
ch_pbc_G .mix( ch_pbc_GR, ch_pbc_H, ch_pbc_BAM, ch_pbc_BAMPE, ch_pbc_BED)
        .map { it -> [it[1],[it[2]]].combinations() }
        .flatten()
        .collate( 2 )
        .set{ ch_pbc_all }


process PerBaseCoverage {

    input:
    file(fileList) from ch_pbc_all

    output:
    file '*coverageStatistics.txt' into ch_coverageStatistics

    maxForks 6

    shell:
    '''
    pre=$( echo !{fileList[0]} | rev | cut -d . -f2- | rev)
    awk '{if ($2>=0) print $0}' !{fileList[0]} > no-negatives.!{fileList[0]}
    bedtools sort -g !{params.genomeInfo} -i no-negatives.!{fileList[0]} > sorted.!{fileList[0]}
    bedtools coverage -a sorted.!{fileList[0]} -b !{fileList[1]} -d -sorted -g !{params.genomeInfo} > ${pre}.per-base-coverage.txt
    rm sorted.!{fileList[0]} no-negatives.!{fileList[0]}
    bases=$(cat ${pre}.per-base-coverage.txt | wc -l)
    bad=$(awk '{if ($(NF) <= 2) { print $(NF) } }' ${pre}.per-base-coverage.txt | wc -l)
    bad_pct=$(echo $bad $bases | awk '{ print $1/$2 }')
    okay=$(awk '{if ($(NF) > 2 && $(NF) <= 6) { print $(NF) } }' ${pre}.per-base-coverage.txt | wc -l)
    okay_pct=$(echo $okay $bases | awk '{ print $1/$2 }')
    good=$(awk '{if ($(NF) >= 7) { print $(NF) } }' ${pre}.per-base-coverage.txt | wc -l)
    good_pct=$(echo $good $bases | awk '{ print $1/$2 }')
    echo -e "${pre}\t${bad}(${bad_pct})\t${okay}(${okay_pct})\t${good}(${good_pct})" >> ${pre}.coverageStatistics.txt
    rm ${pre}.per-base-coverage.txt 
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
        cat $file >> coverageStatistics.txt
        rm $file
    done
    sort -k1,1 coverageStatistics.txt >> coverageStatistics_sorted.txt
    rm coverageStatistics.txt
   '''
}
