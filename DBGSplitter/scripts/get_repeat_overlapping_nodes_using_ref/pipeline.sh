#!/bin/sh
#identifies, using STAR and the RM tracks, the nodes that are due to *genomic* repeats.


#These are not from this script - it should be given, but I am putting here just as a reference
#getting the ref genome
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

#creating star index
#/data/Public/STAR-2.5.0b/bin/Linux_x86_64/STAR --runMode genomeGenerate --runThreadN 8 --genomeDir /data2/leandro/repeats_on_transcriptome/ref_genome/STAR_index --genomeFastaFiles /data2/leandro/repeats_on_transcriptome/ref_genome/hg38.fa --sjdbGTFfile /data2/leandro/repeats_on_transcriptome/ref_annotation/hg38_gencode_v24.gtf --sjdbOverhang 100
#End



set -eux

#get pars
nodesFile=${1}
prefix=${2}
STARlong=${3}
STARGenomeDir=${4}
BEDTools=${5}
RepeatMaskerBedTrack=${6}
RepeatIntersectionThreshold=${7}
threads=${8}



#create the .unitigs file from .nodes
awk '{OFS=""; print ">",$1,"\n",$2}' ${nodesFile} > ${prefix}.unitigs.fasta



#aligning the unitigs to the ref genome using star:
#care for --outFilterMatchNminOverLread 0.9: we are requiring at least 90% of matches to map an unitig
#I saw that with this, unmapped reads are not in the bam... They were supposed to be, I think, since --outSAMunmapped Within
#But I can get them with --outReadsUnmapped Fastx. Anyway, what I wanted is to know which unmapped reads were due to way too many mapped loci, but I can't get this information with any of the two methods...
${STARlong} --runThreadN ${threads} --genomeDir ${STARGenomeDir} --readFilesIn ${prefix}.unitigs.fasta --outFileNamePrefix ${prefix}.unitigs.fasta_mapped --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outSAMunmapped Within --outFilterMatchNminOverLread 0.9

#intersecting tha bam with UCSC's repeat track
#-split : Treat “split” BAM (i.e., having an “N” CIGAR operation) or BED12 entries as distinct BED intervals.
${BEDTools} intersect -abam ${prefix}.unitigs.fasta_mappedAligned.sortedByCoord.out.bam -b ${RepeatMaskerBedTrack} -split -bed > ${prefix}.unitigs.fasta_mappedAligned.intersected_with_repeatmasker_tracks.bed


#getting a list of <unitigId, repeatSize>
awk '{print $4, $3-$2}' ${prefix}.unitigs.fasta_mappedAligned.intersected_with_repeatmasker_tracks.bed | sort -nk1,1  > ${prefix}.unitigs.fasta_mappedAligned.unitigId_repeatSize

#getting the list of unitigs that could have problematic repeats (i.e. unitigs having repeat size >= $RepeatIntersectionThreshold)
awk "{if(\$2>=${RepeatIntersectionThreshold}){print \$1}}" ${prefix}.unitigs.fasta_mappedAligned.unitigId_repeatSize | sort -n | uniq > ${prefix}.unitigs.fasta_mappedAligned.unitigs_containing_problematic_repeats

#we get the unmapped unitigs
#this is not ideal, I would like to get the unmapped unitigs due to way too many multi mapping, but as I can't do this, let's get all unmapped unitigs
#it is not so bad also:
#        Number of reads mapped to too many loci |       1170861
#             % of reads mapped to too many loci |       4.80%
#                                  UNMAPPED READS:
#       % of reads unmapped: too many mismatches |       0.00%
#                 % of reads unmapped: too short |       2.01%
#                     % of reads unmapped: other |       1.48%
#The majority is due to too many loci
grep ">" ${prefix}.unitigs.fasta_mappedUnmapped.out.mate1 | awk '{print substr($1, 2)}' | sort -n | uniq > ${prefix}.unitigs.fasta_mapped.unmapped_unitigs

#put unmapped and repeat unitigs together
cat ${prefix}.unitigs.fasta_mappedAligned.unitigs_containing_problematic_repeats ${prefix}.unitigs.fasta_mapped.unmapped_unitigs | sort -n | uniq > ${prefix}.unitigs.fasta_mapped.unitigs_containing_repeats_or_unmapped








#NOT USED/NEEDED RIGHT NOW, BUT I LEFT HERE BECAUSE IT WAS IN THE OLD SCRIPT
#NOT USED/NEEDED RIGHT NOW, BUT I LEFT HERE BECAUSE IT WAS IN THE OLD SCRIPT
#NOT USED/NEEDED RIGHT NOW, BUT I LEFT HERE BECAUSE IT WAS IN THE OLD SCRIPT
#Generating a repeat-free unitig file:
#python ../scripts/get_repeat_overlapping_nodes_using_ref/generateRepeateFreeUnitigFile.py ${prefix}.unitigs.fasta ${prefix}.unitigs.fasta_mapped.unitigs_containing_repeats_or_unmapped > ${prefix}.unitigs.repeat_free.fasta


#running KS in the repeat-free unitig file just to build the graph and quit
#kissplice-2.4.0-p1_build_graph_only/build/bin/kissplice -r graph_mcf7_k41.nodes.repeat_free.fasta -o ks_repeat_free_graph -t 4 -v -c 0 -C 0

#running KS to get the repeat free BCCs
#ulimit -s unlimited
#kissplice-2.4.0-p1/build/bin/kissplice -b 1 -g ks_repeat_free_graph/graph_graph_mcf7_k41_k41 -c 0 -C 0 -o repeat_free_bccs -t 4 -v --keep-bccs --experimental

#getting the distribution of the size of the repeat free bccs
#ls -1 repeat_free_bccs/bcc/*.nodes | xargs -I '{}' wc -l {} > size_repeat_free_bccs
#awk '{print $1}' size_repeat_free_bccs | sort -n | uniq -c | awk '{print $1, $2}' | cat <(echo "count size_of_bcc") - > distribution_count_size_bcc

#getting the AS events with b=5 to check if we lose or not AS events when filtering out the repeat nodes from the graph
#command line should be inspired by http://kissplice.prabi.fr/pipeline_ks_farline/
#ulimit -s unlimited
#kissplice-2.4.0-p1/build/bin/kissplice -r /data2/leandro/mcf7/raw_data/siCTL_N1_GGCTAC_R1_trim_right_cutadapt_match_trim_to_100_match.fastq.gz -r /data2/leandro/mcf7/raw_data/siCTL_N1_GGCTAC_R2_trim_right_cutadapt_match_trim_to_100_match.fastq.gz -r /data2/leandro/mcf7/raw_data/siCTL_N2_CTTGTA_R1_trim_right_cutadapt_match_trim_to_100_match.fastq.gz -r /data2/leandro/mcf7/raw_data/siCTL_N2_CTTGTA_R2_trim_right_cutadapt_match_trim_to_100_match.fastq.gz -r /data2/leandro/mcf7/raw_data/siDDX5_17_N1_AGTCAA_R1_trim_right_cutadapt_match_trim_to_100_match.fastq.gz -r /data2/leandro/mcf7/raw_data/siDDX5_17_N1_AGTCAA_R2_trim_right_cutadapt_match_trim_to_100_match.fastq.gz -r /data2/leandro/mcf7/raw_data/siDDX5_17_N2_AGTTCC_R1_trim_right_cutadapt_match_trim_to_100_match.fastq.gz -r /data2/leandro/mcf7/raw_data/siDDX5_17_N2_AGTTCC_R2_trim_right_cutadapt_match_trim_to_100_match.fastq.gz -g ks_repeat_free_graph/graph_graph_mcf7_k41_k41 -c 0 -C 0 -o repeat_free_results -t 4 -v --experimental --mismatches 2 --counts 2 --min_overlap 5
#NOT USED/NEEDED RIGHT NOW, BUT I LEFT HERE BECAUSE IT WAS IN THE OLD SCRIPT
#NOT USED/NEEDED RIGHT NOW, BUT I LEFT HERE BECAUSE IT WAS IN THE OLD SCRIPT
#NOT USED/NEEDED RIGHT NOW, BUT I LEFT HERE BECAUSE IT WAS IN THE OLD SCRIPT
