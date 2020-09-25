#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess

"""
All parameter are desbribe here :

    -sr            (1 arg) :    A text file containing the input reads. Each input read is described in a line containing 3 columns: 1) path to a fasta, fastq, or .gz file containing the reads; 2) condition of the input reads (int); 3) replicate number of the input reads (int). This file needs a header. Files can have the same condition and replicate (in case of paired-end reads, for example). Check https://gitlab.inria.fr/lishisoa/DBGSplitter/blob/master/tests/test_map2k2/MAP2K2 for an example.
       -k             (1 arg) :    K-mer size  [default '41']
       -rel-cutoff    (1 arg) :    A threshold. Edges which counts are relatively smaller than this threshold are removed. This is applied only on the short reads graph.  [default '0.02'] = -C
       -min-abundance (1 arg) :    An integer, k-mers present strictly less than thisnumber of times in the dataset will be discarded. This is applied only on the short reads graph.  [default '2'] = -c
       -beta-max      (1 arg) :    Maximum beta layer  [default '100']
       -rep-int-th    (1 arg) :    The minimum size of an intersection of a node and a repeat to consider the node as a repeat  [default '20']
       -rmtrack       (1 arg) :    Path to Repeat Masker bed track  [default '/data2/leandro/repeats_on_transcriptome/ref_repeats_track/hg38_repeatmasker_tracks.bed']
       -bedtools      (1 arg) :    Path to bedtools  [default '/data2/leandro/repeats_on_transcriptome/bedtools/bedtools2/bin/bedtools']
       -stargenome    (1 arg) :    Path to STAR genome dir  [default '/data2/leandro/repeats_on_transcriptome/ref_genome/STAR_index']
       -starlong      (1 arg) :    Path to STARlong binary  [default '/data2/leandro/repeats_on_transcriptome/STAR_2.5.3a/STAR-2.5.3a/bin/Linux_x86_64_static/STARlong']
       -rep-id-mode   (1 arg) :    Repeat identification mode (reference or de-novo)  [default 'reference']
       -prefix        (1 arg) :    Prefix of the name of the built files related to the graph  [default 'graph']
       -simplify      (0 arg) :    Simplify the graph using GATB's error correction procedures (tested in genomics, not in transcriptomics). See https://github.com/GATB/gatb-core/blob/v1.4.1/gatb-core/src/gatb/debruijn/impl/Simplifications.hpp.
       -skip1         (0 arg) :    Skips step 1: building the short reads graph. Uses pre-built short reads graph.
       -skip2         (0 arg) :    Skips steps 1 and 2: identifying nodes overlapping repeats. Uses pre-computed identification.
       -skip3         (0 arg) :    Skips steps 1, 2 and 3: building gene components. Uses pre-computed gene components.
       -nb-cores      (1 arg) :    number of cores  [default '0']
       -verbose       (1 arg) :    verbosity level  [default '1']
       -version       (0 arg) :    version
       -help          (0 arg) :    help
"""



if __name__ == "__main__":

    #PARSE ARGUMENTS :

    parser = argparse.ArgumentParser(description='To run our complex splicing event dectection and differential analysis pipeline in RNAseq data')
    parser.add_argument("--sr", action="store", dest="sr", type=str, help="A text file containing the input reads. Each input read is described in a line containing 3 columns: 1) path to a fasta, fastq, or .gz file containing the reads; 2) condition of the input reads (int); 3) replicate number of the input reads (int). This file needs a header. Files can have the same condition and replicate (in case of paired-end reads, for example). Check https://gitlab.inria.fr/lishisoa/DBGSplitter/blob/master/tests/test_map2k2/MAP2K2 for an example.", required = True)
    #new output folder:
    parser.add_argument("--output_dir", action="store", dest="output_dir", type=str, help="output directory (without the / at the end)", required = False, default=".")

    parser.add_argument("--k", action="store", dest="k", type=int, help="K-mer size  [default '41']", required = False, default=41)
    parser.add_argument("--C", action="store", dest="C", type=float, help="A threshold. Edges which counts are relatively smaller than this threshold are removed. This is applied only on the short reads graph.  [default '0.02'] = -C", required = False, default=0.02)
    parser.add_argument("--c ", action="store", dest="c", type=float, help="An integer, k-mers present strictly less than thisnumber of times in the dataset will be discarded. This is applied only on the short reads graph.  [default '2']", required = False, default=2)
    parser.add_argument("--beta_max", action="store", dest="beta_max", type=int, help="Maximum beta layer  [default '100']", required = False, default=100)
    parser.add_argument("--rep_int_th", action="store", dest="rep_int_th", type=int, help="The minimum size of an intersection of a node and a repeat to consider the node as a repeat  [default '20']", required = False, default=20)
    parser.add_argument("--bedtools", action="store", dest="bedtools", type=str, help="Path to bedtools  [default '/data2/leandro/repeats_on_transcriptome/bedtools/bedtools2/bin/bedtools']", required = False, default="/data2/leandro/repeats_on_transcriptome/bedtools/bedtools2/bin/bedtools")
    parser.add_argument("--stargenome", action="store", dest="stargenome", type=str, help="Path to STAR index genome dir  [default '/data2/leandro/repeats_on_transcriptome/ref_genome/STAR_index']", required = False, default="/data2/leandro/repeats_on_transcriptome/ref_genome/STAR_index")
    parser.add_argument("--starlong", action="store", dest="starlong", type=str, help="Path to STARlong binary  [default '/data2/leandro/repeats_on_transcriptome/STAR_2.5.3a/STAR-2.5.3a/bin/Linux_x86_64_static/STARlong']", required = False, default="/data2/leandro/repeats_on_transcriptome/STAR_2.5.3a/STAR-2.5.3a/bin/Linux_x86_64_static/STARlong")
    parser.add_argument("--rep_id_mode", action="store", dest="rep_id_mode", type=str, help="Repeat identification mode (reference or de-novo)  [default 'reference']", required = False)
    parser.add_argument("--prefix", action="store", dest="prefix", type=str, help="Prefix of the name of the built files related to the graph  [default 'graph']", required = False, default="graph")
    #parser.add_argument("--simplify", action="store_true", dest="simplify", help="Simplify the graph using GATB's error correction procedures (tested in genomics, not in transcriptomics). See https://github.com/GATB/gatb-core/blob/v1.4.1/gatb-core/src/gatb/debruijn/impl/Simplifications.hpp.", required = False)
    parser.add_argument("--skip1", action="store_true", dest="skip1", help="Skips step 1: building the short reads graph. Uses pre-built short reads graph.", required = False)
    parser.add_argument("--skip2", action="store_true", dest="skip2", help="Skips steps 1 and 2: identifying nodes overlapping repeats. Uses pre-computed identification.", required = False)
    parser.add_argument("--skip3", action="store_true", dest="skip3", help="Skips steps 1, 2 and 3: building gene components. Uses pre-computed gene components.", required = False)
    parser.add_argument("--nb_cores", action="store", dest="nb_cores", type=int, help="number of cores, [default '1']", required = False, default=1)
    parser.add_argument("--verbose", action="store_true", dest="verbose", help="verbosity level  [default '1']", required = False)
    parser.add_argument("--version", action="store_true", dest="version", help="print current version", required = False)
    #TODO :  add output folder


    args = parser.parse_args()


    #Create output directory :
    if not os.path.exists(args.output_dir) :
        commmand=["mkdir", str(args.output_dir)]
        subprocess.run(commmand)


    # PATH TO executable :
    scriptDir = os.path.dirname(os.path.realpath(__file__))
    step1="%s/step1"%scriptDir
    step2="%s/pipeline.sh"%scriptDir
    step3="%s/step3"%scriptDir
    step4="%s/run_diff_ana.R"%scriptDir
    repeat_masker_track="/data2/leandro/repeats_on_transcriptome/ref_repeats_track/hg38_repeatmasker_tracks.bed"



    #Create output directory :
    if not os.path.exists(args.output_dir) :
        commmand=["mkdir", str(args.output_dir)]
        subprocess.run(command)


    #RUN STEP 1 : Build DBG and cumpute some metrics
    folder1=args.output_dir+"/output_step1/"
    if not args.skip1 :
        if not os.path.exists(folder1):
            #create output folder:
            command_folder1=["mkdir", str(folder1)]
            subprocess.run(command_folder1)

        #run
        command_line_1=[step1, "-sr", args.sr, "-k", str(args.k), "-rel-cutoff", str(args.C), "-min-abundance", str(args.c), "-step1Folder", folder1, "-prefix", str(args.prefix), "-nb-cores", str(args.nb_cores),\
                                "-beta-max", str(args.beta_max)]
        print("Running step 1 ... ")
        print(command_line_1)
        subprocess.run(command_line_1, check=True)

    #RUN STEP 2 : split the DBG to obtain one component for each gene/gene family
    """
    Exemple for MAP2k2 :
    /data2/camille/complex_event/test_pipeline_complet/DBGSplitter/build/test_package/DBGSplitter-0.0.5-Linux/bin/pipeline.sh graph /data2/leandro/repeats_on_transcriptome/STAR_2.5.3a/STAR-2.5.3a/bin/Linux_x86_64_static/STARlong /data2/leandro/repeats_on_transcriptome/ref_genome/STAR_index /data2/leandro/repeats_on_transcriptome/bedtools/bedtools2/bin/bedtools /data2/leandro/repeats_on_transcriptome/ref_repeats_track/hg38_repeatmasker_tracks.bed

    => pipeline.sh --prefix --path_star_long --index_STAR --path_bedtools --path_track-repeat_masker.bed
    """

    folder2=args.output_dir+"/output_step2/"
    if not args.skip2 :
        if not os.path.exists(folder2):
            #create output folder:
            command_folder2=["mkdir", str(folder2)]
            subprocess.run(command_folder2)

        #run
        command_line_2=[step2, "%s%s.nodes"%(folder1, args.prefix), "%s%s"%(folder2, str(args.prefix)), str(args.starlong), str(args.stargenome), str(args.bedtools), repeat_masker_track, str(args.rep_int_th), str(args.nb_cores), str(folder2)]
        print("Running step 2 ... ")
        print(command_line_2)
        subprocess.run(command_line_2, check=True)


    #RUN STEP 3 : output quantify graph component (or bubbles)
    folder3=args.output_dir+"/output_step3/"
    if not args.skip3 :
        if not os.path.exists(folder3):
            #create output folder:
            command_folder3=["mkdir", str(folder3)]
            subprocess.run(command_folder3)
        #run
        command_line_3=[step3, "-step1Folder", folder1, "-step2Folder", folder2, "-step3Folder", folder3, "-k", str(args.k), "-nb-cores", str(args.nb_cores)]
        print("Running step 3 ... ")
        print(command_line_3)
        subprocess.run(command_line_3, check=True)


    # RUN STEP 4 : differential analysis
    # add output file name as an argument (Rscript + here)
    # run the diff analysis on all the inputR file. they are listed in the diff_analysis_file
    # the number of replicate is store in graph.nb.replicat : this is not
    #the better way to do it, I should modify my Rscript to cumpute it automatically
    #and not have to provide it as an argument.

    """
    exemple for MAP2K2 :
    Rscript /data2/camille/complex_event/test_pipeline_complet/DBGSplitter/build/
    test_package/DBGSplitter-0.0.5-Linux/bin/run_diff_ana.R graph.repeat_free.inputR 2
    """
    print("Running step 4...")

    folder4=args.output_dir+"/output_step4/"
    if not os.path.exists(folder4):
        # create output folder
        command_folder4=["mkdir", str(folder4)]
        subprocess.run(command_folder4, check=True)

    #needed files to open :
    inputR_files_list=open(folder3 + str(args.prefix)+".diff_analysis_files", "r")
    nb_rep = open(folder1 + str(args.prefix)+".nb_replicates","r")
    nb_replicat = nb_rep.readline()

    print("nb_replicat="+str(nb_replicat))

    for inputR_file in inputR_files_list :
        #print(inputR_file[:-1])
        inputR_file_full = str(folder3)+ str(inputR_file.split("/")[-1])
        output_file_name = folder4 + str(args.prefix) + "." + inputR_file.split("/")[-1].split(".")[1]+ "_output_diff_analysis.tab"
        #print(output_file_name)
        Command_line_4=["Rscript", step4, inputR_file_full[:-1], nb_replicat, str(output_file_name)]
        subprocess.run(Command_line_4)
