# About DBGSplitter

We try to split RNA-seq DBG into gene components.

Just works on Linux for the moment.

This is a software tool relying on GATB-CORE library (this forked version: https://github.com/leoisl/gatb-core).

# Downloading the precompiled binaries

v0.0.5: TODO

## Previous versions:

v0.0.4: https://www.dropbox.com/s/26hfrf3oo1vn7x0/DBGSplitter-0.0.4-Linux.tar.gz?dl=1

v0.0.3: https://www.dropbox.com/s/2seyxsuhpikxgjn/DBGSplitter-0.0.3-Linux.tar.gz?dl=1

v0.0.2: https://www.dropbox.com/s/duyz3q89uuvtyxx/DBGSplitter-0.0.2-Linux.tar.gz?dl=1

v0.0.1: https://www.dropbox.com/s/y0p0v7tz6xcl1d1/DBGSplitter-0.0.1-Linux.tar.gz?dl=1

# Dependencies

The following third parties should be already installed:

* cmake 3.1+ (mandatory)
* Boost

# Project build

For building your project, you should do the following (tested only on Linux)
```
    git clone --recursive https://gitlab.inria.fr/lishisoa/DBGSplitter
    cd DBGSplitter && mkdir build && cd build && cmake -DABUNDANCE_TYPE=u_int32_t -DSKIP_DISCRETIZATION=1 .. &&  make -j8 && make package
```
Then, you should get a package `DBGSplitter-<version>-Linux.tar.gz`.

# Running on a sample dataset
To test on reads mapping to a sample human gene (MAP2K2), with 2 conditions and 2 replicates, do:
```
    cd bin; cp ../tests/test_map2k2/* . ; ./step1 -sr MAP2K2 ; ./step2 ; ./step3 ; ./step4
```

To test on reads mapping to another sample human gene (BCS1L), with 2 conditions and 8 replicates, do:
```
    cd bin; cp ../tests/test_BCS1L/* . ; ./step1 -sr BCS1L ; ./step2 ; ./step3 ; ./step4
```

# If running from `build/tools`:

To test on reads mapping to a sample human gene (MAP2K2), with 2 conditions and 2 replicates, do:
```
cd build/tools; cp ../../tests/test_map2k2/* .; ./step1 -sr MAP2K2 ; ./step2 ; ./step3 ; ./step4
```

To test on reads mapping to another sample human gene (BCS1L), with 2 conditions and 8 replicates, do:
```
cd build/tools; cp ../../tests/test_BCS1L/* .; ./step1 -sr BCS1L ; ./step2 ; ./step3 ; ./step4
```

# Parameters
```
       -sr            (1 arg) :    A text file containing the input reads. Each input read is described in a line containing 3 columns: 1) path to a fasta, fastq, or .gz file containing the reads; 2) condition of the input reads (int); 3) replicate number of the input reads (int). This file needs a header. Files can have the same condition and replicate (in case of paired-end reads, for example). Check https://gitlab.inria.fr/lishisoa/DBGSplitter/blob/master/tests/test_map2k2/MAP2K2 for an example.
       -k             (1 arg) :    K-mer size  [default '41']
       -rel-cutoff    (1 arg) :    A threshold. Edges which counts are relatively smaller than this threshold are removed. This is applied only on the short reads graph.  [default '0.02']
       -min-abundance (1 arg) :    An integer, k-mers present strictly less than thisnumber of times in the dataset will be discarded. This is applied only on the short reads graph.  [default '2']
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
```