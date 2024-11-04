# Yarrowia-lipolytica-library-design
The scripts presented here can be used for designing an 6-fold coverage sgRNAs targeting the first 40% of all coding sequences and tRNA regions within _Yarrowia lipolytica_ CLIB89 genome. These scripts can be modified as needed for targeting other loci or other haploid organisms and to design libraries with various coverages. 

# Publication

## Prerequisite
python 2.7.14

scipy 1.2.2

gffutils 0.10.1

mySQL-python 1.2.3

numpy 1.16.6

pandas 0.24.2

scikit-learn 0.18.1

scipy 1.2.2

Laptop or desktop computer that meets the requirements to run python 2.7.14

## Running the code

The first step in designing the CRISPR-Cas9 sgRNA library using our codes is to be able to run CHOPCHOP v3 on the command line (https://doi.org/10.1093/nar/gkz365). All the necessary requirements and step-by-step installation guides for CHOPCHOP v3 are explained by its authors here: https://bitbucket.org/valenlab/chopchop/src/master/. After installation, you should be able to run the examples mentioned on their website properly.

To use our source code to design a genome-wide sgRNA library, clone this repository: 

<code>git clone https://github.com/ianwheeldon/Yarrowia_lipolytica_library_design.git/</code>

Running <code>execute.sh</code> would download _Yarrowia lipolytica_'s genome from NCBI website, clone CHOPCHOP repository, and make the necessary files for running CHOPCHOP.

<code>bash execute.sh</code>

<code>cd library_design/chopchop</code>

<h3>Run example on Yarrowia lipolytica CLIB89 genomic data</h3>

To make sure cCHOPCHOP is running properly, run this command as an example:

<code>./chopchop.py -T 1 -M NGG --maxMismatches 3 -g 20 -G CLIB89 -o Results -Target NC_090770.1:5100-5200 --scoringMethod ALL</code>

You should get these results: 

<code>Rank    Target sequence Genomic location        Strand  GC content (%)  Self-complementarity    MM0     MM1     MM2     MM3     XU_2015 DOENCH_2014     DOENCH_2016       MORENO_MATEOS_2015      CHARI_2015      G_20    ALKAN_2018      ZHANG_2019
1       CCACAGTGCTGCTGAGATGCCGG NC_090770.1:5101        +       60      0       0       0       0       0       0.50    0.03    57.61   0.59    0.00    0.00    24.70     57.43
2       GCAGCACTGTGGAACATGCATGG NC_090770.1:5090        -       55      0       0       0       0       0       0.40    0.50    57.83   0.59    0.00    0.00    24.39     80.07
3       AAGCCATATGATAGGGCTCTCGG NC_090770.1:5156        -       45      0       0       0       0       0       0.39    0.04    57.42   0.66    0.00    0.00    24.17     11.42
4       CTGCCGAGAGCCCTATCATATGG NC_090770.1:5153        +       55      0       0       0       0       0       0.53    0.09    38.45   0.62    0.00    0.00    23.17     30.02
5       CGTATGCAATAACGAGCTCCTGG NC_090770.1:5195        -       50      0       0       0       0       0       0.42    0.05    48.30   0.57    0.00    0.00    22.21     39.66
6       TTGATTGTGGGAATACAACCAGG NC_090770.1:5177        +       40      0       0       0       0       0       0.55    0.21    60.01   0.62    0.00    0.00    21.65     76.95
7       GTTTAAGTTTCACGGCTGACCGG NC_090770.1:5120        -       45      0       0       0       0       0       0.42    0.32    54.41   0.59    0.00    0.00    19.36     38.26
8       CACAATCAAGCCATATGATAGGG NC_090770.1:5163        -       35      0       0       0       0       0       0.52    0.06    43.13   0.58    0.00    0.00    18.51     51.12
9       CCTATCATATGGCTTGATTGTGG NC_090770.1:5164        +       40      0       0       0       0       0       0.43    0.12    55.08   0.54    0.00    1.00    17.99     47.81
10      CAGCACTGTGGAACATGCATGGG NC_090770.1:5089        -       50      0       0       0       0       1       0.59    0.40    63.49   0.69    0.00    0.00    21.82     68.88
11      CTATCATATGGCTTGATTGTGGG NC_090770.1:5165        +       35      0       0       0       0       0       0.41    0.22    43.68   0.61    0.00    0.00    16.49     25.21
12      AGAGTATAGTTTAAGTTTCACGG NC_090770.1:5128        -       25      0       0       0       0       0       0.38    0.08    35.32   0.58    0.00    0.00    11.83     18.29</code>

<h3>Designing the genome-wide CRISPR-Cas9 sgRNA library</h3>

Running the <code>Library_design.py</code> code would start designing the library and will generate three <code>.csv</code> files: 1- <code>BEST_LIBRARY.csv</code> which is the n-fold coverage library; 2- <code>CHOPCHOP_Total.csv</code> containing all of the sgRNAs CHOPCHOP found within the specified region; and 3- <code>non-targeting.csv</code> that is the list of non-targeting sgRNAs designed based on the library size. 

<code>python Library_design.py</code>

