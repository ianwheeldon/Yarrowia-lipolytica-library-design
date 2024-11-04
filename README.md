# Yarrowia-lipolytica-library-design
The scripts presented here can be used for designing an n-fold coverage sgRNA library using CRISPR-Cas9 system for _Yarrowia lipolytica_, capable of targeting specific regions within each feature of the genome (i.e., gene, coding sequence, exon, etc.). These scripts can also be used to design sgRNA libraries for any organism of interest with a haploid genome.
# Publication

## Prerequsite
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

Run this command:

<code>./chopchop.py -T 1 -M NGG --maxMismatches 3 -g 20 -G GS115 -o Results -Target CP014715.1:5164-5200 --scoringMethod ALL</code>

