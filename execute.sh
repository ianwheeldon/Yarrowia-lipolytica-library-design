mkdir library_design
cd library_design
git clone https://bitbucket.org/valenlab/chopchop.git
cd chopchop
cd uCRISPR
g++ -o uCRISPR uCRISPR.cpp -std=c++11
cd ..

#Downloading the CLIB89 fasta file
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/761/485/GCF_001761485.1_ASM176148v1/GCF_001761485.1_ASM176148v1_genomic.fna.gz;
gzip -d GCF_001761485.1_ASM176148v1_genomic.fna.gz;

#Changes the chromosome names in the fasta file to match the names in the gff3 file
sed -E '/^>/ s/(>[^[:space:]]*).*/\1/' GCF_001761485.1_ASM176148v1_genomic.fna > GCF_001761485.1_ASM176148v1_genomic_fix.fna

#Changes all the characters of the fasta file to uppercase
tr '[:lower:]' '[:upper:]' < GCF_001761485.1_ASM176148v1_genomic_fix.fna > GCF_001761485.1_ASM176148v1_genomic_upper.fna;
rm GCF_001761485.1_ASM176148v1_genomic_fix.fna

#Downloading the CLIB89 gff3 file
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/761/485/GCF_001761485.1_ASM176148v1/GCF_001761485.1_ASM176148v1_genomic.gff.gz;
gzip -d GCF_001761485.1_ASM176148v1_genomic.gff.gz;

#making the required files for CHOPCHOP
mkdir CLIB89;
mv GCF_001761485.1_ASM176148v1_genomic_upper.fna CLIB89/CLIB89.fasta;
bowtie/bowtie-build CLIB89/CLIB89.fasta CLIB89/CLIB89

rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToTwoBit ./

./faToTwoBit CLIB89/CLIB89.fasta CLIB89/CLIB89.2bit

rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/gff3ToGenePred ./

./gff3ToGenePred -geneNameAttr=gene_name GCF_001761485.1_ASM176148v1_genomic.gff GCF_001761485.1_ASM176148v1_genomic.genePred
awk 'BEGIN { FS = OFS = "\t" } { gsub(/\.[^.]*$/, "", $1) }1' GCF_001761485.1_ASM176148v1_genomic.genePred > GCF_001761485.1_ASM176148v1_genomic_fixed.genePred
echo -e "name\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames" | cat - GCF_001761485.1_ASM176148v1_genomic_fixed.genePred > GCF_001761485.1_ASM176148v1_genomic_fixed_fix.genePred
mv GCF_001761485.1_ASM176148v1_genomic_fixed_fix.genePred CLIB89/CLIB89.gene_table
rm *.genePred

cat<< EOF > config_local.json
{
  "PATH": {
    "PRIMER3": "./primer3_core",
    "BOWTIE": "bowtie/bowtie",
    "TWOBITTOFA": "./twoBitToFa",
    "TWOBIT_INDEX_DIR": "CLIB89",
    "BOWTIE_INDEX_DIR": "CLIB89",
    "ISOFORMS_INDEX_DIR": "CLIB89",
    "ISOFORMS_MT_DIR": "CLIB89",
    "GENE_TABLE_INDEX_DIR": "CLIB89"
  },
  "THREADS": 1
}
EOF

module load viennarna/2.5.0
