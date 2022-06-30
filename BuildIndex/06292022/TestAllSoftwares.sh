sudo apt install openjdk-11-jdk
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh   ##install rustc
rustc  --version
sudo apt install samblaster
pip3 install    pyopenssl ndg-httpsclient   cython         six   pysam  numpy scipy matplotlib       nose  virtualenv   theano       multiqc   deeptools    macs3   htseq                                         
pip3 install   --upgrade   cython   pytest progressbar2   tqdm  khmer   six  pysam  numpy   scipy  matplotlib     nose  virtualenv   theano    multiqc   deeptools    macs3   htseq
pip3 install   pysam   coverage eta swalign  pyBigWig  bio_utils 
pip3 install   setuptools  pyparsing cython numpy PyYAML alignlib-lite biopython drmaa  hgapi  matplotlib-venn   matplotlib  networkx   openpyxl   pandas pysam  rdflib 
pip3 install    ruffus  scipy  bx-python  sphinx  sphinxcontrib-programoutput  sqlalchemy  threadpool  web.py  weblogo  xlwt  pybedtools  pep8  CGATReport  cgat
pip3 install   --upgrade   setuptools  pyparsing cython numpy  PyYAML alignlib-lite biopython    hgapi   networkx   openpyxl   pandas pysam  rdflib 
pip3 install   --upgrade   ruffus  scipy  bx-python  sphinx  sphinxcontrib-programoutput  sqlalchemy  threadpool  web.py  weblogo  xlwt  pybedtools  pep8  CGATReport  cgat
pip3 install   CrossMap    jellyfish  illuminate   mictools   fastqp   fastools  MultiQC  editdistance  pysam  coverage  eta   swalign   Atropos  cutadapt fastqp  fastools   kPAL 
pip3 install   PePr   toolshed   RSeQC  bx-python  bx pybedtools  Pandas  Seaborn intervene rpy2   deeptools multiqc htseq

sudo apt install  libbam-dev  libhts-dev  libhts3  libseqlib-dev  libseqlib2  tabix  liblzma-dev  libcurl4  libcrypto++-dev  libcrypto++-utils libghc-zlib-dev
sudo apt install  libsparsehash-dev  libsdsl-dev  libsdsl3
sudo apt install genometools
conda    install -c bioconda fqtools  
conda    install -c bioconda "hmmer>=3.0"
conda    install -c bioconda "pychopper>=2.0"
conda    install -c bioconda falco  mrsfast
cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)
cpanm    Getopt::Long   Pod::Usage   File::Temp   Fcntl   Digest::MD5   Cwd   List::Util  JSON  Cairo  Statistics::PCA  MIME::Base64
cpanm    CGI   File::Path   IO::Uncompress::AnyUncompress   LWP::Simple   File::Copy   File::Basename  Statistics::R  Text::Table Text::Levenshtein::XS
if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager") }   
BiocManager::install( c("ngsReports", "seqTools", "Rqc", "systemPipeR", "qrqc", "ShortRead" ) )
BiocManager::install( "alpine" )
pip3 install --upgrade  kb-python  cutadapt  multiqc  fastools  fastqp  umi_tools  pyfastaq    fast2q  ## pycoQC squidpy
pip3 install -U loompy
pip3 install pyfastx  biopython  fastq  miniFasta   fastq-statistic

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install ngs-bits
conda install -c conda-forge -c bioconda busco=5.3.2  superstr
 
cpanm  CGI   File::Path   IO::Uncompress::AnyUncompress   LWP::Simple   File::Copy   File::Basename
sudo apt  install   libmaus2-2  libmaus2-dev
##  python3 setup.py install  --prefix /media/yp/1one/MyProgramFiles/3_BAM_QC/SAMstats
pip3 install midr  RSeQC  rnaseqc
if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager") }   
BiocManager::install( c("BatchQC", "QuasR",   "scater", "DESeq2", "baySeq"  ) )
library(devtools)
install_github("slzhao/MultiRankSeq")
install_github("liuqivandy/scRNABatchQC")
install.packages("spp", dependencies=TRUE) 
devtools::install_github("skyhorsetomoon/Trumpet")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("exomePeak")  ## sudo R CMD INSTALL  exomePeak_2.16.0.tar.gz
install.packages("mashr")
conda install -c bioconda checkqc  mirtrace  mirge3 pheniqs  alevin-fry
pip install interop  jcvi  mirtop  slamdunk
conda install -c conda-forge -c bioconda nanoq=0.9.0
## PacBio Secondary Analysis Tools on Bioconda: https://github.com/PacificBiosciences/pbbioconda 
conda install -c bioconda bam2fastx  extracthifi isoseq3 lima minorseq pbaa pbccs pbbam pbcommand pbcopper  pbcore pbcoretools  pbgcpp pbipa  pblaa   pbmarkdup  pbmm2  pbpigeon pbskera pbsv   recalladapters zmwfilter      
conda install -c bioconda pbccs  pbmarkdup  ## https://ccs.how/
pip install git+https://github.com/yodeng/fsplit.git
pip install barcode_splitter  pod5_format  pod5_format_tools  ont-bonito  ont-remora






  

#################################################################################################   
## 1_Pre-align  33+27
AdapterRemoval   --version
barcode_split_trim.pl  --version
biobloommaker    --version
bustools | grep bustools
csvtk  version
defastq  2>&1  | grep version
falco  --version
FaQCs  --version
fasten_regex  --version
fastp   --version
fastq-multx  2>&1 | grep Version
fastq-scan   -v
fastq_screen    --version
fastq-uniq      --version   ## fastq-tools
fastq_pre_barcodes   ## fastq_utils
fastqc  --version
FASTQuick    2>&1  |  grep Version
fastv  --version
fastx_clipper  -h   | grep "0.0" ## fastx_toolkit
kmc   | grep "ver. "
kmer-db | grep "version"
MinIONQC.R  -h
mutscan   --help
prinseq-lite.pl  --version
prinseq++   --version
qc3.pl  -h  2>&1  | grep Version
seqkit   version
seqtk
seqyclean  -h  | grep "Version"
fastq-dump    --version     ## sratoolkit
trim_galore   --version
trimmomatic.jar  -version  ## java -jar /media/yp/1one/MyProgramFiles/1_Pre-align/Trimmomatic/trimmomatic.jar   -version
unikmer  version

library( "ngsReports" )
library( "seqTools" )
library( "Rqc" )
library( "systemPipeR" )
library( "qrqc" )
library( "ShortRead" )
kb -h
cutadapt --version
multiqc  --version
fastools -v
fastqp -h
fqtools -version
umi_tools  --version
fastaq   version   ##pyfastaq
pycoQC  --version
cdna_classifier.py  -h  ## pychopper
python -m fast2q -h ##
pyfastx  --version
python -m interop --test  ##import fastq as fq  ## python3
fastq-stat --help
checkqc --version
pheniqs  --version
fsplit  -h
barcode_splitter   --version
pod5-inspect  -h
bonito basecaller
remora  infer  -h






#################################################################################################   
## git clone --recursive  https://github.com/smithlabcode/abismal.git
## n a word, if your intent is to clone-only a repo, 
##     use HTTPS URL (https://github.com/{user_name}/{project_name}.git) instead of SSH URL 
##      (git@github.com:{user_name}/{project_name}.git), 
##      which avoids (unnecessary) public key validation.
## zip -r xxx.zip file1 file2 .... , tar -xvzf xxx.tar.gz
## 2_Aligners  30+5
abismal  -about  ## two-letter alphabet
Arioc
bbmap.sh  --version
biscuit  version
bismark  --version
bowtie   --version
bowtie2  --version
bwa
bwa-mem2.avx512bw   version  ##  bwa-mem2.avx2  mem
gem-mapper   --version
gmap --version  ## gsnap  --version
graphmap
hicup   --version
hisat2  --version   ## hisat-3n  --version  , hisat2_extract_snps_haplotypes_VCF.py   -h
kallisto  version
lastal  --version
minigraph  --version
minimap2  --version
test-mwf  ## miniwfa
rsem-calculate-expression  --version
salmon   --version  ## conda install -c bioconda salmon,   which generateDecoyTranscriptome.sh
## seqan3
snap-aligner  2>&1  | grep version
SNPsplit  --version
STAR  --version
subread-align  -v ## featureCounts  -v , subread-fullscan
## UNCALLED
snp2h5  ## WASP
align_benchmark  -h  ##WFA
winnowmap  --version

library( "alpine" )
mrsfast      --version
alevin-fry   --version  ##  framework for single-cell analysis
mirge3  #  miRge3.0 -h   ## microRNA alignment software for small RNA-seq data
superstr    --help   ## A lightweight, alignment-free utility for detecting repeat-containing reads in short-read WGS, WES and RNA-seq data.

conda create -n myenv1 -c bioconda -c conda-forge superstr  
conda activate myenv1
conda deactivate  myenv1
conda create -n myenv2
conda activate myenv2 
which python 
python -m pip install --user cutadapt reportlab==3.5.42 biopython==1.78  scikit-learn==0.23.1  hypothesis==5.15.1 pytest==5.4.2  scipy==1.4.1  matplotlib==3.2.1  joblib==0.15.1  pandas==1.0.3 future==0.18.2
python -m pip install --user mirge3








#################################################################################################
##  3_BAM_QC  26+17
run_bam_qc.py  --version
bamstats   --version
bamtools   --version
bam help |  grep Version  ## bamUtil
bcftools  -v
bamcollate2  --version  ## biobambam2
busco  -v  ## which busco 
which chase.jar
which fgbio.jar
goleft -v
tabix 2>&1  | grep Version
# iSeqQC_shinyapp
mosdepth   --version
which cal_mrin.R
which run_spp.R
which picard.jar
preseq
which QoRTs.jar
qualimap   -v
quast.py  --version
RQC-parallel-qc  -h
rnaseqc --version
sambamba   --version
samstat   --version
which SAMstats
samtools version

IDR ## (IDRfilter in R package ChIPpeakAnno)
library( BatchQC )
library( QuasR )
library( scater )
library( spp )
library( idr2d  )
library( MultiRankSeq  )
library( scRNABatchQC  )
library(Trumpet)
mIDR    ## python3
RSeQC   ## python3
RNA-SeQC   ## python3
jcvi 
mirtop
mirtrace  -v 
slamdunk  --version 
deeptools  --version 









#################################################################################################
##  4_Post-align   8+3
bedops   --version
bedtk version
bedtools --version
cgmaptools  | grep Version
gt -version  # genometools
gfatools  version
homer
bedCoverage  # UCSC_Utilities

macs2 --version
macs3 --version
library(MAnorm2)







#################################################################################################
## 5_For_Figures (1+6)
methstates  ## methpipe

library(ComplexHeatmap)
intervene -v
DANPOS3
library(methylKit)
library(DSS)
library(RnBeads)







#################################################################################################
## git clone --recursive    https://github.com/broadinstitute/gatk.git
## https://github.com/seandavi/awesome-single-cell
## http://bioconductor.org/books/3.15/OSCA/
## 6_XWAS_xQTY (6+5)
echtvar  --version
gatk --version
octopus  --version
which VarScan.jar
vcftools  --version
vg  version

QTLtools 
DeepVariant
whatshap --version  ##WhatsHap is a software for phasing genomic variants using DNA sequencing reads, also called read-based phasing or haplotype assembly. 
library(vcfR)
library(mashr)



