mamba activate chipseq
cd /mnt/s/AG/AG-Scholz-NGS/Daten/Simon/P3026_ChIP-Seq_epiSVF/deeptools
mkdir bigWig

# effective genome size mm39
# 2,654,621,783

gs=2654621783
run="pe";

dirdata="/mnt/s/AG/AG-Scholz-NGS/Daten/ChIP_epiSVF_P3026";
dirs="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/P3026_ChIP-Seq_epiSVF"; 
bl="/mnt/s/AG/AG-Scholz-NGS/Daten/Simon/Genomic_data/Mus/GRCm39/mm39.excluderanges.bed.gz"
bdir="$dirs/Bowtie2";
bamout=$bdir/$run;
bamin=$bamout
name=$run"_";

dtoolsdir="$dirs/deeptools";
dtout=$dtoolsdir/bigWig/$run
mkdir $dtout

# big wig
log=$dtout"/"$name$(date +%Y_%m_%d)"_bigwig.log";

echo -e "########################################## \n\
# Date: $(date +%Y_%m_%d) $(date +%T) \n\
# Bam files: $bamin \n\
# Output folder: $dtout" | tee -a $log;
echo $bamin/pe_P3026_01_ChIP_1580_S23.bam


bamCoverage --bam $bamin/pe_P3026_01_ChIP_1580_S23.bam -o $dtout/pe_P3026_01_ChIP_1580_S23-2.bw \
    --binSize 100 \
    --smoothLength 200 \
    --normalizeUsing BPM \
    --effectiveGenomeSize $gs \
    --extendReads \
    -bl $bl \
    --centerReads \
    --ignoreForNormalization chrX \
    --filterRNAstrand forward \
    --skipNonCoveredRegions \
    --numberOfProcessors max/2
    
bamCompare -b1 $bamin/pe_P3026_01_ChIP_1580_S23.bam \
    -b2 $bamin/pe_P3026_01_ChIP_1578_S22.bam \
    -o $dtout/pe_P3026_01_ChIP_1580_S23_1578-2.bw \
    --binSize 100 \
    --smoothLength 200 \
    --scaleFactorsMethod None \
    -l 1000 \
    --normalizeUsing BPM \
    --effectiveGenomeSize $gs \
    --skipZeroOverZero \
    --centerReads \
    --extendReads \
    --skipNonCoveredRegions \
    --numberOfProcessors max/2

# Maunals
## tidyGenomeBrowser
## https://github.com/MalteThodberg/tidyGenomeBrowser/blob/master/README.Rmd

## trackViewer
## https://bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html

## Generate bigwig
## https://hbctraining.github.io/Intro-to-ChIPseq-flipped/lessons/08_creating_bigwig_files.html

## Normalize bigwig
## https://www.biostars.org/p/394434/

