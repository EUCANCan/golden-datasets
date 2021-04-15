#!/bin/bash
set -e
#EUCANCAN SNV & INDEL vcf handling
KEEP=false

while getopts "t:s:i:v:f:o:n:kh" option; do
    case "${option}" in
        t) truth=${OPTARG};;
        s) snv=${OPTARG};;
        i) indel=${OPTARG};;
        v) sv=${OPTARG};;
        u) truth_sv=${OPTARG};;
        f) FASTA=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
        n) SAMPLE_NAME=${OPTARG};;
        k) KEEP=true;;
        h)  echo "Usage:"
            echo "bash tmb_dragon.sh -t truth_file.vcf"
            echo "                   -s snv.vcf"
            echo "                   -i indel.vcf"
            echo "                   -v sv.vcf"
            echo "                   -u truth_sv.vcf"
            echo "                   -f ref_fasta.fa"
            echo "                   -o /OUTPUT_DIR/PATH"
            echo "                   -n Sample name"
            echo "                   -k (to keep intermediates files)"
            exit
    esac
done

echo " "
echo -e "[General Information]:\n"
echo "Truth File:" $truth
echo "SNV vcf file:" $snv
echo "INDEL vcf file:" $indel
echo "SV vcf file:" $sv
echo "TRUTH SV file:" $truth_sv
echo "Reference fasta file:" $FASTA
echo "output path:" $OUTPUT_DIR
echo "sample Name:" $SAMPLE_NAME
echo "keep intermediate files ?:" $KEEP
echo " "

# Create output dir:

mkdir -p $OUTPUT_DIR/$SAMPLE_NAME
OUTPUT_DIR=$OUTPUT_DIR/$SAMPLE_NAME

# Load conda env:
#conda env create -n eucancan -f golden-datasets/scripts/environment_snv.yml

source activate eucancan

# Check if vcf is mono sample:

if [ `bcftools query -l $snv |wc -l` -gt 1 ] || [ `bcftools query -l $truth |wc -l` -gt 1 ]; then
    echo $snv "is a multisample"
    echo "the vcf must be single sample "
    exit
fi
# bcftools view -c1 -Oz -s COLO829T -o curie_colo829_indels.sample.vcf curie_colo829_indels.vcf


# if [ `$BCF_DIR/bcftools query -l $file |wc -l` -gt 1 ] ; then
#     echo $ID "is a multisample"
#     echo `$BCF_DIR/bcftools query -l $file`
#     for sample_vcf in `$BCF_DIR/bcftools query -l $file`; do
#         echo $sample_vcf
#
#         $BCF_DIR/bcftools view -c1 -Oz -s $sample_vcf -o $OUTPUT_DIR/$dataset/$sample_vcf"_multisample.hg19_multianno.vcf.gz" $file
#         gunzip $OUTPUT_DIR/$dataset/$sample_vcf"_multisample.hg19_multianno.vcf.gz"


# Sorting vcf files
echo -e "[Running Information]: sorting vcf files\n"
bcftools sort $snv -o $OUTPUT_DIR/`basename $snv .vcf.gz`".sort.vcf.gz" -O z
bcftools sort $indel -o $OUTPUT_DIR/`basename $indel .vcf.gz`".sort.vcf.gz" -O z
bcftools sort $truth -o $OUTPUT_DIR/`basename $truth .vcf.gz`".sort.vcf.gz" -O z

snv=`basename $snv .vcf.gz`".sort.vcf.gz"
indel=`basename $indel .vcf.gz`".sort.vcf.gz"
truth=`basename $truth .vcf.gz`".sort.vcf.gz"

# indexing sorted vcf files
echo -e "[Running Information]: indexing vcf files\n"
bcftools index -f -o $OUTPUT_DIR/$snv".csi" $OUTPUT_DIR/$snv
bcftools index -f -o $OUTPUT_DIR/$indel".csi" $OUTPUT_DIR/$indel
bcftools index -f -o $OUTPUT_DIR/$truth".csi" $OUTPUT_DIR/$truth

# Merging SNV.vcf & INDEL.vcf:
echo -e "[Running Information]: concatenate vcf files\n"
bcftools concat -a $OUTPUT_DIR/$snv $OUTPUT_DIR/$indel -O z -o $OUTPUT_DIR/$SAMPLE_NAME"_merge.vcf.gz"

# Preparing for normalization
echo -e "[Running Information]: preparing for normalizing\n"
echo -e "[Running Information]: indexing...]"

export HGREF=$FASTA
#export HGREF=/data/annotations/pipelines/Human/hg19_base/genome/hg19_base.fa

bcftools index -f -o $OUTPUT_DIR/$SAMPLE_NAME"_merge.vcf.gz.csi" $OUTPUT_DIR/$SAMPLE_NAME"_merge.vcf.gz"
#bcftools index -f $truth

echo -e "[Running Information]: multimerge...]"

#echo $truth
#echo $OUTPUT_DIR/`basename $truth .vcf.gz`".prep.vcf.gz"
multimerge $OUTPUT_DIR/$SAMPLE_NAME"_merge.vcf.gz" -r $HGREF -f true -o $OUTPUT_DIR/$SAMPLE_NAME"_merge.prep.vcf" --process-full=1
multimerge $OUTPUT_DIR/$truth -r $HGREF -o $OUTPUT_DIR/`basename $truth .vcf.gz`".prep.vcf.gz" --process-full=1

# Normalizing vcfs:
echo -e "[Running Information]: normalizing test file\n"
merge_file=$OUTPUT_DIR/$SAMPLE_NAME"_merge.prep.vcf"

pre.py $merge_file $OUTPUT_DIR/`basename $merge_file _merge.prep.vcf`".norm.vcf.gz" -L --decompose --somatic --pass-only

echo -e "\n[Running Information]: normalizing truth file\n"

truth=$OUTPUT_DIR/`basename $truth .vcf.gz`".prep.vcf.gz"
pre.py $truth $OUTPUT_DIR/`basename $truth .vcf.gz`".norm.vcf.gz" -L --decompose --somatic --pass-only

# Running ingestion script:
echo -e "[Running Information]: Running ingestion_snv.py script \n"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

python $DIR/ingest_snv.py -samplename "SAMPLE" -o $OUTPUT_DIR/`basename $merge_file .vcf` $OUTPUT_DIR/`basename $merge_file _merge.prep.vcf`".norm.vcf.gz"

python $DIR/ingest_snv.py -samplename "SAMPLE" -o $OUTPUT_DIR/`basename $truth .vcf.gz` $OUTPUT_DIR/`basename $truth .vcf.gz`".norm.vcf.gz"

# Running som.py:
echo -e "[Running Information]: Running som.py evaluation script\n"

som.py $OUTPUT_DIR/`basename $truth .vcf.gz`".filtered.vcf" $OUTPUT_DIR/`basename $merge_file .vcf`".filtered.vcf" -o $OUTPUT_DIR/`basename $merge_file .vcf`"eval.txt" --verbose -N

echo -e "[Running Information]: script ended successfully\n"

# Running SV ingestion script:

echo -e "[Running Information]: Running ingestion SV script\n"

# Cleaning
echo -e "[Running Information]: Cleaning files\n"
rm $OUTPUT_DIR/*{tbi,csi,json,vcf,vcf.gz}
rm $truth".csi"

#run ${repo_path} "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/vcfs/Hartwig/5.19/COLO829v003T.purple.sv.vcf.gz" "COLO829v003T" "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/dataframes/hartwig.csv"

# Get the samplename:
#bcftools query -l test_merge.norm.vcf.gz

# example command:

#bash ../golden-datasets-fork-tom/golden-datasets/scripts/ingest_snv.sh -t test_truth.vcf -s results_dream/CURIE/insilico_1_snvs.vcf -i results_dream/CURIE/insilico_1_indels.vcf -o . -c test_ingestions_sh/config.txt  -n test -f /data/annotations/pipelines/Human/hg19_base/genome/hg19_base.fa

# test with curie colo files:
#bash ../golden-datasets-fork-tom/golden-datasets/scripts/ingest_snv.sh -t ../Benchmark/test_truth.vcf -s curie_colo829_snps.sample.vcf -i curie_colo829_indels.sample.vcf -o . -c ../Benchmark/test_ingestions_sh/config.txt  -n COLO829T -f /data/annotations/pipelines/Human/hg19_base/genome/hg19_base.fa

# bcftools norm -m -any ./COLO829T_merge.vcf > test.multi.vcf
# pre.py test.multi.vcf test.multi.norm.vcf -L --decompose --somatic --pass-only
# pre.py COLO829T_merge.vcf test_output.vcf -L --decompose --somatic --pass-only
#
# hap.py test.multi.vcf test.multi.norm.vcf -L --decompose --somatic --pass-only
# hap.py COLO829T_merge.vcf ../Benchmark/test_truth.vcf -o test.multi.norm.vcf -L --decompose --somatic --pass-only --unhappy
#
# /data/users/tgutman/.conda/envs/test_benchmark/bin/vcfcheck COLO829T_merge.vcf --check-bcf-errors 1
# /data/users/tgutman/.conda/envs/test_benchmark/bin/vcfcheck test.multi.vcf --check-bcf-errors 1

#handle header problems with truth files:
# Add the contig info by hand with vim - with chr (example: chr1,chr2)
# awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' truth.snvs.synthetic.challenge.set2.vcf > test_truth.vcf

#handle mulisample files:

# bcftools view -c1 -Oz -s COLO829T -o curie_colo829_indels.sample.vcf curie_colo829_indels.vcf
# bcftools view -c1 -Oz -s COLO829T -o curie_colo829_snps.sample.vcf curie_colo829_snps.vcf
#
# gunzip $OUTPUT_DIR/$dataset/$sample_vcf"_multisample.hg19_multianno.rec.vcf.gz"
#
# ${GATK_DIR}/gatk --java-options "-Xmx20g" SelectVariants \
#     -O test.multi.vcf \
#     -V COLO829T_merge.vcf \
#     --select-type-to-include SNP --select-type-to-include INDEL --exclude-filtered --select-type-to-exclude MNP

# Running Dream data challenge:

#SCRIPT_DIR=/bioinfo/users/tgutman/Documents/Tom/EUCANCan/golden-datasets-fork-tom/golden-datasets/scripts/
#OUT_DIR=/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/test_ingestions_sh

#bash $SCRIPT_DIR/ingest_snv.sh -t truth_dream/truth.snvs.synthetic.challenge.set1.chr.vcf.gz -s test_indel.vcf.gz -i test_snv.vcf.gz -o $OUT_DIR  -n curie_dream_1 -f /data/annotations/pipelines/Human/hg19_base/genome/hg19_base.fa


#/!\ ajouter etape de split multisample if it is the case
