repo_path=/Users/daphnevanbeek/IdeaProjects/golden-datasets/scripts/
source ${repo_path}/../venv/bin/activate


function run {
    repo=$1
    prefix=$2
    SV_VCF=$3
    samplename=$4
    out_dir=$5
    truth_sv_df=$6

    out_snv_df=${out_dir}"/"${prefix}"_snv_df.csv"
    out_sv_df=${out_dir}"/"${prefix}"_sv_df.csv"
    reclass_snv_df=${out_dir}"/"${prefix}"_snv_reclass.csv"
    reclass_sv_df=${out_dir}"/"${prefix}"_snv_reclass.csv"
    
    python ${repo}/ingest.py ${SV_VCF} -samplename ${samplename} -outputfile ${out_sv_df}

    # For Tom: conditional test in bash script to test if the input dataframes are empty.
    python ${repo}/reclassification.py ${out_snv_df} ${out_sv_df} ${reclass_snv_df} ${reclass_sv_df}
    python ${repo}/compare_node_to_truth.py ${reclass_sv_df} ${truth_sv_df}
    
    

}


# Hartwig
run ${repo_path} "Hartwig" "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/vcfs/Hartwig/5.19/COLO829v003T.purple.sv.vcf.gz" "COLO829v003T" "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/dataframes/" "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/dataframes/truth.csv"

# BSC
#run ${repo_path} "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/vcfs/BSC_SV_COLO.tsv" "none" "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/dataframes/bsc.csv"

# Curie
#run ${repo_path} "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/vcfs/Curie/COLO_curie_SV.vcf" "COLO829T" "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/dataframes/curie.csv"

# Truth file
#run ${repo_path} "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/truth_file/COLO829_truth.tsv" "none" "/Users/daphnevanbeek/Projects/eucancan/COLO_RESULTS/dataframes/truth.csv"
