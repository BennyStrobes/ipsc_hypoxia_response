#######################################################################################################
# This is the driver script that runs EAGLE model on ipsc hypoxia data
#######################################################################################################



###############################################################################
# Input Data
###############################################################################

# Directory created by "time_step_independent_qtl_pipelines" scripts
# Contains 1 file per sample with information on each test (variant, target region)
# Each file (sample) has the same number of lines (tests)
cht_input_file_dir="/work-zfs/abattle4/bstrober/ipsc_hypoxia_response/raw_data/"

# File containing conversions from ensamble ids to gene symbols
gencode_file="/work-zfs/abattle4/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/gencode.v19.annotation.gtf.gz"


###############################################################################
# Output directories (aasume all of these exist prior to starting analysis)
###############################################################################

# Root directory for this of all ipsc data based results
output_root="/work-zfs/abattle4/bstrober/ipsc_hypoxia_response/dynamic_qtl_pipelines/total_expression/"

# Directory containing necessary input files to qtl tests
input_data_dir=$output_root"input_data/"


# Directory containing text files with results from dynamic qtl analysis
qtl_results_dir=$output_root"qtl_results/"

# Directory containing visualization of results found in qtl_results_dir
qtl_visualization_dir=$output_root"qtl_visualization/"






###############################################################################
# Dynamic QTL Calling
###############################################################################


##########################################
# Step 1: Create joint_test_input_file
##########################################
# joint_test_input_file is a table with 1 row per sample
# each row has 3 columns:
#### 1. Sample id
#### 2. environmental variable
#### 3. Absolute directory to CHT input file for that sample
# NOTE: This script is very specific to our data
# Takes less than a minute to run
# $environmental_variable is a parameter that describes how we parameterize the environmental variable. So far, this is done with:
### 1. 'binary': Encode state 'A' as 0 and state 'B' as 1. Ignore other samples
### 2. 'three_state': Encode state 'A' as 0, state 'B' as 1, and state 'C' as 2
environmental_variable_form="three_state"
joint_test_input_file=$input_data_dir"joint_test_input_file_"$environmental_variable_form".txt"
if false; then
python create_joint_test_input_file.py $cht_input_file_dir $joint_test_input_file $environmental_variable_form
fi



##########################################
# Step 2: Learn genome wide hyperparameters
##########################################

# File to contain learned library size correction factors for each sample (samples ordered same as $joint_test_input_file)
correction_factor_file=$input_data_dir"library_size_correction_factor_"$environmental_variable_form".txt"


# takes around 30 minutes to run
if false; then
sh learn_genome_wide_hyperparameters.sh $joint_test_input_file $correction_factor_file
fi
