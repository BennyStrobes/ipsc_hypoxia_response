import numpy as np 
import os
import sys
import pdb

# Encode environmental variable according to binary format:
# 'binary': Encode state 'A' as 0 and state 'B' as 1. Ignore other samples
def extract_binary_environmental_variable(letter):
    if letter == 'A':
        return 0
    elif letter == 'B':
        return 1
    else:
        return 'NaN'

# Encode environmental variable according to three_state format:
# 'three_state': Encode state 'A' as 0, state 'B' as 1, and state 'C' as 2
def extract_three_state_environmental_variable(letter):
    if letter == 'A':
        return 0
    elif letter == 'B':
        return 1
    elif letter == 'C':
        return 2
    else:
        return 'NaN'
input_directory = sys.argv[1]  # Directory containing all CHT input files
output_file = sys.argv[2]  # Output file
environmental_variable_form = sys.argv[3]  # String option describing how to parameterize the environmetal variable


t = open(output_file, 'w')  # Open output file handle
t.write('sample_id\tenvironmental_variable\tcht_input_file\n')

for file_name in sorted(os.listdir(input_directory)):
    if file_name.startswith('haplotype_read_counts') == False:
        continue
    if file_name.endswith('.txt.gz') == False:
        continue
    # Extract sample id from filename
    sample_name = file_name.split('.')[1].split('H')[1]
    # Throw out sample replicates ( can be identified by having two underscores in sample name)
    if len(sample_name.split('_')) == 3:
        continue

    # Extract cell line
    cell_line = sample_name.split('_')[0]

    # Extract encoded environmental variable
    if environmental_variable_form == 'binary':
        env_variable = extract_binary_environmental_variable(sample_name.split('_')[1])
    elif environmental_variable_form == 'three_state':
        env_variable = extract_three_state_environmental_variable(sample_name.split('_')[1])

    # Throw out samples that do no not belong to the environmental variable encoding
    if env_variable == 'NaN':
        continue

    # Print information to output file
    t.write(sample_name + '\t' + str(env_variable) + '\t' + input_directory + file_name + '\n')
t.close()