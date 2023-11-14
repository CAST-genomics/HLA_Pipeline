import regex as re
import pandas as pd
import vcf
import re
import regex as re
import os
import subprocess
from tqdm import tqdm
import io
import allel
import numpy as np
import time
import haptools
from haptools.data import GenotypesPLINK
import math

def Extract_Impact(string):
    for x in ['HIGH', 'LOW', 'MODERATE', 'MODIFIER']:
        if x in string:
            return re.match(string, x).group(0)

        
def HapDF_Generator(df, Impact):
    """
    input: a whole dataframe
    output: a list containing its corresponding haplotypes
    """
    # Create MODIFIER dataframe
    df_MODIFIER = df[df['Impact'].str.contains(Impact)]

    # Store genotype info into a list
    Impact_list = []
    for column in df.iloc[:, 8:-1]:
        Impact_list.append(df[column].to_list())
  
    # Get two haplotypes for each two samples
    hap_list_Impact = []
    for x in Impact_list:
        hap = []
        hap1 = ''
        hap2 = ''
        for i in x:
            hap1+=str(i[0])
            hap2+=str(i[1])
        hap = [hap1, hap2]
        hap_list_Impact.append(hap)  
    hap_list_Impact_ = [0,0,0,0,0,0,0,0]+hap_list_Impact+[0]
    return hap_list_Impact, hap_list_Impact_


def modify_string(input_string):
    # Split the input string by semicolon
    pairs = input_string.split(';')

    # Create an empty list to store modified pairs
    modified_pairs = []

    # Iterate through the pairs
    for pair in pairs:
        pair = pair.strip()
        if ' ' in pair:
            key_value = pair.split(' ')
            key = key_value[0]
            value = key_value[1:]
        else:
            key = pair  # If no space, the whole pair is the key
            value = []

        # Check if the key is 'transcript_id' or 'protein_id'
        if key == 'transcript_id' or key == 'protein_id':
            # Split the value by dot and keep the part before the dot
            value_parts = value[0].split('.')
            value[0] = value_parts[0]
            modified_pair = f'{key} {" ".join(value)}"'
        else:
            modified_pair = f'{key} {" ".join(value)}'

        modified_pairs.append(modified_pair)

    # Join the modified pairs back together with semicolons
    result = '; '.join(modified_pairs)
    return result


def write_fasta(output_filename, fasta_dict):
    with open(output_filename, 'w') as file:
        for header, sequence_list in fasta_dict.items():
            file.write(f'>{header}\n')
            file.write('\n'.join(sequence_list) + '\n')
            
            
def modify_protein_fasta(output_filename, gene_id):
    
    # Define the original FASTA file name
    input_file = "/frazer01/home/ximei/analysis/Ximei/new_dict/protein_old.fa"

    # Create a dictionary to store the original sequences
    sequences = {}

    # Read the input FASTA file
    with open(input_file, "r") as file:
        current_sequence = None
        for line in file:
            line = line.strip()
            #print(line)
            #if line[39:58] == 'ENSG00000206503.13_6':
            if line.startswith(">"):
                # Start of a new sequence
                current_sequence = line[1:]
                sequences.setdefault(current_sequence, [])
            else:
                # Append sequence lines
                sequences[current_sequence].append(line)
                
    # get the aligned keys
    keys_l = list(sequences.keys())
    keys_l_split = []
    for x in keys_l:
        keys_l_split.append(x.split('|'))
    
    # get the edited and index keys
    indice_keep = []
    indice_keep_edited = []
    for x in keys_l_split:
        if x[2][:-5] == gene_id.split('.')[0]: # ENSG00000206503.13_6
            indice_keep.append(x)
            indice_keep_edited.append(x[0][:-2])

            
    # Transform back to original key form
    indice_keep_ori = ['|'.join(x[:8]) + '|' for x in indice_keep]
    
    # Get values according to the key
    values = []
    for key in indice_keep_ori:
        values.append(sequences[key])
        
    # Get the ready to export dictionary
    final_dict = dict(zip(indice_keep_edited, values))
    write_fasta(output_filename, final_dict)
    print(f'FASTA file "{output_filename}" has been created.')
    return
    

def modify_cds_fasta(output_filename, gene_id):
    # Define the input FASTA file name
    input_file = "/frazer01/home/ximei/analysis/Ximei/new_dict/cds_old.fa"

    # Create a dictionary to store the sequences
    sequences = {}

    # Read the input FASTA file
    with open(input_file, "r") as file:
        current_sequence = None
        for line in file:
            line = line.strip()
            #print(line)
            #if line[39:58] == 'ENSG00000206503.13_6':
            if line.startswith(">"):
                # Start of a new sequence
                current_sequence = line[1:]
                sequences.setdefault(current_sequence, [])
            else:
                # Append sequence lines
                sequences[current_sequence].append(line)

    # get the aligned keys
    keys_l = list(sequences.keys())
    keys_l_split = []
    for x in keys_l:
        keys_l_split.append(x.split('|'))
    
    # get the common gene-id
    indice_keep = []
    indice_keep_edited = []
    for x in keys_l_split:
        if x[1][:-5] == gene_id.split('.')[0]:
            indice_keep.append(x)
            indice_keep_edited.append(x[0].split('.')[0])
         
            
    # Transform back to original key form
    indice_keep_ori = ['|'.join(x[:8]) + '|' for x in indice_keep]
    
    # Get the values
    values = []
    for key in indice_keep_ori:
        values.append(sequences[key])
    
    # Get the final dictionry based on values and keys
    final_dict = dict(zip(indice_keep_edited, values))
    write_fasta(output_filename, final_dict)
    print(f'FASTA file "{output_filename}" has been created.')
    return



def pgen_generator(index,df_new,output_directory,vcf_dict,df_vcf_new_HAPLOTYPE_index, low_pos,impact):
    """
    """
    input_matrix = np.tile(df_new['index'].to_list()[index[0]:index[1]], (int(len(df_new.HAPLOTYPE.unique())/5)+1, 1))
    final_input_matrix = list(input_matrix) * 5
    final_input_matrix = final_input_matrix[:len(df_new['HAPLOTYPE'].value_counts().to_list())]
    
    matrix = np.array(final_input_matrix)
    list_to_compare = df_vcf_new_HAPLOTYPE_index
    comparison_matrix = (matrix == np.array(list_to_compare)[:, np.newaxis]) # Create a broadcasted boolean matrix
    result_matrix = comparison_matrix.astype(int)# Convertresult_matrix = result_matrix.reshape(int(result_matrix.shape[0]),int(result_matrix.shape[1]/2),2)
    result_matrix = result_matrix.reshape(int(result_matrix.shape[0]),int(result_matrix.shape[1]/2),2)

    
    df_header = pd.DataFrame(df_new['HAPLOTYPE'].unique())
    df_header['CHR'] = 'chr6'
    if math.isnan(low_pos):
        low_pos = 100
    df_header['POS'] = low_pos
    df_header['REF/ALT'] = df_header.apply(lambda row: ("A", "T"), axis=1)
    df_header_np = df_header.to_numpy()
    df_header_as_tuples = list(map(tuple, df_header_np.tolist()))
    print(df_header_as_tuples[:5])
    
    gts = GenotypesPLINK(os.path.join(output_directory,'{}_{}_{}.pgen'.format(str(index[0]),str(index[1]), impact)))
    gts.data = np.transpose(result_matrix, (1, 0, 2)) # replace this with your numpy array
    gts.variants = np.array(df_header_as_tuples, dtype=gts.variants.dtype)# Provide the ID, CHROM, POS, and (REF, ALT) columns in another numpy array
    gts.samples = [item for item in [*map(lambda x: x.split('_')[0], vcf_dict['samples'][int(index[0]/2):int(index[1]/2)])]] # now, provide the sample IDs:
    gts.write() # finally, we can write the PGEN file
    
    
def interval_generator(unique_hap_count):
    """
    the [ unique_hap_count * interval_size !> 3300000000 ] will not kill memory
    """
    interval_size = 3300000000 // unique_hap_count
    
    start = 0
    end = 974818
    
    if interval_size % 2 == 0:
        step = interval_size
    else:
        step = interval_size - 1

    iter_list = []

    while start < end:
        iter_list.append([start, min(start + step, end)])
        start += step
    return iter_list

def get_flattened_list(array):
    if array.size == 0:
        return np.array([[1] * 584])  # All samples; or any other value you want to assign when the size is 0
        #return np.array([[1] * 1996])  # test samples: 20
    else:
        return array.reshape((array.shape[0], -1))
    