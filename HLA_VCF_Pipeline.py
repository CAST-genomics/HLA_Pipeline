from HLA_Util import Extract_Impact
from HLA_Util import HapDF_Generator
from HLA_Util import modify_string
from HLA_Util import write_fasta
from HLA_Util import modify_protein_fasta
from HLA_Util import modify_cds_fasta
from HLA_Util import pgen_generator
from HLA_Util import interval_generator
from HLA_Util import get_flattened_list
from collections import Counter
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

class HLA_VCF_Pipeline:
    def __init__(self, gene_name, geneIndex_File, output_directory):
        self.gene_name = gene_name
        self.geneIndex_File = geneIndex_File
        self.output_directory = output_directory
        self.FROM = None
        self.TO = None
        self.PREFIX = None
        self.stats = None
        self.out =  None
        self.MODIFIER_MinPos = None
        self.LOW_MinPos = None
        self.HIGH_MinPos = None
        self.MODERATE_MinPos = None
        self.df = None
        self.hap_list_MODIFIER = None
        self.hap_list_MODIFIER_ = None
        self.hap_list_LOW = None
        self.hap_list_LOW_ = None
        self.hap_list_HIGH = None
        self.hap_list_HIGH_ = None
        self.hap_list_MODERATE = None
        self.hap_list_MODERATE_ = None
        self.df_MODIFIER = None
        self.df_LOW = None
        self.df_HIGH = None
        self.df_MODERATE = None
        self.my_dict_MODIFIER = None
        self.my_dict_low = None
        self.my_dict_high = None
        self.my_dict_MODERATE = None
        self.concatenated_df = None
        self.vcf_dict = None
        self.df_gene = None
        self.df_header_annotated = None
        self.filter_out = None
        self.df_header_filtered = None
        self.df_new = None
        self.LOW_POS = None
        self.MODERATE_POS = None
        self.MODIFIER_POS = None
        self.HIGH_POS = None
        self.vcf_dict = None
        self.transformed_LOW_haplotypes = None
        self.transformed_MODERATE_haplotypes = None
        self.transformed_HIGH_haplotypes = None
        self.transformed_MODIFIER_haplotypes = None
        self.pos_list = []
        self.df_vcf_new_HAPLOTYPE_index = None
        self.df_HIGH_MODERATE_LOW = None
        self.df_HIGH_MODERATE = None
        self.df_HIGH = None
    
    def GenerateVCF(self):
        """
        input: pgen file
        steps: 0. get coordinates based on gene_name
               1. convert PGEN into VCF using PLINK
               2. index VCF
               3. Phase VCF
               4. index VCF
               5. Shapeit VCF
        output: phased & indexed VCF
        """
        # Create folder
        self.output_directory = os.path.join(self.output_directory,self.gene_name)
        self.final_output = os.path.join(self.output_directory,'result')
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)
        else:
            print('The folder {} already exists.'.format(self.output_directory))
            
        if not os.path.exists(self.final_output):
            os.makedirs(self.final_output)
        else:
            print('The folder {} already exists.'.format(self.final_output))
        
        #**************************** step 0: Identify position based on geneID *************************************************************
        self.df_gene = pd.read_csv(self.geneIndex_File, sep='\t')
        self.FROM = self.df_gene[self.df_gene['gene_name'] == '{}'.format(self.gene_name)]['start'].values[0]
        self.TO = self.df_gene[self.df_gene['gene_name'] == '{}'.format(self.gene_name)]['end'].values[0]
        print(self.FROM)
        print(self.TO)
        self.pos_list.append(self.FROM)
        self.pos_list.append(self.TO)
        self.geneID = self.df_gene[self.df_gene['gene_name'] == '{}'.format(self.gene_name)].gene_id.values[0] # find gene ID
        if self.df_gene[self.df_gene['gene_name'] == self.gene_name]['strand'].to_list()[0] == '-':
            self.TO += 10000
        if self.df_gene[self.df_gene['gene_name'] == self.gene_name]['strand'].to_list()[0] == '+':
            self.FROM -=  10000
        print(self.FROM)
        print(self.TO)
        self.pos_list.append(self.FROM)
        self.pos_list.append(self.TO)
        #*************************** STEP 1: convert PGEN into VCF using PLINK **************************************************************
        self.PREFIX = os.path.join(self.output_directory,self.gene_name)
        bcftools_command1 = "/frazer01/software/plink-2.3/plink2_64 --indep-pairwise 250000 1000 0.8 --out {} --memory 128000 --threads 16 --pfile /frazer01/projects/CEGS/pipeline/hla_genotypes/mhc.imputed --export vcf bgz --chr 6 --from-bp {} --to-bp {}".format(self.PREFIX, self.FROM, self.TO)
        
        # Step 1: convert PGEN into VCF using PLINK
        if not os.path.exists(f"{self.PREFIX}.vcf.gz") and not os.path.exists(f"{self.PREFIX}.log"):
            print('Step 1: convert PGEN into VCF using PLINK')
            try:
                subprocess.run(bcftools_command1, shell=True, check=True)
                print("Step 1: convert PGEN into VCF using PLINK completed successfully.")
            except subprocess.CalledProcessError as e:
                print("Error occurred in step 1: convert PGEN into VCF using PLINK")
        
        
        #*************************** STEP 1.1: Filter AF > 1% (PLINK) ********************************************************************************
        self.filter_out = os.path.join(self.output_directory, 'Filtered_{}'.format(self.gene_name))
        bcftools_command1_1 = '/frazer01/software/plink-2.3/plink2_64 --memory 128000 --threads 16 --vcf {} --maf 0.01 --recode vcf --out {}'.format(self.PREFIX +'.vcf.gz', self.filter_out)
        # Step 1.1: Filter AF > 1% (PLINK)
        if not os.path.exists(self.filter_out+'.vcf'):
            print('Starting Step 1.1: Filter AF > 1% (PLINK)')
            try:
                subprocess.run(bcftools_command1_1, shell=True, check=True)
                print("Step 1.1: Filter AF > 1% (PLINK) completed successfully.")
            except subprocess.CalledProcessError as e:
                print("Error occurred in step 1.1: Filter AF > 1% (PLINK)")
                
        #*************************** STEP 1.2: Zip filtered File  ********************************************************************************
        bcftools_command1_2 = "bgzip -c {} > {}".format(self.filter_out + '.vcf', self.filter_out + '.vcf.gz')
        # Step 1.2: Filter AF > 1% (PLINK)
        if not os.path.exists(self.filter_out + '.vcf.gz'):
            print('Starting Step 1.2:  Zip filtered File')
            try:
                subprocess.run(bcftools_command1_2, shell=True, check=True)
                print("Step 1.2: Zip filtered File ")
            except subprocess.CalledProcessError as e:
                print("Error occurred in step 1.2: Zip filtered File")
        
        #******************************* STEP 2: index VCF *****************************************************************************************
        bcftools_command2 = "bcftools index -t {}".format(self.filter_out+'.vcf.gz')
        
        # Step 2: index VCF
        if not os.path.exists(f"{self.filter_out}.tbi"):
            print('Starting step 2: index VCF')
            try:
                subprocess.run(bcftools_command2, shell=True, check=True)
                print("Step 2 completed successfully.")
            except subprocess.CalledProcessError as e:
                print("Error occurred in step 2: index VCF")
        
        
        #*********************************  STEP 3: Phase VCF  ****************************************************************************************
        self.output = os.path.join(self.output_directory, self.gene_name +'_filtered_phased.vcf.gz')
        bcftools_command3 = "shapeit4 --input {} --thread 16 --map /frazer01/projects/CEGS/analysis/ukbb_hla_type_gwas/input/phase/chr6.b37.gmap.gz --region 6:{}-{} --log {} --output {}".format(self.filter_out+ '.vcf.gz' , self.FROM, self.TO, self.gene_name + '_filtered.log', self.output)
        
        #Step 3: Phase VCF
        if not os.path.exists(self.output):
            print('Starting step 3: Phase VCF')
            try:
                subprocess.run(bcftools_command3, shell=True, check=True)
                print("Step 3 completed successfully")
            except subprocess.CalledProcessError as e:
                print("Error occurred in step 3: Phase VCF")
        
        
        #********************************* STEP 4: index VCF  *******************************************************************************************
        bcftools_command4 = "bcftools index -t {}".format(self.output) 
        
        # Step 4: index VCF
        if not os.path.exists(f"{self.output}.tbi"):
            print('Starting step 4: index VCF')
            try:
                subprocess.run(bcftools_command4, shell=True, check=True)
                print("Step 4 completed successfully")
            except subprocess.CalledProcessError as e:
                print("Error occurred in step 4: index VCF")
        
        #********************************* STEP 5 Build snpeff database ********************************************************************************
        # 5_1: make genome directory
        bcftools_command_5_1 = 'mkdir /frazer01/home/ximei/snpEff/data/{}'.format(self.gene_name) 
        
        # 5_2: make gtf file
        gtf_out = os.path.join('/frazer01/home/ximei/snpEff/data/' + self.gene_name ,  'genes.gtf')
        bcftools_command_5_2 = 'grep {} /frazer01/reference/public/Gencode.v34lift37/gencode.v34lift37.annotation.gtf > {}'.format(self.geneID, gtf_out)
        
        # 5_3: go to snpEff directory
        bcftools_command_5_3 = 'cd /frazer01/home/ximei/snpEff/data/{}'.format(self.gene_name)
        
        # 5_4: link sequences.fa
        bcftools_command_5_4 = "ln /home/ximei/snpEff/data/hla_genome/sequences.fa /home/ximei/snpEff/data/{}/sequences.fa".format(self.gene_name)
        
        # step 5: build genome database
        if not os.path.exists('/frazer01/home/ximei/snpEff/data/{}'.format(self.gene_name)):
            print('Starting step 5: build genome database')
            try:
                subprocess.run(bcftools_command_5_1, shell=True, check=True)
                print("Step 5_1 completed successfully")
            except subprocess.CalledProcessError as e:
                print("Error occurred in step 5_1")
            try:
                subprocess.run(bcftools_command_5_2, shell=True, check=True)
                print("Step 5_2 completed successfully")
            except subprocess.CalledProcessError as e:
                print("Error occurred in step 5_2")
            try:
                subprocess.run(bcftools_command_5_3, shell=True, check=True)
                print("Step 5_3 completed successfully")
            except subprocess.CalledProcessError as e:
                print("Error occurred in step 5_3")
            try:
                subprocess.run(bcftools_command_5_4, shell=True, check=True)
                print("Step 5_4 completed successfully")
            except subprocess.CalledProcessError as e:
                print("Error occurred in step 5_4")
                
        # 5_5: edit gene file
        self.df_gene = pd.read_csv('/home/ximei/snpEff/data/{}/genes.gtf'.format(self.gene_name),sep = '\t', header=None)
        self.df_gene[8] = self.df_gene[8].apply(modify_string)
        self.df_gene.to_csv('/home/ximei/snpEff/data/{}/genes.gtf'.format(self.gene_name), sep='\t', header=False, index=False, quoting=3, escapechar='\\')
        print("Step 5_5 completed successfully")
        
        # 5_6: edit fasta file
        output_filename = '/frazer01/home/ximei/snpEff/data/{}/protein.fa'.format(self.gene_name)
        print(self.geneID)
        modify_protein_fasta(output_filename, self.geneID)
        print("Step 5_6 completed successfully")
        
        # 5_7: edit cds file
        output_filename = '/frazer01/home/ximei/snpEff/data/{}/cds.fa'.format(self.gene_name)
        modify_cds_fasta(output_filename, self.geneID)
        print("Step 5_7 completed successfully")

        
        # 5_8: build database
        bcftools_command_5_8 = '/software/jdk-16.0.1/bin/java -jar /home/ximei/snpEff/snpEff.jar build -gtf22 -v {} -noCheckCds -noCheckProtein'.format(self.gene_name)

        try:
            subprocess.run(bcftools_command_5_8, shell=True, check=True)
            print("Step 5_8 completed successfully")
        except subprocess.CalledProcessError as e:
            print("Error occurred in step 5_8")
        
        #****************************************** step 6: Run snpeff  *******************************************************************************
        self.stats = os.path.join(self.output_directory, 'report_' + self.gene_name + '.html')
        self.out = os.path.join(self.output_directory, 'annoted_' + self.gene_name + '_phased.vcf')
        bcftools_command6 = "/software/jdk-16.0.1/bin/java -Xmx100g -jar /home/ximei/snpEff/snpEff.jar -v -stats {} {} {} > {}".format(self.stats, self.gene_name, self.output, self.out)
        
        # Step 6: Run snpeff
        if not os.path.exists(self.out):
            print('Starting step 6: Run snpeff')
            try:
                subprocess.run(bcftools_command6, shell=True, check=True)
                print("Step 6 completed successfully")
            except subprocess.CalledProcessError as e:
                print("Error occurred in step 6: Run snpeff")
        else:
            print("Skipping step 6: Run snpeff because the output file already exists")

        
        #****************************************** step 7: generate annotation txt file  *************************************************************
        annotated_vcf = self.out
        annotation_txt = os.path.join(self.output_directory, 'INFO_annoted_' + self.gene_name + '_phased.txt')
        bcftools_command7 = "bcftools query -f '%INFO/ANN\n' {} > {}".format(annotated_vcf, annotation_txt)
            
        # Step 7: generate annotation txt file
        if not os.path.exists(annotation_txt):
            print('Starting step 7: generate annotation txt file')
            try:
                subprocess.run(bcftools_command7, shell=True, check=True)
                print("step 7 completed successfully")
            except subprocess.CalledProcessError as e:
                print("Error occurred in step 7: generate annotation txt file")
        else:
            print("Skipping step step 7: generate annotation txt file because the output file already exists")
        return self
    
    def readVCF(self):
        """
        input: phased & indexed VCF
        output: df read as VCF
        """
        print('step 1: Starting to read VCF')
        start_time = time.time()
        vcf_file_path = self.out
        self.vcf_dict = allel.read_vcf(vcf_file_path)
        elapsed_time = time.time() - start_time
        print(f"Elapsed Time: {elapsed_time:.2f} seconds")
  
        df_header = pd.DataFrame()
        df_header['CHROM'] = self.vcf_dict['variants/CHROM']
        df_header['POS'] = self.vcf_dict['variants/POS']
        df_header['ID'] = self.vcf_dict['variants/ID']
        df_header['REF'] = self.vcf_dict['variants/REF']
        df_header['ALT'] = np.array(self.vcf_dict['variants/ALT']).tolist()
        df_header['QUAL'] = self.vcf_dict['variants/QUAL']
        df_header['FILTER'] = self.vcf_dict['variants/FILTER_PASS']
        print('df_header is shape:{}'.format(df_header.shape))
        
        self.df_INFO = pd.read_fwf('/home/ximei/analysis/Ximei/HLA_result/{}/INFO_annoted_{}_phased.txt'.format(self.gene_name, self.gene_name), names = ['INFO'])
        self.df_INFO['first_rowish'] = self.df_INFO['INFO'].apply(lambda x: x[0:100])
        self.df_INFO['Impact'] = self.df_INFO['first_rowish'].apply(lambda x: Extract_Impact(x))
        print('df_INFO is shape:{}'.format(self.df_INFO.shape))
        
        self.df_header_annotated = pd.concat([df_header,self.df_INFO], axis = 1)
        print('df_header_annotated is shape:{}'.format(self.df_header_annotated.shape))
        
        l = self.pos_list
        print('l:{}'.format(l))
        l_c = Counter(l)
        print('l_c:{}'.format(l_c))
        ll = list(set(l))
        print('ll:{}'.format(ll))
        ll_c = Counter(list(set(l)))
        print('ll_c:{}'.format(ll_c))
        index_list = list(ll_c-(l_c-ll_c))
        print('index_list:{}'.format(index_list))
        print('min(index_list):{}'.format(min(index_list)))
        print('max(index_list):{}'.format(max(index_list)))
        
        # HIGH + MODERATE + MODIFIER + LOW
        self.df_filtered_modifier = self.df_header_annotated[(self.df_header_annotated['Impact'] == 'MODIFIER') & (self.df_header_annotated['POS'] > min(index_list)) &(self.df_header_annotated['POS']< max(index_list))]
        df_Low_High_Moderate = self.df_header_annotated[self.df_header_annotated['Impact'] != 'MODIFIER']
        self.df_header_filtered = pd.concat([self.df_filtered_modifier,df_Low_High_Moderate], axis = 0)
        self.df_header_filtered['index'] = list(range(0, self.df_header_filtered.shape[0]))
        print(self.df_header_filtered.shape)
        print(self.df_header_filtered)
        
        return self
    
    def Haplotypes(self):
        """
        1. Get high/low/moderate/modifier row index
        2. Generate haplotypes for each Impact for each sample
        """
        print('step 2.1: Starting Get high/low/moderate/modifier row index')
        
        LOW_index = np.array(self.df_header_filtered[self.df_header_filtered['Impact'] == 'LOW']['index'])
        MODERATE_index = np.array(self.df_header_filtered[self.df_header_filtered['Impact'] == 'MODERATE']['index'])
        HIGH_index = np.array(self.df_header_filtered[self.df_header_filtered['Impact'] == 'HIGH']['index'])
        MODIFIER_index = np.array(self.df_header_filtered[self.df_header_filtered['Impact'] == 'MODIFIER']['index'])
        
        print('step 2.2: Starting get sub-matrix')
        index_list_ = self.df_header_filtered.index.to_list()
        LOW_list = self.vcf_dict['calldata/GT'][index_list_][LOW_index]
        MODERATE_list = self.vcf_dict['calldata/GT'][index_list_][MODERATE_index]
        HIGH_list = self.vcf_dict['calldata/GT'][index_list_][HIGH_index]
        MODIFIER_list = self.vcf_dict['calldata/GT'][index_list_][MODIFIER_index]
        print(LOW_list.shape)
        print(MODERATE_list.shape)
        print(HIGH_list.shape)
        print(MODIFIER_list.shape)
        
        print('step 2.3: flatten the list based on each row (variant)')
        LOW_list_flattend = get_flattened_list(LOW_list)
        MODERATE_list_flattend = get_flattened_list(MODERATE_list)
        HIGH_list_flattend = get_flattened_list(HIGH_list)
        MODIFIER_list_flattend = get_flattened_list(MODIFIER_list)
        print(LOW_list_flattend.shape)
        print(MODERATE_list_flattend.shape)
        print(HIGH_list_flattend.shape)
        print(MODIFIER_list_flattend.shape)
        
        print('step 2.4: Get the transpose of matrix (the haplotypes)')
        LOW_list_flattend_transposed = np.transpose(LOW_list_flattend)
        MODERATE_list_flattend_transposed = np.transpose(MODERATE_list_flattend)
        HIGH_list_flattend_transposed = np.transpose(HIGH_list_flattend)
        MODIFIER_list_flattend_transposed = np.transpose(MODIFIER_list_flattend)
        print(LOW_list_flattend_transposed.shape)
        print(MODERATE_list_flattend_transposed.shape)
        print(HIGH_list_flattend_transposed.shape)
        print(MODIFIER_list_flattend_transposed.shape)
        
        print('step 2.5: get haplotypes: Convert each row to a binary string')
        #get haplotypes: Convert each row to a binary string
        start_time = time.time()
        LOW_list_flattend_transposed_binary = np.apply_along_axis(lambda row: ''.join(map(str, row)), axis=1, arr=LOW_list_flattend_transposed)
        print('LOW_list_flattend_transposed_binary done')
        MODERATE_list_flattend_transposed_binary = np.apply_along_axis(lambda row: ''.join(map(str, row)), axis=1, arr=MODERATE_list_flattend_transposed)
        print('MODERATE_list_flattend_transposed_binary done')
        HIGH_list_flattend_transposed_binary = np.apply_along_axis(lambda row: ''.join(map(str, row)), axis=1, arr=HIGH_list_flattend_transposed)
        print('HIGH_list_flattend_transposed_binary done')
        MODIFIER_list_flattend_binary = np.apply_along_axis(lambda row: ''.join(map(str, row)), axis=1, arr=MODIFIER_list_flattend_transposed)
        print('MODIFIER_list_flattend_binary done')
        elapsed_time = time.time() - start_time
        print(f"Elapsed Time: {elapsed_time:.2f} seconds")

        print(LOW_list_flattend_transposed_binary.shape)
        print(MODERATE_list_flattend_transposed_binary.shape)
        print(HIGH_list_flattend_transposed_binary.shape)
        print(MODIFIER_list_flattend_binary.shape)
        
        print('step 2.6: get haplotypes unique values')
        
        # get haplotypes unique values
        unique_LOW_haplotypes = set(LOW_list_flattend_transposed_binary)
        print(len(unique_LOW_haplotypes))
        unique_MODERATE_haplotypes = set(MODERATE_list_flattend_transposed_binary)
        print(len(unique_MODERATE_haplotypes))
        unique_HIGH_haplotypes = set(HIGH_list_flattend_transposed_binary)
        print(len(unique_HIGH_haplotypes))
        unique_MODIFIER_haplotypes = set(MODIFIER_list_flattend_binary)
        print(len(unique_MODIFIER_haplotypes))
        
        
        print('step 2.7: Create haplotype dictionary')
        LOW_index_list = list(range(1, len(unique_LOW_haplotypes)+1))
        MODERATE_index_list = list(range(1, len(unique_MODERATE_haplotypes)+1))
        HIGH_index_list = list(range(1, len(unique_HIGH_haplotypes)+1))
        MODIFIER_index_list = list(range(1, len(unique_MODIFIER_haplotypes)+1))

        LOW_dictionary = dict(zip(unique_LOW_haplotypes,LOW_index_list))
        MODERATE_dictionary = dict(zip(unique_MODERATE_haplotypes,MODERATE_index_list))
        HIGH_dictionary = dict(zip(unique_HIGH_haplotypes,HIGH_index_list))
        MODIFIER_dictionary = dict(zip(unique_MODIFIER_haplotypes,MODIFIER_index_list))
        
        
        print('step 2.8: Transform original haplotypes based on dictioanry')
        self.transformed_LOW_haplotypes = list(map(lambda item: LOW_dictionary[item], LOW_list_flattend_transposed_binary))
        self.transformed_MODERATE_haplotypes = list(map(lambda item: MODERATE_dictionary[item], MODERATE_list_flattend_transposed_binary))
        self.transformed_HIGH_haplotypes = list(map(lambda item: HIGH_dictionary[item], HIGH_list_flattend_transposed_binary))
        self.transformed_MODIFIER_haplotypes = list(map(lambda item: MODIFIER_dictionary[item], MODIFIER_list_flattend_binary))
        
        return self
    
    

    def PrepareVCF(self):
        """
        Generate vcf based on above dataframe
        """
   
        
        # get position info
        print('step 3.1: get position info')
        self.LOW_POS = self.df_header_annotated[self.df_header_annotated['Impact'] == 'LOW']['POS'].min()
        self.MODERATE_POS = self.df_header_annotated[self.df_header_annotated['Impact'] == 'MODERATE']['POS'].min()
        self.MODIFIER_POS = self.df_header_annotated[self.df_header_annotated['Impact'] == 'MODIFIER']['POS'].min()
        self.HIGH_POS = self.df_header_annotated[self.df_header_annotated['Impact'] == 'HIGH']['POS'].min()
        
        # creating df_new
        print('step 3.2: creating df_new')
        self.df_new = pd.DataFrame(columns = ['GENE', 'SAMPLE', 'ALLELE', 'HAPLOTYPE', 'HIGH', 'MODERATE', 'LOW', 'MODIFIER'])
        self.df_new['SAMPLE'] = [item for item in [*map(lambda x: x.split('_')[0], self.vcf_dict['samples'])] for _ in range(2)]
        self.df_new['GENE'] = [self.gene_name] * self.df_new.shape[0]
        self.df_new['ALLELE'] = [1,2] * int(self.df_new.shape[0]/2)
        self.df_new['LOW'] = self.transformed_LOW_haplotypes
        self.df_new['HIGH'] = self.transformed_HIGH_haplotypes
        self.df_new['MODERATE'] = self.transformed_MODERATE_haplotypes
        self.df_new['MODIFIER'] = self.transformed_MODIFIER_haplotypes
        self.df_new['HAPLOTYPE'] = self.df_new.apply(lambda x: f"{x['HIGH']}:{x['MODERATE']}:{x['LOW']}:{x['MODIFIER']}", axis = 1)
        df_vcf_new_HAPLOTYPE = self.df_new['HAPLOTYPE'].unique()# Assuming df_vcf_new_HAPLOTYPE is a list or array of sorted unique elements

        # Create a mapping from unique elements to integers
        element_to_index = {element: index for index, element in enumerate(df_vcf_new_HAPLOTYPE)}
        self.df_new['index'] = self.df_new['HAPLOTYPE'].apply(lambda x: element_to_index[x])
        self.df_new.to_csv(os.path.join(self.final_output, 'table.csv'))
        print('df_new shape:{}'.format(self.df_new.shape))
        print(self.df_new)
        
        self.df_vcf_new_HAPLOTYPE_index = list(map(element_to_index.get, df_vcf_new_HAPLOTYPE))
        
        # Creating df_HIGH_MODERATE_LOW
        self.df_HIGH_MODERATE_LOW = self.df_new[['GENE','SAMPLE', 'ALLELE', 'HAPLOTYPE', 'HIGH', 'MODERATE', 'LOW']]
        self.df_HIGH_MODERATE_LOW['HAPLOTYPE'] = self.df_HIGH_MODERATE_LOW.apply(lambda x: f"{x['HIGH']}:{x['MODERATE']}:{x['LOW']}", axis = 1)
        df_HIGH_MODERATE_LOW_HAPLOTYPE = self.df_HIGH_MODERATE_LOW['HAPLOTYPE'].unique()
        element_to_index = {element: index for index, element in enumerate(df_HIGH_MODERATE_LOW_HAPLOTYPE)}
        self.df_HIGH_MODERATE_LOW['index'] = self.df_HIGH_MODERATE_LOW['HAPLOTYPE'].apply(lambda x: element_to_index[x])
        self.df_HIGH_MODERATE_LOW_HAPLOTYPE_index = list(map(element_to_index.get, df_HIGH_MODERATE_LOW_HAPLOTYPE))
        print('df_HIGH_MODERATE_LOW shape:{}'.format(self.df_HIGH_MODERATE_LOW.shape))
        print(self.df_HIGH_MODERATE_LOW)
        
        # Creating df_HIGH_MODERATE
        self.df_HIGH_MODERATE = self.df_new[['GENE','SAMPLE', 'ALLELE', 'HAPLOTYPE', 'HIGH', 'MODERATE']]
        self.df_HIGH_MODERATE['HAPLOTYPE'] = self.df_HIGH_MODERATE.apply(lambda x: f"{x['HIGH']}:{x['MODERATE']}", axis = 1)
        df_HIGH_MODERATE_HAPLOTYPE = self.df_HIGH_MODERATE['HAPLOTYPE'].unique()
        element_to_index = {element: index for index, element in enumerate(df_HIGH_MODERATE_HAPLOTYPE)}
        self.df_HIGH_MODERATE['index'] = self.df_HIGH_MODERATE['HAPLOTYPE'].apply(lambda x: element_to_index[x])
        self.df_HIGH_MODERATE_HAPLOTYPE_index = list(map(element_to_index.get, df_HIGH_MODERATE_HAPLOTYPE))
        print('df_HIGH_MODERATE shape:{}'.format(self.df_HIGH_MODERATE.shape))
        print(self.df_HIGH_MODERATE)
        
        # Creating df_HIGH
        self.df_HIGH = self.df_new[['GENE','SAMPLE', 'ALLELE', 'HAPLOTYPE', 'HIGH']]
        #print('1')
        self.df_HIGH['HAPLOTYPE'] = self.df_HIGH.apply(lambda x: f"{x['HIGH']}", axis = 1)
        #self.df_HIGH['HAPLOTYPE'] = self.df_HIGH['HIGH']
        #print('2')
        df_HIGH_HAPLOTYPE = self.df_HIGH['HAPLOTYPE'].unique()
        #print('3')
        element_to_index = {element: index for index, element in enumerate(df_HIGH_HAPLOTYPE)}
        #print('4')
        self.df_HIGH['index'] = self.df_HIGH['HAPLOTYPE'].apply(lambda x: element_to_index[x])
        #print('5')
        self.df_HIGH_HAPLOTYPE_index = list(map(element_to_index.get, df_HIGH_HAPLOTYPE))
        print('df_HIGH shape:{}'.format(self.df_HIGH.shape))
        print(self.df_HIGH)
        
        return self

    def Generate_Pgen(self):
        """
        Create the PGen
        1. create each individual pgen file
        """
        #*************************************************  Get the iter_list on all IMPACT  **************************************************************************
        print('step 4.1: Starting to Create the pegn for all impact')
        unique_hap_count = len(self.df_new['HAPLOTYPE'].unique())
        iter_list = interval_generator(unique_hap_count)
        
        for index in iter_list:
            # if file exists, don't create again
            if os.path.exists(os.path.join(self.output_directory, '{}_{}_ALL.pgen'.format(str(index[0]),str(index[1])))):
                print("Skipping processing for {}, file already exists.".format(str(index)))
            else:
                print('start processing {}'.format(str(index)))
                pgen_generator(index,self.df_new,self.output_directory,self.vcf_dict,self.df_vcf_new_HAPLOTYPE_index,self.LOW_POS,'ALL')
                
                
        #**************************************************  Get the iter_list on HIGH_MODERATE_LOW  *****************************************************************
        print('step 4.2: Starting to Create the pegn on HIGH_MODERATE_LOW')
        unique_hap_count = len(self.df_HIGH_MODERATE_LOW['HAPLOTYPE'].unique())
        iter_list = interval_generator(unique_hap_count)
        
        for index in iter_list:
            # if file exists, don't create again
            if os.path.exists(os.path.join(self.output_directory, '{}_{}_HIGH_MODERATE_LOW.pgen'.format(str(index[0]),str(index[1])))):
                print("Skipping processing for {}_HIGH_MODERATE_LOW, file already exists.".format(str(index)))
            else:
                print('start processing {}_HIGH_MODERATE_LOW'.format(str(index)))
                pgen_generator(index,self.df_HIGH_MODERATE_LOW,self.output_directory,self.vcf_dict,self.df_HIGH_MODERATE_LOW_HAPLOTYPE_index,self.LOW_POS,'HIGH_MODERATE_LOW')
        
        
        #**************************************************  Get the iter_list on HIGH_MODERATE  ********************************************************************
        print('step 4.3: Starting to Create the pegn on HIGH_MODERATE')
        unique_hap_count = len(self.df_HIGH_MODERATE['HAPLOTYPE'].unique())
        iter_list = interval_generator(unique_hap_count)
        
        for index in iter_list:
            # if file exists, don't create again
            if os.path.exists(os.path.join(self.output_directory, '{}_{}_HIGH_MODERATE.pgen'.format(str(index[0]),str(index[1])))):
                print("Skipping processing for {}_HIGH_MODERATE, file already exists.".format(str(index)))
            else:
                print('start processing {}_HIGH_MODERATE'.format(str(index)))
                pgen_generator(index,self.df_HIGH_MODERATE,self.output_directory,self.vcf_dict,self.df_HIGH_MODERATE_HAPLOTYPE_index,self.LOW_POS,'HIGH_MODERATE')
        
        
        #**************************************************  Get the iter_list on HIGH  ******************************************************************************
        print('step 4.4: Starting to Create the pegn on HIGH')
        unique_hap_count = len(self.df_HIGH['HAPLOTYPE'].unique())
        iter_list = interval_generator(unique_hap_count)
        
        for index in iter_list:
            # if file exists, don't create again
            if os.path.exists(os.path.join(self.output_directory, '{}_{}_HIGH.pgen'.format(str(index[0]),str(index[1])))):
                print("Skipping processing for {}_HIGH, file already exists.".format(str(index)))
            else:
                print('start processing {}_HIGH'.format(str(index)))
                pgen_generator(index,self.df_HIGH,self.output_directory,self.vcf_dict,self.df_HIGH_HAPLOTYPE_index,self.LOW_POS,'HIGH')

        
    def merge(self):
        """
        2. transform pgen into bed files
        3. Do the merging
        4. Transform the bed back to pgen file
        """
        #************************************************   Create the merge_list   **********************************************************************************
        print('step 5: Starting to Create the merge_list for all impact')
        index_list = interval_generator(len(self.df_new['HAPLOTYPE'].unique()))
        index_list = ['{}_{}_ALL'.format(str(interval[0]), str(interval[1])) for interval in index_list]
        my_string = ''
        for interval in index_list:
            my_string += '{}/{}\n'.format(self.output_directory,interval)
        
        with open(os.path.join(self.output_directory,'merge.txt'), 'w') as file:
            file.write(my_string)
            
        #************************************************   Create the merge_list on HIGH_MODERATE_LOW  **************************************************************
        print('step 5: Starting to Create the merge_list on HIGH_MODERATE_LOW')
        index_list_HIGH_MODERATE_LOW = interval_generator(len(self.df_HIGH_MODERATE_LOW['HAPLOTYPE'].unique()))
        index_list_HIGH_MODERATE_LOW = ['{}_{}_HIGH_MODERATE_LOW'.format(str(interval[0]), str(interval[1])) for interval in index_list_HIGH_MODERATE_LOW]
        my_string = ''
        for interval in index_list_HIGH_MODERATE_LOW:
            my_string += '{}/{}\n'.format(self.output_directory,interval)
        
        with open(os.path.join(self.output_directory,'merge_HIGH_MODERATE_LOW.txt'), 'w') as file:
            file.write(my_string)
            
        #*************************************************   Create the merge_list on HIGH_MODERATE   *****************************************************************
        print('step 5: Starting to Create the merge_list on HIGH_MODERATE')
        index_list_HIGH_MODERATE = interval_generator(len(self.df_HIGH_MODERATE['HAPLOTYPE'].unique()))
        index_list_HIGH_MODERATE = ['{}_{}_HIGH_MODERATE'.format(str(interval[0]), str(interval[1])) for interval in index_list_HIGH_MODERATE]
        my_string = ''
        for interval in index_list_HIGH_MODERATE:
            my_string += '{}/{}\n'.format(self.output_directory,interval)
        
        with open(os.path.join(self.output_directory,'merge_HIGH_MODERATE.txt'), 'w') as file:
            file.write(my_string)
        
        #***************************************************   Create the merge_list on HIGH   ************************************************************************
        print('step 5: Starting to Create the merge_list on HIGH')
        index_list_HIGH = interval_generator(len(self.df_HIGH['HAPLOTYPE'].unique()))
        index_list_HIGH = ['{}_{}_HIGH'.format(str(interval[0]), str(interval[1])) for interval in index_list_HIGH]
        my_string = ''
        for interval in index_list_HIGH:
            my_string += '{}/{}\n'.format(self.output_directory,interval)
        
        with open(os.path.join(self.output_directory,'merge_HIGH.txt'), 'w') as file:
            file.write(my_string)
            
        #****************************************************  Create the pgen name list   ***************************************************************************   
        #pgen_list_l = []
        pgen_list_l = [index_list, index_list_HIGH_MODERATE_LOW, index_list_HIGH_MODERATE, index_list_HIGH]
        #for list_ in index_list_l:
        #    pgen_list = []
        #    for interval in list_:
        #        pgen_list.append('{}_{}'.format(str(interval[0]), str(interval[1])))
        #    pgen_list_l.append(pgen_list)
            
        print(pgen_list_l)
            
        #******************************************************  Transform pgen into bed   ****************************************************************************
        for pgen_list in pgen_list_l:
            for pgen in pgen_list:
                command = "/frazer01/home/ximei/plink2 --pfile {} --make-bed --out {}".format(os.path.join(self.output_directory , pgen), os.path.join(self.output_directory , pgen))
                print('Starting step 5.1: transforming {}.pgen into bed files'.format(pgen))
                try:
                    subprocess.run(command, shell=True, check=True)
                    print("{}.pgen transformed into bed files successfully".format(pgen))
                except subprocess.CalledProcessError as e:
                    print("Error occurred in transforming {}.pgen into bed files".format(pgen))
        
        
        #********************************************************  Merge the bed files  ********************************************************************************
        merge_list_l = ['merge.txt', 'merge_HIGH_MODERATE_LOW.txt', 'merge_HIGH_MODERATE.txt', 'merge_HIGH.txt']
        for x in merge_list_l:
            command_2 = "/frazer01/software/plink-1.90b3x/plink --merge-list {} --make-bed --out {}".format(os.path.join(self.output_directory, x), os.path.join(self.output_directory , x.split('.')[0]))
            if not os.path.exists(self.output_directory + '{}.bed'.format(x.split('.')[0])):
                print('Starting step 5.2:: merging {}'.format(x.split('.')[0]))
                try:
                    subprocess.run(command_2, shell=True, check=True)
                    print("Step 5.2: completed successfully")
                except subprocess.CalledProcessError as e:
                    print("Error occurred in step 5.2:: merging on {}".format(x.split('.')[0]))
                                  
        
        
        #************************************************   Transform merge.bed back into merge.pgen   *****************************************************************
        transform_l = ['merge', 'merge_HIGH_MODERATE_LOW', 'merge_HIGH_MODERATE', 'merge_HIGH']  
        
        for x in transform_l:
            command_3 = "/frazer01/home/ximei/plink2 --bfile  {} --make-pgen --out {}".format(os.path.join(self.output_directory,x),os.path.join(self.final_output,x))
            print('Starting step 5.3:: Transform the bed on {} back to pgen file'.format(x))
            try:
                subprocess.run(command_3, shell=True, check=True)
                print("Step 5.3: completed successfully")
            except subprocess.CalledProcessError as e:
                    print("Error occurred in step 5.3:: Transform the bed on {} back to pgen file".format(x))
           
        

            
    def run(self):
        """
        Run GenerateVCF --> readVCF --> Haplotypes --> PrepareVCF --> Generate_Pgen --> merge
        """
        # Step 1
        self.GenerateVCF()

        # Step 2
        self.readVCF()
        
        # Step 3
        self.Haplotypes()

        # Step 5
        self.PrepareVCF()
        
        # Step 6
        self.Generate_Pgen()
        
        # Step 7
        self.merge()
        

# Run the tool in the terminal
if __name__ == "__main__":
    import argparse
    
    try:
        parser = argparse.ArgumentParser(description="HLA Analyzer Tool")
        parser.add_argument("gene_name", help="gene that we want to test")
        parser.add_argument("geneIndex_File", help="gene_info.txt")
        parser.add_argument("output_directory", help="Desired output folder")
        args = parser.parse_args()

        HLA_VCF = HLA_VCF_Pipeline(args.gene_name, args.geneIndex_File, args.output_directory)
        HLA_VCF.run()
        
    except subprocess.CalledProcessError as e:
        print('Error occurred: {}'.format(e))
    except Exception as ex:
        print('Other error occurred: {}'.format(ex))
        