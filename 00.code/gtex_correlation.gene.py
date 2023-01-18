#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 3.8 #############################################################
# Libraries ####################################################################

import os
import pathlib
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import fdrcorrection

################################################################################
# Description/Notes ############################################################
################################################################################
"""
"""

################################################################################
# Base-level Functions #########################################################
################################################################################


################################################################################
# Task-specific Functions ######################################################
################################################################################


################################################################################
# Initiating Variables #########################################################
################################################################################

# name - raw data dn
raw_data_dn = "01.raw_data"
# name - data dn
data_dn = "02.data"
# name - output dir
output_dn = "03.output"

# name - gtex data dn
gtex_data_dn = os.path.join(raw_data_dn, "V8")

# name - gtex data fn
gtex_data_fn = os.path.join(gtex_data_dn, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")

# name - attributes of sample data filename
sample_attributes_data_fn = os.path.join(gtex_data_dn, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

# name - correlations dir
output_correlations_dn = os.path.join(output_dn, "correlations")

# create - correlations dir
pathlib.Path(output_correlations_dn).mkdir(parents=True, exist_ok=True)

# load - attributes of sample data
sample_attributes_data_df = pd.read_csv(sample_attributes_data_fn, sep="\t")

# name - target genes fn
target_genes_fn = "target_genes.txt"

# load - target genes
target_genes_df = pd.read_csv(target_genes_fn, sep="\t", header=None, names=["target_genes"])
target_genes_ls = target_genes_df['target_genes'].tolist() 
target_genes_ls = [x.upper() for x in target_genes_ls]

# set - target vals
fdr = 0.05
                            
# set - GTEx data headers dict
tissue_header_d = {"toptissue":"SMTS", "subtissue":"SMTSD"}

# set - spearman correlation threshold
spearman_correlation_cutoff = 0.25
################################################################################
# Execution ####################################################################
################################################################################
for target_gene_symbol in target_genes_ls:
    for target_tissue_level in tissue_header_d.keys():
        print(target_gene_symbol, target_tissue_level)

        # set all results fn and initiate df to hold results from each tissue
        all_results_fn = os.path.join(output_correlations_dn, ".".join(["00.correlations_combined", target_tissue_level, target_gene_symbol, "tsv"]))

        # test - all tissues fn exists
        if not os.path.exists(all_results_fn):

            # initiate - dataframe to capture tissue correlation results; of this tissue level
            all_results_df = pd.DataFrame([], columns = ["target_tissue", "tissue_n", "target_gene_symbol", "target_gene_id", "correlated_gene_symbol", "correlated_gene_id", "correlation", "pval", "target_gene_expr_avg", "correlated_gene_expr_avg"])

            # set - sample attribute tissue column, based on top (SMTS) or sub-level (SMTSD)
            tissue_col = tissue_header_d[target_tissue_level]

            # name - unique tissues at this tissue level
            target_tissues = sorted(sample_attributes_data_df[tissue_col].unique().tolist())

            # iterate - all tissues in this target tissue level
            for target_tissue in target_tissues:
                print(target_tissue)

                # set - filename for correlation results output
                target_tissue_fn = os.path.join(output_correlations_dn, ".".join([target_tissue, target_gene_symbol, "correlation.tsv"]))

                # test - tissue fn exists
                if not os.path.exists(target_tissue_fn):

                    # select - sample ids corresponding to this tissue                
                    target_tissue_sample_ids = sample_attributes_data_df[(sample_attributes_data_df[tissue_col]==target_tissue) & (sample_attributes_data_df['SMAFRZE']=="RNASEQ")]['SAMPID'].tolist()

                    # test - number of samples is above threshold
                    if len(target_tissue_sample_ids) >= 10:
                        
                        # load - expression data restricted to these sample ids (tissue)
                        target_tissue_df = pd.read_csv(gtex_data_fn, sep="\t", header=2, usecols=['Name', 'Description'] + target_tissue_sample_ids)

                        # select - expression of target gene within this tissue (by gene name)
                        target_gene_expr_df = target_tissue_df[target_tissue_df['Description']==target_gene_symbol]

                        if len(target_gene_expr_df) == 1:

                            # name - Ensembl ID of target gene
                            target_gene_id = target_gene_expr_df.iloc[0]['Name']
                            
                            # select - expression values of target gene
                            target_gene_expr_ls = target_gene_expr_df.values.flatten().tolist()[2:]

                            # calculate - average expression of target gene in tissue
                            target_gene_expr_avg = round(sum(target_gene_expr_ls)/len(target_gene_expr_ls), 3)

                            # set - list to capture correlations
                            tissue_results_ls = []
                            count = 0
                            
                            # iterate - all genes expression rows
                            for r in zip(*target_tissue_df.to_dict("list").values()):
                                other_gene_symbol, other_gene_id = r[:2]

                                # select - expression values of other gene in tissue
                                other_gene_expr_ls = r[2:]

                                # calculate - average expression of other gene in tissue
                                other_gene_expr_avg = round(sum(other_gene_expr_ls)/len(other_gene_expr_ls), 3)

                                # calculate - spearman correlation
                                spearman_correlation, pval = spearmanr(target_gene_expr_ls, other_gene_expr_ls)
                                tissue_results_ls.append([target_tissue, len(target_tissue_sample_ids), target_gene_symbol, target_gene_id, other_gene_symbol, other_gene_id, round(spearman_correlation, 4), pval, target_gene_expr_avg, other_gene_expr_avg])

                                # track - number of correlations calculated (target gene - other genes; in target tissue)
                                count+=1
                                if count%5000==0:
                                    print(count)

                            # initiate - df of correlation results for target tissue
                            tissue_results_df = pd.DataFrame(tissue_results_ls, columns = ["target_tissue", "tissue_n", "target_gene_symbol", "target_gene_id", "correlated_gene_symbol", "correlated_gene_id", "correlation", "pval", "target_gene_expr_avg", "correlated_gene_expr_avg"])

                            # filter - entries where correlation/pval were not determined
                            tissue_results_df = tissue_results_df[tissue_results_df['pval'].notnull()]
                            
                            # sort - correlation values, for this gene-tissue correlation analysis
                            tissue_results_df.sort_values(by='correlation', ascending=False, inplace=True)

                            # calculate - correct for multiple testing by tissue
                            rejected, corrected_pvals = fdrcorrection(tissue_results_df['pval'].tolist(), is_sorted=False, alpha=fdr)

                            # add - corrected pvals
                            tissue_results_df['pval corrected'] = corrected_pvals
                            
                            # filter - tissue results to strong correlations
                            tissue_results_df = tissue_results_df[abs(tissue_results_df['correlation'])>=spearman_correlation_cutoff]                      
                            
                            # export - tissue results df
                            tissue_results_df.to_csv(target_tissue_fn, index=False, sep="\t")

                            # add - tissue results to all results df
                            all_results_df = pd.concat([all_results_df, tissue_results_df])

                elif os.path.exists(target_tissue_fn):
                    # load tissue results df and add to all results df
                    tissue_results_df = pd.read_csv(target_tissue_fn, sep="\t")
                    all_results_df = pd.concat([all_results_df, tissue_results_df])

            # sort by correlation across all tissues all results df
            all_results_df.sort_values("correlation", ascending=False, inplace=True)
            all_results_df.to_csv(all_results_fn, sep="\t", index=False)



