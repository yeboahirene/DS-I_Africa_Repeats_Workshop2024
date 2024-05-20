import csv
import pandas as pd
import numpy as np
import math
from tqdm import tqdm
import random
import statsmodels.api as sm
from sklearn import preprocessing
import statsmodels.stats.multitest
import sys


Pop = sys.argv[1]
chrom = sys.argv[2]
name_file = sys.argv[3]
GB_summary = sys.argv[4]

def getData(filename):
    with open(filename, "r") as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            yield row

def exonstrSLR(chrom,Pop,name_addr,gb_addr,minTs=10,mingt=3,minsPgt=0,search_range=100_000):
    """
    Pop: population 
    minTs: minimum samples required for regression
    mingt: minimum genotypes required for regression
    minsPgt: minimum samples requried for genotypes to be included in regression
    search_range: ranges search for associated STRs
    """

    #load covariates
    cov_df=pd.read_csv(f'/home/grader-dsi-africa/public/project2-association/dataset/covariates_all.csv',sep='\t',index_col='sample_id')
    exp_df=pd.read_csv(f'/home/grader-dsi-africa/public/project2-association/dataset/{Pop}_normalized_and_filtered_hg38_chr{chrom}.csv',sep='\t')
    gt_dir=gb_addr
    names = []
    with open (name_addr) as f:
        for line in f:
            names.append(line.strip())
    
    csv_gen = getData(gt_dir)
    row_count = 0
    p_df = pd.DataFrame(columns =[ "str-gene","str_end","motif","gene_name", "sample_n","GT_n","p_values","slope",'error',"shuffled_p","shuffled_slope","shuffled_error"])
    #the total length just for illustration purpose, can remove tqdm when running in large batches
    col_name = ['CHROM','POS','motif','END'] + names
    for row in csv_gen:
        gt_value=list(filter(lambda a: a!='',row[0].split('\t')))
        gt_df=pd.DataFrame([gt_value],columns=col_name)
        gt_df[gt_df.columns[4:]]=gt_df[gt_df.columns[4:]]. \
                applymap(lambda x: [int(x.split('|')[0]),int(x.split('|')[1])] if '|' in x else [None,None])
        gt_samples=gt_df.columns[4:][gt_df.iloc[0,4:].apply(lambda x: x !=[None,None])].to_list()
        gt_psi_samples=list(set(exp_df.columns[5:]) & set(gt_samples))
        #STR GENOTYPE FILTER
        gt_psi_phased=gt_df[gt_psi_samples].T.rename(columns={0:'GT'}).applymap(lambda x: sum(x))
        gt_sum=gt_psi_phased.groupby('GT').size()

        #get joint data
        gt_ab_3=gt_sum[gt_sum>=minsPgt].index.tolist()
        #check number of genotype and check how many samples remains
        if len(gt_ab_3)<mingt or gt_sum[gt_ab_3].sum()<minTs:
            row_count+=1
            continue
        
        gt_psi_filted_samples=set(gt_psi_phased[gt_psi_phased.GT.isin(gt_ab_3)].index) & set(cov_df.index)
        gt_psi_samples=list(gt_psi_filted_samples)

        gt_df=gt_df[gt_df.columns[0:4].to_list()+gt_psi_samples]
        joint_cov_df=cov_df.loc[gt_psi_samples]
        joint_exp_df=exp_df[exp_df.columns[0:5].tolist()+gt_psi_samples]
        joint_exp_df=joint_exp_df[joint_exp_df.chromosome.values == gt_df['CHROM'].values]

        paired_df=joint_exp_df[(joint_exp_df.start - search_range <= int(gt_df.POS)) \
                         & (joint_exp_df.end + search_range >= (int(gt_df.POS)))]

        if paired_df.empty:
            row_count+=1
            continue

        for index, curr_exon in paired_df.iterrows():
            #contat the PSI, genotype, peer and pc for current exon of all samples
            a=pd.concat([curr_exon[5:].astype('float64'),\
                       gt_df[gt_psi_samples].T.rename(columns={0:'GT'}).applymap(lambda x: sum(x)),\
                       joint_cov_df],\
                       axis=1).rename(columns={index:'exp'})
            a=a[~a.exp.isnull()]
            #standardization
            a_scaled = preprocessing.StandardScaler().fit_transform(a)
            y=a_scaled[:,0]
            x=a_scaled[:,1:]
            x=sm.add_constant(x)
            mod_ols  = sm.OLS(y,x)
            res_ols = mod_ols.fit()
            p_values=res_ols.pvalues[1]
            slope=res_ols.params[1]
            err=res_ols.bse[1]

            shuffled_y = random.sample(list(y),len(y))
            mod_ols_s = sm.OLS(shuffled_y,x)
            res_ols_s = mod_ols_s.fit()
            shuffled_p=res_ols_s.pvalues[1]
            slope_p=res_ols_s.params[1]
            err_p=res_ols_s.bse[1]

            p_df = p_df.append({"str-gene":list(gt_df.CHROM +'_'+ gt_df.POS.str.rstrip()+'-'+curr_exon.gene_id)[0],\
                                "str_end":gt_df.END.tolist()[0],\
                                "motif":gt_df.motif.tolist()[0],\
                                "gene_name":curr_exon.gene_name,"sample_n":len(a),"GT_n":len(gt_ab_3),"p_values":p_values,\
                                "slope":slope,"error":err,"shuffled_p":shuffled_p,"shuffled_slope":slope_p,\
                                "shuffled_error":err_p}, ignore_index=True)
            row_count+=1

    return p_df

#running regression
reg_results=exonstrSLR(chrom,Pop,name_file, GB_summary)
for index, row in reg_results.iterrows():
    print(f'P-value for association between the repeat copy number and {row["gene_name"]} expression: {row["p_values"]}')