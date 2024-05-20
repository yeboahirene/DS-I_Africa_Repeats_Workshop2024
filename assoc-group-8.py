import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns


dfafr= pd.read_csv('/Users/haroun/Downloads/AFR_expr_results.tab',sep="\t")

dfeur= pd.read_csv('/Users/haroun/Downloads/EUR_expr_results.tab',sep="\t")

reject, qvaluesafr, _, _ = multipletests(dfafr['p_values'], method='bonferroni')#Correcting p-values

dfafr['adj.p_values']=qvaluesafr

reject, qvalueseur, _, _ = multipletests(dfeur['p_values'], method='bonferroni')#extract adj.p_values
dfeur['adj.p_values']=qvalueseur
dfafr=dfafr[dfafr['p_values']<=0.05]
dfeur=dfeur[dfeur['p_values']<=0.05]
dfafrf=dfafr[dfafr['adj.p_values']<=0.05]

dfeurf=dfeur[dfeur['adj.p_values']<=0.05]

#print(dfafrf)# After p values adustment for AFR population
#print(dfafrf[['str-gene','motif','gene_name','sample_n','GT_n','q_values']])
#!pip3 install --user dataframe_image
import dataframe_image as dfi

dfi.export(dfafrf[['str-gene','motif',  'gene_name',   'sample_n','GT_n',  'adj.p_values']], 'AFRoutput.png',table_conversion="matplotlib")#To save data summary
dfi.export(dfeurf[['str-gene',  'motif',  'gene_name',  'sample_n','GT_n',  'adj.p_values']], 'EURoutput.png',table_conversion="matplotlib")
#chr11_57528484-ENSG00000134809.4 strs:#CHROM  chr11
#POS     57528484
print(dfeurf)# After p values adustment for EUR population

#chr11_57520862-ENSG00000134809.4 str:AAAAAAAAAAAAAAAAA Human_STR_231162
##CHROM  chr11
#POS     57524347
#ID      Human_STR_231165
#REF     CTTTTTTTTTTTTTTTTTTTAATTG
#ALT     (T)CTTTTTTTTAATCG,CTTTTTTTTTTTTTTTTTTAATTG,CATTTTTTTTTTTTTTTTTTAATTG,CTTTTTTTTTTTTTTTTTTTAATCG,CTTTTTTTTTTTTTTTTTTTAATTT,CTTTTTTTTTTTTTTTTTTTTAATTG
#POS     57528484ID      Human_STR_231168
#REF     TATATATATATATATATATATATATATATATATA
#ALT     (TA)TATATATATATATATATATATATA,TATATATATATATATATATATACATA,TATATATATATATATATATATATACA,TATATATATATATATATATATATATA,TATATATATATATATATATATATATATA,TATATATATATATATATATATATATATATA,TATATATATATATATATATATATATATATATA,TATATATATATATATATATATATATATATATATATA,TATATATATATATATATATATATATATATATATATATA,TATATATATATATATATATATATATATATATATATATATA,TATATATATATATATATATATATATATATATATATATATATA
# Create bar plot
#dfeurf['q_values']=-np.log10(dfeurf['q_values'])
dfeurf['STR/Gene']=['A rpts-TIMM10','A rpts-SMTNL1','T rpts-TIMM10','T rpts-SMTNL1','AT rpts-TIMM10','AT rpts-SMTNL1']
dfafrf['STR/Gene']=['AT rpts-TIMM10']

"""
Association plot pvalues
"""
fig=px.bar(dfafrf, x="adj.p_values", y='STR/Gene', orientation='h',title='pop')
fig.update_layout(xaxis=dict(tickfont=dict(size=12)),  # Adjust size as needed
                  yaxis=dict(tickfont=dict(size=12))) 
fig.update_layout(title_text='Association between STRs and gene expression for AFR  population', title_x=0.3)
fig.show()
"""
Association plot slop
"""
fig=px.bar(dfafrf, x="slope", y='STR/Gene', orientation='h',title='pop')
fig.update_layout(xaxis=dict(tickfont=dict(size=12)),  # Adjust size as needed
                  yaxis=dict(tickfont=dict(size=12))) 
fig.update_layout(title_text='Association between STRs and gene expression for AFR  population', title_x=0.3)
fig.show()
