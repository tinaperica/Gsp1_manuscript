#!/usr/bin/env python

import pandas as pd

library_gene_names = list(pd.read_csv('../../Data/E-MAP/gsp1_pEMAP_avg_merged_gene_names.txt', sep='\t').columns)[1:]
library_gene_ORFs = list(pd.read_csv('../../Data/E-MAP/gsp1_pEMAP_avg_merged.txt', sep='\t').columns)[1:]
gene_names_to_ORF = pd.DataFrame({'name': library_gene_names,'ORF': library_gene_ORFs})

# remove ' - DAmP' and turn 'ORF - ORF' into 'ORF'
gene_names_to_ORF['name'], _ = gene_names_to_ORF['name'].str.split(' - ', 1).str
gene_names_to_ORF['ORF'], _ = gene_names_to_ORF['ORF'].str.split(' - ', 1).str

# get annotations from SGD, code generated from the following URL:
# "https://yeastmine.yeastgenome.org/yeastmine/results.do?trail=%257Cquery"

from intermine.webservice import Service
service = Service("https://yeastmine.yeastgenome.org:443/yeastmine/service")
query = service.new_query("Gene")
query.add_view("secondaryIdentifier", "symbol", "name", "length", "sgdAlias", "description")
query.add_constraint("status", "IS NULL", code = "D")
query.add_constraint("status", "=", "Active", code = "C")
query.add_constraint("dataSets.name", "=", "SGD data set", code = "F")
query.add_constraint("organism.name", "=", "Saccharomyces cerevisiae", code = "E")
query.set_logic("(C or D) and E and F")

# make a dataframe from the query results
SGD_annotations = pd.DataFrame(query.results('dict'))
cols = list(SGD_annotations.columns)
SGD_annotations.columns = [col.split('.')[1] for col in cols]
SGD_annotations = SGD_annotations.drop(['cytoLocation','featAttribute','geneSummary',
                       'id','length','primaryIdentifier','qualifier',
                       'score','scoreType','status'], axis=1)
SGD_to_merge = SGD_annotations[['symbol','secondaryIdentifier','description','name']]
SGD_to_merge.columns = ['name', 'ORF', 'description', 'name_meaning']
SGD_to_merge.to_csv('../../Data/E-MAP/SGD_descriptions_all.txt', sep='\t', na_rep='', index=False)
df = pd.merge(gene_names_to_ORF, SGD_to_merge, how='left', on=['name','ORF'])

# add in hand-curated descriptions copied from SGD that weren't in the yeastmine download
# NOTE: many unknown gene descriptions were left out
df = pd.merge(df, pd.read_csv('../../Data/E-MAP/SGD_missing_descriptions.txt', sep='\t'), how = 'outer')

# drop if description is empty (won't add any information)
SGD_descriptions = df.loc[~pd.isnull(df.description)]

# save dataframe
SGD_descriptions.to_csv('../../Data/E-MAP/SGD_descriptions_library_genes.txt', sep='\t', na_rep='', index=False)
