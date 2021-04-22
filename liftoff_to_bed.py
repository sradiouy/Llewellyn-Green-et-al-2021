## transform liftoff to bedfile

import pandas as pd 

dict_tolerant_sensitive = {"GIM-012":"TOLERANT","MUN-008":"TOLERANT","MUN-020":"TOLERANT","COR-018":"SENSITIVE","JUT-008":"SENSITIVE","AKA-018":"SENSITIVE"}


input_file = "mapped_features_aka-018.gtf"




df = pd.read_table(input_file,sep="\t",header=None,names=["chr","source","feature","start","end","frame","strand","score","attribute"])


rows = []
for index,row in df.iterrows():
    if row.feature == "gene":
        chromosome = row.chr
        start = row.start
        end = row.end
        gene_id = row.attribute.split(";")[0].split("=")[1]
        gene_symbol = row.attribute.split(";")[1].split("=")[1]
        coverage = row.attribute.split(";")[2].split("=")[1]
        sequence_ID = row.attribute.split(";")[3].split("=")[1]
        extra_copy_number = row.attribute.split(";")[4].split("=")[1]
        rows.append([chromosome,start,end,gene_id,gene_symbol,coverage,sequence_ID,extra_copy_number])


df_info = pd.DataFrame(rows,columns=["chromosome","start","end","gene_id","gene_symbol","coverage","sequence_ID","extra_copy_number"])


# load deg information

deg_file = "DEG_Cu_Sen.txt"


df_deg = pd.read_table(deg_file,sep="\t",header=0)

df_deg = df_deg[["Flybase ID","log2FoldChange","padj"]]
df_deg.columns= ["gene_id","log2FoldChange","padj"]


df_info = df_info.merge(df_deg,on="gene_id",how="left").fillna("-")


df_info.to_csv("AKA-018_genes.tsv",sep="\t",header=True,index=False)


