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
        strand = row.strand
        gene_id = row.attribute.split(";")[0].split("=")[1]
        gene_symbol = row.attribute.split(";")[1].split("=")[1]
        coverage = row.attribute.split(";")[2].split("=")[1]
        sequence_ID = row.attribute.split(";")[3].split("=")[1]
        extra_copy_number = row.attribute.split(";")[4].split("=")[1]
        rows.append([chromosome,start,end,strand,gene_id,gene_symbol,coverage,sequence_ID,extra_copy_number])


df_info = pd.DataFrame(rows,columns=["chromosome","start","end","strand","gene_id","gene_symbol","coverage","sequence_ID","extra_copy_number"])


# load deg information

deg_file = "DEG_Cu_Sen.txt"


df_deg = pd.read_table(deg_file,sep="\t",header=0)

df_deg = df_deg[["Flybase ID","log2FoldChange","padj"]]
df_deg.columns= ["gene_id","log2FoldChange","padj"]


df_info = df_info.merge(df_deg,on="gene_id",how="left").fillna("-")


df_info.to_csv("AKA-018_genes.tsv",sep="\t",header=True,index=False)

### Filtro cromosoma 4 y transformo en bed

df_info = df_info[df_info["chromosome"] != "4"]

df_info[["chromosome","start","end","strand","gene_id"]].to_csv("AKA-018_genes.bed",sep="\t",header=False,index=False)

####

te_strain = "AKA-018.GenesAllClosestGenes.bed"
df_te_strain = pd.read_table(te_strain,sep="\t",header=None)[[0,1,2,3]]
df_te_strain.columns = ["te_chr","te_start","te_end","te_name"]





te_info = "geo_location_by_genome_both.tsv"
df_te_info = pd.read_table(te_info,sep="\t",header=0)





te_ids_dict = {}
rows = []
for index,row in df_te_info.iterrows():
    for te_id in row.te_ids.split(";"):
        recomb = row.RecomConclusion
        te_ref = row.te_ref
        g_count = row.g_count
        Class = row.Class
        Order = row.Order
        SuperFamily = row.SuperFamily
        family = row.family
        genomes = row.genomes
        te_ids_dict[te_id] = [te_ref,g_count,Class,Order,SuperFamily,family,genomes,recomb]


rows = []
for index,row in df_te_strain.iterrows():
    recomb      = "-"
    te_ref      = "-"
    g_count     = "-"
    Class       = "-"
    Order       = "-"
    SuperFamily = "-"
    family      = "-"
    genomes     = "-"
    if row.te_name in te_ids_dict:
        rows.append([row.te_chr,row.te_start,row.te_end,row.te_name] + te_ids_dict[row.te_name])
    else:
        rows.append([row.te_chr,row.te_start,row.te_end,row.te_name,te_ref,g_count,Class,Order,SuperFamily,family,genomes,recomb])


df_te_final = pd.DataFrame(rows,columns=["te_chr","te_start","te_end","te_name","te_ref","g_count","Class","Order","SuperFamily","family","genomes","recombRate"])


df_te_final.to_csv("AKA-018_tes.tsv",sep="\t",header=True,index=False)
df_te_final[["te_chr","te_start","te_end","te_name"]].to_csv("AKA-018_tes.bed",sep="\t",header=False,index=False)



####
import pandas as pd 

genes_te = "AKA-018_genes_tes.tsv"


df_genes_te = pd.read_table(genes_te,sep="\t",header=None)
df_genes_te.columns = ["chromosome","start","end","strand","gene_id","te_chr","te_start","te_end","te_name","distance"]


rows = []
for name,grp in df_genes_te.groupby("gene_id"):
    association = 0
    for index,row in grp.iterrows():
        if row.distance > 0:
            if row.strand == "+":
                orientation = "downstream"
            else:
                orientation = "upstream"
        elif row.distance < 0:
            if row.strand == "+":
                orientation = "upstream"
            else:
                orientation = "downstream"
        else:
            orientation = "inside"
        if abs(row.distance) <= 1000:
            association += 1
            rows.append(row.tolist() + [orientation,association])
        if association == 0:
            rows.append(row.tolist() + [orientation,-1])
            break


df = pd.DataFrame(rows,columns=["chromosome","start","end","strand","gene_id","te_chr","te_start","te_end","te_name","distance","orientation","association"])

### Ahora defino la distancia y cuantos te se asocian a genes en el radio de 1000 


rows = []
for name,grp in df.groupby("gene_id"):
    if len(grp) == 1:
        if grp.association.iloc[0] == -1:
            rows.append(grp.iloc[0].tolist() + [0])
        else:
            rows.append(grp.iloc[0].tolist() + [1])
    else:
        for index,row in grp.iterrows():
            rows.append(row.tolist() + [len(grp)])

df = pd.DataFrame(rows,columns=["chromosome","start","end","strand","gene_id","te_chr","te_start","te_end","te_name","distance","orientation","association","nTEsinrange"])



# Agrego info de los genes

gene_info = "AKA-018_genes.tsv"
df_gene_info = pd.read_table(gene_info,sep="\t",header=0)

df = df.merge(df_gene_info,on="gene_id",how="left")


df = df[["chromosome_x","start_x","end_x","strand_x",'gene_id', 'gene_symbol', 'coverage', 'sequence_ID','extra_copy_number', 'te_name', 'distance', 'orientation','association', 'nTEsinrange','log2FoldChange', 'padj']]

df.columns = ["chromosome","start","end","strand",'gene_id', 'gene_symbol', 'coverage', 'sequence_ID','extra_copy_number', 'te_name', 'distance', 'orientation','association', 'nTEsinrange','log2FoldChange', 'padj']



# Agrego info de los TEs

te_info = "AKA-018_tes.tsv"
df_te_info = pd.read_table(te_info,sep="\t",header=0)

df = df.merge(df_te_info,on="te_name",how="inner")


df = df[["chromosome","start","end","strand",'gene_id', 'te_name','distance','association','nTEsinrange','gene_symbol', 'coverage', 'sequence_ID','extra_copy_number',   'orientation', 'log2FoldChange', 'padj','te_ref', 'g_count', 'Class', 'Order','SuperFamily', 'family', 'genomes', 'recombRate']]

df = df.sort_values(["chromosome","start","end"])

df.to_csv("AKA-018_Llew.tsv",sep="\t",header=True,index=False)
df.to_excel("AKA-018_Llew.xlsx",float_format="%.2f")

## Ahora filtro solo por los DEGs

df_degs = df[df.log2FoldChange != "-"]

df_degs.to_csv("AKA-018_DEGs_Llew.tsv",sep="\t",header=True,index=False)
df_degs.to_excel("AKA-018_DEGs_Llew.xlsx",float_format="%.2f")


