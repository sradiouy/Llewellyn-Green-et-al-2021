import pandas as pd 

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
