import pandas as pd
from strainstability import config

cd = pd.read_pickle("metagenomic_data_files/Poyet_collection_dates.pkl")
df = pd.read_csv("metagenomic_data_files/Poyet_SRA_report.txt",sep="\t",index_col=3)

hosts = config.hosts

df=pd.read_csv("metagenomic_data_files/Poyet_SRA_report.txt",index_col=0,sep="\t")
df["host"] = [s[0] for s in df["sample_alias"].str.split("-")]

df_mg = df.loc[df["scientific_name"] == "human gut metagenome"]
df_mg = df_mg.loc[df_mg["library_strategy"] == "WGS"]

df_mg = df_mg.loc[df_mg["host"].isin(hosts)]

df_mg.index = df_mg["run_accession"]

df_mg = df_mg[["host"]]

df_mg["collection_date"] = cd

df_mg_hosts = df_mg.groupby("host")

repeat_dic = {}
for host in hosts:

    df_am = df_mg_hosts.get_group(host)
    
    repeat_dic[host] = df_am.loc[df_am["collection_date"].duplicated(keep=False)].groupby("collection_date").groups
    
pd.Series(repeat_dic).to_pickle("metagenomic_data_files/Poyet_repeat_day_samples.pkl")