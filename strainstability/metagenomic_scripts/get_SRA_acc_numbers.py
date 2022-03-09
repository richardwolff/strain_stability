import pandas as pd

hosts = ["am","an","ae","ao"]

df=pd.read_csv("metagenomic_data_files/Poyet_SRA_report.txt",index_col=0,sep="\t")
df["host"] = [s[0] for s in df["sample_alias"].str.split("-")]

df_mg = df.loc[df["scientific_name"] == "human gut metagenome"]
df_mg = df_mg.loc[df_mg["library_strategy"] == "WGS"]

df_mg = df_mg.loc[df_mg["host"].isin(hosts)]

df_mg.index = df_mg["run_accession"]

df_mg["sra_ftp"].to_csv("metagenomic_data_files/SRA_accession_download_numbers.txt",index=None)