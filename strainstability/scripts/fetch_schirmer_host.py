import pandas as pd
import sys

acc = sys.argv[1]
df = pd.read_csv("/u/home/r/rwolff/strain_stability_revisions/strainstability/metadata/host_dates/Schirmer_dates.txt",index_col=0)
host = df.loc[acc]["subject"]
print(host)
