import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


df = pd.read_csv("Data_P64.txt", header=None, names=["time in log10 scale in seconds"])
df["cat"] = np.tile(np.repeat(["MPI_Send", "MPI_Pack", "MPI DDT"], 5), 7)
df["Data Points Per Process"] = np.repeat(["16^2", "32^2", "64^2", "128^2", "256^2", "512^2", "1024^2"], 15)
#df["Data Points Per Process"] = df["Data Points Per Process"].astype(str)

#df['time in log10 scale'] = np.log10(df['time in log10 scale'].values)


ax = sns.boxplot(data=df, x="Data Points Per Process", y="time in log10 scale in seconds", hue="cat", zorder=0.9, dodge=False)
ax.set_yscale('log')
#add the lines connecting the median values
sns.pointplot(data=df, x="Data Points Per Process", y="time in log10 scale in seconds", hue="cat", estimator=np.median, ci=None, markers="None", ax=ax)

#remove duplicate entries, setting the legend outside the graph, new title
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:3], labels[:3], bbox_to_anchor=(1, 1), title="parallel")
plt.tight_layout()
#plt.show()
plt.savefig("Plot for P=64")
