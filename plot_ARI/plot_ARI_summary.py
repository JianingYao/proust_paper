import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

IF_ari = pd.read_csv('IF.csv')
df_melted = pd.melt(IF_ari, id_vars=['Sample'], var_name='Column')

# Create the bar plot
ax = sns.barplot(x='Column', y='value', hue='Sample', data=df_melted)
ax.set_xlabel('Methods')
ax.set_ylabel('ARI')
plt.savefig('./plot/IF-ari.png', dpi=600)
plt.close()


histology_ari = pd.read_csv('histology.csv')
df_melted = pd.melt(histology_ari, id_vars=['Sample'], var_name='Methods', value_name='ARI')

# Create the box plot
sns.boxplot(x='Methods', y='ARI', data=df_melted, width=0.5)
sns.stripplot(x='Methods', y='ARI', data=df_melted, color='black')
plt.savefig('./plot/histology-ari.png', dpi=600)
plt.close()