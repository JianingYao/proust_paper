import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

#### IFDLPFC ####
df_scores = pd.read_csv('./tables/IFDLPFC_silhouette_scores.csv', index_col=0)
df_melted = df_scores.melt(var_name='Method', value_name='Silhouette Score')

plt.figure(figsize=(8, 6))
sns.boxplot(x='Method', y='Silhouette Score', data=df_melted, showfliers=False)
sns.stripplot(x='Method', y='Silhouette Score', data=df_melted, 
              color='black', size=5, jitter=True, dodge=True)

plt.xlabel('Method')
plt.ylabel('Silhouette Score')
plt.tight_layout()
plt.savefig('./plot/IFDLPFC_silhouette_scores_boxplot.png')
plt.close()



#### IFmouse ####
df_scores = pd.read_csv('./tables/IFmouse_full_silhouette_scores.csv', index_col=0)
df_melted = df_scores.melt(var_name='Method', value_name='Silhouette Score')

plt.figure(figsize=(8, 6))
sns.boxplot(x='Method', y='Silhouette Score', data=df_melted, showfliers=False)
sns.stripplot(x='Method', y='Silhouette Score', data=df_melted, 
              color='black', size=5, jitter=True, dodge=True)

plt.xlabel('Method')
plt.ylabel('Silhouette Score')
plt.tight_layout()
plt.savefig('./plot/IFmouse_full_silhouette_scores_boxplot.png')
plt.close()
