import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np

dataset_samples = {
    'VSPG DLPFC': ['2720',  '6432', '6522', '8667'], 
    'H&E DLPFC': ['151507', '151508', '151509', '151510', '151669', '151670', '151671', '151672', '151673', '151674',
     '151675', '151676'],
    'VSPG mouse': ['Visium_CKp25_rep1',
                    'Visium_CKp25_rep2',
                    'Visium_CKp25_rep3',
                    'Visium_CKp25_rep4'], 
    'VSPG AD': ['VIFAD1_V10A27004_D1_Br3880',
                    'VIFAD2_V10A27106_B1_Br3854',
                    'VIFAD2_V10A27106_C1_Br3873',
                    'VIFAD2_V10A27106_D1_Br3880',
                    'VIFAD3_V10T31036_B1_Br3854',
                    'VIFAD3_V10T31036_C1_Br3873',
                    'VIFAD3_V10T31036_D1_Br3880']
}

results = []

for dataset_name, samples in dataset_samples.items():
    for sample in samples:
        file_path = f'/users/jyao/projects/proust/Analysis/computation_time_memory_usage/table/memory_time_usage_{sample}.csv'
        df = pd.read_csv(file_path)
        df['Sample'] = sample 
        df['Dataset'] = dataset_name  
        results.append(df)

df_combined = pd.concat(results, ignore_index=True)
df_combined = df_combined[df_combined['Stage'].isin(['CNN Image Processing', 'GNN Gene Processing', 'GNN Image Processing'])]

######################## Computational time by study: four boxplots ########################
fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharey=True)
axes = axes.flatten()

for idx, (dataset_name, ax) in enumerate(zip(dataset_samples.keys(), axes)):
    sns.boxplot(
        data=df_combined[df_combined['Dataset'] == dataset_name],
        x='Stage',
        y='Time (s)',
        ax=ax,
        showfliers=False  
    )
    sns.stripplot(
        data=df_combined[df_combined['Dataset'] == dataset_name],
        x='Stage',
        y='Time (s)',
        ax=ax,
        jitter=True,  
        color='black',  
        size=5,  
        alpha=1.0 
    )
    ax.set_title(dataset_name)
    ax.set_xlabel('Stage')
    ax.set_ylabel('Time (s)')

plot_directory = './plot'
plot_path = os.path.join(plot_directory, 'computation_time_four_boxplots.png')
plt.tight_layout()
plt.savefig(plot_path)
plt.close()

######################## Computation time by study: one boxplot ########################
custom_palette = {
    'CNN Image Processing': '#66c2a5',  
    'GNN Gene Processing': '#fc8d62',    
    'GNN Image Processing': '#8da0cb'    
}

plt.figure(figsize=(8, 6))
sns.boxplot(
    data=df_combined,
    x='Stage',
    y='Time (s)',
    color='grey',
    showfliers=False  
)
sns.stripplot(
    data=df_combined,
    x='Stage',
    y='Time (s)',
    hue='Dataset', 
    jitter=True,  
    size=5,  
    alpha=0.7,  
    dodge=True
)

plt.legend(title='Dataset', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.title('Processing Time by Stage')
plt.xlabel('Stage')
plt.ylabel('Time (s)')

plot_path = os.path.join(plot_directory, 'computation_time_one_boxplot.png')
plt.tight_layout()
plt.savefig(plot_path)
plt.close()

