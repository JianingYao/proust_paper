import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np

dataset_samples = {
    'VSPG DLPFC': ['2720', '6432', '6522', '8667'], 
    'H&E DLPFC': ['151507', '151508', '151509', '151510', '151669', '151670', '151671', '151672', '151673', '151674',
                  '151675', '151676'],
    'VSPG mouse': ['Visium_CKp25_rep1', 'Visium_CKp25_rep2', 'Visium_CKp25_rep3', 'Visium_CKp25_rep4'], 
    'VSPG AD': ['VIFAD1_V10A27004_D1_Br3880', 'VIFAD2_V10A27106_B1_Br3854', 'VIFAD2_V10A27106_C1_Br3873',
                'VIFAD2_V10A27106_D1_Br3880', 'VIFAD3_V10T31036_B1_Br3854', 'VIFAD3_V10T31036_C1_Br3873',
                'VIFAD3_V10T31036_D1_Br3880']
}
proust_results = []
for dataset_name, samples in dataset_samples.items():
    for sample in samples:
        file_path = f'/users/jyao/projects/proust/Analysis/computation_time_memory_usage/table/memory_time_usage_{sample}.csv'
        df = pd.read_csv(file_path)
        df['Sample'] = sample
        df['Method'] = 'Proust'
        max_gpu_memory = df.loc[df['Stage'].isin(['CNN Image Processing', 'GNN Gene Processing', 'GNN Image Processing']),'GPU Memory (MB)'].max()
        df.loc[df['Stage'] == 'Total', 'GPU Memory (MB)'] = max_gpu_memory
        proust_results.append(df)

proust_df_combined = pd.concat(proust_results, ignore_index=True)
proust_df_combined = proust_df_combined[proust_df_combined['Stage'] == 'Total']
proust_df_combined = proust_df_combined.drop(columns=['Stage'])
cols = ['Sample'] + [col for col in proust_df_combined.columns if col != 'Sample']
proust_df_combined = proust_df_combined[cols]



datasets = ['IFDLPFC', 'DLPFC', 'IFmouse', 'IFAD']
graphst_results = []
for dataset in datasets:
    file_path = f'/users/jyao/projects/proust/Analysis/computation_time_memory_usage/table/memory_time_usage_graphST_{dataset}.csv'
    df = pd.read_csv(file_path)
    df['Method'] = 'GraphST'
    graphst_results.append(df)

graphst_df_combined = pd.concat(graphst_results, ignore_index=True)

result_df = pd.concat([proust_df_combined, graphst_df_combined], ignore_index=True)

### boxplot for Time and GPU Memory
fig, axes = plt.subplots(1, 2, figsize=(12, 6))
axes = axes.flatten()
palette = {"Proust": "#1f77b4", "GraphST": "#ff7f0e"}

# Time Boxplot
sns.boxplot(
    x="Method", y="Time (s)", data=result_df, ax=axes[0], palette=palette, showfliers=False
)
sns.stripplot(
    x="Method", y="Time (s)", data=result_df, ax=axes[0], 
    color='black', jitter=True, linewidth=0.5, size=4
)
axes[0].set_title('Total Time Usage')
axes[0].set_xlabel('Method')
axes[0].set_ylabel('Time (s)')

# GPU Memory Boxplot
sns.boxplot(
    x="Method", y="GPU Memory (MB)", data=result_df, ax=axes[1], palette=palette, showfliers=False
)
sns.stripplot(
    x="Method", y="GPU Memory (MB)", data=result_df, ax=axes[1], 
    color='black', jitter=True, linewidth=0.5, size=4
)
axes[1].set_title('Peak GPU Memory Comsumption')
axes[1].set_xlabel('Method')
axes[1].set_ylabel('GPU Memory (MB)')

# Adjust layout and overall title
plt.suptitle('Performance Comparison')

plot_directory = './plot'
plot_path = os.path.join(plot_directory, 'performance_comparison.png')
plt.tight_layout()
plt.savefig(plot_path)
plt.close()
