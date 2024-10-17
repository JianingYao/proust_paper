import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np

results = []

for i in range(1, 7):
    file_path = f'/users/jyao/projects/proust/Analysis/computation_time_memory_usage/table/memory_time_usage_{i}channels.csv'
    df = pd.read_csv(file_path)
    df['Channel'] = i 
    results.append(df)

df_combined = pd.concat(results, ignore_index=True)
df_combined = df_combined[df_combined['Stage'].isin(['CNN Image Processing', 'GNN Gene Processing', 'GNN Image Processing'])]

plot_directory = './plot'
os.makedirs(plot_directory, exist_ok=True)

# Create the plot
fig, ax = plt.subplots(1, 2, figsize=(14, 6), sharex=True)

# Plot for Time (s)
sns.lineplot(
    data=df_combined, 
    x='Channel', 
    y='Time (s)', 
    hue='Stage', 
    ax=ax[0],
    marker='o'
)
ax[0].set_title('Computation Time')
ax[0].set_xlabel('Number of Channels')
ax[0].set_ylabel('Time (s)')

# Plot for GPU Memory (MB)
sns.lineplot(
    data=df_combined, 
    x='Channel', 
    y='GPU Memory (MB)', 
    hue='Stage', 
    ax=ax[1],
    marker='o'
)
ax[1].set_title('Peak GPU Memory Usage')
ax[1].set_xlabel('Number of Channels')
ax[1].set_ylabel('GPU Memory (MB)')

# Create a separate legend and place it outside the plot
handles, labels = ax[0].get_legend_handles_labels()
ax[0].get_legend().remove()
ax[1].get_legend().remove()
legend = fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), title='Stage')

# Adjust layout to ensure the legend fits
# plt.tight_layout(rect=[0, 0, 0.9, 1])  # Leave space on the right for the legend
plt.tight_layout()

# Save the plot and show it for verification
plot_path = os.path.join(plot_directory, 'channel_line_chart.png')
plt.savefig(plot_path, bbox_extra_artists=(legend,), bbox_inches='tight')
plt.show()
plt.close()