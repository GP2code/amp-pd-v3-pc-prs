# Input: plink1.9 pca results of the study merged with the 1000 genome project
# Output1: Ancestry labeling determined by +/- 6SD of PC1-5 of the reference panel
# Output2: PCs with the 1000 genome project
# Output3: Plots for the whole cohort and the European cohort. 

# PCA results visualization and population inference
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
import argparse
import requests

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="calculate pca from plink file")

# Define the command-line arguments
parser.add_argument("--input", type=str, help="Description of input pca file")

# Parse the command-line arguments
args = parser.parse_args()
input_file = args.input

# Scree plot
dfpcaVal = pd.read_csv(f'{input_file}.eigenval', header=None)
dfpcaVal.plot.line()
plt.title('Scree plot')
# plt.show()
plt.savefig('scree_plot.png')





# PCA plot
# Read the population labeling
# The URL of the reference label file
url = "https://raw.githubusercontent.com/hirotaka-i/1kg_ref/main/all_hg38_filtered_chrpos_pop.txt"
# Send an HTTP GET request to the URL
response = requests.get(url)
local_file_path='all_hg38_filtered_chrpos_pop.txt'
# Check if the request was successful (status code 200)
if response.status_code == 200:
    # Open the local file in binary write mode and write the content of the response
    with open(local_file_path, 'wb') as file:
        file.write(response.content)
    print(f"File '{local_file_path}' downloaded successfully.")
else:
    print(f"Failed to download file. Status code: {response.status_code}")
# Read the file into a pandas DataFrame
t=pd.read_csv(local_file_path, sep='\t')


# merge it with eigenvec
d=pd.read_csv(f'{input_file}.eigenvec', delim_whitespace=True, 
               names=['FID', 'IID'] + [f'PC{i}' for i in range(1,21)])
df=pd.merge(d, t, on=['IID'], how='left')

df['Population']=df['Population'].fillna('Study')
df['Population2']=df['Population2'].fillna('Study')

colors = ['red', 'green', 'blue', 'cyan', 'magenta', 
        'black', 'purple', 'orange', 'pink', 'brown','yellow',
        'gray', 'navy', 'maroon', 'violet', 'turquoise', 'lime',
        'teal', 'indigo', 'coral', 'gold', 'darkred', 'darkgreen',
        'darkblue', 'lightgray', 'darkgray', 'beige', 'lightgreen', 'lightblue']
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
legend_lines = []

for i, (j, group) in enumerate(df.groupby('Population')):
    if j == 'Study':
        ax.scatter(x='PC1', y='PC2', color=colors[i], label=j, s=10, alpha=1, data=group)
    else:
        sns.kdeplot(x='PC1', y='PC2', n_levels=4, ax=ax, color=colors[i], alpha=0.8, data=group)
        legend_lines.append(mlines.Line2D([], [], color=colors[i], label=j))

# Add custom legend
ax.legend(handles=legend_lines, title='Population', loc='upper right')
plt.title('Whole cohort')
# plt.show()


# Save the figure
fig.savefig('PopulationPlot_w_1kg.png', dpi=300, format='png')



# generate a figure of the population pc plot
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
legend_lines = []

for i, (j, group) in enumerate(df.groupby('Population2')):
    if j == 'Study':
        ax.scatter(x='PC1', y='PC2', color=colors[i], label=j, s=10, alpha=1, data=group)
    elif j in ['AFR', 'EUR', 'EAS']:
        sns.kdeplot(x='PC1', y='PC2', n_levels=4, ax=ax, color=colors[i], alpha=0.8, data=group)
        legend_lines.append(mlines.Line2D([], [], color=colors[i], label=j))
    # else:
    #     print(j, 'not shown')

# Add custom legend
ax.legend(handles=legend_lines, title='Population', loc='upper right')
plt.title('Whole cohort with reference population of AFR, EUR, EAS')
# plt.show()

fig.savefig('PopulationPlot_w_1kg_v2.png', dpi=300, format='png')


# generate the population label based on 5 PCs and save
# Infer ancestry
dfpop = df.pivot_table(index='Population2', values=['PC1', 'PC2', 'PC3', 'PC4', 'PC5'], aggfunc=['mean', 'std'])

# Get the threshold table of mean +/- 6SD
def funcThres(x):
    lwl = x['mean'] - 6 * x['std']
    hgl = x['mean'] + 6 * x['std']
    return pd.Series({'lwl':lwl, 'hgl':hgl})
thres = dfpop.apply(funcThres, axis=1)

# function to infer ancestry for "OTH"
def funcInfPop(x):
    if x.Population2 != 'Study':
        InfPop = 'REF'
    else:
        InfPop = 'OTHER'
        for Population2 in ['EUR', 'EAS', 'AFR']:
            if (thres.loc[Population2, 'lwl']['PC1'] < x.PC1) & \
              (x.PC1 < thres.loc[Population2, 'hgl']['PC1']) & \
              (thres.loc[Population2, 'lwl']['PC2'] < x.PC2) & \
              (x.PC2 < thres.loc[Population2, 'hgl']['PC2']) & \
              (thres.loc[Population2, 'lwl']['PC3'] < x.PC3) & \
              (x.PC3 < thres.loc[Population2, 'hgl']['PC3']) & \
              (thres.loc[Population2, 'lwl']['PC4'] < x.PC4) & \
              (x.PC4 < thres.loc[Population2, 'hgl']['PC4'])& \
              (thres.loc[Population2, 'lwl']['PC5'] < x.PC5) & \
              (x.PC5 < thres.loc[Population2, 'hgl']['PC5']):
                    InfPop = Population2
    return InfPop
df['InfPop'] = df.apply(funcInfPop, axis=1)

# Count for each population
df.InfPop.value_counts()

# save the population
df.loc[df.InfPop!='REF', 
       ['FID', 'IID', 'InfPop'] + [f'PC{i+1}' for i in range(10)]].to_csv('genetic_ancestry_all_pca.csv', index=False)
for continent in ['EUR', 'EAS', 'AFR', 'OTHER']:
    t = df.loc[df.InfPop==continent, ['FID', 'IID']]
    print(continent, t.shape[0])
    t.to_csv(f'genetic_ancestry_{continent}.txt', index=False, sep='\t')


# Eurepeans plot PC1 vs PC2
df_euro = df[(df.Population2=='EUR') | (df.InfPop=='EUR')].copy()
colors = ['green', 'orange',  'pink', 'yellow', 'black', 'purple', 'blue', 'red', 'brown', 'gray']
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
legend_lines = []
df_euro.Population.fillna('Study', inplace=True)
for i, (j, group) in enumerate(df_euro.groupby('Population')):
    if j == 'Study':
        ax.scatter(x='PC1', y='PC2', color=colors[i], label=j, s=10, alpha=1, data=group)
    else:
        sns.kdeplot(x='PC1', y='PC2', n_levels=4, ax=ax, color=colors[i], alpha=0.8, data=group)
        legend_lines.append(mlines.Line2D([], [], color=colors[i], label=j))

# Add custom legend
ax.legend(handles=legend_lines, title='Popluation', loc='upper right')
plt.title('European cohort')
# plt.show()


# Save the figure
fig.savefig('PopulationPlot_w_1kg_EUR.png', dpi=300, format='png')


# Eurepeans plot PC1 vs PC3
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
legend_lines = []
for i, (j, group) in enumerate(df_euro.groupby('Population')):
    if j == 'Study':
        ax.scatter(x='PC1', y='PC3', color=colors[i], label=j, s=10, alpha=1, data=group)
    else:
        sns.kdeplot(x='PC1', y='PC3', n_levels=4, ax=ax, color=colors[i], alpha=0.8, data=group)
        legend_lines.append(mlines.Line2D([], [], color=colors[i], label=j))

# Add custom legend
ax.legend(handles=legend_lines, title='Popluation', loc='upper right')
plt.title('European cohort')
# plt.show()


# Save the figure
fig.savefig('PopulationPlot_w_1kg_EUR2.png', dpi=300, format='png')
