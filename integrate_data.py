"""read protein sequence from .csv and calculate the FCR value"""
import pandas as pd
import numpy as np
from localcider.sequenceParameters import SequenceParameters
data1 = pd.read_csv("/Users/zoey/Desktop/PS/PSdata/293tnpsp.csv")
data2 = pd.read_csv("/Users/zoey/Desktop/PS/PSdata/293tpsp.csv")
# data.rename(columns={'Unnamed: 0':'names','Amino_Acids_Number':'AA_number'}, inplace=True)
sequences1 = data1['Sequences']
sequences2 = data2['Sequences']
test10 = [seq for seq in list(sequences1)]
test2 = [seq for seq in list(sequences2)]

# del test1[704], test1[3031], test1[3988], test1[4322]
AAs = ['R','H','K','D','E','S','T','N','Q','C','G','P','A','I','L','M','F','W','Y','V']
def delete_seq(sequences):
    seq_out = []
    for seq in sequences:
        for c in seq: 
            if c not in AAs:
                seq_out.append(seq)
    return seq_out
seq_out = delete_seq(test10) 
test1 = [x for x in test10 if x not in seq_out]        
fcr1, hydro1 = np.zeros(len(test1)), np.zeros(len(test1))
fcr2, hydro2 = np.zeros(len(test2)), np.zeros(len(test2))
for i in range(len(test1)):
    seq = SequenceParameters(test1[i])
    fcr1[i] = seq.get_FCR()
    hydro1[i] = np.mean(seq.get_linear_hydropathy()[1,:])
dic = {'Sequences':test1, 'FCR':pd.Series(fcr1), 'Hydropathy':pd.Series(hydro1)}
data11 = pd.DataFrame(data=dic)
#data1['FCR'] = pd.Series(fcr1)
#data1['Hydropathy'] = pd.Series(hydro1)
data11['Names'] = np.nan
data11['AA_Number'] = np.nan
for i in range(len(data11['Sequences'])):
    name = data11['Sequences'][i]
    for j in range(len(data1['Sequences'])):
        if data1['Sequences'][j] == name:
            data11['Names'][i], data11['AA_Number'][i] = data1['Names'][j], data1['AA_Number'][j]
data11.to_csv("/Users/zoey/Desktop/PS/localCIDER/293tnonPS.csv")

for i in range(len(test2)):
    seq = SequenceParameters(test2[i])
    fcr2[i] = seq.get_FCR()
    hydro2[i] = np.mean(seq.get_linear_hydropathy()[1,:])
data2['FCR'] = pd.Series(fcr2)
data2['Hydropathy'] = pd.Series(hydro2)
data2.to_csv("/Users/zoey/Desktop/PS/localCIDER/293tPSP-FCR-Hydro.csv")

# add a column of PLAAC to this csv
data = pd.read_csv("/Users/zoey/Desktop/PS/localCIDER/293tPS.csv")
data1 = pd.read_csv("/Users/zoey/Desktop/PS/PLAAC/plaac_ps.tsv", sep = "\t")
data['PLAAC'] = np.nan
for i in range(0,len(data['Names'])):
    name = data['Names'][i]
    for j in range(0,len(data1['SEQid'])):
        if data1['SEQid'][j] == name:
            data['PLAAC'][i] = data1['NLLR'][j]
data.to_csv("/Users/zoey/Desktop/PS/output/293tPSP.csv")

# add column of IDR to .csv
# IDR score of a sequence is the proportion of disordered amino acids in a sequemce
import os
# directory = os.path.join("c:\\","path")
dic = {'names':[],'IDR':[]}
for root,dirs,files in os.walk("/Users/zoey/Desktop/PS/IDR-Espritz/PSbatch"):
    for file in files:
       if file.endswith(".espritz"):
           f=open(root + '/' + file, 'r')
           data = pd.read_csv(root + '/' + file, sep = "\t", header=None)
           isorder = list(data.iloc[:,0])
           count = 0
           for x in isorder:
            if x == 'D':
               count += 1
           score = count/len(isorder)
           name = os.path.basename(f.name)
           dic['names'].append(name)
           dic['IDR'].append(score)
           f.close()
dataIDR = pd.DataFrame(data=dic)
#namels = [x.replace('.fasta.espritz','') for x in list(dataIDR['names'])]
dataIDR['names'] = dataIDR['names'].str.replace('.fasta.espritz','')
for i in range(0,len(dataIDR['names'])):
    x = dataIDR['names'][i]
    x = x[:2] + '|' + x[2:]
    dataIDR['names'][i] = x
dataIDR.to_csv("/Users/zoey/Desktop/PS/IDR-Espritz/293tPS_IDR.csv")
## add
import operator as op
data = pd.read_csv("/Users/zoey/Desktop/PS/output/293tPSP.csv")
data_add = pd.read_csv("/Users/zoey/Desktop/PS/IDR-Espritz/293tPS_IDR.csv")
data['IDR'] = np.nan
for i in range(0,len(data['Names'])):
    name = data['Names'][i]
    for j in range(0,len(data_add['names'])):
        if (op.contains(name, data_add['names'][j])):
            data['IDR'][i] = data_add['IDR'][j]
data.to_csv("/Users/zoey/Desktop/PS/output/293tPSP-IDR.csv")

# add a PScore
data = []
with open("/Users/zoey/Desktop/PS/PScore/293tPS.out", 'r') as data_file:
        for line in data_file:
                data.append(line.split())
table = [[x[2], x[1]] for x in data]
df = pd.DataFrame(table)
df = df.rename(columns = {0:'Name', 1:'PScore'})
df['Name'] = df['Name'].str.replace('>','')
data0 = pd.read_csv("/Users/zoey/Desktop/PS/output/293tPSP-IDR.csv")
data0['PScore'] = np.nan
for i in range(0,len(data0['Names'])):
    name = data0['Names'][i]
    for j in range(0,len(df)):
        if (op.contains(name, df['Name'][j])):
            data0['PScore'][i] = df['PScore'][j]
data0.to_csv("/Users/zoey/Desktop/PS/output/293tPSP5.csv")

# catCIDER 
names, granule = [],[]
for root,dirs,files in os.walk("/Users/zoey/Desktop/PS/catGRANULE/nonPS"):
    for file in files:
       if file.endswith(".txt"):
           f=open(root + '/' + file, 'r')
           data = pd.read_csv(root + '/' + file, sep = " ", header=None)
           names = names + list(data.iloc[:,0])
           granule = granule + list(data.iloc[:,1])
data_cat = {'Names':names,'Score':granule}
df_cat = pd.DataFrame(data=data_cat)

df_cat = pd.read_csv("/Users/zoey/Desktop/PS/catGRANULE/PS.txt", sep = " ", header=None)
df_cat = df_cat.rename(columns = {0:'Names', 1:'Score'})
data = pd.read_csv("/Users/zoey/Desktop/PS/output/293tPSP5.csv")
data['catGRANULE'] = np.nan
for i in range(0,len(data['Names'])):
    name = data['Names'][i]
    for j in range(0,len(df_cat)):
        if (op.contains(name, df_cat['Names'][j])):
            data['catGRANULE'][i] = df_cat['Score'][j]
data.to_csv("/Users/zoey/Desktop/PS/output/293tPSP6.csv")

# phospho
data = pd.read_csv("/Users/zoey/Desktop/PS/output/293tnonPSP6.csv")
seq = data['Sequences']
org = ['Human'] * len(data['Sequences'])
dic = {'organism':org, 'sequence':seq}
df = pd.DataFrame(data=dic)
# df.to_csv("/Users/zoey/Desktop/PS/phospho/293tPS.csv")
with open('/Users/zoey/Desktop/PS/phospho/file2.txt', 'w') as f:
    for i in dic:
        f.write(i + "\t")
        numSubItems = len(dic[i])
    f.write("\n")

    for level in range(numSubItems):
        for i in dic:
            f.write(str(dic[i][level]) + "\t")
        f.write("\n")
# split file
sep = 1000
df1 = df.iloc[:sep,:]
df2 = df.iloc[sep:2*sep,:]
df3 = df.iloc[2*sep:3*sep,:]
df4 = df.iloc[3*sep:4*sep,:]
df5 = df.iloc[4*sep:5*sep,:]
df6 = df.iloc[5*sep:6*sep,:]
df7 = df.iloc[6*sep:7*sep,:]
df8 = df.iloc[7*sep:8*sep,:]
df9 = df.iloc[8*sep:9*sep,:]
df10 = df.iloc[9*sep:,:] 
#np.savetxt("/Users/zoey/Desktop/PS/phospho/file1.txt",df1.values,fmt='%s',header='organism sequence')
df1.to_csv("/Users/zoey/Desktop/PS/phospho/file1.txt", header=True, sep='\t',index=False)
df2.to_csv("/Users/zoey/Desktop/PS/phospho/file2.txt", header=True, sep='\t',index=False)
df3.to_csv("/Users/zoey/Desktop/PS/phospho/file3.txt", header=True, sep='\t',index=False)
df4.to_csv("/Users/zoey/Desktop/PS/phospho/file4.txt", header=True, sep='\t',index=False)
df5.to_csv("/Users/zoey/Desktop/PS/phospho/file5.txt", header=True, sep='\t',index=False)
df6.to_csv("/Users/zoey/Desktop/PS/phospho/file6.txt", header=True, sep='\t',index=False)
df7.to_csv("/Users/zoey/Desktop/PS/phospho/file7.txt", header=True, sep='\t',index=False)
df8.to_csv("/Users/zoey/Desktop/PS/phospho/file8.txt", header=True, sep='\t',index=False)
df9.to_csv("/Users/zoey/Desktop/PS/phospho/file9.txt", header=True, sep='\t',index=False)
df10.to_csv("/Users/zoey/Desktop/PS/phospho/file10.txt", header=True, sep='\t',index=False)

# data interpret
import pandas as pd 
import operator as op
import numpy as np 
# data = pd.read_csv("/Users/zoey/Desktop/PS/phospho/PSP-BulkSequenceSearch.txt", sep = '\t', index_col=False, on_bad_lines='warn')
data1 = pd.read_csv("/Users/zoey/Desktop/PS/phospho/nonPSP-BulkSequenceSearch1.txt", sep = '\t', index_col=False, on_bad_lines='warn')
data2 = pd.read_csv("/Users/zoey/Desktop/PS/phospho/nonPSP-BulkSequenceSearch2.txt", sep = '\t', index_col=False, on_bad_lines='warn')
data3 = pd.read_csv("/Users/zoey/Desktop/PS/phospho/nonPSP-BulkSequenceSearch3.txt", sep = '\t', index_col=False, on_bad_lines='warn')
data4 = pd.read_csv("/Users/zoey/Desktop/PS/phospho/nonPSP-BulkSequenceSearch4.txt", sep = '\t', index_col=False, on_bad_lines='warn')
data5 = pd.read_csv("/Users/zoey/Desktop/PS/phospho/nonPSP-BulkSequenceSearch5.txt", sep = '\t', index_col=False, on_bad_lines='warn')
data6 = pd.read_csv("/Users/zoey/Desktop/PS/phospho/nonPSP-BulkSequenceSearch6.txt", sep = '\t', index_col=False, on_bad_lines='warn')
data7 = pd.read_csv("/Users/zoey/Desktop/PS/phospho/nonPSP-BulkSequenceSearch7.txt", sep = '\t', index_col=False, on_bad_lines='warn')
data8 = pd.read_csv("/Users/zoey/Desktop/PS/phospho/nonPSP-BulkSequenceSearch8.txt", sep = '\t', index_col=False, on_bad_lines='warn')
data9 = pd.read_csv("/Users/zoey/Desktop/PS/phospho/nonPSP-BulkSequenceSearch9.txt", sep = '\t', index_col=False, on_bad_lines='warn')
data10 = pd.read_csv("/Users/zoey/Desktop/PS/phospho/nonPSP-BulkSequenceSearch10.txt", sep = '\t', index_col=False, on_bad_lines='warn')
data = data1.append(data2,ignore_index=True).append(data3,ignore_index=True).append(data4,ignore_index=True).append(data5,ignore_index=True).append(data6,ignore_index=True).append(data7,ignore_index=True).append(data8,ignore_index=True).append(data9,ignore_index=True).append(data10,ignore_index=True)
data['score'] = np.nan
# data.to_csv("/Users/zoey/Desktop/PS/phospho/nonPSP-BulkSequenceSearch.txt")
for i in range(len(data['organism'])):
    if type(data['Site'][i]) is float:
        data['score'][i] = 0
    else:
        site = data['Site'][i].split(',')
        pho = len([x for x in site if op.contains(x,'-ac')])
        score = pho/(data['End position'][i] - data['Start position'][i] + 1)
        data['score'][i] = score
data['Accession'] = data['Accession'].str.replace('UP:','')
data0 = pd.read_csv("/Users/zoey/Desktop/PS/output/293tAll_2.csv")
# data0['Acetylation'] = np.nan
for i in range(0,len(data0['Names'])):
    name = data0['Names'][i]
    for j in range(0,len(data['Accession'])):
        if (op.contains(name, data['Accession'][j])):
            data0['Acetylation'][i] = data['score'][j]
data0.to_csv("/Users/zoey/Desktop/PS/output/293tAll_2.csv")

# deepcoil
from deepcoil import DeepCoil
from deepcoil.utils import plot_preds
from Bio import SeqIO
dc = DeepCoil(use_gpu=True)
inp = {str(entry.id): str(entry.seq) for entry in SeqIO.parse('example/example.fas', 'fasta')}
results = dc.predict(inp)
# plot_preds(results['3WPA_1'], out_file='example/example.png')
# result analyze
import pandas as pd 
import numpy as np 
import operator as op
base_dir = "/Users/zoey/Desktop/PS/DeepCoil/output/"
data = pd.DataFrame()
for k in range(10):
    dir = "nonPS_dc" + str(k) + ".csv"
    data0 = pd.read_csv(base_dir + dir, header=None, index_col=None)
    # data = data0[:3]
    df = data0.T
    df.iloc[0,0] = 'name'
    df.iloc[0,1] = 'cc'
    df.iloc[0,2] = 'hept'
    new_header = df.iloc[0]
    df = df[1:]
    df.columns = new_header
    df['#AA'] = np.nan
    df['score'] = np.nan
    for i in range(len(df['name'])):
        string = df.iloc[i,1]
        cc = list(string.split(" "))
        cc1 = [x.replace("[","").replace("]","").replace("\n","") for x in cc]
        cc2 = [x for x in cc1 if x != "..." and x!= '']
        df.iloc[i,3] = len(cc2)
        thresh = [x for x in cc2 if float(x) > 0.82]
        if thresh:
            df.iloc[i,4] = 1
        else: 
            df.iloc[i,4] = 0
    data = data.append(df)
da = pd.read_csv("/Users/zoey/Desktop/PS/output/293tnonPS8.csv")
da['DeepCoil'] = np.nan
for i in range(len(da['Names'])):
    name = da['Names'][i]
    for j in range(len(data['name'])):
        if (op.contains(name, data.iloc[j,0])):
            da['DeepCoil'][i] = data.iloc[j,4]
da.to_csv("/Users/zoey/Desktop/PS/output/293tnonPS9.csv")
    

# deep phase 
# use pretrained model but need to find images for each protein
import pandas as pd 
data = pd.read_csv("/Users/zoey/Desktop/PS/output/293tPS10.csv")
names = data['Names']
pid = []
for i in range(len(names)):
    pid.append(names[i][3:9])
with open("/Users/zoey/Desktop/PS/output/protID.txt", "w") as output:
    output.write(str(pid))

# calculate score for deep phase
import os 
import pandas as pd
import numpy as np 
dir = "/Volumes/yzy/PSdata/PS/"
df = []
for filename in os.listdir(dir):
    data = pd.read_csv(dir + filename)
    if len(data): 
        score = data.iloc[:,2].max()
        df.append([filename,score])
df = pd.DataFrame(data=df, columns=['name','score'])
df['name'] = df['name'].str[:15]
# df.to_csv("/Users/zoey/Desktop/PS/output/EN-DP.csv")
id = pd.read_csv("/Users/zoey/Desktop/PS/output/ensid.csv")
id = id.rename(columns={'Unnamed: 0':'names'})
df['id'] = np.nan
for ii in range(len(df)):
    for jj in range(len(id)):
        if df['name'][ii] == id['ensid'][jj]:
            df['id'][ii] = id['names'][jj]
# df.to_csv("/Users/zoey/Desktop/PS/output/EN-DP.csv")
tab = pd.read_csv("/Users/zoey/Desktop/PS/output/idmap.csv")
df['names'] = np.nan
for ii in range(len(df)):
    name = df['id'][ii]
    for jj in range(len(tab)):
        if tab['symbol'][jj] == name:
            df['names'][ii] = tab['query'][jj]
for ii in range(len(df)):
    if type(df['names'][ii]) == float:
        df = df.drop([ii])
df.to_csv("/Users/zoey/Desktop/PS/output/EN-DP.csv")

df = pd.read_csv("/Users/zoey/Desktop/PS/output/EN-DP.csv")
df = df.drop(columns=['Unnamed: 0'])
data = pd.read_csv("/Users/zoey/Desktop/PS/output/293tPS9.csv")
data['DeepPhase'] = np.nan
for ii in range(len(data)):
    name = data['Names'][ii]
    for jj in range(len(df)):
        if op.contains(name, df['names'][jj]):
            data['DeepPhase'][ii] = df['score'][jj]
data.to_csv("/Users/zoey/Desktop/PS/output/293tPS10.csv")

## embedding modify
import pandas as pd
data = pd.read_csv("/Users/zoey/Desktop/PS/output/293tAll.csv")
seq = []
names = []
for i in range(len(data)):
    if len(data['Sequences'][i]) <= 1000:
        seq.append(data['Sequences'][i])
        names.append(data['Names'][i])
df = {'id':names, 'sequence':seq}
dataf = pd.DataFrame(df)
dataf.to_csv("/Users/zoey/Desktop/PS/output/293tAll_Embed.csv")
import pandas as pd
import math
import os

# Load the CSV file
csv_file = '/Users/zoey/Desktop/PS/output/293tAll_Embed.csv'  # replace with your CSV file path
df = pd.read_csv(csv_file)

# Ensure the CSV has at least 'id' and 'sequence' columns
assert 'id' in df.columns and 'sequence' in df.columns, "CSV must contain 'id' and 'sequence' columns"

# Convert to FASTA format
fasta_content = ""
for index, row in df.iterrows():
    fasta_content += f">{row['id']}\n{row['sequence']}\n"

# Save the FASTA content to a file
fasta_file = '/Users/zoey/Desktop/PS/output/293tAll_Embed.fasta'
with open(fasta_file, 'w') as f:
    f.write(fasta_content)

# Split the FASTA file into 10 subfiles
with open(fasta_file, 'r') as f:
    lines = f.readlines()

# Count the number of sequences
num_sequences = len([line for line in lines if line.startswith(">")])
sequences_per_file = math.ceil(num_sequences / 10)

# Create directory for subfiles if it doesn't exist
subfiles_dir = '/Users/zoey/Desktop/PS/output/subfiles'
os.makedirs(subfiles_dir, exist_ok=True)

subfile_index = 0
sequence_count = 0
subfile_content = ""

for line in lines:
    if line.startswith(">"):
        if sequence_count == sequences_per_file:
            # Save the current subfile
            subfile_path = os.path.join(subfiles_dir, f'subfile_{subfile_index + 1}.fasta')
            with open(subfile_path, 'w') as subfile:
                subfile.write(subfile_content)
            subfile_index += 1
            sequence_count = 0
            subfile_content = ""

        sequence_count += 1
    subfile_content += line

# Save the last subfile
if subfile_content:
    subfile_path = os.path.join(subfiles_dir, f'subfile_{subfile_index + 1}.fasta')
    with open(subfile_path, 'w') as subfile:
        subfile.write(subfile_content)

print("FASTA file split into 10 subfiles successfully.")
