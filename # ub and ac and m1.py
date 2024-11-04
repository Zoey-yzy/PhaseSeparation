# ub and ac and m1
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
data0 = pd.read_csv("/Users/zoey/Desktop/PS/output/293tAll.csv")
data0['Acetylation'] = np.nan
for i in range(0,len(data0['Names'])):
    name = data0['Names'][i]
    for j in range(0,len(data['Accession'])):
        if (op.contains(name, data['Accession'][j])):
            data0['Acetylation'][i] = data['score'][j]
data0.to_csv("/Users/zoey/Desktop/PS/output/293tAll_2.csv")