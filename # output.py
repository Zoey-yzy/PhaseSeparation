# output
"""output file from result data to .csv"""
# PLAAC 
import pandas as pd 
data = pd.read_csv("/Users/zoey/Desktop/PS/PLAAC/output.txt", sep = "\t")
score = data[["SEQid", "NLLR"]]
score.to_csv('/Users/zoey/Desktop/PS/output/PLAAC.csv')

# PScore
data = []
with open("/Users/zoey/Desktop/PS/PScore/293tPS.out", 'r') as data_file:
        for line in data_file:
                data.append(line.split())
table = [[x[2], x[1]] for x in data]
df = pd.DataFrame(table)
df.rename(columns = {'0':'Name', '1':'PScore'}, inplace=True)

# IDR
all = pd.read_csv("") # read all protein in the .csv
ls_v, ls_n = [], []
for i in range(0,len(all)):
    #if file_name contain 'fasta.stats':
        data = pd.read_csv("/Users/zoey/Desktop/PS/IDR-Espritz/batch/spO00154.fasta.stats", sep=':', header=None)
        IDR_score = data.loc[data[0]=="Total % disorder",1].values
        ls_v.insert(i,float(IDR_score.item()))
        ls_n.insert(i,file_name.proid)
ls = [ls_n, ls_v]

# input
import pandas as pd 
data1 = pd.read_csv("/Users/zoey/Desktop/PS/PSdata/293tall1-2500.csv")
data2 = pd.read_csv("/Users/zoey/Desktop/PS/PSdata/293tall2501-5000.csv")
data3 = pd.read_csv("/Users/zoey/Desktop/PS/PSdata/293tall5001-6302.csv")
data0 = pd.read_csv("/Users/zoey/Desktop/PS/PSdata/293tpsp.csv")
data_all = data1.append(data2,ignore_index=True).append(data3,ignore_index=True)
data_all.rename(columns={'Unnamed: 0':'names','Amino_Acids_Number':'AA_number'}, inplace=True)
data0.rename(columns={'Unnamed: 0':'names','Amino_Acids_Number':'AA_number'}, inplace=True)
ls_ps, ls_all = data0.values.tolist(), data_all.values.tolist()
def find_dup(ls, x):
        count = 0
        for x in ls:
                if x in ls: 
                        count += 1 
        return count 

def remove_same(ls0):
        ls1 = [x for x in ls0 if find_dup(ls,x) == 1]
        ls2 = [x for x in ls0 if find_dup(ls,x) > 1]
        ls3 = []
        for x in ls2:
                if x not in ls3:
                        ls3.append(x)
        return ls1 + ls3

ls_psm, ls_allm = remove_same(ls_ps), remove_same(ls_all)

def complement(ls1,ls2):
        ls = []
        lso = []
        for ele in ls1:
                if ele not in ls2:
                        ls.append(ele)
                else: 
                        lso.append(ele)          
        return ls, lso
def intersection(ls1, ls2):
        ls = []
        for ele in ls1:
                if ele in ls2:
                        ls.append(ele)
        return ls
ls_nps, ls_npsc = complement(ls_allm, ls_psm)
def output(ls):                        
        dic = {'Names':[x[0] for x in ls], 'AA_Number':[x[1] for x in ls], 'Sequences':[x[2] for x in ls]}
        out = pd.DataFrame(data=dic)
        return out
psp = output(ls_npsc)
npsp = output(ls_nps)
all = output(ls_allm)        
psp.to_csv("/Users/zoey/Desktop/PS/PSdata/293tpsp.csv")
npsp.to_csv("/Users/zoey/Desktop/PS/PSdata/293tnpsp.csv")
all.to_csv("/Users/zoey/Desktop/PS/PSdata/293tall.csv")
# check for protein
ls_inter = intersection(ls_npsc, ls_ps)
# >>> return to a list len = 169, which is equal to the length of ls_npsc