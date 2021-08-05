# importation
from Bio import AlignIO
from Bio.Align import AlignInfo
import csv
import itertools
import pandas as pd
# create dictionary
gene={}
f=open("/Users/bee4/Documents/Summer 2021/Sh_topo4A_align.fasta","r")
alignment = AlignIO.read(f,"fasta")
f.close()
# create gap consensus and append mutations to list
summary=AlignInfo.SummaryInfo(alignment)
consensus= summary.gap_consensus(threshold=1)
dictionary={}
for al in alignment:
    place=0
    mutations=[]
    for x in al.seq._data:
        if(x!=consensus._data[place]):
            mutations.append(x)
        if(x==consensus._data[place] and x=="X"):
            mutations.append(x)
        place+=1
    # add entry to dictionary for each sequence
    dictionary[al.id]=mutations
# save amount of mutations
gene["topo4A"]=len(mutations)
f=open("/Users/bee4/Documents/Summer 2021/Sh_topo4B_align.fasta","r")
alignment = AlignIO.read(f,"fasta")
f.close()
# create gap consensus and append mutations to list

summary=AlignInfo.SummaryInfo(alignment)
consensus= summary.gap_consensus(threshold=1)

for al in alignment:
    place=0
    mutations=[]
    for x in al.seq._data:
        if(x!=consensus._data[place]):
            mutations.append(x)
        if(x==consensus._data[place] and x=="X"):
            mutations.append(x)
        place+=1
    # extend entry in dictionary for each sequence
        
    dictionary[al.id].extend(mutations)
# save amount of mutations
gene["topo4B"]=gene["topo4A"]+len(mutations)

f=open("/Users/bee4/Documents/Summer 2021/Sh_gyrA_align.fasta","r")
alignment = AlignIO.read(f,"fasta")
f.close()
# create gap consensus and append mutations to list

summary=AlignInfo.SummaryInfo(alignment)
consensus= summary.gap_consensus(threshold=1)

for al in alignment:
    place=0
    mutations=[]
    for x in al.seq._data:
        if(x!=consensus._data[place]):
            mutations.append(x)
        if(x==consensus._data[place] and x=="X"):
            mutations.append(x)
        place+=1
    # extend entry in dictionary for each sequence
    dictionary[al.id].extend(mutations)
# save amount of mutations
gene["gyrA"]=gene["topo4B"]+len(mutations)

f=open("/Users/bee4/Documents/Summer 2021/Sh_gyrB_align.fasta","r")
alignment = AlignIO.read(f,"fasta")
f.close()
# create gap consensus and append mutations to list
summary=AlignInfo.SummaryInfo(alignment)
consensus= summary.gap_consensus(threshold=1)

for al in alignment:
    place=0
    mutations=[]
    for x in al.seq._data:
        if(x!=consensus._data[place]):
            mutations.append(x)
        if(x==consensus._data[place] and x=="X"):
            mutations.append(x)
        place+=1
    # extend entry in dictionary for each sequence
    dictionary[al.id].extend(mutations)
# save amount of mutations
gene["gyrB"]=gene["gyrA"]+len(mutations)

# find unique alleles
alleles=[]
for entry in dictionary:
    if dictionary[entry] not in alleles:
        alleles.append(dictionary[entry])
outfile=open("/Users/bee4/Documents/Summer 2021/Sh_alleles.txt","w")
# write unique alleles
outfile.write("Number of unique alleles: "+str(len(alleles))+"\n")
for a in alleles:
    outfile.write(str(a)+"\n")
outfile.close()
# create MIC value dictionary
f=open("/Users/bee4/Documents/Summer 2021/Resistance_MIC.csv","r")
mic_values=csv.reader(f)
mut_amount=len(a)
mic_dict={}
#for i in range(0,mut_amount-1):
for mi in mic_values:
    if mi[0][0:2]=="AA" and mi[1].isdigit():
        mic_dict[mi[0]]=mi[1]
# implementing code for correlation
correlation=[]
correlation_values=[]
# take the size of a list within the dictionary as the range
for i in range(1,mut_amount+1):
    #create dictionary
    my_combos={}
    # save all possible combinations of mutations at certain locations
    for a in dictionary:
        position=range(len(dictionary[a]))
        combo_object=itertools.combinations(position,i)
        #my_combos[a]=(combo_object)
        my_combos[a]=list(combo_object)
        #length = sum(1 for ignore in my_combos[a])
    #for c in range(length):
    for c in range(len(my_combos[a])):
        #create empty lists
        l=[]
        m=[]
        # find all mutation letters in a particular location
        for acids in my_combos:
            ord_value=""
            for b in my_combos[acids][c]:
                ord_value+=dictionary[acids][b]
            # append mutations and mic values to lists
            l.append(ord_value)
            m.append(mic_dict[acids])
        # create empty list
        DataValues=[]
        # create a list of lists
        for z in range(len(l)):
            DataValues.append([l[z],m[z]])
        
        #Create the Data Frame
        ColumnNames=['MUTATIONS','MIC']
        LoanData=pd.DataFrame(data=DataValues,columns=ColumnNames)

        CrosstabResult=pd.crosstab(index=LoanData['MUTATIONS'],columns=LoanData['MIC'])         
        # importing the required function
        from scipy.stats import chi2_contingency
         
        # Performing Chi-sq test
        ChiSqResult = chi2_contingency(CrosstabResult)
        correlation.append(ChiSqResult[1])
        correlation_values.append(my_combos[a][c])
    # find best index and print correlation
    # print gene and MIC values
best_index=correlation.index(min(correlation))
outfile=open("/Users/bee4/Documents/Summer 2021/Sh_results.txt","w")
outfile.write(str(correlation[best_index])+" "+str(correlation_values[best_index])+"\n")
for i in correlation_values[best_index]:
    if i<gene["topo4A\n"]:
        outfile.write("topo4A\n")
    elif i<gene["topo4B"]:
        outfile.write("topo4B\n")
    elif i<gene["gyrA"]:
        outfile.write("gyrA\n")
    elif i<gene["gyrB"]:
        outfile.write("gyrB\n")
for am in dictionary:
    outfile.write(str(am)+",")
    for c in correlation_values[best_index]:
        outfile.write(str(dictionary[am][c])+",")
    outfile.write(str(mic_dict[am])+"\n")
outfile.close()

# for i in range(mut_amount):
#     l=[]
#     m=[]
#     for name in dictionary:
#         combinations=[]
#         for r in range(1,len(dictionary[name])+1):
#             combo_object=itertools.combinations(dictionary[name],r)
#             combo_list=list(combo_object)
#             combinations+=combo_list
#         for c in combinations:
#             ord_value=0
#             for acids in c:
#                 ord_value+=ord(acids)
#             l.append(ord_value)
#             m.append(int(mic_dict[name]))
#     z=numpy.corrcoef(numpy.array(l),numpy.array(m))
#     correlation.append(abs(z[0,1]))
# print(correlation.index(max(correlation)))
