import sys
import os


path=sys.argv[1]

ighDict={}
for line in open(path+'/IGH.fasta.table.parsed'):
    a=line.strip().split('\t')
    cell_name=a[0].split('-')[0]
    if cell_name not in ighDict:
        ighDict[cell_name]=[]
    ighDict[cell_name].append(a[0].split('-')[1])


iglDict={}
for line in open(path+'/IGLL.fasta.table.parsed'):
    a=line.strip().split('\t')
    cell_name=a[0].split('-')[0]
    if cell_name not in iglDict:
        iglDict[cell_name]=[]
    iglDict[cell_name].append(a[0].split('-')[1])

igkDict={}
for line in open(path+'/IGLK.fasta.table.parsed'):
    a=line.strip().split('\t')
    cell_name=a[0].split('-')[0]
    if cell_name not in igkDict:
        igkDict[cell_name]=[]
    igkDict[cell_name].append(a[0].split('-')[1])


print('cells with IGH sequences', len(ighDict))
print('cells with IGL sequences', len(iglDict))
print('cells with IGK sequences', len(igkDict))


ighANDigl=0
ighANDigk=0
All=0
for cell in ighDict:
    if cell in iglDict:
        ighANDigl+=1
    if cell in igkDict:
        ighANDigk+=1
    if cell in iglDict and cell in igkDict:
        All+=1


print('cells with IGH and IGL sequences', ighANDigl)
print('cells with IGH and IGK sequences', ighANDigk)
print('cells with IGH and both IGL and IGK  sequences', All)


tcraDict={}
for line in open(path+'/TCRA.fasta.table.parsed'):
    a=line.strip().split('\t')
    cell_name=a[0].split('-')[0]
    if cell_name not in tcraDict:
       tcraDict[cell_name]=[]
    tcraDict[cell_name].append(a[0].split('-')[1])

tcrbDict={}
for line in open(path+'/TCRB.fasta.table.parsed'):
    a=line.strip().split('\t')
    cell_name=a[0].split('-')[0]
    if cell_name not in tcrbDict:
        tcrbDict[cell_name]=[]
    tcrbDict[cell_name].append(a[0].split('-')[1])


print('cells with TCRA sequences', len(tcraDict))
print('cells with TCRB sequences', len(tcrbDict))


tcraANDtcrb=0
for cell in tcraDict:
    if cell in tcrbDict:
        tcraANDtcrb+=1


print('cells with TCRA and TCRB sequences', tcraANDtcrb)
