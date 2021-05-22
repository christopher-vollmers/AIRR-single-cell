import os

#os.system('python3 run_igblast_IGH.py /home/vollmers/data/10X_PBMCs/AIRR/Rep1/Rep1_IGH.fasta /home/vollmers/data/10X_PBMCs/AIRR/Rep1/IGH')
#os.system('python3 run_igblast_IGL.py /home/vollmers/data/10X_PBMCs/AIRR/Rep1/Rep1_IGLL.fasta /home/vollmers/data/10X_PBMCs/AIRR/Rep1/IGL')
#os.system('python3 run_igblast_IGK.py /home/vollmers/data/10X_PBMCs/AIRR/Rep1/Rep1_IGLK.fasta /home/vollmers/data/10X_PBMCs/AIRR/Rep1/IGK')
#os.system('python3 run_igblast_TRA.py /home/vollmers/data/10X_PBMCs/AIRR/Rep1/Rep1_TCRA.fasta /home/vollmers/data/10X_PBMCs/AIRR/Rep1/TRA')
#os.system('python3 run_igblast_TRB.py /home/vollmers/data/10X_PBMCs/AIRR/Rep1/Rep1_TCRB.fasta /home/vollmers/data/10X_PBMCs/AIRR/Rep1/TRB')

IGH_file='/mnt/holocron4/10x_PBMCs/4rep2/analysis_04072021/AIRR/IGH.fasta.table.parsed'
IGL_file='/mnt/holocron4/10x_PBMCs/4rep2/analysis_04072021/AIRR/IGLL.fasta.table.parsed'
IGK_file='/mnt/holocron4/10x_PBMCs/4rep2/analysis_04072021/AIRR/IGLK.fasta.table.parsed'
TRA_file='/mnt/holocron4/10x_PBMCs/4rep2/analysis_04072021/AIRR/TCRA.fasta.table.parsed'
TRB_file='/mnt/holocron4/10x_PBMCs/4rep2/analysis_04072021/AIRR/TCRB.fasta.table.parsed'
out=open('/mnt/holocron4/10x_PBMCs/4rep2/analysis_04072021/AIRR/combined.txt','w')

def add_cell(cell_dict):
    cell_dict[cell]={}
    cell_dict[cell]['IGHV']=[]
    cell_dict[cell]['IGHD']=[]
    cell_dict[cell]['IGHJ']=[]
    cell_dict[cell]['Isotype']=[]
    cell_dict[cell]['IGLV']=[]
    cell_dict[cell]['IGLJ']=[]
    cell_dict[cell]['IGKV']=[]
    cell_dict[cell]['IGKJ']=[]
    cell_dict[cell]['IGH_Mutations']=[]
    cell_dict[cell]['IGH_Abundance']=[]
    cell_dict[cell]['IGL_Mutations']=[]
    cell_dict[cell]['IGL_Abundance']=[]
    cell_dict[cell]['IGK_Mutations']=[]
    cell_dict[cell]['IGK_Abundance']=[]
    cell_dict[cell]['TRA_Mutations']=[]
    cell_dict[cell]['TRA_Abundance']=[]
    cell_dict[cell]['TRB_Mutations']=[]
    cell_dict[cell]['TRB_Abundance']=[]


    cell_dict[cell]['TRAV']=[]
    cell_dict[cell]['TRAJ']=[]
    cell_dict[cell]['TRBV']=[]
    cell_dict[cell]['TRBD']=[]
    cell_dict[cell]['TRBJ']=[]
    return cell_dict


cell_dict={}
for line in open(IGH_file):
    a=line.strip().split('\t')
    cell=a[0].split('_')[2].split('-')[0]
    abundance=a[0].split('_')[-1]
    if cell not in cell_dict:
        cell_dict=add_cell(cell_dict)
    if a[8]!='_':
        cell_dict[cell]['IGHV'].append(a[2])
        cell_dict[cell]['IGHD'].append(a[3])
        cell_dict[cell]['IGHJ'].append(a[4])
        cell_dict[cell]['Isotype'].append(a[8]+'_'+a[7])
        if len(a)>10:
            Mutations=len(a[10].split(','))
        else:
            Mutations=0
        cell_dict[cell]['IGH_Mutations'].append(str(Mutations))
        cell_dict[cell]['IGH_Abundance'].append(abundance)

for line in open(IGL_file):
    a=line.strip().split('\t')
    cell=a[0].split('_')[2].split('-')[0]
    abundance=a[0].split('_')[-1]
    if cell not in cell_dict:
        cell_dict=add_cell(cell_dict)
    cell_dict[cell]['IGLV'].append(a[2])
    cell_dict[cell]['IGLJ'].append(a[4])
    if len(a)>10:
        Mutations=len(a[10].split(','))
    else:
        Mutations=0
    cell_dict[cell]['IGL_Mutations'].append(str(Mutations))
    cell_dict[cell]['IGL_Abundance'].append(abundance)

for line in open(IGK_file):
    a=line.strip().split('\t')
    cell=a[0].split('_')[2].split('-')[0]
    abundance=a[0].split('_')[-1]
    if cell not in cell_dict:
        cell_dict=add_cell(cell_dict)
    cell_dict[cell]['IGKV'].append(a[2])
    cell_dict[cell]['IGKJ'].append(a[4])
    if len(a)>10:
        Mutations=len(a[10].split(','))
    else:
        Mutations=0
    cell_dict[cell]['IGK_Mutations'].append(str(Mutations))
    cell_dict[cell]['IGK_Abundance'].append(abundance)

for line in open(TRA_file):
    a=line.strip().split('\t')
    cell=a[0].split('_')[2].split('-')[0]
    abundance=a[0].split('_')[-1]
    if cell not in cell_dict:
        cell_dict=add_cell(cell_dict)
    cell_dict[cell]['TRAV'].append(a[2])
    cell_dict[cell]['TRAJ'].append(a[4])
    if len(a)>10:
        Mutations=len(a[10].split(','))
    else:
        Mutations=0
    cell_dict[cell]['TRA_Mutations'].append(str(Mutations))
    cell_dict[cell]['TRA_Abundance'].append(abundance)

for line in open(TRB_file):
    a=line.strip().split('\t')
    cell=a[0].split('_')[2].split('-')[0]
    abundance=a[0].split('_')[-1]
    if cell not in cell_dict:
        cell_dict=add_cell(cell_dict)
    cell_dict[cell]['TRBV'].append(a[2])
    cell_dict[cell]['TRBD'].append(a[3])
    cell_dict[cell]['TRBJ'].append(a[4])
    if len(a)>10:
        Mutations=len(a[10].split(','))
    else:
        Mutations=0
    cell_dict[cell]['TRB_Mutations'].append(str(Mutations))
    cell_dict[cell]['TRB_Abundance'].append(abundance)

count={}
count['IGH']=0
count['IGL']=0
count['IGK']=0
count['IG_paired']=0
count['TRA']=0
count['TRB']=0
count['TR_paired']=0


for cell in cell_dict:
    IGH=False
    IGL=False
    IGK=False
    TRA=False
    TRB=False
    target=cell_dict[cell]
    if target['IGHV']:
        count['IGH']+=1
        IGH=True
    if target['IGLV']:
        count['IGL']+=1
        IGL=True
    if target['IGKV']:
        count['IGK']+=1
        IGK=True
    if target['TRAV']:
        count['TRA']+=1
        TRA=True
    if target['TRBV']:
        count['TRB']+=1
        TRB=True
    if IGH:
        if IGL or IGK:
            count['IG_paired']+=1
    if TRA and TRB:
        count['TR_paired']+=1


    out.write(cell+'\t')
    out.write((',').join(target['IGHV'])+'\t')
    out.write((',').join(target['IGHD'])+'\t')
    out.write((',').join(target['IGHJ'])+'\t')
    out.write((',').join(target['Isotype'])+'\t')
    out.write((',').join(target['IGH_Mutations'])+'\t')
    out.write((',').join(target['IGH_Abundance'])+'\t')

    out.write((',').join(target['IGKV'])+'\t')
    out.write((',').join(target['IGKJ'])+'\t')
    out.write((',').join(target['IGK_Mutations'])+'\t')
    out.write((',').join(target['IGK_Abundance'])+'\t')

    out.write((',').join(target['IGLV'])+'\t')
    out.write((',').join(target['IGLJ'])+'\t')
    out.write((',').join(target['IGL_Mutations'])+'\t')
    out.write((',').join(target['IGL_Abundance'])+'\t')

    out.write((',').join(target['TRAV'])+'\t')
    out.write((',').join(target['TRAJ'])+'\t')
    out.write((',').join(target['TRA_Mutations'])+'\t')
    out.write((',').join(target['TRA_Abundance'])+'\t')

    out.write((',').join(target['TRBV'])+'\t')
    out.write((',').join(target['TRBD'])+'\t')
    out.write((',').join(target['TRBJ'])+'\t')
    out.write((',').join(target['TRB_Mutations'])+'\t')
    out.write((',').join(target['TRB_Abundance'])+'\t')



    out.write('\n')

print(count)
