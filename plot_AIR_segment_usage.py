import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.stats as stats
import matplotlib.image as mpimg
import matplotlib.patches as mplpatches
import sys
sys.path.append('~/Downloads')
#import BME163_Plot_functions_Vollmers_Christopher as Custom_plots

colors={}
colors['IGHM']='black'
colors['IGHD']=(1/255,188/255,98/255)
colors['IGHA1']=(200/255,0/255,80/255)
colors['IGHA2']=(200/255,0/255,80/255)
colors['IGHG1']=(56/255,66/255,156/255)
colors['IGHG2']=(56/255,66/255,156/255)
colors['IGHG3']=(56/255,66/255,156/255)
colors['IGHG4']=(56/255,66/255,156/255)
colors['IGHE']='yellow'

def segment_order(infile,position):
    segments={}
    start=position
    for line in open(infile):
        a=line.strip()
        segments[a]=position
        position+=1
    end=position
    rectangle1=mplpatches.Rectangle((start,400),end-start,1,facecolor='black', edgecolor='black',linewidth=0)
    panel2.add_patch(rectangle1)
    return segments,position

def plot_composition(composition,y_pos,panel):

    print(composition)
    if composition[0].split('*')[0] in IGHV:
        x_pos=IGHV[composition[0].split('*')[0]]
        print(y_pos,'V',x_pos)
        rectangle1=mplpatches.Rectangle((x_pos,y_pos-0.25),1,1.5,facecolor='black', edgecolor='black',linewidth=0)
        panel.add_patch(rectangle1)

    if composition[1].split('*')[0] in IGHD:
        x_pos=IGHD[composition[1].split('*')[0]]
        print(y_pos,'D',x_pos)
        rectangle1=mplpatches.Rectangle((x_pos,y_pos-0.25),1,1.5,facecolor='black', edgecolor='black',linewidth=0)
        panel.add_patch(rectangle1)

    if composition[2].split('*')[0] in IGHJ:
        x_pos=IGHJ[composition[2].split('*')[0]]
        print(y_pos,'J',x_pos)

        rectangle1=mplpatches.Rectangle((x_pos,y_pos-0.25),1,1.5,facecolor='black', edgecolor='black',linewidth=0)
        panel.add_patch(rectangle1)
    for iso in composition[3].split(','):
        iso=iso.split('_')[0]
        if iso in Isotype_dict:
            x_pos=Isotype_dict[iso]
            print(y_pos,'Isotype',x_pos)
            rectangle1=mplpatches.Rectangle((x_pos,y_pos-0.25),1,1.5,facecolor=colors[iso], edgecolor=colors[iso],linewidth=0)
            panel.add_patch(rectangle1)



    if composition[4].split('*')[0] in IGKV:
        x_pos=IGKV[composition[4].split('*')[0]]
        print(y_pos,'KV',x_pos)
        rectangle1=mplpatches.Rectangle((x_pos,y_pos-0.25),1,1.5,facecolor='black', edgecolor='black',linewidth=0)
        panel.add_patch(rectangle1)

    if composition[5].split('*')[0] in IGKJ:
        x_pos=IGKJ[composition[5].split('*')[0]]
        print(y_pos,'KJ',x_pos)
        rectangle1=mplpatches.Rectangle((x_pos,y_pos-0.25),1,1.5,facecolor='black', edgecolor='black',linewidth=0)
        panel.add_patch(rectangle1)

    if composition[6].split('*')[0] in IGLV:
        x_pos=IGLV[composition[6].split('*')[0]]
        print(y_pos,'LV',x_pos)
        rectangle1=mplpatches.Rectangle((x_pos,y_pos-0.25),1,1.5,facecolor='black', edgecolor='black',linewidth=0)
        panel.add_patch(rectangle1)

    if composition[7].split('*')[0] in IGLJ:
        x_pos=IGLJ[composition[7].split('*')[0]]
        print(y_pos,'LJ',x_pos)
        rectangle1=mplpatches.Rectangle((x_pos,y_pos-0.25),1,1.5,facecolor='black', edgecolor='black',linewidth=0)
        panel.add_patch(rectangle1)

    HeavyChain=False
    LightChain=False
    if composition[0]!='':
        rectangle1=mplpatches.Rectangle((-8,y_pos-0.25),6,1.5,facecolor='grey', edgecolor='grey',linewidth=0)
        panel.add_patch(rectangle1)
        HeavyChain=True
    if composition[4]!='':
        rectangle1=mplpatches.Rectangle((-1,y_pos-0.25),6,1.5,facecolor=(0/255,128/255,128/255), edgecolor=(0/255,128/255,128/255),linewidth=0)
        panel.add_patch(rectangle1)
        LightChain=True
    if composition[6]!='':
        rectangle1=mplpatches.Rectangle((6,y_pos-0.25),6,1.5,facecolor=(225/255,165/255,0/255), edgecolor=(225/255,165/255,0/255),linewidth=0)
        panel.add_patch(rectangle1)
        LightChain=True
    if HeavyChain and LightChain:
        print('Paired')
        rectangle1=mplpatches.Rectangle((-20,y_pos-0.25),11,1.5,facecolor='black', edgecolor='black',linewidth=0)
        panel.add_patch(rectangle1)

plt.style.use('~/Downloads/BME163.mplstyle')
fig_1 = plt.figure(figsize=(2,6))
panel2=plt.axes([0, 0 , 1 ,1],frameon=False)

position=15
IGHV,position=segment_order('/home/bd1/Downloads/IGH_Pipeline/IGH_locus_V_segments.txt',position)

position+=5
IGHD,position=segment_order('/home/bd1/Downloads/IGH_Pipeline/IGH_locus_D_segments.txt',position)
position+=5
IGHJ,position=segment_order('/home/bd1/Downloads/IGH_Pipeline/IGH_locus_J_segments.txt',position)
position+=5
Isotype_dict,position=segment_order('/home/bd1/Downloads/IGH_Pipeline/Isotypes.txt',position)

position+=10
IGKV,position=segment_order('/home/bd1/Downloads/IGH_Pipeline/IGK_locus_V_segments.txt',position)
position+=5
IGKJ,position=segment_order('/home/bd1/Downloads/IGH_Pipeline/IGK_locus_J_segments.txt',position)

position+=10
IGLV,position=segment_order('/home/bd1/Downloads/IGH_Pipeline/IGL_locus_V_segments.txt',position)
position+=5
IGLJ,position=segment_order('/home/bd1/Downloads/IGH_Pipeline/IGL_locus_J_segments.txt',position)


isotypes=[]
x_list=[]
y_list=[]
xx_list=[]
yy_list=[]
count=0
name_list=[]
bTotal=0
tTotal=0
Hcount=0
Kcount=0
Lcount=0
Acount=0
Bcount=0
bPaired=0
tPaired=0

y_pos=0
files=[['/mnt/holocron4/10x_PBMCs/4rep1/analysis_04072021/AIRR/combined.txt','/mnt/holocron4/10x_PBMCs/4rep1/analysis_04072021/fc/tsne/R2C2_41_cell_to_celltype'],['/mnt/holocron4/10x_PBMCs/4rep2/analysis_04072021/AIRR/combined.txt','/mnt/holocron4/10x_PBMCs/4rep2/analysis_04072021/fc/tsne/R2C2_42_cell_to_celltype']]
file_count=0
Segments={}
for file,assignment in files:
    file_count+=1

    for line in open(file):

        count+=1
        a=line.split('\t')
        IGH_mutation=''
        IGL_mutation=''
        IGK_mutation=''
        name=a[0]
        VH=''
        JH=''
        DH=''
        VL=''
        JL=''
        VK=''
        JK=''
        VA=''
        JA=''
        VB=''
        DB=''
        JB=''
        Isotypes=''
        AbundanceH=''
        AbundanceK=''
        AbundanceL=''
        if a[1]:
            IGH_mutation=a[5].strip().split(',')
#            print(count,'H',IGH_mutation)
            IGH_mutation=int(IGH_mutation[0])
            Isotype=a[4].split('_')[0]
            Isotypes=a[4]
            print(Isotypes)
            VH=a[1].split(',')[0]
            DH=a[2].split(',')[0]
            JH=a[3].split(',')[0]
            AbundanceH=sum(np.array(a[6].split(','),dtype=int))

        if a[7]:
            IGL_mutation=a[9].strip().split(',')
#            print(count,'L',IGL_mutation)
            IGL_mutation=int(IGL_mutation[0])
            VL=a[7].split(',')[0]
            JL=a[8].split(',')[0]
            AbundanceL=sum(np.array(a[10].split(','),dtype=int))


        if a[11]:
            IGK_mutation=a[13].strip().split(',')
#            print(count,'K',IGK_mutation)
            IGK_mutation=int(IGK_mutation[0])
            VK=a[11].split(',')[0]
            JK=a[12].split(',')[0]
            AbundanceK=sum(np.array(a[14].split(','),dtype=int))

        if a[15]:
            VA=a[15].split(',')[0]
            JA=a[16].split(',')[0]

        if a[19]:
            VB=a[19].split(',')[0]
            DB=a[20].split(',')[0]
            JB=a[21].split(',')[0]

        if a[19]:
            VB=a[19].split(',')[0]
            DB=a[20].split(',')[0]
            JB=a[21].split(',')[0]


        Segments['Rep'+str(file_count)+'_'+a[0]]=(VH,DH,JH,Isotypes,VL,JL,VK,JK,VA,JA,VB,DB,JB)
        if IGH_mutation!='':
            if IGK_mutation!='':
                print(count,'H','K',IGH_mutation,IGK_mutation)
                name_list.append(str(file_count)+'_'+name)
                x_list.append(IGH_mutation)
                y_list.append(IGK_mutation)
                xx_list.append(AbundanceH)
                yy_list.append(AbundanceK)
                isotypes.append(Isotype)
            elif IGL_mutation!='':
                name_list.append(str(file_count)+'_'+name)
                print(count,'H','L',IGH_mutation,IGL_mutation)
                x_list.append(IGH_mutation)
                y_list.append(IGL_mutation)
                xx_list.append(AbundanceH)
                yy_list.append(AbundanceL)

                isotypes.append(Isotype)

    line_set=[]
    shuffled_lines=[]


    for line in open(assignment):
        line_set.append(line)
    random_indexes=np.random.choice(range(0,len(line_set),1),size=len(line_set),replace=False)
    for index in random_indexes:
        shuffled_lines.append(line_set[index])
    draw=False
    for line in shuffled_lines:
        a=line.strip().split('\t')

        if a[1]=='bCell':
            bTotal+=1
            y_pos+=2
            if draw:
                print('DRAW',y_pos)
                facecolor=(240/255,240/255,240/255)
#                facecolor='red'
                rectangle1=mplpatches.Rectangle((-20,y_pos-0.5),300,2,facecolor=facecolor, edgecolor=(240/255,240/255,240/255),linewidth=0)
                panel2.add_patch(rectangle1)
                draw=False
            else:
                draw=True
            print('Rep'+str(file_count)+'_'+a[2])
            if 'Rep'+str(file_count)+'_'+a[2] in Segments:
                print('SEGMENTS',y_pos)
                Heavy=False
                Light=False
                Composition=Segments['Rep'+str(file_count)+'_'+a[2]]
                plot_composition(Composition,y_pos,panel2)
                if Composition[0]!='':
                    Hcount+=1
                    Heavy=True
                if Composition[4]!='':
                    Lcount+=1
                    Light=True
                if Composition[6]!='':
                    Kcount+=1
                    Light=True
                if Heavy and Light:
                    bPaired+=1
        if a[1]=='tCell':
            tTotal+=1
            if 'Rep'+str(file_count)+'_'+a[2] in Segments:
                Alpha=False
                Beta=False
                Composition=Segments['Rep'+str(file_count)+'_'+a[2]]
                if Composition[8]!='':
                    Acount+=1
                    Alpha=True
                if Composition[10]!='':
                    Bcount+=1
                    Beta=True
                if Alpha and Beta:
                    tPaired+=1

print(bTotal,\
      tTotal,\
      Hcount,\
      Kcount,\
      Lcount,\
      bPaired,\
      Acount,\
      Bcount,\
      tPaired)


panel2.tick_params(axis='both',which='both', bottom=False, right=False, left=False, top=False,labelbottom=False,labelleft=False,labeltop=False,labelright=False)
panel2.set_xlim(-20,260)
panel2.set_ylim(0,410)
#panel2.set_yticks(np.arange(0,51,10))
#panel2.set_xticks(np.arange(0,51,10))

plt.savefig('V_segment_usage.png',dpi=2400)
plt.close()





fig_1 = plt.figure(figsize=(2,2))
panel2=plt.axes([0.15, 0.15 , 0.8 ,0.8],frameon=True)


for index in range(0,len(x_list),1):
#    print(name_list[index],x_list[index],y_list[index],isotypes[index])
    panel2.plot(x_list[index],y_list[index],marker='o',markersize=2, mew=0,linewidth=0,mfc=colors[isotypes[index]])

print(stats.pearsonr(x_list,y_list))

panel2.tick_params(axis='both',which='both', bottom=True, right=False, left=True, top=False,labelbottom=True,labelleft=True,labeltop=False,labelright=False)
panel2.set_xlim(0,60)
panel2.set_ylim(0,60)
panel2.set_yticks(np.arange(0,61,10))
panel2.set_xticks(np.arange(0,61,10))

#panel3.set_yticklabels(sorted(V_segments,reverse=True),fontsize=4)

plt.savefig('Mutation_correlation.png',dpi=2400)
plt.close()

fig_1 = plt.figure(figsize=(2,2))
panel2=plt.axes([0.15, 0.15 , 0.8 ,0.8],frameon=True)

#xx_list=np.log2(np.array(xx_list))
#yy_list=np.log2(np.array(yy_list))
for index in range(0,len(xx_list),1):
#    print(name_list[index],x_list[index],y_list[index],isotypes[index])
    panel2.plot(xx_list[index],yy_list[index],marker='o',markersize=2, mew=0,linewidth=0,mfc=colors[isotypes[index]])

print(stats.pearsonr(xx_list,yy_list))

panel2.tick_params(axis='both',which='both', bottom=True, right=False, left=True, top=False,labelbottom=True,labelleft=True,labeltop=False,labelright=False)
#panel2.set_xlim(0,50)
#panel2.set_ylim(0,50)
#panel2.set_yticks(np.arange(0,51,10))
#panel2.set_xticks(np.arange(0,51,10))

#panel3.set_yticklabels(sorted(V_segments,reverse=True),fontsize=4)

plt.savefig('Abundance_correlation.png',dpi=2400)
plt.close()
