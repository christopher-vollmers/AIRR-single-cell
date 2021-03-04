#!/usr/bin/env python3
#!/usr/bin/env python3
# Roger Volden

'''
This is going to be like defineAndQuantifyIsoformsWrapper.py
Except this is going to be more about controlling the exact
inputs going into each script because they had to be heavily
modified for taking each cell.

Usage
python3 antibodyIsotypeWrapper.py inDir
'''

import os
import sys
import numpy as np

def configReader(configIn):
    '''Parses the config file.'''
    progs = {}
    for line in open(configIn):
        if line.startswith('#') or not line.rstrip().split():
            continue
        line = line.rstrip().split('\t')
        progs[line[0]] = line[1]
    # should have minimap, poa, racon, water, consensus
    # check for extra programs that shouldn't be there
    possible = set(['poa', 'minimap2', 'gonk', 'consensus', 'racon', 'blat','emtrey','psl2pslx'])
    inConfig = set()
    for key in progs.keys():
        inConfig.add(key)
        if key not in possible:
            raise Exception('Check config file')
    # check for missing programs
    # if missing, default to path
    for missing in possible-inConfig:
        if missing == 'consensus':
            path = 'consensus.py'
        else:
            path = missing
        progs[missing] = path
        sys.stderr.write('Using ' + str(missing)
                         + ' from your path, not the config file.\n')
    return progs




def read_fastq_file(seq_file):
    '''
    Takes a FASTQ file and returns a list of tuples
    In each tuple:
        name : str, read ID
        seed : int, first occurrence of the splint
        seq : str, sequence
        qual : str, quality line
        average_quals : float, average quality of that line
        seq_length : int, length of the sequence
    '''

    lineNum=0
    lastPlus = False
    read_list=[]
    for line in open(seq_file):
        line = line.rstrip()
        if not line:
            continue
        # make an entry as a list and append the header to that list
        if lineNum % 4 == 0 and line[0] == '@':
            if lastPlus==True:
                read_list.append((name, seed, seq, qual, average_quals, seq_length))
            name_root =line[1:]
            name, seed = name_root, lineNum/4# int(name_root[1])
        if lineNum % 4 == 1:
            seq=line
            seq_length = len(seq)
        if lineNum % 4 == 2:
            lastPlus = True
        if lineNum % 4 == 3 and lastPlus:
            qual=line
            quals = []
            for character in qual:
                number = ord(character) - 33
                quals.append(number)
            average_quals = np.average(quals)
        lineNum += 1
    return read_list



def determine_consensus(name, fasta_reads,fastq_reads):
        '''Aligns and returns the consensus'''
#        print('determine')
        corrected_consensus = ''

#        fastq_reads=read_fastq_file(fastq)
        out_Fq=temp_folder + '/subsampled.fastq'
        out=open(out_Fq,'w')

        indexes=np.random.choice(np.arange(0,len(fastq_reads),1),min(len(fastq_reads),subsample),replace=False)
        subsample_fastq_reads=[]
        for index in indexes:
            subsample_fastq_reads.append(fastq_reads[index])
        for read in subsample_fastq_reads:
            out.write('@'+read[0]+'_'+str(read[1])+'\n'+read[2]+'\n+\n'+read[3]+'\n')
        print('subsampled',len(subsample_fastq_reads))



        out.close()


        poa_cons = temp_folder + '/consensus.fasta'
        final = temp_folder + '/corrected_consensus.fasta'
        overlap = temp_folder +'/overlaps.sam'
        pairwise = temp_folder + '/prelim_consensus.fasta'


        reads=fasta_reads
        repeats=0
        qual=[]
        raw=[]
        before=[]

        after=[]
        combined_name=name
        for read in reads:
            best=read


        out_cons_file = open(poa_cons, 'w')
        out_cons_file.write('>' + best + '\n' + reads[best].replace('-', '') + '\n')
        out_cons_file.close()



        final=poa_cons
        for i in np.arange(1,2,1):
            if i==1:
                input_cons=poa_cons
                output_cons=poa_cons.replace('.fasta','_'+str(i)+'.fasta')
            else:
                input_cons=poa_cons.replace('.fasta','_'+str(i-1)+'.fasta')
                output_cons=poa_cons.replace('.fasta','_'+str(i)+'.fasta')

#                print(i,minimap2,input_cons,out_Fq,overlap)
            minimap2_command='%s -t 1 --secondary=no -ax map-ont \
                              %s %s >%s 2>./minimap2_messages' \
                              %(minimap2, input_cons, out_Fq,overlap)
#            minimap2_process=subprocess.run(minimap2_command,shell=True)
            os.system(minimap2_command)
            racon_command='%s -q 5 -t 1 %s %s %s >%s 2>./racon_messages.txt' \
                      %(racon, out_Fq, overlap, input_cons, output_cons)
#            racon_process=subprocess.run(racon_command, shell=True)
            os.system(racon_command)
            final = output_cons




        reads = read_fasta(final)
        if len(reads)==0:
            reads = read_fasta(poa_cons)


        for read in reads:
            corrected_consensus = reads[read]



        return corrected_consensus


def make_consensus(name,fastq_reads,fasta_reads):
     if len(fasta_reads)>0:
        corrected_consensus=determine_consensus(name,fasta_reads,fastq_reads)
        return '>%s\n%s\n' %(name,corrected_consensus)
     else:
         return 'nope'



def read_fasta(infile):
    reads={}
    sequence=''

    for line in open(infile):
      if line:

        a=line.strip()
        if len(a)>0:
            if a[0]=='>':
                if sequence!='':
                    reads[name]=sequence
                name=a[1:]
                sequence=''
            else:
                sequence+=a
    if sequence!='':
        reads[name]=sequence
    return reads

def check_for_complete_air(fasta,fastq,V,D,J):
    current_dir=os.getcwd()
    fasta=fasta
    os.chdir('/home/ig88/Downloads/ncbi-igblast-1.14.0/')
    os.system('bin/igblastn \
              -germline_db_V AIRR_references/%s \
              -germline_db_J AIRR_references/%s \
              -germline_db_D AIRR_references/%s \
              -organism human \
              -query %s \
              -auxiliary_data optional_file/human_gl.aux \
              -show_translation \
              -outfmt 19 \
              > %s' % (V,J,D,fasta,fasta+'.table'))
    os.chdir(current_dir)
    os.system('python3 /mnt/memorycore1/10X_PBMCs/isoform_scripts/convert_igblast_results_to_table_for_R2C2.py \
               %s \
               yes \
               %s' % (fasta+'.table',fasta+'.table.parsed') )


    name_list=[]
    for line in open(fasta+'.table.parsed'):
#        print(line[:70])
        a=line.strip().split('\t')
        name=a[0]
        name_list.append(name)

    print(len(name_list))

    fastq_reads=read_fastq_file(fastq)
    fasta_reads=read_fasta(fasta)

    print(len(name_list),len(fastq_reads))


    new_fastq_reads=[]
    new_fasta_reads={}

    for read in fastq_reads:
#        print(read[0],name_list)
        if read[0] in name_list:
            new_fastq_reads.append(read)

    for read in fasta_reads:
        if read in name_list:
            new_fasta_reads[read]=fasta_reads[read]
    return new_fastq_reads,new_fasta_reads,len(name_list)

temp_folder = sys.argv[2]
subsample = int(sys.argv[3])
config_file= sys.argv[4]

progs = configReader(config_file)
poa = progs['poa']
minimap2 = progs['minimap2']
racon = progs['racon']
consensus = progs['consensus']


def main():
    '''
    spliceSites.py
        psl file
        path
        '0.05'
        annotation
        'g'
        sam file

    defineAndQuantifyIsoforms.py

    createConsensi.py
    '''
    inDir = sys.argv[1]
    V= sys.argv[5]
    D= sys.argv[6]
    J= sys.argv[7]
    type1=sys.argv[8]
    final=open('/mnt/memorycore1/10X_PBMCs/isoform_scripts/'+type1+'.fasta','w')
    fileList = os.listdir(inDir)
#    print(fileList)
#    pslList = sorted([x for x in fileList if x.endswith('.psl') and type1 in x])
    faList = []
    for file in fileList:
        if type1 in file:
            if 'fasta' in file and 'table' not in file:
                faList.append(file)

#    print(len(set(samList)), len(set(faList)), len(set(fqList)))


    cellDict = {}
    for fasta in sorted(faList,key=lambda x: int(x.split('_')[1])):
    
        root = fasta.split('.')[0]
        fastq=root+'.subreads.fastq'
         
#        psl = pslList[i]
        cellDict[root] = [fasta, fastq, root]

    for _, inputs in cellDict.items():
        fasta,fastq, root=inputs[0],inputs[1], inputs[2]
        print(fasta,fastq)
        temp_fastq_reads,temp_fasta_reads,abundance=check_for_complete_air(inDir+'/'+fasta,inDir+'/'+fastq,V,D,J)
#        print(root+'_'+str(abundance),temp_fastq_reads,temp_fasta_reads)
        new_read=make_consensus(root+'_'+str(abundance),temp_fastq_reads,temp_fasta_reads)
#        print(new_read)
        if new_read!='nope':

            final.write(new_read)
    final.close()


if __name__ == '__main__':
    main()
