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
import mappy as mm
import pyabpoa as poa
poa_aligner = poa.msa_aligner(match=5)



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
    possible = set(['poa', 'minimap2', 'gonk', 'consensus', 'racon', 'blat','emtrey','psl2pslx','medaka'])
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
    read_list1 = []
    read_list2 = []
    length = 0
    for line in open(seq_file):
        length += 1
    lineNum = 0
    seq_file_open = open(seq_file, 'r')
    previous = ''
    burn = False
    while lineNum < length:
        name_root = seq_file_open.readline().strip()[1:].split('_')
        name = name_root[0]
        number = name_root[1]
        seq = seq_file_open.readline().strip()
        _ = seq_file_open.readline().strip()
        qual = seq_file_open.readline().strip()

        if previous != name:
            burn = False
            previous = name

        if number == '0':
            burn = True

        seq_length = len(seq)

        if burn or 'I' in number:
            read_list2.append((name, seq, qual, seq_length))
        else:
            read_list1.append((name, seq, qual, seq_length))

        lineNum += 4
    return read_list1, read_list2

def determine_consensus(name, fasta, fastq_reads_full,fastq_reads_partial, counter):
    '''Aligns and returns the consensus'''
    corrected_consensus = ''
    repeats = '0'

    fasta_read_dict=fasta
    fasta_reads = []
    for read, seq in fasta_read_dict.items():
        fasta_reads.append((read, seq))
    repeats = str(len(fasta_reads))

    out_Fq = temp_folder + '/' + counter + '_subsampled.fastq'
    out_F = temp_folder + '/' + counter + '_subsampled.fasta'
    combined_consensus_file = open(temp_folder + '/' + counter + '.fasta', 'w')
    out = open(out_Fq, 'w')

    poa_cons = temp_folder + '/consensus.'+counter+'.fasta'
    output_cons = temp_folder + '/corrected_consensus.'+counter+'.fasta'
    overlap = temp_folder +'/overlaps.'+counter+'.paf'
    overlap_fh=open(overlap,'w')


    fastq_reads = fastq_reads_full + fastq_reads_partial
    if len(fastq_reads) > 0:
        if len(fastq_reads_full) < subsample:
            subsample_fastq_reads = fastq_reads
        else:
            indeces = np.random.choice(
                np.arange(0, len(fastq_reads_full)),
                min(len(fastq_reads_full), subsample), replace=False)
            subsample_fastq_reads = []
            for index in indeces:
                subsample_fastq_reads.append(fastq_reads_full[index])

        subread_counter = 0

        subsample_fastq_reads_numbered=[]
        for read in subsample_fastq_reads:
            subread_counter += 1
            out.write('@' + read[0] + '_' + str(subread_counter) + '\n'
                      + read[1] + '\n+\n' + read[2] + '\n')
            subsample_fastq_reads_numbered.append((read[0] + '_' + str(subread_counter), read[1], read[2], read[3]))
        out.close()
        subsample_fastq_reads=list(subsample_fastq_reads_numbered)


        indeces = np.random.choice(np.arange(0, len(fasta_reads)),
                                   min(len(fasta_reads), 20), replace=False)
        subsample_fasta_reads = []
        for index in indeces:
            subsample_fasta_reads.append(fasta_reads[index])

        first = subsample_fasta_reads[0][1]
        sequences=[]
        mm_align = mm.Aligner(seq=first, preset='map-ont')
        for read,sequence in subsample_fasta_reads:
            for hit in mm_align.map(sequence):
                 if hit.is_primary:
                     if hit.strand==1:
                         sequences.append(sequence)
                     elif hit.strand==-1:
                         sequences.append(mm.revcomp(sequence))


        res = poa_aligner.msa(sequences, out_cons=True, out_msa=False)
        if len(sequences)<=2:
            consensus_sequence = sequences[0]
        elif not res.cons_seq:
            consensus_sequence = sequences[0]
        else:
            consensus_sequence = res.cons_seq[0]

        out_cons_file = open(poa_cons, 'w')
        out_cons_file.write('>Consensus\n' + consensus_sequence + '\n')
        out_cons_file.close()


        final=poa_cons
        mm_align = mm.Aligner(seq=consensus_sequence, preset='map-ont')
        for name,sequence,q,le in subsample_fastq_reads:
            for hit in mm_align.map(sequence):
                if hit.is_primary:
                    overlap_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        name, str(len(sequence)), hit.q_st, hit.q_en,
                        hit.strand, 'Consensus', hit.ctg_len, hit.r_st,
                        hit.r_en, hit.mlen, hit.blen, hit.mapq))

        overlap_fh.close()

        os.system('%s -q 5 -t 1 --no-trimming %s %s %s >%s 2>./racon_messages.txt' \
                   %(racon,out_Fq, overlap, poa_cons, output_cons))
        final=output_cons

        reads = read_fasta(final)
        if len(reads)==0:
            print('racon no')
            reads = read_fasta(poa_cons)

        forMedaka = open(output_cons,'w')
        for read in reads:
            corrected_consensus = reads[read]
            forMedaka.write('>Corrected_Consensus\n'+corrected_consensus+'\n')
        forMedaka.close()

        os.system('mkdir ' + temp_folder + '/' + counter)
        os.system('%s -f -i %s -d %s -o %s > %s_medaka_messages.txt 2>&1'
                  % (medaka, out_Fq, final,
                     temp_folder + '/' + counter, temp_folder + '/' + counter))
        final = temp_folder + '/' + counter + '/consensus.fasta'
        reads = read_fasta(final)
        for read in reads:
            corrected_consensus = reads[read]  # if no read in file, corrected_consensus from racon output is used implicitly
        return corrected_consensus


def make_consensus(name, fastq_reads_full, fastq_reads_partial, fasta_reads):
     if len(fasta_reads)>0:
        corrected_consensus=determine_consensus(name,fasta_reads,fastq_reads_full,fastq_reads_partial,'1')
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

def check_for_complete_air(fasta,fastq,V,D,J,igblast_folder):
    current_dir=os.getcwd()
    fasta=fasta
    os.chdir(igblast_folder)
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
    os.system('python3 convert_igblast_results_to_table_for_R2C2.py \
               %s \
               yes \
               %s' % (fasta+'.table',fasta+'.table.parsed') )


    name_list=[]
    for line in open(fasta+'.table.parsed'):
        a=line.strip().split('\t')
        name=a[0]
        name_list.append(name.split('_')[0])


    fastq_reads_full, fastq_reads_partial =read_fastq_file(fastq)
    fasta_reads=read_fasta(fasta)



    new_fastq_reads_full=[]
    new_fastq_reads_partial=[]
    new_fasta_reads={}

    for read in fastq_reads_full:
        if read[0].split('_')[0] in name_list:
            new_fastq_reads_full.append(read)

    for read in fastq_reads_partial:
        if read[0].split('_')[0] in name_list:
            new_fastq_reads_partial.append(read)

    for read in fasta_reads:
        if read.split('_')[0] in name_list:
            new_fasta_reads[read]=fasta_reads[read]
    return new_fastq_reads_full,new_fastq_reads_partial,new_fasta_reads,len(name_list)

temp_folder = sys.argv[2]
subsample = int(sys.argv[3])
config_file= sys.argv[4]

progs = configReader(config_file)
poa = progs['poa']
minimap2 = progs['minimap2']
racon = progs['racon']
consensus = progs['consensus']
medaka = progs['medaka']
print(medaka)

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
    igblast_folder=sys.argv[9]
    output_folder=sys.argv[10]
    final=open(output_folder+'/'+type1+'.fasta','w')
    fileList = os.listdir(inDir)
#    print(fileList)
#    pslList = sorted([x for x in fileList if x.endswith('.psl') and type1 in x])
    faList = []
    for file in fileList:
        if type1 in file:
            if 'fasta' in file and 'table' not in file:
                faList.append(file)

    cellDict = {}
    for fasta in sorted(faList,key=lambda x: int(x.split('_')[1])):
        root = fasta.split('.')[0]
        fastq=root+'.subreads.fastq'
        cellDict[root] = [fasta, fastq, root]

    for _, inputs in cellDict.items():
        fasta,fastq, root=inputs[0],inputs[1], inputs[2]
        print(fasta,fastq)
        temp_fastq_reads_full,temp_fastq_reads_partial,temp_fasta_reads,abundance=check_for_complete_air(inDir+'/'+fasta,inDir+'/'+fastq,V,D,J,igblast_folder)
        new_read=make_consensus(root+'_'+str(abundance),temp_fastq_reads_full, temp_fastq_reads_partial, temp_fasta_reads)
        if new_read!='nope':

            final.write(new_read)
    final.close()
    check_for_complete_air(output_folder+'/'+type1+'.fasta',output_folder+'/'+type1+'.fasta',V,D,J,igblast_folder)
    os.system('sort '+output_folder+'/'+type1+'.fasta.table.parsed > '+output_folder+'/'+type1+'.fasta.table.parsed.sorted')

if __name__ == '__main__':
    main()
