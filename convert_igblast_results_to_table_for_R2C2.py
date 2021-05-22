import os
import sys
import re
import editdistance

handle = open(sys.argv[1])
isoform_test=sys.argv[2]
print(isoform_test)
out=open(sys.argv[3],'w')
header = handle.readline()
field_to_index = {}


def determine_isoform(C_seq,match,IsoDict):
    typeDict={}
    typeDict['S']=0
    typeDict['M']=0
    for i in range(0,len(C_seq),3):
        kmer=C_seq[i:i+12]
        if kmer in IsoDict:
            matches=IsoDict[kmer]
            for isoform in matches:
                if match in isoform:
                     type1=isoform.split('_')[1]
                     typeDict[type1]+=1
    if typeDict['S']>typeDict['M']:
         isoform='S'

    elif typeDict['S']<typeDict['M']:
         isoform='M'
    else:
         isoform='-'
    return isoform,typeDict['S'],typeDict['M']

def find_constant_region(C_seq,Isotypes):

    match='_'
    min1=10000000

    for Isotype in Isotypes:

        C_limit=min(len(C_seq),len(Isotype[1]),300)
        if C_limit>8:
            distance1=editdistance.eval(C_seq[:C_limit],Isotype[1][:C_limit])
            if distance1/(C_limit+1)<0.25:
                if distance1<min1:
                    min1=distance1
                    match=Isotype[0]
    return match,min1

def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readList = []
    tempSeqs, headers, sequences = [], [], []
    for line in open(inFile):

        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            headers.append(line.split()[0][1:])
        # covers the case where the file ends while reading sequences
        if line.startswith('>'):
            sequences.append(''.join(tempSeqs).upper())
            tempSeqs = []
        else:
            tempSeqs.append(line)
    sequences.append(''.join(tempSeqs).upper())
    sequences = sequences[1:]
    for i in range(len(headers)):
        readList.append((headers[i],sequences[i]))
    return readList

Isotypes=read_fasta('spliced_constant_regions_edited.fasta')

Isoforms=read_fasta('membrane_secreted_only_human_C')

IsoDict={}
for isoform,sequence in Isoforms:
    if isoform not in IsoDict:
        for i in range(0,120,1):
            if sequence[i:i+12] not in IsoDict:
                IsoDict[sequence[i:i+12]]=[]
            IsoDict[sequence[i:i+12]].append(isoform)



header_list = header.strip().split("\t")
for line in handle:
    la = line.strip().split("\t")
    print(len(la))
    if len(la) >= 80:
        ID = la[0]
        seq = la[1]
        V_seg = la[7].split(",")[0]
        D_seg = la[8].split(",")[0]
        J_seg = la[9].split(",")[0]
        if len(la[7].split(","))>1:
            Amb='Ambiguous:'+la[7]
        else:
            Amb='N/A'
        CDR3 = la[42]
        seq_type = la[2]
        germline_mismatch_pos = []
        sequence_V = la[20]
        germline_V = la[22]
        s_index = int(la[60])
        g_index = int(la[62])
        v_support=float(la[54])
        print(V_seg, v_support)

        try:
            C_seq=seq[int(la[69])-1:-20]
        except:
            C_seq=''
        match,distance=find_constant_region(C_seq,Isotypes)
        if isoform_test=='yes':

            isoform,S,M=determine_isoform(C_seq,match,IsoDict)
            print(match,distance,isoform,S,M)
        else:
            isoform='N/A'


        for s,g in zip(sequence_V,germline_V):
            if s != '-' and g != '-' and s != g:
                germline_mismatch_pos.append(str(g_index))
            if s != "-":
                s_index+=1
            if g != '-':
                g_index+=1
#        out.write("\t".join([re.sub("reversed\|","",ID),seq_type,V_seg,D_seg,J_seg,CDR3,"N/A",isoform,match,seq,",".join(germline_mismatch_pos),'\n']))
        string='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ID.replace('reversed|',''),seq_type,V_seg,D_seg,J_seg,CDR3,Amb,isoform,match,seq,",".join(germline_mismatch_pos))
        if v_support<1e-25:
            out.write(string)
