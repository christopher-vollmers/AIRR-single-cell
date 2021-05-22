import os
import sys
import argparse

def argParser():
    '''Parses arguments.'''
    parser = argparse.ArgumentParser(description = 'Makes consensus sequences \
                                                    from R2C2 reads.',
                                     add_help = True,
                                     prefix_chars = '-')
    required = parser.add_argument_group('required arguments')
    required.add_argument('--input_folder', '-i', type=str, action='store', required=True,
                          help='This folder should contain the consensus read fasta and subread fastq files for all the single cells to be analyzed')
    required.add_argument('--sam_folder', '-s', type=str, action='store', required=True,
                          help='This folder should contain the sam files for all the single cells to be analyzed')
    parser.add_argument('--output_folder', '-o', type=str, action='store',
                        help='Output folders and files will be generate here')
    parser.add_argument('--isotype_positions', '-t', type=str, action='store',
                        help='File storing genome coordinates for IGH isotypes/isoforms')
    parser.add_argument('--igblast', '-b', type=str, action='store',
                    help='Path to igblast. This should point to the folder that contains the bin folder, not the bin folder itself')
    parser.add_argument('--config', '-c', type=str, action='store', default='',
                        help='If you want to use a config file to specify paths to\
                              programs, specify them here. Use for poa, racon, gonk,\
                              blat, and minimap2 if they are not in your path.')

    return vars(parser.parse_args())

def makeFolders():
    os.system('rm -r %s/TCR' %(output_folder))
    os.system('rm -r %s/IGH' %(output_folder))
    os.system('rm -r %s/IGL' %(output_folder))
    os.system('mkdir %s/TCR' %(output_folder))
    os.system('mkdir %s/IGH' %(output_folder))
    os.system('mkdir %s/IGH/isotypes' %(output_folder))
    os.system('mkdir %s/IGL' %(output_folder))
    os.system('mkdir %s/IGL/temp' %(output_folder))
    os.system('mkdir %s/IGH/isotypes/temp' %(output_folder))
    os.system('mkdir %s/TCR/temp' %(output_folder))

def processIGL():
    os.system('python3 igl_filter.py %s %s/IGL' %(sam_folder,output_folder))
    os.system('python3 clean.py %s/IGL'  %(output_folder))
    os.system('python3 filter_fasta.py %s %s/IGL' %(input_folder,output_folder))
    os.system('python3 IGLWrapper_simple.py %s/IGL %s/IGL/temp 500 %s IGKV_processed IGHD_processed IGKJ_processed IGLK %s %s ' %(output_folder,output_folder,config_file,igBlastPath,output_folder))
    os.system('python3 IGLWrapper_simple.py %s/IGL %s/IGL/temp 500 %s  IGLV_processed IGHD_processed IGLJ_processed IGLL %s %s' %(output_folder,output_folder,config_file,igBlastPath,output_folder))


def processTCR():
    os.system('python3 tcr_filter.py %s %s/TCR' %(sam_folder,output_folder))
    os.system('python3 clean.py %s/TCR'  %(output_folder))
    os.system('python3 filter_fasta.py %s %s/TCR' %(input_folder,output_folder))
    os.system('python3 TCRWrapper_simple.py %s/TCR %s/TCR/temp 500 %s TRAV_processed TRBD_processed TRAJ_processed TCRA %s ' %(output_folder,output_folder,config_file,igBlastPath))
    os.system('python3 TCRWrapper_updated.py %s/TCR %s/TCR/temp 500 %s TRBV_processed TRBD_processed TRBJ_processed TCRB %s ' %(output_folder,output_folder,config_file,igBlastPath))

def processIGH():
    os.system('python3 antibody_filter.py %s %s/IGH' %(sam_folder,output_folder))
    os.system('python3 clean.py %s/IGH'  %(output_folder))
    os.system('python3 cell_sam2psl.py %s/IGH' %(output_folder))
    os.system('python3 identifyIsotypes.py %s %s/IGH' %(isotype_positions,output_folder))
    os.system('python3 filter_fa_psl.py %s %s/IGH/isotypes' %(input_folder,output_folder))
    os.system('python3 antibodyIsotypeWrapper_simple.py %s/IGH/isotypes %s/IGH/isotypes 500 %s %s' %(output_folder,output_folder,config_file,igBlastPath))


args = argParser()


input_folder=args['input_folder']
sam_folder=args['sam_folder']
output_folder=args['output_folder']
isotype_positions=args['isotype_positions']
config_file=args['config']
igBlastPath=args['igblast']

def main():
    makeFolders()
    processIGL()
    processTCR()
    processIGH()

main()
