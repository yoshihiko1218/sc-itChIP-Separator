import argparse 
import gzip
import pandas as pd

def calc_hmm(a,b):
    dis = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            dis = dis + 1
    return dis

def find_match_with_hmm(long,short):
    minhmd = 100
    for i in range(len(long)-len(short)+1):
        hmd = calc_hmm(long[i:i+len(short)],short)
        if hmd < minhmd:
            minhmd = hmd
    return minhmd

def Separate_itChip(out_prefix,outfolder,f_filename,r_filename,f_ref,r_ref,start_forward,start_reverse,thresh_forward,thresh_reverse):
    
    print(out_prefix+' start...')
    
    # Read input files
    print('Reading input file...')
    f=gzip.open(f_filename,'rb')
    r=gzip.open(r_filename,'rb')
    f_fastq=f.read()
    r_fastq=r.read()
    f_list = f_fastq.decode().strip().split("\n")
    r_list = r_fastq.decode().strip().split("\n")

    # Match each forward strand to one of the forward index.
    print('Matching forward file...')
    f_barcodes = {}
    f_nomatch = []
    for i in f_ref:
        f_barcodes[i] = []
    f_dupmatchcount = 0
    f_nomatchcount = 0
    f_matchcount = 0
    for i in range(int(len(f_list)/4)):
        match = f_list[i*4+1][(start_forward-1):(start_forward+1+len(f_ref[0]))]
        matchcount = 0
        matchref = ''
        for ref in f_ref:
            if find_match_with_hmm(match,ref) <= thresh_forward:
                matchref = ref
                matchcount += 1
        if matchcount == 0:
            f_nomatch.append(f_list[i*4+1])
            f_nomatchcount += 1
        elif matchcount == 1:
            f_barcodes[matchref].append(i)
            f_matchcount += 1
        elif matchcount > 1:
            f_dupmatchcount += 1

    # Match each reverse strand to one of the reverse index.
    print('Matching reverse file...')
    r_barcodes = {}
    r_nomatch = []
    for i in r_ref:
        r_barcodes[i] = []
    r_dupmatchcount = 0
    r_nomatchcount = 0
    r_matchcount = 0
    for i in range(int(len(r_list)/4)):
        match = r_list[i*4+1][(start_reverse-1):(start_reverse+1+len(r_ref[0]))]
        matchcount = 0
        matchref = ''
        for ref in r_ref:
            if find_match_with_hmm(match,ref) <= thresh_reverse:
                matchref = ref
                matchcount += 1
        if matchcount == 0:
            r_nomatch.append(r_list[i*4+1])
            r_nomatchcount += 1
        elif matchcount == 1:
            r_barcodes[matchref].append(i)
            r_matchcount += 1
        elif matchcount > 1:
            r_dupmatchcount += 1

    # Make match of forward and reverse indexes.
    print('Matching forward and inverse and creating report...')
    barcodes = {}
    for i in f_ref:
        for j in r_ref:
            barcodes[i+j] = list(set(f_barcodes[i])&set(r_barcodes[j]))

    # Creat report file
    report = '\t#Match\t#DuplicatedMatch\t#Mismatch\n'
    report = report + 'Forward\t%d(%.1f%%)\t%d(%.1f%%)\t%d(%.1f%%)\n'%(f_matchcount,(f_matchcount/(len(f_list)/4)*100),
                                                                 f_dupmatchcount,(f_dupmatchcount/(len(f_list)/4)*100),
                                                                 f_nomatchcount,(f_nomatchcount/(len(f_list)/4)*100))
    report = report + 'Reverse\t%d(%.1f%%)\t%d(%.1f%%)\t%d(%.1f%%)'%(r_matchcount,(r_matchcount/(len(r_list)/4)*100),
                                                                 r_dupmatchcount,(r_dupmatchcount/(len(r_list)/4)*100),
                                                                 r_nomatchcount,(r_nomatchcount/(len(r_list)/4)*100))


    with(open(outfolder+'/'+out_prefix+'report.txt','w')) as f:
        f.write(report)

    # Write fastq output file
    print('Creating output fastq...')
    for key in barcodes.keys():
        indexes = barcodes[key]
        allindexes = []
        for i in indexes:
            allindexes.append(4*i)
            allindexes.append(4*i+1)
            allindexes.append(4*i+2)
            allindexes.append(4*i+3)

        pd.DataFrame(f_list).iloc[allindexes].to_csv(outfolder+'/'+out_prefix+'_'+key+'_R1.fastq',header=False,index=False,sep='\t')
        pd.DataFrame(r_list).iloc[allindexes].to_csv(outfolder+'/'+out_prefix+'_'+key+'_R2.fastq',header=False,index=False,sep='\t')


def main():
    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    parser._action_groups.append(optional)
    required.add_argument("-f", "--forward", dest="f_filename",
                        help="forward fq.gz file", metavar="F_FILE",required=True)
    required.add_argument("-r", "--reverse", dest="r_filename",
                        help="reverse fq.gz file", metavar="R_FILE", required=True)
    required.add_argument("-p", "--prefix", dest="out_prefix",
                        help="prefix of output file", metavar="PREFIX", required=True)
    required.add_argument("-o", "--outfolder", dest="outfolder",
                        help="output folder", metavar="OUT_FOLDER", required=True)

    optional.add_argument("-if", "--index_forward", dest="index_forward",
                        help="input file for list of forward index", metavar="F_INDEX_FILE")
    optional.add_argument("-ir", "--index_reverse", dest="index_reverse",
                        help="input file for list of reverse index", metavar="R_INDEX_FILE")

    optional.add_argument("-sf", "--start_forward", dest="start_forward",
                        help="start position of forward index(0 based). default: 21", metavar="START_FORWARD")
    optional.add_argument("-sr", "--start_reverse", dest="start_reverse",
                        help="start position of reverse index(0 based). default: 27", metavar="START_REVERSE")

    optional.add_argument("-tf", "--thresh_forward", dest="thresh_forward",
                        help="threshold for forward index matching(inclusive). default: 2", metavar="THRESH_FORWARD")
    optional.add_argument("-tr", "--thresh_reverse", dest="thresh_reverse",
                        help="threshold for reverse index matching(inclusive). default: 2", metavar="THRESH_REVERSE")

    optional.add_argument("-rg", "--range", dest="range",
                        help="size of range that check index match. default: 2", metavar="RANGE")

    args = parser.parse_args()
    args = vars(args)

    if args['index_forward'] != None:
        with open(args['index_forward'],'r') as f:
            f_ref = [line.rstrip() for line in f]
    else:
        f_ref = ['TATAGCCT','ATAGAGGC','CCTATCCT','GGCTCTGA','AGGCGAAG','TAATCTTA','CAGGACGT','GTACTGAC']

    if args['index_reverse'] != None:
        with open(args['index_reverse'],'r') as f:
            r_ref = [line.rstrip() for line in f]
    else:
        r_ref = ['CGAGTAAT','TCTCCGGA','AATGAGCG','GGAATCTC','TTCTGAAT','ACGAATTC','AGCTTCAG','GCGCATTA','CATAGCCG','TTCGCGGA',
                'GCGCGAGA','CTATCGCT']

    if args['start_forward'] != None:
        start_forward = int(args['start_forward'])
    else:
        start_forward = 21
    if args['start_reverse'] != None:
        start_reverse = int(args['start_reverse'])
    else:
        start_reverse = 27

    if args['thresh_forward'] != None:
        thresh_forward = int(args['thresh_forward'])
    else:
        thresh_forward = 2
    if args['thresh_reverse'] != None:
        thresh_reverse = int(args['thresh_reverse'])
    else:
        thresh_reverse = 2

    if args['range'] != None:
        start_forward = int(args['range'])
    else:
        start_forward = 1

    Separate_itChip(args['out_prefix'],args['outfolder'],args['f_filename'],args['r_filename'],f_ref,r_ref,start_forward,start_reverse,thresh_forward,thresh_reverse)