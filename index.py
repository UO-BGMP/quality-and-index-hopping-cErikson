# -*- coding: utf-8 -*-
#!/bin/python
"""
Take dual index paired end reads, that are in 4 seprate FASTQ(read1, index1, index2, read2). And a TSV of used barcodes in col 2.
ASSUMES_ALL_FASTQ_FILES_ARE_SORTED_THE_SAME
Filter out entries whos indexs are below a qual score, then seperates reads into new files by index pair. 
"""
###### Imports #####
import argparse as arg
import os.path
import json
import gzip
import pdb
import sys
#from itertools import izip #Python2 compatibility

testing=True

##### Args #####
if __name__ == "__main__" and testing != True:
    parser = arg.ArgumentParser()
    # Positional mandatory arguments
    parser.add_argument('-1', "--read1_file", help="read 1 File", type=str)
    parser.add_argument('-2', "--index1_file", help="index 1 File", type=str)
    parser.add_argument('-3', "--index2_file", help="index 2 File", type=str)
    parser.add_argument('-4', "--read2_file", help="read 2 File", type=str)
    parser.add_argument('-k', "--known_file", help="Known indexes File", type=str)
    parser.add_argument('-i', "--index_col", help="Column in Known indexes File that contains barcodes ", type=int)
    parser.add_argument('-o', "--out_dir", help="output dir", type=str)
    parser.add_argument("-c", "--cutoff", help="cutoff", type=int)
    parser.add_argument("-g", "--gzip", help="Are files zipped", default=False, action='store_true')
    parser.add_argument("-r", "--rev_comp", help="Do R2 & I2 need to be rev comp?", default=False, action='store_true')
    
    # Parse arguments
    ARGS = parser.parse_args()
    
    ##### Testing #####
else:   # Else test
    sys.stderr.writelines("!!!!!___RUNNING_IN_TESTING_MODE_WITH_TEST_ARGS___!!!!!\n")
    class test_args(object):
        def __init__(self):
            self.read1_file='/home/christian/gdrive/School/bi622/index_hopping/r1.gz'
            self.index1_file='/home/christian/gdrive/School/bi622/index_hopping/r2.gz'
            self.index2_file='/home/christian/gdrive/School/bi622/index_hopping/r3.gz'
            self.read2_file='/home/christian/gdrive/School/bi622/index_hopping/r4.gz'
            self.known_file='/home/christian/gdrive/School/bi622/index_hopping/indexes.txt'
            self.index_col=5
            self.out_dir='/home/christian/gdrive/School/bi622/index_hopping/counts_test.json'
            self.cutoff=1
            self.gzip=True
            self.rev_comp=True
    ARGS = test_args()


##### Defs #####    
def yield_fastq(fastq, gz=False, rev_comp=False):
    '''Make a generator that yields (seq_header, seq, score) for each entry)'''
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X', 'M':'K', 'R':'Y', 'W':'W', 'S':'S', 'Y':'R', 'K':'M','V':'B', 'H':'D', 'D':'H', 'B':'V'}
    if gz==False:
        fhs=open(fastq, 'r')
    else:
        fhs=gzip.open(fastq, 'rt')
    l=fhs.readline().strip()
    while l != '':
        if l.startswith('@'):   #if header
            header=l  #save header
            if rev_comp is False:
                seq=fhs.readline().strip().upper()
            else:
                seq=''
                for base in fhs.readline().strip().upper()[::-1]:
                    seq+=complement[base]
            fhs.readline()    # pass over qual head
            score=[ord(a)-33 for a in fhs.readline().strip()] if rev_comp is False else [ord(a)-33 for a in fhs.readline().strip()[::-1]]# read phred 
            yield header, seq, score # yield data
        l=fhs.readline().strip()

class qual_info():
    def __init__(self):
        self.run_avg=0
        self.rec_num=0
        self.read_avg_dist={}
        self.pos_avg_dist=[0.0]*10000
    
    def add(self, phred):
        read_sum=round(sum(phred)/len(phred)) 
        self.rec_num+=1
        self.read_avg_dist.setdefault(read_sum, 0)
        self.read_avg_dist[read_sum] += 1
        self.pos_avg_dist=[o_pos+((n_pos-o_pos)/self.rec_num) for n_pos, o_pos in zip(phred, self.pos_avg_dist)]
        self.run_avg=self.run_avg+((read_sum-self.run_avg)/self.rec_num)
        
    
###### MAIN ######
#Read in known codes make all possible combos
code_set={'swaped':[], 'known':[]}
codes=set()
with open(ARGS.known_file,'r') as kf:
    for line in kf:
        if line.startswith('#') == False:
            codes.add(line.split('\t')[ARGS.index_col-1].strip())
for c1 in codes:
    for c2 in codes:
        if c1!=c2:
            code_set['swaped'].append(c1+'_'+c2) 
        else:
            code_set['known'].append(c1+'_'+c2) 

# make counts file
count_file={'tag_counts':{}}
counts = count_file['tag_counts']
    
# Create Qual_info classes
qual_book=[qual_info() for x in range(4)]

# Filter, seperate, and write records
for r1, i1, i2, r2 in zip(yield_fastq(ARGS.read1_file, ARGS.gzip), yield_fastq(ARGS.index1_file, ARGS.gzip), yield_fastq(ARGS.index2_file, ARGS.gzip, ARGS.rev_comp), yield_fastq(ARGS.read2_file, ARGS.gzip, ARGS.rev_comp)): # for each input fastq
    if all(x >= ARGS.cutoff for x in i1[2]+i2[2]): # if all the phred scores pass cutoff
        for chap, phred in zip(qual_book, (r1[2],i1[2],i2[2],r2[2])): # for each qual klass 
            chap.add(phred) #add the phred info to chapter
        counts.setdefault(i1[1]+'_'+i2[1],0) # set default tag in json
        counts[i1[1]+'_'+i2[1]] += 1    # increment index count    
        
# write JSON 
count_file['qual_info']=[vars(chap) for chap in qual_book]
count_file['metadata']=(vars(ARGS))
count_file['codes']=code_set
with open(ARGS.out_dir, 'w+') as fp:
        json.dump(count_file, fp, indent=2) 

