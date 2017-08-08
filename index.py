# -*- coding: utf-8 -*-
#!/bin/python
"""
Take dual index paired end reads, that are in 4 seprate FASTQ(read1, index1, index2, read2). And a TSV of used barcodes in col 2.
ASSUMES_ALL_FASTQ_FILES_ARE_SORTED_THE_SAME
Filter out entries whos indexs are below a qual score, then seperates reads into new files by index pair. 
"""
###### Imports #####
import argparse as arg

##### Args #####
if __name__ == "__main__":
    parser = arg.ArgumentParser()
    # Positional mandatory arguments
    parser.add_argument("read1_file", help="read 1 File", type=str)
    parser.add_argument("index1_file", help="index 1 File", type=str)
    parser.add_argument("index2_file", help="index 2 File", type=str)
    parser.add_argument("read2_file", help="read 2 File", type=str)
    parser.add_argument("known_file", help="Known indexes File", type=str)
    parser.add_argument("out_dir", help="output dir", type=str)
    parser.add_argument("-c", "--cutoff", help="cutoff", type=int)
    # Parse arguments
    ARGS = parser.parse_args()

##### Defs #####    
def yield_fastq(fastq):
    '''Make a generator that yields (seq_header, seq, score) for each entry)'''
    fhs=open(fastq, 'r')
    l=fhs.readline().strip()
    while l != '':
        if l.startswith('@'):   #if header
            header=l  #save header
            seq=fhs.readline().strip() # read seq
            fhs.readline()    # pass over qual head
            score=[ord(a)-33 for a in fhs.readline().strip()] # read phred 
            yield header, seq, score # yield data
        l=fhs.readline().strip()
        
###### MAIN ######
#Read in known codes make all possible combos
known=[]
codes=set()
with open(ARGS.known_file,'r') as kf:
    for line in kf:
        codes.add(line.split('\t')[1].strip())
for c1 in codes:
    for c2 in codes:
        known.append((c1,c2))
known.append(('unknown','index'))
del codes

# make output files
files={}
for code in known:
    files[code]=[open(ARGS.out_dir+code[0]+'_'+code[1]+x+'.fastq', 'w') for x in ('_R1','_I1', '_I2', '_R2')]
del known

# Filter, seperate, and write records
counts={}
for r1, i1, i2, r2 in zip(yield_fastq(ARGS.read1_file), yield_fastq(ARGS.index1_file), yield_fastq(ARGS.index2_file), yield_fastq(ARGS.read2_file)): # for each input fastq
        if all(x >= ARGS.cutoff for x in i1[2] and i2[2]): # if all the phred scores pass cutoff
            counts.setdefault((i1[1],i2[1]),0) # set default index count
            counts[(i1[1],i2[1])]+=1    # increment index count            
            for fh, data in zip(files[(i1[1],i2[1]) if (i1[1], i2[1]) in files else ('unknown','index')], [r1, i1, i2, r2]):   # for each output file and its respective input entry, and determine if known          
                fh.write('{0}_{1}\n{2}\n+\n{3}\n'.format(data[0].split()[0], i1[1]+'_'+i2[1], data[1], ''.join((chr(x+33) for x in data[2])))) # write the fastq
        
# close files                     
for idx in files.values(): 
    for fh in idx: 
        fh.close() 
        