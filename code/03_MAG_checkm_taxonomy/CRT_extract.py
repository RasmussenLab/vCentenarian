#!/bin/python
import argparse
import os
import sys
import gzip
import re

parser = argparse.ArgumentParser(description='''
   Extract MAG CRISPR Spacers from CRT prediction files 
   ''')
parser.add_argument('-i', help='CRT file')
parser.add_argument('-o', help='Spacer fileout')

def is_allowed_specific_char(string):
    charRe = re.compile(r'[^a-zA-Z0-9.]')
    string = charRe.search(string)
    return not bool(string)

def parse_CRT(crt_file,outfile):
    binid = crt_file.replace('.crt','')
    contig = None
    spacer_counter = 0
    with open(crt_file,'r') as infile,open(outfile,'w') as outfile:
        dash_flag_start = False
        dash_flag_end = True
        for line in infile:

            if 'ORGANISM' in line:
                contig = line.strip().split('ORGANISM:')[1].replace(' ','')
                spacer_counter = 0

            if '--------' in line and dash_flag_end:
                dash_flag_start = True
                dash_flag_end = False
                continue
            elif '--------' in line and dash_flag_start:
                dash_flag_start = False
                dash_flag_end = True
            
            ### Parse Repeat-Spacer Line
            if dash_flag_start:
                line = line.strip().split()
                if len(line) > 2:
                        spacer_counter += 1
                        position = line[0]
                        repeat = line[1]
                        spacer = line[2]
                        if re.search(r"[^ATGC]", spacer):
                            continue
                        if len(spacer) >= 19:
                            header = '>'+binid + '::' + contig + ':CRT:' + str(spacer_counter)
                            outfile.write("{}\n{}\n".format(header,spacer))
                        
if __name__ == "__main__":
    args = parser.parse_args()

    ### Extract Spacers and Arrays from CRT files
    crt_dictionary = parse_CRT(args.i,args.o)
    print('Done')



