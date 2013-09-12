#!/bin/env python3
'''
Merge VCF files into a big table (or VCF) file.
'''
import sys;
import re;
import argparse;

parser=argparse.ArgumentParser(description="merge VCF files into a big VCF file, preserving all information.");
parser.add_argument('-l','--label',help='The labels of each sample, separated by the comma');
parser.add_argument('-x','--region',help='Only output mutations falling into a specific region, for example chr11:1-10000');
parser.add_argument('-m','--merged',action='store_true',help='Merge forward and reverse fields into one');
parser.add_argument('-v','--output-vcf',action='store_true',help='Output VCF formatted files');
parser.add_argument('-r','--min-read',type=int,default=4,help='Minimum read requirement for ALT reads. Default 4.');
parser.add_argument('vcffiles',nargs=argparse.REMAINDER,help='VCF file names');

args=parser.parse_args();

#filenames=sys.argv[1:];
filenames=args.vcffiles;
nfile=len(filenames);

# process labels
if args.label is None:
  labels=filenames;
else:
  labels=args.label.split(',');
  if len(labels)!=nfile:
    print('Error: the size of the labels must equal the number of files provided.',file=sys.stderr);
    sys.exit(-1);

# process regions
hasregion=False;
regionchr='';
regionstart=0;
regionend=0;
if args.region is not None:
  hasregion=True;
  rgpattern=re.findall('(\w+):(\d+)-(\d+)',args.region);
  if len(rgpattern)==0:
    print('Error: unknown region '+args.region,file=sys.stderr);
    sys.exit(-1);
  regionchr=rgpattern[0][0];
  regionstart=int(rgpattern[0][1]);
  regionend=int(rgpattern[0][2]);

nfilehandle=[];


# printing labels
print('#chrom\tpos\tref\talt',end='');
for nitem in range(nfile):
  vcffile=filenames[nitem];
  label=labels[nitem];
  fhdl=open(vcffile);
  nfilehandle.append(fhdl);
  if args.merged:
    print('\t'+label+'.ref\t'+label+'.alt',end='');
  else:
    print('\t'+label+'.reff\t'+label+'.refv\t'+label+'.altf\t'+label+'.altv',end='');
print('');

while True:
  alltheline=[];
  endoffile=False;
  for fhdl in nfilehandle:
    tline=fhdl.readline();
    if len(tline)==0:
      endoffile=True;
      break;
    tline=tline.strip();
    alltheline.append(tline);
  if endoffile==True:
    break;
  # print header
  flinefld=alltheline[0].split();
  # process regions
  if hasregion:
    if regionchr != flinefld[0]:
      continue;
    if regionstart > int(flinefld[1]) or regionend < int(flinefld[1]):
      continue;
  # check header
  for lines in alltheline:
    fld=lines.split();
    if fld[0:4] != flinefld[0:4]:
      print('Error: unmatched keys ',file=sys.stderr);
      sys.exit(-1);
  # filter for # of reference reads
  ntotalref=0;
  for lines in alltheline:
    fld=lines.split();
    reff=int(fld[6]);
    refv=int(fld[7]);
    ntotalref=ntotalref+reff+refv;
  if ntotalref<args.min_read:
    continue;
  # printing
  if args.output_vcf==True:
    print('\t'.join([flinefld[0],flinefld[1],".",flinefld[2],flinefld[3],"1.0"]),end='\t');
  else:
    print('\t'.join(flinefld[0:4]),end='');
  thislabel=0;
  for lines in alltheline:
    fld=lines.split();
    if args.output_vcf==True:
      print(labels[thislabel]+".DP4="+','.join(fld[4:8]),end='');
      print(';',end='');
    else:
      if args.merged:
        ref=int(fld[4])+int(fld[5]);
        alt=int(fld[6])+int(fld[7]);
        print('\t'+str(ref)+'\t'+str(alt),end='');
      else:
        print('\t'+'\t'.join(fld[4:8]),end='');
    thislabel=thislabel+1;
  print();

for fhd in nfilehandle:
  fhd.close();
