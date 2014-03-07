#!/bin/env python3
'''
Filter mutations based on the various conditions defined by user
'''

from __future__ import print_function
import sys,math;
import argparse;
import re;

parser=argparse.ArgumentParser(description="filter mutations in samples.");

sampledef=parser.add_argument_group('Sample definitions');
sampledef.add_argument('-c','--control',default='',help='The index of control group samples, separated by comma. For example, 0,2,4 defines a control group of 3 samples: the 1st, the 3rd and the 5th sample in the order of the input table. Default: empty (do not use any control samples)');
sampledef.add_argument('-i','--treatment',default='',help='The index of treatment group samples, separated by comma. Default: complement of control samples (if -c/--control option is not specified, use all samples as treatment samples).');
sampledef.add_argument('-t','--min-recurrent',type=int,default=1,help='Print mutations only occuring in at least this number of good treatment samples, defined as those with mutation >=min-recread reads and >=min-recfrac percent frequency. Default 1.');
sampledef.add_argument('-d','--min-recread',type=int,default=10,help='Minimum alt reads defined in treated good samples, default 10.');
sampledef.add_argument('-f','--min-recfrac',type=float,default=0.2,help='Minimum alt reads frequency in treated good samples, default 0.2.');
sampledef.add_argument('-a','--min-ref',type=int,default=4,help='Minimum reference reads in control samples. Default 4.');
sampledef.add_argument('-b','--max-alt',type=int,default=4,help='Maximum alternative reads in control samples. Default 4.');

outputgroup=parser.add_argument_group('Input/output options');
outputgroup.add_argument('-p','--passall',action='store_true',help='Do not do any filtering');
outputgroup.add_argument('-n','--no-header',action='store_true',help='Do not print header and script used');
outputgroup.add_argument('-z','--no-vcf',action='store_true',help='Do not print in vcf format; print as it is.');
outputgroup.add_argument('-l','--labels',default='',help='Labels used for each sample, separated by comma. Default: SAMPLE_x where x is the sample ID. The number of lables MUST be identical to the number of all samples in the original table, not only those defined by the -c/-i parameter.');
outputgroup.add_argument('--DP2',action='store_true',help='DP2 field instad of DP4 field is used in both input/output files');
outputgroup.add_argument('--DP2-out',action='store_true',help='DP2 field instad of DP4 field is used in output files. This option is automatically set true if --DP2 is specified.');
outputgroup.add_argument('-x','--region',help='Only output mutations falling into a specific region, for example chr11:1-10000');


args=parser.parse_args();

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

# process samples
# control and treatment groups
nsample=-1;
ctrlsample=[];
treatmentsample=[];
if args.control != '':
  ctrlsample= [int(x) for x in args.control.split(',')];
if args.treatment != '':
  treatmentsample= [int(x) for x in args.treatment.split(',')];
# sample labels
# allsamplestr="CONTROL CLR4602 CLR4603 CLR4604 CLR4575 CLR4576 CLR4597 CLR4599 CLR4600 CLR4577 CLR4579 CLR4585 CLR4581 CLR4583 CLR4589 CLR4591 CLR4593 CLR4595 HR-RCT-1 HR-RCT-2 HR-RCT-3 HR-RCT-4 HR-RCT-5 HR-RCT-6 HR-RCT-7 HR-RCT-8 HR-RCT-9 HR-RCT-10 HR-RCT-11 HR-RCT-12 CLR4580 CLR4582 CLR4584 CLR4586 CLR4588 CLR4590 CLR4592 CLR4594 CLR4596 CLR4598"
allsamplelabel=[];
if args.labels != '':
  allsamplelabel=args.labels.split(',');

allusedsample=[];

if args.no_header==False:
  print('# parameters:'+' '.join(sys.argv));
  print('# scripts used:');
  for line in open(sys.argv[0]):
    print('# '+line.strip());

nline=0;
for line in sys.stdin:
  if line[0]=='#':
    print(line,end='');
    continue;
  nline=nline+1;
  if nline % 100000 ==1:
    fd=line.strip().split();
    print(str(nline)+' '+fd[0]+' '+fd[1],file=sys.stderr);
  fd=line.strip().split();
  # get the information of control and treatment only using the 1st line
  if nline==1:
    fdi=[int(x) for x in fd[4:]];
    if args.DP2:
      nsample= int(len(fdi)/2);
    else:
      nsample= int(len(fdi)/4);
    print('Number of samples:'+str(nsample),file=sys.stderr);
    # check sample labels
    if args.labels != '':
      if len(allsamplelabel) != nsample:
        print('Error: the number of samples defined in the --labels parameter does not correspond to the columns in the table.',file=sys.stderr);
        sys.exit(-1);
    else:
      allsamplelabel=['SAMPLE_'+str(x) for x in range(nsample)];
    # check control groups and treatment groups
    if args.treatment == '':
      treatmentsample=[x for x in range(nsample) if x not in ctrlsample];
    print('CONTROL group definition:'+','.join([str(x) for x in ctrlsample]),file=sys.stderr);
    print('TREATMENT group definition:'+','.join([str(x) for x in treatmentsample]),file=sys.stderr);
    allusedsample=ctrlsample + treatmentsample;
    # redefine allsamplelabel
    allsamplelabel=[allsamplelabel[x] for x in allusedsample];
  # region
  if hasregion:
    if regionchr != fd[0]:
      continue;
    cpos=int(fd[1]);
    if cpos < regionstart or cpos > regionend:
      continue;
  fdi=[int(x) for x in fd[4:]];
  fdi_b=[];
  fdi2=[];
  for x in allusedsample:
    if args.DP2:
      fdi_b+=[fdi[2*x], fdi[2*x+1]];
      fdi2+=[fdi[2*x],fdi[2*x+1]];
    else:
      fdi2+=[fdi[4*x]+fdi[4*x+1],fdi[4*x+2]+fdi[4*x+3]];
      if args.DP2_out:
        fdi_b+=[fdi[4*x]+fdi[4*x+1],fdi[4*x+2]+fdi[4*x+3]];
      else:
        fdi_b+=[fdi[4*x],fdi[4*x+1],fdi[4*x+2],fdi[4*x+3]];
  fdiref=[fdi2[2*x] for x in range(int(len(fdi2)/2))];
  fdialt=[fdi2[2*x+1] for x in range(int(len(fdi2)/2))];
  fdisum=[fdiref[x]+fdialt[x] for x in range(len(fdiref))];
  for k in range(len(fdisum)):
    if fdisum[k]==0:
      fdisum[k]=1;
  fdifrac=[fdialt[i]*1.0/fdisum[i] for i in range(len(fdisum))];
  fdigood=[( (fdifrac[i]>=args.min_recfrac) and  (fdialt[i]>=args.min_recread))*1 for i in range(len(fdisum))];
  #print(fdiref);
  #print(fdialt);
  #print([int(x*100)/100.0 for x in fdifrac]);
  #print(fdigood);
  # print only if it occurs in recurrent tumors
  if args.passall==False:
    if sum(fdialt[:len(ctrlsample)])>args.max_alt: 
      continue;
    # exclude those with too few read coverages
    if sum(fdiref[:len(ctrlsample)])<args.min_ref:
      continue;
    if sum(fdigood[len(ctrlsample):])<args.min_recurrent:
      continue;
  # print only if at least 2 recurrent samples with read count >=10 and frequency >=10%
  # print in VCF file
  if args.no_vcf==False:
    print('\t'.join([fd[0],fd[1],".",fd[2],fd[3],"1.0"]),end='\t');
    for i in range(len(allsamplelabel)):
      if args.DP2 or args.DP2_out:
        print(allsamplelabel[i]+".DP2="+','.join([str(x) for x in fdi_b[(2*i):(2*i+2)]]),end='');
        print(';',end='');
      else:
        print(allsamplelabel[i]+".DP4="+','.join([str(x) for x in fdi_b[(4*i):(4*i+4)]]),end='');
        print(';',end='');
    print();
  else:
    print('\t'.join([fd[0],fd[1],fd[2],fd[3]]),end='\t');
    print('\t'.join([str(x) for x in fdi_b]));

