#include "parseargs.h"
#include "tclap/CmdLine.h"
#include "refio.h"

#include <iostream>

using namespace std;
using namespace TCLAP;


/* parsing arguments */
int parseArguments(int argc, char* argv[],CallingArgs& args){
  try{
    CmdLine cmd("Calling mutations in BAM file.",' ',"0.1");
    
    ValueArg<int> mutspanarg("m","mut_span","The minimum distance of the mutation to the beginning (end) of the read. Default 4.",false,4,"mut_span");
    cmd.add(mutspanarg);
    ValueArg<int> minreadarg("i","min_read","The minimum read count for the mutation to output. Default 1.",false,1,"min_read");
    cmd.add(minreadarg);
    SwitchArg hasnarg("n","output_n","Treat the character N as substitutions.");
    cmd.add(hasnarg);

    SwitchArg indelarg("d","with_indel","Do not skip indels. By default all indels are skipped as most RNA-Seq are performed by Illumina sequencing, which is prone to indel errors. This option will be ignored if -r/--ref_fasta option is provided.");
    cmd.add(indelarg);

    SwitchArg indelreadarg("k","with_indel_read","Do not skip reads with indels. By default all reads with indels are skipped as most RNA-Seq are performed by Illumina sequencing, which is prone to indel errors.");
    cmd.add(indelreadarg);

    ValueArg<string> mutfilearg("l","mutation_list","The text file of a given, sorted list of mutations. Each line in a file records one mutations, with chromosome, location, reference and alternative sequence (separated by tab). The output will only include mutations within a given mutation list.",false,"","mutation_list");
    cmd.add(mutfilearg);

    ValueArg<string> reffilearg("r","ref_fasta","The (optional) fasta file for the reference genome. When this option is set, -d/--with_indel option will be ignored.",false,"","ref_fasta");
    cmd.add(reffilearg);

    SwitchArg mdtagarg("t","use_mdtag","Use MD Tag to call mutations instead of using reference genome (by -r/--ref_fasta option). This option is automatically set if the reference genome is not provided, and requires the BAM file contains the MD tag.");
    cmd.add(mdtagarg);

    UnlabeledValueArg<string> bamfile("bamfile","The bam file from which mutation will be called",true,"","bam_file");
    cmd.add(bamfile);
    
    // parse
    cmd.parse(argc,argv);
    
    args.mut_span=mutspanarg.getValue();
    args.bamfilename=bamfile.getValue();
    args.mutfile=mutfilearg.getValue();
    args.min_read=minreadarg.getValue();
    args.printn=hasnarg.getValue();
    args.skipindel=!indelarg.getValue();
    args.skipindelread=!indelreadarg.getValue();
    args.usemdtag=mdtagarg.getValue();
    
    args.ref_fasta=reffilearg.getValue();
    if(args.ref_fasta != ""){
      args.has_fasta = true;
      if(refseq_init(args.ref_fasta)==-1) return -1;
      args.skipindel=true;
    }
    else{
      args.has_fasta = false;
      args.usemdtag=true;
    }

    cerr<<"BAM file:"<<args.bamfilename<<endl;
    if(args.mutfile=="") args.mutation_given=false;
    else args.mutation_given=true;
    cerr<<"Mut span:"<<args.mut_span<<endl;
    cerr<<"Mutation list:"<<args.mutfile<<endl;
    cerr<<"Min read:"<<args.min_read<<endl;
    cerr<<"Reference genome:"<<args.ref_fasta<<endl;
  }catch(ArgException &e){
    cerr<<"error: "<<e.error()<<" for arg "<<e.argId();
    return -1;
  }
  return 0;
}


