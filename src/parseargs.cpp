#include "parseargs.h"
#include "tclap/CmdLine.h"

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

    ValueArg<string> mutfilearg("l","mutation_list","The text file of a given, sorted list of mutations. Each line in a file records one mutations, with chromosome, location, reference and alternative sequence (separated by tab). The output will only include mutations within a given mutation list.",false,"","mutation_list");
    cmd.add(mutfilearg);

    UnlabeledValueArg<string> bamfile("bamfile","The bam file from which mutation will be called",true,"","bam_file");
    cmd.add(bamfile);
    
    // parse
    cmd.parse(argc,argv);
    
    args.mut_span=mutspanarg.getValue();
    args.bamfilename=bamfile.getValue();
    args.mutfile=mutfilearg.getValue();
    args.min_read=minreadarg.getValue();

    cerr<<"BAM file:"<<args.bamfilename<<endl;
    if(args.mutfile=="") args.mutation_given=false;
    else args.mutation_given=true;
    cerr<<"Mut span:"<<args.mut_span<<endl;
    cerr<<"Mutation list:"<<args.mutfile<<endl;
    cerr<<"Min read:"<<args.min_read<<endl;
  }catch(ArgException &e){
    cerr<<"error: "<<e.error()<<" for arg "<<e.argId();
    return -1;
  }
  return 0;
}


