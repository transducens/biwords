#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <ctime>
#include <clocale>
#include <cstdlib>

#include "configure.h"
#include "alignment.h"
#include "zfstream.h"

using namespace std;

void help(char *name) {
  cerr<<"USAGE:\n";
  //cerr<<name<<" --align1 alignments [--align2 alignments] [--offsets] [--gzip] [--non-contiguous] [--debug]\n\n";
  cerr<<name<<" --align alignments [--switch] [--offsets [--encode-offsets]] [--non-contiguous] [--freq-th thheshold] [--gzip] [--debug]\n\n";
  cerr<<"ARGUMENTS: \n"
      //<<"   --align1|-a: Specify a file containing the alignments 1:N (or N:1)\n"
      <<"   --align|-a: Specify a file containing the alignments 1:N (or N:1)\n"
      //<<"   --align2|-l: Specify a file containing the alignments N:1 (or 1:N)\n"
      <<"   --offsets|-o: Use offsets instead of discarding crossing alignments\n"
      <<"   --freq-th|-f: Minimum freq. for a biword to be used\n"
      <<"    --switch|-s: Switch SL and TL parts. To process N:1 alignments as 1:N\n"
      <<"   --encode-offsets|-e: Encode all the offsets in base "<<Alignment::code_base<<" to have a single integer\n"
      <<"   --non-contiguous|-n: Allow non-contigous multiwords\n"
      <<"   --gzip|-z: Tell the program that the alignment files are gziped\n"
      <<"   --help|-h: Show this help\n"
      <<"   --debug|-d: Show debug information\n"
      <<"   --version|-v: Show version information\n";
}


void process_input(istream *fal1, istream *fal2, bool with_offsets, bool non_contiguous, bool encode_offsets, bool rafa, double freqth, bool switch_alg) {
  string al1, al2;

  int nal=0;
  int nal_not_used=0;
  while ((!fal1->eof()) && ((fal2==NULL) || (!fal2->eof()))) {
    getline(*fal1,al1);
    if (fal2!=NULL)
      getline(*fal2,al2);

    if ((al1.length()>0) && ((fal2==NULL) || (al2.length()>0))) {
      nal++;
      Alignment *alig1=NULL, *alig2=NULL;

      alig1 = new Alignment(al1);
      if (fal2!=NULL) 
        alig2= new Alignment(al2);

      string output;
      if (with_offsets) {
        //if (non_contiguous)
        //  output=alig1->biwords_with_offsets_non_contiguous(alig2);
        //else
        //  output=alig1->biwords_with_offsets(alig2);

        if (rafa) 
          output=alig1->biwords_with_offsets_new_version_rafa(non_contiguous, encode_offsets, alig2);
        else {
          output=alig1->biwords_with_offsets_new_version(switch_alg, freqth, non_contiguous, encode_offsets, alig2);
        }

        //output=alig1->biwords_with_offsets_alg1_n_discontinuous(alig2);      OLD
        //output=alig1->biwords_with_offsets_alg1_1();                         OLD
      } else
        output=alig1->biwords_without_offsets(alig2);
      //output=alig1->biwords_without_offsets_alg1_1();

      if (output.length()>0)
        cout<<output;
      else 
        nal_not_used++;

      if ((nal%10000)==0)
        cerr<<nal<<" alignments processed\n";

      delete alig1;
      if(alig2!=NULL)
        delete alig2;
    }
  }

  cerr<<nal<<" alignments processed\n";
  cerr<<nal_not_used<<" alignments not used\n";
}

int main(int argc, char* argv[]) {
  int c;
  int option_index=0;

  string align1_file="";
  string align2_file="";

  bool use_zlib=false;
  bool with_offsets=false;
  bool encode_offsets=false;
  bool debug=false;
  bool non_contiguous=false;

  bool switch_alg=false;
  bool rafa=false;

  double freqth=0.0;

  while (true) {
    static struct option long_options[] =
      {
	//{"align1",    required_argument,  0, 'a'},
        {"align",     required_argument,  0, 'a'},
	{"align2",    required_argument,  0, 'l'},
        {"freq-th",   required_argument,  0, 'f'},
	{"offsets",         no_argument,  0, 'o'},
	{"encode-offsets",  no_argument,  0, 'e'},
	{"non-contiguous",  no_argument,  0, 'n'},
	{"gzip",            no_argument,  0, 'z'},
	{"help",            no_argument,  0, 'h'},
        {"switch",          no_argument,  0, 's'},
        {"rafa",            no_argument,  0, 'r'},
	{"version",         no_argument,  0, 'v'},
	{"debug",           no_argument,  0, 'd'},
	{0, 0, 0, 0}
      };

    c=getopt_long(argc, argv, "a:l:f:oenszhrvd",long_options, &option_index);
    if (c==-1)
      break;
      
    switch (c) {
    case 'a':
      align1_file=optarg;
      break;
    case 'l':
      align2_file=optarg;
      break;
    case 'f':
      freqth=atof(optarg);
      break;
    case 'o':
      with_offsets=true;
      break;
    case 'e':
      encode_offsets=true;
      break;
    case 'n':
      non_contiguous=true;
      break;
    case 'z':
      use_zlib=true;
      break;
    case 'd':
      debug=true;
      break;
    case 's':
      switch_alg=true;
      break;
    case 'r':
      rafa=true;
      break;
    case 'h': 
      help(argv[0]);
      exit(EXIT_SUCCESS);
      break;
    case 'v':
      cerr<<PACKAGE_STRING<<"\n";
      exit(EXIT_SUCCESS);
      break;    
    default:
      help(argv[0]);
      exit(EXIT_FAILURE);
      break;
    }
  }

  if (align1_file=="") {
    cerr<<"Error: No alignment file was given. You need to provide it with the --align1 option\n";
    help(argv[0]);
    exit(EXIT_FAILURE);
  }

  istream *fal1=NULL, *fal2=NULL;
  if (use_zlib) {
    fal1 = new gzifstream(align1_file.c_str());
    if (align2_file != "")
      fal2 = new gzifstream(align2_file.c_str());
  }  else {
    fal1 = new ifstream(align1_file.c_str());
    if (align2_file != "")
      fal2 = new ifstream(align2_file.c_str());
  }

  if (fal1->fail()) {
    cerr<<"Error: Cannot open input file '"<<align1_file<<"'\n";
    delete fal1;
    if (fal2!=NULL)
      delete fal2;
    exit(EXIT_FAILURE);
  }

  if ((fal2!=NULL) && (fal2->fail())) {
    cerr<<"Error: Cannot open input file '"<<align2_file<<"'\n";
    delete fal2;
    delete fal1;
    exit(EXIT_FAILURE);
  }

  time_t start_time, end_time;

  Alignment::debug=debug;
  Alignment::init_stats();

  start_time=time(NULL);

  cerr<<"Are alignment crossings allowed? \n";
  if (with_offsets) {
    cerr<<"    YES. ";
    if (non_contiguous) {
      cerr<<"Non-contiguos multiwords will be allowed (more the one offset may be required)\n";
      if (encode_offsets)
        cerr<<"    Offsets will be encoded in base "<<Alignment::code_base<<" to have a single integer\n";
      else
        cerr<<"    Offsets will not be encoded (comma-separated offsets)\n";
    } else
      cerr<<"Non-contiguos multiwords will *not* be allowed (only one offset is required)\n";
  } else {
    cerr<<"    NO.\n";
  }

  cerr<<"Switch? "<<switch_alg<<"\n";
  cerr<<"Minimum frequency for a biword to be used: "<<freqth<<"\n";
  cerr<<"Generation of the text to compress started at: "<<ctime(&start_time)<<"\n";

  if (rafa)
    cerr<<"Rafa version!\n";

  if (freqth==0.0) 
    process_input(fal1, fal2, with_offsets, non_contiguous, encode_offsets, rafa, freqth, switch_alg);
  else {
    cerr<<"First pass to get biword statistics\n";

    process_input(fal1, fal2, with_offsets, non_contiguous, encode_offsets, rafa, -1.0, switch_alg);

    //Alignment::print_biwords_vocabulary();

    cerr<<"Computing biword frequencies\n";
    int nbiwords_first=Alignment::compute_biwords_frequency();

    if (fal1!=NULL) delete fal1;
    if (fal2!=NULL) delete fal2;

    if (use_zlib) {
      fal1 = new gzifstream(align1_file.c_str());
      if (align2_file != "")
        fal2 = new gzifstream(align2_file.c_str());
    }  else {
      fal1 = new ifstream(align1_file.c_str());
      if (align2_file != "")
        fal2 = new ifstream(align2_file.c_str());
    }
    cerr<<"Second pass to get the final biword representation\n";
    process_input(fal1, fal2, with_offsets, non_contiguous, encode_offsets, rafa, freqth, switch_alg);
    int nbiwords_second=Alignment::compute_biwords_frequency();

    cerr<<"# biwords in the first pass:  "<<nbiwords_first<<"\n";
    cerr<<"# biwords in the second pass: "<<nbiwords_second<<"\n";
    cerr<<"% of reduction: "<<((double)(nbiwords_first-nbiwords_second))/((double)nbiwords_first)*100<<" %\n";
  }
  end_time=time(NULL);

  cerr<<"\nFinished at: "<<ctime(&end_time)<<"\n";
  cerr<<"Generation of the text to compress took "<<difftime(end_time, start_time)<<" seconds\n";

  Alignment::print_alig_stats();

  delete fal1;

  if (fal2!=NULL)
    delete fal2;
}
