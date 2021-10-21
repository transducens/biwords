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
  cerr<<name<<" --input|-i alignments-file [--gzip]\n";
}

int main(int argc, char* argv[]) {
  int c;
  int option_index=0;

  string al_file="";

  bool use_zlib=false;

  while (true) {
    static struct option long_options[] =
      {
	{"input",   required_argument,  0, 'i'},
	{"gzip",          no_argument,  0, 'z'},
	{"help",          no_argument,  0, 'h'},
	{"version",       no_argument,  0, 'v'},
	{0, 0, 0, 0}
      };

    c=getopt_long(argc, argv, "i:s:zhv",long_options, &option_index);
    if (c==-1)
      break;
      
    switch (c) {
    case 'i':
      al_file=optarg;
      break;
    case 'z':
      use_zlib=true;
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

  if (al_file=="") {
    cerr<<"Error: No input file was given. You need to provide it with the --input option\n";
    help(argv[0]);
    exit(EXIT_FAILURE);
  }

  istream *fal;
  if (use_zlib) {
    fal = new gzifstream(al_file.c_str());
  }  else {
    fal = new ifstream(al_file.c_str());
  }

  if (fal->fail()) {
    cerr<<"Error: Cannot open input file '"<<al_file<<"'\n";
    delete fal;
    exit(EXIT_FAILURE);
  }

  string al;

  long nal=0;

  ulong num_source_words=0;
  ulong num_source_words_unaligned=0;

  ulong num_target_words=0;
  ulong num_target_words_unaligned=0;

  ulong alg_one_to_many=0;

  while (!fal->eof()) {
    getline(*fal,al);

    if (al.length()>0) {
      nal++;
      Alignment alig(al);

      num_source_words+=alig.source_words();
      num_source_words_unaligned+=alig.unaligned_source_words();

      num_target_words+=alig.target_words();
      num_target_words_unaligned+=alig.unaligned_target_words();

      alg_one_to_many+=alig.one_to_many();

      if ((nal%10000)==0)
	cerr<<nal<<" alignments processed\n";
    }				    
  }
  if (alg_one_to_many > (nal-alg_one_to_many)) 
    cout<<"ONE TO MANY (1:N)"<<endl<<endl;
  else
    cout<<"MANY TO ONE (N:1)"<<endl<<endl;

  cout<<"unaligned source words: "<<((double) num_source_words_unaligned)/((double) num_source_words)*100.0<<" %"<<endl;
  cout<<"unaligned target words: "<<((double) num_target_words_unaligned)/((double) num_target_words)*100.0<<" %"<<endl;
  cout<<endl;
  cout<<"unaligned words: "<<((double) (num_source_words_unaligned+num_target_words_unaligned))/((double) (num_source_words+num_target_words))*100.0<<" %"<<endl;

  delete fal;
}
