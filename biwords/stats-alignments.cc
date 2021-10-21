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
  cerr<<name<<" --align alignments [-read num] [--gzip]\n\n";
  cerr<<"ARGUMENTS: \n"
      <<"   --align|-a: Specify a file containing the alignments\n"
      <<"   --gzip|-z: Tell the program that the alignment file is gziped\n"
      <<"   --read|-r: Tell the program how many alignments to read from the input file\n"
      <<"   --help|-h: Show this help\n"
      <<"   --version|-v: Show version information\n";
}

int main(int argc, char* argv[]) {
  int c;
  int option_index=0;

  string align_file="";

  bool use_zlib=false;

  int read=-1;

  //cerr<<"Command line: ";
  //for(int i=0; i<argc; i++)
  //  cerr<<argv[i]<<" ";
  //cerr<<"\n\n";
  //cerr<<"LOCALE: "<<setlocale(LC_ALL,"")<<"\n";


  while (true) {
    static struct option long_options[] =
      {
	{"align",  required_argument,  0, 'a'},
	{"read",  required_argument,  0, 'r'},
	{"gzip",         no_argument,  0, 'z'},
	{"help",         no_argument,  0, 'h'},
	{"version",      no_argument,  0, 'v'},
	{0, 0, 0, 0}
      };

    c=getopt_long(argc, argv, "a:r:zhv",long_options, &option_index);
    if (c==-1)
      break;
      
    switch (c) {
    case 'a':
      align_file=optarg;
      break;
    case 'r':
      read=atoi(optarg);
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

  if (align_file=="") {
    cerr<<"Error: No alignment file was given. You need to provide it with the --align option\n";
    help(argv[0]);
    exit(EXIT_FAILURE);
  }

  istream *fal;
  if (use_zlib) {
    fal = new gzifstream(align_file.c_str());
  }  else {
    fal = new ifstream(align_file.c_str());
  }

  if (fal->fail()) {
    cerr<<"Error: Cannot open input file '"<<align_file<<"'\n";
    delete fal;
    exit(EXIT_FAILURE);
  }

  time_t start_time, end_time;

  start_time=time(NULL);
  cerr<<"Collecting statistics, start time: "<<ctime(&start_time)<<"\n";

  string al;

  long nal=0;
  Alignment::init_stats();

  int one_to_many=0;

  while (!fal->eof()) {
    getline(*fal,al);

    if (al.length()>0) {
      nal++;
      Alignment alig(al);

      alig.collect_stats();

      if (alig.one_to_many())
        one_to_many++;

      if ((nal%50000)==0)
	cerr<<nal<<" alignments processed\n";
    }

    if ((read>0) && (nal>read))
      break;
  }

  Alignment::print_stats();

  if (one_to_many==nal)
    cerr<<"ONE TO MANY (1:N)\n";
  else 
    cerr<<"MANY TO ONE (N:1)\n";

    end_time=time(NULL);
  cerr<<"Finished at: "<<ctime(&end_time)<<"\n";
  cerr<<"It took "<<difftime(end_time, start_time)<<" seconds\n";

  delete fal;
}
