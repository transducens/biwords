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
  cerr<<name<<" --input1|-i algin1 --input2|-j algin2 --output1|-o algout1 --output2|-p algout2 [--gzip]\n";
}

int main(int argc, char* argv[]) {
  int c;
  int option_index=0;

  string al1_file="";
  string al2_file="";
  string alout1_file="";
  string alout2_file="";

  int al1_one_to_many=0;
  int al1_many_to_one=0;

  int al2_one_to_many=0;
  int al2_many_to_one=0;

  bool use_zlib=false;

  while (true) {
    static struct option long_options[] =
      {
	{"input1",    required_argument,  0, 'i'},
	{"input2",    required_argument,  0, 'j'},
	{"output1",   required_argument,  0, 'o'},
	{"output2",   required_argument,  0, 'p'},
	{"gzip",          no_argument,  0, 'z'},
	{"help",          no_argument,  0, 'h'},
	{"version",       no_argument,  0, 'v'},
	{0, 0, 0, 0}
      };

    c=getopt_long(argc, argv, "i:j:o:p:zhv",long_options, &option_index);
    if (c==-1)
      break;
      
    switch (c) {
    case 'i':
      al1_file=optarg;
      break;
    case 'j':
      al2_file=optarg;
      break;
    case 'o':
      alout1_file=optarg;
      break;
    case 'p':
      alout2_file=optarg;
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

  if (al1_file=="") {
    cerr<<"Error: No input1 file was given. You need to provide it with the --input1 option\n";
    help(argv[0]);
    exit(EXIT_FAILURE);
  }

  if (al2_file=="") {
    cerr<<"Error: No input2 file was given. You need to provide it with the --input2 option\n";
    help(argv[0]);
    exit(EXIT_FAILURE);
  }

  if (alout1_file=="") {
    cerr<<"Error: No output1file was given. You need to provide it with the --output1 option\n";
    help(argv[0]);
    exit(EXIT_FAILURE);
  }

  if (alout2_file=="") {
    cerr<<"Error: No output2 file was given. You need to provide it with the --output2 option\n";
    help(argv[0]);
    exit(EXIT_FAILURE);
  }
  

  istream *fal1;
  if (use_zlib) {
    fal1 = new gzifstream(al1_file.c_str());
  }  else {
    fal1 = new ifstream(al1_file.c_str());
  }

  if (fal1->fail()) {
    cerr<<"Error: Cannot open input file '"<<al1_file<<"'\n";
    delete fal1;
    exit(EXIT_FAILURE);
  }

  istream *fal2;
  if (use_zlib) {
    fal2 = new gzifstream(al2_file.c_str());
  }  else {
    fal2 = new ifstream(al2_file.c_str());
  }

  if (fal2->fail()) {
    cerr<<"Error: Cannot open input file '"<<al2_file<<"'\n";
    delete fal1;
    delete fal2;
    exit(EXIT_FAILURE);
  }

  ostream *fout1;
  if(use_zlib) {
    fout1 = new gzofstream(alout1_file.c_str());
  } else {
    fout1 = new ofstream(alout1_file.c_str());
  }

  ostream *fout2;
  if(use_zlib) {
    fout2 = new gzofstream(alout2_file.c_str());
  } else {
    fout2 = new ofstream(alout2_file.c_str());
  }

  if (fout1->fail()) {
    cerr<<"Error: Cannot open output file '"<<alout1_file<<"'\n";
    delete fal1;
    delete fal2;
    delete fout1;
    delete fout2;
    exit(EXIT_FAILURE);
  }

  if (fout2->fail()) {
    cerr<<"Error: Cannot open output file '"<<alout2_file<<"'\n";
    delete fal1;
    delete fal2;
    delete fout1;
    delete fout2;
    exit(EXIT_FAILURE);
  }

  string al1;
  string al2;

  long nal=0;

  while ((!fal1->eof())&&(!fal2->eof())) {
    getline(*fal1,al1);
    getline(*fal2,al2);

    if ((al1.length()>0) && (al2.length()>0)) {
      nal++;
      Alignment alig1(al1);
      Alignment alig2(al2);


      if (alig1.are_the_same_alignment(alig2)) {
        if (alig1.one_to_many())
          al1_one_to_many++;
        else
          al1_many_to_one++;

        if (alig2.one_to_many())
          al2_one_to_many++;
        else
          al2_many_to_one++;

	(*fout1)<<alig1<<"\n";	
	(*fout2)<<alig2<<"\n";	
      } else {
	cerr<<"Warning: the alignment at line "<<nal<<" will be discarded\n";
        cerr<<"ALG1: "<<alig1<<"\n";
        cerr<<"ALG2: "<<alig2<<"\n";
      }
      if ((nal%10000)==0)
	cerr<<nal<<" alignments processed\n";
    }

    if ((fal1->eof()!=fal2->eof()) || (fal1->good()!=fal2->good())) {
      cerr<<"Error: The two input file do not contain the same number of alignments.\n";
      exit(EXIT_FAILURE);
    }
				    
  }

  if (al1_one_to_many > al1_many_to_one)
    cerr<<al1_file<<" is ONE TO MANY\n";
  else
    cerr<<al1_file<<" is MANY TO ONE\n";

  if (al2_one_to_many > al2_many_to_one)
    cerr<<al1_file<<" is ONE TO MANY\n";
  else
    cerr<<al1_file<<" is MANY TO ONE\n";

  delete fal1;
  delete fal2;
  delete fout1;
  delete fout2;
}
