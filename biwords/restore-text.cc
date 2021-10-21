#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <cstdlib>
#include <list>
#include <vector>

#include "configure.h"
#include "utils.h"
#include "alignment.h"

using namespace std;

bool debug;


void restore_right_text_rafa(ifstream &fin) {
  //This method assumes that offsets are not encoded and multiple offsets occcur
  //It also assumes that offsets refers to the right-hand side of the biwords (one-to-many)

  string str;
  string output[256]; //Maxium number of words per sentence
  unsigned n=0;

  for (unsigned i=0; i<256; i++)
    output[i]="";

  while (!fin.eof()) {
    getline(fin,str);

    if (str.length() == 0)
      continue;

    if (str == "Ɛ|Ɛ") {
      for(unsigned i=0; i<n; i++) {
        if (i!=0) cout<<" ";
        cout<<output[i];
      }
      if (n>0) cout<<"\n";

      n=0;
      for (unsigned i=0; i<256; i++)
        output[i]="";

      continue;
    }

    vector<string> fields = Utils::split_string(str, "|");

    if (fields.size()==0) {
      cerr<<"Error: no fields found in string '"<<str<<"'\n";
      exit(EXIT_FAILURE);
    }

    vector<int> offset;
    vector<string> target;

    //if (fields.size()>0) cerr<<"fields[0]="<<fields[0]<<endl;
    //if (fields.size()>1) cerr<<"fields[1]="<<fields[1]<<endl;
    //if (fields.size()>2) cerr<<"fields[2]="<<fields[2]<<endl;

    if (fields.size()>1) target = Utils::split_string(fields[1], " ");

    if (fields.size()>2) {
      vector<string> vaux=Utils::split_string(fields[2], ",");
      for(unsigned i=0; i<vaux.size(); i++)
        offset.push_back(atoi(vaux[i].c_str()));
    } 

    /*
    cerr<<"TARGET: ";
    for (unsigned i=0; i<target.size(); i++)
      cerr<<target[i]<<" ";
    cerr<<endl;
    
    cerr<<"OFFSETS: ";
    for (unsigned i=0; i<offset.size(); i++)
      cerr<<offset[i]<<" ";
    cerr<<endl;
    */

    while (offset.size()<target.size())
      offset.push_back(0);

    if (offset.size()!=target.size())
      cerr<<"Error: offset.size()= "<<offset.size()<<" target.size()= "<<target.size()<<"\n";

    unsigned k=n-1;
    for(unsigned i=0; i<target.size(); i++) {
      k+=offset[i]+1;
      output[k]=target[i];
      //k++;
    }

    //nextGap
    while ((n<256) && (output[n]!=""))
      n++;    
  }

  for(unsigned i=0; i<n; i++) {
    if (i!=0) cout<<" ";
    cout<<output[i];
  }
  if (n>0) cout<<"\n";
}

// -----------------------------------------------------------------------------------


void update_offsets_and_print (list<pair<string, int> > &buffer, list<string> &output) {
  bool printed = true;

  while (printed) {
    printed = false;
    list<pair<string, int> >::iterator it;

    for (it=buffer.begin(); it!=buffer.end(); it++) {
      if (it->second == 0) {
        cerr<<"Warning: 0 offset found in offsets buffer. This should *NEVER* happen\n";
        exit(EXIT_FAILURE);
      }
      it->second--;
    }

    for (it=buffer.begin(); it!=buffer.end(); it++) {
      if (it->second == 0) {
        output.push_back(it->first);
        it=buffer.erase(it);
        it--;
        printed = true;
      }
    }
  }
}

void restore_text(bool multiple_offsets, bool encoded_offsets, int side, ifstream &fin) {

  list<pair<string, int> > buffer; //first pos = word; second pos = offset
  string str;
  list<string> output;

  int max_offset=0;

  bool one_to_many_alg=true;

  while (!fin.eof()) {
    getline(fin,str);

    if (str.length() == 0)
      continue;

    if ((str == "Ɛ|Ɛ") || (str == "Ɛ|Ɛ|1")) {
      if (output.size()>0) {
        list<string>::iterator it;
        for(it=output.begin(); it!=output.end(); it++) {
          if (it!=output.begin())
            cout<<" ";
          cout<<*it;
        }
        cout<<"\n";
        buffer.clear();
        output.clear();
      }

      if (str == "Ɛ|Ɛ")
        one_to_many_alg=true; //offset will refer to side 2
      else 
        one_to_many_alg=false; //offset will refer to side 1

      continue;
    }

    vector<string> fields = Utils::split_string(str, "|");
    string src, tgt; 
    vector<int> offsets;

    if (fields.size()==0) {
      cerr<<"Error: no fields found in string '"<<str<<"'\n";
      exit(EXIT_FAILURE);
    }

    src = fields[0];
    if (fields.size()>1) tgt = fields[1];
    else tgt = "";

    if (fields.size()>2) {
      int offset_ant=0;

      vector<int> offsets_aux;

      if (multiple_offsets) {
        if (encoded_offsets)
          offsets_aux = Utils::decode(Alignment::code_base, atoi(fields[2].c_str()));
        else {
          vector<string> vaux=Utils:: split_string(fields[2], ",");
          for(unsigned i=0; i<vaux.size(); i++)
            offsets_aux.push_back(atoi(vaux[i].c_str()));
        }
      } else
        offsets_aux.push_back(atoi(fields[2].c_str()));

      for (unsigned i=0; i<offsets_aux.size(); i++) {
        int relative_offset = offsets_aux[i]; //atoi(offsets_str[i].c_str());
        int offset = relative_offset + offset_ant;
        offsets.push_back(offset);
        offset_ant = offset + 1;

        if (relative_offset>max_offset)
          max_offset=relative_offset;
      }
    } 

    if (debug) {
      cerr<<"side="<<side<<"; one_to_many_alg="<<one_to_many_alg<<"\n";
      cerr<<"string="<<str<<"\n";
      cerr<<"src="<<src<<"; tgt="<<tgt<<"; offsets={"<<Utils::vector2string(offsets)<<"}\n";
      cerr<<"Press ENTER ";
      getchar();
    }

    if ( ((side == 1) && (one_to_many_alg)) ||  ((side == 2) && (!one_to_many_alg)) ) {
      list<string>::iterator it=output.end();
      vector<string> words; 
      if (side == 1)
        words = Utils::split_string(src," ");
      else //side 2
        words = Utils::split_string(tgt," ");

      vector<string>::iterator it2;
      for (it2=words.begin(); it2!=words.end(); it2++)
        output.insert(it, *it2);

    } else {
      list<string>::iterator it=output.end();
      vector<string> words;
      if (side == 1)
        words = Utils::split_string(src," ");
      else
        words = Utils::split_string(tgt," ");

      if (offsets.size()==0) {
        vector<string>::iterator it2;
        for (it2=words.begin(); it2!=words.end(); it2++) {
          output.insert(it, *it2);
          update_offsets_and_print(buffer, output);
        }
      } else { 

        pair<string, int> auxpair;

        //Insert remaining offsets
        int prev_offset=0;
        if (offsets.size()>0)
          prev_offset = *(--offsets.end());
        while(offsets.size()<words.size()) {
          offsets.push_back(++prev_offset);
        }

        bool printed=false;
        for(unsigned i=0; i<words.size(); i++) {

          if (offsets[i]==0) {
            output.push_back(words[i]);
            printed=true;
          } else {
            auxpair.first=words[i];
            auxpair.second=offsets[i];
            buffer.push_back(auxpair);
          }
        }
        if(printed)
          update_offsets_and_print(buffer, output);
      }
    }

    if (debug) {
      cerr<<"-------------------------------\nbuffer\n-------------------------------\n";
      list<pair<string, int> >::iterator it;
      for(it=buffer.begin(); it!=buffer.end(); it++)
        cerr<<it->first<<" "<<it->second<<"; ";

      cerr<<"\n";

      cerr<<"-------------------------------\noutput\n-------------------------------\n";
      list<string>::iterator it2;
      for(it2=output.begin(); it2!=output.end(); it2++)
        cerr<<*it2<<" ";
      cerr<<"\n\n";
    }
  }

  if (output.size()>0) {
    list<string>::iterator it;
    for(it=output.begin(); it!=output.end(); it++) {
      if (it!=output.begin())
        cout<<" ";
      cout<<*it;
    }
    cout<<"\n";
    buffer.clear();
    output.clear();
  }

  //cerr<<"MAX OFFSET: "<<max_offset<<"\n";
}

void help(char *name) {
  cerr<<"USAGE:\n";
  cerr<<name<<" --input file [--moffsets [--encoded-offsets]] --side 1|2 [--debug]\n\n";
  cerr<<"ARGUMENTS: \n"
      <<"   --input|-i: Specify a file containing text in the format generated by gen-text-to-compress\n"
      <<"   --moffsets|-m: Specify the presence of multiple offsets in the input file\n"
      <<"   --encoded-offsets|-e: Multiple oofsets have been encoded in a single integer (base: "<<Alignment::code_base<<")\n"
      <<"   --side|-s: Specify the side of the bitext to restore\n"
      <<"   --help|-h: Show this help\n"
      <<"   --debug|-d: Show debug information\n"
      <<"   --version|-v: Show version information\n";
}

int main(int argc, char* argv[]) {
  int c;
  int option_index=0;

  string input_file="";
  int side=0;
  bool multiple_offsets=false; 
  bool encoded_offsets=false;

  bool rafa_version=false;

  debug=false;

  while (true) {
    static struct option long_options[] =
      {
	{"input",     required_argument,  0, 'i'},
	{"side",      required_argument,  0, 's'},
	{"moffsets",        no_argument,  0, 'm'},
	{"encoded-offsets", no_argument,  0, 'e'},
        {"rafa",            no_argument,  0, 'r'},
	{"help",            no_argument,  0, 'h'},
	{"version",         no_argument,  0, 'v'},
	{"debug",           no_argument,  0, 'd'},
	{0, 0, 0, 0}
      };

    c=getopt_long(argc, argv, "i:s:merhvd",long_options, &option_index);
    if (c==-1)
      break;
      
    switch (c) {
    case 'i':
      input_file=optarg;
      break;
    case 's':
      side=atoi(optarg);
      break;
    case 'm':
      multiple_offsets=true;
      break;
    case 'e':
      encoded_offsets=true;
      break;
    case 'd':
      debug=true;
      break;
    case 'r':
      rafa_version=true;
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

  if (input_file=="") {
    cerr<<"Error: No input file was given. You need to provide it with the --input option\n";
    help(argv[0]);
    exit(EXIT_FAILURE);
  }

  if ((side != 1) && (side != 2)) {
    cerr<<"Error: side to restore is not either 1 or 2.\n";
    help(argv[0]);
    exit(EXIT_FAILURE);
  }

  ifstream fin(input_file.c_str());

  if (fin.fail()) {
    cerr<<"Error: Cannot open input file '"<<input_file<<"'\n";
    exit(EXIT_FAILURE);
  }

  if (rafa_version) {
    cerr<<"Rafa version!!!!!\n";
    restore_right_text_rafa(fin);
  } else
    restore_text(multiple_offsets, encoded_offsets, side, fin);

  fin.close();
}
