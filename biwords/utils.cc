#include "utils.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

//Delete white spaces from the end and the begining of the string
string 
Utils::trim(string str) { 
  string::iterator it;
  
  while ((str.length()>0)&&((*(it=str.begin()))==' ')) {
     str.erase(it);
  }
  
  while ((str.length()>0)&&((*(it=(str.end()-1)))==' ')) {
     str.erase(it);
  }

  return str;
}

vector<string>
Utils::split_string(const string& input, const string& delimiter) {
  int pos;
  int new_pos;
  vector<string> result;
  string s="";
  pos=0;

  while (pos<(int)input.size()) {
    new_pos=input.find(delimiter, pos);
    if(new_pos<0)
      new_pos=input.size();
    s=input.substr(pos, new_pos-pos);
    //if (s.length()==0) {
    //cerr<<"Warning in Utils::split_string: After splitting there is an empty string\n";
    //cerr<<"Skipping this empty string\n";
    //} else
    result.push_back(s);
    pos=new_pos+delimiter.size();
  }
  return result;
}

string 
Utils::vector2string(const vector<string>& v) {
  string s="";
  for(unsigned i=0; i<v.size(); i++) {
    if (i>0)
      s+=" ";
    s+=v[i];
  }
  return s;
}

string 
Utils::vector2string(const vector<int>& v) {
  string s="";
  for(unsigned i=0; i<v.size(); i++) {
    if (i>0)
      s+=" ";
    s+=Utils::itoa(v[i]);
  }
  return s;
}

string 
Utils::set2string(const set<int>& s) {
  string str="";
  set<int>::iterator it;

  for(it=s.begin(); it!=s.end(); it++) {
    if (str.length()>0)
      str += " ";
    str += itoa(*it);
  }

  return str;
}

string 
Utils::substitute(const string& source, const string& olds, const string& news) {
  string s=source;

  int p=s.find(olds,0);
  while (p!=(int)string::npos) {
    s.replace(p, olds.length(), news);
    p+=news.length();
    p=s.find(olds,p);
  }

  return s;
}

string
Utils::itoa(int n) {
  char str[32];
  sprintf(str, "%d",n);
  return str;
}

string
Utils::ftoa(double f) {
  char str[64];
  sprintf(str, "%g",f);
  return str;
}

bool 
Utils::case_insensitive_equal(const string& a, const string& b) {
  string alower="";
  string blower="";

  for(unsigned i=0; i<a.length(); i++) {
    alower+=tolower(a[i]);
  }

  for(unsigned i=0; i<b.length(); i++) {
    blower+=tolower(b[i]);
  }

  return (alower==blower);
}

string
Utils::strtolower(const string& s) {
  string l="";
  for(unsigned i=0; i<s.length(); i++)
    l+=tolower(s[i]);
  return l;
}

int 
Utils::code(int base, vector<int> values) {
  int v=0;

  //cerr<<"base "<<base<<" code ";
  for (unsigned i=0; i<values.size(); i++) {
    //cerr<<values[i]<<" ";
    v += values[i]*pow((double)base, (double)i);
  }
  //cerr<<"\n";
  return v;
}

vector<int>
Utils::decode(int base, int v) {
  vector<int> values;

  div_t div_result;

  while (v>=base) {
    div_result = div(v, base);
    values.push_back(div_result.rem);
    v=div_result.quot;
  }

  values.push_back(v);
  return values;
}
