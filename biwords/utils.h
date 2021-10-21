#ifndef __UTILS_H_
#define __UTILS_H_

#include <string>
#include <vector>
#include <set>

using namespace std;

class Utils {
  public:
  
  static string trim(string str);

  static vector<string> split_string(const string& input, const string& delimiter);

  static string vector2string(const vector<string>& v);
  static string vector2string(const vector<int>& v);

  static string set2string(const set<int>& s);

  //Replace each ocurrence of the string 'olds' by the string 'news' in string 'source'
  static string substitute(const string& source, const string& olds, const string& news);

  static string itoa(int n);
  
  static string ftoa(double f);

  static string strtolower(const string& s);

  static bool case_insensitive_equal(const string& a, const string& b);

  //This function transforms all the values into a single integer
  //by changing the numbering base.
  static int code(int base, vector<int> values);

  static vector<int> decode(int base, int v);
};

#endif
