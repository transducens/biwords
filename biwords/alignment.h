#ifndef __ALIGNMENT_H_
#define __ALIGNMENT_H_

using namespace std;

#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <set>

struct stats_struct { 
  int source_nwords;
  int target_nwords;
  map<int, int> alig_source; // alig_source[n] = num. of SL words aligned with n TL words
  map<int, int> alig_target; // alig_target[m] = num. of TL words aligned with m SL words

  int nalignments;
  int nalignments_discarded;
};

class Alignment {
private:
  double score;
  vector<string> source;
  vector<string> target;

  string alignment_str;

  static stats_struct stats;
  static map<string, map<string, int> > count_alignments;

  static map<string, map<string, int> > count_biwords;
  static int nbiwords;
  static map<string, map<string, double> > freq_biwords;

  //Matrix source.size x target.size with the alignments
  map<int, map<int, bool> > alignment;

  //Returns the set of indexes of the target words aligned with the
  //source word received as input (-1 if unaligned)
  set<int> target_aligned (int s);

  //Returns the set of indexes of the source word aligned with the
  //target word received as input (-1 if unaligned)
  set<int> source_aligned (int t);

  //True if the elements in SP are consecutive
  bool consecutive(const set<int>& SP);

  //Return a string with the words at the given indexes (SP) found in
  //words
  string build_str(const set<int>& SP, const vector<string>& words);

  void update_count_alignments();

  void update_count_biwords(const vector<string>& source_output, const vector<string>& target_output);

  //Process the alignment and updates s, t, source_output, target_output and offsets_output accordingly
  //The bahaviour is controlled by flag non_contiguous
  void get_biword_offsets(Alignment *alg, int &s, int &t, const set<int> &sa_set, const set<int> &ta_set, bool non_contiguous, 
                          vector<string> &source_output, vector<string> &target_output, vector<string> &offsets_output, 
                          bool encode_offsets);

  void remove_least_frequent_alignment(int frompos, const set<int>& topos, bool from_source);

  //Change source side by target side, and vice versa, updating the alignment matrix accordingly
  void switch_alignment();

  int move_forward(int s, int t);
public: 

  static bool debug;

  static int const code_base;
       
  Alignment();

  Alignment(string al, int nfields=4);

  Alignment(const Alignment& al);
    
  ~Alignment();


  unsigned unaligned_source_words();
  unsigned unaligned_target_words();

  inline unsigned source_words() {return source.size();};
  inline unsigned target_words() {return target.size();};

  static void init_stats();

  static void print_stats();
  static void print_alig_stats();

  static int compute_biwords_frequency();
  static void print_biwords_vocabulary();

  string to_string();

  bool one_to_many();

  //True if source size and the target size of both alignments (*this
  //and al2) are equal
  bool are_the_same_alignment(const Alignment& al2);

  //Generates biwords using offsets (alignment crossing is allowed)
  string biwords_with_offsets(Alignment* alg2=NULL);                 // ok - old implementation
  string biwords_with_offsets_non_contiguous(Alignment* alg2=NULL);  // ok - old implementation

  string biwords_with_offsets_new_version(bool switch_alg, double freqth, bool non_contiguous=false, bool encode_offsets=false, Alignment* alg2=NULL);  // ok - integrates previous two versions
  string biwords_with_offsets_new_version_rafa(bool non_contiguous=false, bool encode_offsets=false, Alignment* alg2=NULL);  // ok - Rafa's version (comentarios paper)

  string biwords_with_offsets_old(Alignment* alg2=NULL); // old implementation, uses 1:N and N:1 at the same time

  string biwords_with_offsets_alg1_1(); // old implementation
  string biwords_with_offsets_alg1_n_discontinuous(Alignment* alg2=NULL); //old implementation


  void collect_stats();

  //Generates biwords without using offsets (alignment crossing is *not* allowed)
  string biwords_without_offsets(Alignment* alg2=NULL); // ok

  string biwords_without_offsets_alg1_1(); //old implementation

  //(*this) becomes the intersected alignment
  bool intersection(Alignment& al2);

  //(*this) becomes the united alignment
  bool unionn(Alignment& al2);

  //(*this) becomes the refined intersected alignment
  bool refined_intersection(Alignment& al2);

  friend ostream& operator << (ostream& os, Alignment& al);
};

#endif
