#include "alignment.h"
#include "utils.h"

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>

stats_struct Alignment::stats;
map<string, map<string, int> > Alignment::count_alignments;
bool Alignment::debug;
int const Alignment::code_base=30;

map<string, map<string, int> > Alignment::count_biwords;
map<string, map<string, double> > Alignment::freq_biwords;
int Alignment::nbiwords=0;

void 
Alignment::init_stats() {
  stats.source_nwords=0;
  stats.target_nwords=0;
  stats.nalignments=0;
  stats.nalignments_discarded=0;
}

Alignment::Alignment() {
}

Alignment::Alignment(string al, int nfields) {
  vector<string> v;
  vector<string> alig;

  v=Utils::split_string(al, "|");

  if (v.size()==(unsigned)(nfields-1)) { //The alignment information is missing
    cerr<<"Warning: Following alignment has all words unaligned: '"<<al<<"'\n";
    v.push_back("");
  }

  if (v.size()!=(unsigned)nfields) {
    cerr<<"Error in Alignment::Alignment when reading alignment from string '"<<al<<"'\n";
    cerr<<"Unespected number of fields separated by '|'\n";
    cerr<<"Num. fields: "<<v.size()<<"\n";
    for (unsigned i=0; i<v.size(); i++)
      cerr<<"v["<<i<<"]=+"<<v[i]<<"+\n";

    exit(EXIT_FAILURE); 
  }

  score=atof(Utils::trim(v[0]).c_str());
  source=Utils::split_string(Utils::trim(v[1]), " ");
  target=Utils::split_string(Utils::trim(v[2]), " ");

  alignment_str=Utils::trim(v[3]);

  alig=Utils::split_string(alignment_str, " ");

  for(unsigned i=0; i<alig.size(); i++) {
    vector<string> an_alig;

    an_alig=Utils::split_string(alig[i], ":");
    if (an_alig.size()!=2) {
      cerr<<"Error in Alignment::Alignment when reading alignment from string '"<<al<<"'\n";
      cerr<<"Unespected number of alignment values separated by ':'\n";
      exit(EXIT_FAILURE);
    }
    alignment[atoi(an_alig[0].c_str())][atoi(an_alig[1].c_str())]=true;
    stats.nalignments++;
  }
}

Alignment::Alignment(const Alignment& al) {
  source=al.source;
  target=al.target;
  score=al.score;
  alignment=al.alignment;
  alignment_str=al.alignment_str;
}
    
Alignment::~Alignment() {
}

string
Alignment::to_string() {
  string s;
  s=Utils::ftoa(score)+" |";

  for(unsigned i=0; i<source.size(); i++) 
    s+=" "+source[i];
  s+=" |";

  for(unsigned i=0; i<target.size(); i++) 
    s+=" "+target[i];
  s+=" |";

  for (unsigned i=0; i<source.size(); i++) {
    for(unsigned j=0; j<target.size(); j++) {
      if (alignment[i][j])
	s+=" "+Utils::itoa(i)+":"+Utils::itoa(j);
    }
  }

  return s;
}

string 
Alignment::biwords_with_offsets_alg1_1() {

  // This function assumes that the alignments are 1:1
  // This assumption is not checked

  int s=0, t=0;
  string str="";

  while (s<(int)source.size()) {
    set<int> ta_set = target_aligned(s);
    set<int> sa_set = source_aligned(t);

    int ta, sa;

    if (ta_set.size()==0) ta=-1;
    else ta=*(ta_set.begin());

    if (sa_set.size()==0) sa=-1;
    else sa=*(sa_set.begin());

    //cerr<<"      source word = "<<source[s]<<"\n";
    //cerr<<"      target word = "<<target[t]<<"\n";
    //cerr<<"      s="<<s<<"; t="<<t<<"; sa="<<sa<<"; ta="<<ta<<"\n";
    //cerr<<"\n";

    if (ta == -1 ) {
      str += source[s] + "|" + "\n";
      s++;
    } else if ((ta > t) && (sa == -1)) {
      str += "|" + target[t] + "\n";
      t++;
    } else if ((ta > t) && (sa < s)) {
      t++;
    } else if ((ta == t) && (sa == s)) {
      str += source[s] + "|" + target[t] + "\n";
      s++; t++;
    } else if ((ta > t) && (sa != -1)) {
      int offset = ta - t;
      str += source[s] + "|" + target[ta] + "|" + Utils::itoa(offset) + "\n";
      s++;
    } else {
      cerr<<"Warning at print_alignment_offsets: This message should *NEVER* appear.\n";
      cerr<<"Write an e-mail to fsanchez@dlsi.ua.es with the following:\n";
      cerr<<"SL: "<<Utils::vector2string(source)<<"\n";
      cerr<<"TL: "<<Utils::vector2string(target)<<"\n";
      cerr<<"AL: "<<alignment_str<<"\n";
      exit(EXIT_FAILURE);
    }
  }

  for (; t<(int)target.size(); t++) {
    set<int> sa_set = source_aligned(t);
    if (sa_set.size()==0) {
      str += "|" + target[t] + "\n";
    }
  }

  return str;
}

string 
Alignment::biwords_without_offsets_alg1_1() {

  // This function assumes that the alignments are 1:1
  // This assumption is not checked

  string str="";

  int s=0, t=0;

  while (s<(int)source.size()) {
    set<int> ta_set=target_aligned(s);

    int ta;
    if (ta_set.size()==0) ta=-1;
    else ta=*(ta_set.begin());

    if (ta == -1) {
      str += source[s] + "|" + "\n";
      s++;
    } else if (ta == t) {
      str += source[s] + "|" + target[t] + "\n";
      s++; t++;
    } else if (ta<t) {
      //Ignore this alignment. There is a crossing in alignments
      //target[ta] has been already printed out
      str += source[s] + "|" + "\n";
      s++;
    } else {
      for(; t<ta; t++) {
        str += "|" + target[t] + "\n";
        //t++;
      }
    }
  }

  for(; t<(int)target.size(); t++) {
    str += "|" + target[t] + "\n";
  }

  return str;
}



string 
Alignment::biwords_with_offsets_alg1_n_discontinuous(Alignment* alg2) {

  Alignment* alg = this; //Alignment to use

  // This function assumes 1:N (alg1) and N:1 (alg2) alignments

  if (alg2!=NULL) {
    if (!are_the_same_alignment(*alg2)) {
      cerr<<"Warning, following alignments are not equal, and therefore discarded:\n";
      cerr<<to_string()<<"\n";
      cerr<<alg2->to_string()<<"\n";
      return "";
    }

    //We use the alignment with the highest score
    if (alg2->score>score)
      alg=alg2;
  }

  int s=0, t=0;
  string str="";

  //Cuando se ecuentre un alineamientos 1:N en el que las palabras no sean
  //consecutivas, end_alg será igual a la posición de fin en la otra lengua
  int end_alg=-1; bool alg_1_n=true; //true si el alineamiento 1:N, false en otro caso
  int pivot=-1;  //Será igual a la posicion de la que parten N alineamientos no consecutivos
                 //si alg_1_n==true => se refiere a source, sino a target
  set<int> alg_no_cont; //Conjunto de posiciones alineadas con pivot 

  while (s<(int)alg->source.size()) {
    set<int> ta_set = alg->target_aligned(s);
    set<int> sa_set = alg->source_aligned(t);

    bool consecutive_ta = alg->consecutive(ta_set);
    bool consecutive_sa = alg->consecutive(sa_set);

    int ta_size = ta_set.size();
    int sa_size = sa_set.size();

    int ta_first;
    int ta_last;
    if (ta_size > 0) {
      ta_first=*ta_set.begin();
      ta_last=*(--ta_set.end());
    } else {
      ta_first=-1;
      ta_last=-1;
    }

    int sa_first;
    int sa_last;
    if (sa_size > 0) {
      sa_first=*sa_set.begin();
      sa_last=*(--sa_set.end());
    } else {
      sa_first=-1;
      sa_last=-1;
    }

    if (alg_1_n) {
      alg_no_cont.erase(t);
    } else {
      alg_no_cont.erase(s);
    }

    if (debug) {
      cerr<<"s="<<s<<"; t="<<t<<"\n";
      cerr<<"ta_set={"<<Utils::set2string(ta_set)<<"}; consecutive="<<consecutive_ta<<"\n";
      cerr<<"sa_set={"<<Utils::set2string(sa_set)<<"}; consecutive="<<consecutive_sa<<"\n";
      cerr<<"ta_first="<<ta_first<<"; ta_last="<<ta_last<<"\n";
      cerr<<"sa_first="<<sa_first<<"; sa_last="<<sa_last<<"\n";
      cerr<<"end_alg="<<end_alg<<"; alg_1_n="<<alg_1_n<<"; pivot="<<pivot<<"\n"; 
      cerr<<"alg_no_cont={"<<Utils::set2string(alg_no_cont)<<"}\n";
      //cerr<<"Press ENTER "<<flush;
      //getchar();
    }

    if ((!alg_1_n) && (end_alg==s)) {
      if (debug) cerr<<"6\n";
      s++;
      end_alg=-1;
      pivot=-1;
      alg_1_n=true;
      alg_no_cont.clear();
    } else if ((alg_1_n) && (end_alg==t)) {
      if (debug) cerr<<"7\n";
      t++;
      end_alg=-1;
      pivot=-1;
      alg_1_n=true;
      alg_no_cont.clear();
    } else if ((!alg_1_n) && (pivot==t)) {
      //Este alineamiento ya se imprimió
      if (debug) cerr<<"8\n";
      t++;
    } else if ((!alg_1_n) && (ta_size == 1) && (ta_first == pivot)) {
      //Este alineamiento ya se imprimió
      if (debug) cerr<<"9\n";
      s++;
    } else if (ta_size == 0) { //s no alineado
      if (debug) cerr<<"1\n";
      str += alg->source[s] + "|";

      //if (end_alg != -1)
      if ((end_alg != -1) && (!alg_1_n))
        str += "|-" + Utils::itoa(alg_no_cont.size());

      str +=  "\n";
      s++;
    } else if ((ta_size > 0) && (ta_first > t) && (sa_size == 0)) { //t no alineado
      if (debug) cerr<<"2\n";

      if (end_alg != -1) {
        if (debug) cerr<<"2.1\n";
        str += "|" + alg->target[t] + "|-" + Utils::itoa(alg_no_cont.size()) + "\n";
      } else {
        if (debug) cerr<<"2.2\n";
        str += "|" + alg->target[t] + "\n";
      }
      t++;
    } else if ((ta_size > 0) && (sa_size > 0) && (ta_first > t) && (sa_first < s)) {
      if (debug) cerr<<"3\n";
      t++; //Este alineamiento ya se imprimió
    } else if ((ta_size > 0) && (sa_size > 0) && (ta_first == t) && (sa_first == s)) {
      if (debug) cerr<<"4\n";
      if ((consecutive_ta) && (consecutive_sa)) {
        if (debug) cerr<<"4.1\n";

        str += alg->build_str(sa_set,alg->source) + "|" + alg->build_str(ta_set, alg->target);
        if (end_alg != -1) {
          if (debug) cerr<<"4.1.1\n";
          str += "|-" +Utils::itoa(alg_no_cont.size()) + "\n";
        } else {
          if (debug) cerr<<"4.1.2\n";
          str += "\n";
        }

        s+=sa_size;
        t+=ta_size;
      } else {
        if (debug) { 
          cerr<<"4.2\n";
          cerr<<"No consecutive !!!!!\n";
        }
        if (end_alg == -1) { //OK.
          if (debug) cerr<<"4.2.1\n";
          if (ta_size > 1) {
            end_alg = ta_last;
            alg_1_n = true;
            pivot=s;
            alg_no_cont = ta_set;
            alg_no_cont.erase(t);
          } else if (sa_size > 1) {
            end_alg = sa_last;
            alg_1_n = false;
            pivot = t;
            alg_no_cont = sa_set;
            alg_no_cont.erase(s);
          }

          str += alg->build_str(sa_set,alg->source) + "|" + alg->build_str(ta_set, alg->target) + "\n";
          s++;
          t++;
        } else { //No Ok. Hay que eliminar el alineamiento más improbable
          if (debug) {
            cerr<<"4.2.2\n";
            cerr<<"end_alg != -1\n";
          }
          if (ta_size == 1)
            alg->remove_least_frequent_alignment(t, sa_set, false);
          else
            alg->remove_least_frequent_alignment(s, ta_set, true);
        }
      }
    } else if ((ta_size > 0) && (sa_size > 0) && (ta_first > t)) {
      if (debug) cerr<<"5\n";

      if (ta_size == 1) {
        set<int> new_sa_set = alg->source_aligned(ta_first);
        if (consecutive(new_sa_set)) {
          if (debug) cerr<<"5.1\n";
          int offset = ta_first - t - alg_no_cont.size();
          str += alg->build_str(new_sa_set,alg->source) + "|" + alg->build_str(ta_set, alg->target) + "|" + Utils::itoa(offset) + "\n";
          s+=new_sa_set.size();
        } else {
          if (debug) {
            cerr<<"5.2\n";
            cerr<<"No consecutive !!!!!\n";
          }
          if (end_alg == -1) { //OK.
            end_alg = *(--new_sa_set.end());
            alg_1_n = false;
            pivot = ta_first;
            int offset = ta_first - t;
            alg_no_cont = new_sa_set;
            alg_no_cont.erase(s);

            str += alg->build_str(new_sa_set,alg->source) + "|" + alg->target[ta_first] + "|" + Utils::itoa(offset) + "\n";
            s++;
          } else {//No Ok. Hay que eliminar el alineamiento más improbable
            if (debug) cerr<<"end_alg != -1\n";
            alg->remove_least_frequent_alignment(ta_first, new_sa_set, false);
          }
        }
      } else { // ta_size > 1
        if (consecutive_ta) {
          if (debug) cerr<<"5.3\n";
          int offset = ta_first - t;
          str += alg->source[s] + "|" + alg->build_str(ta_set, alg->target) + "|" + Utils::itoa(offset) + "\n";
          s++;
        } else {
          if (debug) {
            cerr<<"5.4\n";
            cerr<<"No consecutive !!!!!\n";
          }
          alg->remove_least_frequent_alignment(s, ta_set, true);
        }
      }
    } else {
      cerr<<"Warning at Alignment::biwords_with_offsets: This message should *NEVER* appear.\n";
      cerr<<"Write an e-mail to fsanchez@dlsi.ua.es with the following:\n";
      cerr<<"SL: "<<Utils::vector2string(alg->source)<<"\n";
      cerr<<"TL: "<<Utils::vector2string(alg->target)<<"\n";
      cerr<<"AL: "<<alg->alignment_str<<"\n";
      exit(EXIT_FAILURE);
    }

    if (debug) cerr<<"output\n-----------\n"<<str<<"\n";
  }

  if (debug) cerr<<"Fin while\n";

  for (; t<(int)alg->target.size(); t++) {
    set<int> sa_set = alg->source_aligned(t);
    if (sa_set.size() == 0) {
      str += "|" + alg->target[t] + "\n";
    }
  }

  alg->update_count_alignments();

  return str;
}

string 
Alignment::biwords_with_offsets_old(Alignment* alg2) {

  Alignment* alg = this; //Alignment to use

  // This function assumes 1:N (alg1) and N:1 (alg2) alignments
  if (alg2!=NULL) {
    if (!are_the_same_alignment(*alg2)) {
      cerr<<"Warning, following alignments are not equal, and therefore discarded:\n";
      cerr<<to_string()<<"\n";
      cerr<<alg2->to_string()<<"\n";
      return "";
    }

    //We use the alignment with the highest score
    //Warning: This comparison could be unfair, scores come from different
    //alignment models
    if (alg2->score>score)
      alg=alg2;
  }

  if (alg->one_to_many()) {
    cerr<<"ONE TO MANY (1:N)\n";
  } else {
    cerr<<"MANY TO ONE (N:1)\n";
  }


  int s=0, t=0;
  string str="";

  while (s<(int)alg->source.size()) {
    set<int> ta_set = alg->target_aligned(s);
    set<int> sa_set = alg->source_aligned(t);

    bool consecutive_ta = alg->consecutive(ta_set);
    bool consecutive_sa = alg->consecutive(sa_set);

    int ta_size = ta_set.size();
    int sa_size = sa_set.size();

    int ta_first;
    int ta_last;
    if (ta_size > 0) {
      ta_first=*ta_set.begin();
      ta_last=*(--ta_set.end());
    } else {
      ta_first=-1;
      ta_last=-1;
    }

    int sa_first;
    int sa_last;
    if (sa_size > 0) {
      sa_first=*sa_set.begin();
      sa_last=*(--sa_set.end());
    } else {
      sa_first=-1;
      sa_last=-1;
    }

    if (debug) {
      cerr<<"s="<<s<<"; t="<<t<<"\n";
      cerr<<"ta_set={"<<Utils::set2string(ta_set)<<"}; consecutive="<<consecutive_ta<<"\n";
      cerr<<"sa_set={"<<Utils::set2string(sa_set)<<"}; consecutive="<<consecutive_sa<<"\n";
      cerr<<"ta_first="<<ta_first<<"; ta_last="<<ta_last<<"\n";
      cerr<<"sa_first="<<sa_first<<"; sa_last="<<sa_last<<"\n";
      //cerr<<"Press ENTER "<<flush;
      //getchar();
    }

    if (ta_size == 0) { //s is not aligned
      if (debug) cerr<<"1\n";
      str += alg->source[s] + "|\n";
      s++;
    } else if ((ta_size > 0) && (ta_first > t) && (sa_size == 0)) { //t is not aligned
      if (debug) cerr<<"2\n";
      str += "|" + alg->target[t] + "\n";
      t++;
    } else if ((ta_size > 0) && (sa_size > 0) && (ta_first > t) && (sa_first < s)) {
      if (debug) cerr<<"3\n";
      t++; // This alignment has been already printed
    } else if ((ta_size > 0) && (sa_size > 0) && (ta_first == t) && (sa_first == s)) {
      if (debug) cerr<<"4\n";
      if ((consecutive_ta) && (consecutive_sa)) {
        if (debug) cerr<<"4.1\n";
        str += alg->build_str(sa_set,alg->source) + "|" + alg->build_str(ta_set, alg->target) + "\n";
        s+=sa_size;
        t+=ta_size;
      } else {
        if (debug) cerr<<"4.2 - No consecutive !!!!!\n";
        if (ta_size == 1)
          alg->remove_least_frequent_alignment(t, sa_set, false);
        else
          alg->remove_least_frequent_alignment(s, ta_set, true);        
      }
    } else if ((ta_size > 0) && (sa_size > 0) && (ta_first > t)) {
      if (debug) cerr<<"5\n";

      if (ta_size == 1) {
        set<int> new_sa_set = alg->source_aligned(ta_first);
        if (alg->consecutive(new_sa_set)) {
          if (debug) cerr<<"5.1\n";
          int offset = ta_first - t;
          str += alg->build_str(new_sa_set,alg->source) + "|" + alg->build_str(ta_set, alg->target) + "|" + Utils::itoa(offset) + "\n";
          s+=new_sa_set.size();
        } else {
          if (debug) cerr<<"5.2 - No consecutive !!!!!\n";
          alg->remove_least_frequent_alignment(ta_first, new_sa_set, false);
        }
      } else { // ta_size > 1
        if (consecutive_ta) {
          if (debug) cerr<<"5.3\n";
          int offset = ta_first - t;
          str += alg->source[s] + "|" + alg->build_str(ta_set, alg->target) + "|" + Utils::itoa(offset) + "\n";
          s++;
        } else {
          if (debug) cerr<<"5.4 -No consecutive !!!!!\n";
          alg->remove_least_frequent_alignment(s, ta_set, true);
        }
      }
    } else {
      cerr<<"Warning at Alignment::biwords_with_offsets: This message should *NEVER* appear.\n";
      cerr<<"Write an e-mail to fsanchez@dlsi.ua.es with the following:\n";
      cerr<<"SL: "<<Utils::vector2string(alg->source)<<"\n";
      cerr<<"TL: "<<Utils::vector2string(alg->target)<<"\n";
      cerr<<"AL: "<<alg->alignment_str<<"\n";
      exit(EXIT_FAILURE);
    }

    if (debug) cerr<<"output\n-----------\n"<<str<<"\n";
  }

  if (debug) cerr<<"Fin while\n";

  for (; t<(int)alg->target.size(); t++) {
    set<int> sa_set = alg->source_aligned(t);
    if (sa_set.size() == 0) {
      str += "|" + alg->target[t] + "\n";
    }
  }

  alg->update_count_alignments();

  str = "Ɛ|Ɛ\n" + str;
  return str;
}


string 
Alignment::biwords_with_offsets(Alignment* alg2) {

  cerr<<"Running Alignment::biwords_with_offsets\n\n";

  Alignment* alg = this; //Alignment to use

  // This function assumes 1:N (alg1) and N:1 (alg2) alignments
  if (alg2!=NULL) {
    if (!are_the_same_alignment(*alg2)) {
      cerr<<"Warning, following alignments are not equal, and therefore discarded:\n";
      cerr<<to_string()<<"\n";
      cerr<<alg2->to_string()<<"\n";
      return "";
    }

    //We use the alignment with the highest score
    //Warning: This comparison could be unfair, scores come from different
    //alignment models
    if (alg2->score>score)
      alg=alg2;
  }

  int s=0, t=0;

  vector<string> source_output, target_output, offsets_output;

  //First we need to know if we are dealing with one-to-many (1:N) or
  //many-to-one (N:1) alignments
  bool one_to_many_alg = alg->one_to_many();

  if (!one_to_many_alg) {
    alg->switch_alignment(); 
    //cerr<<"MANY TO ONE (N:1)\n";
  } else {
    //cerr<<"ONE TO MANY (1:N)\n";
  }

  while (s<(int)alg->source.size()) {
    set<int> ta_set = alg->target_aligned(s);
    set<int> sa_set = alg->source_aligned(t);

    bool consecutive_ta = alg->consecutive(ta_set);
    bool consecutive_sa = alg->consecutive(sa_set);

    int ta_size = ta_set.size();
    int sa_size = sa_set.size();

    int ta_first;
    int ta_last;
    if (ta_size > 0) {
      ta_first=*ta_set.begin();
      ta_last=*(--ta_set.end());
    } else {
      ta_first=-1;
      ta_last=-1;
    }

    int sa_first;
    int sa_last;
    if (sa_size > 0) {
      sa_first=*sa_set.begin();
      sa_last=*(--sa_set.end());
    } else {
      sa_first=-1;
      sa_last=-1;
    }

    if (sa_size > 1) {
      cerr<<"Warning: This is not a 1:N alignment\n";
      cerr<<to_string()<<"\n";
      exit(EXIT_FAILURE);
    }

    if (debug) {
      cerr<<"s="<<s<<"; t="<<t<<"\n";
      cerr<<"ta_set={"<<Utils::set2string(ta_set)<<"}; consecutive="<<consecutive_ta<<"\n";
      cerr<<"sa_set={"<<Utils::set2string(sa_set)<<"}; consecutive="<<consecutive_sa<<"\n";
      cerr<<"ta_first="<<ta_first<<"; ta_last="<<ta_last<<"\n";
      cerr<<"sa_first="<<sa_first<<"; sa_last="<<sa_last<<"\n";
      //cerr<<"Press ENTER "<<flush;
      //getchar();
    }

    if (ta_size == 0) { //s is not aligned
      if (debug) cerr<<"1\n";
      //str += alg->source[s] + "|\n";
      source_output.push_back(alg->source[s]);
      target_output.push_back("");
      offsets_output.push_back("");

      s++;
    } else if ((ta_first > t) && (sa_size == 0)) { //t is not aligned
      if (debug) cerr<<"2\n";
      //str += "|" + alg->target[t] + "\n"
      source_output.push_back("");
      target_output.push_back(alg->target[t]);
      offsets_output.push_back("");

      t++;
    } else if ((sa_size > 0) && (ta_first > t) && (sa_first < s)) {
      if (debug) cerr<<"3\n";
      t++; // This alignment has been already printed
    } else if ((sa_size > 0) && (ta_first == t) && (sa_first == s)) {
      if (debug) cerr<<"4\n";
      if (consecutive_ta) {
        if (debug) cerr<<"4.1\n";
        //str += alg->build_str(sa_set,alg->source) + "|" + alg->build_str(ta_set, alg->target) + "\n";
        source_output.push_back(alg->source[s]); //alg->build_str(sa_set,alg->source));
        target_output.push_back(alg->build_str(ta_set, alg->target));
        offsets_output.push_back("");

        s++;
        t+=ta_size;
      } else {
        if (debug) cerr<<"4.2 - No consecutive !!!!!\n";
        alg->remove_least_frequent_alignment(s, ta_set, true);        
      }
    } else if ((sa_size > 0) && (ta_first > t)) {
      if (debug) cerr<<"5\n";
      //if (ta_size == 1) {
      //  set<int> new_sa_set = alg->source_aligned(ta_first);
      //  if (alg->consecutive(new_sa_set)) {
      //    if (debug) cerr<<"5.1\n";
      //    int offset = ta_first - t;
      //    str += alg->build_str(new_sa_set,alg->source) + "|" + alg->build_str(ta_set, alg->target) + "|" + Utils::itoa(offset) + "\n";
      //    s+=new_sa_set.size();
      //  } else {
      //    if (debug) cerr<<"5.2 - No consecutive !!!!!\n";
      //    alg->remove_least_frequent_alignment(ta_first, new_sa_set, false);
      //  }
      //} else { // ta_size > 1

      if (consecutive_ta) {
        if (debug) cerr<<"5.3\n";
        int offset = ta_first - t;
        //str += alg->source[s] + "|" + alg->build_str(ta_set, alg->target) + "|" + Utils::itoa(offset) + "\n";
        source_output.push_back(alg->source[s]);
        target_output.push_back(alg->build_str(ta_set, alg->target));
        offsets_output.push_back(Utils::itoa(offset));

        s++;
      } else {
        if (debug) cerr<<"5.4 -No consecutive !!!!!\n";
        alg->remove_least_frequent_alignment(s, ta_set, true);
      }
      //}
    } else {
      cerr<<"Warning at Alignment::biwords_with_offsets: This message should *NEVER* appear.\n";
      cerr<<"Write an e-mail to fsanchez@dlsi.ua.es with the following:\n";
      cerr<<"SL: "<<Utils::vector2string(alg->source)<<"\n";
      cerr<<"TL: "<<Utils::vector2string(alg->target)<<"\n";
      cerr<<"AL: "<<alg->alignment_str<<"\n";
      exit(EXIT_FAILURE);
    }

  }

  if (debug) cerr<<"Fin while\n";

  for (; t<(int)alg->target.size(); t++) {
    set<int> sa_set = alg->source_aligned(t);
    if (sa_set.size() == 0) {
      //str += "|" + alg->target[t] + "\n";
      source_output.push_back("");
      target_output.push_back(alg->target[t]);
      offsets_output.push_back("");
    }
  }

  alg->update_count_alignments();

  string str="";
  if (!one_to_many_alg) { //switch the output to reduce the number of different symbols
    vector<string> aux;

    aux=source_output;
    source_output=target_output;
    target_output=aux;
    str+="Ɛ|Ɛ|1\n";
  } else {
    str+="Ɛ|Ɛ\n";
  }

  for(unsigned i=0; i<source_output.size(); i++) {
    str += source_output[i] + "|" + target_output[i];
    if (offsets_output[i].length()>0)
      str += "|" + offsets_output[i];
    str+="\n";
  }

  return str;
}

string 
Alignment::biwords_with_offsets_non_contiguous(Alignment* alg2) {

  cerr<<"Running Alignment::biwords_with_offsets_non_contiguous\n\n";

  Alignment* alg = this; //Alignment to use

  // This function assumes 1:N (alg1) and N:1 (alg2) alignments
  if (alg2!=NULL) {
    if (!are_the_same_alignment(*alg2)) {
      cerr<<"Warning, following alignments are not equal, and therefore discarded:\n";
      cerr<<to_string()<<"\n";
      cerr<<alg2->to_string()<<"\n";
      return "";
    }

    //We use the alignment with the highest score 
    //Warning: This comparison could bu unfair, scores come from different
    //alignment models
    if (alg2->score>score)
      alg=alg2;
  }

  int s=0, t=0;

  vector<string> source_output, target_output, offsets_output;

  //First we need to know if we are dealing with one-to-many (1:N) or
  //many-to-one (N:1) alignments
  bool one_to_many_alg = alg->one_to_many();

  if (!one_to_many_alg) 
    alg->switch_alignment(); 

  while (s<(int)alg->source.size()) {
    set<int> ta_set = alg->target_aligned(s);
    set<int> sa_set = alg->source_aligned(t);

    int ta_size = ta_set.size();
    int sa_size = sa_set.size();

    int ta_first;
    int ta_last;
    if (ta_size > 0) {
      ta_first=*ta_set.begin();
      ta_last=*(--ta_set.end());
    } else {
      ta_first=-1;
      ta_last=-1;
    }

    int sa_first;
    int sa_last;
    if (sa_size > 0) {
      sa_first=*sa_set.begin();
      sa_last=*(--sa_set.end());
    } else {
      sa_first=-1;
      sa_last=-1;
    }

    if (sa_size > 1) {
      cerr<<"Warning: This is not a 1:N alignment\n";
      cerr<<to_string()<<"\n";
      exit(EXIT_FAILURE);
    }

    if (debug) {
      cerr<<"s="<<s<<"; t="<<t<<"\n";
      cerr<<"ta_set={"<<Utils::set2string(ta_set)<<"}\n";//; consecutive="<<consecutive_ta<<"\n";
      cerr<<"sa_set={"<<Utils::set2string(sa_set)<<"}\n";//; consecutive="<<consecutive_sa<<"\n";
      cerr<<"ta_first="<<ta_first<<"; ta_last="<<ta_last<<"\n";
      cerr<<"sa_first="<<sa_first<<"; sa_last="<<sa_last<<"\n";
      //cerr<<"Press ENTER "<<flush;
      //getchar();
    }

    if (ta_size == 0 ) { //s no alineada
      if (debug) cerr<<"1\n";
      //str += alg->source[s] + "|" + "\n";
      source_output.push_back(alg->source[s]);
      target_output.push_back("");
      offsets_output.push_back("");

      s++;
    } else if ((ta_first > t) && (sa_size == 0)) { //t no alineada
      if (debug) cerr<<"2\n";
      //str += "|" + alg->target[t] + "\n";
      source_output.push_back("");
      target_output.push_back(alg->target[t]);
      offsets_output.push_back("");

      t++;
    } else if ((ta_first > t) && (sa_last < s)) { //este alineamiento ya se imprimió
      if (debug) cerr<<"3\n";
      t++;
    } else if ((ta_first >= t) && (sa_size != 0)) {
      if (debug) cerr<<"4\n";
      string str_offsets="";
      set<int>::iterator it;

      int it_ant=t;
      //if (debug) cerr<<"Calculating offsets\n";
      vector<int> offsets;

      bool alignment_removed=false;

      for(it=ta_set.begin(); (it!=ta_set.end()) && (!alignment_removed); it++) {
        int offset=*it-it_ant;
        //if (debug) cerr<<"it_ant="<<it_ant<<"; it="<<*it<<"; offset="<<offset<<"\n";
      
        if (offset>=code_base) {
          //This offset cannot be coded
          if (debug) {
            cerr<<"Warning: Alignment between words '"<<alg->source[s]<<"' ("<<s<<") and '"<<alg->target[*it]<<"' ("<<*it<<")\n";
            cerr<<"Relative offset = "<<offset<<"; code base: "<<code_base<<"\n";
            cerr<<"s="<<s<<"; t="<<t<<"\n";
            cerr<<"ta_set={"<<Utils::set2string(ta_set)<<"}\n";
            cerr<<"sa_set={"<<Utils::set2string(sa_set)<<"}\n";
            cerr<<*alg<<"\n\n";
          }
          alg->alignment[s][*it]=false;
          alignment_removed=true;
        } else {
          offsets.push_back(offset);
          it_ant=*it+1;
        }
  
        //if (str_offsets.length()>0)
        //  str_offsets += "‡";
        //str_offsets+=Utils::itoa(offset);

      }

      if (alignment_removed)
        continue;


      if (debug) cerr<<"\n";
      //str += alg->source[s] + "|" + alg->build_str(ta_set, alg->target) + "|" + str_offsets + "\n";
      source_output.push_back(alg->source[s]);
      target_output.push_back(alg->build_str(ta_set, alg->target));

      int aux = Utils::code(code_base, offsets);

      if (debug) {
        cerr<<"offsets={"<<Utils::vector2string(offsets)<<"}; code="<<aux<<"\n";
      }

      if (aux==0)
        offsets_output.push_back("");
      else 
        offsets_output.push_back(Utils::itoa(aux));
      s++; 
      if ((ta_first==t) && (sa_first == s))
        t++;
    } else {
      cerr<<"Warning at print_alignment_offsets: This message should *NEVER* appear.\n";
      cerr<<"Write an e-mail to fsanchez@dlsi.ua.es with the following:\n";
      cerr<<"SL: "<<Utils::vector2string(source)<<"\n";
      cerr<<"TL: "<<Utils::vector2string(target)<<"\n";
      cerr<<"AL: "<<alignment_str<<"\n";
      exit(EXIT_FAILURE);
    }

    if (debug) {
      cerr<<"output\n-----------\n";
      string str="";
      for(unsigned i=0; i<source_output.size(); i++) {
        str += source_output[i] + "|" + target_output[i];
        if (offsets_output[i].length()>0)
          str += "|" + offsets_output[i];
        str+="\n";
      }
      cerr<<str;
    }
  }

  if (debug) cerr<<"Fin while\n";

  for (; t<(int)alg->target.size(); t++) {
    set<int> sa_set = alg->source_aligned(t);
    if (sa_set.size()==0) {
      //str += "|" + alg->target[t] + "\n";
      source_output.push_back("");
      target_output.push_back(alg->target[t]);
      offsets_output.push_back("");
    }
  }

  string str="";
  if (!one_to_many_alg) { //switch the output to reduce the number of different symbols
    vector<string> aux;

    aux=source_output;
    source_output=target_output;
    target_output=aux;
    str+="Ɛ|Ɛ|1\n";
  } else {
    str+="Ɛ|Ɛ\n";
  }

  for(unsigned i=0; i<source_output.size(); i++) {
    str += source_output[i] + "|" + target_output[i];
    if (offsets_output[i].length()>0)
      str += "|" + offsets_output[i];
    str+="\n";
  }

  return str;
}

void 
Alignment::get_biword_offsets(Alignment *alg, int &s, int &t, const set<int> &sa_set, const set<int> &ta_set, bool non_contiguous, 
                              vector<string> &source_output, vector<string> &target_output, vector<string> &offsets_output,
                              bool encode_offsets) {

  int ta_first=*ta_set.begin(); // ta_set is not empty
  int sa_first=*sa_set.begin(); // sa_set is not empty

  if (non_contiguous) {

    if (debug) cerr<<"get_biword_offsets, non_contiguous=true\n";
    string str_offsets="";

    int it_ant=t;
    vector<int> offsets;

    bool alignment_removed=false;
    set<int>::iterator it;
    for(it=ta_set.begin(); (it!=ta_set.end()) && (!alignment_removed); it++) {
      int offset=*it-it_ant;
      //if (debug) cerr<<"it_ant="<<it_ant<<"; it="<<*it<<"; offset="<<offset<<"\n";
      
      if ((offset>=code_base)&&(encode_offsets)) {
        //This offset cannot be coded
        if (debug) {
          cerr<<"Warning: Alignment between words '"<<alg->source[s]<<"' ("<<s<<") and '"<<alg->target[*it]<<"' ("<<*it<<")\n";
          cerr<<"Relative offset = "<<offset<<"; code base: "<<code_base<<"\n";
          cerr<<"s="<<s<<"; t="<<t<<"\n";
          cerr<<"ta_set={"<<Utils::set2string(ta_set)<<"}\n";
          cerr<<"sa_set={"<<Utils::set2string(sa_set)<<"}\n";
          cerr<<*alg<<"\n\n";
        }
        alg->alignment[s][*it]=false;
        alignment_removed=true;
        stats.nalignments_discarded++;
      } else {
        offsets.push_back(offset);
        it_ant=*it+1;
      }
    }

    if (alignment_removed)
      return;

    if (debug) cerr<<"\n";

    source_output.push_back(alg->source[s]);
    target_output.push_back(alg->build_str(ta_set, alg->target));

    if (encode_offsets) {
      int aux = Utils::code(code_base, offsets);

      if (debug) 
        cerr<<"offsets={"<<Utils::vector2string(offsets)<<"}; code="<<aux<<"\n";

      if (aux==0)
        offsets_output.push_back("");
      else 
        offsets_output.push_back(Utils::itoa(aux));

    } else { //Do not encode offsets
      string str_aux="";
      bool all_zero=true;
      for (unsigned i=0; i<offsets.size(); i++) {
        if (str_aux.length()>0)
          str_aux += ",";
        str_aux += Utils::itoa(offsets[i]);
        if (offsets[i]!=0)
          all_zero=false;
      }
      if (all_zero)
        offsets_output.push_back("");
      else
        offsets_output.push_back(str_aux);
    }

    s++; 

    if ((ta_first==t) && (sa_first == s))
      t++;
  } else {

    if (alg->consecutive(ta_set)) {
      source_output.push_back(alg->source[s]); 
      target_output.push_back(alg->build_str(ta_set, alg->target));

      if (ta_first == t) { //&& (sa_first == s)) {
        offsets_output.push_back("");
        t+=ta_set.size();
      } else { //ta_first > t
        int offset = ta_first - t;
        offsets_output.push_back(Utils::itoa(offset));
      }

      s++;
    } else {
      if (debug) cerr<<"4.2 - No consecutive !!!!!\n";
      alg->remove_least_frequent_alignment(s, ta_set, true);
    }
  }
}


string 
Alignment::biwords_with_offsets_new_version(bool switch_alg, double freqth, bool non_contiguous, bool encode_offsets, Alignment* alg2) {

  //cerr<<"Running Alignment::biwords_with_offsets_new_version\n\n";

  Alignment* alg = this; //Alignment to use

  // This function assumes 1:N (alg1) and N:1 (alg2) alignments
  if (alg2!=NULL) {
    if (!are_the_same_alignment(*alg2)) {
      cerr<<"Warning, following alignments are not equal, and therefore discarded:\n";
      cerr<<to_string()<<"\n";
      cerr<<alg2->to_string()<<"\n";
      return "";
    }

    //We use the alignment with the highest score
    //Warning: This comparison could be unfair, scores come from different
    //alignment models
    if (alg2->score>score)
      alg=alg2;
  }

  int s=0, t=0;

  vector<string> source_output, target_output, offsets_output;

  //First we need to know if we are dealing with one-to-many (1:N) or
  //many-to-one (N:1) alignments
  bool one_to_many_alg = alg->one_to_many();

  if((!one_to_many_alg) && (!switch_alg)) {
    cerr<<"Error: A many-to-one word alignment was found. Please re-run with --switch\n";
    exit(EXIT_FAILURE);
  }

  if (switch_alg) 
    alg->switch_alignment(); 

  while (s<(int)alg->source.size()) {
    set<int> ta_set = alg->target_aligned(s);
    set<int> sa_set = alg->source_aligned(t);

    int ta_size = ta_set.size();
    int sa_size = sa_set.size();

    int ta_first;
    int ta_last;
    if (ta_size > 0) {
      ta_first=*ta_set.begin();
      ta_last=*(--ta_set.end());
    } else {
      ta_first=-1;
      ta_last=-1;
    }

    int sa_first;
    int sa_last;
    if (sa_size > 0) {
      sa_first=*sa_set.begin();
      sa_last=*(--sa_set.end());
    } else {
      sa_first=-1;
      sa_last=-1;
    }

    if (sa_size > 1) {
      cerr<<"Warning: This is not a 1:N alignment\n";
      cerr<<to_string()<<"\n";
      exit(EXIT_FAILURE);
    }

    if (debug) {
      cerr<<"s="<<s<<"; t="<<t<<"\n";
      cerr<<"ta_set={"<<Utils::set2string(ta_set)<<"}\n"; //; consecutive="<<consecutive_ta<<"\n";
      cerr<<"sa_set={"<<Utils::set2string(sa_set)<<"}\n"; //; consecutive="<<consecutive_sa<<"\n";
      cerr<<"ta_first="<<ta_first<<"; ta_last="<<ta_last<<"\n";
      cerr<<"sa_first="<<sa_first<<"; sa_last="<<sa_last<<"\n";
      //cerr<<"Press ENTER "<<flush;
      //getchar();
    }

    if (ta_size == 0) { //s is not aligned
      if (debug) cerr<<"1\n";
      source_output.push_back(alg->source[s]);
      target_output.push_back("");
      offsets_output.push_back("");

      s++;
    } else if ((ta_first > t) && (sa_size == 0)) { //t is not aligned
      if (debug) cerr<<"2\n";
      source_output.push_back("");
      target_output.push_back(alg->target[t]);
      offsets_output.push_back("");

      t++;
    } else if ((sa_size > 0) && (ta_first > t) && (sa_last < s)) { //(sa_size > 0) se puede eliminar?
      if (debug) cerr<<"3\n";
      t++; // This alignment has been already printed
    } else if (ta_first >= t) { //((sa_size > 0) && (ta_first >= t)) {
      if (debug) cerr<<"4\n";

      if(freqth==-1) {
        get_biword_offsets(alg, s, t, sa_set, ta_set, non_contiguous, source_output, target_output, offsets_output, encode_offsets);
      } else {
        int aux_s=s;
        int aux_t=t;
        vector<string> aux_src_out;
        vector<string> aux_tgt_out;
        vector<string> aux_off_out;

        get_biword_offsets(alg, aux_s, aux_t, sa_set, ta_set, non_contiguous, aux_src_out, aux_tgt_out, aux_off_out, encode_offsets);

        if (aux_src_out.size()>0) {
          string src_str=aux_src_out[0];
          string tgt_str=aux_tgt_out[0];
          string off_str=aux_off_out[0];
          if (freq_biwords[src_str][tgt_str]<freqth) {
            alg->remove_least_frequent_alignment(s, ta_set, true);
          } else {
            s=aux_s;
            t=aux_t;
            source_output.push_back(src_str);
            target_output.push_back(tgt_str);
            offsets_output.push_back(off_str);
          }
        }
      }
    } else {
      cerr<<"Warning at Alignment::biwords_with_offsets_new_version: This message should *NEVER* appear.\n";
      cerr<<"Write an e-mail to fsanchez@dlsi.ua.es with the following:\n";
      cerr<<"SL: "<<Utils::vector2string(alg->source)<<"\n";
      cerr<<"TL: "<<Utils::vector2string(alg->target)<<"\n";
      cerr<<"AL: "<<alg->alignment_str<<"\n";
      exit(EXIT_FAILURE);
    }

    if (debug) {
      cerr<<"OUTPUT\n-------------\n";

      for(unsigned i=0; i<source_output.size(); i++) {
        cerr<<source_output[i]<<"|"<<target_output[i];
        if (offsets_output[i].length()>0)
          cerr<<"|"<<offsets_output[i];
        cerr<<"\n";
      }

      cerr<<"-------------\n";
    }
  }

  if (debug) cerr<<"Fin while\n";


  for (; t<(int)alg->target.size(); t++) {
    set<int> sa_set = alg->source_aligned(t);
    if (sa_set.size() == 0) {
      source_output.push_back("");
      target_output.push_back(alg->target[t]);
      offsets_output.push_back("");
    }
  }

  alg->update_count_alignments();
  alg->update_count_biwords(source_output, target_output);

  string str="";
  /*
    if (switch_alg) { //switch the output to reduce the number of different symbols
    vector<string> aux;

    aux=source_output;
    source_output=target_output;
    target_output=aux;
    if (freqth!=-1) str+="Ɛ|Ɛ|1\n";
    } else {
    if (freqth!=-1) str+="Ɛ|Ɛ\n";
    }  
  */
  if (freqth!=-1.0) str+="Ɛ|Ɛ\n";

  if (freqth!=-1.0) {
    for(unsigned i=0; i<source_output.size(); i++) {
      str += source_output[i] + "|" + target_output[i];
      if (offsets_output[i].length()>0)
        str += "|" + offsets_output[i];
      str+="\n";
    }
  }

  return str;
}

void 
Alignment::update_count_biwords(const vector<string>& source_output, const vector<string>& target_output) {

  if (source_output.size() != target_output.size()) {
    cerr<<"Error: In Alignment::update_count_biwords: "<<source_output.size()<<" "<<target_output.size()<<"\n";
  }

  for(unsigned i=0; i<source_output.size(); i++) {
    nbiwords++;
    count_biwords[source_output[i]][target_output[i]]++;
  }
}

int 
Alignment::compute_biwords_frequency() {
  map<string, map<string, int> >::iterator it1;
  map<string,int>::iterator it2;

  int ndiff_biwords=0;
  for (it1=count_biwords.begin(); it1!=count_biwords.end(); it1++) 
    for(it2=it1->second.begin(); it2!=it1->second.end(); it2++) {
      freq_biwords[it1->first][it2->first]=((double)it2->second)/((double)nbiwords);
      ndiff_biwords++;
    }

  nbiwords=0;
  count_biwords.clear();
  count_alignments.clear();

  return ndiff_biwords;
}

void 
Alignment::print_biwords_vocabulary() {
  map<string, map<string, int> >::iterator it1;
  map<string,int>::iterator it2;

  int ndiff_biwords=0;

  cerr<<"BIWORDS VOCABULARY\n-------------------------\n";
  for (it1=count_biwords.begin(); it1!=count_biwords.end(); it1++) 
    for(it2=it1->second.begin(); it2!=it1->second.end(); it2++) {
      cerr<<it2->second<<" "<<((double)it2->second)/((double)nbiwords)<<" "<<it1->first<<" | "<<it2->first<<"\n";
      ndiff_biwords++;
    }
  cerr<<"---------------------------------------\n";
  cerr<<"Size of the biwords vocabulary: "<<ndiff_biwords<<" "<<nbiwords<<"\n";
}



int 
Alignment::move_forward (int s, int t) {
  
  //Move t to the next unligned target word or to the next target word
  //aligned with a source word not yet processed 

  set<int> sa_set_aux;
  do {
    t++;

    sa_set_aux=source_aligned(t);
  } while ((t<(int)target.size()) && (sa_set_aux.size()!=0) && (*(sa_set_aux.begin())<s));

  return t;
}

string 
Alignment::biwords_with_offsets_new_version_rafa(bool non_contiguous, bool encode_offsets, Alignment* alg2) {

  //cerr<<"Running Alignment::biwords_with_offsets_new_version\n\n";

  Alignment* alg = this; //Alignment to use

  // This function assumes 1:N (alg1) and N:1 (alg2) alignments
  if (alg2!=NULL) {
    if (!are_the_same_alignment(*alg2)) {
      cerr<<"Warning, following alignments are not equal, and therefore discarded:\n";
      cerr<<to_string()<<"\n";
      cerr<<alg2->to_string()<<"\n";
      return "";
    }

    //We use the alignment with the highest score
    //Warning: This comparison could be unfair, scores come from different
    //alignment models
    if (alg2->score>score)
      alg=alg2;
  }

  int s=0, t=0;

  vector<string> source_output, target_output, offsets_output;

  //First we need to know if we are dealing with one-to-many (1:N) or
  //many-to-one (N:1) alignments
  bool one_to_many_alg = alg->one_to_many();

  if (!one_to_many_alg) {
    alg->switch_alignment(); 
    //cerr<<"MANY TO ONE (N:1)\n";
  } else {
    //cerr<<"ONE TO MANY (1:N)\n";
  }

  while ((s<(int)alg->source.size()) && ((t<(int)alg->target.size()))) {
    set<int> ta_set = alg->target_aligned(s);
    set<int> sa_set = alg->source_aligned(t);

    int ta_size = ta_set.size();
    int sa_size = sa_set.size();

    int ta_first;
    int ta_last;
    if (ta_size > 0) {
      ta_first=*ta_set.begin();
      ta_last=*(--ta_set.end());
    } else {
      ta_first=-1;
      ta_last=-1;
    }

    int sa_first;
    int sa_last;
    if (sa_size > 0) {
      sa_first=*sa_set.begin();
      sa_last=*(--sa_set.end());
    } else {
      sa_first=-1;
      sa_last=-1;
    }

    if (sa_size > 1) {
      cerr<<"Warning: This is not a 1:N alignment\n";
      cerr<<to_string()<<"\n";
      exit(EXIT_FAILURE);
    }

    if (debug) {
      cerr<<"s="<<s<<"; t="<<t<<"\n";
      cerr<<"ta_set={"<<Utils::set2string(ta_set)<<"}\n"; //; consecutive="<<consecutive_ta<<"\n";
      cerr<<"sa_set={"<<Utils::set2string(sa_set)<<"}\n"; //; consecutive="<<consecutive_sa<<"\n";
      cerr<<"ta_first="<<ta_first<<"; ta_last="<<ta_last<<"\n";
      cerr<<"sa_first="<<sa_first<<"; sa_last="<<sa_last<<"\n";
      //cerr<<"Press ENTER "<<flush;
      //getchar();
    }

    if (ta_size == 0) { //s is not aligned
      if (debug) cerr<<"1\n";
      source_output.push_back(alg->source[s]);
      target_output.push_back("");
      offsets_output.push_back("");

      s++;
    } else if ((sa_size == 0)) { // && (ta_first>t)) { //t is not aligned
      if (debug) cerr<<"2\n";
      source_output.push_back("");
      target_output.push_back(alg->target[t]);
      offsets_output.push_back("");


      if (debug) cerr<<"Actualizando t desde "<<t<<"\n";
      t=move_forward(s,t);
      if (debug) cerr<<"Actualizado t a "<<t<<"\n";
    } else { 
      if (debug) cerr<<"3\n";
      string str_offsets="";
      int it_ant=t;
      vector<int> offsets;
      set<int>::iterator it;
      for(it=ta_set.begin(); it!=ta_set.end(); it++) {
        int offset=*it-it_ant;
        if (debug) cerr<<"it_ant="<<it_ant<<"; it="<<*it<<"; offset="<<offset<<"\n";
      
        offsets.push_back(offset);
        it_ant=*it+1;

        if ((*it) == t) {
          if (debug) cerr<<"4\n";
          if (debug) cerr<<"Actualizando t desde "<<t<<"\n";
          t=move_forward(s,t);
          if (debug) cerr<<"Actualizado t a "<<t<<"\n";
        }
      }

      source_output.push_back(alg->source[s]);
      target_output.push_back(alg->build_str(ta_set, alg->target));

      string str_aux="";
      bool all_zero=true;
      for (unsigned i=0; i<offsets.size(); i++) {
        if (str_aux.length()>0)
          str_aux += ",";
        str_aux += Utils::itoa(offsets[i]);
        if (offsets[i]!=0)
          all_zero=false;
      }

      if (all_zero)
        offsets_output.push_back("");
      else
        offsets_output.push_back(str_aux);

      s++;

      //get_biword_offsets(alg, s, t, sa_set, ta_set, non_contiguous, source_output, target_output, offsets_output, encode_offsets);
    }

    if (debug) {
      cerr<<"OUTPUT\n-------------\n";

      for(unsigned i=0; i<source_output.size(); i++) {
        cerr<<source_output[i]<<"|"<<target_output[i];
        if (offsets_output[i].length()>0)
          cerr<<"|"<<offsets_output[i];
        cerr<<"\n";
      }

      cerr<<"-------------\n\n";
    }
  }

  if (debug) cerr<<"Fin while\n";

  if (debug) cerr<<"s: "<<s<<"; t: "<<t<<endl;
  

  for (; s<(int)alg->source.size(); s++) {
    source_output.push_back(alg->source[s]);
    target_output.push_back("");
    offsets_output.push_back("");
  }

  for (; t<(int)alg->target.size(); t++) {
    set<int> sa_set = alg->source_aligned(t);
    if (sa_set.size() == 0) {
      source_output.push_back("");
      target_output.push_back(alg->target[t]);
      offsets_output.push_back("");
    }
  }

  alg->update_count_alignments();

  string str="";
  if (!one_to_many_alg) { //switch the output to reduce the number of different symbols
    vector<string> aux;

    aux=source_output;
    source_output=target_output;
    target_output=aux;
    str+="Ɛ|Ɛ|1\n";
  } else {
    str+="Ɛ|Ɛ\n";
  }

  for(unsigned i=0; i<source_output.size(); i++) {
    str += source_output[i] + "|" + target_output[i];
    if (offsets_output[i].length()>0)
      str += "|" + offsets_output[i];
    str+="\n";
  }

  return str;
}

string 
Alignment::biwords_without_offsets(Alignment* alg2) {

  // This function assumes that the alignments are 1:1
  // This assumption is not checked

  string str="Ɛ|Ɛ\n";
  int s=0, t=0;

  while (s<(int)source.size()) {
    int ta = *(target_aligned(s).begin());
    if (ta == -1) {
      str += source[s] + "|" + "\n";
      s++;
    } else if (ta == t) {
      str += source[s] + "|" + target[t] + "\n";
      s++; t++;
    } else if (ta<t) {
      //Ignore this alignment. There is a crossing in alignments
      //target[ta] has been already printed out
      str += source[s] + "|" + "\n";
      s++;
      stats.nalignments_discarded++;
    } else {
      for(; t<ta; t++) {
        str += "|" + target[t] + "\n";
        //t++;
      }
    }
  }

  for(; t<(int)target.size(); t++) {
    str += "|" + target[t] + "\n";
  }

  update_count_alignments();
  return str;
}


bool 
Alignment::one_to_many() {

  for(unsigned  t=0; t<target.size(); t++) {
    if (source_aligned(t).size()>1)
      return false;
  }
  return true;
}

void 
Alignment::switch_alignment() {
  vector<string> aux;

  unsigned source_size = source.size();
  unsigned target_size = target.size();

  aux = source;
  source = target;
  target = aux;

  map<int, map<int, bool> > new_alig;

  for(unsigned i=0; i<source_size; i++) {
    for(unsigned j=0; j<target_size; j++) {
      new_alig[j][i]=alignment[i][j];
    }
  }

  alignment=new_alig;
}

void 
Alignment::collect_stats() {
  for(unsigned i=0; i<source.size(); i++) {
    set<int> ta = target_aligned(i);
    stats.source_nwords++;
    if (*(ta.begin())==-1)
      stats.alig_source[0]++;
    else
      stats.alig_source[ta.size()]++;
  }

  for(unsigned i=0; i<target.size(); i++) {
    set<int> sa = source_aligned(i);
    stats.target_nwords++;
    if (*(sa.begin())==-1)
      stats.alig_target[0]++;
    else
      stats.alig_target[sa.size()]++;
  }
}

unsigned Alignment::unaligned_source_words() {
  unsigned unaligned=0;

  for(unsigned i=0; i<source.size(); i++) {
    bool aligned=false;
    for(unsigned j=0; (j<target.size()) && (!aligned); j++) {
      if (alignment[i][j])
        aligned=true;
    }
    if (!aligned)
      unaligned++;
  }

  return unaligned;
}

unsigned Alignment::unaligned_target_words() {
  unsigned unaligned=0;

  for(unsigned j=0; j<target.size(); j++) {
    bool aligned=false;
    for(unsigned i=0; (i<source.size()) && (!aligned); i++) {
      if (alignment[i][j])
        aligned=true;
    }
    if (!aligned)
      unaligned++;
  }

  return unaligned;
}

void 
Alignment::print_stats() {
  cout<<"Source language\n-------------------------------\n";
  map<int, int>::iterator itmap;
  for(itmap=stats.alig_source.begin(); itmap!=stats.alig_source.end(); itmap++)
    cout<<itmap->first<<" "<<itmap->second<<" "
        <<((double)itmap->second)/((double)stats.source_nwords)<<"\n";
  cout<<"\n";

  cout<<"Target language\n-------------------------------\n";
  for(itmap=stats.alig_target.begin(); itmap!=stats.alig_target.end(); itmap++)
    cout<<itmap->first<<" "<<itmap->second<<" "
        <<((double)itmap->second)/((double)stats.target_nwords)<<"\n";
  cout<<"\n";

}

void 
Alignment::print_alig_stats() {
  cerr<<"# of word alignments in total:  "<<stats.nalignments<<"\n";
  cerr<<"# of word alignments discarded: "<<stats.nalignments_discarded<<"\n";
  cerr<<"% of word alignments discarded: "<<(((double)stats.nalignments_discarded)/((double)stats.nalignments))*100<<" %\n";
}

set<int>
Alignment::target_aligned(int s) {
  set<int> ta;
  for(unsigned j=0; j<target.size(); j++) {
    if (alignment[s][j]) 
      ta.insert(j);
  }

  return ta;
}

set<int>
Alignment::source_aligned(int t) {
  set<int> sa;
  for(unsigned i=0; i<source.size(); i++) {
    if (alignment[i][t])
      sa.insert(i);
  }

  return sa;
}

bool 
Alignment::consecutive(const set<int>& SP) {
  set<int>::iterator it;

  if (SP.size()==0)
    return false;

  it=SP.begin();
  int prev=(*it);

  for(it++; it!=SP.end(); it++) {
    if (((*it)-prev)!=1)
      return false;
    prev=(*it);
  }

  return true;
}


string 
Alignment::build_str(const set<int>& SP, const vector<string>& words) {
  string str="";
  set<int>::iterator it;

  for(it=SP.begin(); it!=SP.end(); it++) {
    if (str.length()>0)
      str += " ";
    str += words[*it];
  }

  return str;
}

void 
Alignment::update_count_alignments() {

  for(unsigned i=0; i<source.size(); i++) {
    for (unsigned j=0; j<target.size(); j++) {
      if (alignment[i][j])
        count_alignments[source[i]][target[j]]++;
    }
  }
}

void 
Alignment::remove_least_frequent_alignment(int frompos, const set<int>& topos, bool from_source) {

  int min_count=numeric_limits<int>::max();
  int min_pos=*(--topos.end());

  set<int>::iterator it;
  string from_str;

  if (from_source) {
    from_str = source[frompos];
    string to_str = target[min_pos];
    min_count = count_alignments[from_str][to_str];
  } else {
    from_str = target[frompos];
    string to_str = source[min_pos];
    min_count = count_alignments[to_str][from_str];
  }

  if (debug) cerr<<"min_pos="<<min_pos<<"; min_count="<<min_count<<"\n";

  for(it=topos.begin(); it!=topos.end(); it++) {
    int count;

    if (from_source) {
      string to_str = target[*it];
      count = count_alignments[from_str][to_str];
    } else {
      string to_str = source[*it];
      count = count_alignments[to_str][from_str];
    }

    if (count<min_count) {
      min_count = count;
      min_pos = *it;
    } 
  }

  if (from_source) {
    if (debug) 
      cerr<<"Removing alignment between '"<<source[frompos]<<"' ("<<frompos<<") and '"<<target[min_pos]<<"' ("<<min_pos<<")\n";
    alignment[frompos][min_pos]=false;
  } else {
    if (debug) 
      cerr<<"Removing alignment between '"<<source[min_pos]<<"' ("<<min_pos<<") and '"<<target[frompos]<<"' ("<<frompos<<")\n";
    alignment[min_pos][frompos]=false;
  }
  stats.nalignments_discarded++;
}

ostream& operator << (ostream& os, Alignment& al) {

  os<<al.to_string();
  return os;
}

bool 
Alignment::are_the_same_alignment(const Alignment& al2) {
  bool ok=true;

  if ((source.size()!=al2.source.size())||
      (target.size()!=al2.target.size())) {
    ok=false;
  }

  for(unsigned i=0; (i<source.size()) && (ok); i++) {
    if (source[i]!=al2.source[i])
      ok=false;
  }

  for(unsigned i=0; (i<target.size()) && (ok); i++) {
    if (target[i]!=al2.target[i])
      ok=false;
  }
  return ok;
}

bool 
Alignment::intersection(Alignment& al2) {
  if (!are_the_same_alignment(al2)) {
    cerr<<"Error when intersecting the following two alignments:\n";
    cerr<<to_string()<<"\n";
    cerr<<al2.to_string()<<"\n";
    return false;
  }

  for(unsigned i=0; i<source.size(); i++) {
    for(unsigned j=0; j<target.size(); j++) {
      if ((alignment[i][j]) && (!al2.alignment[i][j])) {
        alignment[i][j]=false;
      }
    }
  }

  return true;
}

bool
Alignment::unionn(Alignment& al2) {
  if (!are_the_same_alignment(al2)) {
    cerr<<"Error when uniting the following two alignments:\n";
    cerr<<to_string()<<"\n";
    cerr<<al2.to_string()<<"\n";
    return false;
  }

  for(unsigned i=0; i<source.size(); i++) {
    for(unsigned j=0; j<target.size(); j++) {
      if (al2.alignment[i][j]) {
        alignment[i][j]=true;
      }
    }
  }

  return true;
}

bool 
Alignment::refined_intersection(Alignment& al2) {
  if (!are_the_same_alignment(al2)) {
    cerr<<"Error when performing the refined intersection of the following two alignments:\n";
    cerr<<to_string()<<"\n";
    cerr<<al2.to_string()<<"\n";
    return false;
  }

  //We save the alignment of (*this) before intersecting
  map<int, map<int, bool> > al1;
  for(unsigned i=0; i<source.size(); i++) {
    for(unsigned j=0; j<target.size(); j++) {
      al1[i][j]=alignment[i][j];
    }
  }

  //First, the intersection
  intersection(al2);
  //cerr<<"(0) "<<to_string()<<"\n";
  bool alignments_changed;
  //int nit=0;
  do {
    alignments_changed=false;
    //nit++;
    for (unsigned i=0; i<source.size(); i++) {
      for(unsigned j=0; j<target.size(); j++) {
        if (!alignment[i][j]) {
          if ((al1[i][j]) || (al2.alignment[i][j])) {
            //cerr<<"   considering ("<<i<<","<<j<<")\n";
            bool add_alignment=true;
            for(unsigned k=0; k<target.size(); k++) {
              if (alignment[i][k])
                add_alignment=false;
            }
            for(unsigned k=0; k<source.size(); k++) {
              if (alignment[k][j])
                add_alignment=false;
            }
            if (!add_alignment) {
              //cerr<<"   ("<<i<<","<<"*) or (*,"<<j<<") are already in A\n";
              //The aligment can still be added if it has an horizontal
              //neighbor or a vertical neighbor already in 'alignment'
              if ( ((alignment[i-1][j])||(alignment[i+1][j])) ||
                   ((alignment[i][j-1])||(alignment[i][j+1])) ) {
                //cerr<<"   ("<<i<<","<<j<<") has an horizontal or a vertical neighbor in A\n";
                alignment[i][j]=true; 
                //Now we that test there is no aligments in 'alignment' with
                //both horizontal and vertical neighbors
                bool neighbors=false;
                for(unsigned ii=0; (ii<source.size()) && (!neighbors); ii++) {
                  for(unsigned jj=0; (jj<target.size()) && (!neighbors); jj++) {
                    if (alignment[ii][jj]) {
                      if ( ((alignment[ii-1][jj])||(alignment[ii+1][jj])) &&
                           ((alignment[ii][jj-1])||(alignment[ii][jj+1])) ) {
                        //cerr<<"   ("<<i<<","<<j<<") has both horizontal and vertical neighbors\n";
                        neighbors=true;
                      }
                    }
                  }
                }
                alignment[i][j]=false;
                if(!neighbors)
                  add_alignment=true;
              }
            }
            if (add_alignment) {
              //cerr<<"   adding ("<<i<<","<<j<<")\n";
              alignment[i][j]=true;
              alignments_changed=true;
            }
          }
        }
      }
    }
    //cerr<<"("<<nit<<") "<<to_string()<<"\n";
  } while(alignments_changed);

  return true;
}
