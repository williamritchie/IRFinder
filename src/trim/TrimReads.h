#ifndef CODE_TRIMREADS
#define CODE_TRIMREADS

#include "includedefine.h"
#include "sequenceTools.h"

class TrimReads {
  char* chunk1;
  char* chunk2;
  char* Seq1;
  char* Seq2;
  char* Compare1;
  char* Compare1Seq;
  char* Compare2;
  char* Compare2Seq;

  int lAdapt1;
  int lAdapt2;
  int lAdapt1prefix;
  int lAdapt2prefix;
  int minSpan;

  unsigned long int countInputReads;
  unsigned long int countOutputTrimmed;
  unsigned long int countOutputUntrimmed;
  unsigned long int countOutputTooShort;

  std::vector<unsigned long int> insertLenDistribution;

  time_t startTime, endTime;

  istream* IF1;
  istream* IF2;
  ostream* OF1;
  ostream* OF2;
  ostream* OFdebug;
  char debug;

  public:
    TrimReads(istream * _IF1, istream * _IF2, ostream * _OF1, ostream * _OF2, ostream * _DebugLog, char * Adapt1, char * Adapt2, char _debug = 0);
    int trimAll();
};

#endif
