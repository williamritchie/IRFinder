#ifndef CODE_SEQUENCETOOLS
#define CODE_SEQUENCETOOLS

#include "includedefine.h"

uint localAlign(const char *, uint, const char *, uint ny, uint minspan, uint maxspan);
void seqToNum(const char*, char*, uint); // do we really need length, or just run until \0?
void seqToNumComp(const char*, char*, uint);
void seqToNumRevComp(const char*, char*, uint);
void seqToNumRev(const char*, char*, uint);

std::string timeMonthDayTime();
std::string timeMonthDayTime(time_t &rawTime);

#endif
