#include "sequenceTools.h"

uint localAlign(const char *x, uint nx, const char *y, uint ny, uint minspan, uint maxspan)
{
  // Expecting x, a numeric Seq string. y, a pre-complemented numeric Seq string. Both are expected in their original direction.
  uint nMatch;
  int nMMcounter;
  double nScore;
  double nScoreBest=0;
  uint spanBest=maxspan;
  uint ixbegin;
  uint ixlimit;
  // Min score = 0.8 therefore, at worst position (integer comparison):
  //uint maxMismatch = (maxspan/10) + 1;
  uint maxMismatch;

  for (uint span=minspan; span<=maxspan; span++ ) {
    nMatch=0;
    ixbegin = max(0,int(span)-int(ny));
    ixlimit = min(span,nx);
    maxMismatch = ((ixlimit-ixbegin)/10) + 1;
    nMMcounter=maxMismatch;
    for (uint ix=ixbegin; ix<ixlimit; ix++) {
      char cy = y[span-ix-1];
      if (x[ix] == 0 || cy == 0) continue;

      if (x[ix] == cy) {
        nMatch++;
      }else{
        nMMcounter--;
      }

      if (nMMcounter < 0) break;
    }
    if (nMMcounter >= 0) {
      nScore = ((double)(nMatch - maxMismatch + nMMcounter))/(ixlimit-ixbegin);
      if (nScore >= 0.8 && nScore > nScoreBest) {
        nScoreBest = nScore;
        spanBest=span;
      }
    }
  }

  return spanBest;
}

void seqToNum(const char* in, char* out, uint nin) // do we really need length, or just run until \0?
{
  for (uint jj=0;jj<nin;jj++) {
    switch(in[jj]){
//       case ('N'): case ('n'): case ('.'):  out[jj]=char(0);break;
       case ('A'): case ('a'):  out[jj]=char(1);break;
       case ('T'): case ('t'): case ('U'): case ('u'):  out[jj]=char(2);break;
       case ('C'): case ('c'):  out[jj]=char(3);break;
       case ('G'): case ('g'):  out[jj]=char(4);break;
       default:  out[jj]=char(0);
    }
/*    switch(in[jj]){
       case ('N'): case ('n'): case ('.'):  out[jj]=char(0);break;
       case ('A'): case ('a'):  out[jj]='A';break;
       case ('T'): case ('t'):  out[jj]='T';break;
       case ('C'): case ('c'):  out[jj]='C';break;
       case ('G'): case ('g'):  out[jj]='G';break;
       default: out[jj]=char(99);
    }*/
  }
}

void seqToNumComp(const char* in, char* out, uint nin)
{
  for (uint jj=0;jj<nin;jj++) {
    switch(in[jj]){
//       case ('N'): case ('n'): case ('.'):  out[jj]=char(0);break;
       case ('A'): case ('a'):  out[jj]=char(2);break;
       case ('T'): case ('t'): case ('U'): case ('u'):  out[jj]=char(1);break;
       case ('C'): case ('c'):  out[jj]=char(4);break;
       case ('G'): case ('g'):  out[jj]=char(3);break;
       default: out[jj]=char(0);
    }
/*    switch(in[jj]){
       case ('N'): case ('n'): case ('.'):  out[jj]=char(0);break;
       case ('A'): case ('a'):  out[jj]='T';break;
       case ('T'): case ('t'):  out[jj]='A';break;
       case ('C'): case ('c'):  out[jj]='G';break;
       case ('G'): case ('g'):  out[jj]='C';break;
       default: out[jj]=char(99);
    }*/
  }
}

void seqToNumRevComp(const char* in, char* out, uint nin)
{
  for (uint jj=0;jj<nin;jj++) {
    switch(in[jj]){
//       case ('N'): case ('n'): case ('.'):  out[nin-jj-1]=char(0);break;
       case ('A'): case ('a'):  out[nin-jj-1]=char(2);break;
       case ('T'): case ('t'): case ('U'): case ('u'):  out[nin-jj-1]=char(1);break;
       case ('C'): case ('c'):  out[nin-jj-1]=char(4);break;
       case ('G'): case ('g'):  out[nin-jj-1]=char(3);break;
       default: out[nin-jj-1]=char(0);
    }
/*    switch(in[jj]){
       case ('N'): case ('n'): case ('.'):  out[nin-jj-1]=char(0);break;
       case ('A'): case ('a'):  out[nin-jj-1]='T';break;
       case ('T'): case ('t'):  out[nin-jj-1]='A';break;
       case ('C'): case ('c'):  out[nin-jj-1]='G';break;
       case ('G'): case ('g'):  out[nin-jj-1]='C';break;
       default: out[nin-jj-1]=char(99);
    }*/
  }
}
void seqToNumRev(const char* in, char* out, uint nin)
{
  for (uint jj=0;jj<nin;jj++) {
    switch(in[jj]){
//       case ('N'): case ('n'): case ('.'):  out[nin-jj-1]=char(0);break;
       case ('A'): case ('a'):  out[nin-jj-1]=char(1);break;
       case ('T'): case ('t'): case ('U'): case ('u'):  out[nin-jj-1]=char(2);break;
       case ('C'): case ('c'):  out[nin-jj-1]=char(3);break;
       case ('G'): case ('g'):  out[nin-jj-1]=char(4);break;
       default: out[nin-jj-1]=char(0);
    }
/*    switch(in[jj]){
       case ('N'): case ('n'): case ('.'):  out[nin-jj-1]=char(0);break;
       case ('A'): case ('a'):  out[nin-jj-1]='A';break;
       case ('T'): case ('t'):  out[nin-jj-1]='T';break;
       case ('C'): case ('c'):  out[nin-jj-1]='C';break;
       case ('G'): case ('g'):  out[nin-jj-1]='G';break;
       default: out[nin-jj-1]=char(99);
    }*/
  }
}

std::string timeMonthDayTime() {
    time_t rawTime;
    char timeChar[100];
    time(&rawTime);
    strftime(timeChar,80,"%b %d %H:%M:%SS",localtime(&rawTime));
    std::string timeString=timeChar;
    timeString.erase(timeString.end()-1,timeString.end());
    return timeString;
};

std::string timeMonthDayTime(time_t &rawTime) {
    char timeChar[100];
    strftime(timeChar,80,"%b %d %H:%M:%SS",localtime(&rawTime));
    std::string timeString=timeChar;
    timeString.erase(timeString.end()-1,timeString.end());
    return timeString;
};

