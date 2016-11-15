#include "includedefine.h"

#include "sequenceTools.h"
#include "TrimReads.h"


int main(int argc, char * argv[])
{
  char debug = 0;
  
  if (argc != 7) {
  	if (argc == 8 && strcmp(argv[7], "debug") == 0) {
  		debug = 1;
  	}else{
    	cerr << "Usage: cmd in_1.fastq in_2.fastq out_1.fastq out_2.fastq adapt1 adapt2 [debug]" << endl;
    	exit(1);
    } 
  }

  ifstream IN1;
  IN1.open (argv[1], ifstream::in);
  ifstream IN2;
  IN2.open (argv[2], ifstream::in);
  ofstream OUT1;
  OUT1.open (argv[3], ifstream::out);
  ofstream OUT2;
  OUT2.open (argv[4], ifstream::out);

  TrimReads * TR = new TrimReads(&IN1, &IN2, &OUT1, &OUT2, &cerr, argv[5], argv[6], debug);
  int success = TR->trimAll();

  OUT1.flush();
  OUT2.flush();

  IN1.close();
  IN2.close();

  OUT1.close();
  OUT2.close();

  exit(success);
//  f1InStream.getline
};
