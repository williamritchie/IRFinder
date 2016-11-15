#include "TrimReads.h"


TrimReads::TrimReads (istream * _IF1, istream * _IF2, ostream * _OF1, ostream * _OF2, ostream * _DebugLog, char * Adapt1, char * Adapt2, char _debug)
{
  IF1 = _IF1;
  IF2 = _IF2;
  OF1 = _OF1;
  OF2 = _OF2;
  OFdebug = _DebugLog;
  debug = _debug;

  time(&startTime);
  cout << "Started trimming at: " << timeMonthDayTime(startTime) << endl;


  Compare1 = new char[DEF_lineLengthMax + DEF_adaptLengthMax];
  Compare2 = new char[DEF_lineLengthMax + DEF_adaptLengthMax];

  lAdapt1 = strlen(Adapt1);
  lAdapt2 = strlen(Adapt2);

  seqToNumRevComp(Adapt2, Compare1, lAdapt2);
  seqToNumRev(Adapt1, Compare2, lAdapt1);
  
  Compare1Seq = Compare1 + lAdapt2; //Set pointer inside Compare1 just after the adapter - copy Sequence into this location each time.
  Compare2Seq = Compare2 + lAdapt1;

  lAdapt1prefix = 0;
  while (lAdapt1prefix <= lAdapt1) {
  	if (Compare2[lAdapt1-1-lAdapt1prefix] != 0) break;
  	lAdapt1prefix++;
  }
  lAdapt2prefix = 0;
  while (lAdapt2prefix <= lAdapt2) {
  	if (Compare1[lAdapt2-1-lAdapt2prefix] != 0) break;
  	lAdapt2prefix++;
  }

  minSpan = lAdapt1 + lAdapt2 - lAdapt1prefix - lAdapt2prefix;

  chunk1 = new char[DEF_lineLengthMax*4];
  chunk2 = new char[DEF_lineLengthMax*4];

  countInputReads = 0;
  countOutputTrimmed = 0;
  countOutputUntrimmed = 0;
  countOutputTooShort = 0;

};

int TrimReads::trimAll()
{
//  *OFdebug << "trimAll" << endl;
  while(IF1->peek() == '@') {

  countInputReads++;

  uint chunkLen1 = 0;
  uint chunkLen2 = 0;

  IF1->getline(chunk1, DEF_lineLengthMax);
  chunkLen1 += IF1->gcount();
  chunk1[chunkLen1-1] = '\n';
  Seq1 = chunk1 + chunkLen1;
  uint lName1 = IF1->gcount() - 1;

  IF1->getline(chunk1 + chunkLen1, DEF_lineLengthMax);
  chunkLen1 += IF1->gcount();
  chunk1[chunkLen1-1] = '\n';
  uint lR1 = IF1->gcount() - 1;

  IF1->getline(chunk1 + chunkLen1, DEF_lineLengthMax);
  chunkLen1 += IF1->gcount();
  chunk1[chunkLen1-1] = '\n';
  uint lQName1 = IF1->gcount() - 1;

  IF1->getline(chunk1 + chunkLen1, DEF_lineLengthMax);
  chunkLen1 += IF1->gcount();
  chunk1[chunkLen1-1] = '\n';
  uint lQual1 = IF1->gcount() - 1;


  IF2->getline(chunk2, DEF_lineLengthMax);
  chunkLen2 += IF2->gcount();
  chunk2[chunkLen2-1] = '\n';
  Seq2 = chunk2 + chunkLen2;
  uint lName2 = IF2->gcount() - 1;

  IF2->getline(chunk2 + chunkLen2, DEF_lineLengthMax);
  chunkLen2 += IF2->gcount();
  chunk2[chunkLen2-1] = '\n';
  uint lR2 = IF2->gcount() - 1;

  IF2->getline(chunk2 + chunkLen2, DEF_lineLengthMax);
  chunkLen2 += IF2->gcount();
  chunk2[chunkLen2-1] = '\n';
  uint lQName2 = IF2->gcount() - 1;

  IF2->getline(chunk2 + chunkLen2, DEF_lineLengthMax);
  chunkLen2 += IF2->gcount();
  chunk2[chunkLen2-1] = '\n';
  uint lQual2 = IF2->gcount() - 1;


  if (lR1 != lQual1 || lR2 != lQual2) {
    //*OFdebug << "FATAL: a quality line has different length to Sequence line - corrupt input file." << endl;
    cout << "FATAL: a quality line has different length to Sequence line - corrupt input file." << endl;
    cout << "Error at read number: " << countInputReads << endl;
    return(1);
  }

  seqToNum(Seq1, Compare1Seq, lR1);
  seqToNumComp(Seq2, Compare2Seq, lR2);

  uint maxOverlap = minSpan+max(lR1+lAdapt2prefix, lR2+lAdapt1prefix);
//  uint overlapPos = localAlign(Compare1, lR1+lAdapt2, Compare2, lR2+lAdapt1, lAdapt1+lAdapt2, maxOverlap);
  uint overlapPos = localAlign(Compare1, lR1+lAdapt2, Compare2, lR2+lAdapt1, minSpan, maxOverlap);

  uint InsertLen = overlapPos-lAdapt1-lAdapt2+lAdapt1prefix+lAdapt2prefix;

  if (insertLenDistribution.size() <= InsertLen) {
    insertLenDistribution.resize( InsertLen+1 ,0);
  }

  insertLenDistribution.at(InsertLen)++;

  if (maxOverlap == overlapPos) {
    //No trimming, output exactly what we read - as a single fast block.
    OF1->write(chunk1, chunkLen1);
    OF2->write(chunk2, chunkLen2);
    countOutputUntrimmed++;
  }else if (InsertLen < 30) {
    // Don't output, too short.
    countOutputTooShort++;
  }else{
  	if (!debug) {
		chunk1[lName1] = '\0';
		uint firstSpace = strcspn(chunk1," ");
		chunk1[lName1] = '\n';
		OF1->write(chunk1, firstSpace);
		*OF1 << "_INS_" << InsertLen;
		OF1->write(chunk1+firstSpace, lName1+1-firstSpace+min(InsertLen-lAdapt2prefix,lR1));
		OF1->write(chunk1+lName1+lR1+1, lQName1+2+min(InsertLen-lAdapt2prefix,lR1));
		OF1->write("\n",1);

		chunk2[lName2] = '\0';
		firstSpace = strcspn(chunk2," ");
		chunk2[lName2] = '\n';
		OF2->write(chunk2, firstSpace);
		*OF2 << "_INS_" << InsertLen;
		OF2->write(chunk2+firstSpace, lName2+1-firstSpace+min(InsertLen-lAdapt1prefix,lR2));
		OF2->write(chunk2+lName2+lR2+1, lQName2+2+min(InsertLen-lAdapt1prefix,lR2));
		OF2->write("\n",1);
	}else{
		chunk1[lName1] = '\0';
		uint firstSpace = strcspn(chunk1," ");
		chunk1[lName1] = '\n';
		OF1->write(chunk1, firstSpace);
		*OF1 << "_INS_" << InsertLen;

		for(int i = lName1+1+min(InsertLen-lAdapt2prefix,lR1); i<(lName1+1+lR1); i++){
		  chunk1[i] = tolower(chunk1[i]);
		}		
		OF1->write(chunk1+firstSpace, lName1+1-firstSpace+lR1);
		OF1->write(chunk1+lName1+lR1+1, lQName1+2+lR1);
		OF1->write("\n",1);

		chunk2[lName2] = '\0';
		firstSpace = strcspn(chunk2," ");
		chunk2[lName2] = '\n';
		OF2->write(chunk2, firstSpace);
		*OF2 << "_INS_" << InsertLen;

		for(int i = lName2+1+min(InsertLen-lAdapt1prefix,lR2); i<(lName2+1+lR2); i++){
		  chunk2[i] = tolower(chunk2[i]);
		}		
		OF2->write(chunk2+firstSpace, lName2+1-firstSpace+lR2);
		OF2->write(chunk2+lName2+lR2+1, lQName2+2+lR2);
		OF2->write("\n",1);	
	}

    countOutputTrimmed++;
  }


  }

  // Output summary statistics.
  time(&endTime);

  cout << "Completed trimming at: " << timeMonthDayTime(endTime) << endl;
  cout << double(countInputReads)/1e6/difftime(endTime,startTime)*3600 << "\tTrimming speed, million reads per hour" << endl;

  ios::fmtflags old_output_settings = cout.flags();
  cout << fixed << setprecision(4);
  cout << (countOutputTrimmed+countOutputTooShort)/double(countInputReads)*100 << "\t% with adaptor" << endl;
  cout.flags(old_output_settings);
  cout << countInputReads << "\tTotal input reads" << endl;
  cout << (countOutputUntrimmed+countOutputTrimmed) << "\tTotal output reads" << endl;
  cout << countOutputUntrimmed << "\tTotal unmodified reads output" << endl;
  cout << countOutputTrimmed << "\tTotal trimmed reads output" << endl;
  cout << countOutputTooShort << "\tTotal trimmed reads too short" << endl;
  cout << endl;
  cout << "------ Insert length distribution ------" << endl;
  cout << "Length\tCount" << endl;

  for (uint i=0; i<insertLenDistribution.size() ; i++) {
    cout << i << "\t" << insertLenDistribution.at(i) << endl;
  }

  return 0;
};
