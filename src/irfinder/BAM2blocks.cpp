// WARNING: code is little endian only!

#include "BAM2blocks.h"
#include "includedefine.h"
// using namespace std;

//const char cigarChar[] = {'M','I','D','N','S','H','P','=','X'};

BAM2blocks::BAM2blocks() {
	oBlocks = FragmentBlocks(); //Right syntax to call the default constructor on an object variable, declared but not initialised?

  	cShortPairs = 0;
  	cIntersectPairs = 0;
  	cLongPairs = 0;
	cSingleReads = 0;
	cPairedReads = 0;
	cErrorReads = 0;
	cSkippedReads = 0;
}

// OK.
void BAM2blocks::readBamHeader() {
  char buffer[1000];
  std::string chrName;
  //std::vector<std::string> chr_names;
  bam_header bamhead;

  IN->read(bamhead.c, BAM_HEADER_BYTES);


  //IN->ignore(bamhead.l_text);
  char headertext[bamhead.l_text+1];
  IN->read(headertext, bamhead.l_text);
  samHeader = string(headertext, bamhead.l_text);

  stream_int32 i32;
  IN->read(i32.c ,4);
  uint n_chr = i32.i;

  for (uint i = 0; i < n_chr; i++) {
    IN->read(i32.c ,4);
    IN->read(buffer , i32.i);
    chrName = string(buffer, i32.i-1);
    chr_names.push_back(chrName);

    IN->read(i32.c ,4);
    chr_lens.push_back(i32.i);
  }
  
	for (auto & callback : callbacksChrMappingChange ) {
		callback(chr_names);
	}
  
}

void BAM2blocks::cigar2block(int32_t * cigar, uint16_t n_cigar_op, std::vector<int> &starts, std::vector<int> &lens, int &ret_genome_len) {
  bool inBlock = true;
  int relpos = 0;
  int curblock = 0;
  starts.resize(1);  // Is this expensive or not -- does this call destroy on further items, or is it a single op, adjusting the end? If expensive we can revert to earlier behaviour where we keep track of how many blocks, just overwriting relevant parts of the vector.
  lens.resize(1);
  starts[curblock] = 0;
  lens[curblock] = 0;

  for (; n_cigar_op > 0; n_cigar_op--) {
    if (inBlock) {
      switch (*cigar & 15) {
        case 0: case 2: case 7: case 8:
          // increment len of last block
          lens[curblock] += (*cigar >> 4);
          relpos += (*cigar >> 4);
          break;
        case 3:
          curblock++;
          relpos += (*cigar >> 4);
          // extend arrays. 
          starts.push_back(relpos);
          lens.push_back(0);
          inBlock = false;
          break;
      }
    }else{
      switch (*cigar & 15) {
        case 0: case 2: case 7: case 8:
          lens[curblock] = (*cigar >> 4);
          relpos += (*cigar >> 4);
          inBlock = true;
          break;
        case 3:
          // push start of next further out
          relpos += (*cigar >> 4);
          starts[curblock] = relpos;
          break;
        }
    }
    cigar++;
  }
  ret_genome_len = relpos;
//  *ret_blocks = curblock+1;  // Unnecessary if we are using vectors in the expected manner - changing their length as needed.
}


//OK - translated - doesn't call the callbacks yet though.
unsigned int BAM2blocks::processPair(bam_read_core * read1, bam_read_core * read2) {
  // R1 is to the left of R2 (or equal starts).
  int r1_genome_len;
  //int r1_blocks;
  int r2_genome_len;

  string debugstate;

  //int r2_blocks;
  //char dir;

  if (read1->flag & 0x40) {
    //this is first of pair.
    if (read1->flag & 0x10) {
      oBlocks.direction = 0;
    }else{
      oBlocks.direction = 1;
    }
  }else{
    if (read1->flag & 0x20) {
      oBlocks.direction = 0;
    }else{
      oBlocks.direction = 1;
    }
  }


  cigar2block(read1->cigar, read1->n_cigar_op, oBlocks.rStarts[0], oBlocks.rLens[0], r1_genome_len);
  cigar2block(read2->cigar, read2->n_cigar_op, oBlocks.rStarts[1], oBlocks.rLens[1], r2_genome_len);

  if (read1->pos + r1_genome_len < read2->pos) {
    cLongPairs++;
    //reads do not intersect
    oBlocks.readCount = 2;
    debugstate.append( "-Long-");
  }else if (read1->pos + r1_genome_len >= read2->pos + r2_genome_len){
    cShortPairs++;
    // Read 2 is a short read & read 1 fully contains it (or perhaps just a trimmed read with two exactly complementary reads remaining).
    oBlocks.readCount = 1;
    debugstate.append( "-Short-");    
  }else{
    debugstate.append( "-Intersect-");
    cIntersectPairs++;
    bool goodPair = true;
    oBlocks.readCount = 1;
    // We have two reads that intersect - construct just one complete fragment.

// Guaranteed assumptions:
//   Read 1 starts to the left of Read 2.
//   Read 2 end extends beyond the end of Read 1 end.
    int r1pos = read1->pos;
    int r2pos = read2->pos;
    for (uint i = 0; i < oBlocks.rStarts[0].size(); i++) {
        if (r1pos + oBlocks.rStarts[0][i] + oBlocks.rLens[0][i] >= r2pos) {
          if (r1pos + oBlocks.rStarts[0][i] <= r2pos) {
            oBlocks.rLens[0][i] = r2pos - r1pos - oBlocks.rStarts[0][i] + oBlocks.rLens[1][0];
            //r1_blocks = i + r2_blocks;
            oBlocks.rStarts[0].resize(i + oBlocks.rStarts[1].size());
            oBlocks.rLens[0].resize(i + oBlocks.rStarts[1].size());
            // Maybe this can be optimised by using push_back below instead of running resize.
            for (uint j = 1; j < oBlocks.rStarts[1].size(); j++) {
              i++;
              oBlocks.rLens[0][i] = oBlocks.rLens[1][j];
              oBlocks.rStarts[0][i] = oBlocks.rStarts[1][j] + r2pos - r1pos;
            }
            r1_genome_len = r2pos - r1pos + r2_genome_len;
            break;
          }else{
            //cerr << "Fault with this synthetic read, outputting each of the overlapping reads as singles: " << read1->read_name << endl;
            // This error is not worth reporting. The current version of STAR outputs a good number of these, concordance would be nice, but it is better to get at least one read illustrating the splice junction.
            goodPair = false;
            oBlocks.readCount = 2;
          }
        }
    }

    if (!goodPair) {
	    oBlocks.readCount = 2;
    }
  }
	oBlocks.chr_id = read1->refID;
	oBlocks.readStart[0] = read1->pos;
	oBlocks.readEnd[0] = read1->pos + r1_genome_len;
	oBlocks.readName.resize(read1->l_read_name - 1);
	oBlocks.readName.replace(0, read1->l_read_name - 1, read1->read_name, read1->l_read_name - 1); // is this memory/speed efficient?

	uint totalBlockLen = 0;
	for (auto blockLen: oBlocks.rLens[0]) {
		totalBlockLen += blockLen;
	}
	if (oBlocks.readCount > 1) {
		oBlocks.readStart[1] = read2->pos;
		oBlocks.readEnd[1] = read2->pos + r2_genome_len;
		for (auto blockLen: oBlocks.rLens[1]) {
			totalBlockLen += blockLen;
		}
	}
	//DEBUG:
	oBlocks.readName.append(debugstate);
	oBlocks.readName.append(to_string(oBlocks.readCount));
// TODO - restructure -- we could instead do the manipulation from 2 reads-> 1 synthetic in a non-const callback.
//        not required until that future flexibility is needed if part of the framework is repurposed.
	for (auto & callback : callbacksProcessBlocks ) {
		callback(oBlocks);
	}
	return totalBlockLen;
}


unsigned int BAM2blocks::processSingle(bam_read_core * read1) {
  int r1_genome_len;

  string debugstate;

  if (read1->flag & 0x10) {
    oBlocks.direction = 0;
  }else{
    oBlocks.direction = 1;
  }

  cigar2block(read1->cigar, read1->n_cigar_op, oBlocks.rStarts[0], oBlocks.rLens[0], r1_genome_len);
  oBlocks.readCount = 1;

	oBlocks.chr_id = read1->refID;
	oBlocks.readStart[0] = read1->pos;
	oBlocks.readEnd[0] = read1->pos + r1_genome_len;
	oBlocks.readName.resize(read1->l_read_name - 1);
	oBlocks.readName.replace(0, read1->l_read_name - 1, read1->read_name, read1->l_read_name - 1); // is this memory/speed efficient?
	//DEBUG:
	oBlocks.readName.append(debugstate);
	oBlocks.readName.append(to_string(oBlocks.readCount));
	//cout << "process pair - callbacks" << endl;  
	for (auto & callback : callbacksProcessBlocks ) {
		callback(oBlocks);
	}
	uint totalBlockLen = 0;
	for (auto blockLen: oBlocks.rLens[0]) {
		totalBlockLen += blockLen;
	}
	return totalBlockLen;
}



int BAM2blocks::processAll() {

	unsigned long long totalNucleotides = 0;
	ulong j = 0;
	int idx = 0;
	int pair = 0;
	//int bytesread = 0;
  
	while(1) {
		j++;
		if (IN->eof() && !(IN->fail()) ) {
			cout << "Total reads processed: " << j-1 << endl;
			cout << "Total nucleotides: " << totalNucleotides << endl;
			cout << "Total singles processed: " << cSingleReads << endl;
			cout << "Total pairs processed: " << cShortPairs+cIntersectPairs+cLongPairs << endl;
			cout << "Short pairs: " << cShortPairs << endl;
			cout << "Intersect pairs: " << cIntersectPairs << endl;
			cout << "Long pairs: " << cLongPairs << endl;
			cout << "Skipped reads: " << cSkippedReads << endl;
			cout << "Error reads: " << cErrorReads << endl;
			return(0);   
		}
		IN->read(reads[idx].c, BAM_READ_CORE_BYTES);
		if (IN->fail()) {
			cerr << "Input error at line:" << j << endl;
			cerr << "Characters read on last read call:" << IN->gcount() << endl;
			cout << "ERR-Total reads processed: " << j-1 << endl;
			cout << "ERR-Total nucleotides: " << totalNucleotides << endl;
			cout << "ERR-Total singles processed: " << cSingleReads << endl;
			cout << "ERR-Total pairs processed: " << cShortPairs+cIntersectPairs+cLongPairs << endl;
			cout << "ERR-Short pairs: " << cShortPairs << endl;
			cout << "ERR-Intersect pairs: " << cIntersectPairs << endl;
			cout << "ERR-Long pairs: " << cLongPairs << endl;
			cout << "ERR-Skipped reads: " << cSkippedReads << endl;
			cout << "ERR-Error reads: " << cErrorReads << endl;
			return(1);
			//This is possibly also just about the end of the file (say an extra null byte).
			//IN->gcount() knows how many characters were actually read last time.
		}
		IN->read(reads[idx].read_name, reads[idx].l_read_name);
		IN->read(reads[idx].cigar_buffer, reads[idx].n_cigar_op*4);    
		IN->ignore(reads[idx].block_size - BAM_READ_CORE_BYTES + 4 - reads[idx].l_read_name - (reads[idx].n_cigar_op*4));

		if (reads[idx].flag & 0x904) {
			/* If is an unmapped / secondary / supplementary alignment -- discard/overwrite */
			cSkippedReads ++;
		}else if (! (reads[idx].flag & 0x1)) {
			/* If is a single read -- process it as a single -- then discard/overwrite */
			cSingleReads ++;
			totalNucleotides += processSingle(&reads[idx]);
		}else{
			/* If it is potentially a paired read, store it in our buffer, process the pair together when it is complete */
			pair++;

			if (pair >= 2) {
				//There are now 2 reads sitting in the buffer. Maybe they are a true pair.
				if (reads[0].refID == reads[1].refID
						&& reads[0].l_read_name == reads[1].l_read_name
						&& strncmp(reads[0].read_name, reads[1].read_name, reads[0].l_read_name)==0) {
					cPairedReads ++;
					if (reads[0].pos <= reads[1].pos) {
						//cout << "procesPair call1" << endl;        
						totalNucleotides += processPair(&reads[0], &reads[1]);
					}else{
						//cout << "procesPair call2" << endl;                
						totalNucleotides += processPair(&reads[1], &reads[0]);
					}
					pair = 0;
				}else{
					cErrorReads++;
					//Bad pair? TODO.
					//The most likely cause is the user has provided a sorted BAM file, rather than one with reads of the pair next to each other. The ErrorReads counter implies this, and can be picked up by the cleanup QC script.
					// Process as a single?
					// Increment error counter?
					// -- for single end read we are expecting to process as a single.
					// -- need to eat the oldest & process it, allow the next to potentially pair.
				}
			}    
			//The index is only moved in case of pair storage.
			idx++;
			if (idx > 1) idx = 0;
		}
	}
	return(0);
}



void BAM2blocks::openFile(istream * _IN) {
	IN = _IN;
  	readBamHeader(); // readBamHeader needs to call the ChrMappingChange callbacks.
}

void BAM2blocks::registerCallbackChrMappingChange( std::function<void(const std::vector<string> &)> callback ) {
	callbacksChrMappingChange.push_back(callback);
}

void BAM2blocks::registerCallbackProcessBlocks( std::function<void(const FragmentBlocks &)> callback ) {	
	callbacksProcessBlocks.push_back(callback);
}
