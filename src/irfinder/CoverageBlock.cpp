#include "CoverageBlock.h"
#include "includedefine.h"
// using namespace std;


CoverageBlock::CoverageBlock(uint start, uint end) {
	blockStart = start;
	blockEnd = end;
	
    firstDepth[0] = 0;
    firstDepth[1] = 0;

    blockExtents = NULL;
    blockExtentsL = NULL;
}

//direction -- 0=False/Neg, 1=True/Pos.
void CoverageBlock::RecordCover(uint readStart, uint readEnd, bool dir) {
	if (readStart <= blockStart && readEnd > blockStart) {
		firstDepth[dir]++;
	}else if (readStart < blockEnd) {
		// Need to increment the starts vector.
		uint inc_index = readStart - blockStart - 1;
		if (blockExtentsL) { //already an int vector
			blockExtentsL->at(inc_index).start[dir]++;
		}else if (!blockExtents) { //don't have a char vector either - create first.
			blockExtents = new std::vector<start_stops>(vectorLen());
			blockExtents->at(inc_index).start[dir]++;
		}else{
			if (blockExtents->at(inc_index).start[dir] == 254) {
				blockExtentsL = new std::vector<start_stopsL>(blockExtents->begin(), blockExtents->end());
				delete blockExtents;
				blockExtents = NULL;
				blockExtentsL->at(inc_index).start[dir]++;
			}else{
				blockExtents->at(inc_index).start[dir]++;
			}
		}
	}else{
		return;
	}

	if (readEnd >= blockEnd) {
		return;
	}else{
		// Need to increment the ends vector.
		uint inc_index = readEnd - blockStart - 1;

		if (blockExtentsL) { //already an int vector
			blockExtentsL->at(inc_index).end[dir]++;
		}else if (!blockExtents) { //don't have a char vector either - create first.
			blockExtents = new std::vector<start_stops>(vectorLen());
			blockExtents->at(inc_index).end[dir]++;
		}else{
			if (blockExtents->at(inc_index).end[dir] == 254) {
				blockExtentsL = new std::vector<start_stopsL>(blockExtents->begin(), blockExtents->end());
				delete blockExtents;
				blockExtents = NULL;
				blockExtentsL->at(inc_index).end[dir]++;
			}else{
				blockExtents->at(inc_index).end[dir]++;
			}
		}
	}	
	// Can Throw: Out of range exception.
}


void CoverageBlock::updateCoverageHist(std::map<uint,uint> &hist, uint start, uint end) const {
	if (!blockExtentsL && !blockExtents) {
		// how many bases in this block?
		hist[firstDepth[0]+firstDepth[1]] += min(blockEnd, end) - max(blockStart,start);
	}else{
		// There are read starts and ends -- need to walk the positions from the start of this block
		//  even if not in the region of interest.

		//special handling for the first base -- the one before the vector starts.
		uint depth = firstDepth[0]+firstDepth[1];
		if (start <= blockStart) {
			// use the first depth, before commencing in the vector.
			hist[depth] ++;
		}

		uint startindex = max(blockStart+1, start) - blockStart - 1;
		uint endindex = min(blockEnd, end) - blockStart - 1;
		if (blockExtents) {
			for (uint i=0; i<endindex; i++) {
				depth += - (*blockExtents)[i].end[0] - (*blockExtents)[i].end[1] + (*blockExtents)[i].start[0] + (*blockExtents)[i].start[1];
				if (i>=startindex) {
					hist[depth] ++;
				}
			}
		}else{
			for (uint i=0; i<endindex; i++) {
				depth += - (*blockExtentsL)[i].end[0] - (*blockExtentsL)[i].end[1] + (*blockExtentsL)[i].start[0] + (*blockExtentsL)[i].start[1];
				if (i>=startindex) {
					hist[depth] ++;
				}
			}
		}
		//  When in the region of interest, update the hist each step.
	}
}

void CoverageBlock::updateCoverageHist(std::map<uint,uint> &hist, uint start, uint end, bool dir) const {
	if (!blockExtentsL && !blockExtents) {
		// how many bases in this block?
		hist[firstDepth[dir]] += min(blockEnd, end) - max(blockStart,start);
	}else{
		//special handling for the first base -- the one before the vector starts.
		uint depth = firstDepth[dir];
		if (start <= blockStart) {
			// use the first depth, before commencing in the vector.
			hist[depth] ++;
		}

		uint startindex = max(blockStart+1, start) - blockStart - 1;
		uint endindex = min(blockEnd, end) - blockStart - 1;
		if (blockExtents) {
			for (uint i=0; i<endindex; i++) {
				depth += - (*blockExtents)[i].end[dir] + (*blockExtents)[i].start[dir];
				if (i>=startindex) {
					hist[depth] ++;
				}
			}
		}else{
			for (uint i=0; i<endindex; i++) {
				depth += - (*blockExtentsL)[i].end[dir] + (*blockExtentsL)[i].start[dir];
				if (i>=startindex) {
					hist[depth] ++;
				}
			}
		}
	}
}


void CoverageBlock::print( std::ostream& os) const {
	os << "Start\t" << firstDepth[0] << "\t" << firstDepth[1] << "\t" << (bool)blockExtents << (bool)blockExtentsL;
/*	os << " Start Depth\t" << firstDepth << "\tStarts" << (bool)readStartsL << (bool)readStarts << "  Ends" << (bool)readEndsL << (bool)readEnds << "  Pos: " << blockStart << "\t";
	if (readStartsL) {
		os << (uint)*max_element(readStartsL->begin(), readStartsL->end());
	}else if(readStarts) {
		os << (uint)*max_element(readStarts->begin(), readStarts->end());
	}else{
		os << "-";
	}
	os << "\t";
	if (readEndsL) {
		os << (uint)*max_element(readEndsL->begin(), readEndsL->end());
	}else if(readEnds) {
		os << (uint)*max_element(readEnds->begin(), readEnds->end());
	}else{
		os << "-";
	}*/
}

std::ostream& operator<<( std::ostream& os, const CoverageBlock& cb) {
	cb.print( os );
	return os;
}