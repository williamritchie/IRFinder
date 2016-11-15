#include "ReadBlockProcessor_CoverageBlocks.h"
#include "includedefine.h"

void CoverageBlocks::loadRef(std::istream &IN) {
	std::string myLine;
	std::string myField;
	myLine.reserve(1000);
	myField.reserve(100);

	//string s_chr;
	string sLengths;
	string sOffsets;
	//uint i_start;
	//uint i_end;
	uint i_block_start;
	uint i_block_end;
	uint i_segments;
	//string s_keydata;
	string s_dir;
	//s_keydata.reserve(400);
	BEDrecord BEDrec;

	std::map<string, std::vector<std::pair<uint, uint>> > temp_segments;

	while (!IN.eof()) {
		getline(IN, myLine);
		std::istringstream lineStream;
		std::istringstream lensStream;
		std::istringstream offsetsStream;
		lineStream.str(myLine);

		getline(lineStream, BEDrec.chrName, '\t');		//Chr
		getline(lineStream, myField, '\t');		//BED Start
		BEDrec.start = stoul(myField);
		getline(lineStream, myField, '\t');		//BED End
		BEDrec.end = stoul(myField);
		getline(lineStream, BEDrec.name, '\t');	//BED entry name.
				/* This is separated by '/' -- do we actually need to split any of this for our calculations? Maybe - final output should probably be BED like, showing start,end pos of the original BED blocks/intron rather than the excl trimmed form. */
		lineStream.ignore( numeric_limits<streamsize>::max(), '\t' ); //Score /* Throw away the next field */
		getline(lineStream, s_dir, '\t');		//Block direction +/-/.
		BEDrec.direction = (s_dir == "+");
		lineStream.ignore( numeric_limits<streamsize>::max(), '\t' );	//Thick BED Start
		lineStream.ignore( numeric_limits<streamsize>::max(), '\t' );	//Thick BED End
		lineStream.ignore( numeric_limits<streamsize>::max(), '\t' );	//BED Colour
		getline(lineStream, myField, '\t');		//BED block count
		i_segments = stoul(myField);
		getline(lineStream, sLengths, '\t');	//Comma separated lengths.  
		lensStream.str(sLengths);

		if (IN.eof()) {
			//ie: we don't have a complete line here -- maybe because the last line was just a "\n".
			break;
			// TODO -- better use eof() / fail().
		}
		getline(lineStream, sOffsets, '\t');	//Comma separated offsets.
		offsetsStream.str(sOffsets);

		// The Blocked BED records - simply in a vector per Chromosome - sorted as input. (this may need a storage class for these records, optionally that storage class could also do the processing at output stage)
		// It would be more efficient to store on reading directly into the final location in the BEDrecords vector ..
		// Optimisation unnecessary, this only runs at startup.

		BEDrec.blocks.clear();
		//Process the discovered blocks to work out the regions we need to measure & store coverage depth.
		for (int i = i_segments; i > 0; i--) {
			getline(offsetsStream, myField, ',');
			i_block_start = BEDrec.start + stoul(myField);
			getline(lensStream, myField, ',');
			i_block_end = i_block_start + stoul(myField);
			BEDrec.blocks.push_back(std::make_pair( i_block_start, i_block_end ));
			temp_segments[BEDrec.chrName].push_back(std::make_pair( i_block_start, i_block_end ) ); 
		}		
		BEDrecords.push_back(BEDrec);
	}
	// Read from file complete.


	// We have a map of Chr->Vectors. Each vector is a pair of start/stop BED coords.
	//  We want to sort and merge these & produce as result, a set of non-overlapping BED
	//  blocks, each of a maximum length, which we can use to create our data structures
	//  for the computation of coverage.

	// Construct chrName_CoverageBlocks(Pos/Neg). The ChrName->Vector map of non-overlapping "CoverageBlock" objects. (this is where depth gets recorded)
	for (std::map<string, std::vector<std::pair<uint, uint>>>::iterator it_chr=temp_segments.begin(); it_chr!=temp_segments.end(); it_chr++) {
		// foreach chromosome.
		std::sort(it_chr->second.begin(), it_chr->second.end());

		uint next_start = 0;
		uint next_end = 0;
		uint max_gap = 50;
		bool first_record = true;
		for (std::vector<std::pair<uint, uint>>::iterator it_blocks=it_chr->second.begin(); it_blocks!=it_chr->second.end(); it_blocks++) {
			if (it_blocks->first > next_end + max_gap) {
				//We are passed the last block, output the old.
				if (first_record) {
					first_record = false;
				}else{
					while (next_start + coverage_block_max_length < next_end) {
						chrName_CoverageBlocks[it_chr->first].push_back(CoverageBlock(next_start, next_start + coverage_block_max_length));
						next_start += coverage_block_max_length;
					}
					chrName_CoverageBlocks[it_chr->first].push_back(CoverageBlock(next_start, next_end));
				}
				next_start = it_blocks->first;
				next_end = it_blocks->second;
			}else if (it_blocks->second > next_end) {
        		next_end = it_blocks->second;
			}
		}
		if (!first_record) {
			while (next_start + coverage_block_max_length < next_end) {
				chrName_CoverageBlocks[it_chr->first].push_back(CoverageBlock(next_start, next_start + coverage_block_max_length));
				next_start += coverage_block_max_length;
			}
			chrName_CoverageBlocks[it_chr->first].push_back(CoverageBlock(next_start, next_end));
		}
		chrName_CoverageBlocks[it_chr->first].shrink_to_fit();
	}
}

void CoverageBlocks::ChrMapUpdate(const std::vector<string> &chrmap) {
	chrID_CoverageBlocks.resize(0);
	for (uint i = 0; i < chrmap.size(); i++) {
		chrID_CoverageBlocks.push_back( &chrName_CoverageBlocks[chrmap.at(i)] );
	}
}


void CoverageBlocks::ProcessBlocks(const FragmentBlocks &blocks) {
	std::vector<CoverageBlock>::iterator it_coverblock;
	uint start;
	uint end;

	//Walk each read within the fragment (1 or 2).
	for (int index = 0; index < blocks.readCount; index ++) {
		//Walk each block within each read.
		for (uint j = 0; j < blocks.rLens[index].size(); j++) {
			// Have a block.
			start = blocks.readStart[index] + blocks.rStarts[index][j];
			end = start + blocks.rLens[index][j];
			// need abs coords & blocks.chr_id.
			// Then do appropriate upper_bound/lower_bound search - find our first block, make the call.

			it_coverblock = std::upper_bound(
				(*chrID_CoverageBlocks.at(blocks.chr_id)).begin(),
				(*chrID_CoverageBlocks.at(blocks.chr_id)).end(),
				start
			);
			while (it_coverblock != (*chrID_CoverageBlocks.at(blocks.chr_id)).end() && it_coverblock->posIsAfterStart(end)) {
				it_coverblock->RecordCover(start, end, blocks.direction);

				it_coverblock++;
			}
		}
	}
}

void CoverageBlocks::fillHist(std::map<uint,uint> &hist, const std::string &chrName, const std::vector<std::pair<uint,uint>> &blocks) const {
	std::vector<CoverageBlock>::const_iterator it_CB;
	for (std::vector<std::pair<uint,uint>>::const_iterator it_blocks=blocks.begin(); it_blocks!=blocks.end(); it_blocks++) {
		it_CB = upper_bound(
			chrName_CoverageBlocks.at(chrName).begin(),
			chrName_CoverageBlocks.at(chrName).end(),
			it_blocks->first
			);
		while (it_CB != chrName_CoverageBlocks.at(chrName).end() && it_CB->posIsAfterStart(it_blocks->second)) {
			it_CB->updateCoverageHist(hist, it_blocks->first, it_blocks->second);  //for non-dir.
			it_CB++;
		}
	}
}

void CoverageBlocks::fillHist(std::map<uint,uint> &hist, const std::string &chrName, const std::vector<std::pair<uint,uint>> &blocks, bool direction) const {
	std::vector<CoverageBlock>::const_iterator it_CB;
	for (std::vector<std::pair<uint,uint>>::const_iterator it_blocks=blocks.begin(); it_blocks!=blocks.end(); it_blocks++) {
		it_CB = upper_bound(
			chrName_CoverageBlocks.at(chrName).begin(),
			chrName_CoverageBlocks.at(chrName).end(),
			it_blocks->first
			);
		while (it_CB != chrName_CoverageBlocks.at(chrName).end() && it_CB->posIsAfterStart(it_blocks->second)) {
			it_CB->updateCoverageHist(hist, it_blocks->first, it_blocks->second, direction);  //directional.
			it_CB++;
		}
	}
}


double CoverageBlocks::meanFromHist(const std::map<uint,uint> &hist) const {
	unsigned long long total = 0;
	uint count = 0;
	
	for (auto h : hist) {
		total += h.first * h.second;
		count += h.second;
	}
	return (total/(double)count);
}

double CoverageBlocks::coverageFromHist(const std::map<uint,uint> &hist) const {
	if (hist.find(0) == hist.end()) {
		return 1.0; //No bases are at zero cover.
	}
	uint count = 0;
	for (auto h : hist) {
		count += h.second;
	}
	return ((count - hist.at(0))/(double)count);
}


double CoverageBlocks::percentileFromHist(const std::map<uint,uint> &hist, uint percentile) const {
	uint size = 0;
	for (auto h : hist) {
		size += h.second;
	}
	double percentile_frac = (size + 1)*(double)percentile/100;
	uint percentile_index = percentile_frac;  //round down
	percentile_frac = percentile_frac - percentile_index;

	uint count = 0;
	for (auto h = hist.begin(); h != hist.end(); h++) {
		count += h->second;
		if (count >= percentile_index) {
			if (count > percentile_index || percentile_frac == 0) {
				return h->first;
			}else{
				double ret = h->first - (percentile_frac * h->first);
				h++;
				ret += (percentile_frac * h->first);
				return ret;
			}
		}
	}
	return NAN;
}

double CoverageBlocks::trimmedMeanFromHist(const std::map<uint,uint> &hist, uint centerPercent) const {
	uint size = 0;
	for (auto h : hist) {
		size += h.second;
	}
	uint skip = size * ((100 - centerPercent)/2) / 100;
	
	unsigned long long total = 0;
	uint count = 0;
	
	for (auto h : hist) {
		if (count + h.second > size - skip) {
			// This bar will enter the max skip section.
			if (count > skip) {
				//already inside target range
				total += h.first * (size - skip - count);
			}else{
				//yet to enter target range
				return h.first; //(all relevant numbers are the same for this mean)
			}
			break;
		}
		if (count > skip) {
			// Start and stop are fully inside the counted section.
			total += h.first * h.second;
		}else if (count + h.second > skip) {
			// We leave the min skip section and use some of the size of this hist bar.
			total += h.first * (count + h.second - skip);
		}
		count += h.second;
	}
	return ((double)total/(size-2*skip));
}



int CoverageBlocks::WriteOutput(std::ostream *os) const {

// This output function will be generic -- outputting Chr/Start/Stop/Name/Dir/ Score - Mean50 (that bit probably cmd line customisable).
// The output we need will be in the extended class.

	for (std::vector<BEDrecord>::const_iterator it_BED=BEDrecords.begin(); it_BED!=BEDrecords.end(); it_BED++) {
		uint len=0;
		for (std::vector<std::pair<uint,uint>>::const_iterator it_blocks=it_BED->blocks.begin(); it_blocks!= it_BED->blocks.end(); it_blocks++) {
			len += (it_blocks->second - it_blocks->first);
		}
		std::map<uint,uint> hist;
		fillHist(hist, it_BED->chrName, it_BED->blocks);

		uint histPositions = 0;
		for (auto h : hist) {
			histPositions += h.second;
			//DEBUGGING
			*os << h.first << "\t" << h.second << "\n";
		}

		//*os << "\n";
		*os << it_BED->chrName << "\t" << it_BED->start << "\t" << it_BED->end << "\t" << (it_BED->end - it_BED->start) << "\t" << histPositions << "\t" << hist.size() << "\t" << trimmedMeanFromHist(hist, 50)  << "\t" << trimmedMeanFromHist(hist, 20) << "\t" << coverageFromHist(hist) << "\t" << meanFromHist(hist) << "\t" << it_BED->direction << "\t" << it_BED->name << "\n";
		*os << percentileFromHist(hist, 25) << "\t" << percentileFromHist(hist, 50) << "\t" << percentileFromHist(hist, 75) << "\t" << "\n";
	}
	
	return 0;
}


int CoverageBlocksIRFinder::WriteOutput(std::ostream *os, const JunctionCount &JC, const SpansPoint &SP, int directionality) const {
	// Custom output function - related to the IRFinder needs
	uint recordNumber = 0;
	
	for (auto BEDrec : BEDrecords) {
		recordNumber++;
		// if name indicates it is a Dir/Non-dir record of interest - output it.
		// We need to separate dir&non-dir by name. (.startswith)
		if ((directionality != 0 && (0 == BEDrec.name.compare(0, 4, "dir/"))) || (directionality == 0 && (0 == BEDrec.name.compare(0, 3, "nd/")))) {
			try {
				uint intronStart;
				uint intronEnd;
				uint exclBases;
				double intronTrimmedMean;
				double coverage;
				bool measureDir;
				uint JCleft;
				uint JCright;
				uint JCexact;
				uint SPleft;
				uint SPright;


				std::string s_buffer;
				std::string s_name;
				std::string s_ID;
				std::string s_clean;

				std::istringstream lineStream;
				lineStream.str(BEDrec.name);
				lineStream.ignore( numeric_limits<streamsize>::max(), '/' );
				getline(lineStream, s_name, '/');
				getline(lineStream, s_ID, '/');
				lineStream.ignore( numeric_limits<streamsize>::max(), '/' );
				lineStream.ignore( numeric_limits<streamsize>::max(), '/' );
				getline(lineStream, s_buffer, '/');
				intronStart = stol(s_buffer);
				getline(lineStream, s_buffer, '/');
				intronEnd = stol(s_buffer);
				lineStream.ignore( numeric_limits<streamsize>::max(), '/' );
				getline(lineStream, s_buffer, '/');
				exclBases = stol(s_buffer);
				getline(lineStream, s_clean, '/');

	//1       860574  861258  nd/SAMD11/ENSG00000187634/+/2/860569/861301/732/121/anti-over   0       +       860574  861258  255,0,0 2       538,73  0,611
	//1       860574  861296  dir/SAMD11/ENSG00000187634/+/2/860569/861301/732/83/clean       0       +       860574  861296  255,0,0 2       538,111 0,611


				//eg: PHF13/ENSG00000116273/+/3/6676918/6679862/2944/10/clean
				*os << BEDrec.chrName << "\t" << intronStart << "\t" << intronEnd << "\t" << s_name << "/" << s_ID << "/" << s_clean << "\t0\t" << ((BEDrec.direction) ?  "+" : "-" ) << "\t";

				measureDir = BEDrec.direction;
				if (directionality == -1) {
					measureDir = !BEDrec.direction;
				}

				std::map<uint,uint> hist;
				if (directionality == 0) {
					fillHist(hist, BEDrec.chrName, BEDrec.blocks);
				}else{
					fillHist(hist, BEDrec.chrName, BEDrec.blocks, measureDir);
				}
				intronTrimmedMean = trimmedMeanFromHist(hist, 40);
				coverage = coverageFromHist(hist);
				*os << exclBases << "\t"
					<< coverage << "\t"
					<< intronTrimmedMean << "\t"
					<< percentileFromHist(hist, 25) << "\t"
					<< percentileFromHist(hist, 50) << "\t"
					<< percentileFromHist(hist, 75) << "\t";

				if (directionality != 0) {
					SPleft = SP.lookup(BEDrec.chrName, intronStart, measureDir);
					SPright = SP.lookup(BEDrec.chrName, intronEnd, measureDir);
					*os << SPleft << "\t"
						<< SPright << "\t";

					hist.clear();
					fillHist(hist, BEDrec.chrName, {{intronStart + 5, intronStart + 55}}, measureDir);
					*os << trimmedMeanFromHist(hist, 40) << "\t";
					hist.clear();
					fillHist(hist, BEDrec.chrName, {{intronEnd - 55, intronEnd - 5}}, measureDir);
					*os << trimmedMeanFromHist(hist, 40) << "\t";
					JCleft = JC.lookupLeft(BEDrec.chrName, intronStart, measureDir);
					JCright = JC.lookupRight(BEDrec.chrName, intronEnd, measureDir);
					JCexact = JC.lookup(BEDrec.chrName, intronStart, intronEnd, measureDir);
					*os << JCleft << "\t"
						<< JCright << "\t"
						<< JCexact << "\t";
				}else{
					SPleft = SP.lookup(BEDrec.chrName, intronStart);
					SPright = SP.lookup(BEDrec.chrName, intronEnd);
					*os << SPleft << "\t"
						<< SPright << "\t";			

					hist.clear();
					fillHist(hist, BEDrec.chrName, {{intronStart + 5, intronStart + 55}});
					*os << trimmedMeanFromHist(hist, 40) << "\t";
					hist.clear();
					fillHist(hist, BEDrec.chrName, {{intronEnd - 55, intronEnd - 5}});
					*os << trimmedMeanFromHist(hist, 40) << "\t";
					JCleft = JC.lookupLeft(BEDrec.chrName, intronStart);
					JCright = JC.lookupRight(BEDrec.chrName, intronEnd);
					JCexact = JC.lookup(BEDrec.chrName, intronStart, intronEnd);
					*os << JCleft << "\t"
						<< JCright << "\t"
						<< JCexact << "\t";
				}
				if (intronTrimmedMean == 0 && JCleft == 0 && JCright == 0) {
					*os << "0" << "\t";
				}else if (intronTrimmedMean < 1) {
					*os << ( coverage / (coverage + max(JCleft, JCright)) ) << "\t";
				}else{
					*os << ( intronTrimmedMean /(intronTrimmedMean + max(JCleft, JCright)) ) << "\t";
				}
				
				
				// Final column -- don't try to be tri-state. Just say if it is "not ok".
				// Not ok due to:
				//	- insufficient spliced depth
				//  - insufficient exact spliced compared to in-exact spliced depth
				//  - too much variation between depths & crossings.  ... hmm, but at low depth, high probability of this failing.
				
				// Can only make a strong exclude call on spliced depth. Describe on the tool website ways to make a call for IR def true / IR def false.
//				if (JCexact < 10 || JCexact*1.33333333 < max(JCleft, JCright) ) {
//					*os << "-" << "\n";
//				}else{
//					*os << "ok" << "\n";
//				}

				if (JCexact + intronTrimmedMean < 10) {
					*os << "LowCover" << "\n";
				}else if (JCexact < 4) {
					*os << "LowSplicing" << "\n";
				}else if (JCexact*1.33333333 < max(JCleft, JCright) ) {
					*os << "MinorIsoform" << "\n";
				// TODO: check, logic below. Crossing should differ by more than 2 & more than 50% before a fault is called.
				}else if (  (max(SPleft, SPright) > intronTrimmedMean+2 && max(SPleft, SPright) > intronTrimmedMean*1.5 )
						|| (min(SPleft, SPright)+2 < intronTrimmedMean && min(SPleft, SPright)*1.5 < intronTrimmedMean ) ){
					*os << "NonUniformIntronCover" << "\n";
				}else{
					*os << "-" << "\n";
				}
				
			}catch (const std::out_of_range& e) {
				std::cerr << "Format error in name attribute - column 4 - of CoverageBlocks reference file. Record/line number: " << recordNumber << "\n";
			}catch (const std::invalid_argument& e) {
				std::cerr << "Format error in name attribute - column 4 - of CoverageBlocks reference file. Record/line number: " << recordNumber << "\n";			
			}
		}
	}
	return 0;
}
