#include "ReadBlockProcessor.h"
#include "includedefine.h"
// using namespace std;



void BED12Output::SetOutputStream(ostream *os) {
	out = os;
}

void BED12Output::ProcessBlocks(const FragmentBlocks &blocks) {
	// First, check how many reads are needed.
	for (int index = 0; index < blocks.readCount; index ++) {
		*out << chr_names.at(blocks.chr_id) << "\t"
			<< blocks.readStart[index] << "\t"
			<< blocks.readEnd[index] << "\t"
			<< blocks.readName << "\t0\t"
			<< ((blocks.direction) ?  "+" : "-" ) << "\t"
			<< blocks.readStart[index] << "\t"
			<< blocks.readEnd[index] << "\t0,0,0\t"
			<< blocks.rStarts[index].size() << "\t";
		copy_n(blocks.rLens[index].begin(), blocks.rStarts[index].size(), infix_ostream_iterator<int>(*out, ","));
		*out << "\t";
		copy_n(blocks.rStarts[index].begin(), blocks.rStarts[index].size(), infix_ostream_iterator<int>(*out, ",")); 
		*out << "\n";
	}
}

void BED12Output::ChrMapUpdate(const std::vector<string> &chrmap) {
	chr_names = chrmap; // copy.
}



//chrName_junc_count holds the data structure -- ChrName(string) -> Junc Start/End -> count.
//chrID_junc_count holds the ChrID -> ...
//  where the ChrID is the ChrID relating to the appropriate ChrName, as understood by the currently processed BAM file.
void JunctionCount::ChrMapUpdate(const std::vector<std::string> &chrmap) {
	chrID_junc_count.resize(0);
	chrID_juncLeft_count.resize(0);
	chrID_juncRight_count.resize(0);
	// Below could be done with an iterator - i is not used except for element access of the single collection.
	for (uint i = 0; i < chrmap.size(); i++) {
		chrID_junc_count.push_back( &chrName_junc_count[chrmap.at(i)] );
		chrID_juncLeft_count.push_back( &chrName_juncLeft_count[chrmap.at(i)] );
		chrID_juncRight_count.push_back( &chrName_juncRight_count[chrmap.at(i)] );
	}
	
	
//		std::map<string, std::map<uint,uint[2]>> chrName_juncLeft_count;
//		std::vector<std::map<uint,uint[2]>*> chrID_juncLeft_count;

//		std::map<string, std::map<uint,uint[2]>> chrName_juncRight_count;
//		std::vector<std::map<uint,uint[2]>*> chrID_juncRight_count;

}

void JunctionCount::loadRef(std::istream &IN) {
	// ChrName, Start, End, direction(+/-/.).
	std::string myLine;
	std::string myField;
	myLine.reserve(1000);
	myField.reserve(100);
	uint start;
	uint end;
	string s_chr;
	s_chr.reserve(30);
	string direction;

	while(!IN.eof() && !IN.fail()) {
		getline(IN, myLine, '\n');

		if (IN.eof() || IN.fail()) {
			if (myLine.length() == 0) {
				// This line is empty - just a blank line at the end of the file.
				// Checking at this stage allows correct handling of files both with and without a trailing \n after the last record.
				break;
			}else{
				// Error line in input, ignore.
				break;
			}
		}

		std::istringstream lineStream;
		lineStream.str(myLine);
		
		getline(lineStream, s_chr, '\t');
		getline(lineStream, myField, '\t');
		start = stol(myField);
		getline(lineStream, myField, '\t');
		end = stol(myField);
		getline(lineStream, direction, '\t');
		
		if (direction == "-")  {
			chrName_junc_count[s_chr][make_pair(start,end)][2] = 1;
		}else if (direction == "+") {
			chrName_junc_count[s_chr][make_pair(start,end)][2] = 2;
		}
	}
}

void JunctionCount::ProcessBlocks(const FragmentBlocks &blocks) {
	for (int index = 0; index < blocks.readCount; index ++) {
		//Walk each *pair* of blocks. ie: ignore a read that is just a single block.
		for (uint j = 1; j < blocks.rLens[index].size(); j++) {
			if ((blocks.rLens[index][j-1] >= 5) && (blocks.rLens[index][j] >= 5)) {
				(*chrID_junc_count[blocks.chr_id])[
					make_pair(
						blocks.readStart[index] + blocks.rStarts[index][j-1] + blocks.rLens[index][j-1],
						blocks.readStart[index] + blocks.rStarts[index][j])
					][blocks.direction]++;
				(*chrID_juncLeft_count[blocks.chr_id])[
						blocks.readStart[index] + blocks.rStarts[index][j-1] + blocks.rLens[index][j-1]
					][blocks.direction]++;
				(*chrID_juncRight_count[blocks.chr_id])[
						blocks.readStart[index] + blocks.rStarts[index][j]
					][blocks.direction]++;
			}
		}
	}
}

int JunctionCount::WriteOutput(std::ostream *os) const {
	for (auto itChr=chrName_junc_count.begin(); itChr!=chrName_junc_count.end(); itChr++) {
		string chr = itChr->first;
		for (auto itJuncs=itChr->second.begin(); itJuncs!=itChr->second.end(); ++itJuncs) {
			*os << chr << "\t" << itJuncs->first.first << "\t" << itJuncs->first.second
				<< "\t" << ( (itJuncs->second)[2] == 1 ? "-" : (itJuncs->second)[2] == 2 ? "+" : "." )
				<< "\t" << ((itJuncs->second)[1] + (itJuncs->second)[0])
				<< "\t" << (itJuncs->second)[1]
				<< "\t" << (itJuncs->second)[0] << "\n";
		}
	}
	return 0;
}

int JunctionCount::Directional() const {
	uint dir_same = 0;
	uint dir_diff = 0;

	uint dir_evidence = 0;
	uint nondir_evidence = 0;
	uint dir_evidence_known = 0;
	uint nondir_evidence_known = 0;

	for (auto itChr=chrName_junc_count.begin(); itChr!=chrName_junc_count.end(); itChr++) {
		for (auto itJuncs=itChr->second.begin(); itJuncs!=itChr->second.end(); ++itJuncs) {
			if (((itJuncs->second)[1] + (itJuncs->second)[0]) > 8) {
				if ((itJuncs->second)[0] > (itJuncs->second)[1] * 4) {
					dir_evidence++;
					if ((itJuncs->second)[2] == 1) { //Ref is "-"
						dir_same++;
					}else if ((itJuncs->second)[2] == 2) {
						dir_diff++;
					}
				}else if ((itJuncs->second)[1] > (itJuncs->second)[0] * 4) {
					dir_evidence++;
					if ((itJuncs->second)[2] == 2) { //Ref is "+"
						dir_same++;
					}else if ((itJuncs->second)[2] == 1) {
						dir_diff++;
					}				
				}else{
					nondir_evidence++;
					if ((itJuncs->second)[2] > 0) {
						nondir_evidence_known++;
					}
				}
			}
		}
	}
	dir_evidence_known = dir_same + dir_diff;
	cout << "Directionality: Dir evidence:\t" << dir_evidence << "\n";
	cout << "Directionality: Nondir evidence:\t" << nondir_evidence << "\n";
	cout << "Directionality: Dir evidence known junctions:\t" << dir_evidence_known << "\n";
	cout << "Directionality: Nondir evidence known junctions:\t" << nondir_evidence_known << "\n";
	cout << "Directionality: Dir matches ref:\t" << dir_same << "\n";
	cout << "Directionality: Dir opposed to ref:\t" << dir_diff << "\n";
	cout << "Directionality: Dir score all (0-10000):\t" << ((long long)dir_evidence * 10000 / (dir_evidence + nondir_evidence + 1)) << "\n"; //+1 to prevent divide by zero errors.
	long dir_score_known = ((long long)dir_evidence_known * 10000 / (dir_evidence_known + nondir_evidence_known + 1));
	cout << "Directionality: Dir score known junctions (0-10000):\t" << dir_score_known << "\n";

	if ((dir_same > dir_diff * 100) && (dir_score_known >= 9000)) {
		return 1;
	}else if ((dir_diff > dir_same * 100) && (dir_score_known >= 9000)) {
		return -1;
	}else{
		return 0;
	}
}

uint JunctionCount::lookup(std::string ChrName, uint left, uint right, bool direction) const {
	try {
		return chrName_junc_count.at(ChrName).at(make_pair(left, right))[direction];
	}catch (const std::out_of_range& e) {
	}
	return 0;
}
uint JunctionCount::lookup(std::string ChrName, uint left, uint right) const {
	try {
		return chrName_junc_count.at(ChrName).at(make_pair(left, right))[0] + chrName_junc_count.at(ChrName).at(make_pair(left, right))[1];
	}catch (const std::out_of_range& e) {
	}
	return 0;
}
uint JunctionCount::lookupLeft(std::string ChrName, uint left, bool direction) const {
	try {
		return chrName_juncLeft_count.at(ChrName).at(left)[direction];
	}catch (const std::out_of_range& e) {
	}
	return 0;
}
uint JunctionCount::lookupLeft(std::string ChrName, uint left) const {
	try {
		return chrName_juncLeft_count.at(ChrName).at(left)[0] + chrName_juncLeft_count.at(ChrName).at(left)[1];
	}catch (const std::out_of_range& e) {
	}
	return 0;
}
uint JunctionCount::lookupRight(std::string ChrName, uint right, bool direction) const {
	try {
		return chrName_juncRight_count.at(ChrName).at(right)[direction];
	}catch (const std::out_of_range& e) {
	}
	return 0;
}
uint JunctionCount::lookupRight(std::string ChrName, uint right) const {
	try {
		return chrName_juncRight_count.at(ChrName).at(right)[0] + chrName_juncRight_count.at(ChrName).at(right)[1];
	}catch (const std::out_of_range& e) {
	}
	return 0;
}





int SpansPoint::WriteOutput(std::ostream *os) const {
	for (auto itChrPos=chrName_pos.begin(); itChrPos!=chrName_pos.end(); itChrPos++) {
		string chr = itChrPos->first;

		auto itCountPos=chrName_count[1].at(chr).begin();
		auto itCountNeg=chrName_count[0].at(chr).begin();

		
		for (auto itPosition=itChrPos->second.begin(); itPosition!=itChrPos->second.end(); ++itPosition) {
			*os << chr << "\t" << *itPosition << "\t" << (*itCountPos + *itCountNeg) << "\t" << *itCountPos << "\t" << *itCountNeg << "\n";
			++itCountPos;
			++itCountNeg;
		}
	}
	return 0;
}


uint SpansPoint::lookup(std::string chrName, uint pos, bool direction) const {
	auto it_pos = std::lower_bound(chrName_pos.at(chrName).begin(), chrName_pos.at(chrName).end(), pos);
	if (it_pos == chrName_pos.at(chrName).end() || *it_pos != pos) {
		// throw not-found/out-of-bounds exception?
		throw std::out_of_range("Pos not found - SpansPoint::lookup");
		return 0;
	}else{
		// Then use that offset into the other vectors.
		return chrName_count[direction].at(chrName).at(it_pos - chrName_pos.at(chrName).begin());
	}
}

uint SpansPoint::lookup(std::string chrName, uint pos) const {
	//	std::map<string, std::vector<int>> chrName_pos;
	auto it_pos = std::lower_bound(chrName_pos.at(chrName).begin(), chrName_pos.at(chrName).end(), pos);
	if (it_pos == chrName_pos.at(chrName).end() || *it_pos != pos) {
		// throw not-found/out-of-bounds exception?
		throw std::out_of_range("Pos not found - SpansPoint::lookup");
		return 0;
	}else{
		// Then use that offset into the other vectors.
		return (
			chrName_count[0].at(chrName).at(it_pos - chrName_pos.at(chrName).begin())
			+ chrName_count[1].at(chrName).at(it_pos - chrName_pos.at(chrName).begin())
			);
	}
}


void SpansPoint::setSpanLength(uint overhang_left, uint overhang_right) {
	overhangLeft = overhang_left;
	overhangRight = overhang_right;
	overhangTotal = overhang_right + overhang_left;
}

void SpansPoint::ProcessBlocks(const FragmentBlocks &blocks) {
	std::vector<uint>::iterator it_position;

	//Walk each read within the fragment (1 or 2).
	for (int index = 0; index < blocks.readCount; index ++) {
		//Walk each block within each read.
		for (uint j = 0; j < blocks.rLens[index].size(); j++) {
			if ( blocks.rLens[index][j] > overhangTotal ) {
				//Block is long enough it may sufficiently overhang a point of interest.
				it_position = std::upper_bound(
						(*chrID_pos.at(blocks.chr_id)).begin(),
						(*chrID_pos.at(blocks.chr_id)).end(),
						blocks.readStart[index] + blocks.rStarts[index][j] + overhangLeft - 1
				);  // -1 --- as the test is > rather than >=.
				while (it_position != (*chrID_pos.at(blocks.chr_id)).end() && *it_position < (blocks.readStart[index] + blocks.rStarts[index][j] + blocks.rLens[index][j] )) {
					//increment corresponding counter.
					(*chrID_count[blocks.direction].at(blocks.chr_id)).at(it_position - (*chrID_pos.at(blocks.chr_id)).begin())++;
					it_position++;
				}
			}
		}
	}
}


void SpansPoint::loadRef(std::istream &IN) {
	// TODO: will we ever want to store some additional info -- eg: String name of each position? Not right now.
	std::string myLine;
	std::string myField;
	myLine.reserve(1000);
	myField.reserve(100);
	int pos;
	string s_chr;
	s_chr.reserve(30);
	string direction;

	while ( !IN.eof() && !IN.fail() ) {
		getline(IN, myLine, '\n');

		if (IN.eof() || IN.fail()) {
			if (myLine.length() == 0) {
				// This line is empty - just a blank line at the end of the file.
				// Checking at this stage allows correct handling of files both with and without a trailing \n after the last record.
				break;
			}else{
				// Error line in input, ignore.
				break;
			}
		}

		std::istringstream lineStream;
		lineStream.str(myLine);
		
		getline(lineStream, s_chr, '\t');
		getline(lineStream, myField, '\t');
		pos = stol(myField);
		
		getline(lineStream, direction, '\t');
		
		chrName_pos[s_chr].push_back(pos);
	}
	
	for (std::map<string, std::vector<uint>>::iterator it_chr=chrName_pos.begin(); it_chr!=chrName_pos.end(); it_chr++) {	
		std::sort( it_chr->second.begin(), it_chr->second.end() );
		// We now have chrName_pos sorted by position.
		chrName_count[0][it_chr->first].resize(it_chr->second.size(), 0);
		chrName_count[1][it_chr->first].resize(it_chr->second.size(), 0);
		// Just created a vector of the same size as the position one to store the counter.
		// initialise the count vector to the same length, with zero start.
	}
}

void SpansPoint::ChrMapUpdate(const std::vector<std::string> &chrmap) {
	chrID_pos.resize(0);
	chrID_count[0].resize(0);
	chrID_count[1].resize(0);
	for (uint i = 0; i < chrmap.size(); i++) {
		chrID_pos.push_back( &chrName_pos[chrmap.at(i)] );
		chrID_count[0].push_back( &chrName_count[0][chrmap.at(i)] );
		chrID_count[1].push_back( &chrName_count[1][chrmap.at(i)] );
	}
}


void FragmentsInROI::ChrMapUpdate(const std::vector<std::string> &chrmap) {
	chrID_ROI.resize(0);
	chrID_count[0].resize(0);
	chrID_count[1].resize(0);
	for (uint i = 0; i < chrmap.size(); i++) {
		chrID_ROI.push_back( &chrName_ROI[chrmap.at(i)] );
		chrID_count[0].push_back( &chrName_count[0][chrmap.at(i)] );
		chrID_count[1].push_back( &chrName_count[1][chrmap.at(i)] );
	}
}

int FragmentsInROI::WriteOutput(std::ostream *os) const {
	for (std::map<string, ulong>::const_iterator itID=RegionID_counter[1].begin(); itID!=RegionID_counter[1].end(); ++itID) {
		*os << itID->first << "\t" << (itID->second + RegionID_counter[0].at(itID->first)) << "\t" << itID->second << "\t" << RegionID_counter[0].at(itID->first) << "\n";
		//Outputs tab separated: ROIname, total hits, positive-strand hits, negative-strand hits.
	}
	return 0;
}

void FragmentsInROI::loadRef(std::istream &IN) {

	std::string myLine;
	std::string myField;
	myLine.reserve(1000);
	myField.reserve(100);
	int start;
	int end;
	string s_chr;
	s_chr.reserve(30);
	string s_name;
	s_name.reserve(200);

	while ( !IN.eof() && !IN.fail() ) {
		// Input ref:  chr - start - end - name - dir(?). (name\tdir could be considered a single variable)
		getline(IN, myLine, '\n');

		if (IN.eof() || IN.fail()) {
			if (myLine.length() == 0) {
				// This line is empty - just a blank line at the end of the file.
				// Checking at this stage allows correct handling of files both with and without a trailing \n after the last record.
				break;
			}else{
				// Error line in input, ignore.
				break;
			}
		}

		std::istringstream lineStream;
		lineStream.str(myLine);
		
		getline(lineStream, s_chr, '\t');
		getline(lineStream, myField, '\t');
		start = stol(myField);
		getline(lineStream, myField, '\t');
		end = stol(myField);

		getline(lineStream, s_name, '\t');

		chrName_ROI[s_chr].push_back(std::make_pair(end, start));
		chrName_count[0][s_chr].push_back(&RegionID_counter[0][s_name]);
		chrName_count[1][s_chr].push_back(&RegionID_counter[1][s_name]);

	}
}

void FragmentsInROI::ProcessBlocks(const FragmentBlocks &blocks) {
	std::vector<std::pair<uint,uint>>::iterator it_ROI;

	uint frag_start = blocks.readStart[0];
	uint frag_end = blocks.readEnd[0];
	if (blocks.readCount > 1 && blocks.readEnd[1] > frag_end) {
		frag_end = blocks.readEnd[1];
	}
	
	// Frag start, Frag end.
	// See if this is fully inside one of the ref-regions.

	it_ROI = std::lower_bound(
			(*chrID_ROI.at(blocks.chr_id)).begin(),
			(*chrID_ROI.at(blocks.chr_id)).end(),
			std::make_pair(frag_end, frag_end)
	);
	
	if (it_ROI != (*chrID_ROI.at(blocks.chr_id)).end() ) {
		if (frag_start >= it_ROI->second && frag_end <= it_ROI->first) {
			(*(*chrID_count[blocks.direction].at(blocks.chr_id)).at(it_ROI - (*chrID_ROI.at(blocks.chr_id)).begin()))++;			
		}
	}
}



void FragmentsInChr::ProcessBlocks(const FragmentBlocks &blocks) {
	(*chrID_count.at(blocks.chr_id))[blocks.direction]++;
}

void FragmentsInChr::ChrMapUpdate(const std::vector<string> &chrmap) {
	chrID_count.resize(0);
	for (uint i = 0; i < chrmap.size(); i++) {
		chrName_count[chrmap.at(i)].resize(2); // This data structure isn't auto initializing - unlike all the other structures. Or maybe just a vector can't access via [] until a position exists? But a map is fine. Makes sense.
		chrID_count.push_back( &chrName_count[chrmap.at(i)] );
	}
}

int FragmentsInChr::WriteOutput(std::ostream *os) const {
	for (auto itChr=chrName_count.begin(); itChr!=chrName_count.end(); itChr++) {
		*os << itChr->first << "\t"
			<< ((itChr->second)[1] + (itChr->second)[0]) << "\t"
			<< (itChr->second)[1] << "\t"
			<< (itChr->second)[0] << "\n";
	}
	return 0;
}
