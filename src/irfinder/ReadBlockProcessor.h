#ifndef CODE_READBLOCKPROCESSOR
#define CODE_READBLOCKPROCESSOR

#include "FragmentBlocks.h"

/*
The code can be finished faster if we force a requirement that all input files are coordinate sorted by the start of each block.
ie: sort -k2,2n (for BED files).
Chromosome sorted or not won't matter, as these get split into different vectors in all cases.
*/



class ReadBlockProcessor {
	public:
		virtual void ProcessBlocks(const FragmentBlocks &fragblock) = 0;
		virtual void ChrMapUpdate(const std::vector<std::string> &chrmap) = 0; //Maybe some of these funcs shouldn't be pure virtual - overloadable if needed, but default often ok.
};


class BED12Output : public ReadBlockProcessor {
	private:
		std::vector<std::string> chr_names;
		std::ostream* out;
	public:
		void ProcessBlocks(const FragmentBlocks &fragblock);
		void ChrMapUpdate(const std::vector<std::string> &chrmap);
		void SetOutputStream(std::ostream *os);
};


class JunctionCount : public ReadBlockProcessor {
	private:
		std::map<string, std::map<std::pair<uint,uint>,uint[3]>> chrName_junc_count;
		std::vector<std::map<std::pair<uint,uint>,uint[3]>*> chrID_junc_count;
		//uint[3] - 0, neg strand count; 1, pos strand count; 2 = expected direction from ref: 0=unknown, 1=neg, 2=pos.

		std::map<string, std::map<uint,uint[2]>> chrName_juncLeft_count;
		std::vector<std::map<uint,uint[2]>*> chrID_juncLeft_count;

		std::map<string, std::map<uint,uint[2]>> chrName_juncRight_count;
		std::vector<std::map<uint,uint[2]>*> chrID_juncRight_count;
		  //chrID_... stores a fast access pointer to the appropriate structure in chrName_... 
	public:
		void ProcessBlocks(const FragmentBlocks &fragblock);
		void ChrMapUpdate(const std::vector<std::string> &chrmap);
		int WriteOutput(std::ostream *os) const;
		void loadRef(std::istream &IN); //loadRef is optional, it allows directional detection to determine not just non-dir vs dir, but also which direction.

		int Directional() const;
		
		uint lookup(std::string ChrName, uint left, uint right, bool direction) const;
		uint lookup(std::string ChrName, uint left, uint right) const;
		uint lookupLeft(std::string ChrName, uint left, bool direction) const;
		uint lookupLeft(std::string ChrName, uint left) const;
		uint lookupRight(std::string ChrName, uint right, bool direction) const;
		uint lookupRight(std::string ChrName, uint right) const;

// Ideally we would read the XS junction strand attribute from the BAM if we want to count junctions from non-directional sequencing.
//   that will require BAM2blocks to be informed it should read the optional attributes looking for that attrib in that case.
// -- or we can just ignore direction -- the splice start/end information effectively determines the XS info (by ref to the reference)
};


class SpansPoint : public ReadBlockProcessor {
	private:
		std::map<string, std::vector<uint>> chrName_pos;
		std::map<string, std::vector<uint>> chrName_count[2];
		std::vector<std::vector<uint>*> chrID_pos;
		std::vector<std::vector<uint>*> chrID_count[2];
		char overhangLeft;
		char overhangRight;
		char overhangTotal;
		//chrID_... stores a fast access pointer to the appropriate structure in chrName_... 
	public:
		void setSpanLength(uint overhang_left, uint overhang_right);
		void loadRef(std::istream &IN);
		void ProcessBlocks(const FragmentBlocks &fragblock);
		void ChrMapUpdate(const std::vector<std::string> &chrmap);
		//void SetOutputStream(std::ostream *os);
		int WriteOutput(std::ostream *os) const;
		uint lookup(std::string ChrName, uint pos, bool direction) const;
		uint lookup(std::string ChrName, uint pos) const;
};

class FragmentsInChr : public ReadBlockProcessor {
	// Counts the number of fragments in each Chromosome. (for both + & - strands).
	private:
		std::map<string, std::vector<uint>> chrName_count; //only expecting 2 items in our vector.
		std::vector<std::vector<uint>*> chrID_count;
	public:
		void ProcessBlocks(const FragmentBlocks &blocks);
		void ChrMapUpdate(const std::vector<string> &chrmap);
		int WriteOutput(std::ostream *os) const;		
};


class FragmentsInROI : public ReadBlockProcessor {
	// Counts the number of fragments fully contained within a ROI.
	//   the ROIs may not overlap. Direction ignored for overlap detect.
	private:
		std::map<string, ulong> RegionID_counter[2];
 
		std::map<string, std::vector<std::pair<uint,uint>>> chrName_ROI;
		std::map<string, std::vector<ulong*>> chrName_count[2];

		std::vector<std::vector<std::pair<uint,uint>>*> chrID_ROI;
		std::vector<std::vector<ulong*>*> chrID_count[2];

		// Perhaps we want to store some text relating to each record too? Easy to do if the input is pre-sorted (at least within each Chr).
		//   if pre-sorted, it may be easier to check for no overlapping blocks on read .. or can do this immediately after read with a single nested-walk.
		std::map<string, std::vector<string>> chrName_ROI_text;
	public:
		void ProcessBlocks(const FragmentBlocks &blocks);
		void ChrMapUpdate(const std::vector<string> &chrmap);
		void loadRef(std::istream &IN);
		int WriteOutput(std::ostream *os) const;		
};


/*
class CoverageBlocks : public ReadBlockProcessor { ... }
// In it's own file -- bigger code.
*/

#endif
