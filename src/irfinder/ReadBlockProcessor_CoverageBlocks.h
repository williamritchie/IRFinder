#ifndef CODE_READBLOCKPROCESSOR_COVERAGEBLOCKS
#define CODE_READBLOCKPROCESSOR_COVERAGEBLOCKS

#include "CoverageBlock.h"
#include "ReadBlockProcessor.h"
#include "FragmentBlocks.h"

struct BEDrecord {
	std::string chrName;
	std::string name;
	uint start;
	uint end;
	bool direction;
	
	std::vector<std::pair<uint,uint>> blocks;
};


class CoverageBlocks : public ReadBlockProcessor {
	//Store the Blocked BED record for each ROI/intron. This won't be referred to again until the end.
	//XX Create the temporary vectors (per Chr) which simply list the blocks sequentially as read.
	//XX Sort the temporary vectors
	//XX Build the final vectors of "blocks of interest"
	//xx Delete the temporary vectors
	//xx Create the parallel vectors with counter objects. (do these as a batch at the end, once vector size is known - for best memory layout)
	//xx Process fragments against the counter structure. (have I already written a class/object for this?)
	
	//Produce summary statistical output for each Blocked BED record, using the counter structure.

	private:

		// Coverage depth data-structures.
		std::map<string, std::vector<CoverageBlock>> chrName_CoverageBlocks;
		// Shortcut pointers to depth data-structures.
		std::vector<std::vector<CoverageBlock>*> chrID_CoverageBlocks;

		// TODO: what is optimal for speed & memory usage?
//		static const uint coverage_block_max_length = 5000;
		static const uint coverage_block_max_length = 500;

	protected:
		std::vector<BEDrecord> BEDrecords;


	public:
		void ProcessBlocks(const FragmentBlocks &fragblock);
		void ChrMapUpdate(const std::vector<string> &chrmap);
		void loadRef(std::istream &IN);
		int WriteOutput(std::ostream *os) const;
		
		void fillHist(std::map<uint,uint> &hist, const std::string &chrName, const std::vector<std::pair<uint,uint>> &blocks) const;
		void fillHist(std::map<uint,uint> &hist, const std::string &chrName, const std::vector<std::pair<uint,uint>> &blocks, bool direction) const;

		double meanFromHist(const std::map<uint,uint> &hist) const;
		double coverageFromHist(const std::map<uint,uint> &hist) const;
		double percentileFromHist(const std::map<uint,uint> &hist, uint percentile) const;
		double trimmedMeanFromHist(const std::map<uint,uint> &hist, uint centerPercent) const;
};

class CoverageBlocksIRFinder : public CoverageBlocks {
	public:
		int WriteOutput(std::ostream *os, const JunctionCount &JC, const SpansPoint &SP, int directionality = 0) const;
};


#endif
