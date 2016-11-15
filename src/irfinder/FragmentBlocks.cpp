#include "FragmentBlocks.h"
#include "includedefine.h"
// using namespace std;


// This class is an information storage container only -- pretty much a struct.
// It allows all the relevant information relating to an interpreted fragment to be passed
// to the variety of callback watchers that require fragment blocks to update their stats.

FragmentBlocks::FragmentBlocks() {
	rStarts[0].reserve(initial_alloc);
	rLens[0].reserve(initial_alloc);
	rStarts[1].reserve(initial_alloc);
	rLens[1].reserve(initial_alloc);
	readName.reserve(max_read_name);
	readCount = 0;
}

// Return a string representation of the Chromosome name.
const std::string FragmentBlocks::chrName() const {
	return chr_names.at(chr_id);
}

// Update the internal data structure with a new mapping between Chromosome ID# and Chromosome name (string).
void FragmentBlocks::ChrMapUpdate(const std::vector<string> &chrmap) {
	chr_names = chrmap;
}
