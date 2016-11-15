#ifndef CODE_FRAGMENTBLOCKS
#define CODE_FRAGMENTBLOCKS

#include "includedefine.h"

/* A class to store up to 2 reads belonging to a single fragment.
 * It is a storage class, almost a struct, it does not perform processing itself.
 * Read1 is always valid.
 * Read2 is only valid if readCount == 2.
 *
 * There may only be a single read if:
 *  - the sequencing is single end rather than paired end..
 *  - the sequencing is paired end, but the two reads overlapped and have been combined
 *    into a single synthetic read / block of coverage.
 */
class FragmentBlocks {
	private:
		static const int initial_alloc = 100;
		static const int max_read_name = 300;
		std::vector<std::string> chr_names; //TODO - this is currently unused??
	public:
		FragmentBlocks();
		const std::string chrName() const;
		void ChrMapUpdate(const std::vector<std::string>& chrmap);

		std::string readName;
		std::vector<int> rStarts[2];
		std::vector<int> rLens[2];
		uint readStart[2];
		uint readEnd[2];
		int readCount;
		uint chr_id; // Assumption that both r1 & r2 are on the same chromosome?
					//   if they aren't we shouldn't process them as a single fragment.
					//   perhaps a sanity check in pairing, only treat them as a pair
					//   if the name of the reads matches and the Chr matches.
		bool direction;
};

#endif