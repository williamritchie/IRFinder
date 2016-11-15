#ifndef CODE_READBLOCKPROCESSOR_OUTPUTBAM
#define CODE_READBLOCKPROCESSOR_OUTPUTBAM

#include "ReadBlockProcessor.h"
#include "FragmentBlocks.h"



class OutputBAM : public ReadBlockProcessor {
	private:
		union stream_uint32 {
			char c[4];
			uint32_t i;
		};

		union stream_int32 {
			char c[4];
			int32_t i;
		};

		union stream_uint16 {
			char c[2];
			uint16_t i;
		};

		struct bam_read_core {
			union {
			  char c[36];
			  struct {
				int32_t block_size;
				int32_t refID;
				int32_t pos;
				uint8_t l_read_name;
				uint8_t mapq;
				uint16_t bin;
				uint16_t n_cigar_op;
				uint16_t flag;
				int32_t l_seq;
				int32_t next_refID;
				int32_t next_pos;
				int32_t tlen;
			  }; // anonymous struct to allow easy access to members.
			};
			char read_name[256];
			union {
			  char cigar_buffer[2000];
			  int32_t cigar[500];
			};
		};

		static const int bamEOFlength = 28;
		static const char bamEOF[bamEOFlength+1];

		static const int bamGzipHeadLength = 16;  // +2 a uint16 with the full block length.
		static const char bamGzipHead[bamGzipHeadLength+1];

		//bam_read_core bam[2];
		bam_read_core bamP;
		bam_read_core bamS;
	
		std::vector<std::string> chrID_Name;
		std::ostream *out;

		static const int bufferMax = 65500; //Some header and footer bytes are required to also fit within the 65535 length buffer.
		char buffer[65535];
		int bufferPos = 0;
		// Should the actual header and footer be written into this buffer too? That would allow a different thread to do the write if desired at a later stage.
		// At present the buffer is for the data part only.
		
		uint16_t reg2bin(int beg, int end);
		void CreateCigar(const std::vector<int> &starts, const std::vector<int> &lens,int32_t *cigar);
		void FeedBuffer(bam_read_core &bamrec);
	public:
		OutputBAM();
		void ProcessBlocks(const FragmentBlocks &fragblock);
		void ChrMapUpdate(const std::vector<string> &chrmap);
		void SetOutputHandle(std::ostream *os);
		void FlushOutput(const int final = 0);
		// output initial BAM header. "BAM\x01"
		void OutputHeader(const std::string &samHeader, const std::vector<std::string> chr_names, const std::vector<int32_t> chr_lens);
};

#endif
