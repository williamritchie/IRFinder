#include "ReadBlockProcessor_OutputBAM.h"
#include "includedefine.h"
#include "crc32.h"

const char OutputBAM::bamEOF[OutputBAM::bamEOFlength+1] =
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00";
const char OutputBAM::bamGzipHead[OutputBAM::bamGzipHeadLength+1] = 
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00";


OutputBAM::OutputBAM() {
	bamP.mapq = 255;
	bamP.l_seq = 0;

	bamS.mapq = 255;
	bamS.l_seq = 0;
	bamS.tlen = 0;
	bamS.next_refID = -1;
	bamS.next_pos = -1;
}

void OutputBAM::CreateCigar(const std::vector<int> &starts, const std::vector<int> &lens,int32_t *cigar) {
	if (lens.size() < 1) return;
	*cigar = lens[0] << 4; //M opcode = 0;
	cigar++;
	// this should probably use iterators, though maybe not, complicates look-back.
	for (uint i = 1; i < lens.size(); i++) {
		*cigar = (starts[i] - starts[i-1] - lens[i-1]) << 4 | 3; //N opcode = 3;
		cigar++;
		*cigar = lens[i] << 4; //M opcode = 0;
		cigar++;
	}
}

void OutputBAM::ProcessBlocks(const FragmentBlocks &fragblock) {
	if (fragblock.readCount == 2) {
		//Pair.

		//Common for both reads in pair.
		bamP.refID = fragblock.chr_id;
		bamP.next_refID = fragblock.chr_id;
		bamP.l_read_name = fragblock.readName.length() + 1;
		memcpy(bamP.read_name, fragblock.readName.c_str(), fragblock.readName.length());
		bamP.read_name[fragblock.readName.length()] = 0; //set a trailing null.

		//Differ between reads. Read 1.
		bamP.tlen = max(fragblock.readEnd[0], fragblock.readEnd[1]) - fragblock.readStart[0];
		//in our model the first read always has a lower or equal coordinate compared to the second read.
		bamP.pos = fragblock.readStart[0];
		bamP.next_pos = fragblock.readStart[1]; //yes, we need to ref the other of the pair.
		bamP.bin = reg2bin(fragblock.readStart[0], fragblock.readEnd[0]);
		bamP.n_cigar_op = (2 * fragblock.rStarts[0].size()) - 1;
		//bam[0].flag
		//bamP.block_size = bamP.n_cigar_op + bamP.l_read_name + 36 - 4;
		CreateCigar(fragblock.rStarts[0], fragblock.rLens[0], bamP.cigar);
		if (fragblock.direction) {
			// +
			bamP.flag = 0x03 | 0x20 | 0x40;
		}else{
			// -
			bamP.flag = 0x03 | 0x10 | 0x40;		
		}
		FeedBuffer(bamP);


// Dir = +
// x03 | x20 | x40
// Dir = -
// x03 | x10 | x40

		//Differ between reads. Read 2.
		bamP.tlen = -bamP.tlen;
		bamP.pos = fragblock.readStart[1];
		bamP.next_pos = fragblock.readStart[0]; //yes, we need to ref the other of the pair.
		bamP.bin = reg2bin(fragblock.readStart[1], fragblock.readEnd[1]);
		bamP.n_cigar_op = (2 * fragblock.rStarts[1].size()) - 1;
		//bamP.block_size = bamP.n_cigar_op + bamP.l_read_name + 36 - 4;
		CreateCigar(fragblock.rStarts[1], fragblock.rLens[1], bamP.cigar);
		if (fragblock.direction) {
			// +
			bamP.flag = 0x03 | 0x10 | 0x80;
		}else{
			// -
			bamP.flag = 0x03 | 0x20 | 0x80;		
		}
		FeedBuffer(bamP);

// Dir = +
// x03 | x10 | x80
// Dir = -
// x03 | x20 | x80


// 		//bam[0].flag
	}else{
		//Single.

		bamS.l_read_name = fragblock.readName.length() + 1;
		memcpy(bamS.read_name, fragblock.readName.c_str(), fragblock.readName.length());
		bamS.read_name[fragblock.readName.length()] = 0; //set a trailing null.

		bamS.refID = fragblock.chr_id;

		bamS.pos = fragblock.readStart[0];
		bamS.bin = reg2bin(fragblock.readStart[0], fragblock.readEnd[0]);
		bamS.n_cigar_op = (2 * fragblock.rStarts[0].size()) - 1;
		//bam[0].flag
		//bamS.block_size = bamP.n_cigar_op + bamP.l_read_name + 36 - 4;
		CreateCigar(fragblock.rStarts[0], fragblock.rLens[0], bamS.cigar);
		if (fragblock.direction) {
			// +
			bamS.flag = 0;
		}else{
			// -
			bamS.flag = 0x10;		
		}
		FeedBuffer(bamS);
// Dir = +
// x00
// Dir = -
// x10
	}

	// Gzip standard header
	// Gzip Xtra: BS 0x02 0xXX (length of overall block)
	// Deflate:
	// 00 - no compression, not final deflate block.
	// 01 - no compression, final deflate block.
	// then 2 bytes length that follows. 2 bytes 1's complement of that length.
	// deflate blocks.
	//
	// Gzip footer: 0xXXXX CRC-32 of uncompressed stream (not of any of the gzip and deflate headers)
	//    0xXXXX length of uncompressed stream.
	//
	// CRC32 can be calculated as the buffer is filled (rather than at the end when we want to output -- ensures fastest possible mass output .. but if output is done on a separate thread, the checksum should be calculated there)
	//   but algo so fast anyway no difference? / or faster just to finish the calc than be in and out of functions.
	
}

void OutputBAM::FeedBuffer(bam_read_core &bamrec) {
	bamrec.block_size = 4*bamrec.n_cigar_op + bamrec.l_read_name + 36 - 4;
	if (bufferPos > bufferMax - bamrec.block_size) {
		FlushOutput(0);
	}
	memcpy(buffer+bufferPos, bamrec.c, 36);
	bufferPos += 36;
	memcpy(buffer+bufferPos, bamrec.read_name, bamrec.l_read_name);
	bufferPos += bamrec.l_read_name;
	memcpy(buffer+bufferPos, bamrec.cigar_buffer, 4*bamrec.n_cigar_op);
	bufferPos += 4*bamrec.n_cigar_op;
}

void OutputBAM::OutputHeader(const std::string &samHeader, const std::vector<std::string> chr_names, const std::vector<int32_t> chr_lens) {
	stream_int32 myInt32;
	myInt32.i = samHeader.length();
	memcpy(buffer+bufferPos, "BAM\x01",4);
	bufferPos += 4;
	memcpy(buffer+bufferPos, myInt32.c,4);
	bufferPos += 4;
	memcpy(buffer+bufferPos, samHeader.c_str(), samHeader.length());
	bufferPos += samHeader.length();

	myInt32.i = chr_names.size();
	memcpy(buffer+bufferPos, myInt32.c,4);
	bufferPos += 4;

	for (uint i = 0; i < chr_names.size(); i++) {
		myInt32.i = chr_names[i].length() + 1;
		memcpy(buffer+bufferPos, myInt32.c,4);
		bufferPos += 4;

		memcpy(buffer+bufferPos, chr_names[i].c_str(), chr_names[i].length());
		bufferPos += chr_names[i].length();
		memcpy(buffer+bufferPos, "\x00", 1);
		bufferPos += 1;
		myInt32.i = chr_lens[i];
		memcpy(buffer+bufferPos, myInt32.c,4);
		bufferPos += 4;
	}

	FlushOutput(0);
}

void OutputBAM::ChrMapUpdate(const std::vector<string> &chrmap) {

}

void OutputBAM::SetOutputHandle(std::ostream *os) {
	out = os;
}

void OutputBAM::FlushOutput(const int final) {
	stream_uint16 myInt16;
	stream_uint32 myInt32;
	CRC32 crc;

	crc.add(buffer, bufferPos);
	
	out->write(bamGzipHead, bamGzipHeadLength);
	myInt16.i = bufferPos + bamGzipHeadLength + 15 - 1;
	out->write(myInt16.c, 2);	
	out->write("\x01", 1);  //Deflate mode: no compression: binary 00. Final block: binary 1.
	myInt16.i = bufferPos;  //Deflate length.
	out->write(myInt16.c, 2);
	myInt16.i = ~myInt16.i; //1's complement deflate length.
	out->write(myInt16.c, 2);

	out->write(buffer, bufferPos); //Write the buffer itself.

	myInt32.i = crc.getRawHash(); //Uncompressed CRC32.
	out->write(myInt32.c, 4);
	myInt32.i = bufferPos; //Uncompressed Length.
	out->write(myInt32.c, 4);
	//Done - that's a complete BAM gzip block with a single block of compression zero deflate within.

	bufferPos = 0;	
	//if final, then output an empty gzip record at the end (BAM EOF marker).

	if (final) {
		out->write(bamEOF, bamEOFlength);	
	}
}

/* calculate bin given an alignment covering [beg,end) (zero-based, half-closed-half-open) */
uint16_t OutputBAM::reg2bin(int beg, int end) {
	--end;
	if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
	if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
	if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
	if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
	if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
	return 0;
}
