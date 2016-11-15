#include "ReadBlockProcessor.h"
#include "ReadBlockProcessor_CoverageBlocks.h"
#include "ReadBlockProcessor_OutputBAM.h"
#include "BAM2blocks.h"
#include "includedefine.h"

// Ask for an older version of "memcpy" so we don't need new GLIBC
// (GLIBC_2.2.5) memcpy
// __asm__(".symver memcpy,memcpy@GLIBC_2.2.5");

int main(int argc, char * argv[])
{

	if (argc != 7) {
    	cerr << "Usage: gzip -cd < input.bam | IRFinder output-directory ref-coverage.bed ref-SJ.ref ref-spans-point.ref ROI-named.bed(or NULL) out.bam(or NULL)" << endl;
    	exit(1);
	}

	std::string outputDir = argv[1];
	std::string s_inCoverageBlocks = argv[2];
	std::string s_inSJ = argv[3];
	std::string s_inSpansPoint = argv[4];
	std::string s_inROI = argv[5];
	std::string s_outBAM = argv[6];
	
	cout << "IRFinder run with options:\nOutput Dir\t" << outputDir << "\n"
		<< "Main intron reference:\t" << s_inCoverageBlocks << "\n"
		<< "Read spans reference:\t" << s_inSpansPoint << "\n"
		<< "Optional ROI reference:\t" << s_inROI << "\n"
		<< "Optional BAM output:\t" << s_outBAM << "\n";
		
//	std::ofstream BED12;
//	BED12Output oBED12Output;
	FragmentsInROI oFragmentsInROI;
	FragmentsInChr oFragmentsInChr;

	JunctionCount oJuncCount;
	std::ifstream inJuncCount;
	inJuncCount.open(s_inSJ, std::ifstream::in);
	oJuncCount.loadRef(inJuncCount);
	inJuncCount.close();

	SpansPoint oSpansPoint;
	oSpansPoint.setSpanLength(5,4);
	std::ifstream inSpansPoint;
	inSpansPoint.open(s_inSpansPoint, std::ifstream::in);
	oSpansPoint.loadRef(inSpansPoint);
	inSpansPoint.close();

	CoverageBlocksIRFinder oCoverageBlocks;
	std::ifstream inCoverageBlocks;
	inCoverageBlocks.open(s_inCoverageBlocks, std::ifstream::in);
	oCoverageBlocks.loadRef(inCoverageBlocks);
	inCoverageBlocks.close();
	

	BAM2blocks BB;

	OutputBAM oOutputBAM;
	std::ofstream BAMout;
	if (s_outBAM != "NULL") {
		BAMout.open(s_outBAM, std::ifstream::out);
		oOutputBAM.SetOutputHandle(&BAMout);
		BB.registerCallbackProcessBlocks( std::bind(&OutputBAM::ProcessBlocks, &oOutputBAM, std::placeholders::_1) );
	}
// 	if (s_outBED12 != "NULL") {
// 	  	BED12.open(s_outBED12, std::ifstream::out);
// 		oBED12Output.SetOutputStream(&BED12);
// 		BB.registerCallbackChrMappingChange( std::bind(&BED12Output::ChrMapUpdate, &oBED12Output, std::placeholders::_1) );
// 		BB.registerCallbackProcessBlocks( std::bind(&BED12Output::ProcessBlocks, &oBED12Output, std::placeholders::_1) );
// 	}



	BB.registerCallbackChrMappingChange( std::bind(&JunctionCount::ChrMapUpdate, &oJuncCount, std::placeholders::_1) );
	BB.registerCallbackProcessBlocks( std::bind(&JunctionCount::ProcessBlocks, &oJuncCount, std::placeholders::_1) );

	BB.registerCallbackChrMappingChange( std::bind(&FragmentsInChr::ChrMapUpdate, &oFragmentsInChr, std::placeholders::_1) );
	BB.registerCallbackProcessBlocks( std::bind(&FragmentsInChr::ProcessBlocks, &oFragmentsInChr, std::placeholders::_1) );

	BB.registerCallbackChrMappingChange( std::bind(&SpansPoint::ChrMapUpdate, &oSpansPoint, std::placeholders::_1) );
	BB.registerCallbackProcessBlocks( std::bind(&SpansPoint::ProcessBlocks, &oSpansPoint, std::placeholders::_1) );

	if (s_inROI != "NULL") {	
		std::ifstream inFragmentsInROI;
		inFragmentsInROI.open(s_inROI, std::ifstream::in);
		oFragmentsInROI.loadRef(inFragmentsInROI);
		inFragmentsInROI.close();

		BB.registerCallbackChrMappingChange( std::bind(&FragmentsInROI::ChrMapUpdate, &oFragmentsInROI, std::placeholders::_1) );
		BB.registerCallbackProcessBlocks( std::bind(&FragmentsInROI::ProcessBlocks, &oFragmentsInROI, std::placeholders::_1) );
	}

	BB.registerCallbackChrMappingChange( std::bind(&CoverageBlocks::ChrMapUpdate, &oCoverageBlocks, std::placeholders::_1) );
	BB.registerCallbackProcessBlocks( std::bind(&CoverageBlocks::ProcessBlocks, &oCoverageBlocks, std::placeholders::_1) );

	BB.openFile(&std::cin); // This file needs to be a decompressed BAM. (setup via fifo / or expect already decompressed via stdin).

	if (s_outBAM != "NULL") {
		oOutputBAM.OutputHeader(BB.samHeader, BB.chr_names, BB.chr_lens);
	}

	BB.processAll();

	if (s_outBAM != "NULL") {
		oOutputBAM.FlushOutput(1);
		BAMout.flush(); BAMout.close();
	}

// 	if (s_outBED12 != "NULL") {
// 		// Output stream cleanup.
// 		BED12.flush(); BED12.close();
// 	}

	if (s_inROI != "NULL") {
		// Output computed statistics from data structures.
		//oFragmentsInROI -- this tells us if the data was directional or not -- if we need to know for other output modules.
		std::ofstream outFragmentsInROI;
		outFragmentsInROI.open(outputDir + "/IRFinder-ROI.txt", std::ifstream::out);
		oFragmentsInROI.WriteOutput(&outFragmentsInROI);
		outFragmentsInROI.flush(); outFragmentsInROI.close();
	}

	std::ofstream outJuncCount;
	outJuncCount.open(outputDir + "/IRFinder-JuncCount.txt", std::ifstream::out);
	oJuncCount.WriteOutput(&outJuncCount);
	outJuncCount.flush(); outJuncCount.close();

	int directionality = oJuncCount.Directional();
	cout << "RNA-Seq directionality -1/0/+1:\t" << directionality << "\n";
	
	
	std::ofstream outSpansPoint;
	outSpansPoint.open(outputDir + "/IRFinder-SpansPoint.txt", std::ifstream::out);
	oSpansPoint.WriteOutput(&outSpansPoint);
	outSpansPoint.flush(); outSpansPoint.close();

	std::ofstream outFragmentsInChr;
	outFragmentsInChr.open(outputDir + "/IRFinder-ChrCoverage.txt", std::ifstream::out);
	oFragmentsInChr.WriteOutput(&outFragmentsInChr);
	outFragmentsInChr.flush(); outFragmentsInChr.close();

	std::ofstream outCoverageBlocks;
	outCoverageBlocks.open(outputDir + "/IRFinder-IR-nondir.txt", std::ifstream::out);
	oCoverageBlocks.WriteOutput(&outCoverageBlocks, oJuncCount, oSpansPoint);
	outCoverageBlocks.flush(); outCoverageBlocks.close();

	if (directionality != 0) {
		outCoverageBlocks.open(outputDir + "/IRFinder-IR-dir.txt", std::ifstream::out);
		oCoverageBlocks.WriteOutput(&outCoverageBlocks, oJuncCount, oSpansPoint, directionality); // Directional.
		outCoverageBlocks.flush(); outCoverageBlocks.close();
	}

	
}

