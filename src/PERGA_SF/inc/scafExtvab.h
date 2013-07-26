/*
 * scafExt.h
 *
 *  Created on: Jun 11, 2011
 *      Author: zhuxiao
 */

#ifndef SCAFEXT_H_
#define SCAFEXT_H_ 1

#include "scafStructure.h"


//=========== global variables ============
extern char inputPathStr[256];
extern char outputPathStr[256];
extern char outputPrefix[256];

extern char contigsFile[256];
extern char contigsFileFiltered[256];
extern char contigsFileUltraShort[256];
extern char contigsFileAfterOverlap[256];
extern char contigsFileAfterGapFilling[256];

extern char readFilesFastqPE[2][256];				// [0]-- xxx_1.fastq, [1]-- xxx_2.fastq
extern char *readFilesInput[256];
extern int readFileNum;
extern int readsFileFormatType;
extern char readMatchFiles[2][256];				// [0]-- xxx_1.fastq.match.bin, [1]-- xxx_2.fastq.match.bin
extern char readSeqFile[256];
extern char monoReadSeqFile[256];

extern char contigIndexFile[256];
extern char readListFiles[4][256];						//[0]-- RL1, [1]-- RL2, [2]-- SRL, [3]-- MRL
extern char readListFilesAfterOverlap[4][256];		//[0]-- RL1, [1]-- RL2, [2]-- SRL, [3]-- MRL
extern char contigListFile[256];
extern char contigListFileAfterOverlap[256];

extern char meanSdevFile[256];
extern char contigOverlapFile[256];
extern char contigOverlapFileAfterFapFilling[256];
extern char linkInfoFile[256];
extern char averLinkNumFile[256];

extern char scafGraphFile[256];
extern char scafSeqFile[256];

extern int averReadLen;
extern int readLen;
extern int entriesPerRead;
extern uint64_t lastEntryMaskRead;
extern int lastEntryBaseNumRead;

extern int contigEndLen;				// contig end length to build contigs
extern int minContigLenThres;			// the minimal contig length threshold
extern int gapFillFlag;
extern int pairedMode;

extern double meanSizeInsert;						// the mean size of the insert fragments
extern double stardardDeviationInsert;				// the standard deviation of the insert fragments


//================== variables for scafContigIndex.c ==================
extern int contigsNum;
extern baseSeqIndex *baseSeqIndexArr;
extern uint64_t *seqArr;		// sequence array
extern int itemNumSeqArr;		// the item number in sequence array
extern int maxItemNumSeqArr;

extern int uniqueSeqNum;
int rightShift_bitsNum;
extern uniqueSeqKmer *uniqueSeqKmerArr;
extern uint64_t *uniqueSeqArr;
extern seqRow *seqRowArr;
extern contigMatchInfo *contigMatchInfoArr;

extern sortkmertype *sortKmerArr;
extern int *sortRowIndexArr;


//===================== variables for scafReadList.c =======================
extern int64_t readNumInRL[3];			// the read numbers of RL
extern int64_t matchItemNumInRP[3];		// the maximal read numbers of RL
extern ReadList *readListArr[3];		// pointer array to the read List
extern ReadPos *readPosArr[3];			// pointer array to the read position array


//===================== variables for scafContigList.c =======================
extern int64_t readItemNumInSRL;		// the number of reads in shared read list (SRL)
extern int64_t matchItemNumInSRP;		// the number of match items in shared read position array(SRP)
extern ReadList *sharedReadListArr;	// pointer to the shared read List array
extern ReadPos *sharedReadPosArr;		// pointer to the shared read position array

extern int64_t contigItemNumInCL;		// the number of contigs in contig list (CL) array
extern int64_t contigReadItemNumInCR;	// the number of contigs in contig list array
extern ContigList *contigListArr;		// pointer to contig list (CL) array
extern ContigRead *contigReadArr;		// pointer to contig read (CR) array


//=================== variables for scafLinking.c =========================
extern contigInfo *contigInfoArr;			// the contig information array

extern arrLoc *arrLocArr;
extern int64_t itemNumArrLocArr;
extern int64_t maxItemNumArrLocArr;

extern ContigGraph *contigGraphArr;
extern int64_t itemNumContigGraphArr;

extern double maxLinksNumContigsThres;		// the maximal number of links between contigs , added 2012-11-24
extern double minLinksNumContigsThres;			// the minimal number of links between contigs, it is also used in scafContigMerge.c
extern double maxRatioSecondFirstLinkNum;	// the maximal ratio of the second link number to the maximal link number between contigs
extern double secondLinkNumFactor;			// the second link number factor

extern double averLinkNum;
extern double maxSecondLinkNumThres;
extern double minBreakNumThres;

extern int itemNumContigLinkArr;			// item number of contig link array
extern int headRowContigLinkArr;			// head row of contig link array
extern int tailRowContigLinkArr;			// tail row of contig link array
extern contigLink *contigLinkArr;			// the contig link array

//================ New Begin ===================
extern linkCandidate *linkCandidateArr;
extern int64_t itemNumLinkCandidateArr;
extern int64_t maxItemNumLinkCandidateArr;
//================ New End ===================

//=================== variables for scafContigMerge.c =========================
//extern int gapFillingFlag;				// YES or NO: YES-- the gap will be filled, No-- the gap will not be filled.

extern int maxOverlapSeqLen;						// the maximal length of sequence in contig ends to be overlapped
extern char *overlapSeq1;							// the first sequence for overlapping
extern char *overlapSeq2;							// the second sequence for overlapping
extern int *scoreArr;								// the score array for pairwise alignment
extern char *alignResultArr[3];						// the alignment result array

//extern double meanSizeInsert;						// the mean size of the insert fragments
//extern double stardardDeviationInsert;				// the standard deviation of the insert fragments
extern double standardDeviationFactor;				// the factor of the standard deviation of the insert fragments

extern int minOverlapThres;							// the minimal overlap threshold
extern int minExactOverlapThres;					// the minimal exact overlap threshold
extern int mismatchThres;							// the mismatch threshold
extern double maxMisMatchRatioThres;				// the maximal reatio threshold of mismatched bases
extern int subMismatchThres;						// the sub mismatch threshold
extern int minAdjustGapSizeThres;					// the minimal adjusted gap size threshold
extern int maxAdjustGapSizeThres;					// the maximal adjusted gap size threshold
extern int minBaseNumInGap;							// the minimal bases in gap region
extern double gapSizeSdevFactorOverlap;				// the standard deviation factor for gap size in contig overlap
extern double gapSizeSdevFactorGapFilling;			// the standard deviation factor for gap size in gap filling
//extern double exactOverlapSdevThres;				// the exact overlap standard deviation range threshold

extern int breakLinkNumThres;						// the break links number threshold

extern int matchScore;								// the match score
extern int mismatchScore;							// the mismatch score
extern int gapScore;								// the gap score

extern int scaffoldsNumInCOI;						// total number of scaffolds
extern contigOverlapIndex *contigOverlapIndexArr;	// pointer to contig overlap index array

extern int itemNumInContigOverlapArr;
extern contigOverlap *contigOverlapInfoArr;			// pointer to contig overlap information array


//=================== variables for scafGraph.c =========================
extern int64_t readItemNumInMRL;					// the number of reads in mono read list (SRL)
extern int64_t matchItemNumInMRP;					// the number of match items in mono read position array(SRP)
extern ReadList *monoReadListArr;					// pointer to the mono read List array
extern ReadPos *monoReadPosArr;						// pointer to the mono read position array

extern scafGraph *scafGrapDeBruijn;					// the de Bruijn graph in scaffolding


//=================== variables for scafGgpFilling.c =========================
extern scafKmer *scafKmers[2];						// scafKmers: [0] -- scafKmer with plus orientation; [1] -- scafKmer with minus orientation
extern int kmer_len;

extern int maxScafAssemblingReadsNumInArr;			// the maximal number of reads in scafAssemblingReadArr
extern int maxScafSuccessReadsNumInArr;				// the maximal number of successful reads in scafSuccessReadArr
//extern int maxScafContigEndSeqLen;					// the maximal length of scafContig end sequence
//extern int maxComparisonSeqLenInScaf;				// the maximal comparison sequence length in scaffolding

extern int successFilledFlag;						// YES or NO

extern scafAssemblingRead *scafAssemblingReadArr;	// the decision table
extern int scafAssemblingReadsNum;					// reads number in scafAssemblingReadArr
extern scafSuccessRead *scafSuccessReadArr;			// the successful reads array in scaffolding
extern int scafSuccessReadsNum;						// the successful reads number in scaffolding

extern scafContig *scafContighead;					// pointer to the head of scafContigs
extern scafContig *scafContigtail;  				// pointer to the tail of scafContigs
extern scafContig *scafContigLastReadLen;			// pointer to the scafContig node having READ_LEN nodes before scafContig tail
extern int scafContigIndex;
extern scafContig *successScafContig;				// pointer to the most recent successful scafContig node
extern char *scafContigSeqLastReadLen;				// the last base sequence with length READ_LEN
extern int scafContigSeqLenLastReadLen;				// the base number of last base sequence

extern char *scafContigEndSeqArr[2];					// the two contig end sequences with the maximal length of READ_LEN
extern int scafContigEndSeqLenArr[2];					// the sequences length of the two contig end sequences with the maximal length of READ_LEN
extern char *comparisonSeqInScaf;					// the comaprison sequence in scaffolding
extern int comparisonSeqLenInScaf;					// the comparison sequence length

extern int prepareAssemblyLenArr[2];				// the assembly length of prepared assembly in contig ends

extern int lockedReadsNumInScaf;					// the number of locked reads in scaffolding


//++++++++++++++++++++++++++++++++++++
extern int PEGivenType, oldPEGivenType;

extern int occsNumSE[4];  //下一个kmer的得分，对应ACGT的顺序
extern int occsNumPE[4];

extern int longKmerSize;
extern int longKmerStepSize;
extern double averKmerOcc;
extern double firstKmerThres;
extern double minKmerOccSE;
extern double minKmerOccPE;
extern double maxSecondOcc;
extern double maxFirstOcc;
extern double minLongKmerOcc;
extern double lockedReadsNumThres;
extern double maxOccNumFaiedPE;

extern int errorRegLenEnd3;				// erroneous region length of 3' end of a read

//++++++++++++++++++++++++++++++++++++
extern int kmerSize;
extern int entriesPerKmer;
extern uint64_t lastEntryMaskKmer;
extern int lastEntryBaseNumKmer;		// the base number of the last entry of kmer sequence array

extern uint64_t hashTableSize;

extern char kmerBaseSeq[1000];

extern uint64_t *kmerSeqInt;
extern uint64_t *kmerSeqIntRev;

extern uint64_t *kmerSeqIntAssembly;
extern uint64_t *kmerSeqIntAssemblyRev;
extern uint64_t *tmpKmerSeqIntAssembly;

extern scafKmer *firstKmer;
extern uint64_t tabooSeqInt[4];

extern int navigationID;			// 0 -- Single End, 1 -- Paired End, -1 -- unknown
extern int navigationNumSE;
extern int maxNavigationNumSE;


#endif /* SCAFEXT_H_ */
