/*
 * scafGlobal.h
 *
 *  Created on: Jun 11, 2011
 *      Author: zhuxiao
 */

#ifndef SCAFGLOBAL_H_
#define SCAFGLOBAL_H_ 1

//#include "scafStdinc.h"
#include "scafStructure.h"


//=========== global variables ============
char inputPathStr[256];
char outputPathStr[256];
char outputPrefix[256];

char contigsFile[256];
char contigsFileFiltered[256];
char contigsFileUltraShort[256];
char contigsFileAfterOverlap[256];
char contigsFileAfterGapFilling[256];

char readFilesFastqPE[2][256];				// [0]-- xxx_1.fastq, [1]-- xxx_2.fastq
char *readFilesInput[256];
int readFileNum;
int readsFileFormatType;
char readMatchFiles[2][256];				// [0]-- xxx_1.fastq.match.bin, [1]-- xxx_2.fastq.match.bin
char readSeqFile[256];
char monoReadSeqFile[256];

char contigIndexFile[256];
char readListFiles[4][256];					//[0]-- RL1, [1]-- RL2, [2]-- SRL, [3]-- MRL
char readListFilesAfterOverlap[4][256];		//[0]-- RL1, [1]-- RL2, [2]-- SRL, [3]-- MRL
char contigListFile[256];
char contigListFileAfterOverlap[256];

char meanSdevFile[256];
char contigOverlapFile[256];
char contigOverlapFileAfterFapFilling[256];
char linkInfoFile[256];
char averLinkNumFile[256];

char scafGraphFile[256];
char scafSeqFile[256];

int averReadLen;
int readLen;
int entriesPerRead;
uint64_t lastEntryMaskRead;
int lastEntryBaseNumRead;

int contigEndLen;				// contig end length to build contigs
int minContigLenThres;			// the minimal contig length threshold
int gapFillFlag;
int pairedMode;

double meanSizeInsert;						// the mean size of the insert fragments
double stardardDeviationInsert;				// the standard deviation of the insert fragments


//=========== variables for scafContigIndex.c and scafMapping.c ============
int contigsNum = 0;				// the total number of contigs, and this variable is also used in scafContigLinking.c and scafContigMerge.c
baseSeqIndex *baseSeqIndexArr;	// the pointer to base sequence index array
uint64_t *seqArr;				// the pointer to sequence array
int itemNumSeqArr;				// the item number in sequence array
int maxItemNumSeqArr;			// the maximal item number of the sequence array

int uniqueSeqNum;				// the number of unique base sequence
int rightShift_bitsNum;			// the bits number of right shift for a k-mer
uniqueSeqKmer *uniqueSeqKmerArr;// the pointer to unique sequence k-mer array
uint64_t *uniqueSeqArr;			// the pointer to unique sequence array
seqRow *seqRowArr;				// the pointer to seqRow array
contigMatchInfo *contigMatchInfoArr;	// the pointer to contigMatchInfo array

sortkmertype *sortKmerArr;		// the sort k-mer array pointer
int *sortRowIndexArr;			// the sort row index array pointer


//=================== variables for scafReadList.c ==========================
int64_t readNumInRL[3];			// the read numbers of RL
int64_t matchItemNumInRP[3];	// the maximal read numbers of RL
ReadList *readListArr[3];		// pointer array to the read List
ReadPos *readPosArr[3];			// pointer array to the read position array


//=========== variables for scafContigList.c and scafContigLinking.c =============
int64_t readItemNumInSRL;		// the number of reads in shared read list (SRL)
int64_t matchItemNumInSRP;		// the number of match items in shared read position array(SRP)
ReadList *sharedReadListArr;	// pointer to the shared read List array
ReadPos *sharedReadPosArr;		// pointer to the shared read position array

int64_t contigItemNumInCL;		// the number of contigs in contig list (CL) array
int64_t contigReadItemNumInCR;	// the number of reads in contig list array
ContigList *contigListArr;		// pointer to contig list (CL) array
ContigRead *contigReadArr;		// pointer to contig read (CR) array


//=================== variables for scafContigLinking.c =========================
contigInfo *contigInfoArr;			// the contig information array, which is also used in scafContigMerge.c

arrLoc *arrLocArr;
int64_t itemNumArrLocArr;
int64_t maxItemNumArrLocArr;

ContigGraph *contigGraphArr;
int64_t itemNumContigGraphArr;

double maxLinksNumContigsThres;		// the maximal number of links between contigs , added 2012-11-24
double minLinksNumContigsThres;		// the minimal number of links between contigs
double maxRatioSecondFirstLinkNum;	// the maximal ratio of the second link number to the maximal link number between contigs
double secondLinkNumFactor;			// the second link number factor

double averLinkNum;
double maxSecondLinkNumThres;
double minBreakNumThres;


int itemNumContigLinkArr;			// item number of contig link array
int headRowContigLinkArr;			// head row of contig link array
int tailRowContigLinkArr;			// tail row of contig link array
contigLink *contigLinkArr;			// the contig link array

//================ New Begin ===================
linkCandidate *linkCandidateArr;
int64_t itemNumLinkCandidateArr;
int64_t maxItemNumLinkCandidateArr;
//================ New End ===================

//=================== variables for scafContigOverlap.c =========================
//int gapFillingFlag;				// YES or NO: YES-- the gap will be filled, No-- the gap will not be filled.

int maxOverlapSeqLen;						// the maximal length of sequence in contig ends to be overlapped
char *overlapSeq1;							// the first sequence for overlapping
char *overlapSeq2;							// the second sequence for overlapping
int *scoreArr;								// the score array for pairwise alignment
char *alignResultArr[3];					// the alignment result array

//double meanSizeInsert;						// the mean size of the insert fragments
//double stardardDeviationInsert;				// the standard deviation of the insert fragments
double standardDeviationFactor;				// the factor of the standard deviation of the insert fragments

int minOverlapThres;						// the minimal overlap threshold
int minExactOverlapThres;					// the minimal exact overlap threshold
int mismatchThres;							// the mismatch threshold
double maxMisMatchRatioThres;				// the maximal reatio threshold of mismatched bases
int subMismatchThres;						// the sub mismatch threshold
int minAdjustGapSizeThres;					// the minimal adjusted gap size threshold
int maxAdjustGapSizeThres;					// the maximal adjusted gap size threshold
int minBaseNumInGap;						// the minimal bases in gap region
double gapSizeSdevFactorOverlap;			// the standard deviation factor for gap size in contig overlap
double gapSizeSdevFactorGapFilling;			// the standard deviation factor for gap size in gap filling
//double exactOverlapSdevThres;				// the exact overlap standard deviation range threshold

int breakLinkNumThres;						// the break links number threshold

int matchScore;								// the match score
int mismatchScore;							// the mismatch score
int gapScore;								// the gap score

int scaffoldsNumInCOI;						// total number of scaffolds
contigOverlapIndex *contigOverlapIndexArr;	// pointer to contig overlap index array

int itemNumInContigOverlapArr;
contigOverlap *contigOverlapInfoArr;		// pointer to contig overlap information array


//=================== variables for scafGraph.c =========================
int64_t readItemNumInMRL;					// the number of reads in mono read list (SRL)
int64_t matchItemNumInMRP;					// the number of match items in mono read position array(SRP)
ReadList *monoReadListArr;					// pointer to the mono read List array
ReadPos *monoReadPosArr;					// pointer to the mono read position array

scafGraph *scafGrapDeBruijn;				// the de Bruijn graph in scaffolding


//=================== variables for scafGgpFilling.c =========================
scafKmer *scafKmers[2];						// scafKmers: [0] -- scafKmer with plus orientation; [1] -- scafKmer with minus orientation
int kmer_len;

int maxScafAssemblingReadsNumInArr;			// the maximal number of reads in scafAssemblingReadArr
int maxScafSuccessReadsNumInArr;			// the maximal number of successful reads in scafSuccessReadArr
//int maxScafContigEndSeqLen;					// the maximal length of scafContig end sequence
//int maxComparisonSeqLenInScaf;				// the maximal comparison sequence length in scaffolding

int successFilledFlag;						// YES or NO

scafAssemblingRead *scafAssemblingReadArr;	// the decision table
int scafAssemblingReadsNum;					// reads number in scafAssemblingReadArr
scafSuccessRead *scafSuccessReadArr;		// the successful reads array in scaffolding
int scafSuccessReadsNum;					// the successful reads number in scaffolding

scafContig *scafContighead;					// pointer to the head of scafContigs
scafContig *scafContigtail;  				// pointer to the tail of scafContigs
scafContig *scafContigLastReadLen;			// pointer to the scafContig node having READ_LEN nodes before scafContig tail
int scafContigIndex;
scafContig *successScafContig;				// pointer to the most recent successful scafContig node
char *scafContigSeqLastReadLen;				// the last base sequence with length READ_LEN
int scafContigSeqLenLastReadLen;			// the base number of last base sequence

char *scafContigEndSeqArr[2];					// the two contig end sequences with the maximal length of READ_LEN
int scafContigEndSeqLenArr[2];					// the sequences length of the two contig end sequences with the maximal length of READ_LEN
char *comparisonSeqInScaf;					// the comaprison sequence in scaffolding
int comparisonSeqLenInScaf;					// the comparison sequence length

int prepareAssemblyLenArr[2];				// the assembly length of prepared assembly in contig ends

int lockedReadsNumInScaf;					// the number of locked reads in scaffolding


//++++++++++++++++++++++++++++++++++++
int PEGivenType, oldPEGivenType;

int occsNumSE[4];  //下一个kmer的得分，对应ACGT的顺序
int occsNumPE[4];

int longKmerSize;
int longKmerStepSize;
double averKmerOcc;
double firstKmerThres;
double minKmerOccSE;
double minKmerOccPE;
double maxSecondOcc;
double maxFirstOcc;
double minLongKmerOcc;
double lockedReadsNumThres;
double maxOccNumFaiedPE;

int errorRegLenEnd3;				// erroneous region length of 3' end of a read

//++++++++++++++++++++++++++++++++++++
int kmerSize;
int entriesPerKmer;
uint64_t lastEntryMaskKmer;
int lastEntryBaseNumKmer;		// the base number of the last entry of kmer sequence array

uint64_t hashTableSize;

char kmerBaseSeq[1000];

uint64_t *kmerSeqInt;
uint64_t *kmerSeqIntRev;

uint64_t *kmerSeqIntAssembly;
uint64_t *kmerSeqIntAssemblyRev;
uint64_t *tmpKmerSeqIntAssembly;

//scafKmer *firstKmer;
uint64_t tabooSeqInt[4];

int navigationID;			// 0 -- Single End, 1 -- Paired End, -1 -- unknown
int navigationNumSE;
int maxNavigationNumSE;


#endif /* SCAFGLOBAL_H_ */
