/*
 * inc/global.h
 *
 *  Created on: May 9, 2011
 *      Author: zhuxiao
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_ 1

#include "svmcons.h"


//=================== Global files or paths ==================
char inputPathStr[256];
char outputPathStr[256];
char outputPrefixStr[256];

char *readFilesInput[256];
int32_t readFileNum;
int32_t readsFileFormatType;
char contigsFileFasta[256];
char contigsFileHanging[256];

char fragmentSizeFile[256];
char graphFile[256];
char sampleContigsFile[256];
char readMatchInfoFile[256];
char scafSeqFile[256];

int32_t operationMode;
int32_t kmerSize;
int32_t readLenCutOff;					// the read length after cutting at the 3' end of reads
int32_t pairedMode;
double meanSizeInsert, standardDev;
int32_t minContigLen;
int32_t contigAlignRegSize;
int32_t gapFillFlag;


//*********************** reads.c **************************
int64_t totalReadNumSample;
int32_t maxReadLenInFileSample;				// the maximal length of read length in file while sampling
int32_t minReadLenInFileSample;				// the minimal length of read length in file while sampling
int32_t averReadLenInFileSample;			// the average length of read length in file while sampling

int32_t entriesPerReadseq;
uint64_t lastEntryMaskReadseq;
int32_t lastEntryBaseNumReadseq;		// the base number of the last entry of read sequence array
int32_t bytesPerReadseq;
uint64_t hashTableSizeReadseq;

int32_t maxUnknownBaseNumPerRead;
int32_t minSuccessiveAppearedBaseNum, maxSuccessiveAppearedBaseNum;

float kmerRegLenRatioEnd5, kmerRegLenRatioEnd3;

readSet_t *readSet;
int32_t reserveHashItemBlocksFlag;

readBlock_t *pReadBlockTmp;
read_t *pReadTmpDoing;

readseqBlock_t *pReadseqBlockTmp;
uint64_t *pReadseqTmpDoing;

readseqHashItemBlock_t *pReadseqHashItemBlockTmp;
readseqHashItem_t *pReadseqHashItemTmpDoing;


//*********************** graph.c **************************
graphtype *deBruijnGraph;

double timeuse_deBruijn;


uint64_t *pKmerSeqTmpDone, *pKmerSeqTmpDoing;
kmerseqBlock_t *pKmerseqBlockTmp;

kmertype *pKmerTmpDoing;
kmerBlock_t *pKmerBlockTmp;


//*********************** contig.c **************************
contigtype *contigArr;
int64_t itemNumContigArr;
int64_t maxItemNumContigArr;
//int64_t validHeadRowContigArr, validTailRowContigArr;

int64_t successContigIndex;
int32_t *countingHangingBucketArr;

int32_t naviSuccessFlag;
int32_t naviTandFlag, newCandBaseNumAfterTandPathPE, newCandBaseNumAfterTandPathSE;

int32_t occsNumSE[4], occsNumPE[4];   // corresponds to A, C, G, T, respectively
int32_t occsNumIndexSE[4], occsNumIndexPE[4];
double svmFeatureArr[1000];

int32_t maxOccIndexSE, maxOccIndexPE, secondOccIndexSE, secondOccIndexPE;
double maxOccSE, maxOccPE, secondOccSE, secondOccPE;
double sumSecondOccPE, sumSecondOccSE;

int32_t singleBaseQualThres;			// single base quality threshold

double readSimilarityThres;
int32_t minMatchNumSuccessRead;

int32_t longKmerSize;
int32_t longKmerStepSize;
double averKmerOcc;
double firstKmerThres;
double minKmerOccSE;
double minKmerOccPE;
double maxSecondOcc;
double maxFirstOcc;
double minLongKmerOcc;
double lockedReadsNumThres;
double minReadsNumPEHashThres;
int32_t navigationFlag;

int32_t maxUnmatchBaseNumPerRead;

int32_t hangingContigOutFlag;		// whether output the hanging contig file: YES / NO (default)

short assemblyRound;		// FRIST_ROUND_ASSEMBLY  or SECOND_ROUND_ASSEMBLY
kmertype *kmers[2];			// [0]- plus k-mer, [1]- reverse complement k-mer
int32_t kmer_len;
double lowerBoundFirstKmer, upperBoundFirstKmer;

assemblingreadtype *decisionTable;	// the decision table
int32_t itemNumDecisionTable;			// item number in decision table
int32_t maxItemNumDecisionTable;		// the maximal number of reads in decision table

dtRowIndex_t **dtRowHashtable;		// the hash table to retrieve the reads in DT quickly, added 2012-12-02
//int maxitemNumInDTRowHashtable;

int64_t localContigID;
int32_t contigsNum;						// total contig count
int64_t basesNum;					// total base count of all contigs
int64_t this_successReadNum;
uint64_t kmerIndex;					// start row of the first k-mer of a contig
int32_t turnContigIndex;

int32_t itemNumSuccessReadsArr;
int32_t maxItemNumSuccessReadsArr;
successRead_t *successReadsArr;


int32_t lockedReadsNum = 0;

int64_t successReadNum = 0;			// total successful reads count
int32_t number_of_corrected_reads = 0;
int32_t number_of_overlap_less_than_threshold = 0;


//++++++++++++++++++++++++++++++++++++
int32_t PEGivenType, oldPEGivenType;
int32_t estimateSuccessFlag;
PERead_t **PEHashArr;
int32_t readsNumInPEHashArr;
int32_t regLenPEHash, maxRegLenPEHash, minRegLenUsingPE, minMarginLenPEHash, maxMarginLenPEHash;// int leftMarginLenPEHash;
int32_t minContigLenUsingPE, shiftLenRound1, validReadOrientPEHash;
//contigtype *hashRegLeftContig, *hashRegRightContig, *shiftedRegLeftContig, *shiftedRegRightContig;  // =================
int32_t leftContigRowHashReg, rightContigRowHashReg, leftContigRowShiftedReg, rightContigRowShiftedReg;
int32_t minContigLenCheckGap;

int32_t shortInsertFlag;
double oldMeanSizeInsert, oldStandardDev, standardDevFactor, meanSizeInsertEst, standardDevEst;
int32_t contigNumEstContigArr;
int32_t minContigLenEst;			// the minimal contig length that can be used to estimate the insert size and standard deviation

readPosTemp_t *readPosTmpArr, *readPosTmpArrBuf;
int64_t readsNumSingleContig, totalReadsNumContigs;

readList_t *readListArr;
readPos_t *readPosArr;
int64_t itemNumInReadListArr, itemNumInReadPosArr;



//++++++++++++++++++++++++++++++++++++
int32_t readLen;
//int kmerSize;
int32_t entriesPerKmer;
uint64_t lastEntryMask;
int32_t lastEntryBaseNum;		// the base number of the last entry of its sequence array
int32_t bytesPerKmerseq;
uint64_t hashTableSize;

char baseSeq[MAX_READ_LEN_IN_BUF+1];

uint64_t *kmerSeqInt;
uint64_t *kmerSeqIntRev;

uint64_t *kmerSeqIntAssembly;
uint64_t *kmerSeqIntAssemblyRev;
uint64_t *tmpKmerSeqIntAssembly;

kmertype *firstKmer;
uint64_t tabooSeqInt[4];


//++++++++++++++++++++++++++++++++++++
int32_t regLenReadsNumReg, maxRegLenReadsNumReg;
int32_t minContigLenCheckingReadsNum;
int32_t leftContigRowReadsNumReg, rightContigRowReadsNumReg;
int32_t readsNumTotal, readsNumReadsNumReg;
double readsNumRatio;
double maxReadsNumRatioThres, minReadsNumRatioThres;
int32_t solvedRepeatsNum;

int32_t errNumArr[11];
int32_t totalErrReadNum;

//=============================================
// SVM variables
char svmKernelFuncFile[2][256], svmSupportVectorFile[2][256], svmAlphaFile[2][256], svmBiasFile[2][256], svmScaleDataFile[2][256];
svmModel_t *svmModel, *svmModelSE, *svmModelPE;
svmSampleVector_t *svmSample, *svmSampleSE, *svmSamplePE;


//=============================================
//candPath variables
int32_t naviBeforeCandPathPE, naviBeforeCandPathSE, naviAfterCandPathPE, naviAfterCandPathSE;
int32_t maxBaseIndexBeforeCandPathPE, maxBaseIndexBeforeCandPathSE;
int32_t correctNumCandPath, wrongNumCandPath;

int32_t naviBeforeTandPathPE, naviBeforeTandPathSE, naviAfterTandPathPE, naviAfterTandPathSE;
int32_t maxBaseIndexBeforeTandPathPE, maxBaseIndexBeforeTandPathSE;
int32_t correctNumTandPath, wrongNumTandPath;

int32_t incorrectBaseNumCandPathPE, incorrectBaseNumTandPathPE;

int32_t maxBaseIndexAfterCandPathPE, maxBaseIndexAfterCandPathSE, incorrectBaseNumCandPathSE;
int32_t maxBaseIndexAfterTandPathPE, maxBaseIndexAfterTandPathSE, incorrectBaseNumTandPathSE;

FILE *fpCandPathNum, *fpTandPathNum;

candPath_t *candPath;
contigPath_t *contigPath;

int32_t maxItemNumContigPath = 0, maxItemNumContigPathAdjusted = 0;

//=============================================
//contigGraph variables
contigGraph_t *contigGraph;

FILE *fpContigHang, *fpReadMatchInfo;

//=============================================
// scaffolding variables
scafContigIndex_t *scafContigIndex;

double maxLinksNumContigsThres;		// the maximal number of links between contigs , added 2012-11-24
double minLinksNumContigsThres;		// the minimal number of links between contigs
double maxRatioSecondFirstLinkNum;	// the maximal ratio of the second link number to the maximal link number between contigs
double secondLinkNumFactor;			// the second link number factor

double maxSecondLinkNumThres;
double minBreakNumThres;

contigLink_t *contigLinkSet;			// the contig link array

scaffoldSet_t *scaffoldSet;

//=================== variables for scafOverlap.c =========================
int32_t maxOverlapSeqLen;						// the maximal length of sequence in contig ends to be overlapped
char *overlapSeq1;							// the first sequence for overlapping
char *overlapSeq2;							// the second sequence for overlapping
int32_t *scoreArr;								// the score array for pairwise alignment
char *alignResultArr[3];					// the alignment result array

//double meanSizeInsert;						// the mean size of the insert fragments
//double stardardDeviationInsert;				// the standard deviation of the insert fragments
double standardDeviationFactor;				// the factor of the standard deviation of the insert fragments

int32_t minOverlapThres;						// the minimal overlap threshold
int32_t minExactOverlapThres;					// the minimal exact overlap threshold
int32_t mismatchThres;							// the mismatch threshold
double maxMisMatchRatioThres;				// the maximal reatio threshold of mismatched bases
int32_t subMismatchThres;						// the sub mismatch threshold
int32_t minAdjustGapSizeThres;					// the minimal adjusted gap size threshold
int32_t maxAdjustGapSizeThres;					// the maximal adjusted gap size threshold
int32_t minBaseNumInGap;						// the minimal bases in gap region
double gapSizeSdevFactorOverlap;			// the standard deviation factor for gap size in contig overlap
double gapSizeSdevFactorGapFilling;			// the standard deviation factor for gap size in gap filling
//double exactOverlapSdevThres;				// the exact overlap standard deviation range threshold

int32_t breakLinkNumThres;						// the break links number threshold

int32_t matchScore;								// the match score
int32_t mismatchScore;							// the mismatch score
int32_t gapScore;								// the gap score


//=================== variables for scafGap.c =========================
char *scafContigEndSeqArr[2];					// the two contig end sequences with the maximal length of READ_LEN
int32_t scafContigEndSeqLenArr[2];					// the sequences length of the two contig end sequences with the maximal length of READ_LEN
char *comparisonSeqInScaf;					// the comaprison sequence in scaffolding
int32_t comparisonSeqLenInScaf;					// the comparison sequence length

int32_t prepareAssemblyLenArr[2];				// the assembly length of prepared assembly in contig ends


#endif /* GLOBAL_H_ */
