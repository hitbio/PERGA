/*
 * inc/extvab.h
 *
 *  Created on: May 16, 2011
 *      Author: zhuxiao
 */

#ifndef EXTVAB_H_
#define EXTVAB_H_

#include "constants.h"


//=================== Global files or paths ==================
extern char inputPathStr[256];
extern char outputPathStr[256];
extern char outputPrefixStr[256];

extern char *readFilesInput[256];
extern int32_t readFileNum;
extern int32_t readsFileFormatType;
extern char contigsFileFasta[256];
extern char contigsFileHanging[256];

extern char fragmentSizeFile[256];
extern char graphFile[256];
extern char sampleContigsFile[256];
extern char readMatchInfoFile[256];
extern char scafSeqFile[256];

extern int32_t operationMode;
extern int32_t kmerSize;
extern int32_t readLenCutOff;					// the read length after cutting at the 3' end of reads
extern int32_t pairedMode;
extern double meanSizeInsert, standardDev;
extern int32_t minContigLen;
extern int32_t contigAlignRegSize;
extern int32_t gapFillFlag;


//*********************** reads.c **************************
extern int64_t totalReadNumSample;
extern int32_t maxReadLenInFileSample;				// the maximal length of read length in file while sampling
extern int32_t minReadLenInFileSample;				// the minimal length of read length in file while sampling
extern int32_t averReadLenInFileSample;				// the average length of read length in file while sampling

extern int32_t entriesPerReadseq;
extern uint64_t lastEntryMaskReadseq;
extern int32_t lastEntryBaseNumReadseq;			// the base number of the last entry of read sequence array
extern int32_t bytesPerReadseq;
extern uint64_t hashTableSizeReadseq;

extern int32_t maxUnknownBaseNumPerRead;
extern int32_t minSuccessiveAppearedBaseNum, maxSuccessiveAppearedBaseNum;

extern float kmerRegLenRatioEnd5, kmerRegLenRatioEnd3;

extern readSet_t *readSet;
extern int32_t reserveHashItemBlocksFlag;

extern readBlock_t *pReadBlockTmp;
extern read_t *pReadTmpDoing;

extern readseqBlock_t *pReadseqBlockTmp;
extern uint64_t *pReadseqTmpDoing;

extern readseqHashItemBlock_t *pReadseqHashItemBlockTmp;
extern readseqHashItem_t *pReadseqHashItemTmpDoing;


//**************** for graph.c *******************
extern graphtype *deBruijnGraph;

extern double timeuse_deBruijn;

extern uint64_t *pKmerSeqTmpDone, *pKmerSeqTmpDoing;
extern kmerseqBlock_t *pKmerseqBlockTmp;

extern kmertype *pKmerTmpDoing;
extern kmerBlock_t *pKmerBlockTmp;


//**************** for contig.c *******************
extern int32_t singleBaseQualThres;			// single base quality threshold

extern double readSimilarityThres;
extern int32_t minMatchNumSuccessRead;

extern int32_t longKmerSize;
extern int32_t longKmerStepSize;
extern double averKmerOcc;
extern double firstKmerThres;
extern double minKmerOccSE;
extern double minKmerOccPE;
extern double maxSecondOcc;
extern double maxFirstOcc;
extern double minLongKmerOcc;
extern double lockedReadsNumThres;
extern double minReadsNumPEHashThres;
extern int32_t navigationFlag;

extern int32_t maxUnmatchBaseNumPerRead;

extern int32_t kmer_len;
extern int32_t occsNumSE[4], occsNumPE[4];  // corresponds to A, C, G, T, respectively
extern int32_t occsNumIndexSE[4], occsNumIndexPE[4];
extern double svmFeatureArr[1000];

extern int32_t naviSuccessFlag;
extern int32_t naviTandFlag, newCandBaseNumAfterTandPathPE, newCandBaseNumAfterTandPathSE;

extern int32_t maxOccIndexSE, maxOccIndexPE, secondOccIndexSE, secondOccIndexPE;
extern double maxOccSE, maxOccPE, secondOccSE, secondOccPE;
extern double sumSecondOccPE, sumSecondOccSE;

extern short assemblyRound;		//FRIST_ROUND_ASSEMBLY  or SECOND_ROUND_ASSEMBLY
extern int32_t lockedReadsNum;
extern kmertype *kmers[2];		//[0]- plus k-mer, [1]- reverse complement k-mer
extern double lowerBoundFirstKmer, upperBoundFirstKmer;

extern contigtype *contigArr;
extern int64_t itemNumContigArr;
extern int64_t maxItemNumContigArr;
//extern int64_t validHeadRowContigArr, validTailRowContigArr;

extern int64_t successContigIndex;
extern int32_t *countingHangingBucketArr;

extern assemblingreadtype *decisionTable;		// the decision table
extern int32_t itemNumDecisionTable;				// item number in decision table
extern int32_t maxItemNumDecisionTable;				// the maximal number of reads in decision table

extern dtRowIndex_t **dtRowHashtable;			// the hash table to retrieve the reads in DT quickly, added 2012-12-02
//extern int maxitemNumInDTRowHashtable;

extern int64_t localContigID;
extern int32_t contigsNum;			// total contig count
extern int64_t basesNum;		// total base count of all contigs
extern int64_t this_successReadNum;
extern uint64_t kmerIndex;  	// start row of the first k-mer of a contig
extern int32_t turnContigIndex;

extern int32_t itemNumSuccessReadsArr;
extern int32_t maxItemNumSuccessReadsArr;
extern successRead_t *successReadsArr;

extern int32_t hangingContigOutFlag;		// whether output the hanging contig file: YES / NO (default)

extern int32_t number_of_overlap_less_than_threshold;
extern int64_t successReadNum;


//++++++++++++++++++++++++++++++++++++
extern int32_t PEGivenType, oldPEGivenType;
extern int32_t estimateSuccessFlag;;
extern PERead_t **PEHashArr;
extern int32_t readsNumInPEHashArr;
extern int32_t regLenPEHash, maxRegLenPEHash, minRegLenUsingPE, minMarginLenPEHash, maxMarginLenPEHash;// int leftMarginLenPEHash;
extern int32_t minContigLenUsingPE, shiftLenRound1, validReadOrientPEHash;
//extern contigtype *hashRegLeftContig, *hashRegRightContig, *shiftedRegLeftContig, *shiftedRegRightContig;  // ================
extern int32_t leftContigRowHashReg, rightContigRowHashReg, leftContigRowShiftedReg, rightContigRowShiftedReg;
extern int32_t minContigLenCheckGap;

extern int32_t shortInsertFlag;
extern double oldMeanSizeInsert, oldStandardDev, standardDevFactor, meanSizeInsertEst, standardDevEst;
extern int32_t contigNumEstContigArr;
extern int32_t minContigLenEst;			// the minimal contig length that can be used to estimate the insert size and standard deviation

extern readPosTemp_t *readPosTmpArr, *readPosTmpArrBuf;
extern int64_t readsNumSingleContig, totalReadsNumContigs;

extern readList_t *readListArr;
extern readPos_t *readPosArr;
extern int64_t itemNumInReadListArr, itemNumInReadPosArr;



//++++++++++++++++++++++++++++++++++++
extern int32_t readLen;
//extern int kmerSize;
extern int32_t entriesPerKmer;
extern uint64_t lastEntryMask;
extern int32_t lastEntryBaseNum;		// the base number of the last entry of its sequence array
extern int32_t bytesPerKmerseq;

extern uint64_t hashTableSize;

extern char baseSeq[MAX_READ_LEN_IN_BUF+1];

extern uint64_t *kmerSeqInt;
extern uint64_t *kmerSeqIntRev;

extern uint64_t *kmerSeqIntAssembly;
extern uint64_t *kmerSeqIntAssemblyRev;
extern uint64_t *tmpKmerSeqIntAssembly;

extern kmertype *firstKmer;
extern uint64_t tabooSeqInt[4];


//++++++++++++++++++++++++++++++++++++
extern int32_t regLenReadsNumReg, maxRegLenReadsNumReg;
extern int32_t minContigLenCheckingReadsNum;
//extern contigtype *leftContigReadsNumReg, *rightContigReadsNumReg; //================
extern int32_t leftContigRowReadsNumReg, rightContigRowReadsNumReg;
extern int32_t readsNumTotal, readsNumReadsNumReg;
extern double readsNumRatio;
extern double maxReadsNumRatioThres, minReadsNumRatioThres;
extern int32_t solvedRepeatsNum;

extern int32_t errNumArr[11];
extern int32_t totalErrReadNum;

//=============================================
// SVM variables
extern char svmKernelFuncFile[2][256], svmSupportVectorFile[2][256], svmAlphaFile[2][256], svmBiasFile[2][256], svmScaleDataFile[2][256];
extern svmModel_t *svmModel, *svmModelSE, *svmModelPE;
extern svmSampleVector_t *svmSample, *svmSampleSE, *svmSamplePE;


//=============================================
//Multi-align variables
extern int32_t naviBeforeCandPathPE, naviBeforeCandPathSE, naviAfterCandPathPE, naviAfterCandPathSE;
extern int32_t maxBaseIndexBeforeCandPathPE, maxBaseIndexBeforeCandPathSE;
extern int32_t correctNumCandPath, wrongNumCandPath;

extern int32_t naviBeforeTandPathPE, naviBeforeTandPathSE, naviAfterTandPathPE, naviAfterTandPathSE;
extern int32_t maxBaseIndexBeforeTandPathPE, maxBaseIndexBeforeTandPathSE;
extern int32_t correctNumTandPath, wrongNumTandPath;

extern int32_t incorrectBaseNumCandPathPE, incorrectBaseNumTandPathPE;

extern int32_t maxBaseIndexAfterCandPathPE, maxBaseIndexAfterCandPathSE, incorrectBaseNumCandPathSE;
extern int32_t maxBaseIndexAfterTandPathPE, maxBaseIndexAfterTandPathSE, incorrectBaseNumTandPathSE;

extern FILE *fpCandPathNum, *fpTandPathNum;

extern candPath_t *candPath;
extern contigPath_t *contigPath;

extern int32_t maxItemNumContigPath, maxItemNumContigPathAdjusted;

//=============================================
//contigGraph variables
extern contigGraph_t *contigGraph;

extern FILE *fpContigHang, *fpReadMatchInfo;


//=============================================
// scaffolding variables
extern scafContigIndex_t *scafContigIndex;

extern double maxLinksNumContigsThres;		// the maximal number of links between contigs , added 2012-11-24
extern double minLinksNumContigsThres;		// the minimal number of links between contigs
extern double maxRatioSecondFirstLinkNum;	// the maximal ratio of the second link number to the maximal link number between contigs
extern double secondLinkNumFactor;			// the second link number factor

extern double maxSecondLinkNumThres;
extern double minBreakNumThres;

extern contigLink_t *contigLinkSet;			// the contig link array

extern scaffoldSet_t *scaffoldSet;

//=================== variables for scafOverlap.c =========================
extern int32_t maxOverlapSeqLen;						// the maximal length of sequence in contig ends to be overlapped
extern char *overlapSeq1;							// the first sequence for overlapping
extern char *overlapSeq2;							// the second sequence for overlapping
extern int32_t *scoreArr;								// the score array for pairwise alignment
extern char *alignResultArr[3];					// the alignment result array

//double meanSizeInsert;						// the mean size of the insert fragments
//double stardardDeviationInsert;				// the standard deviation of the insert fragments
extern double standardDeviationFactor;				// the factor of the standard deviation of the insert fragments

extern int32_t minOverlapThres;						// the minimal overlap threshold
extern int32_t minExactOverlapThres;					// the minimal exact overlap threshold
extern int32_t mismatchThres;							// the mismatch threshold
extern double maxMisMatchRatioThres;				// the maximal reatio threshold of mismatched bases
extern int32_t subMismatchThres;						// the sub mismatch threshold
extern int32_t minAdjustGapSizeThres;					// the minimal adjusted gap size threshold
extern int32_t maxAdjustGapSizeThres;					// the maximal adjusted gap size threshold
extern int32_t minBaseNumInGap;						// the minimal bases in gap region
extern double gapSizeSdevFactorOverlap;			// the standard deviation factor for gap size in contig overlap
extern double gapSizeSdevFactorGapFilling;			// the standard deviation factor for gap size in gap filling
//double exactOverlapSdevThres;				// the exact overlap standard deviation range threshold

extern int32_t breakLinkNumThres;						// the break links number threshold

extern int32_t matchScore;								// the match score
extern int32_t mismatchScore;							// the mismatch score
extern int32_t gapScore;								// the gap score


//=================== variables for scafGap.c =========================
//extern char *scafContigSeqLastReadLen;				// the last base sequence with length READ_LEN
//extern int32_t scafContigSeqLenLastReadLen;			// the base number of last base sequence

extern char *scafContigEndSeqArr[2];					// the two contig end sequences with the maximal length of READ_LEN
extern int32_t scafContigEndSeqLenArr[2];					// the sequences length of the two contig end sequences with the maximal length of READ_LEN
extern char *comparisonSeqInScaf;					// the comaprison sequence in scaffolding
extern int32_t comparisonSeqLenInScaf;					// the comparison sequence length

extern int32_t prepareAssemblyLenArr[2];				// the assembly length of prepared assembly in contig ends




// ======================== variables for SVM ======================
extern double svmSV_PE[rowsNumSupportVector_PE*colsNumSupportVector_PE];
extern double svmAlpha_PE[rowsNumSupportVector_PE];
extern double svmBias_PE;
extern double svmScaleData_PE[2*colsNumSupportVector_PE];
extern double svmSV_SE[rowsNumSupportVector_SE*colsNumSupportVector_SE];
extern double svmAlpha_SE[rowsNumSupportVector_SE];
extern double svmBias_SE;
extern double svmScaleData_SE[2*colsNumSupportVector_SE];



#endif /* EXTVAB_H_ */
