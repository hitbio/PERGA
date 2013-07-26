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
extern int readFileNum;
extern int readsFileFormatType;
extern char contigsFileFasta[256];
extern char contigsFileFastaCorrected[256];
extern char contigsFileHanging[256];

extern char fragmentSizeFile[256];
extern char graphFile[256];
extern char graphFileCorrected[256];
extern char sampleContigsFile[256];
extern char readCorrectedFile[256];

extern int operationMode;
extern int kmerSize;
extern int readLenCutOff;					// the read length after cutting at the 3' end of reads
extern int pairedMode;
extern int errorCorrectionFlag;
extern double meanSizeInsert, standardDev;
extern int minContigLen;
extern int trimReadLenFlag;


//*********************** reads.c **************************
extern int64_t totalReadNumSample;
extern int maxReadLenInFileSample;				// the maximal length of read length in file while sampling
extern int minReadLenInFileSample;				// the minimal length of read length in file while sampling
extern int averReadLenInFileSample;				// the average length of read length in file while sampling

extern int32_t entriesPerReadseq;
extern uint64_t lastEntryMaskReadseq;
extern int32_t lastEntryBaseNumReadseq;			// the base number of the last entry of read sequence array
extern int32_t bytesPerReadseq;
extern uint64_t hashTableSizeReadseq;

extern int32_t maxUnknownBaseNumPerRead;
extern int32_t minSuccessiveAppearedBaseNum, maxSuccessiveAppearedBaseNum;

extern float kmerRegLenRatioEnd5, kmerRegLenRatioEnd3;

extern readSet_t *readSet;
extern int reserveHashItemBlocksFlag;

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
//extern int readLenInFile;
//extern int maxReadLenInFile;				// the maximal length of read length in file
//extern int minReadLenInFile;				// the minimal length of read length in file
//extern int averReadLenInFile;				// the average length of read length in file

extern int singleBaseQualThres;			// single base quality threshold

extern double readSimilarityThres;

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
extern double minReadsNumPEHashThres;
extern double maxOccNumFaiedPE;
extern int navigationFlag;

extern int32_t maxUnmatchBaseNumPerRead;
extern int32_t minUnmatchBaseNumAlign;

extern double *naviOccQueue;
extern int itemNumNaviOccQueue;
extern int maxItemNumNaviOccQueue;
extern int frontRowNaviOccQueue;
extern int rearRowNaviOccQueue;
extern double lowOccThresNaviOccQueue;

extern int kmer_len;
extern int occsNumSE[4], occsNumPE[4];  // corresponds to A, C, G, T, respectively
extern int occsNumIndexSE[4], occsNumIndexPE[4];
extern double svmFeatureArr[1000];

extern int naviSuccessFlag;

extern int maxOccIndexSE, maxOccIndexPE, secondOccIndexSE, secondOccIndexPE;
extern double maxOccSE, maxOccPE, secondOccSE, secondOccPE;
extern double sumSecondOccPE, sumSecondOccSE;

extern short assemblyCycle;		// values: 1 for 0.5x <= firstKmerThres <= 15x, 2 for firstKmerThres > 15x, 3 for 2 <= firstKmerThres < 0.5x
extern short assemblyRound;		//FRIST_ROUND_ASSEMBLY  or SECOND_ROUND_ASSEMBLY
extern int lockedReadsNum;
extern kmertype *kmers[2];		//[0]- plus k-mer, [1]- reverse complement k-mer
extern double lowerBoundFirstKmer, upperBoundFirstKmer;

extern contigtype *contigArr;
extern int64_t itemNumContigArr;
extern int64_t maxItemNumContigArr;
extern int64_t validHeadRowContigArr, validTailRowContigArr;

extern int64_t successContigIndex;
extern int32_t *countingHangingBucketArr;

extern assemblingreadtype *decisionTable;		// the decision table
extern int itemNumDecisionTable;				// item number in decision table
extern int maxItemNumDecisionTable;				// the maximal number of reads in decision table

extern dtRowIndex_t **dtRowHashtable;			// the hash table to retrieve the reads in DT quickly, added 2012-12-02
//extern int maxitemNumInDTRowHashtable;

extern int64_t localContigID;
extern int contigsNum;			// total contig count
extern int64_t basesNum;		// total base count of all contigs
extern int64_t this_successReadNum;
extern uint64_t kmerIndex;  	// start row of the first k-mer of a contig

extern int itemNumSuccessReadsArr;
extern int maxItemNumSuccessReadsArr;
extern successRead_t *successReadsArr;

extern FILE *fpContigsBase, *fpContigsHanging;
extern int hangingContigOutFlag;		// whether output the hanging contig file: YES / NO (default)

extern int number_of_overlap_less_than_threshold;
extern int64_t successReadNum;

extern int32_t *contigsLenArr;
extern int32_t itemNumContigsLenArr;
extern int32_t maxItemNumContigsLenArr;


//++++++++++++++++++++++++++++++++++++
extern int PEGivenType, oldPEGivenType;
extern int estimateSuccessFlag;;
extern PERead_t **PEHashArr;
extern int readsNumInPEHashArr;
extern int allowedUpdatePEHashArrFlag;
extern int regLenPEHash, maxRegLenPEHash, minRegLenUsingPE, minMarginLenPEHash, maxMarginLenPEHash;// int leftMarginLenPEHash;
extern int minContigLenUsingPE, shiftLenRound1, validReadOrientPEHash;
//extern contigtype *hashRegLeftContig, *hashRegRightContig, *shiftedRegLeftContig, *shiftedRegRightContig;  // ================
extern int32_t leftContigRowHashReg, rightContigRowHashReg, leftContigRowShiftedReg, rightContigRowShiftedReg;

extern double oldMeanSizeInsert, oldStandardDev, standardDevFactor, meanSizeInsertEst, standardDevEst;
extern int contigNumEstContigArr;
extern int minContigLenEst;			// the minimal contig length that can be used to estimate the insert size and standard deviation

extern readPosTemp_t *readPosTmpArr, *readPosTmpArrBuf;
extern int64_t readsNumSingleContig, totalReadsNumContigs;

extern readList_t *readListArr;
extern readPos_t *readPosArr;
extern int64_t itemNumInReadListArr, itemNumInReadPosArr;



//++++++++++++++++++++++++++++++++++++
extern int readLen;
//extern int kmerSize;
extern int entriesPerKmer;
extern uint64_t lastEntryMask;
extern int lastEntryBaseNum;		// the base number of the last entry of its sequence array
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
extern int navigationID;		// 0 -- Single End, 1 -- Paired End, -1 -- unknown
extern int navigationNumSE;
extern int maxNavigationNumSE;


//++++++++++++++++++++++++++++++++++++
extern int regLenReadsNumReg, maxRegLenReadsNumReg;
extern int minContigLenCheckingReadsNum;
//extern contigtype *leftContigReadsNumReg, *rightContigReadsNumReg; //================
extern int32_t leftContigRowReadsNumReg, rightContigRowReadsNumReg;
extern int readsNumTotal, readsNumReadsNumReg;
extern double readsNumRatio;
extern double maxReadsNumRatioThres, minReadsNumRatioThres;
extern int solvedRepeatsNum;


//=================== correction.c ======================
extern int32_t maxSeqLenAlign;				// default: 2 * (maxReadLenInReadset + 1)
extern char *readSeqAlign;
extern int32_t readSeqLenAlign;
extern char *contigSeqAlign;
extern int32_t contigSeqLenAlign;
extern char *alignResultSeqArr[3];			// [0]- read bases, [1]- match letters, [2]- contig bases
extern int32_t alignResultSeqLen;

extern char *matchFlagArr;

extern int32_t *alignScoreArr;				// alignment score array, with size: alignScoreArr[maxSeqLenAlign+1][maxSeqLenAlign+1]

extern int32_t maxUnmatchBaseNumAfterAlign;
extern int32_t totalAlignedSuccessReadNum;
extern int32_t maxErrBaseNumInCorrection;

extern int32_t errNumArr[11];
extern int32_t totalErrReadNum;

extern int32_t correctionAllowed;
extern uint64_t *readseqCorrected;
extern int32_t readseqLenCorrected;
extern int64_t totalReadsNumCorreted;

extern FILE *fpReadCorrected;

extern int32_t matchScore;					// the match score
extern int32_t mismatchScore;				// the mismatch score
extern int32_t gapScore;					// the gap score


//============================================
// reference part
extern char refFile[256];
extern char *occPointFileArr1[4][6], *occPointFileArr2[4][6];
extern char occPointFile[256], occExtensionCorrectFile[4][256], occExtensionIncorrectFile[4][256];
extern char occStopCorrectFile[4][256], occStopIncorrectFile[4][256];
extern char occPointFile2[256], occExtensionCorrectFile2[4][256], occExtensionIncorrectFile2[4][256];
extern char occStopCorrectFile2[4][256], occStopIncorrectFile2[4][256];
extern char refPosFile[256], refPosFile2[256];
extern FILE *fpOccPoint, *fpOccExtensionCorrect[4], *fpOccExtensionIncorrect[4];
extern FILE *fpOccStopCorrect[4], *fpOccStopIncorrect[4];
extern FILE *fpRefPos;

extern ref_t *refArr;
extern int32_t itemNumRefArr;
extern int32_t maxItemNumRefArr;

extern int32_t refPosSoildFlag;
extern int32_t refStrandContig;
extern int32_t refPosContig;				// starts from 1
extern int32_t refBaseIntContig;
extern int32_t refPosMatchFlag;
extern int32_t refPosFirstBase[3];		// [0]-- solid flag; [1]-- strand; [2]-- refPos (>=1)


//=============================================
// SVM variables
extern char svmKernelFuncFile[2][256], svmSupportVectorFile[2][256], svmAlphaFile[2][256], svmBiasFile[2][256], svmScaleDataFile[2][256];
extern svmModel_t *svmModel, *svmModelSE, *svmModelPE;
extern svmSampleVector_t *svmSample, *svmSampleSE, *svmSamplePE;


//=============================================
//Multi-align variables
extern char **readseqMalignArr, **alignResultsMalign;
extern int64_t *readIDMalignArr;
extern int32_t maxItemNumReadsMalignArr;


#endif /* EXTVAB_H_ */
