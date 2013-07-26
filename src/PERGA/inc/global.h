/*
 * inc/global.h
 *
 *  Created on: May 9, 2011
 *      Author: zhuxiao
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_ 1




//=================== Global files or paths ==================
char inputPathStr[256];
char outputPathStr[256];
char outputPrefixStr[256];

char *readFilesInput[256];
int readFileNum;
int readsFileFormatType;
char contigsFileFasta[256];
char contigsFileFastaCorrected[256];
char contigsFileHanging[256];

char fragmentSizeFile[256];
char graphFile[256];
char graphFileCorrected[256];
char sampleContigsFile[256];
char readCorrectedFile[256];

int operationMode;
int kmerSize;
int readLenCutOff;					// the read length after cutting at the 3' end of reads
int pairedMode;
int errorCorrectionFlag;
double meanSizeInsert, standardDev;
int minContigLen;
int trimReadLenFlag;


//*********************** reads.c **************************
int64_t totalReadNumSample;
int maxReadLenInFileSample;				// the maximal length of read length in file while sampling
int minReadLenInFileSample;				// the minimal length of read length in file while sampling
int averReadLenInFileSample;			// the average length of read length in file while sampling

int32_t entriesPerReadseq;
uint64_t lastEntryMaskReadseq;
int32_t lastEntryBaseNumReadseq;		// the base number of the last entry of read sequence array
int32_t bytesPerReadseq;
uint64_t hashTableSizeReadseq;

int32_t maxUnknownBaseNumPerRead;
int32_t minSuccessiveAppearedBaseNum, maxSuccessiveAppearedBaseNum;

float kmerRegLenRatioEnd5, kmerRegLenRatioEnd3;

readSet_t *readSet;
int reserveHashItemBlocksFlag;

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
int64_t validHeadRowContigArr, validTailRowContigArr;

int64_t successContigIndex;
int32_t *countingHangingBucketArr;

int32_t naviSuccessFlag;

int32_t occsNumSE[4], occsNumPE[4];   // corresponds to A, C, G, T, respectively
int32_t occsNumIndexSE[4], occsNumIndexPE[4];
double svmFeatureArr[1000];

int32_t maxOccIndexSE, maxOccIndexPE, secondOccIndexSE, secondOccIndexPE;
double maxOccSE, maxOccPE, secondOccSE, secondOccPE;
double sumSecondOccPE, sumSecondOccSE;

//int readLenInFile;					// the read length in the data file
//int maxReadLenInFile;				// the maximal length of read length in file
//int minReadLenInFile;				// the minimal length of read length in file
//int averReadLenInFile;				// the average length of read length in file
//int readLenCutOff;					// the read length after cutting at the 3' end of reads

int singleBaseQualThres;			// single base quality threshold

double readSimilarityThres;

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
double minReadsNumPEHashThres;
double maxOccNumFaiedPE;
int navigationFlag;

int32_t maxUnmatchBaseNumPerRead;
int32_t minUnmatchBaseNumAlign;

double *naviOccQueue;
int itemNumNaviOccQueue;
int maxItemNumNaviOccQueue;
int frontRowNaviOccQueue;
int rearRowNaviOccQueue;
double lowOccThresNaviOccQueue;

FILE *fpContigsBase, *fpContigsHanging;
int hangingContigOutFlag;		// whether output the hanging contig file: YES / NO (default)

short assemblyCycle;		// values: 1 for 0.5x <= firstKmerThres <= 15x, 2 for firstKmerThres > 15x, 3 for 2 <= firstKmerThres < 0.5x
short assemblyRound;		// FRIST_ROUND_ASSEMBLY  or SECOND_ROUND_ASSEMBLY
kmertype *kmers[2];			// [0]- plus k-mer, [1]- reverse complement k-mer
int kmer_len;
double lowerBoundFirstKmer, upperBoundFirstKmer;

assemblingreadtype *decisionTable;	// the decision table
int itemNumDecisionTable;			// item number in decision table
int maxItemNumDecisionTable;		// the maximal number of reads in decision table

dtRowIndex_t **dtRowHashtable;		// the hash table to retrieve the reads in DT quickly, added 2012-12-02
//int maxitemNumInDTRowHashtable;

int64_t localContigID;
int contigsNum;						// total contig count
int64_t basesNum;					// total base count of all contigs
int64_t this_successReadNum;
uint64_t kmerIndex;					// start row of the first k-mer of a contig

int itemNumSuccessReadsArr;
int maxItemNumSuccessReadsArr;
successRead_t *successReadsArr;


int lockedReadsNum = 0;

int64_t successReadNum = 0;			// total successful reads count
int number_of_corrected_reads = 0;
int number_of_overlap_less_than_threshold = 0;

int32_t *contigsLenArr;
int32_t itemNumContigsLenArr;
int32_t maxItemNumContigsLenArr;


//++++++++++++++++++++++++++++++++++++
int PEGivenType, oldPEGivenType;
int estimateSuccessFlag;
PERead_t **PEHashArr;
int readsNumInPEHashArr;
int allowedUpdatePEHashArrFlag;		// YES or NO
int regLenPEHash, maxRegLenPEHash, minRegLenUsingPE, minMarginLenPEHash, maxMarginLenPEHash;// int leftMarginLenPEHash;
int minContigLenUsingPE, shiftLenRound1, validReadOrientPEHash;
//contigtype *hashRegLeftContig, *hashRegRightContig, *shiftedRegLeftContig, *shiftedRegRightContig;  // =================
int32_t leftContigRowHashReg, rightContigRowHashReg, leftContigRowShiftedReg, rightContigRowShiftedReg;

double oldMeanSizeInsert, oldStandardDev, standardDevFactor, meanSizeInsertEst, standardDevEst;
int contigNumEstContigArr;
int minContigLenEst;			// the minimal contig length that can be used to estimate the insert size and standard deviation

readPosTemp_t *readPosTmpArr, *readPosTmpArrBuf;
int64_t readsNumSingleContig, totalReadsNumContigs;

readList_t *readListArr;
readPos_t *readPosArr;
int64_t itemNumInReadListArr, itemNumInReadPosArr;



//++++++++++++++++++++++++++++++++++++
int readLen;
//int kmerSize;
int entriesPerKmer;
uint64_t lastEntryMask;
int lastEntryBaseNum;		// the base number of the last entry of its sequence array
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

int navigationID;			// 0 -- Single End, 1 -- Paired End, -1 -- unknown
int navigationNumSE;
int maxNavigationNumSE;


//++++++++++++++++++++++++++++++++++++
int regLenReadsNumReg, maxRegLenReadsNumReg;
int minContigLenCheckingReadsNum;
//contigtype *leftContigReadsNumReg, *rightContigReadsNumReg;  // ===================
int32_t leftContigRowReadsNumReg, rightContigRowReadsNumReg;
int readsNumTotal, readsNumReadsNumReg;
double readsNumRatio;
double maxReadsNumRatioThres, minReadsNumRatioThres;
int solvedRepeatsNum;


//=================== correction.c ======================
int32_t maxSeqLenAlign;				// default: 2 * (maxReadLenInReadset + 1)
char *readSeqAlign;
int32_t readSeqLenAlign;
char *contigSeqAlign;
int32_t contigSeqLenAlign;
char *alignResultSeqArr[3];			// [0]- read bases, [1]- match letters, [2]- contig bases
int32_t alignResultSeqLen;

char *matchFlagArr;

int32_t *alignScoreArr;				// alignment score array, with size: alignScoreArr[maxSeqLenAlign+1][maxSeqLenAlign+1]

int32_t maxUnmatchBaseNumAfterAlign;
int32_t totalAlignedSuccessReadNum;
int32_t maxErrBaseNumInCorrection;

int32_t errNumArr[11];
int32_t totalErrReadNum;

int32_t correctionAllowed;
uint64_t *readseqCorrected;
int32_t readseqLenCorrected;
int64_t totalReadsNumCorreted;

FILE *fpReadCorrected;

int32_t matchScore;						// the match score
int32_t mismatchScore;					// the mismatch score
int32_t gapScore;						// the gap score


//============================================
// reference part
char refFile[256];
char *occPointFileArr1[4][6], *occPointFileArr2[4][6];
char occPointFile[256], occExtensionCorrectFile[4][256], occExtensionIncorrectFile[4][256];
char occStopCorrectFile[4][256], occStopIncorrectFile[4][256];
char occPointFile2[256], occExtensionCorrectFile2[4][256], occExtensionIncorrectFile2[4][256];
char occStopCorrectFile2[4][256], occStopIncorrectFile2[4][256];
char refPosFile[256], refPosFile2[256];
FILE *fpOccPoint, *fpOccExtensionCorrect[4], *fpOccExtensionIncorrect[4];
FILE *fpOccStopCorrect[4], *fpOccStopIncorrect[4];
FILE *fpRefPos;

ref_t *refArr;
int32_t itemNumRefArr;
int32_t maxItemNumRefArr;

int32_t refPosSoildFlag;
int32_t refStrandContig;
int32_t refPosContig;				// starts from 1
int32_t refBaseIntContig;
int32_t refPosMatchFlag;
int32_t refPosFirstBase[3];		// [0]-- solid flag; [1]-- strand; [2]-- refPos (>=1)


//=============================================
// SVM variables
char svmKernelFuncFile[2][256], svmSupportVectorFile[2][256], svmAlphaFile[2][256], svmBiasFile[2][256], svmScaleDataFile[2][256];
svmModel_t *svmModel, *svmModelSE, *svmModelPE;
svmSampleVector_t *svmSample, *svmSampleSE, *svmSamplePE;


//=============================================
//Multi-align variables
char **readseqMalignArr, **alignResultsMalign;
int64_t *readIDMalignArr;
int32_t maxItemNumReadsMalignArr;



#endif /* GLOBAL_H_ */
