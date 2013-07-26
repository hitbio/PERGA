/*
 * scafMethods.h
 *
 *  Created on: Mar 31, 2011
 *      Author: zhuxiao
 */

#ifndef SCAF_METHODS_H_
#define SCAF_METHODS_H_ 1

#include "scafStructure.h"


//######################## scafMain.c begin ########################//
short parseCommandParasAndExe(int argc, char **argv);
short showUsageInfo();
//######################## scafMain.c end ########################//


//######################## scafSRGA.c begin ########################//
short startSRGA_SF(char *contigsFilePara, int alignRegSizePara, int minContigLenPara, int pairedModePara, char **readFilesPara, int readFileNumPara, double meanSizeInsertPara, double standardDevPara, char *gapFillFlagPara, char *outputPathPara, char *outputPrefixPara);
short setGlobalParas(const char *outputPathName, const char *outputPrefixName, const char *contigFileName, int minContigLenThreshold, int alignRegSizeThreshold, char **readFilesPara, int readFileNumPara, int pairedModePara, double meanSizeInsertPara, double standardDevPara, char *gapFillFlagPara);
short setGlobalPath(const char *outPathStr);
void freeGlobalParas();
short setContigFileNames(const char *outputPathStr, const char *outputPrefixName, const char *contigFileName);
short getReadsFileFormatInScaf(int *readsFileFormatType, char **readFilesInput, int readFileNum);
short getReadLenFromFilesInScaf(int *readLen, int *averReadLen, char **readFilesInput, int readFileNum, int readsFileFormatType);
short getMinReadLenFromFastaFilesInScaf(int *readLenInFile, int *averReadLen, char **readFilesInput, int readFileNum);
short getMinReadLenFromFastqFilesInScaf(int *readLenInFile, int *averReadLen, char **readFilesInput, int readFileNum);
short getReadLenFromFastaInScaf(int *tmpReadLen, const char *fastaFile);
short getReadLenFromFastqInScaf(int *tmpReadLen, const char *fastqFile);
short startScaffolding();
//######################## scafSRGA.c end ########################//


//######################## util.c begin ########################//
short checkSortResult();
short checkContigIndex();
void testSeqSearch();
void checkSharedReadList();
short checkSortedContigReads(ContigList *pContigListArray, int contigItemNumInCL, ContigRead *pContigReadArray);
void outputContigReadsSingleContigEnd(ContigRead *pContigReadArray, int readsNum);
short checkArrLocArray(arrLoc *arrLocArray, int64_t itemNumArrLocArray);
short isSymmetricCartesianProductArr();
short saveCartesianProductArrayToFile(const char *cartesianProductFile, int *pCartesianProductArr, int rowNum);
short outputScoreArr(int *scoreArray, int rowsNum, int colsNum);
short outputScoreArrToFile(int *scoreArray, int rowsNum, int colsNum);

short outputGapSizeInScaf(contigOverlapIndex *pContigOverlapIndexArray, int scaffoldsNum, contigOverlap *pContigOverlapInfoArray, int minBaseNumInGap);
short checkScafKmersInGrapInScaf(const char *monoSeqFile, scafGraph *pScafGrapDeBruijn);
short outputKmerInScaf(uint64_t hashcode, uint64_t *seqInt, scafGraph *graph);
short outputContigInScaf(scafContig *contighead);
short outputScafContigSequence(scafContig *contighead);
void outputSuccessReads(scafSuccessRead *pScafSuccessReadArr, int tablesize);
short outputContigSeqRegionsInScaf(char *contigRegionFile, contigOverlapIndex *pContigOverlapIndexArray, int scaffoldsNum, contigOverlap *pContigOverlapInfoArray, contigInfo *pContigInfoArray);
short getSingleScaffoldSeqLen(int *scaffoldLen, contigOverlap *pContigOverlapInfo, int linkedContigsNum, int rowsNum, int minBaseNumInGap, contigInfo *pContigInfoArray);
//######################## util.c end ########################//


//######################## scafContig.c begin ########################//
short generateContigsFasta(char *outContigFileName, char *srcContigFileName);
short generateSingleContigFasta();
short outputSingleContigToFile(FILE *fpContigs, int64_t contigID, char *contigHeadTitle, char *contigSeq, int64_t contigLen);
short filterShortContigs(char *contigsFileFiltered, char *contigsFileUltraShort, char *contigsFile);
short outputSingleContigToFile(FILE *fpContigs, int64_t contigID, char *contigHeadTitle, char *contigSeq, int64_t contigLen);
//######################## scafContig.c end ########################//

//######################## scafContigIndex.c begin ########################//
short buildContigIndex(const char *contigFileName, const char *contigIndexFile);
short initMemory(const char *contigFileName);
int getContigsNum(const char *contigFileName);
short extractSequencesFromContigEnds();
short getMaxContigLenFromFile(int64_t *maxContigLen, const char *contigFileName);
short getSingleContigFastaFromFile(FILE *fp_contig, char *contigHeadTitle, char *contigSeq, int64_t *contigLength);
short fillBaseSeqIndex(int contigID, char *contigSeq, int contigLen);
short seqHash(uint64_t *hashArray, char *contigSeq, int start);
short seqSort();
short initSortMem();
void freeSortMem();
short fillKmerData();
short initRightShiftBitsNum();
short kmersSort();
short kmerSelectionSort(int *rowArr, int rowsNum);
short updateBaseSeqArr();
short convertToContigIndex();
short initConvertingMem();
void freeMemAfterConveting();
void freeMemContigIndex(uniqueSeqKmer **pUniqueSeqKmerArr, uint64_t **pUniqueSeqArr, seqRow **pSeqRowArr, int *uniqueSeqNum, contigMatchInfo **pContigMatchInfoArr, int *itemNumContigInfo);

int seqSearch(uint64_t *seq);

short saveContigIndexToFile(const char *contigIndexFile);
short loadContigIndex(const char *contigIndexFile, uniqueSeqKmer **pUniqueSeqKmerArr, uint64_t **pUniqueSeqArr, seqRow **pSeqRowArr, int *uniqueSeqNum, contigMatchInfo **pContigMatchInfoArr, int *itemNumContigInfo);
short initGlobalVariables();

//######################## scafContigIndex.c end ########################//


//######################## scafMapping.c begin ########################//
short mapPEs(const char *readSeqFile, const char *readMatchFile1, const char *readMatchFile2, char **readFilesInput, int readFileNum, int readsFileFormatType, int pairedMode, const char *contigIndexFile);
short mapPEFastaSeparate(const char *readSeqFile, const char *readMatchFile1, const char *readMatchFile2, char **readFilesInput, int readFileNum);
short mapPEFastqSeparate(const char *readSeqFile, const char *readMatchFile1, const char *readMatchFile2, char **readFilesInput, int readFileNum);
short mapPEFastaInterleaved(const char *readSeqFile, const char *readMatchFile1, const char *readMatchFile2, char **readFilesInput, int readFileNum);
short mapPEFastqInterleaved(const char *readSeqFile, const char *readMatchFile1, const char *readMatchFile2, char **readFilesInput, int readFileNum);
short initMemReadsBuf(readBuf_t **pBuf);
void freeMemReadsBuf(readBuf_t **pBuf);
short fillReadsToBufFasta(FILE *fpReads, readBuf_t *pBuf, uint64_t *readsNum, int readLenThreshold);
short fillReadsToBufFastq(FILE *fpReads, readBuf_t *pBuf, uint64_t *readsNum, int readLenThreshold);
short mapSingleRead(uint64_t readID, char *readSeq, int *matchFlag, FILE *fpPE_Result, uint64_t *matchedContigNum, uint64_t *matchedReadNum);
short getSingleReadFasta(FILE *fpPE, readBuf_t *pBuf, int readLenThreshold);
short getSingleReadFastq(FILE *fpPE, readBuf_t *pBuf, int readLenThreshold);
short getReverseComplements(uint64_t *reverse_readSeqInt, uint64_t *readSeqInt);
short getResultFileName(char *resultFileName, const char *readFileName);
//######################## scafMapping.c end ########################//

//######################## scafReadList.c begin ########################//
short buildReadLists(const char *readListFile1, const char *readListFile2, const char *sharedReadListFile, const char *matchResultFile1, const char *matchResultFile2);
short fillSingleReadList(const char *matchResultFile, int arrIndex);
short initMemReadList(const char *matchResultFile1, const char *matchResultFile2);
void freeMemReadList();
int64_t getMatchItemNum(const char *matchResultFile);
int64_t getFileSize(const char *filename);
short fillSharedReadList();
short fillMonoReadList();
int64_t getReadRowFromReadList(const int64_t readID, const ReadList *pReadListArr, const int64_t readItemNum);
short saveSharedReadListToFile(const char *sharedReadListFile);
short saveMonoReadListToFile(const char *monoReadListFile);
short saveReadListsToFile(const char *readListFile1, const char *readListFile2);
short loadSingleReadList(const char *sharedReadListFile, ReadList **pReadListArr, int64_t *readItemNum, ReadPos **pReadPosArr, int64_t *matchItemNum);
void freeSingleReadList(ReadList **pReadListArr, int64_t *readItemNum, ReadPos **pReadPosArr, int64_t *matchItemNum);
//######################## scafReadList.c end ########################//


//######################## scafContigList.c begin ########################//
short buildContigList(const char *contigListFile, const char *sharedReadListFile, const char *contigFile);
short initMemContigList();
void freeMemContigList();
short fillContigList();
int getReadOrientContigEnd(int64_t readID, int contigID, int contigEnd, ContigList *pContigListArr, int contigListSize, ContigRead *contigReadArr);
int getReadRowInContigReads(int64_t readID, ContigRead *pContigReads, int arraySize);
short sortContigReads();
short sortContigReadSingleContigEnd(ContigRead *pContigReadArr, int readsNum);
short selectionSortContigReads(int *pSortRowArray, int rowsNum, ContigRead *pContigReadArray);
short saveContigListToFile(const char *contigListFile);
short loadContigList(const char *contigListFile, ContigList **pContigListArr, int64_t *contigItemNum, ContigRead **pContigReadArr, int64_t *contigReadItemNum);
void freeContigList(ContigList **pContigListArr, int64_t *contigItemNum, ContigRead **pContigReadArr, int64_t *contigReadItemNum);
//######################## scafContigList.c end ########################//


//######################## scafContigLinking.c begin ########################//
short contigsLinking(const char *linkResultFile, const char *averLinkNumFile, const char *contigFileName, const char *sharedReadListFile, const char *contigListFile);
short initMemLinking(const char *contigFileName, const char *sharedReadListFile, const char *contigListFile);
void freeMemLinking();
short initContigInfoArray(const char *contigFileName);
void freeMemContigInfo(contigInfo **pContigInfoArr, int *contigsNum);
short initContigLinkArr(const int contigsNum);
void freeMemContigLinkArr();
short constructContigGraph();
short initMemArrLocArr();
void freeMemArrLocArr();
short fillArrLocArr();
short sortArrLocArr(arrLoc *pArrLocArray, int64_t itemNumArrLocArray);
short generateContigGraph(ContigGraph *pContigGraphArray, int64_t itemNumContigGraphArray, arrLoc *pArrLocArray, int64_t itemNumArrLocArray);
short generateContigGraphEdge(ContigGraph *pContigGraphArray, int64_t itemNumContigGraphArray, arrLoc *pArrLocArray, int64_t itemNumArrLocArray);
short getItemNumSingleRow(int64_t *itemNum, int64_t row, arrLoc *pArrLocArray, int64_t itemNumArrLocArray);
short FillContigGraphEdge(ContigGraph *pContigGraphArray, int64_t itemNumContigGraphArray, contigInfo *contigInfoArray, ContigList *contigListArray, int64_t contigItemNumInCLArray, ContigRead *contigReadArray, ReadList *sharedReadListArray, int64_t readItemNumInSRLArray, ReadPos *sharedReadPosArray);
short removeInvalidGraphEdge(ContigGraph *pContigGraphArray, int64_t itemNumContigGraphArray);
short initMemContigGraph(const int contigsNum);
void freeMemContigGraph();
short fillContigInfoArr(const char *contigFileName, contigInfo *pContigInfoArr);
short initParaLinking(const char *averLinkNumFile);
short computeAverPairsEachContigEdge(double *averageLinkNum, ContigGraph *contigGraphArray, int64_t itemNumContigGraphArray, int minLinksNumContigsThres);
short linkContigs(const char *linkResultFile);
short getMaxColsOfSingleRow(ContigGraph *pContigGraphArray, maxRowCol *pMaxRowColNode, contigInfo *contigInfoArray);
short getMaxRowsOfSingleCol(ContigGraph *pContigGraphArray, maxRowCol *pMaxRowColNode, contigInfo *contigInfoArray);
short changeMaxRowCol(maxRowCol *pMaxRowColNode, contigLink *contigLinkArray, int headRowContigLinkArray, contigInfo *contigInfoArray, int linkRound, int turnRoundFlag);
short getMaxLinksSituationArray(double *pArray, int *maxRowID, int *secondRowID, double *maxValue, double *secondValue);
short getFirstLinkedContigs(int *firstContigID, maxRowCol *maxRowColNode, contigInfo *contigInfoArray, int tmpContigsNum, ContigGraph *pContigGraphArray);
short isLinkSingleton(maxRowCol *pMaxRowColNode, ContigGraph *pContigGraphArray, int linkRound);
short fillSituationArray(ContigEdge *pEdgeNode, int contigID1, int contigID2, int endFlag1, int endFlag2, contigInfo *contigInfoArray, ContigList *contigListArray, ContigRead *contigReadArray, ReadList *sharedReadListArray, int64_t itemNumInSRL, ReadPos *sharedReadPosArray);
short addContigToContigLinkArr(int *linkStatus, contigLink *contigLinkArray, int *itemNumContigLinkArray, int *headRowContigLinkArray, int *tailRowContigLinkArray, contigInfo *contigInfoArray, maxRowCol *pMaxRowColNode, int newContigNum, int linkRound);
short markContigGraphEdge(ContigGraph *pContigGraphArray, maxRowCol *pMaxRowColNode);
void saveLinkResultToFile(FILE *fpLinkResult, int linkID);
void saveUnlinkedContigsToFile(FILE *fpLinkResult, int *startLinkID);
int getScaffoldsNum(const char *linkResultFile);
//######################## scafContigLinking.c end ########################//


//######################## scafContigOverlap.c begin ########################//
short overlapContigsInScaffolds(const char *contigOverlapInfoFile, const char *meanSdevFile, const char *newContigFile, const char *newSharedReadListFile, const char *newReadListFile1, const char *newReadListFile2, const char *newContigListFile, const char *linkResultFile, const char *contigFile, const char *sharedReadListFile, const char *readListFile1, const char *readListFile2, const char *contigListFile, const char *averLinkNumFile);
short initMemContigOverlap(const char *linkResultFile, const char *contigFile, const char *sharedReadListFile, const char *readListFile1, const char *readListFile2, const char *contigListFile, const char *averLinkNumFile);
void freeMemContigOverlap();
short initContigOverlapInfoArray(const char *linkResultFile);
int getItemNumContigOverlapInfo(const char *linkResultFile);
short getAverLinkNum(double *averageLinkNum, const char *averLinkNumFile);
short fillContigOverlapData(const char *linkResultFile);
short generateContigOverlapInfo();
short updateContigOverlapLen(contigOverlap *pContigOverlapInfo, contigInfo *pContigInfoArr);
short reverseSeq(char *seq, int seq_len);
short computeSeqOverlapLenByAlignment(const char *seq1, const int seqLen1, const char *seq2, const int seqLen2, int *scoreArray, char **pAlignResultArray, int *overlapLen, int *mismatchNum);
short computeSeqOverlapLenExact(int *overlapLen, const char *seq1, const int seqLen1, const char *seq2, const int seqLen2, int *scoreArray, int gapSize);
short updateContigs(contigInfo *pContigInfoArr1, contigInfo *pContigInfoArr2, const int contigOrient1, const int contigOrient2, char *seq1, const int originalSeqLen1, char *seq2, const int originalSeqLen2);
short adjustOverlapSeq(char *seq1, char *seq2, char **pAlignResultArray, int *overlapLen);
short meanSizeSDEstimate();
short gapSizeEstimateBetweenContigs(contigOverlap *pContigOverlapInfo, int *gapSize, int *validPairedNum, int endCutRound, int *cutOrderArray, int *uncoveredEndLenArray);
short updateOverlapLenByCutUncoveredContigEnds(contigOverlap *pContigOverlapInfo, contigInfo *pContigInfoArr);
short getUncoveredLenAtContigEnds(contigOverlap *pContigOverlapInfo, int *pUncoveredEndLenArray);
short cutContigEnds(contigOverlap *pContigOverlapInfo, contigInfo *pContigInfoArr, int endCutRound, int *cutOrderArray, int *uncoveredEndLenArray);
short updateReadListsAndContigListSingleContigAfterAlignment(contigInfo *pContigInfoArr, int contigOrient, int contigIndex, int contigLenBeforeAlignment);
short updateReadListsAndContigListAfterCut(contigOverlap *pContigOverlapInfo, contigInfo *pContigInfoArr, int endCutRound, int *cutOrderArray, int *contigLenBeforeCut);
short removeReadsFromReadLists(uint64_t readID, int contigIDFrom, int contigPosFrom);
short updateReadsInReadLists(uint64_t readID, int contigIDFrom, int contigPosFrom, int contigLenDecrease);
short outputContigOverlapInfoToFile(const char *contigOverlapInfoFile);
short outputMeanSizeAndSDev(const char *meanSdevFile);
short outputContigInfoArrToFile(const char *contigFileName, contigInfo *pContigInfoArray, int tmpContigsNum);
short rewriteReadListsAndContigListToFiles(const char *sharedReadListFile, const char *readListFile1, const char *readListFile2, const char *contigListFile);
short rewriteSingleReadListToFile(const char *readListFile, ReadList *pReadListArray, int64_t readItemNum, ReadPos *pReadPosArray, int64_t matchItemNum);
short rewriteContigListToFile(const char *contigListFile, ContigList *pContigListArray, int64_t contigItemNum, ContigRead *pContigReadArray, int64_t contigReadItemNum);
short loadMeanSizeAndSDev(const char *meanSdevFile, double *meanSizeFrag, double *sDevFrag);
short loadContigOverlapInfo(const char *contigOverlapInfoFile, contigOverlapIndex **pContigOverlapIndexArray, int *scaffoldsNum, contigOverlap **pContigOverlapInfoArray, int *itemNumInContigOverlapArray);
short getOverlapItemNumInOverlapInfoFile(const char *contigOverlapInfoFile);
short fillDataFromOverlapInfoFile(const char *contigOverlapInfoFile, contigOverlapIndex *pContigOverlapIndexArray, int scaffoldsNum, contigOverlap *pContigOverlapInfoArray, int itemNumInContigOverlapArray);
void freeContigOverlapInfo(contigOverlapIndex **pContigOverlapIndexArray, int *scaffoldsNum, contigOverlap **pContigOverlapInfoArray, int *itemNumInContigOverlapArray);
short updateContigEndsInfo(contigOverlap *pContigOverlapInfo, contigInfo *pContigInfoArray);
//######################## scafContigOverlap.c end ########################//

//######################## scafContigMerge.c begin ########################//
short generateScaffoldSequence(const char *scafSeqFile, const char *contigOverlapInfoFile, const char *newContigFile);
short initMemGeSeq(const char *contigOverlapInfoFile, const char *newContigFile);
void freeMemGeSeq();
short generateScafSeq(const char *scafSeqFile);
short getMaxContigLenInSingleScafflod(contigOverlapIndex *pContigOverlapIndex, int *maxContigLen);
short getApproximateSingleScafflodLen(contigOverlapIndex *pContigOverlapIndex, int *approximateScaffoldLen);
//######################## scafContigMerge.c end ########################//

//######################## scafGraph.c begin ########################//
short buildGraph(const char *graphFile, const char *monoReadSeqFile, const char *monoReadListFile, const char *readListFile1, const char *readListFile2, const char *readSeqFile);
short generateMonoReadList(const char *monoReadListFile, const char *readListFile1, const char *readListFile2);
short initMemGeMonoReadList(const char *readListFile1, const char *readListFile2);
void freeMemGeMonoReadList();
short getMonoReadSequence(const char *monoReadSeqFile, const char *monoReadListFile, const char *readSeqFile);
short getMonoReadSeqFromReadSeqFile(const char *monoReadSeqFile, const char *readSeqFile);
short containUnknownBase(char *seq);
short constructScafGraph(const char *graphFile, const char *monoReadSeqFile);
short getReadLenFromReadSeqFile(const char *monoReadSeqFile);
short countKmersInScaf(const char* pchFileName);
short initGraphInScaf(scafGraph **pScafGrapDeBruijn);
short countKmersInReadInScaf(char *seq);
short generateKmerSeqIntInScaf(uint64_t *seqInt, char *seq);
uint64_t kmerhashInScaf(uint64_t *seqInt);
short countSingleKmerInScaf(int hashcode, uint64_t *kmerSeqInt);
short addKmersInScaf(const char* pchFileName);
short addReadInScaf(char *seq, uint64_t rid);
short addKmerInScaf(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, uint16_t rpos);
char *getKmerBaseByIntInScaf(uint64_t *kmerSeqInt);
scafKmer *getKmerInScaf(uint64_t *kmerSeqInt, scafGraph *graph);
scafKmer *getKmerByHashInScaf(uint64_t hashvalue, uint64_t *kmerSeqInt, scafGraph *graph);
short identicalKmerSeqInScaf(uint64_t *kmerSeqInt1, uint64_t *kmerSeqInt2);
scafKmer *getReverseKmerInScaf(uint64_t *kmerSeqIntRev, uint64_t *kmerSeqInt, scafGraph *graph);
short delScafKmerByHashInScaf(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, uint32_t rpos, scafGraph *graph);
int32_t findDirectScafRidposIndexInScaf(uint64_t rid, uint32_t rpos, scafRidpos *ridpostable, uint32_t posNum);
int32_t findStartScafRidposIndexInScaf(uint64_t rid, scafRidpos *rid_pos_table, uint32_t posNum);
int32_t getExectScafRidposIndexInScaf(uint64_t rid, uint32_t rpos, int32_t startIndex, scafRidpos *ridpostable, uint32_t posNum);

short saveGraphInScaf(const char* gFileName);
short loadGraphInScaf(scafGraph **pScafGrapDeBruijn, const char* gFileName);
void freeGraphInScaf(scafGraph **pScafGrapDeBruijn);
//######################## scafGraph.c end ########################//

//######################## scafGapFilling.c begin ########################//
short gapFilling(const char *graphFile, const char *monoReadSeqFile, const char *monoReadListFile, const char *newContigOverlapInfoFile, const char *newContigFile, const char *readListFile1, const char *readListFile2, const char *readSeqFile, const char *meanSdevFile, const char *contigOverlapInfoFile, const char *contigFile);
short fillGaps(const char *newContigOverlapInfoFile, const char *newContigFile, const char *meanSdevFile, const char *contigOverlapInfoFile, const char *contigFile, const char *monoReadListFile, const char *graphFile);
short initGapFillingParas();
void freeGapFillingParas();
short initMemGapFilling(const char *contigOverlapInfoFile, const char *meanSdevFile, const char *contigFile, const char *monoReadListFile, const char *graphFile);
void freeMemGapFilling();
short computeAverKmerOcc(double *averKmerOcc, scafGraph *pScafGrapDeBruijn);
short localAssemblyInScaf();
short getScafContigEndSeqs(contigOverlap *pContigOverlapInfo);
short getComparisonSeqInScaf(int assemblyRound);
short prepareAssemblyInScaf(int contigID, int contigLen, int contigOrient, int assemblyRound);
short initScafContig(char *seq, int seq_len);
short initScafAssemblingReads(int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound);
short addAssemblingReadsInScaf(int rid, int rpos, int orientation, int matedFlag);
short reallocateDecisionTableInScaf();
short getMaxMinAssemblyLen(int *maxAssemblyLen, int *minAssemblyLen, int assemblyRound, int *localScafContigNodesNumArr, int gapSize);
short getNextKmerBySEInScaf(int contigNodesNum, int contigID, int contigLen, int contigOrient, int assemblyRound);
short validReadPairInScaf(uint64_t readID, int readOrient, int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound);
short computeLongKmerOccNumBySEInScaf(scafKmer *tmp_kmers[2], int *occNum, int length_k, int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound);
short computeKmerOccNumBySEInScaf(scafKmer *tmp_kmers[2], int *occNum, int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound);
short computeKmerOccNumUnlockedBySEInScaf(scafKmer *tmp_kmers[2], int *occNum, int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound);
short computeKmerOccNumLockedBySEInScaf(scafKmer *tmp_kmers[2], int *occNum);
int getProperIndexInScaf(scafAssemblingRead *assemblingread);
int getProperIndexLimitedInScaf(scafAssemblingRead *assemblingread, int limitLastpos);
short appendContigBaseInScaf(unsigned char base, int contigIndex);
void updateLockedReadsInScaf();
short delReadsFromGraphInScaf();
short delRemainedKmersInScaf(char *seq, uint64_t *tmp_kmerseq, uint64_t rid, uint32_t rpos, scafGraph *graph);
short addRidposToContigInScaf(int contigNodesNum);
scafContig *getSuccessContigInScaf(int contigIndex);
short updateContigtailnodesInScaf(int *scafContigIndex);
short updateScafContigEndSeqs(int assemblyRound, int scafContigNodesNum);
short detectOverlapsInScaf(int *successFilledFlag, int *overlapLen, int *newGapSize, int *breakFlag, int gapSize, int *localScafContigNodesNumArray, int assemblyRound);
short computeNewGapSizeInScaf(int *tmpGapSize, int gapSize, int *localScafContigNodesNumArray, int assemblyRound);
short updateContigOverlapInfoInScaf(contigOverlap *pContigOverlapInfo, int successFilledFlag, int overlapLen, int newGapSize, int breakFlag);
short updateContigInfoInScaf(contigOverlap *pContigOverlapInfo, scafContig **localScafContigheadArr, int *localScafContigNodesNumArr, int *oldEndSeqLenArr, int assemblyRound);
short updateSingleContigInfoInScaf(contigInfo *pContigInfo, int contigOrient, int contigIndex, scafContig *scafContighead, int scafContigNodesNum, char *endSeq, int endSeqLen, int oldEndSeqLen, int prepareAssemblyLen);
short updateOtherContigInfoInScaf(contigInfo *pContigInfo, int contigOrient, char *contigEndSeq, int endSeqLen, int oldEndSeqLen);
void releaseScafContigInScaf(scafContig *contighead);
//######################## scafGapFilling.c end ########################//

//######################## scafUpdate.c begin ########################//
short updateScafAssemblingReadsInScaf(int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound);
short updateAssemblingreadsStatusInScaf();
short updateSameReadStatusToFailedInScaf(uint64_t rid, int assemblingreadIndex);
short updateFinisedScafReadsInScaf();
short getScafSuccessReadsInScaf();
short reallocateScafSuccessReadsArrInScaf();
short removeFinisedReadsInScaf();
//######################## scafUpdate.c end ########################//


//######################## scafPEAssembly.c begin ########################//
short getNextKmerByMixInScaf(int contigNodesNum, int contigID, int contigLen, int contigOrient, int assemblyRound);
short getNextKmerByPEInScaf(int contigNodesNum, int contigID, int contigLen, int contigOrient, int assemblyRound);
short computeKmerOccNumByPEInScaf(scafKmer *tmp_kmers[2], int *occNum, int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound);
short computeKmerOccNumLockedByPEInScaf(scafKmer *tmp_kmers[2], int *occNum, int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound);
short computeKmerOccNumUnlockedByPEInScaf(scafKmer *tmp_kmers[2], int *occNum, int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound);
//######################## scafPEAssembly.c end ########################//



#endif /* SCAF_METHODS_H_ */
