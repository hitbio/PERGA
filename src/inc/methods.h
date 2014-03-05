#ifndef METHODS_H_INCLUDED
#define METHODS_H_INCLUDED

#include "structure.h"

//================= pergaMain.c declaration begin ================/
short parseCommandParasAndExe(int argc, char **argv);
short showUsageInfo();
//================= pergaMain.c declaration end ================/

//================= perga.c declaration begin ================/
short startPERGA(int operationModePara, int kmerSizePara, int readLenCutOffPara, int pairedModePara, char **readFilesPara, int readFileNumPara, int singleBaseQualThresPara, char *graphFilePara, char *contigsFilePara, char *readMatchInfoFilePara, double meanSizeInsertPara, double standardDevPara, char *outputPathPara, char *outputPrefixPara, int minContigLenPara, int contigAlignRegLenPara, char *gapFillFlagPara);
short initGlobalParas(int operationModepara, char *outputPathName, char *prefix, char **readFilesPara, int readFileNumPara, int singleBaseQualThresPara, int pairedModePara, int kmerLen, int readLenCut, char *graphFilePara, char *contigsFilePara, char *readMatchInfoFilePara, double meanSizeInsertPara, double standardDevPara, int minContigLenPara, int contigAlignRegLenPara, char *gapFillFlagPara);
short setGlobalPath(const char *outPathStr);
void freeGlobalParas();
short getReadsFileFormat(int *readsFileFormatType, char **readFilesInput, int readFileNum);
short getReadLenAndCountFromFilesFasta(int64_t *tmpTotalReadNum, int *tmpAverReadLenInFile, int *tmpMaxReadLenInFile, int *tmpMinReadLenInFile, char **readsFileNames, int readsFileNum);
short getReadLenAndCountFromFilesFastq(int64_t *tmpTotalReadNum, int *tmpAverReadLenInFile, int *tmpMaxReadLenInFile, int *tmpMinReadLenInFile, char **readsFileNames, int readsFileNum);
short getReadLenAndCountFromSingleFasta(int64_t *tmpReadCount, int64_t *tmpSumReadLen, int64_t *tmpMaxReadLen, int64_t *tmpMinReadLen, char *fastaFile);
short getReadLenAndCountFromSingleFastq(int64_t *tmpReadCount, int64_t *tmpSumReadLen, int64_t *tmpMaxReadLen, int64_t *tmpMinReadLen, char *fastqFile);
int32_t getSysMemorySize();
//================= perga.c declaration end ================/


//================= util.c declaration begin ================/
short outputContig(contigtype *contigArray, int64_t contigNodesNum);
short outputContigEnd3(contigtype *contigArray, int64_t contigNodesNum, int endNum);
short outputContigBaseEnd3(contigtype *contigArray, int64_t contigNodesNum, int endNum);
short outputContigEnd5(contigtype *contigArray, int64_t contigNodesNum, int endNum);
short outputUndelKmerpos(graphtype *graph);
short outputRemainedKmers(graphtype *graph);
short outputReadsInDecisionTable(assemblingreadtype *decisionTable, int readsNum);
short outputReadsInDecisionTableToFile(char *outfile, assemblingreadtype *decisionTable, int readsNum);
short outputMatedReadsInDecisionTableToFile(char *outfile, assemblingreadtype *decisionTable, int readsNum);
short outputFailedReadsInDecisionTable(assemblingreadtype *decisionTable, int itemNumDecisionTable, int contigID, int contigNodesNum);
short outputLockedReadsInDecisionTable(assemblingreadtype *decisionTable, int readsNum);
short outputMatedReadsInDecisionTable(assemblingreadtype *decisionTable, int readsNum);
short outputKmer(graphtype *graph, int hashcode, uint64_t *kmerSeqInt); //输出kmer中的内容
short outputRidpos(ridpostype *ridpos, int posNum);
void outputSuccessReads(successRead_t *successReadArray, int32_t successReadNum);
short checkGraph(graphtype *graph);
short outputContigToTmpFile(char *fileOut, contigtype *contigArray, int64_t contigNodesNum, int outFileType);
short outputPEHashArray(PERead_t **PEHashArray);
short checkReadListArr(readList_t *readListArray, int64_t itemNumInReadListArray);
short outputReadListToFile(char *readListFile, readList_t *readListArray, readPos_t *readPosArray, int64_t itemNumInReadListArray);
short outputMatedReadsInReadListToFile(char *readListFile, readList_t *readListArray, readPos_t *readPosArray, int64_t itemNumInReadListArray);
short convertFragmentSizeFileInText(const char *fragmentSizeFile);
short outputReadPosInGraph(int64_t readID, graphtype *graph);
short outputNaviOccQueue(double *naviOccQueuePara, int itemNumNaviOccQueuePara, int frontRowNaviOccQueuePara);
short outputDtRowHashtable(dtRowIndex_t **pDtRowHashtable);
short visitKmers(graphtype *graph);
short outputReadseqInReadset(char *outfile, readSet_t *readSet);
//================= util.c declaration end ================/

//================= reads.c declaration begin ================/
short constructReadset(readSet_t **readSet, char **readsFileNames, int readsFileNum, int reserveHashItemBlocksFlag);
short constructReadsetBySEFasta(readSet_t **readSet, char **readsFileNames, int readsFileNum);
short constructReadsetByPEFastaSeparate(readSet_t **readSet, char **readsFileNames, int readsFileNum);
short constructReadsetByPEFastaInterleaved(readSet_t **readSet, char **readsFileNames, int readsFileNum);
short constructReadsetBySEFastq(readSet_t **readSet, char **readsFileNames, int readsFileNum);
short constructReadsetByPEFastqSeparate(readSet_t **readSet, char **readsFileNames, int readsFileNum);
short constructReadsetByPEFastqInterleaved(readSet_t **readSet, char **readsFileNames, int readsFileNum);
short addReadToReadset(char *seq_data, char *qual_data, int32_t seqLen, readSet_t *readSet);
short addReadToReadsetWithRefpos(char *seq, char *qual_data, int32_t seqLen, char *headname, int32_t headlen, readSet_t *readSet);
short initReadSet(readSet_t **pReadSet);
short initReadBlockInReadset(readSet_t *pReadSet);
short addNewBlockRead(readSet_t *pReadSet);
short initReadseqBlockInReadset(readSet_t *pReadSet);
short initReadMatchInfoBlockInReadset(readSet_t *pReadSet);
short addNewBlockReadseq(readSet_t *pReadSet);
short initReadseqHashtableInReadset(readSet_t *pReadSet);
short initReadseqHashItemBlockInGraph(readSet_t *pReadSet);
short addNewBlockReadseqHashItem(readSet_t *pReadSet);
short releaseReadset(readSet_t **pReadSet);
void releaseHashItemReadset(readSet_t *readSet);
uint64_t readseqHashInt(uint64_t *seqInt, int32_t baseNum, int32_t entriesNum);
short generateReadseqInt(uint64_t *seqInt, char *seq, int32_t seqLen, int32_t entriesNum);
inline readseqHashItem_t *getReadseqHashItemByHash(uint64_t hashvalue, uint64_t *readseqInt, int32_t seqLen, int32_t entriesNum, readSet_t *readSet);
short identicalReadseq(uint64_t *kmerSeqInt1, uint64_t *kmerSeqInt2, int32_t entriesNum);
char *getReadBaseByInt(uint64_t *readseqInt, int32_t seqLen);
char *getReadBaseFromPosByInt(char *readBaseSeq, uint64_t *readseqInt, int32_t seqLen, int32_t startBasePos, int32_t baseNum);
char *getReverseReadBaseFromPosByInt(char *readBaseSeq, uint64_t *readseqInt, int32_t seqLen, int32_t startRevBasePos, int32_t baseNum);
char *getReverseReadBaseByInt(uint64_t *readseqInt, int32_t seqLen);
short getReverseReadseqInt(uint64_t *readseqIntRev, uint64_t *readseqInt, int32_t seqLen);
short replaceUnknownBasesInReads(char *seq, int32_t nBaseNum);
short computeMaxReadLenInReadset(readSet_t *readSet);
short extractRefPosFromHeadName(int32_t *strandPosArray, char *headname, int32_t headlen);
//================= reads.c declaration end ================/

//================= graph.c declaration begin ================/
short constructGraph(char *graphFileName, char **readsFileNames, int readsFileNum);
short constructGraphFromReadset(char *graphFileName, char **readsFileNames, int readsFileNum);
short countKmerOccsFromReadset(graphtype *deBruijnGraph);
short addKmerRidposFromReadset(graphtype *graph);
short addReadPreFromReadset(read_t *pRead, graphtype *graph);
short addReadFromReadset(int64_t rid, read_t *pRead, graphtype *graph);
short generateKmerSeqIntFromReadset(uint64_t *seqInt, uint64_t *readseq, int32_t startReadPos, int32_t entriesNum, int32_t baseNumLastEntry);
short getMinReadLenFromFastaFiles(int *readLenInFile, char **readFilesInput, int readFileNum);
short getMinReadLenFromFastqFiles(int *readLenInFile, char **readFilesInput, int readFileNum);
short getReadLenFromFasta(int *tmpReadLen, char *fastqFile);
short getReadLenFromFastq(int *tmpReadLen, char *fastqFile);
short fillReadsToBufFasta(FILE *fpReads, readBuf_t *pBuf, uint64_t *readsNum);
short fillReadsToBuf(FILE *fpReads, readBuf_t *pBuf, uint64_t *readsNum);
short getSingleReadFasta(FILE *fpPE, readBuf_t *pReadBuf);
short getSingleReadFastq(FILE *fpPE, readBuf_t *pReadBuf);
short containUnknownBase(char *seq);
int32_t getUnknownBaseNum(char *seq);
float calcAverQual3End(char *qual_data, int32_t seqLen);
float calcAverQual5End(char *qual_data, int32_t seqLen);
short singleQualSatisfied(char *qual_data);
float getRatioBase(char *seq, char targetBase);
short generateKmerSeqInt(uint64_t *seqInt, char *seq);
uint64_t kmerhashInt(uint64_t *seqInt);
short initgraph(graphtype **graph);
short initGraphAfterCorrection(graphtype *graph);
short countKmer(uint64_t hashcode, graphtype *graph);
short addKmer(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, uint16_t rpos, graphtype *graph);
short delKmer(char *str, uint64_t rid, unsigned short rpos, graphtype *graph);
short delKmerByHash(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, uint32_t rpos, graphtype *graph);

int findStartIndex(uint64_t rid, ridpostype *rid_pos_table, int posNum);
int getExectIndex(uint64_t rid, uint16_t rpos, int startIndex, ridpostype *ridpostable, int posNum);
int findDirectIndex(uint64_t rid, uint16_t rpos, ridpostype *ridpostable, int posNum);
inline kmertype *getKmerByBase(char *tmp_kmerseq, graphtype *graph);
inline kmertype *getKmer(uint64_t *kmerSeqInt, graphtype *graph);
inline kmertype *getKmerByHash(uint64_t hashvalue, uint64_t *kmerSeqInt, graphtype *graph);
short identicalKmerSeq(uint64_t *kmerSeqInt1, uint64_t *kmerSeqInt2);
char *getKmerBaseByInt(uint64_t *kmerSeqInt);
kmertype *getReverseKmer(uint64_t *kmerSeqIntRev, uint64_t *kmerSeqInt, graphtype *graph);
unsigned int getReverseKmerseqIntByHash(unsigned int hashcode);
void reverse(char *str);
short outputGraphToFile(char *graphFile, graphtype *graph);
short loadGraph(graphtype **graph, char *graphFile);
short loadReadSetFromGraphFile(readSet_t **readSet, char *graphFile);
short GlobalParasFromGraph(int32_t *readLenPara, int32_t *kmerSizePara, uint64_t *hashTableSizePara, uint64_t *hashTableSizeReadseqPara, int32_t *pairedModePara, int32_t *PEGivenTypePara, double *meanSizeInsertPara, double *standardDevPara, char *graphFileName);
short updateInsertSizeAndSdevInGraph(int64_t PEGivenType, double meanSizeInsert, double standardDev, const char *graphFileName);
short releaseGraph(graphtype **graph);
short cleanKmerInfoInGraph(graphtype **graph);
short resetGraph(graphtype *graph);
short recoverKmersInRead(int64_t rid, uint64_t *readseqInt, int32_t seqLen, graphtype *graph);
//================= graph.c declaration end ================/

//================= kmerseqBlock.c declaration begin ================/
short initKmerseqBlockInGraph(graphtype *graph);
short addNewBlockKmerSeq(graphtype *graph);
//================= kmerseqBlock.c declaration end ================/

//================= ridposBlock.c declaration begin ================/
short initRidposBlocksInGraph(graphtype *graph);
short addNewBlockRidpos(graphtype *graph);
//================= ridposBlock.c declaration end ================/

//================= kmerBlock.c declaration begin ================/
short initKmerBlockInGraph(graphtype *graph);
short addNewBlockKmer(graphtype *graph);
//================= kmerBlock.c declaration end ================/

//================= contig.c declaration begin ================/
short buildContigs(char *contigFile, char *graphFileName);
short initMemory();
void freeMemory();
short initFirstKmerBounder(double *lowerBoundFirstKmer, double *upperBoundFirstKmer, double averKmerOccNum);
short initContig();
struct kmertype *getFirstKmer(graphtype *graph, unsigned int *kmerIndex);
short initFirstKmerThreshold();
short getFirstKmers(uint64_t *kmerIndex, kmertype **firstKmer);
short isValidFirstKmer(kmertype *kmer, uint64_t *seqInt);
short containFirstPos(kmertype *kmer);
short containLastPos(kmertype *kmer);
short initAssemblingReadTable(assemblingreadtype *assemblingreads);
short getNextKmerBySE(int contigNodesNum);
short addContigBase(uint32_t baseInt);
short addRidposToContig(successRead_t *successReadArray, int32_t successReadNum, int32_t contigNodesNum);
void cleanContigArray(contigtype *contigArr, int64_t *contigNodesNum);
short outputContigToFile(FILE *fpContig, int outFileType, int contigID, contigtype *contigArr, int64_t contigNodeNum);
short addFirstKmerToDecisionTable(kmertype **kmers);
short addReadToDecisionTable(uint64_t rid, int32_t rpos, int32_t orientation, int32_t matedFlag, int32_t seqLen, uint64_t *readseq, int32_t mismatchNum, int64_t itemNumContigArray);
short replaceReadInDecisionTable(uint64_t rid, int32_t rpos, int32_t orientation, int32_t matedFlag, int32_t mismatchNum, assemblingreadtype *dtReadOld, int64_t itemNumContigArray, contigPath_t *contigPath);
short reallocateDecisionTable();
short computeLongKmerOccNum(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t length_k, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag);
short computeKmerOccNum(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag);
short computeKmerOccNumUnlocked(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag);
short computeKmerOccNumLocked(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag);
short delReadsFromGraph(successRead_t *successReadArray, int successReadNum);
short reverseReadseq(char *seq);
short updateLockedReads();
short updateContigtailnodes(contigtype *contigArr, int64_t successContigIndex, int64_t *contigNodesNum);
short updateContigNodes(contigtype *contigArr, int64_t *validHeadRowContigArray, int64_t *validTailRowContigArray, int64_t *successContigIndex, int64_t *contigNodesNum);
void updateContigheadnodes(contigtype **contighead, int *contigNodeNum);

short getSecondAssemblyFirstKmers(contigtype *contigArr, int64_t contigNodesNum, graphtype *graph);
short getNewHeadContigIndex(int64_t *contigHeadIndexTmp, contigtype *contigArr, int64_t contigNodesNum);
short trimContigBeforeCycle2(contigtype *contigArr, int64_t *successContigIndex, int64_t *contigNodeNum);
short reverseContig(contigtype *contigArr, int64_t contigNodesNum);
short initAssemblingTableSecondAssembly(contigtype *contigArray, int64_t contigNodeNum, graphtype *graph);

short getReversedSeq(char *reversed_seq, char *seq, int seq_len);
short recoverRemainedKmers(char *seq, uint64_t *tmp_kmerseq, uint64_t rid, uint16_t rpos, graphtype *graph);
short recoverKmerByHash(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, uint16_t rpos, graphtype *graph);
short recoverReadFromGraph(char *seq, uint64_t rid, graphtype *graph);
short recoverDeledReads( contigtype *startContig);

short initSecondAssembly();
short updateReadsNumReg(int itemNumSuccessReadsArr, int contigNodesNum, int assemblyRound);
short initReadsNumRegSecondAssembly(int32_t contigNodesNum);
short addOccsToContig(contigtype *contigArray, int32_t contigNodesNum, int32_t naviTandFlag, int32_t newCandBaseNumAfterTandPathPE, contigPath_t *contigPath);
short getValidSuccessReadFlag(int32_t *validSuccessReadFlag, successRead_t *successReadsArray, int32_t itemNumSuccessReadsArray, int32_t minMatchNum);
short confirmSingleEndNavi(int32_t *thisNaviSuccessFlag, contigtype *contigArray, int32_t itemNumContigArray);
//================= contig.c declaration end ================/

//================= dtRowHash.c declaration begin ================/
short addReadToDTRowHashtable(int64_t rid, int32_t dtRow, dtRowIndex_t **pDtRowHashtable);
short delReadFromDTRowHashtable(int64_t rid, int32_t dtRow, dtRowIndex_t **pDtRowHashtable);
short updateReadInDTRowHashtable(int64_t rid, int32_t dtRowOld, int32_t dtRowNew, dtRowIndex_t **pDtRowHashtable);
short cleanDTRowIndexHashtable(dtRowIndex_t **pDtRowHashtable);
int32_t getProperDtRow(assemblingreadtype *assemblingread, assemblingreadtype *assemblingreads, dtRowIndex_t **pDtRowHashtable);
int32_t getProperDtRowLimited(assemblingreadtype *assemblingread, assemblingreadtype *assemblingreads, dtRowIndex_t **pDtRowHashtable, int limitLastpos);
int32_t existReadInDecisionTable(int64_t rid, dtRowIndex_t **pDtRowHashtable);
short getExistReadInDT(assemblingreadtype **pReadDt, int64_t rid, assemblingreadtype *decisionTable, dtRowIndex_t **pDtRowHashtable);
int32_t existReadWithPosInDecisionTable(int64_t rid, int32_t basePos, int32_t orientation, assemblingreadtype *decisionTable, dtRowIndex_t **pDtRowHashtable);
int32_t getCopyNumOfReadInDecisionTable(int32_t *copyNum, int64_t rid, dtRowIndex_t **pDtRowHashtable);
//================= dtRowHash.c declaration end ================/

//================= update.c declaration begin ================/
short updateDecisionTable(kmertype *tmp_kmers[2], int32_t baseInt_kmer);
ridpostype *getRidpos(assemblingreadtype assemblingRead, ridpostype *rid_pos_table, int posNum);
short reallocateSuccessReadsArr();
short removeFinishedReadsFromDecisionTable();
short updateFinishedReadsInDecisionTable();
short updateAssemblingreadsStatus();
short updateSameReadStatusToFailed(int64_t rid, int32_t dtRow, dtRowIndex_t **pDtRowHashtable);
short getCandPathFromSingleCopyReadInDT(candPath_t *candPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, dtRowIndex_t **dtRowHashtable);
short getCandPathItemsFromSingleCopyReadInDT(candPath_t *candPath,  assemblingreadtype *decisionTable, int32_t readsNumDT, dtRowIndex_t **dtRowHashtable);
short rmIncorrectRowInDecisionTable(int32_t *bestRowInDT, int64_t rid, int32_t readseqLen, candPath_t *candPath, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable);
short initMemMultiCopyReadInDecisionTable(multiCopyReadErrNum_t **errBaseNumArray, int32_t *readCopyNum, int64_t rid, int32_t readseqLen, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable);
void releaseMemMultiCopyReadInDecisionTable(multiCopyReadErrNum_t **errBaseNumArray, int32_t *readCopyNum);
short fillReadseqAndErrNumArray(multiCopyReadErrNum_t *errBaseNumArray, int32_t readCopyNum, int64_t rid, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable);
short getBestRowMultiCopyReadInDT(int32_t *bestRow, multiCopyReadErrNum_t *errBaseNumArray, int32_t readCopyNum, candPath_t *candPath, assemblingreadtype *decisionTable);
short remainBestCopyReadInDT(int32_t bestRow, multiCopyReadErrNum_t *errBaseNumArray, int32_t readCopyNum, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable);
short getMismatchNumWithContigPath(int32_t *matchFlag, int32_t *mismatchNum, uint64_t *readseqInt, int32_t readseqLen, int32_t kmerpos, int32_t orientation, int32_t itemNumDT, contigPath_t *contigPath);
short getReadtailMismatchNumContigPath(int32_t *tailMismatchNum, uint64_t *readseqInt, int32_t readseqLen, int32_t kmerpos, int32_t orientation, contigPath_t *contigPath);
//================= update.c declaration end ================/

//================= hashPE.c declaration begin ================/
short estimateInsertSizeAndSdev(char *graphFileName);
short initPEHashParas();
short updatePEHashTable(int contigNodesNum, int assemblyRound);
short getReadFromPEHashtable(PERead_t **pRead, uint64_t readID);
short addReadToPEHashtable(successRead_t *ridposOrient, int contigPos, int assemblyRound);
short delReadfromPEHashtable(uint64_t readID);
short cleanReadsFromPEHashtable();
short initPEHashtableSecondAssembly(contigtype *contigArray, int contigNodesNum);

short meanSizeInsertAndSdevEstimation(const char *fragmentSizeFile, const char *graphFileName);
short getPairedEndsFromSingleContig(FILE *fpFragSize, contigtype *contigArr, int64_t itemNumContigArr);
short initMemGetPESingleContig(contigtype *contigArr, int64_t itemNumContigArr);
void freeMemGetPESingleContig();
short getTotalReadsNumOfSingleContig(int64_t *totalReadsNum, contigtype *contigArr, int64_t contigNodesNum);
short fillDataReadPosTmpArr(readPosTemp_t *readPosTmpArray, contigtype *contigArr, int64_t contigNodesNum);
short radixSortReadPosTmpArr(readPosTemp_t *readPosTmpArray, readPosTemp_t *readPosTmpArrBuf, int64_t contigNodesNum);
short fillDataReadList(readList_t *readListArray, readPos_t *readPosArray, int64_t *itemNumInReadListArray, int64_t *itemNumInReadPosArray, readPosTemp_t *readPosTmpArray, int64_t itemNumInReadPosTmpArray);
short outputInsertSizeToFile(FILE *fpFragSize, readList_t *readListArray, readPos_t *readPosArray, int64_t itemNumInReadListArray);
short computeInsertSizeAndSdev(double *meanSizeInsert, double *standardDev, const char *fragmentSizeFile);
short addSuccessReadsToPEHashtable(successRead_t *successReadsArray, int32_t itemNumSuccessReadsArray, int32_t assemblyRound);
//================= hashPE.c declaration end ================/

//================= PEAssembly.c declaration begin ================/
short buildEstContigs(char *contigFile);
short getNextKmerByMix(int contigNodesNum, int assemblyRound);
short getNextKmerByPE(int contigNodesNum);
short computeKmerOccNumByPE(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag);
short computeKmerOccNumUnlockedByPE(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag);
short computeKmerOccNumLockedByPE(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag);
short computeLongKmerOccNumByPE(int32_t *occNum, kmertype *tmp_kmers[2],  int32_t baseInt_kmer, int32_t length_k, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag);
short validReadPair(assemblingreadtype **dtReadPaired, uint64_t readID, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable);
short trimContigTail(int64_t *successContigIndex, int64_t *contigNodesNum, int32_t trimLen, int assemblyRound);
short computeGapSizeInContig(int *gapSize, contigtype *contigArr, int64_t contigNodesNum, int assemblyRound);
//================= PEAssembly.c declaration end ================/

//================= lenStatistics.c declaration begin ================/
short contigsLenStatistics(contigGraph_t *contigGraph, int32_t minContigLen);
short radixSortContigsLen(int32_t *contigsLenArray, int32_t *contigsLenBufTmp, int32_t itemNumContigsLenArray);
short computeLenStatistics(int32_t *contigsLenArray, int32_t itemNumContigsLenArray, int32_t minContigLen, char *itemName);
//================= lenStatistics.c declaration end ================/

//================= svm.c declaration begin ================/
short extensionDecisionBySvm(int32_t *naviExtensionFlag, svmSampleVector_t *svmSample, svmModel_t *svmModel);
short loadSvmModel(svmModel_t **svmModel, char *svmKernelFuncName, double *supportVector, double *alphaVector, double biasData, double *scaleData, int32_t rowsNumSV, int32_t colsNumSV, int32_t useScaleDataFlag);
short initMemSvmModel(svmModel_t **svmModel, int32_t rowsNumSV, int32_t colsNumSV, int32_t useScaleDataFlag);
short loadSupportVectorData(svmModel_t *svmModel, double *supportVector);
short loadAlphaVectorData(svmModel_t *svmModel, double *alphaVector);
short loadScaleData(svmModel_t *svmModel, double *scaleData);

short loadSvmModelFromFile(svmModel_t **svmModel, char *kernelFuncName, char *supportVectorFile, char *alphaFile, char *biasFile, char *scaleDataFile);
short initMemSvmModelFromFile(svmModel_t **svmModel, char *supportVectorFile, char *scaleDataFile);
short initSampleMemSvm(int32_t colsNumSupportVector);
void freeSampleMemSvm();
short freeSvmModel(svmModel_t **svmModel);
short getRowsColsNumSupportVectorsFromFile(int32_t *rowsNumSupportVector, int32_t *colsNumSupportVector, char *supportVectorFile);
short loadKernelFuncNameFromFile(svmModel_t *svmModel, char *kernelFuncName);
short loadSupportVectorDataFromFile(svmModel_t *svmModel, char *supportVectorFile);
short loadAlphaVectorDataFromFile(svmModel_t *svmModel, char *alphaVectorFile);
short loadBiasDataFromFile(svmModel_t *svmModel, char *biasFile);
short loadScaleDataFromFile(svmModel_t *svmModel, char *scaleDataFile);
short readLine(char *lineStr, int32_t *lineLen, FILE *fp);
short mySvmClassifySE(int32_t *outClass, svmSampleVector_t *svmSampleSE, svmModel_t *svmModelSE);
double myLinearKernelFuntion(double *x, double *y, int32_t dimNum);
double myRBFKernelFuntion(double *x, double *y, int32_t dimNum);
double myPolyKernelFuntion(double *x, double *y, int32_t dimNum);
short fillSampleDataSVM(svmSampleVector_t *svmSample, double *svmFeatureArray);
//================= svm.c declaration end ================/

//================= contiggraph.c declaration begin ================/
short initContigGraph(contigGraph_t **contigGraph);
short releaseContigGraph(contigGraph_t **contigGraph);
short addContigItemToContigGraph(contigGraph_t *contigGraph, int32_t contigID, int32_t localContigID, contigtype *contigArray, int32_t itemNumContigArray);
short saveReadsMatchInfo(FILE *fpReadMatchInfo, int32_t contigID, contigtype *contigArray, int64_t itemNumContigArray);
short outputContigFromContigGraph(char *contigFile, contigGraph_t *contigGraph, int32_t minContigLen);
//================= contiggraph.c declaration end ================/

//================= candPath.c declaration begin ================/
short decideByCandPathPE(int32_t *naviCandPathPE, int32_t *maxBaseIndexAfterCandPathPE, int32_t *incorrectBaseNumCandPathPE, int32_t *occNumPE, int32_t *occsNumIndexPE, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, contigPath_t *contigPath);
short decideByCandPathSE(int32_t *naviCandPathSE, int32_t *maxBaseIndexAfterCandPathSE, int32_t *incorrectBaseNumCandPathSE, int32_t *occNumSE, int32_t *occsNumIndexSE, int32_t *occsNumPE, int32_t *occsNumIndexPE, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, contigPath_t *contigPath);
short initMemCandPath(candPath_t **candPath);
void releaseMemCandPath(candPath_t **candPath);
short getCandPathPE(candPath_t *candPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable);
short getCandPathSE(candPath_t *candPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable);
short getPathItemsFromDTCandPathPE(candPath_t *candPath, assemblingreadtype *decisionTable, int32_t readsNumDT);
short getPathItemsFromDTCandPathSE(candPath_t *candPath, assemblingreadtype *decisionTable, int32_t readsNumDT);
short getMatchRowCandPathItem(int32_t *matchRow, char *readseqTmp, int32_t readseqLenTmp, candPath_t *candPath);
short addReadseqToCandPathItem(char *readseqTmp, int32_t readseqLenTmp, int32_t supportReadsNum, int32_t matchRow, candPath_t *candPath);
short addNewCandPathItem(candPath_t *candPath, char *readseqTmp, int32_t readseqLenTmp, int32_t supportReadsNum);
short adjustCandPath(candPath_t *candPath);
short computeBaseNumArrayCandPath(int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, candPath_t *candPath);
short removeErrorBaseCandPath(int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, candPath_t *candPath);
short mergeCandPathItem(candPath_t *candPath);
void outputCandPath(candPath_t *candPath);
short computeIncorrectBaseNumCandPath(int32_t *incorrectNum, int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum);
short decideNaviCandPathPE(int32_t *naviCandPath, int32_t *maxIndex, int32_t *occNumPE, int32_t *occsNumIndexPE, int32_t incorrectNum, int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, contigPath_t *contigPath);
short decideNaviCandPathSE(int32_t *naviCandPath, int32_t *maxIndex, int32_t *occNumSE, int32_t *occsNumIndexSE, int32_t incorrectNum, int32_t *occNumPE, int32_t *occNumIndexPE, int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, contigPath_t *contigPath);
short confirmNaviByShiftOpCandPath(int32_t *naviCandPath, int32_t *maxIndex, int32_t *incorrectBaseNum, candPath_t *candPath, contigPath_t *contigPath);
short confirmNaviCandPathByContigPath(int32_t *maxIndex, int32_t *occNumArray, int32_t *occsNumIndexArray, contigPath_t *contigPath);
//================= candPath.c declaration end ================/


//================= tandPath.c declaration begin ================/
short decideByCheckingTandemRepeatPE(int32_t *naviCandPathPE, int32_t *maxBaseIndexAfterTandPathPE, int32_t *incorrectBaseNumTandPath, int32_t *newCandBaseNumAfterTandPathPE, int32_t *occNumArray, int32_t *occIndexArray, int32_t contigPos, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, contigPath_t *contigPath, graphtype *graph);
short decideByCheckingTandemRepeatSE(int32_t *naviTandPathSE, int32_t *maxBaseIndexAfterTandPathSE, int32_t *incorrectBaseNumTandPath, int32_t *newCandBaseNumAfterTandPathSE, int32_t *occNumArray, int32_t *occIndexArray, int32_t contigPos, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, contigPath_t *contigPath, graphtype *graph);
void releaseTandemPathList(tandemPathItem_t **tandemPathList);
short getTandemPathPE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, graphtype *graph);
short getTandemPathSE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, graphtype *graph);
short getTandemPathFromNewReadsPE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, uint64_t *kmerseqInt, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, graphtype *graph);
short getTandemPathFromNewReadsSE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, uint64_t *kmerseqInt, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, graphtype *graph);
short getTandemPathFromSingleKmerTandPathPE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, kmertype *tmp_kmers[2], graphtype *graph, assemblingreadtype *decisionTable, int32_t itemNumDT);
short getTandemPathFromSingleKmerTandPathSE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, kmertype *tmp_kmers[2], graphtype *graph, assemblingreadtype *decisionTable, int32_t itemNumDT);
short getTandemPathFromDecisionTablePE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable);
short getTandemPathFromDecisionTableSE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable);
short getMatchedTandPath(tandemPathItem_t **targetTandPath, int32_t *perfectMatchFlag, tandemPathItem_t *tandemPathList, char *readseq, int32_t readseqLen, int32_t contigtailReadPos);
short addTandemPath(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, char *readseqTmp, int32_t readseqLen, int32_t contigtailReadPos, int32_t contigtailPos, int32_t dtRow);
short addReadToTandPath(tandemPathItem_t *targetTandPath, char *readseq, int32_t readseqLen, int32_t contigtailReadPos, int32_t dtRow);
short adjustTandemPathByShiftOverlap(tandemPathItem_t *tandemPathList, int32_t *itemNumTandemPathList, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, graphtype *graph);
short adjustTandemPath(tandemPathItem_t *tandemPathList, int32_t *itemNumTandemPathList);
short getMaxRowNumBaseNumArrayTandPath(int32_t *maxRowNumBaseNumArray, int32_t *contigtailRowBaseNumArray, tandemPathItem_t *tandemPathList);
short computeBaseNumArrayTandPath(int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, int32_t contigtailRowBaseNumArray, tandemPathItem_t *tandemPathList);
short removeErrorBaseTandPath(int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, int32_t contigtailRowBaseNumArray, tandemPathItem_t *tandemPathList);
short mergeAdjustedTandPath(tandemPathItem_t *tandemPathList, int32_t *itemNumTandemPathList);
short getPerfectMatchFlagTandPath(int32_t *perfectMatchFlag, char *pathseq1, int32_t pathlen1, int32_t contigtailPathPos1, char *pathseq2, int32_t pathlen2, int32_t contigtailPathPos2);
short getMatchFlagTandPath(int32_t *matchFlag, tandemPathItem_t *targetTandPath, tandemPathItem_t *tandPath, int32_t mismatchThreshold);
short mergeTandPathToTargetPath(tandemPathItem_t *targetTandPath, tandemPathItem_t *tandPath);
short outputTandPath(tandemPathItem_t *tandemPathList, int32_t itemNumTandemPathList, assemblingreadtype *decisionTable);
short removeShortFragSizeTandPath(int32_t *shortFragSizeRemovedFlag, tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, contigPath_t *contigPath);
short computeAverFragSizeTandPath(double *averFragSize, tandemPathItem_t *tandPath, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable);
short getContigMatchFlagTandPath(int32_t *matchWithContigFlag, int32_t *mismatchNum, tandemPathItem_t *tandPath, contigPath_t *contigPath);
short removeOverlappedTandPath(int32_t *overlapTandPathRemovedFlag, tandemPathItem_t *tandemPathList, int32_t *itemNumTandemPathList, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable);
short computeOverlapSizeTandPath(int32_t *overlapSize, int32_t *overlapType, tandemPathItem_t *tandPath1, tandemPathItem_t *tandPath2);
short computeShiftSizeTandPathByShiftOp(int32_t *mismatchNum, int32_t *shiftType, int32_t *shiftSize, tandemPathItem_t *tandPath1, tandemPathItem_t *tandPath2);
short removeErrOverlappedReadsInTandPath(int32_t *pathNodeDeleteFlag, tandemPathItem_t *targetTandPath, tandemPathItem_t *tandPath, int32_t overlapSize, assemblingreadtype *decisionTable);
short exchangeNodeInfoTandPath(tandemPathItem_t *tandPath1, tandemPathItem_t *tandPath2);
short decideNaviByTandPath(int32_t *naviTandPath, int32_t *maxBaseIndexAfterTandPath, int32_t *incorrectBaseNumTandPath, int32_t *newCandBaseNumAfterTandPath, int32_t *occNumArray, int32_t *occIndexArray, tandemPathItem_t *tandemPathList, int32_t itemNumTandemPathList);
short getCandPathsTandPath(candPath_t *candPath, tandemPathItem_t *tandemPathList, int32_t itemNumTandemPathList);
short getCandPathsFromTandPathList(candPath_t *candPath, tandemPathItem_t *tandemPathList, int32_t itemNumTandemPathList);
short decideNaviTandPath(int32_t *naviTandPath, int32_t *maxIndex, int32_t *newCandBaseNumAfterTandPath, int32_t incorrectNum, int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum);
short checkTandPath(tandemPathItem_t *tandemPathList, int32_t itemNumTandemPathList);
short removeShortOverlappedTandPath(int32_t *shortFragSizeRemovedFlag, tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable);
short isValidNewCandTandPath(int32_t *validFlag, char *newCandTandPathseq, int32_t newCandTandPathLen, contigPath_t *contigPath);
short getMaxSecRowCandTandPath(int32_t *maxRow, int32_t *secRow, candPath_t *candPath);
//================= tandPath.c declaration end ================/

//================= contigPath.c declaration begin ================/
short updateContigPath(contigPath_t *contigPath, int32_t navigationFlag, kmertype *kmers[2], assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, contigtype *contigArray, int32_t itemNumContigArray);
short initMemContigPath(contigPath_t **contigPath);
void releaseMemContigPath(contigPath_t **contigPath);
short cleanContigPath(contigPath_t *contigPath);
short getMultiReadCopyFlag(int32_t *multiCopyFlag, kmertype *tmp_kmers[2]);
short getContigPathSEUsingCandTandPathPE(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, int32_t itemNumContigArray);
short getContigPathPE(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, int32_t itemNumContigArray);
short getContigPathSE(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, int32_t itemNumContigArray);
short getContigPathFromDecisionTablePE(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable);
short getContigPathFromDecisionTableSE(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable);
short resetMaxOverlapSizeWithContigContigPathItem(contigPath_t *contigPath, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable);
short getMatchedContigPathItem(contigPathItem_t **matchedPathItem, char *readseqTmp, int32_t readseqLenTmp, contigPath_t *contigPath);
short addReadseqToContigPathItem(char *readseqTmp, int32_t readseqLenTmp, int32_t overlapSize, contigPath_t *contigPath, contigPathItem_t *targetPathItem, assemblingreadtype *dtRead);
short addNewContigPathItem(contigPath_t *contigPath, char *readseqTmp, int32_t readseqLenTmp, int32_t overlapSize, assemblingreadtype *dtRead);
short delReadFromContigPath(assemblingreadtype *dtRead, contigPath_t *contigPath);
short adjustContigPath(contigPath_t *contigPath, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable);
short removeNoReadsContigPathItem(contigPath_t *contigPath);
short computeBaseNumArrayContigPath(int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, contigPath_t *contigPath);
short removeErrorBaseContigPath(int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, contigPath_t *contigPath);
short mergeContigPathItem(contigPath_t *contigPath, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable);
short getContigPathSEFromCandTandPathPE(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable);
short removeUnmatchedContigPathItem(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable);
short replacePathReadInfoInDT(contigPathItem_t *newPathItem, contigPathItem_t *oldPathItem, contigPathItemRead_t *pathItemRead, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable);
short removePathReadInfoInDT(contigPathItem_t *pathItem, contigPathItemRead_t *pathItemRead, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable, int32_t delReadFlag);
short sortContigPathItems(contigPath_t *contigPath);
short comparisonSortContigPathItem(contigPathItemSort_t *contigPathItemSortArray, contigPathItemSort_t *contigPathItemSortBufArray, int32_t itemNum);
short radixSortContigPathItem(contigPathItemSort_t *contigPathItemSortArray, contigPathItemSort_t *contigPathItemSortBufArray, int32_t itemNum);
short getMaxesContigPathItems(contigPath_t *contigPath);
short removeLessSupportedContigPathItems(contigPath_t *contigPath, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable);
short removeShortFragSizeContigPath(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable);
short computeAverFragSizeContigPathItem(double *averFragSize, contigPathItem_t *pathItem, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable);
short removeShortOverlappedContigPathItems(contigPath_t *contigPath, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable);
short remainBestContigPathItems(contigPath_t *contigPath, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable);
short shiftStartRowNewBaseContigPath(contigPath_t *contigPath);
void outputContigPath(contigPath_t *contigPath, int32_t allInfoFlag);
short initTailseqContigPath(contigPath_t *contigPath, contigtype *contigArray, int32_t itemNumContigArray);
short appendTailBaseToContigPath(contigPath_t *contigPath, contigtype *contigArray, int32_t itemNumContigArray);
short removeReadsInIdenticalPathItemPart(contigPath_t *contigPath, int32_t useOldNaviPathseqFlag, char *oldNaviPathseq, int32_t oldNaviPathseqLen, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, int32_t itemNumContigArray);
short setNaviContigPathItem(contigPath_t *contigPath, int32_t useOldNaviPathseqFlag, char *oldNaviPathseq, int32_t oldNaviPathseqLen);
short updateNaviContigPathItem(contigPath_t *contigPath, int32_t contigBaseInt);
short decideByContigPath(int32_t *naviContigPath, contigPath_t *contigPath, int32_t *occsNumArray, int32_t *occsNumIndexArray, graphtype *graph,  double occRatioThreshold);
short getMismatchNumByShiftOp(int32_t *mismatchNum, int32_t *shiftSize, int32_t *shiftType, char *seq1, int32_t seqLen1, char *seq2, int32_t seqLen2);
//================= contigPath.c declaration end ================/

//================= scafPerga.c declaration begin ================/
short startScaffolding(char *scafSeqFile, char *readMatchInfoFile, contigGraph_t *contigGraph, readSet_t *readSet);
short initMemScaffolding(scaffoldSet_t **scaffoldSet);
void freeMemScaffolding(scaffoldSet_t **scaffoldSet, scafContigIndex_t **scafContigIndex);
void releaseScaffoldSet(scaffoldSet_t **scaffoldSet);
short loadContigGraph(contigGraph_t **contigGraph, char *contigsFile);
short fillContigItems(contigGraph_t *contigGraph, char *contigsFile);
short getMaxContigLenFromFile(int64_t *maxContigLen, const char *contigFileName);
short getSingleContigFastaFromFile(FILE *fp_contig, char *contigHeadTitle, char *contigSeq, int64_t *contigLen);
short outputScaffoldSetToFile(char *tmpScafFile, scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph);
//================= scafPerga.c declaration begin ================/

//================= scafContigIndex.c declaration begin ================/
short buildScafContigIndex(scafContigIndex_t **scafContigIndex, contigGraph_t *contigGraph);
short initScafContigIndex(scafContigIndex_t **scafContigIndex);
void releaseScafContigIndex(scafContigIndex_t **scafContigIndex);
short initScafKmerBlockInGraph(scafContigIndex_t *scafContigIndex);
short addNewBlockScafKmer(scafContigIndex_t *scafContigIndex);
short initScafKmerseqBlockInGraph(scafContigIndex_t *scafContigIndex);
short addNewBlockScafKmerSeq(scafContigIndex_t *scafContigIndex);
short countScafKmerOccs(scafContigIndex_t *scafContigIndex, contigGraph_t *contigGraph);
short countScafKmer(uint64_t hashcode, uint64_t *kmerSeqInt, scafContigIndex_t *scafContigIndex);
scafKmer_t *getScafKmerByHash(uint64_t hashvalue, uint64_t *kmerSeqInt, scafContigIndex_t *scafContigIndex);
short addScafKmerContigpos(scafContigIndex_t *scafContigIndex, contigGraph_t *contigGraph);
short initContigposBlocksInContigIndex(scafContigIndex_t *scafContigIndex);
short addNewBlockContigpos(scafContigIndex_t *scafContigIndex);
short addScafKmer(uint64_t hashcode, uint64_t *kmerSeqInt, int32_t contigID, int32_t contigPos, int32_t contigEndFlag, scafContigIndex_t *scafContigIndex);
//================= scafContigIndex.c declaration end ================/

//================= scafMap.c declaration begin ================/
short mapReads(contigGraph_t *contigGraph, readSet_t *readSet, scafContigIndex_t *scafContigIndex);
short getMaxArraySizeFromContigIndex(int32_t *maxArraySize, scafContigIndex_t *scafContigIndex);
short mapSingleReadInScaf(int64_t rid, read_t *pRead, readSet_t *readSet, scafContigpos_t *matchResultArray, scafContigpos_t *matchResultArrayRev, scafContigpos_t *matchResultArrayBuf, int32_t *matchItemNum, int32_t *matchItemNumRev, scafContigIndex_t *scafContigIndex);
short getMatchedContigPos(scafContigpos_t *matchResultArray, scafContigpos_t *matchResultArrayBuf, int32_t *matchItemNum, uint64_t *readSeqInt, int32_t seqLen, scafContigIndex_t *scafContigIndex);
short fillReadMatchInfoContigEnds(contigGraph_t *contigGraph, readSet_t *readSet);
short radixSortContigReadArrayInScaf(contigRead_t *contigReadArray, contigRead_t *contigReadArrayBuf, int32_t itemNum);
short outputContigReadArrayInScaf(contigGraph_t *contigGraph);
//================= scafMap.c declaration end ================/

//================= scafLink.c declaration begin ================/
short contigsLinking(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph, readSet_t *readSet);
short initContigItemInfo(contigGraph_t *contigGraph);
short generateContigEdges(contigGraph_t *contigGraph, readSet_t *readSet);
short initContigEdges(contigGraph_t *contigGraph, readSet_t *readSet);
short calcEdgeNumSingleContigEnd(int32_t *contigEdgeNum, int32_t contigID, contigRead_t *contigReadArray, int32_t contigReadsNum, int32_t *contigEndBucketArray, int32_t bucketArraySize, readSet_t *readSet);
short fillContigEdges(contigGraph_t *contigGraph, readSet_t *readSet);
short fillSituationArray(contigEdge_t *pContigEdge, int32_t contigID1, int32_t contigID2, int32_t endFlag1, int32_t endFlag2, contigGraph_t *contigGraph, readSet_t *readSet);
short removeInvalidGraphEdge(contigGraph_t *contigGraph);
short initParaLinking(contigGraph_t *contigGraph);
short computeAverPairsEachContigEdge(double *averageLinkNum, contigGraph_t *contigGraph);
short linkContigs(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph);
short initContigLinkSet(contigLink_t **contigLinkSet, int32_t contigsNum);
void freeMemContigLinkArr(contigLink_t **contigLinkSet);
short getFirstLinkedContigs(int32_t *firstContigID, maxRowCol_t *pMaxRowColNode, contigGraph_t *contigGraph);
short getMaxColsOfSingleRow(maxRowCol_t *pMaxRowColNode, contigGraph_t *contigGraph);
short getMaxRowsOfSingleCol(maxRowCol_t *pMaxRowColNode, contigGraph_t *contigGraph);
short changeMaxRowCol(maxRowCol_t *pMaxRowColNode, contigLink_t *contigLinkSet, contigGraph_t *contigGraph, int32_t linkRound, int32_t turnRoundFlag);
short isLinkSingleton(maxRowCol_t *pMaxRowColNode, contigGraph_t *contigGraph, int32_t linkRound, int32_t strictFlag);
short addContigToContigLinkSet(int32_t *linkStatus, contigLink_t *contigLinkSet, contigGraph_t *contigGraph, maxRowCol_t *pMaxRowColNode, int32_t newContigNum, int32_t linkRound);
short markContigGraphEdge(contigGraph_t *contigGraph, maxRowCol_t *pMaxRowColNode);
short saveLinkResultToScaffoldSet(scaffoldSet_t *scaffoldSet, contigLink_t *contigLinkSet, int32_t linkID);
short saveUnlinkedContigsToScaffoldSet(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph, int32_t *linkID);
//================= scafLink.c declaration end ================/

//================= scafOverlap.c declaration begin ================/
short overlapContigsInScaf(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph, readSet_t *readSet);
short initMemContigOverlap(contigGraph_t *contigGraph);
void freeMemContigOverlap();
short computeContigOverlapInfo(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph, readSet_t *readSet);
short updateContigOverlapLen(contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph, readSet_t *readSet);
short reverseSeq(char *seq, int seq_len);
short gapSizeEstimateBetweenContigs(double *gapSize, int32_t *validPairedNum, contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph, readSet_t *readSet, int32_t endCutRound, int32_t *cutOrderArray, int32_t *uncoveredEndLenArray);
short updateContigEndsInfo(contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph);
short computeSeqOverlapLenExact(int32_t *overlapLen, const char *seq1, const int32_t seqLen1, const char *seq2, const int32_t seqLen2, int32_t *scoreArray, double gapSize);
short computeSeqOverlapLenByAlignment(const char *seq1, const int32_t seqLen1, const char *seq2, const int32_t seqLen2, int32_t *scoreArray, char **pAlignResultArray, int32_t *overlapLen, int32_t *mismatchNum);
short adjustOverlapSeq(char *seq1, char *seq2, char **pAlignResultArray, int *overlapLen);
short updateContigs(contigGraphItem_t *pContigItem1, contigGraphItem_t *pContigItem2, const int contigOrient1, const int contigOrient2, char *seq1, const int originalSeqLen1, char *seq2, const int originalSeqLen2);
short updateOverlapLenByCutUncoveredContigEnds(contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph);
short getUncoveredLenAtContigEnds(int32_t *pUncoveredEndLenArray, contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph);
short cutContigEnds(contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph, int32_t endCutRound, int32_t *cutOrderArray, int32_t *uncoveredEndLenArray);
short splitScaffolds(scaffoldSet_t *scaffoldSet);
//================= scafOverlap.c declaration end ================/

//================= scafGap.c declaration begin ================/
short gapFilling(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph, readSet_t *readSet);
short constructScafGraphFromReadset(graphtype **deBruijnGraph, contigGraph_t *contigGraph, readSet_t *readSet);
short countKmerOccsFromReadsetInScaf(graphtype *deBruijnGraph);
short addKmerRidposFromReadsetInScaf(graphtype *graph);
short fillGaps(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph, readSet_t *readSet);
short initGapFillingParas();
void freeGapFillingParas();
short localAssemblyInScaf(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph, readSet_t *readSet);
short getScafContigEndSeqs(char *endSeqArray[2], int32_t *endSeqLenArray, contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph);
short getComparisonSeqInScaf(char *comparisonSeqInScaf, int32_t *comparisonSeqLenInScaf, char *endSeqArray[2], int32_t *endSeqLenArray, int32_t assemblyCycle);
short prepareAssemblyInScaf(contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph, char *endSeqArray[2], int32_t *endSeqLenArray, int32_t assemblyCycle);
short initScafContig(contigtype *contigArray, int64_t *contigNodesNum, char *seq, int32_t seq_len);
short initPEHashtableInScaf(contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph, int32_t assemblyCycle);
short getMaxMinAssemblyLen(int32_t *maxAssemblyLen, int32_t *minAssemblyLen, int32_t assemblyCycle, int32_t *localScafContigNodesNumArr, int32_t gapSize);
short getNewContigSeq(char *localScafContigSeq, int32_t *localScafContigNodesNum, contigtype *contigArray, int32_t itemNumContigArray);
short updateScafContigEndSeqs(char *localScafContigSeq, int32_t localScafContigNodesNum, int32_t assemblyCycle);
short detectOverlapsInScaf(int32_t *successFilledFlag, int32_t *overlapLen, int32_t *newGapSize, int32_t *breakFlag, int32_t gapSize, int32_t *localScafContigNodesNumArray, int32_t assemblyCycle);
short computeNewGapSizeInScaf(int32_t *tmpGapSize, int32_t gapSize, int32_t *localScafContigNodesNumArray, int32_t assemblyCycle);
short updateContigOverlapInfoInScaf(contigOverlap_t *pContigOverlapInfo, int32_t successFilledFlag, int32_t overlapLen, int32_t newGapSize, int32_t breakFlag);
short updateContigInfoInScaf(contigOverlap_t *pContigOverlapInfo, char **localScafContigSeqArr, int32_t *localScafContigNodesNumArr, int32_t *oldEndSeqLenArr, contigGraph_t *contigGraph, int32_t assemblyCycle);
short updateSingleContigInfoInScaf(contigGraphItem_t *contigItem, int32_t contigOrient, int32_t contigIndex, char *localScafContigSeq, int32_t scafContigNodesNum, char *endSeq, int32_t endSeqLen, int32_t oldEndSeqLen, int32_t prepareAssemblyLen);
short updateOtherContigInfoInScaf(contigGraphItem_t *contigItem, int32_t contigOrient, char *contigEndSeq, int32_t endSeqLen, int32_t oldEndSeqLen);
//================= scafGap.c declaration end ================/

//================= scafSeq.c declaration begin ================/
short generateScafSeq(char *scafSeqFile, scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph);
short getMaxContigLenInSingleScafflod(scaffoldItem_t *scaffoldItem, int32_t *maxContigLen, contigGraph_t *contigGraph);
short getApproximateSingleScafflodLen(scaffoldItem_t *scaffoldItem, int32_t *approximateScaffoldLen, contigGraph_t *contigGraph);
short contigsLenStatisticsScaf(scaffoldSet_t *scaffoldSet, int32_t minContigLen);
//================= scafSeq.c declaration end ================/

#endif // METHODS_H_INCLUDED
