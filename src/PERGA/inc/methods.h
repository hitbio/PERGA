#ifndef METHODS_H_INCLUDED
#define METHODS_H_INCLUDED

#include "structure.h"

//================= srgaMain.c declaration begin ================/
short parseCommandParasAndExe(int argc, char **argv);
short showUsageInfo();
//================= srgaMain.c declaration end ================/

//================= srga.c declaration begin ================/
short startSRGA(int operationModePara, int kmerSizePara, int readLenCutOffPara, int pairedModePara, char **readFilesPara, int readFileNumPara, int singleBaseQualThresPara, int errorCorrectionPara, char *graphFilePara, double meanSizeInsertPara, double standardDevPara, char *outputPathPara, char *outputPrefixPara, int minContigLenPara);
short initGlobalParas(int operationModepara, char *outputPathName, char *prefix, char **readFilesPara, int readFileNumPara, int singleBaseQualThresPara, int errorCorrectionPara, int pairedModePara, int kmerLen, int readLenCut, char *graphFilePara, double meanSizeInsertPara, double standardDevPara, int minContigLenPara);
short setGlobalPath(const char *outPathStr);
void freeGlobalParas();
short getReadsFileFormat(int *readsFileFormatType, char **readFilesInput, int readFileNum);
short getReadLenAndCountFromFilesFasta(int64_t *tmpTotalReadNum, int *tmpAverReadLenInFile, int *tmpMaxReadLenInFile, int *tmpMinReadLenInFile, char **readsFileNames, int readsFileNum);
short getReadLenAndCountFromFilesFastq(int64_t *tmpTotalReadNum, int *tmpAverReadLenInFile, int *tmpMaxReadLenInFile, int *tmpMinReadLenInFile, char **readsFileNames, int readsFileNum);
short getReadLenAndCountFromSingleFasta(int64_t *tmpReadCount, int64_t *tmpSumReadLen, int64_t *tmpMaxReadLen, int64_t *tmpMinReadLen, char *fastaFile);
short getReadLenAndCountFromSingleFastq(int64_t *tmpReadCount, int64_t *tmpSumReadLen, int64_t *tmpMaxReadLen, int64_t *tmpMinReadLen, char *fastqFile);
int32_t getSysMemorySize();
//================= srga.c declaration end ================/


//================= util.c declaration begin ================/
short outputContig(contigtype *contigArray, int64_t contigNodesNum);
short outputContigEnd3(contigtype *contigArray, int64_t contigNodesNum, int endNum);
short outputContigEnd5(contigtype *contigArray, int64_t contigNodesNum, int endNum);
short outputUndelKmerpos(graphtype *graph);
short outputRemainedKmers(graphtype *graph);
short outputReadsInDecisionTable(assemblingreadtype *decisionTable, int readsNum);
short outputReadsInDecisionTableToFile(char *outfile, assemblingreadtype *decisionTable, int readsNum);
short outputFailedReadsInDecisionTable(assemblingreadtype *decisionTable, int itemNumDecisionTable, int contigID, int contigNodesNum);
short outputLockedReadsInDecisionTable(assemblingreadtype *decisionTable, int readsNum);
short outputMatedReadsInDecisionTable(assemblingreadtype *decisionTable, int readsNum);
short outputKmer(graphtype *graph, int hashcode, uint64_t *kmerSeqInt); //输出kmer中的内容
short outputRidpos(ridpostype *ridpos, int posNum);
void outputSuccessReads(successRead_t *successReadArray, int successReadNum);
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
short testSVM();
short outputClassificationResult();
short resultStatistics();
short calcErrNumSingleFile(int32_t *errClassNumArray, char *svmClassFile);
short outputErrPoints();
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
short addNewBlockReadseq(readSet_t *pReadSet);
short initReadseqHashtableInReadset(readSet_t *pReadSet);
short initReadseqHashItemBlockInGraph(readSet_t *pReadSet);
short addNewBlockReadseqHashItem(readSet_t *pReadSet);
short releaseReadset(readSet_t *pReadSet);
void releaseHashItemReadset(readSet_t *readSet);
uint64_t readseqHashInt(uint64_t *seqInt, int32_t baseNum, int32_t entriesNum);
short generateReadseqInt(uint64_t *seqInt, char *seq, int32_t seqLen, int32_t entriesNum);
inline readseqHashItem_t *getReadseqHashItemByHash(uint64_t hashvalue, uint64_t *readseqInt, int32_t seqLen, int32_t entriesNum, readSet_t *readSet);
short identicalReadseq(uint64_t *kmerSeqInt1, uint64_t *kmerSeqInt2, int32_t entriesNum);
char *getReadBaseByInt(uint64_t *readseqInt, int32_t seqLen);
char *getReverseReadBaseByInt(uint64_t *readseqInt, int32_t seqLen);
short getReverseReadseqInt(uint64_t *readseqIntRev, uint64_t *readseqInt, int32_t seqLen);
short replaceUnknownBasesInReads(char *seq, int32_t nBaseNum);
short computeMaxReadLenInReadset(readSet_t *readSet);
short extractRefPosFromHeadName(int32_t *strandPosArray, char *headname, int32_t headlen);
//================= reads.c declaration end ================/

//================= graph.c declaration begin ================/
short constructGraph(char *graphFileName, char **readsFileNames, int readsFileNum);
short constructGraphAfterCorrection(graphtype *graph, char *graphFileCorrected, char *correctedReadsFile);
short constructGraphFromReadset(char *graphFileName, char **readsFileNames, int readsFileNum);
short countKmerOccsFromReadset(graphtype *deBruijnGraph);
short countKmerOccsFromReadsetAfterCorrection(graphtype *deBruijnGraph);
short addKmerRidposFromReadset(graphtype *graph);
short addKmerRidposFromReadsetAfterCorrection(graphtype *graph);
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
float getRatioBaseA(char *seq);
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
short GlobalParasFromGraph(int *readLenPara, int *averReadLenInFilePara, int *kmerSizePara, uint64_t *hashTableSizePara, uint64_t *hashTableSizeReadseqPara, int *pairedModePara, int32_t *PEGivenTypePara, double *meanSizeInsertPara, double *standardDevPara, char *graphFileName);
short updateInsertSizeAndSdevInGraph(int32_t PEGivenType, double meanSizeInsert, double standardDev, const char *graphFileName);
short releaseGraph(graphtype *graph);
short cleanKmerInfoInGraph(graphtype *graph);
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
short buildContigs(char *contigFile, char *graphFileName, char *(*occPointFileArray)[6]);
short initMemory(char *(*occPointFileArray)[6]);
void freeMemory();
short initFirstKmerBounder(double *lowerBoundFirstKmer, double *upperBoundFirstKmer, short assemblyCycle, double averKmerOccNum);
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
short addRidposToContig(successRead_t *successReadArray, int successReadNum, int contigNodesNum);
void cleanContigArray(contigtype *contigArr, int64_t *contigNodesNum);
short outputContigToFile(FILE *fpContig, int outFileType, int contigID, contigtype *contigArr, int64_t contigNodeNum);
short addFirstKmerToDecisionTable(kmertype **kmers);
short addReadToDecisionTable(uint64_t rid, int32_t rpos, int32_t orientation, int32_t matedFlag, int32_t seqLen, uint64_t *readseq);
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
short initAssemblingTableSecondAssembly(contigtype *contigArr, int64_t contigNodeNum, graphtype *graph);

short getReversedSeq(char *reversed_seq, char *seq, int seq_len);
short recoverRemainedKmers(char *seq, uint64_t *tmp_kmerseq, uint64_t rid, uint16_t rpos, graphtype *graph);
short recoverKmerByHash(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, uint16_t rpos, graphtype *graph);
short recoverReadFromGraph(char *seq, uint64_t rid, graphtype *graph);
short recoverDeledReads( contigtype *startContig);

short initSecondAssembly();
short updateReadsNumReg(int itemNumSuccessReadsArr, int contigNodesNum, int assemblyRound);
short initReadsNumRegSecondAssembly(int contigNodesNum);
//================= contig.c declaration end ================/

//================= dtRowHash.c declaration begin ================/
short addReadToDTRowHashtable(int64_t rid, int32_t dtRow, dtRowIndex_t **pDtRowHashtable);
short delReadFromDTRowHashtable(int64_t rid, int32_t dtRow, dtRowIndex_t **pDtRowHashtable);
short updateReadInDTRowHashtable(int64_t rid, int32_t dtRowOld, int32_t dtRowNew, dtRowIndex_t **pDtRowHashtable);
short cleanDTRowIndexHashtable(dtRowIndex_t **pDtRowHashtable);
int32_t getProperDtRow(assemblingreadtype *assemblingread, assemblingreadtype *assemblingreads, dtRowIndex_t **pDtRowHashtable);
int32_t getProperDtRowLimited(assemblingreadtype *assemblingread, assemblingreadtype *assemblingreads, dtRowIndex_t **pDtRowHashtable, int limitLastpos);
int32_t existReadInDecisionTable(int64_t rid, int32_t basePos, int32_t orientation, assemblingreadtype *desisionTable, dtRowIndex_t **pDtRowHashtable);
//================= dtRowHash.c declaration end ================/

//================= update.c declaration begin ================/
short updateDecisionTable(kmertype *tmp_kmers[2], int32_t baseInt_kmer);
ridpostype *getRidpos(assemblingreadtype assemblingRead, ridpostype *rid_pos_table, int posNum);
short reallocateSuccessReadsArr();
short removeFinishedReadsFromDecisionTable();
short updateFinishedReadsInDecisionTable();
short updateAssemblingreadsStatus();
short updateSameReadStatusToFailed(int64_t rid, int32_t dtRow, dtRowIndex_t **pDtRowHashtable);
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
//================= hashPE.c declaration end ================/

//================= PEAssembly.c declaration begin ================/
short buildEstContigs(char *contigFile);
short getNextKmerByMix(int contigNodesNum, int assemblyRound);
short getNextKmerByPE(int contigNodesNum);
short computeKmerOccNumByPE(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag);
short computeKmerOccNumUnlockedByPE(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag);
short computeKmerOccNumLockedByPE(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag);
short computeLongKmerOccNumByPE(int32_t *occNum, kmertype *tmp_kmers[2],  int32_t baseInt_kmer, int32_t length_k, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag);
short validReadPair(uint64_t readID);
short trimContigTail(int64_t *successContigIndex, int64_t *contigNodesNum, int32_t trimLen, int assemblyRound);
short setEmptyNaviOccQueue(double *naviOccQueuePara, int *itemNumNaviOccQueuePara, int *frontRowNaviOccQueuePara, int *rearRowNaviOccQueuePara);
short updateNaviOccQueue(double *naviOccQueuePara, int *itemNumNaviOccQueuePara, int *frontRowNaviOccQueuePara, int *rearRowNaviOccQueuePara, double maxOccNum);
short calcAverOccNaviOccQueue(double *averOccNum, double *naviOccQueuePara, int itemNumNaviOccQueuePara);
short getLowOccLenNaviOccQueue(int *lowLen, double *naviOccQueuePara, int itemNumNaviOccQueuePara, int frontRowNaviOccQueuePara);
short computeGapSizeInContig(int *gapSize, contigtype *contigArr, int64_t contigNodesNum, int assemblyRound);
//================= PEAssembly.c declaration end ================/

//================= align.c declaration begin ================/
short adjustMatchInfoRead(assemblingreadtype *this_assemblingRead, contigtype *contigArray, int32_t contigNodesNum);
short getReadBasesAlign(char *readSeqAlign, int32_t *readSeqLenAlign, assemblingreadtype *this_assemblingRead);
short getContigBasesAlign(char *contigSeqAlign, int32_t *contigSeqLenAlign, assemblingreadtype *this_assemblingRead, contigtype *contigArray, int32_t contigNodesNum);
short adjustMatchInfo(assemblingreadtype *this_assemblingRead, char *readSeqAlign, int32_t readSeqLenAlign, char *contigSeqAlign, int32_t contigSeqLenAlign);
short generateAlignment(char *alignResultSeqArr[3], char *matchFlagArr, int32_t *alignResultSeqLen, char *seq1, int32_t seqLen1, char *seq2, int32_t seqLen2);
short adjustMatchInfoByAlign(assemblingreadtype *this_assemblingRead, char *alignResultSeqArr[3], char *matchFlagArr, int32_t *alignResultSeqLen, char *readSeqAlign, int32_t *readSeqLenAlign, char *contigSeqAlign, int32_t contigSeqLenAlign);
short addNewShiftedReadBases(char *readSeqAlign, int32_t *readSeqLenAlign, int32_t shiftedBaseNum, assemblingreadtype *this_assemblingRead, char *matchFlagArr, int32_t alignResultSeqLen);
short adjustErrBaseInfoByAlign(assemblingreadtype *this_assemblingRead, char *alignResultSeqArr[3], char *matchFlagArr, int32_t alignResultSeqLen, char *readSeqAlign, int32_t readSeqLenAlign, char *contigSeqAlign, int32_t contigSeqLenAlign);
//================= align.c declaration end ================/

//================= correction.c declaration begin ================/
short correctRead(assemblingreadtype *assemblingRead, contigtype *contigArray, graphtype *graph, FILE *fpReadCorrected);
short getOmittedReadBases(char *readSeqAlign, int32_t *readSeqLenAlign, assemblingreadtype *assemblingRead);
short getOmittedContigBases(char *contigSeqAlign, int32_t *contigSeqLenAlign, assemblingreadtype *assemblingRead, contigtype *contigArray);
short generateNewReadseq(uint64_t *readseqCorrected, int32_t *newSeqLen, assemblingreadtype *assemblingRead);
short updateReadInReadset(int32_t newSeqLen, assemblingreadtype *assemblingRead, graphtype *graph);
short updateAndSaveReadToFile(uint64_t *readseqCorrected, int32_t newSeqLen, assemblingreadtype *assemblingRead, graphtype *graph, FILE *fpReadCorrected);
short fillErrBasesSingleRead(assemblingreadtype *assemblingRead, char *matchFlagArr, int32_t alignResultSeqLen, char *contigSeqAlign, int32_t contigSeqLenAlign);
short AdjustReadsetByCorrectedReads(graphtype *graph, char *correctedReadsFile);
//================= correction.c declaration end ================/

//================= lenStatistics.c declaration begin ================/
short contigsLenStatistics(int32_t *contigsLenArray, int32_t itemNumContigsLenArray);
short radixSortContigsLen(int32_t *contigsLenArray, int32_t *contigsLenBufTmp, int32_t itemNumContigsLenArray);
short computeLenStatistics(int32_t *contigsLenArray, int32_t itemNumContigsLenArray);
//================= lenStatistics.c declaration end ================/

//================= ref.c declaration begin ================/
short loadReferenceFromFile(char *refFileName);
short initMemRef(char *refFileName);
short freeReference();
short updateRefPosContig(successRead_t *successReadArray, int32_t successReadNum, ref_t *refArray);
short computeRefPosContig(successRead_t *successReadArray, int32_t successReadNum, ref_t *refArray);
short drawOccPoints(int32_t navigationFlag, int32_t naviSuccessFlag, ref_t *refArr);
//================= ref.c declaration end ================/

//================= svm.c declaration begin ================/
short extensionDecisionBySvm(int32_t *naviExtensionFlag, svmSampleVector_t *svmSample, svmModel_t *svmModel);
short loadSvmModel(svmModel_t **svmModel, char *kernelFuncName, char *supportVectorFile, char *alphaFile, char *biasFile, char *scaleDataFile);
short initMemSvmModel(svmModel_t **svmModel, char *supportVectorFile, char *scaleDataFile);
short initSampleMemSvm(svmModel_t *svmModel);
void freeSampleMemSvm();
short freeSvmModel(svmModel_t **svmModel);
short getRowsColsNumSupportVectors(int32_t *rowsNumSupportVector, int32_t *colsNumSupportVector, char *supportVectorFile);
short loadKernelFuncName(svmModel_t *svmModel, char *kernelFuncName);
short loadSupportVectorData(svmModel_t *svmModel, char *supportVectorFile);
short loadAlphaVectorData(svmModel_t *svmModel, char *alphaVectorFile);
short loadBiasData(svmModel_t *svmModel, char *biasFile);
short loadScaleData(svmModel_t *svmModel, char *scaleDataFile);
short readLine(char *lineStr, int32_t *lineLen, FILE *fp);
short mySvmClassifySE(int32_t *outClass, svmSampleVector_t *svmSampleSE, svmModel_t *svmModelSE);
double myLinearKernelFuntion(double *x, double *y, int32_t dimNum);
double myRBFKernelFuntion(double *x, double *y, int32_t dimNum);
double myPolyKernelFuntion(double *x, double *y, int32_t dimNum);
short fillSampleDataSVM(svmSampleVector_t *svmSample, double *svmFeatureArray);
//================= svm.c declaration end ================/

//================= malign.c declaration begin ================/
short decideByMalignPE(int32_t *naviMalignPE, int32_t *incorrectNumPE, int32_t maxPEIndex, assemblingreadtype *decisionTable, int32_t readsNum, graphtype *graph);
short decideByMalignSE(int32_t *naviMalignSE, int32_t incorrectNumPE, int32_t maxPEIndex, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, graphtype *graph);
short getTotalReadsNumMalignPE(int32_t *totalReadsNum, int32_t readsNumDecisionTable, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable);
short getTotalReadsNumMalignSE(int32_t *totalReadsNum, int32_t readsNumDecisionTable, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable);
short getNewReadsNumMalignPE(int32_t *newReadsNum, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable);
short getNewReadsNumMalignSE(int32_t *newReadsNum, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable);
short getNewReadsNumSingleKmerMalignPE(int32_t *newReadsNum, kmertype *tmp_kmers[2], graphtype *graph, assemblingreadtype *decisionTable);
short getNewReadsNumSingleKmerMalignSE(int32_t *newReadsNum, kmertype *tmp_kmers[2], graphtype *graph, assemblingreadtype *decisionTable);
short addReadseqMalignPE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *itemNumMalignArray, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable, int32_t readsNumDT);
short addReadseqMalignSE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *itemNumMalignArray, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable, int32_t readsNumDT);
short addNewReadseqMalignPE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *newReadsNumTotal, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable);
short addNewReadseqMalignSE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *newReadsNumTotal, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable);
short addNewReadseqSingleKmerMalignPE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *newReadsNumTotal, kmertype *tmp_kmers[2], graphtype *graph, assemblingreadtype *decisionTable);
short addNewReadseqSingleKmerMalignSE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *newReadsNumTotal, kmertype *tmp_kmers[2], graphtype *graph, assemblingreadtype *decisionTable);
short addNewReadseqInDecisionTableMalignPE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *itemNumMalignArray, assemblingreadtype *decisionTable, int32_t readsNumDT);
short addNewReadseqInDecisionTableMalignSE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *itemNumMalignArray, assemblingreadtype *decisionTable, int32_t readsNumDT);
short sortReadseqMalign(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t itemNumMalignArray);
short adjustReadseqMalign(char **readseqMalignArray, int32_t *shiftPosArray, int32_t itemNumMalignArray);
short computeConsensusBaseNumMalign(int32_t **consensusBaseNumArray, int32_t *maxAlignSeqLen, int32_t rowsNum, int32_t colsNum, char **readseqMalignArray, int32_t *shiftPosArray, int32_t itemNumMalignArray);
short computeIncorrectBaseNumMalign(int32_t *incorrectNum, int32_t **consensusBaseNumArray, int32_t colsNum, int32_t maxAlignSeqLen, char *consensusSeq, int32_t itemNumMalignArray);
short decideNaviMalignPE(int32_t *naviMalign, int32_t *maxIndex, int32_t maxPEIndex, int32_t incorrectNum, int32_t **consensusBaseNumArray, int32_t colsNum, int32_t maxAlignSeqLen, char *consensusSeq);
short decideNaviMalignSE(int32_t *naviMalignSE, int32_t incorrectNumSE, int32_t *maxIndex, int32_t maxPEIndex, int32_t **consensusBaseNumArray, int32_t colsNum, int32_t maxAlignSeqLen, char *consensusSeq);
//================= malign.c declaration end ================/


//================= msa.c declaration begin ================/
short generateMSAMalign(char **alignResultsMalign, char **readseqMalignArray, int32_t itemNumMalignArray);
short getCoreIndexMSAMalign(int32_t *coreIndex, char **readseqMalignArray, int32_t itemNumMalignArray);
short getGolbalScoreMalign(int32_t *score, char *seq1, char *seq2);
short computeScorePAMalign(int32_t **scoreArray, char *seq1, char *seq2, int32_t seqlen1, int32_t seqlen2);
short subScoreMalign(int32_t *score, char ch1, char ch2);
short generateMSAResultMalign(char **alignResultsMalign, char **readseqMalignArray, int32_t itemNumMalignArray);
short computeResultPAMalign(char *tmpResult[2], char *seq1, char *seq2, int32_t seqlen1, int32_t seqlen2, int32_t **scoreArray);
short judgeMaxOrientationMalign(int32_t *orient, int32_t **scoreArray, int32_t i, int32_t j, int32_t rowNum, int32_t colNum);
short removeRightSpaces(char **alignResultsMalign, int32_t itemNumMalignArray);
//================= msa.c declaration end ================/

#endif // METHODS_H_INCLUDED
