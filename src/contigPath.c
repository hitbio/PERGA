/*
 * contigPath.c
 *
 *  Created on: Nov 23, 2013
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Initialize the memory for the contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateContigPath(contigPath_t *contigPath, int32_t navigationFlag, kmertype *kmers[2], assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, contigtype *contigArray, int32_t itemNumContigArray)
{
	//int32_t multiCopyFlag, multiCopyUpdateFlag;

	if(contigPath==NULL)
	{
		printf("line=%d, In %s(), contigPath==NULL, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// update the navigation contig path item, 2014-01-04
	if(updateNaviContigPathItem(contigPath, contigArray[itemNumContigArray-1].base)==FAILED)
	{
		printf("line=%d, In %s(), cannot update the navigation contig path item, error!\n", __LINE__, __func__);
		return FAILED;
	}

	contigPath->startRowNewBase ++;
	contigPath->updateInterval ++;

	if(contigPath->updateInterval>=contigPath->updateIntervalThres || contigPath->itemNumPathItemList==0 || (contigPath->itemNumPathItemList>0 && contigPath->maxContigPathLen-contigPath->startRowNewBase<=0) ||
		(contigPath->itemNumPathItemList>=2 && ((contigPath->maxPathItem && contigPath->maxPathItem->contigPathLen-contigPath->startRowNewBase<=0) || (contigPath->secPathItem && contigPath->secPathItem->contigPathLen-contigPath->startRowNewBase<=0))) ||
		(contigPath->itemNumPathItemList>0 && contigPath->maxContigPathLen-contigPath->startRowNewBase<0.5*contigPath->updateIntervalThres && contigPath->updateInterval>=0.5*contigPath->updateIntervalThres))
	{ // update the contig path

		//if((*readsNumDecisionTable)>1000*averKmerOcc && contigPath->validCandPathTandPathFlagPE==YES && contigPath->candPathLenTandPathPE-contigPath->startRowCandPathTandPathPE>10) // 2014-02-07
		//{
		//	if(getContigPathSEUsingCandTandPathPE(contigPath, decisionTable, readsNumDecisionTable, dtRowHashtable, itemNumContigArray)==FAILED)
		//	{
		//		printf("line=%d, In %s(), cannot update the contig path using paired-ends, error!\n", __LINE__, __func__);
		//		return FAILED;
		//	}
		//}else
		{
//			if(navigationFlag==NAVI_PE_FLAG)
//			{ // using paired-ends
//				if(getContigPathPE(contigPath, decisionTable, readsNumDecisionTable, dtRowHashtable, dtRowHashtable)==FAILED)
//				{
//					printf("line=%d, In %s(), cannot update the contig path using paired-ends, error!\n", __LINE__, __func__);
//					return FAILED;
//				}
//			}else
//			{ // using single-ends

				if(getContigPathSE(contigPath, decisionTable, *readsNumDecisionTable, dtRowHashtable, itemNumContigArray)==FAILED)
				{
					printf("line=%d, In %s(), cannot update the contig path using paired-ends, error!\n", __LINE__, __func__);
					return FAILED;
				}
//			}
		}

		if(contigPath->maxContigPathLen>=contigPath->maxContigPathLenThres)
		{
			if(shiftStartRowNewBaseContigPath(contigPath)==FAILED)
			{
				printf("line=%d, In %s(), cannot update the contig path using paired-ends, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		contigPath->updateInterval = 0;


//		if((*readsNumDecisionTable)>1000*averKmerOcc)
//		{
//			printf("localContigID=%ld, itemNumContigArr=%ld, itemNumDecisionTable=%d\n", localContigID, itemNumContigArr, (*readsNumDecisionTable));
//			outputContigPath(contigPath, NO);
//		}

	}

	if(contigPath->contigtailSeqLen==0)
	{
		if(initTailseqContigPath(contigPath, contigArr, itemNumContigArray)==FAILED)
		{
			printf("line=%d, In %s(), cannot initialize the contig tail sequence of contig path, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		// update the contig tail sequence
		if(appendTailBaseToContigPath(contigPath, contigArr, itemNumContigArray)==FAILED)
		{
			printf("line=%d, In %s(), cannot append the tail base tp the contig path, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}


/**
 * Initialize the memory for the contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemContigPath(contigPath_t **contigPath)
{
	*contigPath = (contigPath_t*) calloc (1, sizeof(contigPath_t));
	if((*contigPath)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	(*contigPath)->updateInterval = -1;
	(*contigPath)->maxContigPathLen = 0;
	(*contigPath)->maxContigPathLenThres = MAX_READ_LEN_IN_BUF - readLen;
	(*contigPath)->startRowNewBase = -1;
	(*contigPath)->updateIntervalThres = (int32_t)(MAX_UPDATE_INTERVAL_FACTOR_CONTIGPATH * readLen);
	(*contigPath)->mismatchFactor = MAX_MISMATCHNUM_FACTOR_CONTIGPATH;
	(*contigPath)->maxMismatchNumThres = MAX_MISMATCHNUM_FACTOR_CONTIGPATH * readLen;
	(*contigPath)->overlapWithContigThres = OVERLAP_SIZE_WITH_CONTIG_FACTOR * readLen;
	(*contigPath)->itemNumPathItemList = 0;
	(*contigPath)->bestItemNumPathItemList = BEST_ITEM_NUM_CONTIGPATH;
	(*contigPath)->contigPathItemList = (*contigPath)->tailPathItem = NULL;
	(*contigPath)->maxPathItem = (*contigPath)->secPathItem = NULL;
	(*contigPath)->naviPathItem = NULL;
	(*contigPath)->naviSuccessSize = 0;
	(*contigPath)->preNaviSuccessSize = 0;
	(*contigPath)->preNaviOverlapSize = 0;

	(*contigPath)->maxContigtailSeqLen = MAX_READ_LEN_IN_BUF;
	(*contigPath)->defaultInitContigtailSeqLen = readLen;
	(*contigPath)->contigtailSeq = (char*) calloc ((*contigPath)->maxContigtailSeqLen+1, sizeof(char));
	if((*contigPath)->contigtailSeq==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Release the memory of the contig path.
 */
void releaseMemContigPath(contigPath_t **contigPath)
{
	//free((*contigPath)->contigPathItemArray);
	free((*contigPath)->contigtailSeq);
	free(*contigPath);
	*contigPath = NULL;
}

/**
 * Clean the contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short cleanContigPath(contigPath_t *contigPath)
{
	contigPathItem_t *pathItem, *nextItemTmp;
	contigPathItemRead_t *pathItemRead, *nextItemReadTmp;

	pathItem = contigPath->contigPathItemList;
	while(pathItem)
	{
		nextItemTmp = pathItem->nextPathItem;
		pathItemRead = pathItem->pathItemReadList;
		while(pathItemRead)
		{
			nextItemReadTmp = pathItemRead->nextPathItemRead;
			free(pathItemRead);
			pathItemRead = nextItemReadTmp;
		}
		free(pathItem);
		pathItem = nextItemTmp;
	}

	contigPath->contigPathItemList = NULL;
	contigPath->tailPathItem = contigPath->maxPathItem = contigPath->secPathItem = contigPath->naviPathItem = NULL;
	contigPath->naviSuccessSize = 0;
	contigPath->preNaviSuccessSize = 0;
	contigPath->maxContigPathLen = 0;
	contigPath->updateInterval = -1;
	contigPath->startRowNewBase = -1;
	contigPath->itemNumPathItemList = 0;
	contigPath->validCandPathTandPathFlagPE = NO;
	contigPath->validCandPathTandPathFlagSE = NO;
	contigPath->startRowCandPathTandPathPE = -1;
	contigPath->startRowCandPathTandPathSE = -1;
	contigPath->candPathLenTandPathPE = 0;
	contigPath->candPathLenTandPathSE = 0;
	contigPath->contigtailSeqLen = 0;
	contigPath->contigtailSeq[0] = '\0';

	return SUCCESSFUL;
}

/**
 * Get the multiCopy flag of reads in k-mers.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMultiReadCopyFlag(int32_t *multiCopyFlag, kmertype *tmp_kmers[2])
{
	int32_t i, posNum, pos;
	ridpostype *ridpostable;

	*multiCopyFlag = NO;
	if(tmp_kmers[0])
	{
		ridpostable = tmp_kmers[0]->ppos;
		posNum = tmp_kmers[0]->arraysize;
		for(i=0; i<posNum-1; i++)
		{
			if(ridpostable[i].rid==ridpostable[i+1].rid)
			{
				*multiCopyFlag = YES;
				break;
			}
		}
	}

	if((*multiCopyFlag)==NO)
	{
		if(tmp_kmers[1])
		{
			ridpostable = tmp_kmers[1]->ppos;
			posNum = tmp_kmers[1]->arraysize;
			for(i=0; i<posNum-1; i++)
			{
				if(ridpostable[i].rid==ridpostable[i+1].rid)
				{
					*multiCopyFlag = YES;
					break;
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the contig path using paired-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getContigPathPE(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, int32_t itemNumContigArray)
{
	char oldNaviPathseq[MAX_READ_LEN_IN_BUF+1];
	int32_t oldNaviPathseqLen, useOldNaviPathseqFlag;

	if(contigPath->naviPathItem && contigPath->naviPathItem->supportReadsNum>0)
	{
		useOldNaviPathseqFlag = YES;
		strcpy(oldNaviPathseq, contigPath->naviPathItem->contigPathStr+contigPath->startRowNewBase);
		oldNaviPathseqLen = contigPath->naviPathItem->contigPathLen - contigPath->startRowNewBase;
	}else
	{
		useOldNaviPathseqFlag = NO;
		oldNaviPathseqLen = 0;
	}

	// remove the path item having no reads
	if(removeNoReadsContigPathItem(contigPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove the path item having no reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(contigPath->itemNumPathItemList>=2)
	{
		// merge contig path items
		if(mergeContigPathItem(contigPath, decisionTable, dtRowHashtable)==FAILED)
		{
			printf("line=%d, In %s(), cannot merge contig path items, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// get the contig path using paired-ends
	if(getContigPathFromDecisionTablePE(contigPath, decisionTable, readsNumDecisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot update the contig path using paired-ends, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ####################### Debug information ########################
	if(maxItemNumContigPath<contigPath->itemNumPathItemList)
	{
		maxItemNumContigPath = contigPath->itemNumPathItemList;
	}
	// ####################### Debug information ########################

	// adjust the contig path items
	if(adjustContigPath(contigPath, decisionTable, dtRowHashtable)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust the contig path, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// remain the main path items
	if(removeLessSupportedContigPathItems(contigPath, useOldNaviPathseqFlag, oldNaviPathseq, oldNaviPathseqLen, decisionTable, dtRowHashtable)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove the less supported contig path, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// remove short overlapped contig path items
	if(removeShortOverlappedContigPathItems(contigPath, decisionTable, dtRowHashtable)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove the less supported contig path, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// sort path items
	if(sortContigPathItems(contigPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot sort the contig path, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// update the naviPathItem
	if(setNaviContigPathItem(contigPath, useOldNaviPathseqFlag, oldNaviPathseq, oldNaviPathseqLen)==FAILED)
	{
		printf("line=%d, In %s(), cannot update navigation contig path item, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ####################### Debug information ########################
	if(maxItemNumContigPathAdjusted<contigPath->itemNumPathItemList)
	{
		maxItemNumContigPathAdjusted = contigPath->itemNumPathItemList;
	}
	// ####################### Debug information ########################

	// remain the best path items
//	if(remainBestContigPathItems(contigPath, decisionTable, dtRowHashtable)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot remove the less supported contig path, error!\n", __LINE__, __func__);
//		return FAILED;
//	}

	// get the maximum and second maximum of the contig path items
//	if(getMaxesContigPathItems(contigPath)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot sort the contig path, error!\n", __LINE__, __func__);
//		return FAILED;
//	}

	return SUCCESSFUL;
}

/**
 * Get the contig path using single-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getContigPathSEUsingCandTandPathPE(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, int32_t itemNumContigArray)
{
	char oldNaviPathseq[MAX_READ_LEN_IN_BUF+1];
	int32_t oldNaviPathseqLen, useOldNaviPathseqFlag;

	if(contigPath->naviPathItem && contigPath->naviPathItem->supportReadsNum>0)
	{
		useOldNaviPathseqFlag = YES;
		strcpy(oldNaviPathseq, contigPath->naviPathItem->contigPathStr+contigPath->startRowNewBase);
		oldNaviPathseqLen = contigPath->naviPathItem->contigPathLen - contigPath->startRowNewBase;
	}else
	{
		useOldNaviPathseqFlag = NO;
		oldNaviPathseqLen = 0;
	}

	// remove the path item having no reads
	if(removeNoReadsContigPathItem(contigPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove the path item having no reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(contigPath->itemNumPathItemList>=2)
	{
		// merge contig path items
		if(mergeContigPathItem(contigPath, decisionTable, dtRowHashtable)==FAILED)
		{
			printf("line=%d, In %s(), cannot merge contig path items, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// remove unmatched path items
	if(removeUnmatchedContigPathItem(contigPath, decisionTable, readsNumDecisionTable, dtRowHashtable)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove the unmatched contig path items, error!\n", __LINE__, __func__);
		return FAILED;
	}

//	// reset the maximal overlap size with contig
//	if(resetMaxOverlapSizeWithContigContigPathItem(contigPath, decisionTable, dtRowHashtable)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot reset the maximal overlap size of contig path items, error!\n", __LINE__, __func__);
//		return FAILED;
//	}

	// added the path item according to candTandPathPE
	if(getContigPathSEFromCandTandPathPE(contigPath, decisionTable, *readsNumDecisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot update the contig path using single-ends, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// adjust the contig path items
	if(adjustContigPath(contigPath, decisionTable, dtRowHashtable)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust contig path, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// sort path items
	if(sortContigPathItems(contigPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot sort the contig path, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// remove reads of identical part of path items
	if(removeReadsInIdenticalPathItemPart(contigPath, useOldNaviPathseqFlag, oldNaviPathseq, oldNaviPathseqLen, decisionTable, *readsNumDecisionTable, dtRowHashtable, itemNumContigArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove the reads in identical part of path items, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// sort path items
	if(sortContigPathItems(contigPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot sort the contig path, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// update the naviPathItem
	if(setNaviContigPathItem(contigPath, useOldNaviPathseqFlag, oldNaviPathseq, oldNaviPathseqLen)==FAILED)
	{
		printf("line=%d, In %s(), cannot update navigation contig path item, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get the contig path using single-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getContigPathSE(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, int32_t itemNumContigArray)
{
	char oldNaviPathseq[MAX_READ_LEN_IN_BUF+1];
	int32_t oldNaviPathseqLen, useOldNaviPathseqFlag;

	if(contigPath->naviPathItem && contigPath->naviPathItem->supportReadsNum>0)
	{
		useOldNaviPathseqFlag = YES;
		strcpy(oldNaviPathseq, contigPath->naviPathItem->contigPathStr+contigPath->startRowNewBase);
		oldNaviPathseqLen = contigPath->naviPathItem->contigPathLen - contigPath->startRowNewBase;
	}else
	{
		useOldNaviPathseqFlag = NO;
		oldNaviPathseqLen = 0;
	}

	// remove the path item having no reads
	if(removeNoReadsContigPathItem(contigPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove the path item having no reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(contigPath->itemNumPathItemList>=2)
	{
		// merge contig path items
		if(mergeContigPathItem(contigPath, decisionTable, dtRowHashtable)==FAILED)
		{
			printf("line=%d, In %s(), cannot merge contig path items, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

//	// reset the maximal overlap size with contig
//	if(resetMaxOverlapSizeWithContigContigPathItem(contigPath, decisionTable, dtRowHashtable)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot reset the maximal overlap size of contig path items, error!\n", __LINE__, __func__);
//		return FAILED;
//	}

	// get the path items
	if(getContigPathFromDecisionTableSE(contigPath, decisionTable, readsNumDecisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot update the contig path using single-ends, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ####################### Debug information ########################
	if(maxItemNumContigPath<contigPath->itemNumPathItemList)
	{
		maxItemNumContigPath = contigPath->itemNumPathItemList;
	}
	// ####################### Debug information ########################

	// adjust the contig path items
	if(adjustContigPath(contigPath, decisionTable, dtRowHashtable)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust contig path, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// remain the main path items
	if(removeLessSupportedContigPathItems(contigPath, useOldNaviPathseqFlag, oldNaviPathseq, oldNaviPathseqLen, decisionTable, dtRowHashtable)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove the less supported contig path, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// remove contig path items with short fragment size
//	if(removeShortFragSizeContigPath(contigPath, decisionTable, dtRowHashtable)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot remove reads in tandem paths, error!\n", __LINE__, __func__);
//		return FAILED;
//	}

	// remove short overlapped contig path items, 2014-01-04
	if(removeShortOverlappedContigPathItems(contigPath, decisionTable, dtRowHashtable)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove the less supported contig path, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// sort path items
	if(sortContigPathItems(contigPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot sort the contig path, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// remove reads of identical part of path items
	if(removeReadsInIdenticalPathItemPart(contigPath, useOldNaviPathseqFlag, oldNaviPathseq, oldNaviPathseqLen, decisionTable, readsNumDecisionTable, dtRowHashtable, itemNumContigArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove the reads in identical part of path items, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// sort path items
	if(sortContigPathItems(contigPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot sort the contig path, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// update the naviPathItem
	if(setNaviContigPathItem(contigPath, useOldNaviPathseqFlag, oldNaviPathseq, oldNaviPathseqLen)==FAILED)
	{
		printf("line=%d, In %s(), cannot update navigation contig path item, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ####################### Debug information ########################
	if(maxItemNumContigPathAdjusted<contigPath->itemNumPathItemList)
	{
		maxItemNumContigPathAdjusted = contigPath->itemNumPathItemList;
	}
	// ####################### Debug information ########################

	// remain the best path items
//	if(remainBestContigPathItems(contigPath, decisionTable, dtRowHashtable)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot remove the less supported contig path, error!\n", __LINE__, __func__);
//		return FAILED;
//	}


	// get the maximal and second maximal contig path items
//	if(getMaxesContigPathItems(contigPath)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot sort the contig path, error!\n", __LINE__, __func__);
//		return FAILED;
//	}

	return SUCCESSFUL;
}

/**
 * Get the contig path from decision table using paired-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getContigPathFromDecisionTablePE(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable)
{
	int32_t i, newStartPos, readseqLenTmp, overlapSize;
	char readseqTmp[MAX_READ_LEN_IN_BUF+1];
	contigPathItem_t *targetPathItem;

	for(i=0; i<readsNumDecisionTable; i++)
	{
		if(decisionTable[i].matedFlag==YES
			&& (decisionTable[i].unmatchBaseNum<=maxUnmatchBaseNumPerRead)
			&& (decisionTable[i].successiveAppearBases>=minSuccessiveAppearedBaseNum)
			&& decisionTable[i].basePos==decisionTable[i].lastMatchedBasePos)
		{
			if(decisionTable[i].contigPathItem==NULL)
			{
				// get the read base sequence
				if(decisionTable[i].orientation==ORIENTATION_PLUS)
				{
					newStartPos = decisionTable[i].basePos + 1;
					readseqLenTmp = decisionTable[i].seqlen - newStartPos;
					getReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, newStartPos, readseqLenTmp);
					overlapSize = newStartPos;
				}else
				{
					newStartPos = decisionTable[i].seqlen - decisionTable[i].basePos;
					readseqLenTmp = decisionTable[i].seqlen - newStartPos;
					getReverseReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, newStartPos, readseqLenTmp);
					overlapSize = newStartPos;
				}

				if(readseqLenTmp>0)
				{
					// get the match row
					if(getMatchedContigPathItem(&targetPathItem, readseqTmp, readseqLenTmp, contigPath)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the matched row in contig path, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(targetPathItem>=0)
					{
						// add the read to the matched contig path item
						if(addReadseqToContigPathItem(readseqTmp, readseqLenTmp, overlapSize, contigPath, targetPathItem, decisionTable+i)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the matched row in contig path, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}else
					{
						// add the new path item
						if(addNewContigPathItem(contigPath, readseqTmp, readseqLenTmp, overlapSize, decisionTable+i)==FAILED)
						{
							printf("line=%d, In %s(), cannot add the contig path item, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
				}
			}
		}else
		{
			if(decisionTable[i].contigPathItem)
			{
				if(delReadFromContigPath(decisionTable+i, contigPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot delete read from contig path item, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Reset the maximal overlap size with contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short resetMaxOverlapSizeWithContigContigPathItem(contigPath_t *contigPath, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable)
{
	int32_t i, newStartPos, readseqLenTmp, overlapSize;
	contigPathItem_t *pathItem;
	contigPathItemRead_t *pathItemRead;
	assemblingreadtype *dtRead;

	pathItem = contigPath->contigPathItemList;
	while(pathItem)
	{
		pathItem->maxOverlapWithContig = 0;

		pathItemRead = pathItem->pathItemReadList;
		while(pathItemRead)
		{
			// get the exist read
			if(getExistReadInDT(&dtRead, pathItemRead->rid, decisionTable, dtRowHashtable)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the exist read %ld in decision table, error!\n", __LINE__, __func__, pathItemRead->rid);
				return FAILED;
			}

			if(dtRead->orientation==ORIENTATION_PLUS)
				overlapSize = dtRead->basePos + 1;
			else
				overlapSize = dtRead->seqlen - dtRead->basePos;

			if(overlapSize>pathItem->maxOverlapWithContig)
				pathItem->maxOverlapWithContig = overlapSize;

			pathItemRead = pathItemRead->nextPathItemRead;
		}

		pathItem = pathItem->nextPathItem;
	}

	return SUCCESSFUL;
}

/**
 * Get the contig path from decision table using single-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getContigPathFromDecisionTableSE(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable)
{
	int32_t i, newStartPos, readseqLenTmp, overlapSize;
	char readseqTmp[MAX_READ_LEN_IN_BUF+1];
	contigPathItem_t *targetPathItem;

	for(i=0; i<readsNumDecisionTable; i++)
	{
		if((decisionTable[i].unmatchBaseNum<=maxUnmatchBaseNumPerRead)
			//(decisionTable[i].unmatchBaseNum<=1)
			&& (decisionTable[i].successiveAppearBases>=minSuccessiveAppearedBaseNum)
			&& decisionTable[i].basePos==decisionTable[i].lastMatchedBasePos)
		{
			if(decisionTable[i].contigPathItem==NULL)
			{
//				if(decisionTable[i].rid==487291)
//				{
//					printf("line=%d, In %s(), rid=%ld, firstContigPos=%d\n", __LINE__, __func__, (int64_t)decisionTable[i].rid, decisionTable[i].firstContigPos);
//				}

				// get the read base sequence
				if(decisionTable[i].orientation==ORIENTATION_PLUS)
				{
					newStartPos = decisionTable[i].basePos + 1;
					readseqLenTmp = decisionTable[i].seqlen - newStartPos;
					getReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, newStartPos, readseqLenTmp);
					overlapSize = newStartPos;
				}else
				{
					newStartPos = decisionTable[i].seqlen - decisionTable[i].basePos;
					readseqLenTmp = decisionTable[i].seqlen - newStartPos;
					getReverseReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, newStartPos, readseqLenTmp);
					overlapSize = newStartPos;
				}

				if(readseqLenTmp>0)
				{
					// get the match row
					if(getMatchedContigPathItem(&targetPathItem, readseqTmp, readseqLenTmp, contigPath)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the matched row in contig path, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(targetPathItem)
					{
						// add the read to the matched contig path item
						if(addReadseqToContigPathItem(readseqTmp, readseqLenTmp, overlapSize, contigPath, targetPathItem, decisionTable+i)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the matched row in contig path, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}else
					{
						// add the new path item
						if(addNewContigPathItem(contigPath, readseqTmp, readseqLenTmp, overlapSize, decisionTable+i)==FAILED)
						{
							printf("line=%d, In %s(), cannot add the contig path item, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
				}
			}
		}else
		{
			if(decisionTable[i].contigPathItem)
			{
				if(delReadFromContigPath(decisionTable+i, contigPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot delete read from contig path item, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}

	// remove the path item having no reads
//	if(removeNoReadsContigPathItem(contigPath)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot remove the path item having no reads, error!\n", __LINE__, __func__);
//		return FAILED;
//	}

	return SUCCESSFUL;
}

/**
 * Get the path items using candTandPathPE.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getContigPathSEFromCandTandPathPE(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable)
{
	int32_t i, j, newStartPos, readseqLenTmp, overlapSize, pathLenCandTandPath, shareLen, mismatchNum;
	char readseqTmp[MAX_READ_LEN_IN_BUF+1], *pathseqCandTandPath;
	contigPathItem_t *targetPathItem;

	pathseqCandTandPath = contigPath->candPathseqTandPathPE + contigPath->startRowCandPathTandPathPE;
	pathLenCandTandPath = contigPath->candPathLenTandPathPE - contigPath->startRowCandPathTandPathPE;

	for(i=0; i<readsNumDecisionTable; i++)
	{
		//if(decisionTable[i].rid==17344030)
		//{
		//	printf("rid=%ld\n", (int64_t)decisionTable[i].rid);
		//}

		if((decisionTable[i].unmatchBaseNum<=maxUnmatchBaseNumPerRead)
			//(decisionTable[i].unmatchBaseNum<=1)
			&& (decisionTable[i].successiveAppearBases>=minSuccessiveAppearedBaseNum)
			&& decisionTable[i].basePos==decisionTable[i].lastMatchedBasePos)
		{
			if(decisionTable[i].contigPathItem==NULL)
			{
				// get the read base sequence
				if(decisionTable[i].orientation==ORIENTATION_PLUS)
				{
					newStartPos = decisionTable[i].basePos + 1;
					readseqLenTmp = decisionTable[i].seqlen - newStartPos;
					getReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, newStartPos, readseqLenTmp);
					overlapSize = newStartPos;
				}else
				{
					newStartPos = decisionTable[i].seqlen - decisionTable[i].basePos;
					readseqLenTmp = decisionTable[i].seqlen - newStartPos;
					getReverseReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, newStartPos, readseqLenTmp);
					overlapSize = newStartPos;
				}

				if(readseqLenTmp>0)
				{
					// get the match row
					if(getMatchedContigPathItem(&targetPathItem, readseqTmp, readseqLenTmp, contigPath)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the matched row in contig path, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(targetPathItem)
					{
						// add the read to the matched contig path item
						if(addReadseqToContigPathItem(readseqTmp, readseqLenTmp, overlapSize, contigPath, targetPathItem, decisionTable+i)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the matched row in contig path, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}else
					{
						if(readseqLenTmp<pathLenCandTandPath)
							shareLen = readseqLenTmp;
						else
							shareLen = pathLenCandTandPath;

						if(shareLen>0)
						{
							mismatchNum = 0;
							for(j=0; j<shareLen; j++)
							{
								if(readseqTmp[j]!=pathseqCandTandPath[j])
								{
									mismatchNum ++;
									if(mismatchNum>3)
										break;
								}
							}
						}else
						{
							mismatchNum = INT_MAX;
						}

						if(mismatchNum<=3)
						//if(mismatchNum<=3 && shareLen>10*mismatchNum)
						{
							// add the new path item
							if(addNewContigPathItem(contigPath, readseqTmp, readseqLenTmp, overlapSize, decisionTable+i)==FAILED)
							{
								printf("line=%d, In %s(), cannot add the contig path item, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}
					}
				}
			}
		}else
		{
			if(decisionTable[i].contigPathItem)
			{
				if(delReadFromContigPath(decisionTable+i, contigPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot delete read from contig path item, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}

	// remove the path item having no reads
//	if(removeNoReadsContigPathItem(contigPath)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot remove the path item having no reads, error!\n", __LINE__, __func__);
//		return FAILED;
//	}

	return SUCCESSFUL;
}

/**
 * Get the match row in contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMatchedContigPathItem(contigPathItem_t **matchedPathItem, char *readseqTmp, int32_t readseqLenTmp, contigPath_t *contigPath)
{
	int32_t j, misMatchNum, minMismatchNum, pathLen, shareLen;
	char *pathseq;
	contigPathItem_t *pathItem;

	minMismatchNum = INT_MAX;
	pathItem = contigPath->contigPathItemList;
	while(pathItem)
	{
		if(pathItem->supportReadsNum>0)
		{
			pathseq = pathItem->contigPathStr + contigPath->startRowNewBase;
			pathLen = pathItem->contigPathLen - contigPath->startRowNewBase;
			if(readseqLenTmp<pathLen)
				shareLen = readseqLenTmp;
			else
				shareLen = pathLen;

			if(shareLen>0)
			{
				misMatchNum = 0;
				for(j=0; j<shareLen; j++)
				{
					if(readseqTmp[j]!=pathseq[j])
					{
						misMatchNum ++;
						//if(misMatchNum>contigPath->maxMismatchNumThres/2)
						if(misMatchNum>0)
							break;
					}
				}

				if(misMatchNum<minMismatchNum)
				{
					minMismatchNum = misMatchNum;
					*matchedPathItem = pathItem;
				}
			}

			if(minMismatchNum==0)
				break;
		}

		pathItem = pathItem->nextPathItem;
	}

	//if(minMismatchNum>contigPath->maxMismatchNumThres/2)
	if(minMismatchNum>0)
		*matchedPathItem = NULL;

	return SUCCESSFUL;
}


/**
 * Add the read to the matched contig path item.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addReadseqToContigPathItem(char *readseqTmp, int32_t readseqLenTmp, int32_t overlapSize, contigPath_t *contigPath, contigPathItem_t *targetPathItem, assemblingreadtype *dtRead)
{
	int32_t pathLen, newBaseLen;
	char *pathseq;
	contigPathItemRead_t *pathItemRead;

	pathLen = targetPathItem->contigPathLen - contigPath->startRowNewBase;
	if(pathLen<readseqLenTmp)
	{ // add the new bases to the path
		newBaseLen = readseqLenTmp - pathLen;
		pathseq = targetPathItem->contigPathStr + contigPath->startRowNewBase;
		strcat(pathseq+pathLen, readseqTmp+pathLen);
		targetPathItem->contigPathLen += newBaseLen;
	}

	// add the read information
	pathItemRead = (contigPathItemRead_t*) calloc(1, sizeof(contigPathItemRead_t));
	if(pathItemRead==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	pathItemRead->rid = dtRead->rid;
	pathItemRead->prePathItemRead = targetPathItem->tailPathItemRead;
	pathItemRead->nextPathItemRead = NULL;
	targetPathItem->tailPathItemRead->nextPathItemRead = pathItemRead;
	targetPathItem->tailPathItemRead = pathItemRead;
	targetPathItem->supportReadsNum ++;

	if(targetPathItem->maxOverlapWithContig<overlapSize)
		targetPathItem->maxOverlapWithContig = overlapSize;

	if(contigPath->maxContigPathLen<targetPathItem->contigPathLen)
		contigPath->maxContigPathLen = targetPathItem->contigPathLen;

	dtRead->contigPathItem = targetPathItem;
	dtRead->contigPathItemRead = pathItemRead;

	return SUCCESSFUL;
}

/**
 * Add the read to the new contig path item.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewContigPathItem(contigPath_t *contigPath, char *readseqTmp, int32_t readseqLenTmp, int32_t overlapSize, assemblingreadtype *dtRead)
{
	contigPathItem_t *pathItem;
	contigPathItemRead_t *pathItemRead;

	pathItem = (contigPathItem_t*) calloc(1, sizeof(contigPathItem_t));
	if(pathItem==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	strcpy(pathItem->contigPathStr+contigPath->startRowNewBase, readseqTmp);
	pathItem->contigPathLen = readseqLenTmp + contigPath->startRowNewBase;
	pathItem->prePathItem = contigPath->tailPathItem;
	pathItem->nextPathItem = NULL;

	pathItemRead = (contigPathItemRead_t*) calloc(1, sizeof(contigPathItemRead_t));
	if(pathItemRead==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	pathItemRead->rid = dtRead->rid;
	pathItemRead->prePathItemRead = NULL;
	pathItemRead->nextPathItemRead = NULL;
	pathItem->pathItemReadList = pathItemRead;
	pathItem->tailPathItemRead = pathItemRead;
	pathItem->supportReadsNum = 1;
	pathItem->maxOverlapWithContig = overlapSize;

	if(contigPath->maxContigPathLen<pathItem->contigPathLen)
		contigPath->maxContigPathLen = pathItem->contigPathLen;

	if(contigPath->tailPathItem)
		contigPath->tailPathItem->nextPathItem = pathItem;
	else
		contigPath->contigPathItemList = pathItem;
	contigPath->tailPathItem = pathItem;
	contigPath->itemNumPathItemList ++;

	dtRead->contigPathItem = pathItem;
	dtRead->contigPathItemRead = pathItemRead;

	return SUCCESSFUL;
}

/**
 * Delete a read from the contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short delReadFromContigPath(assemblingreadtype *dtRead, contigPath_t *contigPath)
{
	contigPathItem_t *contigPathItem;
	contigPathItemRead_t *contigPathItemRead;

	contigPathItem = dtRead->contigPathItem;
	contigPathItemRead = dtRead->contigPathItemRead;

	if(contigPathItem)
	{
		if(contigPathItemRead->rid!=dtRead->rid)
		{
			printf("line=%d, In %s(), invalid contigPathItemRead->rid=%ld, rid=%ld, error!\n", __LINE__, __func__, (int64_t)contigPathItemRead->rid, (int64_t)dtRead->rid);
			return FAILED;
		}

		// delete the read information from contig path
		if(contigPathItemRead->prePathItemRead)
		{ // not the head
			contigPathItemRead->prePathItemRead->nextPathItemRead = contigPathItemRead->nextPathItemRead;
		}else
		{ // the head
			contigPathItem->pathItemReadList = contigPathItemRead->nextPathItemRead;
		}

		// process the tail read
		if(contigPathItemRead->nextPathItemRead)
		{ // not the tail
			contigPathItemRead->nextPathItemRead->prePathItemRead = contigPathItemRead->prePathItemRead;
		}else
		{ // the tail
			contigPathItem->tailPathItemRead = contigPathItemRead->prePathItemRead;
		}

		// free the read
		free(contigPathItemRead);
		contigPathItem->supportReadsNum --;

/*
		if(contigPathItem->supportReadsNum==0)
		{ // remove the item node
			if(contigPathItem->prePathItem)
			{ // not the head
				contigPathItem->prePathItem->nextPathItem = contigPathItem->nextPathItem;
			}else
			{ // the head
				contigPath->contigPathItemList = contigPathItem->nextPathItem;
			}

			if(contigPathItem->nextPathItem)
			{ // not the tail
				contigPathItem->nextPathItem->prePathItem = contigPathItem->prePathItem;
			}else
			{ // the tail
				contigPath->tailPathItem = contigPathItem->prePathItem;
			}

			if(contigPath->maxPathItem==contigPathItem || contigPath->secPathItem==contigPathItem)
			{
				// get the maximal and second maximal contig path items
				if(getMaxesContigPathItems(contigPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot sort the contig path, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}

			// free the item node
			free(contigPathItem);
			contigPath->itemNumPathItemList --;
		}
*/

		dtRead->contigPathItem = NULL;
		dtRead->contigPathItemRead = NULL;
	}

	return SUCCESSFUL;
}

/**
 * Adjust the contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short adjustContigPath(contigPath_t *contigPath, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable)
{
	int32_t *baseNumArray, maxRowNumBaseNumArray, colsNum;  // column=8

	if(contigPath->itemNumPathItemList>=2)
	{
		// initialize the memory
		colsNum = 8;
		maxRowNumBaseNumArray = contigPath->maxContigPathLen - contigPath->startRowNewBase;
		baseNumArray = (int32_t *)calloc(maxRowNumBaseNumArray*colsNum, sizeof(int32_t));
		if(baseNumArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// fill the baseNumArray
		if(computeBaseNumArrayContigPath(baseNumArray, maxRowNumBaseNumArray, colsNum, contigPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the base number array for contig path, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// remove error base in contig path and update the path
		if(removeErrorBaseContigPath(baseNumArray, maxRowNumBaseNumArray, colsNum, contigPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot remove error bases for contig path, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// merge adjusted contig path items
		if(mergeContigPathItem(contigPath, decisionTable, dtRowHashtable)==FAILED)
		{
			printf("line=%d, In %s(), cannot merge contig path, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// free memory
		free(baseNumArray);
	}

	return SUCCESSFUL;
}

/**
 * Remove the contig path item having no reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeNoReadsContigPathItem(contigPath_t *contigPath)
{
	contigPathItem_t *pathItem, *pathItemTmp;

	pathItem = contigPath->contigPathItemList;
	while(pathItem)
	{
		if(pathItem->supportReadsNum==0)
		{
			if(pathItem->prePathItem)
			{ // not the head
				pathItem->prePathItem->nextPathItem = pathItem->nextPathItem;
			}else
			{ // the head
				contigPath->contigPathItemList = pathItem->nextPathItem;
			}

			if(pathItem->nextPathItem)
			{ // not the tail
				pathItem->nextPathItem->prePathItem = pathItem->prePathItem;
			}else
			{ // the tail
				contigPath->tailPathItem = pathItem->prePathItem;
			}

			pathItemTmp = pathItem->nextPathItem;
			free(pathItem);
			pathItem = pathItemTmp;

			contigPath->itemNumPathItemList --;
		}else
		{
			pathItem = pathItem->nextPathItem;
		}
	}

	return SUCCESSFUL;
}

/**
 * Compute the base number array for contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeBaseNumArrayContigPath(int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, contigPath_t *contigPath)
{
	int32_t i, j, baseNumTmp, supportReadsNum, baseInt;
	int32_t value, maxValue, secValue, maxIndex, secIndex, subSum;
	char *pathseq;
	contigPathItem_t *pathItem;

	// fill the base number
	pathItem = contigPath->contigPathItemList;
	while(pathItem)
	{
		if(pathItem->supportReadsNum>0)
		{
			supportReadsNum = pathItem->supportReadsNum;
			baseNumTmp = pathItem->contigPathLen - contigPath->startRowNewBase;
			pathseq = pathItem->contigPathStr + contigPath->startRowNewBase;
			for(j=0; j<baseNumTmp; j++)
			{
				switch(pathseq[j])
				{
					case 'A':
					case 'a': baseInt = 0; break;
					case 'C':
					case 'c': baseInt = 1; break;
					case 'G':
					case 'g': baseInt = 2; break;
					case 'T':
					case 't': baseInt = 3; break;
					case '-': baseInt = 4; break;
					default: printf("line=%d, invalid base %d\n, error!\n", __LINE__, pathseq[j]); return FAILED;
				}

				baseNumArray[j*colsNum + baseInt] += supportReadsNum; // for A, C, G, T, -
				baseNumArray[j*colsNum + 5] += supportReadsNum; 	 // for A+C+G+T+-
			}
		}

		pathItem = pathItem->nextPathItem;
	}

	// get the maxOcc and secOcc
	for(i=0; i<maxRowNumBaseNumArray; i++)
	{
		maxValue = secValue = -1;
		maxIndex = secIndex = -1;
		for(j=0; j<5; j++)
		{
			value = baseNumArray[i*colsNum + j];
			if(value>maxValue)
			{
				secValue = maxValue;
				secIndex = maxIndex;
				maxValue = value;
				maxIndex = j;
			}else if(value>secValue)
			{
				secValue = value;
				secIndex = j;
			}
		}

		subSum = baseNumArray[i*colsNum + 5];

		baseNumArray[i*colsNum + 6] = maxIndex;
		baseNumArray[i*colsNum + 7] = secIndex;
	}

	return SUCCESSFUL;
}

/**
 * Remove error bases for contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeErrorBaseContigPath(int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, contigPath_t *contigPath)
{
	int32_t i, k, pathLen;
	int32_t maxBaseIndex, maxValue, secondBaseIndex, secondValue, sum;
	int32_t gapValue, baseValue, errrBaseIndex, pathseqRow, tmpRow;
	char *pathseq, errBase, newBase;;
	contigPathItem_t *pathItem;

	for(i=0; i<maxRowNumBaseNumArray; i++)
	{
		secondBaseIndex = baseNumArray[i*colsNum + 7];
		if(secondBaseIndex<4)
		{
			maxBaseIndex = baseNumArray[i*colsNum + 6];
			maxValue = baseNumArray[i*colsNum + maxBaseIndex];
			secondValue =baseNumArray[i*colsNum + secondBaseIndex];
			gapValue = baseNumArray[i*colsNum + 4];
			sum = baseNumArray[i*colsNum + 5];
			if(secondValue==1 && gapValue==0 && sum>=10)
			{
				for(k=0; k<4; k++)
				{
					baseValue = baseNumArray[i*colsNum + k];
					if(baseValue==1)
					{
						errrBaseIndex = k;
						switch(errrBaseIndex)
						{
							case 0: errBase = 'A'; break;
							case 1: errBase = 'C'; break;
							case 2: errBase = 'G'; break;
							case 3: errBase = 'T'; break;
						}

						switch(maxBaseIndex)
						{
							case 0: newBase = 'A'; break;
							case 1: newBase = 'C'; break;
							case 2: newBase = 'G'; break;
							case 3: newBase = 'T'; break;
						}

						tmpRow = i;

						// get the single-path according to the single base
						pathItem = contigPath->contigPathItemList;
						while(pathItem)
						{
							if(pathItem->supportReadsNum>0)
							{
								pathseq = pathItem->contigPathStr + contigPath->startRowNewBase;
								pathLen = pathItem->contigPathLen - contigPath->startRowNewBase;
								pathseqRow = tmpRow;
								if(pathseqRow>=0 && pathseqRow<pathLen)
								{
									if(pathseq[pathseqRow]==errBase)
									{
										pathseq[pathseqRow] = newBase;
										break;
									}
								}
							}

							pathItem = pathItem->nextPathItem;
						}
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Merge adjusted contig path items.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mergeContigPathItem(contigPath_t *contigPath, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable)
{
	int32_t i, pathseqLen, targetseqLen, similarFlag, minLen;
	char *pathseq, *targetseq;
	contigPathItem_t *targetItem, *pathItem, *pathItemTmp;
	contigPathItemRead_t *pathItemRead;

	targetItem = contigPath->contigPathItemList;
	while(targetItem)
	{
		if(targetItem->supportReadsNum>0)
		{
			pathItem = targetItem->nextPathItem;
			while(pathItem)
			{
				if(pathItem->supportReadsNum>0)
				{
					targetseq = targetItem->contigPathStr + contigPath->startRowNewBase;
					targetseqLen = targetItem->contigPathLen - contigPath->startRowNewBase;
					pathseq = pathItem->contigPathStr + contigPath->startRowNewBase;
					pathseqLen = pathItem->contigPathLen - contigPath->startRowNewBase;

					if(pathseqLen>targetseqLen)
						minLen = targetseqLen;
					else
						minLen = pathseqLen;

					if(minLen>0)
					{
						similarFlag = YES;
						for(i=0; i<minLen; i++)
						{
							if(pathseq[i]!=targetseq[i])
							{
								similarFlag = NO;
								break;
							}
						}
					}else
					{
						similarFlag = NO;
					}

					if(similarFlag==YES)
					{
						// update the reads in decision table
						pathItemRead = pathItem->pathItemReadList;
						while(pathItemRead)
						{
							if(replacePathReadInfoInDT(targetItem, pathItem, pathItemRead, decisionTable, dtRowHashtable)==FAILED)
							{
								printf("line=%d, In %s(), cannot replace the path read information in decision table, error!\n", __LINE__, __func__);
								return FAILED;
							}

							pathItemRead = pathItemRead->nextPathItemRead;
						}
						targetItem->supportReadsNum += pathItem->supportReadsNum;

						// update the pathItem information
						if(targetseqLen<pathseqLen)
						{
							strcat(targetseq, pathseq+targetseqLen);
							targetItem->contigPathLen += pathseqLen - targetseqLen;
						}

						targetItem->tailPathItemRead->nextPathItemRead = pathItem->pathItemReadList;
						pathItem->pathItemReadList->prePathItemRead = targetItem->tailPathItemRead;
						targetItem->tailPathItemRead = pathItem->tailPathItemRead;

						// update pathItem list
						if(pathItem->prePathItem)
							pathItem->prePathItem->nextPathItem = pathItem->nextPathItem;
						else
							contigPath->contigPathItemList = pathItem->nextPathItem;
						if(pathItem->nextPathItem)
							pathItem->nextPathItem->prePathItem = pathItem->prePathItem;
						else
							contigPath->tailPathItem = pathItem->prePathItem;

						// free pathItem
						pathItemTmp = pathItem->nextPathItem;
						free(pathItem);
						pathItem = pathItemTmp;

						contigPath->itemNumPathItemList --;
					}else
						pathItem = pathItem->nextPathItem;
				}else
					pathItem = pathItem->nextPathItem;
			}
		}

		targetItem = targetItem->nextPathItem;
	}

	return SUCCESSFUL;
}

/**
 * Remove unmatched contig path items.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeUnmatchedContigPathItem(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable)
{
	int32_t i, pathLenCandTandPath, pathLen, shareLen, mismatchNum, delDTReadFlag;
	char *pathseqCandTandPath, *pathseq;
	contigPathItem_t *pathItem, *pathItemTmp;
	contigPathItemRead_t *pathItemRead, *pathItemReadTmp;

	pathseqCandTandPath = contigPath->candPathseqTandPathPE + contigPath->startRowCandPathTandPathPE;
	pathLenCandTandPath = contigPath->candPathLenTandPathPE - contigPath->startRowCandPathTandPathPE;

	delDTReadFlag = YES;
	pathItem = contigPath->contigPathItemList;
	while(pathItem)
	{
		pathseq = pathItem->contigPathStr + contigPath->startRowNewBase;
		pathLen = pathItem->contigPathLen - contigPath->startRowNewBase;
		if(pathLen<pathLenCandTandPath)
			shareLen = pathLen;
		else
			shareLen = pathLenCandTandPath;

		if(shareLen>0)
		{
			mismatchNum = 0;
			for(i=0; i<shareLen; i++)
				if(pathseq[i]!=pathseqCandTandPath[i])
					mismatchNum ++;
		}else
		{
			mismatchNum = INT_MAX;
		}

		if(mismatchNum>3)
		{
			// delete the path information of read in decision table
			pathItemRead = pathItem->pathItemReadList;
			while(pathItemRead)
			{
				if(removePathReadInfoInDT(pathItem, pathItemRead, decisionTable, dtRowHashtable, delDTReadFlag)==FAILED)
				{
					printf("line=%d, In %s(), cannot remove the path read information in decision table, error!\n", __LINE__, __func__);
					return FAILED;
				}

				pathItemRead = pathItemRead->nextPathItemRead;
			}

			// free the read nodes in contigPathItem
			pathItemRead = pathItem->pathItemReadList;
			while(pathItemRead)
			{
				pathItemReadTmp = pathItemRead->nextPathItemRead;
				free(pathItemRead);
				pathItemRead = pathItemReadTmp;
			}
			pathItem->pathItemReadList = NULL;
			pathItem->supportReadsNum = 0;

			// update the pathItem list
			if(pathItem->prePathItem)
				pathItem->prePathItem->nextPathItem = pathItem->nextPathItem;
			else
				contigPath->contigPathItemList = pathItem->nextPathItem;
			if(pathItem->nextPathItem)
				pathItem->nextPathItem->prePathItem = pathItem->prePathItem;
			else
				contigPath->tailPathItem = pathItem->prePathItem;

			// free pathItem
			pathItemTmp = pathItem->nextPathItem;
			free(pathItem);
			pathItem = pathItemTmp;

			contigPath->itemNumPathItemList --;
		}else
		{
			pathItem = pathItem->nextPathItem;
		}
	}

	if(delDTReadFlag==YES)
	{
		// delete the failed reads
		if(removeFinishedReadsFromDecisionTable(decisionTable, readsNumDecisionTable)==FAILED)
		{
			printf("line=%d, In %s(), cannot resort decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Replace the contig path read in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short replacePathReadInfoInDT(contigPathItem_t *newPathItem, contigPathItem_t *oldPathItem, contigPathItemRead_t *pathItemRead, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable)
{
	int32_t hashcode, successFlag;
	int64_t rid;
	dtRowIndex_t *pDTRow;
	assemblingreadtype *dtRead;

	successFlag = NO;
	rid = pathItemRead->rid;
	hashcode = rid & RID_LOW_BITS_MASK;
	pDTRow = dtRowHashtable[hashcode];
	while(pDTRow)
	{
		if(pDTRow->rid==rid)
		{
			dtRead = decisionTable + pDTRow->dtRow;
			if(dtRead->contigPathItem==oldPathItem)
			{
				dtRead->contigPathItem = newPathItem;
				successFlag = YES;
				break;
			}else
			{
				printf("line=%d, In %s(), invalid pathItem, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
		pDTRow = pDTRow->next;
	}

	if(successFlag==YES)
		return SUCCESSFUL;
	else
		return FAILED;
}

/**
 * Remove the contig path read information in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removePathReadInfoInDT(contigPathItem_t *pathItem, contigPathItemRead_t *pathItemRead, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable, int32_t delReadFlag)
{
	int32_t hashcode, successFlag;
	int64_t rid;
	dtRowIndex_t *pDTRow;
	assemblingreadtype *dtRead;

	successFlag = NO;
	rid = pathItemRead->rid;
	hashcode = rid & RID_LOW_BITS_MASK;
	pDTRow = dtRowHashtable[hashcode];
	while(pDTRow)
	{
		if(pDTRow->rid==rid)
		{
			dtRead = decisionTable + pDTRow->dtRow;
			if(dtRead->contigPathItem==pathItem && dtRead->contigPathItemRead==pathItemRead)
			{
				dtRead->contigPathItem = NULL;
				dtRead->contigPathItemRead = NULL;
				pathItem->supportReadsNum --;
				successFlag = YES;

				if(delReadFlag==YES)
					dtRead->status = FAILED_STATUS;

				break;
			}else
			{
				printf("line=%d, In %s(), invalid pathItem, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
		pDTRow = pDTRow->next;
	}

	if(successFlag==YES)
		return SUCCESSFUL;
	else
		return FAILED;
}

/**
 * Sort the contig path items.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxesContigPathItems(contigPath_t *contigPath)
{
	int32_t maxNum, secNum;
	contigPathItem_t *pathItem, *maxPathItem, *secPathItem;

	maxNum = secNum = 0;
	maxPathItem = secPathItem = NULL;
	pathItem = contigPath->contigPathItemList;
	while(pathItem)
	{
		if(pathItem->supportReadsNum>0)
		{
			if(maxNum<pathItem->supportReadsNum)
			{
				secNum = maxNum;
				secPathItem = maxPathItem;
				maxNum = pathItem->supportReadsNum;
				maxPathItem = pathItem;
			}else if(secNum<pathItem->supportReadsNum)
			{
				secNum = pathItem->supportReadsNum;
				secPathItem = pathItem;
			}
		}

		pathItem = pathItem->nextPathItem;
	}

	contigPath->maxPathItem = maxPathItem;
	contigPath->secPathItem = secPathItem;

	return SUCCESSFUL;
}

/**
 * Sort the contig path items.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short sortContigPathItems(contigPath_t *contigPath)
{
	int32_t i, itemNum;
	contigPathItem_t *pathItem;
	contigPathItemSort_t *contigPathItemSortArray, *contigPathItemSortBufArray, *contigPathItemSortResultArray;

	if(contigPath->itemNumPathItemList==1)
	{
		contigPath->maxPathItem = contigPath->contigPathItemList;
		contigPath->secPathItem = NULL;
	}else if(contigPath->itemNumPathItemList==2)
	{
		if(contigPath->contigPathItemList->supportReadsNum<contigPath->tailPathItem->supportReadsNum)
		{
			contigPath->contigPathItemList->prePathItem = contigPath->tailPathItem;
			contigPath->tailPathItem->nextPathItem = contigPath->contigPathItemList;
			contigPath->contigPathItemList = contigPath->tailPathItem;
			contigPath->tailPathItem = contigPath->contigPathItemList->nextPathItem;
			contigPath->contigPathItemList->prePathItem = NULL;
			contigPath->tailPathItem->nextPathItem = NULL;
		}

		contigPath->maxPathItem = contigPath->contigPathItemList;
		contigPath->secPathItem = contigPath->contigPathItemList->nextPathItem;

	}else if(contigPath->itemNumPathItemList>=3)
	{ // radix sort
		contigPathItemSortArray = (contigPathItemSort_t *) calloc (contigPath->itemNumPathItemList, sizeof(contigPathItemSort_t));
		if(contigPathItemSortArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		contigPathItemSortBufArray = (contigPathItemSort_t *) calloc (contigPath->itemNumPathItemList, sizeof(contigPathItemSort_t));
		if(contigPathItemSortBufArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// fill the data
		itemNum = 0;
		pathItem = contigPath->contigPathItemList;
		while(pathItem)
		{
			contigPathItemSortArray[itemNum].supportReadsNum = pathItem->supportReadsNum;
			contigPathItemSortArray[itemNum++].pathItem = pathItem;
			pathItem = pathItem->nextPathItem;
		}

		if(contigPath->itemNumPathItemList<=contigPath->bestItemNumPathItemList)
		{
			// comparison sort of contigsLenArray
			if(comparisonSortContigPathItem(contigPathItemSortArray, contigPathItemSortBufArray, itemNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot sort contigs lengths, error!\n", __LINE__, __func__);
				return FAILED;
			}
			contigPathItemSortResultArray = contigPathItemSortBufArray;
		}else
		{
			// radix sort of contigsLenArray
			if(radixSortContigPathItem(contigPathItemSortArray, contigPathItemSortBufArray, itemNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot sort contigs lengths, error!\n", __LINE__, __func__);
				return FAILED;
			}
			contigPathItemSortResultArray = contigPathItemSortArray;
		}

		// update the contig path item list
		pathItem =  contigPathItemSortResultArray[0].pathItem;
		pathItem->prePathItem = NULL;
		pathItem->nextPathItem = contigPathItemSortResultArray[1].pathItem;
		contigPath->contigPathItemList = pathItem;
		for(i=1; i<itemNum-1; i++)
		{
			pathItem =  contigPathItemSortResultArray[i].pathItem;
			pathItem->prePathItem = contigPathItemSortResultArray[i-1].pathItem;
			pathItem->nextPathItem = contigPathItemSortResultArray[i+1].pathItem;
		}
		pathItem =  contigPathItemSortResultArray[itemNum-1].pathItem;
		pathItem->prePathItem = contigPathItemSortResultArray[itemNum-2].pathItem;
		pathItem->nextPathItem = NULL;
		contigPath->tailPathItem = pathItem;

		contigPath->maxPathItem = contigPathItemSortResultArray[0].pathItem;
		contigPath->secPathItem = contigPathItemSortResultArray[1].pathItem;

		free(contigPathItemSortArray);
		free(contigPathItemSortBufArray);
	}else
	{
		contigPath->maxPathItem = NULL;
		contigPath->secPathItem = NULL;
	}

	return SUCCESSFUL;
}

/**
 * Comparison sort the contig path items.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short comparisonSortContigPathItem(contigPathItemSort_t *contigPathItemSortArray, contigPathItemSort_t *contigPathItemSortBufArray, int32_t itemNum)
{
	int32_t i, j, *sortArray;

	sortArray = (int32_t *) calloc (itemNum, sizeof(int32_t));
	if(sortArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<itemNum-1; i++)
	{
		for(j=i+1; j<itemNum; j++)
		{
			if(contigPathItemSortArray[i].supportReadsNum<contigPathItemSortArray[j].supportReadsNum)
				sortArray[i] ++;
			else
				sortArray[j] ++;
		}
	}

	for(i=0; i<itemNum; i++)
	{
		contigPathItemSortBufArray[sortArray[i]].supportReadsNum = contigPathItemSortArray[i].supportReadsNum;
		contigPathItemSortBufArray[sortArray[i]].pathItem = contigPathItemSortArray[i].pathItem;
	}

	// ######################## Debug information ########################
//	for(i=0; i<itemNum-1; i++)
//	{
//		if(contigPathItemSortBufArray[i].supportReadsNum<contigPathItemSortBufArray[i+1].supportReadsNum)
//		{
//			printf("line=%d, In %s(), contigPathItemSortBufArray[i].supportReadsNum=%d, contigPathItemSortBufArray[i+1].supportReadsNum=%d, itemNum=%d, error!\n", __LINE__, __func__, contigPathItemSortBufArray[i].supportReadsNum, contigPathItemSortBufArray[i+1].supportReadsNum, itemNum);
//			return FAILED;
//		}
//	}
	// ######################## Debug information ########################

	free(sortArray);

	return SUCCESSFUL;
}

/**
 * Radix sort the contig path items.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short radixSortContigPathItem(contigPathItemSort_t *contigPathItemSortArray, contigPathItemSort_t *contigPathItemSortBufArray, int32_t itemNum)
{
	struct partNode
	{
		int curItemNum;
		int totalItemNum;
		int firstRow;
	};

	int64_t i, step, total;
	contigPathItemSort_t *data, *buf;
	struct partNode *part;
	int partArrSize, stepBits, maxStepLen;
	uint64_t bitMask;
	uint64_t hashcode, firstRow, curItemNum;

	stepBits = 16;
	maxStepLen = 32;
	partArrSize = 1 << stepBits;
	bitMask = (1 << stepBits) - 1;

	part = (struct partNode *) malloc(partArrSize * sizeof(struct partNode));
	if(part==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Begin to sort
	step = 0;
	while(step!=maxStepLen)
	{
		// set the data and buf
		if(step==stepBits)
		{
			buf = contigPathItemSortArray;
			data = contigPathItemSortBufArray;
		}else
		{
			data = contigPathItemSortArray;
			buf = contigPathItemSortBufArray;
		}

		// count the number
		if(memset(part, 0, partArrSize * sizeof(struct partNode))==NULL)
		{
			printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		for(i=0; i<itemNum; i++)
			part[ bitMask - ((data[i].supportReadsNum >> step) & bitMask) ].totalItemNum ++;

		// initialize the part array
		total = 0;
		for(i=0; i<partArrSize; i++)
		{
			part[i].firstRow = total;
			total += part[i].totalItemNum;
		}

		// copy the data to the right place
		for(i=0; i<itemNum; i++)
		{
			hashcode = bitMask - ((data[i].supportReadsNum >> step) & bitMask);
			firstRow = part[hashcode].firstRow;
			curItemNum = part[hashcode].curItemNum;
			buf[firstRow+curItemNum].supportReadsNum = data[i].supportReadsNum;
			buf[firstRow+curItemNum].pathItem = data[i].pathItem;
//			if(memcpy(buf+firstRow+curItemNum, data+i, sizeof(int32_t))==NULL)
//			{
//				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
//				free(part);
//				return FAILED;
//			}
			part[hashcode].curItemNum ++;
		}

		step += stepBits;

		//######################## Debug information #######################
//		for(i=0; i<partArrSize; i++)
//		{
//			if(part[i].curItemNum!=part[i].totalItemNum)
//			{
//				printf("line=%d, In %s(), in part[%ld], curItemNum=%d != totalItemNum=%d, error!\n", __LINE__, __func__, i, part[i].curItemNum, part[i].totalItemNum);
//				free(part);
//				return FAILED;
//			}
//		}
		//######################## Debug information #######################
	}

	// ######################## Debug information ########################
//	for(i=0; i<itemNum-1; i++)
//	{
//		if(contigPathItemSortArray[i].supportReadsNum<contigPathItemSortArray[i+1].supportReadsNum)
//		{
//			printf("line=%d, In %s(), contigPathItemSortArray[i].supportReadsNum=%d, contigPathItemSortArray[i+1].supportReadsNum=%d, itemNum=%d, error!\n", __LINE__, __func__, contigPathItemSortArray[i].supportReadsNum, contigPathItemSortArray[i+1].supportReadsNum, itemNum);
//			return FAILED;
//		}
//	}
	// ######################## Debug information ########################

	free(part);
	part = NULL;

	return SUCCESSFUL;
}

/**
 * Remove the less supported contig path items.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeLessSupportedContigPathItems(contigPath_t *contigPath, int32_t useOldNaviPathseqFlag, char *oldNaviPathseq, int32_t oldNaviPathseqLen, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable)
{
	int32_t i, maxLen, maxReadsNum, pathLen, matchFlag;
	char *pathseq;
	contigPathItem_t *pathItem, *pathItemTmp;
	contigPathItemRead_t *pathItemRead, *pathItemReadTmp;

	if(contigPath->itemNumPathItemList>=2)
	{
		// get the maximum reads number
		maxReadsNum = 0;
		pathItem = contigPath->contigPathItemList;
		while(pathItem)
		{
			if(maxReadsNum<pathItem->supportReadsNum)
				maxReadsNum = pathItem->supportReadsNum;
			pathItem = pathItem->nextPathItem;
		}

		if(maxReadsNum>=2)
		{
			pathItem = contigPath->contigPathItemList;
			while(pathItem)
			{
				if(pathItem->supportReadsNum<=1)
				{
					if(useOldNaviPathseqFlag==YES && oldNaviPathseqLen>10)
					{
						matchFlag = YES;
						pathseq = pathItem->contigPathStr + contigPath->startRowNewBase;
						pathLen = pathItem->contigPathLen - contigPath->startRowNewBase;
						if(pathLen>=oldNaviPathseqLen)
						{
							for(i=0; i<oldNaviPathseqLen; i++)
							{
								if(pathseq[i]!=oldNaviPathseq[i])
								{
									matchFlag = NO;
									break;
								}
							}
						}
					}else
						matchFlag = NO;

					if(matchFlag==NO)
					{
						// delete the path information of read in decision table
						pathItemRead = pathItem->pathItemReadList;
						while(pathItemRead)
						{
							if(removePathReadInfoInDT(pathItem, pathItemRead, decisionTable, dtRowHashtable, NO)==FAILED)
							{
								printf("line=%d, In %s(), cannot replace the path read information in decision table, error!\n", __LINE__, __func__);
								return FAILED;
							}

							pathItemRead = pathItemRead->nextPathItemRead;
						}

						// free the read nodes in contigPathItem
						pathItemRead = pathItem->pathItemReadList;
						while(pathItemRead)
						{
							pathItemReadTmp = pathItemRead->nextPathItemRead;
							free(pathItemRead);
							pathItemRead = pathItemReadTmp;
						}
						pathItem->pathItemReadList = NULL;
						pathItem->supportReadsNum = 0;

						// update the pathItem list
						if(pathItem->prePathItem)
							pathItem->prePathItem->nextPathItem = pathItem->nextPathItem;
						else
							contigPath->contigPathItemList = pathItem->nextPathItem;
						if(pathItem->nextPathItem)
							pathItem->nextPathItem->prePathItem = pathItem->prePathItem;
						else
							contigPath->tailPathItem = pathItem->prePathItem;

						// free pathItem
						pathItemTmp = pathItem->nextPathItem;
						free(pathItem);
						pathItem = pathItemTmp;

						contigPath->itemNumPathItemList --;
					}else
						pathItem = pathItem->nextPathItem;
				}else
					pathItem = pathItem->nextPathItem;
			}

			maxLen = 0;
			pathItem = contigPath->contigPathItemList;
			while(pathItem)
			{
				if(maxLen<pathItem->contigPathLen)
					maxLen = pathItem->contigPathLen;

				pathItem = pathItem->nextPathItem;
			}
			contigPath->maxContigPathLen = maxLen;
		}
	}

	return SUCCESSFUL;
}

/**
 * Remove contig path items with short fragment size.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeShortFragSizeContigPath(contigPath_t *contigPath, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable)
{
	double averFragSize, averFragSizeDif;
	int32_t shortFragSizeFlag, deledFlag;
	contigPathItem_t *pathItem, *pathItemTmp;
	contigPathItemRead_t *pathItemRead, *pathItemReadTmp;

	if(contigPath->itemNumPathItemList>=2)
	{
		deledFlag = NO;
		pathItem = contigPath->contigPathItemList;
		while(pathItem)
		{
			// get the average fragment size
			if(computeAverFragSizeContigPathItem(&averFragSize, pathItem, decisionTable, dtRowHashtable)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the overlap type of two tandem paths, error!\n", __LINE__, __func__);
				return FAILED;
			}

			averFragSizeDif = meanSizeInsert - averFragSize;
			if(averFragSizeDif<0)
				averFragSizeDif = -averFragSizeDif;
			//if(averFragSizeDif>2*standardDev || averFragSizeDif>0.2*meanSizeInsert)
			if(averFragSizeDif>2*standardDev || averFragSizeDif>0.1*meanSizeInsert)
				shortFragSizeFlag = YES;
			else
				shortFragSizeFlag = NO;

			if(shortFragSizeFlag==YES)
			{ // remove condition

				deledFlag = YES;
				pathItemRead = pathItem->pathItemReadList;
				while(pathItemRead)
				{
					if(removePathReadInfoInDT(pathItem, pathItemRead, decisionTable, dtRowHashtable, YES)==FAILED)
					{
						printf("line=%d, In %s(), cannot replace the path read information in decision table, error!\n", __LINE__, __func__);
						return FAILED;
					}

					pathItemRead = pathItemRead->nextPathItemRead;
				}

				// free the read nodes in contigPathItem
				pathItemRead = pathItem->pathItemReadList;
				while(pathItemRead)
				{
					pathItemReadTmp = pathItemRead->nextPathItemRead;
					free(pathItemRead);
					pathItemRead = pathItemReadTmp;
				}
				pathItem->pathItemReadList = NULL;
				pathItem->supportReadsNum = 0;

				// update the pathItem list
				if(pathItem->prePathItem)
					pathItem->prePathItem->nextPathItem = pathItem->nextPathItem;
				else
					contigPath->contigPathItemList = pathItem->nextPathItem;
				if(pathItem->nextPathItem)
					pathItem->nextPathItem->prePathItem = pathItem->prePathItem;
				else
					contigPath->tailPathItem = pathItem->prePathItem;

				// free pathItem
				pathItemTmp = pathItem->nextPathItem;
				free(pathItem);
				pathItem = pathItemTmp;

				contigPath->itemNumPathItemList --;
			}else
				pathItem = pathItem->nextPathItem;

			pathItem = pathItem->nextPathItem;
		}

		if(deledFlag==YES)
		{
			// delete the failed reads
			if(removeFinishedReadsFromDecisionTable(decisionTable, readsNumDecisionTable)==FAILED)
			{
				printf("line=%d, In %s(), cannot resort decision table, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Remove reads in second path of two overlapped tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeAverFragSizeContigPathItem(double *averFragSize, contigPathItem_t *pathItem, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable)
{
	uint64_t readID, readID_paired, pairedNum;
	int32_t validReadOrient, validReadOrient_paired;
	PERead_t *pPERead;
	assemblingreadtype *dtRead, *dtReadPaired;
	contigPathItemRead_t *pathItemRead, *pathItemReadTmp;

	*averFragSize = 0;
	pairedNum = 0;
	pathItemRead = pathItem->pathItemReadList;
	while(pathItemRead)
	{
		if(getExistReadInDT(&dtRead, pathItemRead->rid, decisionTable, dtRowHashtable)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the read in decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(dtRead==NULL)
		{
			printf("line=%d, In %s(), cannot get the read in decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}

		readID = pathItemRead->rid;

		// generate the paired readID
		if(readID%2==1)
		{ // odd --> even
			readID_paired = readID + 1;
		}else
		{ // even --> odd
			readID_paired = readID - 1;
		}

		if(getReadFromPEHashtable(&pPERead, readID_paired)==FAILED)
		{
			printf("In %s(), cannot get the paired read %lu, error!\n", __func__, readID_paired);
			return FAILED;
		}

		if(pPERead)
		{
			*averFragSize += dtRead->firstContigPos + dtRead->seqlen - pPERead->cpos;
			pairedNum ++;
		}
		else if(shortInsertFlag==YES)
		{
			// get the exist read
			if(getExistReadInDT(&dtReadPaired, readID_paired, decisionTable, dtRowHashtable)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the exist read %lu in decision table, error!\n", __LINE__, __func__, readID_paired);
				return FAILED;
			}

			//if((*dtReadPaired) && (*dtReadPaired)->orientation==ORIENTATION_PLUS && (*dtReadPaired)->unmatchBaseNum==0)
			if(dtReadPaired && dtReadPaired->orientation==ORIENTATION_PLUS && dtReadPaired->unmatchBaseNum==0)
			{
				*averFragSize += dtRead->firstContigPos + dtRead->seqlen - dtReadPaired->firstContigPos;
				pairedNum ++;
			}
		}

		pathItemRead = pathItemRead->nextPathItemRead;
	}

	if(pairedNum>0)
		*averFragSize /= pairedNum;

	return SUCCESSFUL;
}

/**
 * Remove the short overlapped contig path items.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeShortOverlappedContigPathItems(contigPath_t *contigPath, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable)
{
	int32_t i, j, maxLen, shareLen, mismatchNum, firstDelFlag;
	contigPathItem_t *pathItem, *pathItem2, *pathItemDel, *pathItemTmp;
	contigPathItemRead_t *pathItemRead, *pathItemReadTmp;
	char *pseq1, *pseq2;
	int32_t maxOverlapSizeWithContig;
	double overlapFactor;

	if(contigPath->itemNumPathItemList>=2)
	{
		pathItem = contigPath->contigPathItemList;
		while(pathItem)
		{
			firstDelFlag = NO;
			pathItem2 = pathItem->nextPathItem;
			while(pathItem2)
			{
				// compute the mismatch base number
				if(pathItem->contigPathLen-contigPath->startRowNewBase<pathItem2->contigPathLen-contigPath->startRowNewBase)
					shareLen = pathItem->contigPathLen - contigPath->startRowNewBase;
				else
					shareLen = pathItem2->contigPathLen - contigPath->startRowNewBase;

				pseq1 = pathItem->contigPathStr + contigPath->startRowNewBase;
				pseq2 = pathItem2->contigPathStr + contigPath->startRowNewBase;
				mismatchNum = 0;
				for(i=0; i<shareLen; i++)
					if(pseq1[i]!=pseq2[i])
						mismatchNum ++;

				firstDelFlag = NO;
				if(mismatchNum>contigPath->maxMismatchNumThres)
				{
					pathItemDel = NULL;
//					if(pathItem->maxOverlapWithContig>0.9*readLen && pathItem2->maxOverlapWithContig<0.7*readLen)
//					{
//						pathItemDel = pathItem2;
//					}else
					if(pathItem2->maxOverlapWithContig<contigPath->overlapWithContigThres)
						pathItemDel = pathItem2;
					else  if(pathItem->maxOverlapWithContig<contigPath->overlapWithContigThres)
					{
						firstDelFlag = YES;
						pathItemDel = pathItem;
					}

					if(pathItemDel)
					{
						// delete the path information of read in decision table
						pathItemRead = pathItemDel->pathItemReadList;
						while(pathItemRead)
						{
							if(removePathReadInfoInDT(pathItemDel, pathItemRead, decisionTable, dtRowHashtable, NO)==FAILED)
							{
								printf("line=%d, In %s(), cannot replace the path read information in decision table, error!\n", __LINE__, __func__);
								return FAILED;
							}

							pathItemRead = pathItemRead->nextPathItemRead;
						}

						// free the read nodes in contigPathItem
						pathItemRead = pathItemDel->pathItemReadList;
						while(pathItemRead)
						{
							pathItemReadTmp = pathItemRead->nextPathItemRead;
							free(pathItemRead);
							pathItemRead = pathItemReadTmp;
						}
						pathItemDel->pathItemReadList = NULL;
						pathItemDel->supportReadsNum = 0;

						// update the pathItem list
						if(pathItemDel->prePathItem)
							pathItemDel->prePathItem->nextPathItem = pathItemDel->nextPathItem;
						else
							contigPath->contigPathItemList = pathItemDel->nextPathItem;
						if(pathItemDel->nextPathItem)
							pathItemDel->nextPathItem->prePathItem = pathItemDel->prePathItem;
						else
							contigPath->tailPathItem = pathItemDel->prePathItem;
						contigPath->itemNumPathItemList --;

						// free pathItem
						pathItemTmp = pathItemDel->nextPathItem;
						free(pathItemDel);

						if(firstDelFlag==YES)
						{
							pathItem = pathItemTmp;
							break;
						}else
						{
							pathItem2 = pathItemTmp;
						}
					}else
					{
						pathItem2 = pathItem2->nextPathItem;
					}
				}else
				{
					pathItem2 = pathItem2->nextPathItem;
				}
			}

			if(firstDelFlag==NO)
				pathItem = pathItem->nextPathItem;

			if(contigPath->itemNumPathItemList<=1)
				break;
		}

		maxLen = 0;
		pathItem = contigPath->contigPathItemList;
		while(pathItem)
		{
			if(maxLen<pathItem->contigPathLen)
				maxLen = pathItem->contigPathLen;

			pathItem = pathItem->nextPathItem;
		}
		contigPath->maxContigPathLen = maxLen;
	}

/*
	// remove short contig path item with overlap size
	if(contigPath->itemNumPathItemList>=2)
	{
		// get the maxmial overlap size
		maxOverlapSizeWithContig = 0;
		pathItem = contigPath->contigPathItemList;
		while(pathItem)
		{
			if(maxOverlapSizeWithContig<pathItem->maxOverlapWithContig)
				maxOverlapSizeWithContig = pathItem->maxOverlapWithContig;

			pathItem = pathItem->nextPathItem;
		}

		if(maxOverlapSizeWithContig>0.6*readLen)
		{
			//overlapFactor = (double)maxOverlapSizeWithContig / readLen;
			pathItem = contigPath->contigPathItemList;
			while(pathItem)
			{
				//if(pathItem->maxOverlapWithContig<(overlapFactor-0.2)*readLen)
				if(pathItem->maxOverlapWithContig<0.3*readLen)
				{
					// delete the path information of read in decision table
					pathItemRead = pathItem->pathItemReadList;
					while(pathItemRead)
					{
						if(removePathReadInfoInDT(pathItem, pathItemRead, decisionTable, dtRowHashtable, NO)==FAILED)
						{
							printf("line=%d, In %s(), cannot replace the path read information in decision table, error!\n", __LINE__, __func__);
							return FAILED;
						}

						pathItemRead = pathItemRead->nextPathItemRead;
					}

					// free the read nodes in contigPathItem
					pathItemRead = pathItem->pathItemReadList;
					while(pathItemRead)
					{
						pathItemReadTmp = pathItemRead->nextPathItemRead;
						free(pathItemRead);
						pathItemRead = pathItemReadTmp;
					}
					pathItem->pathItemReadList = NULL;
					pathItem->supportReadsNum = 0;

					// update the pathItem list
					if(pathItem->prePathItem)
						pathItem->prePathItem->nextPathItem = pathItem->nextPathItem;
					else
						contigPath->contigPathItemList = pathItem->nextPathItem;
					if(pathItem->nextPathItem)
						pathItem->nextPathItem->prePathItem = pathItem->prePathItem;
					else
						contigPath->tailPathItem = pathItem->prePathItem;
					contigPath->itemNumPathItemList --;

					// free pathItem
					pathItemTmp = pathItem->nextPathItem;
					free(pathItem);

					pathItem = pathItemTmp;
				}else
				{
					pathItem = pathItem->nextPathItem;
				}
			}
		}

		maxLen = 0;
		pathItem = contigPath->contigPathItemList;
		while(pathItem)
		{
			if(maxLen<pathItem->contigPathLen)
				maxLen = pathItem->contigPathLen;

			pathItem = pathItem->nextPathItem;
		}
		contigPath->maxContigPathLen = maxLen;
	}
*/

	return SUCCESSFUL;
}

/**
 * Remain the best path item.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short remainBestContigPathItems(contigPath_t *contigPath, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable)
{
	int32_t itemNum;
	contigPathItem_t *pathItem, *pathItemTmp;
	contigPathItemRead_t *pathItemRead, *pathItemReadTmp;

	if(contigPath->itemNumPathItemList>contigPath->bestItemNumPathItemList)
	{
		itemNum = 1;
		pathItem = contigPath->contigPathItemList;
		while(pathItem)
		{
			if(itemNum>contigPath->bestItemNumPathItemList)
			{
				pathItem->prePathItem->nextPathItem = NULL;
				contigPath->tailPathItem = pathItem->prePathItem;
				break;
			}
			itemNum ++;
			pathItem = pathItem->nextPathItem;
		}

		while(pathItem)
		{
			// update the reads in decision table
			pathItemRead = pathItem->pathItemReadList;
			while(pathItemRead)
			{
				if(removePathReadInfoInDT(pathItem, pathItemRead, decisionTable, dtRowHashtable, NO)==FAILED)
				{
					printf("line=%d, In %s(), cannot replace the path read information in decision table, error!\n", __LINE__, __func__);
					return FAILED;
				}

				pathItemRead = pathItemRead->nextPathItemRead;
			}

			// delete the reads in path item
			pathItemRead = pathItem->pathItemReadList;
			while(pathItemRead)
			{
				pathItemReadTmp = pathItemRead->nextPathItemRead;
				free(pathItemRead);
				pathItemRead = pathItemReadTmp;
			}
			pathItem->pathItemReadList = NULL;
			pathItem->supportReadsNum = 0;
			contigPath->itemNumPathItemList --;

			// delete the path item node
			pathItemTmp = pathItem->nextPathItem;
			free(pathItem);
			pathItem = pathItemTmp;
		}
	}

	return SUCCESSFUL;
}

/**
 * Shift the start row of new base in contig path to zero.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short shiftStartRowNewBaseContigPath(contigPath_t *contigPath)
{
	int32_t i, pathLen, startRow;
	char *pathseq;
	contigPathItem_t *pathItem;

	startRow = contigPath->startRowNewBase;
	pathItem = contigPath->contigPathItemList;
	while(pathItem)
	{
		pathseq = pathItem->contigPathStr;
		pathLen = pathItem->contigPathLen - startRow;
		for(i=0; i<pathLen; i++)
			pathseq[i] = pathseq[startRow+i];
		pathseq[pathLen] = '\0';
		pathItem->contigPathLen = pathLen;

		pathItem = pathItem->nextPathItem;
	}
	contigPath->startRowNewBase = 0;


	return SUCCESSFUL;
}

/**
 * Output the contig paths.
 */
void outputContigPath(contigPath_t *contigPath, int32_t allInfoFlag)
{
	int32_t startRow, readNumPerLine, readNum, outNum;
	contigPathItem_t *pathItem;
	contigPathItemRead_t *pathItemRead;

	if(contigPath->itemNumPathItemList==0)
	{
		printf("There are no contig paths.\n");
		return;
	}

	printf("---- updateInterval=%d, startRowNewBase=%d, maxContigPathLen=%d, itemNumPathItemList=%d, naviSuccessSize=%d, preNaviSuccessSize=%d, preNaviOverSize=%d\n", contigPath->updateInterval, contigPath->startRowNewBase, contigPath->maxContigPathLen-contigPath->startRowNewBase, contigPath->itemNumPathItemList, contigPath->naviSuccessSize, contigPath->preNaviSuccessSize, contigPath->preNaviOverlapSize);
	if(contigPath->validCandPathTandPathFlagPE==YES)
		printf("---- CTPathseq= %s, len=%d, %d\n", contigPath->candPathseqTandPathPE+contigPath->startRowCandPathTandPathPE, contigPath->candPathLenTandPathPE-contigPath->startRowCandPathTandPathPE, contigPath->appearTimesCandTandPathPE);

	outNum = 0;
	readNumPerLine = 10;
	pathItem = contigPath->contigPathItemList;
	while(pathItem)
	{
		outNum ++;
		if(outNum<=10 && pathItem->contigPathLen-contigPath->startRowNewBase>0 && (pathItem->supportReadsNum>0 || allInfoFlag==YES))
		{
			printf("\t%d\t%s\t%d\t%d\n", pathItem->supportReadsNum, pathItem->contigPathStr+contigPath->startRowNewBase, pathItem->maxOverlapWithContig, pathItem==contigPath->naviPathItem);
/*
			pathItemRead = pathItem->pathItemReadList;
			readNum = 1;
			while(pathItemRead)
			{
				if(readNum==readNumPerLine)
				{
					printf("\t%ld\n", pathItemRead->rid);
					readNum = 0;
				}else
				{
					printf("\t%ld", pathItemRead->rid);
				}

				readNum ++;
				pathItemRead = pathItemRead->nextPathItemRead;
			}
			printf("\n");
*/
		}

		pathItem = pathItem->nextPathItem;
	}

	// output the tail sequence
	startRow = contigPath->contigtailSeqLen-readLen;
	if(startRow<0)
		startRow = 0;
	printf("\tcontigtailSeqLen=%d, contigtail sequence=%s\n", contigPath->contigtailSeqLen, contigPath->contigtailSeq+startRow);
}

/**
 * Initialize the tail sequence of contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initTailseqContigPath(contigPath_t *contigPath, contigtype *contigArray, int32_t itemNumContigArray)
{
	int32_t i, baseNum, startRow;
	char baseTmp;

	if(itemNumContigArray>contigPath->defaultInitContigtailSeqLen)
		baseNum = contigPath->defaultInitContigtailSeqLen;
	else
		baseNum = itemNumContigArray;

	startRow = itemNumContigArray - baseNum;
	for(i=0; i<baseNum; i++)
	{
		switch(contigArray[startRow+i].base)
		{
			case 0: baseTmp = 'A'; break;
			case 1: baseTmp = 'C'; break;
			case 2: baseTmp = 'G'; break;
			case 3: baseTmp = 'T'; break;
			default: printf("line=%d, In %s(), invalid baseInt=%d, error!\n", __LINE__, __func__, contigArray[startRow+i].base); return FAILED;
		}

		contigPath->contigtailSeq[i] = baseTmp;
	}
	contigPath->contigtailSeq[baseNum] = '\0';
	contigPath->contigtailSeqLen = baseNum;

	return SUCCESSFUL;
}

/**
 * Append the tail base to the contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short appendTailBaseToContigPath(contigPath_t *contigPath, contigtype *contigArray, int32_t itemNumContigArray)
{
	int32_t i, startRowOld, startRowNew, baseNum;
	char baseTmp;

	switch(contigArray[itemNumContigArray-1].base)
	{
		case 0: baseTmp = 'A'; break;
		case 1: baseTmp = 'C'; break;
		case 2: baseTmp = 'G'; break;
		case 3: baseTmp = 'T'; break;
		default: printf("line=%d, In %s(), invalid baseInt=%d, error!\n", __LINE__, __func__, contigArray[itemNumContigArray-1].base); return FAILED;
	}

	contigPath->contigtailSeq[contigPath->contigtailSeqLen++] = baseTmp;
	contigPath->contigtailSeq[contigPath->contigtailSeqLen] = '\0';

	if(contigPath->contigtailSeqLen>=contigPath->maxContigtailSeqLen)
	{
		if(itemNumContigArray>contigPath->defaultInitContigtailSeqLen)
			baseNum = contigPath->defaultInitContigtailSeqLen;
		else
			baseNum = itemNumContigArray;

		startRowNew = 0;
		startRowOld = contigPath->contigtailSeqLen - baseNum;
		for(i=0; i<baseNum; i++)
			contigPath->contigtailSeq[startRowNew+i] = contigPath->contigtailSeq[startRowOld+i];
		contigPath->contigtailSeq[baseNum] = '\0';
		contigPath->contigtailSeqLen = baseNum;
	}

	return SUCCESSFUL;
}

/**
 * Remove reads of identical part of path items.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeReadsInIdenticalPathItemPart(contigPath_t *contigPath, int32_t useOldNaviPathseqFlag, char *oldNaviPathseq, int32_t oldNaviPathseqLen, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, int32_t itemNumContigArray)
{
	int32_t i, j, mismatchNum, identicalSize, pathLen1, pathLen2, shareLen, useMateFlag;
	int32_t uniqueReadNum1, uniqueReadNum2, remainBaseNum, newStartPos, readseqLenTmp;
	char *pathseq1, *pathseq2, readseqTmp[MAX_READ_LEN_IN_BUF+1];
	contigPathItemRead_t *pathItemRead, *pathItemReadTmp;
	assemblingreadtype *dtRead;


	if(useOldNaviPathseqFlag==YES && contigPath->itemNumPathItemList>=2)
	{
		pathseq1 = contigPath->maxPathItem->contigPathStr + contigPath->startRowNewBase;
		pathseq2 = contigPath->secPathItem->contigPathStr + contigPath->startRowNewBase;
		pathLen1 = contigPath->maxPathItem->contigPathLen - contigPath->startRowNewBase;
		pathLen2 = contigPath->secPathItem->contigPathLen - contigPath->startRowNewBase;
		if(pathLen1<pathLen2)
			shareLen = pathLen1;
		else
			shareLen = pathLen2;

		identicalSize = 0;
		for(i=0; i<shareLen; i++)
		{
			if(pathseq1[i]==pathseq2[i])
				identicalSize ++;
			else
				break;
		}

		if(oldNaviPathseqLen<=identicalSize)  // contained in the old navi path
		{
			if(meanSizeInsert>0 && itemNumContigArray>2*meanSizeInsert)
				useMateFlag = YES;
			else
				useMateFlag = NO;

			uniqueReadNum1 = 0;
			pathItemRead = contigPath->maxPathItem->pathItemReadList;
			while(pathItemRead)
			{
				// get the exist read in decision table
				if(getExistReadInDT(&dtRead, pathItemRead->rid, decisionTable, dtRowHashtable)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the exist read %lu in decision table, error!\n", __LINE__, __func__, pathItemRead->rid);
					return FAILED;
				}

				if(dtRead)
				{
					if(useMateFlag==YES)
					{
						if(dtRead->matedFlag==YES)
						{
							if(dtRead->orientation==ORIENTATION_PLUS)
								remainBaseNum = dtRead->seqlen - dtRead->basePos - 1;
							else
								remainBaseNum = dtRead->basePos;
							if(remainBaseNum>identicalSize)
								uniqueReadNum1 ++;
						}
					}else
					{
						if(dtRead->orientation==ORIENTATION_PLUS)
							remainBaseNum = dtRead->seqlen - dtRead->basePos - 1;
						else
							remainBaseNum = dtRead->basePos;
						if(remainBaseNum>identicalSize)
							uniqueReadNum1 ++;
					}
				}else
				{
					printf("line=%d, In %s(), cannot get the exist read %lu in decision table, error!\n", __LINE__, __func__, pathItemRead->rid);
					return FAILED;
				}

				pathItemRead = pathItemRead->nextPathItemRead;
			}

			// secPathItem
			uniqueReadNum2 = 0;
			pathItemRead = contigPath->secPathItem->pathItemReadList;
			while(pathItemRead)
			{
				// get the exist read in decision table
				if(getExistReadInDT(&dtRead, pathItemRead->rid, decisionTable, dtRowHashtable)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the exist read %lu in decision table, error!\n", __LINE__, __func__, pathItemRead->rid);
					return FAILED;
				}

				if(dtRead)
				{
					if(useMateFlag==YES)
					{
						if(dtRead->matedFlag==YES)
						{
							if(dtRead->orientation==ORIENTATION_PLUS)
								remainBaseNum = dtRead->seqlen - dtRead->basePos - 1;
							else
								remainBaseNum = dtRead->basePos;
							if(remainBaseNum>identicalSize)
								uniqueReadNum2 ++;
						}
					}else
					{
						if(dtRead->orientation==ORIENTATION_PLUS)
							remainBaseNum = dtRead->seqlen - dtRead->basePos - 1;
						else
							remainBaseNum = dtRead->basePos;
						if(remainBaseNum>identicalSize)
							uniqueReadNum2 ++;
					}
				}else
				{
					printf("line=%d, In %s(), cannot get the exist read %lu in decision table, error!\n", __LINE__, __func__, pathItemRead->rid);
					return FAILED;
				}

				pathItemRead = pathItemRead->nextPathItemRead;
			}

			//printf("uniqueReadNum1=%d, uniqueReadNum2=%d\n", uniqueReadNum1, uniqueReadNum2);

			// update the reads in path item
			if(3*uniqueReadNum1<uniqueReadNum2)
			{
				pathItemRead = contigPath->maxPathItem->pathItemReadList;
				while(pathItemRead)
				{
					pathItemReadTmp = pathItemRead->nextPathItemRead;

					// get the exist read in decision table
					if(getExistReadInDT(&dtRead, pathItemRead->rid, decisionTable, dtRowHashtable)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the exist read %lu in decision table, error!\n", __LINE__, __func__, pathItemRead->rid);
						return FAILED;
					}

					if(dtRead)
					{
						if(dtRead->orientation==ORIENTATION_PLUS)
							remainBaseNum = dtRead->seqlen - dtRead->basePos - 1;
						else
							remainBaseNum = dtRead->basePos;
						if(remainBaseNum<=identicalSize)
						{ // move to secPathItem
							if(pathItemRead->prePathItemRead)
								pathItemRead->prePathItemRead->nextPathItemRead = pathItemRead->nextPathItemRead;
							else
								contigPath->maxPathItem->pathItemReadList = pathItemRead->nextPathItemRead;

							if(pathItemRead->nextPathItemRead)
								pathItemRead->nextPathItemRead->prePathItemRead = pathItemRead->prePathItemRead;
							else
								contigPath->maxPathItem->tailPathItemRead = pathItemRead->prePathItemRead;
							contigPath->maxPathItem->supportReadsNum --;

							pathItemRead->prePathItemRead = contigPath->secPathItem->tailPathItemRead;
							pathItemRead->nextPathItemRead = NULL;
							contigPath->secPathItem->tailPathItemRead->nextPathItemRead = pathItemRead;
							contigPath->secPathItem->tailPathItemRead = pathItemRead;
							contigPath->secPathItem->supportReadsNum ++;

							if(dtRead->contigPathItem==contigPath->maxPathItem && dtRead->contigPathItemRead==pathItemRead)
								dtRead->contigPathItem = contigPath->secPathItem;
							else
							{
								printf("line=%d, In %s(), cannot update the contigPathItem of read %lu in decision table, error!\n", __LINE__, __func__, pathItemRead->rid);
								return FAILED;
							}
						}
					}else
					{
						printf("line=%d, In %s(), cannot get the exist read %lu in decision table, error!\n", __LINE__, __func__, pathItemRead->rid);
						return FAILED;
					}

					pathItemRead = pathItemReadTmp;
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Set the naviPathItem of contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short setNaviContigPathItem(contigPath_t *contigPath, int32_t useOldNaviPathseqFlag, char *oldNaviPathseq, int32_t oldNaviPathseqLen)
{
	int32_t i, matchFlag;
	char *pathseq;
	contigPathItem_t *contigPathItem;

	contigPath->naviPathItem = NULL;

	if(useOldNaviPathseqFlag==YES)
	{
		contigPathItem = contigPath->contigPathItemList;
		while(contigPathItem)
		{
			if(contigPathItem->contigPathLen-contigPath->startRowNewBase>=oldNaviPathseqLen)
			{
				matchFlag = YES;
				pathseq = contigPathItem->contigPathStr + contigPath->startRowNewBase;
				for(i=0; i<oldNaviPathseqLen; i++)
				{
					if(pathseq[i]!=oldNaviPathseq[i])
					{
						matchFlag = NO;
						break;
					}
				}

				if(matchFlag==YES)
				{
					contigPath->naviPathItem = contigPathItem;
					break;
				}
			}

			contigPathItem = contigPathItem->nextPathItem;
		}

		//if(contigPath->naviPathItem==NULL || contigPath->naviPathItem->supportReadsNum<0.6*contigPath->maxPathItem->supportReadsNum)
		if(contigPath->naviPathItem==NULL)
		{
			contigPath->preNaviOverlapSize = 0;
			contigPath->preNaviSuccessSize = contigPath->naviSuccessSize;
			contigPath->naviPathItem = contigPath->maxPathItem;
			contigPath->naviSuccessSize = 0;
		}else
		{
			contigPath->preNaviOverlapSize = oldNaviPathseqLen;
		}
	}else
	{
		contigPath->preNaviOverlapSize = 0;
		contigPath->preNaviSuccessSize = contigPath->naviSuccessSize;
		contigPath->naviPathItem = contigPath->maxPathItem;
		contigPath->naviSuccessSize = 0;
	}

	return SUCCESSFUL;
}

/**
 * Navigation by contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateNaviContigPathItem(contigPath_t *contigPath, int32_t contigBaseInt)
{
	int32_t baseIndex, naviPathseqLen;
	char *naviPathseq;
	contigPathItem_t *naviPathItem;

	if(contigPath->naviPathItem)
	{
		naviPathItem = contigPath->naviPathItem;
		naviPathseq = naviPathItem->contigPathStr + contigPath->startRowNewBase;
		naviPathseqLen = naviPathItem->contigPathLen - contigPath->startRowNewBase;

		if(naviPathseqLen>0)
		{
			switch(naviPathseq[0])
			{
				case 'A': baseIndex = 0; break;
				case 'C': baseIndex = 1; break;
				case 'G': baseIndex = 2; break;
				case 'T': baseIndex = 3; break;
				default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, naviPathseq[0]); return FAILED;
			}


			if(baseIndex==contigBaseInt)
			{
				contigPath->naviSuccessSize ++;
			}else
			{
				contigPath->preNaviSuccessSize = contigPath->naviSuccessSize;
				contigPath->naviPathItem = NULL;
				contigPath->naviSuccessSize = 0;
			}
		}else
		{
			contigPath->preNaviSuccessSize = contigPath->naviSuccessSize;
			contigPath->naviPathItem = NULL;
			contigPath->naviSuccessSize = 0;
		}
	}


	if(contigPath->validCandPathTandPathFlagPE==YES)
	{
		naviPathseq = contigPath->candPathseqTandPathPE + contigPath->startRowCandPathTandPathPE;
		naviPathseqLen = contigPath->candPathLenTandPathPE - contigPath->startRowCandPathTandPathPE;
		if(naviPathseqLen>0)
		{
			switch(naviPathseq[0])
			{
				case 'A': baseIndex = 0; break;
				case 'C': baseIndex = 1; break;
				case 'G': baseIndex = 2; break;
				case 'T': baseIndex = 3; break;
				default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, naviPathseq[0]); return FAILED;
			}

			if(baseIndex==contigBaseInt)
			{
				contigPath->startRowCandPathTandPathPE ++;
			}else
			{
				contigPath->validCandPathTandPathFlagPE = NO;
				contigPath->startRowCandPathTandPathPE = -1;
				contigPath->appearTimesCandTandPathPE = 0;
			}
		}else
		{
			contigPath->validCandPathTandPathFlagPE = NO;
			contigPath->startRowCandPathTandPathPE = -1;
			contigPath->appearTimesCandTandPathPE = 0;
		}
	}

	if(contigPath->validCandPathTandPathFlagSE==YES)
	{
		naviPathseq = contigPath->candPathseqTandPathSE + contigPath->startRowCandPathTandPathSE;
		naviPathseqLen = contigPath->candPathLenTandPathSE - contigPath->startRowCandPathTandPathSE;
		if(naviPathseqLen>0)
		{
			switch(naviPathseq[0])
			{
				case 'A': baseIndex = 0; break;
				case 'C': baseIndex = 1; break;
				case 'G': baseIndex = 2; break;
				case 'T': baseIndex = 3; break;
				default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, naviPathseq[0]); return FAILED;
			}

			if(baseIndex==contigBaseInt)
			{
				contigPath->startRowCandPathTandPathSE ++;
			}else
			{
				contigPath->validCandPathTandPathFlagSE = NO;
				contigPath->startRowCandPathTandPathSE = -1;
				contigPath->appearTimesCandTandPathSE = 0;
			}
		}else
		{
			contigPath->validCandPathTandPathFlagSE = NO;
			contigPath->startRowCandPathTandPathSE = -1;
			contigPath->appearTimesCandTandPathSE = 0;
		}
	}

	return SUCCESSFUL;
}


/**
 * Navigation by contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short decideByContigPath(int32_t *naviContigPath, contigPath_t *contigPath, int32_t *occsNumArray, int32_t *occsNumIndexArray, graphtype *graph, double occRatioThreshold)
{
	int32_t i, baseIndex, baseIndexMax, validFlag, naviPathseqLen, maxPathseqLen, secPathseqLen, candPathLen, shareLen, mismatchNum, shiftSize, shiftType;
	char *naviPathseq, *maxPathseq, *secPathseq, *candPathseq;
	contigPathItem_t *naviPathItem;


	if(contigPath->updateInterval>20 && contigPath->naviPathItem && ((contigPath->naviSuccessSize>=200 && (double)occsNumArray[occsNumIndexArray[1]]/occsNumArray[occsNumIndexArray[0]]<0.3) ||
		(contigPath->naviSuccessSize>=2000 && (double)occsNumArray[occsNumIndexArray[1]]/occsNumArray[occsNumIndexArray[0]]<occRatioThreshold)))  // added 2014-01-16
	{
		//if(contigPath->updateInterval>=20)
		if((contigPath->updateInterval>=20 && (double)occsNumArray[occsNumIndexArray[1]]/occsNumArray[occsNumIndexArray[0]]<0.85*occRatioThreshold) || (contigPath->preNaviOverlapSize>15 && contigPath->updateInterval<contigPath->preNaviOverlapSize) || contigPath->secPathItem==NULL) // 2014-04-21
		{
			validFlag = YES;
		}else if(contigPath->secPathItem)
		{
			maxPathseq = contigPath->maxPathItem->contigPathStr + contigPath->startRowNewBase;
			secPathseq = contigPath->secPathItem->contigPathStr + contigPath->startRowNewBase;
			maxPathseqLen = contigPath->maxPathItem->contigPathLen - contigPath->startRowNewBase;
			secPathseqLen = contigPath->secPathItem->contigPathLen - contigPath->startRowNewBase;
			if(maxPathseqLen<secPathseqLen)
				shareLen = maxPathseqLen;
			else
				shareLen = secPathseqLen;

			mismatchNum = 0;
			if(shareLen>0)
				for(i=0; i<shareLen; i++)
					if(maxPathseq[i]!=secPathseq[i])
						mismatchNum ++;

			if(mismatchNum<=contigPath->maxMismatchNumThres)
				validFlag = YES;
			else
				validFlag = NO;
		}

		if(validFlag==YES)
		{
			naviPathItem = contigPath->naviPathItem;
			naviPathseq = naviPathItem->contigPathStr + contigPath->startRowNewBase;
			naviPathseqLen = naviPathItem->contigPathLen - contigPath->startRowNewBase;
			if(naviPathseqLen>0)
			{
				switch(naviPathseq[0])
				{
					case 'A': baseIndex = 0; break;
					case 'C': baseIndex = 1; break;
					case 'G': baseIndex = 2; break;
					case 'T': baseIndex = 3; break;
					default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, naviPathseq[0]); return FAILED;
				}

				//printf("************* localContigID=%ld, contigNodesNum=%ld, itemNumDT=%d, baseIndex=%d, kmerBase=%d, kmer_len=%d, occNumArray:(%d,%d,%d,%d).\n", localContigID, itemNumContigArr, itemNumDecisionTable, baseIndex, kmerSeqIntAssembly[entriesPerKmer-1]&3, kmer_len, occsNumArray[0], occsNumArray[1], occsNumArray[2], occsNumArray[3]);
				//outputContigPath(contigPath, YES);

				if(baseIndex==occsNumIndexArray[0])
				{
					//printf("************* localContigID=%ld, contigNodesNum=%ld, itemNumDT=%d, baseIndex=%d, kmerBase=%d, kmer_len=%d, occNumArray:(%d,%d,%d,%d).\n", localContigID, itemNumContigArr, itemNumDecisionTable, baseIndex, kmerSeqIntAssembly[entriesPerKmer-1]&3, kmer_len, occsNumArray[0], occsNumArray[1], occsNumArray[2], occsNumArray[3]);
					//outputContigPath(contigPath, YES);

					*naviContigPath = NAVI_SUCCESS;

					kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] >> 2) << 2) | baseIndex;

					kmers[0] = getKmer(kmerSeqIntAssembly, graph);
					kmers[1] = getReverseKmer(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, graph);
				}else
				{
					*naviContigPath = NAVI_FAILED;
				}
			}else
			{
				*naviContigPath = NAVI_FAILED;
			}
		}else
		{
			*naviContigPath = NAVI_FAILED;
		}
	}else
	{
		*naviContigPath = NAVI_FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get the mismatch size by shift operation.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMismatchNumByShiftOp(int32_t *mismatchNum, int32_t *shiftSize, int32_t *shiftType, char *seq1, int32_t seqLen1, char *seq2, int32_t seqLen2)
{
	int32_t i, j, startPos, startPos2, shareLen, mismatchNumTmp, seqLenTmp, validFlag;
	char *seqTmp;

	*shiftType = -1;
	*mismatchNum = INT_MAX;
	*shiftSize = 0;

	mismatchNumTmp = INT_MAX;
	validFlag = NO;
	for(startPos=0; startPos<=10; startPos++)
	{
		seqTmp = seq1 + startPos;
		seqLenTmp = seqLen1 - startPos;
		if(seqLenTmp<seqLen2)
			shareLen = seqLenTmp;
		else
			shareLen = seqLen2;

		if(shareLen>0)
		{
			mismatchNumTmp = 0;
			for(i=0; i<shareLen; i++)
				if(seqTmp[i]!=seq2[i])
					mismatchNumTmp ++;
		}else
		{
			mismatchNumTmp = INT_MAX;
		}

		if(mismatchNumTmp<=1)
		{
			validFlag = YES;
			*mismatchNum = mismatchNumTmp;
			*shiftSize = startPos;
			*shiftType = 1;
			break;
		}
	}

	if(validFlag==YES)
	{
		return SUCCESSFUL;
	}else
	{
		mismatchNumTmp = INT_MAX;
		validFlag = NO;
		for(startPos=0; startPos<=10; startPos++)
		{
			seqTmp = seq2 + startPos;
			seqLenTmp = seqLen2 - startPos;
			if(seqLenTmp<seqLen1)
				shareLen = seqLenTmp;
			else
				shareLen = seqLen1;

			if(shareLen>0)
			{
				mismatchNumTmp = 0;
				for(i=0; i<shareLen; i++)
					if(seq1[i]!=seqTmp[i])
						mismatchNumTmp ++;
			}else
			{
				mismatchNumTmp = INT_MAX;
			}

			if(mismatchNumTmp<=1)
			{
				validFlag = YES;
				*mismatchNum = mismatchNumTmp;
				*shiftSize = startPos;
				*shiftType = 2;
				break;
			}
		}
	}

	return SUCCESSFUL;
}

