/*
 * candPath2.c
 *
 *  Created on: Nov 15, 2013
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Decide the extension by the candidate paths for paired-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short decideByCandPathPE(int32_t *naviCandPathPE, int32_t *maxBaseIndexAfterCandPathPE, int32_t *incorrectBaseNumCandPathPE, int32_t *occNumPE, int32_t *occsNumIndexPE, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, contigPath_t *contigPath)
{
	int32_t *baseNumArray, maxRowNumBaseNumArray, colsNum, maxIndex;
	char base;


	// get the candidate paths
	if(getCandPathPE(candPath, decisionTable, readsNumDecisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the candidate paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(candPath->itemNumCandPathItemArray>1)
	{
		// initialize the memory
		colsNum = 8;
		maxRowNumBaseNumArray = candPath->maxPathLen;
		baseNumArray = (int32_t *)calloc(maxRowNumBaseNumArray*colsNum, sizeof(int32_t));
		if(baseNumArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute the base number array
		if(computeBaseNumArrayCandPath(baseNumArray, maxRowNumBaseNumArray, 8, candPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot adjust the reads for candidate paths, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute the number of incorrect bases
		if(computeIncorrectBaseNumCandPath(incorrectBaseNumCandPathPE, baseNumArray, maxRowNumBaseNumArray, 8)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the number of incorrect bases for candidate paths, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// decide navigation by incorrect base number
		if(decideNaviCandPathPE(naviCandPathPE, &maxIndex, occNumPE, occsNumIndexPE, *incorrectBaseNumCandPathPE, baseNumArray, maxRowNumBaseNumArray, 8, contigPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot navigation from candidate paths, error!\n", __LINE__, __func__);
			return FAILED;
		}


		if((*naviCandPathPE)==NAVI_FAILED) // 2014-01-29
		{
			if(candPath->itemNumCandPathItemArray==2)
			{
				if(confirmNaviByShiftOpCandPath(naviCandPathPE, &maxIndex, incorrectBaseNumCandPathPE, candPath, contigPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot navigation from shifted candidate path, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}

		// chech the contig path
		if((*naviCandPathPE)==NAVI_SUCCESS)
		{
			if(confirmNaviCandPathByContigPath(&maxIndex, occNumPE, occsNumIndexPE, contigPath)==FAILED)
			{
				printf("line=%d, In %s(), cannot confirm navigation from contig path, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		free(baseNumArray);
	}else
	{
		*naviCandPathPE = NAVI_SUCCESS;
		*incorrectBaseNumCandPathPE = 0;
		switch(candPath->candPathItemArray[0].candPathStr[0])
		{
			case 'A': maxIndex = 0; break;
			case 'C': maxIndex = 1; break;
			case 'G': maxIndex = 2; break;
			case 'T': maxIndex = 3; break;
			default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, candPath->candPathItemArray[0].candPathStr[0]); return FAILED;
		}
	}

	if((*naviCandPathPE)==NAVI_SUCCESS)
	{
		kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] >> 2) << 2) | maxIndex;

		kmers[0] = getKmer(kmerSeqIntAssembly, deBruijnGraph);
		kmers[1] = getReverseKmer(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, deBruijnGraph);

		*maxBaseIndexAfterCandPathPE = maxIndex;
	}else
	{
		*maxBaseIndexAfterCandPathPE = -1;
	}


#if(CANDPATH_OUTPUT==YES)
	// output the paths
	outputCandPath(candPath);
	printf("******** localContigID=%ld, contigNodesNum=%ld, naviCandPathPE=%d, incorrectBaseNumCandPathPE=%d\n", localContigID, itemNumContigArr, *naviCandPathPE, *incorrectBaseNumCandPathPE);
#endif

	return SUCCESSFUL;
}

/**
 * Decide the extension by the candidate paths for single-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short decideByCandPathSE(int32_t *naviCandPathSE, int32_t *maxBaseIndexAfterCandPathSE, int32_t *incorrectBaseNumCandPathSE, int32_t *occNumSE, int32_t *occsNumIndexSE, int32_t *occsNumPE, int32_t *occsNumIndexPE, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, contigPath_t *contigPath)
{
	int32_t *baseNumArray, maxRowNumBaseNumArray, colsNum, maxIndex;

	// get the candPath
	if(getCandPathSE(candPath, decisionTable, readsNumDecisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the candidate paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(candPath->itemNumCandPathItemArray>1)
	{
		// initialize the memory
		colsNum = 8;
		maxRowNumBaseNumArray = candPath->maxPathLen;
		baseNumArray = (int32_t *)calloc(maxRowNumBaseNumArray*colsNum, sizeof(int32_t));
		if(baseNumArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute the base number array
		if(computeBaseNumArrayCandPath(baseNumArray, maxRowNumBaseNumArray, 8, candPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot adjust the reads for candidate paths, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute the number of incorrect bases
		if(computeIncorrectBaseNumCandPath(incorrectBaseNumCandPathSE, baseNumArray, maxRowNumBaseNumArray, 8)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the number of incorrect bases for candidate paths, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// decide navigation by incorrect base number
		if(decideNaviCandPathSE(naviCandPathSE, &maxIndex, occNumSE, occsNumIndexSE, *incorrectBaseNumCandPathSE, occsNumPE, occsNumIndexPE, baseNumArray, maxRowNumBaseNumArray, 8, contigPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot navigation from candidate paths, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if((*naviCandPathSE)==NAVI_FAILED) // 2014-03-01
		{
			if(candPath->itemNumCandPathItemArray==2)
			{
				if(confirmNaviByShiftOpCandPath(naviCandPathSE, &maxIndex, incorrectBaseNumCandPathSE, candPath, contigPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot navigation from shifted candidate path, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}

		// chech the contig path
		if((*naviCandPathSE)==NAVI_SUCCESS)
		{
			if(confirmNaviCandPathByContigPath(&maxIndex, occNumSE, occsNumIndexSE, contigPath)==FAILED)
			{
				printf("line=%d, In %s(), cannot confirm navigation from contig path, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		free(baseNumArray);
	}else
	{
		*naviCandPathSE = NAVI_SUCCESS;
		*incorrectBaseNumCandPathSE = 0;
		switch(candPath->candPathItemArray[0].candPathStr[0])
		{
			case 'A': maxIndex = 0; break;
			case 'C': maxIndex = 1; break;
			case 'G': maxIndex = 2; break;
			case 'T': maxIndex = 3; break;
			default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, candPath->candPathItemArray[0].candPathStr[0]); return FAILED;
		}
	}

	if((*naviCandPathSE)==NAVI_SUCCESS)
	{
		kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] >> 2) << 2) | maxIndex;

		kmers[0] = getKmer(kmerSeqIntAssembly, deBruijnGraph);
		kmers[1] = getReverseKmer(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, deBruijnGraph);

		*maxBaseIndexAfterCandPathSE = maxIndex;
	}else
	{
		*maxBaseIndexAfterCandPathSE = -1;
	}


#if(CANDPATH_OUTPUT==YES)
	// output the paths
	outputCandPath(candPath);
	printf("******** localContigID=%ld, contigNodesNum=%ld, naviCandPathSE=%d, incorrectBaseNumCandPathSE=%d\n", localContigID, itemNumContigArr, *naviCandPathSE, *incorrectBaseNumCandPathSE);
#endif

	return SUCCESSFUL;
}

/**
 * Initialize the global memory for candPath.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemCandPath(candPath_t **candPath)
{
	int32_t i;

	*candPath = (candPath_t*) calloc(1, sizeof(candPath_t));
	if((*candPath)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	(*candPath)->maxPathLen = 0;
	(*candPath)->maxCandPathItemLen = MAX_CANDIDATE_PATH_LEN;
	(*candPath)->maxItemNumCandPathItemArray = MAX_CANDIDATE_PATH_NUM;
	(*candPath)->itemNumCandPathItemArray = 0;

	// candPath
	(*candPath)->candPathItemArray = (candPathItem_t *) calloc ((*candPath)->maxItemNumCandPathItemArray, sizeof(candPathItem_t));
	if((*candPath)->candPathItemArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	for(i=0; i<(*candPath)->maxItemNumCandPathItemArray; i++)
	{
		(*candPath)->candPathItemArray[i].candPathStr = (char *) calloc ((*candPath)->maxCandPathItemLen+1, sizeof(char));
		if((*candPath)->candPathItemArray[i].candPathStr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Release the global memory for candPath.
 */
void releaseMemCandPath(candPath_t **candPath)
{
	int32_t i;

	for(i=0; i<(*candPath)->maxItemNumCandPathItemArray; i++)
		free((*candPath)->candPathItemArray[i].candPathStr);
	free((*candPath)->candPathItemArray);
	free(*candPath);
	*candPath = NULL;
}

/**
 * Get the candPath for paired-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getCandPathPE(candPath_t *candPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable)
{
	candPath->itemNumCandPathItemArray = 0;
	candPath->maxPathLen = 0;

	// get candPath
	if(getPathItemsFromDTCandPathPE(candPath, decisionTable, readsNumDecisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the candPath items from decision table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// adjust
	if(adjustCandPath(candPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust the candPath items, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get the candPath for single-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getCandPathSE(candPath_t *candPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable)
{
	candPath->itemNumCandPathItemArray = 0;
	candPath->maxPathLen = 0;

	// get candPath
	if(getPathItemsFromDTCandPathSE(candPath, decisionTable, readsNumDecisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the candPath items from decision table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// adjust
	if(adjustCandPath(candPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust the candPath items, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get path items from decision table for candPath for paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getPathItemsFromDTCandPathPE(candPath_t *candPath, assemblingreadtype *decisionTable, int32_t readsNumDT)
{
	int32_t i, newStartPos, readseqLenTmp, matchRow;
	char readseqTmp[MAX_CANDIDATE_PATH_LEN+1];

	// get the sequences from decision table
	for(i=0; i<readsNumDT; i++)
	{
		if(decisionTable[i].matedFlag==YES
			&& (decisionTable[i].unmatchBaseNum<=maxUnmatchBaseNumPerRead)
			&& (decisionTable[i].successiveAppearBases>=minSuccessiveAppearedBaseNum)
			&& decisionTable[i].basePos==decisionTable[i].lastMatchedBasePos)
		{
			if(decisionTable[i].orientation==ORIENTATION_PLUS)
			{
				newStartPos = decisionTable[i].basePos + 1;
				readseqLenTmp = decisionTable[i].seqlen - newStartPos;
				getReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, newStartPos, readseqLenTmp);
			}else
			{
				newStartPos = decisionTable[i].seqlen - decisionTable[i].basePos;
				readseqLenTmp = decisionTable[i].seqlen - newStartPos;
				getReverseReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, newStartPos, readseqLenTmp);
			}

			if(readseqLenTmp>0)
			{
				// get the matched candPathItem
				if(getMatchRowCandPathItem(&matchRow, readseqTmp, readseqLenTmp, candPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the matched row in candPath, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(matchRow>=0)
				{ // add the read if it matches
					if(addReadseqToCandPathItem(readseqTmp, readseqLenTmp, 1, matchRow, candPath)==FAILED)
					{
						printf("line=%d, In %s(), cannot add the read to candPath, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}else
				{ // generate a new candPath item
					if(addNewCandPathItem(candPath, readseqTmp, readseqLenTmp, 1)==FAILED)
					{
						printf("line=%d, In %s(), cannot add the contig path item, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Get path items from decision table for candPath for single ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getPathItemsFromDTCandPathSE(candPath_t *candPath, assemblingreadtype *decisionTable, int32_t readsNumDT)
{
	int32_t i, newStartPos, readseqLenTmp, matchRow;
	char readseqTmp[MAX_CANDIDATE_PATH_LEN+1];

	// get the sequences from decision table
	for(i=0; i<readsNumDT; i++)
	{
		if((decisionTable[i].unmatchBaseNum<=maxUnmatchBaseNumPerRead)
			&& (decisionTable[i].successiveAppearBases>=minSuccessiveAppearedBaseNum)
			&& decisionTable[i].basePos==decisionTable[i].lastMatchedBasePos)
		{
			if(decisionTable[i].orientation==ORIENTATION_PLUS)
			{
				newStartPos = decisionTable[i].basePos + 1;
				readseqLenTmp = decisionTable[i].seqlen - newStartPos;
				getReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, newStartPos, readseqLenTmp);
			}else
			{
				newStartPos = decisionTable[i].seqlen - decisionTable[i].basePos;
				readseqLenTmp = decisionTable[i].seqlen - newStartPos;
				getReverseReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, newStartPos, readseqLenTmp);
			}

			if(readseqLenTmp>0)
			{
				// get the matched candPathItem
				if(getMatchRowCandPathItem(&matchRow, readseqTmp, readseqLenTmp, candPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the matched row in candPath, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(matchRow>=0)
				{ // add the read if it matches
					if(addReadseqToCandPathItem(readseqTmp, readseqLenTmp, 1, matchRow, candPath)==FAILED)
					{
						printf("line=%d, In %s(), cannot add the read to candPath, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}else
				{ // generate a new candPath item
					if(addNewCandPathItem(candPath, readseqTmp, readseqLenTmp, 1)==FAILED)
					{
						printf("line=%d, In %s(), cannot add the contig path item, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the match row in contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMatchRowCandPathItem(int32_t *matchRow, char *readseqTmp, int32_t readseqLenTmp, candPath_t *candPath)
{
	int32_t i, j, itemNum, misMatchNum, minMismatchNum, pathLen, shareLen;
	char *pathseq;
	candPathItem_t *pathItemArray;

	minMismatchNum = INT_MAX;
	pathItemArray = candPath->candPathItemArray;
	itemNum = candPath->itemNumCandPathItemArray;
	for(i=0; i<itemNum; i++)
	{
		pathseq = pathItemArray[i].candPathStr;
		pathLen = pathItemArray[i].pathLen;
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
					if(misMatchNum>0)
						break;
				}
			}

			if(misMatchNum<minMismatchNum)
			{
				minMismatchNum = misMatchNum;
				*matchRow = i;
			}
		}

		if(minMismatchNum==0)
			break;
	}

	if(minMismatchNum>0)
		*matchRow = -1;

	return SUCCESSFUL;
}

/**
 * Add the read to the matched candPath item.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addReadseqToCandPathItem(char *readseqTmp, int32_t readseqLenTmp, int32_t supportReadsNum, int32_t matchRow, candPath_t *candPath)
{
	int32_t pathLen, newBaseLen;
	char *pathseq;
	candPathItem_t *pathItem;

	pathItem = candPath->candPathItemArray + matchRow;
	pathLen = pathItem->pathLen;
	if(pathLen<readseqLenTmp)
	{ // add the new bases to the path
		newBaseLen = readseqLenTmp - pathLen;
		pathseq = pathItem->candPathStr;
		strcat(pathseq+pathLen, readseqTmp+pathLen);
		pathItem->pathLen += newBaseLen;
	}
	pathItem->supportReadsNum += supportReadsNum;

	// ############################## Debug information ##############################
	if(pathItem->pathLen!=strlen(candPath->candPathItemArray[matchRow].candPathStr))
	{
		printf("line=%d, In %s(), pathLen=%d, pathseq=%s, error!\n", __LINE__, __func__, pathItem->pathLen, candPath->candPathItemArray[matchRow].candPathStr);
		return FAILED;
	}
	// ############################## Debug information ##############################

	if(candPath->maxPathLen<pathItem->pathLen)
		candPath->maxPathLen = pathItem->pathLen;

	return SUCCESSFUL;
}

/**
 * Add the read to the new candPath item.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewCandPathItem(candPath_t *candPath, char *readseqTmp, int32_t readseqLenTmp, int32_t supportReadsNum)
{
	candPathItem_t *pathItemArray;

	if(candPath->itemNumCandPathItemArray<candPath->maxItemNumCandPathItemArray)
	{
		pathItemArray = candPath->candPathItemArray;
		strcpy(pathItemArray[candPath->itemNumCandPathItemArray].candPathStr, readseqTmp);
		pathItemArray[candPath->itemNumCandPathItemArray].pathLen = readseqLenTmp;
		pathItemArray[candPath->itemNumCandPathItemArray].supportReadsNum = supportReadsNum;
		candPath->itemNumCandPathItemArray ++;

		if(readseqLenTmp>candPath->maxPathLen)
			candPath->maxPathLen = readseqLenTmp;
	}

	return SUCCESSFUL;
}

/**
 * Adjust the candPath items.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short adjustCandPath(candPath_t *candPath)
{
	int32_t *baseNumArray, maxRowNumBaseNumArray, colsNum;  // column=8

	if(contigPath->itemNumPathItemList>=2)
	{
		// initialize the memory
		colsNum = 8;
		maxRowNumBaseNumArray = candPath->maxPathLen;
		baseNumArray = (int32_t *)calloc(maxRowNumBaseNumArray*colsNum, sizeof(int32_t));
		if(baseNumArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// fill the baseNumArray
		if(computeBaseNumArrayCandPath(baseNumArray, maxRowNumBaseNumArray, colsNum, candPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the base number array for candPath, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// remove error base in candPath and update the path
		if(removeErrorBaseCandPath(baseNumArray, maxRowNumBaseNumArray, colsNum, candPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot remove error bases for candPath, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// merge adjusted contig path items
		if(mergeCandPathItem(candPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot merge candPath, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// free memory
		free(baseNumArray);
	}

	return SUCCESSFUL;
}

/**
 * Compute the base number array for candPath.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeBaseNumArrayCandPath(int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, candPath_t *candPath)
{
	int32_t i, j, itemNum, baseNumTmp, supportReadsNum, baseInt;
	int32_t value, maxValue, secValue, maxIndex, secIndex, subSum;
	char *pathseq;
	candPathItem_t *pathItemArray;

	// fill the base number
	pathItemArray = candPath->candPathItemArray;
	itemNum = candPath->itemNumCandPathItemArray;
	for(i=0; i<itemNum; i++)
	{
		supportReadsNum = pathItemArray[i].supportReadsNum;
		baseNumTmp = pathItemArray[i].pathLen;
		pathseq = pathItemArray[i].candPathStr;
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
 * Remove error bases for candPath.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeErrorBaseCandPath(int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, candPath_t *candPath)
{
	int32_t i, j, k, itemNum, pathLen;
	int32_t maxBaseIndex, maxValue, secondBaseIndex, secondValue, sum;
	int32_t gapValue, baseValue, errrBaseIndex, pathseqRow, tmpRow;
	char *pathseq, errBase, newBase;;
	candPathItem_t *pathItemArray;

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
						pathItemArray = candPath->candPathItemArray;
						itemNum = candPath->itemNumCandPathItemArray;
						for(j=0; j<itemNum; j++)
						{
							pathseq = pathItemArray[j].candPathStr;
							pathLen = pathItemArray[j].pathLen;
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
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Merge candPath items.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mergeCandPathItem(candPath_t *candPath)
{
	int32_t i, j, k, pathseqLen, targetseqLen, similarFlag, minLen, mismatchNum;
	char *pathseq, *targetseq;
	candPathItem_t *pathItemArray;

	pathItemArray = candPath->candPathItemArray;
	for(i=0; i<candPath->itemNumCandPathItemArray-1; i++)
	{
		for(j=candPath->itemNumCandPathItemArray-1; j>i; j--)
		{
			targetseq = pathItemArray[i].candPathStr;
			targetseqLen = pathItemArray[i].pathLen;
			pathseq = pathItemArray[j].candPathStr;
			pathseqLen = pathItemArray[j].pathLen;

			if(pathseqLen>targetseqLen)
				minLen = targetseqLen;
			else
				minLen = pathseqLen;

			if(minLen>0)
			{
				mismatchNum = 0;
				similarFlag = YES;
				for(k=0; k<minLen; k++)
				{
					if(pathseq[k]!=targetseq[k])
					{
						mismatchNum ++;

						if(mismatchNum>0)
						{
							similarFlag = NO;
							break;
						}
					}
				}
			}else
			{
				similarFlag = NO;
			}

			if(similarFlag==YES && pathItemArray[i].supportReadsNum!=pathItemArray[j].supportReadsNum)
			{
				if(pathItemArray[i].supportReadsNum>pathItemArray[j].supportReadsNum)
				{
					if(targetseqLen<pathseqLen)
					{
						strcat(targetseq, pathseq+targetseqLen);
						pathItemArray[i].pathLen += pathseqLen - targetseqLen;
					}
				}else
				{
					if(targetseqLen<pathseqLen)
					{
						strcpy(targetseq, pathseq);
						pathItemArray[i].pathLen = pathseqLen;
					}else
					{
						strncpy(targetseq, pathseq, pathseqLen);
					}
				}
				pathItemArray[i].supportReadsNum += pathItemArray[j].supportReadsNum;

				k = j + 1;
				while(k<candPath->itemNumCandPathItemArray)
				{
					strcpy(pathItemArray[k-1].candPathStr, pathItemArray[k].candPathStr);
					pathItemArray[k-1].pathLen = pathItemArray[k].pathLen;
					pathItemArray[k-1].supportReadsNum = pathItemArray[k].supportReadsNum;
					k ++;
				}
				candPath->itemNumCandPathItemArray --;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Output the candPath.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
void outputCandPath(candPath_t *candPath)
{
	int32_t i;

	if(candPath->itemNumCandPathItemArray>0)
	{
		printf("%d candidate paths, maxPathLen=%d:\n", candPath->itemNumCandPathItemArray, candPath->maxPathLen);
		for(i=0; i<candPath->itemNumCandPathItemArray; i++)
			printf("%d\t%s, len=%d\n", candPath->candPathItemArray[i].supportReadsNum, candPath->candPathItemArray[i].candPathStr, candPath->candPathItemArray[i].pathLen);
	}
	else
	{
		printf("There are no paths.\n");
	}
}

/**
 * Compute the incorrect base numbers for candidate path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeIncorrectBaseNumCandPath(int32_t *incorrectNum, int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum)
{
	int32_t i, maxValue, secValue, maxIndex, secIndex, subSum;

	// compute the incorrect base number
	*incorrectNum = 0;
	for(i=0; i<maxRowNumBaseNumArray; i++)
	{
		maxIndex = baseNumArray[i*colsNum + 6];
		secIndex = baseNumArray[i*colsNum + 7];

		maxValue = baseNumArray[i*colsNum + maxIndex];
		secValue = baseNumArray[i*colsNum + secIndex];
		subSum = baseNumArray[i*colsNum + 5];

		//if(itemNumMalignArray>10)
		if(subSum>=10)
		{
			subSum = baseNumArray[i*colsNum + 5];
			if((double)maxValue/subSum<INCOR_RATIO_CANDPATH)
				(*incorrectNum) ++;
		}else
		{
			if(secValue>0 && (double)maxValue/secValue<4)
				(*incorrectNum) ++;
		}
	}

	return SUCCESSFUL;
}


/**
 * Decide the navigation from the incorrect base numbers for candidate paths for paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short decideNaviCandPathPE(int32_t *naviCandPath, int32_t *maxIndex, int32_t *occNumPE, int32_t *occsNumIndexPE, int32_t incorrectNum, int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, contigPath_t *contigPath)
{
	int32_t i, firstBaseRow, maxBaseIndex, maxValue, secondBaseIndex, secondValue, sum;
	int32_t lowNum, baseIndex;
	char *pathseq;

	if(incorrectNum>MAX_INCOR_BASE_NUM_CANDPATH)
	{
		*maxIndex = -1;
		*naviCandPath = NAVI_FAILED;

/*
		if(incorrectNum>=2*MAX_INCOR_BASE_NUM_CANDPATH)
		{
			*maxIndex = -1;
			*naviCandPath = NAVI_FAILED;
			//printf("line=%d, In %s(), naviCandPathPE=stop!\n", __LINE__, __func__);
		}else
		{
			// skip the gaps
			for(i=0; i<maxRowNumBaseNumArray; i++)
			{
				maxBaseIndex = baseNumArray[i*colsNum + 6];
				secondBaseIndex = baseNumArray[i*colsNum + 7];

				if(maxBaseIndex<4)
				{
					firstBaseRow = i;
					break;
				}
			}

			lowNum = 0;
			for(i=firstBaseRow; i<maxRowNumBaseNumArray; i++)
			{
				maxBaseIndex = baseNumArray[i*colsNum + 6];
				secondBaseIndex = baseNumArray[i*colsNum + 7];
				maxValue = baseNumArray[i*colsNum + maxBaseIndex];
				secondValue = baseNumArray[i*colsNum + secondBaseIndex];
				sum = baseNumArray[i*colsNum + 5];

				if((double)maxValue/sum<0.6)
					lowNum ++;
			}

			if(lowNum<2)
			{
				maxBaseIndex = baseNumArray[firstBaseRow*colsNum + 6];
				*maxIndex = maxBaseIndex;
				*naviCandPath = NAVI_SUCCESS;
				//printf("line=%d, In %s(), naviCandPathPE=continue!\n", __LINE__, __func__);
			}else
			{
				*maxIndex = -1;
				*naviCandPath = NAVI_FAILED;
				//printf("line=%d, In %s(), naviCandPathPE=stop!\n", __LINE__, __func__);
			}
		}
*/
	}else
	{
		// skip the gaps
		for(i=0; i<maxRowNumBaseNumArray; i++)
		{
			maxBaseIndex = baseNumArray[i*colsNum + 6];
			secondBaseIndex = baseNumArray[i*colsNum + 7];

			if(maxBaseIndex<4)
			{
				firstBaseRow = i;
				break;
			}
		}

		maxValue = baseNumArray[firstBaseRow*colsNum + maxBaseIndex];
		secondValue = baseNumArray[firstBaseRow*colsNum + secondBaseIndex];
		sum = baseNumArray[firstBaseRow*colsNum + 5];

		if(occsNumIndexPE[0]!=-1 && maxBaseIndex!=occsNumIndexPE[0])
		{
			*maxIndex = -1;
			*naviCandPath = NAVI_FAILED;
			//printf("line=%d, In %s(), naviCandPathPE=stop!\n", __LINE__, __func__);
		}else
		{
			//if((double)secondValue/maxValue>0.8 && incorrectNum>1 && sum>10)
			if((double)secondValue/maxValue>0.7 && incorrectNum>1 && sum>10)  // best
			//if((double)secondValue/maxValue>0.85 && incorrectNum>1 && sum>10)  // 2014-01-19
			{
				*maxIndex = -1;
				*naviCandPath = NAVI_FAILED;
				//printf("line=%d, In %s(), naviCandPathPE=stop!\n", __LINE__, __func__);
			}
			else if(((double)maxValue/sum<MAX_OCC_RATIO_CANDPATH) && incorrectNum>2)
			//else if(((double)maxValue/sum<MAX_OCC_RATIO_CANDPATH) && incorrectNum>3)
			{
				*maxIndex = -1;
				*naviCandPath = NAVI_FAILED;
				//printf("line=%d, In %s(), naviCandPathPE=stop!\n", __LINE__, __func__);
			}
/*
			else if((double)secondValue/maxValue>=MAX_OCC_RATIO_CANDPATH)  // added 2013-10-22
			{
				lowNum = 0;
				for(j=0; j<maxRowNumBaseNumArray; j++)
				{
					maxBaseIndex = baseNumArray[j*colsNum + 6];
					//secondBaseIndex = baseNumArray[j*colsNum + 7];
					maxValue = baseNumArray[j*colsNum + maxBaseIndex];
					//secondValue = baseNumArray[j*colsNum + secondBaseIndex];
					sum = baseNumArray[j*colsNum + 5];
					if((double)maxValue/sum<0.7)
						lowNum ++;
				}

				if(lowNum>=2)
				{
					*maxIndex = -1;
					*naviCandPath = NAVI_FAILED;
				}else
				{
					maxBaseIndex = baseNumArray[firstBaseRow*colsNum + 6];
					*maxIndex = maxBaseIndex;
					*naviCandPath = NAVI_SUCCESS;
				}
			}
*/
			else
			{
				*maxIndex = maxBaseIndex;
				*naviCandPath = NAVI_SUCCESS;
				//printf("line=%d, In %s(), naviCandPathPE=continue!\n", __LINE__, __func__);
			}
		}

		// adjust the maxIndex, 2014-02-04
		if((*naviCandPath)==NAVI_SUCCESS && secondValue>0 && (secondValue>0.8*maxValue || maxValue-secondValue<2))
		{
//			if((occsNumIndexPE[1]!=-1 && occNumPE[occsNumIndexPE[0]]<5) && contigPath->naviPathItem==NULL && contigPath->validCandPathTandPathFlagPE==YES && contigPath->appearTimesCandTandPathPE>=3 && contigPath->candPathLenTandPathPE-contigPath->startRowCandPathTandPathPE>0) // 2014-02-24, contig 476
//			{
//				pathseq = contigPath->candPathseqTandPathPE + contigPath->startRowCandPathTandPathPE;
//				switch(pathseq[0])
//				{
//					case 'A': baseIndex = 0; break;
//					case 'C': baseIndex = 1; break;
//					case 'G': baseIndex = 2; break;
//					case 'T': baseIndex = 3; break;
//					default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, pathseq[0]); return FAILED;
//				}
//				*maxIndex = baseIndex;
//			}
//			else
			//if(contigPath->naviPathItem && contigPath->naviPathItem->contigPathLen-contigPath->startRowNewBase>0)
			if(contigPath->naviPathItem && contigPath->naviPathItem->contigPathLen-contigPath->startRowNewBase>1) // 2014-03-12
			{
				pathseq = contigPath->naviPathItem->contigPathStr + contigPath->startRowNewBase;
				switch(pathseq[0])
				{
					case 'A': baseIndex = 0; break;
					case 'C': baseIndex = 1; break;
					case 'G': baseIndex = 2; break;
					case 'T': baseIndex = 3; break;
					default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, pathseq[0]); return FAILED;
				}
				if(baseIndex==secondBaseIndex)
					*maxIndex = baseIndex;
			}
		}

	}

	return SUCCESSFUL;
}

/**
 * Decide the navigation from the incorrect base numbers for candidate paths for single-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short decideNaviCandPathSE(int32_t *naviCandPath, int32_t *maxIndex, int32_t *occNumSE, int32_t *occsNumIndexSE, int32_t incorrectNum, int32_t *occNumPE, int32_t *occNumIndexPE, int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, contigPath_t *contigPath)
{
	int32_t i, j, firstBaseRow, maxBaseIndex, maxValue, secondBaseIndex, secondValue, sum;
	int32_t lowNum;
	int32_t maxPEIndex, secPEIndex, maxPEValue, secPEValue, baseIndex;
	char *pathseq;

	if(incorrectNum>MAX_INCOR_BASE_NUM_CANDPATH)
	{
		//*maxIndex = -1;
		//*naviCandPath = NAVI_FAILED;
		//printf("line=%d, In %s(), naviCandPathSE=stop!\n", __LINE__, __func__);

		if(incorrectNum>=2*MAX_INCOR_BASE_NUM_CANDPATH)
		{
			*maxIndex = -1;
			*naviCandPath = NAVI_FAILED;
			//printf("line=%d, In %s(), naviCandPathPE=stop!\n", __LINE__, __func__);
		}else
		{
			// skip the gaps
			for(i=0; i<maxRowNumBaseNumArray; i++)
			{
				maxBaseIndex = baseNumArray[i*colsNum + 6];
				secondBaseIndex = baseNumArray[i*colsNum + 7];

				if(maxBaseIndex<4)
				{
					firstBaseRow = i;
					break;
				}
			}

			lowNum = 0;
			for(i=firstBaseRow; i<maxRowNumBaseNumArray; i++)
			{
				maxBaseIndex = baseNumArray[i*colsNum + 6];
				secondBaseIndex = baseNumArray[i*colsNum + 7];
				maxValue = baseNumArray[i*colsNum + maxBaseIndex];
				secondValue = baseNumArray[i*colsNum + secondBaseIndex];
				sum = baseNumArray[i*colsNum + 5];

				if((double)maxValue/sum<0.6)
					lowNum ++;
			}

			if(lowNum<2)
			{
				maxBaseIndex = baseNumArray[firstBaseRow*colsNum + 6];
				*maxIndex = maxBaseIndex;
				*naviCandPath = NAVI_SUCCESS;
				//printf("line=%d, In %s(), naviCandPathPE=continue!\n", __LINE__, __func__);
			}else
			{
				*maxIndex = -1;
				*naviCandPath = NAVI_FAILED;
				//printf("line=%d, In %s(), naviCandPathPE=stop!\n", __LINE__, __func__);
			}
		}
	}else
	{
		// skip the gaps
		for(i=0; i<maxRowNumBaseNumArray; i++)
		{
			maxBaseIndex = baseNumArray[i*colsNum + 6];
			secondBaseIndex = baseNumArray[i*colsNum + 7];

			if(maxBaseIndex<4)
			{
				firstBaseRow = i;
				break;
			}
		}

		maxValue = baseNumArray[firstBaseRow*colsNum + maxBaseIndex];
		secondValue = baseNumArray[firstBaseRow*colsNum + secondBaseIndex];
		sum = baseNumArray[firstBaseRow*colsNum + 5];

		if(occNumIndexPE!=NULL)
		{
			maxPEIndex = occNumIndexPE[0];
			secPEIndex = occNumIndexPE[1];
			maxPEValue = occNumPE[maxPEIndex];
			secPEValue = occNumPE[secPEIndex];
		}else
		{
			maxPEIndex = secPEIndex = -1;
			maxPEValue = secPEValue = 0;
		}

		if(occNumSE[occsNumIndexSE[1]]/occNumSE[occsNumIndexSE[0]]>=MAX_OCC_RATIO_CANDPATH && (double)secondValue/maxValue>=MAX_OCC_RATIO_CANDPATH && maxBaseIndex!=occsNumIndexSE[0])
		{ // added 2013-10-21
			*maxIndex = -1;
			*naviCandPath = NAVI_FAILED;
			//printf("line=%d, In %s(), naviCandPathSE=stop!\n", __LINE__, __func__);
		}
		//else if((maxPEIndex!=-1 && (maxBaseIndex!=maxPEIndex || (maxPEValue>0 && (double)secPEValue/maxPEValue>0.8))) || maxBaseIndex!=occsNumIndexSE[0])
		//else if((maxPEIndex!=-1 && ((maxBaseIndex!=maxPEIndex && maxPEValue-secPEValue>2 && maxBaseIndex!=occsNumIndexSE[0]) || (maxPEValue>0 && (double)secPEValue/maxPEValue>0.8))) || maxBaseIndex!=occsNumIndexSE[0])
		else if((maxPEIndex!=-1 && (incorrectNum>=2 && maxValue<5*secondValue) && ((maxBaseIndex!=maxPEIndex && maxPEValue-secPEValue>2 && maxBaseIndex!=occsNumIndexSE[0]) || (maxPEValue>0 && (double)secPEValue/maxPEValue>0.8))) || maxBaseIndex!=occsNumIndexSE[0])
		{
			*maxIndex = -1;
			*naviCandPath = NAVI_FAILED;
			//printf("line=%d, In %s(), naviCandPathSE=stop!\n", __LINE__, __func__);
		}
		else
		{
			//if((double)secondValue/maxValue>0.8 && incorrectNum>1 && sum>10)
			if((double)secondValue/maxValue>0.7 && incorrectNum>1 && sum>10)  // best
			//if((double)secondValue/maxValue>0.85 && incorrectNum>1 && sum>10)  // 2014-01-19
			{
				*maxIndex = -1;
				*naviCandPath = NAVI_FAILED;
				//printf("line=%d, In %s(), naviCandPathSE=stop!\n", __LINE__, __func__);
			}
			else if(((double)maxValue/sum<MAX_OCC_RATIO_CANDPATH) && incorrectNum>2)
			//else if(((double)maxValue/sum<MAX_OCC_RATIO_CANDPATH) && incorrectNum>3)
			{
				*maxIndex = -1;
				*naviCandPath = NAVI_FAILED;
				//printf("line=%d, In %s(), naviCandPathSE=stop!\n", __LINE__, __func__);
			}
			else if((double)secondValue/maxValue>=MAX_OCC_RATIO_CANDPATH)  // added 2013-10-22
			{
				lowNum = 0;
				for(j=0; j<maxRowNumBaseNumArray; j++)
				{
					maxBaseIndex = baseNumArray[j*colsNum + 6];
					//secondBaseIndex = baseNumArray[j*colsNum + 7];
					maxValue = baseNumArray[j*colsNum + maxBaseIndex];
					//secondValue = baseNumArray[j*colsNum + secondBaseIndex];
					sum = baseNumArray[j*colsNum + 5];
					if((double)maxValue/sum<0.7)
						lowNum ++;
				}

				if(lowNum>=2)
				{
					*maxIndex = -1;
					*naviCandPath = NAVI_FAILED;
					//printf("line=%d, In %s(), naviCandPathSE=stop!\n", __LINE__, __func__);
				}else
				{
					maxBaseIndex = baseNumArray[firstBaseRow*colsNum + 6];
					*maxIndex = maxBaseIndex;
					*naviCandPath = NAVI_SUCCESS;
					//printf("line=%d, In %s(), naviCandPathSE=continue!\n", __LINE__, __func__);
				}
			}
			else
			{
				*maxIndex = maxBaseIndex;
				*naviCandPath = NAVI_SUCCESS;
				//printf("line=%d, In %s(), naviCandPathSE=continue!\n", __LINE__, __func__);
			}
		}

		// adjust the maxIndex, 2014-02-04
		if((*naviCandPath)==NAVI_SUCCESS && secondValue>0 && (secondValue>0.8*maxValue || maxValue-secondValue<2))
		{
			//if(contigPath->naviPathItem && contigPath->naviPathItem->contigPathLen-contigPath->startRowNewBase>0)
			if(contigPath->naviPathItem && contigPath->naviPathItem->contigPathLen-contigPath->startRowNewBase>1) // 2014-03-12
			{
				pathseq = contigPath->naviPathItem->contigPathStr + contigPath->startRowNewBase;
				switch(pathseq[0])
				{
					case 'A': baseIndex = 0; break;
					case 'C': baseIndex = 1; break;
					case 'G': baseIndex = 2; break;
					case 'T': baseIndex = 3; break;
					default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, pathseq[0]); return FAILED;
				}
				if(baseIndex==secondBaseIndex)
					*maxIndex = baseIndex;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Confirm the candPath navigation by shift operation.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short confirmNaviByShiftOpCandPath(int32_t *naviCandPath, int32_t *maxIndex, int32_t *incorrectBaseNum, candPath_t *candPath, contigPath_t *contigPath)
{
	char *seq1, *seq2, *pathseq, *naviPathseq;
	int32_t pathLen, seqLen1, seqLen2, mismatchNum, shiftSize, shiftType, supReadNum1, supReadNum2, baseIndex;
	int32_t i, j, mismatchNum2, naviPathLen, shareLen;

	if(candPath->itemNumCandPathItemArray==2)
	{
		seq1 = candPath->candPathItemArray[0].candPathStr;
		seq2 = candPath->candPathItemArray[1].candPathStr;
		seqLen1 = candPath->candPathItemArray[0].pathLen;
		seqLen2 = candPath->candPathItemArray[1].pathLen;
		supReadNum1 = candPath->candPathItemArray[0].supportReadsNum;
		supReadNum2 = candPath->candPathItemArray[1].supportReadsNum;

		if(getMismatchNumByShiftOp(&mismatchNum, &shiftSize, &shiftType, seq1, seqLen1, seq2, seqLen2)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the mismatch size by shift operation, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(shiftType!=-1 && mismatchNum<=3)
		{
			pathseq = NULL;
			if(contigPath->naviSuccessSize>500 && contigPath->updateInterval>20 && contigPath->naviPathItem==contigPath->maxPathItem && contigPath->naviPathItem->maxOverlapWithContig>0.7*readLen) // 2014-02-20, contig 276
			{
				naviPathseq = contigPath->naviPathItem->contigPathStr + contigPath->startRowNewBase;
				naviPathLen = contigPath->naviPathItem->contigPathLen + contigPath->startRowNewBase;
				if(naviPathLen>0)
				{
					mismatchNum2 = 0;
					for(i=0; i<2; i++)
					{
						if(candPath->candPathItemArray[i].pathLen<naviPathLen)
							shareLen = candPath->candPathItemArray[i].pathLen;
						else
							shareLen = naviPathLen;

						for(j=0; j<shareLen; j++)
						{
							if(naviPathseq[j]!=candPath->candPathItemArray[i].candPathStr[j])
							{
								mismatchNum2 ++;
								if(mismatchNum2>=1)
									break;
							}
						}

						if(mismatchNum2<=0)
						{
							pathseq = candPath->candPathItemArray[i].candPathStr;
							break;
						}
					}
				}
			}

			if(pathseq==NULL)
			{
				if(supReadNum1>supReadNum2)
					pathseq = seq1;
				else
					pathseq = seq2;
			}

			switch(pathseq[0])
			{
				case 'A': baseIndex = 0; break;
				case 'C': baseIndex = 1; break;
				case 'G': baseIndex = 2; break;
				case 'T': baseIndex = 3; break;
				default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, pathseq[0]); return FAILED;
			}

			*naviCandPath = NAVI_SUCCESS;
			*maxIndex = baseIndex;
			*incorrectBaseNum = mismatchNum;
		}else
		{
			*naviCandPath = NAVI_FAILED;
			*maxIndex = -1;
			*incorrectBaseNum = INT_MAX;
		}
	}

	return SUCCESSFUL;
}

/**
 * Confirm the candPath navigation by contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short confirmNaviCandPathByContigPath(int32_t *maxIndex, int32_t *occNumArray, int32_t *occsNumIndexArray, contigPath_t *contigPath)
{
	int32_t i, baseIndex, naviPathLen, pathLen, maxPathLen, secPathLen;
	char *naviPathseq, *pathseq, *maxPathseq, *secPathseq;
	contigPathItem_t *naviPathItem;

	if(contigPath->naviPathItem && contigPath->naviSuccessSize>2000 && contigPath->itemNumPathItemList>=2)
	{
		naviPathItem = contigPath->naviPathItem;
		naviPathseq = naviPathItem->contigPathStr + contigPath->startRowNewBase;
		naviPathLen = naviPathItem->contigPathLen - contigPath->startRowNewBase;

		if(naviPathLen>0)
		{
			switch(naviPathseq[0])
			{
				case 'A': baseIndex = 0; break;
				case 'C': baseIndex = 1; break;
				case 'G': baseIndex = 2; break;
				case 'T': baseIndex = 3; break;
				default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, naviPathseq[0]); return FAILED;
			}

			if(baseIndex!=(*maxIndex) && occNumArray[baseIndex]==occNumArray[*maxIndex] && occNumArray[occsNumIndexArray[0]]==occNumArray[occsNumIndexArray[1]])
			{
				*maxIndex = baseIndex;
			}
		}
	}
	//else if(contigPath->naviPathItem && contigPath->naviSuccessSize>100)
	else if((contigPath->naviPathItem && contigPath->naviSuccessSize>100) || contigPath->updateInterval>20) // 2014-02-21, contig 330
	{
		if(contigPath->itemNumPathItemList==1)
		{
			pathseq = contigPath->maxPathItem->contigPathStr + contigPath->startRowNewBase;
			pathLen = contigPath->maxPathItem->contigPathLen - contigPath->startRowNewBase;
			if(pathLen>0)
			{
				switch(pathseq[0])
				{
					case 'A': baseIndex = 0; break;
					case 'C': baseIndex = 1; break;
					case 'G': baseIndex = 2; break;
					case 'T': baseIndex = 3; break;
					default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, pathseq[0]); return FAILED;
				}

				if(baseIndex!=(*maxIndex))
					*maxIndex = baseIndex;
			}
		}else if(contigPath->itemNumPathItemList>=2)
		{
			if(contigPath->maxPathItem->contigPathLen-contigPath->startRowNewBase>0 && contigPath->secPathItem->contigPathLen-contigPath->startRowNewBase>0
				&& contigPath->maxPathItem->contigPathStr[contigPath->startRowNewBase]==contigPath->secPathItem->contigPathStr[contigPath->startRowNewBase])
			{
				if(contigPath->naviPathItem)
				{
					pathseq = contigPath->naviPathItem->contigPathStr + contigPath->startRowNewBase;
					pathLen = contigPath->naviPathItem->contigPathLen - contigPath->startRowNewBase;
					if(pathLen>0)
					{
						switch(pathseq[0])
						{
							case 'A': baseIndex = 0; break;
							case 'C': baseIndex = 1; break;
							case 'G': baseIndex = 2; break;
							case 'T': baseIndex = 3; break;
							default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, pathseq[0]); return FAILED;
						}

						if(baseIndex!=(*maxIndex))
							*maxIndex = baseIndex;
					}
				}else
				{
					maxPathseq = contigPath->maxPathItem->contigPathStr + contigPath->startRowNewBase;
					maxPathLen = contigPath->maxPathItem->contigPathLen - contigPath->startRowNewBase;
					secPathseq = contigPath->secPathItem->contigPathStr + contigPath->startRowNewBase;
					secPathLen = contigPath->secPathItem->contigPathLen - contigPath->startRowNewBase;
					if(maxPathLen>0 && secPathLen>0 && maxPathseq[0]==secPathseq[0])
					{
						switch(maxPathseq[0])
						{
							case 'A': baseIndex = 0; break;
							case 'C': baseIndex = 1; break;
							case 'G': baseIndex = 2; break;
							case 'T': baseIndex = 3; break;
							default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, maxPathseq[0]); return FAILED;
						}

						if(baseIndex!=(*maxIndex))
							*maxIndex = baseIndex;
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}
