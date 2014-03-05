/*
 * scafLink.c
 *
 *  Created on: Dec 18, 2013
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Link contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short contigsLinking(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph, readSet_t *readSet)
{
	printf("Constructing the scaffolds ...\n");

	// initialize the variables in contig graph
	if(initContigItemInfo(contigGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate contig edges, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// generate contig edges
	if(generateContigEdges(contigGraph, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate contig edges, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize the parameters for contigs linking
	if(initParaLinking(contigGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the parameters for contigs linking, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// link contigs
	if(linkContigs(scaffoldSet, contigGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot link contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Link contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initContigItemInfo(contigGraph_t *contigGraph)
{
	int32_t i;

	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
		contigGraph->contigItemArray[i].used = 0;
		contigGraph->contigItemArray[i].usedTimeEnd5 = 0;
		contigGraph->contigItemArray[i].usedTimeEnd3 = 0;
		contigGraph->contigItemArray[i].overlapKindEnd5 = 0;
		contigGraph->contigItemArray[i].maxOverlapLenEnd5 = 0;
		contigGraph->contigItemArray[i].overlapKindEnd3 = 0;
		contigGraph->contigItemArray[i].maxOverlapLenEnd3 = 0;

		if(contigGraph->contigItemArray[i].contigLen<=contigAlignRegSize)
			contigGraph->contigItemArray[i].onlyEnd5 = YES;
		else
			contigGraph->contigItemArray[i].onlyEnd5 = NO;

		//if(contigGraph->contigItemArray[i].contigLen<=1.5*contigAlignRegSize)
		if(contigGraph->contigItemArray[i].contigLen<=MAX_SHORT_LEN_THRES)		// added 2012-11-21
			contigGraph->contigItemArray[i].shortFlag = YES;
		else
			contigGraph->contigItemArray[i].shortFlag = NO;
	}

	return SUCCESSFUL;
}

/**
 * Generate contig edges.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateContigEdges(contigGraph_t *contigGraph, readSet_t *readSet)
{
	// initialize contig edges
	if(initContigEdges(contigGraph, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize contig edges, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the edges
	if(fillContigEdges(contigGraph, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill contig edges, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// remove invalid edges
	if(removeInvalidGraphEdge(contigGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove invalid contigGraphEdges, error!\n", __LINE__, __func__);
		return FAILED;
	}


	return SUCCESSFUL;
}

/**
 * Generate contig edges.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initContigEdges(contigGraph_t *contigGraph, readSet_t *readSet)
{
	contigGraphItem_t *contigItem;
	contigRead_t *contigReadArray;
	int32_t i, j, num, bucketArraySize, contigEdgeNum, *contigEndBucketArray;

	bucketArraySize = 2 * contigGraph->itemNumContigItemArray;
	contigEndBucketArray = (int32_t *) calloc (bucketArraySize, sizeof(int32_t));
	if(contigEndBucketArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
		contigItem = contigGraph->contigItemArray + i;

		// ######################## Debug information ########################
		//if(contigItem->contigID==94)
		//{
		//	printf("contigID=%d, contigLen=%d, endNum5=%d, endNum3=%d\n", contigItem->contigID, contigItem->contigLen, contigItem->contigReadNumEnd5, contigItem->contigReadNumEnd3);
		//}
		// ######################## Debug information ########################

		// 5' end
		if(contigItem->contigReadNumEnd5>0)
		{
			if(memset(contigEndBucketArray, 0, bucketArraySize * sizeof(int32_t))==NULL)
			{
				printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(calcEdgeNumSingleContigEnd(&contigEdgeNum, contigItem->contigID, contigItem->contigReadArrayEnd5, contigItem->contigReadNumEnd5, contigEndBucketArray, bucketArraySize, readSet)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the edge number, error!\n", __LINE__, __func__);
				return FAILED;
			}

			contigItem->itemNumContigEdgeArrayEnd5 = contigItem->contigEdgeArraySizeEnd5 = contigEdgeNum;
			if(contigEdgeNum>0)
			{
				contigItem->contigEdgeArrayEnd5 = (contigEdge_t *) calloc (contigItem->contigEdgeArraySizeEnd5, sizeof(contigEdge_t));
				if(contigItem->contigEdgeArrayEnd5==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}

				num = 0;
				for(j=0; j<bucketArraySize; j++)
				{
					if(contigEndBucketArray[j]>0)
					{
						contigItem->contigEdgeArrayEnd5[num].col = j;
						contigItem->contigEdgeArrayEnd5[num].used = NO;
						contigItem->contigEdgeArrayEnd5[num].pairedNum = contigEndBucketArray[j];
						num ++;
					}
				}
			}
		}

		// 3' end
		if(contigItem->contigReadNumEnd3>0)
		{
			if(memset(contigEndBucketArray, 0L, bucketArraySize * sizeof(int32_t))==NULL)
			{
				printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(calcEdgeNumSingleContigEnd(&contigEdgeNum, contigItem->contigID, contigItem->contigReadArrayEnd3, contigItem->contigReadNumEnd3, contigEndBucketArray, bucketArraySize, readSet)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the edge number, error!\n", __LINE__, __func__);
				return FAILED;
			}

			contigItem->itemNumContigEdgeArrayEnd3 = contigItem->contigEdgeArraySizeEnd3 = contigEdgeNum;
			if(contigEdgeNum>0)
			{
				contigItem->contigEdgeArrayEnd3 = (contigEdge_t *) calloc (contigItem->contigEdgeArraySizeEnd3, sizeof(contigEdge_t));
				if(contigItem->contigEdgeArrayEnd3==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}

				num = 0;
				for(j=0; j<bucketArraySize; j++)
				{
					if(contigEndBucketArray[j]>0)
					{
						contigItem->contigEdgeArrayEnd3[num].col = j;
						contigItem->contigEdgeArrayEnd3[num].used = NO;
						contigItem->contigEdgeArrayEnd3[num].pairedNum = contigEndBucketArray[j];
						num ++;
					}
				}
			}
		}
	}

	free(contigEndBucketArray);

	return SUCCESSFUL;
}

/**
 * Calculate the edge number for single contig end.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short calcEdgeNumSingleContigEnd(int32_t *contigEdgeNum, int32_t contigID, contigRead_t *contigReadArray, int32_t contigReadsNum, int32_t *contigEndBucketArray, int32_t bucketArraySize, readSet_t *readSet)
{
	int32_t i, endFlag, row;
	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;
	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo;
	int64_t readID, paired_readID;

	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

	for(i=0; i<contigReadsNum; i++)
	{
		readID = contigReadArray[i].readID;
		if(readID%2==1)
			paired_readID = readID + 1;
		else
			paired_readID = readID - 1;

		readMatchInfoBlockID = (paired_readID - 1) / maxItemNumPerReadMatchInfoBlock;
		rowNumInReadMatchInfoBlock = (paired_readID - 1) % maxItemNumPerReadMatchInfoBlock;
		pReadMatchInfo = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;
		endFlag = pReadMatchInfo->contigEnd;

		if(endFlag!=-1 && pReadMatchInfo->contigID!=contigID)
		{
			if(pReadMatchInfo->contigID>0)
			{
				if(endFlag==1)
					row = 2 * pReadMatchInfo->contigID - 1;
				else
					row = 2 * pReadMatchInfo->contigID - 2;

				if(row>=bucketArraySize || row<0)
				{
					printf("line=%d, In %s(), invalid row=%d, bucketArraySize=%d, error!\n", __LINE__, __func__, row, bucketArraySize);
					return FAILED;
				}

				contigEndBucketArray[row] ++;
			}
		}
	}

	*contigEdgeNum = 0;
	for(i=0; i<bucketArraySize; i++)
		if(contigEndBucketArray[i]>0)
			(*contigEdgeNum) ++;

	return SUCCESSFUL;
}

/**
 * Fill contig edges.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillContigEdges(contigGraph_t *contigGraph, readSet_t *readSet)
{
	contigEdge_t *pEdgeArray;
	int32_t i, j, col, contigID1, contigID2, endFlag1, endFlag2, edgeNum;

	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
		// 5' end
		if(contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd5>0)
		{
			contigID1 = contigGraph->contigItemArray[i].contigID;
			if(contigGraph->contigItemArray[i].onlyEnd5==NO)
				endFlag1 = 0;
			else
				endFlag1 = 2;

			pEdgeArray = contigGraph->contigItemArray[i].contigEdgeArrayEnd5;
			edgeNum = contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd5;
			for(j=0; j<edgeNum; j++)
			{
				col = pEdgeArray[j].col;
				contigID2 = col / 2 + 1;
				if((col&1)==0)
				{ // 5' contig end
					if(contigGraph->contigItemArray[contigID2-1].onlyEnd5==NO)
						endFlag2 = 0;
					else
						endFlag2 = 2;
				}else
				{ // 3' contig end
					endFlag2 = 1;
				}

				// fill the validNum[4] of contigEdge
				if(fillSituationArray(pEdgeArray+j, contigID1, contigID2, endFlag1, endFlag2, contigGraph, readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot fill situation array, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}

		// 3' end
		if(contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd3>0)
		{
			contigID1 = contigGraph->contigItemArray[i].contigID;
			endFlag1 = 1;

			pEdgeArray = contigGraph->contigItemArray[i].contigEdgeArrayEnd3;
			edgeNum = contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd3;
			for(j=0; j<edgeNum; j++)
			{
				col = pEdgeArray[j].col;
				contigID2 = col / 2 + 1;
				if((col&1)==0)
				{ // 5' contig end
					if(contigGraph->contigItemArray[contigID2-1].onlyEnd5==NO)
						endFlag2 = 0;
					else
						endFlag2 = 2;
				}else
				{ // 3' contig end
					endFlag2 = 1;
				}

				// fill the validNum[4] of contigEdge
				if(fillSituationArray(pEdgeArray+j, contigID1, contigID2, endFlag1, endFlag2, contigGraph, readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot fill situation array, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Fill the situationArray.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillSituationArray(contigEdge_t *pContigEdge, int32_t contigID1, int32_t contigID2, int32_t endFlag1, int32_t endFlag2, contigGraph_t *contigGraph, readSet_t *readSet)
{
	int64_t i, j, readID1, readID2;
	contigRead_t *contigReadArray2;
	int32_t rowNum2, orientation1, orientation2, endRead1, endRead2;
	int32_t tmp_contigID1, tmp_contigPos1;
	int32_t countArray[4];

	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;
	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo1, *pReadMatchInfo2;

	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

	for(i=0; i<4; i++) pContigEdge->validNums[i] = 0;  // reset the situation array

	// check each paired read between the two contig ends
	if(endFlag2==1)
	{
		contigReadArray2 = contigGraph->contigItemArray[contigID2-1].contigReadArrayEnd3;
		rowNum2 = contigGraph->contigItemArray[contigID2-1].contigReadNumEnd3;
	}else
	{
		contigReadArray2 = contigGraph->contigItemArray[contigID2-1].contigReadArrayEnd5;
		rowNum2 = contigGraph->contigItemArray[contigID2-1].contigReadNumEnd5;
	}

	for(i=0; i<rowNum2; i++)
	{
		readID2 = contigReadArray2[i].readID;
		orientation2 = contigReadArray2[i].orientation;
		endRead2 = contigReadArray2[i].contigEnd;

		if(endRead2==endFlag2)
		{
			readMatchInfoBlockID = (readID2 - 1) / maxItemNumPerReadMatchInfoBlock;
			rowNumInReadMatchInfoBlock = (readID2 - 1) % maxItemNumPerReadMatchInfoBlock;
			pReadMatchInfo2 = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;
			if(pReadMatchInfo2->contigID==contigID2 && pReadMatchInfo2->contigEnd==endFlag2)
			{
				// get its paired end read
				if((readID2&1)==1)
				{ // odd number
					readID1 = readID2 + 1;
				}else
				{ // even number
					readID1 = readID2 - 1;
				}

				readMatchInfoBlockID = (readID1 - 1) / maxItemNumPerReadMatchInfoBlock;
				rowNumInReadMatchInfoBlock = (readID1 - 1) % maxItemNumPerReadMatchInfoBlock;
				pReadMatchInfo1 = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;
				if(pReadMatchInfo1->contigID==contigID1 && pReadMatchInfo1->contigEnd==endFlag1)
				{
					orientation1 = pReadMatchInfo1->readOrientation;
					endRead1 = pReadMatchInfo1->contigEnd;

					tmp_contigPos1 = pReadMatchInfo1->contigPos;
					if(contigGraph->contigItemArray[contigID1-1].used==NO)
					{
						// fill the situations array
						if(endFlag1!=0 && orientation1==ORIENTATION_PLUS && endFlag2!=1 && orientation2==ORIENTATION_MINUS)
						{ // (3', +) , (5', -) ==> (+, +)
							pContigEdge->validNums[0] ++;
						}else if(endFlag1!=0 && orientation1==ORIENTATION_PLUS && endFlag2!=0 && orientation2==ORIENTATION_PLUS)
						{ // (3', +) , (3', +) ==> (+, -)
							pContigEdge->validNums[1] ++;
						}else if(endFlag1!=1 && orientation1==ORIENTATION_MINUS && endFlag2!=0 && orientation2==ORIENTATION_PLUS)
						{ // (5', -) , (3', +) ==> (-, -)
							pContigEdge->validNums[2] ++;
						}else if(endFlag1!=1 && orientation1==ORIENTATION_MINUS && endFlag2!=1 && orientation2==ORIENTATION_MINUS)
						{ // (5', -) , (5', -) ==> (-, +)
							pContigEdge->validNums[3] ++;
						}else
						{
#if DEBUG_OUT_FLAG
							printf("line=%d, In %s(), endRead1=%d, endRead2=%d, readID1=%lu, readID2=%lu, (endFlag1=%d, orientation1=%d, endFlag2=%d, orientation2=%d)\n", __LINE__, __func__, endRead1, endRead2, readID1, readID2, endFlag1, orientation1, endFlag2, orientation2);
#endif
						}
					}
				}
			}
		}
	}

	// compute the maxIndexes[4] by counting sort method
	for(i=0; i<4; i++) countArray[i] = 0;
	for(i=0; i<3; i++)
		for(j=i+1; j<4; j++)
			if(pContigEdge->validNums[i] >= pContigEdge->validNums[j])
				countArray[j] ++;
			else if(pContigEdge->validNums[i] < pContigEdge->validNums[j])
				countArray[i] ++;

	for(i=0; i<4; i++)
		pContigEdge->maxIndexes[ countArray[i] ] = i;

	for(i=0; i<4; i++)
		if(pContigEdge->validNums[ pContigEdge->maxIndexes[i] ] > 0)
			pContigEdge->totalSituationNum ++;
		else
			pContigEdge->maxIndexes[i] = -1;

	return SUCCESSFUL;
}

/**
 * Remove invalid graph edges.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeInvalidGraphEdge(contigGraph_t *contigGraph)
{
	int32_t i, top, bot;
	contigEdge_t *pEdgeArray;

	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
		// 5' end
		if(contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd5>0)
		{
			pEdgeArray = contigGraph->contigItemArray[i].contigEdgeArrayEnd5;
			top = 0;
			bot = contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd5 - 1;
			while(top<=bot)
			{
				if(pEdgeArray[top].totalSituationNum==0)
				{
					if(top<bot)
					{
						if(memcpy(pEdgeArray+top, pEdgeArray+bot, sizeof(contigEdge_t))==NULL)
						{
							printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
					bot --;
					contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd5 --;
				}else
				{
					top ++;
				}
			}

			// free the invalid memory of invalid graphEdge nodes
			if(contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd5==0)
			{
				free(contigGraph->contigItemArray[i].contigEdgeArrayEnd5);
				contigGraph->contigItemArray[i].contigEdgeArrayEnd5 = NULL;
				contigGraph->contigItemArray[i].contigEdgeArraySizeEnd5 = 0;
			}else if(contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd5<0)
			{
				printf("line=%d, In %s(), contigItemArray[%d].itemNumContigEdgeArrayEnd5=%d, error!\n", __LINE__, __func__, i, contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd5);
				return FAILED;
			}
		}

		// 3' end
		if(contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd3>0)
		{
			pEdgeArray = contigGraph->contigItemArray[i].contigEdgeArrayEnd3;
			top = 0;
			bot = contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd3 - 1;
			while(top<=bot)
			{
				if(pEdgeArray[top].totalSituationNum==0)
				{
					if(top<bot)
					{
						if(memcpy(pEdgeArray+top, pEdgeArray+bot, sizeof(contigEdge_t))==NULL)
						{
							printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
					bot --;
					contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd3 --;
				}else
				{
					top ++;
				}
			}

			// free the invalid memory of invalid graphEdge nodes
			if(contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd3==0)
			{
				free(contigGraph->contigItemArray[i].contigEdgeArrayEnd3);
				contigGraph->contigItemArray[i].contigEdgeArrayEnd3 = NULL;
				contigGraph->contigItemArray[i].contigEdgeArraySizeEnd3 = 0;
			}else if(contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd3<0)
			{
				printf("line=%d, In %s(), contigItemArray[%d].itemNumContigEdgeArrayEnd3=%d, error!\n", __LINE__, __func__, i, contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd3);
				return FAILED;
			}
		}
	}


	return SUCCESSFUL;
}

/**
 * Initialize the parameters for contigs linking.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initParaLinking(contigGraph_t *contigGraph)
{
	double averLinkNum;

	// get the average linked pairs per contigEdge
	if(computeAverPairsEachContigEdge(&averLinkNum, contigGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the average linked pairs per contigEdge, error!\n", __LINE__, __func__);
		return FAILED;
	}
	contigGraph->averLinkNum = averLinkNum;

	maxLinksNumContigsThres = averLinkNum * MAX_FIRST_LINKNUM_FACTOR;
	minLinksNumContigsThres = averLinkNum * MIN_FIRST_LINKNUM_FACTOR;
	maxRatioSecondFirstLinkNum = SECOND_FIRST_SCORE_RATIO_LINKING;
	secondLinkNumFactor = SECOND_LINKNUM_FACTOR;

	if(minLinksNumContigsThres>MIN_FIRST_LINKNUM_THRES)
	{
		minLinksNumContigsThres = MIN_FIRST_LINKNUM_THRES;
	}
	//else if(minLinksNumContigsThres<MIN_FIRST_LINKNUM_THRES)
	//{
	//	minLinksNumContigsThres = MIN_FIRST_LINKNUM_THRES;
	//}

	maxSecondLinkNumThres = minLinksNumContigsThres * secondLinkNumFactor;
	if(maxSecondLinkNumThres > MAX_SECOND_LINKNUM_THRES)
	{
		maxSecondLinkNumThres = MAX_SECOND_LINKNUM_THRES;
	}

#if (DEBUG_FLAG==YES)
	printf("maxLinksNumContigsThres=%.2f\n", maxLinksNumContigsThres);
	printf("minLinksNumContigsThres=%.2f\n", minLinksNumContigsThres);
	printf("maxSecondLinkNumThres=%.2f\n", maxSecondLinkNumThres);
	printf("maxRatioSecondFirstLinkNum=%.2f\n", maxRatioSecondFirstLinkNum);
#endif

	return SUCCESSFUL;
}

/**
 * Compute average pairs each contigEdge.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeAverPairsEachContigEdge(double *averageLinkNum, contigGraph_t *contigGraph)
{
	int64_t i, j, k;
	contigEdge_t *pEdgeArray;
	int64_t totalPairs, totalEdges;

	totalPairs = 0;
	totalEdges = 0;
	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
		// 5' end
		if(contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd5>0)
		{
			pEdgeArray = contigGraph->contigItemArray[i].contigEdgeArrayEnd5;
			for(j=0; j<contigGraph->contigItemArray[i].itemNumContigEdgeArrayEnd5; j++)
			{
				for(k=0; k<pEdgeArray[j].totalSituationNum; k++)
				{
					//if(pEdgeArray[j].validNums[ pEdgeArray[j].maxIndexes[k] ] >= minLinksNumContigsThreshold)
					if(pEdgeArray[j].validNums[ pEdgeArray[j].maxIndexes[k] ] >= 1)
					{
						totalPairs += pEdgeArray[j].validNums[ pEdgeArray[j].maxIndexes[k] ];
					}
				}

				totalEdges ++;
			}
		}
	}

	totalPairs /= 2;
	totalEdges /= 2;

	*averageLinkNum = (double)totalPairs / totalEdges;

#if (DEBUG_FLAG==YES)
	printf("totalPairs=%ld, totalEdges=%ld, averageLinkNum=%.2f\n", totalPairs, totalEdges, *averageLinkNum);
#endif

	return SUCCESSFUL;
}

/**
 * Link contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short linkContigs(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph)
{
	int32_t firstContigID;
	int32_t linkID, linkRound, linkStatus, returnCode;
	maxRowCol_t *maxRowColNode;

	maxRowColNode = (maxRowCol_t*) malloc (sizeof(maxRowCol_t));
	if(maxRowColNode==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize the contig link array
	if(initContigLinkSet(&contigLinkSet, contigGraph->itemNumContigItemArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the contig link array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// begin linking
	firstContigID = 1;
	linkID = 1;
	while(firstContigID<=contigGraph->itemNumContigItemArray)
	{
		//############################ Debug information ######################
#if DEBUG_SCAF_FLAG
		printf("============ Begin linking scaffolds: %d ============\n", linkID);
		if(linkID==3)
		{
			printf("$$$$$$$$$$$$$$$$$$$$$ linkID=%d!\n", linkID);
		}
#endif
		//############################ Debug information ######################

		linkRound = FIRST_LINK_ROUND;

		// get first contigID
		if(getFirstLinkedContigs(&firstContigID, maxRowColNode, contigGraph)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the first linked contigs in single scaffold, error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(maxRowColNode->contigID1==-1 || maxRowColNode->contigID2==-1)
			break;

		// link the first two contigs
		if(addContigToContigLinkSet(&linkStatus, contigLinkSet, contigGraph, maxRowColNode, 2, linkRound)==FAILED)
		{
			printf("line=%d, In %s(), cannot add contigs to contigLinkArr, error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(linkStatus==FAILED)
		{
			continue;
		}

		// reset the corresponding rows and columns
		if(markContigGraphEdge(contigGraph, maxRowColNode)==FAILED)
		{
			printf("line=%d, In %s(), cannot mark contig graph edge, error!\n", __LINE__, __func__);
			return FAILED;
		}


		// link contigs given a contig
		while(linkRound<=SECOND_LINK_ROUND)
		{
			// given a contig with end
			if(linkRound==FIRST_LINK_ROUND)
			{ // the first link round
				if(changeMaxRowCol(maxRowColNode, contigLinkSet, contigGraph, linkRound, NO)==FAILED)
				{
					printf("line=%d, In %s(), cannot change maxRowCol, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// get the column having the maximal value
				if(getMaxColsOfSingleRow(maxRowColNode, contigGraph)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the maximal values for single row %d, error!\n", __LINE__, __func__, maxRowColNode->maxRow);
					return FAILED;
				}

				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || maxRowColNode->secondMaxValue>minLinksNumContigsThres*secondLinkNumFactor)))
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres)))
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || (maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3))))
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || maxRowColNode->secondMaxValue>0.5*averLinkNum))))  // added 2012-11-19
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres)))))  // added 2012-11-22
				//if(maxRowColNode->maxValue==0 || (maxRowColNode->maxValue<minLinksNumContigsThres && maxRowColNode->secondMaxValue>0) || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres)))))  // added 2012-11-22, deleted 2012-11-24
				//if(maxRowColNode->maxValue==0 || maxRowColNode->maxValue>maxLinksNumContigsThres || (maxRowColNode->maxValue<minLinksNumContigsThres && maxRowColNode->secondMaxValue>0) || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*contigGraph->averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres)))))  // added 2012-11-24
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || maxRowColNode->maxValue>maxLinksNumContigsThres || (maxRowColNode->maxValue<minLinksNumContigsThres && maxRowColNode->secondMaxValue>0) || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*contigGraph->averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres)))))  // added 2014-01-19
				if(maxRowColNode->maxValue<2*minLinksNumContigsThres || maxRowColNode->maxValue>maxLinksNumContigsThres || (maxRowColNode->maxValue<minLinksNumContigsThres && maxRowColNode->secondMaxValue>0) || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*contigGraph->averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres)) || maxRowColNode->maxValue-maxRowColNode->secondMaxValue<2)))  // added 2014-01-24
				{
					if(changeMaxRowCol(maxRowColNode, contigLinkSet, contigGraph, linkRound, YES)==FAILED)
					{
						printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
						return FAILED;
					}

					linkRound ++;
					continue;
					//break;
				//}else if(contigGraph->contigItemArray[maxRowColNode->contigID2-1].onlyEnd5==NO)
				}else if(contigGraph->contigItemArray[maxRowColNode->contigID2-1].shortFlag==NO)
				{
					if(contigGraph->contigItemArray[maxRowColNode->contigID2-1].used==NO)
					{
						returnCode = isLinkSingleton(maxRowColNode, contigGraph, linkRound, NO);
						if(returnCode==NO)
						{
							if(changeMaxRowCol(maxRowColNode, contigLinkSet, contigGraph, linkRound, YES)==FAILED)
							{
								printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
								return FAILED;
							}

							linkRound ++;
							continue;
						}else if(returnCode==ERROR)
						{
							printf("line=%d, In %s(), cannot check link singleton, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}else
					{
						if(changeMaxRowCol(maxRowColNode, contigLinkSet, contigGraph, linkRound, YES)==FAILED)
						{
							printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
							return FAILED;
						}

						linkRound ++;
						continue;
					}
				}else
				{
					//if(maxRowColNode->secondMaxValue>0.1*averLinkNum && contigGraph->contigItemArray[maxRowColNode->contigID1-1].onlyEnd5==YES)
					if(maxRowColNode->secondMaxValue>0.1*contigGraph->averLinkNum && contigGraph->contigItemArray[maxRowColNode->contigID1-1].shortFlag==YES)
					{
						if(changeMaxRowCol(maxRowColNode, contigLinkSet, contigGraph, linkRound, YES)==FAILED)
						{
							printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
							return FAILED;
						}

						linkRound ++;
						continue;
					}else
					{
						returnCode = isLinkSingleton(maxRowColNode, contigGraph, linkRound, NO);
						if(returnCode==NO)
						{
							if(changeMaxRowCol(maxRowColNode, contigLinkSet, contigGraph, linkRound, YES)==FAILED)
							{
								printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
								return FAILED;
							}

							linkRound ++;
							continue;
						}else if(returnCode==ERROR)
						{
							printf("line=%d, In %s(), cannot check link singleton, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}

//					if(maxRowColNode->maxValue<15)
//					{
//						if(changeMaxRowCol(maxRowColNode, contigLinkSet, contigGraph, linkRound, YES)==FAILED)
//						{
//							printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
//							return FAILED;
//						}
//
//						linkRound ++;
//						continue;
//					}
				}

				// orient the contigs of the first round
				if(addContigToContigLinkSet(&linkStatus, contigLinkSet, contigGraph, maxRowColNode, 1, linkRound)==FAILED)
				{
					printf("line=%d, In %s(), cannot add contigs to contigLinkArr, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(linkStatus==FAILED)
				{
					if(changeMaxRowCol(maxRowColNode, contigLinkSet, contigGraph, linkRound, YES)==FAILED)
					{
						printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
						return FAILED;
					}

					linkRound ++;
					continue;
				}

				// reset the corresponding rows and columns
				if(markContigGraphEdge(contigGraph, maxRowColNode)==FAILED)
				{
					printf("line=%d, In %s(), cannot mark contigGraph node, error!\n", __LINE__, __func__);
					return FAILED;
				}

			}else
			{ // the second link round

				if(changeMaxRowCol(maxRowColNode, contigLinkSet, contigGraph, linkRound, NO)==FAILED)
				{
					printf("line=%d, In %s(), cannot change maxRowCol, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// get the column having the maximal value
				if(getMaxRowsOfSingleCol(maxRowColNode, contigGraph)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the maximal rows for column %u, error!\n", __LINE__, __func__, maxRowColNode->maxCol);
					return FAILED;
				}

				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || maxRowColNode->secondMaxValue>minLinksNumContigsThres*secondLinkNumFactor)))
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres)))
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || (maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3))))
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || maxRowColNode->secondMaxValue>0.5*averLinkNum)))) // added 2012-11-19
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres))))) // added 2012-11-22
				//if(maxRowColNode->maxValue==0 || (maxRowColNode->maxValue<minLinksNumContigsThres && maxRowColNode->secondMaxValue>0) || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres))))) // added 2012-11-22, deleted 2012-11-24
				//if(maxRowColNode->maxValue==0 || maxRowColNode->maxValue>maxLinksNumContigsThres || (maxRowColNode->maxValue<minLinksNumContigsThres && maxRowColNode->secondMaxValue>0) || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*contigGraph->averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres))))) // added 2012-11-24
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || maxRowColNode->maxValue>maxLinksNumContigsThres || (maxRowColNode->maxValue<minLinksNumContigsThres && maxRowColNode->secondMaxValue>0) || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*contigGraph->averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres))))) // added 2014-01-19
				if(maxRowColNode->maxValue<2*minLinksNumContigsThres || maxRowColNode->maxValue>maxLinksNumContigsThres || (maxRowColNode->maxValue<minLinksNumContigsThres && maxRowColNode->secondMaxValue>0) || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*contigGraph->averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres)) || maxRowColNode->maxValue-maxRowColNode->secondMaxValue<2))) // added 2014-01-24
				{
					linkRound ++;
					break;
				//}else if(contigGraph->contigItemArray[maxRowColNode->contigID1-1].onlyEnd5==NO)
				}else if(contigGraph->contigItemArray[maxRowColNode->contigID1-1].shortFlag==NO)
				{
					if(contigGraph->contigItemArray[maxRowColNode->contigID1-1].used==NO)
					{
						returnCode = isLinkSingleton(maxRowColNode, contigGraph, linkRound, NO);
						if(returnCode==NO)
						{
							linkRound ++;
							break;
						}else if(returnCode==ERROR)
						{
							printf("line=%d, In %s(), cannot check link singleton, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}else
					{
						linkRound ++;
						break;
					}
				}else
				{
					//if(maxRowColNode->secondMaxValue>0.1*contigGraph->averLinkNum && contigGraph->contigItemArray[maxRowColNode->contigID2-1].onlyEnd5==YES)
					if(maxRowColNode->secondMaxValue>0.1*contigGraph->averLinkNum && contigGraph->contigItemArray[maxRowColNode->contigID2-1].shortFlag==YES)
					{
						linkRound ++;
						break;
					}else
					{
						returnCode = isLinkSingleton(maxRowColNode, contigGraph, linkRound, NO);
						if(returnCode==NO)
						{
							linkRound ++;
							break;
						}else if(returnCode==ERROR)
						{
							printf("line=%d, In %s(), cannot check link singleton, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}


//					if(maxRowColNode->maxValue<15)
//					{
//						if(changeMaxRowCol(maxRowColNode, contigLinkSet, contigGraph, linkRound, YES)==FAILED)
//						{
//							printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
//							return FAILED;
//						}
//
//						linkRound ++;
//						break;
//					}
				}

				// orient the contigs of the second round
				if(addContigToContigLinkSet(&linkStatus, contigLinkSet, contigGraph, maxRowColNode, 1, linkRound)==FAILED)
				{
					printf("line=%d, In %s(), cannot add contigs to contigLinkArr, error!\n", __LINE__, __func__);
					return FAILED;
				}
				if(linkStatus==FAILED)
				{
					linkRound ++;
					break;
				}

				// reset the corresponding rows and columns
				if(markContigGraphEdge(contigGraph, maxRowColNode)==FAILED)
				{
					printf("line=%d, In %s(), cannot mark contigGraph node, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		} // while(1)

		// save the linked contigs
		if(saveLinkResultToScaffoldSet(scaffoldSet, contigLinkSet, linkID)==FAILED)
		{
			printf("line=%d, In %s(), cannot save linked contigs to scaffolds, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// save the linked contigs
		//saveLinkResultToFile(fpLinkResult, contigLinkSet, linkID);

		linkID ++;
	}

	if(saveUnlinkedContigsToScaffoldSet(scaffoldSet, contigGraph, &linkID)==FAILED)
	{
		printf("line=%d, In %s(), cannot save unlinked contigs to scaffolds, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//saveUnlinkedContigsToFile(fpLinkResult, contigGraph, &linkID);

	free(maxRowColNode);
	freeMemContigLinkArr(&contigLinkSet);

	printf("The scaffolds number is: %d\n", linkID-1);

	return SUCCESSFUL;
}

/**
 * Initialize contig link array and its memory.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initContigLinkSet(contigLink_t **contigLinkSet, int32_t contigsNum)
{
	*contigLinkSet = (contigLink_t *) calloc(1, sizeof(contigLink_t));
	if((*contigLinkSet)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	(*contigLinkSet)->itemNumContigLinkItemArray = 0;
	(*contigLinkSet)->maxItemNumContigLinkItemArray = contigsNum;
	(*contigLinkSet)->headRowContigLinkItemArray = 0;
	(*contigLinkSet)->tailRowContigLinkItemArray = 0;

	(*contigLinkSet)->contigLinkItemArray = (contigLinkItem_t *) calloc(contigsNum, sizeof(contigLinkItem_t));
	if((*contigLinkSet)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Release memory contig link array.
 */
void freeMemContigLinkArr(contigLink_t **contigLinkSet)
{
	(*contigLinkSet)->itemNumContigLinkItemArray = 0;
	(*contigLinkSet)->headRowContigLinkItemArray = 0;
	(*contigLinkSet)->tailRowContigLinkItemArray = 0;

	free((*contigLinkSet)->contigLinkItemArray);
	(*contigLinkSet)->contigLinkItemArray = NULL;
	free(*contigLinkSet);
	*contigLinkSet = NULL;
}

/**
 * Get first two linked contigs given the first contigID.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getFirstLinkedContigs(int32_t *firstContigID, maxRowCol_t *pMaxRowColNode, contigGraph_t *contigGraph)
{
	int32_t i, satisfiedFlag, returnCode;

	// get first contigID
	satisfiedFlag = NO;
	while((*firstContigID) <= contigGraph->itemNumContigItemArray)
	{
		//if(contigInfoArray[(*firstContigID)-1].used==NO && contigInfoArray[(*firstContigID)-1].onlyEnd5==NO)
		//if(contigInfoArray[(*firstContigID)-1].used==NO && contigInfoArray[(*firstContigID)-1].shortFlag==NO)
		if(contigGraph->contigItemArray[(*firstContigID)-1].used==NO && contigGraph->contigItemArray[(*firstContigID)-1].shortFlag==NO && contigGraph->contigItemArray[(*firstContigID)-1].onlyEnd5==NO)
		{ // the first contig is unused and not too short

			pMaxRowColNode->maxRow = (*firstContigID) * 2 - 1;
			pMaxRowColNode->contigID1 = pMaxRowColNode->maxRow / 2 + 1;
			pMaxRowColNode->endFlag1 = pMaxRowColNode->maxRow % 2;
			if(contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].onlyEnd5==YES)
				pMaxRowColNode->endFlag1 = 2;

			for(i=0; i<2; i++)
			{
				if(getMaxColsOfSingleRow(pMaxRowColNode, contigGraph)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the maximal values for single row %d, error!\n", __LINE__, __func__, pMaxRowColNode->maxRow);
					return FAILED;
				}

				//if(maxValue>=minLinksNumContigsThres || (secondMaxValue>0 && (secondMaxValue/maxValue<maxRatioSecondFirstLinkNum || secondMaxValue<minLinksNumContigsThres*secondLinkNumFactor)))
				//if(maxValue>=minLinksNumContigsThres && (secondMaxValue/maxValue<maxRatioSecondFirstLinkNum && secondMaxValue<minLinksNumContigsThres*secondLinkNumFactor))
				//if(pMaxRowColNode->maxValue>=minLinksNumContigsThres && ((double)pMaxRowColNode->secondMaxValue/pMaxRowColNode->maxValue<maxRatioSecondFirstLinkNum && pMaxRowColNode->secondMaxValue<minLinksNumContigsThres*secondLinkNumFactor))
				//if(pMaxRowColNode->maxValue>=minLinksNumContigsThres && ((double)pMaxRowColNode->secondMaxValue/pMaxRowColNode->maxValue<maxRatioSecondFirstLinkNum && pMaxRowColNode->secondMaxValue<maxSecondLinkNumThres))
				//if(pMaxRowColNode->maxValue>=minLinksNumContigsThres && ((double)pMaxRowColNode->secondMaxValue/pMaxRowColNode->maxValue<maxRatioSecondFirstLinkNum && pMaxRowColNode->secondMaxValue==0))  // deleted 2012-11-24
				//if(pMaxRowColNode->maxValue>=minLinksNumContigsThres && pMaxRowColNode->maxValue<=maxLinksNumContigsThres && ((double)pMaxRowColNode->secondMaxValue/pMaxRowColNode->maxValue<maxRatioSecondFirstLinkNum && pMaxRowColNode->secondMaxValue==0))  // added 2012-11-24
				if(pMaxRowColNode->maxValue>=2*minLinksNumContigsThres && pMaxRowColNode->maxValue<=maxLinksNumContigsThres && ((double)pMaxRowColNode->secondMaxValue/pMaxRowColNode->maxValue<maxRatioSecondFirstLinkNum && pMaxRowColNode->secondMaxValue==0))  // added 2014-01-24
				{
					//if(contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].onlyEnd5==NO)
					if(contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].shortFlag==NO)
					{
						returnCode = isLinkSingleton(pMaxRowColNode, contigGraph, FIRST_LINK_ROUND, YES);
						if(returnCode==YES)
						{
							if((contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigLen<5*contigAlignRegSize || contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigLen<5*contigAlignRegSize) && pMaxRowColNode->maxValue>2*contigGraph->averLinkNum)
								satisfiedFlag = NO;
							else
								satisfiedFlag = YES;
						}else if(returnCode==NO)
							satisfiedFlag = NO;
						else
						{
							printf("line=%d, In %s(), cannot check link singleton, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}else
					{
						//satisfiedFlag = YES;
						satisfiedFlag = NO;
					}
				}

				if(satisfiedFlag==YES)
					break;

				pMaxRowColNode->maxRow --;
			}
		}

		(*firstContigID) ++;

		if(satisfiedFlag==YES)
			break;
	}

	if(satisfiedFlag==NO)
	{
		pMaxRowColNode->contigID1 = pMaxRowColNode->contigID2 = -1;
		pMaxRowColNode->endFlag1 = pMaxRowColNode->endFlag2 = -1;
	}

	return SUCCESSFUL;
}

/**
 * Get the maximal columns and its value of a single row.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxColsOfSingleRow(maxRowCol_t *pMaxRowColNode, contigGraph_t *contigGraph)
{
	int32_t i, j, edgeNum, situationNum, contigID, shortFlag;
	contigEdge_t *pEdgeArray;

	pMaxRowColNode->maxValue = 0;
	pMaxRowColNode->secondMaxValue = 0;
	pMaxRowColNode->maxCol = -1;
	pMaxRowColNode->secondMaxCol = -1;

	if(pMaxRowColNode->endFlag1==1)
	{
		pEdgeArray = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigEdgeArrayEnd3;
		edgeNum = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].itemNumContigEdgeArrayEnd3;
	}else
	{
		pEdgeArray = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigEdgeArrayEnd5;
		edgeNum = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].itemNumContigEdgeArrayEnd5;
	}

	for(i=0; i<edgeNum; i++)
	{
		contigID = pEdgeArray[i].col / 2 + 1;
		shortFlag = contigGraph->contigItemArray[contigID-1].shortFlag;

		if(pEdgeArray[i].used==NO)
		{
			situationNum = pEdgeArray[i].totalSituationNum;
			for(j=0; j<situationNum; j++)
			{
				if(shortFlag==NO)
				{
					if(pEdgeArray[i].maxIndexes[j]>=0 && pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ]> pMaxRowColNode->maxValue)
					{
						pMaxRowColNode->secondMaxValue = pMaxRowColNode->maxValue;
						pMaxRowColNode->secondMaxArrIndex = pMaxRowColNode->maxArrIndex;
						pMaxRowColNode->secondMaxCol = pMaxRowColNode->maxCol;
						pMaxRowColNode->maxValue = pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ];
						pMaxRowColNode->maxArrIndex = pEdgeArray[i].maxIndexes[j];
						pMaxRowColNode->maxCol = pEdgeArray[i].col;
					}else if(pEdgeArray[i].maxIndexes[j]>=0 && pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ]> pMaxRowColNode->secondMaxValue)
					{
						pMaxRowColNode->secondMaxValue = pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ];
						pMaxRowColNode->secondMaxArrIndex = pEdgeArray[i].maxIndexes[j];
						pMaxRowColNode->secondMaxCol = pEdgeArray[i].col;
					}
				}

				// #################### Debug information #####################
#if DEBUG_SCAF_FLAG
				if(pEdgeArray[i].maxIndexes[j]>=0 && pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ] > 0)
				{
					printf("line=%d, In %s(), contigs(%d, %d), row=%d, col=%d, value=%d, situationID=%d, shortFlag=%d\n", __LINE__, __func__, pMaxRowColNode->maxRow/2+1, pEdgeArray[i].col/2+1, pMaxRowColNode->maxRow, pEdgeArray[i].col, pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ], pEdgeArray[i].maxIndexes[j], shortFlag);
				}
#endif
				// #################### Debug information #####################
			}
		}
	}

	for(i=0; i<edgeNum; i++)  // 2014-01-24
	{
		contigID = pEdgeArray[i].col / 2 + 1;
		shortFlag = contigGraph->contigItemArray[contigID-1].shortFlag;

		if(pEdgeArray[i].used==NO)
		{
			situationNum = pEdgeArray[i].totalSituationNum;
			for(j=0; j<situationNum; j++)
			{
				if(shortFlag==YES && (pEdgeArray[i].maxIndexes[j]>=0 && pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ]> 0.3*pMaxRowColNode->maxValue))
				{
					pMaxRowColNode->maxValue = 0;
					pMaxRowColNode->maxArrIndex = -1;
					pMaxRowColNode->maxCol = -1;
					pMaxRowColNode->secondMaxValue = 0;
					pMaxRowColNode->secondMaxArrIndex = -1;
					pMaxRowColNode->secondMaxCol = -1;
				}
			}
		}
	}


	if(pMaxRowColNode->maxCol>=0)
	{
		pMaxRowColNode->contigID2 = pMaxRowColNode->maxCol / 2 + 1;
		pMaxRowColNode->endFlag2 = pMaxRowColNode->maxCol % 2;
		if(contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].onlyEnd5==YES)
			pMaxRowColNode->endFlag2 = 2;

		// #################### Debug information #####################
#if DEBUG_OUT_FLAG
		printf("line=%d, In %s(), contigID1=%d, contigID2=%d, endFlag1=%d, endFlag2=%d, maxValue=%d, secondMaxValue=%d\n", __LINE__, __func__, pMaxRowColNode->contigID1, pMaxRowColNode->contigID2, pMaxRowColNode->endFlag1, pMaxRowColNode->endFlag2, pMaxRowColNode->maxValue, pMaxRowColNode->secondMaxValue);
#endif
		// #################### Debug information #####################
	}else
	{
		pMaxRowColNode->contigID2 = -1;
		pMaxRowColNode->endFlag2 = -1;
	}

	return SUCCESSFUL;
}

/**
 * Get the maximal columns and its value of a single row.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxRowsOfSingleCol(maxRowCol_t *pMaxRowColNode, contigGraph_t *contigGraph)
{
	int i, j, edgeNum, situationNum, newSituationID, contigID, shortFlag;
	contigEdge_t *pEdgeArray;

	pMaxRowColNode->maxValue = 0;
	pMaxRowColNode->secondMaxValue = 0;
	pMaxRowColNode->maxRow = -1;
	pMaxRowColNode->secondMaxRow = -1;

	if(pMaxRowColNode->endFlag2==1)
	{
		pEdgeArray = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigEdgeArrayEnd3;
		edgeNum = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].itemNumContigEdgeArrayEnd3;
	}else
	{
		pEdgeArray = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigEdgeArrayEnd5;
		edgeNum = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].itemNumContigEdgeArrayEnd5;
	}

	for(i=0; i<edgeNum; i++)
	{
		contigID = pEdgeArray[i].col / 2 + 1;
		shortFlag = contigGraph->contigItemArray[contigID-1].shortFlag;

		if(pEdgeArray[i].used==NO)
		{
			situationNum = pEdgeArray[i].totalSituationNum;
			for(j=0; j<situationNum; j++)
			{
				switch(pEdgeArray[i].maxIndexes[j])
				{
					case 0: newSituationID = 2; break;
					case 1: newSituationID = 1; break;
					case 2: newSituationID = 0; break;
					case 3: newSituationID = 3; break;
					default: printf("line=%d, In %s(), maxIndexes[%d]=%d, error!\n", __LINE__, __func__, j, pEdgeArray[i].maxIndexes[j]); return FAILED;
				}

				if(shortFlag==NO)
				{
					if(pEdgeArray[i].maxIndexes[j]>=0 && pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ]> pMaxRowColNode->maxValue)
					{
						pMaxRowColNode->secondMaxValue = pMaxRowColNode->maxValue;
						pMaxRowColNode->secondMaxArrIndex = pMaxRowColNode->maxArrIndex;
						pMaxRowColNode->secondMaxRow = pMaxRowColNode->maxRow;
						pMaxRowColNode->maxValue = pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ];
						pMaxRowColNode->maxArrIndex = newSituationID;
						pMaxRowColNode->maxRow = pEdgeArray[i].col;
					}else if(pEdgeArray[i].maxIndexes[j]>=0 && pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ]> pMaxRowColNode->secondMaxValue)
					{
						pMaxRowColNode->secondMaxValue = pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ];
						pMaxRowColNode->secondMaxArrIndex = newSituationID;
						pMaxRowColNode->secondMaxRow = pEdgeArray[i].col;
					}
				}

				// #################### Debug information #####################
#if DEBUG_SCAF_FLAG
				if(pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ] > 0)
				{
					printf("line=%d, In %s(), contigs(%d, %d), row=%d, col=%d, value=%d, situationID=%d, shortFlag=%d\n", __LINE__, __func__, pEdgeArray[i].col/2+1, pMaxRowColNode->maxCol/2+1, pEdgeArray[i].col, pMaxRowColNode->maxCol, pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ], newSituationID, shortFlag);
				}
#endif
				// #################### Debug information #####################

			}
		}
	}


	for(i=0; i<edgeNum; i++)  // 2014-01-24
	{
		contigID = pEdgeArray[i].col / 2 + 1;
		shortFlag = contigGraph->contigItemArray[contigID-1].shortFlag;

		if(pEdgeArray[i].used==NO)
		{
			situationNum = pEdgeArray[i].totalSituationNum;
			for(j=0; j<situationNum; j++)
			{
				if(shortFlag==YES && (pEdgeArray[i].maxIndexes[j]>=0 && pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ]>= pMaxRowColNode->maxValue))
				{
					pMaxRowColNode->maxValue = 0;
					pMaxRowColNode->maxArrIndex = -1;
					pMaxRowColNode->maxCol = -1;
					pMaxRowColNode->secondMaxValue = 0;
					pMaxRowColNode->secondMaxArrIndex = -1;
					pMaxRowColNode->secondMaxCol = -1;
				}
			}
		}
	}


	if(pMaxRowColNode->maxRow>=0)
	{
		pMaxRowColNode->contigID1 = pMaxRowColNode->maxRow / 2 + 1;
		pMaxRowColNode->endFlag1 = pMaxRowColNode->maxRow % 2;
		if(contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].onlyEnd5==YES)
			pMaxRowColNode->endFlag1 = 2;

		// #################### Debug information #####################
#if DEBUG_OUT_FLAG
		printf("line=%d, In %s(), contigID1=%d, contigID2=%d, endFlag1=%d, endFlag2=%d, maxValue=%d, secondMaxValue=%d\n", __LINE__, __func__, pMaxRowColNode->contigID1, pMaxRowColNode->contigID2, pMaxRowColNode->endFlag1, pMaxRowColNode->endFlag2, pMaxRowColNode->maxValue, pMaxRowColNode->secondMaxValue);
#endif
		// #################### Debug information #####################
	}else
	{
		pMaxRowColNode->contigID1 = -1;
		pMaxRowColNode->endFlag1 = -1;
	}

	return SUCCESSFUL;
}

/**
 * Change the maximal row and columns.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short changeMaxRowCol(maxRowCol_t *pMaxRowColNode, contigLink_t *contigLinkSet, contigGraph_t *contigGraph, int32_t linkRound, int32_t turnRoundFlag)
{
	if(turnRoundFlag==NO)
	{
		if(linkRound==1)
		{
			if(pMaxRowColNode->endFlag2!=2)
			{
				if(pMaxRowColNode->maxCol%2==0)
				{ // 5' end -> 3 'end
					pMaxRowColNode->maxRow = pMaxRowColNode->maxCol + 1;
					pMaxRowColNode->endFlag1 = 1;
				}else
				{ // 3' end -> 5 'end
					pMaxRowColNode->maxRow = pMaxRowColNode->maxCol - 1;
					pMaxRowColNode->endFlag1 = 0;
				}
				pMaxRowColNode->contigID1 = pMaxRowColNode->contigID2;

			}else
			{
				pMaxRowColNode->maxRow = pMaxRowColNode->maxCol;
				pMaxRowColNode->endFlag1 = 2;
				pMaxRowColNode->contigID1 = pMaxRowColNode->contigID2;
			}
		}else
		{
			if(pMaxRowColNode->endFlag1!=2)
			{
				if(pMaxRowColNode->maxRow%2==0)
				{ // 5' end -> 3 'end
					pMaxRowColNode->maxCol = pMaxRowColNode->maxRow + 1;
					pMaxRowColNode->endFlag2 = 1;
				}else
				{
					pMaxRowColNode->maxCol = pMaxRowColNode->maxRow - 1;
					pMaxRowColNode->endFlag2 = 0;
				}
				pMaxRowColNode->contigID2 = pMaxRowColNode->contigID1;

			}else
			{
				pMaxRowColNode->maxCol = pMaxRowColNode->maxRow;
				pMaxRowColNode->endFlag2 = 2;
				pMaxRowColNode->contigID2 = pMaxRowColNode->contigID1;
			}
		}
	}else
	{ // when turn the first round to the second round
		pMaxRowColNode->contigID1 = contigLinkSet->contigLinkItemArray[contigLinkSet->headRowContigLinkItemArray].contigID;
		if(contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].onlyEnd5==YES)
		{
			printf("line=%d, In %s(), the first contig in the linked list should not be short contigs, error!\n", __LINE__, __func__);
			return FAILED;

			//pMaxRowColNode->maxRow = (pMaxRowColNode->contigID1 - 1) * 2;
			//pMaxRowColNode->endFlag1 = 2;
		}else
		{
			if(contigLinkSet->contigLinkItemArray[contigLinkSet->headRowContigLinkItemArray].orientation==ORIENTATION_PLUS)
			{ // plus orientation
				pMaxRowColNode->maxRow = (pMaxRowColNode->contigID1 - 1) * 2 + 1;
				pMaxRowColNode->endFlag1 = 1;
			}else
			{ // minus orienataion
				pMaxRowColNode->maxRow = (pMaxRowColNode->contigID1 - 1) * 2;
				pMaxRowColNode->endFlag1 = 0;
			}
		}
	}

	return SUCCESSFUL;
}


/**
 * Check whether the link is singleton or not.
 *  @return:
 *   (1) If singleton, return YES; otherwise, return NO.
 *   (2) If errors occurred, return ERROR.
 */
short isLinkSingleton(maxRowCol_t *pMaxRowColNode, contigGraph_t *contigGraph, int32_t linkRound, int32_t strictFlag)
{
	int32_t j, maxRow, maxCol, edgeNum1, edgeNum2, singletonFlag, contigID, shortFlag;
	contigEdge_t *pEdgeArray1, *pEdgeArray2;

	singletonFlag = YES;

	if(linkRound==1)
	{ // the first link round

		if(pMaxRowColNode->endFlag2==1)
		{
			pEdgeArray2 = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigEdgeArrayEnd3;
			edgeNum2 = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].itemNumContigEdgeArrayEnd3;
		}else
		{
			pEdgeArray2 = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigEdgeArrayEnd5;
			edgeNum2 = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].itemNumContigEdgeArrayEnd5;
		}

		if(edgeNum2==1)
		{
			singletonFlag = YES;
		}else if(edgeNum2>1)
		{
			for(j=0; j<edgeNum2; j++)
			{

#if (DEBUG_SCAF_FLAG==YES)
				printf("line=%d, In %s(), linkRound=%d, contigs(%d, %d), row=%d, col=%d, value=%d, situationID=%d, shortFlag=%d\n", __LINE__, __func__, linkRound, pMaxRowColNode->maxCol/2+1, pEdgeArray2[j].col/2+1, pMaxRowColNode->maxCol, pEdgeArray2[j].col, pEdgeArray2[j].validNums[ pEdgeArray2[j].maxIndexes[0] ], pEdgeArray2[j].maxIndexes[0], shortFlag);
#endif

				contigID = pEdgeArray2[j].col / 2 + 1;
				shortFlag = contigGraph->contigItemArray[contigID-1].shortFlag;

				//if(singletonFlag==YES && pEdgeArray2[j].col != pMaxRowColNode->maxRow && (pEdgeArray2[j].validNums[pEdgeArray2[j].maxIndexes[0]] > maxSecondLinkNumThres || pEdgeArray2[j].validNums[pEdgeArray2[j].maxIndexes[0]]>0.2*pMaxRowColNode->maxValue))  // added 2014-01-13, deleted 2014-01-16
				//if(singletonFlag==YES && shortFlag==NO && pEdgeArray2[j].col != pMaxRowColNode->maxRow && (pEdgeArray2[j].validNums[pEdgeArray2[j].maxIndexes[0]] > maxSecondLinkNumThres || pEdgeArray2[j].validNums[pEdgeArray2[j].maxIndexes[0]]>0.2*pMaxRowColNode->maxValue))  // added 2014-01-16
				if(singletonFlag==YES && shortFlag==NO && pEdgeArray2[j].col != pMaxRowColNode->maxRow && ((pEdgeArray2[j].validNums[pEdgeArray2[j].maxIndexes[0]] > maxSecondLinkNumThres && pEdgeArray2[j].validNums[pEdgeArray2[j].maxIndexes[0]]>0.1*pMaxRowColNode->maxValue) || pEdgeArray2[j].validNums[pEdgeArray2[j].maxIndexes[0]]>0.1*pMaxRowColNode->maxValue))  // added 2014-01-24
				{
					singletonFlag = NO;
				}
				else if(singletonFlag==YES && shortFlag==YES && pEdgeArray2[j].col != pMaxRowColNode->maxRow && pEdgeArray2[j].validNums[pEdgeArray2[j].maxIndexes[0]] >= 0.3*pMaxRowColNode->maxValue)  // 2014-01-24
				{
					singletonFlag = NO;
				}
				else if(strictFlag==YES && singletonFlag==YES && shortFlag==NO && pEdgeArray2[j].col != pMaxRowColNode->maxRow && pEdgeArray2[j].validNums[pEdgeArray2[j].maxIndexes[0]]>=1 && pMaxRowColNode->maxValue<2*contigGraph->averLinkNum)  // 2014-01-24
				{
					singletonFlag = NO;
				}
			}

/*
			if(pMaxRowColNode->endFlag1==1)
			{
				pEdgeArray1 = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigEdgeArrayEnd3;
				edgeNum1 = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].itemNumContigEdgeArrayEnd3;
			}else
			{
				pEdgeArray1 = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigEdgeArrayEnd5;
				edgeNum1 = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].itemNumContigEdgeArrayEnd5;
			}

			for(j=0; j<edgeNum1; j++)
			{

#if (DEBUG_SCAF_FLAG==YES)
				//printf("num=%d\n", pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]]);
				printf("line=%d, In %s(), contigs(%d, %d), row=%d, col=%d, value=%d, situationID=%d\n", __LINE__, __func__, pMaxRowColNode->maxCol/2+1, pEdgeArray1[j].col/2+1, pMaxRowColNode->maxCol, pEdgeArray1[j].col, pEdgeArray1[j].validNums[ pEdgeArray1[j].maxIndexes[0] ], pEdgeArray1[j].maxIndexes[0]);
#endif

				//if(pEdgeArray[j].col != maxRow && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] >= minLinksNumContigsThres*secondLinkNumFactor)
				//if(pEdgeArray[j].col != maxRow && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] >= maxSecondLinkNumThres)
				//if(pEdgeArray[j].col != maxRow && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] >= minLinksNumContigsThres)
				//if(singletonFlag==YES && pEdgeArray[j].col != maxRow && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] > 0)  // added 2012-11-21
				if(singletonFlag==YES && pEdgeArray1[j].col != pMaxRowColNode->maxRow && pEdgeArray1[j].validNums[pEdgeArray1[j].maxIndexes[0]] > maxSecondLinkNumThres)  // added 2012-11-21
				{
					singletonFlag = NO;
				}
			}
*/
		}else
		{
			printf("line=%d, In %s(), linkRound=%d, maxCol=%d, edgeNum2=%d, error!\n", __LINE__, __func__, linkRound, pMaxRowColNode->maxCol, edgeNum2);
			return ERROR;
		}
	}else
	{ // the second link round

		if(pMaxRowColNode->endFlag1==1)
		{
			pEdgeArray1 = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigEdgeArrayEnd3;
			edgeNum1 = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].itemNumContigEdgeArrayEnd3;
		}else
		{
			pEdgeArray1 = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigEdgeArrayEnd5;
			edgeNum1 = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].itemNumContigEdgeArrayEnd5;
		}

		if(edgeNum1==1)
		{
			singletonFlag = YES;
		}else if(edgeNum1>1)
		{
			for(j=0; j<edgeNum1; j++)
			{

#if (DEBUG_SCAF_FLAG==YES)
				printf("line=%d, In %s(), linkRound=%d, contigs(%d, %d), row=%d, col=%d, value=%d, situationID=%d, shortFlag=%d\n", __LINE__, __func__, linkRound, pMaxRowColNode->maxRow/2+1, pEdgeArray1[j].col/2+1, pMaxRowColNode->maxRow, pEdgeArray1[j].col, pEdgeArray1[j].validNums[ pEdgeArray1[j].maxIndexes[0] ], pEdgeArray1[j].maxIndexes[0], shortFlag);
#endif

				contigID = pEdgeArray1[j].col / 2 + 1;
				shortFlag = contigGraph->contigItemArray[contigID-1].shortFlag;

				//if(singletonFlag==YES && pEdgeArray1[j].col != pMaxRowColNode->maxCol && (pEdgeArray1[j].validNums[pEdgeArray1[j].maxIndexes[0]] > maxSecondLinkNumThres || pEdgeArray1[j].validNums[pEdgeArray1[j].maxIndexes[0]]>0.2*pMaxRowColNode->maxValue))  // added 2014-01-13, deleted 2014-01-16
				//if(singletonFlag==YES && shortFlag==NO && pEdgeArray1[j].col != pMaxRowColNode->maxCol && (pEdgeArray1[j].validNums[pEdgeArray1[j].maxIndexes[0]] > maxSecondLinkNumThres || pEdgeArray1[j].validNums[pEdgeArray1[j].maxIndexes[0]]>0.2*pMaxRowColNode->maxValue))  // added 2014-01-16
				if(singletonFlag==YES && shortFlag==NO && pEdgeArray1[j].col != pMaxRowColNode->maxCol && ((pEdgeArray1[j].validNums[pEdgeArray1[j].maxIndexes[0]] > maxSecondLinkNumThres && pEdgeArray1[j].validNums[pEdgeArray1[j].maxIndexes[0]]>0.1*pMaxRowColNode->maxValue) || pEdgeArray1[j].validNums[pEdgeArray1[j].maxIndexes[0]]>0.1*pMaxRowColNode->maxValue))  // added 2014-01-24
				{
					singletonFlag = NO;
				}
				else if(singletonFlag==YES && shortFlag==YES && pEdgeArray1[j].col != pMaxRowColNode->maxCol && pEdgeArray1[j].validNums[pEdgeArray1[j].maxIndexes[0]] >= 0.3*pMaxRowColNode->maxValue)  // 2014-01-24
				{
					singletonFlag = NO;
				}
				else if(strictFlag==YES && singletonFlag==YES && shortFlag==NO && pEdgeArray1[j].col != pMaxRowColNode->maxCol && pEdgeArray1[j].validNums[pEdgeArray1[j].maxIndexes[0]]>=1 && pMaxRowColNode->maxValue<2*contigGraph->averLinkNum)  // 2014-01-24
				{
					singletonFlag = NO;
				}
			}

/*
			if(pMaxRowColNode->endFlag2==1)
			{
				pEdgeArray2 = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigEdgeArrayEnd3;
				edgeNum2 = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].itemNumContigEdgeArrayEnd3;
			}else
			{
				pEdgeArray2 = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigEdgeArrayEnd5;
				edgeNum2 = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].itemNumContigEdgeArrayEnd5;
			}

			for(j=0; j<edgeNum2; j++)
			{

#if (DEBUG_SCAF_FLAG==YES)
				//printf("num=%d\n", pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]]);
				printf("line=%d, In %s(), contigs(%d, %d), row=%d, col=%d, value=%d, situationID=%d\n", __LINE__, __func__, pMaxRowColNode->maxRow/2+1, pEdgeArray2[j].col/2+1, pMaxRowColNode->maxRow, pEdgeArray2[j].col, pEdgeArray2[j].validNums[ pEdgeArray2[j].maxIndexes[0] ], pEdgeArray2[j].maxIndexes[0]);
#endif

				//if(pEdgeArray[j].col != maxCol && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] >= minLinksNumContigsThres*secondLinkNumFactor)
				//if(pEdgeArray[j].col != maxCol && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] >= maxSecondLinkNumThres)
				//if(pEdgeArray[j].col != maxCol && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] >= minLinksNumContigsThres)
				//if(singletonFlag==YES && pEdgeArray[j].col != maxCol && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] > 0)  // added 2012-11-21
				if(singletonFlag==YES && pEdgeArray2[j].col != pMaxRowColNode->maxCol && pEdgeArray2[j].validNums[pEdgeArray2[j].maxIndexes[0]] > maxSecondLinkNumThres)  // added 2012-11-21
				{
					singletonFlag = NO;
				}
			}
*/
		}else
		{
			printf("line=%d, In %s(), linkRound=%d, maxRow=%d, edgeNum1=%d, error!\n", __LINE__, __func__, linkRound, pMaxRowColNode->maxRow, edgeNum1);
			return ERROR;
		}
	}

	return singletonFlag;
}

/**
 * Add linked contig to contigLinkSet.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addContigToContigLinkSet(int32_t *linkStatus, contigLink_t *contigLinkSet, contigGraph_t *contigGraph, maxRowCol_t *pMaxRowColNode, int32_t newContigNum, int32_t linkRound)
{
	contigLinkItem_t *linkItemArray;

	linkItemArray = contigLinkSet->contigLinkItemArray;

	if(linkRound==1)
	{ // the first link round
		if(newContigNum==2)
		{ // orient the first two contigs

			contigLinkSet->itemNumContigLinkItemArray = 0;
			contigLinkSet->headRowContigLinkItemArray = -1;
			contigLinkSet->tailRowContigLinkItemArray = -1;

			if(pMaxRowColNode->maxArrIndex==0)
			{
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_PLUS;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigLen;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = -1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = -1;
				contigLinkSet->headRowContigLinkItemArray = contigLinkSet->tailRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->itemNumContigLinkItemArray ++;

				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID2;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_PLUS;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigLen;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = -1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = -1;
				contigLinkSet->tailRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->itemNumContigLinkItemArray ++;

				linkItemArray[contigLinkSet->headRowContigLinkItemArray].next = contigLinkSet->tailRowContigLinkItemArray;
				linkItemArray[contigLinkSet->tailRowContigLinkItemArray].previous = contigLinkSet->headRowContigLinkItemArray;

				contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].used = YES;
				contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].used = YES;

			}else if(pMaxRowColNode->maxArrIndex==1)
			{
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_PLUS;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigLen;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = -1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = -1;
				contigLinkSet->headRowContigLinkItemArray = contigLinkSet->tailRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->itemNumContigLinkItemArray ++;

				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID2;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_MINUS;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigLen;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = -1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = -1;
				contigLinkSet->tailRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->itemNumContigLinkItemArray ++;

				linkItemArray[contigLinkSet->headRowContigLinkItemArray].next = contigLinkSet->tailRowContigLinkItemArray;
				linkItemArray[contigLinkSet->tailRowContigLinkItemArray].previous = contigLinkSet->headRowContigLinkItemArray;

				contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].used = YES;
				contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].used = YES;

			}else if(pMaxRowColNode->maxArrIndex==2)
			{
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_MINUS;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigLen;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = -1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = -1;
				contigLinkSet->headRowContigLinkItemArray = contigLinkSet->tailRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->itemNumContigLinkItemArray ++;

				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID2;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_MINUS;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigLen;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = -1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = -1;
				contigLinkSet->tailRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->itemNumContigLinkItemArray ++;

				linkItemArray[contigLinkSet->headRowContigLinkItemArray].next = contigLinkSet->tailRowContigLinkItemArray;
				linkItemArray[contigLinkSet->tailRowContigLinkItemArray].previous = contigLinkSet->headRowContigLinkItemArray;

				contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].used = YES;
				contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].used = YES;

			}else if(pMaxRowColNode->maxArrIndex==3)
			{
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_MINUS;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigLen;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = -1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = -1;
				contigLinkSet->headRowContigLinkItemArray = contigLinkSet->tailRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->itemNumContigLinkItemArray ++;

				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID2;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_PLUS;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigLen;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = -1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = -1;
				contigLinkSet->tailRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->itemNumContigLinkItemArray ++;

				linkItemArray[contigLinkSet->headRowContigLinkItemArray].next = contigLinkSet->tailRowContigLinkItemArray;
				linkItemArray[contigLinkSet->tailRowContigLinkItemArray].previous = contigLinkSet->headRowContigLinkItemArray;

				contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].used = YES;
				contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].used = YES;

			}else
			{
				printf("line=%d, In %s(), linkRound=%d, situationID=%d\n", __LINE__, __func__, linkRound, pMaxRowColNode->maxArrIndex);
				return FAILED;
			}

			*linkStatus = SUCCESSFUL;

		}else
		{ // link the second contig
			if(pMaxRowColNode->maxArrIndex==0)
			{
				if(linkItemArray[contigLinkSet->tailRowContigLinkItemArray].orientation==ORIENTATION_PLUS)
				{ // the orientation of linked contig is changed
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID2;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_PLUS;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigLen;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = contigLinkSet->tailRowContigLinkItemArray;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = -1;

					linkItemArray[contigLinkSet->tailRowContigLinkItemArray].next = contigLinkSet->itemNumContigLinkItemArray;
					contigLinkSet->tailRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
					contigLinkSet->itemNumContigLinkItemArray ++;

					contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].used = YES;

					*linkStatus = SUCCESSFUL;

				}else
				{
					*linkStatus = FAILED;
				}

			}else if(pMaxRowColNode->maxArrIndex==1)
			{
				if(linkItemArray[contigLinkSet->tailRowContigLinkItemArray].orientation==ORIENTATION_PLUS)
				{ // the orientation of linked contig is changed
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID2;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_MINUS;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigLen;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = contigLinkSet->tailRowContigLinkItemArray;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = -1;

					linkItemArray[contigLinkSet->tailRowContigLinkItemArray].next = contigLinkSet->itemNumContigLinkItemArray;
					contigLinkSet->tailRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
					contigLinkSet->itemNumContigLinkItemArray ++;

					contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].used = YES;

					*linkStatus = SUCCESSFUL;

				}else
				{
					*linkStatus = FAILED;
				}

			}else if(pMaxRowColNode->maxArrIndex==2)
			{
				if(linkItemArray[contigLinkSet->tailRowContigLinkItemArray].orientation==ORIENTATION_MINUS)
				{ // the orientation of linked contig is changed
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID2;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_MINUS;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigLen;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = contigLinkSet->tailRowContigLinkItemArray;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = -1;

					linkItemArray[contigLinkSet->tailRowContigLinkItemArray].next = contigLinkSet->itemNumContigLinkItemArray;
					contigLinkSet->tailRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
					contigLinkSet->itemNumContigLinkItemArray ++;

					contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].used = YES;

					*linkStatus = SUCCESSFUL;

				}else
				{
					*linkStatus = FAILED;
				}

			}else if(pMaxRowColNode->maxArrIndex==3)
			{
				if(linkItemArray[contigLinkSet->tailRowContigLinkItemArray].orientation==ORIENTATION_MINUS)
				{ // the orientation of linked contig is changed
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID2;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_PLUS;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigLen;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = contigLinkSet->tailRowContigLinkItemArray;
					linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = -1;

					linkItemArray[contigLinkSet->tailRowContigLinkItemArray].next = contigLinkSet->itemNumContigLinkItemArray;
					contigLinkSet->tailRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
					contigLinkSet->itemNumContigLinkItemArray ++;

					contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].used = YES;

					*linkStatus = SUCCESSFUL;

				}else
				{
					*linkStatus = FAILED;
				}

			}else
			{
				printf("line=%d, In %s(), linkRound=%d, situationID=%d\n", __LINE__, __func__, linkRound, pMaxRowColNode->maxArrIndex);
				return FAILED;
			}
		}
	}else
	{ // the second link round

		if(pMaxRowColNode->maxArrIndex==0)
		{
			if(linkItemArray[contigLinkSet->headRowContigLinkItemArray].orientation==ORIENTATION_PLUS)
			{ // the orientation of linked contig is changed
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_PLUS;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigLen;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = -1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = contigLinkSet->headRowContigLinkItemArray;

				linkItemArray[contigLinkSet->headRowContigLinkItemArray].previous = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->headRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->itemNumContigLinkItemArray ++;

				contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].used = YES;

				*linkStatus = SUCCESSFUL;

			}else
			{
				*linkStatus = FAILED;
			}

		}else if(pMaxRowColNode->maxArrIndex==1)
		{
			if(linkItemArray[contigLinkSet->headRowContigLinkItemArray].orientation==ORIENTATION_MINUS)
			{ // the orientation of linked contig is changed
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_PLUS;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigLen;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = -1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = contigLinkSet->headRowContigLinkItemArray;

				linkItemArray[contigLinkSet->headRowContigLinkItemArray].previous = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->headRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->itemNumContigLinkItemArray ++;

				contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].used = YES;

				*linkStatus = SUCCESSFUL;

			}else
			{
				*linkStatus = FAILED;
			}

		}else if(pMaxRowColNode->maxArrIndex==2)
		{
			if(linkItemArray[contigLinkSet->headRowContigLinkItemArray].orientation==ORIENTATION_MINUS)
			{ // the orientation of linked contig is changed
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_MINUS;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigLen;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = -1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = contigLinkSet->headRowContigLinkItemArray;

				linkItemArray[contigLinkSet->headRowContigLinkItemArray].previous = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->headRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->itemNumContigLinkItemArray ++;

				contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].used = YES;

				*linkStatus = SUCCESSFUL;

			}else
			{
				*linkStatus = FAILED;
			}

		}else if(pMaxRowColNode->maxArrIndex==3)
		{
			if(linkItemArray[contigLinkSet->headRowContigLinkItemArray].orientation==ORIENTATION_PLUS)
			{ // the orientation of linked contig is changed
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigID = pMaxRowColNode->contigID1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].orientation = ORIENTATION_MINUS;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].contigLen = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigLen;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].previous = -1;
				linkItemArray[contigLinkSet->itemNumContigLinkItemArray].next = contigLinkSet->headRowContigLinkItemArray;

				linkItemArray[contigLinkSet->headRowContigLinkItemArray].previous = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->headRowContigLinkItemArray = contigLinkSet->itemNumContigLinkItemArray;
				contigLinkSet->itemNumContigLinkItemArray ++;

				contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].used = YES;

				*linkStatus = SUCCESSFUL;

			}else
			{
				*linkStatus = FAILED;
			}

		}else
		{
			printf("line=%d, In %s(), linkRound=%d, situationID=%d\n", __LINE__, __func__, linkRound, pMaxRowColNode->maxArrIndex);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Mark contigGraph node.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short markContigGraphEdge(contigGraph_t *contigGraph, maxRowCol_t *pMaxRowColNode)
{
	int32_t i, edgeNum;
	contigEdge_t *pEdgeArray;

	// first contig
	if(pMaxRowColNode->endFlag1==1)
	{
		pEdgeArray = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigEdgeArrayEnd3;
		edgeNum = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].itemNumContigEdgeArrayEnd3;
	}else
	{
		pEdgeArray = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].contigEdgeArrayEnd5;
		edgeNum = contigGraph->contigItemArray[pMaxRowColNode->contigID1-1].itemNumContigEdgeArrayEnd5;
	}

	for(i=0; i<edgeNum; i++)
	{
		if(pEdgeArray[i].col==pMaxRowColNode->maxCol && pEdgeArray[i].used==NO)
		{
			pEdgeArray[i].used = YES;
		}
	}

	// second contig
	if(pMaxRowColNode->endFlag2==1)
	{
		pEdgeArray = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigEdgeArrayEnd3;
		edgeNum = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].itemNumContigEdgeArrayEnd3;
	}else
	{
		pEdgeArray = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].contigEdgeArrayEnd5;
		edgeNum = contigGraph->contigItemArray[pMaxRowColNode->contigID2-1].itemNumContigEdgeArrayEnd5;
	}

	for(i=0; i<edgeNum; i++)
	{
		if(pEdgeArray[i].col==pMaxRowColNode->maxRow && pEdgeArray[i].used==NO)
		{
			pEdgeArray[i].used = YES;
		}
	}

	return SUCCESSFUL;
}

/**
 * Save the linked contigs to scaffold set.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short saveLinkResultToScaffoldSet(scaffoldSet_t *scaffoldSet, contigLink_t *contigLinkSet, int32_t linkID)
{
	int32_t tmpRow, itemNum, contigID1, readOrient1, contigLen1, contigID2, readOrient2, contigLen2;
	contigLinkItem_t *linkItemArray;
	scaffoldItem_t *scaffoldItem;

	// allocate scaffold item node
	scaffoldItem = (scaffoldItem_t *) calloc (1, sizeof(scaffoldItem_t));
	if(scaffoldItem==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	scaffoldItem->scaffoldID = linkID;
	scaffoldItem->linkedContigsNum = contigLinkSet->itemNumContigLinkItemArray;
	scaffoldItem->next = NULL;
	if(scaffoldSet->scaffoldItemList==NULL)
	{
		scaffoldSet->scaffoldItemList = scaffoldSet->tailScaffoldItem = scaffoldItem;
		scaffoldItem->previous = NULL;
	}else
	{
		scaffoldSet->tailScaffoldItem->next = scaffoldItem;
		scaffoldItem->previous = scaffoldSet->tailScaffoldItem;
		scaffoldSet->tailScaffoldItem = scaffoldItem;
	}
	scaffoldSet->scaffoldNum ++;

	if(contigLinkSet->itemNumContigLinkItemArray>=2)
		scaffoldItem->itemNumContigOverlapArray = contigLinkSet->itemNumContigLinkItemArray - 1;
	else
		scaffoldItem->itemNumContigOverlapArray = 1;

	scaffoldItem->contigOverlapArray = (contigOverlap_t *) calloc (scaffoldItem->itemNumContigOverlapArray, sizeof(contigOverlap_t));
	if(scaffoldItem->contigOverlapArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}


	linkItemArray = contigLinkSet->contigLinkItemArray;
	tmpRow = contigLinkSet->headRowContigLinkItemArray;

	contigID1 = linkItemArray[tmpRow].contigID;
	readOrient1 = linkItemArray[tmpRow].orientation;

	if(contigLinkSet->itemNumContigLinkItemArray>=2)
	{
		itemNum = 0;
		while(1)
		{
			tmpRow = linkItemArray[tmpRow].next;

			contigID2 = linkItemArray[tmpRow].contigID;
			readOrient2 = linkItemArray[tmpRow].orientation;

			scaffoldItem->contigOverlapArray[itemNum].contigID1 = contigID1;
			scaffoldItem->contigOverlapArray[itemNum].orientation1 = readOrient1;
			scaffoldItem->contigOverlapArray[itemNum].contigID2 = contigID2;
			scaffoldItem->contigOverlapArray[itemNum].orientation2 = readOrient2;
			itemNum ++;

			contigID1 = contigID2;
			readOrient1 = readOrient2;

			if(tmpRow==contigLinkSet->tailRowContigLinkItemArray)
				break;
		}
	}else
	{
		scaffoldItem->contigOverlapArray[0].contigID1 = contigID1;
		scaffoldItem->contigOverlapArray[0].orientation1 = readOrient1;
	}

	return SUCCESSFUL;
}

/**
 * Save the unlinked contigs to scaffold set.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short saveUnlinkedContigsToScaffoldSet(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph, int32_t *linkID)
{
	int32_t i, j;
	scaffoldItem_t *scaffoldItem;

	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
		if(contigGraph->contigItemArray[i].used==NO)
		{
			// allocate scaffold item node
			scaffoldItem = (scaffoldItem_t *) calloc (1, sizeof(scaffoldItem_t));
			if(scaffoldItem==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			scaffoldItem->scaffoldID = (*linkID) ++ ;
			scaffoldItem->linkedContigsNum = 1;
			scaffoldItem->next = NULL;
			if(scaffoldSet->scaffoldItemList==NULL)
			{
				scaffoldSet->scaffoldItemList = scaffoldSet->tailScaffoldItem = scaffoldItem;
				scaffoldItem->previous = NULL;
			}else
			{
				scaffoldSet->tailScaffoldItem->next = scaffoldItem;
				scaffoldItem->previous = scaffoldSet->tailScaffoldItem;
				scaffoldSet->tailScaffoldItem = scaffoldItem;
			}
			scaffoldSet->scaffoldNum ++;

			scaffoldItem->itemNumContigOverlapArray = 1;
			scaffoldItem->contigOverlapArray = (contigOverlap_t *) calloc (scaffoldItem->itemNumContigOverlapArray, sizeof(contigOverlap_t));
			if(scaffoldItem->contigOverlapArray==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			scaffoldItem->contigOverlapArray[0].contigID1 = contigGraph->contigItemArray[i].contigID;
			scaffoldItem->contigOverlapArray[0].orientation1 = ORIENTATION_PLUS;
		}
	}

	return SUCCESSFUL;
}


