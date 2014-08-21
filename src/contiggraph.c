/*
 * contiggraph.c
 *
 *  Created on: Oct 28, 2013
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Initialize the contigGraph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initContigGraph(contigGraph_t **contigGraph)
{
	// allocate the contigGraph node
	*contigGraph = (contigGraph_t*) calloc(1, sizeof(contigGraph_t));
	if((*contigGraph)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	(*contigGraph)->itemNumContigItemArray = 0;
	(*contigGraph)->maxItemNumContigItemArray = MAX_ITEMNUM_CONTIGGRAPH_DEFAULT;

	(*contigGraph)->contigItemArray = (contigGraphItem_t*) calloc((*contigGraph)->maxItemNumContigItemArray, sizeof(contigGraphItem_t));
	if((*contigGraph)->contigItemArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}


/**
 * Release the contigGraph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short releaseContigGraph(contigGraph_t **contigGraph)
{
	int32_t i, contigItemNum;
	contigGraphItem_t *contigItemArray;

	contigItemNum = (*contigGraph)->itemNumContigItemArray;
	contigItemArray = (*contigGraph)->contigItemArray;

	for(i=0; i<contigItemNum; i++)
	{
		if(contigItemArray[i].contigTitle)
			free(contigItemArray[i].contigTitle);
		if(contigItemArray[i].contigSeq)
			free(contigItemArray[i].contigSeq);
		contigItemArray[i].contigTitle = NULL;
		contigItemArray[i].contigSeq = NULL;
		if(contigItemArray[i].contigReadArrayEnd5)
			free(contigItemArray[i].contigReadArrayEnd5);
		if(contigItemArray[i].contigReadArrayEnd3)
			free(contigItemArray[i].contigReadArrayEnd3);
		contigItemArray[i].contigReadNumEnd5 = 0;
		contigItemArray[i].contigReadNumEnd3 = 0;
		if(contigItemArray[i].contigEdgeArrayEnd5)
			free(contigItemArray[i].contigEdgeArrayEnd5);
		if(contigItemArray[i].contigEdgeArrayEnd3)
			free(contigItemArray[i].contigEdgeArrayEnd3);
		contigItemArray[i].itemNumContigEdgeArrayEnd5 = 0;
		contigItemArray[i].itemNumContigEdgeArrayEnd3 = 0;
	}

	free((*contigGraph)->contigItemArray);
	(*contigGraph)->itemNumContigItemArray = 0;
	(*contigGraph)->maxItemNumContigItemArray = 0;
	free(*contigGraph);
	*contigGraph = NULL;

	return SUCCESSFUL;
}


/**
 * Copy the public contigArr to private contigGraph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addContigItemToContigGraph(contigGraph_t *contigGraph, int32_t contigID, int32_t localContigID, contigtype *contigArray, int32_t itemNumContigArray)
{
	int32_t i;
	contigGraphItem_t *newContigItem;
	char baseChar, *contigSeq;

	if(contigGraph->itemNumContigItemArray>=contigGraph->maxItemNumContigItemArray)
	{
		contigGraph->maxItemNumContigItemArray *= 2;
		contigGraph->contigItemArray = realloc(contigGraph->contigItemArray, contigGraph->maxItemNumContigItemArray * sizeof(contigGraphItem_t));
		if(contigGraph->contigItemArray==NULL)
		{
			printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	if(itemNumContigArray>0)
	{
		newContigItem = contigGraph->contigItemArray + contigGraph->itemNumContigItemArray;
		contigGraph->itemNumContigItemArray ++;

		newContigItem->contigID = contigID;
		newContigItem->localContigID = localContigID;
		newContigItem->contigTitle = NULL;
		newContigItem->contigLen = itemNumContigArray;
		newContigItem->validFlag = YES;

		for(i=0; i<itemNumContigArray; i++)
		{
			if(contigArray[i].ridposnum>0)
				newContigItem->totalAlignedReadsNum += contigArray[i].ridposnum;
		}

		contigSeq = (char *) malloc ((itemNumContigArray+1) * sizeof(char));
		if(contigSeq==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// copy the base sequence
		for(i=0; i<itemNumContigArray; i++)
		{
			switch(contigArray[i].base)
			{
				case 0: baseChar = 'A'; break;
				case 1: baseChar = 'C'; break;
				case 2: baseChar = 'G'; break;
				case 3: baseChar = 'T'; break;
			}
			contigSeq[i] = baseChar;
		}
		contigSeq[i] = '\0';
		newContigItem->contigSeq = contigSeq;
	}else
	{
		printf("line=%d, In %s(), itemNumContigArray=%d, error!\n", __LINE__, __func__, itemNumContigArray);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Save the reads match information.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short saveReadsMatchInfo(FILE *fpReadMatchInfo, int32_t contigID, contigtype *contigArray, int64_t itemNumContigArray, int32_t alignRegSize)
{
	int64_t i, j;
	int32_t posNum, outputFlagEnd5, outputFlagEnd3, startRowEnd5, endRowEnd5, startRowEnd3, endRowEnd3, midRow, tmpBaseNum;
	successRead_t *contigReadArray;
	readMatchInfo_t readMatchInfo;

	if(itemNumContigArray>=2*alignRegSize)
	{
		outputFlagEnd5 = YES;
		outputFlagEnd3 = YES;

		startRowEnd5 = 0;
		endRowEnd5 = alignRegSize - 1;
		startRowEnd3 = itemNumContigArray - alignRegSize;
		endRowEnd3 = itemNumContigArray - 1;
	}else if(itemNumContigArray>alignRegSize)
	{
		outputFlagEnd5 = YES;
		outputFlagEnd3 = YES;

		midRow = itemNumContigArray / 2;
		startRowEnd5 = 0;
		endRowEnd5 = midRow - 1;
		startRowEnd3 = midRow;
		endRowEnd3 = itemNumContigArray - 1;
	}else
	{
		outputFlagEnd5 = YES;
		outputFlagEnd3 = NO;

		startRowEnd5 = 0;
		endRowEnd5 = itemNumContigArray - 1;
	}
/*
	if(outputFlagEnd5==YES)
	{
		tmpBaseNum = 2 * readLen;
		if(tmpBaseNum>itemNumContigArray)
			tmpBaseNum = itemNumContigArray;
		for(i=0; i<tmpBaseNum; i++)
		{
			if(contigArray[i].itemNumContigPath>=2)
			{
				outputFlagEnd5 = NO;
				break;
			}
		}
	}

	if(outputFlagEnd3==YES)
	{
		tmpBaseNum = 2 * readLen;
		if(tmpBaseNum>itemNumContigArray)
			tmpBaseNum = itemNumContigArray;
		for(i=itemNumContigArray-tmpBaseNum; i<itemNumContigArray; i++)
		{
			if(contigArray[i].itemNumContigPath>=2)
			{
				outputFlagEnd3 = NO;
				break;
			}
		}
	}
*/
	if(outputFlagEnd5==YES)
	{
		for(i=startRowEnd5; i<=endRowEnd5; i++)
		{
			if(contigArray[i].ridposnum>0)
			{
				readMatchInfo.contigID = contigID;
				readMatchInfo.contigPos = contigArray[i].index;
				contigReadArray = contigArray[i].pridposorientation;
				posNum = contigArray[i].ridposnum;
				for(j=0; j<posNum; j++)
				{
					readMatchInfo.readID = contigReadArray[j].rid;
					readMatchInfo.seqlen = contigReadArray[j].seqlen;
					readMatchInfo.matchlen = contigReadArray[j].matchlen;
					readMatchInfo.readOrientation = contigReadArray[j].orientation;
					readMatchInfo.contigEnd = -1;

					if(fwrite(&readMatchInfo, sizeof(readMatchInfo_t), 1, fpReadMatchInfo)!=1)
					{
						printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
						return FAILED;
					}

					//fprintf(fpReadMatchInfo, "%d\t%d\t%ld\t%d\n", readMatchInfo.contigID, readMatchInfo.contigPos, (int64_t)readMatchInfo.readID, readMatchInfo.readOrientation);
				}
			}
		}
	}

	if(outputFlagEnd3==YES)
	{
		for(i=startRowEnd3; i<=endRowEnd3; i++)
		{
			if(contigArray[i].ridposnum>0)
			{
				readMatchInfo.contigID = contigID;
				readMatchInfo.contigPos = contigArray[i].index;
				contigReadArray = contigArray[i].pridposorientation;
				posNum = contigArray[i].ridposnum;
				for(j=0; j<posNum; j++)
				{
					readMatchInfo.readID = contigReadArray[j].rid;
					readMatchInfo.seqlen = contigReadArray[j].seqlen;
					readMatchInfo.matchlen = contigReadArray[j].matchlen;
					readMatchInfo.readOrientation = contigReadArray[j].orientation;
					readMatchInfo.contigEnd = -1;

					if(fwrite(&readMatchInfo, sizeof(readMatchInfo_t), 1, fpReadMatchInfo)!=1)
					{
						printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
						return FAILED;
					}

					//fprintf(fpReadMatchInfo, "%d\t%d\t%ld\t%d\n", readMatchInfo.contigID, readMatchInfo.contigPos, (int64_t)readMatchInfo.readID, readMatchInfo.readOrientation);
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Copy the public contigArr to private contigGraph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputContigFromContigGraph(char *contigFile, contigGraph_t *contigGraph, int32_t minContigLen)
{
	int32_t i, contigNum;
	FILE *fpContig;

	fpContig = fopen(contigFile, "w");
	if(fpContig==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigFile);
		return FAILED;
	}

	// output the base sequences
	contigNum = contigGraph->itemNumContigItemArray;
	for(i=0; i<contigNum; i++)
	{
		if(contigGraph->contigItemArray[i].contigLen>=minContigLen)
		{
			fprintf(fpContig, ">%d length: %d\n", contigGraph->contigItemArray[i].contigID, contigGraph->contigItemArray[i].contigLen);
			fprintf(fpContig, "%s\n", contigGraph->contigItemArray[i].contigSeq);
		}
	}
	fclose(fpContig);


	// get the statistics of contig lengths
	if(contigsLenStatistics(contigGraph, minContigLen)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute contig length statistics, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}


