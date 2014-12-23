/*
 * merge.c
 *
 *  Created on: Mar 1, 2014
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Merge overlapped contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mergeOverlappedContigs(char *readMatchInfoFile, contigGraph_t *contigGraph, graphtype *graph)
{
	// initialize scaffoldSet
	if(initMemMergeContigs(contigGraph, &scaffoldSet, graph->readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory for scaffolding, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the reads information to readSet
	if(fillReadInfoToReadSet(contigGraph, graph->readSet, readMatchInfoFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the read information to readSet, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill read information at contig ends
	if(fillReadMatchInfoContigEnds(contigGraph, graph->readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill read match information to read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// link contigs
	if(contigsLinking(scaffoldSet, contigGraph, graph->readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot link contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// overlap
	if(overlapContigsInScaf(scaffoldSet, contigGraph, graph->readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot overlap contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ############################## Debug information ##########################
//	if(outputScaffoldSetToFile("../tmpLinkContigs_merge2.txt", scaffoldSet, contigGraph)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot output overlapped contigs, error!\n", __LINE__, __func__);
//		return FAILED;
//	}
	// ############################## Debug information ##########################

	// split gaped contigs
	if(splitUnoverlappedContigs(scaffoldSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot merge contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ############################## Debug information ##########################
//	if(outputScaffoldSetToFile("../tmpOverlappedContigs_merge2.txt", scaffoldSet, contigGraph)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot output overlapped contigs, error!\n", __LINE__, __func__);
//		return FAILED;
//	}
	// ############################## Debug information ##########################

	// merge overlapped contigs
	if(mergeOverContigs(scaffoldSet, contigGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot merge contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// remove merged contigs
	if(removeMergedContigs(contigGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot merge contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}

	freeMemMergeContigs(&scaffoldSet, contigGraph, graph->readSet);

	remove(readMatchInfoFile); // remove the tmp file

	return SUCCESSFUL;
}


/**
 * Initialize the memory for merging overlapped contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemMergeContigs(contigGraph_t *contigGraph, scaffoldSet_t **scaffoldSet, readSet_t *readSet)
{
	if(initMemScaffolding(scaffoldSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory for scaffolding, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize read blocks
	if(initReadMatchInfoBlockInReadset(readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize read blocks, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(initContigGraphMergeContigs(contigGraph, contigAlignRegSize)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize contig graph, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}


/*
 * Free the memory for merging overlapped contigs.
 */
void freeMemMergeContigs(scaffoldSet_t **scaffoldSet, contigGraph_t *contigGraph, readSet_t *readSet)
{
	releaseScaffoldSet(scaffoldSet);
	releaseReadMatchInfoBlockInReadset(readSet);
	releaseReadAndEdgeInfoContigGraph(contigGraph);
}

void releaseReadAndEdgeInfoContigGraph(contigGraph_t *contigGraph)
{
	int32_t i, contigItemNum;
	contigGraphItem_t *contigItemArray;

	contigItemNum = contigGraph->itemNumContigItemArray;
	contigItemArray = contigGraph->contigItemArray;

	for(i=0; i<contigItemNum; i++)
	{
		if(contigItemArray[i].contigReadArrayEnd5)
		{
			free(contigItemArray[i].contigReadArrayEnd5);
			contigItemArray[i].contigReadArrayEnd5 = NULL;
		}
		if(contigItemArray[i].contigReadArrayEnd3)
		{
			free(contigItemArray[i].contigReadArrayEnd3);
			contigItemArray[i].contigReadArrayEnd3 = NULL;
		}
		contigItemArray[i].contigReadNumEnd5 = 0;
		contigItemArray[i].contigReadNumEnd3 = 0;
		if(contigItemArray[i].contigEdgeArrayEnd5)
		{
			free(contigItemArray[i].contigEdgeArrayEnd5);
			contigItemArray[i].contigEdgeArrayEnd5 = NULL;
		}
		if(contigItemArray[i].contigEdgeArrayEnd3)
		{
			free(contigItemArray[i].contigEdgeArrayEnd3);
			contigItemArray[i].contigEdgeArrayEnd3 = NULL;
		}
		contigItemArray[i].itemNumContigEdgeArrayEnd5 = 0;
		contigItemArray[i].itemNumContigEdgeArrayEnd3 = 0;

		contigItemArray[i].alignRegSizeEnd5 = 0;
		contigItemArray[i].alignRegSizeEnd3 = 0;
		contigItemArray[i].totalAlignedReadsNum = 0;
		contigItemArray[i].onlyEnd5 = NO;
		contigItemArray[i].shortFlag = NO;
		contigItemArray[i].used = NO;
		contigItemArray[i].usedTimeEnd5 = 0;
		contigItemArray[i].usedTimeEnd3 = 0;
		contigItemArray[i].maxOverlapLenEnd5 = 0;
		contigItemArray[i].overlapKindEnd5 = 0;
		contigItemArray[i].maxOverlapLenEnd3 = 0;
		contigItemArray[i].overlapKindEnd3 = 0;
		contigItemArray[i].delReadFlagEnd5 = NO;
		contigItemArray[i].delReadFlagEnd3 = NO;
		contigItemArray[i].averCovNumEnd5 = 0;
		contigItemArray[i].averCovNumEnd3 = 0;
	}
}

/**
 * Initialize contig graph for merging contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initContigGraphMergeContigs(contigGraph_t *contigGraph, int32_t contigAlignRegSize)
{
	int32_t i, j, contigLen;
	contigGraphItem_t *contigItem;

	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
		contigItem = contigGraph->contigItemArray + i;
		contigLen = contigItem->contigLen;

		if(contigLen>2*contigAlignRegSize)
		{
			contigItem->alignRegSizeEnd5 = contigAlignRegSize;
			contigItem->alignRegSizeEnd3 = contigAlignRegSize;
		}
		else if(contigLen>contigAlignRegSize)
		{
			contigItem->alignRegSizeEnd5 = (contigLen + 1) / 2;
			contigItem->alignRegSizeEnd3 = contigItem->contigLen - contigItem->alignRegSizeEnd5;
		}else
		{
			contigItem->alignRegSizeEnd5 = contigLen;
			contigItem->alignRegSizeEnd3 = 0;
		}
	}

	return SUCCESSFUL;
}

/**
 * Fill the read information to readSet.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillReadInfoToReadSet(contigGraph_t *contigGraph, readSet_t *readSet, char *readMatchInfoFile)
{
	FILE *fpMatchInfo;
	readMatchInfo_t readMatchInfoTmp;
	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo;
	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;
	contigGraphItem_t *contigItem;


	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;


	fpMatchInfo = fopen(readMatchInfoFile, "r");
	if(fpMatchInfo==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readMatchInfoFile);
		return FAILED;
	}

	while(1)
	{
		if(fread(&readMatchInfoTmp, sizeof(readMatchInfo_t), 1, fpMatchInfo)!=1)
		{
			if(feof(fpMatchInfo))
			{
				break;
			}else
			{
				printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		contigItem = contigGraph->contigItemArray + readMatchInfoTmp.contigID - 1;

		readMatchInfoTmp.contigEnd = -1;
		if(contigItem->alignRegSizeEnd3>0)
		{
			if(readMatchInfoTmp.contigPos>=1 && readMatchInfoTmp.contigPos-1+readMatchInfoTmp.seqlen<=contigItem->alignRegSizeEnd5)
			{
				readMatchInfoTmp.contigEnd = 0;
			}else if(readMatchInfoTmp.contigPos>=contigItem->contigLen-contigItem->alignRegSizeEnd3+1 && readMatchInfoTmp.contigPos-1+readMatchInfoTmp.seqlen<=contigItem->contigLen)
			{
				readMatchInfoTmp.contigEnd = 1;
			}
		}else
		{
			if(readMatchInfoTmp.contigPos>=1 && readMatchInfoTmp.contigPos-1+readMatchInfoTmp.seqlen<=contigItem->alignRegSizeEnd5)
			{
				readMatchInfoTmp.contigEnd = 2;
			}
		}

		readMatchInfoBlockID = (readMatchInfoTmp.readID - 1) / maxItemNumPerReadMatchInfoBlock;
		rowNumInReadMatchInfoBlock = (readMatchInfoTmp.readID - 1) % maxItemNumPerReadMatchInfoBlock;
		pReadMatchInfo = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;

		if(readMatchInfoTmp.contigEnd!=-1)
		{
			// save the match result
			readMatchInfoBlockArr[readMatchInfoBlockID].itemNum ++;
			readSet->totalValidItemNumReadMatchInfo ++;

			pReadMatchInfo->contigID = readMatchInfoTmp.contigID;
			pReadMatchInfo->contigPos = readMatchInfoTmp.contigPos;
			pReadMatchInfo->matchlen = readMatchInfoTmp.matchlen;
			pReadMatchInfo->readID = readMatchInfoTmp.readID;
			pReadMatchInfo->seqlen = readMatchInfoTmp.seqlen;
			pReadMatchInfo->readOrientation = readMatchInfoTmp.readOrientation;
			pReadMatchInfo->contigEnd = readMatchInfoTmp.contigEnd;

			if(pReadMatchInfo->contigEnd==1)
				contigGraph->contigItemArray[pReadMatchInfo->contigID-1].contigReadNumEnd3 ++;
			else
				contigGraph->contigItemArray[pReadMatchInfo->contigID-1].contigReadNumEnd5 ++;
		}else
		{
			pReadMatchInfo->contigID = -1;
			pReadMatchInfo->contigPos = -1;
		}
	}

	fclose(fpMatchInfo);

	return SUCCESSFUL;
}


/**
 * Merge the overlapped contigs by scaffoldSet.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short splitUnoverlappedContigs(scaffoldSet_t *scaffoldSet)
{
	int32_t i, j, itemNumOverlapArray;
	scaffoldItem_t *scafItem;
	contigOverlap_t *overlapArray, *overlapInfo;
	int32_t contigID1, contigID2, contigLen1, contigLen2, contigNum;

	scafItem = scaffoldSet->scaffoldItemList;
	while(scafItem)
	{
		if(scafItem->linkedContigsNum>=2)
		{
			overlapArray = scafItem->contigOverlapArray;
			itemNumOverlapArray = scafItem->itemNumContigOverlapArray;

			for(i=0; i<itemNumOverlapArray; i++)
			{
				if((overlapArray[i].breakFlag==NO && overlapArray[i].mergeFlag==NO) || (overlapArray[i].mergeFlag==YES && overlapArray[i].overlapLen<8))
				{
					overlapArray[i].breakFlag = YES;
				}
			}
		}

		scafItem = scafItem->next;
	}


	// split broken scaffolds
	if(splitScaffolds(scaffoldSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate contig overlap information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Merge the overlapped contigs by scaffoldSet.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mergeOverContigs(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph)
{
	int32_t i, j, itemNumOverlapArray;
	scaffoldItem_t *scafItem;
	contigOverlap_t *overlapArray, *overlapInfo;
	int32_t contigID1, contigID2, contigLen1, contigLen2, totalLenTmp, newContigLen, overlapLen, targetContigID;
	char *newContigseq, *contigseq1, *contigseq2, *contigseqBuf;

	scafItem = scaffoldSet->scaffoldItemList;
	while(scafItem)
	{
		if(scafItem->linkedContigsNum>=2)
		{
			overlapArray = scafItem->contigOverlapArray;
			itemNumOverlapArray = scafItem->itemNumContigOverlapArray;

			// get the total length
			totalLenTmp = 0;
			contigID1 = overlapArray[0].contigID1;
			totalLenTmp = contigLen1 = contigGraph->contigItemArray[contigID1-1].contigLen;

			for(i=0; i<itemNumOverlapArray; i++)
			{
				contigID2 = overlapArray[i].contigID2;
				contigLen2 = contigGraph->contigItemArray[contigID2-1].contigLen;
				totalLenTmp += contigLen2;
			}

			// allocate the memory
			newContigseq = (char*) calloc (totalLenTmp+1, sizeof(char));
			if(newContigseq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			contigseqBuf = (char*) calloc (totalLenTmp+1, sizeof(char));
			if(contigseqBuf==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// generate the new sequence
			contigID1 = overlapArray[0].contigID1;
			contigseq1 = contigGraph->contigItemArray[contigID1-1].contigSeq;
			contigLen1 = contigGraph->contigItemArray[contigID1-1].contigLen;
			newContigLen = 0;
			if(overlapArray[0].orientation1==ORIENTATION_PLUS)
			{
				strcpy(newContigseq, contigseq1);
				newContigLen += contigLen1;
			}
			else
			{
				strcpy(newContigseq, contigseq1);
				if(reverseSeq(newContigseq, contigLen1)==FAILED) // reverse the sequence
				{
					printf("line=%d, In %s(), cannot reverse sequence, error!\n", __LINE__, __func__);
					return FAILED;
				}
				newContigLen += contigLen1;
			}

			for(i=0; i<itemNumOverlapArray; i++)
			{
				contigID2 = overlapArray[i].contigID2;
				contigseq1 = contigGraph->contigItemArray[contigID2-1].contigSeq;
				contigLen2 = contigGraph->contigItemArray[contigID2-1].contigLen;
				overlapLen = overlapArray[i].overlapLen;

				if(overlapArray[i].orientation2==ORIENTATION_PLUS)
				{
					strcpy(newContigseq+newContigLen, contigGraph->contigItemArray[contigID2-1].contigSeq+overlapLen);
					newContigLen += contigLen2 - overlapLen;
				}else
				{
					strcpy(contigseqBuf, contigGraph->contigItemArray[contigID2-1].contigSeq); // copy sequence
					if(reverseSeq(contigseqBuf, contigLen2)==FAILED) // reverse the sequence
					{
						printf("line=%d, In %s(), cannot reverse sequence, error!\n", __LINE__, __func__);
						return FAILED;
					}

					strcpy(newContigseq+newContigLen, contigseqBuf+overlapLen); // copy to newContigseq
					newContigLen += contigLen2 - overlapLen;
				}
			}

			// mark the merged contigs
			targetContigID = overlapArray[0].contigID1;
			for(i=0; i<itemNumOverlapArray; i++)
			{
				contigID2 = overlapArray[i].contigID2;
				if(contigID2<targetContigID)
					targetContigID = contigID2;
			}

			contigID1 = overlapArray[0].contigID1;
			if(contigID1!=targetContigID)
				contigGraph->contigItemArray[contigID1-1].validFlag = NO;
			for(i=0; i<itemNumOverlapArray; i++)
			{
				contigID2 = overlapArray[i].contigID2;
				if(contigID2!=targetContigID)
					contigGraph->contigItemArray[contigID2-1].validFlag = NO;
			}

			free(contigGraph->contigItemArray[targetContigID-1].contigSeq);
			contigGraph->contigItemArray[targetContigID-1].contigSeq = newContigseq;
			contigGraph->contigItemArray[targetContigID-1].contigLen = newContigLen;

			free(contigseqBuf);
			contigseqBuf = NULL;
		}

		scafItem = scafItem->next;
	}

	return SUCCESSFUL;
}

/**
 * Remove merged contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeMergedContigs(contigGraph_t *contigGraph)
{
	int32_t i, j;
	contigGraphItem_t *contigItemArray;

	i = 0;
	contigItemArray = contigGraph->contigItemArray;
	while(i<contigGraph->itemNumContigItemArray)
	{
		contigItemArray[i].contigID = i + 1;

		if(contigItemArray[i].validFlag==NO)
		{
			if(contigItemArray[i].contigTitle)
			{
				free(contigItemArray[i].contigTitle);
				contigItemArray[i].contigTitle = NULL;
			}
			if(contigItemArray[i].contigSeq)
			{
				free(contigItemArray[i].contigSeq);
				contigItemArray[i].contigSeq = NULL;
			}
			if(contigItemArray[i].contigReadArrayEnd5)
			{
				free(contigItemArray[i].contigReadArrayEnd5);
				contigItemArray[i].contigReadArrayEnd5 = NULL;
			}
			if(contigItemArray[i].contigReadArrayEnd3)
			{
				free(contigItemArray[i].contigReadArrayEnd3);
				contigItemArray[i].contigReadArrayEnd3 = NULL;
			}
			contigItemArray[i].contigReadNumEnd5 = 0;
			contigItemArray[i].contigReadNumEnd3 = 0;
			if(contigItemArray[i].contigEdgeArrayEnd5)
			{
				free(contigItemArray[i].contigEdgeArrayEnd5);
				contigItemArray[i].contigEdgeArrayEnd5 = NULL;
			}
			if(contigItemArray[i].contigEdgeArrayEnd3)
			{
				free(contigItemArray[i].contigEdgeArrayEnd3);
				contigItemArray[i].contigEdgeArrayEnd3 = NULL;
			}
			contigItemArray[i].itemNumContigEdgeArrayEnd5 = 0;
			contigItemArray[i].itemNumContigEdgeArrayEnd3 = 0;

			for(j=i; j<contigGraph->itemNumContigItemArray-1; j++)
			{
				if(memcpy(contigItemArray+j, contigItemArray+j+1, sizeof(contigGraphItem_t))==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
			contigGraph->itemNumContigItemArray --;
		}
		else
		{
			i ++;
		}
	}

	return SUCCESSFUL;
}


