/*
 * scafSeq.c
 *
 *  Created on: Dec 22, 2013
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Generate scaffold sequence.
 *  File format:
 *  	(1) Head field: >scaffoldID, contigsNum, sequence length, which are separated by tab character;
 *  	(2) Body field: base sequence including unknown bases, which are separated by tab character.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateScafSeq(char *scafSeqFile, scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph)
{
	FILE *fpScafSeq;
	int32_t i, j, k;
	char *newScafSeq, *tmpContigSeq;
	int32_t maxContigLen, approximateScaffoldLen, scaffoldLen;
	scaffoldItem_t *scaffoldItem;
	contigOverlap_t *pContigOverlapInfo, *pContigOverlapArray;
	int32_t scaffoldID, linkedContigsNum, rowsNum, startRow;
	int32_t contigID1, contigOrient1, contigLen1, contigID2, contigOrient2, contigLen2;
	int32_t mergeFlag, overlapLen, gapSize, breakFlag, newGapSize;

	fpScafSeq = fopen(scafSeqFile, "w");
	if(fpScafSeq==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, scafSeqFile);
		return FAILED;
	}

	scaffoldItem = scaffoldSet->scaffoldItemList;
	while(scaffoldItem)
	{
		scaffoldID = scaffoldItem->scaffoldID;
		linkedContigsNum = scaffoldItem->linkedContigsNum;
		pContigOverlapArray = scaffoldItem->contigOverlapArray;
		rowsNum = scaffoldItem->itemNumContigOverlapArray;

		// ####################### Debug information #####################
		//printf("Begin generate sequence for scaffold %d ...\n", scaffoldID);
		// ####################### Debug information #####################

		// generate the consensus sequence
		scaffoldLen = 0;

		// deal with the head link
		if(linkedContigsNum==1)
		{ // only one contig in the scaffold, then output the contig directly
			contigID1 = pContigOverlapArray[0].contigID1;
			scaffoldLen = contigGraph->contigItemArray[contigID1-1].contigLen;
			scaffoldItem->scaffoldLen = scaffoldLen;
			fprintf(fpScafSeq, ">%d\t%d\t%d\n", scaffoldID, linkedContigsNum, scaffoldLen);
			fprintf(fpScafSeq, "%s\n", contigGraph->contigItemArray[contigID1-1].contigSeq);
		}else
		{
			// get the maximal contig length
			if(getMaxContigLenInSingleScafflod(scaffoldItem, &maxContigLen, contigGraph)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the maximal contig length for scafflod %d, error!\n", __LINE__, __func__, scaffoldID);
				return FAILED;
			}

			// get the approximate single scaffold sequence length
			if(getApproximateSingleScafflodLen(scaffoldItem, &approximateScaffoldLen, contigGraph)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the maximal contig length for scafflod %d, error!\n", __LINE__, __func__, scaffoldID);
				return FAILED;
			}

			// allocate the auxiliary memory
			tmpContigSeq = (char *) calloc (maxContigLen+1, sizeof(char));
			if(tmpContigSeq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			newScafSeq = (char *) calloc (approximateScaffoldLen+1, sizeof(char));
			if(newScafSeq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// deal with the two contigs in head link
			contigID1 = pContigOverlapArray[0].contigID1;
			contigOrient1 = pContigOverlapArray[0].orientation1;
			contigLen1 = contigGraph->contigItemArray[contigID1-1].contigLen;
			contigID2 = pContigOverlapArray[0].contigID2;
			contigOrient2 = pContigOverlapArray[0].orientation2;
			contigLen2 = contigGraph->contigItemArray[contigID2-1].contigLen;
			mergeFlag = pContigOverlapArray[0].mergeFlag;
			overlapLen = pContigOverlapArray[0].overlapLen;
			gapSize = pContigOverlapArray[0].gapSize;
			breakFlag = pContigOverlapArray[0].breakFlag;

			if(breakFlag==YES)
			{
				printf("line=%d, In %s(), the scaffold should be broken, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// generate the scaffold sequence of the two contigs
			if(contigOrient1==ORIENTATION_PLUS)
			{ // plus orientation of the first contig, then copy it to newScafSeq
				strcpy(newScafSeq, contigGraph->contigItemArray[contigID1-1].contigSeq);
				scaffoldLen += contigLen1;
			}else
			{ // minus orientation of the first contig, reverse the contig sequence, and then copy it to newScafSeq
				strcpy(tmpContigSeq, contigGraph->contigItemArray[contigID1-1].contigSeq); // copy sequence

				if(reverseSeq(tmpContigSeq, contigLen1)==FAILED) // reverse the sequence
				{
					printf("line=%d, In %s(), cannot reverse sequence, error!\n", __LINE__, __func__);
					return FAILED;
				}

				strcpy(newScafSeq, tmpContigSeq); // copy to newScafSeq
				scaffoldLen += contigLen1;
			}

			// deal with the overlap or gap
			if(mergeFlag==YES)
			{ // an overlap happens
				if(contigOrient2==ORIENTATION_PLUS)
				{ // plus orientation of the second contig, then copy it to newScafSeq excluding the overlapped bases
					strcpy(newScafSeq+scaffoldLen, contigGraph->contigItemArray[contigID2-1].contigSeq+overlapLen);
					scaffoldLen += contigLen2 - overlapLen;
				}else
				{ // minus orientation of the second contig, reverse the contig sequence, and then copy it to newScafSeq
					strcpy(tmpContigSeq, contigGraph->contigItemArray[contigID2-1].contigSeq); // copy sequence

					if(reverseSeq(tmpContigSeq, contigLen2)==FAILED) // reverse the sequence
					{
						printf("line=%d, In %s(), cannot reverse sequence, error!\n", __LINE__, __func__);
						return FAILED;
					}

					strcpy(newScafSeq+scaffoldLen, tmpContigSeq+overlapLen); // copy to newScafSeq
					scaffoldLen += contigLen2 - overlapLen;
				}
			}else
			{ // a gap happens

				if(gapSize<minBaseNumInGap)
					newGapSize = minBaseNumInGap;
				else
					newGapSize = gapSize;

				// fill unknown bases
				for(k=0; k<newGapSize; k++) newScafSeq[scaffoldLen+k] = 'N';
				newScafSeq[scaffoldLen+newGapSize] = '\0';
				scaffoldLen += newGapSize;

				// process the second contig
				if(contigOrient2==ORIENTATION_PLUS)
				{ // plus orientation of the second contig, then copy it to newScafSeq
					strcpy(newScafSeq+scaffoldLen, contigGraph->contigItemArray[contigID2-1].contigSeq);
					scaffoldLen += contigLen2;
				}else
				{ // minus orientation of the second contig, reverse the contig sequence, and then copy it to newScafSeq
					strcpy(tmpContigSeq, contigGraph->contigItemArray[contigID2-1].contigSeq); // copy sequence

					if(reverseSeq(tmpContigSeq, contigLen2)==FAILED) // reverse the sequence
					{
						printf("line=%d, In %s(), cannot reverse sequence, error!\n", __LINE__, __func__);
						return FAILED;
					}

					strcpy(newScafSeq+scaffoldLen, tmpContigSeq); // copy to newScafSeq
					scaffoldLen += contigLen2;
				}
			}

			// deal with other contigs
			for(j=1; j<rowsNum; j++)
			{
				contigID2 = pContigOverlapArray[j].contigID2;
				contigOrient2 = pContigOverlapArray[j].orientation2;
				contigLen2 = contigGraph->contigItemArray[contigID2-1].contigLen;
				mergeFlag = pContigOverlapArray[j].mergeFlag;
				overlapLen = pContigOverlapArray[j].overlapLen;
				gapSize = pContigOverlapArray[j].gapSize;
				breakFlag = pContigOverlapArray[j].breakFlag;

				if(breakFlag==YES)
				{
					printf("line=%d, In %s(), the scaffold should be broken, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// deal with the overlap or gap
				if(mergeFlag==YES)
				{ // an overlap happens
					if(contigOrient2==ORIENTATION_PLUS)
					{ // plus orientation of the second contig, then copy it to newScafSeq excluding the overlapped bases
						strcpy(newScafSeq+scaffoldLen, contigGraph->contigItemArray[contigID2-1].contigSeq+overlapLen);
						scaffoldLen += contigLen2 - overlapLen;
					}else
					{ // minus orientation of the second contig, reverse the contig sequence, and then copy it to newScafSeq
						strcpy(tmpContigSeq, contigGraph->contigItemArray[contigID2-1].contigSeq); // copy sequence

						if(reverseSeq(tmpContigSeq, contigLen2)==FAILED) // reverse the sequence
						{
							printf("line=%d, In %s(), cannot reverse sequence, error!\n", __LINE__, __func__);
							return FAILED;
						}

						strcpy(newScafSeq+scaffoldLen, tmpContigSeq+overlapLen); // copy to newScafSeq
						scaffoldLen += contigLen2 - overlapLen;
					}
				}else
				{ // a gap happens

					if(gapSize<minBaseNumInGap)
						newGapSize = minBaseNumInGap;
					else
						newGapSize = gapSize;

					// fill unknown bases
					for(k=0; k<newGapSize; k++) newScafSeq[scaffoldLen+k] = 'N';
					newScafSeq[scaffoldLen+newGapSize] = '\0';
					scaffoldLen += newGapSize;

					// process the second contig
					if(contigOrient2==ORIENTATION_PLUS)
					{ // plus orientation of the second contig, then copy it to newScafSeq
						strcpy(newScafSeq+scaffoldLen, contigGraph->contigItemArray[contigID2-1].contigSeq);
						scaffoldLen += contigLen2;
					}else
					{ // minus orientation of the second contig, reverse the contig sequence, and then copy it to newScafSeq
						strcpy(tmpContigSeq, contigGraph->contigItemArray[contigID2-1].contigSeq); // copy sequence

						if(reverseSeq(tmpContigSeq, contigLen2)==FAILED) // reverse the sequence
						{
							printf("line=%d, In %s(), cannot reverse sequence, error!\n", __LINE__, __func__);
							return FAILED;
						}

						strcpy(newScafSeq+scaffoldLen, tmpContigSeq); // copy to newScafSeq
						scaffoldLen += contigLen2;
					}
				}
			}

			scaffoldItem->scaffoldLen = scaffoldLen;

			// output the scaffold sequence
			fprintf(fpScafSeq, ">%d\t%d\t%d\n", scaffoldID, linkedContigsNum, scaffoldLen);
			fprintf(fpScafSeq, "%s\n", newScafSeq);
			//fflush(fpScafSeq);

			// free the auxiliary memory
			free(newScafSeq);
			newScafSeq = NULL;
			free(tmpContigSeq);
			tmpContigSeq = NULL;
		}

		scaffoldItem = scaffoldItem->next;
	}

	fclose(fpScafSeq);
	fpScafSeq = NULL;

	return SUCCESSFUL;
}

/**
 * Get the maximal contig length in single scaffold.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxContigLenInSingleScafflod(scaffoldItem_t *scaffoldItem, int32_t *maxContigLen, contigGraph_t *contigGraph)
{
	int32_t i, linkedContigsNum, rowsNum, startRow;
	contigOverlap_t *pContigOverlapArray;
	int32_t contigID1, contigLen1, contigID2, contigLen2, breakFlag;

	linkedContigsNum = scaffoldItem->linkedContigsNum;
	rowsNum = scaffoldItem->itemNumContigOverlapArray;
	pContigOverlapArray = scaffoldItem->contigOverlapArray;

	*maxContigLen = 0;
	// process the first two contigs in head link
	contigID1 = pContigOverlapArray[0].contigID1;
	contigLen1 = contigGraph->contigItemArray[contigID1-1].contigLen;
	contigID2 = pContigOverlapArray[0].contigID2;
	contigLen2 = contigGraph->contigItemArray[contigID2-1].contigLen;
	breakFlag = pContigOverlapArray[0].breakFlag;

	if(breakFlag==YES)
	{
		printf("line=%d, In %s(), the scaffold should be broken, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(contigLen1>(*maxContigLen))
		*maxContigLen = contigLen1;

	if(contigLen2>(*maxContigLen))
		*maxContigLen = contigLen2;

	// process remained contigs
	for(i=1; i<rowsNum; i++)
	{
		contigID2 = pContigOverlapArray[i].contigID2;
		contigLen2 = contigGraph->contigItemArray[contigID2-1].contigLen;
		breakFlag = pContigOverlapArray[i].breakFlag;

		if(breakFlag==YES)
		{
			printf("line=%d, In %s(), the scaffold is broken, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(contigLen2>(*maxContigLen))
			*maxContigLen = contigLen2;
	}

	//################### Debug information ##################
	//printf("maxContigLen=%d\n", *maxContigLen);
	//################### Debug information ##################

	if((*maxContigLen)<=0)
	{
		printf("line=%d, In %s(), maxContigLen=%d, error!\n", __LINE__, __func__, *maxContigLen);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get approximate scaffold length for single scaffold.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getApproximateSingleScafflodLen(scaffoldItem_t *scaffoldItem, int32_t *approximateScaffoldLen, contigGraph_t *contigGraph)
{
	int32_t i, linkedContigsNum, rowsNum, startRow;
	contigOverlap_t *pContigOverlapArray;
	int32_t contigID1, contigLen1, contigID2, contigLen2, mergeFlag, gapSize, breakFlag, newGapSize;

	linkedContigsNum = scaffoldItem->linkedContigsNum;
	rowsNum = scaffoldItem->itemNumContigOverlapArray;
	pContigOverlapArray = scaffoldItem->contigOverlapArray;

	*approximateScaffoldLen = 0;
	// process the first two contigs in head link
	contigID1 = pContigOverlapArray[0].contigID1;
	contigLen1 = contigGraph->contigItemArray[contigID1-1].contigLen;
	contigID2 = pContigOverlapArray[0].contigID2;
	contigLen2 = contigGraph->contigItemArray[contigID2-1].contigLen;
	mergeFlag = pContigOverlapArray[0].mergeFlag;
	gapSize = pContigOverlapArray[0].gapSize;
	breakFlag = pContigOverlapArray[0].breakFlag;

	if(breakFlag==YES)
	{
		printf("line=%d, In %s(), the scaffold should be broken, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(mergeFlag==YES)
	{ // a overlap happens, then only consider the lengths of two contigs
		*approximateScaffoldLen += contigLen1 + contigLen2;
	}else
	{ // a gap happens, then consider the length of the two contigs and the gap size
		if(gapSize<minBaseNumInGap)
			newGapSize = minBaseNumInGap;
		else
			newGapSize = gapSize;

		*approximateScaffoldLen += contigLen1 + contigLen2 + newGapSize;
	}

	// process remained contigs
	for(i=1; i<rowsNum; i++)
	{
		contigID2 = pContigOverlapArray[i].contigID2;
		contigLen2 = contigGraph->contigItemArray[contigID2-1].contigLen;
		mergeFlag = pContigOverlapArray[i].mergeFlag;
		gapSize = pContigOverlapArray[i].gapSize;
		breakFlag = pContigOverlapArray[i].breakFlag;

		if(breakFlag==YES)
		{
			printf("line=%d, In %s(), the scaffold is broken, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(mergeFlag==YES)
		{ // a overlap happens, consider the length of the second contig
			*approximateScaffoldLen += contigLen2;
		}else
		{
			if(gapSize<minBaseNumInGap)
				newGapSize = minBaseNumInGap;
			else
				newGapSize = gapSize;

			*approximateScaffoldLen += contigLen2 + newGapSize;
		}
	}

	//################### Debug information ##################
	//printf("approximateScaffoldLen=%d\n", *approximateScaffoldLen);
	//################### Debug information ##################

	if((*approximateScaffoldLen)<=0)
	{
		printf("line=%d, In %s(), approximateScaffoldLen=%d, error!\n", __LINE__, __func__, *approximateScaffoldLen);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Statistics for contigs lengths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short contigsLenStatisticsScaf(scaffoldSet_t *scaffoldSet, int32_t minContigLen)
{
	int32_t i, *scaffoldsLenArray, *scaffoldsLenBufTmp, itemNumScaffoldsLenArray, newContigNum;
	scaffoldItem_t *scaffoldItem;

	itemNumScaffoldsLenArray = scaffoldSet->scaffoldNum;

	scaffoldsLenArray = (int32_t *) calloc (itemNumScaffoldsLenArray, sizeof(int32_t));
	if(scaffoldsLenArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	scaffoldsLenBufTmp = (int32_t *) calloc (itemNumScaffoldsLenArray, sizeof(int32_t));
	if(scaffoldsLenBufTmp==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	newContigNum = 0;
	scaffoldItem = scaffoldSet->scaffoldItemList;
	while(scaffoldItem)
	{
		if(scaffoldItem->scaffoldLen>=minContigLen)
		{
			scaffoldsLenArray[newContigNum] = scaffoldItem->scaffoldLen;
			newContigNum ++;
		}

		scaffoldItem = scaffoldItem->next;
	}

	// radix sort of contigsLenArray
	if(radixSortContigsLen(scaffoldsLenArray, scaffoldsLenBufTmp, newContigNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot sort contigs lengths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get lengths statistics
	if(computeLenStatistics(scaffoldsLenArray, newContigNum, minContigLen, "scaffolds")==FAILED)
	{
		printf("line=%d, In %s(), cannot compute contigs length statistics, error!\n", __LINE__, __func__);
		return FAILED;
	}

	free(scaffoldsLenArray);
	free(scaffoldsLenBufTmp);

	return SUCCESSFUL;
}
