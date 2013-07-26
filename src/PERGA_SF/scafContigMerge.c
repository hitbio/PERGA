/*
 * scafContigMerge.c
 *
 *  Created on: Jul 5, 2011
 *      Author: zhuxiao
 */

#include "inc/scafStdinc.h"
#include "inc/scafExtvab.h"


/**
 * Generate the scaffold sequence based on the overlap information.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateScaffoldSequence(const char *scafSeqFile, const char *contigOverlapInfoFile, const char *newContigFile)
{
	printf("=========== Begin generating scaffold sequences, please wait ...\n");

	// initialize the memory for generating scaffold sequence
	if(initMemGeSeq(contigOverlapInfoFile, newContigFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// generate the scaffold sequence
	if(generateScafSeq(scafSeqFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}


	// ######################## Debug information #####################
#if DEBUG_FLAG
	// Output the regions of contig sequences in scaffolds to text file
	//  Format:
	//  	(1) Header: > scaffoldID, linked contigs number, scaffold length;
	// 		(2) Body: contigID, contigLen, contig orientation, start region position, end region position.
	if(outputContigSeqRegionsInScaf("contigRegPosInScaf.txt", contigOverlapIndexArr, scaffoldsNumInCOI, contigOverlapInfoArr, contigInfoArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot output regions of contig sequences in scaffolds, error!\n", __LINE__, __func__);
		return FAILED;
	}
#endif
	// ######################## Debug information #####################

	freeMemGeSeq();

	printf("=========== End generated scaffold sequences.\n");

	return SUCCESSFUL;
}

/**
 * Initialize the memory for generating scaffold sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemGeSeq(const char *contigOverlapInfoFile, const char *newContigFile)
{
	// initialize the contig information array
	if(initContigInfoArray(newContigFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the contig information array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(loadContigOverlapInfo(contigOverlapInfoFile, &contigOverlapIndexArr, &scaffoldsNumInCOI, &contigOverlapInfoArr, &itemNumInContigOverlapArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the contig overlap information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	minBaseNumInGap = MIN_BASENUM_IN_GAP;

	return SUCCESSFUL;
}

/**
 * Free the memory for generating scaffold sequence.
 */
void freeMemGeSeq()
{
	// free contig information array
	freeMemContigInfo(&contigInfoArr, &contigsNum);

	// free the contig overlap information
	freeContigOverlapInfo(&contigOverlapIndexArr, &scaffoldsNumInCOI, &contigOverlapInfoArr, &itemNumInContigOverlapArr);

	minBaseNumInGap = 0;
}


/**
 * Generate scaffold sequence.
 *  File format:
 *  	(1) Head field: >scaffoldID, contigsNum, sequence length, which are separated by tab character;
 *  	(2) Body field: base sequence including unknown bases, which are separated by tab character.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateScafSeq(const char *scafSeqFile)
{
	FILE *fpScafSeq;
	int i, j, k;
	char *newScafSeq, *tmpContigSeq;
	int maxContigLen, approximateScaffoldLen, scaffoldLen;
	contigOverlap *pContigOverlapInfo;
	int scaffoldID, linkedContigsNum, rowsNum, startRow;
	int contigID1, contigOrient1, contigLen1, contigID2, contigOrient2, contigLen2;
	int mergeFlag, overlapLen, gapSize, breakFlag, newGapSize;

	fpScafSeq = fopen(scafSeqFile, "w");
	if(fpScafSeq==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, scafSeqFile);
		return FAILED;
	}

	for(i=0; i<scaffoldsNumInCOI; i++)
	{
		// ####################### Debug information #####################
		//printf("Begin generate sequence for scaffold %d ...\n", contigOverlapIndexArr[i].scaffoldID);
		// ####################### Debug information #####################

		scaffoldID = contigOverlapIndexArr[i].scaffoldID;
		linkedContigsNum = contigOverlapIndexArr[i].linkedNum;
		rowsNum = contigOverlapIndexArr[i].rowsNum;
		startRow = contigOverlapIndexArr[i].startRow;

		// ####################### Debug information #####################
		//if(scaffoldID==30)
		//{
		//	printf("line=%d, In %s(), scaffoldID=%d, linkedContigsNum=%d, rowsNum=%d, startRow=%d\n", __LINE__, __func__, scaffoldID, linkedContigsNum, rowsNum, startRow);
		//}
		// ####################### Debug information #####################

		// generate the consensus sequence
		pContigOverlapInfo = contigOverlapInfoArr + startRow;
		scaffoldLen = 0;

		// deal with the head link
		if(linkedContigsNum==1)
		{ // only one contig in the scaffold, then output the contig directly
			contigID1 = pContigOverlapInfo->contigID1;
			scaffoldLen = contigInfoArr[contigID1-1].contigLen;
			fprintf(fpScafSeq, ">%d\t%d\t%d\n", scaffoldID, linkedContigsNum, scaffoldLen);
			fprintf(fpScafSeq, "%s\n", contigInfoArr[contigID1-1].contigSeq);
		}else
		{
			// get the maximal contig length
			if(getMaxContigLenInSingleScafflod(contigOverlapIndexArr+i, &maxContigLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the maximal contig length for scafflod %d, error!\n", __LINE__, __func__, scaffoldID);
				return FAILED;
			}

			// get the approximate single scaffold sequence length
			if(getApproximateSingleScafflodLen(contigOverlapIndexArr+i, &approximateScaffoldLen)==FAILED)
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
			contigID1 = pContigOverlapInfo[0].contigID1;
			contigOrient1 = pContigOverlapInfo[0].orientation1;
			contigLen1 = contigInfoArr[contigID1-1].contigLen;
			contigID2 = pContigOverlapInfo[0].contigID2;
			contigOrient2 = pContigOverlapInfo[0].orientation2;
			contigLen2 = contigInfoArr[contigID2-1].contigLen;
			mergeFlag = pContigOverlapInfo[0].mergeFlag;
			overlapLen = pContigOverlapInfo[0].overlapLen;
			gapSize = pContigOverlapInfo[0].gapSize;
			breakFlag = pContigOverlapInfo[0].breakFlag;

			if(breakFlag==YES)
			{
				printf("line=%d, In %s(), the scaffold should be broken, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// generate the scaffold sequence of the two contigs
			if(contigOrient1==ORIENTATION_PLUS)
			{ // plus orientation of the first contig, then copy it to newScafSeq
				strcpy(newScafSeq, contigInfoArr[contigID1-1].contigSeq);
				scaffoldLen += contigLen1;
			}else
			{ // minus orientation of the first contig, reverse the contig sequence, and then copy it to newScafSeq
				strcpy(tmpContigSeq, contigInfoArr[contigID1-1].contigSeq); // copy sequence

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
					strcpy(newScafSeq+scaffoldLen, contigInfoArr[contigID2-1].contigSeq+overlapLen);
					scaffoldLen += contigLen2 - overlapLen;
				}else
				{ // minus orientation of the second contig, reverse the contig sequence, and then copy it to newScafSeq
					strcpy(tmpContigSeq, contigInfoArr[contigID2-1].contigSeq); // copy sequence

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
					strcpy(newScafSeq+scaffoldLen, contigInfoArr[contigID2-1].contigSeq);
					scaffoldLen += contigLen2;
				}else
				{ // minus orientation of the second contig, reverse the contig sequence, and then copy it to newScafSeq
					strcpy(tmpContigSeq, contigInfoArr[contigID2-1].contigSeq); // copy sequence

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
				contigID2 = pContigOverlapInfo[j].contigID2;
				contigOrient2 = pContigOverlapInfo[j].orientation2;
				contigLen2 = contigInfoArr[contigID2-1].contigLen;
				mergeFlag = pContigOverlapInfo[j].mergeFlag;
				overlapLen = pContigOverlapInfo[j].overlapLen;
				gapSize = pContigOverlapInfo[j].gapSize;
				breakFlag = pContigOverlapInfo[j].breakFlag;

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
						strcpy(newScafSeq+scaffoldLen, contigInfoArr[contigID2-1].contigSeq+overlapLen);
						scaffoldLen += contigLen2 - overlapLen;
					}else
					{ // minus orientation of the second contig, reverse the contig sequence, and then copy it to newScafSeq
						strcpy(tmpContigSeq, contigInfoArr[contigID2-1].contigSeq); // copy sequence

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
						strcpy(newScafSeq+scaffoldLen, contigInfoArr[contigID2-1].contigSeq);
						scaffoldLen += contigLen2;
					}else
					{ // minus orientation of the second contig, reverse the contig sequence, and then copy it to newScafSeq
						strcpy(tmpContigSeq, contigInfoArr[contigID2-1].contigSeq); // copy sequence

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

			// output the scaffold sequence
			fprintf(fpScafSeq, ">%d\t%d\t%d\n", scaffoldID, linkedContigsNum, scaffoldLen);
			fprintf(fpScafSeq, "%s\n", newScafSeq);
			fflush(fpScafSeq);

			// free the auxiliary memory
			free(newScafSeq);
			newScafSeq = NULL;
			free(tmpContigSeq);
			tmpContigSeq = NULL;
		}

		// ####################### Debug information #####################
		//printf("End generating sequence for scaffold %d.\n", contigOverlapIndexArr[i].scaffoldID);
		// ####################### Debug information #####################
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
short getMaxContigLenInSingleScafflod(contigOverlapIndex *pContigOverlapIndex, int *maxContigLen)
{
	int i;
	int linkedContigsNum, rowsNum, startRow;
	contigOverlap *pContigOverlapInfo;
	int contigID1, contigLen1, contigID2, contigLen2, breakFlag;

	linkedContigsNum = pContigOverlapIndex->linkedNum;
	rowsNum = pContigOverlapIndex->rowsNum;
	startRow = pContigOverlapIndex->startRow;

	*maxContigLen = 0;
	// process the first two contigs in head link
	pContigOverlapInfo = contigOverlapInfoArr + startRow;
	contigID1 = pContigOverlapInfo[0].contigID1;
	contigLen1 = contigInfoArr[contigID1-1].contigLen;
	contigID2 = pContigOverlapInfo[0].contigID2;
	contigLen2 = contigInfoArr[contigID2-1].contigLen;
	breakFlag = pContigOverlapInfo[0].breakFlag;

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
		contigID2 = pContigOverlapInfo[i].contigID2;
		contigLen2 = contigInfoArr[contigID2-1].contigLen;
		breakFlag = pContigOverlapInfo[i].breakFlag;

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
short getApproximateSingleScafflodLen(contigOverlapIndex *pContigOverlapIndex, int *approximateScaffoldLen)
{
	int i;
	int linkedContigsNum, rowsNum, startRow;
	contigOverlap *pContigOverlapInfo;
	int contigID1, contigLen1, contigID2, contigLen2, mergeFlag, gapSize, breakFlag, newGapSize;

	linkedContigsNum = pContigOverlapIndex->linkedNum;
	rowsNum = pContigOverlapIndex->rowsNum;
	startRow = pContigOverlapIndex->startRow;

	*approximateScaffoldLen = 0;
	// process the first two contigs in head link
	pContigOverlapInfo = contigOverlapInfoArr + startRow;
	contigID1 = pContigOverlapInfo[0].contigID1;
	contigLen1 = contigInfoArr[contigID1-1].contigLen;
	contigID2 = pContigOverlapInfo[0].contigID2;
	contigLen2 = contigInfoArr[contigID2-1].contigLen;
	mergeFlag = pContigOverlapInfo[0].mergeFlag;
	gapSize = pContigOverlapInfo[0].gapSize;
	breakFlag = pContigOverlapInfo[0].breakFlag;

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
		contigID2 = pContigOverlapInfo[i].contigID2;
		contigLen2 = contigInfoArr[contigID2-1].contigLen;
		mergeFlag = pContigOverlapInfo[i].mergeFlag;
		gapSize = pContigOverlapInfo[i].gapSize;
		breakFlag = pContigOverlapInfo[i].breakFlag;

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


