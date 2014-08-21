/*
 * scafOverlap.c
 *
 *  Created on: Dec 19, 2013
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Generate overlap information between contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short overlapContigsInScaf(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph, readSet_t *readSet)
{
	// initialize the memory
	if(initMemContigOverlap(contigGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// compute the overlaps
	if(computeContigOverlapInfo(scaffoldSet, contigGraph, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate contig overlap information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// confirm the overlapped contigs
	if(confirmOverlapInfoInScaf(scaffoldSet, contigGraph, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot confirm contig overlap information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// split broken scaffolds
	if(splitScaffolds(scaffoldSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate contig overlap information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// free memory
	freeMemContigOverlap();

	return SUCCESSFUL;
}

/**
 * Initialize the memory for contig overlap.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemContigOverlap(contigGraph_t *contigGraph)
{
	int32_t i;

	// the auxiliary memory for overlapping
	minLinksNumContigsThres = contigGraph->averLinkNum * MIN_FIRST_LINKNUM_FACTOR;
	if(minLinksNumContigsThres>MIN_FIRST_LINKNUM_THRES)
	{
		minLinksNumContigsThres = MIN_FIRST_LINKNUM_THRES;
	}
	//else if(minLinksNumContigsThres<MIN_FIRST_LINKNUM_THRES)
	//{
	//	minLinksNumContigsThres = MIN_FIRST_LINKNUM_THRES;
	//}

	maxOverlapSeqLen = 2 * readLen;
	minOverlapThres = MIN_OVERLAP_THRESHOLD;
	minExactOverlapThres = MIN_EXACT_OVERLAP_THRESHOLD;
	mismatchThres = MAX_MISMATCH_THRESHOLD;
	maxMisMatchRatioThres = MAX_MISMATCH_RATIO_THRESHOLD;
	subMismatchThres = SUB_MISMATCH_FACTOR * MAX_MISMATCH_THRESHOLD;
	minAdjustGapSizeThres = MIN_ADJUST_GAP_SIZE_THRESHOLD;
	maxAdjustGapSizeThres = MAX_ADJUST_GAP_SIZE_THRESHOLD;
	minBaseNumInGap = MIN_BASENUM_IN_GAP;
	gapSizeSdevFactorOverlap = GAP_SIZE_SDEV_FACTOR_OVERLAP;
	//exactOverlapSdevThres = EXACT_OVERLAP_SDEV_THRESHOLD;
	matchScore = MATCH_SCORE;
	mismatchScore = MISMATCH_SCORE;
	gapScore = GAP_SCORE;

	breakLinkNumThres = BREAK_LINKNUM_FACTOR * minLinksNumContigsThres;
	if(breakLinkNumThres>MIN_BREAK_LINKNUM_THRES)
		breakLinkNumThres = MIN_BREAK_LINKNUM_THRES;

	overlapSeq1 = (char *) calloc((maxOverlapSeqLen+1)*2, sizeof(char));
	if(overlapSeq1==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	overlapSeq2 = (char *) calloc((maxOverlapSeqLen+1)*2, sizeof(char));
	if(overlapSeq2==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	scoreArr = (int *) calloc((maxOverlapSeqLen+1)*(maxOverlapSeqLen+1), sizeof(int));
	if(scoreArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<3; i++)
	{
		alignResultArr[i] = (char *) calloc (2*(maxOverlapSeqLen+1), sizeof(char));
		if(alignResultArr[i]==NULL)
		{
			printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Free the memory of contigs overlap.
 */
void freeMemContigOverlap()
{
	int i;

	free(overlapSeq1);
	overlapSeq1 = NULL;
	free(overlapSeq2);
	overlapSeq2 = NULL;
	free(scoreArr);
	scoreArr = NULL;

	for(i=0; i<3; i++)
	{
		free(alignResultArr[i]);
		alignResultArr[i] = NULL;
	}
}

/**
 * Compute the contig overlap information.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeContigOverlapInfo(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph, readSet_t *readSet)
{
	int32_t scaffoldID, rowsNum;
	contigOverlap_t *contigOverlapArray;
	scaffoldItem_t *scaffoldItem;
	int32_t j;

	// generate contig overlap information of scaffolds
	scaffoldItem = scaffoldSet->scaffoldItemList;
	while(scaffoldItem)
	{
		//####################### Deubg information #####################
#if DEBUG_SCAF_OVERLAP_FLAG
		printf("\n=================== scaffold: %d, linkedContigsNum: %d, itemNumContigOverlapArray:  ================\n", scaffoldItem->scaffoldID, scaffoldItem->linkedContigsNum, scaffoldItem->itemNumContigOverlapArray);
#endif
		//####################### Deubg information #####################

		scaffoldID = scaffoldItem->scaffoldID;
		rowsNum = scaffoldItem->itemNumContigOverlapArray;
		contigOverlapArray = scaffoldItem->contigOverlapArray;

		if(scaffoldItem->linkedContigsNum>=2)
		{ // just consider the scaffolds having more than 1 linked contigs
			for(j=0; j<rowsNum; j++)
			{
				//####################### Debug information ###################
#if DEBUG_SCAF_OVERLAP_FLAG
				if(contigOverlapArray[j].contigID1==763 && contigOverlapArray[j].contigID2==8914)
				{
					printf("contig1: [%d,%d,%d]; contig2: [%d,%d,%d]\n", contigOverlapArray[j].contigID1, contigOverlapArray[j].orientation1, contigGraph->contigItemArray[contigOverlapArray[j].contigID1-1].contigLen, contigOverlapArray[j].contigID2, contigOverlapArray[j].orientation2, contigGraph->contigItemArray[contigOverlapArray[j].contigID2-1].contigLen);
				}

				printf("*****contig1: [%d,%d,%d]; contig2: [%d,%d,%d]\n", contigOverlapArray[j].contigID1, contigOverlapArray[j].orientation1, contigGraph->contigItemArray[contigOverlapArray[j].contigID1-1].contigLen, contigOverlapArray[j].contigID2, contigOverlapArray[j].orientation2, contigGraph->contigItemArray[contigOverlapArray[j].contigID2-1].contigLen);
#endif
				//####################### Debug information ###################

				if(updateContigOverlapLen(contigOverlapArray+j, contigGraph, readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute overlap length for scaffold: %d, error!\n", __LINE__, __func__, scaffoldID);
					return FAILED;
				}
			}
		}else if(scaffoldItem->linkedContigsNum==1)
		{
			//####################### Debug information ###################
//#if DEBUG_SCAF_FLAG
//			printf("*****contig1: [%d,0,%d]; contig2: [0,0,0]\n", contigOverlapArray[0].contigID1, contigGraph->contigItemArray[contigOverlapArray[0].contigID1-1].contigLen);
//#endif
			//####################### Debug information ###################
		}else
		{
			printf("line=%d, In %s(), contig overlap information error!\n", __LINE__, __func__);
			return FAILED;
		}

		scaffoldItem = scaffoldItem->next;
	}

	return SUCCESSFUL;
}

/**
 * Compute the overlap length of two contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateContigOverlapLen(contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph, readSet_t *readSet)
{
	int32_t contigID1, orientation1, contigLen1, contigID2, orientation2, contigLen2;
	int32_t seq_len1, seq_len2, overlapLenExact, overlapLenAlignment, overlapLenAdjust, mismatchNum;
	double gapSize;
	int32_t validPairedNum_gapEstimate;
	int32_t contigLenBeforeAlignment[2];
	int32_t contigEndFlag1, contigEndFlag2;

	// get the information of contg1 and contig2
	contigID1 = pContigOverlapInfo->contigID1;
	orientation1 = pContigOverlapInfo->orientation1;
	contigLen1 = contigGraph->contigItemArray[contigID1-1].contigLen;
	contigID2 = pContigOverlapInfo->contigID2;
	orientation2 = pContigOverlapInfo->orientation2;
	contigLen2 = contigGraph->contigItemArray[contigID2-1].contigLen;

	// get the sequences from the two contigs according to their orientations
	if(contigLen1>=maxOverlapSeqLen)
	{
		seq_len1 = maxOverlapSeqLen;
	}else
	{
		seq_len1 = contigLen1;
	}

	if(orientation1==ORIENTATION_PLUS)
	{
		strncpy(overlapSeq1, contigGraph->contigItemArray[contigID1-1].contigSeq+contigLen1-seq_len1, seq_len1);
		overlapSeq1[seq_len1] = '\0';
	}else
	{
		strncpy(overlapSeq1, contigGraph->contigItemArray[contigID1-1].contigSeq, seq_len1);
		overlapSeq1[seq_len1] = '\0';
		if(reverseSeq(overlapSeq1, seq_len1)==FAILED) // reverse the sequence
		{
			printf("line=%d, In %s(), cannot reverse the sequence [ %s ], error!\n", __LINE__, __func__, overlapSeq1);
			return FAILED;
		}
	}

	if(contigLen2>=maxOverlapSeqLen)
	{
		seq_len2 = maxOverlapSeqLen;
	}else
	{
		seq_len2 = contigLen2;
	}

	if(orientation2==ORIENTATION_PLUS)
	{
		strncpy(overlapSeq2, contigGraph->contigItemArray[contigID2-1].contigSeq, seq_len2);
		overlapSeq2[seq_len2] = '\0';
	}else
	{
		strncpy(overlapSeq2, contigGraph->contigItemArray[contigID2-1].contigSeq+contigLen2-seq_len2, seq_len2);
		overlapSeq2[seq_len2] = '\0';
		if(reverseSeq(overlapSeq2, seq_len2)==FAILED) // reverse the sequence
		{
			printf("line=%d, In %s(), cannot reverse the sequence [ %s ], error!\n", __LINE__, __func__, overlapSeq2);
			return FAILED;
		}
	}

	//############################# Debug information #########################
#if DEBUG_SCAF_OVERLAP_FLAG
	printf("overlapSeq1=%s, len=%d\noverlapSeq2=%s, len=%d\n", overlapSeq1, (int32_t)strlen(overlapSeq1), overlapSeq2, (int32_t)strlen(overlapSeq2));
#endif
	//############################# Debug information #########################

	contigLenBeforeAlignment[0] = contigLen1;
	contigLenBeforeAlignment[1] = contigLen2;

	// estimate gap size between contigs
	if(gapSizeEstimateBetweenContigs(&gapSize, &validPairedNum_gapEstimate, pContigOverlapInfo, contigGraph, readSet, -1, NULL, NULL)==FAILED)
	{
		printf("line=%d, In %s(), cannot estimate the gap size between contigs [%d,%d], error!\n", __LINE__, __func__, contigID1, contigID2);
		return FAILED;
	}

	//printf("line=%d, In %s(), gapSize=%d.\n", __LINE__, __func__, gapSize);

	// closed on 2014-01-02
	//if(validPairedNum_gapEstimate<=0)
	//if(validPairedNum_gapEstimate<minLinksNumContigsThres)
	//if(validPairedNum_gapEstimate<minLinksNumContigsThres || (gapSize>0.6*meanSizeInsert && gapSize>2*averReadLen)) // 2012-11-19
	//if(gapSize<-500 && validPairedNum_gapEstimate>=3)  // 2014-04-05
	if(gapSize<-500 || (gapSize<-300 && validPairedNum_gapEstimate<2))  // 2014-04-19
	{
#if (DEBUG_FLAG==YES)
		printf("line=%d, In %s(), broken.\n", __LINE__, __func__);
#endif
		pContigOverlapInfo->breakFlag = YES;
		pContigOverlapInfo->gapSize = gapSize;
		pContigOverlapInfo->mergeFlag = NO;

		//printf("line=%d, In %s(), validPairedNum_gapEstimate=%d, error!\n", __LINE__, __func__, validPairedNum_gapEstimate);
		return SUCCESSFUL;
	}


	// update the contig overlap information
	//if(gapSize<=stardardDeviationInsert)
	if(gapSize<=gapSizeSdevFactorOverlap*standardDev)
	{
		// compute the overlap length by exact alignment
		if(computeSeqOverlapLenExact(&overlapLenExact, overlapSeq1, seq_len1, overlapSeq2, seq_len2, scoreArr, gapSize)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute sequence overlap length, error!\n", __LINE__, __func__);
			return FAILED;
		}

		//if(overlapLenExact>=minOverlapThres)
		//if(overlapLenExact>=minOverlapThres && (overlapLenExact>=-gapSize-standardDev && overlapLenExact<=-gapSize+standardDev))
		if(overlapLenExact>=minOverlapThres && (overlapLenExact>=-gapSize-2*standardDev && overlapLenExact<=-gapSize+2*standardDev))
		{
			pContigOverlapInfo->mergeFlag = YES;
			pContigOverlapInfo->overlapLen = overlapLenExact;
		}else
		{ // alignment the two sequences
			// pairwise alignment
			if(computeSeqOverlapLenByAlignment(overlapSeq1, seq_len1, overlapSeq2, seq_len2, scoreArr, alignResultArr, &overlapLenAlignment, &mismatchNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot compute sequence overlap length, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(overlapLenAlignment>=minOverlapThres && mismatchNum <= mismatchThres)
			{ // there is an overlap between the two contigs
				overlapLenAdjust = overlapLenAlignment;
				// adjust overlapped sequences
				if(adjustOverlapSeq(overlapSeq1, overlapSeq2, alignResultArr, &overlapLenAdjust)==FAILED)
				{
					printf("line=%d, In %s(), cannot adjust the overlapped sequences, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// update the contig sequences
				if(overlapLenAdjust>0)
				{ // valid overlap, then update contig end sequences
					//############################### Debug information #########################
#if DEBUG_OUT_FLAG
					//printf("line=%d, In %s(), before updateContigs(), contigLen1=%d, contigLen2=%d\n", __LINE__, __func__, pContigInfoArr[contigID1-1].contigLen, pContigInfoArr[contigID2-1].contigLen);
#endif
					//############################### Debug information #########################

					if(updateContigs(contigGraph->contigItemArray+contigID1-1, contigGraph->contigItemArray+contigID2-1, orientation1, orientation2, overlapSeq1, seq_len1, overlapSeq2, seq_len2)==FAILED)
					{
						printf("line=%d, In %s(), cannot update contig end sequences, error!\n", __LINE__, __func__);
						return FAILED;
					}

					//############################### Debug information #########################
#if DEBUG_OUT_FLAG
					//printf("line=%d, In %s(), after updateContigs(), contigLen1=%d, contigLen2=%d\n", __LINE__, __func__, contigGraph->contigItemArray[contigID1-1].contigLen, contigGraph->contigItemArray[contigID2-1].contigLen);
#endif
					//############################### Debug information #########################

					pContigOverlapInfo->mergeFlag = YES;
					pContigOverlapInfo->overlapLen = overlapLenAdjust;

					if(updateContigEndsInfo(pContigOverlapInfo, contigGraph)==FAILED)
					{
						printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
						return FAILED;
					}

				}else // if(overlapLenAdjust==0)
				{ // there is a gap between the two contigs with high probability

					if(validPairedNum_gapEstimate>=minLinksNumContigsThres)
					{
						if(gapSize<minAdjustGapSizeThres)
						{ // gapSize < -10
							pContigOverlapInfo->gapSize = gapSize;
							pContigOverlapInfo->mergeFlag = NO;

							// update the update length by cutting uncovered contig ends
							if(updateOverlapLenByCutUncoveredContigEnds(pContigOverlapInfo, contigGraph, readSet)==FAILED)
							{
								printf("line=%d, In %s(), cannot update the overlap size by cutting uncovered contig ends between contigs [%d,%d], error!\n", __LINE__, __func__, contigID1, contigID2);
								return FAILED;
							}
						}else // gapSize >= -10
						{ // update the overlap information

							if(gapSize<maxAdjustGapSizeThres)
							{ // -10 =< gapSize < 10
								if(overlapLenAlignment>0 && mismatchNum==0)
								{
									pContigOverlapInfo->mergeFlag = YES;
									pContigOverlapInfo->overlapLen = overlapLenAlignment;
								//}else if(overlapLenExact>0)
								}else if(overlapLenExact>=minExactOverlapThres)
								{
									pContigOverlapInfo->mergeFlag = YES;
									pContigOverlapInfo->overlapLen = overlapLenExact;
								}else
								{
									//if(gapSize<minBaseNumInGap) // default minBaseNumInGap: 2
									//{ // -10 =< gapSize < 2
									//	pContigOverlapInfo->gapSize = minBaseNumInGap;
									//}else
									//{ // 2 =< gapSize < 10
									//	pContigOverlapInfo->gapSize = gapSize;
									//}

									pContigOverlapInfo->gapSize = gapSize;
									pContigOverlapInfo->mergeFlag = NO;
								}
							}else
							{ // gapSize >= 10
								pContigOverlapInfo->mergeFlag = NO;
								pContigOverlapInfo->gapSize = gapSize;
							}

							if(updateContigEndsInfo(pContigOverlapInfo, contigGraph)==FAILED)
							{
								printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}
					}
					else
					{ // errors
						//printf("line=%d, In %s(), gapSize=%.2f, validPairedNum_gapEstimate=%d < %.2f, error!\n", __LINE__, __func__, gapSize, validPairedNum_gapEstimate, minLinksNumContigsThres);
						//return FAILED;
						pContigOverlapInfo->mergeFlag = NO;
						pContigOverlapInfo->gapSize = gapSize;
					}
				}
			}else // if(overlapLenAlignment<minOverlapThres || mismatchNum > mismatchThres)
			{ // there is a gap between the two contigs with high probability

				if(overlapLenAlignment==0)
				{
					overlapLenAlignment = strlen(alignResultArr[0]);
				}

				if(validPairedNum_gapEstimate>=minLinksNumContigsThres)
				{
					if(gapSize<minAdjustGapSizeThres)
					{
						// check the overlaps
						if((mismatchNum<=subMismatchThres || (double)mismatchNum/overlapLenAlignment<=maxMisMatchRatioThres) && (overlapLenAlignment<-gapSize+maxAdjustGapSizeThres && overlapLenAlignment>-gapSize-maxAdjustGapSizeThres))
						{
							overlapLenAdjust = overlapLenAlignment;
							// adjust overlapped sequences
							if(adjustOverlapSeq(overlapSeq1, overlapSeq2, alignResultArr, &overlapLenAdjust)==FAILED)
							{
								printf("line=%d, In %s(), cannot adjust the overlapped sequences, error!\n", __LINE__, __func__);
								return FAILED;
							}

							// update the contig sequences
							if(overlapLenAdjust>0)
							{ // valid overlap, then update contig end sequences
								if(updateContigs(contigGraph->contigItemArray+contigID1-1, contigGraph->contigItemArray+contigID2-1, orientation1, orientation2, overlapSeq1, seq_len1, overlapSeq2, seq_len2)==FAILED)
								{
									printf("line=%d, In %s(), cannot update contig end sequences, error!\n", __LINE__, __func__);
									return FAILED;
								}
								pContigOverlapInfo->mergeFlag = YES;
								pContigOverlapInfo->overlapLen = overlapLenAdjust;

								if(updateContigEndsInfo(pContigOverlapInfo, contigGraph)==FAILED)
								{
									printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
									return FAILED;
								}
							}else
							{
								pContigOverlapInfo->gapSize = gapSize;
								pContigOverlapInfo->mergeFlag = NO;

								// update the update length by cutting uncovered contig ends
								if(updateOverlapLenByCutUncoveredContigEnds(pContigOverlapInfo, contigGraph, readSet)==FAILED)
								{
									printf("line=%d, In %s(), cannot update the overlap size by cutting uncovered contig ends between contigs [%d,%d], error!\n", __LINE__, __func__, contigID1, contigID2);
									return FAILED;
								}
							}
						}else
						{ // the overlap length is not similar with the gap size, may be there are some errors at contig ends
							pContigOverlapInfo->gapSize = gapSize;
							pContigOverlapInfo->mergeFlag = NO;

							// update the update length by cutting uncovered contig ends
							if(updateOverlapLenByCutUncoveredContigEnds(pContigOverlapInfo, contigGraph, readSet)==FAILED)
							{
								printf("line=%d, In %s(), cannot update the overlap size by cutting uncovered contig ends between contigs [%d,%d], error!\n", __LINE__, __func__, contigID1, contigID2);
								return FAILED;
							}
						}
					}else // gapSize >= -10
					{ // update the overlap information
						if(gapSize<maxAdjustGapSizeThres)
						{ // -10 =< gapSize < 10
							if(overlapLenAlignment>0 && mismatchNum==0)
							{
								pContigOverlapInfo->mergeFlag = YES;
								pContigOverlapInfo->overlapLen = overlapLenAlignment;
							//}else if(overlapLenExact>0)
							}else if(overlapLenExact>=minExactOverlapThres)
							{
								pContigOverlapInfo->mergeFlag = YES;
								pContigOverlapInfo->overlapLen = overlapLenExact;
							}else
							{
								//if(gapSize<minBaseNumInGap) // default minBaseNumInGap: 2
								//{ // -10 =< gapSize < 2
								//	pContigOverlapInfo->gapSize = minBaseNumInGap;
								//}else
								//{ // 2 =< gapSize < 10
								//	pContigOverlapInfo->gapSize = gapSize;
								//}

								pContigOverlapInfo->gapSize = gapSize;
								pContigOverlapInfo->mergeFlag = NO;
							}
						}else
						{ // gapSize >= 10
							pContigOverlapInfo->mergeFlag = NO;
							pContigOverlapInfo->gapSize = gapSize;
						}

						if(updateContigEndsInfo(pContigOverlapInfo, contigGraph)==FAILED)
						{
							printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
				}else
				{ // errors
					//printf("line=%d, In %s(), gapSize=%.2f, validPairedNum_gapEstimate=%d < %.2f, error!\n", __LINE__, __func__, gapSize, validPairedNum_gapEstimate, minLinksNumContigsThres);
					//return FAILED;
					pContigOverlapInfo->mergeFlag = NO;
					pContigOverlapInfo->gapSize = gapSize;
				}
			}
		}
	}else
	{
		pContigOverlapInfo->mergeFlag = NO;
		pContigOverlapInfo->gapSize = gapSize;

		if(updateContigEndsInfo(pContigOverlapInfo, contigGraph)==FAILED)
		{
			printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Reverse the sequence.
 *  @return:
 *   If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short reverseSeq(char *seq, int seq_len)
{
	int i;
	char ch;

	for(i=0; i<seq_len/2; i++)
	{
		ch = seq[i];
		seq[i] = seq[seq_len-i-1];
		seq[seq_len-i-1] = ch;
	}

	for(i=0; i<seq_len; i++)
	{
		switch(seq[i])
		{
			case 'A': seq[i] = 'T'; break;
			case 'C': seq[i] = 'G'; break;
			case 'G': seq[i] = 'C'; break;
			case 'T': seq[i] = 'A'; break;
			default: printf("line=%d, In %s(), error sequence [ %s ], error!\n", __LINE__, __func__, seq);
					return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Estimate gap size between contigs.
 *  Note:
 *  	(endCutRound!=-1 and cutOrderArray!=NULL and uncoveredEndLenArray!=NULL) means the cuts operation will be done.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short gapSizeEstimateBetweenContigs(double *gapSize, int32_t *validPairedNum, contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph, readSet_t *readSet, int32_t endCutRound, int32_t *cutOrderArray, int32_t *uncoveredEndLenArray)
{
	int32_t i;
	int32_t contigID1, contigID2, contigOrient1, contigOrient2, contigEndFlag1, contigEndFlag2;
	int32_t readOrient1, readOrient2, contigPos1, contigPos2;
	int64_t readID1, readID2;
	int32_t baseNum1, baseNum2, contigLen1, contigLen2, seqLen1, seqLen2, tmpGapLen;
	int32_t tmp_contigID2;
	double totalGapSize;
	contigRead_t *contigReadArray1, *contigReadArray2;
	int32_t contigReadNum1, contigReadNum2;
	readMatchInfo_t *pReadMatchInfo;
	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;
	readMatchInfoBlock_t *readMatchInfoBlockArr;

	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

	contigID1 = pContigOverlapInfo->contigID1;
	contigID2 = pContigOverlapInfo->contigID2;

	contigOrient1 = pContigOverlapInfo->orientation1;
	contigOrient2 = pContigOverlapInfo->orientation2;

	// get the contig end flag
	if(contigGraph->contigItemArray[contigID1-1].onlyEnd5==YES)
		contigEndFlag1 = 2;
	else
	{
		if(contigOrient1==ORIENTATION_PLUS)
			contigEndFlag1 = 1;
		else
			contigEndFlag1 = 0;
	}

	if(contigGraph->contigItemArray[contigID2-1].onlyEnd5==YES)
		contigEndFlag2 = 2;
	else
	{
		if(contigOrient2==ORIENTATION_PLUS)
			contigEndFlag2 = 0;
		else
			contigEndFlag2 = 1;
	}

	// get the number of reads at the ends of the two contigs
	if(contigEndFlag1==1)
	{ // 3' end
		contigReadArray1 = contigGraph->contigItemArray[contigID1-1].contigReadArrayEnd3;
		contigReadNum1 = contigGraph->contigItemArray[contigID1-1].contigReadNumEnd3;
	}else
	{ // 5' end
		contigReadArray1 = contigGraph->contigItemArray[contigID1-1].contigReadArrayEnd5;
		contigReadNum1 = contigGraph->contigItemArray[contigID1-1].contigReadNumEnd5;
	}

	if(contigEndFlag2==1)
	{ // 3' end
		contigReadArray2 = contigGraph->contigItemArray[contigID2-1].contigReadArrayEnd3;
		contigReadNum2 = contigGraph->contigItemArray[contigID2-1].contigReadNumEnd3;
	}else
	{ // 5' end
		contigReadArray2 = contigGraph->contigItemArray[contigID2-1].contigReadArrayEnd5;
		contigReadNum2 = contigGraph->contigItemArray[contigID2-1].contigReadNumEnd5;
	}

	contigLen1 = contigGraph->contigItemArray[contigID1-1].contigLen;
	contigLen2 = contigGraph->contigItemArray[contigID2-1].contigLen;

	// estimate the gap size between contigs
	totalGapSize = 0;
	*validPairedNum = 0;
	for(i=0; i<contigReadNum1; i++)
	{
		readID1 = contigReadArray1[i].readID;
		readOrient1 = contigReadArray1[i].orientation;
		contigPos1 = contigReadArray1[i].contigPos;
		seqLen1 = contigReadArray1[i].seqlen;

		// get its paired end read
		if(readID1%2==1)
		{ // odd number
			readID2 = readID1 + 1;
		}else
		{ // even number
			readID2 = readID1 - 1;
		}

		readMatchInfoBlockID = (readID2 - 1) / maxItemNumPerReadMatchInfoBlock;
		rowNumInReadMatchInfoBlock = (readID2 - 1) % maxItemNumPerReadMatchInfoBlock;
		pReadMatchInfo = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;

		tmp_contigID2 = pReadMatchInfo->contigID;
		readOrient2 = pReadMatchInfo->readOrientation;
		contigPos2 = pReadMatchInfo->contigPos;
		seqLen2 = pReadMatchInfo->seqlen;

		if(tmp_contigID2==contigID2)
		{
			// ########################## Debug information ########################
			//if(contigID1==635)
			//{
			//	printf("contigID1=%d, contigID2=%d, contigPos1=%d, contigPos2=%d\n", contigID1, contigID2, contigPos1, contigPos2);
			//}
			// ########################## Debug information ########################

			// estimate the gap size according to a single read pair
			if((contigEndFlag1!=0 && contigOrient1==ORIENTATION_PLUS && readOrient1==ORIENTATION_PLUS) && (contigEndFlag2!=1 && contigOrient2==ORIENTATION_PLUS && readOrient2==ORIENTATION_MINUS))
			{ // (3', +, +) , (5', +, -)
				baseNum1 = contigLen1 - contigPos1 + 1;
				baseNum2 = contigPos2 + seqLen2 - 1;
				tmpGapLen = meanSizeInsert - baseNum1 - baseNum2;

				totalGapSize += tmpGapLen;
				(*validPairedNum) ++;

			}else if((contigEndFlag1!=0 && contigOrient1==ORIENTATION_PLUS && readOrient1==ORIENTATION_PLUS) && (contigEndFlag2!=0 && contigOrient2==ORIENTATION_MINUS && readOrient2==ORIENTATION_PLUS))
			{ // (3', +, +) , (3', -, +)
				baseNum1 = contigLen1 - contigPos1 + 1;
				baseNum2 = contigLen2 - contigPos2 + 1;
				tmpGapLen = meanSizeInsert - baseNum1 - baseNum2;

				totalGapSize += tmpGapLen;
				(*validPairedNum) ++;

			}else if((contigEndFlag1!=1 && contigOrient1==ORIENTATION_MINUS && readOrient1==ORIENTATION_MINUS) && (contigEndFlag2!=0 && contigOrient2==ORIENTATION_MINUS && readOrient2==ORIENTATION_PLUS))
			{ // (5', -, -) , (3', -, +)
				baseNum1 = contigPos1 + seqLen1 - 1;
				baseNum2 = contigLen2 - contigPos2 + 1;
				tmpGapLen = meanSizeInsert - baseNum1 - baseNum2;

				totalGapSize += tmpGapLen;
				(*validPairedNum) ++;

			}else if((contigEndFlag1!=1 && contigOrient1==ORIENTATION_MINUS && readOrient1==ORIENTATION_MINUS) && (contigEndFlag2!=1 && contigOrient2==ORIENTATION_PLUS && readOrient2==ORIENTATION_MINUS))
			{ // (5', -, -) , (5', +, -)
				baseNum1 = contigPos1 + seqLen1 - 1;
				baseNum2 = contigPos2 + seqLen2 - 1;
				tmpGapLen = meanSizeInsert - baseNum1 - baseNum2;

				totalGapSize += tmpGapLen;
				(*validPairedNum) ++;

			}else
			{
#if DEBUG_OUT_FLAG
				printf("line=%d, In %s(), (contigEndFlag1=%d, contigOrient1=%d, readOrient1=%d), (contigEndFlag2=%d, contigOrient2=%d, readOrient2=%d), readID1=%ld, readID2=%ld\n", __LINE__, __func__, contigEndFlag1, contigOrient1, readOrient1, contigEndFlag2, contigOrient2, readOrient2, readID1, readID2);
#endif
			}

			//printf("tmpGapLen=%d\n", tmpGapLen);
		}
	}

	if(*validPairedNum>0)
		totalGapSize /= (*validPairedNum);
	else
		totalGapSize = 0;

	// cut the contig ends
	if(endCutRound!=-1 && cutOrderArray!=NULL && uncoveredEndLenArray!=NULL)
	{
		if(endCutRound==0 || endCutRound==1)
		{ // the first round or the second round
			if(cutOrderArray[endCutRound]==0)
			{ // cut the end of the first contig
				totalGapSize += uncoveredEndLenArray[0];
			}else
			{ // cut the end of the second contig
				totalGapSize += uncoveredEndLenArray[1];
			}
		}else
		{ // the third round
			totalGapSize += uncoveredEndLenArray[0];
			totalGapSize += uncoveredEndLenArray[1];
		}
	}

	*gapSize = totalGapSize;

	return SUCCESSFUL;
}

/**
 * Update the used time of contig ends.
 */
short updateContigEndsInfo(contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph)
{
	int32_t contigID1, contigID2, contigOrient1, contigOrient2;

	contigID1 = pContigOverlapInfo->contigID1;
	contigID2 = pContigOverlapInfo->contigID2;
	contigOrient1 = pContigOverlapInfo->orientation1;
	contigOrient2 = pContigOverlapInfo->orientation2;

	if(pContigOverlapInfo->breakFlag==NO)
	{
		if(contigOrient1==ORIENTATION_PLUS)
		{ // plus orientation, then update the 3' end of the first contig
			contigGraph->contigItemArray[contigID1-1].usedTimeEnd3 ++;
			if(pContigOverlapInfo->mergeFlag==YES && pContigOverlapInfo->overlapLen>contigGraph->contigItemArray[contigID1-1].maxOverlapLenEnd3)
				contigGraph->contigItemArray[contigID1-1].maxOverlapLenEnd3 = pContigOverlapInfo->overlapLen;
		}else
		{ // minus orientation, then update the 5' end of the first contig
			contigGraph->contigItemArray[contigID1-1].usedTimeEnd5 ++;
			if(pContigOverlapInfo->mergeFlag==YES && pContigOverlapInfo->overlapLen>contigGraph->contigItemArray[contigID1-1].maxOverlapLenEnd5)
				contigGraph->contigItemArray[contigID1-1].maxOverlapLenEnd5 = pContigOverlapInfo->overlapLen;
		}

		if(contigOrient2==ORIENTATION_PLUS)
		{ // plus orientation, then update the 5' end of the second contig
			contigGraph->contigItemArray[contigID1-1].usedTimeEnd5 ++;
			if(pContigOverlapInfo->mergeFlag==YES && pContigOverlapInfo->overlapLen>contigGraph->contigItemArray[contigID1-1].maxOverlapLenEnd5)
				contigGraph->contigItemArray[contigID1-1].maxOverlapLenEnd5 = pContigOverlapInfo->overlapLen;
		}else
		{ // minus orientation, then update the 3' end of the second contig
			contigGraph->contigItemArray[contigID1-1].usedTimeEnd3 ++;
			if(pContigOverlapInfo->mergeFlag==YES && pContigOverlapInfo->overlapLen>contigGraph->contigItemArray[contigID1-1].maxOverlapLenEnd3)
				contigGraph->contigItemArray[contigID1-1].maxOverlapLenEnd3 = pContigOverlapInfo->overlapLen;
		}
	}else
	{
		printf("line=%d, In %s(), cannot update cotnig end info for breaked links of contigs [%d, %d], error!\n", __LINE__, __func__, contigID1, contigID2);
		return FAILED;
	}

	return SUCCESSFUL;
}


/**
 * Compute sequence overlap length.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short computeSeqOverlapLenExact(int32_t *overlapLen, const char *seq1, const int32_t seqLen1, const char *seq2, const int32_t seqLen2, int32_t *scoreArray, double gapSize)
{
	int32_t i, j, k;
	int32_t rowsNum, colsNum, minDiff;
	int32_t tmpOverlap;

	// ####################### Debug information ########################
#if DEBUG_SCAF_OVERLAP_FLAG
	printf("seq1=%s, len=%d\nseq2=%s, len=%d\n", seq1, (int32_t)strlen(seq1), seq2, (int32_t)strlen(seq2));
#endif
	// ####################### Debug information ########################

	// reset value of each item in the score array to zero
	if(memset(scoreArray, 0L, (maxOverlapSeqLen+1)*(maxOverlapSeqLen+1)*sizeof(int32_t))==NULL)
	{
		printf("line=%d, In %s(), cannot reset the memory to zero, error!\n", __LINE__, __func__);
		return FAILED;
	}

	rowsNum = seqLen1;
	colsNum = seqLen2;

	// set the elements
	for(i=0; i<rowsNum; i++)
	{
		for(j=0; j<colsNum; j++)
		{
			if(seq1[i]==seq2[j])
				scoreArray[i*rowsNum+j] = 1;
		}
	}

	// check the diagonal entries
	minDiff = INT_MAX;
	*overlapLen = 0;
	for(k=0; k<rowsNum; k++)
	{
		i = k;
		j = 0;
		tmpOverlap = 0;
		while(i<rowsNum && j<colsNum)
		{
			if(scoreArray[i*rowsNum+j]==0)
			{
				tmpOverlap = 0;
				break;
			}
			i ++;
			j ++;
			tmpOverlap ++;
		}

		if(tmpOverlap>0)
		{
			// #################### Debug information ####################
#if DEBUG_SCAF_OVERLAP_FLAG
			printf("Exact overlapLen=%d\n", tmpOverlap);
#endif
			// #################### Debug information ####################

			if(abs(tmpOverlap+gapSize)<minDiff)
			{
				minDiff = abs(tmpOverlap + gapSize);
				*overlapLen = tmpOverlap;
			}else
			{
				break;
			}
		}
	}

	//*overlapLen = tmpOverlap;

	//printf("Exact overlapLen=%d\n", *overlapLen);


/*
	//############################ Debug information #########################
	if(*overlapLen==0)
	{
		// output the score array to file
		if(outputScoreArrToFile(scoreArray, rowsNum, colsNum)==FAILED)
		{
			printf("In computeSeqOverlapLen(), cannot output the score array, error!\n");
			return FAILED;
		}
	}
	//############################ Debug information #########################
*/

#if DEBUG_SCAF_OVERLAP_FLAG
	printf("Final exact overlapLen=%d\n", *overlapLen);
#endif

	return SUCCESSFUL;
}

/**
 * Compute sequence overlap length.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short computeSeqOverlapLenByAlignment(const char *seq1, const int32_t seqLen1, const char *seq2, const int32_t seqLen2, int32_t *scoreArray, char **pAlignResultArray, int32_t *overlapLen, int32_t *mismatchNum)
{
	int32_t i, j, tmp;
	int32_t rowsNum, colsNum;
	int32_t maxValue, scoreIJ, maxRow, maxCol;
	int32_t itemNumInAlignArray;

	// ####################### Debug information ########################
#if DEBUG_SCAF_OVERLAP_FLAG
	printf("seq1=%s, len=%d\nseq2=%s, len=%d\n", seq1, (int)strlen(seq1), seq2, (int)strlen(seq2));
#endif
	// ####################### Debug information ########################

	// reset value of each item in the score array to zero
	if(memset(scoreArray, 0L, (maxOverlapSeqLen+1)*(maxOverlapSeqLen+1)*sizeof(int))==NULL)
	{
		printf("line=%d, In %s(), cannot reset the memory to zero, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<3; i++)
	{
		if(memset(pAlignResultArray[i], 0L, 2*(maxOverlapSeqLen+1)*sizeof(char))==NULL)
		{
			printf("line=%d, In %s(), cannot reset the memory to zero, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	rowsNum = seqLen1 + 1;
	colsNum = seqLen2 + 1;

	// initial values
	for(i=0; i<colsNum; i++)
		scoreArray[i] = 0;
	for(i=1; i<rowsNum; i++)
		scoreArray[colsNum*i] = 0;

	// compute the scores of each element
	for(i=1; i<rowsNum; i++)
	{
		for(j=1; j<colsNum; j++)
		{
			if(seq1[i-1]==seq2[j-1])
				scoreIJ = matchScore;
			else
				scoreIJ = mismatchScore;

			maxValue = INT_MIN;
			// compute the maximal score
			if(scoreArray[(i-1)*colsNum+j-1]+scoreIJ>maxValue)
			{ // from (i-1, j-1)
				maxValue = scoreArray[(i-1)*colsNum+j-1]+scoreIJ;
			}
			if(scoreArray[(i-1)*colsNum+j]+gapScore>maxValue)
			{ // from (i-1, j)
				maxValue = scoreArray[(i-1)*colsNum+j]+gapScore;
			}
			if(scoreArray[i*colsNum+j-1]+gapScore>maxValue)
			{ // from (i, j-1)
				maxValue = scoreArray[i*colsNum+j-1]+gapScore;
			}

			scoreArray[i*colsNum+j] = maxValue;
		}
	}

	// get the column of the maximal element in last row
	maxValue = 0;
	maxRow = 0;
	maxCol = 0;
	for(i=rowsNum-1; i<rowsNum; i++)
	{
		for(j=0; j<colsNum; j++)
		{
			if(scoreArray[i*colsNum+j]>maxValue)
			{
				maxValue = scoreArray[i*colsNum+j];
				maxRow = i;
				maxCol = j;
			}
		}
	}

#if DEBUG_SCAF_OVERLAP_FLAG
	printf("maxRow=%d, maxCol=%d, maxValue=%d\n", maxRow, maxCol, maxValue);
#endif

	itemNumInAlignArray = 0;
	*mismatchNum = 0;
	i = maxRow;
	j = maxCol;
	while(i>0 && j>0)
	{
		if(seq1[i-1]==seq2[j-1])
			scoreIJ = matchScore;
		else
			scoreIJ = mismatchScore;

		if(scoreArray[(i-1)*colsNum+j-1]+scoreIJ==scoreArray[i*colsNum+j])
		{ // from (i-1, j-1)
			if(seq1[i-1]!=seq2[j-1])
			{
				(*mismatchNum) ++;
				//printf("0:(%d,%d)\n", i, j);

				pAlignResultArray[0][itemNumInAlignArray] = seq1[i-1];
				pAlignResultArray[1][itemNumInAlignArray] = ' ';
				pAlignResultArray[2][itemNumInAlignArray] = seq2[j-1];
				itemNumInAlignArray ++;
			}else
			{
				pAlignResultArray[0][itemNumInAlignArray] = seq1[i-1];
				pAlignResultArray[1][itemNumInAlignArray] = '|';
				pAlignResultArray[2][itemNumInAlignArray] = seq2[j-1];
				itemNumInAlignArray ++;
			}

			i --;
			j --;
		}else if(scoreArray[(i-1)*colsNum+j]+gapScore==scoreArray[i*colsNum+j])
		{ // from (i-1, j)
			(*mismatchNum) ++;
			//printf("1:(%d,%d)\n", i, j);

			pAlignResultArray[0][itemNumInAlignArray] = seq1[i-1];
			pAlignResultArray[1][itemNumInAlignArray] = ' ';
			pAlignResultArray[2][itemNumInAlignArray] = '-';
			itemNumInAlignArray ++;

			i --;
		}else
		{ // from (i, j-1)
			(*mismatchNum) ++;
			//printf("2:(%d,%d)\n", i, j);

			pAlignResultArray[0][itemNumInAlignArray] = '-';
			pAlignResultArray[1][itemNumInAlignArray] = ' ';
			pAlignResultArray[2][itemNumInAlignArray] = seq2[j-1];
			itemNumInAlignArray ++;

			j --;
		}
	}

	*overlapLen = itemNumInAlignArray;

#if DEBUG_SCAF_OVERLAP_FLAG
	printf("i=%d, j=%d, mismatchNum=%d\n", i, j, *mismatchNum);
	printf("maxRow-i=%d, maxCol-j=%d, itemNumInAlignArray=%d\n", maxRow-i, maxCol-j, itemNumInAlignArray);
	printf("Alignment overlapLen=%d\n", *overlapLen);
#endif

	//if(*mismatchNum > 2 || (i==0 && j>0))
	if(i==0 && j>0)
	{
		*overlapLen = 0;
	}

#if DEBUG_SCAF_OVERLAP_FLAG
	printf("Final alignment overlapLen=%d\n", *overlapLen);
#endif

	// recover the alignment result
	for(i=0; i<3; i++)
	{
		for(j=0; j<itemNumInAlignArray/2; j++)
		{
			tmp = pAlignResultArray[i][j];
			pAlignResultArray[i][j] = pAlignResultArray[i][itemNumInAlignArray-1-j];
			pAlignResultArray[i][itemNumInAlignArray-1-j] = tmp;
		}
		pAlignResultArray[i][itemNumInAlignArray] = '\0';
	}

#if DEBUG_SCAF_OVERLAP_FLAG
	// print the alignment result
	for(i=0 ;i<3; i++)
	{
		printf("%s\n", pAlignResultArray[i]);
	}
#endif

	return SUCCESSFUL;
}

/**
 * Adjust overlapped sequences.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short adjustOverlapSeq(char *seq1, char *seq2, char **pAlignResultArray, int *overlapLen)
{
	double endFraction;
	int startPos[3];  // start from 0
	int endNum, tmpBaseNum1, tmpBaseNum2;
	int i;
	char *pTmpEndSeq[3], *accordanceSeq, *tmpStr;
	int baseNumTmpEndSeq[3], accordanceSeqLen;

	//############################### Debug information ########################
#if DEBUG_SCAF_OVERLAP_FLAG
	printf("Alignment overlapLen=%d\n", *overlapLen);
#endif
	//############################### Debug information ########################

	if((*overlapLen)<3)
		return SUCCESSFUL;

	for(i=0; i<3; i++)
	{
		pTmpEndSeq[i] = (char *) calloc((*overlapLen)+1, sizeof(char));
		if(pTmpEndSeq[i]==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	accordanceSeq = (char *) calloc((*overlapLen)+1, sizeof(char));
	if(accordanceSeq==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	tmpStr = (char *) calloc(strlen(seq1)+strlen(seq2)+1+1, sizeof(char));
	if(tmpStr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}


	endFraction = 1.0/3;
	endNum = ceil((*overlapLen) * endFraction);

	startPos[0] = 0;
	startPos[1] = endNum;
	startPos[2] = (*overlapLen) - endNum;

#if DEBUG_SCAF_OVERLAP_FLAG
	printf("endNum=%d\n", endNum);
	printf("overlapLen=%d, pos1=%d, pos2=%d, pos3=%d\n", *overlapLen, startPos[0], startPos[1], startPos[2]);
#endif

	// check the middle part
	for(i=startPos[1]; i<startPos[2]; i++)
	{
		if(pAlignResultArray[1][i]==' ')
		{
#if DEBUG_SCAF_OVERLAP_FLAG
			printf("errors in middle part\n");
#endif

			(*overlapLen) = 0;

			// free the memory
			for(i=0; i<3; i++)
			{
				free(pTmpEndSeq[i]); pTmpEndSeq[i] = NULL;
			}

			free(accordanceSeq); accordanceSeq = NULL;
			free(tmpStr); tmpStr = NULL;

			return SUCCESSFUL;
		}
	}

	// begin adjust the sequence
	baseNumTmpEndSeq[0] = 0;
	if(startPos[1]>startPos[0])
	{
		for(i=startPos[0]; i<startPos[1]; i++)
		{ // replace the second sequence with the first sequence
			if(pAlignResultArray[0][i]!='-')
			{
				pTmpEndSeq[0][ baseNumTmpEndSeq[0] ] = pAlignResultArray[0][i];
				baseNumTmpEndSeq[0] ++;
			}
		}
		pTmpEndSeq[0][ baseNumTmpEndSeq[0] ] = '\0';
	}

	baseNumTmpEndSeq[1] = startPos[2] - startPos[1];
	if(baseNumTmpEndSeq[1]>0)
		strncpy(pTmpEndSeq[1], pAlignResultArray[0]+startPos[1], baseNumTmpEndSeq[1]);

	baseNumTmpEndSeq[2] = 0;
	if((*overlapLen)>startPos[2])
	{
		for(i=startPos[2]; i<(*overlapLen); i++)
		{ // replace the first sequence with the second sequence
			if(pAlignResultArray[2][i]!='-')
			{
				pTmpEndSeq[2][ baseNumTmpEndSeq[2] ] = pAlignResultArray[2][i];
				baseNumTmpEndSeq[2] ++;
			}
		}
		pTmpEndSeq[2][ baseNumTmpEndSeq[2] ] = '\0';
	}

#if DEBUG_SCAF_OVERLAP_FLAG
	for(i=0; i<3; i++)
	{
		printf("tmpEndSeq[%d]=%s, len=%d\n", i, pTmpEndSeq[i], baseNumTmpEndSeq[i]);
	}
#endif

	// concatenate the new sub sequences
	if(baseNumTmpEndSeq[0]>0)
		strcpy(accordanceSeq, pTmpEndSeq[0]);
	if(baseNumTmpEndSeq[1]>0)
		strcat(accordanceSeq, pTmpEndSeq[1]);
	if(baseNumTmpEndSeq[2]>0)
		strcat(accordanceSeq, pTmpEndSeq[2]);
	accordanceSeqLen = baseNumTmpEndSeq[0] + baseNumTmpEndSeq[1] + baseNumTmpEndSeq[2];

	if(accordanceSeqLen!=strlen(accordanceSeq))
	{
		printf("accordanceSeqLen=%d != strlen(accordanceSeq)=%d, error!\n", accordanceSeqLen, (int32_t)strlen(accordanceSeq));
		return FAILED;
	}

#if DEBUG_SCAF_OVERLAP_FLAG
	printf("accordanceSeq=%s, len=%d\n", accordanceSeq, accordanceSeqLen);
#endif

	// get the number of bases of the two sequences
	tmpBaseNum1 = 0;
	for(i=0; i<(*overlapLen); i++)
	{
		if(pAlignResultArray[0][i]!='-')
			tmpBaseNum1 ++;
	}

	tmpBaseNum2 = 0;
	for(i=0; i<(*overlapLen); i++)
	{
		if(pAlignResultArray[2][i]!='-')
			tmpBaseNum2 ++;
	}

	// update the overlapped sequences
	strcpy(seq1+strlen(seq1)-tmpBaseNum1, accordanceSeq);	// sequence 1
	strcpy(tmpStr, accordanceSeq);			// sequence 2
	strcat(tmpStr, seq2+tmpBaseNum2);
	strcpy(seq2, tmpStr);
	*overlapLen = accordanceSeqLen;  // update the overlap length

#if DEBUG_SCAF_OVERLAP_FLAG
	// print the two sequences
	printf("new seq1=%s, len=%d\nnew seq2=%s, len=%d\n", seq1, (int)strlen(seq1), seq2, (int)strlen(seq2));
#endif

	// free the memory
	for(i=0; i<3; i++)
	{
		free(pTmpEndSeq[i]); pTmpEndSeq[i] = NULL;
	}

	free(accordanceSeq); accordanceSeq = NULL;
	free(tmpStr); tmpStr = NULL;

	return SUCCESSFUL;
}

/**
 * Update the two contigs, including their sequence and the length of the sequence.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short updateContigs(contigGraphItem_t *pContigItem1, contigGraphItem_t *pContigItem2, const int contigOrient1, const int contigOrient2, char *seq1, const int originalSeqLen1, char *seq2, const int originalSeqLen2)
{
	int i;
	int seqLen1, seqLen2;
	char *newSeq;

	seqLen1 = strlen(seq1);
	seqLen2 = strlen(seq2);

	// get the contig end sequences
	if(contigOrient1==ORIENTATION_PLUS)
	{ // update the sequence at the 3' end of contig 1
		if(seqLen1>originalSeqLen1)
		{
			newSeq = (char *) malloc((pContigItem1->contigLen+1 + seqLen1-originalSeqLen1) * sizeof(char));
			if(newSeq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// copy the sequence
			strcpy(newSeq, pContigItem1->contigSeq);

			free(pContigItem1->contigSeq);
			pContigItem1->contigSeq = newSeq;
		}

		strcpy(pContigItem1->contigSeq+pContigItem1->contigLen-originalSeqLen1, seq1);
		pContigItem1->contigLen += seqLen1 - originalSeqLen1;

	}else
	{ // update the sequence at the 5' end of contig 1
		if(reverseSeq(seq1, seqLen1)==FAILED)
		{
			printf("line=%d, In %s(), cannot reverse the sequence [ %s ], error!\n", __LINE__, __func__, seq1);
			return FAILED;
		}

		if(seqLen1!=originalSeqLen1)
		{
			newSeq = (char *) malloc((pContigItem1->contigLen+1 + seqLen1-originalSeqLen1) * sizeof(char));
			if(newSeq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			strcpy(newSeq, seq1);
			strcat(newSeq, pContigItem1->contigSeq+originalSeqLen1);

			free(pContigItem1->contigSeq);
			pContigItem1->contigSeq = newSeq;

			pContigItem1->contigLen += seqLen1 - originalSeqLen1;
		}else
		{
			for(i=0; i<seqLen1; i++) pContigItem1->contigSeq[i] = seq1[i];
		}
	}

	if(contigOrient2==ORIENTATION_PLUS)
	{ // update the sequence at the 5' end of contig 2
		if(seqLen2!=originalSeqLen2)
		{
			newSeq = (char *) malloc((pContigItem2->contigLen+1 + seqLen2-originalSeqLen2) * sizeof(char));
			if(newSeq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			strcpy(newSeq, seq2);
			strcat(newSeq, pContigItem2->contigSeq+originalSeqLen2);

			free(pContigItem2->contigSeq);
			pContigItem2->contigSeq = newSeq;

			pContigItem2->contigLen += seqLen2 - originalSeqLen2;
		}else
		{
			for(i=0; i<seqLen2; i++) pContigItem2->contigSeq[i] = seq2[i];
		}

	}else
	{ // update the sequence at the 3' end of contig 2
		if(reverseSeq(seq2, seqLen2)==FAILED)
		{
			printf("line=%d, In %s(), cannot reverse the sequence [ %s ], error!\n", __LINE__, __func__, seq2);
			return FAILED;
		}

		if(seqLen2>originalSeqLen2)
		{
			newSeq = (char *) malloc((pContigItem2->contigLen+1 + seqLen2-originalSeqLen2) * sizeof(char));
			if(newSeq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// copy the sequence
			strcpy(newSeq, pContigItem2->contigSeq);

			free(pContigItem2->contigSeq);
			pContigItem2->contigSeq = newSeq;

		}

		strcpy(pContigItem2->contigSeq+pContigItem2->contigLen-originalSeqLen2, seq2);
		pContigItem2->contigLen += seqLen2 - originalSeqLen2;
	}

	return SUCCESSFUL;
}

/**
 * Update overlap length between contigs by cutting the contig ends.
 * Notes:
 *  	(1) The contig with shorter uncovered length will be first considered.
 *  	(2) Three round of cuts will be applied until the detection of a valid overlap.
 *  		a) The first round will consider the contigs with shorter uncovered length;
 *  		b) The second round will consider the contigs with longer uncovered length if the overlap of the first round is not detected;
 *  		c) The third round will both consider the two contigs if the overlap of the second round is not detected.
 *  	(3) if the overlap is detected, the contig information and the overlap information will be updated; otherwise, the linked contigs will be broken.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateOverlapLenByCutUncoveredContigEnds(contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph, readSet_t *readSet)
{
	int32_t contigID1, orientation1, contigLen1, contigID2, orientation2, contigLen2;
	int32_t seq_len1, seq_len2, overlapLenExact, overlapLenAlignment, overlapLenAdjust, mismatchNum;
	int32_t validPairedNum_gapEstimate;
	double gapSize;

	int32_t uncoveredEndLenArray[2], cutOrderArray[2];
	int32_t endCutRound, maxEndCutRoundNum;
	int32_t contigLenBeforeCut[2];

	// get the information of contg1 and contig2
	contigID1 = pContigOverlapInfo->contigID1;
	orientation1 = pContigOverlapInfo->orientation1;
	contigLen1 = contigGraph->contigItemArray[contigID1-1].contigLen;
	contigID2 = pContigOverlapInfo->contigID2;
	orientation2 = pContigOverlapInfo->orientation2;
	contigLen2 = contigGraph->contigItemArray[contigID2-1].contigLen;

	contigLenBeforeCut[0] = contigLen1;
	contigLenBeforeCut[1] = contigLen2;

	// get the uncovered region length at contig ends
	if(getUncoveredLenAtContigEnds(uncoveredEndLenArray, pContigOverlapInfo, contigGraph, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot check the uncovered length at the ends of contigs [%d,%d], error!\n", __LINE__, __func__, contigID1, contigID2);
		return FAILED;
	}

	if((uncoveredEndLenArray[0]==0 && uncoveredEndLenArray[1]==0) || (contigLen1-uncoveredEndLenArray[0]<=0 && contigLen2-uncoveredEndLenArray[1]<=0))
	{ // No uncovered regions at contig ends, then break the contig links, and return
		// break the links
#if (DEBUG_SCAF_OVERLAP_FLAG==YES)
		printf("line=%d, In %s(), broken.\n", __LINE__, __func__);
#endif
		//pContigOverlapInfo->breakFlag = YES;
		return SUCCESSFUL;
	}

	if((uncoveredEndLenArray[0]>0 && uncoveredEndLenArray[1]>0) && (contigLen1-uncoveredEndLenArray[0]>0 && contigLen2-uncoveredEndLenArray[1]>0))
	{ // the end cuts of the two contigs are both valid
		maxEndCutRoundNum = 3;

		if(uncoveredEndLenArray[0]<=uncoveredEndLenArray[1])
		{ // 0 < len1 <= len2
			cutOrderArray[0] = 0;
			cutOrderArray[1] = 1;
		}else
		{ // len1 >= len2 > 0
			cutOrderArray[0] = 1;
			cutOrderArray[1] = 0;
		}

	}else
	{
		maxEndCutRoundNum = 1;

		if(uncoveredEndLenArray[0]>0 && contigLen1-uncoveredEndLenArray[0]>0)
		{ // the end cut of the first contig is valid
			cutOrderArray[0] = 0;
		}else if(uncoveredEndLenArray[1]>0 && contigLen2-uncoveredEndLenArray[1]>0)
		{
			cutOrderArray[0] = 1;
		}else
		{
			printf("line=%d, In %s(), uncoveredEndLen1=%d, contigLen1=%d, uncoveredEndLen2=%d, contigLen2=%d, error!\n", __LINE__, __func__, uncoveredEndLenArray[0], contigLen1, uncoveredEndLenArray[1], contigLen2);
			return FAILED;
		}
	}

	// To do the cut and sequence overlap
	for(endCutRound=0; endCutRound<maxEndCutRoundNum; endCutRound++)
	{
		if(endCutRound==0 || endCutRound==1)
		{ // the first or the second round
			// get the sequences from the two contigs according to their orientations
			if(cutOrderArray[endCutRound]==0)
			{ // cut the end of the first contig
				if(contigLen1-uncoveredEndLenArray[0] >= maxOverlapSeqLen)
				{
					seq_len1 = maxOverlapSeqLen;
				}else
				{
					seq_len1 = contigLen1 - uncoveredEndLenArray[0];
				}

				if(orientation1==ORIENTATION_PLUS)
				{
					strncpy(overlapSeq1, contigGraph->contigItemArray[contigID1-1].contigSeq+contigLen1-uncoveredEndLenArray[0]-seq_len1, seq_len1);
					overlapSeq1[seq_len1] = '\0';
				}else
				{
					strncpy(overlapSeq1, contigGraph->contigItemArray[contigID1-1].contigSeq+uncoveredEndLenArray[0], seq_len1);
					overlapSeq1[seq_len1] = '\0';
					if(reverseSeq(overlapSeq1, seq_len1)==FAILED) // reverse the sequence
					{
						printf("line=%d, In %s(), cannot reverse the sequence [ %s ], error!\n", __LINE__, __func__, overlapSeq1);
						return FAILED;
					}
				}

				if(contigLen2>=maxOverlapSeqLen)
				{
					seq_len2 = maxOverlapSeqLen;
				}else
				{
					seq_len2 = contigLen2;
				}

				if(orientation2==ORIENTATION_PLUS)
				{
					strncpy(overlapSeq2, contigGraph->contigItemArray[contigID2-1].contigSeq, seq_len2);
					overlapSeq2[seq_len2] = '\0';
				}else
				{
					strncpy(overlapSeq2, contigGraph->contigItemArray[contigID2-1].contigSeq+contigLen2-seq_len2, seq_len2);
					overlapSeq2[seq_len2] = '\0';
					if(reverseSeq(overlapSeq2, seq_len2)==FAILED) // reverse the sequence
					{
						printf("line=%d, In %s(), cannot reverse the sequence [ %s ], error!\n", __LINE__, __func__, overlapSeq2);
						return FAILED;
					}
				}
			}else
			{ // cut the end of the second contig
				if(contigLen1>=maxOverlapSeqLen)
				{
					seq_len1 = maxOverlapSeqLen;
				}else
				{
					seq_len1 = contigLen1;
				}

				if(orientation1==ORIENTATION_PLUS)
				{
					strncpy(overlapSeq1, contigGraph->contigItemArray[contigID1-1].contigSeq+contigLen1-seq_len1, seq_len1);
					overlapSeq1[seq_len1] = '\0';
				}else
				{
					strncpy(overlapSeq1, contigGraph->contigItemArray[contigID1-1].contigSeq, seq_len1);
					overlapSeq1[seq_len1] = '\0';
					if(reverseSeq(overlapSeq1, seq_len1)==FAILED) // reverse the sequence
					{
						printf("line=%d, In %s(), cannot reverse the sequence [ %s ], error!\n", __LINE__, __func__, overlapSeq1);
						return FAILED;
					}
				}

				if(contigLen2-uncoveredEndLenArray[1] >= maxOverlapSeqLen)
				{
					seq_len2 = maxOverlapSeqLen;
				}else
				{
					seq_len2 = contigLen2 - uncoveredEndLenArray[1];
				}

				if(orientation2==ORIENTATION_PLUS)
				{
					strncpy(overlapSeq2, contigGraph->contigItemArray[contigID2-1].contigSeq+uncoveredEndLenArray[1], seq_len2);
					overlapSeq2[seq_len2] = '\0';
				}else
				{
					strncpy(overlapSeq2, contigGraph->contigItemArray[contigID2-1].contigSeq+contigLen2-uncoveredEndLenArray[1]-seq_len2, seq_len2);
					overlapSeq2[seq_len2] = '\0';
					if(reverseSeq(overlapSeq2, seq_len2)==FAILED) // reverse the sequence
					{
						printf("line=%d, In %s(), cannot reverse the sequence [ %s ], error!\n", __LINE__, __func__, overlapSeq2);
						return FAILED;
					}
				}
			}
		}else
		{ // the third round
			// get the sequences from the two contigs according to their orientations
			if(contigLen1-uncoveredEndLenArray[0] >= maxOverlapSeqLen)
			{
				seq_len1 = maxOverlapSeqLen;
			}else
			{
				seq_len1 = contigLen1 - uncoveredEndLenArray[0];
			}

			if(orientation1==ORIENTATION_PLUS)
			{
				strncpy(overlapSeq1, contigGraph->contigItemArray[contigID1-1].contigSeq+contigLen1-uncoveredEndLenArray[0]-seq_len1, seq_len1);
				overlapSeq1[seq_len1] = '\0';
			}else
			{
				strncpy(overlapSeq1, contigGraph->contigItemArray[contigID1-1].contigSeq+uncoveredEndLenArray[0], seq_len1);
				overlapSeq1[seq_len1] = '\0';
				if(reverseSeq(overlapSeq1, seq_len1)==FAILED) // reverse the sequence
				{
					printf("line=%d, In %s(), cannot reverse the sequence [ %s ], error!\n", __LINE__, __func__, overlapSeq1);
					return FAILED;
				}
			}

			if(contigLen2-uncoveredEndLenArray[1] >= maxOverlapSeqLen)
			{
				seq_len2 = maxOverlapSeqLen;
			}else
			{
				seq_len2 = contigLen2 - uncoveredEndLenArray[1];
			}

			if(orientation2==ORIENTATION_PLUS)
			{
				strncpy(overlapSeq2, contigGraph->contigItemArray[contigID2-1].contigSeq+uncoveredEndLenArray[1], seq_len2);
				overlapSeq2[seq_len2] = '\0';
			}else
			{
				strncpy(overlapSeq2, contigGraph->contigItemArray[contigID2-1].contigSeq+contigLen2-uncoveredEndLenArray[1]-seq_len2, seq_len2);
				overlapSeq2[seq_len2] = '\0';
				if(reverseSeq(overlapSeq2, seq_len2)==FAILED) // reverse the sequence
				{
					printf("line=%d, In %s(), cannot reverse the sequence [ %s ], error!\n", __LINE__, __func__, overlapSeq2);
					return FAILED;
				}
			}
		}

		//############################# Debug information #########################
#if DEBUG_SCAF_OVERLAP_FLAG
		printf("overlapSeq1=%s, len=%d\noverlapSeq2=%s, len=%d\n", overlapSeq1, (int)strlen(overlapSeq1), overlapSeq2, (int)strlen(overlapSeq2));
#endif
		//############################# Debug information #########################

		// estimate gap size between contigs
		if(gapSizeEstimateBetweenContigs(&gapSize, &validPairedNum_gapEstimate, pContigOverlapInfo, contigGraph, readSet, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
		{
			printf("line=%d, In %s(), cannot estimate the gap size between contigs [%d,%d], error!\n", __LINE__, __func__, contigID1, contigID2);
			return FAILED;
		}

		if(validPairedNum_gapEstimate<minLinksNumContigsThres)
		{
			continue;
		}

		//if(gapSize<=stardardDeviationInsert)
		if(gapSize<=gapSizeSdevFactorOverlap*standardDev)
		{
			// compute the overlap length by exact alignment
			if(computeSeqOverlapLenExact(&overlapLenExact, overlapSeq1, seq_len1, overlapSeq2, seq_len2, scoreArr, gapSize)==FAILED)
			{
				printf("line=%d, In %s(), cannot compute sequence overlap length, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// update the contig overlap information
			if(overlapLenExact>=minOverlapThres && (overlapLenExact>=-gapSize-standardDev && overlapLenExact<=-gapSize+standardDev))
			{ // update the overlap information, and break the iteration
				if(cutContigEnds(pContigOverlapInfo, contigGraph, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
				{
					printf("line=%d, In %s(), cannot cut the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
					return FAILED;
				}

				pContigOverlapInfo->mergeFlag = YES;
				pContigOverlapInfo->overlapLen = overlapLenExact;

				if(updateContigEndsInfo(pContigOverlapInfo, contigGraph)==FAILED)
				{
					printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
					return FAILED;
				}

				//break;
				return SUCCESSFUL;
			}else
			{ // alignment the two sequences
				// pairwise alignment
				if(computeSeqOverlapLenByAlignment(overlapSeq1, seq_len1, overlapSeq2, seq_len2, scoreArr, alignResultArr, &overlapLenAlignment, &mismatchNum)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute sequence overlap length, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(overlapLenAlignment>=minOverlapThres && mismatchNum <= mismatchThres)
				{ // there is an overlap between the two contigs
					overlapLenAdjust = overlapLenAlignment;
					// adjust overlapped sequences
					if(adjustOverlapSeq(overlapSeq1, overlapSeq2, alignResultArr, &overlapLenAdjust)==FAILED)
					{
						printf("line=%d, In %s(), cannot adjust the overlapped sequences, error!\n", __LINE__, __func__);
						return FAILED;
					}

					// update the contig sequences
					if(overlapLenAdjust>0)
					{ // valid overlap, then update contig end sequences

						//############################### Debug information #########################
#if DEBUG_SCAF_OVERLAP_FLAG
						printf("line=%d, In %s(), before updateContigs(), contigLen1=%d, contigLen2=%d\n", __LINE__, __func__, contigGraph->contigItemArray[contigID1-1].contigLen, contigGraph->contigItemArray[contigID2-1].contigLen);
#endif
						//############################### Debug information #########################

						if(cutContigEnds(pContigOverlapInfo, contigGraph, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
						{
							printf("line=%d, In %s(), cannot cut the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
							return FAILED;
						}

						if(updateContigs(contigGraph->contigItemArray+contigID1-1, contigGraph->contigItemArray+contigID2-1, orientation1, orientation2, overlapSeq1, seq_len1, overlapSeq2, seq_len2)==FAILED)
						{
							printf("line=%d, In %s(), cannot update contig end sequences, error!\n", __LINE__, __func__);
							return FAILED;
						}

						//############################### Debug information #########################
#if DEBUG_SCAF_OVERLAP_FLAG
						printf("line=%d, In %s(), after updateContigs(), contigLen1=%d, contigLen2=%d\n", __LINE__, __func__, contigGraph->contigItemArray[contigID1-1].contigLen, contigGraph->contigItemArray[contigID2-1].contigLen);
#endif
						//############################### Debug information #########################

						pContigOverlapInfo->mergeFlag = YES;
						pContigOverlapInfo->overlapLen = overlapLenAdjust;

						if(updateContigEndsInfo(pContigOverlapInfo, contigGraph)==FAILED)
						{
							printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
							return FAILED;
						}

						//break;
						return SUCCESSFUL;
					}else // if(overlapLenAdjust==0)
					{ // there is a gap between the two contigs with high probability
						if(validPairedNum_gapEstimate>=minLinksNumContigsThres)
						{
							if(gapSize<minAdjustGapSizeThres)
							{ // link errors, then break the links, need next cut round
								//cutContigs(contigID1, contigID2, pContigInfoArr, endCutRound, cutOrderArray, uncoveredEndLenArray);
								//pContigOverlapInfo->breakFlag = YES;
								//break;
								continue; // go to next cut round
							}else // gapSize >= -10
							{ // update the overlap information

								if(cutContigEnds(pContigOverlapInfo, contigGraph, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
								{
									printf("line=%d, In %s(), cannot cut the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
									return FAILED;
								}

								if(gapSize<maxAdjustGapSizeThres)
								{ // -10 =< gapSize < 10
									if(overlapLenAlignment>0 && mismatchNum==0)
									{
										pContigOverlapInfo->mergeFlag = YES;
										pContigOverlapInfo->overlapLen = overlapLenAlignment;
									//}else if(overlapLenExact>0)
									}else if(overlapLenExact>=minExactOverlapThres)
									{
										pContigOverlapInfo->mergeFlag = YES;
										pContigOverlapInfo->overlapLen = overlapLenExact;
									}else
									{
										if(gapSize<minBaseNumInGap) // default minBaseNumInGap: 2
										{ // -10 =< gapSize < 2
											pContigOverlapInfo->gapSize = minBaseNumInGap;
										}else
										{ // 2 =< gapSize < 10
											pContigOverlapInfo->gapSize = gapSize;
										}
										pContigOverlapInfo->mergeFlag = NO;
									}
								}else
								{ // gapSize >= 10
									pContigOverlapInfo->mergeFlag = NO;
									pContigOverlapInfo->gapSize = gapSize;
								}

								if(updateContigEndsInfo(pContigOverlapInfo, contigGraph)==FAILED)
								{
									printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
									return FAILED;
								}

								//break;
								return SUCCESSFUL;
							}
						}else
						{ // errors
							printf("line=%d, In %s(), validPairedNum_gapEstimate=%d < %.2f, error!\n", __LINE__, __func__, validPairedNum_gapEstimate, minLinksNumContigsThres);
							return FAILED;
							//pContigOverlapInfo->mergeFlag = NO;
							//pContigOverlapInfo->gapSize = gapSize;
						}
					}
				}else // if(overlapLenAlignment<minOverlapThres || mismatchNum > mismatchThres)
				{ // there is a gap between the two contigs with high probability

					if(overlapLenAlignment==0)
						overlapLenAlignment = strlen(alignResultArr[0]);

					if(validPairedNum_gapEstimate>=minLinksNumContigsThres)
					{
						if(gapSize<minAdjustGapSizeThres)
						{ // gapSize < -10
							// check the overlaps
							if((mismatchNum<=subMismatchThres || (double)mismatchNum/overlapLenAlignment<=maxMisMatchRatioThres) && (overlapLenAlignment<-gapSize+maxAdjustGapSizeThres && overlapLenAlignment>-gapSize-maxAdjustGapSizeThres))
							{
								overlapLenAdjust = overlapLenAlignment;
								// adjust overlapped sequences
								if(adjustOverlapSeq(overlapSeq1, overlapSeq2, alignResultArr, &overlapLenAdjust)==FAILED)
								{
									printf("line=%d, In %s(), cannot adjust the overlapped sequences, error!\n", __LINE__, __func__);
									return FAILED;
								}

								// update the contig sequences
								if(overlapLenAdjust>0)
								{ // valid overlap, then update contig end sequences

									if(cutContigEnds(pContigOverlapInfo, contigGraph, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
									{
										printf("line=%d, In %s(), cannot cut the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
										return FAILED;
									}

									if(updateContigs(contigGraph->contigItemArray+contigID1-1, contigGraph->contigItemArray+contigID2-1, orientation1, orientation2, overlapSeq1, seq_len1, overlapSeq2, seq_len2)==FAILED)
									{
										printf("line=%d, In %s(), cannot update contig end sequences, error!\n", __LINE__, __func__);
										return FAILED;
									}
									pContigOverlapInfo->mergeFlag = YES;
									pContigOverlapInfo->overlapLen = overlapLenAdjust;

									if(updateContigEndsInfo(pContigOverlapInfo, contigGraph)==FAILED)
									{
										printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
										return FAILED;
									}

									//break;
									return SUCCESSFUL;
								}else
								{ // link errors, then break the links
									//cutContigs(contigID1, contigID2, pContigInfoArr, endCutRound, cutOrderArray, uncoveredEndLenArray);
									//pContigOverlapInfo->breakFlag = YES;
									//return SUCCESSFUL;
									continue; // go to next cut round
								}
							}else
							{
								continue;
							}
						}else // gapSize >= -10
						{ // update the overlap information

							if(gapSize<maxAdjustGapSizeThres)
							{ // -10 =< gapSize < 10
								//if((contigGraph->contigItemArray[contigID1-1].onlyEnd5==YES || contigGraph->contigItemArray[contigID2-1].onlyEnd5==YES) && validPairedNum_gapEstimate<breakLinkNumThres)
								if((contigGraph->contigItemArray[contigID1-1].shortFlag==YES || contigGraph->contigItemArray[contigID2-1].shortFlag==YES) && validPairedNum_gapEstimate<breakLinkNumThres)
								{
#if (DEBUG_SCAF_OVERLAP_FLAG==YES)
									printf("line=%d, In %s(), broken.\n", __LINE__, __func__);
#endif
									pContigOverlapInfo->breakFlag = YES;
								}else
								{
									if(cutContigEnds(pContigOverlapInfo, contigGraph, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
									{
										printf("line=%d, In %s(), cannot cut the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
										return FAILED;
									}

									if(overlapLenAlignment>0 && mismatchNum==0)
									{
										pContigOverlapInfo->mergeFlag = YES;
										pContigOverlapInfo->overlapLen = overlapLenAlignment;
									//}else if(overlapLenExact>0)
									}else if(overlapLenExact>=minExactOverlapThres)
									{
										pContigOverlapInfo->mergeFlag = YES;
										pContigOverlapInfo->overlapLen = overlapLenExact;
									}else
									{
										if(gapSize<minBaseNumInGap) // default minBaseNumInGap: 2
										{ // -10 =< gapSize < 2
											pContigOverlapInfo->gapSize = minBaseNumInGap;
										}else
										{ // 2 =< gapSize < 10
											pContigOverlapInfo->gapSize = gapSize;
										}
										pContigOverlapInfo->mergeFlag = NO;
									}

									if(updateContigEndsInfo(pContigOverlapInfo, contigGraph)==FAILED)
									{
										printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
										return FAILED;
									}
								}

								return SUCCESSFUL;
							}else
							{ // gapSize >= 10

								pContigOverlapInfo->mergeFlag = NO;
								pContigOverlapInfo->gapSize = gapSize;

								if(cutContigEnds(pContigOverlapInfo, contigGraph, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
								{
									printf("line=%d, In %s(), cannot cut the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
									return FAILED;
								}

								if(updateContigEndsInfo(pContigOverlapInfo, contigGraph)==FAILED)
								{
									printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
									return FAILED;
								}

								return SUCCESSFUL;
							}
						}
					}else
					{ // errors
						printf("line=%d, In %s(), validPairedNum_gapEstimate=%d < %.2f, error!\n", __LINE__, __func__, validPairedNum_gapEstimate, minLinksNumContigsThres);
						return FAILED;
						//pContigOverlapInfo->mergeFlag = NO;
						//pContigOverlapInfo->gapSize = gapSize;
					}
				}
			}
		}else
		{
			pContigOverlapInfo->mergeFlag = NO;
			pContigOverlapInfo->gapSize = gapSize;

			if(cutContigEnds(pContigOverlapInfo, contigGraph, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
			{
				printf("line=%d, In %s(), cannot cut the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
				return FAILED;
			}

			if(updateContigEndsInfo(pContigOverlapInfo, contigGraph)==FAILED)
			{
				printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
				return FAILED;
			}

			return SUCCESSFUL;
		}
	} // end for(endCutRound=0; endCutRound<maxEndCutRoundNum; endCutRound++)

	// No valid overlaps, break the links without cutting contig ends
	//pContigOverlapInfo->breakFlag = YES;

	return SUCCESSFUL;
}

/**
 * Get the length of uncovered sub sequences at the ends of the two contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getUncoveredLenAtContigEnds(int32_t *pUncoveredEndLenArray, contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph, readSet_t *readSet)
{
	int32_t i, hitRow;
	int32_t contigID1, contigID2, contigOrient1, contigOrient2, seqLen1, seqLen2, contigEndFlag1, contigEndFlag2;
	int32_t readsNum1, readsNum2, firstRow1, firstRow2, readOrient1, readOrient2, contigPos1, contigPos2;
	contigRead_t *contigReadArray1, *contigReadArray2;
	int64_t readID1, readID2;
	int32_t contigLen1, contigLen2;
	int32_t tmp_contigID, validReadOrient1, validReadOrient2;
	int32_t cutAllowedContig1, cutAllowedContig2;
	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;
	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo;

	contigID1 = pContigOverlapInfo->contigID1;
	contigID2 = pContigOverlapInfo->contigID2;
	contigOrient1 = pContigOverlapInfo->orientation1;
	contigOrient2 = pContigOverlapInfo->orientation2;
	contigLen1 = contigGraph->contigItemArray[contigID1-1].contigLen;
	contigLen2 = contigGraph->contigItemArray[contigID2-1].contigLen;

	// get the contig end flag
	if(contigGraph->contigItemArray[contigID1-1].onlyEnd5==YES)
		contigEndFlag1 = 2;
	else
	{
		if(contigOrient1==ORIENTATION_PLUS)
			contigEndFlag1 = 1;
		else
			contigEndFlag1 = 0;
	}

	if(contigGraph->contigItemArray[contigID2-1].onlyEnd5==YES)
		contigEndFlag2 = 2;
	else
	{
		if(contigOrient2==ORIENTATION_PLUS)
			contigEndFlag2 = 0;
		else
			contigEndFlag2 = 1;
	}

	cutAllowedContig1 = cutAllowedContig2 = YES;
	if(contigOrient1==ORIENTATION_PLUS)
	{ // to cut the 3' end
		if(contigGraph->contigItemArray[contigID1-1].usedTimeEnd3>0)
			cutAllowedContig1 = NO;
	}else
	{ // to cut the 5' end
		if(contigGraph->contigItemArray[contigID1-1].usedTimeEnd5>0)
			cutAllowedContig1 = NO;
	}

	if(contigOrient2==ORIENTATION_PLUS)
	{ // to cut the 5' end
		if(contigGraph->contigItemArray[contigID2-1].usedTimeEnd5>0)
			cutAllowedContig2 = NO;
	}else
	{ // to cut the 3' end
		if(contigGraph->contigItemArray[contigID2-1].usedTimeEnd3>0)
			cutAllowedContig2 = NO;
	}

	if(contigEndFlag1==1)
	{ // 3' end
		contigReadArray1 = contigGraph->contigItemArray[contigID1-1].contigReadArrayEnd3;
		readsNum1 = contigGraph->contigItemArray[contigID1-1].contigReadNumEnd3;
	}else
	{ // 5' end
		contigReadArray1 = contigGraph->contigItemArray[contigID1-1].contigReadArrayEnd5;
		readsNum1 = contigGraph->contigItemArray[contigID1-1].contigReadNumEnd5;
	}

	if(contigEndFlag2==1)
	{ // 3' end
		contigReadArray2 = contigGraph->contigItemArray[contigID2-1].contigReadArrayEnd3;
		readsNum2 = contigGraph->contigItemArray[contigID2-1].contigReadNumEnd3;
	}else
	{ // 5' end
		contigReadArray2 = contigGraph->contigItemArray[contigID2-1].contigReadArrayEnd5;
		readsNum2 = contigGraph->contigItemArray[contigID2-1].contigReadNumEnd5;
	}

	// get the valid read orientations of the two paired reads in the contigs
	if(contigOrient1==ORIENTATION_PLUS)
		validReadOrient1 = ORIENTATION_PLUS;
	else
		validReadOrient1 = ORIENTATION_MINUS;

	if(contigOrient2==ORIENTATION_PLUS)
		validReadOrient2 = ORIENTATION_MINUS;
	else
		validReadOrient2 = ORIENTATION_PLUS;

	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

	// get the number of reads at the end of the first contig
	if(cutAllowedContig1==YES)
	{
		// get the uncovered length at the end of contigs
		pUncoveredEndLenArray[0] = 0;
		if(contigOrient1==ORIENTATION_PLUS)
		{ // (1, +): the plus orientation of the first contig
			for(i=readsNum1-1; i>=0; i--)
			{
				readID1 = contigReadArray1[i].readID;
				readOrient1 = contigReadArray1[i].orientation;
				contigPos1 = contigReadArray1[i].contigPos;
				seqLen1 = contigReadArray1[i].seqlen;

				if(readOrient1==validReadOrient1)
				{
					// get its paired end read
					if(readID1%2==1)
					{ // odd number
						readID2 = readID1 + 1;
					}else
					{ // even number
						readID2 = readID1 - 1;
					}

					readMatchInfoBlockID = (readID2 - 1) / maxItemNumPerReadMatchInfoBlock;
					rowNumInReadMatchInfoBlock = (readID2 - 1) % maxItemNumPerReadMatchInfoBlock;
					pReadMatchInfo = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;

					tmp_contigID = pReadMatchInfo->contigID;
					readOrient2 = pReadMatchInfo->readOrientation;

					if(tmp_contigID==contigID2 && readOrient2==validReadOrient2)
					{
						pUncoveredEndLenArray[0] = contigLen1 - contigPos1 - seqLen1 + 1;
						break;
					}
				}
			}

			if(i<0)
			{
				pUncoveredEndLenArray[0] = contigLen1;
			}
		}else
		{ // (1, -): the minus orientation of the first contig
			for(i=0; i<readsNum1; i++)
			{
				readID1 = contigReadArray1[i].readID;
				readOrient1 = contigReadArray1[i].orientation;
				contigPos1 = contigReadArray1[i].contigPos;
				seqLen1 = contigReadArray1[i].seqlen;

				if(readOrient1==validReadOrient1)
				{
					// get its paired end read
					if(readID1%2==1)
					{ // odd number
						readID2 = readID1 + 1;
					}else
					{ // even number
						readID2 = readID1 - 1;
					}

					readMatchInfoBlockID = (readID2 - 1) / maxItemNumPerReadMatchInfoBlock;
					rowNumInReadMatchInfoBlock = (readID2 - 1) % maxItemNumPerReadMatchInfoBlock;
					pReadMatchInfo = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;

					tmp_contigID = pReadMatchInfo->contigID;
					readOrient2 = pReadMatchInfo->readOrientation;

					if(tmp_contigID==contigID2 && readOrient2==validReadOrient2)
					{
						pUncoveredEndLenArray[0] = contigPos1 - 1;
						break;
					}
				}
			}

			if(i>=readsNum1)
			{
				pUncoveredEndLenArray[0] = contigLen1;
			}
		}
	}else
	{
		pUncoveredEndLenArray[0] = 0;
	}

	// get the number of reads at the end of the second contig
	if(cutAllowedContig2==YES)
	{
		pUncoveredEndLenArray[1] = 0;
		if(contigOrient2==ORIENTATION_PLUS)
		{ // (2, +): the plus orientation of the second contig
			for(i=0; i<readsNum2; i++)
			{
				readID2 = contigReadArray2[i].readID;
				readOrient2 = contigReadArray2[i].orientation;
				contigPos2 = contigReadArray2[i].contigPos;
				seqLen2 = contigReadArray2[i].seqlen;

				if(readOrient2==validReadOrient2)
				{
					// get its paired end read
					if(readID2%2==1)
					{ // odd number
						readID1 = readID2 + 1;
					}else
					{ // even number
						readID1 = readID2 - 1;
					}

					readMatchInfoBlockID = (readID1 - 1) / maxItemNumPerReadMatchInfoBlock;
					rowNumInReadMatchInfoBlock = (readID1 - 1) % maxItemNumPerReadMatchInfoBlock;
					pReadMatchInfo = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;

					tmp_contigID = pReadMatchInfo->contigID;
					readOrient1 = pReadMatchInfo->readOrientation;

					if(tmp_contigID==contigID1 && readOrient1==validReadOrient1)
					{
						pUncoveredEndLenArray[1] = contigPos2 - 1;
						break;
					}
				}
			}

			if(i>=readsNum2)
			{
				pUncoveredEndLenArray[1] = contigLen2;
			}
		}else
		{ // (2, -): the minus orientation of the second contig
			for(i=readsNum2-1; i>=0; i--)
			{
				readID2 = contigReadArray2[i].readID;
				readOrient2 = contigReadArray2[i].orientation;
				contigPos2 = contigReadArray2[i].contigPos;
				seqLen2 = contigReadArray2[i].seqlen;

				if(readOrient2==validReadOrient2)
				{
					// get its paired end read
					if(readID2%2==1)
					{ // odd number
						readID1 = readID2 + 1;
					}else
					{ // even number
						readID1 = readID2 - 1;
					}

					readMatchInfoBlockID = (readID1 - 1) / maxItemNumPerReadMatchInfoBlock;
					rowNumInReadMatchInfoBlock = (readID1 - 1) % maxItemNumPerReadMatchInfoBlock;
					pReadMatchInfo = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;

					tmp_contigID = pReadMatchInfo->contigID;
					readOrient1 = pReadMatchInfo->readOrientation;

					if(tmp_contigID==contigID1 && readOrient1==validReadOrient1)
					{
						pUncoveredEndLenArray[1] = contigLen2 - contigPos2 - seqLen2 + 1;
						break;
					}
				}
			}

			if(i<0)
			{
				pUncoveredEndLenArray[1] = contigLen2;
			}
		}
	}else
	{
		pUncoveredEndLenArray[1] = 0;
	}

	//########################### Debug information ##########################
#if DEBUG_OUT_FLAG
	printf("uncoveredEndLen1=%d, uncoveredEndLen2=%d\n", pUncoveredEndLenArray[0], pUncoveredEndLenArray[1]);
#endif
	//########################### Debug information ##########################

	// adjust the values of uncoveredEndLen1 and uncoveredEndLen2
	if(pUncoveredEndLenArray[0]<0)
		pUncoveredEndLenArray[0] = 0;

	if(pUncoveredEndLenArray[1]<=0)
		pUncoveredEndLenArray[1] = 0;

	//########################### Debug information ##########################
#if DEBUG_OUT_FLAG
	printf("Final uncoveredEndLen1=%d, uncoveredEndLen2=%d\n", pUncoveredEndLenArray[0], pUncoveredEndLenArray[1]);
#endif
	//########################### Debug information ##########################

	return SUCCESSFUL;
}

/**
 * Cut contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short cutContigEnds(contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph, int32_t endCutRound, int32_t *cutOrderArray, int32_t *uncoveredEndLenArray)
{
	int32_t contigID1, contigID2, contigOrient1, contigOrient2;
	int32_t contigLen1, contigLen2;
	char *tmpStr;

	contigID1 = pContigOverlapInfo->contigID1;
	contigID2 = pContigOverlapInfo->contigID2;
	contigOrient1 = pContigOverlapInfo->orientation1;
	contigOrient2 = pContigOverlapInfo->orientation2;

	contigLen1 = contigGraph->contigItemArray[contigID1-1].contigLen;
	contigLen2 = contigGraph->contigItemArray[contigID2-1].contigLen;

	// update contigs information
	if(endCutRound==0 || endCutRound==1)
	{ // the first round or the second round
		if(cutOrderArray[endCutRound]==0)
		{ // cut the end of the first contig
			if(contigOrient1==ORIENTATION_PLUS)
			{ // the plus orientation, cut the 3' end
				contigGraph->contigItemArray[contigID1-1].contigLen -= uncoveredEndLenArray[0];
				contigGraph->contigItemArray[contigID1-1].contigSeq[ contigGraph->contigItemArray[contigID1-1].contigLen ] = '\0';
			}else
			{ // the minus orientation, cut the 5' end
				tmpStr = (char *) calloc(contigLen1+1, sizeof(char));
				if(tmpStr==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// copy the remained bases to tmpStr
				strcpy(tmpStr, contigGraph->contigItemArray[contigID1-1].contigSeq+uncoveredEndLenArray[0]);
				// copy back
				strcpy(contigGraph->contigItemArray[contigID1-1].contigSeq, tmpStr);
				contigGraph->contigItemArray[contigID1-1].contigLen -= uncoveredEndLenArray[0];

				// free tmpStr
				free(tmpStr);
				tmpStr = NULL;
			}

			// ########################### Debug information ####################
#if DEBUG_FLAG
			if(strlen(contigGraph->contigItemArray[contigID1-1].contigSeq)!=contigGraph->contigItemArray[contigID1-1].contigLen)
			{
				printf("line=%d, In %s(), contigID1=%d, contigSeq=%s, strlen=%d, contigLen=%d, error!\n", __LINE__, __func__, contigID1, contigGraph->contigItemArray[contigID1-1].contigSeq, (int)strlen(contigGraph->contigItemArray[contigID1-1].contigSeq), contigGraph->contigItemArray[contigID1-1].contigLen);
				return FAILED;
			}
#endif
			// ########################### Debug information ####################

		}else
		{ // cut the end of the second contig
			if(contigOrient2==ORIENTATION_PLUS)
			{ // the plus orientation, cut the 5' end
				tmpStr = (char *) calloc(contigLen2+1, sizeof(char));
				if(tmpStr==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// copy the remained bases to tmpStr
				strcpy(tmpStr, contigGraph->contigItemArray[contigID2-1].contigSeq+uncoveredEndLenArray[1]);
				// copy back
				strcpy(contigGraph->contigItemArray[contigID2-1].contigSeq, tmpStr);
				contigGraph->contigItemArray[contigID2-1].contigLen -= uncoveredEndLenArray[1];

				// free tmpStr
				free(tmpStr);
				tmpStr = NULL;
			}else
			{ // the minus orientation, cut the 3' end
				contigGraph->contigItemArray[contigID2-1].contigLen -= uncoveredEndLenArray[1];
				contigGraph->contigItemArray[contigID2-1].contigSeq[ contigGraph->contigItemArray[contigID2-1].contigLen ] = '\0';
			}

			// ########################### Debug information ####################
#if DEBUG_FLAG
			if(strlen(contigGraph->contigItemArray[contigID2-1].contigSeq)!=contigGraph->contigItemArray[contigID2-1].contigLen)
			{
				printf("line=%d, In %s(), contigID2=%d, contigSeq=%s, strlen=%d, contigLen=%d, error!\n", __LINE__, __func__, contigID2, contigGraph->contigItemArray[contigID2-1].contigSeq, (int)strlen(contigGraph->contigItemArray[contigID2-1].contigSeq), contigGraph->contigItemArray[contigID2-1].contigLen);
				return FAILED;
			}
#endif
			// ########################### Debug information ####################
		}
	}else
	{ // the third round
		if(contigOrient1==ORIENTATION_PLUS && contigOrient2==ORIENTATION_MINUS)
		{ // A(+), B(-)
			contigGraph->contigItemArray[contigID1-1].contigLen -= uncoveredEndLenArray[0];
			contigGraph->contigItemArray[contigID1-1].contigSeq[ contigGraph->contigItemArray[contigID1-1].contigLen ] = '\0';

			contigGraph->contigItemArray[contigID2-1].contigLen -= uncoveredEndLenArray[1];
			contigGraph->contigItemArray[contigID2-1].contigSeq[ contigGraph->contigItemArray[contigID2-1].contigLen ] = '\0';
		}else
		{
			tmpStr = (char *) calloc(contigLen1+contigLen2+1, sizeof(char));
			if(tmpStr==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(contigOrient1==ORIENTATION_PLUS)
			{ // A(+)
				contigGraph->contigItemArray[contigID1-1].contigLen -= uncoveredEndLenArray[0];
				contigGraph->contigItemArray[contigID1-1].contigSeq[ contigGraph->contigItemArray[contigID1-1].contigLen ] = '\0';
			}else
			{ // A(-)
				// copy the remained bases to tmpStr
				strcpy(tmpStr, contigGraph->contigItemArray[contigID1-1].contigSeq+uncoveredEndLenArray[0]);
				// copy back
				strcpy(contigGraph->contigItemArray[contigID1-1].contigSeq, tmpStr);

				contigGraph->contigItemArray[contigID1-1].contigLen -= uncoveredEndLenArray[0];
			}

			if(contigOrient2==ORIENTATION_PLUS)
			{ // B(+)
				// copy the remained bases to tmpStr
				strcpy(tmpStr, contigGraph->contigItemArray[contigID2-1].contigSeq+uncoveredEndLenArray[1]);
				// copy back
				strcpy(contigGraph->contigItemArray[contigID2-1].contigSeq, tmpStr);

				contigGraph->contigItemArray[contigID2-1].contigLen -= uncoveredEndLenArray[1];
			}else
			{ // B(-)
				contigGraph->contigItemArray[contigID2-1].contigLen -= uncoveredEndLenArray[1];
				contigGraph->contigItemArray[contigID2-1].contigSeq[ contigGraph->contigItemArray[contigID2-1].contigLen ] = '\0';
			}

			free(tmpStr);
		}

		// ########################### Debug information ####################
#if DEBUG_FLAG
		if(strlen(contigGraph->contigItemArray[contigID1-1].contigSeq)!=contigGraph->contigItemArray[contigID1-1].contigLen)
		{
			printf("line=%d, In %s(), contigID1=%d, contigSeq=%s, strlen=%d, contigLen=%d, error!\n", __LINE__, __func__, contigID1, contigGraph->contigItemArray[contigID1-1].contigSeq, (int)strlen(contigGraph->contigItemArray[contigID1-1].contigSeq), contigGraph->contigItemArray[contigID1-1].contigLen);
			return FAILED;
		}

		if(strlen(contigGraph->contigItemArray[contigID2-1].contigSeq)!=contigGraph->contigItemArray[contigID2-1].contigLen)
		{
			printf("line=%d, In %s(), contigID2=%d, contigSeq=%s, strlen=%d, contigLen=%d, error!\n", __LINE__, __func__, contigID2, contigGraph->contigItemArray[contigID2-1].contigSeq, (int)strlen(contigGraph->contigItemArray[contigID2-1].contigSeq), contigGraph->contigItemArray[contigID2-1].contigLen);
			return FAILED;
		}
#endif
		// ########################### Debug information ####################
	}

	return SUCCESSFUL;
}

/**
 * Confirm the link information of overlapped contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short confirmOverlapInfoInScaf(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph, readSet_t *readSet)
{
	int32_t i, j;
	scaffoldItem_t *scaffoldItem;
	contigOverlap_t *contigOverlapArray, *pContigOverlapInfo1, *pContigOverlapInfo2;
	int32_t linkedContigNum, rowsNum, startRow;
	//int32_t subRowsNum, subStartRow, subEndRow, subLinkedContigNum;
	int32_t contigID1, contigID2, contigID3, contigLen1, contigLen2, contigLen3, overlapLen1, overlapLen2, pairedNum_1_3;
	int32_t endFlag1, endFlag2, endFlag3;
	double gapSize1, gapSize2;
	int32_t shortRowsNum1, shortRowsNum2;
	int32_t minPairNumThreshold, validPairNum, invalidPairNumAtEnds;

	scaffoldItem = scaffoldSet->scaffoldItemList;
	while(scaffoldItem)
	{
		contigOverlapArray = scaffoldItem->contigOverlapArray;
		rowsNum = scaffoldItem->itemNumContigOverlapArray;
		linkedContigNum = scaffoldItem->linkedContigsNum;

//		printf("scaffoldID=%d\n", scaffoldItem->scaffoldID);
//		if(scaffoldItem->scaffoldID==73 || scaffoldItem->scaffoldID==318)
//		{
//			printf("scaffoldID=%d\n", scaffoldItem->scaffoldID);
//		}

		if(linkedContigNum>=3)
		{
			startRow = 0;
			for(i=0; i<rowsNum-1; i++)
			{
				pContigOverlapInfo1 = contigOverlapArray + i;
				pContigOverlapInfo2 = contigOverlapArray + i + 1;
				if(pContigOverlapInfo1->breakFlag==NO && pContigOverlapInfo2->breakFlag==NO)
				{
					contigID1 = pContigOverlapInfo1->contigID1;
					contigID2 = pContigOverlapInfo1->contigID2;
					contigID3 = pContigOverlapInfo2->contigID2;
					contigLen1 = contigGraph->contigItemArray[contigID1-1].contigLen;
					contigLen2 = contigGraph->contigItemArray[contigID2-1].contigLen;
					contigLen3 = contigGraph->contigItemArray[contigID3-1].contigLen;
					gapSize1 = pContigOverlapInfo1->gapSize;
					gapSize2 = pContigOverlapInfo2->gapSize;
					overlapLen1 =  pContigOverlapInfo1->overlapLen;
					overlapLen2 =  pContigOverlapInfo2->overlapLen;

					// break unreliable links
					if(contigLen2<300 && (overlapLen1>60 || overlapLen2>60))
					{
						// get the paired count between contig 1 and contig 3
						if(getPairedNumBetweenContigsInScaf(&pairedNum_1_3, contigID1, contigID3, pContigOverlapInfo1->orientation1, pContigOverlapInfo2->orientation2, contigGraph)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the paired reads count between contigs, error!\n", __LINE__, __func__);
							return FAILED;
						}

						if(pairedNum_1_3==0)
						{
							//printf("#############+++++++++++######## line=%d, In %s(), contigID1=%d, contigID2=%d, contigID3=%d, unreliable links was broken.\n", __LINE__, __func__, contigID1, contigID2, contigID3);
							pContigOverlapInfo1->breakFlag = YES;
							pContigOverlapInfo2->breakFlag = YES;
						}
					}
				}
			}
		}

		// remove the short contigs at end if the pair count<3
		if(rowsNum>=2)
		{
			// 5' end
			for(i=0; i<rowsNum; i++)
			{
				contigID1 = contigOverlapArray[i].contigID1;
				if(contigOverlapArray[i].breakFlag==NO && contigGraph->contigItemArray[contigID1-1].onlyEnd5==YES)
				{
					//if(contigOverlapArray[i].pairNum<3)
						contigOverlapArray[i].breakFlag = YES;
				}else
				{
					break;
				}
			}

			// 3' end
			for(i=rowsNum-1; i>=0; i--)
			{
				contigID2 = contigOverlapArray[i].contigID2;
				if(contigOverlapArray[i].breakFlag==NO && contigGraph->contigItemArray[contigID2-1].onlyEnd5==YES)
				{
					//if(contigOverlapArray[i].pairNum<3)
						contigOverlapArray[i].breakFlag = YES;
				}else
				{
					break;
				}
			}
		}

		scaffoldItem = scaffoldItem->next;
	}


	if(computeMinPairNumThresInScaf(&minPairNumThreshold, scaffoldSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the minimal pair count threshold, error!\n", __LINE__, __func__);
		return FAILED;
	}

	scaffoldItem = scaffoldSet->scaffoldItemList;
	while(scaffoldItem)
	{
		contigOverlapArray = scaffoldItem->contigOverlapArray;
		rowsNum = scaffoldItem->itemNumContigOverlapArray;
		linkedContigNum = scaffoldItem->linkedContigsNum;

		// break the links if the reads at the one of the very contig ends
		if(linkedContigNum>=2)
		{
			for(i=0; i<rowsNum; i++)
			{
				if(contigOverlapArray[i].breakFlag==NO && contigOverlapArray[i].pairNum<3)
				{
					if(getPairNumBesidesVeryContigEnds(&validPairNum, &invalidPairNumAtEnds, contigOverlapArray+i, contigGraph, readSet)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the paired reads count between contigs, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(minPairNumThreshold==2 && validPairNum<minPairNumThreshold)
					{
						//printf("--=--=--=--=--=-- minPairNumThreshold=%d, oldPairNum=%d, validPairNum=%d, invalidPairNumAtEnds=%d\n", minPairNumThreshold, contigOverlapArray[i].pairNum, validPairNum, invalidPairNumAtEnds);

						contigOverlapArray[i].breakFlag = YES;
					}else if(minPairNumThreshold==1 && validPairNum<minPairNumThreshold && invalidPairNumAtEnds>=1)
					{
						//printf("++=++=++=++=++=++ minPairNumThreshold=%d, oldPairNum=%d, validPairNum=%d, invalidPairNumAtEnds=%d\n", minPairNumThreshold, contigOverlapArray[i].pairNum, validPairNum, invalidPairNumAtEnds);

						contigOverlapArray[i].breakFlag = YES;
					}
				}
			}
		}

		scaffoldItem = scaffoldItem->next;
	}

	return SUCCESSFUL;
}

/**
 * Compute the new minimal paired count threshold between contigs in scaffolds.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeMinPairNumThresInScaf(int32_t *minPairNumThreshold, scaffoldSet_t *scaffoldSet)
{
	int32_t i, j;
	scaffoldItem_t *scaffoldItem;
	contigOverlap_t *contigOverlapArray;
	int32_t linkedContigNum, rowsNum, startRow;
	int64_t pairNumArray[20], totalLinkNum, totalPairNum;

	for(i=0; i<20; i++) pairNumArray[i] = 0;
	totalLinkNum = totalPairNum = 0;

	scaffoldItem = scaffoldSet->scaffoldItemList;
	while(scaffoldItem)
	{
		contigOverlapArray = scaffoldItem->contigOverlapArray;
		rowsNum = scaffoldItem->itemNumContigOverlapArray;
		linkedContigNum = scaffoldItem->linkedContigsNum;

		if(linkedContigNum>=2)
		{
			for(i=0; i<rowsNum; i++)
			{
				if(contigOverlapArray[i].breakFlag==NO)
				{
					totalLinkNum ++;
					totalPairNum += contigOverlapArray[i].pairNum;

					if(contigOverlapArray[i].pairNum>0 && contigOverlapArray[i].pairNum<=20)
						pairNumArray[contigOverlapArray[i].pairNum-1] ++;
				}
			}
		}

		scaffoldItem = scaffoldItem->next;
	}

	if(totalLinkNum>0)
	{
//		for(i=0; i<20; i++)
//		{
//			printf("i=%d, num=%ld, ratio=%.4f\n", i, pairNumArray[i], (double)pairNumArray[i]/totalLinkNum);
//		}

		if((double)pairNumArray[0]/totalLinkNum>0.1)
			*minPairNumThreshold = 1;
		else
			*minPairNumThreshold = 2;
	}else
	{
		*minPairNumThreshold = 1;
//		printf("totalLinkNum=%ld\n", totalLinkNum);
	}

	return SUCCESSFUL;
}

/**
 * Get the paired count between contigs in scaffolds.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getPairedNumBetweenContigsInScaf(int32_t *pairedNum, int32_t contigID1, int32_t contigID2, int32_t contigOrient1, int32_t contigOrient2, contigGraph_t *contigGraph)
{
	int32_t i, contigID_tmp, endFlag1, endFlag2, edgeNum1, edgeNum2;
	contigEdge_t *pEdgeArray1, *pEdgeArray2;

	if(contigOrient1==ORIENTATION_PLUS)
	{ // 3' end
		if(contigGraph->contigItemArray[contigID1-1].alignRegSizeEnd3>0)
			endFlag1 = 1;
		else
			endFlag1 = 2;
	}else
	{ // 5' end
		if(contigGraph->contigItemArray[contigID1-1].alignRegSizeEnd3>0)
			endFlag1 = 0;
		else
			endFlag1 = 2;
	}
/*
	if(contigOrient2==ORIENTATION_PLUS)
	{ // 5' end
		if(contigGraph->contigItemArray[contigID2-1].alignRegSizeEnd3>0)
			endFlag2 = 0;
		else
			endFlag2 = 2;
	}else
	{ // 3' end
		if(contigGraph->contigItemArray[contigID2-1].alignRegSizeEnd3>0)
			endFlag2 = 1;
		else
			endFlag2 = 2;
	}
*/
	// get the edge and count the number
	if(endFlag1==1)
	{
		pEdgeArray1 = contigGraph->contigItemArray[contigID1-1].contigEdgeArrayEnd3;
		edgeNum1 = contigGraph->contigItemArray[contigID1-1].itemNumContigEdgeArrayEnd3;
	}else
	{
		pEdgeArray1 = contigGraph->contigItemArray[contigID1-1].contigEdgeArrayEnd5;
		edgeNum1 = contigGraph->contigItemArray[contigID1-1].itemNumContigEdgeArrayEnd5;
	}
/*
	if(endFlag2==1)
	{
		pEdgeArray1 = contigGraph->contigItemArray[contigID2-1].contigEdgeArrayEnd3;
		edgeNum1 = contigGraph->contigItemArray[contigID2-1].itemNumContigEdgeArrayEnd3;
	}else
	{
		pEdgeArray2 = contigGraph->contigItemArray[contigID2-1].contigEdgeArrayEnd5;
		edgeNum2 = contigGraph->contigItemArray[contigID2-1].itemNumContigEdgeArrayEnd5;
	}
*/
	*pairedNum = 0;
	for(i=0; i<edgeNum1; i++)
	{
		contigID_tmp = pEdgeArray1[i].col / 2 + 1;
		if(contigID_tmp==contigID2)
			(*pairedNum) += pEdgeArray1[i].validNums[pEdgeArray1[i].maxIndexes[0]];
	}

	return SUCCESSFUL;
}


short getPairNumBesidesVeryContigEnds(int32_t *validPairNum, int32_t *invalidPairNumAtEnds, contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph, readSet_t *readSet)
{
	int32_t i, j, contigID1, contigID2, contigEnd1, contigEnd2, contigReadNum;
	int32_t validMatchFlag1, validMatchFlag2, validStartContigPos1, validEndContigPos1, validStartContigPos2, validEndContigPos2;
	int64_t readID1, readID2, seqLen1, seqLen2, contigPos1, contigPos2;
	int32_t readOrient1, readOrient2, validReadOrient1, validReadOrient2;
	contigRead_t *contigReadArray;

	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo;
	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;


	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

	contigID1 = pContigOverlapInfo->contigID1;
	contigID2 = pContigOverlapInfo->contigID2;

	if(pContigOverlapInfo->orientation1==ORIENTATION_PLUS)
	{
		if(contigGraph->contigItemArray[contigID1-1].alignRegSizeEnd3>0)
			contigEnd1 = 1;
		else
			contigEnd1 = 2;
		validReadOrient1 = ORIENTATION_PLUS;
	}else
	{
		if(contigGraph->contigItemArray[contigID1-1].alignRegSizeEnd3>0)
			contigEnd1 = 0;
		else
			contigEnd1 = 2;
		validReadOrient1 = ORIENTATION_MINUS;
	}

	if(pContigOverlapInfo->orientation2==ORIENTATION_PLUS)
	{
		if(contigGraph->contigItemArray[contigID2-1].alignRegSizeEnd3>0)
			contigEnd2 = 0;
		else
			contigEnd2 = 2;
		validReadOrient2 = ORIENTATION_MINUS;
	}else
	{
		if(contigGraph->contigItemArray[contigID2-1].alignRegSizeEnd3>0)
			contigEnd2 = 1;
		else
			contigEnd2 = 2;
		validReadOrient2 = ORIENTATION_PLUS;
	}

	if(contigEnd1==1)
	{
		contigReadArray = contigGraph->contigItemArray[contigID1-1].contigReadArrayEnd3;
		contigReadNum = contigGraph->contigItemArray[contigID1-1].contigReadNumEnd3;
	}else
	{
		contigReadArray = contigGraph->contigItemArray[contigID1-1].contigReadArrayEnd5;
		contigReadNum = contigGraph->contigItemArray[contigID1-1].contigReadNumEnd5;
	}

	if(contigEnd1==1)
	{
		validStartContigPos1 = contigGraph->contigItemArray[contigID1-1].contigLen - contigGraph->contigItemArray[contigID1-1].alignRegSizeEnd3 + 1;
		validEndContigPos1 = contigGraph->contigItemArray[contigID1-1].contigLen - 1.5 * kmerSize;
	}else
	{
		validStartContigPos1 = 1.5 * kmerSize + 1;
		validEndContigPos1 = contigGraph->contigItemArray[contigID1-1].alignRegSizeEnd5;
	}

	if(contigEnd2==1)
	{
		validStartContigPos2 = contigGraph->contigItemArray[contigID2-1].contigLen - contigGraph->contigItemArray[contigID2-1].alignRegSizeEnd3 + 1;
		validEndContigPos2 = contigGraph->contigItemArray[contigID2-1].contigLen - 1.5 * kmerSize;
	}else
	{
		validStartContigPos2 = 1.5 * kmerSize + 1;
		validEndContigPos2 = contigGraph->contigItemArray[contigID2-1].alignRegSizeEnd5;
	}

	*validPairNum = 0;
	*invalidPairNumAtEnds = 0;
	for(i=0; i<contigReadNum; i++)
	{
		readID1 = contigReadArray[i].readID;
		contigPos1 = contigReadArray[i].contigPos;
		seqLen1 = contigReadArray[i].seqlen;
		readOrient1 = contigReadArray[i].orientation;

		if(readID1%2==1)
			readID2 = readID1 + 1;
		else
			readID2 = readID1 - 1;

		readMatchInfoBlockID = (readID2 - 1) / maxItemNumPerReadMatchInfoBlock;
		rowNumInReadMatchInfoBlock = (readID2 - 1) % maxItemNumPerReadMatchInfoBlock;
		pReadMatchInfo = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;

		contigPos2 = pReadMatchInfo->contigPos;
		seqLen2 = pReadMatchInfo->seqlen;
		readOrient2 = pReadMatchInfo->readOrientation;

		if(pReadMatchInfo->contigID==contigID2 && pReadMatchInfo->contigEnd==contigEnd2 && readOrient1==validReadOrient1 && readOrient2==validReadOrient2)
		{
			if((contigPos1>=validStartContigPos1 && contigPos1+seqLen1-1<=validEndContigPos1) && (contigPos2>=validStartContigPos2 && contigPos2+seqLen2-1<=validEndContigPos2))
			{
				(*validPairNum) ++;
			}else if((seqLen1<0.4*readLen || seqLen2<0.4*readLen) || (seqLen1<0.5*readLen && seqLen2<0.5*readLen))
			{
				(*invalidPairNumAtEnds) ++;
			}
		}
	}

	return SUCCESSFUL;
}


/**
 * Split the broken scaffolds.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short splitScaffolds(scaffoldSet_t *scaffoldSet)
{
	int32_t i, j;
	scaffoldItem_t *scaffoldItem, *newScaffoldItem;
	contigOverlap_t *contigOverlapArray, *pContigOverlapInfo, *newContigOverlapArray;
	int32_t localScaffoldID, linkedContigNum, rowsNum, startRow;
	int32_t subRowsNum, subStartRow, subEndRow, subLinkedContigNum;

	// break the erroneous links
	localScaffoldID = 0;
	scaffoldItem = scaffoldSet->scaffoldItemList;
	while(scaffoldItem)
	{
		localScaffoldID ++;
		scaffoldItem->scaffoldID = localScaffoldID;

		//if(scaffoldItem->scaffoldID==4)
		//{
		//	printf("line=%d, In %s(), scaffoldID=%d\n", __LINE__, __func__, scaffoldItem->scaffoldID);
		//}

		contigOverlapArray = scaffoldItem->contigOverlapArray;
		rowsNum = scaffoldItem->itemNumContigOverlapArray;
		linkedContigNum = scaffoldItem->linkedContigsNum;

		if(linkedContigNum<=1)
		{ // only one contig, then do nothing
			if(linkedContigNum<1)
			{ // error linkedContigNum
				printf("line=%d, In %s(), linkedContigNum=%d, error!\n", __LINE__, __func__, linkedContigNum);
				return FAILED;
			}
		}else
		{
			// generate new linked information
			for(subEndRow=0; subEndRow<rowsNum; subEndRow++)
			{
				if(contigOverlapArray[subEndRow].breakFlag==YES)
					break;
			}

			if(subEndRow<rowsNum)
			{ // A break happens

				subStartRow = subEndRow + 1;
				subRowsNum = rowsNum - subEndRow - 1;

				// allocate memory to store the new sub-item
				newScaffoldItem = (scaffoldItem_t*) calloc (1, sizeof(scaffoldItem_t));
				if(newScaffoldItem==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				newScaffoldItem->scaffoldID = localScaffoldID + 1;
				newScaffoldItem->scaffoldLen = -1;
				newScaffoldItem->previous = scaffoldItem;
				newScaffoldItem->next = scaffoldItem->next;
				if(scaffoldItem->next)
					scaffoldItem->next->previous = newScaffoldItem;
				scaffoldItem->next = newScaffoldItem;
				scaffoldSet->scaffoldNum ++;

				if(subEndRow==rowsNum-1)
				{ // break at the tail row
					newScaffoldItem->itemNumContigOverlapArray = 1;
					newScaffoldItem->linkedContigsNum = 1;
				}else
				{ // break not at the tail row
					newScaffoldItem->itemNumContigOverlapArray = subRowsNum;
					newScaffoldItem->linkedContigsNum = subRowsNum + 1;
				}

				newContigOverlapArray = (contigOverlap_t*) calloc (newScaffoldItem->itemNumContigOverlapArray, sizeof(contigOverlap_t));
				if(newContigOverlapArray==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				newScaffoldItem->contigOverlapArray = newContigOverlapArray;

				if(subEndRow==rowsNum-1)
				{ // break at the tail row
					newContigOverlapArray[0].contigID1 = contigOverlapArray[subEndRow].contigID2;
					newContigOverlapArray[0].orientation2 = ORIENTATION_PLUS;
				}else
				{
					if(memcpy(newContigOverlapArray, contigOverlapArray+subStartRow, newScaffoldItem->itemNumContigOverlapArray*sizeof(contigOverlap_t))==NULL)
					{
						printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

				// deal with the old item
				if(subEndRow==0)
				{ // break at the head row
					scaffoldItem->itemNumContigOverlapArray = 1;
					scaffoldItem->scaffoldLen = 0;
					scaffoldItem->linkedContigsNum = 1;
					contigOverlapArray[0].orientation1 = ORIENTATION_PLUS;
					contigOverlapArray[0].contigID2 = 0;
					contigOverlapArray[0].orientation2 = 0;
					contigOverlapArray[0].mergeFlag = NO;
					contigOverlapArray[0].breakFlag = NO;
					contigOverlapArray[0].overlapLen = 0;
					contigOverlapArray[0].gapSize = 0;
				}else
				{ // break not at the head row
					scaffoldItem->itemNumContigOverlapArray = subEndRow;
					scaffoldItem->linkedContigsNum = subEndRow + 1;
				}
			}
		}

		scaffoldItem = scaffoldItem->next;
	}

	return SUCCESSFUL;
}
