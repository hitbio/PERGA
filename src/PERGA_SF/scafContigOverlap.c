/*
 * scafContigOverlap.c
 *
 *  Created on: Jul 4, 2011
 *      Author: zhuxiao
 */

#include "inc/scafStdinc.h"
#include "inc/scafExtvab.h"


/**
 * Generate overlap information between contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short overlapContigsInScaffolds(const char *contigOverlapInfoFile, const char *meanSdevFile, const char *newContigFile, const char *newSharedReadListFile, const char *newReadListFile1, const char *newReadListFile2, const char *newContigListFile, const char *linkResultFile, const char *contigFile, const char *sharedReadListFile, const char *readListFile1, const char *readListFile2, const char *contigListFile, const char *averLinkNumFile)
{

	printf("=========== Begin generating overlap information of contigs, please wait ...\n");

	// initialize the memory and variables, including scaffold information,
	// and its corresponding contig overlap information
	if(initMemContigOverlap(linkResultFile, contigFile, sharedReadListFile, readListFile1, readListFile2, contigListFile, averLinkNumFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// estimate mean size and standard deviation of paired ends
	if(meanSizeSDEstimate()==FAILED)
	{
		printf("line=%d, In %s(), cannot estimate the gap size, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// generate contig overlap information of a scaffold
	if(generateContigOverlapInfo()==FAILED)
	{
		printf("line=%d, In %s(), cannot generate contig overlap information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// output mean size and standard deviation of fragments to text file
	//  format: mean size \t standard deviation
	if(outputMeanSizeAndSDev(meanSdevFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot output mean size and standard deviation of fragments to file [ %s ], error!\n", __LINE__, __func__, meanSdevFile);
		return FAILED;
	}

	// output the contigs
	if(outputContigInfoArrToFile(newContigFile, contigInfoArr, contigsNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot output the contig information array to file [ %s ], error!\n", __LINE__, __func__, newContigFile);
		return FAILED;
	}

	// output contig overlap information to text file.
	//  Format:
	//  	(1) Header fields: >scaffoldID, linkedNum, rowsNum, which are separated by tab character;
	//  	(2)   Body fields: contigID1, orientation1, contigLen1, contigID2, orientation2, contigLen2, mergeFlag, overlapLen, gapSize, breakFlag, which are separated by tab character.
	if(outputContigOverlapInfoToFile(contigOverlapInfoFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot output the contig overlap information to file [%s ], error!\n", __LINE__, __func__, contigOverlapInfoFile);
		return FAILED;
	}

	// rewrite the RL1, RL2, SRL and CL to binary files
	if(rewriteReadListsAndContigListToFiles(newSharedReadListFile, newReadListFile1, newReadListFile2, newContigListFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot rewrite the read lists and contig list to file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// free the memory
	freeMemContigOverlap();

	printf("=========== End generated overlap information of contigs.\n");

	return SUCCESSFUL;
}

/**
 * Initialize the memory for contig overlap.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemContigOverlap(const char *linkResultFile, const char *contigFile, const char *sharedReadListFile, const char *readListFile1, const char *readListFile2, const char *contigListFile, const char *averLinkNumFile)
{
	int i;
	// initialize the contig information array
	if(initContigInfoArray(contigFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the contig information array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize the contig overlap information array
	if(initContigOverlapInfoArray(linkResultFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the contig overlap information array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Allocate memory and load the data of shared read list (SRL)
	if(loadSingleReadList(sharedReadListFile, &sharedReadListArr, &readItemNumInSRL, &sharedReadPosArr, &matchItemNumInSRP)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the data of shared read list (SRL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Allocate memory and load the data of read list 1 (RL1) and read list 2 (RL2)
	if(loadSingleReadList(readListFile1, readListArr+1, readNumInRL+1, readPosArr+1, matchItemNumInRP+1)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the data of a read list (RL), error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(loadSingleReadList(readListFile2, readListArr+2, readNumInRL+2, readPosArr+2, matchItemNumInRP+2)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the data of a read list (RL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Allocate memory and load the data of contig list (CL)
	if(loadContigList(contigListFile, &contigListArr, &contigItemNumInCL, &contigReadArr, &contigReadItemNumInCR)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the data of contig list (CL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the averLinkNum from averLinkNumFile
	if(getAverLinkNum(&averLinkNum, averLinkNumFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the averLinkNum, error!\n", __LINE__, __func__);
		return FAILED;
	}

	meanSizeInsert = 0;
	stardardDeviationInsert = 0;
	standardDeviationFactor = SDEV_FACTOR;


	// the auxiliary memory for overlapping
	minLinksNumContigsThres = averLinkNum * MIN_FIRST_LINKNUM_FACTOR;
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

	// free contig information array
	freeMemContigInfo(&contigInfoArr, &contigsNum);

	scaffoldsNumInCOI = 0;
	free(contigOverlapIndexArr);
	contigOverlapIndexArr = NULL;

	itemNumInContigOverlapArr = 0;
	free(contigOverlapInfoArr);
	contigOverlapInfoArr = NULL;

	// Release the memory of shared read list (SRL)
	freeSingleReadList(&sharedReadListArr, &readItemNumInSRL, &sharedReadPosArr, &matchItemNumInSRP);

	// Release the memory of read list 1 (RL1) and read list 2 (RL2)
	freeSingleReadList(readListArr+1, readNumInRL+1, readPosArr+1, matchItemNumInRP+1);
	freeSingleReadList(readListArr+2, readNumInRL+2, readPosArr+2, matchItemNumInRP+2);

	// free contig list (CL)
	freeContigList(&contigListArr, &contigItemNumInCL, &contigReadArr, &contigReadItemNumInCR);

	meanSizeInsert = 0;
	stardardDeviationInsert = 0;
	standardDeviationFactor = 0;

	minLinksNumContigsThres = 0;
	maxOverlapSeqLen = 0;
	minOverlapThres = 0;
	minExactOverlapThres = 0;
	mismatchThres = 0;
	maxMisMatchRatioThres = 0;
	subMismatchThres = 0;
	minAdjustGapSizeThres = 0;
	maxAdjustGapSizeThres = 0;
	minBaseNumInGap = 0;
	gapSizeSdevFactorOverlap = 0;
	//exactOverlapSdevThres = 0;
	matchScore = 0;
	mismatchScore = 0;
	gapScore = 0;

	breakLinkNumThres = 0;

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
 * Initialize the contig overlap information array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initContigOverlapInfoArray(const char *linkResultFile)
{
	// get the total number of scaffolds
	scaffoldsNumInCOI = getScaffoldsNum(linkResultFile);
	if(scaffoldsNumInCOI<=0)
	{
		printf("line=%d, In %s(), cannot get the total number of scaffolds, error\n", __LINE__, __func__);
		return FAILED;
	}

	itemNumInContigOverlapArr = getItemNumContigOverlapInfo(linkResultFile);
	if(itemNumInContigOverlapArr<=0)
	{
		printf("line=%d, In %s(), cannot get the total item number of contig overlap information, error\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate the memory
	contigOverlapIndexArr = (contigOverlapIndex*) calloc(scaffoldsNumInCOI, sizeof(contigOverlapIndex));
	if(contigOverlapIndexArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error\n", __LINE__, __func__);
		return FAILED;
	}

	contigOverlapInfoArr = (contigOverlap*) calloc(itemNumInContigOverlapArr, sizeof(contigOverlap));
	if(contigOverlapInfoArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the data of contig overlap index array and contig overlap information array
	if(fillContigOverlapData(linkResultFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the contig overlap data, error\n", __LINE__, __func__);
		return FAILED;
	}


	return SUCCESSFUL;
}

/**
 * Get total item number of contig overlap information from linked contigs file.
 *  @return:
 *  	If succeeds, return the item number; otherwise, return ERROR.
 */
int getItemNumContigOverlapInfo(const char *linkResultFile)
{
	FILE *fpLinkResult;
	int totalItemNum, scaffoldID, num;
	char ch;

	fpLinkResult = fopen(linkResultFile, "r");
	if(fpLinkResult==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, linkResultFile);
		return ERROR;
	}

	// compute the total number
	totalItemNum = 0;
	ch = fgetc(fpLinkResult);
	while(ch!=EOF)
	{
		if(ch=='>')
		{
			fscanf(fpLinkResult, "%d\t%d\n", &scaffoldID, &num);
			if(num>=2)
			{
				totalItemNum += num - 1;
			}else
			{
				totalItemNum += 1;
			}
		}

		ch = fgetc(fpLinkResult);
	}

	fclose(fpLinkResult);

	return totalItemNum;
}

/**
 * Get the average linked read pairs per contigEdge.
 *  @return:
 *   If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getAverLinkNum(double *averageLinkNum, const char *averLinkNumFile)
{
	FILE *fpAverLinkNum;

	fpAverLinkNum = fopen(averLinkNumFile, "r");
	if(fpAverLinkNum==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, averLinkNumFile);
		return FAILED;
	}

	if(fscanf(fpAverLinkNum, "%lf\n", averageLinkNum)!=1)
	{
		printf("line=%d, In %s(), cannot read file [ %s ], error!\n", __LINE__, __func__, averLinkNumFile);
		return FAILED;
	}

	fclose(fpAverLinkNum);
	fpAverLinkNum = NULL;

	return SUCCESSFUL;
}

/**
 * Fill the contig overlap data.
 *  @return:
 *   If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillContigOverlapData(const char *linkResultFile)
{
	FILE *fpLinkResult;
	char ch;
	int i, j, scaffoldID, num/*, tmp_scaffoldItemNum*/, tmp_startRow;
	int tmp_contigID1, tmp_orient1, tmp_contigLen1, tmp_contigID2, tmp_orient2, tmp_contigLen2;

	fpLinkResult = fopen(linkResultFile, "r");
	if(fpLinkResult==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, linkResultFile);
		return FAILED;
	}

	// fill the data
	//tmp_scaffoldItemNum = 0;
	tmp_startRow = 0;
	for(i=0; i<scaffoldsNumInCOI; i++)
	{
		ch = fgetc(fpLinkResult);
		if(ch=='>')
		{
			fscanf(fpLinkResult, "%d\t%d\n", &scaffoldID, &num);
			contigOverlapIndexArr[scaffoldID-1].scaffoldID = scaffoldID;
			contigOverlapIndexArr[scaffoldID-1].linkedNum = num;
			contigOverlapIndexArr[scaffoldID-1].startRow = tmp_startRow;

			if(num>=2)
				contigOverlapIndexArr[scaffoldID-1].rowsNum = num - 1;
			else
				contigOverlapIndexArr[scaffoldID-1].rowsNum = 1;

			fscanf(fpLinkResult, "%d\t%d\t%d\n", &tmp_contigID1, &tmp_orient1, &tmp_contigLen1);

			if(num>=2)
			{
				for(j=0; j<num-1; j++)
				{
					fscanf(fpLinkResult, "%d\t%d\t%d\n", &tmp_contigID2, &tmp_orient2, &tmp_contigLen2);
					contigOverlapInfoArr[tmp_startRow+j].contigID1 = tmp_contigID1;
					contigOverlapInfoArr[tmp_startRow+j].orientation1 = tmp_orient1;
					//contigOverlapInfoArr[tmp_startRow+j].contigLen1 = tmp_contigLen1;
					contigOverlapInfoArr[tmp_startRow+j].contigID2 = tmp_contigID2;
					contigOverlapInfoArr[tmp_startRow+j].orientation2 = tmp_orient2;
					//contigOverlapInfoArr[tmp_startRow+j].contigLen2 = tmp_contigLen2;

					// set contig1 to contig2
					tmp_contigID1 = tmp_contigID2;
					tmp_orient1 = tmp_orient2;
					//tmp_contigLen1 = tmp_contigLen2;

				}
			}else
			{
				contigOverlapInfoArr[tmp_startRow].contigID1 = tmp_contigID1;
				contigOverlapInfoArr[tmp_startRow].orientation1 = tmp_orient1;
				//contigOverlapInfoArr[tmp_startRow].contigLen1 = tmp_contigLen1;
			}

			tmp_startRow += contigOverlapIndexArr[scaffoldID-1].rowsNum;
		}else
		{
			printf("line=%d, In %s(), cannot read file [ %s ], error!\n", __LINE__, __func__, linkResultFile);
			return FAILED;
		}
	}

	fclose(fpLinkResult);
	fpLinkResult = NULL;

	return SUCCESSFUL;
}


/**
 * Generate the contig overlap information.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateContigOverlapInfo()
{
	int scaffoldID, rowsNum;
	contigOverlap *pContigOverlapInfo;
	int i, j;

	// generate contig overlap information of scaffolds
	for(i=0; i<scaffoldsNumInCOI; i++)
	{
		//####################### Deubg information #####################
#if DEBUG_FLAG
		printf("\n=================== scaffold: %d, itemsNum: %d ================\n", contigOverlapIndexArr[i].scaffoldID, contigOverlapIndexArr[i].rowsNum);
#endif
		//####################### Deubg information #####################

		scaffoldID = contigOverlapIndexArr[i].scaffoldID;
		rowsNum = contigOverlapIndexArr[i].rowsNum;
		pContigOverlapInfo = contigOverlapInfoArr + contigOverlapIndexArr[i].startRow;

		if(contigOverlapIndexArr[i].linkedNum>=2)
		{ // just consider the scaffolds having more than 1 linked contigs
			for(j=0; j<rowsNum; j++)
			{
				//####################### Debug information ###################
#if DEBUG_FLAG
				if(pContigOverlapInfo[j].contigID1==78 && pContigOverlapInfo[j].contigID2==81)
				{
					printf("contig1: [%d,%d,%d]; contig2: [%d,%d,%d]\n", pContigOverlapInfo[j].contigID1, pContigOverlapInfo[j].orientation1, contigInfoArr[pContigOverlapInfo[j].contigID1-1].contigLen, pContigOverlapInfo[j].contigID2, pContigOverlapInfo[j].orientation2, contigInfoArr[pContigOverlapInfo[j].contigID2-1].contigLen);
				}

				printf("*****contig1: [%d,%d,%d]; contig2: [%d,%d,%d]\n", pContigOverlapInfo[j].contigID1, pContigOverlapInfo[j].orientation1, contigInfoArr[pContigOverlapInfo[j].contigID1-1].contigLen, pContigOverlapInfo[j].contigID2, pContigOverlapInfo[j].orientation2, contigInfoArr[pContigOverlapInfo[j].contigID2-1].contigLen);
#endif
				//####################### Debug information ###################

				if(updateContigOverlapLen(pContigOverlapInfo+j, contigInfoArr)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute overlap length for scaffold: %d, error!\n", __LINE__, __func__, scaffoldID);
					return FAILED;
				}
			}
		}else if(contigOverlapIndexArr[i].linkedNum==1)
		{
			//####################### Debug information ###################
#if DEBUG_FLAG
			printf("*****contig1: [%d,0,%d]; contig2: [0,0,0]\n", pContigOverlapInfo[0].contigID1, contigInfoArr[pContigOverlapInfo[0].contigID1-1].contigLen);
#endif
			//####################### Debug information ###################
		}else
		{
			printf("line=%d, In %s(), contig overlap information error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Compute the overlap length of two contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateContigOverlapLen(contigOverlap *pContigOverlapInfo, contigInfo *pContigInfoArr)
{
	int contigID1, orientation1, contigLen1, contigID2, orientation2, contigLen2;
	int seq_len1, seq_len2, overlapLenExact, overlapLenAlignment, overlapLenAdjust, mismatchNum;
	int gapSize, validPairedNum_gapEstimate;
	int contigLenBeforeAlignment[2];

	// get the information of contg1 and contig2
	contigID1 = pContigOverlapInfo->contigID1;
	orientation1 = pContigOverlapInfo->orientation1;
	contigLen1 = pContigInfoArr[contigID1-1].contigLen;
	contigID2 = pContigOverlapInfo->contigID2;
	orientation2 = pContigOverlapInfo->orientation2;
	contigLen2 = pContigInfoArr[contigID2-1].contigLen;

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
		strncpy(overlapSeq1, pContigInfoArr[contigID1-1].contigSeq+contigLen1-seq_len1, seq_len1);
		overlapSeq1[seq_len1] = '\0';
	}else
	{
		strncpy(overlapSeq1, pContigInfoArr[contigID1-1].contigSeq, seq_len1);
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
		strncpy(overlapSeq2, pContigInfoArr[contigID2-1].contigSeq, seq_len2);
		overlapSeq2[seq_len2] = '\0';
	}else
	{
		strncpy(overlapSeq2, pContigInfoArr[contigID2-1].contigSeq+contigLen2-seq_len2, seq_len2);
		overlapSeq2[seq_len2] = '\0';
		if(reverseSeq(overlapSeq2, seq_len2)==FAILED) // reverse the sequence
		{
			printf("line=%d, In %s(), cannot reverse the sequence [ %s ], error!\n", __LINE__, __func__, overlapSeq2);
			return FAILED;
		}
	}

	//############################# Debug information #########################
#if DEBUG_OUT_FLAG
	printf("overlapSeq1=%s, len=%d\noverlapSeq2=%s, len=%d\n", overlapSeq1, (int)strlen(overlapSeq1), overlapSeq2, (int)strlen(overlapSeq2));
#endif
	//############################# Debug information #########################

	contigLenBeforeAlignment[0] = contigLen1;
	contigLenBeforeAlignment[1] = contigLen2;

	// estimate gap size between contigs
	if(gapSizeEstimateBetweenContigs(pContigOverlapInfo, &gapSize, &validPairedNum_gapEstimate, -1, NULL, NULL)==FAILED)
	{
		printf("line=%d, In %s(), cannot estimate the gap size between contigs [%d,%d], error!\n", __LINE__, __func__, contigID1, contigID2);
		return FAILED;
	}

	//printf("line=%d, In %s(), gapSize=%d.\n", __LINE__, __func__, gapSize);

	//if(validPairedNum_gapEstimate<=0)
	if(validPairedNum_gapEstimate<minLinksNumContigsThres)
	//if(validPairedNum_gapEstimate<minLinksNumContigsThres || (gapSize>0.6*meanSizeInsert && gapSize>2*averReadLen)) // 2012-11-19
	{
#if (DEBUG_FLAG==YES)
		printf("line=%d, In %s(), broken.\n", __LINE__, __func__);
#endif
		pContigOverlapInfo->breakFlag = YES;
		pContigOverlapInfo->gapSize = 0;
		pContigOverlapInfo->mergeFlag = NO;

		//printf("line=%d, In %s(), validPairedNum_gapEstimate=%d, error!\n", __LINE__, __func__, validPairedNum_gapEstimate);
		return SUCCESSFUL;
	}

	// update the contig overlap information
	//if(gapSize<=stardardDeviationInsert)
	if(gapSize<=gapSizeSdevFactorOverlap*stardardDeviationInsert)
	{
		// compute the overlap length by exact alignment
		if(computeSeqOverlapLenExact(&overlapLenExact, overlapSeq1, seq_len1, overlapSeq2, seq_len2, scoreArr, gapSize)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute sequence overlap length, error!\n", __LINE__, __func__);
			return FAILED;
		}

		//if(overlapLenExact>=minOverlapThres)
		if(overlapLenExact>=minOverlapThres && (overlapLenExact>=-gapSize-stardardDeviationInsert && overlapLenExact<=-gapSize+stardardDeviationInsert))
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

					if(updateContigs(pContigInfoArr+contigID1-1, pContigInfoArr+contigID2-1, orientation1, orientation2, overlapSeq1, seq_len1, overlapSeq2, seq_len2)==FAILED)
					{
						printf("line=%d, In %s(), cannot update contig end sequences, error!\n", __LINE__, __func__);
						return FAILED;
					}

					//############################### Debug information #########################
#if DEBUG_OUT_FLAG
					//printf("line=%d, In %s(), after updateContigs(), contigLen1=%d, contigLen2=%d\n", __LINE__, __func__, pContigInfoArr[contigID1-1].contigLen, pContigInfoArr[contigID2-1].contigLen);
#endif
					//############################### Debug information #########################

					pContigOverlapInfo->mergeFlag = YES;
					pContigOverlapInfo->overlapLen = overlapLenAdjust;

					if(updateContigEndsInfo(pContigOverlapInfo, pContigInfoArr)==FAILED)
					{
						printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
						return FAILED;
					}

					//=====================================================================
					// update the read list 1 (RL1), read list 2 (RL2), shared read list (SRL) and contig list (CL)
					if(pContigInfoArr[contigID1-1].contigLen != contigLenBeforeAlignment[0])
					{
						if(updateReadListsAndContigListSingleContigAfterAlignment(pContigInfoArr+contigID1-1, orientation1, 1, contigLenBeforeAlignment[0])==FAILED)
						{
							printf("line=%d, In %s(), cannot update read lists and contig list after alignment, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}

					if(pContigInfoArr[contigID2-1].contigLen != contigLenBeforeAlignment[1])
					{
						if(updateReadListsAndContigListSingleContigAfterAlignment(pContigInfoArr+contigID2-1, orientation2, 2, contigLenBeforeAlignment[1])==FAILED)
						{
							printf("line=%d, In %s(), cannot update read lists and contig list after alignment, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
					//=====================================================================

				}else // if(overlapLenAdjust==0)
				{ // there is a gap between the two contigs with high probability
					if(validPairedNum_gapEstimate>=minLinksNumContigsThres)
					{
						if(gapSize<minAdjustGapSizeThres)
						{ // gapSize < -10
							pContigOverlapInfo->gapSize = gapSize;
							pContigOverlapInfo->mergeFlag = NO;

							// update the update length by cutting uncovered contig ends
							if(updateOverlapLenByCutUncoveredContigEnds(pContigOverlapInfo, pContigInfoArr)==FAILED)
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

							if(updateContigEndsInfo(pContigOverlapInfo, pContigInfoArr)==FAILED)
							{
								printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
								return FAILED;
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
								if(updateContigs(pContigInfoArr+contigID1-1, pContigInfoArr+contigID2-1, orientation1, orientation2, overlapSeq1, seq_len1, overlapSeq2, seq_len2)==FAILED)
								{
									printf("line=%d, In %s(), cannot update contig end sequences, error!\n", __LINE__, __func__);
									return FAILED;
								}
								pContigOverlapInfo->mergeFlag = YES;
								pContigOverlapInfo->overlapLen = overlapLenAdjust;

								if(updateContigEndsInfo(pContigOverlapInfo, pContigInfoArr)==FAILED)
								{
									printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
									return FAILED;
								}

								//=====================================================================
								// update the read list 1 (RL1), read list 2 (RL2), shared read list (SRL) and contig list (CL)
								if(pContigInfoArr[contigID1-1].contigLen != contigLenBeforeAlignment[0])
								{
									if(updateReadListsAndContigListSingleContigAfterAlignment(pContigInfoArr+contigID1-1, orientation1, 1, contigLenBeforeAlignment[0])==FAILED)
									{
										printf("line=%d, In %s(), cannot update read lists and contig list after alignment, error!\n", __LINE__, __func__);
										return FAILED;
									}
								}

								if(pContigInfoArr[contigID2-1].contigLen != contigLenBeforeAlignment[1])
								{
									if(updateReadListsAndContigListSingleContigAfterAlignment(pContigInfoArr+contigID2-1, orientation2, 2, contigLenBeforeAlignment[1])==FAILED)
									{
										printf("line=%d, In %s(), cannot update read lists and contig list after alignment, error!\n", __LINE__, __func__);
										return FAILED;
									}
								}
								//=====================================================================

							}else
							{
								pContigOverlapInfo->gapSize = gapSize;
								pContigOverlapInfo->mergeFlag = NO;

								// update the update length by cutting uncovered contig ends
								if(updateOverlapLenByCutUncoveredContigEnds(pContigOverlapInfo, pContigInfoArr)==FAILED)
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
							if(updateOverlapLenByCutUncoveredContigEnds(pContigOverlapInfo, pContigInfoArr)==FAILED)
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

						if(updateContigEndsInfo(pContigOverlapInfo, pContigInfoArr)==FAILED)
						{
							printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
							return FAILED;
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

		if(updateContigEndsInfo(pContigOverlapInfo, pContigInfoArr)==FAILED)
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
 * Compute sequence overlap length.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short computeSeqOverlapLenByAlignment(const char *seq1, const int seqLen1, const char *seq2, const int seqLen2, int *scoreArray, char **pAlignResultArray, int *overlapLen, int *mismatchNum)
{
	int i, j, tmp;
	int rowsNum, colsNum;
	int maxValue, scoreIJ, maxRow, maxCol;
	int itemNumInAlignArray;

	// ####################### Debug information ########################
#if DEBUG_OUT_FLAG
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

#if DEBUG_OUT_FLAG
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

#if DEBUG_OUT_FLAG
	printf("i=%d, j=%d, mismatchNum=%d\n", i, j, *mismatchNum);
	printf("maxRow-i=%d, maxCol-j=%d, itemNumInAlignArray=%d\n", maxRow-i, maxCol-j, itemNumInAlignArray);
	printf("Alignment overlapLen=%d\n", *overlapLen);
#endif

	//if(*mismatchNum > 2 || (i==0 && j>0))
	if(i==0 && j>0)
	{
		*overlapLen = 0;
	}

#if DEBUG_OUT_FLAG
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

#if DEBUG_OUT_FLAG
	// print the alignment result
	for(i=0 ;i<3; i++)
	{
		printf("%s\n", pAlignResultArray[i]);
	}
#endif

	return SUCCESSFUL;
}

/**
 * Compute sequence overlap length.
 *  @return:
 *   If succeeds, return FAILED; otherwise, return FAILED.
 */
short computeSeqOverlapLenExact(int *overlapLen, const char *seq1, const int seqLen1, const char *seq2, const int seqLen2, int *scoreArray, int gapSize)
{
	int i, j, k;
	int rowsNum, colsNum, minDiff;
	int tmpOverlap;

	// ####################### Debug information ########################
#if DEBUG_OUT_FLAG
	printf("seq1=%s, len=%d\nseq2=%s, len=%d\n", seq1, (int)strlen(seq1), seq2, (int)strlen(seq2));
#endif
	// ####################### Debug information ########################

	// reset value of each item in the score array to zero
	if(memset(scoreArray, 0L, (maxOverlapSeqLen+1)*(maxOverlapSeqLen+1)*sizeof(int))==NULL)
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
#if DEBUG_OUT_FLAG
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

#if DEBUG_OUT_FLAG
	printf("Final exact overlapLen=%d\n", *overlapLen);
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
#if DEBUG_OUT_FLAG
	printf("Alignment overlapLen=%d\n", *overlapLen);
#endif
	//############################### Debug information ########################

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

#if DEBUG_OUT_FLAG
	printf("endNum=%d\n", endNum);
	printf("overlapLen=%d, pos1=%d, pos2=%d, pos3=%d\n", *overlapLen, startPos[0], startPos[1], startPos[2]);
#endif

	// check the middle part
	for(i=startPos[1]; i<startPos[2]; i++)
	{
		if(pAlignResultArray[1][i]==' ')
		{
#if DEBUG_OUT_FLAG
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

#if DEBUG_OUT_FLAG
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
		printf("accordanceSeqLen=%d != strlen(accordanceSeq)=%d, error!\n", accordanceSeqLen, (int)strlen(accordanceSeq));
		return FAILED;
	}

#if DEBUG_OUT_FLAG
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

#if DEBUG_OUT_FLAG
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
short updateContigs(contigInfo *pContigInfoArr1, contigInfo *pContigInfoArr2, const int contigOrient1, const int contigOrient2, char *seq1, const int originalSeqLen1, char *seq2, const int originalSeqLen2)
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
			newSeq = (char *) malloc((pContigInfoArr1->contigLen+1 + seqLen1-originalSeqLen1) * sizeof(char));
			if(newSeq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// copy the sequence
			strcpy(newSeq, pContigInfoArr1->contigSeq);

			free(pContigInfoArr1->contigSeq);
			pContigInfoArr1->contigSeq = newSeq;
		}

		strcpy(pContigInfoArr1->contigSeq+pContigInfoArr1->contigLen-originalSeqLen1, seq1);
		pContigInfoArr1->contigLen += seqLen1 - originalSeqLen1;

	}else
	{ // update the sequence at the 5' end of contig 1
		if(reverseSeq(seq1, seqLen1)==FAILED)
		{
			printf("line=%d, In %s(), cannot reverse the sequence [ %s ], error!\n", __LINE__, __func__, seq1);
			return FAILED;
		}

		if(seqLen1!=originalSeqLen1)
		{
			newSeq = (char *) malloc((pContigInfoArr1->contigLen+1 + seqLen1-originalSeqLen1) * sizeof(char));
			if(newSeq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			strcpy(newSeq, seq1);
			strcat(newSeq, pContigInfoArr1->contigSeq+originalSeqLen1);

			free(pContigInfoArr1->contigSeq);
			pContigInfoArr1->contigSeq = newSeq;

			pContigInfoArr1->contigLen += seqLen1 - originalSeqLen1;
		}else
		{
			for(i=0; i<seqLen1; i++) pContigInfoArr1->contigSeq[i] = seq1[i];
		}
	}

	if(contigOrient2==ORIENTATION_PLUS)
	{ // update the sequence at the 5' end of contig 2
		if(seqLen2!=originalSeqLen2)
		{
			newSeq = (char *) malloc((pContigInfoArr2->contigLen+1 + seqLen2-originalSeqLen2) * sizeof(char));
			if(newSeq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			strcpy(newSeq, seq2);
			strcat(newSeq, pContigInfoArr2->contigSeq+originalSeqLen2);

			free(pContigInfoArr2->contigSeq);
			pContigInfoArr2->contigSeq = newSeq;

			pContigInfoArr2->contigLen += seqLen2 - originalSeqLen2;
		}else
		{
			for(i=0; i<seqLen2; i++) pContigInfoArr2->contigSeq[i] = seq2[i];
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
			newSeq = (char *) malloc((pContigInfoArr2->contigLen+1 + seqLen2-originalSeqLen2) * sizeof(char));
			if(newSeq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// copy the sequence
			strcpy(newSeq, pContigInfoArr2->contigSeq);

			free(pContigInfoArr2->contigSeq);
			pContigInfoArr2->contigSeq = newSeq;

		}

		strcpy(pContigInfoArr2->contigSeq+pContigInfoArr2->contigLen-originalSeqLen2, seq2);
		pContigInfoArr2->contigLen += seqLen2 - originalSeqLen2;

	}

	return SUCCESSFUL;
}

/**
 * Estimate mean size and standard deviation of paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short meanSizeSDEstimate()
{
	int64_t i;
	int firstRow_1, firstRow_2, contigID_1, contigID_2, contigPos_1, contigPos_2, left_contigPos, right_contigPos;
	int left_readOrient, right_readOrient;
	int64_t totalPairedNum;
	double quadraticSum;

	// compute the mean size of fragments in library
	totalPairedNum = 0;
	for(i=0; i<readItemNumInSRL; i+=2)
	{
		// Only the paired reads both having just one occurrence will be considered.
		if(sharedReadListArr[i].curNum==1 && sharedReadListArr[i+1].curNum==1)
		{// only the unique reads are taken consideration

			firstRow_1 = sharedReadListArr[i].firstRow;
			firstRow_2 = sharedReadListArr[i+1].firstRow;

			contigID_1 = sharedReadPosArr[firstRow_1].contigID;
			contigID_2 = sharedReadPosArr[firstRow_2].contigID;

			if(contigID_1==contigID_2 && contigID_1>0)
			{ // only the paired reads in the same contig are considered to get the left read and the right read
				contigPos_1 = sharedReadPosArr[firstRow_1].contigPos;
				contigPos_2 = sharedReadPosArr[firstRow_2].contigPos;

				if(contigPos_1<=contigPos_2)
				{
					left_contigPos = contigPos_1;
					right_contigPos = contigPos_2 + readLen - 1;
					left_readOrient = sharedReadPosArr[firstRow_1].orientation;
					right_readOrient = sharedReadPosArr[firstRow_2].orientation;
				}else
				{
					left_contigPos = contigPos_2;
					right_contigPos = contigPos_1 + readLen - 1;
					left_readOrient = sharedReadPosArr[firstRow_2].orientation;
					right_readOrient = sharedReadPosArr[firstRow_1].orientation;
				}

				if(left_readOrient==ORIENTATION_PLUS && right_readOrient==ORIENTATION_MINUS)
				{ // only the inner oriented paired reads are considered to compute the mean size and standard deviation of fragments
					meanSizeInsert += right_contigPos - left_contigPos + 1;
					totalPairedNum ++;
				}
			}
		}
	}

	meanSizeInsert /= totalPairedNum;

	// compute the standard deviation size of fragments in library
	quadraticSum = 0;
	totalPairedNum = 0;
	for(i=0; i<readItemNumInSRL; i+=2)
	{
		// Only the paired reads both having just one occurrence will be considered.
		if(sharedReadListArr[i].curNum==1 && sharedReadListArr[i+1].curNum==1)
		{// only the unique reads are taken consideration

			firstRow_1 = sharedReadListArr[i].firstRow;
			firstRow_2 = sharedReadListArr[i+1].firstRow;

			contigID_1 = sharedReadPosArr[firstRow_1].contigID;
			contigID_2 = sharedReadPosArr[firstRow_2].contigID;

			if(contigID_1==contigID_2 && contigID_1>0)
			{ // only the paired reads in the same contig are considered to get the left read and the right read
				contigPos_1 = sharedReadPosArr[firstRow_1].contigPos;
				contigPos_2 = sharedReadPosArr[firstRow_2].contigPos;

				if(contigPos_1<=contigPos_2)
				{
					left_contigPos = contigPos_1;
					right_contigPos = contigPos_2 + readLen - 1;
					left_readOrient = sharedReadPosArr[firstRow_1].orientation;
					right_readOrient = sharedReadPosArr[firstRow_2].orientation;
				}else
				{
					left_contigPos = contigPos_2;
					right_contigPos = contigPos_1 + readLen - 1;
					left_readOrient = sharedReadPosArr[firstRow_2].orientation;
					right_readOrient = sharedReadPosArr[firstRow_1].orientation;
				}

				if(left_readOrient==ORIENTATION_PLUS && right_readOrient==ORIENTATION_MINUS)
				{ // only the inner oriented paired reads are considered to compute the mean size and standard deviation of fragments
					quadraticSum += (right_contigPos - left_contigPos + 1 - meanSizeInsert) * (right_contigPos - left_contigPos + 1 - meanSizeInsert);
					totalPairedNum ++;
				}
			}
		}
	}

	stardardDeviationInsert = sqrt(quadraticSum / (totalPairedNum - 1));

	//########################### Debug information #######################
	printf("Estimated insert size  : %.2f\n", meanSizeInsert);
	printf("Estimated standard dev : %.2f\n", stardardDeviationInsert);
	printf("Used paired reads num  : %ld\n", totalPairedNum);
	//########################### Debug information #######################

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
short gapSizeEstimateBetweenContigs(contigOverlap *pContigOverlapInfo, int *gapSize, int *validPairedNum, int endCutRound, int *cutOrderArray, int *uncoveredEndLenArray)
{
	int i, hitRow;
	int contigID1, contigID2, contigOrient1, contigOrient2, contigEndFlag1, contigEndFlag2;
	int readsNum1, readsNum2, firstRow1, firstRow2, readOrient1, readOrient2, contigPos1, contigPos2;
	int64_t readID1, readID2;
	int baseNum1, baseNum2, contigLen1, contigLen2;
	ReadList *pReadList1, *pReadList2;
	int tmp_contigID2;
	double totalGapSize;

	contigID1 = pContigOverlapInfo->contigID1;
	contigID2 = pContigOverlapInfo->contigID2;

	contigOrient1 = pContigOverlapInfo->orientation1;
	contigOrient2 = pContigOverlapInfo->orientation2;

	// get the contig end flag
	if(contigInfoArr[contigID1-1].onlyEnd5==YES)
		contigEndFlag1 = 2;
	else
	{
		if(contigOrient1==ORIENTATION_PLUS)
			contigEndFlag1 = 1;
		else
			contigEndFlag1 = 0;
	}

	if(contigInfoArr[contigID2-1].onlyEnd5==YES)
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
		readsNum1 = contigListArr[contigID1-1].EndNum3;
		firstRow1 = contigListArr[contigID1-1].firstRow3;
	}else
	{ // 5' end
		readsNum1 = contigListArr[contigID1-1].EndNum5;
		firstRow1 = contigListArr[contigID1-1].firstRow5;
	}

	if(contigEndFlag2==1)
	{ // 3' end
		readsNum2 = contigListArr[contigID2-1].EndNum3;
		firstRow2 = contigListArr[contigID2-1].firstRow3;
	}else
	{ // 5' end
		readsNum2 = contigListArr[contigID2-1].EndNum5;
		firstRow2 = contigListArr[contigID2-1].firstRow5;
	}

	contigLen1 = contigInfoArr[contigID1-1].contigLen;
	contigLen2 = contigInfoArr[contigID2-1].contigLen;

	// estimate the gap size between contigs
	totalGapSize = 0;
	*validPairedNum = 0;
	for(i=0; i<readsNum1; i++)
	{
		readID1 = contigReadArr[firstRow1+i].readID;
		readOrient1 = contigReadArr[firstRow1+i].orientation;
		contigPos1 = contigReadArr[firstRow1+i].contigPos;

		if(readID1>0)
		{
			hitRow = getReadRowFromReadList(readID1, sharedReadListArr, readItemNumInSRL);
			if(hitRow>=0)
			{
				pReadList1 = sharedReadListArr + hitRow;
				if(pReadList1->curNum==1)
				{
					// get its paired end read
					if(readID1%2==1)
					{ // odd number
						readID2 = readID1 + 1;
					}else
					{ // even number
						readID2 = readID1 - 1;
					}

					hitRow = getReadRowFromReadList(readID2, sharedReadListArr, readItemNumInSRL);
					if(hitRow>=0)
					{
						pReadList2 = sharedReadListArr + hitRow;
						if(pReadList2->curNum==1)
						{
							tmp_contigID2 = sharedReadPosArr[ pReadList2->firstRow ].contigID;
							readOrient2 = sharedReadPosArr[ pReadList2->firstRow ].orientation;
							contigPos2 = sharedReadPosArr[ pReadList2->firstRow ].contigPos;

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
									baseNum2 = contigPos2 + readLen - 1;

									totalGapSize += meanSizeInsert - baseNum1 - baseNum2;
									(*validPairedNum) ++;

								}else if((contigEndFlag1!=0 && contigOrient1==ORIENTATION_PLUS && readOrient1==ORIENTATION_PLUS) && (contigEndFlag2!=0 && contigOrient2==ORIENTATION_MINUS && readOrient2==ORIENTATION_PLUS))
								{ // (3', +, +) , (3', -, +)
									baseNum1 = contigLen1 - contigPos1 + 1;
									baseNum2 = contigLen2 - contigPos2 + 1;

									totalGapSize += meanSizeInsert - baseNum1 - baseNum2;
									(*validPairedNum) ++;

								}else if((contigEndFlag1!=1 && contigOrient1==ORIENTATION_MINUS && readOrient1==ORIENTATION_MINUS) && (contigEndFlag2!=0 && contigOrient2==ORIENTATION_MINUS && readOrient2==ORIENTATION_PLUS))
								{ // (5', -, -) , (3', -, +)
									baseNum1 = contigPos1 + readLen - 1;
									baseNum2 = contigLen2 - contigPos2 + 1;

									totalGapSize += meanSizeInsert - baseNum1 - baseNum2;
									(*validPairedNum) ++;

								}else if((contigEndFlag1!=1 && contigOrient1==ORIENTATION_MINUS && readOrient1==ORIENTATION_MINUS) && (contigEndFlag2!=1 && contigOrient2==ORIENTATION_PLUS && readOrient2==ORIENTATION_MINUS))
								{ // (5', -, -) , (5', +, -)
									baseNum1 = contigPos1 + readLen - 1;
									baseNum2 = contigPos2 + readLen - 1;

									totalGapSize += meanSizeInsert - baseNum1 - baseNum2;
									(*validPairedNum) ++;

								}else
								{
#if DEBUG_OUT_FLAG
									printf("line=%d, In %s(), (contigEndFlag1=%d, contigOrient1=%d, readOrient1=%d), (contigEndFlag2=%d, contigOrient2=%d, readOrient2=%d), readID1=%ld, readID2=%ld\n", __LINE__, __func__, contigEndFlag1, contigOrient1, readOrient1, contigEndFlag2, contigOrient2, readOrient2, readID1, readID2);
#endif
								}
							}
						}
					}else
					{
						printf("line=%d, In %s(), cannot get the read %ld from read list (RL), error!\n", __LINE__, __func__, readID2);
						return FAILED;
					}
				}
			}else
			{
				printf("line=%d, In %s(), cannot get the read %ld from read list (RL), error!\n", __LINE__, __func__, readID1);
				return FAILED;
			}
		}
	}

	if(*validPairedNum>0)
		totalGapSize /= (*validPairedNum);

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

	*gapSize = round(totalGapSize);

#if DEBUG_OUT_FLAG
	printf("Final gapSize=%d, validPairedNum=%d\n", *gapSize, (*validPairedNum));
#endif

	return SUCCESSFUL;
}

/**
 * Get the length of uncovered sub sequences at the ends of the two contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getUncoveredLenAtContigEnds(contigOverlap *pContigOverlapInfo, int *pUncoveredEndLenArray)
{
	int i, hitRow;
	int contigID1, contigID2, contigOrient1, contigOrient2, contigEndFlag1, contigEndFlag2;
	int readsNum1, readsNum2, firstRow1, firstRow2, readOrient1, readOrient2, contigPos1, contigPos2;
	int64_t readID1, readID2;
	int contigLen1, contigLen2;
	ReadList *pReadList1, *pReadList2;
	int tmp_contigID, validReadOrient1, validReadOrient2;
	int cutAllowedContig1, cutAllowedContig2;

	contigID1 = pContigOverlapInfo->contigID1;
	contigID2 = pContigOverlapInfo->contigID2;

	contigOrient1 = pContigOverlapInfo->orientation1;
	contigOrient2 = pContigOverlapInfo->orientation2;

	// get the contig end flag
	if(contigInfoArr[contigID1-1].onlyEnd5==YES)
		contigEndFlag1 = 2;
	else
	{
		if(contigOrient1==ORIENTATION_PLUS)
			contigEndFlag1 = 1;
		else
			contigEndFlag1 = 0;
	}

	if(contigInfoArr[contigID2-1].onlyEnd5==YES)
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
		if(contigInfoArr[contigID1-1].usedTimeEnd3>0)
			cutAllowedContig1 = NO;
	}else
	{ // to cut the 5' end
		if(contigInfoArr[contigID1-1].usedTimeEnd5>0)
			cutAllowedContig1 = NO;
	}

	if(contigOrient2==ORIENTATION_PLUS)
	{ // to cut the 5' end
		if(contigInfoArr[contigID2-1].usedTimeEnd5>0)
			cutAllowedContig2 = NO;
	}else
	{ // to cut the 3' end
		if(contigInfoArr[contigID2-1].usedTimeEnd3>0)
			cutAllowedContig2 = NO;
	}

	if(contigEndFlag1==1)
	{ // 3' end
		readsNum1 = contigListArr[contigID1-1].EndNum3;
		firstRow1 = contigListArr[contigID1-1].firstRow3;
	}else
	{ // 5' end
		readsNum1 = contigListArr[contigID1-1].EndNum5;
		firstRow1 = contigListArr[contigID1-1].firstRow5;
	}

	if(contigEndFlag2==1)
	{ // 3' end
		readsNum2 = contigListArr[contigID2-1].EndNum3;
		firstRow2 = contigListArr[contigID2-1].firstRow3;
	}else
	{ // 5' end
		readsNum2 = contigListArr[contigID2-1].EndNum5;
		firstRow2 = contigListArr[contigID2-1].firstRow5;
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


	// get the number of reads at the end of the first contig
	if(cutAllowedContig1==YES)
	{
		// get the uncovered length at the end of contigs
		pUncoveredEndLenArray[0] = 0;
		if(contigOrient1==ORIENTATION_PLUS)
		{ // (1, +): the plus orientation of the first contig
			for(i=readsNum1-1; i>=0; i--)
			{
				readID1 = contigReadArr[firstRow1+i].readID;
				readOrient1 = contigReadArr[firstRow1+i].orientation;
				contigPos1 = contigReadArr[firstRow1+i].contigPos;
				contigLen1 = contigInfoArr[contigID1-1].contigLen;

				if(readID1>0)
				{
					if(readOrient1==validReadOrient1)
					{
						hitRow = getReadRowFromReadList(readID1, sharedReadListArr, readItemNumInSRL);
						if(hitRow>=0)
						{
							pReadList1 = sharedReadListArr + hitRow;
							if(pReadList1->curNum==1)
							{
								// get its paired end read
								if(readID1%2==1)
								{ // odd number
									readID2 = readID1 + 1;
								}else
								{ // even number
									readID2 = readID1 - 1;
								}

								hitRow = getReadRowFromReadList(readID2, sharedReadListArr, readItemNumInSRL);
								if(hitRow>=0)
								{
									pReadList2 = sharedReadListArr + hitRow;
									if(pReadList2->curNum==1)
									{
										tmp_contigID = sharedReadPosArr[ pReadList2->firstRow ].contigID;
										readOrient2 = sharedReadPosArr[ pReadList2->firstRow ].orientation;
										contigPos2 = sharedReadPosArr[ pReadList2->firstRow ].contigPos;
										contigLen2 = contigInfoArr[contigID2-1].contigLen;

										if(tmp_contigID==contigID2 && readOrient2==validReadOrient2)
										{
											pUncoveredEndLenArray[0] = contigLen1 - contigPos1 - readLen + 1;
											break;
										}
									}
								}else
								{
									printf("line=%d, In %s(), cannot get the read %ld from read list (RL), error!\n", __LINE__, __func__, readID2);
									return FAILED;
								}
							}
						}else
						{
							printf("line=%d, In %s(), cannot get the read %ld from read list (RL), error!\n", __LINE__, __func__, readID1);
							return FAILED;
						}
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
				readID1 = contigReadArr[firstRow1+i].readID;
				readOrient1 = contigReadArr[firstRow1+i].orientation;
				contigPos1 = contigReadArr[firstRow1+i].contigPos;
				contigLen1 = contigInfoArr[contigID1-1].contigLen;

				if(readID1>0)
				{
					if(readOrient1==validReadOrient1)
					{
						hitRow = getReadRowFromReadList(readID1, sharedReadListArr, readItemNumInSRL);
						if(hitRow>=0)
						{
							pReadList1 = sharedReadListArr + hitRow;
							if(pReadList1->curNum==1)
							{
								// get its paired end read
								if(readID1%2==1)
								{ // odd number
									readID2 = readID1 + 1;
								}else
								{ // even number
									readID2 = readID1 - 1;
								}

								hitRow = getReadRowFromReadList(readID2, sharedReadListArr, readItemNumInSRL);
								if(hitRow>=0)
								{
									pReadList2 = sharedReadListArr + hitRow;
									if(pReadList2->curNum==1)
									{
										tmp_contigID = sharedReadPosArr[ pReadList2->firstRow ].contigID;
										readOrient2 = sharedReadPosArr[ pReadList2->firstRow ].orientation;
										contigPos2 = sharedReadPosArr[ pReadList2->firstRow ].contigPos;
										contigLen2 = contigInfoArr[contigID2-1].contigLen;

										if(tmp_contigID==contigID2 && readOrient2==validReadOrient2)
										{
											pUncoveredEndLenArray[0] = contigPos1 - 1;
											break;
										}
									}
								}else
								{
									printf("line=%d, In %s(), cannot get the read %ld from read list (RL), error!\n", __LINE__, __func__, readID2);
									return FAILED;
								}
							}
						}else
						{
							printf("line=%d, In %s(), cannot get the read %ld from read list (RL), error!\n", __LINE__, __func__, readID1);
							return FAILED;
						}
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
				readID2 = contigReadArr[firstRow2+i].readID;
				readOrient2 = contigReadArr[firstRow2+i].orientation;
				contigPos2 = contigReadArr[firstRow2+i].contigPos;
				contigLen2 = contigInfoArr[contigID2-1].contigLen;

				if(readID2>0)
				{
					if(readOrient2==validReadOrient2)
					{
						hitRow = getReadRowFromReadList(readID2, sharedReadListArr, readItemNumInSRL);
						if(hitRow>=0)
						{
							pReadList2 = sharedReadListArr + hitRow;
							if(pReadList2->curNum==1)
							{
								// get its paired end read
								if(readID2%2==1)
								{ // odd number
									readID1 = readID2 + 1;
								}else
								{ // even number
									readID1 = readID2 - 1;
								}

								hitRow = getReadRowFromReadList(readID1, sharedReadListArr, readItemNumInSRL);
								if(hitRow>=0)
								{
									pReadList1 = sharedReadListArr + hitRow;
									if(pReadList1->curNum==1)
									{
										tmp_contigID = sharedReadPosArr[ pReadList1->firstRow ].contigID;
										readOrient1 = sharedReadPosArr[ pReadList1->firstRow ].orientation;
										contigPos1 = sharedReadPosArr[ pReadList1->firstRow ].contigPos;
										contigLen1 = contigInfoArr[contigID1-1].contigLen;

										if(tmp_contigID==contigID1 && readOrient1==validReadOrient1)
										{
											pUncoveredEndLenArray[1] = contigPos2 - 1;
											break;
										}
									}
								}else
								{
									printf("line=%d, In %s(), cannot get the read %ld from read list (RL), error!\n", __LINE__, __func__, readID1);
									return FAILED;
								}
							}
						}else
						{
							printf("line=%d, In %s(), cannot get the read %ld from read list (RL), error!\n", __LINE__, __func__, readID2);
							return FAILED;
						}
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
				readID2 = contigReadArr[firstRow2+i].readID;
				readOrient2 = contigReadArr[firstRow2+i].orientation;
				contigPos2 = contigReadArr[firstRow2+i].contigPos;
				contigLen2 = contigInfoArr[contigID2-1].contigLen;

				if(readID2>0)
				{
					if(readOrient2==validReadOrient2)
					{
						hitRow = getReadRowFromReadList(readID2, sharedReadListArr, readItemNumInSRL);
						if(hitRow>=0)
						{
							pReadList2 = sharedReadListArr + hitRow;
							if(pReadList2->curNum==1)
							{
								// get its paired end read
								if(readID2%2==1)
								{ // odd number
									readID1 = readID2 + 1;
								}else
								{ // even number
									readID1 = readID2 - 1;
								}

								hitRow = getReadRowFromReadList(readID1, sharedReadListArr, readItemNumInSRL);
								if(hitRow>=0)
								{
									pReadList1 = sharedReadListArr + hitRow;
									if(pReadList1->curNum==1)
									{
										tmp_contigID = sharedReadPosArr[ pReadList1->firstRow ].contigID;
										readOrient1 = sharedReadPosArr[ pReadList1->firstRow ].orientation;
										contigPos1 = sharedReadPosArr[ pReadList1->firstRow ].contigPos;
										contigLen1 = contigInfoArr[contigID1-1].contigLen;

										if(tmp_contigID==contigID1 && readOrient1==validReadOrient1)
										{
											pUncoveredEndLenArray[1] = contigLen2 - contigPos2 - readLen + 1;
											break;
										}
									}
								}else
								{
									printf("line=%d, In %s(), cannot get the read %ld from read list (RL), error!\n", __LINE__, __func__, readID1);
									return FAILED;
								}
							}
						}else
						{
							printf("line=%d, In %s(), cannot get the read %ld from read list (RL), error!\n", __LINE__, __func__, readID2);
							return FAILED;
						}
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
short updateOverlapLenByCutUncoveredContigEnds(contigOverlap *pContigOverlapInfo, contigInfo *pContigInfoArr)
{
	int contigID1, orientation1, contigLen1, contigID2, orientation2, contigLen2;
	int seq_len1, seq_len2, overlapLenExact, overlapLenAlignment, overlapLenAdjust, mismatchNum;
	int gapSize, validPairedNum_gapEstimate;

	int uncoveredEndLenArray[2], cutOrderArray[2];
	int endCutRound, maxEndCutRoundNum;
	int contigLenBeforeCut[2];

	// get the information of contg1 and contig2
	contigID1 = pContigOverlapInfo->contigID1;
	orientation1 = pContigOverlapInfo->orientation1;
	contigLen1 = pContigInfoArr[contigID1-1].contigLen;
	contigID2 = pContigOverlapInfo->contigID2;
	orientation2 = pContigOverlapInfo->orientation2;
	contigLen2 = pContigInfoArr[contigID2-1].contigLen;

	contigLenBeforeCut[0] = contigLen1;
	contigLenBeforeCut[1] = contigLen2;

	// get the uncovered region length at contig ends
	if(getUncoveredLenAtContigEnds(pContigOverlapInfo, uncoveredEndLenArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot check the uncovered length at the ends of contigs [%d,%d], error!\n", __LINE__, __func__, contigID1, contigID2);
		return FAILED;
	}

	if((uncoveredEndLenArray[0]==0 && uncoveredEndLenArray[1]==0) || (contigLen1-uncoveredEndLenArray[0]<=0 && contigLen2-uncoveredEndLenArray[1]<=0))
	{ // No uncovered regions at contig ends, then break the contig links, and return
		// break the links
#if (DEBUG_FLAG==YES)
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
					strncpy(overlapSeq1, pContigInfoArr[contigID1-1].contigSeq+contigLen1-uncoveredEndLenArray[0]-seq_len1, seq_len1);
					overlapSeq1[seq_len1] = '\0';
				}else
				{
					strncpy(overlapSeq1, pContigInfoArr[contigID1-1].contigSeq+uncoveredEndLenArray[0], seq_len1);
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
					strncpy(overlapSeq2, pContigInfoArr[contigID2-1].contigSeq, seq_len2);
					overlapSeq2[seq_len2] = '\0';
				}else
				{
					strncpy(overlapSeq2, pContigInfoArr[contigID2-1].contigSeq+contigLen2-seq_len2, seq_len2);
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
					strncpy(overlapSeq1, pContigInfoArr[contigID1-1].contigSeq+contigLen1-seq_len1, seq_len1);
					overlapSeq1[seq_len1] = '\0';
				}else
				{
					strncpy(overlapSeq1, pContigInfoArr[contigID1-1].contigSeq, seq_len1);
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
					strncpy(overlapSeq2, pContigInfoArr[contigID2-1].contigSeq+uncoveredEndLenArray[1], seq_len2);
					overlapSeq2[seq_len2] = '\0';
				}else
				{
					strncpy(overlapSeq2, pContigInfoArr[contigID2-1].contigSeq+contigLen2-uncoveredEndLenArray[1]-seq_len2, seq_len2);
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
				strncpy(overlapSeq1, pContigInfoArr[contigID1-1].contigSeq+contigLen1-uncoveredEndLenArray[0]-seq_len1, seq_len1);
				overlapSeq1[seq_len1] = '\0';
			}else
			{
				strncpy(overlapSeq1, pContigInfoArr[contigID1-1].contigSeq+uncoveredEndLenArray[0], seq_len1);
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
				strncpy(overlapSeq2, pContigInfoArr[contigID2-1].contigSeq+uncoveredEndLenArray[1], seq_len2);
				overlapSeq2[seq_len2] = '\0';
			}else
			{
				strncpy(overlapSeq2, pContigInfoArr[contigID2-1].contigSeq+contigLen2-uncoveredEndLenArray[1]-seq_len2, seq_len2);
				overlapSeq2[seq_len2] = '\0';
				if(reverseSeq(overlapSeq2, seq_len2)==FAILED) // reverse the sequence
				{
					printf("line=%d, In %s(), cannot reverse the sequence [ %s ], error!\n", __LINE__, __func__, overlapSeq2);
					return FAILED;
				}
			}
		}

		//############################# Debug information #########################
#if DEBUG_OUT_FLAG
		printf("overlapSeq1=%s, len=%d\noverlapSeq2=%s, len=%d\n", overlapSeq1, (int)strlen(overlapSeq1), overlapSeq2, (int)strlen(overlapSeq2));
#endif
		//############################# Debug information #########################


		// estimate gap size between contigs
		if(gapSizeEstimateBetweenContigs(pContigOverlapInfo, &gapSize, &validPairedNum_gapEstimate, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
		{
			printf("line=%d, In %s(), cannot estimate the gap size between contigs [%d,%d], error!\n", __LINE__, __func__, contigID1, contigID2);
			return FAILED;
		}

		if(validPairedNum_gapEstimate<minLinksNumContigsThres)
		{
			continue;
		}

		//if(gapSize<=stardardDeviationInsert)
		if(gapSize<=gapSizeSdevFactorOverlap*stardardDeviationInsert)
		{
			// compute the overlap length by exact alignment
			if(computeSeqOverlapLenExact(&overlapLenExact, overlapSeq1, seq_len1, overlapSeq2, seq_len2, scoreArr, gapSize)==FAILED)
			{
				printf("line=%d, In %s(), cannot compute sequence overlap length, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// update the contig overlap information
			if(overlapLenExact>=minOverlapThres && (overlapLenExact>=-gapSize-stardardDeviationInsert && overlapLenExact<=-gapSize+stardardDeviationInsert))
			{ // update the overlap information, and break the iteration
				if(cutContigEnds(pContigOverlapInfo, pContigInfoArr, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
				{
					printf("line=%d, In %s(), cannot cut the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
					return FAILED;
				}

				pContigOverlapInfo->mergeFlag = YES;
				pContigOverlapInfo->overlapLen = overlapLenExact;

				if(updateContigEndsInfo(pContigOverlapInfo, pContigInfoArr)==FAILED)
				{
					printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
					return FAILED;
				}

				//==============================================================
				// update the read list 1 (RL1), read list 2 (RL2), shared read list (SRL) and contig list (CL)
				if(updateReadListsAndContigListAfterCut(pContigOverlapInfo, pContigInfoArr, endCutRound, cutOrderArray, contigLenBeforeCut)==FAILED)
				{
					printf("line=%d, In %s(), cannot update reads at the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
					return FAILED;
				}
				//==============================================================

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
#if DEBUG_OUT_FLAG
						//printf("line=%d, In %s(), before updateContigs(), contigLen1=%d, contigLen2=%d\n", __LINE__, __func__, pContigInfoArr[contigID1-1].contigLen, pContigInfoArr[contigID2-1].contigLen);
#endif
						//############################### Debug information #########################

						if(cutContigEnds(pContigOverlapInfo, pContigInfoArr, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
						{
							printf("line=%d, In %s(), cannot cut the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
							return FAILED;
						}

						if(updateContigs(pContigInfoArr+contigID1-1, pContigInfoArr+contigID2-1, orientation1, orientation2, overlapSeq1, seq_len1, overlapSeq2, seq_len2)==FAILED)
						{
							printf("line=%d, In %s(), cannot update contig end sequences, error!\n", __LINE__, __func__);
							return FAILED;
						}

						//############################### Debug information #########################
#if DEBUG_OUT_FLAG
						//printf("line=%d, In %s(), after updateContigs(), contigLen1=%d, contigLen2=%d\n", __LINE__, __func__, pContigInfoArr[contigID1-1].contigLen, pContigInfoArr[contigID2-1].contigLen);
#endif
						//############################### Debug information #########################

						pContigOverlapInfo->mergeFlag = YES;
						pContigOverlapInfo->overlapLen = overlapLenAdjust;

						if(updateContigEndsInfo(pContigOverlapInfo, pContigInfoArr)==FAILED)
						{
							printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
							return FAILED;
						}

						//==============================================================
						// update the read list 1 (RL1), read list 2 (RL2), shared read list (SRL) and contig list (CL)
						if(updateReadListsAndContigListAfterCut(pContigOverlapInfo, pContigInfoArr, endCutRound, cutOrderArray, contigLenBeforeCut)==FAILED)
						{
							printf("line=%d, In %s(), cannot update reads at the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
							return FAILED;
						}
						//==============================================================

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

								if(cutContigEnds(pContigOverlapInfo, pContigInfoArr, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
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

								if(updateContigEndsInfo(pContigOverlapInfo, pContigInfoArr)==FAILED)
								{
									printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
									return FAILED;
								}

								//==============================================================
								// update the read list 1 (RL1), read list 2 (RL2), shared read list (SRL) and contig list (CL)
								if(updateReadListsAndContigListAfterCut(pContigOverlapInfo, pContigInfoArr, endCutRound, cutOrderArray, contigLenBeforeCut)==FAILED)
								{
									printf("line=%d, In %s(), cannot update reads at the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
									return FAILED;
								}
								//==============================================================

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

									if(cutContigEnds(pContigOverlapInfo, pContigInfoArr, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
									{
										printf("line=%d, In %s(), cannot cut the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
										return FAILED;
									}

									if(updateContigs(pContigInfoArr+contigID1-1, pContigInfoArr+contigID2-1, orientation1, orientation2, overlapSeq1, seq_len1, overlapSeq2, seq_len2)==FAILED)
									{
										printf("line=%d, In %s(), cannot update contig end sequences, error!\n", __LINE__, __func__);
										return FAILED;
									}
									pContigOverlapInfo->mergeFlag = YES;
									pContigOverlapInfo->overlapLen = overlapLenAdjust;

									if(updateContigEndsInfo(pContigOverlapInfo, pContigInfoArr)==FAILED)
									{
										printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
										return FAILED;
									}

									//==============================================================
									// update the read list 1 (RL1), read list 2 (RL2), shared read list (SRL) and contig list (CL)
									if(updateReadListsAndContigListAfterCut(pContigOverlapInfo, pContigInfoArr, endCutRound, cutOrderArray, contigLenBeforeCut)==FAILED)
									{
										printf("line=%d, In %s(), cannot update reads at the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
										return FAILED;
									}
									//==============================================================

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
								//if((pContigInfoArr[contigID1-1].onlyEnd5==YES || pContigInfoArr[contigID2-1].onlyEnd5==YES) && validPairedNum_gapEstimate<breakLinkNumThres)
								if((pContigInfoArr[contigID1-1].shortFlag==YES || pContigInfoArr[contigID2-1].shortFlag==YES) && validPairedNum_gapEstimate<breakLinkNumThres)
								{
#if (DEBUG_FLAG==YES)
		printf("line=%d, In %s(), broken.\n", __LINE__, __func__);
#endif
									pContigOverlapInfo->breakFlag = YES;
								}else
								{
									if(cutContigEnds(pContigOverlapInfo, pContigInfoArr, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
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

									if(updateContigEndsInfo(pContigOverlapInfo, pContigInfoArr)==FAILED)
									{
										printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
										return FAILED;
									}

									//==============================================================
									// update the read list 1 (RL1), read list 2 (RL2), shared read list (SRL) and contig list (CL)
									if(updateReadListsAndContigListAfterCut(pContigOverlapInfo, pContigInfoArr, endCutRound, cutOrderArray, contigLenBeforeCut)==FAILED)
									{
										printf("line=%d, In %s(), cannot update reads at the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
										return FAILED;
									}
									//==============================================================

								}

								return SUCCESSFUL;

							}else
							{ // gapSize >= 10
								pContigOverlapInfo->mergeFlag = NO;
								pContigOverlapInfo->gapSize = gapSize;

								if(cutContigEnds(pContigOverlapInfo, pContigInfoArr, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
								{
									printf("line=%d, In %s(), cannot cut the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
									return FAILED;
								}

								if(updateContigEndsInfo(pContigOverlapInfo, pContigInfoArr)==FAILED)
								{
									printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
									return FAILED;
								}

								//==============================================================
								// update the read list 1 (RL1), read list 2 (RL2), shared read list (SRL) and contig list (CL)
								if(updateReadListsAndContigListAfterCut(pContigOverlapInfo, pContigInfoArr, endCutRound, cutOrderArray, contigLenBeforeCut)==FAILED)
								{
									printf("line=%d, In %s(), cannot update reads at the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
									return FAILED;
								}
								//==============================================================

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

			if(cutContigEnds(pContigOverlapInfo, pContigInfoArr, endCutRound, cutOrderArray, uncoveredEndLenArray)==FAILED)
			{
				printf("line=%d, In %s(), cannot cut the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
				return FAILED;
			}

			if(updateContigEndsInfo(pContigOverlapInfo, pContigInfoArr)==FAILED)
			{
				printf("line=%d, In %s(), cannot update contig end info, error!\n", __LINE__, __func__);
				return FAILED;
			}

			//==============================================================
			// update the read list 1 (RL1), read list 2 (RL2), shared read list (SRL) and contig list (CL)
			if(updateReadListsAndContigListAfterCut(pContigOverlapInfo, pContigInfoArr, endCutRound, cutOrderArray, contigLenBeforeCut)==FAILED)
			{
				printf("line=%d, In %s(), cannot update reads at the ends of contigs: [%d, %d], error!\n", __LINE__, __func__, pContigOverlapInfo->contigID1, pContigOverlapInfo->contigID2);
				return FAILED;
			}
			//==============================================================

			return SUCCESSFUL;
		}
	} // end for(endCutRound=0; endCutRound<maxEndCutRoundNum; endCutRound++)

	// No valid overlaps, break the links without cutting contig ends
	//pContigOverlapInfo->breakFlag = YES;

	return SUCCESSFUL;
}

/**
 * Cut contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short cutContigEnds(contigOverlap *pContigOverlapInfo, contigInfo *pContigInfoArr, int endCutRound, int *cutOrderArray, int *uncoveredEndLenArray)
{
	int contigID1, contigID2, contigOrient1, contigOrient2;
	int contigLen1, contigLen2;
	char *tmpStr;

	contigID1 = pContigOverlapInfo->contigID1;
	contigID2 = pContigOverlapInfo->contigID2;
	contigOrient1 = pContigOverlapInfo->orientation1;
	contigOrient2 = pContigOverlapInfo->orientation2;

	contigLen1 = pContigInfoArr[contigID1-1].contigLen;
	contigLen2 = pContigInfoArr[contigID2-1].contigLen;

	// update contigs information
	if(endCutRound==0 || endCutRound==1)
	{ // the first round or the second round
		if(cutOrderArray[endCutRound]==0)
		{ // cut the end of the first contig
			if(contigOrient1==ORIENTATION_PLUS)
			{ // the plus orientation, cut the 3' end
				pContigInfoArr[contigID1-1].contigLen -= uncoveredEndLenArray[0];
				pContigInfoArr[contigID1-1].contigSeq[ pContigInfoArr[contigID1-1].contigLen ] = '\0';
			}else
			{ // the minus orientation, cut the 5' end
				tmpStr = (char *) calloc(contigLen1+1, sizeof(char));
				if(tmpStr==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// copy the remained bases to tmpStr
				strcpy(tmpStr, pContigInfoArr[contigID1-1].contigSeq+uncoveredEndLenArray[0]);
				// copy back
				strcpy(pContigInfoArr[contigID1-1].contigSeq, tmpStr);
				pContigInfoArr[contigID1-1].contigLen -= uncoveredEndLenArray[0];

				// free tmpStr
				free(tmpStr);
				tmpStr = NULL;
			}

			// ########################### Debug information ####################
#if DEBUG_FLAG
			if(strlen(pContigInfoArr[contigID1-1].contigSeq)!=pContigInfoArr[contigID1-1].contigLen)
			{
				printf("line=%d, In %s(), contigID1=%d, contigSeq=%s, strlen=%d, contigLen=%d, error!\n", __LINE__, __func__, contigID1, pContigInfoArr[contigID1-1].contigSeq, (int)strlen(pContigInfoArr[contigID1-1].contigSeq), pContigInfoArr[contigID1-1].contigLen);
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
				strcpy(tmpStr, pContigInfoArr[contigID2-1].contigSeq+uncoveredEndLenArray[1]);
				// copy back
				strcpy(pContigInfoArr[contigID2-1].contigSeq, tmpStr);
				pContigInfoArr[contigID2-1].contigLen -= uncoveredEndLenArray[1];

				// free tmpStr
				free(tmpStr);
				tmpStr = NULL;
			}else
			{ // the minus orientation, cut the 3' end
				pContigInfoArr[contigID2-1].contigLen -= uncoveredEndLenArray[1];
				pContigInfoArr[contigID2-1].contigSeq[ pContigInfoArr[contigID2-1].contigLen ] = '\0';
			}

			// ########################### Debug information ####################
#if DEBUG_FLAG
			if(strlen(pContigInfoArr[contigID2-1].contigSeq)!=pContigInfoArr[contigID2-1].contigLen)
			{
				printf("line=%d, In %s(), contigID2=%d, contigSeq=%s, strlen=%d, contigLen=%d, error!\n", __LINE__, __func__, contigID2, pContigInfoArr[contigID2-1].contigSeq, (int)strlen(pContigInfoArr[contigID2-1].contigSeq), pContigInfoArr[contigID2-1].contigLen);
				return FAILED;
			}
#endif
			// ########################### Debug information ####################
		}
	}else
	{ // the third round
		if(contigOrient1==ORIENTATION_PLUS && contigOrient2==ORIENTATION_MINUS)
		{ // A(+), B(-)
			pContigInfoArr[contigID1-1].contigLen -= uncoveredEndLenArray[0];
			pContigInfoArr[contigID1-1].contigSeq[ pContigInfoArr[contigID1-1].contigLen ] = '\0';

			pContigInfoArr[contigID2-1].contigLen -= uncoveredEndLenArray[1];
			pContigInfoArr[contigID2-1].contigSeq[ pContigInfoArr[contigID2-1].contigLen ] = '\0';
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
				pContigInfoArr[contigID1-1].contigLen -= uncoveredEndLenArray[0];
				pContigInfoArr[contigID1-1].contigSeq[ pContigInfoArr[contigID1-1].contigLen ] = '\0';
			}else
			{ // A(-)
				// copy the remained bases to tmpStr
				strcpy(tmpStr, pContigInfoArr[contigID1-1].contigSeq+uncoveredEndLenArray[0]);
				// copy back
				strcpy(pContigInfoArr[contigID1-1].contigSeq, tmpStr);

				pContigInfoArr[contigID1-1].contigLen -= uncoveredEndLenArray[0];
			}

			if(contigOrient2==ORIENTATION_PLUS)
			{ // B(+)
				// copy the remained bases to tmpStr
				strcpy(tmpStr, pContigInfoArr[contigID2-1].contigSeq+uncoveredEndLenArray[1]);
				// copy back
				strcpy(pContigInfoArr[contigID2-1].contigSeq, tmpStr);

				pContigInfoArr[contigID2-1].contigLen -= uncoveredEndLenArray[1];
			}else
			{ // B(-)
				pContigInfoArr[contigID2-1].contigLen -= uncoveredEndLenArray[1];
				pContigInfoArr[contigID2-1].contigSeq[ pContigInfoArr[contigID2-1].contigLen ] = '\0';
			}

			free(tmpStr);
		}

		// ########################### Debug information ####################
#if DEBUG_FLAG
		if(strlen(pContigInfoArr[contigID1-1].contigSeq)!=pContigInfoArr[contigID1-1].contigLen)
		{
			printf("line=%d, In %s(), contigID1=%d, contigSeq=%s, strlen=%d, contigLen=%d, error!\n", __LINE__, __func__, contigID1, pContigInfoArr[contigID1-1].contigSeq, (int)strlen(pContigInfoArr[contigID1-1].contigSeq), pContigInfoArr[contigID1-1].contigLen);
			return FAILED;
		}

		if(strlen(pContigInfoArr[contigID2-1].contigSeq)!=pContigInfoArr[contigID2-1].contigLen)
		{
			printf("line=%d, In %s(), contigID2=%d, contigSeq=%s, strlen=%d, contigLen=%d, error!\n", __LINE__, __func__, contigID2, pContigInfoArr[contigID2-1].contigSeq, (int)strlen(pContigInfoArr[contigID2-1].contigSeq), pContigInfoArr[contigID2-1].contigLen);
			return FAILED;
		}
#endif
		// ########################### Debug information ####################
	}

	return SUCCESSFUL;
}

/**
 * Update the read list 1 (RL1), read list 2 (RL2), shared read list (SRL) and contig list (CL)
 * after contig end alignment.
 *
 *  Updates:
 *  	(1) remove invalid reads in read list (RL);
 *  	(2) update contigPos of a read in read list (RL);
 *  	(3) remove invalid reads in contig list (CL);
 *  	(4) update contig reads information in contig list (CL).
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateReadListsAndContigListSingleContigAfterAlignment(contigInfo *pContigInfoArr, int contigOrient, int contigIndex, int contigLenBeforeAlignment)
{
	uint32_t i;
	uint32_t contigID, contigEnd;
	uint32_t firstRow, endNum;
	ContigRead *pContigReadInfo;
	int diff;

	if(contigIndex == 1)
	{ // the first contig
		if(contigOrient == ORIENTATION_PLUS)
			contigEnd = 1;
		else
			contigEnd = 0;
	}
	else //if(contigIndex == 2)
	{ // the second contig
		if(contigOrient == ORIENTATION_PLUS)
			contigEnd = 0;
		else
			contigEnd = 1;
	}


	contigID = pContigInfoArr->contigID;
	diff = contigLenBeforeAlignment - pContigInfoArr->contigLen;

	if(contigEnd == 0) // 5' end
	{
		//update 5' end
		firstRow = contigListArr[contigID-1].firstRow5;
		endNum = contigListArr[contigID-1].EndNum5;
		pContigReadInfo = contigReadArr + firstRow;

		for(i=0; i<endNum; i++)
		{
			if(diff > 0)
			{ //the contig becomes shorter than before
				if(pContigReadInfo[i].contigPos <= diff && pContigReadInfo[i].readID != 0)
				{
					if(removeReadsFromReadLists(pContigReadInfo[i].readID, contigID, pContigReadInfo[i].contigPos) == FAILED)
					{
						printf("line=%d, In %s(), cannot remove reads from read lists, error!\n", __LINE__, __func__);
						return FAILED;
					}
					pContigReadInfo[i].readID = 0;
					contigListArr[contigID-1].curNum5--;
				}
			}

			if(pContigReadInfo[i].readID != 0)
			{
				if(updateReadsInReadLists(pContigReadInfo[i].readID, contigID, pContigReadInfo[i].contigPos, diff) == FAILED)
				{
					printf("line=%d, In %s(), cannot update reads from read lists, error!\n", __LINE__, __func__);
					return FAILED;
				}
				pContigReadInfo[i].contigPos -= diff;
			}
		}

		//update 3' end of the first contig
		if(pContigInfoArr->onlyEnd5==NO)
		{
			firstRow = contigListArr[contigID-1].firstRow3;
			endNum = contigListArr[contigID-1].EndNum3;
			pContigReadInfo = contigReadArr + firstRow;
			for(i=0; i<endNum; i++)
			{
				if(pContigReadInfo[i].readID != 0)
				{
					if(updateReadsInReadLists(pContigReadInfo[i].readID, contigID, pContigReadInfo[i].contigPos, diff) == FAILED)
					{
						printf("line=%d, In %s(), cannot update reads from read lists, error!\n", __LINE__, __func__);
						return FAILED;
					}
					pContigReadInfo[i].contigPos -= diff;
				}
			}
		}
	}else if(contigEnd == 1) // 3' end
	{
		if(diff > 0)
		{ //the contig becomes shorter than before
			if(pContigInfoArr->onlyEnd5==NO)
			{
				firstRow = contigListArr[contigID-1].firstRow3;
				endNum = contigListArr[contigID-1].EndNum3;
				pContigReadInfo = contigReadArr + firstRow;

				for(i=endNum-1; i>=0; i--)
				{
					if(pContigReadInfo[i].contigPos+readLen-1 > pContigInfoArr->contigLen)
					{
						if(pContigReadInfo[i].readID != 0)
						{
							if(removeReadsFromReadLists(pContigReadInfo[i].readID, contigID, pContigReadInfo[i].contigPos) == FAILED)
							{
								printf("line=%d, In %s(), cannot remove reads from read lists, error!\n", __LINE__, __func__);
								return FAILED;
							}
							pContigReadInfo[i].readID = 0;
							contigListArr[contigID-1].curNum3--;
						}
					}else
					{
						break;
					}
				}
			}else
			{
				firstRow = contigListArr[contigID-1].firstRow5;
				endNum = contigListArr[contigID-1].EndNum5;
				pContigReadInfo = contigReadArr + firstRow;

				for(i=endNum-1; i>=0; i--)
				{
					if(pContigReadInfo[i].contigPos+readLen-1 > pContigInfoArr->contigLen)
					{
						if(pContigReadInfo[i].readID != 0)
						{
							if(removeReadsFromReadLists(pContigReadInfo[i].readID, contigID, pContigReadInfo[i].contigPos) == FAILED)
							{
								printf("line=%d, In %s(), cannot remove reads from read lists, error!\n", __LINE__, __func__);
								return FAILED;
							}
							pContigReadInfo[i].readID = 0;
							contigListArr[contigID-1].curNum5--;
						}
					}else
					{
						break;
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Update the read list 1 (RL1), read list 2 (RL2), shared read list (SRL) and contig list (CL)
 * after cutting contig ends.
 *
 *  Updates:
 *  	(1) remove invalid reads in read list (RL);
 *  	(2) update contigPos of a read in read list (RL);
 *  	(3) remove invalid reads in contig list (CL);
 *  	(4) update contig reads information in contig list (CL).
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateReadListsAndContigListAfterCut(contigOverlap *pContigOverlapInfo, contigInfo *pContigInfoArr, int endCutRound, int *cutOrderArray, int *contigLenBeforeCut)
{
	int i;
	int contigID1, contigID2, contigOrient1, contigOrient2;
	int contigLen1, contigLen2;
	int readsNum, contigPosDivision;
	ContigRead *pContigReadInfo;
	uint64_t readID;

	contigID1 = pContigOverlapInfo->contigID1;
	contigID2 = pContigOverlapInfo->contigID2;
	contigOrient1 = pContigOverlapInfo->orientation1;
	contigOrient2 = pContigOverlapInfo->orientation2;

	contigLen1 = pContigInfoArr[contigID1-1].contigLen;
	contigLen2 = pContigInfoArr[contigID2-1].contigLen;


	// update read lists and contig list
	if(endCutRound==2 || ((endCutRound==0 || endCutRound==1) && cutOrderArray[endCutRound]==0))
	{ // get the cut end of the first contig
		if(contigOrient1==ORIENTATION_PLUS)
		{ // the plus orientation, the 3' end has been cut

			contigPosDivision = contigLen1; // the valid last contigPos

			// get the reads in contig list
			if(pContigInfoArr[contigID1-1].onlyEnd5==YES)
			{
				readsNum = contigListArr[contigID1-1].EndNum5;
				pContigReadInfo = contigReadArr + contigListArr[contigID1-1].firstRow5;
			}else
			{
				readsNum = contigListArr[contigID1-1].EndNum3;
				pContigReadInfo = contigReadArr + contigListArr[contigID1-1].firstRow3;
			}

			for(i=readsNum-1; i>=0; i--)
			{
				if(pContigReadInfo[i].contigPos+readLen-1 > contigPosDivision && pContigReadInfo[i].readID!=0)
				{
					// remove the reads of this contig in read lists
					// 	take care that the matchNum was not accurate now,
					// 	it must be recompute before outputting the read lists to file
					readID = pContigReadInfo[i].readID;
					if(removeReadsFromReadLists(readID, contigID1, pContigReadInfo[i].contigPos)==FAILED)
					{
						printf("line=%d, In %s(), cannot remove the read [ %ld ] from RLs, error!\n", __LINE__, __func__, readID);
						return FAILED;
					}

					// remove the read from contig list
					pContigReadInfo[i].readID = 0;
					contigListArr[contigID1-1].curNum3 --;
				}else
				{
					break;
				}
			}
		}else
		{ // the minus orientation, the 5' end has been cut

			contigPosDivision = contigLenBeforeCut[0] - contigLen1 + 1; // the valid last contigPos

			// get the reads in contig list
			readsNum = contigListArr[contigID1-1].EndNum5;
			pContigReadInfo = contigReadArr + contigListArr[contigID1-1].firstRow5;

			for(i=0; i<readsNum; i++)
			{
				readID = pContigReadInfo[i].readID;
				if(pContigReadInfo[i].contigPos < contigPosDivision && pContigReadInfo[i].readID!=0)
				{
					// remove the reads of this contig in read lists
					// 	take care that the matchNum was not accurate now,
					// 	it must be recompute before outputting the read lists to file
					if(removeReadsFromReadLists(readID, contigID1, pContigReadInfo[i].contigPos)==FAILED)
					{
						printf("line=%d, In %s(), cannot remove the read [ %ld ] from RLs, error!\n", __LINE__, __func__, readID);
						return FAILED;
					}

					// remove the read from contig list
					pContigReadInfo[i].readID = 0;
					contigListArr[contigID1-1].curNum5 --;
				}else
				{
					// update the contigPos in read lists
					if(updateReadsInReadLists(readID, contigID1, pContigReadInfo[i].contigPos, contigLenBeforeCut[0]-contigLen1)==FAILED)
					{
						printf("line=%d, In %s(), cannot remove the read [ %ld ] in RLs, error!\n", __LINE__, __func__, readID);
						return FAILED;
					}

					// update the contigPos at 5' end in contig list
					pContigReadInfo[i].contigPos -= contigLenBeforeCut[0] - contigLen1;
				}
			}

			// update the reads at the 3' end of the first contig
			readsNum = contigListArr[contigID1-1].EndNum3;
			pContigReadInfo = contigReadArr + contigListArr[contigID1-1].firstRow3;

			for(i=0; i<readsNum; i++)
			{
				// update the contigPos at 3' end in contig list
				pContigReadInfo[i].contigPos -= contigLenBeforeCut[0] - contigLen1;
			}
		}
	}

	if(endCutRound==2 || ((endCutRound==0 || endCutRound==1) && cutOrderArray[endCutRound]==1))
	{ // get the cut end of the second contig
		if(contigOrient2==ORIENTATION_PLUS)
		{ // the plus orientation, the 5' end has been cut
			contigPosDivision = contigLenBeforeCut[1] - contigLen2 + 1; // the valid last contigPos

			// get the reads in contig list
			readsNum = contigListArr[contigID2-1].EndNum5;
			pContigReadInfo = contigReadArr + contigListArr[contigID2-1].firstRow5;

			for(i=0; i<readsNum; i++)
			{
				readID = pContigReadInfo[i].readID;
				if(pContigReadInfo[i].contigPos < contigPosDivision && pContigReadInfo[i].readID!=0)
				{
					// remove the reads of this contig in read lists
					// 	take care that the matchNum was not accurate now,
					// 	it must be recompute before outputting the read lists to file
					if(removeReadsFromReadLists(readID, contigID2, pContigReadInfo[i].contigPos)==FAILED)
					{
						printf("line=%d, In %s(), cannot remove the read [ %ld ] from RLs, error!\n", __LINE__, __func__, readID);
						return FAILED;
					}

					// remove the read from contig list
					pContigReadInfo[i].readID = 0;
					contigListArr[contigID2-1].curNum5 --;
				}else
				{
					// update the contigPos in read lists
					if(updateReadsInReadLists(readID, contigID2, pContigReadInfo[i].contigPos, contigLenBeforeCut[1]-contigLen2)==FAILED)
					{
						printf("line=%d, In %s(), cannot remove the read [ %ld ] in RLs, error!\n", __LINE__, __func__, readID);
						return FAILED;
					}

					// update the contigPos at 5' end in contig list
					pContigReadInfo[i].contigPos -= contigLenBeforeCut[1] - contigLen2;
				}
			}

			// update the reads at the 3' end of the second contig
			readsNum = contigListArr[contigID2-1].EndNum3;
			pContigReadInfo = contigReadArr + contigListArr[contigID2-1].firstRow3;

			for(i=0; i<readsNum; i++)
			{
				// update the contigPos at 3' end in contig list
				pContigReadInfo[i].contigPos -= contigLenBeforeCut[1] - contigLen2;
			}

		}else
		{ // the minus orientation, the 3' end has been cut

			contigPosDivision = contigLen2; // the valid last contigPos

			// get the reads in contig list
			if(pContigInfoArr[contigID2-1].onlyEnd5==YES)
			{
				readsNum = contigListArr[contigID2-1].EndNum5;
				pContigReadInfo = contigReadArr + contigListArr[contigID2-1].firstRow5;
			}else
			{
				readsNum = contigListArr[contigID2-1].EndNum3;
				pContigReadInfo = contigReadArr + contigListArr[contigID2-1].firstRow3;
			}

			for(i=readsNum-1; i>=0; i--)
			{
				if(pContigReadInfo[i].contigPos+readLen-1 > contigPosDivision && pContigReadInfo[i].readID!=0)
				{
					// remove the reads of this contig in read lists
					// 	take care that the matchNum was not accurate now,
					// 	it must be recompute before outputting the read lists to file
					readID = pContigReadInfo[i].readID;
					if(removeReadsFromReadLists(readID, contigID2, pContigReadInfo[i].contigPos)==FAILED)
					{
						printf("line=%d, In %s(), cannot remove the read [ %ld ] from RLs, error!\n", __LINE__, __func__, readID);
						return FAILED;
					}

					// remove the read from contig list
					pContigReadInfo[i].readID = 0;
					contigListArr[contigID2-1].curNum3 --;
				}else
				{
					break;
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Remove reads from read lists (RLs), including RL1, RL2 and SRL.
 *
 *  Take care:
 *  	The matchNum was not accurate after reads removal,
 *  	it must be recompute before outputting the read lists to file.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeReadsFromReadLists(uint64_t readID, int contigIDFrom, int contigPosFrom)
{
	int j, hitRow, tmp_matchNum;
	ReadList *pReadList;
	ReadPos *pReadPos;

	//readID = pContigReadInfo[i].readID;
	if(readID%2==1)
	{ // odd number
		hitRow = getReadRowFromReadList(readID, readListArr[1], readNumInRL[1]);
		if(hitRow>=0)
		{
			//get the match information in read lists
			pReadList = readListArr[1] + hitRow;
			pReadPos = readPosArr[1] + pReadList->firstRow;
			tmp_matchNum = pReadList->matchNum;
			for(j=0; j<tmp_matchNum; j++)
			{
				if(pReadPos->contigID==contigIDFrom && pReadPos->contigPos==contigPosFrom)
				{
					pReadPos->contigID = 0;
					pReadList->curNum --;
				}
			}
		}else
		{ // errors
			printf("line=%d, In %s(), cannot get the read [ %ld ] from RL1, error!\n", __LINE__, __func__, readID);
			return FAILED;
		}
	}else
	{ // even number
		hitRow = getReadRowFromReadList(readID, readListArr[2], readNumInRL[2]);
		if(hitRow>=0)
		{
			//get the match information in read lists
			pReadList = readListArr[2] + hitRow;
			pReadPos = readPosArr[2] + pReadList->firstRow;
			tmp_matchNum = pReadList->matchNum;
			for(j=0; j<tmp_matchNum; j++)
			{
				if(pReadPos->contigID==contigIDFrom && pReadPos->contigPos==contigPosFrom)
				{
					pReadPos->contigID = 0;
					pReadList->curNum --;
				}
			}
		}else
		{ // errors
			printf("line=%d, In %s(), cannot get the read [ %ld ] from RL2, error!\n", __LINE__, __func__, readID);
			return FAILED;
		}
	}

	hitRow = getReadRowFromReadList(readID, sharedReadListArr, readItemNumInSRL);
	if(hitRow>=0)
	{
		//get the match information in read lists
		pReadList = sharedReadListArr + hitRow;
		pReadPos = sharedReadPosArr + pReadList->firstRow;
		tmp_matchNum = pReadList->matchNum;
		for(j=0; j<tmp_matchNum; j++)
		{
			if(pReadPos->contigID==contigIDFrom && pReadPos->contigPos==contigPosFrom)
			{
				pReadPos->contigID = 0;
				pReadList->curNum --;
			}
		}
	}else
	{ // errors
		printf("line=%d, In %s(), cannot get the read [ %ld ] from SRL, error!\n", __LINE__, __func__, readID);
		return FAILED;
	}

	return SUCCESSFUL;
}


/**
 * Update the reads information in read lists, including RL1, RL2 and SRL.
 *  Notes:
 *  	Only update the contigPos of the reads.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateReadsInReadLists(uint64_t readID, int contigIDFrom, int contigPosFrom, int contigLenDecrease)
{
	int j, hitRow, tmp_matchNum;
	ReadList *pReadList;
	ReadPos *pReadPos;

	if(readID%2==1)
	{ // odd number
		hitRow = getReadRowFromReadList(readID, readListArr[1], readNumInRL[1]);
		if(hitRow>=0)
		{
			// get the match information in read lists
			pReadList = readListArr[1] + hitRow;
			pReadPos = readPosArr[1] + pReadList->firstRow;
			tmp_matchNum = pReadList->matchNum;
			for(j=0; j<tmp_matchNum; j++)
			{
				if(pReadPos->contigID==contigIDFrom && pReadPos->contigPos==contigPosFrom)
				{
					pReadPos->contigPos -= contigLenDecrease;
				}
			}
		}else
		{ // errors
			printf("line=%d, In %s(), cannot get the read [ %ld ] from RL1, error!\n", __LINE__, __func__, readID);
			return FAILED;
		}
	}else
	{ // even number
		hitRow = getReadRowFromReadList(readID, readListArr[2], readNumInRL[2]);
		if(hitRow>=0)
		{
			//get the match information in read lists
			pReadList = readListArr[2] + hitRow;
			pReadPos = readPosArr[2] + pReadList->firstRow;
			tmp_matchNum = pReadList->matchNum;
			for(j=0; j<tmp_matchNum; j++)
			{
				if(pReadPos->contigID==contigIDFrom && pReadPos->contigPos==contigPosFrom)
				{
					pReadPos->contigPos -= contigLenDecrease;
				}
			}
		}else
		{ // errors
			printf("line=%d, In %s(), cannot get the read [ %ld ] from RL2, error!\n", __LINE__, __func__, readID);
			return FAILED;
		}
	}

	hitRow = getReadRowFromReadList(readID, sharedReadListArr, readItemNumInSRL);
	if(hitRow>=0)
	{
		//get the match information in read lists
		pReadList = sharedReadListArr + hitRow;
		pReadPos = sharedReadPosArr + pReadList->firstRow;
		tmp_matchNum = pReadList->matchNum;
		for(j=0; j<tmp_matchNum; j++)
		{
			if(pReadPos->contigID==contigIDFrom && pReadPos->contigPos==contigPosFrom)
			{
				pReadPos->contigPos -= contigLenDecrease;
			}
		}
	}else
	{ // errors
		printf("line=%d, In %s(), cannot get the read [ %ld ], from SRL, error!\n", __LINE__, __func__, readID);
		return FAILED;
	}

	return SUCCESSFUL;
}


/**
 * Output contig overlap information to text file.
 *  Format:
 *  	(1) Header fields: >scaffoldID, linkedNum, rowsNum, which are separated by tab character;
 *  	(2)   Body fields: contigID1, orientation1, contigLen1, contigID2, orientation2, contigLen2, mergeFlag, overlapLen, gapSize, breakFlag, which are separated by tab character.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputContigOverlapInfoToFile(const char *contigOverlapInfoFile)
{
	int i, j;
	contigOverlap *pContigOverlapInfo;
	FILE *fpOverlapInfo;
	int scaffoldID, linkedContigNum, rowsNum, startRow;
	int subRowsNum, subStartRow, subEndRow, subLinkedContigNum;

	fpOverlapInfo = fopen(contigOverlapInfoFile, "w");
	if(fpOverlapInfo==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigOverlapInfoFile);
		return FAILED;
	}

	// break the erroneous links and output the overlap information
	scaffoldID = 0;
	for(i=0; i<scaffoldsNumInCOI; i++)
	{
		// generate sub scaffolds
		startRow = contigOverlapIndexArr[i].startRow;
		rowsNum = contigOverlapIndexArr[i].rowsNum;
		linkedContigNum = contigOverlapIndexArr[i].linkedNum;
		pContigOverlapInfo = contigOverlapInfoArr + startRow;

		if(linkedContigNum<=1)
		{ // only one contig, then do nothing
			if(linkedContigNum<1)
			{ // error linkedContigNum
				printf("line=%d, In %s(), linkedContigNum=%d, error!\n", __LINE__, __func__, linkedContigNum);
				return FAILED;
			}else
			{
				fprintf(fpOverlapInfo, ">%d\t%d\t%d\n", ++scaffoldID, 1, 1);
				fprintf(fpOverlapInfo, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
						pContigOverlapInfo->contigID1, ORIENTATION_PLUS, contigInfoArr[pContigOverlapInfo->contigID1-1].contigLen,
						0, 0, 0, 0, 0, 0, 0);

				continue;
			}
		}

		// generate new linked information, and output it to file
		for(j=0; j<rowsNum; j++)
		{
			for(subEndRow=j; subEndRow<rowsNum; subEndRow++)
			{
				if(pContigOverlapInfo[subEndRow].breakFlag==YES)
					break;
			}

			if(subEndRow==rowsNum)
			{
				subRowsNum = subEndRow - j;
				subLinkedContigNum = subRowsNum + 1;

				fprintf(fpOverlapInfo, ">%d\t%d\t%d\n", ++scaffoldID, subLinkedContigNum, subRowsNum);
				for(subStartRow=j; subStartRow<subEndRow; subStartRow++)
				{
					fprintf(fpOverlapInfo, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
							pContigOverlapInfo[subStartRow].contigID1, pContigOverlapInfo[subStartRow].orientation1, contigInfoArr[pContigOverlapInfo[subStartRow].contigID1-1].contigLen,
							pContigOverlapInfo[subStartRow].contigID2, pContigOverlapInfo[subStartRow].orientation2, contigInfoArr[pContigOverlapInfo[subStartRow].contigID2-1].contigLen,
							pContigOverlapInfo[subStartRow].mergeFlag, pContigOverlapInfo[subStartRow].overlapLen, pContigOverlapInfo[subStartRow].gapSize, pContigOverlapInfo[subStartRow].breakFlag);
				}

			}else if(subEndRow<rowsNum)
			{ // A break happens
				subRowsNum = subEndRow - j;

				if(subRowsNum>0)
				{ // There are some valid rows in the sub scaffold
					subLinkedContigNum = subRowsNum + 1;

					fprintf(fpOverlapInfo, ">%d\t%d\t%d\n", ++scaffoldID, subLinkedContigNum, subRowsNum);
					for(subStartRow=j; subStartRow<subEndRow; subStartRow++)
					{
						fprintf(fpOverlapInfo, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
								pContigOverlapInfo[subStartRow].contigID1, pContigOverlapInfo[subStartRow].orientation1, contigInfoArr[pContigOverlapInfo[subStartRow].contigID1-1].contigLen,
								pContigOverlapInfo[subStartRow].contigID2, pContigOverlapInfo[subStartRow].orientation2, contigInfoArr[pContigOverlapInfo[subStartRow].contigID2-1].contigLen,
								pContigOverlapInfo[subStartRow].mergeFlag, pContigOverlapInfo[subStartRow].overlapLen, pContigOverlapInfo[subStartRow].gapSize, pContigOverlapInfo[subStartRow].breakFlag);
					}

					//=========================== Modification begin ==================================
					if(subEndRow==rowsNum-1)
					{ // break the last link
						fprintf(fpOverlapInfo, ">%d\t%d\t%d\n", ++scaffoldID, 1, 1);
						fprintf(fpOverlapInfo, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
								pContigOverlapInfo[subEndRow].contigID2, ORIENTATION_PLUS, contigInfoArr[pContigOverlapInfo[subEndRow].contigID2-1].contigLen,
								0, 0, 0, 0, 0, 0, 0);
					}
					//=========================== Modification end ==================================

				}else // subRowsNum==0
				{ // No valid rows in the sub scaffold
					if(rowsNum==1)
					{ // only two contigs in the scaffold, then break it into two separate contigs
						fprintf(fpOverlapInfo, ">%d\t%d\t%d\n", ++scaffoldID, 1, 1);
						fprintf(fpOverlapInfo, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
								pContigOverlapInfo[subEndRow].contigID1, ORIENTATION_PLUS, contigInfoArr[pContigOverlapInfo[subEndRow].contigID1-1].contigLen,
								0, 0, 0, 0, 0, 0, 0);

						fprintf(fpOverlapInfo, ">%d\t%d\t%d\n", ++scaffoldID, 1, 1);
						fprintf(fpOverlapInfo, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
								pContigOverlapInfo[subEndRow].contigID2, ORIENTATION_PLUS, contigInfoArr[pContigOverlapInfo[subEndRow].contigID2-1].contigLen,
								0, 0, 0, 0, 0, 0, 0);
					}else
					{ // more than two contigs in the scaffold
						if(subEndRow==rowsNum-1)
						{ // break the last link
							if(pContigOverlapInfo[subEndRow-1].breakFlag==YES)
							{ // the previous link has been broken, then break the current link
								fprintf(fpOverlapInfo, ">%d\t%d\t%d\n", ++scaffoldID, 1, 1);
								fprintf(fpOverlapInfo, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
										pContigOverlapInfo[subEndRow].contigID1, ORIENTATION_PLUS, contigInfoArr[pContigOverlapInfo[subEndRow].contigID1-1].contigLen,
										0, 0, 0, 0, 0, 0, 0);

								fprintf(fpOverlapInfo, ">%d\t%d\t%d\n", ++scaffoldID, 1, 1);
								fprintf(fpOverlapInfo, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
										pContigOverlapInfo[subEndRow].contigID2, ORIENTATION_PLUS, contigInfoArr[pContigOverlapInfo[subEndRow].contigID2-1].contigLen,
										0, 0, 0, 0, 0, 0, 0);
							}else
							{
								fprintf(fpOverlapInfo, ">%d\t%d\t%d\n", ++scaffoldID, 1, 1);
								fprintf(fpOverlapInfo, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
										pContigOverlapInfo[subEndRow].contigID2, ORIENTATION_PLUS, contigInfoArr[pContigOverlapInfo[subEndRow].contigID2-1].contigLen,
										0, 0, 0, 0, 0, 0, 0);
							}

						}else if(subEndRow==0)
						{ // break the first link

							fprintf(fpOverlapInfo, ">%d\t%d\t%d\n", ++scaffoldID, 1, 1);
							fprintf(fpOverlapInfo, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
									pContigOverlapInfo[subEndRow].contigID1, ORIENTATION_PLUS, contigInfoArr[pContigOverlapInfo[subEndRow].contigID1-1].contigLen,
									0, 0, 0, 0, 0, 0, 0);
						}else
						{ // breaks in middle part
							if(pContigOverlapInfo[subEndRow-1].breakFlag==YES)
							{ // the previous link has been broken, then break the current link

								fprintf(fpOverlapInfo, ">%d\t%d\t%d\n", ++scaffoldID, 1, 1);
								fprintf(fpOverlapInfo, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
										pContigOverlapInfo[subEndRow].contigID1, ORIENTATION_PLUS, contigInfoArr[pContigOverlapInfo[subEndRow].contigID1-1].contigLen,
										0, 0, 0, 0, 0, 0, 0);
							}else
							{ // the previous link has not been broken, then nothing will be done for the current link

							}
						}
					}
				}
			}else
			{ // error
				printf("line=%d, In %s(), subEndRow=%d > rowsNum=%d, error!\n", __LINE__, __func__, subEndRow, rowsNum);
				return FAILED;
			}

			j = subEndRow;
		}
	}

	printf("There are %d new scaffolds.\n", scaffoldID);

	fclose(fpOverlapInfo);

	return SUCCESSFUL;
}

/**
 * Output mean size and standard deviation of fragments to text file.
 *  Format: mean size \t standard deviation.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputMeanSizeAndSDev(const char *meanSdevFile)
{
	FILE *fpMeanSdev;

	fpMeanSdev = fopen(meanSdevFile, "w");
	if(fpMeanSdev==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, meanSdevFile);
		return FAILED;
	}

	fprintf(fpMeanSdev, "%lf\t%lf\n", meanSizeInsert, stardardDeviationInsert);

	fclose(fpMeanSdev);
	fpMeanSdev = NULL;

	return SUCCESSFUL;
}

/**
 * Output contig information array to text file.
 *  Format:
 *  	(1) Header format: >contigID length: contigLen
 *  	(2)   Body format: base sequence
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputContigInfoArrToFile(const char *contigFileName, contigInfo *pContigInfoArray, int tmpContigsNum)
{
	int i;
	FILE *fpNewContig;

	fpNewContig = fopen(contigFileName, "w");
	if(fpNewContig==NULL)
	{
		printf("line=%d, In %s(), cannot open the new contig file [ %s ], error!\n", __LINE__, __func__, contigFileName);
		return FAILED;
	}

	for(i=0; i<tmpContigsNum; i++)
	{
		fprintf(fpNewContig, ">%d length: %d\n%s\n", pContigInfoArray[i].contigID, pContigInfoArray[i].contigLen, pContigInfoArray[i].contigSeq);
	}

	fclose(fpNewContig);

	return SUCCESSFUL;
}

/**
 * Rewrite the RL1, RL2, SRL and CL to binary files.
 *  @return:
 *   If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short rewriteReadListsAndContigListToFiles(const char *sharedReadListFile, const char *readListFile1, const char *readListFile2, const char *contigListFile)
{
	// rewrite RL1 to file
	if(rewriteSingleReadListToFile(readListFile1, readListArr[1], readNumInRL[1], readPosArr[1], matchItemNumInRP[1])==FAILED)
	{
		printf("line=%d, In %s(), cannot rewrite read list 1 (RL1) to file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// rewrite RL1 to file
	if(rewriteSingleReadListToFile(readListFile2, readListArr[2], readNumInRL[2], readPosArr[2], matchItemNumInRP[2])==FAILED)
	{
		printf("line=%d, In %s(), cannot rewrite read list 2 (RL2) to file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// rewrite SRL to file
	if(rewriteSingleReadListToFile(sharedReadListFile, sharedReadListArr, readItemNumInSRL, sharedReadPosArr, matchItemNumInSRP)==FAILED)
	{
		printf("line=%d, In %s(), cannot rewrite shared read list (SRL) to file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// rewrite contig list to file
	if(rewriteContigListToFile(contigListFile, contigListArr, contigItemNumInCL, contigReadArr, contigReadItemNumInCR)==FAILED)
	{
		printf("line=%d, In %s(), cannot rewrite contig list (CL) to file, error!\n", __LINE__, __func__);
		return FAILED;
	}


	return SUCCESSFUL;
}

/**
 * Rewrite single read list to binary file.
 *  @return:
 *   If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short rewriteSingleReadListToFile(const char *readListFile, ReadList *pReadListArray, int64_t readItemNum, ReadPos *pReadPosArray, int64_t matchItemNum)
{
	int64_t i;
	FILE *fpReadList;
	int64_t total_readItemNum, total_matchItemNum, tmp_totalFirstRow;

	fpReadList = fopen(readListFile, "wb");
	if(fpReadList==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readListFile);
		return FAILED;
	}

	// get the item number of read list array, and the item number of read position array
	total_readItemNum = 0;
	total_matchItemNum = 0;
	for(i=0; i<readItemNum; i++)
	{
		if(pReadListArray[i].curNum>0)
		{
			total_readItemNum ++;
			total_matchItemNum += pReadListArray[i].curNum;
		}
	}

	// output the updated read list to binary file
	// save the item numbers
	if(fwrite(&total_readItemNum, sizeof(int64_t), 1, fpReadList)!=1)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(fwrite(&total_matchItemNum, sizeof(int64_t), 1, fpReadList)!=1)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	// save the read list array
	tmp_totalFirstRow = 0;
	for(i=0; i<readItemNum; i++)
	{
		if(pReadListArray[i].curNum>0)
		{
			pReadListArray[i].matchNum = pReadListArray[i].curNum;
			pReadListArray[i].firstRow = tmp_totalFirstRow;

			tmp_totalFirstRow += pReadListArray[i].curNum;

			if(fwrite(pReadListArray+i, sizeof(ReadList), 1, fpReadList)!=1)
			{
				printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	// save the read position array
	for(i=0; i<matchItemNum; i++)
	{
		if(pReadPosArray[i].contigID>0)
		{
			if(fwrite(pReadPosArray+i, sizeof(ReadPos), 1, fpReadList)!=1)
			{
				printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	fclose(fpReadList);
	fpReadList = NULL;

	return SUCCESSFUL;
}

/**
 * Rewrite contig list to binary file.
 *  @return:
 *   If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short rewriteContigListToFile(const char *contigListFile, ContigList *pContigListArray, int64_t contigItemNum, ContigRead *pContigReadArray, int64_t contigReadItemNum)
{
	int64_t i;
	FILE *fpContigList;
	int64_t total_contigReadItemNum, tmp_totalFirstRow;

	fpContigList = fopen(contigListFile, "wb");
	if(fpContigList==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigListFile);
		return FAILED;
	}

	// get the item number of contig read array
	total_contigReadItemNum = 0;
	for(i=0; i<contigItemNum; i++)
		total_contigReadItemNum += pContigListArray[i].curNum5 + pContigListArray[i].curNum3;

	// save the contig item number in contig list array and the contig read item number in contig read array
	if(fwrite(&contigItemNum, sizeof(int64_t), 1, fpContigList)!=1)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(fwrite(&total_contigReadItemNum, sizeof(int64_t), 1, fpContigList)!=1)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<contigItemNum; i++)
	{
		pContigListArray[i].EndNum5 = pContigListArray[i].curNum5;
		pContigListArray[i].EndNum3 = pContigListArray[i].curNum3;
	}

	tmp_totalFirstRow = 0;
	for(i=0; i<contigItemNum; i++)
	{
		if(pContigListArray[i].EndNum5>0)
		{
			pContigListArray[i].firstRow5 = tmp_totalFirstRow;
			tmp_totalFirstRow += pContigListArray[i].EndNum5;
		}

		if(pContigListArray[i].EndNum3>0)
		{
			pContigListArray[i].firstRow3 = tmp_totalFirstRow;
			tmp_totalFirstRow += pContigListArray[i].EndNum3;
		}
	}

	// save the contig list array data
	if(fwrite(pContigListArray, sizeof(ContigList), contigItemNum, fpContigList)!=contigItemNum)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	// save the contig read array data
	for(i=0; i<contigReadItemNum; i++)
	{
		if(pContigReadArray[i].readID>0)
		{
			if(fwrite(pContigReadArray+i, sizeof(ContigRead), 1, fpContigList)!=1)
			{
				printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	fclose(fpContigList);
	fpContigList = NULL;

	return SUCCESSFUL;
}

/**
 * Load mean size and standard deviation of fragments from text file.
 *  Format: mean size \t standard deviation.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short loadMeanSizeAndSDev(const char *meanSdevFile, double *meanSizeFrag, double *sDevFrag)
{
	FILE *fpMeanSdev;

	fpMeanSdev = fopen(meanSdevFile, "r");
	if(fpMeanSdev==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, meanSdevFile);
		return FAILED;
	}

	fscanf(fpMeanSdev, "%lf\t%lf\n", meanSizeFrag, sDevFrag);

	fclose(fpMeanSdev);
	fpMeanSdev = NULL;

	return SUCCESSFUL;
}

/**
 * Load contig overlap information.
 *  @return:
 *   If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short loadContigOverlapInfo(const char *contigOverlapInfoFile, contigOverlapIndex **pContigOverlapIndexArray, int *scaffoldsNum, contigOverlap **pContigOverlapInfoArray, int *itemNumInContigOverlapArray)
{
	// get the total number of scaffolds
	*scaffoldsNum = getScaffoldsNum(contigOverlapInfoFile);
	if((*scaffoldsNum)<=0)
	{
		printf("line=%d, In %s(), cannot get the total number of scaffolds, error\n", __LINE__, __func__);
		return FAILED;
	}

	*itemNumInContigOverlapArray = getOverlapItemNumInOverlapInfoFile(contigOverlapInfoFile);
	if((*itemNumInContigOverlapArray)<=0)
	{
		printf("line=%d, In %s(), cannot get the total item number of contig overlap information, error\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate the memory
	*pContigOverlapIndexArray = (contigOverlapIndex*) calloc(*scaffoldsNum, sizeof(contigOverlapIndex));
	if((*pContigOverlapIndexArray)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error\n", __LINE__, __func__);
		return FAILED;
	}

	*pContigOverlapInfoArray = (contigOverlap*) calloc(*itemNumInContigOverlapArray, sizeof(contigOverlap));
	if((*pContigOverlapInfoArray)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the data of contig overlap index array and contig overlap information array
	if(fillDataFromOverlapInfoFile(contigOverlapInfoFile, *pContigOverlapIndexArray, *scaffoldsNum, *pContigOverlapInfoArray, *itemNumInContigOverlapArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the contig overlap data, error\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get total item number of contig overlap information from contig overlap information file.
 *  @return:
 *  	If succeeds, return the item number; otherwise, return ERROR.
 */
short getOverlapItemNumInOverlapInfoFile(const char *contigOverlapInfoFile)
{
	FILE *fpOverlapInfo;
	int totalItemNum, scaffoldID, linkedContigNum, rowsNum;
	char ch;

	fpOverlapInfo = fopen(contigOverlapInfoFile, "r");
	if(fpOverlapInfo==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigOverlapInfoFile);
		return ERROR;
	}

	// compute the total number
	totalItemNum = 0;
	ch = fgetc(fpOverlapInfo);
	while(ch!=EOF)
	{
		if(ch=='>')
		{
			fscanf(fpOverlapInfo, "%d\t%d\t%d\n", &scaffoldID, &linkedContigNum, &rowsNum);
			totalItemNum += rowsNum;
		}

		ch = fgetc(fpOverlapInfo);
	}

	fclose(fpOverlapInfo);

	return totalItemNum;
}

/**
 * Fill the contig overlap information data from contig overlap information file.
 *  @return:
 *   If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillDataFromOverlapInfoFile(const char *contigOverlapInfoFile, contigOverlapIndex *pContigOverlapIndexArray, int scaffoldsNum, contigOverlap *pContigOverlapInfoArray, int itemNumInContigOverlapArray)
{
	FILE *fpOverlapInfo;
	char ch;
	int i, j, scaffoldID, linkedContigsNum, rowsNum, tmp_startRow;
	int tmp_contigID1, tmp_orient1, tmp_contigLen1, tmp_contigID2, tmp_orient2, tmp_contigLen2, tmp_mergeFlag, tmp_overlapLen, tmp_gapSize, tmp_breakFlag;

	fpOverlapInfo = fopen(contigOverlapInfoFile, "r");
	if(fpOverlapInfo==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigOverlapInfoFile);
		return FAILED;
	}

	// fill the data
	tmp_startRow = 0;
	for(i=0; i<scaffoldsNum; i++)
	{
		ch = fgetc(fpOverlapInfo);
		if(ch=='>')
		{
			fscanf(fpOverlapInfo, "%d\t%d\t%d\n", &scaffoldID, &linkedContigsNum, &rowsNum);
			pContigOverlapIndexArray[scaffoldID-1].scaffoldID = scaffoldID;
			pContigOverlapIndexArray[scaffoldID-1].linkedNum = linkedContigsNum;
			pContigOverlapIndexArray[scaffoldID-1].rowsNum = rowsNum;
			pContigOverlapIndexArray[scaffoldID-1].startRow = tmp_startRow;

			for(j=0; j<rowsNum; j++)
			{
				fscanf(fpOverlapInfo, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
						&tmp_contigID1, &tmp_orient1, &tmp_contigLen1,
						&tmp_contigID2, &tmp_orient2, &tmp_contigLen2,
						&tmp_mergeFlag, &tmp_overlapLen, &tmp_gapSize, &tmp_breakFlag);

				pContigOverlapInfoArray[tmp_startRow+j].contigID1 = tmp_contigID1;
				pContigOverlapInfoArray[tmp_startRow+j].orientation1 = tmp_orient1;
				pContigOverlapInfoArray[tmp_startRow+j].contigID2 = tmp_contigID2;
				pContigOverlapInfoArray[tmp_startRow+j].orientation2 = tmp_orient2;
				pContigOverlapInfoArray[tmp_startRow+j].mergeFlag = tmp_mergeFlag;
				pContigOverlapInfoArray[tmp_startRow+j].overlapLen = tmp_overlapLen;
				pContigOverlapInfoArray[tmp_startRow+j].gapSize = tmp_gapSize;
				pContigOverlapInfoArray[tmp_startRow+j].breakFlag = tmp_breakFlag;
			}

			tmp_startRow += rowsNum;
		}else
		{
			printf("line=%d, In %s(), cannot read file [ %s ], error!\n", __LINE__, __func__, contigOverlapInfoFile);
			return FAILED;
		}
	}

#if DEBUG_FLAG
	if(tmp_startRow!=itemNumInContigOverlapArray)
	{
		printf("line=%d, In %s(), tmp_startRow=%d != itemNumInContigOverlapArray=%d, error!\n", __LINE__, __func__, tmp_startRow, itemNumInContigOverlapArray);
		return FAILED;
	}
#endif

	fclose(fpOverlapInfo);
	fpOverlapInfo = NULL;

	return SUCCESSFUL;
}

/**
 * Free the contig overlap information array.
 */
void freeContigOverlapInfo(contigOverlapIndex **pContigOverlapIndexArray, int *scaffoldsNum, contigOverlap **pContigOverlapInfoArray, int *itemNumInContigOverlapArray)
{
	*scaffoldsNum = 0;
	*itemNumInContigOverlapArray = 0;

	free(*pContigOverlapIndexArray);
	*pContigOverlapIndexArray = NULL;
	free(*pContigOverlapInfoArray);
	*pContigOverlapInfoArray = NULL;

}

/**
 * Update the used time of contig ends.
 */
short updateContigEndsInfo(contigOverlap *pContigOverlapInfo, contigInfo *pContigInfoArray)
{
	int contigID1, contigID2, contigOrient1, contigOrient2;

	contigID1 = pContigOverlapInfo->contigID1;
	contigID2 = pContigOverlapInfo->contigID2;
	contigOrient1 = pContigOverlapInfo->orientation1;
	contigOrient2 = pContigOverlapInfo->orientation2;

	if(pContigOverlapInfo->breakFlag==NO)
	{
		if(contigOrient1==ORIENTATION_PLUS)
		{ // plus orientation, then update the 3' end of the first contig
			pContigInfoArray[contigID1-1].usedTimeEnd3 ++;
			if(pContigOverlapInfo->mergeFlag==YES && pContigOverlapInfo->overlapLen>pContigInfoArray[contigID1-1].maxOverlapLenEnd3)
				pContigInfoArray[contigID1-1].maxOverlapLenEnd3 = pContigOverlapInfo->overlapLen;
		}else
		{ // minus orientation, then update the 5' end of the first contig
			pContigInfoArray[contigID1-1].usedTimeEnd5 ++;
			if(pContigOverlapInfo->mergeFlag==YES && pContigOverlapInfo->overlapLen>pContigInfoArray[contigID1-1].maxOverlapLenEnd5)
				pContigInfoArray[contigID1-1].maxOverlapLenEnd5 = pContigOverlapInfo->overlapLen;
		}

		if(contigOrient2==ORIENTATION_PLUS)
		{ // plus orientation, then update the 5' end of the second contig
			pContigInfoArray[contigID1-1].usedTimeEnd5 ++;
			if(pContigOverlapInfo->mergeFlag==YES && pContigOverlapInfo->overlapLen>pContigInfoArray[contigID1-1].maxOverlapLenEnd5)
				pContigInfoArray[contigID1-1].maxOverlapLenEnd5 = pContigOverlapInfo->overlapLen;
		}else
		{ // minus orientation, then update the 3' end of the second contig
			pContigInfoArray[contigID1-1].usedTimeEnd3 ++;
			if(pContigOverlapInfo->mergeFlag==YES && pContigOverlapInfo->overlapLen>pContigInfoArray[contigID1-1].maxOverlapLenEnd3)
				pContigInfoArray[contigID1-1].maxOverlapLenEnd3 = pContigOverlapInfo->overlapLen;
		}
	}else
	{
		printf("line=%d, In %s(), cannot update cotnig end info for breaked links of contigs [%d, %d], error!\n", __LINE__, __func__, contigID1, contigID2);
		return FAILED;
	}

	return SUCCESSFUL;
}
