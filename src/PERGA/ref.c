/*
 * ref.c
 *
 *  Created on: Jan 29, 2013
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"

/**
 * Load the reference from file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short loadReferenceFromFile(char *refFileName)
{
	// initialize the memory for the refernece
	if(initMemRef(refFileName)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Initialize the memory for reference.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemRef(char *refFileName)
{
	int32_t refLenTmp;
	FILE *fpRef;
	char ch, *refSeqTmp;

	fpRef = fopen(refFileName, "r");
	if(fpRef==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, refFileName);
		return FAILED;
	}

	itemNumRefArr = 1;
	refArr = (ref_t *) calloc (1, sizeof(ref_t));
	if(refArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get ref length
	ch = fgetc(fpRef);
	while(ch!='\n' && ch!=EOF) ch = fgetc(fpRef);

	refLenTmp = 0;
	ch = fgetc(fpRef);
	while(ch!='>' && ch!=EOF)
	{
		switch(ch)
		{
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case 'a':
			case 'c':
			case 'g':
			case 't':
				refLenTmp ++;;
		}
		ch = fgetc(fpRef);
	}
	refArr[0].reflen = refLenTmp;

	refSeqTmp = (char*) calloc(refLenTmp+1, sizeof(char));
	if(refSeqTmp==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	refArr[0].refseq = refSeqTmp;
	refArr[0].circularFlag = YES;

	// fill ref sequence
	rewind(fpRef);
	ch = fgetc(fpRef);
	while(ch!='\n' && ch!=EOF) ch = fgetc(fpRef);

	refLenTmp = 0;
	ch = fgetc(fpRef);
	while(ch!='>' && ch!=EOF)
	{
		switch(ch)
		{
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case 'a':
			case 'c':
			case 'g':
			case 't':
				refSeqTmp[refLenTmp++] = ch;
		}
		ch = fgetc(fpRef);
	}
	refSeqTmp[refLenTmp] = '\0';

	fclose(fpRef);

	return SUCCESSFUL;
}

/**
 * Free the reference from file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short freeReference()
{
	int32_t i;

	for(i=0; i<itemNumRefArr; i++)
		free(refArr[i].refseq);

	free(refArr);
	refArr = NULL;

	return SUCCESSFUL;
}

/**
 * Update the ref position of contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateRefPosContig(successRead_t *successReadArray, int32_t successReadNum, ref_t *refArray)
{
	char refBase;
	int32_t newRefPosContig, refBaseInt;

	if(refPosContig>0)
	{
		switch(refStrandContig)
		{
			case STRAND_PLUS:
				newRefPosContig = refPosContig + 1;
				if(refArray[0].circularFlag==YES && newRefPosContig>refArray[0].reflen)
					newRefPosContig -= refArray[0].reflen;
				refBase = refArray[0].refseq[newRefPosContig-1];
				break;
			case STRAND_MINUS:
				newRefPosContig = refPosContig - 1;
				if(refArray[0].circularFlag==YES && newRefPosContig<1)
					newRefPosContig += refArray[0].reflen;
				refBase = refArray[0].refseq[newRefPosContig-1];
				break;
			default: printf("line=%d, In %s(), invalid contig ref strand: %d, error!\n", __LINE__, __func__, refStrandContig);
				return FAILED;
		}

		switch(refBase)
		{
			case 'A':
			case 'a':
				refBaseInt = 0;
				break;
			case 'C':
			case 'c':
				refBaseInt = 1;
				break;
			case 'G':
			case 'g':
				refBaseInt = 2;
				break;
			case 'T':
			case 't':
				refBaseInt = 3;
				break;
			default:
				printf("line=%d, In %s(), invalid base in reference %d, error!\n", __LINE__, __func__, refBase);
				return FAILED;
		}

		if(refStrandContig==STRAND_MINUS)
		{
			refBaseInt = (~refBaseInt) & 3;
		}

		refBaseIntContig = refBaseInt;

		if((kmerSeqIntAssembly[entriesPerKmer-1]&3)!=refBaseInt)
		{
			refPosMatchFlag = NO;
			refPosSoildFlag = NO;
		}else
			refPosMatchFlag = YES;
		refPosContig = newRefPosContig;
	}else
	{
		printf("line=%d, In %s(), refPosContig=%d, error!\n", __LINE__, __func__, refPosContig);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Update the ref position of contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeRefPosContig(successRead_t *successReadArray, int32_t successReadNum, ref_t *refArray)
{
	int32_t i, j, refPosContigTmp, allReadsMatchFlag, validContigMatchFlag;
	int32_t refPosReadArray[successReadNum][2], refPosContigArray[successReadNum][2], refPosContigArrayMono[successReadNum][3];
	int32_t readBlockID, rowNumInReadBlock, maxItemNumPerReadBlock; // block id starts from 0
	int32_t maxRow, maxValue, monoRowsNum, tmp, tmpItemNum, minDistance, rowMinDistance, nearestRefPosValidFlag;

	readBlock_t *readBlockArr;
	read_t *pRead;
	readseqBlock_t *readseqBlockArr;

	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;
	readseqBlockArr = deBruijnGraph->readSet->readseqBlockArr;

	// get the ref position of reads
	for(i=0; i<successReadNum; i++)
	{
		readBlockID = (successReadArray[i].rid - 1) / maxItemNumPerReadBlock;
		rowNumInReadBlock = (successReadArray[i].rid - 1) % maxItemNumPerReadBlock;
		pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;

#if(DRAW_CURVE_FLAG==YES)
		refPosReadArray[i][0] = pRead->refStrand;
		refPosReadArray[i][1] = pRead->refPos;
#endif
	}

	// generate the ref position of the contig according to the reads as well as the ref pos match kind
	for(i=0; i<successReadNum; i++)
	{
		if(refPosReadArray[i][0]==STRAND_PLUS && successReadArray[i].orientation==ORIENTATION_PLUS)
		{
			refPosContigArray[i][0] = STRAND_PLUS;
			refPosContigTmp = refPosReadArray[i][1] + successReadArray[i].seqlen - 1;
			if(refArray[0].circularFlag==YES && refPosContigTmp>refArray[0].reflen)
				refPosContigTmp -= refArray[0].reflen;
			refPosContigArray[i][1] = refPosContigTmp;
		}else if(refPosReadArray[i][0]==STRAND_MINUS && successReadArray[i].orientation==ORIENTATION_PLUS)
		{
			refPosContigArray[i][0] = STRAND_MINUS;
			refPosContigTmp = refPosReadArray[i][1] - successReadArray[i].seqlen + 1;
			if(refArray[0].circularFlag==YES && refPosContigTmp<0)
				refPosContigTmp += refArray[0].reflen;
			refPosContigArray[i][1] = refPosContigTmp;
		}else if(refPosReadArray[i][0]==STRAND_PLUS && successReadArray[i].orientation==ORIENTATION_MINUS)
		{
			refPosContigArray[i][0] = STRAND_MINUS;
			refPosContigArray[i][1] = refPosReadArray[i][1];
		}else if(refPosReadArray[i][0]==STRAND_MINUS && successReadArray[i].orientation==ORIENTATION_MINUS)
		{
			refPosContigArray[i][0] = STRAND_PLUS;
			refPosContigArray[i][1] = refPosReadArray[i][1];
		}else
		{
			printf("line=%d, In %s(), invalid situation, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// sort the array
	for(i=0; i<successReadNum-1; i++)
	{
		for(j=i+1; j<successReadNum; j++)
		{
			if(refPosContigArray[i][0]<refPosContigArray[j][0] || refPosContigArray[i][1]<refPosContigArray[j][1])
			{ // sort according to strand, refPos
				// exchange the strand
				tmp = refPosContigArray[i][0];
				refPosContigArray[i][0] = refPosContigArray[j][0];
				refPosContigArray[j][0] = tmp;
				// exchange the refPos
				tmp = refPosContigArray[i][1];
				refPosContigArray[i][1] = refPosContigArray[j][1];
				refPosContigArray[j][1] = tmp;
			}
		}
	}

	// compute the occurrence number of each element
	i = 0; monoRowsNum = 0;
	while(i<successReadNum)
	{
		tmpItemNum = 1;
		for(j=i+1; j<successReadNum; j++)
		{
			if(refPosContigArray[i][0]==refPosContigArray[j][0] && refPosContigArray[i][1]==refPosContigArray[j][1])
				tmpItemNum ++;
			else
				break;
		}

		refPosContigArrayMono[monoRowsNum][0] = refPosContigArray[i][0];
		refPosContigArrayMono[monoRowsNum][1] = refPosContigArray[i][1];
		refPosContigArrayMono[monoRowsNum][2] = tmpItemNum;
		monoRowsNum ++;
		i += tmpItemNum;
	}

	// get the nearest refPos with previous refPos
	if(monoRowsNum>1)
	{
		nearestRefPosValidFlag = NO;
		if(refPosContig>0)
		{
			minDistance = INT_MAX;
			rowMinDistance = -1;
			for(i=1; i<monoRowsNum; i++)
			{
				if(refPosContigArrayMono[i][0]==refStrandContig && refPosContigArrayMono[i][1]==refPosContig)
				{
					tmp = refPosContigArrayMono[i][1] - refPosContig;
					if(tmp<0)
						tmp = -tmp;
					if(tmp<minDistance)
					{
						minDistance = tmp;
						rowMinDistance = i;
					}
				}
			}

			if(rowMinDistance>=0 && minDistance<=MAX_REF_POS_NEAR_THRES)
				nearestRefPosValidFlag = YES;
		}

		if(nearestRefPosValidFlag==YES)
		{
			// set the refPos
			refPosMatchFlag = YES;
			refPosSoildFlag = YES;
			refStrandContig = refPosContigArrayMono[rowMinDistance][0];
			refPosContig = refPosContigArrayMono[rowMinDistance][1];
		}else
		{
			// get the maximal occurrence
			maxRow = 0;
			maxValue = refPosContigArrayMono[maxRow][2];
			for(i=1; i<monoRowsNum; i++)
			{
				if(refPosContigArrayMono[i][2]>maxValue)
				{
					maxValue = refPosContigArrayMono[i][2];
					maxRow = i;
				}
			}

			// set the refPos
			refPosMatchFlag = YES;
			refPosSoildFlag = YES;
			refStrandContig = refPosContigArrayMono[maxRow][0];
			refPosContig = refPosContigArrayMono[maxRow][1];
		}
	}else
	{
		// set the refPos
		refPosMatchFlag = YES;
		refPosSoildFlag = YES;
		refStrandContig = refPosContigArrayMono[0][0];
		refPosContig = refPosContigArrayMono[0][1];
	}

	if(refPosFirstBase[0]==NO)
	{
		refPosFirstBase[0] = YES;
		refPosFirstBase[1] = refStrandContig;
		if(refStrandContig==STRAND_PLUS)
			refPosFirstBase[2] = refPosContig - itemNumContigArr + 1;
		else
			refPosFirstBase[2] = refPosContig + itemNumContigArr - 1;
	}

/*
	allReadsMatchFlag = YES;
	for(i=0; i<successReadNum; i++)
	{
		if(successReadArray[i].matchnum!=successReadArray[i].seqlen)
			allReadsMatchFlag = NO;
	}

	if(allReadsMatchFlag==YES)
	{
		// get the ref position of reads
		for(i=0; i<successReadNum; i++)
		{
			readBlockID = (successReadArray[i].rid - 1) / maxItemNumPerReadBlock;
			rowNumInReadBlock = (successReadArray[i].rid - 1) % maxItemNumPerReadBlock;
			pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;

#if(DRAW_CURVE_FLAG==YES)
			refPosReadArray[i][0] = pRead->refStrand;
			refPosReadArray[i][1] = pRead->refPos;
#endif
		}

		// generate the ref position of the contig according to the reads as well as the ref pos match kind
		for(i=0; i<successReadNum; i++)
		{
			if(refPosReadArray[i][0]==STRAND_PLUS && successReadArray[i].orientation==ORIENTATION_PLUS)
			{
				refPosContigArray[i][0] = STRAND_PLUS;
				refPosContigTmp = refPosReadArray[i][1] + successReadArray[i].seqlen - 1;
				if(refArray[0].circularFlag==YES && refPosContigTmp>refArray[0].reflen)
					refPosContigTmp -= refArray[0].reflen;
				refPosContigArray[i][1] = refPosContigTmp;
			}else if(refPosReadArray[i][0]==STRAND_MINUS && successReadArray[i].orientation==ORIENTATION_PLUS)
			{
				refPosContigArray[i][0] = STRAND_MINUS;
				refPosContigTmp = refPosReadArray[i][1] - successReadArray[i].seqlen + 1;
				if(refArray[0].circularFlag==YES && refPosContigTmp<0)
					refPosContigTmp += refArray[0].reflen;
				refPosContigArray[i][1] = refPosContigTmp;
			}else if(refPosReadArray[i][0]==STRAND_PLUS && successReadArray[i].orientation==ORIENTATION_MINUS)
			{
				refPosContigArray[i][0] = STRAND_MINUS;
				refPosContigArray[i][1] = refPosReadArray[i][1];
			}else if(refPosReadArray[i][0]==STRAND_MINUS && successReadArray[i].orientation==ORIENTATION_MINUS)
			{
				refPosContigArray[i][0] = STRAND_PLUS;
				refPosContigArray[i][1] = refPosReadArray[i][1];
			}else
			{
				printf("line=%d, In %s(), invalid situation, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		// determine the ref position of the contig as well as its ref match kind,
		// and the contig match kind should be the same and the ref pos should also be the same
		validContigMatchFlag = YES;
		if(successReadNum>1)
		{
			for(i=1; i<successReadNum; i++)
			{
				if(refPosContigArray[0][0]!=refPosContigArray[i][0] || refPosContigArray[0][1]!=refPosContigArray[i][1])
				{
					validContigMatchFlag = NO;
					break;
				}
			}
		}

		if(validContigMatchFlag==YES)
		{
			refStrandContig = refPosContigArray[0][0];
			refPosContig = refPosContigArray[0][1];

			if(refPosInfoFirstBase[0]==-1)
			{
				refPosInfoFirstBase[0] = refStrandContig;
				if(refStrandContig==STRAND_PLUS)
					refPosInfoFirstBase[1] = refPosContig - itemNumContigArr + 1;
				else
					refPosInfoFirstBase[1] = refPosContig + itemNumContigArr - 1;
			}
		}else
		{
			refStrandContig = -1;
			refPosContig = -1;
		}
	}else
	{
		refStrandContig = -1;
		refPosContig = -1;
	}
*/

#if(DRAW_CURVE_FLAG==YES)
	fprintf(fpRefPos, "%ld\t%d\t%ld\t%d\t%d\t%d\n", localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);
#endif

	return SUCCESSFUL;
}

/**
 * Draw the points for maxOcc and secondOcc when branches.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short drawOccPoints(int32_t navigationFlag, int32_t naviSuccessFlag, ref_t *refArr)
{
	int32_t extensionContinueFlag, extensionCorrectFlag, tmp_gapSize, sameMaxIndex, extensionContinueFlagPE_Mix, extensionCorrectFlagPE_Mix;

	// check the correctness of the extension
	if(naviSuccessFlag==NAVI_SUCCESS)
	{
		extensionContinueFlag = YES;
	}else
	{
		extensionContinueFlag = NO;
	}

	if(successContigIndex<=0)
	{
		if(extensionContinueFlag==YES)
			extensionCorrectFlag = YES;
		else
			extensionCorrectFlag = NO;
	}else
	{
		//if(refPosContig>0)		// deleted 2013-03-01
		if(refPosMatchFlag==YES)	// added 2013-03-01
		{
			if(extensionContinueFlag==YES)
				extensionCorrectFlag = YES;
			else
				extensionCorrectFlag = NO;
		}else
		{
			if(extensionContinueFlag==YES)
				extensionCorrectFlag = NO;
			else
				extensionCorrectFlag = YES;
		}
	}

	if(successContigIndex>0)
	{
//		if(navigationFlag==NAVI_PE_FLAG)
//		{
			// compute the maximal gap size in contig tail region
			if(computeGapSizeInContig(&tmp_gapSize, contigArr, itemNumContigArr, assemblyRound)==FAILED)
			{
				printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
				return FAILED;
			}
//		}else
//		{
//			tmp_gapSize = itemNumContigArr - successContigIndex;
//		}
	}else
	{
		tmp_gapSize = 0;
	}

	// get the corresponding maxOcc and secondOcc
	if(navigationFlag==NAVI_PE_FLAG)
	{ // paired ends are used for assembling

		//if(maxOccPE==secondOccPE)							// deleted 2013-03-04
		if((maxOccPE==secondOccPE) || (secondOccPE>=2 && secondOccPE>=SEC_MAX_OCC_RATIO_SVM*maxOccPE))		// deleted 2013-03-04
		{
			if(extensionContinueFlag==YES)
			{
				extensionCorrectFlag = NO;
			}else
			{
				extensionCorrectFlag = YES;
			}
		}

//		sumSecondOccPE = 0;
//		for(i=0; i<4; i++) if(i!=occsNumIndexPE[0]) sumSecondOccPE += occsNumPE[i];

		//fprintf(fpOccPoint, "%d\t%d\t%.4f\t%d\t%d\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), refPosContig, extensionContinueFlag, extensionCorrectFlag);
		//fprintf(fpOccPoint, "%d\t%d\t%.4f\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, extensionCorrectFlag, extensionContinueFlag, localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);
		//fprintf(fpOccPoint, "%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag, extensionContinueFlag, localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);
		fprintf(fpOccPoint, "%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, tmp_gapSize, extensionCorrectFlag, extensionContinueFlag, localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);

		if(extensionContinueFlag==YES)
		{ // extension situation
			if(extensionCorrectFlag==YES)
			{ // extension correct
				//fprintf(fpOccExtensionCorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				//fprintf(fpOccExtensionCorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, extensionCorrectFlag);
				//fprintf(fpOccExtensionCorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				fprintf(fpOccExtensionCorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, tmp_gapSize, extensionCorrectFlag);
				fprintf(fpOccExtensionCorrect[1], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, tmp_gapSize, extensionCorrectFlag);
			}else
			{ // extension incorrect
				//fprintf(fpOccExtensionIncorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				//fprintf(fpOccExtensionIncorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, extensionCorrectFlag);
				//fprintf(fpOccExtensionIncorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				fprintf(fpOccExtensionIncorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, tmp_gapSize, extensionCorrectFlag);
				fprintf(fpOccExtensionIncorrect[1], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, tmp_gapSize, extensionCorrectFlag);
			}
		}else
		{ // stop situation
			if(extensionCorrectFlag==YES)
			{ // stop correct
				//fprintf(fpOccStopCorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				//fprintf(fpOccStopCorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, extensionCorrectFlag);
				//fprintf(fpOccStopCorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				fprintf(fpOccStopCorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, tmp_gapSize, extensionCorrectFlag);
				fprintf(fpOccStopCorrect[1], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, tmp_gapSize, extensionCorrectFlag);
			}else
			{ // stop incorrect
				//fprintf(fpOccStopIncorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				//fprintf(fpOccStopIncorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, extensionCorrectFlag);
				//fprintf(fpOccStopIncorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				fprintf(fpOccStopIncorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, tmp_gapSize, extensionCorrectFlag);
				fprintf(fpOccStopIncorrect[1], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, tmp_gapSize, extensionCorrectFlag);
			}
		}
	}else if(navigationFlag==NAVI_MIX_FLAG)
	{ // paired ends and single ends are valid for assembling
		if(secondOccPE>0)
		{
			extensionContinueFlagPE_Mix = NO;

			if(refBaseIntContig==occsNumIndexPE[0])
			{
				extensionCorrectFlagPE_Mix = NO;
			}else
			{
				extensionCorrectFlagPE_Mix = YES;
			}


			//if(maxOccPE==secondOccPE)							// deleted 2013-03-04
			if((maxOccPE==secondOccPE) || (secondOccPE>=2 && secondOccPE>=SEC_MAX_OCC_RATIO_SVM*maxOccPE))		// deleted 2013-03-04
			{
				extensionCorrectFlagPE_Mix = YES;
			}

//			sumSecondOccPE = 0;
//			for(i=0; i<4; i++) if(i!=occsNumIndexPE[0]) sumSecondOccPE += occsNumPE[i];

			//fprintf(fpOccPoint, "%d\t%d\t%.4f\t%d\t%d\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), refPosContig, extensionContinueFlag, extensionCorrectFlag);
			//fprintf(fpOccPoint, "%d\t%d\t%.4f\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, extensionCorrectFlag, extensionContinueFlag, localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);
			//fprintf(fpOccPoint, "%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag, extensionContinueFlag, localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);
			fprintf(fpOccPoint, "%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, tmp_gapSize, extensionCorrectFlagPE_Mix, extensionContinueFlagPE_Mix, localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);

			if(extensionCorrectFlagPE_Mix==YES)
			{ // stop correct
				//fprintf(fpOccStopCorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				//fprintf(fpOccStopCorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, extensionCorrectFlag);
				//fprintf(fpOccStopCorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				fprintf(fpOccStopCorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, tmp_gapSize, extensionCorrectFlagPE_Mix);
				fprintf(fpOccStopCorrect[1], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, tmp_gapSize, extensionCorrectFlagPE_Mix);
			}else
			{ // stop incorrect
				//fprintf(fpOccStopIncorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				//fprintf(fpOccStopIncorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, extensionCorrectFlag);
				//fprintf(fpOccStopIncorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				fprintf(fpOccStopIncorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, tmp_gapSize, extensionCorrectFlagPE_Mix);
				fprintf(fpOccStopIncorrect[1], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccPE, (int32_t)secondOccPE, readsNumRatio, maxOccPE/secondOccPE, tmp_gapSize, extensionCorrectFlagPE_Mix);
			}
		}

		if(secondOccSE>0)
		{
			//if(maxOccSE==secondOccSE)							// deleted 2013-03-04
			if((maxOccSE==secondOccSE) || (secondOccSE>=2 && secondOccSE>=SEC_MAX_OCC_RATIO_SVM*maxOccSE))		// deleted 2013-03-04
			{
				if(extensionContinueFlag==YES)
				{
					extensionCorrectFlag = NO;
				}else
				{
					extensionCorrectFlag = YES;
				}
			}

//			sumSecondOccSE = 0;
//			for(i=0; i<4; i++) if(i!=occsNumIndexSE[0]) sumSecondOccSE += occsNumSE[i];

			if(maxOccPE>0)
			{ // 3+2

				if(maxOccIndexPE==maxOccIndexSE)
				{
					sameMaxIndex = YES;
				}else
				{
					if(maxOccPE==secondOccPE && secondOccIndexPE==maxOccIndexSE)
						sameMaxIndex = YES;
					else
						sameMaxIndex = NO;
				}

				//fprintf(fpOccPoint, "%d\t%d\t%.4f\t%d\t%d\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), refPosContig, extensionContinueFlag, extensionCorrectFlag);
				//fprintf(fpOccPoint, "%d\t%d\t%.4f\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, extensionCorrectFlag, extensionContinueFlag, localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);
				//fprintf(fpOccPoint, "%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag, extensionContinueFlag, localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);
				fprintf(fpOccPoint, "%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag, extensionContinueFlag, localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);

				if(extensionContinueFlag==YES)
				{ // extension situation
					if(extensionCorrectFlag==YES)
					{ // extension correct
						//fprintf(fpOccExtensionCorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						//fprintf(fpOccExtensionCorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, extensionCorrectFlag);
						//fprintf(fpOccExtensionCorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						fprintf(fpOccExtensionCorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
						fprintf(fpOccExtensionCorrect[2], "%d\t%d\t%.4f\t%.4f\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, maxOccPE/secondOccPE, sameMaxIndex, extensionCorrectFlag);
					}else
					{ // extension incorrect
						//fprintf(fpOccExtensionIncorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						//fprintf(fpOccExtensionIncorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, extensionCorrectFlag);
						//fprintf(fpOccExtensionIncorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						fprintf(fpOccExtensionIncorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
						fprintf(fpOccExtensionIncorrect[2], "%d\t%d\t%.4f\t%.4f\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, maxOccPE/secondOccPE, sameMaxIndex, extensionCorrectFlag);
					}
				}else
				{ // stop situation
					if(extensionCorrectFlag==YES)
					{ // stop correct
						//fprintf(fpOccStopCorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						//fprintf(fpOccStopCorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, extensionCorrectFlag);
						//fprintf(fpOccStopCorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						fprintf(fpOccStopCorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
						fprintf(fpOccStopCorrect[2], "%d\t%d\t%.4f\t%.4f\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, maxOccPE/secondOccPE, sameMaxIndex, extensionCorrectFlag);
					}else
					{ // stop incorrect
						//fprintf(fpOccStopIncorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						//fprintf(fpOccStopIncorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, extensionCorrectFlag);
						//fprintf(fpOccStopIncorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						fprintf(fpOccStopIncorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
						fprintf(fpOccStopIncorrect[2], "%d\t%d\t%.4f\t%.4f\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, maxOccPE/secondOccPE, sameMaxIndex, extensionCorrectFlag);
					}
				}
			}else
			{ // 3
				//fprintf(fpOccPoint, "%d\t%d\t%.4f\t%d\t%d\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), refPosContig, extensionContinueFlag, extensionCorrectFlag);
				//fprintf(fpOccPoint, "%d\t%d\t%.4f\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, extensionCorrectFlag, extensionContinueFlag, localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);
				//fprintf(fpOccPoint, "%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag, extensionContinueFlag, localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);
				fprintf(fpOccPoint, "%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag, extensionContinueFlag, localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);

				if(extensionContinueFlag==YES)
				{ // extension situation
					if(extensionCorrectFlag==YES)
					{ // extension correct
						//fprintf(fpOccExtensionCorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						//fprintf(fpOccExtensionCorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, extensionCorrectFlag);
						//fprintf(fpOccExtensionCorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						fprintf(fpOccExtensionCorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
						fprintf(fpOccExtensionCorrect[3], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
					}else
					{ // extension incorrect
						//fprintf(fpOccExtensionIncorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						//fprintf(fpOccExtensionIncorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, extensionCorrectFlag);
						//fprintf(fpOccExtensionIncorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						fprintf(fpOccExtensionIncorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
						fprintf(fpOccExtensionIncorrect[3], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
					}
				}else
				{ // stop situation
					if(extensionCorrectFlag==YES)
					{ // stop correct
						//fprintf(fpOccStopCorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						//fprintf(fpOccStopCorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, extensionCorrectFlag);
						//fprintf(fpOccStopCorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						fprintf(fpOccStopCorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
						fprintf(fpOccStopCorrect[3], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
					}else
					{ // stop incorrect
						//fprintf(fpOccStopIncorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						//fprintf(fpOccStopIncorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, extensionCorrectFlag);
						//fprintf(fpOccStopIncorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
						fprintf(fpOccStopIncorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
						fprintf(fpOccStopIncorrect[3], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
					}
				}
			}
		}

	}else
	{ // single ends are valid and used for assembling

		//if(maxOccSE==secondOccSE)							// deleted 2013-03-04
		if((maxOccSE==secondOccSE) || (secondOccSE>=2 && secondOccSE>=SEC_MAX_OCC_RATIO_SVM*maxOccSE))		// deleted 2013-03-04
		{
			if(extensionContinueFlag==YES)
			{
				extensionCorrectFlag = NO;
			}else
			{
				extensionCorrectFlag = YES;
			}
		}

//		sumSecondOccSE = 0;
//		for(i=0; i<4; i++) if(i!=occsNumIndexSE[0]) sumSecondOccSE += occsNumSE[i];

		//fprintf(fpOccPoint, "%d\t%d\t%.4f\t%d\t%d\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), refPosContig, extensionContinueFlag, extensionCorrectFlag);
		//fprintf(fpOccPoint, "%d\t%d\t%.4f\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, extensionCorrectFlag, extensionContinueFlag, localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);
		//fprintf(fpOccPoint, "%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag, extensionContinueFlag, localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);
		fprintf(fpOccPoint, "%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag, extensionContinueFlag, localContigID, contigsNum+1, itemNumContigArr, assemblyRound, refStrandContig, refPosContig);

		if(extensionContinueFlag==YES)
		{ // extension situation
			if(extensionCorrectFlag==YES)
			{ // extension correct
				//fprintf(fpOccExtensionCorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				//fprintf(fpOccExtensionCorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, extensionCorrectFlag);
				//fprintf(fpOccExtensionCorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				fprintf(fpOccExtensionCorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
				fprintf(fpOccExtensionCorrect[3], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
			}else
			{ // extension incorrect
				//fprintf(fpOccExtensionIncorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				//fprintf(fpOccExtensionIncorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, extensionCorrectFlag);
				//fprintf(fpOccExtensionIncorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				fprintf(fpOccExtensionIncorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
				fprintf(fpOccExtensionIncorrect[3], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
			}
		}else
		{ // stop situation
			if(extensionCorrectFlag==YES)
			{ // stop correct
				//fprintf(fpOccStopCorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				//fprintf(fpOccStopCorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, extensionCorrectFlag);
				//fprintf(fpOccStopCorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				fprintf(fpOccStopCorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
				fprintf(fpOccStopCorrect[3], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
			}else
			{ // stop incorrect
				//fprintf(fpOccStopIncorrect, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				//fprintf(fpOccStopIncorrect, "%d\t%d\t%.4f\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, extensionCorrectFlag);
				//fprintf(fpOccStopIncorrect, "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, (int32_t)(itemNumContigArr-successContigIndex), extensionCorrectFlag);
				fprintf(fpOccStopIncorrect[0], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
				fprintf(fpOccStopIncorrect[3], "%d\t%d\t%.4f\t%.4f\t%d\t%d\n", (int32_t)maxOccSE, (int32_t)secondOccSE, readsNumRatio, maxOccSE/secondOccSE, tmp_gapSize, extensionCorrectFlag);
			}
		}
	}

	return SUCCESSFUL;
}
