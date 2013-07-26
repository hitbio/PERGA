/*
 * align.c
 *
 *  Created on: Dec 15, 2012
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Adjust the match information of a read.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short adjustMatchInfoRead(assemblingreadtype *this_assemblingRead, contigtype *contigArray, int32_t contigNodesNum)
{
	// get read bases according to read orientation and basePos
	if(getReadBasesAlign(readSeqAlign, &readSeqLenAlign, this_assemblingRead)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the read bases for alignment, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get contig bases according to read orientation and basePos
	if(getContigBasesAlign(contigSeqAlign, &contigSeqLenAlign, this_assemblingRead, contigArray, contigNodesNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the contig bases for alignment, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// alignment to adjust the match information, including basePos, lastMatchBasePos, matchBaseNum, unmatchBaseNum, firstpos
	if(adjustMatchInfo(this_assemblingRead, readSeqAlign, readSeqLenAlign, contigSeqAlign, contigSeqLenAlign)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust the bases match information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}


/**
 * Get the align read bases. For plus orientation read, get the bases, otherwise, get the reverse complementary bases.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadBasesAlign(char *readSeqAlign, int32_t *readSeqLenAlign, assemblingreadtype *this_assemblingRead)
{
	int32_t i, baseInt, basePos, seqPos, entriesNum, baseNumLastEntry, entryRow, entryPos;
	uint64_t *readseqInt;

	readseqInt = this_assemblingRead->readseq;
	entriesNum = this_assemblingRead->entriesNumReadseq;
	baseNumLastEntry = this_assemblingRead->baseNumLastEentryReadseq;
	basePos = this_assemblingRead->basePos;

	if(this_assemblingRead->orientation==ORIENTATION_PLUS)
	{ // plus orientation
		// get the read bases from positions [0 ... basePos]
		entryRow = entryPos = 0;
		for(i=0; i<=basePos; i++)
		{
			if(entryRow==entriesNum-1)
				baseInt = (readseqInt[entryRow] >> (2*(baseNumLastEntry-1)-2*entryPos)) & 3;
			else
				baseInt = (readseqInt[entryRow] >> (62-2*entryPos)) & 3;

			entryPos ++;
			if(entryPos==32)
			{
				entryRow ++;
				entryPos = 0;
			}

			switch(baseInt)
			{
				case 0: readSeqAlign[i] = 'A'; break;
				case 1: readSeqAlign[i] = 'C'; break;
				case 2: readSeqAlign[i] = 'G'; break;
				case 3: readSeqAlign[i] = 'T'; break;
			}
		}
		readSeqAlign[i] = '\0';

		*readSeqLenAlign = i;

	}else
	{ // minus orientation
		// get the read bases from positions [seqLen-1 ... basePos]
		seqPos = 0;
		entryRow = (this_assemblingRead->seqlen - 1) >> 5;
		entryPos = (this_assemblingRead->seqlen - 1) % 32;
		for(i=this_assemblingRead->seqlen-1; i>=basePos; i--, seqPos++)
		{
			if(entryRow==entriesNum-1)
				baseInt = (readseqInt[entryRow] >> (2*(baseNumLastEntry-1)-2*entryPos)) & 3;
			else
				baseInt = (readseqInt[entryRow] >> (62-2*entryPos)) & 3;

			entryPos --;
			if(entryPos<0)
			{
				entryRow --;
				entryPos = 31;
			}

			switch(baseInt)
			{
				case 0: readSeqAlign[seqPos] = 'T'; break;
				case 1: readSeqAlign[seqPos] = 'G'; break;
				case 2: readSeqAlign[seqPos] = 'C'; break;
				case 3: readSeqAlign[seqPos] = 'A'; break;
			}
		}

		readSeqAlign[seqPos] = '\0';

		*readSeqLenAlign = seqPos;
	}

	return SUCCESSFUL;
}

/**
 * Get the align contig bases.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getContigBasesAlign(char *contigSeqAlign, int32_t *contigSeqLenAlign, assemblingreadtype *this_assemblingRead, contigtype *contigArray, int32_t contigNodesNum)
{
	int32_t i, startPos, endPos, seqPos;

	// get the contig bases from positions [startPos ... endPos]
	if(this_assemblingRead->orientation==ORIENTATION_PLUS)
		startPos = contigNodesNum - 1 - this_assemblingRead->basePos;
	else
		startPos = contigNodesNum - (this_assemblingRead->seqlen - this_assemblingRead->basePos);
	endPos = contigNodesNum - 1;

	seqPos = 0;
	for(i=startPos; i<=endPos; i++, seqPos++)
	{
		switch(contigArray[i].base)
		{
			case 0: contigSeqAlign[seqPos] = 'A'; break;
			case 1: contigSeqAlign[seqPos] = 'C'; break;
			case 2: contigSeqAlign[seqPos] = 'G'; break;
			case 3: contigSeqAlign[seqPos] = 'T'; break;
		}
	}
	contigSeqAlign[seqPos] = '\0';
	*contigSeqLenAlign = seqPos;

	return SUCCESSFUL;
}

/**
 * Adjust the match information.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short adjustMatchInfo(assemblingreadtype *this_assemblingRead, char *readSeqAlign, int32_t readSeqLenAlign, char *contigSeqAlign, int32_t contigSeqLenAlign)
{
	// generate the alignment between read and contig
	if(generateAlignment(alignResultSeqArr, matchFlagArr, &alignResultSeqLen, readSeqAlign, readSeqLenAlign, contigSeqAlign, contigSeqLenAlign)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate the alignment between read bases and contig bases, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// adjust read match information by alignment
	if(adjustMatchInfoByAlign(this_assemblingRead, alignResultSeqArr, matchFlagArr, &alignResultSeqLen, readSeqAlign, &readSeqLenAlign, contigSeqAlign, contigSeqLenAlign)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust match information between read bases and contig bases, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// adjust erroneous bases information
	if(adjustErrBaseInfoByAlign(this_assemblingRead, alignResultSeqArr, matchFlagArr, alignResultSeqLen, readSeqAlign, readSeqLenAlign, contigSeqAlign, contigSeqLenAlign)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust match information between read bases and contig bases, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Generate alignment between read and contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateAlignment(char *alignResultSeqArr[3], char *matchFlagArr, int32_t *alignResultSeqLen, char *seq1, int32_t seqLen1, char *seq2, int32_t seqLen2)
{
	int32_t i, j, scoreIJ, itemNumInAlignArray, rowsNum, colsNum, maxValue;
	char tmp;


	// reset value of each item in the score array to zero
	if(memset(alignScoreArr, 0L, (maxSeqLenAlign+1)*(maxSeqLenAlign+1)*sizeof(int32_t))==NULL)
	{
		printf("line=%d, In %s(), cannot reset the memory to zero, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<3; i++)
	{
		if(memset(alignResultSeqArr[i], 0L, 2*(maxSeqLenAlign+1)*sizeof(char))==NULL)
		{
			printf("line=%d, In %s(), cannot reset the memory to zero, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	rowsNum = seqLen1 + 1;
	colsNum = seqLen2 + 1;

	// initial values
	for(i=0; i<colsNum; i++)
		alignScoreArr[i] = i * gapScore;
	for(i=1; i<rowsNum; i++)
		alignScoreArr[colsNum*i] = i * gapScore;

	// compute the scores of each element
	for(i=1; i<rowsNum; i++)
	{
		for(j=1; j<colsNum; j++)
		{
			if(seq1[i-1]==seq2[j-1])
			//if(seq1[i]==seq2[j])
				scoreIJ = matchScore;
			else
				scoreIJ = mismatchScore;

			maxValue = INT_MIN;
			// compute the maximal score
			if(alignScoreArr[(i-1)*colsNum+j-1]+scoreIJ>maxValue)
			{ // from (i-1, j-1)
				maxValue = alignScoreArr[(i-1)*colsNum+j-1]+scoreIJ;
			}
			if(alignScoreArr[(i-1)*colsNum+j]+gapScore>maxValue)
			{ // from (i-1, j)
				maxValue = alignScoreArr[(i-1)*colsNum+j]+gapScore;
			}
			if(alignScoreArr[i*colsNum+j-1]+gapScore>maxValue)
			{ // from (i, j-1)
				maxValue = alignScoreArr[i*colsNum+j-1]+gapScore;
			}

			alignScoreArr[i*colsNum+j] = maxValue;
		}
	}


	itemNumInAlignArray = 0;
	//*mismatchNum = 0;
	i = rowsNum - 1;
	j = colsNum - 1;
	while(i>0 && j>0)
	{
		if(seq1[i-1]==seq2[j-1])
		//if(seq1[i]==seq2[j])
			scoreIJ = matchScore;
		else
			scoreIJ = mismatchScore;

		if(alignScoreArr[(i-1)*colsNum+j-1]+scoreIJ==alignScoreArr[i*colsNum+j])
		{ // from (i-1, j-1)
			//if(seq1[i-1]!=seq2[j-1])
			if(scoreIJ==mismatchScore)
			{
				matchFlagArr[itemNumInAlignArray] = '1';

				//(*mismatchNum) ++;
				//printf("0:(%d,%d)\n", i, j);
#if (DEBUG_ALIGNMENT==YES)
				alignResultSeqArr[0][itemNumInAlignArray] = seq1[i-1];
				alignResultSeqArr[1][itemNumInAlignArray] = ' ';
				alignResultSeqArr[2][itemNumInAlignArray] = seq2[j-1];
#endif

				itemNumInAlignArray ++;
			}else
			{
				matchFlagArr[itemNumInAlignArray] = '0';

#if (DEBUG_ALIGNMENT==YES)
				alignResultSeqArr[0][itemNumInAlignArray] = seq1[i-1];
				alignResultSeqArr[1][itemNumInAlignArray] = '|';
				alignResultSeqArr[2][itemNumInAlignArray] = seq2[j-1];
#endif

				itemNumInAlignArray ++;
			}

			i --;
			j --;
		}else if(alignScoreArr[(i-1)*colsNum+j]+gapScore==alignScoreArr[i*colsNum+j])
		{ // from (i-1, j), gap in seq2
			matchFlagArr[itemNumInAlignArray] = '3';

			//(*mismatchNum) ++;
			//printf("1:(%d,%d)\n", i, j);

#if (DEBUG_ALIGNMENT==YES)
			alignResultSeqArr[0][itemNumInAlignArray] = seq1[i-1];
			alignResultSeqArr[1][itemNumInAlignArray] = ' ';
			alignResultSeqArr[2][itemNumInAlignArray] = '-';
#endif

			itemNumInAlignArray ++;

			i --;
		}else
		{ // from (i, j-1), gap in seq1
			matchFlagArr[itemNumInAlignArray] = '2';

			//(*mismatchNum) ++;
			//printf("2:(%d,%d)\n", i, j);

#if (DEBUG_ALIGNMENT==YES)
			alignResultSeqArr[0][itemNumInAlignArray] = '-';
			alignResultSeqArr[1][itemNumInAlignArray] = ' ';
			alignResultSeqArr[2][itemNumInAlignArray] = seq2[j-1];
#endif

			itemNumInAlignArray ++;

			j --;
		}
	}

	*alignResultSeqLen = itemNumInAlignArray;

	// recover the alignment result
#if (DEBUG_ALIGNMENT==YES)
	for(i=0; i<3; i++)
	{
		for(j=0; j<itemNumInAlignArray/2; j++)
		{
			tmp = alignResultSeqArr[i][j];
			alignResultSeqArr[i][j] = alignResultSeqArr[i][itemNumInAlignArray-1-j];
			alignResultSeqArr[i][itemNumInAlignArray-1-j] = tmp;
		}
		alignResultSeqArr[i][itemNumInAlignArray] = '\0';
	}
#endif

	for(j=0; j<itemNumInAlignArray/2; j++)
	{
		tmp = matchFlagArr[j];
		matchFlagArr[j] = matchFlagArr[itemNumInAlignArray-1-j];
		matchFlagArr[itemNumInAlignArray-1-j] = tmp;
	}
	matchFlagArr[itemNumInAlignArray] = '\0';


//#if DEBUG_OUTPUT
//	// print the alignment result
//	for(i=0;i<3; i++)
//	{
//		printf("%s\n", alignResultSeqArr[i]);
//	}
//	printf("%s\n", matchFlagArr);
//#endif


	return SUCCESSFUL;
}

/**
 * Adjust the match information by alignment between read and contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short adjustMatchInfoByAlign(assemblingreadtype *this_assemblingRead, char *alignResultSeqArr[3], char *matchFlagArr, int32_t *alignResultSeqLen, char *readSeqAlign, int32_t *readSeqLenAlign, char *contigSeqAlign, int32_t contigSeqLenAlign)
{
	int32_t i, unmatchNum, rightShiftNumRead, oldReadseqLen, newContigPos, alignResultSeqLenTmp, lastMatchedBasePosTmp;

	unmatchNum = 0;
	for(i=0; i<*alignResultSeqLen; i++)
		if(matchFlagArr[i]!='0')
			unmatchNum ++;

	if(unmatchNum<=maxUnmatchBaseNumAfterAlign)
	{
		this_assemblingRead->alignSuccessTimes ++;
		this_assemblingRead->newErrNumAfterCorrection = 0;

		rightShiftNumRead = 0;
		for(i=(*alignResultSeqLen)-1; i>=0; i--)
		{
			if(matchFlagArr[i]=='0' || matchFlagArr[i]=='1')
			{
				break;
			}else if(matchFlagArr[i]=='2')
			{  // deletion in read, insertion in contig
				rightShiftNumRead --;
			}else
			{ // deletion in contig, insertion in read
				rightShiftNumRead ++;
			}
		}

		if(this_assemblingRead->orientation==ORIENTATION_PLUS)
		{ // plus orientation
			if(rightShiftNumRead>=0)
			{ // delete shifted bases in read
				this_assemblingRead->basePos -= rightShiftNumRead;

				lastMatchedBasePosTmp = this_assemblingRead->basePos;
				for(i=(*alignResultSeqLen)-1-rightShiftNumRead; i>=0; i--)
				{
					if(matchFlagArr[i]!='2')
					{ // not deletion (gap) in contig
						if(matchFlagArr[i]=='0')
							break;

						lastMatchedBasePosTmp --;
					}
				}
				this_assemblingRead->lastMatchedBasePos = lastMatchedBasePosTmp;

				this_assemblingRead->matchBaseNum = (*alignResultSeqLen) - unmatchNum;
				this_assemblingRead->unmatchBaseNum = unmatchNum - rightShiftNumRead;

				// update the read sequence
				*readSeqLenAlign -= rightShiftNumRead;
				readSeqAlign[*readSeqLenAlign] = '\0';

				// update the alignment results
				*alignResultSeqLen -= rightShiftNumRead;
				matchFlagArr[*alignResultSeqLen] = '\0';

#if (DEBUG_ALIGNMENT==YES)
				for(i=0; i<3; i++)
					alignResultSeqArr[i][*alignResultSeqLen] = '\0';
#endif

			}else
			{  // add new bases to read
				oldReadseqLen = *readSeqLenAlign;
				newContigPos = contigSeqLenAlign + rightShiftNumRead;

				lastMatchedBasePosTmp = this_assemblingRead->basePos;
				for(i=(*alignResultSeqLen)-1+rightShiftNumRead; i>=0; i--)
				{
					if(matchFlagArr[i]!='2')
					{
						if(matchFlagArr[i]=='0')
							break;

						lastMatchedBasePosTmp --;
					}
				}
				this_assemblingRead->lastMatchedBasePos = lastMatchedBasePosTmp;

				this_assemblingRead->matchBaseNum = (*alignResultSeqLen) - unmatchNum;
				this_assemblingRead->unmatchBaseNum = unmatchNum + rightShiftNumRead;

				// get new shifted bases
				if(addNewShiftedReadBases(readSeqAlign, readSeqLenAlign, -rightShiftNumRead, this_assemblingRead, matchFlagArr, *alignResultSeqLen)==FAILED)
				{
					printf("line=%d, In %s(), cannot add new bases in reads, error!\n", __LINE__, __func__);
					return FAILED;
				}

				alignResultSeqLenTmp = (*alignResultSeqLen) + rightShiftNumRead;
				for(i=oldReadseqLen; i<*readSeqLenAlign; i++)
				{
					this_assemblingRead->basePos ++;
					if(readSeqAlign[i]==contigSeqAlign[newContigPos])
					{ // match
						this_assemblingRead->lastMatchedBasePos = this_assemblingRead->basePos;
						this_assemblingRead->matchBaseNum ++;

						// update the alignment results
						matchFlagArr[alignResultSeqLenTmp] = '0';

#if (DEBUG_ALIGNMENT==YES)
						alignResultSeqArr[1][alignResultSeqLenTmp] = '|';
						alignResultSeqArr[2][alignResultSeqLenTmp] = readSeqAlign[i];
#endif

						alignResultSeqLenTmp ++;

					}else
					{ // unmatch
						this_assemblingRead->unmatchBaseNum ++;

						// update the alignment results
						matchFlagArr[alignResultSeqLenTmp] = '1';

#if (DEBUG_ALIGNMENT==YES)
						alignResultSeqArr[1][alignResultSeqLenTmp] = ' ';
						alignResultSeqArr[2][alignResultSeqLenTmp] = readSeqAlign[i];
#endif

						alignResultSeqLenTmp ++;
					}
				}

				// ########################## Debug information ##########################
				if(alignResultSeqLenTmp!=*alignResultSeqLen)
				{
					printf("line=%d, In %s(), alignResultSeqLenTmp=%d != alignResultSeqLen=%d, error!\n", __LINE__, __func__, alignResultSeqLenTmp, *alignResultSeqLen);
					return FAILED;
				}
				// ########################## Debug information ##########################
			}
		}else
		{ // minus orientation
			if(rightShiftNumRead>=0)
			{ // delete shifted bases in read

				this_assemblingRead->basePos += rightShiftNumRead;

				lastMatchedBasePosTmp = this_assemblingRead->basePos;
				for(i=(*alignResultSeqLen)-1-rightShiftNumRead; i>=0; i--)
				{
					if(matchFlagArr[i]!='2')
					{
						if(matchFlagArr[i]=='0')
							break;

						lastMatchedBasePosTmp ++;
					}
				}
				this_assemblingRead->lastMatchedBasePos = lastMatchedBasePosTmp;

				this_assemblingRead->matchBaseNum = (*alignResultSeqLen) - unmatchNum;
				this_assemblingRead->unmatchBaseNum = unmatchNum - rightShiftNumRead;

				// update the read sequence
				*readSeqLenAlign -= rightShiftNumRead;
				readSeqAlign[*readSeqLenAlign] = '\0';

				// update the alignment results
				*alignResultSeqLen -= rightShiftNumRead;
				matchFlagArr[*alignResultSeqLen] = '\0';

#if (DEBUG_ALIGNMENT==YES)
				for(i=0; i<3; i++)
					alignResultSeqArr[i][*alignResultSeqLen] = '\0';
#endif

			}else
			{  // add new bases to read
				oldReadseqLen = *readSeqLenAlign;
				newContigPos = contigSeqLenAlign + rightShiftNumRead;

				lastMatchedBasePosTmp = this_assemblingRead->basePos;
				for(i=(*alignResultSeqLen)-1+rightShiftNumRead; i>=0; i--)
				{
					if(matchFlagArr[i]!='2')
					{
						if(matchFlagArr[i]=='0')
							break;

						lastMatchedBasePosTmp ++;
					}
				}
				this_assemblingRead->lastMatchedBasePos = lastMatchedBasePosTmp;

				// get new shifted bases
				if(addNewShiftedReadBases(readSeqAlign, readSeqLenAlign, -rightShiftNumRead, this_assemblingRead, matchFlagArr, *alignResultSeqLen)==FAILED)
				{
					printf("line=%d, In %s(), cannot add new bases in reads, error!\n", __LINE__, __func__);
					return FAILED;
				}

				alignResultSeqLenTmp = (*alignResultSeqLen) + rightShiftNumRead;
				for(i=oldReadseqLen; i<*readSeqLenAlign; i++)
				{
					this_assemblingRead->basePos --;
					if(readSeqAlign[i]==contigSeqAlign[newContigPos])
					{ // match
						this_assemblingRead->lastMatchedBasePos = this_assemblingRead->basePos;
						this_assemblingRead->matchBaseNum ++;

						// update the alignment results
						matchFlagArr[alignResultSeqLenTmp] = '0';

#if (DEBUG_ALIGNMENT==YES)
						alignResultSeqArr[1][alignResultSeqLenTmp] = '|';
						alignResultSeqArr[2][alignResultSeqLenTmp] = readSeqAlign[i];
#endif

						alignResultSeqLenTmp ++;
					}else
					{ // unmatch
						this_assemblingRead->unmatchBaseNum ++;

						// update the alignment results
						matchFlagArr[alignResultSeqLenTmp] = '1';

#if (DEBUG_ALIGNMENT==YES)
						alignResultSeqArr[1][alignResultSeqLenTmp] = ' ';
						alignResultSeqArr[2][alignResultSeqLenTmp] = readSeqAlign[i];
#endif

						alignResultSeqLenTmp ++;
					}
				}

				// ########################## Debug information ##########################
				if(alignResultSeqLenTmp!=*alignResultSeqLen)
				{
					printf("line=%d, In %s(), alignResultSeqLenTmp=%d != alignResultSeqLen=%d, error!\n", __LINE__, __func__, alignResultSeqLenTmp, *alignResultSeqLen);
					return FAILED;
				}
				// ########################## Debug information ##########################
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Add the new shifted read bases after alignment between read and contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewShiftedReadBases(char *readSeqAlign, int32_t *readSeqLenAlign, int32_t shiftedBaseNum, assemblingreadtype *this_assemblingRead, char *matchFlagArr, int32_t alignResultSeqLen)
{
	int32_t i, baseInt, basePos, seqPos, entriesNum, baseNumLastEntry, entryRow, entryPos;
	uint64_t *readseqInt;

	readseqInt = this_assemblingRead->readseq;
	entriesNum = this_assemblingRead->entriesNumReadseq;
	baseNumLastEntry = this_assemblingRead->baseNumLastEentryReadseq;

	if(this_assemblingRead->orientation==ORIENTATION_PLUS)
	{ // plus orientation
		// get the read bases from positions [basePos+1, ..., basePos+shiftedBaseNum]

		basePos = this_assemblingRead->basePos + 1;
		entryRow = basePos >> 5;
		entryPos = basePos % 32;

		for(i=basePos; i<basePos+shiftedBaseNum; i++)
		{
			if(entryRow==entriesNum-1)
				baseInt = (readseqInt[entryRow] >> (2*(baseNumLastEntry-1)-2*entryPos)) & 3;
			else
				baseInt = (readseqInt[entryRow] >> (62-2*entryPos)) & 3;

			entryPos ++;
			if(entryPos==32)
			{
				entryRow ++;
				entryPos = 0;
			}

			switch(baseInt)
			{
				case 0: readSeqAlign[i] = 'A'; break;
				case 1: readSeqAlign[i] = 'C'; break;
				case 2: readSeqAlign[i] = 'G'; break;
				case 3: readSeqAlign[i] = 'T'; break;
			}
		}
		readSeqAlign[i] = '\0';

		*readSeqLenAlign = i;

	}else
	{ // minus orientation
		// get the read bases from positions [basePos-1 ... basePos-shiftedBaseNum]
		basePos = this_assemblingRead->basePos - 1;
		seqPos = *readSeqLenAlign;
		entryRow = basePos >> 5;
		entryPos = basePos % 32;
		for(i=basePos; i>basePos-shiftedBaseNum; i--, seqPos++)
		{
			if(entryRow==entriesNum-1)
				baseInt = (readseqInt[entryRow] >> (2*(baseNumLastEntry-1)-2*entryPos)) & 3;
			else
				baseInt = (readseqInt[entryRow] >> (62-2*entryPos)) & 3;

			entryPos --;
			if(entryPos<0)
			{
				entryRow --;
				entryPos = 31;
			}

			switch(baseInt)
			{
				case 0: readSeqAlign[seqPos] = 'T'; break;
				case 1: readSeqAlign[seqPos] = 'G'; break;
				case 2: readSeqAlign[seqPos] = 'C'; break;
				case 3: readSeqAlign[seqPos] = 'A'; break;
			}
		}

		readSeqAlign[seqPos] = '\0';

		*readSeqLenAlign = seqPos;
	}

	return SUCCESSFUL;
}


/**
 * Adjust the erroneous bases information of reads by alignment between read and contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short adjustErrBaseInfoByAlign(assemblingreadtype *this_assemblingRead, char *alignResultSeqArr[3], char *matchFlagArr, int32_t alignResultSeqLen, char *readSeqAlign, int32_t readSeqLenAlign, char *contigSeqAlign, int32_t contigSeqLenAlign)
{
	int32_t i, readPosTmp, contigPosTmp, baseInt;

	if(this_assemblingRead->alignSuccessTimes==this_assemblingRead->alignSuccessTimesOld)
		return SUCCESSFUL;

	if(this_assemblingRead->orientation==ORIENTATION_PLUS)
	{ // plus orientation
		readPosTmp = 0;
		contigPosTmp = 0;
		this_assemblingRead->errBaseNum = 0;
		for(i=0; i<alignResultSeqLen; i++)
		{
			if(matchFlagArr[i]=='0')
			{
				readPosTmp ++;
				contigPosTmp ++;
			}else if(matchFlagArr[i]=='1')
			{ // mismatch

				//if(this_assemblingRead->errBaseNum<maxErrBaseNumInCorrection)
				//{
					switch(contigSeqAlign[contigPosTmp])
					{
						case 'A': baseInt = 0; break;
						case 'C': baseInt = 1; break;
						case 'G': baseInt = 2; break;
						case 'T': baseInt = 3; break;
						default: printf("line=%d, In %s(), error base %c\n", __LINE__, __func__, contigSeqAlign[contigPosTmp]); return FAILED;
					}

					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].basePos = readPosTmp;
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].newBase = baseInt;
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].errKind = 1;
				//}
				this_assemblingRead->errBaseNum ++;

#if (DEBUG_CONTIG_CHECK==YES)
				// ########################## Debug information ######################
				if(this_assemblingRead->errBaseNum>2*maxErrBaseNumInCorrection)
				{
					printf("line=%d, In %s(), errBaseNum=%d > %d, error!\n", __LINE__, __func__, this_assemblingRead->errBaseNum, 2*maxErrBaseNumInCorrection);
					return FAILED;
				}
				// ########################## Debug information ######################
#endif

				readPosTmp ++;
				contigPosTmp ++;

			}else if(matchFlagArr[i]=='2')
			{ // deletion in read, gap in read

				switch(contigSeqAlign[contigPosTmp])
				{
					case 'A': baseInt = 0; break;
					case 'C': baseInt = 1; break;
					case 'G': baseInt = 2; break;
					case 'T': baseInt = 3; break;
					default: printf("line=%d, In %s(), error base %c\n", __LINE__, __func__, contigSeqAlign[contigPosTmp]); return FAILED;
				}

				//if(this_assemblingRead->errBaseNum<maxErrBaseNumInCorrection)
				//{
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].basePos = readPosTmp;
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].newBase = baseInt;
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].errKind = 2;
				//}
				this_assemblingRead->errBaseNum ++;

#if (DEBUG_CONTIG_CHECK==YES)
				// ########################## Debug information ######################
				if(this_assemblingRead->errBaseNum>2*maxErrBaseNumInCorrection)
				{
					printf("line=%d, In %s(), errBaseNum=%d > %d, error!\n", __LINE__, __func__, this_assemblingRead->errBaseNum, 2*maxErrBaseNumInCorrection);
					return FAILED;
				}
				// ########################## Debug information ######################
#endif

				contigPosTmp ++;

			}else if(matchFlagArr[i]=='3')
			{ // insertion in read, gap in contig

				//if(this_assemblingRead->errBaseNum<maxErrBaseNumInCorrection)
				//{
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].basePos = readPosTmp;
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].newBase = '-';
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].errKind = 3;
				//}
				this_assemblingRead->errBaseNum ++;

#if (DEBUG_CONTIG_CHECK==YES)
				// ########################## Debug information ######################
				if(this_assemblingRead->errBaseNum>2*maxErrBaseNumInCorrection)
				{
					printf("line=%d, In %s(), errBaseNum=%d > %d, error!\n", __LINE__, __func__, this_assemblingRead->errBaseNum, 2*maxErrBaseNumInCorrection);
					return FAILED;
				}
				// ########################## Debug information ######################
#endif

				readPosTmp ++;
			}else
			{
				printf("line=%d, In %s(), error matchFlag %c\n", __LINE__, __func__, matchFlagArr[i]);
				return FAILED;
			}
		}
	}else
	{ // minus orientation
		readPosTmp = this_assemblingRead->seqlen - 1;
		contigPosTmp = 0;
		this_assemblingRead->errBaseNum = 0;
		for(i=0; i<alignResultSeqLen; i++)
		{
			if(matchFlagArr[i]=='0')
			{
				readPosTmp --;
				contigPosTmp ++;
			}else if(matchFlagArr[i]=='1')
			{ // mismatch

				//if(this_assemblingRead->errBaseNum<maxErrBaseNumInCorrection)
				//{
					switch(contigSeqAlign[contigPosTmp])
					{
						case 'A': baseInt = 3; break;
						case 'C': baseInt = 2; break;
						case 'G': baseInt = 1; break;
						case 'T': baseInt = 0; break;
						default: printf("line=%d, In %s(), error base %c\n", __LINE__, __func__, contigSeqAlign[contigPosTmp]); return FAILED;
					}

					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].basePos = readPosTmp;
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].newBase = baseInt;
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].errKind = 1;
				//}
				this_assemblingRead->errBaseNum ++;

#if (DEBUG_CONTIG_CHECK==YES)
				// ########################## Debug information ######################
				if(this_assemblingRead->errBaseNum>2*maxErrBaseNumInCorrection)
				{
					printf("line=%d, In %s(), errBaseNum=%d > %d, error!\n", __LINE__, __func__, this_assemblingRead->errBaseNum, 2*maxErrBaseNumInCorrection);
					return FAILED;
				}
				// ########################## Debug information ######################
#endif

				readPosTmp --;
				contigPosTmp ++;

			}else if(matchFlagArr[i]=='2')
			{ // deletion in read, gap in read

				switch(contigSeqAlign[contigPosTmp])
				{
					case 'A': baseInt = 3; break;
					case 'C': baseInt = 2; break;
					case 'G': baseInt = 1; break;
					case 'T': baseInt = 0; break;
					default: printf("line=%d, In %s(), error base %c\n", __LINE__, __func__, contigSeqAlign[contigPosTmp]); return FAILED;
				}

				//if(this_assemblingRead->errBaseNum<maxErrBaseNumInCorrection)
				//{
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].basePos = readPosTmp;
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].newBase = baseInt;
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].errKind = 2;
				//}
				this_assemblingRead->errBaseNum ++;

#if (DEBUG_CONTIG_CHECK==YES)
				// ########################## Debug information ######################
				if(this_assemblingRead->errBaseNum>2*maxErrBaseNumInCorrection)
				{
					printf("line=%d, In %s(), errBaseNum=%d > %d, error!\n", __LINE__, __func__, this_assemblingRead->errBaseNum, 2*maxErrBaseNumInCorrection);
					return FAILED;
				}
				// ########################## Debug information ######################
#endif

				contigPosTmp ++;

			}else if(matchFlagArr[i]=='3')
			{ // insertion in read, gap in contig

				//if(this_assemblingRead->errBaseNum<maxErrBaseNumInCorrection)
				//{
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].basePos = readPosTmp;
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].newBase = '-';
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].errKind = 3;
				//}
				this_assemblingRead->errBaseNum ++;

#if (DEBUG_CONTIG_CHECK==YES)
				// ########################## Debug information ######################
				if(this_assemblingRead->errBaseNum>2*maxErrBaseNumInCorrection)
				{
					printf("line=%d, In %s(), errBaseNum=%d > %d, error!\n", __LINE__, __func__, this_assemblingRead->errBaseNum, 2*maxErrBaseNumInCorrection);
					return FAILED;
				}
				// ########################## Debug information ######################
#endif

				readPosTmp --;
			}
		}
	}

	return SUCCESSFUL;
}
