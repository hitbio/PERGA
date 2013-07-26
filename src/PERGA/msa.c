/*
 * msa.c
 *
 *  Created on: Apr 11, 2013
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Generate alignment of reads for multi-align.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateMSAMalign(char **alignResultsMalign, char **readseqMalignArray, int32_t itemNumMalignArray)
{
	int32_t coreIndex;
	char tmpStr[MAX_SEQ_LEN_MALIGN+1];

	// get the core sequence index
	if(getCoreIndexMSAMalign(&coreIndex, readseqMalignArray, itemNumMalignArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the core sequence index of multiple sequence alignment, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// exchange the seq[0] and seq[core]
	if(coreIndex!=0)
	{
		strcpy(tmpStr, readseqMalignArray[0]);
		strcpy(readseqMalignArray[0], readseqMalignArray[coreIndex]);
		strcpy(readseqMalignArray[coreIndex], tmpStr);
	}

	// generate the alignment
	if(generateMSAResultMalign(alignResultsMalign, readseqMalignArray, itemNumMalignArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the multiple sequence alignment, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// exchange the seq[0] and seq[core]
	if(coreIndex!=0)
	{
		strcpy(tmpStr, readseqMalignArray[0]);
		strcpy(readseqMalignArray[0], readseqMalignArray[coreIndex]);
		strcpy(readseqMalignArray[coreIndex], tmpStr);

		strcpy(tmpStr, alignResultsMalign[0]);
		strcpy(alignResultsMalign[0], alignResultsMalign[coreIndex]);
		strcpy(alignResultsMalign[coreIndex], tmpStr);
	}

	// remove the right end spaces
	if(removeRightSpaces(alignResultsMalign, itemNumMalignArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove the right spaces of the multiple sequence alignment, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Generate alignment of reads for multi-align.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getCoreIndexMSAMalign(int32_t *coreIndex, char **readseqMalignArray, int32_t itemNumMalignArray)
{
	int32_t i, j, totalScoreArray[itemNumMalignArray][itemNumMalignArray];
	int32_t maxScore, maxIndex, sum, score;

	for(i=0; i<itemNumMalignArray; i++)
	{
		for(j=i+1; j<itemNumMalignArray; j++)
		{
			if(getGolbalScoreMalign(&score, readseqMalignArray[i], readseqMalignArray[j])==FAILED)
			{
				printf("line=%d, In %s(), cannot get the global alignment score, error!\n", __LINE__, __func__);
				return FAILED;
			}
			totalScoreArray[i][j] = score;
			totalScoreArray[j][i] = score;
		}
	}

	maxScore = -INT_MAX, maxIndex = -1;
	for(i=0; i<itemNumMalignArray; i++)
	{
		sum = 0;
		for(j=0; j<itemNumMalignArray; j++)
			if(j!=i)
				sum += totalScoreArray[i][j];

		if(maxScore<sum)
		{
			maxScore = sum;
			maxIndex = i;
		}
	}

	*coreIndex = maxIndex;

	return SUCCESSFUL;
}

/**
 * Get the global alignment score of reads for pairwise alignment.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getGolbalScoreMalign(int32_t *globalScore, char *seq1, char *seq2)
{
	int32_t seqlen1 = strlen(seq1), seqlen2 = strlen(seq2);
	int32_t scoreArray[seqlen1+1][seqlen2+1];

	if(computeScorePAMalign((int32_t**)scoreArray, seq1, seq2, seqlen1, seqlen2)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute pairwise alignment score, error!\n", __LINE__, __func__);
		return FAILED;
	}

	*globalScore = scoreArray[seqlen1][seqlen2];

	return SUCCESSFUL;
}

short computeScorePAMalign(int32_t **scoreArray, char *seq1, char *seq2, int32_t seqlen1, int32_t seqlen2)
{
	int32_t i, j, lenRow, lenCol;
	int32_t maxScore, subScore1, subScore2, subScore3;
	int32_t score1, score2, score3;

	lenRow = seqlen1 + 1;
	lenCol = seqlen2 + 1;

	for(j=0; j<lenCol; j++) {*((int32_t*)scoreArray + j) = -j;}			//初始化第0行元素值为-j
	for(i=0; i<lenRow; i++) {*((int32_t*)scoreArray + i*lenCol) = -i;}	//初始化第0列元素值为-i

	for(i=1; i<lenRow; i++)
	{
		for(j=1; j<lenCol; j++)
		{
			if(subScoreMalign(&subScore1, seq1[i-1],seq2[j-1])==FAILED)
			{
				printf("line=%d, In %s(), cannot get the score, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(j==lenCol-1)
				subScore2 = 0;
			else
			{
				if(subScoreMalign(&subScore2, seq1[i-1],'-')==FAILED)
				{
					printf("line=%d, In %s(), cannot get the score, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
			if(i==lenRow-1)
				subScore3 = 0;
			else
			{
				if(subScoreMalign(&subScore3, '-',seq2[j-1])==FAILED)
				{
					printf("line=%d, In %s(), cannot get the score, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}

			score1 = (*((int32_t*)scoreArray + (i-1)*lenCol + j-1)) + subScore1;
			score2 = (*((int32_t*)scoreArray + (i-1)*lenCol + j)) + subScore2;
			score3 = (*((int32_t*)scoreArray + i*lenCol + j-1)) + subScore3;

//			score1 = scoreArray[i-1][j-1] + subScore1;
//			score2 = scoreArray[i-1][j] + subScore2;
//			score3 = scoreArray[i][j-1] + subScore3;

			maxScore = score1;
			if(maxScore<score2)
			{
				maxScore = score2;
			}
			if(maxScore<score3)
			{
				maxScore = score3;
			}
			*((int32_t*)scoreArray + i*lenCol + j) = maxScore;

//			scoreArray[i][j] = maxScore;
		}
	}

	// output the score array
//	for(i=0; i<lenRow; i++)
//	{
//		for(j=0; j<lenCol-1; j++)
//			printf("%d\t", *((int32_t*)scoreArray + i*lenCol + j));
//		printf("%d\n", *((int32_t*)scoreArray + i*lenCol + lenCol-1));
//	}

	return SUCCESSFUL;
}

/**
 * Get the score of two characters.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short subScoreMalign(int32_t *score, char ch1, char ch2)
{
	*score = 0;

	if(ch1 == ch2) *score = matchScore;
	else if(ch1=='-'|| ch2=='-') *score = gapScore;
	else if(ch1!=ch2) *score = mismatchScore;
	else
	{
		printf("line=%d, In %s(), cannot get the sub score, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Generate the multiple alignment.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateMSAResultMalign(char **alignResultsMalign, char **readseqMalignArray, int32_t itemNumMalignArray)
{
	int32_t i, j, k, t, h, corelen, seqlen2, resultlen0, tmplen, tmplenTmp;
	char *coreSeq, *tmpResult[2], tmpStr1[MAX_SEQ_LEN_MALIGN+1], tmpStr2[MAX_SEQ_LEN_MALIGN+1];

	tmpResult[0] = tmpStr1;
	tmpResult[1] = tmpStr2;

	coreSeq = readseqMalignArray[0];
	corelen = strlen(coreSeq);
	for(k=1; k<itemNumMalignArray; k++)
	{
		seqlen2 = strlen(readseqMalignArray[k]);
		int32_t scoreArray[corelen+1][seqlen2+1];

		if(computeScorePAMalign((int32_t**)scoreArray, coreSeq, readseqMalignArray[k], corelen, seqlen2)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute pairwise alignment score, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// get the pairwise alignment
		if(computeResultPAMalign(tmpResult, coreSeq, readseqMalignArray[k], corelen, seqlen2, (int32_t**)scoreArray)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the result of pairwise alignment, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(k==1)
		{
			strcpy(alignResultsMalign[0], tmpResult[0]);
			strcpy(alignResultsMalign[1], tmpResult[1]);
		}else
		{
			resultlen0 = strlen(alignResultsMalign[0]);
			tmplen = strlen(tmpResult[0]);
			i = 0; j = 0;
			while(i<resultlen0 && j<tmplen)
			{
				if(alignResultsMalign[0][i]=='-' && tmpResult[0][j]!='-')
				{ // [i]='-', [j]!='-', then add '-' before [j] into tmpResult
					for(h=tmplen; h>=j; h--)
					{
						tmpResult[0][h+1] = tmpResult[0][h];
						tmpResult[1][h+1] = tmpResult[1][h];
					}
					tmpResult[0][j] = '-';
					tmpResult[1][j] = '-';
					tmplen ++;
				}else if(alignResultsMalign[0][i]!='-' && tmpResult[0][j]=='-')
				{ // [i]!='-', [j]='-], then add '-' before [i] into the aligned results[0...(k-1)]
					for(t=0; t<k; t++)
					{
						tmplenTmp = strlen(alignResultsMalign[t]);
						for(h=tmplenTmp; h>=i; h--)
							alignResultsMalign[t][h+1] = alignResultsMalign[t][h];
						alignResultsMalign[t][i] = '-';
					}
					resultlen0 ++;
				}
				i++; j++;
			}

			while(i<resultlen0)  //在结果序列中存在未添加的列
				tmpResult[1][i++] = '-';
			tmpResult[1][resultlen0] = '\0';
			while(j<tmplen)  //在新序列中存在未添加的列
			{
				for(t=0; t<k; t++)
					alignResultsMalign[t][j] = '-';
				j++;
			}
			for(t=0; t<k; t++)
				alignResultsMalign[t][tmplen] = '\0';

			strcpy(alignResultsMalign[k], tmpResult[1]);
		}
	}

	return SUCCESSFUL;
}

/**
 * Generate the pairwise alignment result.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeResultPAMalign(char *tmpResult[2], char *seq1, char *seq2, int32_t seqlen1, int32_t seqlen2, int32_t **scoreArray)
{
	int i, j, orient, tmplen;
	char ch;

	tmplen = 0;
	i = seqlen1;
	j = seqlen2;
	while(i>0 && j>0)
	{
		if(judgeMaxOrientationMalign(&orient, scoreArray, i, j, seqlen1+1, seqlen2+1)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the orientation of pairwise alignment, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(orient==0)
		{
			tmpResult[0][tmplen] = seq1[i-1];
			tmpResult[1][tmplen++] = seq2[j-1];
			i--; j--;
		}else if(orient==1)
		{
			tmpResult[0][tmplen] = seq1[i-1];
			tmpResult[1][tmplen++] = '-';
			i--;
		}else
		{
			tmpResult[0][tmplen] = '-';
			tmpResult[1][tmplen++] = seq2[j-1];
			j--;
		}
	}

	while(i>0)
	{
		tmpResult[0][tmplen] = seq1[i-1];
		tmpResult[1][tmplen++] = '-';
		i--;
	}
	while(j>0)
	{
		tmpResult[0][tmplen] = '-';
		tmpResult[1][tmplen++] = seq2[j-1];
		j--;
	}
	tmpResult[0][tmplen] = '\0';
	tmpResult[1][tmplen] = '\0';

	for(i=0; i<tmplen/2; i++)
	{
		ch = tmpResult[0][i];
		tmpResult[0][i] = tmpResult[0][tmplen-i-1];
		tmpResult[0][tmplen-i-1] = ch;

		ch = tmpResult[1][i];
		tmpResult[1][i] = tmpResult[1][tmplen-i-1];
		tmpResult[1][tmplen-i-1] = ch;
	}

	return SUCCESSFUL;
}

/**
 * Get the orientation of pairwise alignment result.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short judgeMaxOrientationMalign(int32_t *orient, int32_t **scoreArray, int32_t i, int32_t j, int32_t rowNum, int32_t colNum)
{
	int tmp0, tmp1, tmp2, tmp;

	tmp0 = *((int32_t*)scoreArray + (i-1)*colNum + j-1);
	tmp1 = *((int32_t*)scoreArray + (i-1)*colNum + j);
	tmp2 = *((int32_t*)scoreArray + i*colNum + j-1);
	tmp = *((int32_t*)scoreArray + i*colNum + j);

	*orient = 0;

	if(i==rowNum-1 && tmp==tmp2)
		*orient = 2;
	else if(i<rowNum-1 && tmp==tmp2+gapScore)
		*orient = 2;
	else if(j==colNum-1 && tmp==tmp1)
		*orient = 1;
	else if(j<colNum-1 && tmp==tmp1+gapScore)
		*orient = 1;
	else if(tmp==tmp0+matchScore || tmp==tmp0+mismatchScore)
		*orient = 0;
	else
	{
		printf("line=%d, In %s(), cannot get maximal orientation, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Remove the right end spaces of the multiple alignment result.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeRightSpaces(char **alignResultsMalign, int32_t itemNumMalignArray)
{
	int32_t i, j, len;

	len = strlen(alignResultsMalign[0]);
	for(i=0; i<itemNumMalignArray; i++)
	{
		for(j=len-1; j>=0; j--)
		{
			if(alignResultsMalign[i][j]!='-')
			{
				alignResultsMalign[i][j+1] = '\0';
				break;
			}
		}
	}

	return SUCCESSFUL;
}
