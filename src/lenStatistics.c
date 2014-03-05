/*
 * lenStatistics.c
 *
 *  Created on: Jan 6, 2013
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Statistics for contigs lengths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short contigsLenStatistics(contigGraph_t *contigGraph, int32_t minContigLen)
{
	int32_t i, *contigsLenArray, *contigsLenBufTmp, itemNumContigsLenArray, newContigNum;

	itemNumContigsLenArray = contigGraph->itemNumContigItemArray;

	contigsLenArray = (int32_t *) calloc (itemNumContigsLenArray, sizeof(int32_t));
	if(contigsLenArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	contigsLenBufTmp = (int32_t *) calloc (itemNumContigsLenArray, sizeof(int32_t));
	if(contigsLenBufTmp==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	newContigNum = 0;
	for(i=0; i<itemNumContigsLenArray; i++)
	{
		if(contigGraph->contigItemArray[i].contigLen>=minContigLen)
		{
			contigsLenArray[newContigNum] = contigGraph->contigItemArray[i].contigLen;
			newContigNum ++;
		}
	}

	// radix sort of contigsLenArray
	if(radixSortContigsLen(contigsLenArray, contigsLenBufTmp, newContigNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot sort contigs lengths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get lengths statistics
	if(computeLenStatistics(contigsLenArray, newContigNum, minContigLen, "contigs")==FAILED)
	{
		printf("line=%d, In %s(), cannot compute contigs length statistics, error!\n", __LINE__, __func__);
		return FAILED;
	}

	free(contigsLenArray);
	free(contigsLenBufTmp);

	return SUCCESSFUL;
}

/**
 * Radix sort of contigs lengths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short radixSortContigsLen(int32_t *contigsLenArray, int32_t *contigsLenBufTmp, int32_t itemNumContigsLenArray)
{
	struct partNode
	{
		int curItemNum;
		int totalItemNum;
		int firstRow;
	};

	int64_t i, step, total;
	int32_t *data, *buf;
	struct partNode *part;
	int partArrSize, stepBits, maxStepLen;
	uint64_t bitMask;
	uint64_t hashcode, firstRow, curItemNum;

	stepBits = 16;
	maxStepLen = 32;
	partArrSize = 1 << stepBits;
	bitMask = (1 << stepBits) - 1;

	part = (struct partNode *) malloc(partArrSize * sizeof(struct partNode));
	if(part==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Begin to sort
	step = 0;
	while(step!=maxStepLen)
	{
		// set the data and buf
		if(step==stepBits)
		{
			buf = contigsLenArray;
			data = contigsLenBufTmp;
		}else
		{
			data = contigsLenArray;
			buf = contigsLenBufTmp;
		}

		// count the number
		if(memset(part, 0, partArrSize * sizeof(struct partNode))==NULL)
		{
			printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
			free(part);
			return FAILED;
		}
		for(i=0; i<itemNumContigsLenArray; i++)
			//part[ (data[i].queryLen >> step) & bitMask ].totalItemNum ++;
			part[ bitMask - ((data[i] >> step) & bitMask) ].totalItemNum ++;

		// initialize the part array
		total = 0;
		for(i=0; i<partArrSize; i++)
		{
			part[i].firstRow = total;
			total += part[i].totalItemNum;
		}

		// copy the data to the right place
		for(i=0; i<itemNumContigsLenArray; i++)
		{
			hashcode = bitMask - ((data[i] >> step) & bitMask);
			firstRow = part[hashcode].firstRow;
			curItemNum = part[hashcode].curItemNum;
			buf[firstRow+curItemNum] = data[i];
//			if(memcpy(buf+firstRow+curItemNum, data+i, sizeof(int32_t))==NULL)
//			{
//				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
//				free(part);
//				return FAILED;
//			}
			part[hashcode].curItemNum ++;
		}

		step += stepBits;

		//######################## Debug information #######################
		for(i=0; i<partArrSize; i++)
		{
			if(part[i].curItemNum!=part[i].totalItemNum)
			{
				printf("line=%d, In %s(), in part[%ld], curItemNum=%d != totalItemNum=%d, error!\n", __LINE__, __func__, i, part[i].curItemNum, part[i].totalItemNum);
				free(part);
				return FAILED;
			}
		}
		//######################## Debug information #######################
	}

	free(part);
	part = NULL;

	return SUCCESSFUL;
}

/**
 * Compute statistics of contigs lengths.
 *  Including: # contigs, SumLen, maxSize, N50 Size, mean Size, median Size.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeLenStatistics(int32_t *contigsLenArray, int32_t itemNumContigsLenArray, int32_t minContigLen, char *itemName)
{
	int32_t i;
	int64_t N50_ArrIndex, totalQueryLen, tmpTotalLen, N50_TotalLen, medianSize;
	double meanSize;

	if(itemNumContigsLenArray==0)
		return SUCCESSFUL;

	totalQueryLen = 0;
	for(i=0; i<itemNumContigsLenArray; i++)
	{
		totalQueryLen += contigsLenArray[i];
	}

	N50_ArrIndex = -1;
	N50_TotalLen = totalQueryLen * 0.5;
	tmpTotalLen = 0;
	for(i=0; i<itemNumContigsLenArray; i++)
	{
		tmpTotalLen += contigsLenArray[i];
		if(tmpTotalLen>=N50_TotalLen)
		{
			N50_ArrIndex = i;
			break;
		}
	}

	if(N50_ArrIndex==-1)
	{
		printf("line=%d, In %s(), cannot compute N50 size of contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(itemNumContigsLenArray>0)
	{
		meanSize = (double)totalQueryLen / itemNumContigsLenArray;
		medianSize = contigsLenArray[itemNumContigsLenArray/2];
	}else
	{
		meanSize = 0;
		medianSize = 0;
	}

	// print the statistics
	printf("There are %d %s with length >= %d bp, and their total length %ld, \n"
			"maximal Size is %d, N50 size %d, mean Size %.2f, median size %ld.\n",
			itemNumContigsLenArray, itemName, minContigLen, totalQueryLen, contigsLenArray[0],
			contigsLenArray[N50_ArrIndex], meanSize, medianSize);

	return SUCCESSFUL;
}
