/*
 * scafContigList.c
 *
 *  Created on: Jun 18, 2011
 *      Author: zhuxiao
 */

#include "inc/scafStdinc.h"
#include "inc/scafExtvab.h"


/**
 * Fill the contig list (CL) from shared read list (SRL).
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short buildContigList(const char *contigListFile, const char *sharedReadListFile, const char *contigFile)
{
	printf("=========== Begin building Contig List (CL), please wait ...\n");

	// Allocate memory and load the data of shared read list (SRL)
	if(loadSingleReadList(sharedReadListFile, &sharedReadListArr, &readItemNumInSRL, &sharedReadPosArr, &matchItemNumInSRP)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the data of shared read list (SRL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate the memory of contig list (CL)
	if(initMemContigList(contigFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the contig list (CL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// begin filling contig list (CL)
	if(fillContigList()==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the contig list (CL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// sort the contig reads to ascending order by contigPos
	if(sortContigReads()==FAILED)
	{
		printf("line=%d, In %s(), cannot sort the contig reads, error!\n", __LINE__, __func__);
		return FAILED;
	}


	//############################### Debug information ########################
#if DEBUG_FLAG
	// check the sorted contig reads
	if(checkSortedContigReads(contigListArr, contigItemNumInCL, contigReadArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot check the sorted the contig reads, error!\n", __LINE__, __func__);
		return FAILED;
	}
#endif
	//############################### Debug information ########################


	// Save contig list (CL) to file.
	// 	Format:
	// 		(1) first two int64_t integers:
	// 			item number of contig list array, and
	// 			the item number of contig read array, respectively;
	// 		(2) contig list array data;
	// 		(3) contig read array data.
	if(saveContigListToFile(contigListFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot save the contig list (CL) to file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Release the memory of shared read list
	freeSingleReadList(&sharedReadListArr, &readItemNumInSRL, &sharedReadPosArr, &matchItemNumInSRP);

	// release the memory of contig list (CL)
	freeMemContigList();

	printf("=========== End building Contig List (CL).\n");

	return SUCCESSFUL;
}

/**
 * Initialize the memory of contig list (CL).
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemContigList(const char *contigFile)
{
	// initialize variable values
	contigReadItemNumInCR = matchItemNumInSRP;
	//contigItemNumInCL = getDistinctContigIDNumFromRP(sharedReadPosArr, matchItemNumInSRP);
	contigItemNumInCL = getContigsNum(contigFile);
	if(contigItemNumInCL<=0)
	{
		printf("line=%d, In %s(), contigItemNumInCL=%ld, error!\n", __LINE__, __func__, contigItemNumInCL);
		return FAILED;
	}

	//============== allocate the memory of contig list (CL)===============
	contigListArr = (ContigList*) calloc(contigItemNumInCL, sizeof(ContigList));
	if(contigListArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	contigReadArr = (ContigRead*) calloc(contigReadItemNumInCR, sizeof(ContigRead));
	if(contigReadArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Release the memory of contig list (CL).
 */
void freeMemContigList()
{
	contigItemNumInCL = 0;
	contigReadItemNumInCR = 0;

	free(contigListArr);
	contigListArr = NULL;

	free(contigReadArr);
	contigReadArr = NULL;
}


/**
 * Fill data of contig list (CL).
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillContigList()
{
	int64_t i, j, total, tmp_itemNumInContigReadArr;
	int endFlag, tmp_firstRow, tmp_curNum, tmp_contigID;
	int64_t tmp_readID, tmp_matchNum, tmp_firstRowInRL;


	// counting the number of reads in each contig
	for(i=0; i<matchItemNumInSRP; i++)
	{
		if(sharedReadPosArr[i].contigEnd==1) // 3' end for long contigs
		{
			contigListArr[sharedReadPosArr[i].contigID-1].EndNum3 ++;
		}else
		{
			contigListArr[sharedReadPosArr[i].contigID-1].EndNum5 ++;
		}
	}

	// compute the firstRow of each contig list node
	total = 0;
	for(i=0; i<contigItemNumInCL; i++)
	{
		contigListArr[i].contigID = i + 1;

		if(contigListArr[i].EndNum5>0)
		{
			contigListArr[i].firstRow5 = total;
			total += contigListArr[i].EndNum5;
		}

		if(contigListArr[i].EndNum3>0)
		{
			contigListArr[i].firstRow3 = total;
			total += contigListArr[i].EndNum3;
		}
	}

	//############################## Debug information ########################
#if DEBUG_FLAG
	if(total!=matchItemNumInSRP)
	{
		printf("line=%d, In %s(), total=%ld is not equal to matchItemNumInSRP=%ld, error!\n", __LINE__, __func__, total, matchItemNumInSRP);
		return FAILED;
	}else if(total!=contigReadItemNumInCR)
	{
		printf("line=%d, In %s(), total=%ld is not equal to contigReadItemNumInCR=%ld, error!\n", __LINE__, __func__, total, contigReadItemNumInCR);
		return FAILED;
	}
#endif
	//############################## Debug information ########################

	// fill contig reads data
	tmp_itemNumInContigReadArr = 0;
	for(i=0; i<readItemNumInSRL; i++)
	{
		tmp_readID = sharedReadListArr[i].readID;
		tmp_matchNum = sharedReadListArr[i].matchNum;
		tmp_firstRowInRL = sharedReadListArr[i].firstRow;

		for(j=0; j<tmp_matchNum; j++)
		{
			tmp_contigID = sharedReadPosArr[tmp_firstRowInRL+j].contigID;
			endFlag = sharedReadPosArr[tmp_firstRowInRL+j].contigEnd;
			if(endFlag==1)
			{
				tmp_firstRow = contigListArr[tmp_contigID-1].firstRow3;
				tmp_curNum = contigListArr[tmp_contigID-1].curNum3;
			}else
			{
				tmp_firstRow = contigListArr[tmp_contigID-1].firstRow5;
				tmp_curNum = contigListArr[tmp_contigID-1].curNum5;
			}

			contigReadArr[tmp_firstRow + tmp_curNum].readID = tmp_readID;
			contigReadArr[tmp_firstRow + tmp_curNum].contigPos = sharedReadPosArr[tmp_firstRowInRL+j].contigPos;
			contigReadArr[tmp_firstRow + tmp_curNum].contigEnd = sharedReadPosArr[tmp_firstRowInRL+j].contigEnd;
			contigReadArr[tmp_firstRow + tmp_curNum].orientation = sharedReadPosArr[tmp_firstRowInRL+j].orientation;

			if(endFlag==1)
			{
				contigListArr[tmp_contigID-1].curNum3 ++;
			}else
			{
				contigListArr[tmp_contigID-1].curNum5 ++;
			}

			tmp_itemNumInContigReadArr ++;
		}
	}

	//############################## Debug information ########################
#if DEBUG_FLAG
	if(tmp_itemNumInContigReadArr!=contigReadItemNumInCR)
	{
		printf("line=%d, In %s(), tmp_itemNumInContigReadArr=%ld is not equal to contigReadItemNumInCR=%ld, error!\n", __LINE__, __func__, tmp_itemNumInContigReadArr, contigReadItemNumInCR);
		return FAILED;
	}
#endif
	//############################## Debug information ########################


	return SUCCESSFUL;
}

/**
 * Get the orientation of a read at end of a contig.
 *  @return:
 *  	If succeeds, return the orientation; otherwise, return ERROR.
 */
int getReadOrientContigEnd(int64_t readID, int contigID, int contigEnd, ContigList *pContigListArr, int contigListSize, ContigRead *contigReadArr)
{
	int hitReadRow;
	ContigRead *pContigReads;
	int readsNum;

	if(contigEnd==0) // 5' end
	{
		pContigReads = contigReadArr + pContigListArr[readID-1].firstRow5;
		readsNum = pContigListArr[readID-1].EndNum5;
	}else
	{
		pContigReads = contigReadArr + pContigListArr[readID-1].firstRow3;
		readsNum = pContigListArr[readID-1].EndNum3;
	}

	hitReadRow = getReadRowInContigReads(readID, pContigReads, readsNum);
	if(hitReadRow>=0)
	{
		return pContigReads[hitReadRow].orientation;
	}else
	{
		printf("line=%d, In %s(), cannot get the read (%ld) from contig read array (CR), error!\n", __LINE__, __func__, readID);
		return ERROR;
	}
}

/**
 * Get the row of the read in contig reads array.
 *  @return:
 *  	If succeeds, return the row; otherwise, return -1.
 */
int getReadRowInContigReads(int64_t readID, ContigRead *pContigReads, int arraySize)
{
	int mid, low, high;
	low = 0;
	high = arraySize - 1;
	while(low<=high)
	{
		mid = (low + high) / 2;
		if(pContigReads[mid].readID==readID)
		{
			return mid;
		}else if(pContigReads[mid].readID<readID)
		{
			low = mid + 1;
		}else
		{
			high = mid - 1;
		}
	}

	return -1;
}

/**
 * Sort the contig reads to ascending order by contigPos.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short sortContigReads()
{
	int i;
	for(i=0; i<contigItemNumInCL; i++)
	{
		//########################## Debug information #######################
		//printf("\n================ contig %d ================\n", contigListArr[i].contigID);
		//########################## Debug information #######################

		// sort the reads at the 5' end of a contig
		if(contigListArr[i].EndNum5>1)
		{
			//########################## Debug information #######################
			//printf("The reads at the 5' end of contig %d before sorting:\n", contigListArr[i].contigID);
			//outputContigReadsSingleContigEnd(contigReadArr+contigListArr[i].firstRow5, contigListArr[i].EndNum5);
			//########################## Debug information #######################

			if(sortContigReadSingleContigEnd(contigReadArr+contigListArr[i].firstRow5, contigListArr[i].EndNum5)==FAILED)
			{
				printf("line=%d, In %s(), cannot sort the reads at the 5' end of contig %d, error!\n", __LINE__, __func__, contigListArr[i].contigID);
				return FAILED;
			}

			//########################## Debug information #######################
			//printf("The reads at the 5' end of contig %d after sorting:\n", contigListArr[i].contigID);
			//outputContigReadsSingleContigEnd(contigReadArr+contigListArr[i].firstRow5, contigListArr[i].EndNum5);
			//########################## Debug information #######################
		}

		// sort the reads at the 3' end of a contig
		if(contigListArr[i].EndNum3>1)
		{
			//########################## Debug information #######################
			//printf("The reads at the 3' end of contig %d before sorting:\n", contigListArr[i].contigID);
			//outputContigReadsSingleContigEnd(contigReadArr+contigListArr[i].firstRow3, contigListArr[i].EndNum3);
			//########################## Debug information #######################

			if(sortContigReadSingleContigEnd(contigReadArr+contigListArr[i].firstRow3, contigListArr[i].EndNum3)==FAILED)
			{
				printf("line=%d, In %s(), cannot sort the reads at the 3' end of contig %d, error!\n", __LINE__, __func__, contigListArr[i].contigID);
				return FAILED;
			}

			//########################## Debug information #######################
			//printf("The reads at the 3' end of contig %d after sorting:\n", contigListArr[i].contigID);
			//outputContigReadsSingleContigEnd(contigReadArr+contigListArr[i].firstRow3, contigListArr[i].EndNum3);
			//########################## Debug information #######################
		}
	}


	return SUCCESSFUL;
}

/**
 * Sort the reads at the the end of a contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short sortContigReadSingleContigEnd(ContigRead *pContigReadArr, int readsNum)
{
	struct sortCounter
	{
		int startRow;
		int rowsNum;
		int curNum;
	};

	int i, j, minPos, maxPos, startRow, curNum, rowsNum;
	int itemNumCounter;
	int *sortRowArray;

	struct sortCounter *sortCounterArray;
	ContigRead *contigReadBuf;


	// get the maxPos and minPos
	minPos = INT_MAX;
	maxPos = 0;
	for(i=0; i<readsNum; i++)
	{
		if(pContigReadArr[i].contigPos>maxPos)
			maxPos = pContigReadArr[i].contigPos;

		if(pContigReadArr[i].contigPos<minPos)
			minPos = pContigReadArr[i].contigPos;
	}

	// allocate the memory
	itemNumCounter = maxPos - minPos + 1;
	sortCounterArray = (struct sortCounter*) calloc(itemNumCounter, sizeof(struct sortCounter));
	if(sortCounterArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	sortRowArray = (int *) calloc(readsNum, sizeof(int));
	if(sortRowArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	contigReadBuf = (ContigRead *) calloc(readsNum, sizeof(ContigRead));
	if(contigReadBuf==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}


	// count the read number of each element in counter array
	for(i=0; i<readsNum; i++)
		sortCounterArray[ pContigReadArr[i].contigPos - minPos ].rowsNum ++;

	// initialize the start row of each element in counter array
	startRow = 0;
	for(i=0; i<itemNumCounter; i++)
	{
		if(sortCounterArray[i].rowsNum>0)
		{
			sortCounterArray[i].startRow = startRow;
			startRow += sortCounterArray[i].rowsNum;
		}
	}


	// fill the data of sort row array
	for(i=0; i<readsNum; i++)
	{
		startRow = sortCounterArray[ pContigReadArr[i].contigPos - minPos ].startRow;
		curNum = sortCounterArray[ pContigReadArr[i].contigPos - minPos ].curNum;
		sortRowArray[ startRow + curNum ] = i;
		sortCounterArray[ pContigReadArr[i].contigPos - minPos ].curNum ++;
	}

	// sort the reads in the counter array
	for(i=0; i<itemNumCounter; i++)
	{
		if(sortCounterArray[i].rowsNum>1)
		{
			if(selectionSortContigReads(sortRowArray+sortCounterArray[i].startRow, sortCounterArray[i].rowsNum, pContigReadArr)==FAILED)
			{
				printf("line=%d, In %s(), cannot sort the reads, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	// copy the reads data to the right place in contig reads buffer array
	for(i=0; i<itemNumCounter; i++)
	{
		if(sortCounterArray[i].rowsNum>0)
		{
			startRow = sortCounterArray[i].startRow;
			rowsNum = sortCounterArray[i].rowsNum;
			for(j=0; j<rowsNum; j++)
			{
				if(memcpy(contigReadBuf+startRow+j, pContigReadArr+sortRowArray[startRow+j], sizeof(ContigRead))==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}

	// copy back the reads data
	if(memcpy(pContigReadArr, contigReadBuf, readsNum*sizeof(ContigRead))==NULL)
	{
		printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
		return FAILED;
	}


	free(sortCounterArray);
	free(sortRowArray);
	free(contigReadBuf);

	return SUCCESSFUL;
}

/**
 * Selection sort the reads to ascending order to their readID.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short selectionSortContigReads(int *pSortRowArray, int rowsNum, ContigRead *pContigReadArray)
{
	int i, j, minRow, tmp;
	uint64_t minReadID;

	for(i=0; i<rowsNum-1; i++)
	{
		// get the minimal contigPos
		minReadID = pContigReadArray[ pSortRowArray[i] ].readID;
		minRow = i;
		for(j=i+1; j<rowsNum; j++)
		{
			if(pContigReadArray[ pSortRowArray[j] ].readID < minReadID)
			{
				minReadID = pContigReadArray[ pSortRowArray[j] ].readID;
				minRow = j;
			}
		}

		if(i!=minRow)
		{ // exchange the item i and j
			tmp = pSortRowArray[i];
			pSortRowArray[i] = pSortRowArray[minRow];
			pSortRowArray[minRow] = tmp;
		}
	}

	return SUCCESSFUL;
}


/**
 * Save contig list (CL) to a binary file.
 * 	File format:
 * 		(1) first two int64_t integers:
 * 			item number of contig list array, and
 * 			the item number of contig read array, respectively;
 * 		(2) contig list array data;
 * 		(3) contig read array data.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short saveContigListToFile(const char *contigListFile)
{
	FILE *fpCL;
	uint64_t i, curContigReadItemNumInCR, totalItemNum;
	ContigList tmp_contigList;

	fpCL = fopen(contigListFile, "wb");
	if(fpCL==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigListFile);
		return FAILED;
	}

	curContigReadItemNumInCR = 0;
	for(i=0; i<contigItemNumInCL; i++)
		curContigReadItemNumInCR += contigListArr[i].curNum3 + contigListArr[i].curNum5;

	// save the contig item number in contig list array and the contig read item number in contig read array
	if(fwrite(&contigItemNumInCL, sizeof(int64_t), 1, fpCL)!=1)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(fwrite(&curContigReadItemNumInCR, sizeof(int64_t), 1, fpCL)!=1)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	// save the contig list array data
	totalItemNum = 0;
	for(i=0; i<contigItemNumInCL; i++)
	{
		if(memcpy(&tmp_contigList, contigListArr+i, sizeof(ContigList))==NULL)
		{
			printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		tmp_contigList.EndNum5 = tmp_contigList.curNum5;
		tmp_contigList.EndNum3 = tmp_contigList.curNum3;

		tmp_contigList.firstRow5 = totalItemNum;
		totalItemNum += tmp_contigList.curNum5;

		tmp_contigList.firstRow3 = totalItemNum;
		totalItemNum += tmp_contigList.curNum3;

		if(fwrite(&tmp_contigList, sizeof(ContigList), 1, fpCL)!=1)
		{
			printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// save the contig read array data
	for(i=0; i<contigReadItemNumInCR; i++)
	{
		if(contigReadArr[i].readID>0)
		{
			if(fwrite(contigReadArr+i, sizeof(ContigRead), 1, fpCL)!=1)
			{
				printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	fclose(fpCL);
	fpCL = NULL;

	return SUCCESSFUL;
}


/**
 * Allocate memory and load the data of contig list (CL).
 * 	File format:
 * 		(1) first two int64_t integers:
 * 			item number of contig list array, and
 * 			the item number of contig read array, respectively;
 * 		(2) contig list array data;
 * 		(3) contig read array data.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short loadContigList(const char *contigListFile, ContigList **pContigListArr, int64_t *contigItemNum, ContigRead **pContigReadArr, int64_t *contigReadItemNum)
{
	int64_t tmp[2];
	FILE *fpCL;

	fpCL = fopen(contigListFile, "rb");
	if(fpCL==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigListFile);
		return FAILED;
	}

	// get the item numbers
	if(fread(tmp, sizeof(int64_t), 2, fpCL)!=2)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}
	*contigItemNum = tmp[0];
	*contigReadItemNum = tmp[1];

	// allocate the memory
	*pContigListArr = (ContigList*) malloc((*contigItemNum) * sizeof(ContigList));
	if(*pContigListArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	*pContigReadArr = (ContigRead*) malloc((*contigReadItemNum) * sizeof(ContigRead));
	if(*pContigReadArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}


	// get the contig list array
	if(fread(*pContigListArr, sizeof(ContigList), *contigItemNum, fpCL)!=(*contigItemNum))
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the contig read array
	if(fread(*pContigReadArr, sizeof(ContigRead), *contigReadItemNum, fpCL)!=(*contigReadItemNum))
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}

	fclose(fpCL);
	fpCL = NULL;

	return SUCCESSFUL;
}

/**
 * Release the memory of contig list.
 */
void freeContigList(ContigList **pContigListArr, int64_t *contigItemNum, ContigRead **pContigReadArr, int64_t *contigReadItemNum)
{
	*contigItemNum = 0;
	*contigReadItemNum = 0;

	free(*pContigListArr);
	*pContigListArr = NULL;

	free(*pContigReadArr);
	*pContigReadArr = NULL;
}
