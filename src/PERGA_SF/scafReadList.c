/*
 * scafReadList.c
 *
 *  Created on: Jun 16, 2011
 *      Author: zhuxiao
 */


#include "inc/scafStdinc.h"
#include "inc/scafExtvab.h"


/**
 * Fill the Read Lists (RL1, RL2, RL3).
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short buildReadLists(const char *readListFile1, const char *readListFile2, const char *sharedReadListFile, const char *matchResultFile1, const char *matchResultFile2)
{
	printf("=========== Begin building Read Lists, please wait ...\n");

	if(initMemReadList(matchResultFile1, matchResultFile2)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the first read list
	if(fillSingleReadList(matchResultFile1, 1)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill read list for file [ %s ], error!\n", __LINE__, __func__, matchResultFile1);
		return FAILED;
	}

	// fill the second read list
	if(fillSingleReadList(matchResultFile2, 2)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill read list for file [ %s ], error!\n", __LINE__, __func__, matchResultFile2);
		return FAILED;
	}

	//################################ Debug information ######################
#if DEBUG_FLAG
	printf("readNum0=%ld, readNum1=%ld, readNum2=%ld\n", readNumInRL[0], readNumInRL[1], readNumInRL[2]);
	printf("matchItemNum0=%ld, matchItemNum1=%ld, matchItemNum2=%ld\n", matchItemNumInRP[0], matchItemNumInRP[1], matchItemNumInRP[2]);
#endif
	//################################ Debug information ######################

	// fill the third read list
	if(fillSharedReadList()==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the shared read list, error!\n", __LINE__, __func__);
		return FAILED;
	}


	//Save shared read list to a binary file.
	// 	Format:
	// 		(1) first two int64_t integers:
	// 			item number of read list array, and
	// 			the item number of read position array, respectively;
	// 		(2) read list array data;
	// 		(3) read position array data.
	if(saveSharedReadListToFile(sharedReadListFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot save the shared read list to file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//checkSharedReadList();

	if(saveReadListsToFile(readListFile1, readListFile2)==FAILED)
	{
		printf("line=%d, In %s(), cannot save the read list to file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// release the memory of read lists (RLs)
	freeMemReadList();

	printf("=========== End building Read Lists.\n");

	return SUCCESSFUL;
}

/**
 * Fill single read list given the corresponding matched result file.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillSingleReadList(const char *matchResultFile, int arrIndex)
{
	FILE *fpResult;
	readMatchInfoTemp *tmp_readMatchInfo, *pMatchInfo_i, *pMatchInfo_j;
	int64_t i, j, readNum, matchItemNum, cur_matchItemNum;
	int tmp_rowsNum, sameItemNum;
	ReadList *this_readListArr;
	ReadPos *this_readPosArr;

	fpResult = fopen(matchResultFile, "r");
	if(fpResult==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, matchResultFile);
		return FAILED;
	}

	matchItemNum = matchItemNumInRP[arrIndex];
	this_readListArr = readListArr[arrIndex];
	this_readPosArr = readPosArr[arrIndex];

	// allocate the memory
	tmp_readMatchInfo = (readMatchInfoTemp*) calloc (matchItemNum, sizeof(readMatchInfoTemp));
	if(tmp_readMatchInfo==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// load the data
	if(fread(tmp_readMatchInfo, sizeof(readMatchInfoTemp), matchItemNum, fpResult)!=matchItemNum)
	{
		if(feof(fpResult))
		{
			printf("line=%d, In %s(), matchItemNumInRP[%d]=%ld, error!\n", __LINE__, __func__, arrIndex, matchItemNum);
			return FAILED;
		}else
		{
			printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// record the match information to its corresponding read list
	i = 0;
	readNum = 0;
	cur_matchItemNum = 0;
	while(i<matchItemNum)
	{
		pMatchInfo_i = tmp_readMatchInfo + i;
		sameItemNum = 1;

		for(j=i+1; j<matchItemNum; j++)
		{
			pMatchInfo_j = tmp_readMatchInfo + j;
			if(pMatchInfo_i->readID != pMatchInfo_j->readID)
			{
				break;
			}

			sameItemNum ++;
		}

		// save the data
		tmp_rowsNum = sameItemNum;
		this_readListArr[readNum].readID = pMatchInfo_i->readID;
		this_readListArr[readNum].firstRow = cur_matchItemNum;
		this_readListArr[readNum].matchNum = tmp_rowsNum;
		this_readListArr[readNum].curNum = tmp_rowsNum;
		readNum ++;

		for(j=0; j<tmp_rowsNum; j++)
		{
			this_readPosArr[cur_matchItemNum+j].contigID = pMatchInfo_i[j].contigID;
			this_readPosArr[cur_matchItemNum+j].contigPos = pMatchInfo_i[j].contigPos;
			this_readPosArr[cur_matchItemNum+j].contigEnd = pMatchInfo_i[j].contigEnd;
			this_readPosArr[cur_matchItemNum+j].orientation = pMatchInfo_i[j].orientation;
		}
		cur_matchItemNum += tmp_rowsNum;

		i += tmp_rowsNum;
	}

	// save the read number
	readNumInRL[arrIndex] = readNum;

	//######################## Debug information #####################
	//printf("this_readListArr[%ld]: readID=%lu, firstRow=%u, matchNum=%u\n",
	//		readNum-1, this_readListArr[readNum-1].readID, this_readListArr[readNum-1].firstRow, this_readListArr[readNum-1].matchNum);
	//######################## Debug information #####################

	// free the memory
	free(tmp_readMatchInfo);
	tmp_readMatchInfo = NULL;

	// close the file
	fclose(fpResult);
	fpResult = NULL;

	return SUCCESSFUL;
}

/**
 * Initialize the memory of read lists.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemReadList(const char *matchResultFile1, const char *matchResultFile2)
{
	int i;
	for(i=0; i<3; i++)
	{
		readNumInRL[i] = 0;
	}

	matchItemNumInRP[1] = getMatchItemNum(matchResultFile1);
	if(matchItemNumInRP[1]<=0)
	{
		printf("line=%d, In %s(), invalid item number of file [ %s ], error!\n", __LINE__, __func__, matchResultFile1);
		return FAILED;
	}
	matchItemNumInRP[2] = getMatchItemNum(matchResultFile2);
	if(matchItemNumInRP[2]<=0)
	{
		printf("line=%d, In %s(), invalid item number of file [ %s ], error!\n", __LINE__, __func__, matchResultFile2);
		return FAILED;
	}

	matchItemNumInRP[0] = matchItemNumInRP[1] + matchItemNumInRP[2];

	// begin to allocate the memory
	for(i=0; i<3; i++)
	{
		readListArr[i] = (ReadList*) calloc(matchItemNumInRP[i], sizeof(ReadList));
		if(readListArr[i]==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		readPosArr[i] = (ReadPos*) calloc(matchItemNumInRP[i], sizeof(ReadPos));
		if(readPosArr[i]==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Free memory of read list.
 */
void freeMemReadList()
{
	int i;
	for(i=0; i<3; i++)
	{
		free(readListArr[i]);
		readListArr[i] = NULL;

		free(readPosArr[i]);
		readPosArr[i] = NULL;
	}
}

/**
 * Get the item number of a matched result file.
 *  @return:
 *  	If succeed, return the maximal item number; otherwise, return ERROR.
 */
int64_t getMatchItemNum(const char *matchResultFile)
{
	int64_t fileBytes;
	fileBytes = getFileSize(matchResultFile);
	if(fileBytes==ERROR)
	{
		printf("line=%d, In %s(), cannot get the size of file [ %s ], error!\n", __LINE__, __func__, matchResultFile);
		return ERROR;
	}else if(fileBytes==0)
	{
		printf("line=%d, In %s(), there are no matched reads in file [ %s ]!\n", __LINE__, __func__, matchResultFile);
		return ERROR;
	}

	return fileBytes / sizeof(readMatchInfoTemp);
}


/**
 * Get the file size.
 *  @return:
 *  	If succeed, return the bytes of the file; otherwise, return ERROR.
 */
int64_t getFileSize(const char *filename)
{
	struct stat buf;
	int intStat;
	intStat = stat(filename, &buf);
	if(intStat==-1)
	{
		printf("line=%d, In %s(), Cannot get the size of file [ %s ], please make sure that the file is exist.\n", __LINE__, __func__, filename);
		return ERROR;
	}
	return (int64_t)buf.st_size;
}


/**
 * Fill the shared read list (or the third read list).
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillSharedReadList()
{
	int64_t i, j, readItemNum1, readItemNum2;
	ReadList *this_readListArr1, *this_readListArr2, *shared_readListArr;
	ReadPos *shared_readPosArr, *this_readPosArr1, *this_readPosArr2;
	int64_t readID, paired_readID, hitRow;
	int64_t sharedReadNum, sharedMatchItemNum;

	// get the global variable values
	readItemNum1 = readNumInRL[1];
	readItemNum2 = readNumInRL[2];
	this_readListArr1 = readListArr[1];
	this_readListArr2 = readListArr[2];
	shared_readListArr = readListArr[0];
	this_readPosArr1 = readPosArr[1];
	this_readPosArr2 = readPosArr[2];
	shared_readPosArr = readPosArr[0];

	sharedReadNum = 0;
	sharedMatchItemNum = 0;

	for(i=0; i<readItemNum1; i++)
	{
		if(this_readListArr1[i].curNum>0)
		{
			readID = this_readListArr1[i].readID;
			if(readID%2==1) // odd number -> even number
			{
				paired_readID = readID + 1;
			}else	//even number -> odd number
			{
				paired_readID = readID - 1;
			}

			hitRow = getReadRowFromReadList(paired_readID, this_readListArr2, readItemNum2);
			if(hitRow>=0)
			{ // find the item
				if(this_readListArr2[hitRow].curNum>0)
				{
					// check the two readIDs
					if(readID<paired_readID)
					{
						// save the match information of the first paired reads
						if(memcpy(shared_readListArr+sharedReadNum, this_readListArr1+i, sizeof(ReadList))==NULL)
						{
							printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
						shared_readListArr[sharedReadNum].firstRow = sharedMatchItemNum;
						shared_readListArr[sharedReadNum].matchNum = this_readListArr1[i].curNum;

						for(j=0; j<this_readListArr1[i].matchNum; j++)
						{
							if(this_readPosArr1[this_readListArr1[i].firstRow+j].contigID>0)
							{
								if(memcpy(shared_readPosArr+sharedMatchItemNum+j, this_readPosArr1+this_readListArr1[i].firstRow+j, sizeof(ReadPos))==NULL)
								{
									printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
									return FAILED;
								}
							}
						}
						sharedReadNum ++;
						sharedMatchItemNum += this_readListArr1[i].curNum;

						// save the match information of the second paired reads
						if(memcpy(shared_readListArr+sharedReadNum, this_readListArr2+hitRow, sizeof(ReadList))==NULL)
						{
							printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
						shared_readListArr[sharedReadNum].firstRow = sharedMatchItemNum;
						shared_readListArr[sharedReadNum].matchNum = this_readListArr2[hitRow].curNum;

						for(j=0; j<this_readListArr2[hitRow].matchNum; j++)
						{
							if(this_readPosArr2[this_readListArr2[hitRow].firstRow+j].contigID>0)
							{
								if(memcpy(shared_readPosArr+sharedMatchItemNum+j, this_readPosArr2+this_readListArr2[hitRow].firstRow+j, sizeof(ReadPos))==NULL)
								{
									printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
									return FAILED;
								}
							}
						}
						sharedReadNum ++;
						sharedMatchItemNum += this_readListArr2[hitRow].curNum;

					}else
					{
						// save the match information of the second paired reads
						if(memcpy(shared_readListArr+sharedReadNum, this_readListArr2+hitRow, sizeof(ReadList))==NULL)
						{
							printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
						shared_readListArr[sharedReadNum].firstRow = sharedMatchItemNum;
						shared_readListArr[sharedReadNum].matchNum = this_readListArr2[hitRow].curNum;

						for(j=0; j<this_readListArr2[hitRow].matchNum; j++)
						{
							if(this_readPosArr2[this_readListArr2[hitRow].firstRow+j].contigID>0)
							{
								if(memcpy(shared_readPosArr+sharedMatchItemNum+j, this_readPosArr2+this_readListArr2[hitRow].firstRow+j, sizeof(ReadPos))==NULL)
								{
									printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
									return FAILED;
								}
							}
						}
						sharedReadNum ++;
						sharedMatchItemNum += this_readListArr2[hitRow].curNum;

						// save the match information of the first paired reads
						if(memcpy(shared_readListArr+sharedReadNum, this_readListArr1+i, sizeof(ReadList))==NULL)
						{
							printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
						shared_readListArr[sharedReadNum].firstRow = sharedMatchItemNum;
						shared_readListArr[sharedReadNum].matchNum = this_readListArr1[i].curNum;

						for(j=0; j<this_readListArr1[i].matchNum; j++)
						{
							if(this_readPosArr1[this_readListArr1[i].firstRow+j].contigID>0)
							{
								if(memcpy(shared_readPosArr+sharedMatchItemNum+j, this_readPosArr1+this_readListArr1[i].firstRow+j, sizeof(ReadPos))==NULL)
								{
									printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
									return FAILED;
								}
							}
						}
						sharedReadNum ++;
						sharedMatchItemNum += this_readListArr1[i].curNum;

					}
				}
			}
		}
	}
	readNumInRL[0] = sharedReadNum;
	matchItemNumInRP[0] = sharedMatchItemNum;

	//########################## Debug information ##########################
#if DEBUG_OUT_FLAG
	printf("sharedReadNum=%ld, sharedMatchItemNum=%ld\n", sharedReadNum, sharedMatchItemNum);
#endif
	//########################## Debug information ##########################

	return SUCCESSFUL;
}

/**
 * Fill the mono read list.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillMonoReadList()
{
	int i, j, k, readItemNum1, readItemNum2;
	ReadList *this_readListArr1, *this_readListArr2, *mono_readListArr;
	ReadPos *mono_readPosArr, *this_readPosArr1, *this_readPosArr2;
	int64_t monoReadNum, monoMatchItemNum;

	// get the global variable values
	readItemNum1 = readNumInRL[1];
	readItemNum2 = readNumInRL[2];
	this_readListArr1 = readListArr[1];
	this_readListArr2 = readListArr[2];
	mono_readListArr = readListArr[0];
	this_readPosArr1 = readPosArr[1];
	this_readPosArr2 = readPosArr[2];
	mono_readPosArr = readPosArr[0];

	monoReadNum = 0;
	monoMatchItemNum = 0;

	i = 0;
	j = 0;
	while(i<readItemNum1 && j< readItemNum2)
	{
		if(this_readListArr1[i].curNum<=0) { i++; continue; }
		if(this_readListArr2[j].curNum<=0) { j++; continue; }

		if(this_readListArr1[i].readID+1 < this_readListArr2[j].readID)
		{

			if(memcpy(mono_readListArr+monoReadNum, this_readListArr1+i, sizeof(ReadList))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			mono_readListArr[monoReadNum].firstRow = monoMatchItemNum;
			mono_readListArr[monoReadNum].matchNum = this_readListArr1[i].curNum;

			for(k=0; k<this_readListArr1[i].matchNum; k++)
			{
				if(this_readPosArr1[this_readListArr1[i].firstRow+k].contigID>0)
				{
					if(memcpy(mono_readPosArr+monoMatchItemNum+k, this_readPosArr1+this_readListArr1[i].firstRow+k, sizeof(ReadPos))==NULL)
					{
						printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
			}
			monoReadNum ++;
			monoMatchItemNum += this_readListArr1[i].curNum;
			i++;

		}else if(this_readListArr1[i].readID+1 > this_readListArr2[j].readID)
		{

			if(memcpy(mono_readListArr+monoReadNum, this_readListArr2+j, sizeof(ReadList))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			mono_readListArr[monoReadNum].firstRow = monoMatchItemNum;
			mono_readListArr[monoReadNum].matchNum = mono_readListArr[monoReadNum].curNum;

			for(k=0; k<this_readListArr2[j].matchNum; k++)
			if(memcpy(mono_readPosArr+monoMatchItemNum+k, this_readPosArr2+this_readListArr2[j].firstRow+k, sizeof(ReadPos))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			monoReadNum ++;
			monoMatchItemNum += this_readListArr2[j].curNum;
			j++;

		}else
		{
			i++;
			j++;
		}
	}

	// process the rest reads
	if(i<readItemNum1)
	{
		for(; i<readItemNum1; i++)
		{
			if(this_readListArr1[i].curNum>0)
			{
				if(memcpy(mono_readListArr+monoReadNum, this_readListArr1+i, sizeof(ReadList))==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				mono_readListArr[monoReadNum].firstRow = monoMatchItemNum;
				mono_readListArr[monoReadNum].matchNum = this_readListArr1[i].curNum;

				for(k=0; k<this_readListArr1[i].matchNum; k++)
				{
					if(this_readPosArr1[this_readListArr1[i].firstRow+k].contigID>0)
					{
						if(memcpy(mono_readPosArr+monoMatchItemNum+k, this_readPosArr1+this_readListArr1[i].firstRow+k, sizeof(ReadPos))==NULL)
						{
							printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
				}
				monoReadNum ++;
				monoMatchItemNum += this_readListArr1[i].curNum;
			}
		}
	}else if(j < readItemNum2)
	{
		for(; j<readItemNum2; j++)
		{
			if(this_readListArr2[j].curNum>0)
			{
				if(memcpy(mono_readListArr+monoReadNum, this_readListArr2+j, sizeof(ReadList))==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				mono_readListArr[monoReadNum].firstRow = monoMatchItemNum;
				mono_readListArr[monoReadNum].matchNum = this_readListArr2[j].curNum;

				for(k=0; k<this_readListArr2[j].matchNum; k++)
				{
					if(this_readPosArr2[this_readListArr2[j].firstRow+k].contigID>0)
					{
						if(memcpy(mono_readPosArr+monoMatchItemNum+k, this_readPosArr2+this_readListArr2[j].firstRow+k, sizeof(ReadPos))==NULL)
						{
							printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
				}
				monoReadNum ++;
				monoMatchItemNum += this_readListArr2[j].curNum;
			}
		}
	}

	readNumInRL[0] = monoReadNum;
	matchItemNumInRP[0] = monoMatchItemNum;


	//########################## Debug information ##########################
#if DEBUG_OUT_FLAG
	printf("monoReadNum=%ld, monoMatchItemNum=%ld\n", monoReadNum, monoMatchItemNum);
#endif
	//########################## Debug information ##########################

	return SUCCESSFUL;
}

/**
 * Get the readID from the read list.
 *  @return:
 *  	If find, return the row number; else, return -1.
 */
int64_t getReadRowFromReadList(const int64_t readID, const ReadList *pReadListArr, const int64_t readItemNum)
{
	int64_t mid, low, high;
	low = 0;
	high = readItemNum - 1;
	while(low<=high)
	{
		mid = (low + high) / 2;
		if(pReadListArr[mid].readID==readID)
		{
			return mid;
		}else if(pReadListArr[mid].readID<readID)
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
 * Save shared read list to a binary file.
 * 	Format:
 * 		(1) first two int64_t integers:
 * 			item number of read list array, and
 * 			the item number of read position array, respectively;
 * 		(2) read list array data;
 * 		(3) read position array data.
 *
 *  @return:
 *   If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short saveSharedReadListToFile(const char *sharedReadListFile)
{
	FILE *fpShared;
	int64_t i, totalMatchItemNum;
	int64_t curReadNumInRL, curMatchItemNumInRP;

	fpShared = fopen(sharedReadListFile, "wb");
	if(fpShared==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, sharedReadListFile);
		return FAILED;
	}

	curReadNumInRL = 0;
	curMatchItemNumInRP = 0;

	for(i=0; i<readNumInRL[0]; i++)
	{
		if(readListArr[0][i].curNum > 0)
		{
			curReadNumInRL ++;
			curMatchItemNumInRP += readListArr[0][i].curNum;
		}
	}

	// ########################## Debug information #########################
	if(curReadNumInRL > readNumInRL[0])
	{
		printf("line=%d, In %s(), curReadNum=%lu > readNum=%lu, error!\n", __LINE__, __func__, curReadNumInRL, readNumInRL[0]);
		return FAILED;
	}
	if(curMatchItemNumInRP > matchItemNumInRP[0])
	{
		printf("line=%d, In %s(), curMatchItemNum=%lu > matchItemNum=%lu, error!\n", __LINE__, __func__, curMatchItemNumInRP, matchItemNumInRP[0]);
		return FAILED;
	}
	// ########################## Debug information #########################

	// save the item numbers
	if(fwrite(&curReadNumInRL, sizeof(int64_t), 1, fpShared)!=1)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(fwrite(&curMatchItemNumInRP, sizeof(int64_t), 1, fpShared)!=1)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	// save the read list array
	totalMatchItemNum = 0;
	for(i=0; i<readNumInRL[0]; i++)
	{
		if(readListArr[0][i].curNum>0)
		{
			readListArr[0][i].matchNum = readListArr[0][i].curNum;
			readListArr[0][i].firstRow = totalMatchItemNum;

			totalMatchItemNum += readListArr[0][i].curNum;

			if(fwrite(readListArr[0]+i, sizeof(ReadList), 1, fpShared)!=1)
			{
				printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	// save the read position array
	for(i=0; i<matchItemNumInRP[0]; i++)
	{
		if(readPosArr[0][i].contigID > 0)
		{
			if(fwrite(readPosArr[0]+i, sizeof(ReadPos), 1, fpShared)!=1)
			{
				printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	fclose(fpShared);
	fpShared = NULL;

	return SUCCESSFUL;
}

/**
 * Save mono read list to a binary file.
 * 	Format:
 * 		(1) first two int64_t integers:
 * 			item number of read list array, and
 * 			the item number of read position array, respectively;
 * 		(2) read list array data;
 * 		(3) read position array data.
 *
 *  @return:
 *   If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short saveMonoReadListToFile(const char *monoReadListFile)
{
	FILE *fpMonoRL;
	int64_t i, totalMatchItemNum;
	int64_t curReadNumInRL, curMatchItemNumInRP;

	fpMonoRL = fopen(monoReadListFile, "wb");
	if(fpMonoRL==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, monoReadListFile);
		return FAILED;
	}

	curReadNumInRL = 0;
	curMatchItemNumInRP = 0;

	for(i=0; i<readNumInRL[0]; i++)
	{
		if(readListArr[0][i].curNum > 0)
		{
			curReadNumInRL ++;
			curMatchItemNumInRP += readListArr[0][i].curNum;
		}
	}

	// ########################## Debug information #########################
	if(curReadNumInRL > readNumInRL[0])
	{
		printf("line=%d, In %s(), curReadNum=%lu > readNum=%lu, error!\n", __LINE__, __func__, curReadNumInRL, readNumInRL[0]);
		return FAILED;
	}
	if(curMatchItemNumInRP > matchItemNumInRP[0])
	{
		printf("line=%d, In %s(), curMatchItemNum=%lu > matchItemNum=%lu, error!\n", __LINE__, __func__, curMatchItemNumInRP, matchItemNumInRP[0]);
		return FAILED;
	}
	// ########################## Debug information #########################

	// save the item numbers
	if(fwrite(&curReadNumInRL, sizeof(int64_t), 1, fpMonoRL)!=1)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(fwrite(&curMatchItemNumInRP, sizeof(int64_t), 1, fpMonoRL)!=1)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	// save the read list array
	totalMatchItemNum = 0;
	for(i=0; i<readNumInRL[0]; i++)
	{
		if(readListArr[0][i].curNum>0)
		{
			readListArr[0][i].matchNum = readListArr[0][i].curNum;
			readListArr[0][i].firstRow = totalMatchItemNum;

			totalMatchItemNum += readListArr[0][i].curNum;

			if(fwrite(readListArr[0]+i, sizeof(ReadList), 1, fpMonoRL)!=1)
			{
				printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	// save the read position array
	for(i=0; i<matchItemNumInRP[0]; i++)
	{
		if(readPosArr[0][i].contigID > 0)
		{
			if(fwrite(readPosArr[0]+i, sizeof(ReadPos), 1, fpMonoRL)!=1)
			{
				printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	fclose(fpMonoRL);
	fpMonoRL = NULL;

	return SUCCESSFUL;
}

/**
 * Save the read list to file.
 * 	Format:
 * 		(1) first two int64_t integers:
 * 			item number of read list array, and
 * 			the item number of read position array, respectively;
 * 		(2) read list array data;
 * 		(3) read position array data.
 *
 *  @return:
 *   If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short saveReadListsToFile(const char *readListFile1, const char *readListFile2)
{
	FILE *fpReadList[2];
	int64_t curReadNumInRL, curMatchItemNumInRP;
	int64_t i, j, totalMatchItemNum;

	fpReadList[0] = fopen(readListFile1, "wb");
	if(fpReadList[0]==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readListFile1);
		return FAILED;
	}

	fpReadList[1] = fopen(readListFile2, "wb");
	if(fpReadList[1]==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readListFile2);
		return FAILED;
	}

	for(i=0; i<2; i++)
	{
		curReadNumInRL = 0;
		curMatchItemNumInRP = 0;
		for(j=0; j<readNumInRL[i+1]; j++)
		{
			if(readListArr[i+1][j].curNum > 0)
			{
				curReadNumInRL ++;
				curMatchItemNumInRP += readListArr[i+1][j].curNum;
			}
		}

		// ########################## Debug information #########################
		if(curReadNumInRL > readNumInRL[i+1])
		{
			printf("line=%d, In %s(), curReadNum=%lu > readNum=%lu, error!\n", __LINE__, __func__, curReadNumInRL, readNumInRL[i+1]);
			return FAILED;
		}
		if(curMatchItemNumInRP > matchItemNumInRP[i+1])
		{
			printf("line=%d, In %s(), curMatchItemNum=%lu > matchItemNum=%lu, error!\n", __LINE__, __func__, curMatchItemNumInRP, matchItemNumInRP[i+1]);
			return FAILED;
		}
		// ########################## Debug information #########################

		// save the item numbers
		if(fwrite(&curReadNumInRL, sizeof(int64_t), 1, fpReadList[i])!=1)
		{
			printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(fwrite(&curMatchItemNumInRP, sizeof(int64_t), 1, fpReadList[i])!=1)
		{
			printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
			return FAILED;
		}

		// save the read list array
		totalMatchItemNum = 0;
		for(j=0; j<readNumInRL[i+1]; j++)
		{
			if(readListArr[i+1][j].curNum>0)
			{
				readListArr[i+1][j].matchNum = readListArr[i+1][j].curNum;
				readListArr[i+1][j].firstRow = totalMatchItemNum;

				totalMatchItemNum += readListArr[i+1][j].curNum;

				if(fwrite(readListArr[i+1]+j, sizeof(ReadList), 1, fpReadList[i])!=1)
				{
					printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}

		// save the read position array
		for(j=0; j<matchItemNumInRP[i+1]; j++)
		{
			if(readPosArr[i+1][j].contigID>0)
			{
				if(fwrite(readPosArr[i+1]+j, sizeof(ReadPos), 1, fpReadList[i])!=1)
				{
					printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}


	for(i=0; i<2; i++)
	{
		fclose(fpReadList[i]);
		fpReadList[i] = NULL;
	}

	return SUCCESSFUL;
}

/**
 * Allocate memory and load the data of a read list (RL) from a binary file.
 * 	File format:
 * 		(1) first two int64_t integers:
 * 			item number of read list array, and
 * 			the item number of read position array, respectively;
 * 		(2) read list array data;
 * 		(3) read position array data.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short loadSingleReadList(const char *readListFile, ReadList **pReadListArr, int64_t *readItemNum, ReadPos **pReadPosArr, int64_t *matchItemNum)
{
	int64_t tmp[2];
	FILE *fpReadList;

	fpReadList = fopen(readListFile, "rb");
	if(fpReadList==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readListFile);
		return FAILED;
	}

	// get the item numbers
	if(fread(tmp, sizeof(int64_t), 2, fpReadList)!=2)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}
	*readItemNum = tmp[0];
	*matchItemNum = tmp[1];

	// allocate the memory
	*pReadListArr = (ReadList*) malloc((*readItemNum) * sizeof(ReadList));
	if(*pReadListArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	*pReadPosArr = (ReadPos*) malloc((*matchItemNum) * sizeof(ReadPos));
	if(*pReadPosArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}


	// get the read list array
	if(fread(*pReadListArr, sizeof(ReadList), *readItemNum, fpReadList)!=(*readItemNum))
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the read position array
	if(fread(*pReadPosArr, sizeof(ReadPos), *matchItemNum, fpReadList)!=(*matchItemNum))
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}

	fclose(fpReadList);
	fpReadList = NULL;

	return SUCCESSFUL;
}

/**
 * Release the memory of a read list.
 */
void freeSingleReadList(ReadList **pReadListArr, int64_t *readItemNum, ReadPos **pReadPosArr, int64_t *matchItemNum)
{
	*readItemNum = 0;
	*matchItemNum = 0;

	free(*pReadListArr);
	*pReadListArr = NULL;

	free(*pReadPosArr);
	*pReadPosArr = NULL;
}
