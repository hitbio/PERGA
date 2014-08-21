/*
 * scafMap.c
 *
 *  Created on: Jan 11, 2014
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Start the scaffolding.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mapReads(contigGraph_t *contigGraph, readSet_t *readSet, scafContigIndex_t *scafContigIndex)
{
	int32_t i, j, seqLen, maxArraySize, matchItemNum, matchItemNumRev;
	readBlock_t *pReadBlockArray;
	read_t *pRead, *pReadArray;
	uint64_t rid, *readSeqInt, *readSeqIntRev;
	scafContigpos_t *matchResultArray, *matchResultArrayRev, *matchResultArrayBuf1, *matchResultArrayBuf2;
	char readseqTmp[MAX_READ_LEN_IN_BUF+1];
	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo;
	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;

	printf("Aligning the reads ...\n");

	// initialize read blocks
	if(initReadMatchInfoBlockInReadset(readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize read blocks, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// the match result array
	if(getMaxArraySizeFromContigIndex(&maxArraySize, scafContigIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the maximal array size from contig index, error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(maxArraySize<=0)
	{
		printf("line=%d, In %s(), invalid array size %d, error!\n", __LINE__, __func__, maxArraySize);
		return FAILED;
	}

	matchResultArray = (scafContigpos_t*) calloc (maxArraySize, sizeof(scafContigpos_t));
	matchResultArrayRev = (scafContigpos_t*) calloc (maxArraySize, sizeof(scafContigpos_t));
	matchResultArrayBuf1 = (scafContigpos_t*) calloc (maxArraySize, sizeof(scafContigpos_t));
	matchResultArrayBuf2 = (scafContigpos_t*) calloc (maxArraySize, sizeof(scafContigpos_t));
	if(matchResultArray==NULL || matchResultArrayRev==NULL || matchResultArrayBuf1==NULL || matchResultArrayBuf2==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}


	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

	rid = 0;
	for(i=0; i<readSet->blocksNumRead; i++)
	{
		pReadBlockArray = readSet->readBlockArr + i;
		pReadArray = pReadBlockArray->readArr;

		for(j=0; j<pReadBlockArray->itemNum; j++)
		{
			rid ++;
			pRead = pReadArray + j;
			seqLen = pRead->seqlen;

			if(pRead->validFlag==YES)
			//if(pRead->validFlag==YES && pRead->seqlen>=readLen)
			{
				// map single read
				if(mapSingleReadInScaf(rid, pRead, readSet, matchResultArray, matchResultArrayRev, matchResultArrayBuf1, matchResultArrayBuf2, &matchItemNum, &matchItemNumRev, scafContigIndex)==FAILED)
				{
					printf("line=%d, In %s(), cannot map the read %lu, error!\n", __LINE__, __func__, rid);
					return FAILED;
				}

				// add unique alignment
				if(matchItemNum+matchItemNumRev==1)
				{
					// save the match result
					readMatchInfoBlockID = (rid - 1) / maxItemNumPerReadMatchInfoBlock;
					rowNumInReadMatchInfoBlock = (rid - 1) % maxItemNumPerReadMatchInfoBlock;
					pReadMatchInfo = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;
					readMatchInfoBlockArr[readMatchInfoBlockID].itemNum ++;
					readSet->totalValidItemNumReadMatchInfo ++;

					if(matchItemNum>0)
					{
						pReadMatchInfo->contigID = matchResultArray[0].contigID;
						pReadMatchInfo->contigPos = matchResultArray[0].contigpos;
						pReadMatchInfo->matchlen = seqLen;
						pReadMatchInfo->readID = rid;
						pReadMatchInfo->seqlen = seqLen;
						pReadMatchInfo->readOrientation = ORIENTATION_PLUS;
						pReadMatchInfo->contigEnd = matchResultArray[0].contigEndFlag;
					}else if(matchItemNumRev>0)
					{
						pReadMatchInfo->contigID = matchResultArrayRev[0].contigID;
						pReadMatchInfo->contigPos = matchResultArrayRev[0].contigpos;
						pReadMatchInfo->matchlen = seqLen;
						pReadMatchInfo->readID = rid;
						pReadMatchInfo->seqlen = seqLen;
						pReadMatchInfo->readOrientation = ORIENTATION_MINUS;
						pReadMatchInfo->contigEnd = matchResultArrayRev[0].contigEndFlag;
					}

					if(pReadMatchInfo->contigEnd==1)
						contigGraph->contigItemArray[pReadMatchInfo->contigID-1].contigReadNumEnd3 ++;
					else
						contigGraph->contigItemArray[pReadMatchInfo->contigID-1].contigReadNumEnd5 ++;

/*
					// output the match result to file
					//if(pReadMatchInfo->contigID==15 || pReadMatchInfo->contigID==247 || pReadMatchInfo->contigID==285 || pReadMatchInfo->contigID==656 || pReadMatchInfo->contigID==646 || pReadMatchInfo->contigID==598 || pReadMatchInfo->contigID==560 || pReadMatchInfo->contigID==191 || pReadMatchInfo->contigID==323)
					if(pReadMatchInfo->contigID==52 || pReadMatchInfo->contigID==1046)
					{
						readSeqInt = readSet->readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
						seqLen = pRead->seqlen;

						if(matchItemNum>0)
						{
							getReadBaseFromPosByInt(readseqTmp, readSeqInt, seqLen, 0, seqLen);
							printf("%ld\t%d\t%d\t%d\t+\t%d\t%s\n", rid, pReadMatchInfo->contigID, pReadMatchInfo->contigEnd, pReadMatchInfo->contigPos, seqLen, readseqTmp);
						}
						else if(matchItemNumRev>0)
						{
							readSeqIntRev = (uint64_t*) calloc (((seqLen-1)>>5)+1, sizeof(uint64_t));
							if(readSeqIntRev==NULL)
							{
								printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
								return FAILED;
							}

							// generate the reverse complements
							if(getReverseReadseqInt(readSeqIntRev, readSeqInt, seqLen)==FAILED)
							{
								printf("line=%d, In %s(), cannot get the reverse complements of a read, error!\n", __LINE__, __func__);
								return FAILED;
							}

							getReadBaseFromPosByInt(readseqTmp, readSeqIntRev, seqLen, 0, seqLen);
							printf("%ld\t%d\t%d\t%d\t-\t%d\t%s\n", rid, pReadMatchInfo->contigID, pReadMatchInfo->contigEnd, pReadMatchInfo->contigPos, seqLen, readseqTmp);

							free(readSeqIntRev);
						}
					}
*/

				}
			}
		}
	}

	free(matchResultArray);
	free(matchResultArrayRev);
	free(matchResultArrayBuf1);
	free(matchResultArrayBuf2);

	return SUCCESSFUL;
}

/**
 * Get the maximal array size for the match result array for scaffolding.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxArraySizeFromContigIndex(int32_t *maxArraySize, scafContigIndex_t *scafContigIndex)
{
	int32_t i, j, maxValue;
	scafKmerBlock_t *pKmerBlock;
	scafKmer_t *kmer;

	maxValue = 0;
	for(i=0; i<scafContigIndex->blocksNumKmer; i++)
	{
		pKmerBlock = scafContigIndex->kmerBlockArr + i;
		for(j=0; j<pKmerBlock->itemNum; j++)
		{
			kmer = pKmerBlock->kmerArr + j;
			if(maxValue<kmer->arraysize)
				maxValue = kmer->arraysize;
		}
	}

	*maxArraySize = maxValue;

	return SUCCESSFUL;
}

/**
 * Map single read for the scaffolding.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mapSingleReadInScaf(int64_t rid, read_t *pRead, readSet_t *readSet, scafContigpos_t *matchResultArray, scafContigpos_t *matchResultArrayRev, scafContigpos_t *matchResultArrayBuf1, scafContigpos_t *matchResultArrayBuf2, int32_t *matchItemNum, int32_t *matchItemNumRev, scafContigIndex_t *scafContigIndex)
{
	int32_t seqLen, entriesNumRead;
	uint64_t *readSeqInt, *readSeqIntRev;

	readSeqInt = readSet->readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
	seqLen = pRead->seqlen;

	entriesNumRead = ((seqLen - 1) >> 5) + 1;

	if(seqLen>=kmerSize)
	{
		if(getMatchedContigPos(matchResultArray, matchResultArrayBuf1, matchResultArrayBuf2, matchItemNum, readSeqInt, seqLen, scafContigIndex)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the match contig position, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// the reverse complements
		if((*matchItemNum)<=1)
		{
			readSeqIntRev = (uint64_t*) calloc (entriesNumRead, sizeof(uint64_t));
			if(readSeqIntRev==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// generate the reverse complements
			if(getReverseReadseqInt(readSeqIntRev, readSeqInt, seqLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the reverse complements of a read, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// align the read
			if(getMatchedContigPos(matchResultArrayRev, matchResultArrayBuf1, matchResultArrayBuf2, matchItemNumRev, readSeqIntRev, seqLen, scafContigIndex)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the match contig position, error!\n", __LINE__, __func__);
				return FAILED;
			}

			free(readSeqIntRev);
		}
	}else
	{
		*matchItemNum = 0;
	}

	return SUCCESSFUL;
}

/**
 * Map matched contig position of a read sequence for scaffolding.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMatchedContigPos(scafContigpos_t *matchResultArray, scafContigpos_t *matchResultArrayBuf1, scafContigpos_t *matchResultArrayBuf2, int32_t *matchItemNum, uint64_t *readSeqInt, int32_t seqLen, scafContigIndex_t *scafContigIndex)
{
	int32_t i, j, k, entriesNumRead, baseNumLastEntryRead, matchItemNum1, matchItemNum2, newMatchItemNum1, validItemNum, mapTimes, startBasePos;
	uint64_t hashcode, kmerSeqInt[scafContigIndex->entriesPerKmer];
	scafKmer_t *kmer;
	int32_t *matchValidFlagArray, matchDistance, matchFlag, startRow, endRow;

	entriesNumRead = ((seqLen - 1) >> 5) + 1;
	baseNumLastEntryRead = ((seqLen - 1) % 32) + 1;

	mapTimes = (seqLen - 1) / kmerSize + 1;
	matchItemNum1 = 0;
	newMatchItemNum1 = 0;
	if(seqLen>=kmerSize)
	{
		startBasePos = 0;

		// generate the kmer integer sequence
		if(generateKmerSeqIntFromReadset(kmerSeqInt, readSeqInt, startBasePos, entriesNumRead, baseNumLastEntryRead)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// get the kmer
		hashcode = kmerhashInt(kmerSeqInt);
		kmer = getScafKmerByHash(hashcode, kmerSeqInt, scafContigIndex);

		// record the match information
		if(kmer==NULL)
		{
			matchItemNum1 = 0;
		}else
		{
			matchItemNum1 = kmer->arraysize;
			if(memcpy(matchResultArrayBuf1, kmer->ppos, matchItemNum1*sizeof(scafContigpos_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			matchValidFlagArray = (int32_t *) malloc (matchItemNum1*sizeof(int32_t));
			if(matchValidFlagArray==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			for(j=0; j<matchItemNum1; j++) matchValidFlagArray[j] = YES;

			// process the other bases
			validItemNum = matchItemNum1;
			for(i=1; i<mapTimes && validItemNum>0; i++)
			{
				if(i<mapTimes-1)
				{
					startBasePos += kmerSize;
					matchDistance = startBasePos;
				}else
				{
					startBasePos = seqLen - kmerSize;
					matchDistance = startBasePos;
				}

				// generate the kmer integer sequence
				if(generateKmerSeqIntFromReadset(kmerSeqInt, readSeqInt, startBasePos, entriesNumRead, baseNumLastEntryRead)==FAILED)
				{
					printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// get the kmer
				hashcode = kmerhashInt(kmerSeqInt);
				kmer = getScafKmerByHash(hashcode, kmerSeqInt, scafContigIndex);

				// record the match information
				if(kmer==NULL)
				{
					matchItemNum2 = 0;
				}else
				{
					matchItemNum2 = kmer->arraysize;
					if(memcpy(matchResultArrayBuf2, kmer->ppos, matchItemNum2*sizeof(scafContigpos_t))==NULL)
					{
						printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

				// get the intersection
				if(matchItemNum2==0)
				{
					matchItemNum1 = 0;
					validItemNum = 0;
					break;
				}else if(matchItemNum2>0)
				{
					for(j=0; j<matchItemNum1; j++)
					{
						if(matchValidFlagArray[j]==YES)
						{
							matchFlag = NO;

							// get the startRow and endRow of Buf2
							if(matchItemNum2>20)
							{
								if(getRowRangeMatchArray(&startRow, &endRow, matchResultArrayBuf2, matchItemNum2, matchResultArrayBuf1[j].contigID)==FAILED)
								{
									printf("line=%d, In %s(), cannot get the array range, error!\n", __LINE__, __func__);
									return FAILED;
								}
							}
							else
							{
								startRow = 0;
								endRow = matchItemNum2 - 1;
							}

							if(startRow>=0 && endRow>=0)
							{
								for(k=startRow; k<=endRow; k++)
								{
									if(matchResultArrayBuf1[j].contigID==matchResultArrayBuf2[k].contigID && matchResultArrayBuf1[j].contigpos+matchDistance==matchResultArrayBuf2[k].contigpos)
									{
										matchFlag = YES;
										break;
									}
								}
							}

							if(matchFlag==NO)
							{
								matchValidFlagArray[j] = NO;
								validItemNum --;
							}
						}
					}
				}
			} // for(i=1; i<mapTimes; i++)

			// copy the valid items
			if(validItemNum>0)
			{
				newMatchItemNum1 = 0;
				for(j=0; j<matchItemNum1; j++)
				{
					if(matchValidFlagArray[j]==YES)
					{
						if(memcpy(matchResultArray+newMatchItemNum1, matchResultArrayBuf1+j, sizeof(scafContigpos_t))==NULL)
						{
							printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
						newMatchItemNum1 ++;
					}
				}

				if(newMatchItemNum1!=validItemNum)
				{
					printf("newMatchItemNum1=%d, validItemNum=%d, error!\n", newMatchItemNum1, validItemNum);
				}
			}

			free(matchValidFlagArray);
		}
	}

	*matchItemNum = newMatchItemNum1;

	return SUCCESSFUL;
}

/**
 * Get the row range of the reads match information array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getRowRangeMatchArray(int32_t *startRow, int32_t *endRow, scafContigpos_t *matchResultArray, int32_t arraySize, int32_t contigID)
{
	int32_t left, right, middle, existFlag;

	left = 0;
	right = arraySize - 1;
	existFlag = NO;
	while(left<=right)
	{
		middle = (left + right) / 2;
		if(matchResultArray[middle].contigID==contigID)
		{
			existFlag = YES;
			break;
		}
		if(matchResultArray[middle].contigID<contigID)
			left = middle + 1;
		else
			right = middle - 1;
	}

	if(existFlag==YES)
	{
		*startRow = middle;
		while((*startRow)>0 && matchResultArray[(*startRow)-1].contigID==contigID)
			(*startRow) --;
		*endRow = middle;
		while((*endRow)<arraySize-1 && matchResultArray[(*endRow)+1].contigID==contigID)
			(*endRow) ++;
	}else
	{
		*startRow = -1;
		*endRow = -1;
	}

	return SUCCESSFUL;
}

/**
 * Fill the reads match information to contig ends for scaffolding.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillReadMatchInfoContigEnds(contigGraph_t *contigGraph, readSet_t *readSet)
{
	int32_t i, j, contigLen, readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock, maxReadsNum, readsCount;
	readBlock_t *readBlock;
	readMatchInfoBlock_t *readMatchInfoBlock;
	readMatchInfo_t *readMatchInfoArrayTmp, *pReadMatchInfo;
	contigGraphItem_t *contigItem;
	contigRead_t *contigReadArrayBuf;
	int32_t contigReadNumTmp;

	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
		if(contigGraph->contigItemArray[i].alignRegSizeEnd3>0)
			contigGraph->contigItemArray[i].onlyEnd5 = NO;
		else
			contigGraph->contigItemArray[i].onlyEnd5 = YES;
	}

	// allocate memory
	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
		contigItem = contigGraph->contigItemArray + i;
		if(contigItem->contigReadNumEnd5>0)
		{
			contigItem->contigReadArrayEnd5 = (contigRead_t *) calloc (contigItem->contigReadNumEnd5, sizeof(contigRead_t));
			if(contigItem->contigReadArrayEnd5==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
		if(contigItem->contigReadNumEnd3>0)
		{
			contigItem->contigReadArrayEnd3 = (contigRead_t *) calloc (contigItem->contigReadNumEnd3, sizeof(contigRead_t));
			if(contigItem->contigReadArrayEnd3==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
		contigItem->contigReadNumEnd5 = 0;
		contigItem->contigReadNumEnd3 = 0;
	}

	// fill the match information
	contigReadNumTmp = 0;
	for(i=0; i<readSet->blocksNumRead; i++)
	{
		readBlock = readSet->readBlockArr + i;
		readMatchInfoBlock = readSet->readMatchInfoBlockArr + i;
		for(j=0; j<readBlock->itemNum; j++)
		{
			pReadMatchInfo = readMatchInfoBlock->readMatchInfoArr + j;

//			if(pReadMatchInfo->contigID==858)
//			{
//				contigReadNumTmp ++;
//				printf("contigID=%d, contigReadNumTmp=%d, rid=%ld, contigEnd=%d\n", pReadMatchInfo->contigID, contigReadNumTmp, pReadMatchInfo->readID, pReadMatchInfo->contigEnd);
//			}

			if(pReadMatchInfo->contigID>0)
			{
				contigItem = contigGraph->contigItemArray + pReadMatchInfo->contigID - 1;
				if(pReadMatchInfo->contigEnd==1)
				{
					contigItem->contigReadArrayEnd3[contigItem->contigReadNumEnd3].readID = pReadMatchInfo->readID;
					contigItem->contigReadArrayEnd3[contigItem->contigReadNumEnd3].contigPos = pReadMatchInfo->contigPos;
					contigItem->contigReadArrayEnd3[contigItem->contigReadNumEnd3].seqlen = pReadMatchInfo->seqlen;
					contigItem->contigReadArrayEnd3[contigItem->contigReadNumEnd3].orientation = pReadMatchInfo->readOrientation;
					contigItem->contigReadArrayEnd3[contigItem->contigReadNumEnd3].contigEnd = pReadMatchInfo->contigEnd;
					contigItem->contigReadNumEnd3 ++;
				}else
				{
					contigItem->contigReadArrayEnd5[contigItem->contigReadNumEnd5].readID = pReadMatchInfo->readID;
					contigItem->contigReadArrayEnd5[contigItem->contigReadNumEnd5].contigPos = pReadMatchInfo->contigPos;
					contigItem->contigReadArrayEnd5[contigItem->contigReadNumEnd5].seqlen = pReadMatchInfo->seqlen;
					contigItem->contigReadArrayEnd5[contigItem->contigReadNumEnd5].orientation = pReadMatchInfo->readOrientation;
					contigItem->contigReadArrayEnd5[contigItem->contigReadNumEnd5].contigEnd = pReadMatchInfo->contigEnd;
					contigItem->contigReadNumEnd5 ++;
				}
			}
		}
	}

	// get the maximal array size and allocate the buffer memory
	maxReadsNum = 0;
	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
		contigItem = contigGraph->contigItemArray + i;
		if(maxReadsNum<contigItem->contigReadNumEnd5)
			maxReadsNum = contigItem->contigReadNumEnd5;
		if(maxReadsNum<contigItem->contigReadNumEnd3)
			maxReadsNum = contigItem->contigReadNumEnd3;
	}
	if(maxReadsNum<=0)
	{
		printf("line=%d, In %s(), invalid reads count %d, error!\n", __LINE__, __func__, maxReadsNum);
		return FAILED;
	}

	contigReadArrayBuf = (contigRead_t*) calloc(maxReadsNum, sizeof(contigRead_t));
	if(contigReadArrayBuf==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// radix sort the match information
	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
		contigItem = contigGraph->contigItemArray + i;
		if(contigItem->contigReadNumEnd5>0)
		{
			if(radixSortContigReadArrayInScaf(contigItem->contigReadArrayEnd5, contigReadArrayBuf, contigItem->contigReadNumEnd5)==FAILED)
			{
				printf("line=%d, In %s(), cannot sort the contig read array, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		if(contigItem->contigReadNumEnd3>0)
		{
			if(radixSortContigReadArrayInScaf(contigItem->contigReadArrayEnd3, contigReadArrayBuf, contigItem->contigReadNumEnd3)==FAILED)
			{
				printf("line=%d, In %s(), cannot sort the contig read array, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	// free contigReadArrayBuf
	free(contigReadArrayBuf);


	// remove the reads information in highly coverage regions
	if(removeMatchInfoInHighCovRegs(contigGraph, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot sort the contig read array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}


/**
 * Radix sort the reads match information at contig ends for scaffolding.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short radixSortContigReadArrayInScaf(contigRead_t *contigReadArray, contigRead_t *contigReadArrayBuf, int32_t itemNum)
{
	struct partNode
	{
		int curItemNum;
		int totalItemNum;
		int firstRow;
	};

	int64_t i, step, total;
	contigRead_t *data, *buf;
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
			buf = contigReadArray;
			data = contigReadArrayBuf;
		}else
		{
			data = contigReadArray;
			buf = contigReadArrayBuf;
		}

		// count the number
		if(memset(part, 0, partArrSize * sizeof(struct partNode))==NULL)
		{
			printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		for(i=0; i<itemNum; i++)
		{
			//part[ bitMask - ((data[i].contigPos >> step) & bitMask) ].totalItemNum ++;  // from big to small
			part[ (data[i].contigPos >> step) & bitMask ].totalItemNum ++;  // from small to big
		}

		// initialize the part array
		total = 0;
		for(i=0; i<partArrSize; i++)
		{
			part[i].firstRow = total;
			total += part[i].totalItemNum;
		}

		// copy the data to the right place
		for(i=0; i<itemNum; i++)
		{
			//hashcode = bitMask - ((data[i].contigPos >> step) & bitMask);  // from big to small
			hashcode = (data[i].contigPos >> step) & bitMask;  // from small to big
			firstRow = part[hashcode].firstRow;
			curItemNum = part[hashcode].curItemNum;
			//buf[firstRow+curItemNum].supportReadsNum = data[i].supportReadsNum;
			//buf[firstRow+curItemNum].pathItem = data[i].pathItem;
			if(memcpy(buf+firstRow+curItemNum, data+i, sizeof(contigRead_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			part[hashcode].curItemNum ++;
		}

		step += stepBits;

		//######################## Debug information #######################
//		for(i=0; i<partArrSize; i++)
//		{
//			if(part[i].curItemNum!=part[i].totalItemNum)
//			{
//				printf("line=%d, In %s(), in part[%ld], curItemNum=%d != totalItemNum=%d, error!\n", __LINE__, __func__, i, part[i].curItemNum, part[i].totalItemNum);
//				free(part);
//				return FAILED;
//			}
//		}
		//######################## Debug information #######################
	}

	// ######################## Debug information ########################
//	for(i=0; i<itemNum-1; i++)
//	{
//		if(contigReadArray[i].contigPos>contigReadArray[i+1].contigPos)
//		{
//			printf("line=%d, In %s(), contigReadArray[i].contigPos=%d, contigReadArray[i+1].contigPos=%d, itemNum=%d, error!\n", __LINE__, __func__, contigReadArray[i].contigPos, contigReadArray[i+1].contigPos, itemNum);
//			return FAILED;
//		}
//	}
	// ######################## Debug information ########################

	free(part);
	part = NULL;

	return SUCCESSFUL;
}

/**
 * Remove the reads match information at highly coverage regions.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeMatchInfoInHighCovRegs(contigGraph_t *contigGraph, readSet_t *readSet)
{
	double ratioThres;

	// compute the average coverage
	if(computeAverCovNumGlobal(contigGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute average coverage reads count, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// remove the match information in contig ends and readSet
	ratioThres = 3;
	if(removeReadsInHighCovRegs(ratioThres, contigGraph, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove the reads in high coverage regions, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Compute the average coverage reads count at contig ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeAverCovNumGlobal(contigGraph_t *contigGraph)
{
	int64_t i, alignLenSum, alingReadsNumSum, alignLenTmp, alingReadsNumTmp;

	alignLenSum = 0;
	alingReadsNumSum = 0;

	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
		if(contigGraph->contigItemArray[i].alignRegSizeEnd3>0)
		{
			alignLenSum += contigGraph->contigItemArray[i].alignRegSizeEnd5;
			alingReadsNumSum += contigGraph->contigItemArray[i].contigReadNumEnd5;

			alignLenSum += contigGraph->contigItemArray[i].alignRegSizeEnd3;
			alingReadsNumSum += contigGraph->contigItemArray[i].contigReadNumEnd3;
		}
	}

	if(alignLenSum>0)
		contigGraph->averCovNum = (double)alingReadsNumSum / alignLenSum;
	else
	{
		contigGraph->averCovNum = 0;
		printf("line=%d, In %s(), cannot averCovNum=%.2f, error!\n", __LINE__, __func__, contigGraph->averCovNum);
		return FAILED;
	}


	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{

		if(contigGraph->contigItemArray[i].alignRegSizeEnd5>0)
		{
			contigGraph->contigItemArray[i].averCovNumEnd5 = (double)contigGraph->contigItemArray[i].contigReadNumEnd5 / contigGraph->contigItemArray[i].alignRegSizeEnd5 / contigGraph->averCovNum;
		}else
		{
			contigGraph->contigItemArray[i].averCovNumEnd5 = 0;
		}
		if(contigGraph->contigItemArray[i].alignRegSizeEnd3>0)
		{
			contigGraph->contigItemArray[i].averCovNumEnd3 = (double)contigGraph->contigItemArray[i].contigReadNumEnd3 / contigGraph->contigItemArray[i].alignRegSizeEnd3 / contigGraph->averCovNum;
		}else
		{
			contigGraph->contigItemArray[i].averCovNumEnd3 = 0;
		}
	}

	return SUCCESSFUL;
}

/**
 * Remove the match information in contig ends and readSet.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeReadsInHighCovRegs(double ratioThres, contigGraph_t *contigGraph, readSet_t *readSet)
{
	int32_t i, j, k, startContigPosAlignReg, endContigPosAlignReg, startRowSubReg, endRowSubReg, subRegSizeDefault, startContigPosSubReg, endContigPosSubReg;
	int32_t rowsNumSubReg, subRegSize, contigReadNum;
	int32_t preStartRowSubReg, preEndRowSubReg, preStartContigPosSubReg, preEndContigPosSubReg, preCovNumRatio, preDeleteFlag;
	contigRead_t *contigReadArray, *newContigReadArray;
	double averCovNum, averCovNumSubReg;
	int32_t *linkArray, linkArraySize, contigID, contigLen, multiLinkFlag, deleteFlag, outputFlag;
	int32_t maxLinkedContigIDBeforeDel, secLinkedContigIDBeforeDel, maxValueLinkedContigBeforeDel, secValueLinkedContigBeforeDel, linkedContigsNumBeforeDel;
	int32_t maxLinkedContigIDAfterDel, secLinkedContigIDAfterDel, maxValueLinkedContigAfterDel, secValueLinkedContigAfterDel, linkedContigsNumAfterDel;


	linkArraySize = contigGraph->itemNumContigItemArray;
	linkArray = (int32_t *) malloc(linkArraySize*sizeof(int32_t));
	if(linkArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	subRegSizeDefault = 30;
	averCovNum = contigGraph->averCovNum;
	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
//		outputFlag = NO;
//		if(contigGraph->contigItemArray[i].contigID==1 || contigGraph->contigItemArray[i].contigID==52 || contigGraph->contigItemArray[i].contigID==1046)
//		{
//			printf("contigID=%d, contigLen=%d\n", contigGraph->contigItemArray[i].contigID, contigGraph->contigItemArray[i].contigLen);
//			outputFlag = YES;
//		}

		if(contigGraph->contigItemArray[i].alignRegSizeEnd3>0)
		{
			contigID = contigGraph->contigItemArray[i].contigID;
			contigLen = contigGraph->contigItemArray[i].contigLen;

			// 5' end
			contigReadArray = contigGraph->contigItemArray[i].contigReadArrayEnd5;
			contigReadNum = contigGraph->contigItemArray[i].contigReadNumEnd5;

			// compute the link status before removing reads
			if(getMaxSecLinkedContigs(&maxLinkedContigIDBeforeDel, &secLinkedContigIDBeforeDel, &maxValueLinkedContigBeforeDel, &secValueLinkedContigBeforeDel, &linkedContigsNumBeforeDel, linkArray, linkArraySize, sizeof(int32_t), contigID, contigReadArray, 0, contigReadNum-1, readSet)==FAILED)
			{
				printf("line=%d, In %s(), cannot compute the status for linked contigs, error!\n", __LINE__, __func__);
				return FAILED;
			}


			preStartRowSubReg = preEndRowSubReg = -1;
			preStartContigPosSubReg = preEndContigPosSubReg = -1;
			preCovNumRatio = -1;
			preDeleteFlag = -1;

			startContigPosAlignReg = 1;
			endContigPosAlignReg = contigGraph->contigItemArray[i].alignRegSizeEnd5;
			startContigPosSubReg = startContigPosAlignReg;
			startRowSubReg = 0;
			while(startRowSubReg<contigReadNum)
			{
				endContigPosSubReg = startContigPosSubReg + subRegSizeDefault - 1;
				if(endContigPosSubReg>endContigPosAlignReg)
					endContigPosSubReg = endContigPosAlignReg;
				subRegSize = endContigPosSubReg - startContigPosSubReg + 1;

				endRowSubReg = -1;
				for(j=startRowSubReg; j<contigReadNum; j++)
				{
					if(contigReadArray[j].contigPos>endContigPosSubReg)
					{
						endRowSubReg = j - 1;
						break;
					}
				}
				if(j>=contigReadNum)
					endRowSubReg = contigReadNum - 1;
				rowsNumSubReg = endRowSubReg - startRowSubReg + 1;

				// compute the average coverage of the sub-region
				if(subRegSize>0)
					averCovNumSubReg = (double)rowsNumSubReg / subRegSize;
				else
					averCovNumSubReg = 0;

//				if(outputFlag==YES)
//					printf("5' end: contigID=%d, contigLen=%d, averCovNum=%.4f, averCovNumSubReg=%.4f, covNumRatio=%.4f, rowsNumSubReg=%d:[%d,%d], subRegSize=%d:[%d,%d] \n", contigGraph->contigItemArray[i].contigID, contigGraph->contigItemArray[i].contigLen, averCovNum, averCovNumSubReg, averCovNumSubReg/averCovNum, rowsNumSubReg, startRowSubReg, endRowSubReg, subRegSize, startContigPosSubReg, endContigPosSubReg);

				// remove the reads in high coverage sub-regions
				deleteFlag = NO;
				if(averCovNumSubReg/averCovNum>ratioThres)
				{
					// check multiple links
					if(computeMultiLinkFlagOfSubReg(&multiLinkFlag, linkArray, linkArraySize, sizeof(int32_t), contigID, contigReadArray, startRowSubReg, endRowSubReg, readSet)==FAILED)
					{
						printf("line=%d, In %s(), cannot compute the multiple link flag, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(multiLinkFlag==YES)
						deleteFlag = YES;
				}
				else if((averCovNumSubReg/averCovNum>1 && endContigPosSubReg<=3*subRegSizeDefault) || (preDeleteFlag==YES && averCovNumSubReg/averCovNum>1 && startContigPosSubReg>3*subRegSizeDefault))
				{
					// check multiple links
					if(computeMultiLinkFlagOfSubReg(&multiLinkFlag, linkArray, linkArraySize, sizeof(int32_t), contigID, contigReadArray, startRowSubReg, endRowSubReg, readSet)==FAILED)
					{
						printf("line=%d, In %s(), cannot compute the multiple link flag, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(multiLinkFlag==YES)
						deleteFlag = YES;
				}

				if(deleteFlag==YES)
				{
					// mark invalid reads information in readSet
					if(markInvalidReadInfoSubReg(contigReadArray, startRowSubReg, endRowSubReg, readSet)==FAILED)
					{
						printf("line=%d, In %s(), cannot delete reads information of sub-region, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}


				// deal with the previous sub-region
				if(deleteFlag==YES && preDeleteFlag==NO && preCovNumRatio>1 && preStartContigPosSubReg>3*subRegSizeDefault)
				{
					// check multiple links
					if(computeMultiLinkFlagOfSubReg(&multiLinkFlag, linkArray, linkArraySize, sizeof(int32_t), contigID, contigReadArray, preStartRowSubReg, preEndRowSubReg, readSet)==FAILED)
					{
						printf("line=%d, In %s(), cannot compute the multiple link flag, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(multiLinkFlag==YES)
					{
						// mark invalid reads information in readSet
						if(markInvalidReadInfoSubReg(contigReadArray, preStartRowSubReg, preEndRowSubReg, readSet)==FAILED)
						{
							printf("line=%d, In %s(), cannot delete reads information of sub-region, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
				}

				preStartRowSubReg = startRowSubReg;
				preEndRowSubReg = endRowSubReg;
				preStartContigPosSubReg = startContigPosSubReg;
				preEndContigPosSubReg = endContigPosSubReg;
				preCovNumRatio = averCovNumSubReg / averCovNum;
				preDeleteFlag = deleteFlag;

				startRowSubReg += rowsNumSubReg;
				startContigPosSubReg = endContigPosSubReg + 1;
			}

			// remove the invalid items
			if(removeInvalidReadsInfoAlignReg(contigGraph->contigItemArray+i, 0)==FAILED)
			{
				printf("line=%d, In %s(), cannot remove reads information of align region, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// compute the link status after removing reads
			if(getMaxSecLinkedContigs(&maxLinkedContigIDAfterDel, &secLinkedContigIDAfterDel, &maxValueLinkedContigAfterDel, &secValueLinkedContigAfterDel, &linkedContigsNumAfterDel, linkArray, linkArraySize, sizeof(int32_t), contigID, contigGraph->contigItemArray[i].contigReadArrayEnd5, 0, contigGraph->contigItemArray[i].contigReadNumEnd5-1, readSet)==FAILED)
			{
				printf("line=%d, In %s(), cannot compute the status for linked contigs, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// confirm the reads information of linked contigs after reads removal
			if(confirmReadsInfoLinkedContigs(maxLinkedContigIDBeforeDel, secLinkedContigIDBeforeDel, maxValueLinkedContigBeforeDel, secValueLinkedContigBeforeDel, linkedContigsNumBeforeDel, maxLinkedContigIDAfterDel, secLinkedContigIDAfterDel, maxValueLinkedContigAfterDel, secValueLinkedContigAfterDel, linkedContigsNumAfterDel, contigGraph->contigItemArray+i, 0, readSet)==FAILED)
			{
				printf("line=%d, In %s(), cannot confirm the reads information for linked contigs, error!\n", __LINE__, __func__);
				return FAILED;
			}


			// 3' end
			if(contigGraph->contigItemArray[i].alignRegSizeEnd3>0 && contigGraph->contigItemArray[i].contigReadNumEnd3>0)
			{
				contigReadArray = contigGraph->contigItemArray[i].contigReadArrayEnd3;
				contigReadNum = contigGraph->contigItemArray[i].contigReadNumEnd3;

				// compute the link status before removing reads
				if(getMaxSecLinkedContigs(&maxLinkedContigIDBeforeDel, &secLinkedContigIDBeforeDel, &maxValueLinkedContigBeforeDel, &secValueLinkedContigBeforeDel, &linkedContigsNumBeforeDel, linkArray, linkArraySize, sizeof(int32_t), contigID, contigReadArray, 0, contigReadNum-1, readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute the status for linked contigs, error!\n", __LINE__, __func__);
					return FAILED;
				}

				preStartRowSubReg = preEndRowSubReg = -1;
				preStartContigPosSubReg = preEndContigPosSubReg = -1;
				preCovNumRatio = -1;
				preDeleteFlag = -1;

				startContigPosAlignReg = contigLen - contigGraph->contigItemArray[i].alignRegSizeEnd3 + 1;
				endContigPosAlignReg = contigLen;
				startContigPosSubReg = startContigPosAlignReg;
				startRowSubReg = 0;
				while(startRowSubReg<contigReadNum)
				{
					endContigPosSubReg = startContigPosSubReg + subRegSizeDefault - 1;
					if(endContigPosSubReg>endContigPosAlignReg)
						endContigPosSubReg = endContigPosAlignReg;
					subRegSize = endContigPosSubReg - startContigPosSubReg + 1;

					endRowSubReg = -1;
					for(j=startRowSubReg; j<contigReadNum; j++)
					{
						if(contigReadArray[j].contigPos>endContigPosSubReg)
						{
							endRowSubReg = j - 1;
							break;
						}
					}
					if(j>=contigReadNum)
						endRowSubReg = contigReadNum - 1;
					rowsNumSubReg = endRowSubReg - startRowSubReg + 1;

					// compute the average coverage of the sub-region
					if(subRegSize>0)
						averCovNumSubReg = (double)rowsNumSubReg / subRegSize;
					else
						averCovNumSubReg = 0;

//					if(outputFlag==YES)
//						printf("3' end: contigID=%d, contigLen=%d, averCovNum=%.4f, averCovNumSubReg=%.4f, covNumRatio=%.4f, rowsNumSubReg=%d:[%d,%d], subRegSize=%d:[%d,%d] \n", contigGraph->contigItemArray[i].contigID, contigGraph->contigItemArray[i].contigLen, averCovNum, averCovNumSubReg, averCovNumSubReg/averCovNum, rowsNumSubReg, startRowSubReg, endRowSubReg, subRegSize, startContigPosSubReg, endContigPosSubReg);

					// remove the reads in high coverage sub-regions
					deleteFlag = NO;
					if(averCovNumSubReg/averCovNum>ratioThres)
					{
						// check multiple links
						if(computeMultiLinkFlagOfSubReg(&multiLinkFlag, linkArray, linkArraySize, sizeof(int32_t), contigID, contigReadArray, startRowSubReg, endRowSubReg, readSet)==FAILED)
						{
							printf("line=%d, In %s(), cannot compute the multiple link flag, error!\n", __LINE__, __func__);
							return FAILED;
						}

						if(multiLinkFlag==YES)
							deleteFlag = YES;
					}
					else if((averCovNumSubReg/averCovNum>1 && startContigPosSubReg>=endContigPosAlignReg-readLen-3*subRegSizeDefault) || (preDeleteFlag==YES && averCovNumSubReg/averCovNum>1 && endContigPosSubReg<endContigPosAlignReg-readLen-3*subRegSizeDefault))
					{
						// check multiple links
						if(computeMultiLinkFlagOfSubReg(&multiLinkFlag, linkArray, linkArraySize, sizeof(int32_t), contigID, contigReadArray, startRowSubReg, endRowSubReg, readSet)==FAILED)
						{
							printf("line=%d, In %s(), cannot compute the multiple link flag, error!\n", __LINE__, __func__);
							return FAILED;
						}

						if(multiLinkFlag==YES)
							deleteFlag = YES;
					}

					if(deleteFlag==YES)
					{
						// mark invalid reads information at contig end and remove them in readSet
						if(markInvalidReadInfoSubReg(contigReadArray, startRowSubReg, endRowSubReg, readSet)==FAILED)
						{
							printf("line=%d, In %s(), cannot delete reads information of sub-region, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}


					// deal with the previous sub-region
					if(deleteFlag==YES && preDeleteFlag==NO && preCovNumRatio>1 && preEndContigPosSubReg<endContigPosAlignReg-readLen-3*subRegSizeDefault)
					{
						// check multiple links
						if(computeMultiLinkFlagOfSubReg(&multiLinkFlag, linkArray, linkArraySize, sizeof(int32_t), contigID, contigReadArray, preStartRowSubReg, preEndRowSubReg, readSet)==FAILED)
						{
							printf("line=%d, In %s(), cannot compute the multiple link flag, error!\n", __LINE__, __func__);
							return FAILED;
						}

						if(multiLinkFlag==YES)
						{
							// mark invalid reads information in readSet
							if(markInvalidReadInfoSubReg(contigReadArray, preStartRowSubReg, preEndRowSubReg, readSet)==FAILED)
							{
								printf("line=%d, In %s(), cannot delete reads information of sub-region, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}
					}

					preStartRowSubReg = startRowSubReg;
					preEndRowSubReg = endRowSubReg;
					preStartContigPosSubReg = startContigPosSubReg;
					preEndContigPosSubReg = endContigPosSubReg;
					preCovNumRatio = averCovNumSubReg / averCovNum;
					preDeleteFlag = deleteFlag;

					startRowSubReg += rowsNumSubReg;
					startContigPosSubReg = endContigPosSubReg + 1;
				}

				// remove the invalid items
				if(removeInvalidReadsInfoAlignReg(contigGraph->contigItemArray+i, 1)==FAILED)
				{
					printf("line=%d, In %s(), cannot remove reads information of align region, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// compute the link status after removing reads
				if(getMaxSecLinkedContigs(&maxLinkedContigIDAfterDel, &secLinkedContigIDAfterDel, &maxValueLinkedContigAfterDel, &secValueLinkedContigAfterDel, &linkedContigsNumAfterDel, linkArray, linkArraySize, sizeof(int32_t), contigID, contigGraph->contigItemArray[i].contigReadArrayEnd3, 0, contigGraph->contigItemArray[i].contigReadNumEnd3-1, readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute the status for linked contigs, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// confirm the reads information of linked contigs after reads removal
				if(confirmReadsInfoLinkedContigs(maxLinkedContigIDBeforeDel, secLinkedContigIDBeforeDel, maxValueLinkedContigBeforeDel, secValueLinkedContigBeforeDel, linkedContigsNumBeforeDel, maxLinkedContigIDAfterDel, secLinkedContigIDAfterDel, maxValueLinkedContigAfterDel, secValueLinkedContigAfterDel, linkedContigsNumAfterDel, contigGraph->contigItemArray+i, 1, readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot confirm the reads information for linked contigs, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}

	free(linkArray);

	return SUCCESSFUL;
}

/**
 * Compute the link maximal and second maximum linked contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxSecLinkedContigs(int32_t *maxLinkedContigID, int32_t *secLinkedContigID, int32_t *maxValueLinkedContig, int32_t *secValueLinkedContig, int32_t *linkedContigsNum, int32_t *linkArray, int32_t linkArraySize, int32_t elementSize, int32_t contigID, contigRead_t *contigReadArray, int32_t startRowSubReg, int32_t endRowSubReg, readSet_t *readSet)
{
	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo;
	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;
	int32_t i, linkNum, maxValue, secValue, maxContigID, secContigID;
	int64_t readID_paired;

	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

	if(memset(linkArray, 0, linkArraySize*elementSize)==NULL)
	{
		printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=startRowSubReg; i<=endRowSubReg; i++)
	{
		if(contigReadArray[i].readID%2==1)
			readID_paired = contigReadArray[i].readID + 1;
		else
			readID_paired = contigReadArray[i].readID - 1;

		readMatchInfoBlockID = (readID_paired - 1) / maxItemNumPerReadMatchInfoBlock;
		rowNumInReadMatchInfoBlock = (readID_paired - 1) % maxItemNumPerReadMatchInfoBlock;
		pReadMatchInfo = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;
		if(pReadMatchInfo->contigID>0)
		{
			linkArray[pReadMatchInfo->contigID-1] ++;
		}
	}

	linkNum = 0;
	maxValue = secValue = 0;
	maxContigID = secContigID = -1;
	for(i=0; i<linkArraySize; i++)
	{
		if(linkArray[i]>0 && i+1!=contigID)
		{
			linkNum ++;

			if(maxValue<linkArray[i])
			{
				secValue = maxValue;
				secContigID = maxContigID;
				maxValue = linkArray[i];
				maxContigID = i + 1;
			}else if(secValue<linkArray[i])
			{
				secValue = linkArray[i];
				secContigID = i + 1;
			}
		}
	}

	*maxLinkedContigID = maxContigID;
	*secLinkedContigID = secContigID;
	*maxValueLinkedContig = maxValue;
	*secValueLinkedContig = secValue;
	*linkedContigsNum = linkNum;

	return SUCCESSFUL;
}

/**
 * Confirm the reads information of linked contigs after reads removal.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short confirmReadsInfoLinkedContigs(int32_t maxLinkedContigIDBeforeDel, int32_t secLinkedContigIDBeforeDel, int32_t maxValueLinkedContigBeforeDel, int32_t secValueLinkedContigBeforeDel, int32_t linkedContigsNumBeforeDel, int32_t maxLinkedContigIDAfterDel, int32_t secLinkedContigIDAfterDel, int32_t maxValueLinkedContigAfterDel, int32_t secValueLinkedContigAfterDel, int32_t linkedContigsNumAfterDel, contigGraphItem_t *contigItem, int32_t contigEndFlag, readSet_t *readSet)
{
	int32_t i, j, delFlag;

	delFlag = 0;
	if(contigEndFlag==0)
	{ // 5' end
		if(contigItem->delReadFlagEnd5==YES && contigItem->contigReadNumEnd5>0)
		{
			if(maxValueLinkedContigAfterDel<3 && maxLinkedContigIDBeforeDel!=maxLinkedContigIDAfterDel && linkedContigsNumBeforeDel>=2 && linkedContigsNumBeforeDel!=linkedContigsNumAfterDel)
			{
				delFlag = YES;
			}

			if(delFlag==YES)
			{
				// mark invalid reads information in readSet
				if(markInvalidReadInfoSubReg(contigItem->contigReadArrayEnd5, 0, contigItem->contigReadNumEnd5-1, readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot delete reads information of align region, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// remove the invalid items
				if(removeInvalidReadsInfoAlignReg(contigItem, contigEndFlag)==FAILED)
				{
					printf("line=%d, In %s(), cannot remove reads information of align region, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}else
	{ // 3' end
		if(contigItem->delReadFlagEnd3==YES && contigItem->contigReadNumEnd3>0)
		{
			if(maxValueLinkedContigAfterDel<3 && maxLinkedContigIDBeforeDel!=maxLinkedContigIDAfterDel && linkedContigsNumBeforeDel>=2 && linkedContigsNumBeforeDel!=linkedContigsNumAfterDel)
			{
				delFlag = YES;
			}

			if(delFlag==YES)
			{
				// mark invalid reads information in readSet
				if(markInvalidReadInfoSubReg(contigItem->contigReadArrayEnd3, 0, contigItem->contigReadNumEnd3-1, readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot delete reads information of align region, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// remove the invalid items
				if(removeInvalidReadsInfoAlignReg(contigItem, contigEndFlag)==FAILED)
				{
					printf("line=%d, In %s(), cannot remove reads information of align region, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Compute the multiple link flag of a sub-region.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeMultiLinkFlagOfSubReg(int32_t *multiLinkFlag, int32_t *linkArray, int32_t linkArraySize, int32_t elementSize, int32_t contigID, contigRead_t *contigReadArray, int32_t startRowSubReg, int32_t endRowSubReg, readSet_t *readSet)
{
	int32_t linkNum, maxValue, secValue, maxContigID, secContigID;

	if(getMaxSecLinkedContigs(&maxContigID, &secContigID, &maxValue, &secValue, &linkNum, linkArray, linkArraySize, elementSize, contigID, contigReadArray, startRowSubReg, endRowSubReg, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the linked contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}

	*multiLinkFlag = NO;
	if(linkNum>5 || (linkNum>2 && maxValue>0 && secValue>=2 && (double)secValue/maxValue>0.2))
	{
		//printf("\tlinkNum=%d, maxValue=%d, secValue=%d\n", linkNum, maxValue, secValue);
		*multiLinkFlag = YES;
	}

	return SUCCESSFUL;
}

/**
 * Mark invalid reads information of a sub-region.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short markInvalidReadInfoSubReg(contigRead_t *contigReadArray, int32_t startRowSubReg, int32_t endRowSubReg, readSet_t *readSet)
{
	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo;
	int32_t i, readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;

	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

	for(i=startRowSubReg; i<=endRowSubReg; i++)
	{
		readMatchInfoBlockID = (contigReadArray[i].readID - 1) / maxItemNumPerReadMatchInfoBlock;
		rowNumInReadMatchInfoBlock = (contigReadArray[i].readID - 1) % maxItemNumPerReadMatchInfoBlock;
		pReadMatchInfo = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;
		if(pReadMatchInfo->contigID>0)
		{
			pReadMatchInfo->contigID = -1;
			pReadMatchInfo->contigPos = -1;

			readMatchInfoBlockArr[readMatchInfoBlockID].itemNum --;
			readSet->totalValidItemNumReadMatchInfo --;
		}

		//printf("\treadID=%ld, contigPos=%d, seqLen=%d, orient=%d\n", contigReadArray[i].readID, contigReadArray[i].contigPos, contigReadArray[i].seqlen, contigReadArray[i].orientation);

		// remove reads information at contig ends
		contigReadArray[i].readID = -1;
		contigReadArray[i].contigPos = -1;
	}

	return SUCCESSFUL;
}

/**
 * Output the reads match information at contig ends while scaffolding.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeInvalidReadsInfoAlignReg(contigGraphItem_t *contigItem, int32_t endFlag)
{
	int32_t i, j, contigReadNum, newItemNum;
	contigRead_t *contigReadArray, *newContigReadArray;

	if(endFlag==0)
	{
		contigReadArray = contigItem->contigReadArrayEnd5;
		contigReadNum = contigItem->contigReadNumEnd5;
	}else
	{
		contigReadArray = contigItem->contigReadArrayEnd3;
		contigReadNum = contigItem->contigReadNumEnd3;
	}

	newItemNum = 0;
	for(i=0; i<contigReadNum; i++)
	{
		if(contigReadArray[i].readID>=0)
			newItemNum ++;
	}

	if(newItemNum==0)
	{
		if(endFlag==0)
		{
			free(contigItem->contigReadArrayEnd5);
			contigItem->contigReadArrayEnd5 = NULL;
			contigItem->contigReadNumEnd5 = 0;
			contigItem->delReadFlagEnd5 = YES;
		}else
		{
			free(contigItem->contigReadArrayEnd3);
			contigItem->contigReadArrayEnd3 = NULL;
			contigItem->contigReadNumEnd3 = 0;
			contigItem->delReadFlagEnd3 = YES;
		}
	}
	else if(newItemNum<contigReadNum)
	{
		newContigReadArray = (contigRead_t*) calloc(newItemNum, sizeof(contigRead_t));
		if(newContigReadArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		j = 0;
		for(i=0; i<contigReadNum; i++)
		{
			if(contigReadArray[i].readID>0)
			{
				if(memcpy(newContigReadArray+j, contigReadArray+i, sizeof(contigRead_t))==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				j ++;
			}
		}
		if(j!=newItemNum)
		{
			printf("line=%d, In %s(), j=%d, newItemNum=%d, error!\n", __LINE__, __func__, j, newItemNum);
			return FAILED;
		}

		if(endFlag==0)
		{
			free(contigItem->contigReadArrayEnd5);
			contigItem->contigReadArrayEnd5 = newContigReadArray;
			contigItem->contigReadNumEnd5 = newItemNum;
			contigItem->delReadFlagEnd5 = YES;
		}else
		{
			free(contigItem->contigReadArrayEnd3);
			contigItem->contigReadArrayEnd3 = newContigReadArray;
			contigItem->contigReadNumEnd3 = newItemNum;
			contigItem->delReadFlagEnd3 = YES;
		}
	}

	return SUCCESSFUL;
}

/**
 * Output the reads match information at contig ends while scaffolding.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputContigReadArrayInScaf(contigGraph_t *contigGraph)
{
	int32_t i, j, contigID, contigLen, itemNum;
	char orientation;
	contigGraphItem_t *contigItem;
	contigRead_t *contigReadArray;

	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
		contigItem = contigGraph->contigItemArray + i;
		contigID = contigItem->contigID;
		contigLen = contigItem->contigLen;
		if(contigItem->contigReadNumEnd5>0)
		{
			contigReadArray = contigItem->contigReadArrayEnd5;
			itemNum = contigItem->contigReadNumEnd5;
			printf(">%d\t%d\t0\t%d\n", contigID, contigLen, itemNum);
			for(j=0; j<itemNum; j++)
			{
				if(contigReadArray[j].orientation==ORIENTATION_PLUS)
					orientation = '+';
				else
					orientation = '-';
				printf("%d\t%d\t%ld\t%c\t%d\n", contigReadArray[j].contigPos, contigReadArray[j].contigEnd, contigReadArray[j].readID, orientation, contigReadArray[j].seqlen);
			}
		}

		if(contigItem->contigReadNumEnd3>0)
		{
			contigReadArray = contigItem->contigReadArrayEnd3;
			itemNum = contigItem->contigReadNumEnd3;
			printf(">%d\t%d\t0\t%d\n", contigID, contigLen, itemNum);
			for(j=0; j<itemNum; j++)
			{
				if(contigReadArray[j].orientation==ORIENTATION_PLUS)
					orientation = '+';
				else
					orientation = '-';
				printf("%d\t%d\t%ld\t%c\t%d\n", contigReadArray[j].contigPos, contigReadArray[j].contigEnd, contigReadArray[j].readID, orientation, contigReadArray[j].seqlen);
			}
		}
	}

	return SUCCESSFUL;
}

