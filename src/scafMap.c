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
	scafContigpos_t *matchResultArray, *matchResultArrayRev, *matchResultArrayBuf;
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
	matchResultArrayBuf = (scafContigpos_t*) calloc (maxArraySize, sizeof(scafContigpos_t));
	if(matchResultArray==NULL || matchResultArrayRev==NULL || matchResultArrayBuf==NULL)
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

			//if(pRead->validFlag==YES)
			if(pRead->validFlag==YES && pRead->seqlen>=readLen)
			{
				// map single read
				if(mapSingleReadInScaf(rid, pRead, readSet, matchResultArray, matchResultArrayRev, matchResultArrayBuf, &matchItemNum, &matchItemNumRev, scafContigIndex)==FAILED)
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
					if(pReadMatchInfo->contigID==21 || pReadMatchInfo->contigID==42 || pReadMatchInfo->contigID==58 || pReadMatchInfo->contigID==87)
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
	free(matchResultArrayBuf);

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
short mapSingleReadInScaf(int64_t rid, read_t *pRead, readSet_t *readSet, scafContigpos_t *matchResultArray, scafContigpos_t *matchResultArrayRev, scafContigpos_t *matchResultArrayBuf, int32_t *matchItemNum, int32_t *matchItemNumRev, scafContigIndex_t *scafContigIndex)
{
	int32_t i, j, seqLen, entriesNumRead;
	uint64_t *readSeqInt, *readSeqIntRev;

	readSeqInt = readSet->readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
	seqLen = pRead->seqlen;

	entriesNumRead = ((seqLen - 1) >> 5) + 1;

	if(seqLen>=kmerSize)
	{
		if(getMatchedContigPos(matchResultArray, matchResultArrayBuf, matchItemNum, readSeqInt, seqLen, scafContigIndex)==FAILED)
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
			if(getMatchedContigPos(matchResultArrayRev, matchResultArrayBuf, matchItemNumRev, readSeqIntRev, seqLen, scafContigIndex)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the match contig position, error!\n", __LINE__, __func__);
				return FAILED;
			}

			free(readSeqIntRev);
		}

	}

	return SUCCESSFUL;
}

/**
 * Map matched contig position of a read sequence for scaffolding.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMatchedContigPos(scafContigpos_t *matchResultArray, scafContigpos_t *matchResultArrayBuf, int32_t *matchItemNum, uint64_t *readSeqInt, int32_t seqLen, scafContigIndex_t *scafContigIndex)
{
	int32_t i, j, k, entriesNumRead, baseNumLastEntryRead, matchItemNum1, matchItemNum2, mapTimes, startBasePos;
	uint64_t hashcode, kmerSeqInt[scafContigIndex->entriesPerKmer];
	scafKmer_t *kmer;
	int32_t *matchValidFlagArray, matchDistance;

	entriesNumRead = ((seqLen - 1) >> 5) + 1;
	baseNumLastEntryRead = ((seqLen - 1) % 32) + 1;

	mapTimes = (seqLen - 1) / kmerSize + 1;

	matchItemNum1 = 0;
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
			if(memcpy(matchResultArray, kmer->ppos, matchItemNum1*sizeof(scafContigpos_t))==NULL)
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

			// process the other bases
			for(i=1; i<mapTimes; i++)
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

				for(j=0; j<matchItemNum1; j++) matchValidFlagArray[j] = NO;

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
					if(memcpy(matchResultArrayBuf, kmer->ppos, matchItemNum2*sizeof(scafContigpos_t))==NULL)
					{
						printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

				// get the intersection
				if(matchItemNum2==0)
				{
					matchItemNum1 = 0;
					break;
				}else if(matchItemNum2>0)
				{
					for(j=0; j<matchItemNum1; j++)
					{
						for(k=0; k<matchItemNum2; k++)
						{
							if(matchResultArray[j].contigID==matchResultArrayBuf[k].contigID && matchResultArray[j].contigpos+matchDistance==matchResultArrayBuf[k].contigpos)
							{
								matchValidFlagArray[j] = YES;
							}
						}
					}

					j = 0;
					while(j<matchItemNum1)
					{
						if(matchValidFlagArray[j]==NO)
						{
							for(k=j; k<matchItemNum1-1; k++)
							{
								if(memcpy(matchResultArray+k, matchResultArray+k+1, sizeof(scafContigpos_t))==NULL)
								{
									printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
									return FAILED;
								}
								matchValidFlagArray[k] = matchValidFlagArray[k+1];
							}

							matchItemNum1 --;
						}else
						{
							j ++;
						}
					}
				}
			} // for(i=1; i<mapTimes; i++)
			free(matchValidFlagArray);
		}
	}

	*matchItemNum = matchItemNum1;

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

//			if(pReadMatchInfo->contigID==1)
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









