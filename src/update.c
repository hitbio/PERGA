#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Update decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateDecisionTable(kmertype *tmp_kmers[2], int32_t baseInt_kmer)
{
	int32_t i, posNum, basePos, baseInt_read;
	int32_t entryIDReadseq, entryPosReadseq, baseNumEntryReadseq;
	int32_t readBlockID, rowNumInReadBlock, maxItemNumPerReadBlock, seqLen, errorRegLenEnd3; // block id starts from 0
	ridpostype *ridpostable;
	assemblingreadtype *this_assemblingRead, *dtRead;

	readBlock_t *readBlockArr;
	read_t *pRead;
	readseqBlock_t *readseqBlockArr;
	uint64_t *pReadseq;
	int32_t bestRowInDT;

	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;
	readseqBlockArr = deBruijnGraph->readSet->readseqBlockArr;

	this_assemblingRead = decisionTable;
	for(i=0; i<itemNumDecisionTable; i++)
	{

		// ######################### Debug information ###########################
//		if(this_assemblingRead->rid==1112397410218)
//		{
//			printf("line=%d, In %s(), rid=%lu\n", __LINE__, __func__, (uint64_t)this_assemblingRead->rid);
//		}

//		if(this_assemblingRead->matchBaseNum+this_assemblingRead->unmatchBaseNum>=50
//			&& (double)this_assemblingRead->unmatchBaseNum/this_assemblingRead->matchBaseNum>=0.5)
//		{
//			printf("line=%d, In %s(), rid=%lu, matchBaseNum=%d, unmatchBaseNum=%d\n", __LINE__, __func__, (uint64_t)this_assemblingRead->rid, this_assemblingRead->matchBaseNum, this_assemblingRead->unmatchBaseNum);
//		}
		// ######################### Debug information ###########################


		if(this_assemblingRead->orientation==ORIENTATION_PLUS)
			basePos = this_assemblingRead->basePos + 1;
		else
			basePos = this_assemblingRead->basePos - 1;

		entryIDReadseq = basePos >> 5;
		entryPosReadseq = basePos % 32;

		if(entryIDReadseq==this_assemblingRead->entriesNumReadseq-1)
			baseNumEntryReadseq = this_assemblingRead->baseNumLastEentryReadseq;
		else
			baseNumEntryReadseq = 32;

		baseInt_read = (this_assemblingRead->readseq[entryIDReadseq] >> (2*(baseNumEntryReadseq-entryPosReadseq-1))) & 3;

		if(this_assemblingRead->orientation==ORIENTATION_PLUS)
		{ // plus orientation
			if(baseInt_read==baseInt_kmer)
			{ // base match
				this_assemblingRead->matchBaseNum ++;
				this_assemblingRead->lastMatchedBasePos = basePos;

				if(this_assemblingRead->successiveUnappearBases>0)
				{
					this_assemblingRead->successiveAppearBases = 1;
					this_assemblingRead->successiveUnappearBases = 0;
				}else
				{
					this_assemblingRead->successiveAppearBases ++;
				}
			}else
			{
				this_assemblingRead->unmatchBaseNum ++;

				if(this_assemblingRead->successiveUnappearBases==0)
				{
					this_assemblingRead->unappearBlocksNum ++;
					this_assemblingRead->successiveAppearBases = 0;
				}
				this_assemblingRead->successiveUnappearBases ++;

			}
			this_assemblingRead->basePos ++;
		}else
		{ // minus orientation
			if(((~baseInt_read)&3)==baseInt_kmer)
			{ // base match
				this_assemblingRead->matchBaseNum ++;
				this_assemblingRead->lastMatchedBasePos = basePos;

				if(this_assemblingRead->successiveUnappearBases>0)
				{
					this_assemblingRead->successiveAppearBases = 1;
					this_assemblingRead->successiveUnappearBases = 0;
				}else
				{
					this_assemblingRead->successiveAppearBases ++;
				}
			}else
			{
				this_assemblingRead->unmatchBaseNum ++;

				if(this_assemblingRead->successiveUnappearBases==0)
				{
					this_assemblingRead->unappearBlocksNum ++;
					this_assemblingRead->successiveAppearBases = 0;
				}
				this_assemblingRead->successiveUnappearBases ++;

			}
			this_assemblingRead->basePos --;
		}

		this_assemblingRead ++;
	}

	int pos;
	int returnCode, matedFlag;
	int32_t matchFlag, mismatchNum, multiCopyReadNum;
	assemblingreadtype *dtReadPaired;
	multiCopyReadNum = 0;
	//if(PEGivenType>NONE_PE_GIVEN_TYPE && readsNumInPEHashArr>=MIN_READ_NUM_PE_HASH)
	if(PEGivenType>=NONE_PE_GIVEN_TYPE)
	{
		if(tmp_kmers[0])
		{
			ridpostable = tmp_kmers[0]->ppos;
			posNum = tmp_kmers[0]->arraysize;
			for(i=0; i<posNum; i++)
			{
				pos = ridpostable[i].pos;
				if(ridpostable[i].delsign==0 && ridpostable[i].reserved==0 && pos==1)
				{

					// ############################ Debug information ##############################
					//if(ridpostable[i].rid==5687137)
					//{
					//	printf("line=%d, In %s(), i=%d, itemNumContigArr=%ld, rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, i, itemNumContigArr, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
					//}
					// ############################ Debug information ##############################

					readBlockID = (ridpostable[i].rid - 1) / maxItemNumPerReadBlock;
					rowNumInReadBlock = (ridpostable[i].rid - 1) % maxItemNumPerReadBlock;
					pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;
					pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;

					// get the exist read
					if(getExistReadInDT(&dtRead, ridpostable[i].rid, decisionTable, dtRowHashtable)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the exist read %lu in decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
						return FAILED;
					}

					if(dtRead==NULL)
					{
						// compute the mismatch base number
						if(getMismatchNumWithContigPath(&matchFlag, &mismatchNum, pReadseq, pRead->seqlen, pos, ORIENTATION_PLUS, itemNumDecisionTable, contigPath)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
							return FAILED;
						}

						// add the new read
						if(matchFlag==YES)
						{
							if(addReadToDecisionTable(ridpostable[i].rid, pos, ORIENTATION_PLUS, NO, pRead->seqlen, pReadseq, mismatchNum, itemNumContigArr)==FAILED)
							{
								printf("line=%d, In %s(), cannot add read %lu to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
								return FAILED;
							}
						}
					}else if(dtRead->mismatchNumWithContigPath>0 || dtRead->unmatchBaseNum>0)
					{
						// compute the mismatch base number
						if(getMismatchNumWithContigPath(&matchFlag, &mismatchNum, pReadseq, pRead->seqlen, ridpostable[i].pos, ORIENTATION_PLUS, itemNumDecisionTable, contigPath)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
							return FAILED;
						}

						// replace the read if its mismatch base number is smaller than the exist one
						if(mismatchNum<dtRead->mismatchNumWithContigPath)
						{
							if(replaceReadInDecisionTable(ridpostable[i].rid, pos, ORIENTATION_PLUS, NO, mismatchNum, dtRead, itemNumContigArr, contigPath)==FAILED)
							{
								printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
								return FAILED;
							}
						}

						multiCopyReadNum ++;
					}
				}//end if
				ridpostable[i].reserved = 0;
			}//end for
		}

		if(tmp_kmers[1])
		{
			ridpostable = tmp_kmers[1]->ppos;
			posNum = tmp_kmers[1]->arraysize;
			for(i=0; i<posNum; i++)
			{
				pos = ridpostable[i].pos;

				readBlockID = (ridpostable[i].rid - 1) / maxItemNumPerReadBlock;
				rowNumInReadBlock = (ridpostable[i].rid - 1) % maxItemNumPerReadBlock;
				pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;
				pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
				seqLen = pRead->seqlen;
				errorRegLenEnd3 = ceil(seqLen * ERROR_REGION_LEN_3End_FACTOR);

				if(ridpostable[i].delsign==0 && ridpostable[i].reserved==0 && pos>=seqLen-kmerSize+1-errorRegLenEnd3)
				{

					// ############################ Debug information ##############################
					//if(ridpostable[i].rid==27160337)
					//{
					//	printf("line=%d, In %s(), i=%d, itemNumContigArr=%ld, rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, i, itemNumContigArr, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
					//}
					// ############################ Debug information ##############################

					if(existReadWithPosInDecisionTable(ridpostable[i].rid, pos-1, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
					{
						returnCode = validReadPair(&dtReadPaired, ridpostable[i].rid, kmerSize, seqLen, itemNumContigArr, decisionTable, dtRowHashtable);
						if(returnCode==YES)
						{
							if(dtReadPaired)
								dtReadPaired->matedFlag = YES;
							matedFlag = YES;
						}else if(returnCode==NO)
						{
							matedFlag = NO;
						}else
						{
							printf("line=%d, In %s(), cannot valid read pair, error!\n", __LINE__, __func__);
							return FAILED;
						}

						// get the exist read
						if(getExistReadInDT(&dtRead, ridpostable[i].rid, decisionTable, dtRowHashtable)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the exist read %lu in decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
							return FAILED;
						}


						if(dtRead==NULL)
						{
							// compute the mismatch base number
							if(getMismatchNumWithContigPath(&matchFlag, &mismatchNum, pReadseq, seqLen, pos, ORIENTATION_MINUS, itemNumDecisionTable, contigPath)==FAILED)
							{
								printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
								return FAILED;
							}

							// add the new read
							if(matchFlag==YES)
							{
								if(addReadToDecisionTable(ridpostable[i].rid, pos, ORIENTATION_MINUS, matedFlag, seqLen, pReadseq, mismatchNum, itemNumContigArr)==FAILED)
								{
									printf("line=%d, In %s(), cannot add read %lu to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
									return FAILED;
								}
							}
						}else if(dtRead->mismatchNumWithContigPath>0 || dtRead->unmatchBaseNum>0)
						{
							// compute the mismatch base number
							if(getMismatchNumWithContigPath(&matchFlag, &mismatchNum, pReadseq, seqLen, pos, ORIENTATION_MINUS, itemNumDecisionTable, contigPath)==FAILED)
							{
								printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
								return FAILED;
							}

							// replace the read if its mismatch base number is smaller than the exist one
							if(mismatchNum<dtRead->mismatchNumWithContigPath)
							{
								if(replaceReadInDecisionTable(ridpostable[i].rid, pos, ORIENTATION_MINUS, matedFlag, mismatchNum, dtRead, itemNumContigArr, contigPath)==FAILED)
								{
									printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
									return FAILED;
								}
							}

							multiCopyReadNum ++;
						}
					}
				}//end if
				ridpostable[i].reserved = 0;
			}//end for
		}
	}else
	{  // single end assembly
		if(tmp_kmers[0])
		{
			ridpostable = tmp_kmers[0]->ppos;
			posNum = tmp_kmers[0]->arraysize;
			for(i=0; i<posNum; i++)
			{
				pos = ridpostable[i].pos;
				if(ridpostable[i].delsign==0 && ridpostable[i].reserved==0 && pos==1)
				{

					// ############################ Debug information ##############################
					//if(ridpostable[i].rid==5687137)
					//{
					//	printf("line=%d, In %s(), i=%d, itemNumContigArr=%ld, rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, i, itemNumContigArr, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
					//}
					// ############################ Debug information ##############################

					readBlockID = (ridpostable[i].rid - 1) / maxItemNumPerReadBlock;
					rowNumInReadBlock = (ridpostable[i].rid - 1) % maxItemNumPerReadBlock;
					pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;
					pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;

					// get the exist read
					if(getExistReadInDT(&dtRead, ridpostable[i].rid, decisionTable, dtRowHashtable)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the exist read %lu in decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
						return FAILED;
					}

					if(dtRead==NULL)
					{
						// compute the mismatch base number
						if(getMismatchNumWithContigPath(&matchFlag, &mismatchNum, pReadseq, pRead->seqlen, pos, ORIENTATION_PLUS, itemNumDecisionTable, contigPath)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
							return FAILED;
						}

						// add the new read
						if(matchFlag==YES)
						{
							if(addReadToDecisionTable(ridpostable[i].rid, pos, ORIENTATION_PLUS, NO, pRead->seqlen, pReadseq, mismatchNum, itemNumContigArr)==FAILED)
							{
								printf("line=%d, In %s(), cannot add read %lu to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
								return FAILED;
							}
						}
					}else if(dtRead->mismatchNumWithContigPath>0 || dtRead->unmatchBaseNum>0)
					{
						// compute the mismatch base number
						if(getMismatchNumWithContigPath(&matchFlag, &mismatchNum, pReadseq, pRead->seqlen, ridpostable[i].pos, ORIENTATION_PLUS, itemNumDecisionTable, contigPath)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
							return FAILED;
						}

						// replace the read if its mismatch base number is smaller than the exist one
						if(mismatchNum<dtRead->mismatchNumWithContigPath)
						{
							if(replaceReadInDecisionTable(ridpostable[i].rid, pos, ORIENTATION_PLUS, NO, mismatchNum, dtRead, itemNumContigArr, contigPath)==FAILED)
							{
								printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
								return FAILED;
							}
						}

						multiCopyReadNum ++;
					}
				}//end if
				ridpostable[i].reserved = 0;
			}//end for
		}

		if(tmp_kmers[1])
		{
			ridpostable = tmp_kmers[1]->ppos;
			posNum = tmp_kmers[1]->arraysize;
			for(i=0; i<posNum; i++)
			{
				pos = ridpostable[i].pos;

				readBlockID = (ridpostable[i].rid - 1) / maxItemNumPerReadBlock;
				rowNumInReadBlock = (ridpostable[i].rid - 1) % maxItemNumPerReadBlock;
				pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;
				pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
				seqLen = pRead->seqlen;
				errorRegLenEnd3 = ceil(seqLen * ERROR_REGION_LEN_3End_FACTOR);

				if(ridpostable[i].delsign==0 && ridpostable[i].reserved==0 && pos>=seqLen-kmerSize+1-errorRegLenEnd3)
				{

					// ############################ Debug information ##############################
					//if(ridpostable[i].rid==5687137)
					//{
					//	printf("line=%d, In %s(), i=%d, itemNumContigArr=%ld, rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, i, itemNumContigArr, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
					//}
					// ############################ Debug information ##############################

					if(existReadWithPosInDecisionTable(ridpostable[i].rid, pos-1, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
					{
						// get the exist read
						if(getExistReadInDT(&dtRead, ridpostable[i].rid, decisionTable, dtRowHashtable)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the exist read %lu in decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
							return FAILED;
						}

						if(dtRead==NULL)
						{
							// compute the mismatch base number
							if(getMismatchNumWithContigPath(&matchFlag, &mismatchNum, pReadseq, seqLen, pos, ORIENTATION_MINUS, itemNumDecisionTable, contigPath)==FAILED)
							{
								printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
								return FAILED;
							}

							// add the new read
							if(matchFlag==YES)
							{
								if(addReadToDecisionTable(ridpostable[i].rid, pos, ORIENTATION_MINUS, NO, seqLen, pReadseq, mismatchNum, itemNumContigArr)==FAILED)
								{
									printf("line=%d, In %s(), cannot add read %lu to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
									return FAILED;
								}
							}
						}else if(dtRead->mismatchNumWithContigPath>0 || dtRead->unmatchBaseNum>0)
						{
							// compute the mismatch base number
							if(getMismatchNumWithContigPath(&matchFlag, &mismatchNum, pReadseq, seqLen, pos, ORIENTATION_MINUS, itemNumDecisionTable, contigPath)==FAILED)
							{
								printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
								return FAILED;
							}

							// replace the read if its mismatch base number is smaller than the exist one
							if(mismatchNum<dtRead->mismatchNumWithContigPath)
							{
								if(replaceReadInDecisionTable(ridpostable[i].rid, pos, ORIENTATION_MINUS, NO, mismatchNum, dtRead, itemNumContigArr, contigPath)==FAILED)
								{
									printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
									return FAILED;
								}
							}

							multiCopyReadNum ++;
						}
					}
				}//end if
				ridpostable[i].reserved = 0;
			}//end for
		}
	}

	// ############################ Debug information ##############################
	//if(multiCopyReadNum>0)
	//	printf("line=%d, In %s(), multiCopyReadNum=%d, itemNumPathItemList=%d, itemNumcontigArr=%ld, itemNumDecisionTable=%d\n", __LINE__, __func__, multiCopyReadNum, contigPath->itemNumPathItemList, itemNumContigArr, itemNumDecisionTable);
	// ############################ Debug information ##############################

	return SUCCESSFUL;
}


/**
 * Update finished reads and record the succefful reads to successReadsArr.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateFinishedReadsInDecisionTable()
{
	int32_t i;

	itemNumSuccessReadsArr = 0;
	for(i = 0; i<itemNumDecisionTable; i++)
	{
		if(decisionTable[ i ].status==SUCCESSFUL_STATUS)
		{
			// ###################### Debug information #####################
			//if(decisionTable[i].rid==1707556)
			//{
			//	printf("line=%d, In %s(), rid=%ld\n", __LINE__, __func__, (int64_t)decisionTable[i].rid);
			//}
			// ###################### Debug information #####################

			successReadsArr[ itemNumSuccessReadsArr ].rid = decisionTable[i].rid;
			successReadsArr[ itemNumSuccessReadsArr ].pos = decisionTable[i].firstBasePos + 1;
			successReadsArr[ itemNumSuccessReadsArr ].matchnum = decisionTable[i].matchBaseNum;
			successReadsArr[ itemNumSuccessReadsArr ].matchlen = decisionTable[i].matchBaseNum + decisionTable[i].unmatchBaseNum;
			successReadsArr[ itemNumSuccessReadsArr ].orientation = decisionTable[i].orientation;
			successReadsArr[ itemNumSuccessReadsArr ].seqlen = decisionTable[i].seqlen;
			successReadsArr[ itemNumSuccessReadsArr ].readseq = decisionTable[i].readseq;
			successReadsArr[ itemNumSuccessReadsArr ].entriesNumReadseq = decisionTable[i].entriesNumReadseq;
			successReadsArr[ itemNumSuccessReadsArr ].baseNumLastEentryReadseq = decisionTable[i].baseNumLastEentryReadseq;
			successReadsArr[ itemNumSuccessReadsArr ].kmerNumEnd5 = decisionTable[i].kmerNumEnd5;
			successReadsArr[ itemNumSuccessReadsArr ].kmerNumEnd3 = decisionTable[i].kmerNumEnd3;

			itemNumSuccessReadsArr ++;

			if(itemNumSuccessReadsArr>=maxItemNumSuccessReadsArr)
			{
				if(reallocateSuccessReadsArr()==FAILED)
				{
					printf("line=%d, In %s(), cannot reallocate memory for success reads array, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}else if(decisionTable[ i ].status==FAILED_STATUS)
		{
		}

	}

	if(removeFinishedReadsFromDecisionTable(decisionTable, &itemNumDecisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot resort decision table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Reallocate memory for successful reads array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short reallocateSuccessReadsArr()
{
	successRead_t *pSuccessReadsArr;

	pSuccessReadsArr = (successRead_t *) malloc(2 * maxItemNumSuccessReadsArr * sizeof(successRead_t));
	if(pSuccessReadsArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory for successful reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(memcpy(pSuccessReadsArr, successReadsArr, maxItemNumSuccessReadsArr * sizeof(successRead_t))==NULL)
	{
		printf("line=%d, In %s(), cannot copy memory for successful reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	free(successReadsArr);
	successReadsArr = pSuccessReadsArr;
	maxItemNumSuccessReadsArr *= 2;

	return SUCCESSFUL;
}

/**
 * Remove the finished element in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeFinishedReadsFromDecisionTable(assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable)
{
	int32_t i , j ;
	int64_t readID_paired;
	assemblingreadtype *dtReadPaired;

	for(i=0; i<(*readsNumDecisionTable); i++)
	{
		// update the mate flag of the paired end
		if(decisionTable[i].status==FAILED_STATUS && decisionTable[i].matedFlag==YES)
		{
			// generate the paired readID
			if(decisionTable[i].rid%2==1)
			{ // odd --> even
				readID_paired = decisionTable[i].rid + 1;
			}else
			{ // even --> odd
				readID_paired = decisionTable[i].rid - 1;
			}

			// get the exist read
			if(getExistReadInDT(&dtReadPaired, readID_paired, decisionTable, dtRowHashtable)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the exist read %lu in decision table, error!\n", __LINE__, __func__, readID_paired);
				return FAILED;
			}

			if(dtReadPaired)
				dtReadPaired->matedFlag = NO;
		}
	}

	i = 0;
	j= (*readsNumDecisionTable) - 1;
	while(i <= j)
	{
		if(decisionTable[i].status !=  ASSEMBLING_STATUS && decisionTable[j].status == ASSEMBLING_STATUS)
		{

			//if(decisionTable[i].rid==32141401)
			//{
			//	printf("====== rid=%ld\n", (int64_t)decisionTable[i].rid);
			//}

			if(decisionTable[i].locked==YES)
				lockedReadsNum --;

			// delete read from the contig path item
			if(delReadFromContigPath(decisionTable+i, contigPath)==FAILED)
			{
				printf("line=%d, In %s(), cannot delete read from contig path item, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// delete the read from dtRowHashtable
			if(delReadFromDTRowHashtable(decisionTable[i].rid, i, dtRowHashtable)==FAILED)
			{
				printf("line=%d, In %s(), cannot delete read rid=%lu, dtRow=%d, error!\n", __LINE__, __func__, (int64_t)decisionTable[i].rid, i);
				return FAILED;
			}

			// update the read in dtRowHashtable
			if(updateReadInDTRowHashtable(decisionTable[j].rid, j, i, dtRowHashtable)==FAILED)
			{
				printf("line=%d, In %s(), cannot update read rid=%lu, dtRowOld=%d, dtRowNew=%d, error!\n", __LINE__, __func__, (int64_t)decisionTable[i].rid, j, i);
				return FAILED;
			}

			if(memcpy(decisionTable+i, decisionTable+j, sizeof(assemblingreadtype))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			decisionTable[i].status = ASSEMBLING_STATUS;
			(*readsNumDecisionTable) --;
			j--;
			i++;
		}

		if(decisionTable[i].status == ASSEMBLING_STATUS)
			i++;

		if(decisionTable[j].status != ASSEMBLING_STATUS)
		{
//			if(decisionTable[j].rid==487291)
//			{
//				printf("rid=%ld\n", (int64_t)decisionTable[j].rid);
//			}

			if(decisionTable[j].locked==YES)
				lockedReadsNum --;

			// delete read from the contig path item
			if(delReadFromContigPath(decisionTable+j, contigPath)==FAILED)
			{
				printf("line=%d, In %s(), cannot delete read from contig path item, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// delete the read from dtRowHashtable
			if(delReadFromDTRowHashtable(decisionTable[j].rid, j, dtRowHashtable)==FAILED)
			{
				printf("line=%d, In %s(), cannot delete read rid=%lu, dtRow=%d, error!\n", __LINE__, __func__, (int64_t)decisionTable[i].rid, i);
				return FAILED;
			}

			j--;
			(*readsNumDecisionTable) --;
		}
	}

	return SUCCESSFUL;
}

/**
 * Update reads status in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateAssemblingreadsStatus()
{
	assemblingreadtype *this_assemblingRead;
	int32_t i;

	this_assemblingRead = decisionTable;

	//int successNum = 0, failedNum = 0;
	for(i=0; i<itemNumDecisionTable; i++)
	{
		// ######################### Debug information ###########################
		//if(this_assemblingRead->rid==212805)
		//{
		//	printf("line=%d, In %s(), rid=%lu\n", __LINE__, __func__, (uint64_t)this_assemblingRead->rid);
		//}
		// ######################### Debug information ###########################

		if(this_assemblingRead->delsign==YES)
		{
			this_assemblingRead->status = FAILED_STATUS;
			//failedNum++;
		}
		else if(this_assemblingRead->unmatchBaseNum>maxUnmatchBaseNumPerRead)  // added 2013-12-16
		{
			this_assemblingRead->status = FAILED_STATUS;
		}
		else
		{
			if((this_assemblingRead->basePos==this_assemblingRead->seqlen-1 && this_assemblingRead->orientation==ORIENTATION_PLUS)
				|| (this_assemblingRead->basePos==0 && this_assemblingRead->orientation==ORIENTATION_MINUS))
			{
				if((double)this_assemblingRead->matchBaseNum/this_assemblingRead->seqlen>=readSimilarityThres)
				{


#if (DEBUG_PARA_PRINT==YES)
					// ################################ Debug information ##############################
					if(this_assemblingRead->unmatchBaseNum<11)
					{
						totalErrReadNum ++;
						errNumArr[this_assemblingRead->unmatchBaseNum] ++;
					}
					// ################################ Debug information ##############################
#endif

					if(readsNumRatio>1.5 && itemNumContigArr>10*minContigLenUsingPE)
					{
						if(this_assemblingRead->matchBaseNum>=this_assemblingRead->seqlen-1)
						{
							this_assemblingRead->status = SUCCESSFUL_STATUS;
							updateSameReadStatusToFailed(this_assemblingRead->rid, i, dtRowHashtable);
						}else
						{
							this_assemblingRead->status = FAILED_STATUS;
						}
					}else
					{
						this_assemblingRead->status = SUCCESSFUL_STATUS;
						updateSameReadStatusToFailed(this_assemblingRead->rid, i, dtRowHashtable);
					}
				}
				else
				{
					this_assemblingRead->status = FAILED_STATUS;
				}
			}
		}

		this_assemblingRead ++;
	}

	return SUCCESSFUL;
}

/**
 * Update reads status of the reads except the one having dtRow in decision table to be failed.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateSameReadStatusToFailed(int64_t rid, int32_t dtRow, dtRowIndex_t **pDtRowHashtable)
{
	dtRowIndex_t *pDtRow;

	pDtRow = pDtRowHashtable[rid & RID_LOW_BITS_MASK];
	while(pDtRow)
	{
		if(pDtRow->rid==rid && pDtRow->dtRow!=dtRow)
			decisionTable[pDtRow->dtRow].status = FAILED_STATUS;

		pDtRow = pDtRow->next;
	}

	return SUCCESSFUL;
}

/**
 * Get the candPath items from unique reads in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getCandPathFromSingleCopyReadInDT(candPath_t *candPath, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, dtRowIndex_t **dtRowHashtable)
{
	// add the read sequences from unique reads in decision table
	if(getCandPathItemsFromSingleCopyReadInDT(candPath, decisionTable, readsNumDecisionTable, dtRowHashtable)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// adjust
	if(adjustCandPath(candPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust the candPath items, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(MAX_MISMATCH_NUM_CANDIDATE_PATH>0)
	{
		// merge path items that due to sequencing errors
		if(mergeCandPathItem(candPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot merge similar candPath items, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}


/**
 * Add reads in decision table for multi-align for single ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getCandPathItemsFromSingleCopyReadInDT(candPath_t *candPath, assemblingreadtype *decisionTable, int32_t readsNumDT, dtRowIndex_t **dtRowHashtable)
{
	int32_t i, newStartPos, copyNum, readseqLenTmp, matchRow;
	char readseqTmp[MAX_CANDIDATE_PATH_LEN+1];

	// get the sequences from decision table
	candPath->itemNumCandPathItemArray = 0;
	candPath->maxPathLen = 0;
	for(i=0; i<readsNumDT; i++)
	{
		if((decisionTable[i].unmatchBaseNum<=maxUnmatchBaseNumPerRead)
			&& (decisionTable[i].successiveAppearBases>=minSuccessiveAppearedBaseNum)
			&& decisionTable[i].basePos==decisionTable[i].lastMatchedBasePos)
		{
			if(getCopyNumOfReadInDecisionTable(&copyNum, decisionTable[i].rid, dtRowHashtable)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the copy number of a read in decision table, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(copyNum==1)
			{
				if(decisionTable[i].orientation==ORIENTATION_PLUS)
				{
					newStartPos = decisionTable[i].basePos + 1;
					readseqLenTmp = decisionTable[i].seqlen - newStartPos;
					getReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, newStartPos, readseqLenTmp);
				}else
				{
					newStartPos = decisionTable[i].seqlen - decisionTable[i].basePos;
					readseqLenTmp = decisionTable[i].seqlen - newStartPos;
					getReverseReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, newStartPos, readseqLenTmp);
				}

				if(readseqLenTmp>0)
				{
					// get the matched candPathItem
					if(getMatchRowCandPathItem(&matchRow, readseqTmp, readseqLenTmp, candPath)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the matched row in candPath, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(matchRow>=0)
					{ // add the read if it matches
						if(addReadseqToCandPathItem(readseqTmp, readseqLenTmp, 1, matchRow, candPath)==FAILED)
						{
							printf("line=%d, In %s(), cannot add the read to candPath, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}else
					{ // generate a new candPath item
						if(addNewCandPathItem(candPath, readseqTmp, readseqLenTmp, 1)==FAILED)
						{
							printf("line=%d, In %s(), cannot add the contig path item, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Remove the incorrect row in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short rmIncorrectRowInDecisionTable(int32_t *bestRowInDT, int64_t rid, int32_t readseqLen, candPath_t *candPath, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable)
{
	int32_t readCopyNum;
	multiCopyReadErrNum_t *errBaseNumArray;

	//if(rid==5687137)
	//{
	//	printf("line=%d, In %s(), rid=%ld\n", __LINE__, __func__, rid);
	//}

	// initialize the memory
	if(initMemMultiCopyReadInDecisionTable(&errBaseNumArray, &readCopyNum, rid, readseqLen, decisionTable, dtRowHashtable)==FAILED)
	{
		printf("line=%d, In %s(), cannot the memory for multiple-copy reads in decision table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the read sequences
	if(fillReadseqAndErrNumArray(errBaseNumArray, readCopyNum, rid, decisionTable, dtRowHashtable)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the multiple-copy read sequences, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the correct row in decision table
	if(getBestRowMultiCopyReadInDT(bestRowInDT, errBaseNumArray, readCopyNum, candPath, decisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the correct row for the multiple-copy reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// update the decision table
	if(remainBestCopyReadInDT(*bestRowInDT, errBaseNumArray, readCopyNum, decisionTable, readsNumDecisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove the incorrect row for the multiple-copy reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	releaseMemMultiCopyReadInDecisionTable(&errBaseNumArray, &readCopyNum);

	return SUCCESSFUL;
}


/**
 * Initialize the memory for multiple-copy reads in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemMultiCopyReadInDecisionTable(multiCopyReadErrNum_t **errBaseNumArray, int32_t *readCopyNum, int64_t rid, int32_t readseqLen, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable)
{
	int32_t i;

	if(getCopyNumOfReadInDecisionTable(readCopyNum, rid, dtRowHashtable)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the read copy number and the maximal sequence length, error!\n", __LINE__, __func__);
		return FAILED;
	}

	*errBaseNumArray = (multiCopyReadErrNum_t *) calloc(*readCopyNum, sizeof(multiCopyReadErrNum_t));
	if((*errBaseNumArray)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<(*readCopyNum); i++)
	{
		(*errBaseNumArray)[i].hangingSeq= (char*) calloc(readseqLen+1, sizeof(char));
		if((*errBaseNumArray)[i].hangingSeq==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Release the memory for multiple-copy reads in decision table.
 */
void releaseMemMultiCopyReadInDecisionTable(multiCopyReadErrNum_t **errBaseNumArray, int32_t *readCopyNum)
{
	int32_t i;

	for(i=0; i<(*readCopyNum); i++)
		free((*errBaseNumArray)[i].hangingSeq);
	free(*errBaseNumArray);
	*errBaseNumArray = NULL;
	*readCopyNum = 0;
}

/**
 * Fill the multiple-copy reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillReadseqAndErrNumArray(multiCopyReadErrNum_t *errBaseNumArray, int32_t readCopyNum, int64_t rid, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable)
{
	int32_t hashcode, itemNum, newStartPos;
	dtRowIndex_t *pDtRow;
	assemblingreadtype *pReadDt;

	hashcode = rid & RID_LOW_BITS_MASK;
	pDtRow = dtRowHashtable[hashcode];

	itemNum = 0;
	while(pDtRow)
	{
		if(pDtRow->rid==rid)
		{
			// get the read sequence
			pReadDt = decisionTable + pDtRow->dtRow;
			if(pReadDt->orientation==ORIENTATION_PLUS)
			{ // plus orientation
				newStartPos = pReadDt->basePos + 1;
				getReadBaseFromPosByInt(errBaseNumArray[itemNum].hangingSeq, pReadDt->readseq, pReadDt->seqlen, newStartPos, pReadDt->seqlen-newStartPos);
			}else
			{ // minus orientation
				newStartPos = pReadDt->seqlen - pReadDt->basePos;
				getReverseReadBaseFromPosByInt(errBaseNumArray[itemNum].hangingSeq, pReadDt->readseq, pReadDt->seqlen, newStartPos, pReadDt->seqlen-newStartPos);
			}

			errBaseNumArray[itemNum].readID = rid;
			errBaseNumArray[itemNum].errBaseNum = 0;
			errBaseNumArray[itemNum].seqLen = pReadDt->seqlen - newStartPos;
			errBaseNumArray[itemNum].rowInDecisionTable = pDtRow->dtRow;
			itemNum ++;
		}

		pDtRow = pDtRow->next;
	}

	if(itemNum!=readCopyNum)
	{
		printf("line=%d, In %s(), itemNum=%d != readCopyNum=%d, error!\n", __LINE__, __func__, itemNum, readCopyNum);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get the correct row of the multiple-copy reads in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getBestRowMultiCopyReadInDT(int32_t *bestRow, multiCopyReadErrNum_t *errBaseNumArray, int32_t readCopyNum, candPath_t *candPath, assemblingreadtype *decisionTable)
{
	int32_t i, j, k, seqLen, pathLen, shareLen, *errNumArray, minErrNum, rowMinErrNum, arrayNum;
	char *pSeq, *pPathseq;
	assemblingreadtype *pReadDt;

	errNumArray = (int32_t *) calloc(candPath->itemNumCandPathItemArray, sizeof(int32_t));
	if(errNumArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<readCopyNum; i++)
	{
		pSeq = errBaseNumArray[i].hangingSeq;
		seqLen = errBaseNumArray[i].seqLen;
		pReadDt = decisionTable + errBaseNumArray[i].rowInDecisionTable;
		if(candPath->itemNumCandPathItemArray<=20)
			arrayNum = candPath->itemNumCandPathItemArray;
		else
			arrayNum = 20;
		//for(j=0; j<candPath->itemNumCandPathItemArray; j++)
		for(j=0; j<arrayNum; j++)
		{
			pPathseq = candPath->candPathItemArray[j].candPathStr;
			pathLen = candPath->candPathItemArray[j].pathLen;
			if(seqLen<pathLen)
				shareLen = seqLen;
			else
				shareLen = pathLen;

			errNumArray[j] = 0;
			for(k=0; k<shareLen; k++)
			{
				if(pSeq[k]!=pPathseq[k])
					errNumArray[j] ++;
			}
			errNumArray[j] += pReadDt->unmatchBaseNum;

			if(pReadDt->orientation==ORIENTATION_PLUS)
				errNumArray[j] += pReadDt->firstBasePos;
			else
				errNumArray[j] += pReadDt->seqlen - 1 - pReadDt->firstBasePos;
		}

		// get the minimal error number
		minErrNum = INT_MAX;
		//for(j=0; j<candPath->itemNumCandPathItemArray; j++)
		for(j=0; j<arrayNum; j++)
			if(minErrNum>errNumArray[j])
				minErrNum = errNumArray[j];

		errBaseNumArray[i].errBaseNum = minErrNum;
	}

	minErrNum = INT_MAX;
	rowMinErrNum = -1;
	for(i=0; i<readCopyNum; i++)
	{
		if(minErrNum>errBaseNumArray[i].errBaseNum)
		{
			minErrNum = errBaseNumArray[i].errBaseNum;
			rowMinErrNum = i;
		}else if(minErrNum==errBaseNumArray[i].errBaseNum)
		{
			if(errBaseNumArray[rowMinErrNum].seqLen<errBaseNumArray[i].seqLen)
				rowMinErrNum = i;
		}
	}

	*bestRow = rowMinErrNum;

	free(errNumArray);

	return SUCCESSFUL;
}

/**
 * Remain the best row of the multiple-copy reads in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short remainBestCopyReadInDT(int32_t bestRow, multiCopyReadErrNum_t *errBaseNumArray, int32_t readCopyNum, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable)
{
	int32_t i;

	if(bestRow>=0)
	{
		for(i=0; i<readCopyNum; i++)
			if(i!=bestRow)
				decisionTable[errBaseNumArray[i].rowInDecisionTable].status = FAILED_STATUS;

		if(removeFinishedReadsFromDecisionTable(decisionTable, readsNumDecisionTable)==FAILED)
		{
			printf("line=%d, In %s(), cannot remove the error read copy in decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		printf("line=%d, In %s(), invalid bestRow=%d, error!\n", __LINE__, __func__, bestRow);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get the mismatched base number of a read compared to the contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMismatchNumWithContigPath(int32_t *matchFlag, int32_t *mismatchNum, uint64_t *readseqInt, int32_t readseqLen, int32_t kmerpos, int32_t orientation, int32_t itemNumDT, contigPath_t *contigPath)
{
	int32_t j, newStartPos, tailMismatchNum, mismatchNumTmp, minMismatchNum, seqLen, pathLen, shareLen, checkNum;
	char readseq[MAX_READ_LEN_IN_BUF+1], *pathseq;
	contigPathItem_t *pathItem;
	int32_t maxMismatchThres;

	if(contigPath->itemNumPathItemList==0)
	{
		*mismatchNum = 0;
		*matchFlag = YES;
		return SUCCESSFUL;
	}

	if(orientation==ORIENTATION_PLUS)
	{
		newStartPos = kmerpos + kmerSize - 2;
		seqLen = readseqLen - newStartPos;
		getReadBaseFromPosByInt(readseq, readseqInt, readseqLen, newStartPos, seqLen);
	}else
	{
		newStartPos = readseqLen - kmerpos;
		seqLen = readseqLen - newStartPos;
		getReverseReadBaseFromPosByInt(readseq, readseqInt, readseqLen, newStartPos, seqLen);
	}
	tailMismatchNum = newStartPos - kmerSize + 1;

//	if(itemNumDT>100*averKmerOcc)
//		maxMismatchThres = contigPath->maxMismatchNumThres / 2;
//	else
		maxMismatchThres = contigPath->maxMismatchNumThres;

	*matchFlag = NO;
	checkNum = 0;
	minMismatchNum = INT_MAX;
	pathItem = contigPath->contigPathItemList;
	while(pathItem)
	{
		if(pathItem->contigPathLen-contigPath->startRowNewBase>0)
		//if(pathItem->supportReadsNum>0 && pathItem->contigPathLen-contigPath->startRowNewBase>0)
		{
			mismatchNumTmp = 0;
			pathseq = pathItem->contigPathStr + contigPath->startRowNewBase;
			pathLen = pathItem->contigPathLen - contigPath->startRowNewBase;
			if(seqLen<pathLen)
				shareLen = seqLen;
			else
				shareLen = pathLen;

			for(j=0; j<shareLen; j++)
			{
				if(readseq[j]!=pathseq[j])
				{
					mismatchNumTmp ++;
					//if(mismatchNumTmp>contigPath->mismatchFactor*pathLen)
					if(mismatchNumTmp>maxMismatchThres)
						break;
				}
			}

			if(minMismatchNum>mismatchNumTmp)
				minMismatchNum = mismatchNumTmp;
			//if(minMismatchNum==0)
			//if(minMismatchNum+tailMismatchNum<contigPath->mismatchFactor*pathLen)
			//if(minMismatchNum+tailMismatchNum<maxMismatchThres)
			if(minMismatchNum+tailMismatchNum<maxMismatchThres) // deleted 2014-03-18
			//if(minMismatchNum+tailMismatchNum<maxMismatchThres && minMismatchNum+tailMismatchNum<contigPath->mismatchFactor*readseqLen) // added 2014-03-18
			{
				*matchFlag = YES;
				break;
			}
		}

		checkNum ++;
		if(checkNum>contigPath->bestItemNumPathItemList)
			break;

		pathItem = pathItem->nextPathItem;
	}

	//if((*matchFlag)==NO && tailMismatchNum>0 && minMismatchNum<=contigPath->mismatchFactor*(contigPath->maxContigPathLen-contigPath->startRowNewBase))
	if((*matchFlag)==NO && tailMismatchNum>0 && minMismatchNum<=maxMismatchThres) // deleted 2013-03-18
	//if((*matchFlag)==NO && tailMismatchNum>0 && minMismatchNum<=maxMismatchThres && minMismatchNum<contigPath->mismatchFactor*readseqLen) // added 2013-03-18
	{
		// get the tail mismatch base number of read
		if(getReadtailMismatchNumContigPath(&tailMismatchNum, readseqInt, readseqLen, kmerpos, orientation, contigPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the tail mismatch base number of read, error!\n", __LINE__, __func__);
			return FAILED;
		}

		//if(minMismatchNum+tailMismatchNum<=contigPath->mismatchFactor*(contigPath->maxContigPathLen-contigPath->startRowNewBase))
		if(minMismatchNum+tailMismatchNum<=maxMismatchThres) // deleted 2013-03-18
		//if(minMismatchNum+tailMismatchNum<=maxMismatchThres && minMismatchNum+tailMismatchNum<contigPath->mismatchFactor*readseqLen) // added 2013-03-18
			*matchFlag = YES;
	}

	*mismatchNum = minMismatchNum + tailMismatchNum;

//	if(readseqLen<40 && (*mismatchNum)>contigPath->mismatchFactor*readseqLen && (*matchFlag)==YES)
//	{
//		*matchFlag = NO;
//		printf("******** localContigID=%ld, contigNodesNum=%ld, itemNumDT=%d, readseqLen=%d, mismatchNum=%d\n", localContigID, itemNumContigArr, itemNumDecisionTable, readseqLen, *mismatchNum);
//	}

	return SUCCESSFUL;
}

/**
 * Get the tail mismatched base number of a read compared to the contig path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadtailMismatchNumContigPath(int32_t *tailMismatchNum, uint64_t *readseqInt, int32_t readseqLen, int32_t kmerpos, int32_t orientation, contigPath_t *contigPath)
{
	int32_t i, tailReadBaseNum, tailBaseNum, startRowReadtail, startRowContigtail;
	char readtailSeq[MAX_READ_LEN_IN_BUF+1], *contigtailSeq;

	if(orientation==ORIENTATION_PLUS)
		tailReadBaseNum = kmerpos - 1;
	else
		tailReadBaseNum = readseqLen - (kmerpos + kmerSize - 1);

	if(tailReadBaseNum>0)
	{
		startRowContigtail = contigPath->contigtailSeqLen - kmerSize + 1 - tailReadBaseNum;
		if(startRowContigtail>=0)
		{
			startRowReadtail = 0;
			tailBaseNum = tailReadBaseNum;
		}else
		{
			startRowReadtail = -startRowContigtail;
			tailBaseNum = tailReadBaseNum - startRowReadtail;
			startRowContigtail = 0;
		}
		if(tailBaseNum<0)
		{
			printf("line=%d, In %s(), startRowReadtail=%d, startRowContigtail=%d, tailBaseNum=%d\n", __LINE__, __func__, startRowReadtail, startRowContigtail, tailBaseNum);
			return FAILED;
		}

		if(orientation==ORIENTATION_PLUS)
			getReadBaseFromPosByInt(readtailSeq, readseqInt, readseqLen, startRowReadtail, tailBaseNum);
		else
			getReverseReadBaseFromPosByInt(readtailSeq, readseqInt, readseqLen, startRowReadtail, tailBaseNum);

		contigtailSeq = contigPath->contigtailSeq + startRowContigtail;

		*tailMismatchNum = 0;
		for(i=0; i<tailBaseNum; i++)
		{
			if(readtailSeq[i]!=contigtailSeq[i])
			{
				(*tailMismatchNum) ++;
				//if((*tailMismatchNum)>contigPath->mismatchFactor*(contigPath->maxContigPathLen-contigPath->startRowNewBase))
				if((*tailMismatchNum)>contigPath->maxMismatchNumThres)
					break;
			}
		}
	}else
	{
		*tailMismatchNum = 0;
	}

	return SUCCESSFUL;
}
