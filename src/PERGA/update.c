#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Update decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateDecisionTable(kmertype *tmp_kmers[2], int32_t baseInt_kmer)
{
	int32_t i, posNum, basePos, baseInt_read, adjustFlag;
	int32_t entryIDReadseq, entryPosReadseq, baseNumEntryReadseq;
	int32_t readBlockID, rowNumInReadBlock, maxItemNumPerReadBlock, seqLen, errorRegLenEnd3; // block id starts from 0
	ridpostype *ridpostable;
	assemblingreadtype *this_assemblingRead;

	readBlock_t *readBlockArr;
	read_t *pRead;
	readseqBlock_t *readseqBlockArr;
	uint64_t *pReadseq;

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

				//if(this_assemblingRead->errBaseNum<maxErrBaseNumInCorrection)
				//{
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].basePos = basePos;
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].newBase = baseInt_kmer;
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].errKind = 1;
				//}
				this_assemblingRead->errBaseNum ++;
				this_assemblingRead->newErrNumAfterCorrection ++;

				if(this_assemblingRead->successiveUnappearBases==0)
				{
					this_assemblingRead->unappearBlocksNum ++;
					this_assemblingRead->successiveAppearBases = 0;
				}
				this_assemblingRead->successiveUnappearBases ++;

#if (DEBUG_CONTIG_CHECK==YES)
				// ########################## Debug information ######################
				if(this_assemblingRead->errBaseNum>(MIN_UNMATCH_BASE_NUM_ALIGN_FACTOR+1)*MAX_UNMATCH_BASE_NUM_PER_READ)
				{
					printf("line=%d, In %s(), errBaseNum=%d > %d, error!\n", __LINE__, __func__, this_assemblingRead->errBaseNum, 2*maxErrBaseNumInCorrection);
					return FAILED;
				}
				// ########################## Debug information ######################
#endif

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

				//if(this_assemblingRead->errBaseNum<maxErrBaseNumInCorrection)
				//{
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].basePos = basePos;
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].newBase = (~baseInt_kmer) & 3;
					this_assemblingRead->errBaseArr[this_assemblingRead->errBaseNum].errKind = 1;
				//}
				this_assemblingRead->errBaseNum ++;
				this_assemblingRead->newErrNumAfterCorrection ++;

				if(this_assemblingRead->successiveUnappearBases==0)
				{
					this_assemblingRead->unappearBlocksNum ++;
					this_assemblingRead->successiveAppearBases = 0;
				}
				this_assemblingRead->successiveUnappearBases ++;

#if (DEBUG_CONTIG_CHECK==YES)
				// ########################## Debug information ######################
				if(this_assemblingRead->errBaseNum>(MIN_UNMATCH_BASE_NUM_ALIGN_FACTOR+1)*MAX_UNMATCH_BASE_NUM_PER_READ)
				{
					printf("line=%d, In %s(), errBaseNum=%d > %d, error!\n", __LINE__, __func__, this_assemblingRead->errBaseNum, 2*maxErrBaseNumInCorrection);
					return FAILED;
				}
				// ########################## Debug information ######################
#endif

			}
			this_assemblingRead->basePos --;
		}

		adjustFlag = NO;
		if(this_assemblingRead->unmatchBaseNum>minUnmatchBaseNumAlign && this_assemblingRead->alignNum==0)
			adjustFlag = YES;
//		else if(this_assemblingRead->orientation==ORIENTATION_PLUS && this_assemblingRead->basePos==this_assemblingRead->seqlen-1 && this_assemblingRead->newErrNumAfterCorrection>=2)
//			adjustFlag = YES;
//		else if(this_assemblingRead->orientation==ORIENTATION_MINUS && this_assemblingRead->basePos==0 && this_assemblingRead->newErrNumAfterCorrection>=2)
//			adjustFlag = YES;

		// adjust the bases by alignment
		if(adjustFlag==YES)
		{

			// ################################## Debug information ###############################
			//if(this_assemblingRead->rid==1371453)
			//{
			//	printf("line=%d, In %s(), rid=%ld\n", __LINE__, __func__, (int64_t)this_assemblingRead->rid);
			//}
			// ################################## Debug information ###############################

			// adjust the match information of a read
			if(adjustMatchInfoRead(this_assemblingRead, contigArr, itemNumContigArr)==FAILED)
			{
				printf("line=%d, In %s(), cannot adjust the match information of a read, error!\n", __LINE__, __func__);
				return FAILED;
			}

			this_assemblingRead->alignNum ++;
			this_assemblingRead->alignSuccessTimesOld = this_assemblingRead->alignSuccessTimes;
		}

		this_assemblingRead ++;
	}

	int pos;
	int returnCode, matedFlag;
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
					//if(ridpostable[i].rid==212805)
					//{
					//	printf("line=%d, In %s(), rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
					//}
					// ############################ Debug information ##############################

					readBlockID = (ridpostable[i].rid - 1) / maxItemNumPerReadBlock;
					rowNumInReadBlock = (ridpostable[i].rid - 1) % maxItemNumPerReadBlock;
					pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;
					pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;

					if(addReadToDecisionTable(ridpostable[i].rid, pos, ORIENTATION_PLUS, NO, pRead->seqlen, pReadseq)==FAILED)
					{
						printf("line=%d, In %s(), cannot add read %lu to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
						return FAILED;
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
					//if(ridpostable[i].rid==1707556)
					//{
					//	printf("line=%d, In %s(), rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
					//}
					// ############################ Debug information ##############################

					if(existReadInDecisionTable(ridpostable[i].rid, pos-1, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
					{
						returnCode = validReadPair(ridpostable[i].rid);
						if(returnCode==YES)
						{
							matedFlag = YES;
						}else if(returnCode==NO)
						{
							matedFlag = NO;
						}else
						{
							printf("line=%d, In %s(), cannot valid read pair, error!\n", __LINE__, __func__);
							return FAILED;
						}


						if(addReadToDecisionTable(ridpostable[i].rid, pos, ORIENTATION_MINUS, matedFlag, seqLen, pReadseq)==FAILED)
						{
							printf("line=%d, In %s(), cannot add read %lu to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
							return FAILED;
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
					//if(ridpostable[i].rid==212805)
					//{
					//	printf("line=%d, In %s(), rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
					//}
					// ############################ Debug information ##############################

					readBlockID = (ridpostable[i].rid - 1) / maxItemNumPerReadBlock;
					rowNumInReadBlock = (ridpostable[i].rid - 1) % maxItemNumPerReadBlock;
					pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;
					pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;

					if(addReadToDecisionTable(ridpostable[i].rid, pos, ORIENTATION_PLUS, NO, pRead->seqlen, pReadseq)==FAILED)
					{
						printf("line=%d, In %s(), cannot add read %lu to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
						return FAILED;
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
					//if(ridpostable[i].rid==612747)
					//{
					//	printf("line=%d, In %s(), rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
					//}
					// ############################ Debug information ##############################

					if(existReadInDecisionTable(ridpostable[i].rid, pos-1, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
					{
						if(addReadToDecisionTable(ridpostable[i].rid, pos, ORIENTATION_MINUS, NO, pRead->seqlen, pReadseq)==FAILED)
						{
							printf("line=%d, In %s(), cannot add read %lu to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
							return FAILED;
						}
					}
				}//end if
				ridpostable[i].reserved = 0;
			}//end for
		}
	}

	// ############################ Debug information ##############################
	//printf("In %s(), new reads Num=%d, numassemblingreads=%d\n", __func__, newNum, itemNumDecisionTable);
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

	if(removeFinishedReadsFromDecisionTable()==FAILED)
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
short removeFinishedReadsFromDecisionTable()
{
	int32_t i , j ;
	i = 0;
	j= itemNumDecisionTable - 1;
	while(i <= j)
	{
		if(decisionTable[i].status !=  ASSEMBLING_STATUS && decisionTable[j].status == ASSEMBLING_STATUS)
		{

			//if(decisionTable[i].rid==1707556)
			//{
			//	printf("rid=%ld\n", (int64_t)decisionTable[i].rid);
			//}

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
			itemNumDecisionTable --;
			j--;
			i++;
		}

		if(decisionTable[i].status == ASSEMBLING_STATUS)
			i++;

		if(decisionTable[j].status != ASSEMBLING_STATUS)
		{
			//if(decisionTable[j].rid==1707556)
			//{
			//	printf("rid=%ld\n", (int64_t)decisionTable[j].rid);
			//}

			// delete the read from dtRowHashtable
			if(delReadFromDTRowHashtable(decisionTable[j].rid, j, dtRowHashtable)==FAILED)
			{
				printf("line=%d, In %s(), cannot delete read rid=%lu, dtRow=%d, error!\n", __LINE__, __func__, (int64_t)decisionTable[i].rid, i);
				return FAILED;
			}

			j--;
			itemNumDecisionTable --;
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
	assemblingreadtype *this_assemblingRead = decisionTable;
	int32_t i;

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
		}else
		{
			if((this_assemblingRead->basePos==this_assemblingRead->seqlen-1 && this_assemblingRead->orientation==ORIENTATION_PLUS)
				|| (this_assemblingRead->basePos==0 && this_assemblingRead->orientation==ORIENTATION_MINUS))
			{
				if((double)this_assemblingRead->matchBaseNum/this_assemblingRead->seqlen>=readSimilarityThres)
				{


#if (DEBUG_PARA_PRINT==YES)
					// ################################ Debug information ##############################
					if(this_assemblingRead->alignSuccessTimes>0)
					{
						//printf("rid=%ld, aligned success\n", (int64_t)this_assemblingRead->rid);
						totalAlignedSuccessReadNum ++;
					}

					if(this_assemblingRead->unmatchBaseNum<11)
					{
						totalErrReadNum ++;
						errNumArr[this_assemblingRead->unmatchBaseNum] ++;
					}
					// ################################ Debug information ##############################
#endif

//					//if(this_assemblingRead->unmatchBaseNum!=this_assemblingRead->errBaseNum)
//					{
//						// adjust the match information of a read
//						if(adjustMatchInfoRead(this_assemblingRead, contigArr, itemNumContigArr)==FAILED)
//						{
//							printf("line=%d, In %s(), cannot adjust the match information of a read, error!\n", __LINE__, __func__);
//							return FAILED;
//						}
//
//						this_assemblingRead->alignNum ++;
//						this_assemblingRead->alignSuccessTimesOld = this_assemblingRead->alignSuccessTimes;
//					}

					// correct errors
					//if(correctionAllowed==YES && (this_assemblingRead->matchBaseNum!=this_assemblingRead->seqlen || this_assemblingRead->unmatchBaseNum>0))
					//if(correctionAllowed==YES && (this_assemblingRead->unmatchBaseNum>0 && this_assemblingRead->matchBaseNum>=this_assemblingRead->seqlen-2))
					if(correctionAllowed==YES && (this_assemblingRead->unmatchBaseNum>0 && this_assemblingRead->matchBaseNum>=this_assemblingRead->seqlen-maxErrBaseNumInCorrection))
					{
						if(correctRead(this_assemblingRead, contigArr, deBruijnGraph, fpReadCorrected)==FAILED)
						{
							printf("line=%d, In %s(), cannot correct the read bases, error!\n", __LINE__, __func__);
							return FAILED;
						}

						totalReadsNumCorreted ++;
					}

					this_assemblingRead->status = SUCCESSFUL_STATUS;
					updateSameReadStatusToFailed(this_assemblingRead->rid, i, dtRowHashtable);
				}
//				//else if((double)this_assemblingRead->unmatchBaseNum<maxUnmatchBaseNumAfterAlign && this_assemblingRead->alignNum>0)
//				else if((double)this_assemblingRead->unmatchBaseNum>maxUnmatchBaseNumAfterAlign && this_assemblingRead->alignNum>0)
//				//else if((double)this_assemblingRead->unmatchBaseNum>10)
//				{
//					this_assemblingRead->status = FAILED_STATUS;
//				}
				else
				{
					// To do ...
					// further check the similarity by alignment

					this_assemblingRead->status = FAILED_STATUS;
				}
			}
			else if(this_assemblingRead->unmatchBaseNum>maxErrBaseNumInCorrection && this_assemblingRead->alignNum>0)
			//else if(this_assemblingRead->unmatchBaseNum>6 && this_assemblingRead->alignNum>0)
			{
				this_assemblingRead->status = FAILED_STATUS;
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
