/*
 * hashPE.c
 *
 *  Created on: Dec 1, 2011
 *      Author: xiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Estimate insertSize and Sdev.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short estimateInsertSizeAndSdev(char *graphFileName)
{
	printf("Begin estimating the insert size and standard deviation of paired end fragment library ...\n");

	if(PEGivenType>=INSERT_PE_GIVEN_TYPE)
		minContigLenEst = meanSizeInsert * MIN_CONTIG_LEN_EST_FACTOR;
	else
		minContigLenEst = MIN_CONTIG_LEN_EST;

	if(initPEHashParas()==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the PEHash table parameters, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// build estimated contigs
	estimateSuccessFlag = buildEstContigs(sampleContigsFile);
	if(estimateSuccessFlag!=FAILED)
	{
		//printf("line=%d, In %s(), can not build estimated contigs. Error!\n", __LINE__, __func__);
		//return FAILED;
		// estimation
		if(meanSizeInsertAndSdevEstimation(fragmentSizeFile, graphFileName)==FAILED)
		{
			printf("line=%d, In %s(), cannot estimate the insert size and standard deviation of library fragments, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{

	}

	// reset k-mer hash table
	if(resetGraph(deBruijnGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot reset the De Bruijn graph, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ############################ Debug information ##############################
	//if(checkGraph(deBruijnGraph)==FAILED)
	//{
	//	printf("line=%d, In %s(), checking graph error!\n", __LINE__, __func__);
	//	return FAILED;
	//}
	// ############################ Debug information ##############################

	printf("End estimating the insert size and standard deviation of paired end fragment library.\n");

	return SUCCESSFUL;
}

/**
 * Initialize the PEHash table parameters.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initPEHashParas()
{
	int32_t minContigLenUsingPEThresholdTmp;

	if(PEGivenType==INSERT_PE_GIVEN_TYPE)
	{
		standardDev = meanSizeInsert * DRAFT_SDEV_FACTOR;
	}

	if(PEGivenType>NONE_PE_GIVEN_TYPE)
	{
		if(standardDev>=0.1*meanSizeInsert)
			standardDevFactor = SDEV_FACTOR;
		else
			standardDevFactor = SDEV_FACTOR + 1;

		minContigLenUsingPEThresholdTmp = 0.5 * readLen;
		if(minContigLenUsingPEThresholdTmp<kmerSize)
			minContigLenUsingPEThresholdTmp = kmerSize;

		minContigLenUsingPE = minMarginLenPEHash = (int32_t)(meanSizeInsert - standardDevFactor * standardDev - (readLen-kmerSize+1) + 1);
		//if(minContigLenUsingPE<KMER_SIZE)
//		if(minContigLenUsingPE<readLen)  // added 2014-01-15, deleted 2014-12-19
		if(minContigLenUsingPE<minContigLenUsingPEThresholdTmp)  // added 2014-12-19
		{
			//minContigLenUsingPE = minMarginLenPEHash = KMER_SIZE;
			//minContigLenUsingPE = minMarginLenPEHash = readLen + kmerSize; // deleted 2014-01-15
			//minContigLenUsingPE = minMarginLenPEHash = readLen; // added 2014-01-15, deleted 2014-12-19
			minContigLenUsingPE = minMarginLenPEHash = minContigLenUsingPEThresholdTmp; // added 2014-12-19
		}
		maxMarginLenPEHash = (int32_t)(meanSizeInsert + standardDevFactor * standardDev - (readLen-kmerSize+1));
		//maxMarginLenPEHash = (int32_t)(meanSizeInsert + standardDevFactor * standardDev);
		if(maxMarginLenPEHash<=minMarginLenPEHash)
		{
			printf("line=%d, In %s(), maxMarginLenPEHash=%d <= minMarginLenPEHash=%d, incorrect given insert size, error!\n", __LINE__, __func__, maxMarginLenPEHash, minMarginLenPEHash);
			return FAILED;
		}
		maxRegLenPEHash = maxMarginLenPEHash - minMarginLenPEHash + 1;
		if(maxRegLenPEHash<MAX_REG_LEN_PE_HASH_THRES)  // added 2014-01-15
		{
			minMarginLenPEHash -= (MAX_REG_LEN_PE_HASH_THRES - maxRegLenPEHash - 1) / 2 + 1;
			minContigLenUsingPE = minMarginLenPEHash;
			maxMarginLenPEHash += (MAX_REG_LEN_PE_HASH_THRES - maxRegLenPEHash - 1) / 2 + 1;
			maxRegLenPEHash = maxMarginLenPEHash - minMarginLenPEHash + 1;
		}

//		if(minMarginLenPEHash<readLen)
//			shiftLenRound1 = minMarginLenPEHash;
//		else
			shiftLenRound1 = readLen;

		minRegLenUsingPE = maxRegLenPEHash * REG_LEN_PE_HASH_FACTOR;

		minContigLenCheckGap = maxRegLenPEHash + 1.5 * readLen;

		if(minContigLenUsingPE<readLen)
			shortInsertFlag = YES;

#if (DEBUG_PARA_PRINT==YES)
		printf("maxMarginLenPEHash=%d\n", maxMarginLenPEHash);
		printf("minMarginLenPEHash=%d\n", minMarginLenPEHash);
		printf("minContigLenUsingPE=%d\n", minContigLenUsingPE);
		printf("maxRegLenPEHash=%d\n", maxRegLenPEHash);
		printf("minRegLenUsingPE=%d\n", minRegLenUsingPE);
		printf("minContigLenCheckGap=%d\n", minContigLenCheckGap);
#endif
	}

	// ######################### Debug information ##########################
//	if(PEGivenType==INSERT_PE_GIVEN_TYPE)
//	{
//		printf("ReadLen=%d, meanSizeInsert=%.2f\n", readLen, meanSizeInsert);
//		printf("The standard deviation is set to be %.2f by default.\n", standardDev);
//	}else if(PEGivenType==BOTH_PE_GIVEN_TYPE)
//		printf("readLen=%d, meanSizeInsert=%.2f, standardDev=%.2f\n", readLen, meanSizeInsert, standardDev);
	// ######################### Debug information ##########################

	return SUCCESSFUL;
}

/**
 * Initialize the PE hash table before second round assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initPEHashtableSecondAssembly(contigtype *contigArray, int32_t contigNodesNum, int32_t cleanFlag)
{
	int32_t i, j, tmpLeftContigIndex, tmpRightContigIndex;
	//contigtype *tmpContig;
	int32_t ridposnum;
	successRead_t *tmpRidposorient;

	if(cleanFlag==YES)
	{
		if(cleanReadsFromPEHashtable()==FAILED)
		{
			printf("line=%d, In %s, cannot clean PE hash table, error!\n", __LINE__, __func__);
			return ERROR;
		}

		//###################### Debug information ########################
		if(readsNumInPEHashArr!=0)
		{
			printf("line=%d, In %s(), readsNumInPEHashArr=%d != 0, error!\n", __LINE__, __func__, readsNumInPEHashArr);
			return FAILED;
		}
		//###################### Debug information ########################
	}

	validReadOrientPEHash = -1;
	if(contigNodesNum>=maxMarginLenPEHash)
	{ // full length of the hash region
		tmpRightContigIndex = contigNodesNum - minMarginLenPEHash + 1;
		tmpLeftContigIndex = contigNodesNum - maxMarginLenPEHash + 1;

		regLenPEHash = tmpRightContigIndex - tmpLeftContigIndex + 1;
		if(regLenPEHash!=maxRegLenPEHash)
		{
			printf("line=%d, In %s(), regLenPEHash=%d, != maxRegLenPEHash=%d, error!\n", __LINE__, __func__, regLenPEHash, maxRegLenPEHash);
			return FAILED;
		}

		//hashRegLeftContig = contigArray + tmpLeftContigIndex - 1;
		//hashRegRightContig = contigArray + tmpRightContigIndex - 1;
		leftContigRowHashReg = tmpLeftContigIndex - 1;
		rightContigRowHashReg = tmpRightContigIndex - 1;
		validReadOrientPEHash = ORIENTATION_PLUS;

		// add reads into PE hash table
		//tmpContig = contigArray + leftContigRowHashReg;
		//while(tmpContig)
		for(j=leftContigRowHashReg; j<=rightContigRowHashReg; j++)
		{
			ridposnum = contigArray[j].ridposnum;
			tmpRidposorient = contigArray[j].pridposorientation;
			for(i=0; i<ridposnum; i++)
			{
				if(tmpRidposorient[i].orientation==validReadOrientPEHash && tmpRidposorient[i].matchnum==tmpRidposorient[i].seqlen)
				{
					if(addReadToPEHashtable(tmpRidposorient+i, contigArray[j].index, SECOND_ROUND_ASSEMBLY)==FAILED)
					{
						printf("In %s(), cannot add read %lu to PE hash table, error!\n", __func__, (uint64_t)tmpRidposorient[i].rid);
						return FAILED;
					}
				}
			}

			//if(tmpContig==hashRegRightContig)
			//	break;

			//tmpContig ++;
		}

		//########################## Debug information ###########################
//		if(tmpContig==NULL || tmpContig!=hashRegRightContig)
//		{
//			printf("line=%d, In %s(), PE hash region error!\n", __LINE__, __func__);
//			return FAILED;
//		}
		//########################## Debug information ###########################

	}else if(contigNodesNum>=minMarginLenPEHash)
	{ // some reads in the hash region
		tmpRightContigIndex = contigNodesNum - minMarginLenPEHash + 1;
		tmpLeftContigIndex = 1;

		regLenPEHash = tmpRightContigIndex - tmpLeftContigIndex + 1;
		if(regLenPEHash<=0)
		{
			printf("line=%d, In %s(), regLenPEHash=%d, error!\n", __LINE__, __func__, regLenPEHash);
			return FAILED;
		}

		//hashRegLeftContig = contigArray;
		//hashRegRightContig = contigArray + tmpRightContigIndex - 1;
		leftContigRowHashReg = 0;
		rightContigRowHashReg = tmpRightContigIndex - 1;
		validReadOrientPEHash = ORIENTATION_PLUS;

		// add reads into PE hash table
		//tmpContig = hashRegLeftContig;
		//while(tmpContig)
		for(j=leftContigRowHashReg; j<=rightContigRowHashReg; j++)
		{
			ridposnum = contigArray[j].ridposnum;
			tmpRidposorient = contigArray[j].pridposorientation;
			for(i=0; i<ridposnum; i++)
			{
				if(tmpRidposorient[i].orientation==validReadOrientPEHash && tmpRidposorient[i].matchnum==tmpRidposorient[i].seqlen)
				{
					if(addReadToPEHashtable(tmpRidposorient+i, contigArray[j].index, SECOND_ROUND_ASSEMBLY)==FAILED)
					{
						printf("In %s(), cannot add read %lu to PE hash table, error!\n", __func__, (uint64_t)tmpRidposorient[i].rid);
						return FAILED;
					}
				}
			}

			//if(tmpContig==hashRegRightContig)
			//	break;

			//tmpContig ++;
		}

		//########################## Debug information ###########################
//		if(tmpContig==NULL || tmpContig!=hashRegRightContig)
//		{
//			printf("line=%d, In %s(), hash region error!\n", __LINE__, __func__);
//			return FAILED;
//		}
		//########################## Debug information ###########################

	}else
	{ // do nothing

	}

	return SUCCESSFUL;
}

/**
 * Update the PE hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updatePEHashTable(int32_t contigNodesNum, int32_t assemblyRound)
{
	int32_t i, newContigIndex, ridposnum;
	successRead_t *tmpRidposOrient;
	int32_t tmpContigRowLeftReg, tmpContigRowRightReg;

	if(contigNodesNum>=maxMarginLenPEHash)
	{ // slide the window by one base pair
		if(assemblyRound==FIRST_ROUND_ASSEMBLY)
		{ // the first round
			tmpContigRowLeftReg = leftContigRowShiftedReg;
			tmpContigRowRightReg = rightContigRowShiftedReg;
		}else
		{ // the second round
			tmpContigRowLeftReg = leftContigRowHashReg;
			tmpContigRowRightReg = rightContigRowHashReg;
		}

		// remove the reads hanging on the left margin
		ridposnum = contigArr[tmpContigRowLeftReg].ridposnum;
		if(ridposnum>0)
		{
			tmpRidposOrient = contigArr[tmpContigRowLeftReg].pridposorientation;
			for(i=0; i<ridposnum; i++)
			{
				if(tmpRidposOrient[i].orientation==validReadOrientPEHash && tmpRidposOrient[i].matchnum==tmpRidposOrient[i].seqlen)
				{
					if(delReadfromPEHashtable(tmpRidposOrient[i].rid)==FAILED)
					{
						printf("line=%d, In %s(), cannot remove read %lu from PE hash table, error!\n", __LINE__, __func__, (uint64_t)tmpRidposOrient[i].rid);
						return FAILED;
					}
				}
			}
		}

		// add new reads
		if(contigArr[tmpContigRowRightReg].index==contigNodesNum)
		{
			printf("line=%d, In %s(), pContig->next==NULL, error!\n", __LINE__, __func__);
			return FAILED;
		}
		tmpContigRowRightReg ++;

		ridposnum = contigArr[tmpContigRowRightReg].ridposnum;
		if(ridposnum>0)
		{
			tmpRidposOrient = contigArr[tmpContigRowRightReg].pridposorientation;
			for(i=0; i<ridposnum; i++)
			{
				if(tmpRidposOrient[i].orientation==validReadOrientPEHash && tmpRidposOrient[i].matchnum==tmpRidposOrient[i].seqlen)
				{
					//printf("line=%d, In %s(), assemblyRound=%d, rid=%ld, hangingIndex=%d, leftContigRowHashReg=%d, rightContigRowHashReg=%d\n", __LINE__, __func__, assemblyRound, (int64_t)tmpRidposOrient[i].rid, tmpRidposOrient[i].hangingIndex, tmpContigRowLeftReg, tmpContigRowRightReg);

					if(addReadToPEHashtable(tmpRidposOrient+i, contigArr[tmpContigRowRightReg].index, assemblyRound)==FAILED)
					{
						printf("line=%d, In %s(), cannot add read %lu to PE hash table, error!\n", __LINE__, __func__, (uint64_t)tmpRidposOrient[i].rid);
						return FAILED;
					}
				}
			}
		}

		if(assemblyRound==FIRST_ROUND_ASSEMBLY)
		{ // the first round
			leftContigRowShiftedReg ++;
			rightContigRowShiftedReg ++;
		}else
		{ // the second round
			leftContigRowHashReg ++;
			rightContigRowHashReg ++;
		}
	}else if(contigNodesNum>=minMarginLenPEHash)
	{
		//if(contigNodesNum==minMarginLenPEHash)
		if(validReadOrientPEHash==-1)
		{ // initialize the PE hash table
//			if(readsNumInPEHashArr>0)
//			{
//				if(cleanReadsFromPEHashtable()==FAILED)
//				{
//					printf("line=%d, In %s(), cannot clean PE hash table, error!\n", __LINE__, __func__);
//					return FAILED;
//				}
//			}

			readsNumInPEHashArr = 0;
			regLenPEHash = 1;
			if(assemblyRound==FIRST_ROUND_ASSEMBLY)
			{ // the first round
				newContigIndex = shiftLenRound1;
				if(newContigIndex>contigNodesNum)
					newContigIndex = contigNodesNum;
				leftContigRowShiftedReg = rightContigRowShiftedReg = newContigIndex - 1;
				validReadOrientPEHash = ORIENTATION_MINUS;

				// add reads to PE hash table
				ridposnum = contigArr[rightContigRowShiftedReg].ridposnum;
				if(ridposnum>0)
				{
					tmpRidposOrient = contigArr[rightContigRowShiftedReg].pridposorientation;
					for(i=0; i<ridposnum; i++)
					{
						if(tmpRidposOrient[i].orientation==validReadOrientPEHash && tmpRidposOrient[i].matchnum==tmpRidposOrient[i].seqlen)
						{
							//printf("line=%d, In %s(), assemblyRound=%d, rid=%ld, hangingIndex=%d, leftContigRowHashReg=%d, rightContigRowHashReg=%d\n", __LINE__, __func__, assemblyRound, (int64_t)tmpRidposOrient[i].rid, tmpRidposOrient[i].hangingIndex, leftContigRowShiftedReg, rightContigRowShiftedReg);

							if(addReadToPEHashtable(tmpRidposOrient+i, contigArr[rightContigRowShiftedReg].index, assemblyRound)==FAILED)
							{
								printf("line=%d, In %s(), cannot add read %lu to PE hash table, error!\n", __LINE__, __func__, (uint64_t)tmpRidposOrient[i].rid);
								return FAILED;
							}
						}
					}
				}
			}else
			{ // the second round
				leftContigRowHashReg = rightContigRowHashReg = 0;
				validReadOrientPEHash = ORIENTATION_PLUS;

				// add reads to PE hash table
				ridposnum = contigArr[rightContigRowHashReg].ridposnum;
				if(ridposnum>0)
				{
					tmpRidposOrient = contigArr[rightContigRowHashReg].pridposorientation;
					for(i=0; i<ridposnum; i++)
					{
						if(tmpRidposOrient[i].orientation==validReadOrientPEHash && tmpRidposOrient[i].matchnum==tmpRidposOrient[i].seqlen)
						{
							//printf("line=%d, In %s(), assemblyRound=%d, rid=%ld, hangingIndex=%d, leftContigRowHashReg=%d, rightContigRowHashReg=%d\n", __LINE__, __func__, assemblyRound, (int64_t)tmpRidposOrient[i].rid, tmpRidposOrient[i].hangingIndex, leftContigRowHashReg, rightContigRowHashReg);

							if(addReadToPEHashtable(tmpRidposOrient+i, contigArr[rightContigRowHashReg].index, assemblyRound)==FAILED)
							{
								printf("line=%d, In %s(), cannot add read %lu to PE hash table, error!\n", __LINE__, __func__, (uint64_t)tmpRidposOrient[i].rid);
								return FAILED;
							}
						}
					}
				}
			}
		}else
		{ // enlarge the window by one base pair
			if(assemblyRound==FIRST_ROUND_ASSEMBLY)
			{ // the first round
				tmpContigRowLeftReg = leftContigRowShiftedReg;
				tmpContigRowRightReg = rightContigRowShiftedReg;
			}else
			{ // the second round
				tmpContigRowLeftReg = leftContigRowHashReg;
				tmpContigRowRightReg = rightContigRowHashReg;
			}

			regLenPEHash ++;

			// add new reads
			if(contigArr[tmpContigRowRightReg].index==contigNodesNum)
			{
				printf("line=%d, In %s(), pContig->next==NULL, error!\n", __LINE__, __func__);
				return FAILED;
			}
			tmpContigRowRightReg ++;

			ridposnum = contigArr[tmpContigRowRightReg].ridposnum;
			if(ridposnum>0)
			{
				tmpRidposOrient = contigArr[tmpContigRowRightReg].pridposorientation;
				for(i=0; i<ridposnum; i++)
				{
					if(tmpRidposOrient[i].orientation==validReadOrientPEHash && tmpRidposOrient[i].matchnum==tmpRidposOrient[i].seqlen)
					{
						//printf("line=%d, In %s(), assemblyRound=%d, rid=%ld, hangingIndex=%d, leftContigRowHashReg=%d, rightContigRowHashReg=%d\n", __LINE__, __func__, assemblyRound, (int64_t)tmpRidposOrient[i].rid, tmpRidposOrient[i].hangingIndex, tmpContigRowLeftReg, tmpContigRowRightReg);

						if(addReadToPEHashtable(tmpRidposOrient+i, contigArr[tmpContigRowRightReg].index, assemblyRound)==FAILED)
						{
							printf("line=%d, In %s(), cannot add read %lu to PE hash table, error!\n", __LINE__, __func__, (uint64_t)tmpRidposOrient[i].rid);
							return FAILED;
						}
					}
				}
			}

			if(assemblyRound==FIRST_ROUND_ASSEMBLY)
			{ // the first round
				rightContigRowShiftedReg ++;
			}else
			{ // the second round
				rightContigRowHashReg ++;
			}
		}
	}else
	{ // do nothing

	}

	return SUCCESSFUL;
}

/**
 * Get the paired read by the given readID.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadFromPEHashtable(PERead_t **pRead, uint64_t readID)
{
	uint64_t hashcode;

	hashcode = readID & RID_LOW_BITS_MASK;
	if(PEHashArr[hashcode])
	{
		*pRead = PEHashArr[hashcode];
		while(*pRead)
		{
			if((*pRead)->rid==readID)
			{
				break;
			}
			*pRead = (*pRead)->next;
		}
	}else
	{
		*pRead = NULL;
	}

	return SUCCESSFUL;
}

/**
 * Add a read into PE hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addReadToPEHashtable(successRead_t *ridposOrient, int32_t contigPos, int32_t assemblyRound)
{
	uint64_t hashcode;
	PERead_t *tmpRead;

	// ######################## Debug information ############################
	//if(ridposOrient->rid==134048)
	//{
	//	printf("line=%d, In %s(), rid=%ld, matchnum=%d, orient=%c, firstpos=%d\n", __LINE__, __func__, (int64_t)ridposOrient->rid, ridposOrient->matchnum, ridposOrient->orientation, ridposOrient->pos);
	//}
	// ######################## Debug information ############################

	hashcode = ridposOrient->rid & RID_LOW_BITS_MASK;

	tmpRead = (PERead_t*) malloc(sizeof(PERead_t));
	if(!tmpRead)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{ // the first round
		tmpRead->cpos = contigPos - ridposOrient->seqlen + 1;
		if(tmpRead->cpos<=0)
		{
			printf("line=%d, In %s(), cpos=%d, error!\n", __LINE__, __func__, tmpRead->cpos);
			return FAILED;
		}
	}else
	{ // the second round
		tmpRead->cpos = contigPos;
	}

	tmpRead->rid = ridposOrient->rid;
	tmpRead->orient = ORIENTATION_PLUS;
	tmpRead->seqlen = ridposOrient->seqlen;
	//insert the new node from head
	tmpRead->next = PEHashArr[hashcode];
	PEHashArr[hashcode] = tmpRead;

	readsNumInPEHashArr ++;

	return SUCCESSFUL;
}

/**
 * Remove a read into PE hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short delReadfromPEHashtable(uint64_t readID)
{
	uint64_t hashcode;
	PERead_t *tmpRead, *preRead;

	// ######################## Debug information ############################
//	if(readID==16798181)
//	{
//		printf("line=%d, In %s(), readID=%lu\n", __LINE__, __func__, readID);
//	}
	// ######################## Debug information ############################

	hashcode = readID & RID_LOW_BITS_MASK;
	tmpRead = PEHashArr[hashcode];
	if(!tmpRead)
	{
		printf("line=%d, In %s(), cannot get the read %lu from PE hash table, error!\n", __LINE__, __func__, readID);
		return FAILED;
	}

	preRead = NULL;
	while(tmpRead)
	{
		if(tmpRead->rid==readID)
			break;

		preRead = tmpRead;
		tmpRead = tmpRead->next;
	}

	if(!tmpRead)
	{
		printf("line=%d, In %s(), cannot get read %lu from PE hash table, error!\n", __LINE__, __func__, readID);
		return FAILED;
	}

	if(preRead==NULL)
	{
		PEHashArr[hashcode] = tmpRead->next;
	}else
	{
		preRead->next = tmpRead->next;
	}
	free(tmpRead);

	readsNumInPEHashArr --;

	return SUCCESSFUL;
}

/**
 * Clean the PE hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short cleanReadsFromPEHashtable()
{
	int i;
	PERead_t *tmpRead, *head;

	//###################### Debug information ########################
	//printf("line=%d, In %s(), before clean, readsNumInPEHashArr=%d\n", __LINE__, __func__, readsNumInPEHashArr);
	//###################### Debug information ########################

	for(i=0; i<TABLE_SIZE_HASH_PE; i++)
	{
		if(PEHashArr[i])
		{
			head = PEHashArr[i];
			while(head)
			{
				tmpRead = head->next;
				free(head);
				head = tmpRead;
				readsNumInPEHashArr --;
			}
			PEHashArr[i] = NULL;
		}
	}

	//###################### Debug information ########################
	//printf("line=%d, In %s(), after clean, readsNumInPEHashArr=%d\n", __LINE__, __func__, readsNumInPEHashArr);
	//###################### Debug information ########################

	regLenPEHash = 0;
	readsNumInPEHashArr = 0;

	// ######################### Debug information ##############################
//	if(readsNumInPEHashArr!=0)
//	{
//		printf("line=%d, In %s(), readsNumInPEHashArr=%d != 0, error!\n", __LINE__, __func__, readsNumInPEHashArr);
//		return FAILED;
//	}
	// ######################### Debug information ##############################

	return SUCCESSFUL;
}

// ============================== Estimation ================================
/**
 * Estimate the insert size and standard deviation of library fragments.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short meanSizeInsertAndSdevEstimation(const char *fragmentSizeFile, const char *graphFileName)
{

/*
	// ################################# Debug information ###############################
	if(convertFragmentSizeFileInText(fragmentSizeFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot convert fragment size file in text mode, error!\n", __LINE__, __func__);
		return FAILED;
	}
	// ################################# Debug information ###############################
*/

	// reload the file and estimate the insert size and standard deviation
	if(computeInsertSizeAndSdev(&meanSizeInsertEst, &standardDevEst, fragmentSizeFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the insert size and standard deviation, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(PEGivenType>NONE_PE_GIVEN_TYPE)
	{
		oldMeanSizeInsert = meanSizeInsert;
		oldStandardDev = standardDev;
	}

	meanSizeInsert = meanSizeInsertEst;
	standardDev = standardDevEst;

	oldPEGivenType = PEGivenType;
	PEGivenType = BOTH_PE_GIVEN_TYPE;

	// update the insert size and sdev in k-mer hash table
	if(updateInsertSizeAndSdevInGraph(PEGivenType, meanSizeInsert, standardDev, graphFileName)==FAILED)
	{
		printf("line=%d, In %s(), cannot update the insert size and standard deviation in k-mer hash table, error!\n", __LINE__, __func__);
		return FAILED;
	}

#if	(RESERVE_EST_FILES==YES)
	remove(fragmentSizeFile);
#endif

	return SUCCESSFUL;
}

/**
 * Get the paired ends from single contig, and output them to a binary file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getPairedEndsFromSingleContig(FILE *fpFragSize, contigtype *contigArr, int64_t itemNumContigArr)
{
	// initialize the memory
	if(initMemGetPESingleContig(contigArr, itemNumContigArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize memory for getting paired ends of a contig, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill data
	if(fillDataReadPosTmpArr(readPosTmpArr, contigArr, itemNumContigArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill data of readPosTmp array of a contig, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// sort the data in readPosTmp array
	if(radixSortReadPosTmpArr(readPosTmpArr, readPosTmpArrBuf, readsNumSingleContig)==FAILED)
	{
		printf("line=%d, In %s(), cannot sort the reads in readPosTmp Array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Fill data of the readPosTmp array to read list (RL)
	if(fillDataReadList(readListArr, readPosArr, &itemNumInReadListArr, &itemNumInReadPosArr, readPosTmpArr, readsNumSingleContig)==FAILED)
	{
		printf("line=%d, In %s(), cannot convert the readPosTmp array to read list (RL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ############################### Debug information ###############################
	//if(outputMatedReadsInReadListToFile("../matedReads", readListArr, readPosArr, itemNumInReadListArr)==FAILED)
	//{
	//	printf("line=%d, In %s(), cannot mated reads in read list (RL) to file, error!\n", __LINE__, __func__);
	//	return FAILED;
	//}
	// ############################### Debug information ###############################

	// get the valid insert size, and output their fragment size to binary file
	if(outputInsertSizeToFile(fpFragSize, readListArr, readPosArr, itemNumInReadListArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot output valid insert size to file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// free memory
	freeMemGetPESingleContig();

	return SUCCESSFUL;
}

/**
 * Initialize the memory for getting paired ends of single contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemGetPESingleContig(contigtype *contigArr, int64_t itemNumContigArr)
{
	// get the total reads number in single contig
	if(getTotalReadsNumOfSingleContig(&readsNumSingleContig, contigArr, itemNumContigArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot get total reads number of a contig, error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(readsNumSingleContig<=0)
	{
		printf("line=%d, In %s(), readsNumSingleContig=%ld, error!\n", __LINE__, __func__, readsNumSingleContig);
		return FAILED;
	}

	// allocate memory
	readPosTmpArr = (readPosTemp_t*) malloc(readsNumSingleContig * sizeof(readPosTemp_t));
	if(readPosTmpArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	readPosTmpArrBuf = (readPosTemp_t*) malloc(readsNumSingleContig * sizeof(readPosTemp_t));
	if(readPosTmpArrBuf==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	readListArr = (readList_t *) calloc(readsNumSingleContig, sizeof(readList_t));
	if(readListArr==NULL)
	{
		printf("line=%d, In %s(), cannot callocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	readPosArr = (readPos_t *) calloc(readsNumSingleContig, sizeof(readPos_t));
	if(readPosArr==NULL)
	{
		printf("line=%d, In %s(), cannot callocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Free memory of getting paired ends of single contig.
 */
void freeMemGetPESingleContig()
{
	readsNumSingleContig = 0;

	free(readPosTmpArr);
	readPosTmpArr = NULL;
	free(readPosTmpArrBuf);
	readPosTmpArrBuf = NULL;

	free(readListArr);
	readListArr = NULL;
	free(readPosArr);
	readPosArr = NULL;
}

/**
 * Get total reads number of a contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getTotalReadsNumOfSingleContig(int64_t *totalReadsNum, contigtype *contigArr, int64_t contigNodesNum)
{
	int64_t i;

	*totalReadsNum = 0;
	for(i=0; i<contigNodesNum; i++)
		if(contigArr[i].ridposnum>0)
			*totalReadsNum += contigArr[i].ridposnum;

	return SUCCESSFUL;
}

/**
 * Fill the data of readPosTmp Array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillDataReadPosTmpArr(readPosTemp_t *readPosTmpArray, contigtype *contigArr, int64_t contigNodesNum)
{
	int64_t tmpReadsNum;
	int i, j, ridposnum;
	successRead_t *pReadsArray;

	tmpReadsNum = 0;
	for(i=0; i<contigNodesNum; i++)
	{
		if(contigArr[i].ridposnum>0)
		{
			ridposnum = contigArr[i].ridposnum;
			pReadsArray = contigArr[i].pridposorientation;
			for(j=0; j<ridposnum; j++)
			{
				readPosTmpArray[tmpReadsNum].readID = pReadsArray[j].rid;
				readPosTmpArray[tmpReadsNum].contigPos = i + 1;
				readPosTmpArray[tmpReadsNum].orientation = pReadsArray[j].orientation;
				readPosTmpArray[tmpReadsNum].matchBaseNum = pReadsArray[j].matchnum;
				readPosTmpArray[tmpReadsNum].seqlen = pReadsArray[j].seqlen;
				tmpReadsNum ++;

				// ################### Debug information ##################
				if(tmpReadsNum>readsNumSingleContig)
				{
					printf("line=%d, In %s(), tmpReadsNum=%ld > readsNumSingleContig=%ld, error!\n", __LINE__, __func__, tmpReadsNum, readsNumSingleContig);
					return FAILED;
				}
				// ################### Debug information ##################
			}
		}
	}

	// ################### Debug information ##################
	if(tmpReadsNum!=readsNumSingleContig)
	{
		printf("line=%d, In %s(), tmpReadsNum=%ld != readsNumSingleContig=%ld, error!\n", __LINE__, __func__, tmpReadsNum, readsNumSingleContig);
		return FAILED;
	}
	// ################### Debug information ##################

	return SUCCESSFUL;
}

/**
 * Sort the reads in readPosTmp array to ascending order according to their readID.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short radixSortReadPosTmpArr(readPosTemp_t *readPosTmpArray, readPosTemp_t *readPosTmpArrayBuf, int64_t readsNumInArray)
{
	struct partNode
	{
		uint32_t curItemNum;
		uint32_t totalItemNum;
		uint64_t firstRow;
	};

	int64_t i, step, total;
	readPosTemp_t *data, *buf;
	struct partNode *part;
	int partArrSize, stepBits, maxStepLen;
	unsigned int bitMask;
	unsigned int hashcode, firstRow, curItemNum;

	stepBits = 16;
	maxStepLen = 64;
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
		if((step/stepBits)%2==1)
		{
			buf = readPosTmpArray;
			data = readPosTmpArrayBuf;
		}else
		{
			data = readPosTmpArray;
			buf = readPosTmpArrayBuf;
		}

		// count the number
		if(memset(part, 0, partArrSize * sizeof(struct partNode))==NULL)
		{
			printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
			free(part);
			return FAILED;
		}
		for(i=0; i<readsNumInArray; i++)
			part[ (data[i].readID >> step) & bitMask ].totalItemNum ++;

		// initialize the part array
		total = 0;
		for(i=0; i<partArrSize; i++)
		{
			part[i].firstRow = total;
			total += part[i].totalItemNum;
		}

		// copy the data to the right place
		for(i=0; i<readsNumInArray; i++)
		{
			hashcode = (data[i].readID >> step) & bitMask;
			firstRow = part[hashcode].firstRow;
			curItemNum = part[hashcode].curItemNum;
			if(memcpy(buf+firstRow+curItemNum, data+i, sizeof(readPosTemp_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				free(part);
				return FAILED;
			}
			part[hashcode].curItemNum ++;
		}

		step += stepBits;

		//######################## Debug information #######################
		for(i=0; i<partArrSize; i++)
		{
			if(part[i].curItemNum!=part[i].totalItemNum)
			{
				printf("line=%d, In %s(), in part[%ld], curItemNum=%u != totalItemNum=%u, error!\n", __LINE__, __func__, i, part[i].curItemNum, part[i].totalItemNum);
				free(part);
				return FAILED;
			}
		}
		//######################## Debug information #######################
	}

	free(part);

	return SUCCESSFUL;
}

/**
 * Fill data of the readPosTmp array to read list (RL).
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillDataReadList(readList_t *readListArray, readPos_t *readPosArray, int64_t *itemNumInReadListArray, int64_t *itemNumInReadPosArray, readPosTemp_t *readPosTmpArray, int64_t itemNumInReadPosTmpArray)
{
	int64_t i, j;
	int sameItemNum;

	*itemNumInReadListArray = *itemNumInReadPosArray = 0;

	i = 0;
	while(i<itemNumInReadPosTmpArray)
	{
		sameItemNum = 1;
		for(j=i+1; j<itemNumInReadPosTmpArray; j++)
		{
			if(readPosTmpArray[i].readID!=readPosTmpArray[j].readID)
				break;
			sameItemNum ++;
		}

		// fill data of the item
		readListArray[*itemNumInReadListArray].readID = readPosTmpArray[i].readID;
		readListArray[*itemNumInReadListArray].seqlen = readPosTmpArray[i].seqlen;
		readListArray[*itemNumInReadListArray].matchNum = sameItemNum;
		readListArray[*itemNumInReadListArray].firstRow = *itemNumInReadPosArray;
		(*itemNumInReadListArray) ++;

		for(j=0; j<sameItemNum; j++)
		{
			readPosArray[(*itemNumInReadPosArray)+j].contigPos = readPosTmpArray[i+j].contigPos;
			readPosArray[(*itemNumInReadPosArray)+j].matchBaseNum = readPosTmpArray[i+j].matchBaseNum;
			readPosArray[(*itemNumInReadPosArray)+j].orientation = readPosTmpArray[i+j].orientation;
		}
		(*itemNumInReadPosArray) += sameItemNum;

		i += sameItemNum;
	}

	// ###################### Debug information ########################
	if((*itemNumInReadPosArray)!=itemNumInReadPosTmpArray)
	{
		printf("line=%d, In %s(), itemNumInReadPosArray=%ld != itemNumInReadPosTmpArray=%ld, error!\n", __LINE__, __func__, *itemNumInReadPosArray, itemNumInReadPosTmpArray);
		return FAILED;
	}
	// ###################### Debug information ########################

	// ###################### Debug information ########################
	//printf("itemNumInReadListArray=%ld, itemNumInReadPosArray=%ld\n", *itemNumInReadListArray, *itemNumInReadPosArray);
	// ###################### Debug information ########################

	return SUCCESSFUL;
}

/**
 * Get the valid insert size, and output their fragment size to binary file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputInsertSizeToFile(FILE *fpFragSize, readList_t *readListArray, readPos_t *readPosArray, int64_t itemNumInReadListArray)
{
	int64_t i;
	int32_t fragmentSize;
	readList_t *pLeftRead, *pRightRead;
	readPos_t *pLeftReadPos, *pRightReadPos;

	i = 0;
	while(i<itemNumInReadListArray-1)
	{
		if(readListArray[i].readID % 2 == 1 && readListArray[i].readID+1==readListArray[i+1].readID)
		{ // the odd readID
			//if(readListArray[i].matchNum==1 && readListArray[i+1].matchNum==1)
			if(readListArray[i].matchNum==1 && readListArray[i+1].matchNum==1
				&& readPosArray[readListArray[i].firstRow].matchBaseNum==readListArray[i].seqlen
				&& readPosArray[readListArray[i+1].firstRow].matchBaseNum==readListArray[i+1].seqlen)
			{
				if(readPosArray[ readListArray[i].firstRow ].contigPos < readPosArray[ readListArray[i+1].firstRow ].contigPos)
				{
					pLeftRead = readListArray + i;
					pRightRead = readListArray + i + 1;
					pLeftReadPos = readPosArray + readListArray[i].firstRow;
					pRightReadPos = readPosArray + readListArray[i+1].firstRow;
				}else
				{
					pLeftRead = readListArray + i + 1;
					pRightRead = readListArray + i;
					pLeftReadPos = readPosArray + readListArray[i+1].firstRow;
					pRightReadPos = readPosArray + readListArray[i].firstRow;
				}

				if(pLeftReadPos->orientation==ORIENTATION_PLUS && pRightReadPos->orientation==ORIENTATION_MINUS)
				{
					fragmentSize = pRightReadPos->contigPos + pRightReadPos->matchBaseNum - 1 - pLeftReadPos->contigPos + 1;
					if(PEGivenType>=INSERT_PE_GIVEN_TYPE && fragmentSize<=MAX_INSERT_SIZE_FACTOR*meanSizeInsert)
					{
						if(fwrite(&fragmentSize, sizeof(int32_t), 1, fpFragSize)!=1)
						{
							printf("line=%d, In %s(), (%lu,%u,%u,%u,%c), (%lu,%u,%u,%u,%c), fwrite error!\n", __LINE__, __func__,
								(int64_t)pLeftRead->readID, (int32_t)pLeftRead->seqlen, (int32_t)pLeftReadPos->contigPos, (int32_t)pLeftReadPos->matchBaseNum, (int32_t)pLeftReadPos->orientation,
								(int64_t)pRightRead->readID, (int32_t)pRightRead->seqlen, (int32_t)pRightReadPos->contigPos, (int32_t)pRightReadPos->matchBaseNum, (int32_t)pRightReadPos->orientation);
							return FAILED;
						}
					}else
					{
						if(fwrite(&fragmentSize, sizeof(int32_t), 1, fpFragSize)!=1)
						{
							printf("line=%d, In %s(), (%lu,%u,%u,%u,%c), (%lu,%u,%u,%u,%c), fwrite error!\n", __LINE__, __func__,
								(int64_t)pLeftRead->readID, (int32_t)pLeftRead->seqlen, (int32_t)pLeftReadPos->contigPos, (int32_t)pLeftReadPos->matchBaseNum, (int32_t)pLeftReadPos->orientation,
								(int64_t)pRightRead->readID, (int32_t)pRightRead->seqlen, (int32_t)pRightReadPos->contigPos, (int32_t)pRightReadPos->matchBaseNum, (int32_t)pRightReadPos->orientation);
							return FAILED;
						}
					}
				}else
				{
					//printf("line=%d, In %s(), (%lu,%u,%u,%c), (%lu,%u,%u,%c), invalid paired end!\n", __LINE__, __func__, pLeftRead->readID, pLeftReadPos->contigPos, pLeftReadPos->matchBaseNum, pLeftReadPos->orientation, pRightRead->readID, pRightReadPos->contigPos, pRightReadPos->matchBaseNum, pRightReadPos->orientation);
				}
			}

			i += 2;
		}else
		{
			i ++;
		}
	}

	return SUCCESSFUL;
}

/**
 * Reload the file and estimate the insert size and standard deviation.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeInsertSizeAndSdev(double *meanSizeInsert, double *standardDev, const char *fragmentSizeFile)
{
	int64_t fragNum, fragNum2;
	FILE *fpFragSize;
	int32_t fragmentSize;
	double tmpInsertSize, tmpSdev;
	int32_t longFragNum;

	fpFragSize = fopen(fragmentSizeFile, "rb");
	if(fpFragSize==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fragmentSizeFile);
		return FAILED;
	}

	// estimate the insert size
	longFragNum = 0;
	fragNum = 0;
	tmpInsertSize = 0;
	while(1)
	{
		if(fread(&fragmentSize, sizeof(int32_t), 1, fpFragSize)!=1)
		{
			if(feof(fpFragSize))
				break;
			else
			{
				printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		// ######################## Debug information #####################
		if(fragmentSize<=MAX_INSERT_SIZE_THRES)
		{
			tmpInsertSize += fragmentSize;
			fragNum ++;
		}else
		{
			longFragNum ++;
			//printf("fragmentSize=%d\n", fragmentSize);
		}
		// ######################## Debug information #####################

	}
	*meanSizeInsert = tmpInsertSize / fragNum;

//	if(longFragNum>0)
//		printf("longFragNum=%d\n", longFragNum);

	if(fragNum<=1)
	{
		printf("line=%d, In %s(), fragNum=%ld, error!\n", __LINE__, __func__, fragNum);
		return FAILED;
	}

	// estimate the standard deviation
	fseek(fpFragSize, 0, SEEK_SET);   // seek to the beginning of the file
	tmpSdev = 0;
	fragNum2 = 0;
	while(1)
	{
		if(fread(&fragmentSize, sizeof(int32_t), 1, fpFragSize)!=1)
		{
			if(feof(fpFragSize))
				break;
			else
			{
				printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		if(fragmentSize<=MAX_INSERT_SIZE_THRES)
		{
			tmpSdev += (fragmentSize - *meanSizeInsert) * (fragmentSize - *meanSizeInsert);
			fragNum2 ++;
		}
	}
	*standardDev = sqrt(tmpSdev / (fragNum-1));

	fclose(fpFragSize);
	fpFragSize = NULL;

	remove(fragmentSizeFile);

	// ########################## Debug information ###########################
	if(fragNum!=fragNum2)
	{
		printf("line=%d, In %s(), fragNum=%ld != fragNum2=%ld, error!\n", __LINE__, __func__, fragNum, fragNum2);
		return FAILED;
	}
	// ########################## Debug information ###########################

	printf("The estimated insert size  : %.2f\n", *meanSizeInsert);
	printf("The estimated standard dev : %.2f\n", *standardDev);
	printf("Pairs of used paired ends  : %ld\n", fragNum);
	printf("Number of estimated contigs: %d\n", contigNumEstContigArr);

	return SUCCESSFUL;
}

/**
 * Add success reads to PEHashtable.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addSuccessReadsToPEHashtable(successRead_t *successReadsArray, int32_t itemNumSuccessReadsArray, int32_t assemblyRound)
{
	int32_t i, validReadOrient;
	int32_t tmpLeftContigRowHashReg, tmpRightContigRowHashReg;

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{
		tmpLeftContigRowHashReg = leftContigRowShiftedReg;
		tmpRightContigRowHashReg = rightContigRowShiftedReg;
	}else
	{
		tmpLeftContigRowHashReg = leftContigRowHashReg;
		tmpRightContigRowHashReg = rightContigRowHashReg;
	}
	validReadOrient = ORIENTATION_PLUS;

	if(tmpLeftContigRowHashReg!=tmpRightContigRowHashReg)
	{
		for(i=0; i<itemNumSuccessReadsArray; i++)
		{
			//if(successReadsArray[i].rid==9743836)
			//{
			//	printf("line=%d, In %s(), rid=%ld, hangingIndex=%d\n", __LINE__, __func__, (int64_t)successReadsArray[i].rid, successReadsArray[i].hangingIndex);
			//}

			if(successReadsArray[i].hangingIndex-1>=tmpLeftContigRowHashReg && successReadsArray[i].hangingIndex-1<=tmpRightContigRowHashReg)
			{
				//printf("line=%d, In %s(), assemblyRound=%d, rid=%ld, hangingIndex=%d, leftContigRowHashReg=%d, rightContigRowHashReg=%d\n", __LINE__, __func__, assemblyRound, (int64_t)successReadsArray[i].rid, successReadsArray[i].hangingIndex, tmpLeftContigRowHashReg, tmpRightContigRowHashReg);

				if(successReadsArray[i].orientation==validReadOrient && successReadsArray[i].matchnum==successReadsArray[i].seqlen)
				{
					if(addReadToPEHashtable(successReadsArray+i, successReadsArray[i].hangingIndex, assemblyRound)==FAILED)
					{
						printf("line=%d, In %s(), cannot add read %ld to PE hash table, error!\n", __LINE__, __func__, (int64_t)successReadsArray[i].rid);
						return FAILED;
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}
