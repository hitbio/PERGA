/*
 * tandPath.c
 *
 *  Created on: Nov 20, 2013
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Decide the navigation by checking the tandem paths for paired-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short decideByCheckingTandemRepeatPE(int32_t *naviTandPathPE, int32_t *maxBaseIndexAfterTandPathPE, int32_t *incorrectBaseNumTandPath, int32_t *newCandBaseNumAfterTandPath, int32_t *occNumArray, int32_t *occIndexArray, int32_t contigPos, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, contigPath_t *contigPath, graphtype *graph)
{
	tandemPathItem_t *tandemPathList;
	int32_t itemNumTandemPathList, shortFragSizeRemovedFlag, overlapTandPathRemovedFlag;

	// get the tandem paths
	if(getTandemPathPE(&tandemPathList, &itemNumTandemPathList, contigPos, decisionTable, *readsNumDecisionTable, dtRowHashtable, graph)==FAILED)
	{
		printf("line=%d, In %s(), cannot get tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

#if(TANDPATH_OUTPUT==YES)
	// output the paths
	printf("********* contigPos=%d\n", contigPos);
	printf("********* Before remove short tandPath:\n");
	if(outputTandPath(tandemPathList, itemNumTandemPathList, decisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot output tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}
#endif

	// remove short tandem paths with short fragment size
	if(removeShortFragSizeTandPath(&shortFragSizeRemovedFlag, &tandemPathList, &itemNumTandemPathList, decisionTable, readsNumDecisionTable, dtRowHashtable, contigPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove reads in tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

#if(TANDPATH_OUTPUT==YES)
	// output the paths
	printf("********* After remove short tandPath:\n");
	if(outputTandPath(tandemPathList, itemNumTandemPathList, decisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot output tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}
#endif

	// remove reads in second path of two overlapped tandem paths
	if(removeOverlappedTandPath(&overlapTandPathRemovedFlag, tandemPathList, &itemNumTandemPathList, decisionTable, readsNumDecisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove reads in tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

#if(TANDPATH_OUTPUT==YES)
	// output the paths
	printf("********* After remove overlapped tandPath:\n");
	if(outputTandPath(tandemPathList, itemNumTandemPathList, decisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot output tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}
#endif

//	if(removeShortOverlappedTandPath(&shortFragSizeRemovedFlag, &tandemPathList, &itemNumTandemPathList, decisionTable, dtRowHashtable)==FAILED)  // 2014-01-29
//	{
//		printf("line=%d, In %s(), cannot remove short overlapped tandem paths, error!\n", __LINE__, __func__);
//		return FAILED;
//	}

#if(TANDPATH_OUTPUT==YES)
	// output the paths
	printf("********* After remove short overlapped tandPath:\n");
	if(outputTandPath(tandemPathList, itemNumTandemPathList, decisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot output tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}
#endif

	// decide the navigation from the paths
	if(shortFragSizeRemovedFlag==YES || overlapTandPathRemovedFlag==YES || itemNumTandemPathList==2)
	{
		// decide the navigation
		if(decideNaviByTandPath(naviTandPathPE, maxBaseIndexAfterTandPathPE, incorrectBaseNumTandPath, newCandBaseNumAfterTandPath, occNumArray, occIndexArray, tandemPathList, itemNumTandemPathList)==FAILED)
		{
			printf("line=%d, In %s(), cannot navigate by tandem paths, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		*naviTandPathPE = NAVI_FAILED;
		*maxBaseIndexAfterTandPathPE = INT_MAX;
		*incorrectBaseNumTandPath = INT_MAX;
		*newCandBaseNumAfterTandPath = itemNumTandemPathList;
	}

	// release the memory
	releaseTandemPathList(&tandemPathList);

	return SUCCESSFUL;
}

/**
 * Decide the navigation by checking the tandem paths for single-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short decideByCheckingTandemRepeatSE(int32_t *naviTandPathSE, int32_t *maxBaseIndexAfterTandPathSE, int32_t *incorrectBaseNumTandPath, int32_t *newCandBaseNumAfterTandPathSE, int32_t *occNumArray, int32_t *occIndexArray, int32_t contigPos, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, contigPath_t *contigPath, graphtype *graph)
{
	tandemPathItem_t *tandemPathList;
	int32_t itemNumTandemPathList, overlapTandPathRemovedFlag;

	// get the tandem paths
	if(getTandemPathSE(&tandemPathList, &itemNumTandemPathList, contigPos, decisionTable, *readsNumDecisionTable, graph)==FAILED)
	{
		printf("line=%d, In %s(), cannot get tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

#if(TANDPATH_OUTPUT==YES)
	// output the paths
	if(outputTandPath(tandemPathList, itemNumTandemPathList, decisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot output tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}
#endif

	// remove reads in second path of two overlapped tandem paths
	if(removeOverlappedTandPath(&overlapTandPathRemovedFlag, tandemPathList, &itemNumTandemPathList, decisionTable, readsNumDecisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove reads in tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// decide the navigation from the paths
	if(overlapTandPathRemovedFlag==YES)
	{
		// decide the navigation
		if(decideNaviByTandPath(naviTandPathSE, maxBaseIndexAfterTandPathSE, incorrectBaseNumTandPath, newCandBaseNumAfterTandPathSE, occNumArray, occIndexArray, tandemPathList, itemNumTandemPathList)==FAILED)
		{
			printf("line=%d, In %s(), cannot navigate by tandem paths, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		*naviTandPathSE = NAVI_FAILED;
		*maxBaseIndexAfterTandPathSE = INT_MAX;
		*incorrectBaseNumTandPath = INT_MAX;
	}

#if(TANDPATH_OUTPUT==YES)
	// output the paths
	if(outputTandPath(tandemPathList, itemNumTandemPathList, decisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot output tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}
#endif

	// release the memory
	releaseTandemPathList(&tandemPathList);

	return SUCCESSFUL;
}

/**
 * Release the tandem paths.
 */
void releaseTandemPathList(tandemPathItem_t **tandemPathList)
{
	tandemPathItem_t *pTandPath, *pTandPathHead;
	dtReadInTandemPath_t *pDtRead, *pDtReadHead;

	pTandPathHead = *tandemPathList;
	while(pTandPathHead)
	{
		pDtReadHead = pTandPathHead->dtReadsList;
		while(pDtReadHead)
		{
			pDtRead = pDtReadHead->next;
			free(pDtReadHead);
			pDtReadHead = pDtRead;
		}
		pTandPathHead->dtReadsList = NULL;

		pTandPath = pTandPathHead->next;
		free(pTandPathHead);
		pTandPathHead = pTandPath;
	}
	*tandemPathList = NULL;
}

/**
 * Get the tandem paths from decision table for paired-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getTandemPathPE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, graphtype *graph)
{
	*tandemPathList = NULL;
	*itemNumTandemPathList = 0;

/*
	if(getTandemPathFromNewReadsPE(tandemPathList, itemNumTandemPathList, contigPos, kmerSeqIntAssembly, decisionTable, readsNumDecisionTable, graph)==FAILED)
	{
		printf("line=%d, In %s(), cannot get tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}
*/

	// get the tandem paths and their corresponding reads
	if(getTandemPathFromDecisionTablePE(tandemPathList, itemNumTandemPathList, contigPos, decisionTable, readsNumDecisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot get tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ############################ Debug information #########################
	// check the tandPath sequence length
//	if(checkTandPath(*tandemPathList, *itemNumTandemPathList)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot check tandem paths, error!\n", __LINE__, __func__);
//		return FAILED;
//	}
	// ############################ Debug information #########################

	// adjust the tandem paths by shift overlap
//	if(adjustTandemPathByShiftOverlap(*tandemPathList, itemNumTandemPathList, decisionTable, readsNumDecisionTable, dtRowHashtable, graph)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot adjust tandem paths, error!\n", __LINE__, __func__);
//		return FAILED;
//	}

	// adjust the tandem paths
	if(adjustTandemPath(*tandemPathList, itemNumTandemPathList)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ############################ Debug information #########################
	// check the tandPath sequence length
//	if(checkTandPath(*tandemPathList, *itemNumTandemPathList)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot check tandem paths, error!\n", __LINE__, __func__);
//		return FAILED;
//	}
	// ############################ Debug information #########################

	return SUCCESSFUL;
}

/**
 * Get the tandem paths from decision table for single-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getTandemPathSE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, graphtype *graph)
{
	*tandemPathList = NULL;
	*itemNumTandemPathList = 0;

/*
	if(getTandemPathFromNewReadsSE(tandemPathList, itemNumTandemPathList, contigPos, kmerSeqIntAssembly, decisionTable, readsNumDecisionTable, graph)==FAILED)
	{
		printf("line=%d, In %s(), cannot get tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}
*/

	// get the tandem paths and their corresponding reads
	if(getTandemPathFromDecisionTableSE(tandemPathList, itemNumTandemPathList, contigPos, decisionTable, readsNumDecisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot get tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ############################ Debug information #########################
	// check the tandPath sequence length
	if(checkTandPath(*tandemPathList, *itemNumTandemPathList)==FAILED)
	{
		printf("line=%d, In %s(), cannot check tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}
	// ############################ Debug information #########################

	// adjust the tandem paths
	if(adjustTandemPath(*tandemPathList, itemNumTandemPathList)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ############################ Debug information #########################
	// check the tandPath sequence length
	if(checkTandPath(*tandemPathList, *itemNumTandemPathList)==FAILED)
	{
		printf("line=%d, In %s(), cannot check tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}
	// ############################ Debug information #########################


	return SUCCESSFUL;
}

/**
 * Get the tandem paths from new paired-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getTandemPathFromNewReadsPE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, uint64_t *kmerseqInt, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, graphtype *graph)
{
	int32_t i, j;
	uint64_t tmp_kmerseq[entriesPerKmer], tmp_kmerseqRev[entriesPerKmer];
	kmertype *tmp_kmers[4][2];

	for(i=0; i<4; i++)
	{
		for(j=0; j<entriesPerKmer; j++)
			tmp_kmerseq[j] = kmerseqInt[j];

		tmp_kmerseq[entriesPerKmer-1] = ((kmerseqInt[entriesPerKmer-1] >> 2) << 2) | i;

		tmp_kmers[i][0] = getKmer(tmp_kmerseq, graph);
		tmp_kmers[i][1] = getReverseKmer(tmp_kmerseqRev, tmp_kmerseq, graph);
	}

	for(i=0; i<4; i++)
	{
		if(getTandemPathFromSingleKmerTandPathPE(tandemPathList, itemNumTandemPathList, contigPos, tmp_kmers[i], graph, decisionTable, readsNumDecisionTable)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the new number of reads of single kmer for multi-align, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the tandem paths from new single-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getTandemPathFromNewReadsSE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, uint64_t *kmerseqInt, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, graphtype *graph)
{
	int32_t i, j;
	uint64_t tmp_kmerseq[entriesPerKmer], tmp_kmerseqRev[entriesPerKmer];
	kmertype *tmp_kmers[4][2];

	for(i=0; i<4; i++)
	{
		for(j=0; j<entriesPerKmer; j++)
			tmp_kmerseq[j] = kmerseqInt[j];

		tmp_kmerseq[entriesPerKmer-1] = ((kmerseqInt[entriesPerKmer-1] >> 2) << 2) | i;

		tmp_kmers[i][0] = getKmer(tmp_kmerseq, graph);
		tmp_kmers[i][1] = getReverseKmer(tmp_kmerseqRev, tmp_kmerseq, graph);
	}

	for(i=0; i<4; i++)
	{
		if(getTandemPathFromSingleKmerTandPathSE(tandemPathList, itemNumTandemPathList, contigPos, tmp_kmers[i], graph, decisionTable, readsNumDecisionTable)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the new number of reads of single kmer for multi-align, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the tandem paths from new paired-ends using paired k-mer.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getTandemPathFromSingleKmerTandPathPE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, kmertype *tmp_kmers[2], graphtype *graph, assemblingreadtype *decisionTable, int32_t itemNumDT)
{
	int32_t returnCode;
	ridpostype *rid_pos_table;
	readBlock_t *readBlockArr;
	read_t *pRead;
	readseqBlock_t *readseqBlockArr;
	uint64_t *pReadseq, rid;

	int32_t i, posNum;
	int32_t maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3, matchFlag, mismatchNum;

	int32_t contigtailReadPos, perfectMatchFlag;
	char readseqTmp[MAX_READ_LEN_IN_BUF+1];
	tandemPathItem_t *targetTandPath;
	assemblingreadtype *dtReadPaired;


	readBlockArr = graph->readSet->readBlockArr;
	maxItemNumPerReadBlock = graph->readSet->maxItemNumPerReadBlock;
	readseqBlockArr = graph->readSet->readseqBlockArr;

	if(tmp_kmers[1])
	{
		rid_pos_table = tmp_kmers[1]->ppos;
		posNum = tmp_kmers[1]->arraysize;
		for(i=0; i<posNum; i++)
		{
			// get seqLen and errorRegLenEnd3
			rid = rid_pos_table[i].rid;
			blockRowRead = (rid - 1) / maxItemNumPerReadBlock;
			itemRowBlockRead = (rid - 1) % maxItemNumPerReadBlock;
			pRead = readBlockArr[blockRowRead].readArr + itemRowBlockRead;
			pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
			seqLen = pRead->seqlen;
			errorRegLenEnd3 = ceil(seqLen * ERROR_REGION_LEN_3End_FACTOR);

			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos>seqLen-kmerSize+1-errorRegLenEnd3)
			{
				if(existReadWithPosInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
				{
					returnCode = validReadPair(&dtReadPaired, rid_pos_table[i].rid, decisionTable, dtRowHashtable);
					if(returnCode==YES)
					{
						if(getMismatchNumWithContigPath(&matchFlag, &mismatchNum, pReadseq, seqLen, rid_pos_table[i].pos, ORIENTATION_MINUS, itemNumDT, contigPath)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)rid_pos_table[i].rid);
							return FAILED;
						}

						if(matchFlag==YES)
						{
							// get the read base sequence
							contigtailReadPos = seqLen - rid_pos_table[i].pos - 1;
							getReverseReadBaseFromPosByInt(readseqTmp, pReadseq, seqLen, 0, seqLen);

							perfectMatchFlag = YES;
							if((*itemNumTandemPathList)==0)
							{
								perfectMatchFlag = NO;
							}else
							{ // get the path it can be matched
								if(getMatchedTandPath(&targetTandPath, &perfectMatchFlag, *tandemPathList, readseqTmp, seqLen, contigtailReadPos)==FAILED)
								{
									printf("line=%d, In %s(), cannot get the matched tandem path, error!\n", __LINE__, __func__);
									return FAILED;
								}
							}

							if(perfectMatchFlag==YES)
							{ // add the read to the path
								if(addReadToTandPath(targetTandPath, readseqTmp, seqLen, contigtailReadPos, -1)==FAILED)
								{
									printf("line=%d, In %s(), cannot add new read to tandem path, error!\n", __LINE__, __func__);
									return FAILED;
								}
							}
							else
							{ // add tandem path
								if(addTandemPath(tandemPathList, itemNumTandemPathList, readseqTmp, seqLen, contigtailReadPos, contigPos, -1)==FAILED)
								{
									printf("line=%d, In %s(), cannot add new tandem path, error!\n", __LINE__, __func__);
									return FAILED;
								}
							}
						}
					}else if(returnCode==ERROR)
					{
						printf("In %s(), cannot valid read pair, error!\n", __func__);
						return FAILED;
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the tandem paths from new paired-ends using single k-mer.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getTandemPathFromSingleKmerTandPathSE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, kmertype *tmp_kmers[2], graphtype *graph, assemblingreadtype *decisionTable, int32_t itemNumDT)
{
	ridpostype *rid_pos_table;
	readBlock_t *readBlockArr;
	read_t *pRead;
	readseqBlock_t *readseqBlockArr;
	uint64_t *pReadseq, rid;

	int32_t i, posNum;
	int32_t maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3, matchFlag, mismatchNum;

	int32_t contigtailReadPos, perfectMatchFlag;
	char readseqTmp[MAX_READ_LEN_IN_BUF+1];
	tandemPathItem_t *targetTandPath;


	readBlockArr = graph->readSet->readBlockArr;
	maxItemNumPerReadBlock = graph->readSet->maxItemNumPerReadBlock;
	readseqBlockArr = graph->readSet->readseqBlockArr;

	if(tmp_kmers[0])
	{
		rid_pos_table = tmp_kmers[0]->ppos;
		posNum = tmp_kmers[0]->arraysize;
		for(i=0; i<posNum; i++)
		{
			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
			{
				rid = rid_pos_table[i].rid;
				blockRowRead = (rid - 1) / maxItemNumPerReadBlock;
				itemRowBlockRead = (rid - 1) % maxItemNumPerReadBlock;
				pRead = readBlockArr[blockRowRead].readArr + itemRowBlockRead;
				pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
				seqLen = pRead->seqlen;

				if(getMismatchNumWithContigPath(&matchFlag, &mismatchNum, pReadseq, seqLen, rid_pos_table[i].pos, ORIENTATION_PLUS, itemNumDT, contigPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)rid);
					return FAILED;
				}

				if(matchFlag==YES)
				{
					// get the read base sequence
					contigtailReadPos = rid_pos_table[i].pos + kmerSize - 3;
					getReadBaseFromPosByInt(readseqTmp, pReadseq, seqLen, 0, seqLen);

					perfectMatchFlag = YES;
					if((*itemNumTandemPathList)==0)
					{
						perfectMatchFlag = NO;
					}else
					{ // get the path it can be matched
						if(getMatchedTandPath(&targetTandPath, &perfectMatchFlag, *tandemPathList, readseqTmp, seqLen, contigtailReadPos)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the matched tandem path, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}

					if(perfectMatchFlag==YES)
					{ // add the read to the path
						if(addReadToTandPath(targetTandPath, readseqTmp, seqLen, contigtailReadPos, -1)==FAILED)
						{
							printf("line=%d, In %s(), cannot add new read to tandem path, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
					else
					{ // add tandem path
						if(addTandemPath(tandemPathList, itemNumTandemPathList, readseqTmp, seqLen, contigtailReadPos, contigPos, -1)==FAILED)
						{
							printf("line=%d, In %s(), cannot add new tandem path, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
				}
			}
		}
	}

	if(tmp_kmers[1])
	{
		rid_pos_table = tmp_kmers[1]->ppos;
		posNum = tmp_kmers[1]->arraysize;
		for(i=0; i<posNum; i++)
		{
			// get seqLen and errorRegLenEnd3
			rid = rid_pos_table[i].rid;
			blockRowRead = (rid - 1) / maxItemNumPerReadBlock;
			itemRowBlockRead = (rid - 1) % maxItemNumPerReadBlock;
			pRead = readBlockArr[blockRowRead].readArr + itemRowBlockRead;
			pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
			seqLen = pRead->seqlen;
			errorRegLenEnd3 = ceil(seqLen * ERROR_REGION_LEN_3End_FACTOR);

			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos>seqLen-kmerSize+1-errorRegLenEnd3)
			{
				if(existReadWithPosInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
				{
					if(getMismatchNumWithContigPath(&matchFlag, &mismatchNum, pReadseq, seqLen, rid_pos_table[i].pos, ORIENTATION_MINUS, itemNumDT, contigPath)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)rid);
						return FAILED;
					}

					if(matchFlag==YES)
					{
						// get the read base sequence
						contigtailReadPos = seqLen - rid_pos_table[i].pos - 1;
						getReverseReadBaseFromPosByInt(readseqTmp, pReadseq, seqLen, 0, seqLen);

						perfectMatchFlag = YES;
						if((*itemNumTandemPathList)==0)
						{
							perfectMatchFlag = NO;
						}else
						{ // get the path it can be matched
							if(getMatchedTandPath(&targetTandPath, &perfectMatchFlag, *tandemPathList, readseqTmp, seqLen, contigtailReadPos)==FAILED)
							{
								printf("line=%d, In %s(), cannot get the matched tandem path, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}

						if(perfectMatchFlag==YES)
						{ // add the read to the path
							if(addReadToTandPath(targetTandPath, readseqTmp, seqLen, contigtailReadPos, -1)==FAILED)
							{
								printf("line=%d, In %s(), cannot add new read to tandem path, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}
						else
						{ // add tandem path
							if(addTandemPath(tandemPathList, itemNumTandemPathList, readseqTmp, seqLen, contigtailReadPos, contigPos, -1)==FAILED)
							{
								printf("line=%d, In %s(), cannot add new tandem path, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the tandem paths from decision table for paired-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getTandemPathFromDecisionTablePE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable)
{
	int32_t i, perfectMatchFlag, readseqLen, contigtailReadPos;
	char readseqTmp[MAX_READ_LEN_IN_BUF+1];
	tandemPathItem_t *targetTandPath;

	for(i=0; i<readsNumDecisionTable; i++)
	{
		if(decisionTable[i].matedFlag==YES
			&& (decisionTable[i].unmatchBaseNum<=maxUnmatchBaseNumPerRead)
			&& (decisionTable[i].successiveAppearBases>=minSuccessiveAppearedBaseNum)
			&& decisionTable[i].basePos==decisionTable[i].lastMatchedBasePos)
		{
			// get the read base sequence
			if(decisionTable[i].orientation==ORIENTATION_PLUS)
			{
				contigtailReadPos = decisionTable[i].basePos;
				getReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, 0, decisionTable[i].seqlen);
			}else
			{
				contigtailReadPos = decisionTable[i].seqlen - decisionTable[i].basePos - 1;
				getReverseReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, 0, decisionTable[i].seqlen);
			}

			perfectMatchFlag = YES;
			readseqLen = decisionTable[i].seqlen;
			if((*itemNumTandemPathList)==0)
			{
				perfectMatchFlag = NO;
			}else
			{ // get the path it can be matched
				if(getMatchedTandPath(&targetTandPath, &perfectMatchFlag, *tandemPathList, readseqTmp, readseqLen, contigtailReadPos)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the matched tandem path, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}

			if(perfectMatchFlag==YES)
			{ // add the read to the path
				if(addReadToTandPath(targetTandPath, readseqTmp, readseqLen, contigtailReadPos, i)==FAILED)
				{
					printf("line=%d, In %s(), cannot add new read to tandem path, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
			else if((*itemNumTandemPathList)<20)
			{ // add tandem path
				if(addTandemPath(tandemPathList, itemNumTandemPathList, readseqTmp, readseqLen, contigtailReadPos, contigPos, i)==FAILED)
				{
					printf("line=%d, In %s(), cannot add new tandem path, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the tandem paths from decision table for single-ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getTandemPathFromDecisionTableSE(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, int32_t contigPos, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable)
{
	int32_t i, perfectMatchFlag, readseqLen, contigtailReadPos;
	char readseqTmp[MAX_READ_LEN_IN_BUF+1];
	tandemPathItem_t *targetTandPath;

	for(i=0; i<readsNumDecisionTable; i++)
	{
		if((decisionTable[i].unmatchBaseNum<=maxUnmatchBaseNumPerRead)
			&& (decisionTable[i].successiveAppearBases>=minSuccessiveAppearedBaseNum)
			&& decisionTable[i].basePos==decisionTable[i].lastMatchedBasePos)
		{
			// get the read base sequence
			if(decisionTable[i].orientation==ORIENTATION_PLUS)
			{
				contigtailReadPos = decisionTable[i].basePos;
				getReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, 0, decisionTable[i].seqlen);
			}else
			{
				contigtailReadPos = decisionTable[i].seqlen - decisionTable[i].basePos - 1;
				getReverseReadBaseFromPosByInt(readseqTmp, decisionTable[i].readseq, decisionTable[i].seqlen, 0, decisionTable[i].seqlen);
			}

			perfectMatchFlag = YES;
			readseqLen = decisionTable[i].seqlen;
			if((*itemNumTandemPathList)==0)
			{
				perfectMatchFlag = NO;
			}else
			{ // get the path it can be matched
				if(getMatchedTandPath(&targetTandPath, &perfectMatchFlag, *tandemPathList, readseqTmp, readseqLen, contigtailReadPos)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the matched tandem path, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}

			if(perfectMatchFlag==YES)
			{ // add the read to the path
				if(addReadToTandPath(targetTandPath, readseqTmp, readseqLen, contigtailReadPos, i)==FAILED)
				{
					printf("line=%d, In %s(), cannot add new read to tandem path, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
			else if((*itemNumTandemPathList)<20)
			{ // add tandem path
				if(addTandemPath(tandemPathList, itemNumTandemPathList, readseqTmp, readseqLen, contigtailReadPos, contigPos, i)==FAILED)
				{
					printf("line=%d, In %s(), cannot add new tandem path, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * get the matched tandem path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMatchedTandPath(tandemPathItem_t **targetTandPath, int32_t *perfectMatchFlag, tandemPathItem_t *tandemPathList, char *readseq, int32_t readseqLen, int32_t contigtailReadPos)
{
	tandemPathItem_t *tandPath;

	if(tandemPathList==NULL)
	{
		*targetTandPath = NULL;
		*perfectMatchFlag = NO;
		return SUCCESSFUL;
	}

	// check each tandem path
	tandPath = tandemPathList;
	while(tandPath)
	{
		if(getPerfectMatchFlagTandPath(perfectMatchFlag, tandPath->tandemPathStr, tandPath->tandemPathLen, tandPath->contigtailPathPos, readseq, readseqLen, contigtailReadPos)==FAILED)
		{
			printf("line=%d, In %s(), cannot determine the match flag of two tandem paths, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if((*perfectMatchFlag)==YES)
		{
			*targetTandPath = tandPath;
			break;
		}else
		{
			*targetTandPath = NULL;
		}

		tandPath = tandPath->next;
	}

	return SUCCESSFUL;
}

/**
 * Add tandem path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addTandemPath(tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, char *readseqTmp, int32_t readseqLen, int32_t contigtailReadPos, int32_t contigtailPos, int32_t dtRow)
{
	tandemPathItem_t *tandPath, *tandPathTail;
	dtReadInTandemPath_t *dtRead;

	if((*itemNumTandemPathList)>0)
	{
		// get the tail node
		tandPath = *tandemPathList;
		while(tandPath)
		{
			if(tandPath->next==NULL)
			{
				tandPathTail = tandPath;
				break;
			}
			tandPath = tandPath->next;
		}
	}else
		tandPathTail = NULL;

	// generate the new path node
	tandPath = (tandemPathItem_t*) malloc (sizeof(tandemPathItem_t));
	if(tandPath==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	strcpy(tandPath->tandemPathStr, readseqTmp);
	tandPath->tandemPathLen = readseqLen;
	tandPath->contigtailPathPos = contigtailReadPos;
	tandPath->contigtailPos = contigtailPos;
	tandPath->startContigPos = contigtailPos - contigtailReadPos;
	tandPath->endContigPos = contigtailPos + readseqLen - contigtailReadPos;
	tandPath->itemNumDtReadsList = 1;
	tandPath->fragSize = tandPath->fragSizeDif = 0;
	tandPath->shortFragFlag = NO;
	tandPath->validFlag = YES;
	tandPath->shiftMergeFlag = NO;
	tandPath->next = NULL;

	if(tandPathTail)
		tandPathTail->next = tandPath;
	else
		*tandemPathList = tandPath;
	(*itemNumTandemPathList) ++;

	// generate the new read node
	dtRead = (dtReadInTandemPath_t*) malloc(sizeof(dtReadInTandemPath_t));
	if(dtRead==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	tandPath->dtReadsList = dtRead;

	dtRead->dtRow = dtRow;
	dtRead->readseqLen = readseqLen;
	dtRead->contigtailReadPos = contigtailReadPos;
	dtRead->startRowInTandemPath = 0;
	dtRead->startContigPos = contigtailPos - contigtailReadPos;
	dtRead->endContigPos = contigtailPos + readseqLen - contigtailReadPos;
	dtRead->next = NULL;

	return SUCCESSFUL;
}

/**
 * Add read to tandem path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addReadToTandPath(tandemPathItem_t *targetTandPath, char *readseq, int32_t readseqLen, int32_t contigtailReadPos, int32_t dtRow)
{
	int32_t i, leftNewBaseLen, rightNewBaseLen, tandemPathLen;
	char *pathseq;
	dtReadInTandemPath_t *dtRead, *dtReadTail;

	// get the new increased bases of the path
	leftNewBaseLen = contigtailReadPos - targetTandPath->contigtailPathPos;
	rightNewBaseLen = (readseqLen - contigtailReadPos - 1) - (targetTandPath->tandemPathLen - targetTandPath->contigtailPathPos - 1);
	if(leftNewBaseLen<0)
		leftNewBaseLen = 0;
	if(rightNewBaseLen<0)
		rightNewBaseLen = 0;

	// update the path node
	if(leftNewBaseLen>0)
	{
		pathseq = targetTandPath->tandemPathStr;
		tandemPathLen = targetTandPath->tandemPathLen;
		for(i=tandemPathLen; i>=0; i--)	pathseq[i+leftNewBaseLen] = pathseq[i];
		for(i=0; i<leftNewBaseLen; i++) pathseq[i] = readseq[i];

		targetTandPath->tandemPathLen += leftNewBaseLen;
		targetTandPath->contigtailPathPos += leftNewBaseLen;
		targetTandPath->startContigPos -= leftNewBaseLen;
		if(targetTandPath->startContigPos<0)
			targetTandPath->startContigPos = 0;
	}
	if(rightNewBaseLen>0)
	{
		strcat(targetTandPath->tandemPathStr, readseq+readseqLen-rightNewBaseLen);
		targetTandPath->tandemPathLen += rightNewBaseLen;
		targetTandPath->endContigPos += rightNewBaseLen;
	}

	// update the dtReads
	dtRead = targetTandPath->dtReadsList;
	while(dtRead)
	{
		if(leftNewBaseLen>0)
			dtRead->startRowInTandemPath += leftNewBaseLen;
		dtRead = dtRead->next;
	}

	// add the new read
	dtRead = targetTandPath->dtReadsList;
	while(dtRead)
	{
		if(dtRead->next==NULL)
		{
			dtReadTail = dtRead;
			break;
		}
		dtRead = dtRead->next;
	}
	if(dtReadTail==NULL)
	{
		printf("line=%d, In %s(), dtReadTail==NULL, error!\n", __LINE__, __func__);
		return FAILED;
	}

	dtRead = (dtReadInTandemPath_t*) malloc(sizeof(dtReadInTandemPath_t));
	if(dtRead==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	dtRead->readseqLen = readseqLen;
	dtRead->dtRow = dtRow;
	dtRead->contigtailReadPos = contigtailReadPos;
	dtRead->startRowInTandemPath = targetTandPath->contigtailPathPos - contigtailReadPos;
	dtRead->startContigPos = targetTandPath->contigtailPos - contigtailReadPos;
	dtRead->endContigPos = dtRead->startContigPos + readseqLen;
	dtRead->next = NULL;

	dtReadTail->next = dtRead;
	targetTandPath->itemNumDtReadsList ++;

	return SUCCESSFUL;
}

/**
 * Adjust tandem path by shift overlap.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short adjustTandemPathByShiftOverlap(tandemPathItem_t *tandemPathList, int32_t *itemNumTandemPathList, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, graphtype *graph)
{
	int32_t shiftType;   // -1: No overlap, 1: shift in path1; 2: shift in path2
	int32_t mismatchNum, shiftSize, pathNodeDeleteFlag;
	tandemPathItem_t *targetTandPath, *preTandPath, *tandPath;
	dtReadInTandemPath_t *dtReadInTandemPath;
	assemblingreadtype *dtRead;

	pathNodeDeleteFlag = NO;
	targetTandPath = tandemPathList;
	while(targetTandPath)
	{
		preTandPath = targetTandPath;
		tandPath = preTandPath->next;
		while(tandPath)
		{
			// get the overlap size
			if(computeShiftSizeTandPathByShiftOp(&mismatchNum, &shiftType, &shiftSize, targetTandPath, tandPath)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the overlap type of two tandem paths, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(shiftType!=-1)
			{
				pathNodeDeleteFlag = YES;

				if(targetTandPath->itemNumDtReadsList<tandPath->itemNumDtReadsList)
				{
					if(exchangeNodeInfoTandPath(targetTandPath, tandPath)==FAILED)
					{
						printf("line=%d, In %s(), cannot exchange two tandem path nodes, error!\n", __LINE__, __func__);
						return FAILED;
					}

					targetTandPath->shiftMergeFlag = YES;
					shiftSize = -shiftSize;
					if(shiftType==1)
						shiftType = 2;
					else
						shiftType = 1;
				}

				// update the information of tandPath
				tandPath->contigtailPathPos += shiftSize;
				tandPath->startContigPos -= shiftSize;
				tandPath->endContigPos -= shiftSize;
				dtReadInTandemPath = tandPath->dtReadsList;
				while(dtReadInTandemPath)
				{
					dtReadInTandemPath->startRowInTandemPath += shiftSize;
					dtReadInTandemPath->contigtailReadPos += shiftSize;
					dtReadInTandemPath->startContigPos -= shiftSize;
					dtReadInTandemPath->endContigPos -= shiftSize;
					dtReadInTandemPath = dtReadInTandemPath->next;
				}

				// update the read in decision table
				dtReadInTandemPath = tandPath->dtReadsList;
				while(dtReadInTandemPath)
				{
					dtRead = decisionTable + dtReadInTandemPath->dtRow;
					if(dtRead->orientation==ORIENTATION_PLUS)
					{
						dtRead->firstBasePos = 0;
						dtRead->basePos += shiftSize;
						dtRead->successiveAppearBases = dtRead->basePos + 1;
					}else
					{
						dtRead->firstBasePos = dtRead->seqlen - 1;
						dtRead->basePos -= shiftSize;
						dtRead->successiveAppearBases = dtRead->seqlen - dtRead->basePos;
					}
					dtRead->firstContigPos = tandPath->contigtailPos + 2 - dtRead->successiveAppearBases;
					dtRead->successiveUnappearBases = 0;
					dtRead->unappearBlocksNum = 0;
					dtRead->matchBaseNum = dtRead->successiveAppearBases;
					dtRead->unmatchBaseNum = 0;
					dtRead->lastMatchedBasePos = dtRead->basePos;

					dtReadInTandemPath = dtReadInTandemPath->next;
				}

				// merge tandPath to targetTandPath
				if(mergeTandPathToTargetPath(targetTandPath, tandPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot merge two tandem paths, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// delete tandPath
				preTandPath->next = tandPath->next;
				free(tandPath);
				tandPath = preTandPath->next;
				(*itemNumTandemPathList) --;
			}else
			{
				preTandPath = tandPath;
				tandPath = tandPath->next;
			}
		}

		targetTandPath = targetTandPath->next;
	}

//	if(pathNodeDeleteFlag==YES)
//	{
//		// update the reads in decision table by shift overlap
//		if(updateReadsInDTByShiftOpTandPath(decisionTable, readsNumDecisionTable, tandemPathItem_t *tandemPathList, dtRowIndex_t **dtRowHashtable, graphtype *graph)==FAILED)
//		{
//			printf("line=%d, In %s(), cannot resort decision table, error!\n", __LINE__, __func__);
//			return FAILED;
//		}
//	}

	return SUCCESSFUL;
}

/**
 * Adjust tandem path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short adjustTandemPath(tandemPathItem_t *tandemPathList, int32_t *itemNumTandemPathList)
{
	int32_t *baseNumArray, maxRowNumBaseNumArray, contigtailRowBaseNumArray, colsNum;  // column=8

	if((*itemNumTandemPathList)<=1)
		return SUCCESSFUL;

	// get the maximal row number and the contigtailRowBaseNumArray
	if(getMaxRowNumBaseNumArrayTandPath(&maxRowNumBaseNumArray, &contigtailRowBaseNumArray, tandemPathList)==FAILED)
	{
		printf("line=%d, In %s(), get the maximal row number of the baseNumArray, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize the memory
	colsNum = 8;
	baseNumArray = (int32_t *)calloc(maxRowNumBaseNumArray*colsNum, sizeof(int32_t));
	if(baseNumArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the baseNumArray
	if(computeBaseNumArrayTandPath(baseNumArray, maxRowNumBaseNumArray, colsNum, contigtailRowBaseNumArray, tandemPathList)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the base number array for tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// remove error base in tandem path and update the path
	if(removeErrorBaseTandPath(baseNumArray, maxRowNumBaseNumArray, colsNum, contigtailRowBaseNumArray, tandemPathList)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove error bases for tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// merge tandem paths
	if(mergeAdjustedTandPath(tandemPathList, itemNumTandemPathList)==FAILED)
	{
		printf("line=%d, In %s(), cannot merge tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// free memory
	free(baseNumArray);

	return SUCCESSFUL;
}

/**
 * Get the maximal row number for tandem path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxRowNumBaseNumArrayTandPath(int32_t *maxRowNumBaseNumArray, int32_t *contigtailRowBaseNumArray, tandemPathItem_t *tandemPathList)
{
	int32_t maxLeftBaseNum, maxRightBaseNum, tmpBaseNum;
	tandemPathItem_t *tandPath;

	maxLeftBaseNum = maxRightBaseNum = 0;
	tandPath = tandemPathList;
	while(tandPath)
	{
		// the left side
		tmpBaseNum = tandPath->contigtailPathPos + 1;
		if(maxLeftBaseNum<tmpBaseNum)
			maxLeftBaseNum = tmpBaseNum;

		// the right side
		tmpBaseNum = tandPath->tandemPathLen - 1 - tandPath->contigtailPathPos;
		if(maxRightBaseNum<tmpBaseNum)
			maxRightBaseNum = tmpBaseNum;

		tandPath = tandPath->next;
	}

	*maxRowNumBaseNumArray = maxLeftBaseNum + maxRightBaseNum;
	*contigtailRowBaseNumArray = maxLeftBaseNum - 1;

	if((*maxRowNumBaseNumArray)<0 || (*contigtailRowBaseNumArray)<0)
	{
		printf("line=%d, In %s(), maxRowNumBaseNumArray=%d, contigtailRowBaseNumArray=%d, error!\n", __LINE__, __func__, *maxRowNumBaseNumArray, *contigtailRowBaseNumArray);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Compute the base number array for tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeBaseNumArrayTandPath(int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, int32_t contigtailRowBaseNumArray, tandemPathItem_t *tandemPathList)
{
	int32_t i, j, baseNumTmp, supportReadsNum, startRowBaseNumArray, baseInt;
	int32_t value, maxValue, secValue, maxIndex, secIndex, subSum;
	char *pathseq;
	tandemPathItem_t *tandPath;

	// fill the base number
	tandPath = tandemPathList;
	while(tandPath)
	{
		pathseq = tandPath->tandemPathStr;
		supportReadsNum = tandPath->itemNumDtReadsList;

		startRowBaseNumArray = contigtailRowBaseNumArray - tandPath->contigtailPathPos;
		baseNumTmp = tandPath->tandemPathLen;
		for(i=0; i<baseNumTmp; i++)
		{
			switch(pathseq[i])
			{
				case 'A':
				case 'a': baseInt = 0; break;
				case 'C':
				case 'c': baseInt = 1; break;
				case 'G':
				case 'g': baseInt = 2; break;
				case 'T':
				case 't': baseInt = 3; break;
				case '-': baseInt = 4; break;
				default: printf("line=%d, invalid base %d\n, error!\n", __LINE__, pathseq[i]); return FAILED;
			}

			baseNumArray[(startRowBaseNumArray+i) * colsNum + baseInt] += supportReadsNum; // for A, C, G, T, -
			baseNumArray[(startRowBaseNumArray+i) * colsNum + 5] += supportReadsNum; // // for A+C+G+T+-
		}

		tandPath = tandPath->next;
	}

	// get the maxOcc and secOcc
	for(i=0; i<maxRowNumBaseNumArray; i++)
	{
		maxValue = secValue = -1;
		maxIndex = secIndex = -1;
		for(j=0; j<5; j++)
		{
			value = baseNumArray[i*colsNum + j];
			if(value>maxValue)
			{
				secValue = maxValue;
				secIndex = maxIndex;
				maxValue = value;
				maxIndex = j;
			}else if(value>secValue)
			{
				secValue = value;
				secIndex = j;
			}
		}

		subSum = baseNumArray[i*colsNum + 5];

		baseNumArray[i*colsNum + 6] = maxIndex;
		baseNumArray[i*colsNum + 7] = secIndex;
	}

	return SUCCESSFUL;
}

/**
 * Remove error bases for tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeErrorBaseTandPath(int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum, int32_t contigtailRowBaseNumArray, tandemPathItem_t *tandemPathList)
{
	int32_t i, k;
	int32_t maxBaseIndex, maxValue, secondBaseIndex, secondValue, sum;
	int32_t gapValue, baseValue, errrBaseIndex, pathseqRow, tmpRow;
	char *pathseq, errBase, newBase;
	tandemPathItem_t *tandPath;

	for(i=0; i<maxRowNumBaseNumArray; i++)
	{
		secondBaseIndex = baseNumArray[i*colsNum + 7];
		if(secondBaseIndex<4)
		{
			maxBaseIndex = baseNumArray[i*colsNum + 6];
			maxValue = baseNumArray[i*colsNum + maxBaseIndex];
			secondValue =baseNumArray[i*colsNum + secondBaseIndex];
			gapValue = baseNumArray[i*colsNum + 4];
			sum = baseNumArray[i*colsNum + 5];
			if(secondValue==1 && gapValue==0 && sum>=10)
			{
				for(k=0; k<4; k++)
				{
					baseValue = baseNumArray[i*colsNum + k];
					if(baseValue==1)
					{
						errrBaseIndex = k;
						switch(errrBaseIndex)
						{
							case 0: errBase = 'A'; break;
							case 1: errBase = 'C'; break;
							case 2: errBase = 'G'; break;
							case 3: errBase = 'T'; break;
						}

						switch(maxBaseIndex)
						{
							case 0: newBase = 'A'; break;
							case 1: newBase = 'C'; break;
							case 2: newBase = 'G'; break;
							case 3: newBase = 'T'; break;
						}

						tmpRow = i;

						// get the single-path according to the single base
						tandPath = tandemPathList;
						while(tandPath)
						{
							pathseq = tandPath->tandemPathStr;
							pathseqRow = tmpRow - (contigtailRowBaseNumArray - tandPath->contigtailPathPos);
							if(pathseqRow>=0 && pathseqRow<tandPath->tandemPathLen)
							{
								if(pathseq[pathseqRow]==errBase)
								{
									pathseq[pathseqRow] = newBase;
									break;
								}
							}

							tandPath = tandPath->next;
						}
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Merge adjusted tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mergeAdjustedTandPath(tandemPathItem_t *tandemPathList, int32_t *itemNumTandemPathList)
{
	int32_t matchFlag;
	tandemPathItem_t *targetTandPath, *preTandPath, *tandPath;
	dtReadInTandemPath_t *dtRead, *dtReadHead;

	targetTandPath = tandemPathList;
	while(targetTandPath)
	{
		preTandPath = targetTandPath;
		tandPath = preTandPath->next;
		while(tandPath)
		{
			// get the match flag
			if(getMatchFlagTandPath(&matchFlag, targetTandPath, tandPath, 1)==FAILED)
			{
				printf("line=%d, In %s(), cannot determine the match flag of two tandem paths, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(matchFlag==YES)
			{
				// merge the two paths
				if(mergeTandPathToTargetPath(targetTandPath, tandPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot merge two tandem paths, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// delete path node
				dtReadHead = tandPath->dtReadsList;
				while(dtReadHead)
				{
					dtRead = dtReadHead->next;
					free(dtReadHead);
					dtReadHead = dtRead;
				}
				tandPath->dtReadsList = NULL;

				preTandPath->next = tandPath->next;
				free(tandPath);
				tandPath = preTandPath->next;
				(*itemNumTandemPathList) --;
			}else
			{
				preTandPath = tandPath;
				tandPath = tandPath->next;
			}
		}

		targetTandPath = targetTandPath->next;
	}

	return SUCCESSFUL;
}

/**
 * Determine whether the two tandem path sequences perfectly matched.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getPerfectMatchFlagTandPath(int32_t *perfectMatchFlag, char *pathseq1, int32_t pathlen1, int32_t contigtailPathPos1, char *pathseq2, int32_t pathlen2, int32_t contigtailPathPos2)
{
	int32_t i, startRow1, startRow2, seqLenTmp1, seqLenTmp2, shareLen;

	if(contigtailPathPos1<contigtailPathPos2)
	{
		startRow1 = 0;
		startRow2 = contigtailPathPos2 - contigtailPathPos1;
	}else
	{
		startRow1 = contigtailPathPos1 - contigtailPathPos2;
		startRow2 = 0;
	}

	seqLenTmp1 = pathlen1 - startRow1;
	seqLenTmp2 = pathlen2 - startRow2;
	if(seqLenTmp1<seqLenTmp2)
		shareLen = seqLenTmp1;
	else
		shareLen = seqLenTmp2;

	if(shareLen>0)
	{
		*perfectMatchFlag = YES;
		for(i=0; i<shareLen; i++)
		{
			if(pathseq1[startRow1+i]!=pathseq2[startRow2+i])
			{
				*perfectMatchFlag = NO;
				break;
			}
		}
	}else
		*perfectMatchFlag = NO;

	return SUCCESSFUL;
}

/**
 * Determine whether the two tandem path sequences matched.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMatchFlagTandPath(int32_t *matchFlag, tandemPathItem_t *targetTandPath, tandemPathItem_t *tandPath, int32_t mismatchThreshold)
{
	int32_t i, startRow1, startRow2, seqLenTmp1, seqLenTmp2, shareLen, mismatchNum, correctedNum;

	if(targetTandPath->contigtailPathPos<tandPath->contigtailPathPos)
	{
		startRow1 = 0;
		startRow2 = tandPath->contigtailPathPos - targetTandPath->contigtailPathPos;
	}else
	{
		startRow1 = targetTandPath->contigtailPathPos - tandPath->contigtailPathPos;
		startRow2 = 0;
	}

	seqLenTmp1 = targetTandPath->tandemPathLen - startRow1;
	seqLenTmp2 = tandPath->tandemPathLen - startRow2;
	if(seqLenTmp1<seqLenTmp2)
		shareLen = seqLenTmp1;
	else
		shareLen = seqLenTmp2;

	if(shareLen>0)
	{
		mismatchNum = 0;
		for(i=0; i<shareLen; i++)
		{
			if(targetTandPath->tandemPathStr[startRow1+i]!=tandPath->tandemPathStr[startRow2+i])
			{
				mismatchNum ++;
				//if(mismatchNum>MAX_INCOR_BASE_NUM_TANDPATH)
				if(mismatchNum>mismatchThreshold)
					break;
			}
		}

		if(mismatchNum<=mismatchThreshold)
		{
			*matchFlag = YES;
		}
		//else if(mismatchNum>MAX_INCOR_BASE_NUM_TANDPATH)
		else
		{
			*matchFlag = NO;
		}
/*
		else
		{
			correctedNum = 0;
			for(i=0; i<shareLen; i++)
			{
				if(targetTandPath->tandemPathStr[startRow1+i]!=tandPath->tandemPathStr[startRow2+i])
				{
					if(targetTandPath->itemNumDtReadsList>tandPath->itemNumDtReadsList)
						tandPath->tandemPathStr[startRow2+i] = targetTandPath->tandemPathStr[startRow1+i];
					else
						targetTandPath->tandemPathStr[startRow1+i] = tandPath->tandemPathStr[startRow2+i];

					correctedNum ++;
					if(correctedNum>=mismatchNum)
						break;
				}
			}

			*matchFlag = YES;
		}
*/
	}else
		*matchFlag = NO;

	return SUCCESSFUL;
}

/**
 * Merge the tandem path to the target tandem path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mergeTandPathToTargetPath(tandemPathItem_t *targetTandPath, tandemPathItem_t *tandPath)
{
	char *readseq;
	dtReadInTandemPath_t *dtRead;

	readseq = (char*) calloc(tandPath->tandemPathLen+1, sizeof(char));
	if(readseq==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	dtRead = tandPath->dtReadsList;
	while(dtRead)
	{
		strcpy(readseq, tandPath->tandemPathStr+dtRead->startRowInTandemPath);
		readseq[dtRead->readseqLen] = '\0';

		if(addReadToTandPath(targetTandPath, readseq, dtRead->readseqLen, dtRead->contigtailReadPos, dtRead->dtRow)==FAILED)
		{
			printf("line=%d, In %s(), cannot add the read to tandem path, error!\n", __LINE__, __func__);
			return FAILED;
		}

		dtRead = dtRead->next;
	}

	free(readseq);

	return SUCCESSFUL;
}

/**
 * Update the reads in decision table by shift overlap tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
/*
short updateReadsInDTByShiftOpTandPath(tandemPathItem_t *tandemPathList, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable, graphtype *graph)
{
	int32_t startRowKmerTandPath, pathLen, rightSize, entriesNumPathseq, baseNumLastEntry, ridposnum;
	char *pathseq;
	uint64_t pathseqInt[200], kmerSeqIntTmp[entriesPerKmer], hashcode;
	kmertype *kmerTmp;
	ridpostype *ridpostable;
	tandemPathItem_t *tandPath;

	tandPath = tandemPathList;
	while(tandPath)
	{
		if(tandPath->shiftMergeFlag==YES)
		{
			pathseq = tandPath->tandemPathStr;
			pathLen = tandPath->tandemPathLen;

			if(tandPath->tandemPathLen-tandPath->contigtailPathPos-1>=10)
				rightSize = 10;
			else
				rightSize = tandPath->tandemPathLen - tandPath->contigtailPathPos - 1;
			startRowKmerTandPath = tandPath->contigtailPathPos - kmerSize + rightSize;
			if(startRowKmerTandPath<0)
			{
				tandPath = tandPath->next;
				continue;
			}

			// generate the readseq
			entriesNumPathseq = ((pathLen - 1) >> 5) + 1;
			baseNumLastEntry = ((pathLen - 1) % 32) + 1;
			if(generateReadseqInt(pathseqInt, pathseq, pathLen, entriesNumPathseq)==FAILED)
			{
				printf("line=%d, In %s(), unknownBaseNum=%d, cannot generate the readseq integer sequence, error!\n", __LINE__, __func__, unknownBaseNum);
				return FAILED;
			}

			// generate the kmer integer sequence
			if(generateKmerSeqIntFromReadset(kmerSeqIntTmp, pathseqInt, startRowKmerTandPath, entriesNumPathseq, baseNumLastEntry)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			hashcode = kmerhashInt(kmerSeqIntTmp);
			kmerTmp = getKmerByHash(hashcode, kmerSeqIntTmp, graph);
			ridpostable = kmerTmp->ppos;
			ridposnum = kmerTmp->arraysize;



		}

		tandPath = tandPath->next;
	}



	return SUCCESSFUL;
}
*/

/**
 * Merge adjusted tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputTandPath(tandemPathItem_t *tandemPathList, int32_t itemNumTandemPathList, assemblingreadtype *decisionTable)
{
	int64_t readID;
	char readseqTmp[MAX_READ_LEN_IN_BUF+1];
	tandemPathItem_t *tandPath;
	dtReadInTandemPath_t *dtRead;

	printf("==== itemNumTandemPathList=%d\n", itemNumTandemPathList);
	tandPath = tandemPathList;
	while(tandPath)
	{
		printf("%d\t%s\t%d\n", tandPath->itemNumDtReadsList, tandPath->tandemPathStr, tandPath->contigtailPathPos);
		printf("\tvalidFlag=%d, shortFragFlag=%d, fragSize=%.2f, fragSizeDif=%.2f, matchWithContigFlag=%d, mismatchNum=%d, "
				"\n\tcontigtailPos=%d, startContigPos=%d, endContigPos=%d, tandemPathLen=%d\n",
				tandPath->validFlag, tandPath->shortFragFlag, tandPath->fragSize, tandPath->fragSizeDif, tandPath->matchWithContigFlag, tandPath->mismatchNum,
				tandPath->contigtailPos, tandPath->startContigPos, tandPath->endContigPos, tandPath->tandemPathLen);

		dtRead = tandPath->dtReadsList;
		while(dtRead)
		{
			strcpy(readseqTmp, tandPath->tandemPathStr+dtRead->startRowInTandemPath);
			readseqTmp[dtRead->readseqLen] = '\0';

			if(dtRead->dtRow>=0)
			{
				readID = decisionTable[dtRead->dtRow].rid;
			}else
				readID = -1;
			printf("\t%ld\t%s\t%d\n", readID, readseqTmp, dtRead->contigtailReadPos);

			dtRead = dtRead->next;
		}

		tandPath = tandPath->next;
	}

	return SUCCESSFUL;
}

/**
 * Remove tandem paths with short fragment size.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeShortFragSizeTandPath(int32_t *shortFragSizeRemovedFlag, tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable, contigPath_t *contigPath)
{
	double averFragSize, averFragSizeDif, overlapFactorTmp;
	int32_t shortFragSizeFlag, matchWithContigFlag, maxSupportReadsNum, mismatchNum, maxOverlapSizeWithContig, validTandPathFlag;
	tandemPathItem_t *preTandPath, *tandPath, *tandPathTmp;
	dtReadInTandemPath_t *dtRead, *dtReadHead;
	int32_t i, j, k, startPos, endPos, baseInt, matchBeforeFlag, maxOverlapLen;
	char *pseq;

	if((*itemNumTandemPathList)>=2)
	{
		preTandPath = NULL;
		tandPath = *tandemPathList;
		while(tandPath)
		{
			// get the average fragment size
			if(computeAverFragSizeTandPath(&averFragSize, tandPath, decisionTable, dtRowHashtable)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the overlap type of two tandem paths, error!\n", __LINE__, __func__);
				return FAILED;
			}

			averFragSizeDif = meanSizeInsert - averFragSize;
			if(averFragSizeDif<0)
				averFragSizeDif = -averFragSizeDif;
			tandPath->fragSize = averFragSize;
			tandPath->fragSizeDif = averFragSizeDif;
			//if(averFragSizeDif>2*standardDev || averFragSizeDif>0.2*meanSizeInsert)
			if(averFragSizeDif>2*standardDev || averFragSizeDif>0.1*meanSizeInsert)
				tandPath->shortFragFlag = YES;
			else
				tandPath->shortFragFlag = NO;

			// get the contig match flag
			if(getContigMatchFlagTandPath(&tandPath->matchWithContigFlag, &tandPath->mismatchNum, tandPath, contigPath)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the overlap type of two tandem paths, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if((tandPath->shortFragFlag==YES && tandPath->matchWithContigFlag==NO) || tandPath->mismatchNum>=2)
			{ // remove condition

				*shortFragSizeRemovedFlag = YES;

				// remove the error tandem path and its reads
				dtReadHead = tandPath->dtReadsList;
				while(dtReadHead)
				{
					dtRead = dtReadHead->next;
					if(dtReadHead->dtRow>=0)
						decisionTable[dtReadHead->dtRow].status = FAILED_STATUS;
					free(dtReadHead);
					tandPath->itemNumDtReadsList --;
					dtReadHead = dtRead;
				}
				tandPath->dtReadsList = NULL;

				if(tandPath->itemNumDtReadsList!=0)
				{
					printf("line=%d, In %s(), invalid itemNumDtReadsList=%d, error!\n", __LINE__, __func__, tandPath->itemNumDtReadsList);
					return FAILED;
				}

				// delete the nodes
				tandPathTmp = tandPath->next;
				if(preTandPath)
					preTandPath->next = tandPathTmp;
				else
					*tandemPathList = tandPathTmp;
				free(tandPath);
				(*itemNumTandemPathList) --;
				tandPath = tandPathTmp;
			}else
			{
				preTandPath = tandPath;
				tandPath = tandPath->next;
			}

			if((*itemNumTandemPathList)<2)
				break;
		}
	}

	if((*itemNumTandemPathList)>=2)
	{
		// get the maximum read count
		maxSupportReadsNum = 0;
		tandPath = *tandemPathList;
		while(tandPath)
		{
			if(maxSupportReadsNum<tandPath->itemNumDtReadsList && tandPath->shortFragFlag==NO)
				maxSupportReadsNum = tandPath->itemNumDtReadsList;
			tandPath = tandPath->next;
		}

		// check the second round
		preTandPath = NULL;
		tandPath = *tandemPathList;
		while(tandPath)
		{
			if(tandPath->shortFragFlag==YES && tandPath->itemNumDtReadsList<0.5*maxSupportReadsNum && tandPath->fragSizeDif>0.2*meanSizeInsert)
			{ // remove condition

				*shortFragSizeRemovedFlag = YES;

				// remove the error tandem path and its reads
				dtReadHead = tandPath->dtReadsList;
				while(dtReadHead)
				{
					dtRead = dtReadHead->next;
					if(dtReadHead->dtRow>=0)
						decisionTable[dtReadHead->dtRow].status = FAILED_STATUS;
					free(dtReadHead);
					tandPath->itemNumDtReadsList --;
					dtReadHead = dtRead;
				}
				tandPath->dtReadsList = NULL;

				if(tandPath->itemNumDtReadsList!=0)
				{
					printf("line=%d, In %s(), invalid itemNumDtReadsList=%d, error!\n", __LINE__, __func__, tandPath->itemNumDtReadsList);
					return FAILED;
				}

				// delete the nodes
				tandPathTmp = tandPath->next;
				if(preTandPath)
					preTandPath->next = tandPathTmp;
				else
					*tandemPathList = tandPathTmp;
				free(tandPath);
				(*itemNumTandemPathList) --;
				tandPath = tandPathTmp;
			}else
			{
				preTandPath = tandPath;
				tandPath = tandPath->next;
			}
		}
	}


	if((*itemNumTandemPathList)>=2)
	{
		// check the second round
		preTandPath = NULL;
		tandPath = *tandemPathList;
		while(tandPath)
		{
			pseq = tandPath->tandemPathStr;
			startPos = itemNumContigArr - 2 * meanSizeInsert;
			if(startPos<0)
				startPos = 0;
			endPos = itemNumContigArr - tandPath->tandemPathLen;
			matchBeforeFlag = NO;
			for(i=startPos; i<=endPos; i++)
			{
				mismatchNum = 0;
				for(j=i, k=0; k<tandPath->tandemPathLen; j++, k++)
				{
					switch(pseq[k])
					{
						case 'A': baseInt = 0; break;
						case 'C': baseInt = 1; break;
						case 'G': baseInt = 2; break;
						case 'T': baseInt = 3; break;
						default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, pseq[k]); return FAILED;
					}

					if(baseInt!=contigArr[j].base)
					{
						mismatchNum ++;
						if(mismatchNum>3)
							break;
					}
				}

				if(k==tandPath->tandemPathLen && mismatchNum<=3)
				{
					matchBeforeFlag = YES;
					break;
				}
			}

			if(matchBeforeFlag==YES)
			{
				*shortFragSizeRemovedFlag = YES;

				// remove the error tandem path and its reads
				dtReadHead = tandPath->dtReadsList;
				while(dtReadHead)
				{
					dtRead = dtReadHead->next;
					if(dtReadHead->dtRow>=0)
						decisionTable[dtReadHead->dtRow].status = FAILED_STATUS;
					free(dtReadHead);
					tandPath->itemNumDtReadsList --;
					dtReadHead = dtRead;
				}
				tandPath->dtReadsList = NULL;

				if(tandPath->itemNumDtReadsList!=0)
				{
					printf("line=%d, In %s(), invalid itemNumDtReadsList=%d, error!\n", __LINE__, __func__, tandPath->itemNumDtReadsList);
					return FAILED;
				}

				// delete the nodes
				tandPathTmp = tandPath->next;
				if(preTandPath)
					preTandPath->next = tandPathTmp;
				else
					*tandemPathList = tandPathTmp;
				free(tandPath);
				(*itemNumTandemPathList) --;
				tandPath = tandPathTmp;
			}else
			{
				preTandPath = tandPath;
				tandPath = tandPath->next;
			}
		}
	}


/*
	if((*itemNumTandemPathList)>=2)
	{
		preTandPath = NULL;
		tandPath = *tandemPathList;
		while(tandPath)
		{
			// get the average fragment size
			if(computeAverFragSizeTandPath(&averFragSize, tandPath, decisionTable, dtRowHashtable)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the overlap type of two tandem paths, error!\n", __LINE__, __func__);
				return FAILED;
			}

			averFragSizeDif = meanSizeInsert - averFragSize;
			if(averFragSizeDif<0)
				averFragSizeDif = -averFragSizeDif;
			tandPath->fragSize = averFragSize;
			tandPath->fragSizeDif = averFragSizeDif;
			//if(averFragSizeDif>2*standardDev || averFragSizeDif>0.2*meanSizeInsert)
			if(averFragSizeDif>2*standardDev || averFragSizeDif>0.1*meanSizeInsert)
				tandPath->shortFragFlag = YES;
			else
				tandPath->shortFragFlag = NO;

			// get the contig match flag
			if(getContigMatchFlagTandPath(&tandPath->matchWithContigFlag, &tandPath->mismatchNum, tandPath, contigPath)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the overlap type of two tandem paths, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if((tandPath->shortFragFlag==YES && tandPath->matchWithContigFlag==NO) || tandPath->mismatchNum>=2)
			{ // remove condition
				tandPath->validFlag = NO;
			}

			tandPath = tandPath->next;

			if((*itemNumTandemPathList)<2)
				break;
		}
	}

	if((*itemNumTandemPathList)>=2)
	{
		// get the maximum read count
		maxSupportReadsNum = 0;
		tandPath = *tandemPathList;
		while(tandPath)
		{
			if(maxSupportReadsNum<tandPath->itemNumDtReadsList && tandPath->shortFragFlag==NO)
				maxSupportReadsNum = tandPath->itemNumDtReadsList;
			tandPath = tandPath->next;
		}

		// check the second round
		preTandPath = NULL;
		tandPath = *tandemPathList;
		while(tandPath)
		{
			if(tandPath->shortFragFlag==YES && tandPath->itemNumDtReadsList<0.5*maxSupportReadsNum && tandPath->fragSizeDif>0.2*meanSizeInsert)
			{ // remove condition
				tandPath->validFlag = NO;
			}
			tandPath = tandPath->next;
		}
	}


	if((*itemNumTandemPathList)>=2)
	{
		// check the second round
		preTandPath = NULL;
		tandPath = *tandemPathList;
		while(tandPath)
		{
			pseq = tandPath->tandemPathStr;
			startPos = itemNumContigArr - 2 * meanSizeInsert;
			if(startPos<0)
				startPos = 0;
			endPos = itemNumContigArr - tandPath->tandemPathLen;
			matchBeforeFlag = NO;
			for(i=startPos; i<=endPos; i++)
			{
				mismatchNum = 0;
				for(j=i, k=0; k<tandPath->tandemPathLen; j++, k++)
				{
					switch(pseq[k])
					{
						case 'A': baseInt = 0; break;
						case 'C': baseInt = 1; break;
						case 'G': baseInt = 2; break;
						case 'T': baseInt = 3; break;
						default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, pseq[k]); return FAILED;
					}

					if(baseInt!=contigArr[j].base)
					{
						mismatchNum ++;
						if(mismatchNum>3)
							break;
					}
				}

				if(k==tandPath->tandemPathLen && mismatchNum<=3)
				{
					matchBeforeFlag = YES;
					break;
				}
			}

			if(matchBeforeFlag==YES)
			{
				tandPath->validFlag = NO;
			}
			tandPath = tandPath->next;
		}
	}


	if((*itemNumTandemPathList)>=2)
	{
		validTandPathFlag = NO;
		tandPath = *tandemPathList;
		while(tandPath)
		{
			if(tandPath->validFlag==YES && tandPath->shortFragFlag==NO)
			{
				validTandPathFlag = YES;
				break;
			}
			tandPath = tandPath->next;
		}

		if(validTandPathFlag==YES)
		{
			tandPath = *tandemPathList;
			while(tandPath)
			{
				if(tandPath->shortFragFlag==YES)
					tandPath->validFlag = NO;
				tandPath = tandPath->next;
			}
		}
	}


	// check the first round
	*shortFragSizeRemovedFlag = NO;

	// check the average contigtailPathPos
	if((*itemNumTandemPathList)>=2)
	{
		preTandPath = NULL;
		tandPath = *tandemPathList;
		while(tandPath)
		{
			if(tandPath->validFlag==NO)
			{ // remove condition
				*shortFragSizeRemovedFlag = YES;

				// remove the error tandem path and its reads
				dtReadHead = tandPath->dtReadsList;
				while(dtReadHead)
				{
					dtRead = dtReadHead->next;
					if(dtReadHead->dtRow>=0)
						decisionTable[dtReadHead->dtRow].status = FAILED_STATUS;
					free(dtReadHead);
					tandPath->itemNumDtReadsList --;
					dtReadHead = dtRead;
				}
				tandPath->dtReadsList = NULL;

				if(tandPath->itemNumDtReadsList!=0)
				{
					printf("line=%d, In %s(), invalid itemNumDtReadsList=%d, error!\n", __LINE__, __func__, tandPath->itemNumDtReadsList);
					return FAILED;
				}

				// delete the nodes
				tandPathTmp = tandPath->next;
				if(preTandPath)
					preTandPath->next = tandPathTmp;
				else
					*tandemPathList = tandPathTmp;
				free(tandPath);
				(*itemNumTandemPathList) --;
				tandPath = tandPathTmp;
			}else
			{
				preTandPath = tandPath;
				tandPath = tandPath->next;
			}
		}
	}


	if((*itemNumTandemPathList)>=2)
	{
		overlapFactorTmp = 0.9;
		for(i=0; i<3; i++)
		{
			validTandPathFlag = NO;
			tandPath = *tandemPathList;
			while(tandPath)
			{
				if(tandPath->validFlag==YES && tandPath->fragSizeDif<3*standardDev && tandPath->contigtailPathPos>overlapFactorTmp*readLen)
				{
					validTandPathFlag = YES;
					break;
				}
				tandPath = tandPath->next;
			}
			if(validTandPathFlag==YES)
			{
				tandPath = *tandemPathList;
				while(tandPath)
				{
					if(tandPath->validFlag==YES && tandPath->shortFragFlag==YES && tandPath->fragSizeDif>2*standardDev && tandPath->contigtailPathPos<(overlapFactorTmp-0.2)*readLen)
						tandPath->validFlag = NO;
					tandPath = tandPath->next;
				}
			}
		}

		validTandPathFlag = NO;
		tandPath = *tandemPathList;
		while(tandPath)
		{
			if(tandPath->validFlag==YES && tandPath->fragSizeDif<2*standardDev && tandPath->contigtailPathPos>0.7*readLen)
			{
				validTandPathFlag = YES;
				break;
			}
			tandPath = tandPath->next;
		}
		if(validTandPathFlag==YES)
		{
			tandPath = *tandemPathList;
			while(tandPath)
			{
				if(tandPath->validFlag==YES && tandPath->contigtailPathPos<0.3*readLen)
					tandPath->validFlag = NO;
				tandPath = tandPath->next;
			}
		}


		validTandPathFlag = NO;
		tandPath = *tandemPathList;
		while(tandPath)
		{
			if(tandPath->validFlag==YES && tandPath->fragSizeDif<standardDev && tandPath->contigtailPathPos>0.9*readLen)
			{
				validTandPathFlag = YES;
				break;
			}
			tandPath = tandPath->next;
		}

		if(validTandPathFlag==YES)
		{
			tandPath = *tandemPathList;
			while(tandPath)
			{
				if(tandPath->validFlag==YES && tandPath->fragSizeDif<standardDev && tandPath->contigtailPathPos<0.5*readLen)
					tandPath->validFlag = NO;
				tandPath = tandPath->next;
			}
		}

	}

	// check the average contigtailPathPos
	if((*itemNumTandemPathList)>=2)
	{
		preTandPath = NULL;
		tandPath = *tandemPathList;
		while(tandPath)
		{
			if(tandPath->validFlag==NO)
			{ // remove condition
				*shortFragSizeRemovedFlag = YES;

				// remove the error tandem path and its reads
				dtReadHead = tandPath->dtReadsList;
				while(dtReadHead)
				{
					dtRead = dtReadHead->next;
					//if(dtReadHead->dtRow>=0)
					//	decisionTable[dtReadHead->dtRow].status = FAILED_STATUS;
					free(dtReadHead);
					tandPath->itemNumDtReadsList --;
					dtReadHead = dtRead;
				}
				tandPath->dtReadsList = NULL;

				if(tandPath->itemNumDtReadsList!=0)
				{
					printf("line=%d, In %s(), invalid itemNumDtReadsList=%d, error!\n", __LINE__, __func__, tandPath->itemNumDtReadsList);
					return FAILED;
				}

				// delete the nodes
				tandPathTmp = tandPath->next;
				if(preTandPath)
					preTandPath->next = tandPathTmp;
				else
					*tandemPathList = tandPathTmp;
				free(tandPath);
				(*itemNumTandemPathList) --;
				tandPath = tandPathTmp;
			}else
			{
				preTandPath = tandPath;
				tandPath = tandPath->next;
			}
		}
	}
*/

	// remove the reads in decision table
	if((*shortFragSizeRemovedFlag)==YES)
	{
		// delete the failed reads
		if(removeFinishedReadsFromDecisionTable(decisionTable, readsNumDecisionTable)==FAILED)
		{
			printf("line=%d, In %s(), cannot resort decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Remove reads in second path of two overlapped tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeAverFragSizeTandPath(double *averFragSize, tandemPathItem_t *tandPath, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable)
{
	uint64_t readID, readID_paired, pairedNum;
	int32_t validReadOrient, validReadOrient_paired;
	PERead_t *pPERead;
	assemblingreadtype *dtRead, *dtReadPaired;
	dtReadInTandemPath_t *dtReadInTandPath;

	*averFragSize = 0;
	pairedNum = 0;
	dtReadInTandPath = tandPath->dtReadsList;
	while(dtReadInTandPath)
	{
		dtRead = decisionTable + dtReadInTandPath->dtRow;
		readID = dtRead->rid;

		// generate the paired readID
		if(readID%2==1)
		{ // odd --> even
			readID_paired = readID + 1;
		}else
		{ // even --> odd
			readID_paired = readID - 1;
		}

		if(getReadFromPEHashtable(&pPERead, readID_paired)==FAILED)
		{
			printf("In %s(), cannot get the paired read %lu, error!\n", __func__, readID_paired);
			return FAILED;
		}

		if(pPERead)
		{
			*averFragSize += dtRead->firstContigPos + dtRead->seqlen - pPERead->cpos;
			pairedNum ++;
		}else if(shortInsertFlag==YES)
		{
			// get the exist read
			if(getExistReadInDT(&dtReadPaired, readID_paired, decisionTable, dtRowHashtable)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the exist read %lu in decision table, error!\n", __LINE__, __func__, readID_paired);
				return FAILED;
			}

			if(dtReadPaired && dtReadPaired->orientation==ORIENTATION_PLUS && dtReadPaired->unmatchBaseNum==0)
			{
				*averFragSize += dtRead->firstContigPos + dtRead->seqlen - dtReadPaired->firstContigPos;
				pairedNum ++;
			}
		}

		dtReadInTandPath = dtReadInTandPath->next;
	}

	if(pairedNum>0)
		*averFragSize /= pairedNum;

	return SUCCESSFUL;
}

/**
 * Get the contig match flag for tandem path.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getContigMatchFlagTandPath(int32_t *matchWithContigFlag, int32_t *mismatchNum, tandemPathItem_t *tandPath, contigPath_t *contigPath)
{
	int32_t i, pathLen, contigseqLen, shareLen, startRowPathseq, startRowContigseq;
	char *pathseq, *contigseq;

	pathLen = tandPath->contigtailPathPos + 1;
	contigseqLen = contigPath->contigtailSeqLen;
	if(pathLen<=contigseqLen)
	{
		shareLen = pathLen;
		startRowPathseq = 0;
		startRowContigseq = contigseqLen - shareLen;
	}else
	{
		shareLen = contigseqLen;
		startRowPathseq = pathLen - shareLen;
		startRowContigseq = 0;
	}

	pathseq = tandPath->tandemPathStr + startRowPathseq;
	contigseq = contigPath->contigtailSeq + startRowContigseq;

	*mismatchNum = 0;
	for(i=0; i<shareLen; i++)
		if(pathseq[i]!=contigseq[i])
			(*mismatchNum) ++;

	if((*mismatchNum)<=1)
		*matchWithContigFlag = YES;
	else
		*matchWithContigFlag = NO;

	return SUCCESSFUL;
}

/**
 * Remove reads in second path of two overlapped tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeOverlappedTandPath(int32_t *overlapTandPathRemovedFlag, tandemPathItem_t *tandemPathList, int32_t *itemNumTandemPathList, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable)
{
	int32_t overlapType;   // 0: No overlap, 1: path1 => path2; 2: path2 => path1
	int32_t overlapSize, pathNodeDeleteFlag;
	tandemPathItem_t *targetTandPath, *preTandPath, *tandPath;

	*overlapTandPathRemovedFlag = NO;
	targetTandPath = tandemPathList;
	while(targetTandPath)
	{
		preTandPath = targetTandPath;
		tandPath = preTandPath->next;
		while(tandPath)
		{
			// get the overlap size
			if(computeOverlapSizeTandPath(&overlapSize, &overlapType, targetTandPath, tandPath)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the overlap type of two tandem paths, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(overlapSize>0)
			{
				*overlapTandPathRemovedFlag = YES;

				if(overlapType==2)
				{
					if(exchangeNodeInfoTandPath(targetTandPath, tandPath)==FAILED)
					{
						printf("line=%d, In %s(), cannot exchange two tandem path nodes, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

				// remove the error overlapped reads
				if(removeErrOverlappedReadsInTandPath(&pathNodeDeleteFlag, targetTandPath, tandPath, overlapSize, decisionTable)==FAILED)
				{
					printf("line=%d, In %s(), cannot remove the overlap type of two tandem paths, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// delete the nodes
				if(pathNodeDeleteFlag==YES)
				{
					preTandPath->next = tandPath->next;
					free(tandPath);
					(*itemNumTandemPathList) --;
				}else
				{
					preTandPath = tandPath;
				}
			}else
			{
				preTandPath = tandPath;
			}

			tandPath = preTandPath->next;
		}

		targetTandPath = targetTandPath->next;
	}

	if((*overlapTandPathRemovedFlag)==YES)
	{
		// delete the failed reads
		if(removeFinishedReadsFromDecisionTable(decisionTable, readsNumDecisionTable)==FAILED)
		{
			printf("line=%d, In %s(), cannot resort decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Compute the overlap size of two tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeOverlapSizeTandPath(int32_t *overlapSize, int32_t *overlapType, tandemPathItem_t *tandPath1, tandemPathItem_t *tandPath2)
{
	int32_t i, pathlen1, pathlen2, maxOverlapSize, overSize, startRow1, startRow2, mismatchNum, matchFlag;
	char *pathseq1, *pathseq2;

	pathseq1 = tandPath1->tandemPathStr;
	pathseq2 = tandPath2->tandemPathStr;
	pathlen1 = tandPath1->tandemPathLen;
	pathlen2 = tandPath2->tandemPathLen;

	if(pathlen1>pathlen2)
		maxOverlapSize = pathlen2;
	else
		maxOverlapSize = pathlen1;

	// detect the overlap: seq1 => seq2
	for(overSize=maxOverlapSize; overSize>=MIN_OVERLAP_SIZE_TANDPATH; overSize--)
	{
		startRow1 = pathlen1 - overSize;
		startRow2 = 0;

		matchFlag = YES;
		mismatchNum = 0;
		for(i=0; i<overSize; i++)
		{
			if(pathseq1[startRow1+i]!=pathseq2[startRow2+i])
			{
				mismatchNum ++;
				if(mismatchNum>1)
				{
					matchFlag = NO;
					break;
				}
			}
		}

		if(matchFlag==YES)
		{
			*overlapSize = overSize;
			*overlapType = 1;
			break;
		}
	}

	// detect the overlap: seq2 => seq1
	if(matchFlag==NO)
	{
		for(overSize=maxOverlapSize; overSize>=MIN_OVERLAP_SIZE_TANDPATH; overSize--)
		{
			startRow1 = 0;
			startRow2 = pathlen2 - overSize;

			matchFlag = YES;
			mismatchNum = 0;
			for(i=0; i<overSize; i++)
			{
				if(pathseq1[startRow1+i]!=pathseq2[startRow2+i])
				{
					mismatchNum ++;
					if(mismatchNum>1)
					{
						matchFlag = NO;
						break;
					}
				}
			}

			if(matchFlag==YES)
			{
				*overlapSize = overSize;
				*overlapType = 2;
				break;
			}
		}
	}

	if(matchFlag==NO)
	{
		*overlapSize = -1;
		*overlapType = 0;
	}

	return SUCCESSFUL;
}

/**
 * Compute the overlap size of two tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeShiftSizeTandPathByShiftOp(int32_t *mismatchNum, int32_t *shiftType, int32_t *shiftSize, tandemPathItem_t *tandPath1, tandemPathItem_t *tandPath2)
{
	int32_t i, j, k, pathLenTmp1, pathLenTmp2, shareLen, mismatchNumTmp, shiftSizeFactor, shiftSizeTmp, shiftTypeTmp, startRow1, startRow2, matchFlag;
	char *pathseqTmp1, *pathseqTmp2;

	*mismatchNum = INT_MAX;
	*shiftType = -1;
	*shiftSize = 0;


	shiftTypeTmp = -1;
	matchFlag = NO;
	for(i=1; i<=2; i++)
	{
		if(i==1)
			shiftSizeFactor = -1;
		else
			shiftSizeFactor = 1;

		for(j=1; j<=10; j++)
		{
			shiftSizeTmp = shiftSizeFactor * j;

			if(tandPath1->contigtailPathPos<=tandPath2->contigtailPathPos+shiftSizeTmp)
			{
				startRow1 = 0;
				startRow2 = (tandPath2->contigtailPathPos + shiftSizeTmp) - tandPath1->contigtailPathPos;
				shiftTypeTmp = 2;
			}else
			{
				startRow1 = tandPath1->contigtailPathPos - (tandPath2->contigtailPathPos + shiftSizeTmp);
				startRow2 = 0;
				shiftTypeTmp = 1;
			}

			pathseqTmp1 = tandPath1->tandemPathStr + startRow1;
			pathseqTmp2 = tandPath2->tandemPathStr + startRow2;
			pathLenTmp1 = tandPath1->tandemPathLen - startRow1;
			pathLenTmp2 = tandPath2->tandemPathLen - startRow2;
			if(pathLenTmp1<pathLenTmp2)
				shareLen = pathLenTmp1;
			else
				shareLen = pathLenTmp2;

			if(shareLen>0)
			{
				mismatchNumTmp = 0;
				for(k=0; k<shareLen; k++)
				{
					if(pathseqTmp1[startRow1+k]!=pathseqTmp2[startRow2+k])
					{
						mismatchNumTmp ++;
						if(mismatchNumTmp>1)
							break;
					}
				}
			}else
			{
				mismatchNumTmp = INT_MAX;
			}

			if(mismatchNumTmp<=1)
			{
				matchFlag = YES;
				*mismatchNum = mismatchNumTmp;
				*shiftType = shiftTypeTmp;
				*shiftSize = shiftSizeTmp;
				break;
			}
		}

		if(matchFlag==YES)
			break;
	}

	return SUCCESSFUL;
}

/**
 * Remove the error overlapped reads of two tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeErrOverlappedReadsInTandPath(int32_t *pathNodeDeleteFlag, tandemPathItem_t *targetTandPath, tandemPathItem_t *tandPath, int32_t overlapSize, assemblingreadtype *decisionTable)
{
	int32_t baseShiftNum;
	dtReadInTandemPath_t *dtRead, *dtReadHead, *preDtRead, *nextDtRead;
	assemblingreadtype *pRead;

	if(overlapSize>0)
	{
		baseShiftNum = tandPath->tandemPathLen - overlapSize - (tandPath->endContigPos - targetTandPath->endContigPos);
		tandPath->startContigPos += baseShiftNum;
		tandPath->endContigPos += baseShiftNum;
		tandPath->contigtailPathPos -= baseShiftNum;

		// delete the reads that do not overlapped with the contigs
		if(tandPath->contigtailPathPos>=kmerSize-1)
		{ // select the reads to delete
			*pathNodeDeleteFlag = NO;

			// update the the read status in decision table
			preDtRead = nextDtRead = NULL;
			dtRead = tandPath->dtReadsList;
			while(dtRead)
			{
				nextDtRead = dtRead->next;
				dtRead->endContigPos += baseShiftNum;
				dtRead->startContigPos += baseShiftNum;
				dtRead->contigtailReadPos -= baseShiftNum;
				if(dtRead->contigtailReadPos>=kmerSize-1)
				{ // update the read information in decision table
					if(dtRead->dtRow>=0)
					{
						pRead = decisionTable + dtRead->dtRow;
						if(pRead->orientation==ORIENTATION_PLUS)
						{
							pRead->firstBasePos = 0;
							pRead->basePos = tandPath->contigtailPos - dtRead->startContigPos;
						}else
						{
							pRead->firstBasePos = pRead->seqlen - 1;
							pRead->basePos = pRead->seqlen - 1 - (tandPath->contigtailPos - dtRead->startContigPos);
						}
						pRead->lastMatchedBasePos = pRead->basePos;
						pRead->firstContigPos = dtRead->startContigPos;
						pRead->matchBaseNum = tandPath->contigtailPos - pRead->firstContigPos + 1;
						pRead->unmatchBaseNum = 0;

						pRead->successiveAppearBases = pRead->matchBaseNum;
						pRead->successiveUnappearBases = 0;
						pRead->unappearBlocksNum = 0;
					}

					preDtRead = dtRead;
					dtRead = nextDtRead;
				}else
				{ // delete the read
					if(dtRead->dtRow>=0)
						decisionTable[dtRead->dtRow].status = FAILED_STATUS;

					if(preDtRead==NULL)
						tandPath->dtReadsList = nextDtRead;
					else
						preDtRead->next = nextDtRead;
					free(dtRead);
					dtRead = nextDtRead;

					tandPath->itemNumDtReadsList --;
				}
			}
		}else
		{ // delete all the reads
			*pathNodeDeleteFlag = YES;

			// update the the read status in decision table, and delete the reads
			dtReadHead = tandPath->dtReadsList;
			while(dtReadHead)
			{
				dtRead = dtReadHead->next;
				if(dtReadHead->dtRow>=0)
					decisionTable[dtReadHead->dtRow].status = FAILED_STATUS;
				free(dtReadHead);
				tandPath->itemNumDtReadsList --;
				dtReadHead = dtRead;
			}
			tandPath->dtReadsList = NULL;

			if(tandPath->itemNumDtReadsList!=0)
			{
				printf("line=%d, In %s(), invalid itemNumDtReadsList=%d, error!\n", __LINE__, __func__, tandPath->itemNumDtReadsList);
				return FAILED;
			}
		}
	}else
	{
		printf("line=%d, In %s(), invalid overlap size for tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Exchange the information of two tandem path nodes.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short exchangeNodeInfoTandPath(tandemPathItem_t *tandPath1, tandemPathItem_t *tandPath2)
{
	tandemPathItem_t tandPathTmp;

	if(memcpy(&tandPathTmp, tandPath1, sizeof(tandemPathItem_t))==NULL)
	{
		printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(memcpy(tandPath1, tandPath2, sizeof(tandemPathItem_t))==NULL)
	{
		printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(memcpy(tandPath2, &tandPathTmp, sizeof(tandemPathItem_t))==NULL)
	{
		printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	tandPath2->next = tandPath1->next;
	tandPath1->next = tandPathTmp.next;

	return SUCCESSFUL;
}

/**
 * Determine the navigation by tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short decideNaviByTandPath(int32_t *naviTandPath, int32_t *maxBaseIndexAfterTandPath, int32_t *incorrectBaseNumTandPath, int32_t *newCandBaseNumAfterTandPath, int32_t *occNumArray, int32_t *occIndexArray, tandemPathItem_t *tandemPathList, int32_t itemNumTandemPathList)
{
	int32_t *baseNumArray, maxRowNumBaseNumArray, colsNum, maxIndex;

	// get the candidate paths
	if(getCandPathsTandPath(candPath, tandemPathList, itemNumTandemPathList)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the candPath from tandem paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(candPath->itemNumCandPathItemArray>=2)
	{
		// initialize the memory
		colsNum = 8;
		maxRowNumBaseNumArray = candPath->maxPathLen;
		baseNumArray = (int32_t *)calloc(maxRowNumBaseNumArray*colsNum, sizeof(int32_t));
		if(baseNumArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute the base number array
		if(computeBaseNumArrayCandPath(baseNumArray, maxRowNumBaseNumArray, colsNum, candPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the base number for candPath, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// compute the number of incorrect bases
		if(computeIncorrectBaseNumCandPath(incorrectBaseNumTandPath, baseNumArray, maxRowNumBaseNumArray, colsNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the number of incorrect bases for candidate paths, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// decide navigation by incorrect base number
		if(decideNaviTandPath(naviTandPath, &maxIndex, newCandBaseNumAfterTandPath, *incorrectBaseNumTandPath, baseNumArray, maxRowNumBaseNumArray, colsNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot navigation from candidate paths, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if((*naviTandPath)==NAVI_FAILED) // 2014-01-29
		{
			if(candPath->itemNumCandPathItemArray==2)
			{
				if(confirmNaviByShiftOpCandPath(naviTandPath, &maxIndex, incorrectBaseNumTandPath, candPath, contigPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot navigation from shifted candidate path, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}

		free(baseNumArray);
	}
	else if(candPath->itemNumCandPathItemArray==1)
	{
		*naviTandPath = NAVI_SUCCESS;
		*incorrectBaseNumTandPath = 0;
		*newCandBaseNumAfterTandPath = 1;
		switch(candPath->candPathItemArray[0].candPathStr[0])
		{
			case 'A': maxIndex = 0; break;
			case 'C': maxIndex = 1; break;
			case 'G': maxIndex = 2; break;
			case 'T': maxIndex = 3; break;
			default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, candPath->candPathItemArray[0].candPathStr[0]); return FAILED;
		}
	}
	else
	{
		*naviTandPath = NAVI_FAILED;
		*incorrectBaseNumTandPath = INT_MAX;
		*newCandBaseNumAfterTandPath = 0;
	}

	if((*naviTandPath)==NAVI_SUCCESS)
	{
		kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] >> 2) << 2) | maxIndex;

		kmers[0] = getKmer(kmerSeqIntAssembly, deBruijnGraph);
		kmers[1] = getReverseKmer(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, deBruijnGraph);

		*maxBaseIndexAfterTandPath = maxIndex;
	}else
	{
		*maxBaseIndexAfterTandPath = -1;
	}


#if(CANDPATH_OUTPUT==YES)
	// output the paths
	outputCandPath(candPath);
	printf("=*=*=*=* localContigID=%ld, contigNodesNum=%ld, naviTandPath=%d, incorrectBaseNumTandPath=%d\n", localContigID, itemNumContigArr, *naviTandPath, *incorrectBaseNumTandPath);
#endif

	return SUCCESSFUL;
}

/**
 * Get the candidate paths from tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getCandPathsTandPath(candPath_t *candPath, tandemPathItem_t *tandemPathList, int32_t itemNumTandemPathList)
{
	// get the candPath from tandem paths
	if(getCandPathsFromTandPathList(candPath, tandemPathList, itemNumTandemPathList)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the candidate paths from tandem paths, error!\n", __LINE__, __func__);
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
 * Get the candidate paths from tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getCandPathsFromTandPathList(candPath_t *candPath, tandemPathItem_t *tandemPathList, int32_t itemNumTandemPathList)
{
	int32_t i, j, k, pathLen, targetLen, shareLen, perfectMatchFlag, matchRow;
	char *pathseq, *targetseq;
	tandemPathItem_t *tandPath;

	candPath->itemNumCandPathItemArray = 0;
	candPath->maxPathLen = 0;

/*
	tandPath = tandemPathList;
	while(tandPath)
	{
		pathLen = tandPath->tandemPathLen - tandPath->contigtailPathPos - 1;
		if(pathLen>0)
		{
			strcpy(candPath->candPathItemArray[candPath->itemNumCandPathItemArray].candPathStr, tandPath->tandemPathStr+tandPath->contigtailPathPos+1);
			candPath->candPathItemArray[candPath->itemNumCandPathItemArray].pathLen = pathLen;
			candPath->candPathItemArray[candPath->itemNumCandPathItemArray].supportReadsNum = tandPath->itemNumDtReadsList;
			if(candPath->maxPathLen<pathLen)
				candPath->maxPathLen = pathLen;
			candPath->itemNumCandPathItemArray ++;

			if(candPath->itemNumCandPathItemArray>=candPath->maxItemNumCandPathItemArray)
				break;
		}

		tandPath = tandPath->next;
	}

	// merge the same candidate paths
	for(i=0; i<candPath->itemNumCandPathItemArray; i++)
	{
		targetseq = candPath->candPathItemArray[i].candPathStr;
		targetLen = candPath->candPathItemArray[i].pathLen;
		for(j=candPath->itemNumCandPathItemArray-1; j>i; j--)
		{
			pathseq = candPath->candPathItemArray[j].candPathStr;
			pathLen = candPath->candPathItemArray[j].pathLen;
			if(targetLen<pathLen)
				shareLen = targetLen;
			else
				shareLen = pathLen;

			if(shareLen>0)
			{
				perfectMatchFlag = YES;
				for(k=0; k<shareLen; k++)
				{
					if(targetseq[k]!=pathseq[k])
					{
						perfectMatchFlag = NO;
						break;
					}
				}
			}else
			{
				perfectMatchFlag = NO;
			}

			if(perfectMatchFlag==YES)
			{ // merge the two paths
				if(targetLen<pathLen)
				{
					strcat(candPath->candPathItemArray[i].candPathStr, pathseq+targetLen);
					candPath->candPathItemArray[i].pathLen = pathLen;
				}
				candPath->candPathItemArray[i].supportReadsNum += candPath->candPathItemArray[j].supportReadsNum;

				k = j + 1;
				while(k<candPath->itemNumCandPathItemArray)
				{
					strcpy(candPath->candPathItemArray[k-1].candPathStr, candPath->candPathItemArray[k].candPathStr);
					candPath->candPathItemArray[k-1].pathLen = candPath->candPathItemArray[k].pathLen;
					candPath->candPathItemArray[k-1].supportReadsNum = candPath->candPathItemArray[k].supportReadsNum;
					k ++;
				}

				candPath->itemNumCandPathItemArray --;
			}
		}
	}
*/


	tandPath = tandemPathList;
	while(tandPath)
	{
		pathseq = tandPath->tandemPathStr + tandPath->contigtailPathPos + 1;
		pathLen = tandPath->tandemPathLen - tandPath->contigtailPathPos - 1;
		if(pathLen>0)
		{
			// get the matched candPathItem
			if(getMatchRowCandPathItem(&matchRow, pathseq, pathLen, candPath)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the matched row in candPath, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(matchRow>=0)
			{ // add the read if it matches
				if(addReadseqToCandPathItem(pathseq, pathLen, tandPath->itemNumDtReadsList, matchRow, candPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot add the read to candPath, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}else
			{ // generate a new candPath item
				if(addNewCandPathItem(candPath, pathseq, pathLen, tandPath->itemNumDtReadsList)==FAILED)
				{
					printf("line=%d, In %s(), cannot add the contig path item, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}

		tandPath = tandPath->next;
	}


	return SUCCESSFUL;
}

/**
 * Decide the navigation by candidate paths from tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short decideNaviTandPath(int32_t *naviTandPath, int32_t *maxIndex, int32_t *newCandBaseNumAfterTandPath, int32_t incorrectNum, int32_t *baseNumArray, int32_t maxRowNumBaseNumArray, int32_t colsNum)
{
	int32_t i, firstBaseRow, maxBaseIndex, maxValue, secondBaseIndex, secondValue, sum;
	int32_t lowNum;

	// 2014-01-28
	for(i=0; i<maxRowNumBaseNumArray; i++)
	{
		maxBaseIndex = baseNumArray[i*colsNum + 6];
		secondBaseIndex = baseNumArray[i*colsNum + 7];

		if(maxBaseIndex<4)
		{
			firstBaseRow = i;
			break;
		}
	}

	*newCandBaseNumAfterTandPath = 0;
	for(i=0; i<4; i++)
		if(baseNumArray[firstBaseRow*colsNum+i]>0)
			(*newCandBaseNumAfterTandPath) ++;


	if(incorrectNum>MAX_INCOR_BASE_NUM_TANDPATH)
	{
		*maxIndex = -1;
		*naviTandPath = NAVI_FAILED;
	}else
	{
		// skip the gaps
		for(i=0; i<maxRowNumBaseNumArray; i++)
		{
			maxBaseIndex = baseNumArray[i*colsNum + 6];
			secondBaseIndex = baseNumArray[i*colsNum + 7];

			if(maxBaseIndex<4)
			{
				firstBaseRow = i;
				break;
			}
		}

		maxValue = baseNumArray[firstBaseRow*colsNum + maxBaseIndex];
		secondValue = baseNumArray[firstBaseRow*colsNum + secondBaseIndex];
		sum = baseNumArray[firstBaseRow*colsNum + 5];

		if((double)secondValue/maxValue>0.8 && incorrectNum>1 && sum>10)
		{
			*maxIndex = -1;
			*naviTandPath = NAVI_FAILED;
			//printf("line=%d, In %s(), naviTandPath=stop!\n", __LINE__, __func__);
		}
		else if(((double)maxValue/maxValue>MAX_OCC_RATIO_TANDPATH) && incorrectNum>2)
		{
			*maxIndex = -1;
			*naviTandPath = NAVI_FAILED;
			//printf("line=%d, In %s(), naviTandPath=stop!\n", __LINE__, __func__);
		}
		else
		{
			*maxIndex = maxBaseIndex;
			*naviTandPath = NAVI_SUCCESS;
			//printf("line=%d, In %s(), naviTandPath=continue!\n", __LINE__, __func__);
		}
	}

	return SUCCESSFUL;
}


/**
 * Decide the navigation by candidate paths from tandem paths.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short checkTandPath(tandemPathItem_t *tandemPathList, int32_t itemNumTandemPathList)
{
	int32_t i, j, pathLen;
	tandemPathItem_t *tandPath;

	tandPath = tandemPathList;
	while(tandPath)
	{
		pathLen = strlen(tandPath->tandemPathStr);
		if(pathLen==0 || tandPath->tandemPathLen==0)
		{
			printf("line=%d, In %s(), invalid pathLen=%d, tandemPathLen=%d, tandemPathStr=%s, error!\n", __LINE__, __func__, pathLen, tandPath->tandemPathLen, tandPath->tandemPathStr);
			return FAILED;
		}else if(pathLen!=tandPath->tandemPathLen)
		{
			printf("line=%d, In %s(), invalid pathLen=%d, tandemPathLen=%d, tandemPathStr=%s, error!\n", __LINE__, __func__, pathLen, tandPath->tandemPathLen, tandPath->tandemPathStr);
			return FAILED;
		}

		tandPath = tandPath->next;
	}

	return SUCCESSFUL;
}


/**
 * Remove short overlap size tandem paths with contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeShortOverlappedTandPath(int32_t *shortFragSizeRemovedFlag, tandemPathItem_t **tandemPathList, int32_t *itemNumTandemPathList, assemblingreadtype *decisionTable, int32_t *readsNumDecisionTable, dtRowIndex_t **dtRowHashtable)
{
	int32_t maxSupportReadsNum, mismatchNum, maxOverlapSizeWithContig, validTandPathFlag, removedFlag;
	tandemPathItem_t *preTandPath, *tandPath, *tandPathTmp;
	dtReadInTandemPath_t *dtRead, *dtReadHead;

	if((*itemNumTandemPathList)>=2)
	{
		validTandPathFlag = NO;
		maxSupportReadsNum = 0;
		tandPath = *tandemPathList;
		while(tandPath)
		{
			if(tandPath->validFlag==YES && tandPath->shortFragFlag==NO && tandPath->matchWithContigFlag==YES && tandPath->contigtailPathPos>0.9*readLen)
			{
				validTandPathFlag = YES;
				maxSupportReadsNum = tandPath->itemNumDtReadsList;
				break;
			}
			tandPath = tandPath->next;
		}

		if(validTandPathFlag==YES)
		{
			tandPath = *tandemPathList;
			while(tandPath)
			{
				if(tandPath->contigtailPathPos<0.6*readLen && tandPath->itemNumDtReadsList<0.6*maxSupportReadsNum)
					tandPath->validFlag = NO;
				tandPath = tandPath->next;
			}
		}
	}

	removedFlag = NO;

	// check the average contigtailPathPos
	if((*itemNumTandemPathList)>=2)
	{
		preTandPath = NULL;
		tandPath = *tandemPathList;
		while(tandPath)
		{
			if(tandPath->validFlag==NO)
			{ // remove condition
				removedFlag = YES;

				// remove the error tandem path and its reads
				dtReadHead = tandPath->dtReadsList;
				while(dtReadHead)
				{
					dtRead = dtReadHead->next;
					if(dtReadHead->dtRow>=0)
						decisionTable[dtReadHead->dtRow].status = FAILED_STATUS;
					free(dtReadHead);
					tandPath->itemNumDtReadsList --;
					dtReadHead = dtRead;
				}
				tandPath->dtReadsList = NULL;

				if(tandPath->itemNumDtReadsList!=0)
				{
					printf("line=%d, In %s(), invalid itemNumDtReadsList=%d, error!\n", __LINE__, __func__, tandPath->itemNumDtReadsList);
					return FAILED;
				}

				// delete the nodes
				tandPathTmp = tandPath->next;
				if(preTandPath)
					preTandPath->next = tandPathTmp;
				else
					*tandemPathList = tandPathTmp;
				free(tandPath);
				(*itemNumTandemPathList) --;
				tandPath = tandPathTmp;
			}else
			{
				preTandPath = tandPath;
				tandPath = tandPath->next;
			}
		}
	}

	// remove the reads in decision table
	if(removedFlag==YES)
	{
		*shortFragSizeRemovedFlag = removedFlag;

		// delete the failed reads
		if(removeFinishedReadsFromDecisionTable(decisionTable, readsNumDecisionTable)==FAILED)
		{
			printf("line=%d, In %s(), cannot resort decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Check the valid flag of new candTandPath.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short isValidNewCandTandPath(int32_t *validFlag, char *newCandTandPathseq, int32_t newCandTandPathLen, contigPath_t *contigPath)
{
	char *preCandPathseq, *maxPathseq;
	int32_t i, preCandPathLen, maxPathLen, shareLen, mismatchNum1, mismatchNum2, shiftSize, shiftType;

	if(contigPath->validCandPathTandPathFlagPE==YES && contigPath->naviPathItem==contigPath->maxPathItem && contigPath->maxPathItem->contigPathLen-contigPath->startRowNewBase>5)
	{
		preCandPathseq = contigPath->candPathseqTandPathPE + contigPath->startRowCandPathTandPathPE;
		preCandPathLen = contigPath->candPathLenTandPathPE - contigPath->startRowCandPathTandPathPE;
		maxPathseq = contigPath->maxPathItem->contigPathStr + contigPath->startRowNewBase;
		maxPathLen = contigPath->maxPathItem->contigPathLen - contigPath->startRowNewBase;
		if(preCandPathLen>maxPathLen)
			shareLen = maxPathLen;
		else
			shareLen = preCandPathLen;

		mismatchNum1 = 0;
		for(i=0; i<shareLen; i++)
		{
			if(preCandPathseq[i]!=maxPathseq[i])
			{
				mismatchNum1 ++;
				if(mismatchNum1>1)
					break;
			}
		}


		if(newCandTandPathLen>maxPathLen)
			shareLen = maxPathLen;
		else
			shareLen = newCandTandPathLen;

		mismatchNum2 = 0;
		for(i=0; i<shareLen; i++)
		{
			if(newCandTandPathseq[i]!=maxPathseq[i])
			{
				mismatchNum2 ++;
				if(mismatchNum2>1)
					break;
			}
		}

		if(mismatchNum1>=mismatchNum2 && mismatchNum2<=1)
			*validFlag = YES;
		else if(mismatchNum1<mismatchNum2 && mismatchNum1<=1)
			*validFlag = NO;
		else if(mismatchNum1>1 && mismatchNum2>1)
			*validFlag = YES;
		else
			*validFlag = NO;

/*
		if(getMismatchNumByShiftOp(&mismatchNum1, &shiftSize, &shiftType, preCandPathseq, preCandPathLen, maxPathseq, maxPathLen)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the mismatch size by shift operation, error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(shiftType==-1)
			mismatchNum1 = INT_MAX;

		if(getMismatchNumByShiftOp(&mismatchNum2, &shiftSize, &shiftType, newCandTandPathseq, newCandTandPathLen, maxPathseq, maxPathLen)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the mismatch size by shift operation, error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(shiftType==-1)
			mismatchNum2 = INT_MAX;

		if(mismatchNum1>=mismatchNum2 && mismatchNum1!=INT_MAX)
			*validFlag = YES;
		else if(mismatchNum1==INT_MAX && mismatchNum2==INT_MAX)
			*validFlag = YES;
		else
			*validFlag = NO;
*/
	}else
	{
		*validFlag = YES;
	}

	return SUCCESSFUL;
}

/**
 * Check the valid flag of new candTandPath.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxSecRowCandTandPath(int32_t *maxRow, int32_t *secRow, candPath_t *candPath)
{
	int32_t i, j, maxValue, secValue;

	maxValue = secValue = 0;
	*maxRow = *secRow = -1;
	for(i=0; i<candPath->itemNumCandPathItemArray; i++)
	{
		if(candPath->candPathItemArray[i].supportReadsNum>maxValue)
		{
			secValue = maxValue;
			*secRow = *maxRow;
			maxValue = candPath->candPathItemArray[i].supportReadsNum;
			*maxRow = i;
		}else if(candPath->candPathItemArray[i].supportReadsNum>secValue)
		{
			secValue = candPath->candPathItemArray[i].supportReadsNum;
			*secRow = i;
		}
	}

	return SUCCESSFUL;
}
