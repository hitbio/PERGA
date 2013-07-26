/*
 * malign.c
 *
 *  Created on: Apr 4, 2013
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Use the multi-align to generate the extension by paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short decideByMalignPE(int32_t *naviMalignPE, int32_t *incorrectNumPE, int32_t maxPEIndex, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, graphtype *graph)
{
	int32_t totalReadsNum, itemNumMalignArray, maxAlignSeqLen, maxIndex;
	char consensusSeq[MAX_SEQ_LEN_MALIGN+1];
	int32_t *shiftPosArray, consensusBaseNum[MAX_SEQ_LEN_MALIGN+1][8];  // [0-3]: --A, C, G, T, -; [5]-- A+C+G+T; [6]--maxIndex, [7]-- secIndex


	// get the number of total reads for multi-align
	if(getTotalReadsNumMalignPE(&totalReadsNum, readsNumDecisionTable, kmerSeqIntAssembly, graph, decisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot get total the number of new reads for multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate the shift pos array
	shiftPosArray = (int32_t *) calloc(totalReadsNum, sizeof(int32_t));
	if(shiftPosArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// add the read sequences for multi-align
	if(addReadseqMalignPE(readseqMalignArr, readIDMalignArr, &itemNumMalignArray, kmerSeqIntAssembly, graph, decisionTable, readsNumDecisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the reads for multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(itemNumMalignArray>totalReadsNum)
	{
		printf("line=%d, In %s(), itemNumMalignArray=%d > totalReadsNum=%d, error!\n", __LINE__, __func__, itemNumMalignArray, totalReadsNum);
		return FAILED;
	}

	// adjust the read sequences
	if(adjustReadseqMalign(readseqMalignArr, shiftPosArray, itemNumMalignArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust the reads for multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}

//	if(generateMSAMalign(alignResultsMalign, readseqMalignArr, itemNumMalignArray)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot adjust the reads for multi-align, error!\n", __LINE__, __func__);
//		return FAILED;
//	}

	// compute the consensus base numbers for multi-align.
	//if(computeConsensusBaseNumMalign((int32_t**)consensusBaseNum, &maxAlignSeqLen, MAX_SEQ_LEN_MALIGN, 8, alignResultsMalign, itemNumMalignArray)==FAILED)
	if(computeConsensusBaseNumMalign((int32_t**)consensusBaseNum, &maxAlignSeqLen, MAX_SEQ_LEN_MALIGN, 8, readseqMalignArr, shiftPosArray, itemNumMalignArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust the reads for multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// compute the number of incorrect bases
	if(computeIncorrectBaseNumMalign(incorrectNumPE, (int32_t**)consensusBaseNum, 8, maxAlignSeqLen, consensusSeq, itemNumMalignArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the number of incorrect bases for multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// decide navigation by incorrect base number
	if(decideNaviMalignPE(naviMalignPE, &maxIndex, maxPEIndex, *incorrectNumPE, (int32_t**)consensusBaseNum, 8, maxAlignSeqLen, consensusSeq)==FAILED)
	{
		printf("line=%d, In %s(), cannot navigation from multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(*naviMalignPE==NAVI_SUCCESS)
	{
		kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] >> 2) << 2) | maxIndex;

		kmers[0] = getKmer(kmerSeqIntAssembly, deBruijnGraph);
		kmers[1] = getReverseKmer(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, deBruijnGraph);
	}


#if(MALIGN_OUTPUT==YES)
	int32_t i, j;
	if(*naviMalignPE==NAVI_SUCCESS)
		printf("*********** line=%d, itemNumContigArr=%ld, incorrectNumPE=%d, continue!\n", __LINE__, itemNumContigArr, *incorrectNumPE);
	else
		printf("*********** line=%d, itemNumContigArr=%ld, incorrectNumPE=%d, stop!\n", __LINE__, itemNumContigArr, *incorrectNumPE);

	printf("occsNumPE: (%d,%d,%d,%d), readsNumRatio=%.2f\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3], readsNumRatio);
	for(i=0; i<itemNumMalignArray; i++)
	{
		printf("%ld\t", readIDMalignArr[i]);
		for(j=0; j<shiftPosArray[i]; j++)
			printf("-");
		printf("%s\n", readseqMalignArr[i]);
	}
//	printf("align result:\n");
//	for(i=0; i<itemNumMalignArray; i++)
//	{
//		printf("%ld\t", readIDMalignArr[i]);
////			for(j=0; j<shiftPosArray[i]; j++)
////				printf("-");
//		printf("%s\n", alignResultsMalign[i]);
//	}
	for(i=0; i<maxAlignSeqLen; i++)
		printf("(%d,%d,%d,%d,%d), total=%d\n", consensusBaseNum[i][0], consensusBaseNum[i][1], consensusBaseNum[i][2], consensusBaseNum[i][3], consensusBaseNum[i][4], consensusBaseNum[i][5]);
#endif

	// release memory
	free(shiftPosArray);

	return SUCCESSFUL;
}

/**
 * Use the multi-align to generate the extension by single ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short decideByMalignSE(int32_t *naviMalignSE, int32_t incorrectNumPE, int32_t maxPEIndex, assemblingreadtype *decisionTable, int32_t readsNumDecisionTable, graphtype *graph)
{
	int32_t totalReadsNum, itemNumMalignArray, maxAlignSeqLen, maxIndex;
	char consensusSeq[MAX_SEQ_LEN_MALIGN+1];
	int32_t *shiftPosArray, consensusBaseNum[MAX_SEQ_LEN_MALIGN+1][8];  // [0-4]: --A, C, G, T, -; [5]-- A+C+G+T+-; [6]--maxIndex
	int32_t incorrectBaseNumSE;

	// get the number of total reads for multi-align
	if(getTotalReadsNumMalignSE(&totalReadsNum, readsNumDecisionTable, kmerSeqIntAssembly, graph, decisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot get total the number of new reads for multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate the shift pos array
	shiftPosArray = (int32_t *) calloc(totalReadsNum, sizeof(int32_t));
	if(shiftPosArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// add the read sequences for multi-align
	if(addReadseqMalignSE(readseqMalignArr, readIDMalignArr, &itemNumMalignArray, kmerSeqIntAssembly, graph, decisionTable, readsNumDecisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the reads for multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(itemNumMalignArray>totalReadsNum)
	{
		printf("line=%d, In %s(), itemNumMalignArray=%d > totalReadsNum=%d, error!\n", __LINE__, __func__, itemNumMalignArray, totalReadsNum);
		return FAILED;
	}

	// adjust the read sequences
	if(adjustReadseqMalign(readseqMalignArr, shiftPosArray, itemNumMalignArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust the reads for multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}

//	if(generateMSAMalign(alignResultsMalign, readseqMalignArr, itemNumMalignArray)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot adjust the reads for multi-align, error!\n", __LINE__, __func__);
//		return FAILED;
//	}

	// compute the consensus base numbers for multi-align.
	if(computeConsensusBaseNumMalign((int32_t**)consensusBaseNum, &maxAlignSeqLen, MAX_SEQ_LEN_MALIGN, 8, readseqMalignArr, shiftPosArray, itemNumMalignArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot adjust the reads for multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// compute the number of incorrect bases
	if(computeIncorrectBaseNumMalign(&incorrectBaseNumSE, (int32_t**)consensusBaseNum, 8, maxAlignSeqLen, consensusSeq, itemNumMalignArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the number of incorrect bases for multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// decide navigation by incorrect base number
	if(decideNaviMalignSE(naviMalignSE, incorrectBaseNumSE, &maxIndex, maxPEIndex, (int32_t**)consensusBaseNum, 8, maxAlignSeqLen, consensusSeq)==FAILED)
	{
		printf("line=%d, In %s(), cannot navigation from multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(*naviMalignSE==NAVI_SUCCESS)
	{
		kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] >> 2) << 2) | maxIndex;

		kmers[0] = getKmer(kmerSeqIntAssembly, deBruijnGraph);
		kmers[1] = getReverseKmer(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, deBruijnGraph);
	}

#if(MALIGN_OUTPUT==YES)
	int32_t i, j;
	if(*naviMalignSE==NAVI_SUCCESS)
		printf("*********** line=%d, itemNumContigArr=%ld, incorrectNumSE=%d, continue!\n", __LINE__, itemNumContigArr, incorrectBaseNumSE);
	else
		printf("*********** line=%d, itemNumContigArr=%ld, incorrectNumSE=%d, stop!\n", __LINE__, itemNumContigArr, incorrectBaseNumSE);

	printf("occsNumSE: (%d,%d,%d,%d), readsNumRatio=%.2f\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3], readsNumRatio);
	for(i=0; i<itemNumMalignArray; i++)
	{
		printf("%ld\t", readIDMalignArr[i]);
		for(j=0; j<shiftPosArray[i]; j++)
			printf("-");
		printf("%s\n", readseqMalignArr[i]);
	}
//	printf("align result:\n");
//	for(i=0; i<itemNumMalignArray; i++)
//	{
//		printf("%ld\t", readIDMalignArr[i]);
////			for(j=0; j<shiftPosArray[i]; j++)
////				printf("-");
////			printf("%s\n", readseqMalignArr[i]);
//		printf("%s\n", alignResultsMalign[i]);
//	}

	for(i=0; i<maxAlignSeqLen; i++)
		printf("(%d,%d,%d,%d,%d), total=%d\n", consensusBaseNum[i][0], consensusBaseNum[i][1], consensusBaseNum[i][2], consensusBaseNum[i][3], consensusBaseNum[i][4], consensusBaseNum[i][5]);
#endif

	// release memory
	free(shiftPosArray);

	return SUCCESSFUL;
}

/**
 * Get the total number of reads for multi-align for paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getTotalReadsNumMalignPE(int32_t *totalReadsNum, int32_t readsNumDecisionTable, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable)
{
	int32_t i, newReadsNumArray[4], tmpNum;

	// compute the new reads number
	if(getNewReadsNumMalignPE(newReadsNumArray, kmerSeqIntAssembly, graph, decisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the number of new reads for multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the total reads number
	*totalReadsNum = readsNumDecisionTable;
	for(i=0; i<4; i++)
		(*totalReadsNum) += newReadsNumArray[i];

	if((*totalReadsNum)>maxItemNumReadsMalignArr)
	{
		tmpNum = maxItemNumReadsMalignArr;
		while((*totalReadsNum)>maxItemNumReadsMalignArr)
			maxItemNumReadsMalignArr *= 2;
		readseqMalignArr = (char **) realloc(readseqMalignArr, maxItemNumReadsMalignArr*sizeof(char*));
		if(readseqMalignArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
//		alignResultsMalign = (char **) realloc(alignResultsMalign, maxItemNumReadsMalignArr*sizeof(char*));
//		if(alignResultsMalign==NULL)
//		{
//			printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
//			return FAILED;
//		}

		for(i=tmpNum; i<maxItemNumReadsMalignArr; i++)
		{
			readseqMalignArr[i] = (char *) calloc(MAX_SEQ_LEN_MALIGN+1, sizeof(char));
			if(readseqMalignArr[i]==NULL)
			{
				printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
//			alignResultsMalign[i] = (char *) calloc(MAX_SEQ_LEN_MALIGN+1, sizeof(char));
//			if(alignResultsMalign[i]==NULL)
//			{
//				printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
//				return FAILED;
//			}
		}

		readIDMalignArr = (int64_t *) realloc(readIDMalignArr, maxItemNumReadsMalignArr*sizeof(int64_t));
		if(readIDMalignArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the total number of reads for multi-align for single ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getTotalReadsNumMalignSE(int32_t *totalReadsNum, int32_t readsNumDecisionTable, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable)
{
	int32_t i, newReadsNumArray[4], tmpNum;

	// compute the new reads number
	if(getNewReadsNumMalignSE(newReadsNumArray, kmerSeqIntAssembly, graph, decisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the number of new reads for multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the total reads number
	*totalReadsNum = readsNumDecisionTable;
	for(i=0; i<4; i++)
		(*totalReadsNum) += newReadsNumArray[i];

	if((*totalReadsNum)>maxItemNumReadsMalignArr)
	{
		tmpNum = maxItemNumReadsMalignArr;
		while((*totalReadsNum)>maxItemNumReadsMalignArr)
			maxItemNumReadsMalignArr *= 2;
		readseqMalignArr = (char **) realloc(readseqMalignArr, maxItemNumReadsMalignArr*sizeof(char*));
		if(readseqMalignArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
//		alignResultsMalign = (char **) realloc(alignResultsMalign, maxItemNumReadsMalignArr*sizeof(char*));
//		if(alignResultsMalign==NULL)
//		{
//			printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
//			return FAILED;
//		}

		for(i=tmpNum; i<maxItemNumReadsMalignArr; i++)
		{
			readseqMalignArr[i] = (char *) calloc(MAX_SEQ_LEN_MALIGN+1, sizeof(char));
			if(readseqMalignArr[i]==NULL)
			{
				printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
//			alignResultsMalign[i] = (char *) calloc(MAX_SEQ_LEN_MALIGN+1, sizeof(char));
//			if(alignResultsMalign[i]==NULL)
//			{
//				printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
//				return FAILED;
//			}
		}

		readIDMalignArr = (int64_t *) realloc(readIDMalignArr, maxItemNumReadsMalignArr*sizeof(int64_t));
		if(readIDMalignArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the number of new reads for multi-align for paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNewReadsNumMalignPE(int32_t *newReadsNum, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable)
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
		if(getNewReadsNumSingleKmerMalignPE(newReadsNum+i, tmp_kmers[i], graph, decisionTable)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the new number of reads of single kmer for multi-align, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the number of new reads for multi-align for single ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNewReadsNumMalignSE(int32_t *newReadsNum, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable)
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
		if(getNewReadsNumSingleKmerMalignSE(newReadsNum+i, tmp_kmers[i], graph, decisionTable)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the new number of reads of single kmer for multi-align, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the number of new reads of for multi-align for paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNewReadsNumSingleKmerMalignPE(int32_t *newReadsNum, kmertype *tmp_kmers[2], graphtype *graph, assemblingreadtype *decisionTable)
{
	int returnCode;
	ridpostype *rid_pos_table;
	readBlock_t *readBlockArr;
	read_t *pRead;

	int32_t i, posNum;
	int32_t rid, maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3;

	readBlockArr = graph->readSet->readBlockArr;
	maxItemNumPerReadBlock = graph->readSet->maxItemNumPerReadBlock;

	*newReadsNum = 0;
//	if(tmp_kmers[0])
//	{
//		rid_pos_table = tmp_kmers[0]->ppos;
//		posNum = tmp_kmers[0]->arraysize;
//		for(i=0; i<posNum; i++)
//		{
//			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
//			{
//				returnCode = validReadPair(rid_pos_table[i].rid);
//				if(returnCode==YES)
//				{
//					(*newReadsNum) ++;
//				}else if(returnCode==ERROR)
//				{
//					printf("In %s(), cannot valid read pair, error!\n", __func__);
//					return FAILED;
//				}
//			}
//		}
//	}

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
			seqLen = pRead->seqlen;
			errorRegLenEnd3 = ceil(seqLen * ERROR_REGION_LEN_3End_FACTOR);

			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos>seqLen-kmerSize+1-errorRegLenEnd3)
			{
				if(existReadInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
				{
					returnCode = validReadPair(rid_pos_table[i].rid);
					if(returnCode==YES)
					{
						(*newReadsNum) ++;
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
 * Get the number of new reads for multi-align for single ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNewReadsNumSingleKmerMalignSE(int32_t *newReadsNum, kmertype *tmp_kmers[2], graphtype *graph, assemblingreadtype *decisionTable)
{
	ridpostype *rid_pos_table;
	readBlock_t *readBlockArr;
	read_t *pRead;

	int32_t i, posNum;
	int32_t rid, maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3;

	readBlockArr = graph->readSet->readBlockArr;
	maxItemNumPerReadBlock = graph->readSet->maxItemNumPerReadBlock;

	*newReadsNum = 0;
	if(tmp_kmers[0])
	{
		rid_pos_table = tmp_kmers[0]->ppos;
		posNum = tmp_kmers[0]->arraysize;
		for(i=0; i<posNum; i++)
			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
				(*newReadsNum) ++;
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
			seqLen = pRead->seqlen;
			errorRegLenEnd3 = ceil(seqLen * ERROR_REGION_LEN_3End_FACTOR);

			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos>seqLen-kmerSize+1-errorRegLenEnd3)
				if(existReadInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
					(*newReadsNum) ++;
		}
	}

	return SUCCESSFUL;
}

/**
 * Add reads for multi-align for paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addReadseqMalignPE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *itemNumMalignArray, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable, int32_t readsNumDT)
{
	*itemNumMalignArray = 0;

	// get new sequences
	if(addNewReadseqMalignPE(readseqMalignArray, readIDMalignArray, itemNumMalignArray, kmerSeqIntAssembly, graph, decisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot add the new reads for multi-align for paired ends, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the sequences from decision table
	if(addNewReadseqInDecisionTableMalignPE(readseqMalignArray, readIDMalignArray, itemNumMalignArray, decisionTable, readsNumDT)==FAILED)
	{
		printf("line=%d, In %s(), cannot add the reads in decision table for multi-align for paired ends, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// sort the read seqs
	if(sortReadseqMalign(readseqMalignArray, readIDMalignArray, *itemNumMalignArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot get sort read sequences for multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Add reads for multi-align for single ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addReadseqMalignSE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *itemNumMalignArray, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable, int32_t readsNumDT)
{
	*itemNumMalignArray = 0;

	// add new sequences
	if(addNewReadseqMalignSE(readseqMalignArray, readIDMalignArray, itemNumMalignArray, kmerSeqIntAssembly, graph, decisionTable)==FAILED)
	{
		printf("line=%d, In %s(), cannot add the new reads for multi-align for paired ends, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the sequences from decision table
	if(addNewReadseqInDecisionTableMalignSE(readseqMalignArray, readIDMalignArray, itemNumMalignArray, decisionTable, readsNumDT)==FAILED)
	{
		printf("line=%d, In %s(), cannot add the reads in decision table for multi-align for paired ends, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// sort the read seqs
	if(sortReadseqMalign(readseqMalignArray, readIDMalignArray, *itemNumMalignArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot get sort read sequences for multi-align, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

short sortReadseqMalign(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t itemNumMalignArray)
{
	int32_t i, j;
	int64_t tmpRid;
	char tmpStr[MAX_SEQ_LEN_MALIGN+1];

	// sort the read sequences
	for(i=0; i<itemNumMalignArray-1; i++)
	{
		for(j=i+1; j<itemNumMalignArray; j++)
		{
			if(strlen(readseqMalignArray[i])<strlen(readseqMalignArray[j]))
			{
				strcpy(tmpStr, readseqMalignArray[i]);
				strcpy(readseqMalignArray[i], readseqMalignArray[j]);
				strcpy(readseqMalignArray[j], tmpStr);

				tmpRid = readIDMalignArray[i];
				readIDMalignArray[i] = readIDMalignArray[j];
				readIDMalignArray[j] = tmpRid;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Add new reads for multi-align for paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewReadseqMalignPE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *newReadsNumTotal, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable)
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

	*newReadsNumTotal = 0;
	for(i=0; i<4; i++)
	{
		if(addNewReadseqSingleKmerMalignPE(readseqMalignArray, readIDMalignArray, newReadsNumTotal, tmp_kmers[i], graph, decisionTable)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the new number of reads of single kmer for multi-align, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Add new reads for multi-align for single ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewReadseqMalignSE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *newReadsNumTotal, uint64_t *kmerseqInt, graphtype *graph, assemblingreadtype *decisionTable)
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

	*newReadsNumTotal = 0;
	for(i=0; i<4; i++)
	{
		if(addNewReadseqSingleKmerMalignSE(readseqMalignArray, readIDMalignArray, newReadsNumTotal, tmp_kmers[i], graph, decisionTable)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the new number of reads of single kmer for multi-align, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Add reads in decision table for multi-align for paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewReadseqInDecisionTableMalignPE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *itemNumMalignArray, assemblingreadtype *decisionTable, int32_t readsNumDT)
{
	int32_t i, newStartPos;

	// get the sequences from decision table
	for(i=0; i<readsNumDT; i++)
	{
		if(decisionTable[i].matedFlag==YES
			&& (decisionTable[i].unmatchBaseNum<=maxUnmatchBaseNumPerRead)
			&& (decisionTable[i].successiveAppearBases>=minSuccessiveAppearedBaseNum)
			&& decisionTable[i].basePos==decisionTable[i].lastMatchedBasePos)
		{
			if(decisionTable[i].orientation==ORIENTATION_PLUS)
			{
				newStartPos = decisionTable[i].basePos + 1;
				strcpy(readseqMalignArray[*itemNumMalignArray], getReadBaseByInt(decisionTable[i].readseq, decisionTable[i].seqlen)+newStartPos);
			}else
			{
				newStartPos = decisionTable[i].seqlen - decisionTable[i].basePos;
				strcpy(readseqMalignArray[*itemNumMalignArray], getReverseReadBaseByInt(decisionTable[i].readseq, decisionTable[i].seqlen)+newStartPos);
			}
			readIDMalignArray[*itemNumMalignArray] = decisionTable[i].rid;
			(*itemNumMalignArray) ++;
		}
	}

	return SUCCESSFUL;
}

/**
 * Add reads in decision table for multi-align for single ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewReadseqInDecisionTableMalignSE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *itemNumMalignArray, assemblingreadtype *decisionTable, int32_t readsNumDT)
{
	int32_t i, newStartPos;

	// get the sequences from decision table
	for(i=0; i<readsNumDT; i++)
	{
		if((decisionTable[i].unmatchBaseNum<=maxUnmatchBaseNumPerRead)
			&& (decisionTable[i].successiveAppearBases>=minSuccessiveAppearedBaseNum)
			&& decisionTable[i].basePos==decisionTable[i].lastMatchedBasePos)
		{
			if(decisionTable[i].orientation==ORIENTATION_PLUS)
			{
				newStartPos = decisionTable[i].basePos + 1;
				strcpy(readseqMalignArray[*itemNumMalignArray], getReadBaseByInt(decisionTable[i].readseq, decisionTable[i].seqlen)+newStartPos);
			}else
			{
				newStartPos = decisionTable[i].seqlen - decisionTable[i].basePos;
				strcpy(readseqMalignArray[*itemNumMalignArray], getReverseReadBaseByInt(decisionTable[i].readseq, decisionTable[i].seqlen)+newStartPos);
			}
			readIDMalignArray[*itemNumMalignArray] = decisionTable[i].rid;
			(*itemNumMalignArray) ++;
		}
	}

	return SUCCESSFUL;
}

/**
 * Add the number of new reads for multi-align for paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewReadseqSingleKmerMalignPE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *newReadsNumTotal, kmertype *tmp_kmers[2], graphtype *graph, assemblingreadtype *decisionTable)
{
	int32_t returnCode;
	ridpostype *rid_pos_table;
	readBlock_t *readBlockArr;
	read_t *pRead;
	readseqBlock_t *readseqBlockArr;
	uint64_t *pReadseq;

	int32_t i, posNum;
	int32_t rid, maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3;

	readBlockArr = graph->readSet->readBlockArr;
	maxItemNumPerReadBlock = graph->readSet->maxItemNumPerReadBlock;
	readseqBlockArr = deBruijnGraph->readSet->readseqBlockArr;

//	if(tmp_kmers[0])
//	{
//		rid_pos_table = tmp_kmers[0]->ppos;
//		posNum = tmp_kmers[0]->arraysize;
//		for(i=0; i<posNum; i++)
//		{
//			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
//			{
//				rid = rid_pos_table[i].rid;
//				blockRowRead = (rid - 1) / maxItemNumPerReadBlock;
//				itemRowBlockRead = (rid - 1) % maxItemNumPerReadBlock;
//				pRead = readBlockArr[blockRowRead].readArr + itemRowBlockRead;
//				pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
//
//				// add a read sequence
//				strcpy(readseqMalignArr[*newReadsNumTotal], getReadBaseByInt(pReadseq, pRead->seqlen)+kmerSize-1);
//				(*newReadsNumTotal) ++;
//			}
//		}
//	}

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
				if(existReadInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
				{
					returnCode = validReadPair(rid_pos_table[i].rid);
					if(returnCode==YES)
					{
						// add a read sequence
						strcpy(readseqMalignArray[*newReadsNumTotal], getReverseReadBaseByInt(pReadseq, seqLen)+rid_pos_table[i].pos-1);
						readIDMalignArray[*newReadsNumTotal] = rid_pos_table[i].rid;
						(*newReadsNumTotal) ++;
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
 * Add the number of new reads of for multi-align for single ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewReadseqSingleKmerMalignSE(char **readseqMalignArray, int64_t *readIDMalignArray, int32_t *newReadsNumTotal, kmertype *tmp_kmers[2], graphtype *graph, assemblingreadtype *decisionTable)
{
	ridpostype *rid_pos_table;
	readBlock_t *readBlockArr;
	read_t *pRead;
	readseqBlock_t *readseqBlockArr;
	uint64_t *pReadseq;

	int32_t i, posNum;
	int32_t rid, maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3;

	readBlockArr = graph->readSet->readBlockArr;
	maxItemNumPerReadBlock = graph->readSet->maxItemNumPerReadBlock;
	readseqBlockArr = deBruijnGraph->readSet->readseqBlockArr;

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

				// add a read sequence
				strcpy(readseqMalignArr[*newReadsNumTotal], getReadBaseByInt(pReadseq, pRead->seqlen)+kmerSize-1);
				readIDMalignArray[*newReadsNumTotal] = rid_pos_table[i].rid;
				(*newReadsNumTotal) ++;
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
				if(existReadInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
				{
					// add a read sequence
					strcpy(readseqMalignArray[*newReadsNumTotal], getReverseReadBaseByInt(pReadseq, seqLen)+rid_pos_table[i].pos-1);
					readIDMalignArray[*newReadsNumTotal] = rid_pos_table[i].rid;
					(*newReadsNumTotal) ++;
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Adjust reads of read sequences for multi-align.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short adjustReadseqMalign(char **readseqMalignArray, int32_t *shiftPosArray, int32_t itemNumMalignArray)
{
	int32_t i, j, k, seqlen_0, minValue;
	int32_t len, alignOK, startPos1, startPos2;


	// adjust the sequences
	seqlen_0 = strlen(readseqMalignArray[0]);
	for(i=1; i<itemNumMalignArray; i++)
	{
		len = strlen(readseqMalignArray[i]);
		if(seqlen_0<len)
			len = seqlen_0;
		if(len>15)
			len = 15;
//		else if(len<=3)
//			continue;

		//get the startPos1, startPos2
		startPos1 = startPos2 = 0;
		for(k=0; k<2; k++)
		{
			alignOK = YES;
			j = 0;
			if(j<len && j+startPos1<len)
			{
				while(j<len && j+startPos1<len)
				{
					if(readseqMalignArray[0][j+startPos1]!=readseqMalignArray[i][j])
					{
						alignOK = NO;
						break;
					}

					j++;
				}
			}else
			{
				alignOK = NO;
			}

			if(alignOK==YES)
				break;

			startPos1 ++;
		}

		//get the startPos2
		if(alignOK==NO)
		{
			startPos2 = 0;
			for(k=0; k<2; k++)
			{
				alignOK = YES;
				j = 0;
				if(j<len && j+startPos2<len)
				{
					while(j<len && j+startPos2<len)
					{
						if(readseqMalignArray[0][j]!=readseqMalignArray[i][j+startPos2])
						{
							alignOK = NO;
							break;
						}

						j++;
					}
				}else
				{
					alignOK = NO;
				}

				if(alignOK==YES)
					break;

				startPos2 ++;
			}
		}

		if(alignOK==NO)
			startPos1 = startPos2 = 0;

		if(startPos1>0)
		{
			shiftPosArray[i] = startPos1;
		}else if(startPos2>0)
		{
			shiftPosArray[i] = -startPos2;
		}
	}

	minValue = shiftPosArray[0];
	for(i=1; i<itemNumMalignArray; i++)
		if(minValue>shiftPosArray[i])
			minValue = shiftPosArray[i];

	if(minValue!=0)
		for(i=0; i<itemNumMalignArray; i++)
			shiftPosArray[i] -= minValue;

	return SUCCESSFUL;
}

/**
 * Compute the Consensus base numbers for multi-align.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeConsensusBaseNumMalign(int32_t **consensusBaseNumArray, int32_t *maxAlignSeqLen, int32_t rowsNum, int32_t colsNum, char **readseqMalignArray, int32_t *shiftPosArray, int32_t itemNumMalignArray)
{
	int32_t i, j, chNum, maxAlignSeqLenTmp, seqlen, baseInt;
	char *readseq;

	// compute the numbers
	if(memset(consensusBaseNumArray, 0, rowsNum * colsNum * sizeof(int32_t))==NULL)
	{
		printf("line=%d, In %s(), cannot memset memory to zero, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<itemNumMalignArray; i++)
	{
		readseq = readseqMalignArray[i];
		//readseq = alignResultsMalignArray[i];
		seqlen = strlen(readseq);
		chNum = shiftPosArray[i];
		for(j=0; j<chNum; j++)
		{
			(*((int32_t*)consensusBaseNumArray + j*colsNum + 4)) ++;  	// '-'
			(*((int32_t*)consensusBaseNumArray + j*colsNum + 5)) ++;  	 // for A+C+G+T

			//consensusBaseNumArray[j][4] ++;  	// '-'
			//consensusBaseNumArray[j][5] ++;  	 // for A+C+G+T
		}

		for(j=0; j<seqlen; j++)
		{
			switch(readseq[j])
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
				default: printf("line=%d, invalid base %d\n, error!\n", __LINE__, readseq[j]); return FAILED;
			}
//			(*((int32_t*)consensusBaseNumArray + j*colsNum + baseInt)) ++;  // for A, C, G, T
//			(*((int32_t*)consensusBaseNumArray + j*colsNum + 5)) ++;  	 // for A+C+G+T

			(*((int32_t*)consensusBaseNumArray + (j+chNum)*colsNum + baseInt)) ++;  // for A, C, G, T
			(*((int32_t*)consensusBaseNumArray + (j+chNum)*colsNum + 5)) ++;  	 // for A+C+G+T
		}
	}

	(*maxAlignSeqLen) = 0;
	for(i=0; i<itemNumMalignArray; i++)
	{
		maxAlignSeqLenTmp = strlen(readseqMalignArray[i]) + shiftPosArray[i];
//		maxAlignSeqLenTmp = strlen(alignResultsMalignArray[i]);
		if((*maxAlignSeqLen)<maxAlignSeqLenTmp)
			(*maxAlignSeqLen) = maxAlignSeqLenTmp;
	}

	if((*maxAlignSeqLen)>MAX_ALIGN_LEN_MALIGN)
		(*maxAlignSeqLen) = MAX_ALIGN_LEN_MALIGN;

	return SUCCESSFUL;
}

/**
 * Compute the incorrect base numbers for multi-align.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeIncorrectBaseNumMalign(int32_t *incorrectNum, int32_t **consensusBaseNumArray, int32_t colsNum, int32_t maxAlignSeqLen, char *consensusSeq, int32_t itemNumMalignArray)
{
	int32_t i, j, value, maxValue, secValue, maxIndex, secIndex, subSum;

	// compute the incorrect base number
	*incorrectNum = 0;
	for(i=0; i<maxAlignSeqLen; i++)
	{
		maxValue = secValue = -1;
		maxIndex = secIndex = -1;
		for(j=0; j<5; j++)
		{
			//if(consensusBaseNum[i][j]>maxValue)
			value = *((int32_t*)consensusBaseNumArray + i*colsNum + j);
			if(value>maxValue)
			{
				secValue = maxValue;
				secIndex = maxIndex;
				maxValue = value;
				maxIndex = j;
			}else if(value>secValue)
			{
				secValue = value;
				secIndex = maxIndex;
			}
		}

		subSum = *((int32_t*)consensusBaseNumArray + i*colsNum + 5);

		//if(itemNumMalignArray>10)
		if(subSum>10)
		{
			subSum = *((int32_t*)consensusBaseNumArray + i*colsNum + 5);
			if((double)maxValue/subSum<INCOR_RATIO_MALIGN)
				(*incorrectNum) ++;
		}else
		{
			if(secValue>0 && (double)maxValue/secValue<4)
				(*incorrectNum) ++;
		}

		*((int32_t*)consensusBaseNumArray + i*colsNum + 6) = maxIndex;
		*((int32_t*)consensusBaseNumArray + i*colsNum + 7) = secIndex;

		//consensusBaseNumArray[i][6] = maxIndex;
		//consensusBaseNumArray[i][7] = secIndex;

		switch(maxIndex)
		{
			case 0: consensusSeq[i] = 'A'; break;
			case 1: consensusSeq[i] = 'C'; break;
			case 2: consensusSeq[i] = 'G'; break;
			case 3: consensusSeq[i] = 'T'; break;
			case 4: consensusSeq[i] = '-'; break;
		}
	}
	consensusSeq[maxAlignSeqLen] = '\0';

	return SUCCESSFUL;
}

/**
 * Decide the navigation from the incorrect base numbers for multi-align for paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short decideNaviMalignPE(int32_t *naviMalign, int32_t *maxIndex, int32_t maxPEIndex, int32_t incorrectNum, int32_t **consensusBaseNumArray, int32_t colsNum, int32_t maxAlignSeqLen, char *consensusSeq)
{
	int32_t i, j, firstBaseRow, maxBaseIndex, maxValue, secondBaseIndex, secondValue, sum, invalidNum;

	if(incorrectNum>MAX_INCOR_BASE_NUM_MALIGN)
	{
//		if(incorrectNum>2*MAX_INCOR_BASE_NUM_MALIGN)
//		{
			*maxIndex = -1;
			*naviMalign = NAVI_FAILED;
//		}else
//		{
//			// skip the gaps
//			for(i=0; i<maxAlignSeqLen; i++)
//			{
//				if(consensusSeq[i]!='-')
//				{
//					//maxIndex = consensusBaseNum[i][6];
//					maxBaseIndex = *((int32_t*)consensusBaseNumArray + i*colsNum + 6);
//					secondBaseIndex = *((int32_t*)consensusBaseNumArray + i*colsNum + 7);
//					firstBaseRow = i;
//					break;
//				}
//			}
//
//			maxValue = *((int32_t*)consensusBaseNumArray + firstBaseRow*colsNum + maxBaseIndex);
//			secondValue = *((int32_t*)consensusBaseNumArray + firstBaseRow*colsNum + secondBaseIndex);
//
//			if(maxValue/secondValue>=4)
//			{
//				*maxIndex = maxBaseIndex;
//				*naviMalign = NAVI_SUCCESS;
//			}else
//			{
//				*maxIndex = -1;
//				*naviMalign = NAVI_FAILED;
//			}
//		}
	}else
	{
		// skip the gaps
		for(i=0; i<maxAlignSeqLen; i++)
		{
			if(consensusSeq[i]!='-')
			{
				//maxIndex = consensusBaseNum[i][6];
				maxBaseIndex = *((int32_t*)consensusBaseNumArray + i*colsNum + 6);
				secondBaseIndex = *((int32_t*)consensusBaseNumArray + i*colsNum + 7);
				firstBaseRow = i;
				break;
			}
		}

		maxValue = *((int32_t*)consensusBaseNumArray + firstBaseRow*colsNum + maxBaseIndex);
		secondValue = *((int32_t*)consensusBaseNumArray + firstBaseRow*colsNum + secondBaseIndex);
		sum = *((int32_t*)consensusBaseNumArray + firstBaseRow*colsNum + 5);

		if(maxPEIndex!=-1 && maxBaseIndex!=maxPEIndex)
		{
//			if(secondValue>=0.5*maxValue)
//			{
				*maxIndex = -1;
				*naviMalign = NAVI_FAILED;
//			}else
//			{
//				*maxIndex = maxBaseIndex;
//				*naviMalignSE = NAVI_SUCCESS;
//			}
		}else
		{
			if(maxValue==secondValue && incorrectNum>1 && sum>10)
			{
				*maxIndex = -1;
				*naviMalign = NAVI_FAILED;
			}
			else if(((double)maxValue/sum<MXA_OCC_RATIO_MALIGN) && incorrectNum>2)
			{
				*maxIndex = -1;
				*naviMalign = NAVI_FAILED;
			}
//			else if(incorrectNum>1)
//			{
//				invalidNum = 0;
//				for(j=0; j<maxAlignSeqLen; j++)
//				{
//					maxBaseIndex = *((int32_t*)consensusBaseNumArray + j*colsNum + 6);
//					secondBaseIndex = *((int32_t*)consensusBaseNumArray + j*colsNum + 7);
//					maxValue = *((int32_t*)consensusBaseNumArray + j*colsNum + maxBaseIndex);
//					secondValue = *((int32_t*)consensusBaseNumArray + j*colsNum + secondBaseIndex);
//					sum = *((int32_t*)consensusBaseNumArray + j*colsNum + 5);
//					if(secondValue>1 && (double)maxValue/sum<0.7)
//						invalidNum ++;
//				}
//
//				if(invalidNum==MAX_INCOR_BASE_NUM_MALIGN)
//				{
//					printf("invalidNum=%d\n", invalidNum);
//					*maxIndex = -1;
//					*naviMalign = NAVI_FAILED;
//				}else
//				{
//					printf("invalidNum=%d\n", invalidNum);
//					*maxIndex = maxBaseIndex;
//					*naviMalign = NAVI_SUCCESS;
//				}
//			}
			else
			{
				*maxIndex = maxBaseIndex;
				*naviMalign = NAVI_SUCCESS;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Decide the navigation from the incorrect base numbers for multi-align for single ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short decideNaviMalignSE(int32_t *naviMalignSE, int32_t incorrectNumSE, int32_t *maxIndex, int32_t maxPEIndex, int32_t **consensusBaseNumArray, int32_t colsNum, int32_t maxAlignSeqLen, char *consensusSeq)
{
	int32_t i, j, firstBaseRow, maxBaseIndex, maxValue, secondBaseIndex, secondValue, sum, invalidNum;

	if(incorrectNumSE>MAX_INCOR_BASE_NUM_MALIGN)
	{
//		if(incorrectNumSE>2*MAX_INCOR_BASE_NUM_MALIGN)
//		{
			*maxIndex = -1;
			*naviMalignSE = NAVI_FAILED;
//		}else
//		{
//			// skip the gaps
//			for(i=0; i<maxAlignSeqLen; i++)
//			{
//				if(consensusSeq[i]!='-')
//				{
//					//maxIndex = consensusBaseNum[i][6];
//					maxBaseIndex = *((int32_t*)consensusBaseNumArray + i*colsNum + 6);
//					secondBaseIndex = *((int32_t*)consensusBaseNumArray + i*colsNum + 7);
//					firstBaseRow = i;
//					break;
//				}
//			}
//
//			maxValue = *((int32_t*)consensusBaseNumArray + firstBaseRow*colsNum + maxBaseIndex);
//			secondValue = *((int32_t*)consensusBaseNumArray + firstBaseRow*colsNum + secondBaseIndex);
//			//sum = *((int32_t*)consensusBaseNumArray + firstBaseRow*colsNum + 5);
//
//			if(maxValue/secondValue>=4)
//			{
//				*maxIndex = maxBaseIndex;
//				*naviMalignSE = NAVI_SUCCESS;
//			}else
//			{
//				*maxIndex = -1;
//				*naviMalignSE = NAVI_FAILED;
//			}
//		}
	}else
	{
		// skip the gaps
		for(i=0; i<maxAlignSeqLen; i++)
		{
			if(consensusSeq[i]!='-')
			{
				//maxIndex = consensusBaseNum[i][6];
				maxBaseIndex = *((int32_t*)consensusBaseNumArray + i*colsNum + 6);
				secondBaseIndex = *((int32_t*)consensusBaseNumArray + i*colsNum + 7);
				firstBaseRow = i;
				break;
			}
		}

		maxValue = *((int32_t*)consensusBaseNumArray + firstBaseRow*colsNum + maxBaseIndex);
		secondValue = *((int32_t*)consensusBaseNumArray + firstBaseRow*colsNum + secondBaseIndex);
		sum = *((int32_t*)consensusBaseNumArray + firstBaseRow*colsNum + 5);

		if(maxPEIndex!=-1 && maxBaseIndex!=maxPEIndex)
		{
//			if(secondValue>=0.5*maxValue)
//			{
				*maxIndex = -1;
				*naviMalignSE = NAVI_FAILED;
//			}else
//			{
//				*maxIndex = maxBaseIndex;
//				*naviMalignSE = NAVI_SUCCESS;
//			}
		}else
		{
			if(maxValue==secondValue && incorrectNumSE>1 && sum>10)
			{
				*maxIndex = -1;
				*naviMalignSE = NAVI_FAILED;
			}
			else if(((double)maxValue/sum<MXA_OCC_RATIO_MALIGN) && incorrectNumSE>2)
			{
				*maxIndex = -1;
				*naviMalignSE = NAVI_FAILED;
			}
//			else if(incorrectNumSE>1)
//			{
//				invalidNum = 0;
//				for(j=0; j<maxAlignSeqLen; j++)
//				{
//					maxBaseIndex = *((int32_t*)consensusBaseNumArray + j*colsNum + 6);
//					secondBaseIndex = *((int32_t*)consensusBaseNumArray + j*colsNum + 7);
//					maxValue = *((int32_t*)consensusBaseNumArray + j*colsNum + maxBaseIndex);
//					secondValue = *((int32_t*)consensusBaseNumArray + j*colsNum + secondBaseIndex);
//					sum = *((int32_t*)consensusBaseNumArray + j*colsNum + 5);
//					if(secondValue>1 && (double)maxValue/sum<0.7)
//						invalidNum ++;
//				}
//
//				if(invalidNum==MAX_INCOR_BASE_NUM_MALIGN)
//				{
//					printf("invalidNum=%d\n", invalidNum);
//					*maxIndex = -1;
//					*naviMalign = NAVI_FAILED;
//				}else
//				{
//					printf("invalidNum=%d\n", invalidNum);
//					*maxIndex = maxBaseIndex;
//					*naviMalign = NAVI_SUCCESS;
//				}
//			}
			else
			{
				*maxIndex = maxBaseIndex;
				*naviMalignSE = NAVI_SUCCESS;
			}
		}
	}

	return SUCCESSFUL;
}
