/*
 * util.c
 *
 *  Created on: Jun 14, 2011
 *      Author: zhuxiao
 */

#include "inc/scafStdinc.h"
#include "inc/scafExtvab.h"


/**
 * Check the sorted result.
 *  @return:
 *  	if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short checkSortResult()
{
	printf("Begin checking sort result ...\n");

	uint64_t *pSeqArr_1, *pSeqArr_2;
	int i, k;
	for(i=0; i<itemNumSeqArr-1; i++)
	{
		pSeqArr_1 = seqArr + baseSeqIndexArr[i].startRow;
		pSeqArr_2 = seqArr + baseSeqIndexArr[i+1].startRow;

		for(k=0; k<entriesPerRead; k++)
		{
			if(pSeqArr_1[k] > pSeqArr_2[k])
			{
				printf("line=%d, In %s(), Sort error, i=%d, k=%d\n", __LINE__, __func__, i, k);
				return FAILED;
			}else if(pSeqArr_1[k] < pSeqArr_2[k])
			{
				break;
			}
		}
	}

	printf("End checking sort result, and congratulations.\n");

	return SUCCESSFUL;
}

/**
 * Check the contig index (CI).
 *  @return:
 *  	if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short checkContigIndex()
{
	printf("Begin checking contig index (CI) ...\n");

	int i;
	uint64_t *pUniqueSeq;
	for(i=0; i<SORT_KMER_ARR_SIZE; i++)
	{
		if(uniqueSeqKmerArr[i].itemNum>0)
		{
			pUniqueSeq = uniqueSeqArr + uniqueSeqKmerArr[i].itemRow * entriesPerRead;
			if((pUniqueSeq[0]>>rightShift_bitsNum)!=i)
			{
				printf("line=%d, In %s(), contig index error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	printf("End checking contig index (CI), and congratulations.\n");

	return SUCCESSFUL;
}

/**
 * Test sequence search.
 */
void testSeqSearch()
{
	// initialize the right shift bits number
	if(initRightShiftBitsNum()==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the right shift bits number, error!\n", __LINE__, __func__);
		exit(1);
	}

	char *ch_seq = "AGAAGCCACGAATACGATCAACAATGCTCGCCGCTTCGC";

	printf("len=%d\n", (int) strlen(ch_seq));

	uint64_t *seq_int;
	contigMatchInfo *pContigMatchInfo;
	int i, j, contigItemNum, hitRow;
	seq_int = (uint64_t*) calloc(entriesPerRead, sizeof(uint64_t));
	if(seq_int==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		exit(1);
	}

	for(i=0; i<10000000; i++)
	{
		if(seqHash(seq_int, ch_seq, 1)==FAILED)
		{
			printf("line=%d, In %s(), cannot hash the sequence [ %s ], error!\n", __LINE__, __func__, ch_seq);
			exit(1);
		}

		hitRow = seqSearch(seq_int);
		if(hitRow>=0)
		{
			pContigMatchInfo = contigMatchInfoArr + seqRowArr[hitRow].startRow;
			contigItemNum = seqRowArr[hitRow].rowsNum;

			printf("contigItemNum=%d\n", contigItemNum);
			for(j=0; j<contigItemNum; j++)
			{
				printf("\thitRow=%d, contigID=%d, contigPos=%d, ends=%d\n", hitRow, pContigMatchInfo[j].contigID, pContigMatchInfo[j].contigPos, pContigMatchInfo[j].contigEnd);
			}

		}else
		{
			printf("No hits.\n");
		}
	}


	free(seq_int);
	seq_int = NULL;

}

/**
 * Check the shared read list.
 */
void checkSharedReadList()
{
	printf("Begin checking shared read list ...\n");

	int64_t i, readsNum;
	ReadList *pReadListArr, *pReadListArr_1, *pReadListArr_2;

	readsNum = readNumInRL[0];
	pReadListArr = readListArr[0];

	for(i=0; i<readsNum-1; i++)
	{
		pReadListArr_1 = pReadListArr + i;
		pReadListArr_2 = pReadListArr + i + 1;

		if(pReadListArr_1->readID >= pReadListArr_2->readID)
		{
			printf("line=%d, In %s(), i=%ld, readID1=%lu, >=, readID2=%lu, error in shared read list.\n", __LINE__, __func__, i, pReadListArr_1->readID, pReadListArr_2->readID);
			exit(1);
		}
	}

	for(i=0; i<readsNum; i+=2)
	{
		pReadListArr_1 = pReadListArr + i;
		pReadListArr_2 = pReadListArr + i + 1;

		if((pReadListArr_1->readID != pReadListArr_2->readID - 1) || (pReadListArr_1->readID%2==0) || (pReadListArr_2->readID%2==1))
		{
			printf("line=%d, In %s(), i=%ld, readID1=%lu, readID2=%lu, error in shared read list.\n", __LINE__, __func__, i, pReadListArr_1->readID, pReadListArr_2->readID);
			exit(1);
		}
	}

	printf("End checking shared read list, and congratulations.\n");
}

/**
 * Check the sorted contig reads in contig list (CL) array.
 * If the reads are correctly sorted, return SUCCESSFUL; otherwise, return FAILED.
 */
short checkSortedContigReads(ContigList *pContigListArray, int contigItemNumInCL, ContigRead *pContigReadArray)
{
	printf("Begin checking the sorted contig reads in Contig List (CL) ...\n");

	int i, j;
	int contigID, contigPos1, contigPos2;
	int rowsNum;
	for(i=0; i<contigItemNumInCL; i++)
	{
		contigID = pContigListArray[i].contigID;

		// check the 5' end of a contig
		if(pContigListArray[i].EndNum5>1)
		{
			contigPos1 = pContigReadArray[ pContigListArray[i].firstRow5 ].contigPos;
			rowsNum = pContigListArray[i].EndNum5;
			for(j=1; j<rowsNum; j++)
			{
				contigPos2 = pContigReadArray[ pContigListArray[i].firstRow5 + j ].contigPos;
				if(contigPos1>contigPos2)
				{
					printf("line=%d, In %s(), at the 5' end of contig %d, contigPos1=%d is larger than contigPos2=%d, error!\n", __LINE__, __func__, contigID, contigPos1, contigPos2);
					return FAILED;
				}

				contigPos1 = contigPos2;
			}
		}

		// check the 3' end of a contig
		if(pContigListArray[i].EndNum3>1)
		{
			contigPos1 = pContigReadArray[ pContigListArray[i].firstRow3 ].contigPos;
			rowsNum = pContigListArray[i].EndNum3;
			for(j=1; j<rowsNum; j++)
			{
				contigPos2 = pContigReadArray[ pContigListArray[i].firstRow3 + j ].contigPos;
				if(contigPos1>contigPos2)
				{
					printf("line=%d, In %s(), at the 3' end of contig %d, contigPos1=%d is larger than contigPos2=%d, error!\n", __LINE__, __func__, contigID, contigPos1, contigPos2);
					return FAILED;
				}

				contigPos1 = contigPos2;
			}
		}
	}

	printf("End checking the sorted contig reads in Contig List (CL), and congratulations.\n");

	return SUCCESSFUL;
}

/**
 * Output the reads information at the end of a contig.
 */
void outputContigReadsSingleContigEnd(ContigRead *pContigReadArray, int readsNum)
{
	int i, totalNum;

	printf("readsNum=%d\n", readsNum);

	totalNum = 0;
	for(i=0; i<readsNum; i++)
	{
		totalNum ++;
		printf("contigPos=%d, readID=%ld, orientation=%d, end=%d\n", pContigReadArray[i].contigPos, pContigReadArray[i].readID, pContigReadArray[i].orientation, pContigReadArray[i].contigEnd);
	}

	printf("Total reads number=%d\n", totalNum);
}

/**
 * Check the arrLocArr.
 *  @return:
 *  	If correctly sorted, return YES; otherwise, return NO.
 */
short checkArrLocArray(arrLoc *arrLocArray, int64_t itemNumArrLocArray)
{
	int64_t i;

	for(i=0; i<itemNumArrLocArray-1; i++)
	{
		if(arrLocArray[i].row > arrLocArray[i+1].row)
		{
			printf("line=%d, In %s(), sort error!\n", __LINE__, __func__);
			return NO;
		}
	}

	printf("Congratulations!\n");

	return YES;
}


/**
 * Output the score array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputScoreArr(int *scoreArray, int rowsNum, int colsNum)
{
	int i, j;
	// print the score array
	for(i=0; i<rowsNum; i++)
	{
		for(j=0; j<colsNum-1; j++)
		{
			printf("%d\t", scoreArray[i*colsNum+j]);
		}
		printf("%d\n", scoreArray[i*colsNum+colsNum-1]);
	}

	return SUCCESSFUL;
}

/**
 * Output the score array to file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputScoreArrToFile(int *scoreArray, int rowsNum, int colsNum)
{
	FILE *fp;

	fp = fopen("../scoreArr.txt", "w");
	if(fp==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ ../scoreArr.txt ], error!\n", __LINE__, __func__);
		return FAILED;
	}

	int i, j;
	// print the score array
	for(i=0; i<rowsNum; i++)
	{
		for(j=0; j<colsNum-1; j++)
		{
			fprintf(fp, "%d\t", scoreArray[i*colsNum+j]);
		}
		fprintf(fp, "%d\n", scoreArray[i*colsNum+colsNum-1]);
	}

	fclose(fp);

	return SUCCESSFUL;
}

/**
 * Output the total number of gap regions, and their total length.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputGapSizeInScaf(contigOverlapIndex *pContigOverlapIndexArray, int scaffoldsNum, contigOverlap *pContigOverlapInfoArray, int minBaseNumInGap)
{
	int i, j, linkedContigsNum, rowsNum, startRow, gapSize, newGapSize;
	int64_t totalGapSize, gapNum;
	contigOverlap *pContigOverlapInfo;

	totalGapSize = 0;
	gapNum = 0;
	for(i=0; i<scaffoldsNum; i++)
	{
		linkedContigsNum = pContigOverlapIndexArray[i].linkedNum;
		if(linkedContigsNum>=2)
		{
			rowsNum = pContigOverlapIndexArray[i].rowsNum;
			startRow = pContigOverlapIndexArray[i].startRow;
			pContigOverlapInfo = pContigOverlapInfoArray + startRow;

			for(j=0; j<rowsNum; j++)
			{
				if(pContigOverlapInfo[j].mergeFlag==NO)
				{
					gapSize = pContigOverlapInfo[j].gapSize;
					if(gapSize<minBaseNumInGap)
					{
						newGapSize = minBaseNumInGap;
					}else
					{
						newGapSize = gapSize;
					}

					totalGapSize += newGapSize;
					gapNum ++;
				}
			}
		}else if(linkedContigsNum<=0)
		{
			printf("line=%d, In %s(), linkedContigsNum=%d, error!\n", __LINE__, __func__, linkedContigsNum);
			return FAILED;
		}
	}

	printf("Number of gaps: %ld, total size: %ld, mean: %.2f\n", gapNum, totalGapSize, (double)totalGapSize/gapNum);

	return SUCCESSFUL;
}

/**
 * Check the scafKmers in scafGrap in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short checkScafKmersInGrapInScaf(const char *monoSeqFile, scafGraph *pScafGrapDeBruijn)
{
	int i, pos;
	FILE *fpMonoSeq;
	uint64_t readID, hashcode;
	scafKmer *kmer;
	scafRidpos *ridpostable;
	int entryIndex;


	for(hashcode=0; hashcode<hashTableSize; hashcode++)
	{
		if(pScafGrapDeBruijn->pkmers[hashcode])
		{
			if(pScafGrapDeBruijn->pkmers[hashcode]->arraysize != pScafGrapDeBruijn->pkmers[hashcode]->multiplicity)
			{
				printf("line=%d, In %s(), arraysize (%d) != multiplicity (%d), error!\n", __LINE__, __func__, pScafGrapDeBruijn->pkmers[hashcode]->arraysize, pScafGrapDeBruijn->pkmers[hashcode]->multiplicity);
				return FAILED;
			}
		}
	}


	fpMonoSeq = fopen(monoSeqFile, "r");
	if(fpMonoSeq==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, monoSeqFile);
		return FAILED;
	}

	char *readSeq = (char*) calloc(readLen+1, sizeof(char));
	if(readSeq==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}


	while(fscanf(fpMonoSeq, "%ld\t%s\n", &readID, readSeq) != EOF)
	{
		// generate the kmer integer sequence
		if(generateKmerSeqIntInScaf(kmerSeqInt, readSeq)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}

		hashcode = kmerhashInScaf(kmerSeqInt);
		kmer = getKmerByHashInScaf(hashcode, kmerSeqInt, pScafGrapDeBruijn);
		if(kmer)
		{
			ridpostable = kmer->ppos;
			entryIndex = findDirectScafRidposIndexInScaf(readID, 1, ridpostable, kmer->arraysize);//二分查找rid_pos_table
			if(entryIndex>=0)
			{  //find the scafRidpos

			}else
			{
				printf("line=%d, In %s(), hashcode=%lu, cannot get the ridpos: (%lu,%d) memory, error!\n", __LINE__, __func__, hashcode, readID, i+1);
				if(outputKmerInScaf(hashcode, kmerSeqInt, pScafGrapDeBruijn)==FAILED)
				{
					printf("line=%d, In %s(), cannot output k-mer %s, error!\n", __LINE__, __func__, getKmerBaseByIntInScaf(kmerSeqInt));
				}
				return FAILED;
			}
		}else
		{
			printf("line=%d, In %s(), hashcode=%lu, cannot get the ridpos: (%lu,%d) memory, error!\n", __LINE__, __func__, hashcode, readID, i+1);
			if(outputKmerInScaf(hashcode, kmerSeqInt, pScafGrapDeBruijn)==FAILED)
			{
				printf("line=%d, In %s(), cannot output k-mer %s, error!\n", __LINE__, __func__, getKmerBaseByIntInScaf(kmerSeqInt));
			}
			return FAILED;
		}

		pos = 2;
		for(i=kmerSize; i<readLen-kmerSize+1; i++, pos++)
		{
			// generate the kmer integer sequence
			if(generateKmerSeqIntInScaf(kmerSeqInt, readSeq+pos-1)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			hashcode = kmerhashInScaf(kmerSeqInt);
			kmer = getKmerByHashInScaf(hashcode, kmerSeqInt, pScafGrapDeBruijn);
			if(kmer)
			{
				ridpostable = kmer->ppos;
				entryIndex = findDirectScafRidposIndexInScaf(readID, pos, ridpostable, kmer->arraysize);//二分查找rid_pos_table
				if(entryIndex>=0)
				{  //find the scafRidpos

				}else
				{
					printf("line=%d, In %s(), hashcode=%lu, cannot get the ridpos: (%lu,%d) memory, error!\n", __LINE__, __func__, hashcode, readID, i+1);
					if(outputKmerInScaf(hashcode, kmerSeqInt, pScafGrapDeBruijn)==FAILED)
					{
						printf("line=%d, In %s(), cannot output k-mer %s, error!\n", __LINE__, __func__, getKmerBaseByIntInScaf(kmerSeqInt));
					}
					return FAILED;
				}
			}else
			{
				printf("line=%d, In %s(), hashcode=%lu, cannot get the ridpos: (%lu,%d) memory, error!\n", __LINE__, __func__, hashcode, readID, i+1);
				if(outputKmerInScaf(hashcode, kmerSeqInt, pScafGrapDeBruijn)==FAILED)
				{
					printf("line=%d, In %s(), cannot output k-mer %s, error!\n", __LINE__, __func__, getKmerBaseByIntInScaf(kmerSeqInt));
				}
				return FAILED;
			}
		}
	}

	free(readSeq);
	readSeq = NULL;

	fclose(fpMonoSeq);
	fpMonoSeq = NULL;

	return SUCCESSFUL;
}

/**
 * 输出kmer中的内容.
 */
short outputKmerInScaf(uint64_t hashcode, uint64_t *seqInt, scafGraph *graph)
{
	scafKmer *kmer;
	scafRidpos *rid_pos;
	int i, posNum;

	kmer = getKmerByHashInScaf(hashcode, kmerSeqInt, graph);
	if(!kmer)
	{
		printf("line=%d, In %s(), the kmer==NULL.\n", __LINE__, __func__);
		return FAILED;
	}
	printf("kmerseq=%s, multi=%d, arraysize=%d\n", getKmerBaseByIntInScaf(seqInt), kmer->multiplicity, kmer->arraysize);

	rid_pos = kmer->ppos;
	posNum = kmer->arraysize;
	for(i=0; i<posNum; i++)  //输出ridpos表
	{
		printf("\trid=%lu, pos=%d, used=%d, reserved=%d\n", (uint64_t)rid_pos->rid, (int32_t)rid_pos->pos, (int32_t)rid_pos->used, (int32_t)rid_pos->reserved);
		rid_pos ++;
	}

	return SUCCESSFUL;
}

/**
 * Output scafContig nodes in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputContigInScaf(scafContig *contighead)
{
	scafContig *contig = contighead;
	scafSuccessRead *ridposorientation = NULL;
	int i = 0, j = 0, num = 0;
	char base, orientCh;
	while(contig)
	{
		switch(contig->base)
		{
			case 0: base = 'A'; break;
			case 1: base = 'C'; break;
			case 2: base = 'G'; break;
			case 3: base = 'T'; break;
			default: printf("line=%d, In %s(), error base integer form: %d, error!\n", __LINE__, __func__, contig->base); return FAILED;
		}

		printf("%d\t%c\t%d", contig->index, base, contig->ridposnum);
		ridposorientation = contig->pridposorientation;
		num = contig->ridposnum;
		for(i=0; i<num; i++)
		{
			switch(ridposorientation[i].orientation)
			{
			case ORIENTATION_PLUS: orientCh = '+'; break;
			case ORIENTATION_MINUS: orientCh = '-'; break;
			default: printf("line=%d, In %s(), error read orientation: %d, error!\n", __LINE__, __func__, ridposorientation[i].orientation); return FAILED;
			}

			printf("\t(%lu,%u,%u,%c)", (uint64_t)ridposorientation[i].rid, ridposorientation[i].startmatchpos, ridposorientation[i].matchnum, orientCh);
		}
		printf("\n");
		contig = contig->next;
		j++;
	}
	printf("There are %d nodes.\n", j);

	return SUCCESSFUL;
}

/**
 * Output scafContig sequence in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputScafContigSequence(scafContig *contighead)
{
	char base;
	scafContig *tmp_contig = contighead;
	while(tmp_contig)
	{
		switch(tmp_contig->base)
		{
			case 0: base = 'A'; break;
			case 1: base = 'C'; break;
			case 2: base = 'G'; break;
			case 3: base = 'T'; break;
			default: printf("line=%d, In %s(), unknown base %d, error!\n", __LINE__, __func__, tmp_contig->base); return FAILED;
		}
		printf("%c", base);
		tmp_contig = tmp_contig->next;
	}

	printf("\n");

	return SUCCESSFUL;
}

/**
 * 输出成功的reads信息.
 */
void outputSuccessReads(scafSuccessRead *pScafSuccessReadArr, int tablesize)
{
	int i;
	for(i=0; i<tablesize; i++)
	{
		printf("Success Reads[%d]: rid=%lu, startmatchpos=%u, matchnum=%u, orientation=%lu\n",
				i+1, (uint64_t)pScafSuccessReadArr[i].rid, pScafSuccessReadArr[i].startmatchpos, pScafSuccessReadArr[i].matchnum, (uint64_t)pScafSuccessReadArr[i].orientation);
	}
	printf("There are %d success reads\n", i);
}


/**
 * Output the regions of contig sequences in scaffolds to text file.
 *  Format:
 *  	(1) Header: > scaffoldID, linked contigs number, scaffold length;
 *  	(2) Body: contigID, contigLen, contig orientation, start region position, end region position.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputContigSeqRegionsInScaf(char *contigRegionFile, contigOverlapIndex *pContigOverlapIndexArray, int scaffoldsNum, contigOverlap *pContigOverlapInfoArray, contigInfo *pContigInfoArray)
{
	int i, j, scaffoldID, linkedContigsNum, rowsNum, startRow;
	int contigID, contigLen, contigOrient;
	int mergeFlag, overlapLen, gapSize, breakFlag, newGapSize;
	int scaffoldLen, totalScaffoldLen, startRegPos, endRegPos;
	contigOverlap *pContigOverlapInfo;
	FILE *fpRegion;

	char fullContigRegionFile[256];

	strcpy(fullContigRegionFile, outputPathStr);
	if(strlen(outputPrefix)>0)
	{
		strcat(fullContigRegionFile, outputPrefix);
		strcat(fullContigRegionFile, "_");
	}
	strcat(fullContigRegionFile, contigRegionFile);


	// check variables
	if(minBaseNumInGap!=MIN_BASENUM_IN_GAP)
	{
		printf("line=%d, In %s(), minBaseNumInGap!=%d, error!\n", __LINE__, __func__, MIN_BASENUM_IN_GAP);
		return FAILED;
	}
	if(pContigOverlapIndexArray==NULL || pContigInfoArray==NULL)
	{
		printf("line=%d, In %s(), invalid pointer, error!\n", __LINE__, __func__);
		return FAILED;
	}

	fpRegion = fopen(fullContigRegionFile, "w");
	if(fpRegion==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fullContigRegionFile);
		return FAILED;
	}

	// get the regions and output them to file
	for(i=0; i<scaffoldsNum; i++)
	{
		scaffoldID = pContigOverlapIndexArray[i].scaffoldID;
		linkedContigsNum = pContigOverlapIndexArray[i].linkedNum;
		rowsNum = pContigOverlapIndexArray[i].rowsNum;
		startRow = pContigOverlapIndexArray[i].startRow;
		pContigOverlapInfo = pContigOverlapInfoArray + startRow;

		// deal with the head link
		if(linkedContigsNum==1)
		{ // only one contig in the scaffold, then output the contig directly
			contigID = pContigOverlapInfo->contigID1;
			contigLen = pContigInfoArray[contigID-1].contigLen;
			contigOrient = pContigOverlapInfo->orientation1;
			scaffoldLen = contigLen;
			startRegPos = 1;
			endRegPos = scaffoldLen;
			totalScaffoldLen = scaffoldLen;
			fprintf(fpRegion, ">%d\t%d\t%d\n", scaffoldID, linkedContigsNum, totalScaffoldLen);
			fprintf(fpRegion, "%d\t%d\t%d\t%d\t%d\n", contigID, contigLen, contigOrient, startRegPos, endRegPos);
		}else
		{
			if(getSingleScaffoldSeqLen(&totalScaffoldLen, pContigOverlapInfo, linkedContigsNum, rowsNum, minBaseNumInGap, pContigInfoArray)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the scaffold sequence length, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(totalScaffoldLen<=0)
			{
				printf("line=%d, In %s(), totalScaffoldLen=%d, error!\n", __LINE__, __func__, totalScaffoldLen);
				return FAILED;
			}

			contigID = pContigOverlapInfo[0].contigID1;
			contigLen = pContigInfoArray[contigID-1].contigLen;
			contigOrient = pContigOverlapInfo[0].orientation1;
			scaffoldLen = contigLen;
			startRegPos = 1;
			endRegPos = scaffoldLen;
			fprintf(fpRegion, ">%d\t%d\t%d\n", scaffoldID, linkedContigsNum, totalScaffoldLen);
			fprintf(fpRegion, "%d\t%d\t%d\t%d\t%d\n", contigID, contigLen, contigOrient, startRegPos, endRegPos);


			for(j=0; j<rowsNum; j++)
			{
				contigID = pContigOverlapInfo[j].contigID2;
				contigOrient = pContigOverlapInfo[j].orientation2;
				contigLen = pContigInfoArray[contigID-1].contigLen;
				mergeFlag = pContigOverlapInfo[j].mergeFlag;
				overlapLen = pContigOverlapInfo[j].overlapLen;
				gapSize = pContigOverlapInfo[j].gapSize;
				breakFlag = pContigOverlapInfo[j].breakFlag;

				if(breakFlag==YES)
				{
					printf("line=%d, In %s(), the scaffold should be broken, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// deal with the overlap or gap
				if(mergeFlag==YES)
				{ // an overlap happens
					startRegPos = scaffoldLen - overlapLen + 1;
					scaffoldLen += contigLen - overlapLen;
					endRegPos = scaffoldLen;
				}else
				{ // a gap happens

					if(gapSize<minBaseNumInGap)
						newGapSize = minBaseNumInGap;
					else
						newGapSize = gapSize;

					startRegPos = scaffoldLen + newGapSize + 1;
					scaffoldLen += contigLen + newGapSize;
					endRegPos = scaffoldLen;
				}

				fprintf(fpRegion, "%d\t%d\t%d\t%d\t%d\n", contigID, contigLen, contigOrient, startRegPos, endRegPos);
			}
		}

		// ############################ Debug information ###########################
		if(scaffoldLen!=totalScaffoldLen)
		{
			printf("line=%d, In %s(), scaffoldLen=%d, totalScaffoldLen=%d, they are not equal, error!\n", __LINE__, __func__, scaffoldLen, totalScaffoldLen);
			return FAILED;
		}
		// ############################ Debug information ###########################
	}

	fclose(fpRegion);
	fpRegion = NULL;

	return SUCCESSFUL;
}

/**
 * Get single scaffold sequence length.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getSingleScaffoldSeqLen(int *scaffoldLen, contigOverlap *pContigOverlapInfo, int linkedContigsNum, int rowsNum, int minBaseNumInGap, contigInfo *pContigInfoArray)
{
	int j, contigID, contigLen;
	int mergeFlag, overlapLen, gapSize, breakFlag, newGapSize;

	// check variables
	if(minBaseNumInGap!=MIN_BASENUM_IN_GAP)
	{
		printf("line=%d, In %s(), minBaseNumInGap!=%d, error!\n", __LINE__, __func__, MIN_BASENUM_IN_GAP);
		return FAILED;
	}
	if(pContigOverlapInfo==NULL)
	{
		printf("line=%d, In %s(), invalid pointer, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// deal with the head link
	if(linkedContigsNum==1)
	{ // only one contig in the scaffold, then output the contig directly
		contigID = pContigOverlapInfo->contigID1;
		contigLen = pContigInfoArray[contigID-1].contigLen;
		*scaffoldLen = contigLen;
	}else
	{
		contigID = pContigOverlapInfo[0].contigID1;
		contigLen = pContigInfoArray[contigID-1].contigLen;
		*scaffoldLen = contigLen;

		for(j=0; j<rowsNum; j++)
		{
			contigID = pContigOverlapInfo[j].contigID2;
			contigLen = pContigInfoArray[contigID-1].contigLen;
			mergeFlag = pContigOverlapInfo[j].mergeFlag;
			overlapLen = pContigOverlapInfo[j].overlapLen;
			gapSize = pContigOverlapInfo[j].gapSize;
			breakFlag = pContigOverlapInfo[j].breakFlag;

			if(breakFlag==YES)
			{
				printf("line=%d, In %s(), the scaffold should be broken, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// deal with the overlap or gap
			if(mergeFlag==YES)
			{ // an overlap happens
				(*scaffoldLen) += contigLen - overlapLen;
			}else
			{ // a gap happens

				if(gapSize<minBaseNumInGap)
					newGapSize = minBaseNumInGap;
				else
					newGapSize = gapSize;

				(*scaffoldLen) += contigLen + newGapSize;
			}
		}
	}

	return SUCCESSFUL;
}

