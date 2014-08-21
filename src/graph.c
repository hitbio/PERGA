/*
 * graph.c
 *
 *  Created on: Jul 2, 2010
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Build k-mer hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructGraph(char *graphFileName, char **readsFileNames, int readsFileNum)
{
	printf("\n============= Begin to construct k-mer hash table, please wait ... =============\n");

	struct timeval tpstart, tpend;
	double timeused_graph;
	gettimeofday(&tpstart, NULL);


	if(constructGraphFromReadset(graphFileName, readsFileNames, readsFileNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot construct the graph from read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// output the graph to file
	if(outputGraphToFile(graphFileName, deBruijnGraph)==FAILED)
	{
		printf("line=%d, In %s(), can not output graph to file. Error!\n", __LINE__, __func__);
		return FAILED;
	}

	gettimeofday(&tpend, NULL);
	timeused_graph = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Total used time: %.2f seconds.\n", timeused_graph);

	printf("============= End constructed k-mer hash table. =============\n");

	return SUCCESSFUL;
}

/**
 * Build k-mer hash table from read set.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructGraphFromReadset(char *graphFileName, char **readsFileNames, int readsFileNum)
{
	double kmerOccTmp;
	int32_t endKmerNumTmp, kmerNumTmp;

	if(constructReadset(&readSet, readsFileNames, readsFileNum, reserveHashItemBlocksFlag)==FAILED)
	{
		printf("line=%d, In %s(), cannot construct the readset, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize the k-mer hash table
	if(initgraph(&deBruijnGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the graph, error!\n", __LINE__, __func__);
		return FAILED;
	}
	deBruijnGraph->readSet = readSet;

	// count the k-mer occurrences
	if(countKmerOccsFromReadset(deBruijnGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot count the kmers from read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// add the k-mer ridpos information
	if(addKmerRidposFromReadset(deBruijnGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot count the kmers from read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

#if (DEBUG_PARA_PRINT==YES)
	kmerOccTmp = (double)deBruijnGraph->totalItemNumRidpos / deBruijnGraph->totalItemNumKmer;
	kmerNumTmp = averReadLenInFileSample - kmerSize + 1;

	endKmerNumTmp = ceil(readLen * kmerRegLenRatioEnd5) + ceil(readLen * kmerRegLenRatioEnd3);
	kmerOccTmp *= (double)kmerNumTmp / endKmerNumTmp;

/*
	endKmerNumTmp = ceil(averReadLenInFileSample * kmerRegLenRatioEnd5) + ceil(averReadLenInFileSample * kmerRegLenRatioEnd3);
	//endKmerNumTmp = (ceil(averReadLenInFileSample * kmerRegLenRatioEnd5)/deBruijnGraph->kmerSampleInterval + 1) + (ceil(averReadLenInFileSample * kmerRegLenRatioEnd3)/deBruijnGraph->kmerSampleInterval + 1);
	kmerOccTmp *= (double)kmerNumTmp / endKmerNumTmp;
	kmerOccTmp *= (double)readLen / averReadLenInFileSample;
*/
	//kmerOccTmp /= 1.4;
	kmerOccTmp /= 1.6;
	printf("totalReads=%ld, validReads=%ld, validRatio=%.4f, kmerNum=%ld, ridposNum=%ld, averKmerOcc=%.2f, adjustKmerOcc=%.2f\n", readSet->totalItemNumRead, readSet->totalValidItemNumRead, (double)readSet->totalValidItemNumRead/readSet->totalItemNumRead, deBruijnGraph->totalItemNumKmer, deBruijnGraph->totalItemNumRidpos, (double)deBruijnGraph->totalItemNumRidpos/deBruijnGraph->totalItemNumKmer, kmerOccTmp);
#endif

	return SUCCESSFUL;
}

/**
 * Count the k-mer occurrences from read set.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short countKmerOccsFromReadset(graphtype *deBruijnGraph)
{
	int32_t i, j, percent;
	readSet_t *readSet;
	readBlock_t *pReadBlockArray;
	read_t *pRead, *pReadArray;
	int64_t rid, totalReadNumTmp;

	struct timeval tpstart, tpend;
	double timeused_counting;
	gettimeofday(&tpstart, NULL);


	printf("Filling k-mer information ...\n");

	percent = 0;
	readSet = deBruijnGraph->readSet;
	totalReadNumTmp = readSet->totalItemNumRead;

	pKmerBlockTmp = deBruijnGraph->kmerBlockArr + deBruijnGraph->blocksNumKmer - 1;
	pKmerTmpDoing = pKmerBlockTmp->kmerArr;

	pKmerseqBlockTmp = deBruijnGraph->kmerSeqBlockArr + deBruijnGraph->blocksNumKmerSeq - 1;
	pKmerSeqTmpDoing = pKmerseqBlockTmp->kmerSeqArr;

	// get the reads from read set
	rid = 0;
	for(i=0; i<readSet->blocksNumRead; i++)
	{
		pReadBlockArray = readSet->readBlockArr + i;
		pReadArray = pReadBlockArray->readArr;
		for(j=0; j<pReadBlockArray->itemNum; j++)
		{
			rid ++;
			pRead = pReadArray + j;

			if(pRead->validFlag==YES)
			{
				//count the kmers
				if(addReadPreFromReadset(pRead, deBruijnGraph)==FAILED)
				{
					printf("line=%d, In %s(), cannot count the kmer, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}

			// update the processing percentage
			if((int)((double)rid/totalReadNumTmp*100) > percent)
			{
				percent = (int)((double)rid/totalReadNumTmp*100);
				printf("%d%%", percent);
				if(percent%10==0) printf("\n");
				else printf("\t");
				fflush(stdout);
			}

			//fprintf(fpOut, "rid=%ld, seqlen=%d, nBaseNum=%d, validFlag=%d, successFlag=%d, seq=%s\n", rid, (int32_t)pRead->seqlen, (int32_t)pRead->nBaseNum, (int32_t)pRead->validFlag, (int32_t)pRead->successFlag, getReadBaseByInt(readseqInt, pRead->seqlen));
		}
	}

	// the last k-mer block is empty, remove it
	if(pKmerBlockTmp->itemNum==0)
	{
		free(pKmerBlockTmp->kmerArr);
		deBruijnGraph->kmerBlockArr[deBruijnGraph->blocksNumKmer-1].kmerArr = NULL;
		deBruijnGraph->blocksNumKmer --;
	}

	// the last kmerseq block is empty, remove it
	if(pKmerseqBlockTmp->itemNum==0)
	{
		free(pKmerseqBlockTmp->kmerSeqArr);
		deBruijnGraph->kmerSeqBlockArr[deBruijnGraph->blocksNumKmerSeq-1].kmerSeqArr = NULL;
		deBruijnGraph->blocksNumKmerSeq --;
	}


	gettimeofday(&tpend, NULL);
	timeused_counting = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Filling k-mers used time: %.2f seconds.\n", timeused_counting);

	return SUCCESSFUL;
}

/**
 * Add the k-mer ridpos information from read set.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addKmerRidposFromReadset(graphtype *graph)
{
	int32_t i, j, percent;
	int64_t rid, totalReadNumTmp;
	readSet_t *readSet;
	readBlock_t *pReadBlockArray;
	read_t *pRead, *pReadArray;

	struct timeval tpstart, tpend;
	double timeused_adding;
	gettimeofday(&tpstart, NULL);


	printf("Filling reads information ...\n");

	rid = 0;
	totalReadNumTmp = graph->readSet->totalItemNumRead;
	percent = 0;

	// initialize ridpos blocks and the ridpos regions in that blocks
	if(initRidposBlocksInGraph(graph)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the reads from read set
	rid = 0;
	readSet = graph->readSet;
	for(i=0; i<readSet->blocksNumRead; i++)
	{
		pReadBlockArray = readSet->readBlockArr + i;
		pReadArray = pReadBlockArray->readArr;
		for(j=0; j<pReadBlockArray->itemNum; j++)
		{
			rid ++;
			pRead = pReadArray + j;

			if(pRead->validFlag==YES)
			{
				//count the kmers
				if(addReadFromReadset(rid, pRead, graph)==FAILED)
				{
					printf("line=%d, In %s(), cannot count the kmer, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}

			// update the processing percentage
			if((int)((double)rid/totalReadNumTmp*100) > percent)
			{
				percent = (int)((double)rid/totalReadNumTmp*100);
				printf("%d%%", percent);
				if(percent%10==0) printf("\n");
				else printf("\t");
				fflush(stdout);
			}

			//fprintf(fpOut, "rid=%ld, seqlen=%d, nBaseNum=%d, validFlag=%d, successFlag=%d, seq=%s\n", rid, (int32_t)pRead->seqlen, (int32_t)pRead->nBaseNum, (int32_t)pRead->validFlag, (int32_t)pRead->successFlag, getReadBaseByInt(readseqInt, pRead->seqlen));
		}
	}

	gettimeofday(&tpend, NULL);
	timeused_adding = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Filling reads information used time: %.2f seconds.\n", timeused_adding);

	return SUCCESSFUL;
}

/**
 * Count the k-mer occurrences of a read from read set.
 *  @return:
 *  	If succeeds, return SUCESSFUL; otherwise, return FAILED.
 */
short addReadPreFromReadset(read_t *pRead, graphtype *graph)
{
	int32_t i, j, basePos, startKmerPos, startBasePos, endBasePos, endKmerNum;
	int32_t baseInt, seqLen, entriesNum, baseNumLastEntry, entryRow, entryPos, kmerIntervalTmp;
	uint64_t hashcode, *readseqInt;

	readseqInt = graph->readSet->readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
	seqLen = pRead->seqlen;
	entriesNum = ((seqLen - 1) >> 5) + 1;
	baseNumLastEntry = ((seqLen - 1) % 32) + 1;


	if(ceil(seqLen * kmerRegLenRatioEnd5) + ceil(seqLen * kmerRegLenRatioEnd3) >= seqLen - kmerSize + 1)
	{ // the 5' + 3' >= kmerNum
		// generate the kmer integer sequence
		if(generateKmerSeqIntFromReadset(pKmerSeqTmpDoing, readseqInt, 0, entriesNum, baseNumLastEntry)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}
		pKmerSeqTmpDone = pKmerSeqTmpDoing;

		// get the hashcode and count the kmer occurrence
		kmerIntervalTmp = 0;
		hashcode = kmerhashInt(pKmerSeqTmpDoing);
		if(countKmer(hashcode, graph)==FAILED)
		{
			printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
			return FAILED;
		}

		entryRow = kmerSize >> 5;
		entryPos = kmerSize % 32;
		for(basePos=kmerSize; basePos<seqLen; basePos++)
		{
			// get the baseInt
			if(entryRow<entriesNum-1)
				baseInt = (readseqInt[entryRow] >> (62-2*entryPos)) & 3;
			else
				baseInt = (readseqInt[entryRow] >> (2*(baseNumLastEntry-1)-2*entryPos)) & 3;

			entryPos ++;
			if(entryPos==32)
			{
				entryRow ++;
				entryPos = 0;
			}

			// generate the kmer integer sequence
			if(entriesPerKmer>=2)
			{
				for(j=0; j<entriesPerKmer-2; j++)
				{
					pKmerSeqTmpDoing[j] = (pKmerSeqTmpDone[j] << 2) | (pKmerSeqTmpDone[j+1] >> 62);
				}

				pKmerSeqTmpDoing[entriesPerKmer-2] = (pKmerSeqTmpDone[entriesPerKmer-2] << 2) | (pKmerSeqTmpDone[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
			}
			pKmerSeqTmpDoing[entriesPerKmer-1] = ((pKmerSeqTmpDone[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

			pKmerSeqTmpDone = pKmerSeqTmpDoing;

			// get the hashcode and count the kmer occurrence
			kmerIntervalTmp ++;
			if(kmerIntervalTmp==graph->kmerSampleInterval || basePos==seqLen-1)
			{
				hashcode = kmerhashInt(pKmerSeqTmpDoing);
				if(countKmer(hashcode, graph)==FAILED)
				{
					printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
					return FAILED;
				}
				kmerIntervalTmp = 0;
			}
		}
	}
	else
	{
		for(i=0; i<2; i++)
		{
			if(i==0)
			{
				endKmerNum = ceil(seqLen * kmerRegLenRatioEnd5);
				startKmerPos = 0;
				startBasePos = kmerSize;
				endBasePos = startBasePos + endKmerNum - 2;
			}
			else
			{
				endKmerNum = ceil(seqLen * kmerRegLenRatioEnd3);
				//startKmerPos = seqLen - endKmerNum - kmerSize + 1;
				//startBasePos = seqLen - endKmerNum + 1;
				startKmerPos = seqLen - kmerSize - 1 - ((endKmerNum-1)/graph->kmerSampleInterval)*graph->kmerSampleInterval + 1;
				startBasePos = startKmerPos + kmerSize;
				endBasePos = seqLen - 1;
			}

			// generate the kmer integer sequence
			if(generateKmerSeqIntFromReadset(pKmerSeqTmpDoing, readseqInt, startKmerPos, entriesNum, baseNumLastEntry)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}
			pKmerSeqTmpDone = pKmerSeqTmpDoing;

			// get the hashcode and count the kmer occurrence
			kmerIntervalTmp = 0;
			hashcode = kmerhashInt(pKmerSeqTmpDoing);
			if(countKmer(hashcode, graph)==FAILED)
			{
				printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			entryRow = startBasePos >> 5;
			entryPos = startBasePos % 32;
			for(basePos=startBasePos; basePos<=endBasePos; basePos++)
			{
				// get the baseInt
				if(entryRow<entriesNum-1)
					baseInt = (readseqInt[entryRow] >> (62-2*entryPos)) & 3;
				else
					baseInt = (readseqInt[entryRow] >> (2*(baseNumLastEntry-1)-2*entryPos)) & 3;

				entryPos ++;
				if(entryPos==32)
				{
					entryRow ++;
					entryPos = 0;
				}

				// generate the kmer integer sequence
				if(entriesPerKmer>=2)
				{
					for(j=0; j<entriesPerKmer-2; j++)
					{
						pKmerSeqTmpDoing[j] = (pKmerSeqTmpDone[j] << 2) | (pKmerSeqTmpDone[j+1] >> 62);
					}

					pKmerSeqTmpDoing[entriesPerKmer-2] = (pKmerSeqTmpDone[entriesPerKmer-2] << 2) | (pKmerSeqTmpDone[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
				}
				pKmerSeqTmpDoing[entriesPerKmer-1] = ((pKmerSeqTmpDone[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

				pKmerSeqTmpDone = pKmerSeqTmpDoing;

				// get the hashcode and count the kmer occurrence
				kmerIntervalTmp ++;
				if(kmerIntervalTmp==graph->kmerSampleInterval)
				{
					hashcode = kmerhashInt(pKmerSeqTmpDoing);
					if(countKmer(hashcode, graph)==FAILED)
					{
						printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
						return FAILED;
					}
					kmerIntervalTmp = 0;
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Add a read to graph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addReadFromReadset(int64_t rid, read_t *pRead, graphtype *graph)
{
	int32_t i, j, pos, basePos, startKmerPos, startBasePos, endBasePos, endKmerNum;
	int32_t baseInt, seqLen, entriesNum, baseNumLastEntry, entryRow, entryPos, kmerIntervalTmp;
	uint64_t hashcode, *readseqInt;


	readseqInt = graph->readSet->readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
	seqLen = pRead->seqlen;
	entriesNum = ((seqLen - 1) >> 5) + 1;
	baseNumLastEntry = ((seqLen - 1) % 32) + 1;

	if(ceil(seqLen * kmerRegLenRatioEnd5) + ceil(seqLen * kmerRegLenRatioEnd3) >= seqLen - kmerSize + 1)
	{ // the 5' + 3' >= kmerNum
		// generate the kmer integer sequence
		if(generateKmerSeqIntFromReadset(kmerSeqInt, readseqInt, 0, entriesNum, baseNumLastEntry)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// get the hashcode and count the kmer occurrence
		kmerIntervalTmp = 0;
		hashcode = kmerhashInt(kmerSeqInt);
		if(addKmer(hashcode, kmerSeqInt, rid, 1, graph)==FAILED)
		{
			printf("line=%d, In %s(), cannot add kmer for read %lu, error!\n", __LINE__, __func__, rid);
			return FAILED;
		}

		pos = 2;
		entryRow = kmerSize >> 5;
		entryPos = kmerSize % 32;
		for(basePos=kmerSize; basePos<seqLen; basePos++, pos++)
		{
			// get the baseInt
			if(entryRow<entriesNum-1)
				baseInt = (readseqInt[entryRow] >> (62-2*entryPos)) & 3;
			else
				baseInt = (readseqInt[entryRow] >> (2*(baseNumLastEntry-1)-2*entryPos)) & 3;

			entryPos ++;
			if(entryPos==32)
			{
				entryRow ++;
				entryPos = 0;
			}

			// generate the kmer integer sequence
			if(entriesPerKmer>=2)
			{
				for(j=0; j<entriesPerKmer-2; j++)
				{
					kmerSeqInt[j] = (kmerSeqInt[j] << 2) | (kmerSeqInt[j+1] >> 62);
				}
				kmerSeqInt[entriesPerKmer-2] = (kmerSeqInt[entriesPerKmer-2] << 2) | (kmerSeqInt[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
			}
			kmerSeqInt[entriesPerKmer-1] = ((kmerSeqInt[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

			// get the hashcode and count the kmer occurrence
			kmerIntervalTmp ++;
			if(kmerIntervalTmp==graph->kmerSampleInterval || basePos==seqLen-1)
			{
				hashcode = kmerhashInt(kmerSeqInt);
				if(addKmer(hashcode, kmerSeqInt, rid, pos, graph)==FAILED)
				{
					printf("line=%d, In %s(), cannot add kmer, error!\n", __LINE__, __func__);
					return FAILED;
				}
				kmerIntervalTmp = 0;
			}
		}
	}
	else
	{
		for(i=0; i<2; i++)
		{
			if(i==0)
			{
				endKmerNum = ceil(seqLen * kmerRegLenRatioEnd5);
				startKmerPos = 0;
				startBasePos = kmerSize;
				endBasePos = startBasePos + endKmerNum - 2;
			}else
			{
				endKmerNum = ceil(seqLen * kmerRegLenRatioEnd3);
				//startKmerPos = seqLen - endKmerNum - kmerSize + 1;
				//startBasePos = seqLen - endKmerNum + 1;
				startKmerPos = seqLen - kmerSize - 1 - ((endKmerNum-1)/graph->kmerSampleInterval)*graph->kmerSampleInterval + 1;
				startBasePos = startKmerPos + kmerSize;
				endBasePos = seqLen - 1;
			}


			// generate the kmer integer sequence
			if(generateKmerSeqIntFromReadset(kmerSeqInt, readseqInt, startKmerPos, entriesNum, baseNumLastEntry)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// get the hashcode and count the kmer occurrence
			kmerIntervalTmp = 0;
			pos = startKmerPos + 1;
			hashcode = kmerhashInt(kmerSeqInt);
			if(addKmer(hashcode, kmerSeqInt, rid, pos, graph)==FAILED)
			{
				printf("line=%d, In %s(), cannot add kmer for read %lu, error!\n", __LINE__, __func__, rid);
				return FAILED;
			}

			pos ++;
			entryRow = startBasePos >> 5;
			entryPos = startBasePos % 32;
			for(basePos=startBasePos; basePos<=endBasePos; basePos++, pos++)
			{
				// get the baseInt
				if(entryRow<entriesNum-1)
					baseInt = (readseqInt[entryRow] >> (62-2*entryPos)) & 3;
				else
					baseInt = (readseqInt[entryRow] >> (2*(baseNumLastEntry-1)-2*entryPos)) & 3;

				entryPos ++;
				if(entryPos==32)
				{
					entryRow ++;
					entryPos = 0;
				}

				// generate the kmer integer sequence
				if(entriesPerKmer>=2)
				{
					for(j=0; j<entriesPerKmer-2; j++)
					{
						kmerSeqInt[j] = (kmerSeqInt[j] << 2) | (kmerSeqInt[j+1] >> 62);
					}
					kmerSeqInt[entriesPerKmer-2] = (kmerSeqInt[entriesPerKmer-2] << 2) | (kmerSeqInt[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
				}
				kmerSeqInt[entriesPerKmer-1] = ((kmerSeqInt[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

				// get the hashcode and count the kmer occurrence
				kmerIntervalTmp ++;
				if(kmerIntervalTmp==graph->kmerSampleInterval)
				{
					hashcode = kmerhashInt(kmerSeqInt);
					if(addKmer(hashcode, kmerSeqInt, rid, pos, graph)==FAILED)
					{
						printf("line=%d, In %s(), cannot add kmer, error!\n", __LINE__, __func__);
						return FAILED;
					}
					kmerIntervalTmp = 0;
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Generate the kmer integer sequence from a read specified by a startPos (>=0).
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateKmerSeqIntFromReadset(uint64_t *seqInt, uint64_t *readseq, int32_t startReadPos, int32_t entriesNum, int32_t baseNumLastEntry)
{
	int32_t i, j, startEntryPos, remainedBaseNum, rightRemainedNum, entryRow, baseNumInEntry;
	uint64_t *readseqStart;

	for(i=0; i<entriesPerKmer; i++) seqInt[i] = 0;

	entryRow = startReadPos >> 5;
	readseqStart = readseq + entryRow;
	startEntryPos = startReadPos % 32;

	i = 0;
	j = 0;
	remainedBaseNum = kmerSize;
	while(remainedBaseNum>0)
	{
		// process first entry
		if(entryRow!=entriesNum-1)
			baseNumInEntry = 32;
		else
			baseNumInEntry = baseNumLastEntry;

		rightRemainedNum = baseNumInEntry - remainedBaseNum - startEntryPos;
		if(rightRemainedNum<0)
			rightRemainedNum = 0;

		seqInt[i] = (readseqStart[j] << (2*startEntryPos)) >> (2*(startEntryPos+rightRemainedNum));
		remainedBaseNum -= baseNumInEntry - rightRemainedNum - startEntryPos;

		if(startEntryPos>remainedBaseNum)
			startEntryPos = remainedBaseNum;

		j++;
		entryRow ++;

		// process second entry
		if(startEntryPos>0)
		{
			if(entryRow!=entriesNum-1)
				baseNumInEntry = 32;
			else
				baseNumInEntry = baseNumLastEntry;

			seqInt[i] = (seqInt[i] << (2*startEntryPos)) | (readseqStart[j] >> (2*(baseNumInEntry-startEntryPos)));
			remainedBaseNum -= startEntryPos;
		}

		i++;
	}

	return SUCCESSFUL;
}

/**
 * Get the read length from fasta file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMinReadLenFromFastaFiles(int *readLenInFile, char **readFilesInput, int readFileNum)
{
	int i, tmpReadLen;

	*readLenInFile = INT_MAX;
	for(i=0; i<readFileNum; i++)
	{
		if(getReadLenFromFasta(&tmpReadLen, readFilesInput[i])==FAILED)
		{
			//printf("line=%d, In %s(), cannot get readlenInFile, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(tmpReadLen<*readLenInFile)
			*readLenInFile = tmpReadLen;
	}

	return SUCCESSFUL;
}

/**
 * Get the read length from fastq file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMinReadLenFromFastqFiles(int *readLenInFile, char **readFilesInput, int readFileNum)
{
	int i, tmpReadLen;

	*readLenInFile = INT_MAX;
	for(i=0; i<readFileNum; i++)
	{
		if(getReadLenFromFastq(&tmpReadLen, readFilesInput[i])==FAILED)
		{
			//printf("line=%d, In %s(), cannot get readlenInFile, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(tmpReadLen<*readLenInFile)
			*readLenInFile = tmpReadLen;
	}

	return SUCCESSFUL;
}

/**
 * Get the read length from fasta file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadLenFromFasta(int *tmpReadLen, char *fastqFile)
{
	FILE *fpFasta;
	readBuf_t readBuf;
	char seq_data[5000];
	int i, returnCode, tmpLen, validFlag;

	fpFasta = fopen(fastqFile, "r");
	if(fpFasta==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fastqFile);
		return FAILED;
	}

	readBuf.seq = seq_data;

	tmpLen = 0;
	validFlag = YES;
	for(i=0; i<100; i++)
	{
		returnCode = getSingleReadFasta(fpFasta, &readBuf);
		if(returnCode==ERROR)
		{
			printf("line=%d, In %s(), cannot get single read, error!\n", __LINE__, __func__);
			return FAILED;
		}else if(returnCode==FAILED)
		{
			if(tmpLen>0 && tmpLen!=strlen(readBuf.seq))
			{
				validFlag = NO;
				break;
			}

			break;
		}

		if(i==0)
			tmpLen = strlen(readBuf.seq);
		else if(tmpLen!=strlen(readBuf.seq))
		{
			validFlag = NO;
			break;
		}
	}

	if(validFlag==YES)
		*tmpReadLen = tmpLen;
	else
	{
		*tmpReadLen = 0;

		fclose(fpFasta);
		fpFasta = NULL;

		printf("Reads in data sets should be in equal size.\n");
		return FAILED;
	}

	if(*tmpReadLen<kmerSize)
	{
		printf("readLen=%d < kmerSize=%d, error!\n", *tmpReadLen, kmerSize);
		return FAILED;
	}

	fclose(fpFasta);
	fpFasta = NULL;

	return SUCCESSFUL;
}

/**
 * Get the read length from fastq file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadLenFromFastq(int *tmpReadLen, char *fastqFile)
{
	FILE *fpFastq;
	readBuf_t readBuf;
	char seq_data[5000], qual_data[5000];
	int i, returnCode, tmpLen, validFlag;

	fpFastq = fopen(fastqFile, "r");
	if(fpFastq==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fastqFile);
		return FAILED;
	}

	readBuf.seq = seq_data;
	readBuf.qual = qual_data;

	tmpLen = 0;
	validFlag = YES;
	for(i=0; i<100; i++)
	{
		returnCode = getSingleReadFastq(fpFastq, &readBuf);
		if(returnCode==ERROR)
		{
			printf("line=%d, In %s(), cannot get single read, error!\n", __LINE__, __func__);
			return FAILED;
		}else if(returnCode==FAILED)
		{
			if(tmpLen>0 && tmpLen!=strlen(readBuf.seq))
			{
				validFlag = NO;
				break;
			}

			break;
		}

		if(i==0)
			tmpLen = strlen(readBuf.seq);
		else if(tmpLen!=strlen(readBuf.seq))
		{
			validFlag = NO;
			break;
		}
	}

	if(validFlag==YES)
		*tmpReadLen = tmpLen;
	else
	{
		*tmpReadLen = 0;

		fclose(fpFastq);
		fpFastq = NULL;

		printf("Reads in data sets should be in equal size.\n");
		return FAILED;
	}

	if(*tmpReadLen<kmerSize)
	{
		printf("readLen=%d < kmerSize=%d, error!\n", *tmpReadLen, kmerSize);
		return FAILED;
	}

	fclose(fpFastq);
	fpFastq = NULL;

	return SUCCESSFUL;
}

/**
 * Fill reads to reads buffer.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillReadsToBufFasta(FILE *fpReads, readBuf_t *pBuf, uint64_t *readsNum)
{
	uint64_t i;
	int returnCode;

	*readsNum = 0;
	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		returnCode = getSingleReadFasta(fpReads, pBuf+i);
		if(returnCode==FAILED)
			break;
		else if(returnCode==ERROR)
		{
			printf("line=%d, In %s(), cannot get single read, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	*readsNum = i;

	return SUCCESSFUL;
}

/**
 * Fill reads to reads buffer.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillReadsToBuf(FILE *fpReads, readBuf_t *pBuf, uint64_t *readsNum)
{
	uint64_t i;
	int returnCode;

	*readsNum = 0;
	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		returnCode = getSingleReadFastq(fpReads, pBuf+i);
		if(returnCode==FAILED)
			break;
		else if(returnCode==ERROR)
		{
			printf("line=%d, In %s(), cannot get single read, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	*readsNum = i;

	return SUCCESSFUL;
}

/**
 * get the single read from the given file.
 *  @return:
 *  	If succeed, return SUCCESSFUL;
 *  	else if the file end is reached, return FAILED;
 *  	otherwise, return ERROR.
 */
short getSingleReadFasta(FILE *fpPE, readBuf_t *pReadBuf)
{
	int i;
	//char qual_data[5000];	// read quality data which is encoded in ASCII
	char *readSeq;
	readSeq = pReadBuf->seq;

	char ch = fgetc(fpPE);
	if(feof(fpPE))// the file end is reached.
	{
		return FAILED;
	}

	while(ch!='\n') ch = fgetc(fpPE);

	while(!feof(fpPE))
	{
		i = 0;
		ch = fgetc(fpPE);
		while(ch!='>' && ch!=-1)
		{
			if(ch!='\n')
				readSeq[i++] = ch;
			ch = fgetc(fpPE);
		}
		readSeq[i] = '\0';

		pReadBuf->len = i;

		if(ch=='>')  //a read is read finished
		{
			break;
		}
	}

	return SUCCESSFUL;
}

/**
 * get the single read from the given file.
 *  @return:
 *  	If succeed, return SUCCESSFUL;
 *  	else if the file end is reached, return FAILED;
 *  	otherwise, return ERROR.
 */
short getSingleReadFastq(FILE *fpPE, readBuf_t *pReadBuf)
{
	int32_t i, line_index, tmpLen;
	//char qual_data[5000];	// read quality data which is encoded in ASCII
	char ch, *readSeq, *qual_data, *head_name;

	readSeq = pReadBuf->seq;
	qual_data = pReadBuf->qual;

	ch = fgetc(fpPE);
	if(feof(fpPE))// the file end is reached.
	{
		return FAILED;
	}

	line_index = 0;
	while(!feof(fpPE))
	{
		if(line_index==0)  //the sequence name line
		{
			ch = fgetc(fpPE);
			while(ch!='\n' && ch!=-1)
				ch = fgetc(fpPE);
		}else if(line_index==1)  //the sequence line
		{
			tmpLen = 0;
			ch = fgetc(fpPE);
			while(ch!='\n' && ch!=-1)
			{
				readSeq[tmpLen++] = ch;
				ch = fgetc(fpPE);
			}
			readSeq[tmpLen] = '\0';
		}else if(line_index==2)  //the sequence name line
		{
			ch = fgetc(fpPE);
			while(ch!='\n' && ch!=-1)
			{
				ch = fgetc(fpPE);
			}
		}else
		{
			i = 0;
			ch = fgetc(fpPE);
			while(ch!='\n'  && ch!=-1)
			{
				qual_data[i++] = ch;
				ch = fgetc(fpPE);
			}
			qual_data[i] = '\0';
		}
		line_index++;

		if(line_index==4)  //a read is read finished
		{
			pReadBuf->len = tmpLen;
			break;
		}
	}


	// check the return code
	if(line_index==4)
	{
		return SUCCESSFUL;
	}else
	{
		return ERROR;
	}
}


/**
 * Check whether the read contain unknown bases.
 *  @return:
 *  	If exists unknown bases, return YES; otherwise, return NO.
 */
short containUnknownBase(char *seq)
{
	char *p_str = seq;
	while(*p_str)
	{
		if(*p_str=='N' || *p_str=='.' || *p_str=='n')
			return YES;
		p_str++;
	}
	return NO;
}

/**
 * Get the unknown base count.
 *  @return:
 *  	Return the count.
 */
int32_t getUnknownBaseNum(char *seq)
{
	int32_t unknownBaseNum;
	char *p_str;

	p_str = seq;
	unknownBaseNum = 0;
	while(*p_str)
	{
		if(*p_str=='N' || *p_str=='.' || *p_str=='n')
			unknownBaseNum ++;
		p_str++;
	}

	return unknownBaseNum;
}

short singleQualSatisfied(char *qual_data)
{
	char *p_qual = qual_data;
	while(*p_qual)
	{
		if((*p_qual)-33 < singleBaseQualThres)
			return NO;
		p_qual ++;
	}

	return YES;
}

/**
 * Compute the average quality value at 3' end of a read.
 *  @return:
 *  	Return the average quality value.
 */
float calcAverQual3End(char *qual_data, int32_t seqLen)
{
	float sumQual;
	char *p_qual;
	int32_t baseNumEnd3;

	sumQual = 0;
	baseNumEnd3 = ceil(seqLen * QUAL_BASE_NUM_3End_FACTOR);
	p_qual = qual_data + seqLen - baseNumEnd3;
	while(*p_qual)
	{
		sumQual += (*p_qual) - 33;
		p_qual ++;
	}

	return sumQual/baseNumEnd3;
}

/**
 * Compute the average quality value at 5' end of a read.
 *  @return:
 *  	Return the average quality value.
 */
float calcAverQual5End(char *qual_data, int32_t seqLen)
{
	float sumQual;
	char *p_qual;
	int32_t i, baseNumEnd5;

	sumQual = 0;
	baseNumEnd5 = seqLen - ceil(seqLen * QUAL_BASE_NUM_3End_FACTOR);
	p_qual = qual_data;
	for(i=0; i<baseNumEnd5; i++)
	{
		sumQual += (*p_qual) - 33;
		p_qual ++;
	}

	return sumQual/baseNumEnd5;
}

/**
 * Compute the ratio of base 'A' of a read.
 *  @return:
 *  	Return the ratio.
 */
float getRatioBase(char *seq, char targetBase)
{
	int ratioA = 0, baseNum = 0;
	char *p_seq = seq;
	while(*p_seq)
	{
		if(*p_seq==targetBase)
			ratioA ++;
		baseNum ++;
		p_seq ++;
	}
	return (float)ratioA/baseNum;
}

/**
 * Generate the kmer interger sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateKmerSeqInt(uint64_t *seqInt, char *seq)
{
	int i, j, baseInt;

	for(i=0; i<entriesPerKmer; i++) seqInt[i] = 0;

	j = 0;
	i = 0;
	while(i<kmerSize)
	{
		switch(seq[i])
		{
			case 'A':
			case 'a': baseInt = 0; break;
			case 'C':
			case 'c': baseInt = 1; break;
			case 'G':
			case 'g': baseInt = 2; break;
			case 'T':
			case 't': baseInt = 3; break;
			default: printf("line=%d, In %s(), base %c, error!\n", __LINE__, __func__, seq[i]); return FAILED;
		}

		seqInt[j] = (seqInt[j] << 2) | baseInt;
		i ++;

		if(i%32==0)
			j ++;
	}

	return SUCCESSFUL;
}

/**
 * Compute the kmer hash code.
 *  @return:
 *    The hash code.
 */
uint64_t kmerhashInt(uint64_t *seqInt)
{
	//return *seqInt;

	uint64_t hashcode;
	int i, j;

	hashcode = 5381;
	for(i=0; i<entriesPerKmer-1; i++)
	{
		for(j=0; j<32; j++)
		{
			hashcode += (hashcode << 5) | ((seqInt[i] >> (62-2*j)) & 3);
		}
	}

	for(j=0; j<lastEntryBaseNum; j++)
	{
		hashcode += (hashcode << 5) | ((seqInt[entriesPerKmer-1] >> (2*(lastEntryBaseNum-j-1))) & 3);
	}

	//return (hashcode & 0x7FFFFFFF) % hashTableSize;
	return hashcode % hashTableSize;

}

/**
 * Initialize k-mer hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initgraph(graphtype **graph)
{
	*graph = (graphtype *) calloc(1, sizeof(graphtype));
	if((*graph)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize k-mer block
	if(initKmerBlockInGraph(*graph)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the k-mer blocks in k-mer hash table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize kmerseq block
	if(initKmerseqBlockInGraph(*graph)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the kmerseq blocks in k-mer hash table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Count the k-mer occurrences.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short countKmer(uint64_t hashcode, graphtype *graph)
{
	kmertype *kmer;

	kmer = getKmerByHash(hashcode, pKmerSeqTmpDoing, graph);
	if(!kmer)
	{
		// process k-mer blocks
		kmer = pKmerTmpDoing;

		kmer->kmerseqBlockID = pKmerBlockTmp->blockID;
		kmer->itemRowKmerseqBlock = pKmerBlockTmp->itemNum;
		kmer->arraysize = 1;
		kmer->nextKmerBlockID = graph->kmerHashtable[hashcode].kmerBlockID;
		kmer->nextItemRowKmerBlock = graph->kmerHashtable[hashcode].itemRowKmerBlock;
		graph->kmerHashtable[hashcode].kmerBlockID = kmer->kmerseqBlockID;
		graph->kmerHashtable[hashcode].itemRowKmerBlock = kmer->itemRowKmerseqBlock;

		pKmerBlockTmp->itemNum ++;
		graph->totalItemNumKmer ++;
		pKmerTmpDoing ++;

		if(pKmerBlockTmp->itemNum >= graph->maxItemNumPerKmerBlock)
		{
			// add new kmer block
			if(addNewBlockKmer(graph)==FAILED)
			{
				printf("line=%d, In %s(), cannot add new k-mer block, Error!\n", __LINE__, __func__);
				return FAILED;
			}
			pKmerBlockTmp = graph->kmerBlockArr + graph->blocksNumKmer - 1;
			pKmerTmpDoing = pKmerBlockTmp->kmerArr;
		}

		// process kmerseq blocks
		kmer->kmerseqBlockID = pKmerseqBlockTmp->blockID;
		kmer->itemRowKmerseqBlock = pKmerseqBlockTmp->itemNum;
		pKmerseqBlockTmp->itemNum ++;
		graph->totalItemNumKmerSeq ++;
		pKmerSeqTmpDoing += graph->entriesPerKmer;

		if(pKmerseqBlockTmp->itemNum >= graph->maxItemNumPerKmerSeqBlock)
		{
			// add new kmerseq block
			if(addNewBlockKmerSeq(graph)==FAILED)
			{
				printf("line=%d, In %s(), cannot add new kmerseq block, Error!\n", __LINE__, __func__);
				return FAILED;
			}
			pKmerseqBlockTmp = graph->kmerSeqBlockArr + graph->blocksNumKmerSeq - 1;
			pKmerSeqTmpDoing = pKmerseqBlockTmp->kmerSeqArr;
		}

	}else
	{
		kmer ->arraysize ++;
	}

	return SUCCESSFUL;
}

/**
 * Add a kmer to graph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addKmer(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, uint16_t rpos, graphtype *graph)
{
	kmertype *kmer;
	ridpostype *ridpos;

	kmer = getKmerByHash(hashcode, kmerSeqInt, graph);
	if(kmer)
	{
		ridpos = kmer->ppos + kmer->multiplicity;
		//ridpos->delsign = NO;
		ridpos->pos = rpos;
		ridpos->rid = rid;
		//ridpos->reserved = 0;

		kmer->multiplicity ++;

	}else
	{
		printf("line=%d, In %s(), can not get the k-mer %s. Error and exit.\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt));
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Delete the ridpos of a k-mer in a read.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short delKmer(char *str, uint64_t rid, unsigned short rpos, graphtype *graph)
{
	kmertype *kmer;
	ridpostype *ridpostable = NULL;
	uint64_t hashcode;

	// generate the kmer integer sequence
	if(generateKmerSeqInt(kmerSeqInt, str)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	hashcode = kmerhashInt(kmerSeqInt);
	kmer = getKmerByHash(hashcode, kmerSeqInt, graph);
	if(!kmer)
	{
		//printf("line=%d, In %s(), can not delete the kmer %s of (%d,%d), it does not exist in the graph. Error!\n", __LINE__, __func__, str, rid, rpos);
		return FAILED;
	}

	ridpostable = kmer->ppos;
	if(!ridpostable)
	{
		printf("line=%d, In %s(), can not delete the kmer %s of (%lu, %u), the position array does not exist in the graph. Error!\n", __LINE__, __func__, str, rid, rpos);
		return FAILED;
	}

	if(kmer->multiplicity==0)
	{
		//printf("line=%d, In %s(), can not delete the kmer %s of (%lu,%u), all of its positions have been deleted.\n", __LINE__, __func__, str, rid, rpos);
		return FAILED;
	}

	int index = findDirectIndex(rid, rpos, ridpostable, kmer->arraysize);
	if(index>=0)
	{

		if(ridpostable[index].delsign==0)
		{
			ridpostable[index].delsign = 1;
			kmer->multiplicity--;
		}else
		{
			printf("line=%d, In %s(), can not delete the kmer %s (%lu, %u), it has been deleted in the graph.\n", __LINE__, __func__, str, rid, rpos);
			return FAILED;
		}

	}else
	{
		//printf("line=%d, In %s(), can not delete the kmer %s s(%d, %d), it does not exist in the graph.\n", __LINE__, __func__, str, rid, rpos);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Delete the ridpos of a k-mer in a read using k-mer hash code.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short delKmerByHash(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, uint32_t rpos, graphtype *graph)
{
	ridpostype *ridpostable = NULL;
	kmertype *kmer;
	int32_t row;

	kmer = getKmerByHash(hashcode, kmerSeqInt, graph);
	if(!kmer)
	{
		printf("line=%d, In %s(), can not delete the kmer %s of (%lu,%u), it does not exist in the graph. Error!\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	ridpostable = kmer->ppos;
	if(!ridpostable)
	{
		printf("line=%d, In %s(), can not delete the kmer %s of (%lu,%u), the position array does not exist in the graph. Error!\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	if(kmer->multiplicity==0)
	{
		printf("line=%d, In %s(), can not delete the kmer %s of (%lu, %u), all of its positions have been deleted.\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	row = findDirectIndex(rid, rpos, ridpostable, kmer->arraysize);
	if(row>=0)
	{
		if(ridpostable[row].delsign==NO)
		{
			ridpostable[row].delsign = YES;
			kmer->multiplicity--;
		}else
		{
			printf("line=%d, In %s(), can not delete the kmer %s of (%lu, %u), it has been deleted in the graph.\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
			return FAILED;
		}

	}else
	{
		printf("line=%d, In %s(), can not delete the kmer %s of (%lu, %u), it does not exist in the graph.\n", __LINE__, __func__,  getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Recover the k-mer in a read.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short recoverKmerByHash(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, uint16_t rpos, graphtype *graph)
{
	ridpostype *ridpostable;
	kmertype *kmer;
	int32_t row;

	hashcode = kmerhashInt(kmerSeqInt);
	kmer = getKmerByHash(hashcode, kmerSeqInt, graph);
	if(!kmer)
	{
		printf("line=%d, In %s(), can not recover the kmer %s of (%lu,%u), it does not exist in the graph, kmer=NULL. Error!\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	ridpostable = kmer->ppos;
	if(!ridpostable)
	{
		printf("line=%d, In %s(), can not recover the kmer %s of (%lu,%u), the position array does not exist in the graph. Error!\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	if(kmer->multiplicity==kmer->arraysize)
	{
		printf("line=%d, In %s(), can not recover the kmer %s of (%lu,%u), multiplicity==arraysize. Error!\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	row = findDirectIndex(rid, rpos, ridpostable, kmer->arraysize);
	if(row>=0)
	{
		if(ridpostable[row].delsign==NO)
		{ // it has already been recovered
			printf("line=%d, In %s(), can not recover the kmer %s of (%lu, %u), it has been recovered in the graph.\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
			return FAILED;
		}
		ridpostable[row].delsign = NO;
		ridpostable[row].reserved = 0;
		kmer->multiplicity++;

	}else
	{ // does not find the ridpos
		printf("line=%d, In %s(), can not recover the kmer %s of (%lu,%u), it does not exist in the graph.\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get the row (starts from 0) of a ridpos in the ridpos table.
 *  @return:
 *  	If finds, return the row; otherwise, return -1.
 */
int findDirectIndex(uint64_t rid, uint16_t rpos, ridpostype *ridpostable, int posNum)
{
	int index = -1;
	int startIndex = findStartIndex(rid, ridpostable, posNum);
	if(startIndex>=0)
	{
		index = getExectIndex(rid, rpos, startIndex, ridpostable, posNum);
	}
	return index;
}

/**
 * Get the start row (starts from 0) of a read in the ridpos table.
 *  @return:
 *  	If finds, return the row; otherwise, return -1.
 */
int findStartIndex(uint64_t rid, ridpostype *rid_pos_table, int posNum)
{
	int left, right, middle, existFlag = 0;
	left = 0;
	right = posNum - 1;
	while(left<=right)
	{
		middle = (left+right)/2;
		if(rid==rid_pos_table[middle].rid)
		{
			existFlag = 1;
			break;
		}
		if(rid>rid_pos_table[middle].rid)
			left = middle+1;
		else
			right = middle-1;
	}

	if(existFlag)
	{
		while(middle>0 && rid_pos_table[middle-1].rid==rid)
			middle--;
		return middle;
	}

	return -1;
}

/**
 * Get the exact row (starts from 0) of a ridpos from a start row in the ridpos table.
 *  @return:
 *  	If finds, return the row; otherwise, return -1.
 */
int getExectIndex(uint64_t rid, uint16_t rpos, int startIndex, ridpostype *ridpostable, int posNum)
{
	int i = startIndex, exectIndex = -1;
	while(i<posNum && ridpostable[i].rid==rid)
	{
		if(ridpostable[i].pos==rpos)
		{
			exectIndex = i;
			break;
		}
		i++;
	}

	return exectIndex;
}


inline kmertype *getKmerByBase(char *str, graphtype *graph)
{
	uint64_t hashcode;
	kmertype *kmer;

	// generate the kmer integer sequence
	if(generateKmerSeqInt(kmerSeqInt, str)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
		return NULL;
	}

	hashcode = kmerhashInt(kmerSeqInt);
	kmer = getKmerByHash(hashcode, kmerSeqInt, graph);

	return kmer;
}

inline kmertype *getKmer(uint64_t *kmerSeqInt, graphtype *graph)
{
	uint64_t hashcode;

	hashcode = kmerhashInt(kmerSeqInt);

	return getKmerByHash(hashcode, kmerSeqInt, graph);
}

inline kmertype *getKmerByHash(uint64_t hashvalue, uint64_t *kmerSeqInt, graphtype *graph)
{
	kmertype *kmer;
	uint64_t *kmerseq;
	kmerHashBucket_t *pKmerBucket;

	pKmerBucket = graph->kmerHashtable + hashvalue;
	if(pKmerBucket->kmerBlockID>0)
	{
		kmer = graph->kmerBlockArr[pKmerBucket->kmerBlockID-1].kmerArr + pKmerBucket->itemRowKmerBlock;
		while(kmer)
		{
			kmerseq = graph->kmerSeqBlockArr[kmer->kmerseqBlockID-1].kmerSeqArr + kmer->itemRowKmerseqBlock * graph->entriesPerKmer;
			if(identicalKmerSeq(kmerSeqInt, kmerseq)==YES)
				break;

			if(kmer->nextKmerBlockID>0)
				kmer = graph->kmerBlockArr[kmer->nextKmerBlockID-1].kmerArr + kmer->nextItemRowKmerBlock;
			else
				kmer = NULL;
		}
	}else
	{
		kmer = NULL;
	}

	return kmer;
}

/**
 * Check whether the two sequence is identical.
 *  @return:
 *  	If identical, return YES; otherwise, return NO.
 */
short identicalKmerSeq(uint64_t *kmerSeqInt1, uint64_t *kmerSeqInt2)
{
	int i;
	for(i=0; i<entriesPerKmer; i++)
	{
		if(kmerSeqInt1[i] != kmerSeqInt2[i])
			return NO;
	}

	return YES;
}

/**
 * Get the kmer bases from integer.
 */
char *getKmerBaseByInt(uint64_t *kmerSeqInt)
{
	int i, j, k, baseInt;

	k = 0;
	for(i=0; i<entriesPerKmer-1; i++)
	{
		for(j=0; j<32; j++)
		{
			baseInt = (kmerSeqInt[i] >> (62-2*j)) & 3;
			switch(baseInt)
			{
				case 0: baseSeq[k]='A'; break;
				case 1: baseSeq[k]='C'; break;
				case 2: baseSeq[k]='G'; break;
				case 3: baseSeq[k]='T'; break;
			}
			k ++;
		}
	}

	for(j=0; j<lastEntryBaseNum; j++)
	{
		baseInt = (kmerSeqInt[entriesPerKmer-1] >> (2*(lastEntryBaseNum-j-1))) & 3;
		switch(baseInt)
		{
			case 0: baseSeq[k]='A'; break;
			case 1: baseSeq[k]='C'; break;
			case 2: baseSeq[k]='G'; break;
			case 3: baseSeq[k]='T'; break;
		}
		k ++;
	}

	baseSeq[k] = '\0';

	return baseSeq;
}


/**
 * Get the kmer integer sequence of a kmer by its kmer integer sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReverseKmerseqInt(uint64_t *kmerSeqIntRev, uint64_t *kmerSeqInt)
{
	int i, j, baseInt;
	int entryIndexRev, baseNum;

	for(i=0; i<entriesPerKmer; i++) kmerSeqIntRev[i] = 0;

	entryIndexRev = baseNum = 0;
	for(j=0; j<lastEntryBaseNum; j++)
	{
		baseInt = (~(kmerSeqInt[entriesPerKmer-1] >> (j*2))) & 3;
		kmerSeqIntRev[entryIndexRev] = (kmerSeqIntRev[entryIndexRev] << 2) | baseInt;
		baseNum ++;
		if(baseNum==32)
		{
			entryIndexRev ++;
			baseNum = 0;
		}
	}
	for(i=entriesPerKmer-2; i>=0; i--)
	{
		for(j=0; j<32; j++)
		{
			baseInt = (~(kmerSeqInt[i] >> (2*j))) & 3;
			kmerSeqIntRev[entryIndexRev] = (kmerSeqIntRev[entryIndexRev] << 2) | baseInt;
			baseNum ++;
			if(baseNum==32)
			{
				entryIndexRev ++;
				baseNum = 0;
			}
		}
	}

	// ########################### Debug information ########################
	if(kmerSize<32 && baseNum!=kmerSize)
	{
		printf("line=%d, In %s(), baseNum=%d != kmerSize=%d, error!\n", __LINE__, __func__, baseNum, kmerSize);
		return FAILED;
	}
	// ########################### Debug information ########################

	return SUCCESSFUL;
}

/**
 * The fast version of getReverseKmer().
 */
kmertype *getReverseKmer(uint64_t *kmerSeqIntRev, uint64_t *kmerSeqInt, graphtype *graph)
{
	uint64_t hashcode;

	if(getReverseKmerseqInt(kmerSeqIntRev, kmerSeqInt)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the reverse kmer, error!\n", __LINE__, __func__);
		exit(1);
	}

	hashcode = kmerhashInt(kmerSeqIntRev);

	return getKmerByHash(hashcode, kmerSeqIntRev, graph);
}


/**
 *  Output the k-mer hash table to file.
 *   File format:
 *   	(1) readLen, kmerSize, hashTableSize, hashTableSizeReadseq, pairedMode,
 *   		PEGivenType, meanSizeInsert, standardDev;
 *   	(2) graph node;
 *   	(3) readset node;
 *   	(4) read block head nodes;
 *   	(5) read blocks;
 *   	(6) readseq block head nodes;
 *   	(7) readseq blocks;
 *   	(8) readseq hash table;
 *   	(9) readseq hash item head nodes;
 *   	(10) readseq hash item blocks;
 *   	(11) kmer hash table;
 *   	(12) kmer block head nodes;
 *   	(13) kmer blocks;
 *   	(14) kmerseq block head nodes;
 *   	(15) kmerseq blocks;
 *   	(16) ridpos block head nodes;
 *   	(17) ridposBlocks;
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputGraphToFile(char *graphFile, graphtype *graph)
{
	printf("Begin to output the k-mer hash table to file ...\n");

	uint64_t i, tmp[6];
	FILE *fpGraph;
	double insSdev[2];
	readBlock_t *readBlockArr;
	readseqBlock_t *readseqBlockArr;
	readseqHashItemBlock_t *readseqHashItemBlockArr;
	kmerBlock_t *kmerBlockArr;
	kmerseqBlock_t *kmerseqBlockArr;
	ridposBlock_t *ridposBlockArr;

	struct timeval tpstart, tpend;
	double timeused_file;
	gettimeofday(&tpstart, NULL);


	fpGraph = fopen(graphFile, "wb");
	if(fpGraph==NULL)
	{
		printf("line=%d, In %s(), can not open the file [ %s ], Error!\n", __LINE__, __func__, graphFile);
		return FAILED;
	}

	tmp[0] = readLen;
	tmp[1] = kmerSize;
	tmp[2] = hashTableSize;
	tmp[3] = hashTableSizeReadseq;
	tmp[4] = pairedMode;
	tmp[5] = PEGivenType;

	// part (1)-1: readLen, kmerSize, hashTableSize, hashTableSizeReadseq, pairedMode, PEGivenType
	if(fwrite(tmp, sizeof(uint64_t), 6, fpGraph)!=6)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (1)-2: meanSizeInsert, standardDev
	insSdev[0] = meanSizeInsert;
	insSdev[1] = standardDev;
	if(fwrite(insSdev, sizeof(double), 2, fpGraph)!=2)
	{
		printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
		return FAILED;
	}


	// part (2): graph node;
	if(fwrite(graph, sizeof(graphtype), 1, fpGraph)!=1)
	{
		printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (3): readset node
	if(fwrite(graph->readSet, sizeof(readSet_t), 1, fpGraph)!=1)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (4): read block head nodes
	if(fwrite(graph->readSet->readBlockArr, sizeof(readBlock_t), graph->readSet->blocksNumRead, fpGraph)!=graph->readSet->blocksNumRead)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (5): read blocks
#if (DEBUG_OUTPUT==YES)
	printf("outputting read blocks ...\n");
#endif
	readBlockArr = graph->readSet->readBlockArr;
	for(i=0; i<graph->readSet->blocksNumRead; i++)
	{
		if(fwrite(readBlockArr[i].readArr, graph->readSet->bytesPerRead, readBlockArr[i].itemNum, fpGraph)!=readBlockArr[i].itemNum)
		{
			printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// part (6): readseq block head nodes
	if(fwrite(graph->readSet->readseqBlockArr, sizeof(readseqBlock_t), graph->readSet->blocksNumReadseq, fpGraph)!=graph->readSet->blocksNumReadseq)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (7): readseq blocks
#if (DEBUG_OUTPUT==YES)
	printf("outputting readseq blocks ...\n");
#endif
	readseqBlockArr = graph->readSet->readseqBlockArr;
	for(i=0; i<graph->readSet->blocksNumReadseq; i++)
	{
		if(fwrite(readseqBlockArr[i].readseqArr, graph->readSet->bytesPerEntryReadseq, readseqBlockArr[i].rowsNum, fpGraph)!=readseqBlockArr[i].rowsNum)
		{
			printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
			return FAILED;
		}
	}


	if(graph->readSet->readseqHashtable)
	{
		// part (8): readseq hash table
#if (DEBUG_OUTPUT==YES)
		printf("outputting readseq hash table ...\n");
#endif
		if(fwrite(graph->readSet->readseqHashtable, sizeof(readseqHashBucket_t), graph->readSet->hashTableSizeReadseq, fpGraph)!=graph->readSet->hashTableSizeReadseq)
		{
			printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
			return FAILED;
		}

		// part (9): readseq hash item head nodes
		if(fwrite(graph->readSet->readseqHashItemBlockArr, sizeof(readseqHashItemBlock_t), graph->readSet->blocksNumReadseqHashItem, fpGraph)!=graph->readSet->blocksNumReadseqHashItem)
		{
			printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
			return FAILED;
		}

		// part (10): readseq hash item blocks
#if (DEBUG_OUTPUT==YES)
		printf("outputting readseq hash item blocks ...\n");
#endif
		readseqHashItemBlockArr = graph->readSet->readseqHashItemBlockArr;
		for(i=0; i<graph->readSet->blocksNumReadseqHashItem; i++)
		{
			if(fwrite(readseqHashItemBlockArr[i].readseqHashItemArr, graph->readSet->bytesPerReadseqHashItem, readseqHashItemBlockArr[i].itemNum, fpGraph)!=readseqHashItemBlockArr[i].itemNum)
			{
				printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	// part (11): kmer hash table
#if (DEBUG_OUTPUT==YES)
	printf("outputting k-mer hash table ...\n");
#endif
	if(fwrite(graph->kmerHashtable, sizeof(kmerHashBucket_t), graph->hashTableSize, fpGraph)!=graph->hashTableSize)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (12): kmer block head nodes
	if(fwrite(graph->kmerBlockArr, sizeof(kmerBlock_t), graph->blocksNumKmer, fpGraph)!=graph->blocksNumKmer)
	{
		printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (13): kmer blocks
#if (DEBUG_OUTPUT==YES)
	printf("outputting kmer blocks ...\n");
#endif
	kmerBlockArr = graph->kmerBlockArr;
	for(i=0; i<graph->blocksNumKmer; i++)
	{
		if(fwrite(kmerBlockArr[i].kmerArr, graph->bytesPerKmer, kmerBlockArr[i].itemNum, fpGraph)!=kmerBlockArr[i].itemNum)
		{
			printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
			return FAILED;
		}
	}


	// part (14): kmerseq block head nodes
	if(fwrite(graph->kmerSeqBlockArr, sizeof(kmerseqBlock_t), graph->blocksNumKmerSeq, fpGraph)!=graph->blocksNumKmerSeq)
	{
		printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (15): kmerseq blocks
#if (DEBUG_OUTPUT==YES)
	printf("outputting kmerseq blocks ...\n");
#endif
	kmerseqBlockArr = graph->kmerSeqBlockArr;
	for(i=0; i<graph->blocksNumKmerSeq; i++)
	{
		if(fwrite(kmerseqBlockArr[i].kmerSeqArr, graph->bytesPerKmerseq, kmerseqBlockArr[i].itemNum, fpGraph)!=kmerseqBlockArr[i].itemNum)
		{
			printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// part (16): ridpos block head nodes
	if(fwrite(graph->ridposBlockArr, sizeof(ridposBlock_t), graph->blocksNumRidpos, fpGraph)!=graph->blocksNumRidpos)
	{
		printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (17): ridpos blocks
#if (DEBUG_OUTPUT==YES)
	printf("outputting ridpos blocks ...\n");
#endif
	ridposBlockArr = graph->ridposBlockArr;
	for(i=0; i<graph->blocksNumRidpos; i++)
	{
		if(fwrite(ridposBlockArr[i].ridposArr, graph->bytesPerRidpos, ridposBlockArr[i].itemNum, fpGraph)!=ridposBlockArr[i].itemNum)
		{
			printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	fclose(fpGraph);
	fpGraph = NULL;

	printf("End outputting the k-mer hash table to file.\n");

	gettimeofday(&tpend, NULL);
	timeused_file = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Outputting file used time: %.2f seconds.\n", timeused_file);

	return SUCCESSFUL;
}

/**
 *  Load the graph to memory.
 *   File format:
 *   	(1) readLen, kmerSize, hashTableSize, hashTableSizeReadseq, pairedMode,
 *   		PEGivenType, meanSizeInsert, standardDev;
 *   	(2) graph node;
 *   	(3) readset node;
 *   	(4) read block head nodes;
 *   	(5) read blocks;
 *   	(6) readseq block head nodes;
 *   	(7) readseq blocks;
 *   	(8) readseq hash table;
 *   	(9) readseq hash item head nodes;
 *   	(10) readseq hash item blocks;
 *   	(11) kmer hash table;
 *   	(12) kmer block head nodes;
 *   	(13) kmer blocks;
 *   	(14) kmerseq block head nodes;
 *   	(15) kmerseq blocks;
 *   	(16) ridpos block head nodes;
 *   	(17) ridposBlocks;
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short loadGraph(graphtype **graph, char *graphFile)
{
	printf("\nBegin loading the k-mer hash table ...\n");

	FILE *fpGraph;
	uint64_t tmp[6], i, j, sum, itemNumKmerBlock;
	double insSdev[2];
	kmertype *kmer;

	ridposBlock_t *ridposBlockTmp;
	ridpostype *ridposTmp;

	struct timeval tpstart, tpend;
	double timeused_file;
	gettimeofday(&tpstart, NULL);


	fpGraph = fopen(graphFile, "rb");
	if(fpGraph==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, graphFile);
		return FAILED;
	}

	// part (1)-1: readLen, kmerSize, hashTableSize, hashTableSizeReadseq, pairedMode, PEGivenType
	if(fread(tmp, sizeof(uint64_t), 6, fpGraph)!=6)
	{
		printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
		return FAILED;
	}

	readLen = tmp[0];
	kmerSize = tmp[1];
	hashTableSize = tmp[2];
	hashTableSizeReadseq = tmp[3];
	pairedMode = tmp[4];
	PEGivenType = tmp[5];

	// part (1)-2: meanSizeInsert, standardDev
	if(fread(insSdev, sizeof(double), 2, fpGraph)!=2)
	{
		printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
		return FAILED;
	}
	meanSizeInsert = insSdev[0];
	standardDev = insSdev[1];


	// part (2): graph node
	*graph = (graphtype *) malloc(sizeof(graphtype));
	if((*graph)==NULL)
	{
		printf("line=%d, In %s(), cannot malloc the graph, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fread(*graph, sizeof(graphtype), 1, fpGraph)!=1)
	{
		printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
		return FAILED;
	}
	(*graph)->kmerHashtable = NULL;
	(*graph)->kmerSeqBlockArr = NULL;
	//(*graph)->maxBlocksNumKmer = (*graph)->blocksNumKmer;
	//(*graph)->maxBlocksNumKmerSeq = (*graph)->blocksNumKmerSeq;
	//(*graph)->maxBlocksNumRidpos = (*graph)->blocksNumRidpos;

	// part (3): readset node
	(*graph)->readSet = (readSet_t *) malloc (sizeof(readSet_t));
	if((*graph)->readSet==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fread((*graph)->readSet, sizeof(readSet_t), 1, fpGraph)!=1)
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}
	//(*graph)->readSet->maxBlocksNumRead = (*graph)->readSet->blocksNumRead;
	//(*graph)->readSet->maxBlocksNumReadseq = (*graph)->readSet->blocksNumReadseq;
	//(*graph)->readSet->maxBlocksNumReadseqHashItem = (*graph)->readSet->blocksNumReadseqHashItem;

	// part (4): read block head nodes
	//(*graph)->readSet->readBlockArr = (readBlock_t *) malloc ((*graph)->readSet->blocksNumRead * sizeof(readBlock_t));
	(*graph)->readSet->readBlockArr = (readBlock_t *) malloc ((*graph)->readSet->maxBlocksNumRead * sizeof(readBlock_t));
	if((*graph)->readSet->readBlockArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fread((*graph)->readSet->readBlockArr, sizeof(readBlock_t), (*graph)->readSet->blocksNumRead, fpGraph)!=(*graph)->readSet->blocksNumRead)
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (5): read blocks
#if (DEBUG_OUTPUT==YES)
	printf("loading read blocks ...\n");
#endif
	for(i=0; i<(*graph)->readSet->blocksNumRead; i++)
	{
		//(*graph)->readSet->readBlockArr[i].readArr = (read_t *) malloc((*graph)->readSet->readBlockArr[i].itemNum * (*graph)->readSet->bytesPerRead);
		(*graph)->readSet->readBlockArr[i].readArr = (read_t *) malloc((*graph)->readSet->maxItemNumPerReadBlock * (*graph)->readSet->bytesPerRead);
		if((*graph)->readSet->readBlockArr[i].readArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(fread((*graph)->readSet->readBlockArr[i].readArr, (*graph)->readSet->bytesPerRead, (*graph)->readSet->readBlockArr[i].itemNum, fpGraph)!=(*graph)->readSet->readBlockArr[i].itemNum)
		{
			printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// part (6): readseq block head nodes
	//(*graph)->readSet->readseqBlockArr = (readseqBlock_t *) malloc ((*graph)->readSet->blocksNumReadseq * sizeof(readseqBlock_t));
	(*graph)->readSet->readseqBlockArr = (readseqBlock_t *) malloc ((*graph)->readSet->maxBlocksNumReadseq * sizeof(readseqBlock_t));
	if((*graph)->readSet->readseqBlockArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fread((*graph)->readSet->readseqBlockArr, sizeof(readseqBlock_t), (*graph)->readSet->blocksNumReadseq, fpGraph)!=(*graph)->readSet->blocksNumReadseq)
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (7): readseq blocks
#if (DEBUG_OUTPUT==YES)
	printf("loading readseq blocks ...\n");
#endif
	for(i=0; i<(*graph)->readSet->blocksNumReadseq; i++)
	{
		//(*graph)->readSet->readseqBlockArr[i].readseqArr = (uint64_t *) malloc((*graph)->readSet->readseqBlockArr[i].rowsNum * (*graph)->readSet->bytesPerEntryReadseq);
		(*graph)->readSet->readseqBlockArr[i].readseqArr = (uint64_t *) malloc((*graph)->readSet->maxEntryNumReadseqBlock * (*graph)->readSet->bytesPerEntryReadseq);
		if((*graph)->readSet->readseqBlockArr[i].readseqArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(fread((*graph)->readSet->readseqBlockArr[i].readseqArr, (*graph)->readSet->bytesPerEntryReadseq, (*graph)->readSet->readseqBlockArr[i].rowsNum, fpGraph)!=(*graph)->readSet->readseqBlockArr[i].rowsNum)
		{
			printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	if((*graph)->readSet->readseqHashtable)
	{
		// part (8): readseq hash table
#if (DEBUG_OUTPUT==YES)
		printf("loading readseq hash table ...\n");
#endif
		(*graph)->readSet->readseqHashtable = (readseqHashBucket_t *) malloc ((*graph)->readSet->hashTableSizeReadseq * sizeof(readseqHashBucket_t));
		if((*graph)->readSet->readseqHashtable==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(fread((*graph)->readSet->readseqHashtable, sizeof(readseqHashBucket_t), (*graph)->readSet->hashTableSizeReadseq, fpGraph)!=(*graph)->readSet->hashTableSizeReadseq)
		{
			printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
			return FAILED;
		}

		// part (9): readseq hash item head nodes
		//(*graph)->readSet->readseqHashItemBlockArr = (readseqHashItemBlock_t *) malloc ((*graph)->readSet->blocksNumReadseqHashItem * sizeof(readseqHashItemBlock_t));
		(*graph)->readSet->readseqHashItemBlockArr = (readseqHashItemBlock_t *) malloc ((*graph)->readSet->maxBlocksNumReadseqHashItem * sizeof(readseqHashItemBlock_t));
		if((*graph)->readSet->readseqHashItemBlockArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(fread((*graph)->readSet->readseqHashItemBlockArr, sizeof(readseqHashItemBlock_t), (*graph)->readSet->blocksNumReadseqHashItem, fpGraph)!=(*graph)->readSet->blocksNumReadseqHashItem)
		{
			printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
			return FAILED;
		}

		// part (10): readseq hash item blocks
#if (DEBUG_OUTPUT==YES)
		printf("loading readseq hash item blocks ...\n");
#endif
		for(i=0; i<(*graph)->readSet->blocksNumReadseqHashItem; i++)
		{
			//(*graph)->readSet->readseqHashItemBlockArr[i].readseqHashItemArr = (readseqHashItem_t *) malloc((*graph)->readSet->readseqHashItemBlockArr[i].itemNum * (*graph)->readSet->bytesPerReadseqHashItem);
			(*graph)->readSet->readseqHashItemBlockArr[i].readseqHashItemArr = (readseqHashItem_t *) malloc((*graph)->readSet->maxItemNumPerReadseqHashItemBlock * (*graph)->readSet->bytesPerReadseqHashItem);
			if((*graph)->readSet->readseqHashItemBlockArr[i].readseqHashItemArr==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(fread((*graph)->readSet->readseqHashItemBlockArr[i].readseqHashItemArr, (*graph)->readSet->bytesPerReadseqHashItem, (*graph)->readSet->readseqHashItemBlockArr[i].itemNum, fpGraph)!=(*graph)->readSet->readseqHashItemBlockArr[i].itemNum)
			{
				printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	// part (11): kmer hash table
#if (DEBUG_OUTPUT==YES)
	printf("loading k-mer hash table ...\n");
#endif
	(*graph)->kmerHashtable = (kmerHashBucket_t *) malloc ((*graph)->hashTableSize * sizeof(kmerHashBucket_t));
	if((*graph)->kmerHashtable==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fread((*graph)->kmerHashtable, sizeof(kmerHashBucket_t), (*graph)->hashTableSize, fpGraph)!=(*graph)->hashTableSize)
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (12): kmer block head nodes
	//(*graph)->kmerBlockArr = (kmerBlock_t *) malloc ((*graph)->blocksNumKmer * sizeof(kmerBlock_t));
	(*graph)->kmerBlockArr = (kmerBlock_t *) malloc ((*graph)->maxBlocksNumKmer * sizeof(kmerBlock_t));
	if((*graph)->kmerBlockArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fread((*graph)->kmerBlockArr, sizeof(kmerBlock_t), (*graph)->blocksNumKmer, fpGraph)!=(*graph)->blocksNumKmer)
	{
		printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (13): kmer blocks
#if (DEBUG_OUTPUT==YES)
	printf("loading kmer blocks ...\n");
#endif
	for(i=0; i<(*graph)->blocksNumKmer; i++)
	{
		//(*graph)->kmerBlockArr[i].kmerArr = (kmertype *) malloc((*graph)->kmerBlockArr[i].itemNum * (*graph)->bytesPerKmer);
		(*graph)->kmerBlockArr[i].kmerArr = (kmertype *) malloc((*graph)->maxItemNumPerKmerBlock * (*graph)->bytesPerKmer);
		if((*graph)->kmerBlockArr[i].kmerArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(fread((*graph)->kmerBlockArr[i].kmerArr, (*graph)->bytesPerKmer, (*graph)->kmerBlockArr[i].itemNum, fpGraph)!=(*graph)->kmerBlockArr[i].itemNum)
		{
			printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
			return FAILED;
		}
	}


	// part (14): kmerseq block head nodes
	//(*graph)->kmerSeqBlockArr = (kmerseqBlock_t *) malloc ((*graph)->blocksNumKmerSeq * sizeof(kmerseqBlock_t));
	(*graph)->kmerSeqBlockArr = (kmerseqBlock_t *) malloc ((*graph)->maxBlocksNumKmerSeq * sizeof(kmerseqBlock_t));
	if((*graph)->kmerSeqBlockArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fread((*graph)->kmerSeqBlockArr, sizeof(kmerseqBlock_t), (*graph)->blocksNumKmerSeq, fpGraph)!=(*graph)->blocksNumKmerSeq)
	{
		printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (15): kmerseq blocks
#if (DEBUG_OUTPUT==YES)
	printf("loading kmerseq blocks ...\n");
#endif
	for(i=0; i<(*graph)->blocksNumKmerSeq; i++)
	{
		//(*graph)->kmerSeqBlockArr[i].kmerSeqArr = (uint64_t *) malloc((*graph)->kmerSeqBlockArr[i].itemNum * (*graph)->bytesPerKmerseq);
		(*graph)->kmerSeqBlockArr[i].kmerSeqArr = (uint64_t *) malloc((*graph)->maxItemNumPerKmerSeqBlock * (*graph)->bytesPerKmerseq);
		if((*graph)->kmerSeqBlockArr[i].kmerSeqArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(fread((*graph)->kmerSeqBlockArr[i].kmerSeqArr, (*graph)->bytesPerKmerseq, (*graph)->kmerSeqBlockArr[i].itemNum, fpGraph)!=(*graph)->kmerSeqBlockArr[i].itemNum)
		{
			printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
			return FAILED;
		}
	}


	// part (16): ridpos block head nodes
	//(*graph)->ridposBlockArr = (ridposBlock_t *) malloc ((*graph)->blocksNumRidpos * sizeof(ridposBlock_t));
	(*graph)->ridposBlockArr = (ridposBlock_t *) malloc ((*graph)->maxBlocksNumRidpos * sizeof(ridposBlock_t));
	if((*graph)->ridposBlockArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fread((*graph)->ridposBlockArr, sizeof(ridposBlock_t), (*graph)->blocksNumRidpos, fpGraph)!=(*graph)->blocksNumRidpos)
	{
		printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (17): k-mer ridpos blocks
#if (DEBUG_OUTPUT==YES)
	printf("loading k-mer ridpos blocks ...\n");
#endif
	for(i=0; i<(*graph)->blocksNumRidpos; i++)
	{
		//(*graph)->ridposBlockArr[i].ridposArr = (ridpostype *) malloc((*graph)->ridposBlockArr[i].itemNum * (*graph)->bytesPerRidpos);
		(*graph)->ridposBlockArr[i].ridposArr = (ridpostype *) malloc((*graph)->maxItemNumPerRidposBlock * (*graph)->bytesPerRidpos);
		if((*graph)->ridposBlockArr[i].ridposArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(fread((*graph)->ridposBlockArr[i].ridposArr, (*graph)->bytesPerRidpos, (*graph)->ridposBlockArr[i].itemNum, fpGraph)!=(*graph)->ridposBlockArr[i].itemNum)
		{
			printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// part (18): initialize ridpos regions
	// set start ridpos regions of k-mers
#if (DEBUG_OUTPUT==YES)
	printf("initialize k-mer ridpos regions ...\n");
#endif
	ridposBlockTmp = (*graph)->ridposBlockArr;
	ridposTmp = ridposBlockTmp->ridposArr;
	sum = 0;
	for(i=0; i<(*graph)->blocksNumKmer; i++)
	{
		kmer = (*graph)->kmerBlockArr[i].kmerArr;
		itemNumKmerBlock = (*graph)->kmerBlockArr[i].itemNum;
		for(j=0; j<itemNumKmerBlock; j++)
		{
			kmer->ppos = ridposTmp;
			ridposTmp += kmer->arraysize;

			sum += kmer->arraysize;
			if(sum == ridposBlockTmp->itemNum)
			{
				ridposBlockTmp ++;
				ridposTmp = ridposBlockTmp->ridposArr;
				sum = 0;
			}

			kmer ++;
		}
	}

	fclose(fpGraph);
	fpGraph = NULL;

	printf("End loading the k-mer hash table.\n");

	gettimeofday(&tpend, NULL);
	timeused_file = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Loading file used time: %.2f seconds.\n", timeused_file);

	return SUCCESSFUL;
}

/**
 *  Load the readSet from graph file.
 *   File format:
 *   	(1) readLen, kmerSize, hashTableSize, hashTableSizeReadseq, pairedMode,
 *   		PEGivenType, meanSizeInsert, standardDev;
 *   	(2) graph node;
 *   	(3) readset node;
 *   	(4) read block head nodes;
 *   	(5) read blocks;
 *   	(6) readseq block head nodes;
 *   	(7) readseq blocks;
 *   	(8) readseq hash table;
 *   	(9) readseq hash item head nodes;
 *   	(10) readseq hash item blocks;
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short loadReadSetFromGraphFile(readSet_t **readSet, char *graphFile)
{
	printf("\nBegin loading the readSet from the k-mer hash table ...\n");

	FILE *fpGraph;
	uint64_t tmp[6], i, j, sum, itemNumKmerBlock;
	double insSdev[2];
	graphtype graphTmp;

	struct timeval tpstart, tpend;
	double timeused_file;
	gettimeofday(&tpstart, NULL);


	fpGraph = fopen(graphFile, "rb");
	if(fpGraph==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, graphFile);
		return FAILED;
	}

	// part (1)-1: readLen, kmerSize, hashTableSize, hashTableSizeReadseq, pairedMode, PEGivenType
	if(fread(tmp, sizeof(uint64_t), 6, fpGraph)!=6)
	{
		printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (1)-2: meanSizeInsert, standardDev
	if(fread(insSdev, sizeof(double), 2, fpGraph)!=2)
	{
		printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (2): graph node
	if(fread(&graphTmp, sizeof(graphtype), 1, fpGraph)!=1)
	{
		printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (3): readset node
	*readSet = (readSet_t *) malloc (sizeof(readSet_t));
	if((*readSet)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fread(*readSet, sizeof(readSet_t), 1, fpGraph)!=1)
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (4): read block head nodes
	(*readSet)->readBlockArr = (readBlock_t *) malloc ((*readSet)->maxBlocksNumRead * sizeof(readBlock_t));
	if((*readSet)->readBlockArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fread((*readSet)->readBlockArr, sizeof(readBlock_t), (*readSet)->blocksNumRead, fpGraph)!=(*readSet)->blocksNumRead)
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (5): read blocks
#if (DEBUG_OUTPUT==YES)
	printf("loading read blocks ...\n");
#endif
	for(i=0; i<(*readSet)->blocksNumRead; i++)
	{
		(*readSet)->readBlockArr[i].readArr = (read_t *) malloc((*readSet)->maxItemNumPerReadBlock * (*readSet)->bytesPerRead);
		if((*readSet)->readBlockArr[i].readArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(fread((*readSet)->readBlockArr[i].readArr, (*readSet)->bytesPerRead, (*readSet)->readBlockArr[i].itemNum, fpGraph)!=(*readSet)->readBlockArr[i].itemNum)
		{
			printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// part (6): readseq block head nodes
	(*readSet)->readseqBlockArr = (readseqBlock_t *) malloc ((*readSet)->maxBlocksNumReadseq * sizeof(readseqBlock_t));
	if((*readSet)->readseqBlockArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fread((*readSet)->readseqBlockArr, sizeof(readseqBlock_t), (*readSet)->blocksNumReadseq, fpGraph)!=(*readSet)->blocksNumReadseq)
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}

	// part (7): readseq blocks
#if (DEBUG_OUTPUT==YES)
	printf("loading readseq blocks ...\n");
#endif
	for(i=0; i<(*readSet)->blocksNumReadseq; i++)
	{
		(*readSet)->readseqBlockArr[i].readseqArr = (uint64_t *) malloc((*readSet)->maxEntryNumReadseqBlock * (*readSet)->bytesPerEntryReadseq);
		if((*readSet)->readseqBlockArr[i].readseqArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(fread((*readSet)->readseqBlockArr[i].readseqArr, (*readSet)->bytesPerEntryReadseq, (*readSet)->readseqBlockArr[i].rowsNum, fpGraph)!=(*readSet)->readseqBlockArr[i].rowsNum)
		{
			printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	if((*readSet)->readseqHashtable)
	{
		// part (8): readseq hash table
#if (DEBUG_OUTPUT==YES)
		printf("loading readseq hash table ...\n");
#endif
		(*readSet)->readseqHashtable = (readseqHashBucket_t *) malloc ((*readSet)->hashTableSizeReadseq * sizeof(readseqHashBucket_t));
		if((*readSet)->readseqHashtable==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(fread((*readSet)->readseqHashtable, sizeof(readseqHashBucket_t), (*readSet)->hashTableSizeReadseq, fpGraph)!=(*readSet)->hashTableSizeReadseq)
		{
			printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
			return FAILED;
		}

		// part (9): readseq hash item head nodes
		(*readSet)->readseqHashItemBlockArr = (readseqHashItemBlock_t *) malloc ((*readSet)->maxBlocksNumReadseqHashItem * sizeof(readseqHashItemBlock_t));
		if((*readSet)->readseqHashItemBlockArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(fread((*readSet)->readseqHashItemBlockArr, sizeof(readseqHashItemBlock_t), (*readSet)->blocksNumReadseqHashItem, fpGraph)!=(*readSet)->blocksNumReadseqHashItem)
		{
			printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
			return FAILED;
		}

		// part (10): readseq hash item blocks
#if (DEBUG_OUTPUT==YES)
		printf("loading readseq hash item blocks ...\n");
#endif
		for(i=0; i<(*readSet)->blocksNumReadseqHashItem; i++)
		{
			(*readSet)->readseqHashItemBlockArr[i].readseqHashItemArr = (readseqHashItem_t *) malloc((*readSet)->maxItemNumPerReadseqHashItemBlock * (*readSet)->bytesPerReadseqHashItem);
			if((*readSet)->readseqHashItemBlockArr[i].readseqHashItemArr==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(fread((*readSet)->readseqHashItemBlockArr[i].readseqHashItemArr, (*readSet)->bytesPerReadseqHashItem, (*readSet)->readseqHashItemBlockArr[i].itemNum, fpGraph)!=(*readSet)->readseqHashItemBlockArr[i].itemNum)
			{
				printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	fclose(fpGraph);
	fpGraph = NULL;

	printf("End loading the readSet from the k-mer hash table.\n");

	gettimeofday(&tpend, NULL);
	timeused_file = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Loading file used time: %.2f seconds.\n", timeused_file);

	return SUCCESSFUL;
}


/**
 * Get the global parameters from graph.
 * 		(1) readLen, averReadLenInFileSample, kmerSize, hashTableSize, hashTableSizeReadseq, pairedMode,
 * 			PEGivenType, meanSizeInsert, standardDev;
 *
 *  @return:
 *  	If succeeds, return SUCCESFUL; otherwise, return FAILED.
 */
short GlobalParasFromGraph(int32_t *readLenPara, int32_t *kmerSizePara, uint64_t *hashTableSizePara,
		uint64_t *hashTableSizeReadseqPara, int32_t *pairedModePara, int32_t *PEGivenTypePara,
		double *meanSizeInsertPara, double *standardDevPara, char *graphFileName)
{
	FILE *fpGraph;
	uint64_t tmp[6];
	double insSdev[2];

	fpGraph = fopen(graphFileName, "rb");
	if(fpGraph==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, graphFileName);
		return FAILED;
	}

	// get the item number of the kmer array and ridpos array, respectively
	if(fread(tmp, sizeof(uint64_t), 6, fpGraph)!=6)
	{
		printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
		return FAILED;
	}

	*readLenPara = tmp[0];
	*kmerSizePara = tmp[1];
	*hashTableSizePara = tmp[2];
	*hashTableSizeReadseqPara = tmp[3];
	*pairedModePara = tmp[4];
	*PEGivenTypePara = tmp[5];

	// get meanSizeInsert, standardDev
	if(fread(insSdev, sizeof(double), 2, fpGraph)!=2)
	{
		printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
		return FAILED;
	}

	*meanSizeInsertPara = insSdev[0];
	*standardDevPara = insSdev[1];

	fclose(fpGraph);
	fpGraph = NULL;

	return SUCCESSFUL;
}

/**
 * Update the insert size and sdev in k-mer hash table.
 * 		(1) PEGivenType, meanSizeInsert, standardDev.
 *
 *  @return:
 *  	If succeeds, return SUCCESFUL; otherwise, return FAILED.
 */
short updateInsertSizeAndSdevInGraph(int64_t PEGivenType, double meanSizeInsert, double standardDev, const char *graphFileName)
{
	FILE *fpGraph;
	double insSdev[2];

	fpGraph = fopen(graphFile, "rb+");
	if(fpGraph==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, graphFile);
		return FAILED;
	}

	// seek to the right position in the file
	if(fseek(fpGraph, 5*sizeof(uint64_t), SEEK_SET)==-1)
	{
		printf("line=%d, In %s(), fseek error!\n", __LINE__, __func__);
		return FAILED;
	}

	// write: PEGivenType
	if(fwrite(&PEGivenType, sizeof(uint64_t), 1, fpGraph)!=1)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	// write: meanSizeInsert, standardDev
	insSdev[0] = meanSizeInsert;
	insSdev[1] = standardDev;
	if(fwrite(insSdev, sizeof(double), 2, fpGraph)!=2)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	fclose(fpGraph);
	fpGraph = NULL;

	return SUCCESSFUL;
}

/**
 * Release the k-mer hash table.
 */
short releaseGraph(graphtype **graph)
{
	uint64_t i;

	if((*graph)==NULL)
		return SUCCESSFUL;

	// free read set
	releaseReadset(&((*graph)->readSet));

	// free k-mer blocks
	if((*graph)->kmerBlockArr)
	{
		for(i=0; i<(*graph)->blocksNumKmer; i++)
			free((*graph)->kmerBlockArr[i].kmerArr);
		(*graph)->blocksNumKmer = 0;
		free((*graph)->kmerBlockArr);
		(*graph)->kmerBlockArr = NULL;
	}
	if((*graph)->kmerHashtable)
	{
		free((*graph)->kmerHashtable);
		(*graph)->kmerHashtable = NULL;
	}

	// free kmerseq blocks
	if((*graph)->kmerSeqBlockArr)
	{
		for(i=0; i<(*graph)->blocksNumKmerSeq; i++)
			free((*graph)->kmerSeqBlockArr[i].kmerSeqArr);
		(*graph)->blocksNumKmerSeq = 0;
		free((*graph)->kmerSeqBlockArr);
		(*graph)->kmerSeqBlockArr = NULL;
	}

	// free ridpos blocks
	if((*graph)->ridposBlockArr)
	{
		for(i=0; i<(*graph)->blocksNumRidpos; i++)
			free((*graph)->ridposBlockArr[i].ridposArr);
		(*graph)->blocksNumRidpos = 0;
		free((*graph)->ridposBlockArr);
		(*graph)->ridposBlockArr = NULL;
	}

	// free graph node
	free(*graph);
	*graph = NULL;

	return SUCCESSFUL;
}

/**
 * Clean k-mer information in k-mer hash table, and only remains the read set information.
 */
short cleanKmerInfoInGraph(graphtype **graph)
{
	uint64_t i;

	if((*graph)==NULL)
	{
		printf("line=%d, In %s(), graph=NULL, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// free k-mer blocks
	if((*graph)->kmerBlockArr)
	{
		for(i=0; i<(*graph)->blocksNumKmer; i++)
			free((*graph)->kmerBlockArr[i].kmerArr);
		(*graph)->blocksNumKmer = 0;
		free((*graph)->kmerBlockArr);
		(*graph)->kmerBlockArr = NULL;
	}
	if((*graph)->kmerHashtable)
	{
		free((*graph)->kmerHashtable);
		(*graph)->kmerHashtable = NULL;
	}

	// free kmerseq blocks
	if((*graph)->kmerSeqBlockArr)
	{
		for(i=0; i<(*graph)->blocksNumKmerSeq; i++)
			free((*graph)->kmerSeqBlockArr[i].kmerSeqArr);
		(*graph)->blocksNumKmerSeq = 0;
		free((*graph)->kmerSeqBlockArr);
		(*graph)->kmerSeqBlockArr = NULL;
	}

	// free ridpos blocks
	if((*graph)->ridposBlockArr)
	{
		for(i=0; i<(*graph)->blocksNumRidpos; i++)
			free((*graph)->ridposBlockArr[i].ridposArr);
		(*graph)->blocksNumRidpos = 0;
		free((*graph)->ridposBlockArr);
		(*graph)->ridposBlockArr = NULL;
	}

	// free graph node
	free(*graph);

	return SUCCESSFUL;
}

/**
 * Reset the de Bruijn graph.
 *  @return:
 *  	If succeeds, return SUCCESFUL; otherwise, return FAILED.
 */
short resetGraph(graphtype *graph)
{
	readBlock_t *readBlockArr;
	read_t *pRead;
	readseqBlock_t *readseqBlockArr;
	uint64_t rid, *readseq;
	int32_t i, j, blocksNumRead, itemNumReadBlock;

	readBlockArr = graph->readSet->readBlockArr;
	readseqBlockArr = graph->readSet->readseqBlockArr;
	blocksNumRead = graph->readSet->blocksNumRead;


	// reset the read nodes
	rid = 0;
	for(i=0; i<blocksNumRead; i++)
	{
		pRead = readBlockArr[i].readArr;
		itemNumReadBlock = readBlockArr[i].itemNum;
		for(j=0; j<itemNumReadBlock; j++)
		{
			rid ++;
			if(pRead->successFlag==YES)
			{
				readseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
				if(recoverKmersInRead(rid, readseq, pRead->seqlen, graph)==FAILED)
				{
					printf("line=%d, In %s(), cannot recover read %ld, error!\n", __LINE__, __func__, rid);
					return FAILED;
				}

				pRead->successFlag = NO;
			}

			pRead ++;
		}
	}

	return SUCCESSFUL;
}

/**
 * Recover k-mers in a read.
 *  @return:
 *  	If succeeds, return SUCCESFUL; otherwise, return FAILED.
 */
short recoverKmersInRead(int64_t rid, uint64_t *readseqInt, int32_t seqLen, graphtype *graph)
{
	int32_t i, j, pos, basePos, startKmerPos, startBasePos, endBasePos, endKmerNum;
	int32_t baseInt, entriesNum, baseNumLastEntry, entryRow, entryPos, kmerIntervalTmp;
	uint64_t hashcode;

	entriesNum = ((seqLen - 1) >> 5) + 1;
	baseNumLastEntry = ((seqLen - 1) % 32) + 1;

	if(ceil(seqLen * kmerRegLenRatioEnd5) + ceil(seqLen * kmerRegLenRatioEnd3) >= seqLen - kmerSize + 1)
	{ // the 5' + 3' >= kmerNum
		// generate the k-mer integer sequence
		if(generateKmerSeqIntFromReadset(kmerSeqInt, readseqInt, 0, entriesNum, baseNumLastEntry)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// get the hashcode and count the kmer occurrence
		kmerIntervalTmp = 0;
		hashcode = kmerhashInt(kmerSeqInt);
		if(recoverKmerByHash(hashcode, kmerSeqInt, rid, 1, graph)==FAILED)
		{
			printf("line=%d, In %s(), cannot recover kmer for read %lu, error!\n", __LINE__, __func__, rid);
			return FAILED;
		}

		pos = 2;
		entryRow = kmerSize >> 5;
		entryPos = kmerSize % 32;
		for(basePos=kmerSize; basePos<seqLen; basePos++, pos++)
		{
			// get the baseInt
			if(entryRow<entriesNum-1)
				baseInt = (readseqInt[entryRow] >> (62-2*entryPos)) & 3;
			else
				baseInt = (readseqInt[entryRow] >> (2*(baseNumLastEntry-1)-2*entryPos)) & 3;

			entryPos ++;
			if(entryPos==32)
			{
				entryRow ++;
				entryPos = 0;
			}

			// generate the kmer integer sequence
			if(entriesPerKmer>=2)
			{
				for(j=0; j<entriesPerKmer-2; j++)
				{
					kmerSeqInt[j] = (kmerSeqInt[j] << 2) | (kmerSeqInt[j+1] >> 62);
				}
				kmerSeqInt[entriesPerKmer-2] = (kmerSeqInt[entriesPerKmer-2] << 2) | (kmerSeqInt[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
			}
			kmerSeqInt[entriesPerKmer-1] = ((kmerSeqInt[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

			// get the hash code and recover the k-mer occurrence
			kmerIntervalTmp ++;
			if(kmerIntervalTmp==graph->kmerSampleInterval || basePos==seqLen-1)
			{
				hashcode = kmerhashInt(kmerSeqInt);
				if(recoverKmerByHash(hashcode, kmerSeqInt, rid, pos, graph)==FAILED)
				{
					printf("line=%d, In %s(), cannot recover kmer, error!\n", __LINE__, __func__);
					return FAILED;
				}
				kmerIntervalTmp = 0;
			}
		}
	}else
	{
		for(i=0; i<2; i++)
		{
			if(i==0)
			{
				endKmerNum = ceil(seqLen * kmerRegLenRatioEnd5);
				startKmerPos = 0;
				startBasePos = kmerSize;
				endBasePos = startBasePos + endKmerNum - 2;
			}else
			{
				endKmerNum = ceil(seqLen * kmerRegLenRatioEnd3);
				//startKmerPos = seqLen - endKmerNum - kmerSize + 1;
				//startBasePos = seqLen - endKmerNum + 1;
				startKmerPos = seqLen - kmerSize - 1 - ((endKmerNum-1)/graph->kmerSampleInterval)*graph->kmerSampleInterval + 1;
				startBasePos = startKmerPos + kmerSize;
				endBasePos = seqLen - 1;
			}

			// generate the kmer integer sequence
			if(generateKmerSeqIntFromReadset(kmerSeqInt, readseqInt, startKmerPos, entriesNum, baseNumLastEntry)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// get the hash code and recover the k-mer occurrence
			kmerIntervalTmp = 0;
			pos = startKmerPos + 1;
			hashcode = kmerhashInt(kmerSeqInt);
			if(recoverKmerByHash(hashcode, kmerSeqInt, rid, pos, graph)==FAILED)
			{
				printf("line=%d, In %s(), cannot recover kmer for read %lu, error!\n", __LINE__, __func__, rid);
				return FAILED;
			}

			pos ++;
			entryRow = startBasePos >> 5;
			entryPos = startBasePos % 32;
			for(basePos=startBasePos; basePos<=endBasePos; basePos++, pos++)
			{
				// get the baseInt
				if(entryRow<entriesNum-1)
					baseInt = (readseqInt[entryRow] >> (62-2*entryPos)) & 3;
				else
					baseInt = (readseqInt[entryRow] >> (2*(baseNumLastEntry-1)-2*entryPos)) & 3;

				entryPos ++;
				if(entryPos==32)
				{
					entryRow ++;
					entryPos = 0;
				}

				// generate the kmer integer sequence
				if(entriesPerKmer>=2)
				{
					for(j=0; j<entriesPerKmer-2; j++)
					{
						kmerSeqInt[j] = (kmerSeqInt[j] << 2) | (kmerSeqInt[j+1] >> 62);
					}
					kmerSeqInt[entriesPerKmer-2] = (kmerSeqInt[entriesPerKmer-2] << 2) | (kmerSeqInt[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
				}
				kmerSeqInt[entriesPerKmer-1] = ((kmerSeqInt[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

				// get the hashcode and recover the kmer occurrence
				kmerIntervalTmp ++;
				if(kmerIntervalTmp==graph->kmerSampleInterval)
				{
					hashcode = kmerhashInt(kmerSeqInt);
					if(recoverKmerByHash(hashcode, kmerSeqInt, rid, pos, graph)==FAILED)
					{
						printf("line=%d, In %s(), cannot add kmer, error!\n", __LINE__, __func__);
						return FAILED;
					}
					kmerIntervalTmp = 0;
				}
			}
		}
	}

	return SUCCESSFUL;
}
