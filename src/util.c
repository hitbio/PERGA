/*
 * util.c
 *
 *  Created on: May 29, 2010
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Output all contig nodes.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputContig(contigtype *contigArray, int64_t contigNodesNum)
{
	successRead_t *ridposorientation = NULL;
	int32_t i, j, num;
	char tmp_base, orient;

	for(i=0; i<contigNodesNum; i++)
	{
		switch(contigArray[i].base)
		{
			case 0: tmp_base = 'A'; break;
			case 1: tmp_base = 'C'; break;
			case 2: tmp_base = 'G'; break;
			case 3: tmp_base = 'T'; break;
			default: printf("line=%d, In %s(), error base int: %d\n", __LINE__, __func__, contigArray[i].base); return FAILED;
		}

		printf("%d\t%c\t%d", contigArray[i].index, tmp_base, contigArray[i].ridposnum);
		ridposorientation = contigArray[i].pridposorientation;
		num = contigArray[i].ridposnum;
		for(j=0; j<num; j++)
		{
			switch(ridposorientation[j].orientation)
			{
				case ORIENTATION_PLUS: orient = '+'; break;
				case ORIENTATION_MINUS: orient = '-'; break;
				default: printf("line=%d, In %s(), j=%d, invalid read orientation=%d, error!\n", __LINE__, __func__, j, (int32_t)ridposorientation[j].orientation);
			}
			printf("\t(%lu,%d,%d,%c)", (int64_t)ridposorientation[j].rid, (int32_t)ridposorientation[j].pos, ridposorientation[j].matchnum, orient);
		}
		printf("\n");
	}
	printf("There are %ld nodes.\n", contigNodesNum);

	return SUCCESSFUL;
}

/**
 * Output the contig nodes near 5' end nodeNum bases.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputContigEnd5(contigtype *contigArray, int64_t contigNodesNum, int endNum)
{
	successRead_t *ridposorientation;
	int32_t i, j, num;
	char tmp_base, orient;

	for(i=0; i<contigNodesNum; i++)
	{
		switch(contigArray[i].base)
		{
			case 0: tmp_base = 'A'; break;
			case 1: tmp_base = 'C'; break;
			case 2: tmp_base = 'G'; break;
			case 3: tmp_base = 'T'; break;
			default: printf("line=%d, In %s(), error base int: %d\n", __LINE__, __func__, contigArray[i].base); return FAILED;
		}

		printf("%d\t%c\t%d", contigArray[i].index, tmp_base, contigArray[i].ridposnum);
		ridposorientation = contigArray[i].pridposorientation;
		num = contigArray[i].ridposnum;
		for(j=0; j<num; j++)
		{
			switch(ridposorientation[j].orientation)
			{
				case ORIENTATION_PLUS: orient = '+'; break;
				case ORIENTATION_MINUS: orient = '-'; break;
				default: printf("line=%d, In %s(), j=%d, invalid read orientation=%d, error!\n", __LINE__, __func__, j, (int32_t)ridposorientation[j].orientation);
			}
			printf("\t(%lu,%d,%d,%c)", (uint64_t)ridposorientation[i].rid, (int32_t)ridposorientation[i].pos, ridposorientation[i].matchnum, orient);
		}
		printf("\n");

		if(i+1>=endNum)
			break;
	}
	printf("There are %d nodes.\n", i+1);

	return SUCCESSFUL;
}

/**
 * Output the contig nodes near 3' end nodeNum bases.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputContigEnd3(contigtype *contigArray, int64_t contigNodesNum, int endNum)
{
	successRead_t *ridposorientation;
	int32_t i, j, num, nodeNum;
	char tmp_base, orient;

	nodeNum = 0;
	for(i=contigNodesNum-endNum; i<contigNodesNum; i++)
	{
		switch(contigArray[i].base)
		{
			case 0: tmp_base = 'A'; break;
			case 1: tmp_base = 'C'; break;
			case 2: tmp_base = 'G'; break;
			case 3: tmp_base = 'T'; break;
			default: printf("line=%d, In %s(), error base int: %d\n", __LINE__, __func__, contigArray[i].base); return FAILED;
		}

		printf("%d\t%c\t%d", contigArray[i].index, tmp_base, contigArray[i].ridposnum);
		ridposorientation = contigArray[i].pridposorientation;
		num = contigArray[i].ridposnum;
		for(j=0; j<num; j++)
		{
			switch(ridposorientation[j].orientation)
			{
				case ORIENTATION_PLUS: orient = '+'; break;
				case ORIENTATION_MINUS: orient = '-'; break;
				default: printf("line=%d, In %s(), j=%d, invalid read orientation=%d, error!\n", __LINE__, __func__, j, (int32_t)ridposorientation[j].orientation);
			}
			printf("\t(%lu,%d,%d,%c)", (uint64_t)ridposorientation[j].rid, (int32_t)ridposorientation[j].pos, ridposorientation[j].matchnum, orient);
		}
		printf("\n");

		nodeNum ++;
	}
	printf("There are %d nodes.\n", nodeNum);

	return SUCCESSFUL;
}

/**
 * Output the contig bases near 3' end nodeNum bases.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputContigBaseEnd3(contigtype *contigArray, int64_t contigNodesNum, int endNum)
{
	int i, nodeNum;
	char tmp_base, baseArray[endNum+1];

	nodeNum = 0;
	for(i=contigNodesNum-endNum; i<contigNodesNum; i++)
	{
		switch(contigArray[i].base)
		{
			case 0: tmp_base = 'A'; break;
			case 1: tmp_base = 'C'; break;
			case 2: tmp_base = 'G'; break;
			case 3: tmp_base = 'T'; break;
			default: printf("line=%d, In %s(), error base int: %d\n", __LINE__, __func__, contigArray[i].base); return FAILED;
		}
		baseArray[nodeNum++] = tmp_base;
	}
	baseArray[nodeNum] = '\0';
	printf("%s, len=%d\n", baseArray, nodeNum);

	return SUCCESSFUL;
}

/**
 * Output the undeleted k-mer rispos.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputUndelKmerpos(graphtype *graph)
{
	uint64_t i, j, k, count, posNum;
	kmertype *kmer;
	ridpostype *ridpostable;
	int64_t itemNumKmer;


	printf("Begin checking undelKmerpos:\n");

	count = 0;
	for(i=0; i<graph->blocksNumKmer; i++)
	{
		kmer = graph->kmerBlockArr[i].kmerArr;
		itemNumKmer = graph->kmerBlockArr[i].itemNum;
		for(j=0; j<itemNumKmer; j++)
		{
			posNum = kmer->arraysize;
			ridpostable = kmer->ppos;
			if(kmer->multiplicity>0)
			{
				for(k=0; k<posNum; k++)
				{
					if(ridpostable->delsign==0)
					{
						//printf("(%d,%d) ", ridpostable->rid, ridpostable->pos);
						count++;
					}
					ridpostable++;
				}
				//printf("\n");
			}
		}
	}

	printf("checking undelKmerpos finished, the count=%lu\n", count);

	return SUCCESSFUL;
}

/*
 * Output the remained Kmer and its multiplicity and arraySize.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputRemainedKmers(graphtype *graph)
{
	FILE *fpRemainedKmer;
	int64_t i, count;
	kmertype *kmer;
	kmerHashBucket_t *pKmerBucket;
	uint64_t *kmerseq;

	printf("Begin outputting  remained k-mers:\n");

	fpRemainedKmer = fopen("../remainedKmers.txt", "w");
	if(fpRemainedKmer==NULL)
	{
		printf("line=%d, In %s(), cannot open file [../remainedKmers.txt], error!\n", __LINE__, __func__);
		return FAILED;
	}

	count = 0;
	for(i=0; i<hashTableSize; i++)
	{
		pKmerBucket = graph->kmerHashtable + i;
		if(pKmerBucket->kmerBlockID>0)
		{
			kmer = graph->kmerBlockArr[pKmerBucket->kmerBlockID-1].kmerArr + pKmerBucket->itemRowKmerBlock;
			while(kmer)
			{
				if(kmer->multiplicity>0)
				{
					kmerseq = graph->kmerSeqBlockArr[kmer->kmerseqBlockID-1].kmerSeqArr + kmer->itemRowKmerseqBlock * graph->entriesPerKmer;
					fprintf(fpRemainedKmer, "%ld\t%s\t%u\t%u\n", i, getKmerBaseByInt(kmerseq), kmer->multiplicity, kmer->arraysize);

					count ++;
				}

				if(kmer->nextKmerBlockID>0)
					kmer = graph->kmerBlockArr[kmer->nextKmerBlockID-1].kmerArr + kmer->nextItemRowKmerBlock;
				else
					kmer = NULL;
			}
		}
	}

	fclose(fpRemainedKmer);

	printf("checking undelKmerpos finished, the count=%lu\n", count);

	return SUCCESSFUL;
}

/**
 * Output the reads in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputReadsInDecisionTable(assemblingreadtype *decisionTable, int readsNum)
{
	assemblingreadtype *this_assemblingRead;
	int32_t i, lockedNum, appearNum, matedNum;
	char orient;

	lockedNum = 0;
	appearNum = 0;
	matedNum = 0;
	this_assemblingRead = decisionTable;
	for(i=0; i<readsNum; i++)
	{
		switch(this_assemblingRead->orientation)
		{
			case ORIENTATION_PLUS: orient = '+'; break;
			case ORIENTATION_MINUS: orient = '-'; break;
			default: printf("line=%d, In %s(), i=%d, invalid read orientation=%d, error!\n", __LINE__, __func__, i, (int32_t)this_assemblingRead->orientation);
		}

		printf(	"rid=%lu, firstBasePos=%d, firstContigPos=%d, orientation=%c, status=%d, "
				"matchBaseNum=%d, unmatchBaseNum=%d, basePos=%d, lastMatchedBasePos=%d, "
				"seqlen=%d, entriesNumReadseq=%d, baseNumLastEentryReadseq=%d, "
				"kmerNumEnd5=%d, kmerNumEnd3=%d, "
				"locked=%d, matedFlag=%d\n",
				(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstBasePos, this_assemblingRead->firstContigPos, orient, this_assemblingRead->status,
				this_assemblingRead->matchBaseNum, this_assemblingRead->unmatchBaseNum, this_assemblingRead->basePos, this_assemblingRead->lastMatchedBasePos,
				this_assemblingRead->seqlen, this_assemblingRead->entriesNumReadseq, this_assemblingRead->baseNumLastEentryReadseq,
				this_assemblingRead->kmerNumEnd5, this_assemblingRead->kmerNumEnd3,
				this_assemblingRead->locked, this_assemblingRead->matedFlag);
		printf("\t+: %s\n", getReadBaseByInt(this_assemblingRead->readseq, this_assemblingRead->seqlen));
		printf("\t-: %s\n", getReverseReadBaseByInt(this_assemblingRead->readseq, this_assemblingRead->seqlen));


		if(this_assemblingRead->locked==1)
			lockedNum++;

		if(this_assemblingRead->basePos==this_assemblingRead->lastMatchedBasePos)
			appearNum ++;

		if(this_assemblingRead->matedFlag==YES)
			matedNum ++;

		this_assemblingRead++;
	}
	printf("The readsNum=%d, appearNum=%d, locked=%d, matedNum=%d\n", readsNum, appearNum, lockedNum, matedNum);
	return SUCCESSFUL;
}

short outputReadsInDecisionTableToFile(char *outfile, assemblingreadtype *decisionTable, int readsNum)
{
	FILE *fpReads;
	assemblingreadtype *this_assemblingRead;
	int32_t i, lockedNum, appearNum, matedNum;
	char orient;

	fpReads = fopen(outfile, "w");
	if(fpReads==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, outfile);
		return FAILED;
	}

	lockedNum = 0;
	appearNum = 0;
	matedNum = 0;
	this_assemblingRead = decisionTable;
	for(i=0; i<readsNum; i++)
	{
		switch(this_assemblingRead->orientation)
		{
			case ORIENTATION_PLUS: orient = '+'; break;
			case ORIENTATION_MINUS: orient = '-'; break;
			default: printf("line=%d, In %s(), i=%d, invalid read orientation=%d, error!\n", __LINE__, __func__, i, (int32_t)this_assemblingRead->orientation);
		}

		fprintf(fpReads, "rid=%lu, firstBasePos=%d, firstContigPos=%d, orientation=%c, status=%d, "
				"matchBaseNum=%d, unmatchBaseNum=%d, basePos=%d, lastMatchedBasePos=%d, "
				"seqlen=%d, entriesNumReadseq=%d, baseNumLastEentryReadseq=%d, "
				"kmerNumEnd5=%d, kmerNumEnd3=%d, "
				"locked=%d, matedFlag=%d\n",
				(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstBasePos, this_assemblingRead->firstContigPos, orient, this_assemblingRead->status,
				this_assemblingRead->matchBaseNum, this_assemblingRead->unmatchBaseNum, this_assemblingRead->basePos, this_assemblingRead->lastMatchedBasePos,
				this_assemblingRead->seqlen, this_assemblingRead->entriesNumReadseq, this_assemblingRead->baseNumLastEentryReadseq,
				this_assemblingRead->kmerNumEnd5, this_assemblingRead->kmerNumEnd3,
				this_assemblingRead->locked, this_assemblingRead->matedFlag);
		fprintf(fpReads, "\t+: %s\n", getReadBaseByInt(this_assemblingRead->readseq, this_assemblingRead->seqlen));
		fprintf(fpReads, "\t-: %s\n", getReverseReadBaseByInt(this_assemblingRead->readseq, this_assemblingRead->seqlen));

		if(this_assemblingRead->locked==1)
			lockedNum++;

		if(this_assemblingRead->basePos==this_assemblingRead->lastMatchedBasePos)
			appearNum ++;

		if(this_assemblingRead->matedFlag==YES)
			matedNum ++;

		this_assemblingRead++;
	}
	fprintf(fpReads, "readsNum=%d, appearNum=%d, locked=%d, matedNum=%d\n", readsNum, appearNum, lockedNum, matedNum);

	fclose(fpReads);
	fpReads = NULL;

	return SUCCESSFUL;
}

short outputMatedReadsInDecisionTableToFile(char *outfile, assemblingreadtype *decisionTable, int readsNum)
{
	FILE *fpReads;
	assemblingreadtype *this_assemblingRead;
	int32_t i, lockedNum, appearNum, matedNum;
	char orient;

	fpReads = fopen(outfile, "w");
	if(fpReads==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, outfile);
		return FAILED;
	}

	lockedNum = 0;
	appearNum = 0;
	matedNum = 0;
	this_assemblingRead = decisionTable;
	for(i=0; i<readsNum; i++)
	{
		if(this_assemblingRead->matedFlag==YES)
		{
			switch(this_assemblingRead->orientation)
			{
				case ORIENTATION_PLUS: orient = '+'; break;
				case ORIENTATION_MINUS: orient = '-'; break;
				default: printf("line=%d, In %s(), i=%d, invalid read orientation=%d, error!\n", __LINE__, __func__, i, (int32_t)this_assemblingRead->orientation);
			}

			fprintf(fpReads, "rid=%lu, firstBasePos=%d, firstContigPos=%d, orientation=%c, status=%d, "
					"matchBaseNum=%d, unmatchBaseNum=%d, basePos=%d, lastMatchedBasePos=%d, "
					"seqlen=%d, entriesNumReadseq=%d, baseNumLastEentryReadseq=%d, "
					"kmerNumEnd5=%d, kmerNumEnd3=%d, "
					"locked=%d, matedFlag=%d\n",
					(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstBasePos, this_assemblingRead->firstContigPos, orient, this_assemblingRead->status,
					this_assemblingRead->matchBaseNum, this_assemblingRead->unmatchBaseNum, this_assemblingRead->basePos, this_assemblingRead->lastMatchedBasePos,
					this_assemblingRead->seqlen, this_assemblingRead->entriesNumReadseq, this_assemblingRead->baseNumLastEentryReadseq,
					this_assemblingRead->kmerNumEnd5, this_assemblingRead->kmerNumEnd3,
					this_assemblingRead->locked, this_assemblingRead->matedFlag);
			fprintf(fpReads, "\t+: %s\n", getReadBaseByInt(this_assemblingRead->readseq, this_assemblingRead->seqlen));
			fprintf(fpReads, "\t-: %s\n", getReverseReadBaseByInt(this_assemblingRead->readseq, this_assemblingRead->seqlen));

			if(this_assemblingRead->locked==1)
				lockedNum++;

			if(this_assemblingRead->basePos==this_assemblingRead->lastMatchedBasePos)
				appearNum ++;

			if(this_assemblingRead->matedFlag==YES)
				matedNum ++;
		}

		this_assemblingRead++;
	}
	fprintf(fpReads, "readsNum=%d, appearNum=%d, locked=%d, matedNum=%d\n", readsNum, appearNum, lockedNum, matedNum);

	fclose(fpReads);
	fpReads = NULL;

	return SUCCESSFUL;
}

short outputFailedReadsInDecisionTable(assemblingreadtype *decisionTable, int itemNumDecisionTable, int contigID, int contigNodesNum)
{
	assemblingreadtype *this_assemblingRead;
	int32_t i, failedNum;
	char orient;

	failedNum = 0;
	this_assemblingRead = decisionTable;
	for(i=0; i<itemNumDecisionTable; i++)
	{
		if(this_assemblingRead->status==FAILED_STATUS)
		{
			switch(this_assemblingRead->orientation)
			{
				case ORIENTATION_PLUS: orient = '+'; break;
				case ORIENTATION_MINUS: orient = '-'; break;
				default: printf("line=%d, In %s(), i=%d, invalid read orientation=%d, error!\n", __LINE__, __func__, i, (int32_t)this_assemblingRead->orientation);
			}

			printf("\trid=%lu, firstBasePos=%d, firstContigPos=%d, orientation=%c, basePos=%d, lastMatchedBasePos=%d, locked=%d, matedFlag=%d, status=%d\n",
					(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstBasePos, this_assemblingRead->firstContigPos, orient,
					this_assemblingRead->basePos, this_assemblingRead->lastMatchedBasePos, this_assemblingRead->locked, this_assemblingRead->matedFlag, this_assemblingRead->status);

			failedNum ++;
		}

		this_assemblingRead ++;
	}
	if(failedNum>0)
		printf("contigID=%d, contigNodesNum=%d, failedNum=%d\n", contigID, contigNodesNum, failedNum);

	return SUCCESSFUL;
}

/**
 * Output the locked reads in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputLockedReadsInDecisionTable(assemblingreadtype *decisionTable, int readsNum)
{
	assemblingreadtype *this_assemblingRead;
	int32_t i, lockedNum;
	char orient;

	printf("Output locked reads:\n");

	lockedNum = 0;
	this_assemblingRead = decisionTable;
	for(i=0; i<readsNum; i++)
	{
		if(this_assemblingRead->locked==YES)
		{
			switch(this_assemblingRead->orientation)
			{
				case ORIENTATION_PLUS: orient = '+'; break;
				case ORIENTATION_MINUS: orient = '-'; break;
				default: printf("line=%d, In %s(), i=%d, invalid read orientation=%d, error!\n", __LINE__, __func__, i, (int32_t)this_assemblingRead->orientation);
			}

			printf("rid=%lu, firstBasePos=%d, firstContigPos=%d, orientation=%c, status=%d, basePos=%d, lastMatchedBasePos=%d, matedFlag=%d\n",
					(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstBasePos, this_assemblingRead->firstContigPos, orient, this_assemblingRead->status,
					this_assemblingRead->basePos,	this_assemblingRead->lastMatchedBasePos, this_assemblingRead->matedFlag);
			lockedNum++;
		}

		this_assemblingRead++;
	}
	printf("readsNum=%d, LockedNum=%d\n", readsNum, lockedNum);
	return SUCCESSFUL;
}

/**
 * Output the mated reads in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputMatedReadsInDecisionTable(assemblingreadtype *decisionTable, int readsNum)
{
	assemblingreadtype *this_assemblingRead;
	int32_t i, matedNum;
	char orient;

	printf("Output mated reads in decision table:\n");

	matedNum = 0;
	this_assemblingRead = decisionTable;
	for(i=0; i<readsNum; i++)
	{
		if(this_assemblingRead->matedFlag==YES)
		{
			switch(this_assemblingRead->orientation)
			{
				case ORIENTATION_PLUS: orient = '+'; break;
				case ORIENTATION_MINUS: orient = '-'; break;
				default: printf("line=%d, In %s(), i=%d, invalid read orientation=%d, error!\n", __LINE__, __func__, i, (int32_t)this_assemblingRead->orientation);
			}

			printf("rid=%lu, firstBasePos=%d, firstContigPos=%d, orientation=%c, status=%d, basePos=%d, lastMatchedBasePos=%d, matedFlag=%d\n",
					(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstBasePos, this_assemblingRead->firstContigPos, orient, this_assemblingRead->status,
					this_assemblingRead->basePos,	this_assemblingRead->lastMatchedBasePos, this_assemblingRead->matedFlag);
			matedNum++;
		}

		this_assemblingRead ++;
	}
	printf("readsNum=%d, matedNum=%d\n", readsNum, matedNum);
	return SUCCESSFUL;
}

/**
 * Output the kmer in k-mer hash table by specifying the hash code and integer sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputKmer(graphtype *graph, int hashcode, uint64_t *kmerSeqInt)
{
	kmertype *kmer;

	kmer = getKmerByHash(hashcode, kmerSeqInt, graph);;
	if(!kmer)
	{
		printf("In outputKmer(), the kmer==NULL.\n");
		return FAILED;
	}
	printf("kmerseq=%s, multi=%u, arraysize=%u\n", getKmerBaseByInt(kmerSeqInt), kmer->multiplicity, kmer->arraysize);
	ridpostype *rid_pos = kmer->ppos;
	int i = 0, posNum = kmer->arraysize;
	for(; i<posNum; i++)
	{
		printf("\trid=%lu, pos=%u, delsign=%u, reserved=%u\n", (uint64_t)rid_pos->rid, rid_pos->pos, rid_pos->delsign, rid_pos->reserved);
		rid_pos++;
	}
	return SUCCESSFUL;
}

/**
 * Output the kmer in k-mer hash table by specifying the hash code and integer sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputReadsOfKmer(char *readseq, int32_t startPos, graphtype *graph)
{
	uint64_t hashcode, kmerSeqIntTmp[entriesPerKmer], kmerSeqIntRevTmp[entriesPerKmer];
	kmertype *kmer;
	int32_t readseqLen, remainSize;
	char *kmerseq;

	readseqLen = strlen(readseq);
	remainSize = readseqLen - startPos;
	if(readseqLen<kmerSize)
	{
		printf("seq length should not less than %d.\n", kmerSize);
		return FAILED;
	}
	else if(remainSize<kmerSize)
	{
		printf("start pos should not larger than %d.\n", readseqLen-kmerSize);
		return FAILED;
	}

	// generate the kmer integer sequence
	kmerseq = readseq + startPos;
	if(generateKmerSeqInt(kmerSeqIntTmp, kmerseq)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	hashcode = kmerhashInt(kmerSeqIntTmp);
	kmer = getKmerByHash(hashcode, kmerSeqIntTmp, graph);

	if(kmer)
	{
		printf("========= ridpos information for plus orientation kmer:\n");
		if(outputRidposReadseq(kmer->ppos, kmer->arraysize, ORIENTATION_PLUS, graph)==FAILED)
		{
			printf("line=%d, In %s(), cannot output the kmer information, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	kmer = getReverseKmer(kmerSeqIntRevTmp, kmerSeqIntTmp, graph);
	if(kmer)
	{
		printf("========= ridpos information for minus orientation kmer:\n");
		if(outputRidposReadseq(kmer->ppos, kmer->arraysize, ORIENTATION_MINUS, graph)==FAILED)
		{
			printf("line=%d, In %s(), cannot output the kmer information, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}


	return SUCCESSFUL;
}

short outputRidpos(ridpostype *ridpos, int posNum)
{
	int32_t i;

	printf("arraysize=%d\n", posNum);
	for(i=0; i<posNum; i++)
	{
		printf("\trid=%lu, pos=%u, delsign=%u, reserved=%u\n", (uint64_t)ridpos->rid, ridpos->pos, ridpos->delsign, ridpos->reserved);
		ridpos++;
	}

	return SUCCESSFUL;
}

short outputRidposReadseq(ridpostype *ridpos, int32_t posNum, int32_t orient, graphtype *graph)
{
	char readseq[1000];
	int32_t i, readBlockID, rowNumInReadBlock, maxItemNumPerReadBlock, seqLen;

	readBlock_t *readBlockArr;
	read_t *pRead;
	readseqBlock_t *readseqBlockArr;
	uint64_t *pReadseq;

	if(graph->readSet)
	{
		readBlockArr = graph->readSet->readBlockArr;
		maxItemNumPerReadBlock = graph->readSet->maxItemNumPerReadBlock;
		readseqBlockArr = graph->readSet->readseqBlockArr;
	}else
	{
		printf("readSet is NULL, error!\n");
		return FAILED;
	}


	printf("arraysize=%d\n", posNum);
	for(i=0; i<posNum; i++)
	{

		if(ridpos[i].rid==5551836)
		{
			printf("rid=%ld", (int64_t)ridpos[i].rid);
		}


		readBlockID = (ridpos[i].rid - 1) / maxItemNumPerReadBlock;
		rowNumInReadBlock = (ridpos[i].rid - 1) % maxItemNumPerReadBlock;
		pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;
		pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
		seqLen = pRead->seqlen;

		if(orient==ORIENTATION_PLUS)
			getReadBaseFromPosByInt(readseq, pReadseq, seqLen, 0, seqLen);
		else
			getReverseReadBaseFromPosByInt(readseq, pReadseq, seqLen, 0, seqLen);

		printf("\trid=%lu, pos=%u, delsign=%u, seqLen=%u, seq=%s\n", (uint64_t)ridpos[i].rid, ridpos[i].pos, ridpos[i].delsign, seqLen, readseq);
	}

	return SUCCESSFUL;
}


/**
 * Output the successful reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
void outputSuccessReads(successRead_t *successReadArray, int32_t successReadNum)
{
	int32_t i;
	char orient;

	for(i=0; i<successReadNum; i++)
	{
		switch(successReadArray[i].orientation)
		{
			case ORIENTATION_PLUS: orient = '+'; break;
			case ORIENTATION_MINUS: orient = '-'; break;
			default: printf("line=%d, In %s(), i=%d, invalid read orientation=%d, error!\n", __LINE__, __func__, i, (int32_t)successReadArray[i].orientation);
		}

		printf("Success Reads[%d]: rid=%ld, pos=%d, matchnum=%u, orientation=%c, seqlen=%d, hangingIndex=%d\n",
				i+1, (int64_t)successReadArray[i].rid, (int32_t)successReadArray[i].pos, successReadArray[i].matchnum, orient, successReadArray[i].seqlen, successReadArray[i].hangingIndex);
	}
	printf("There are %d success reads\n", i);
}

short outputReadseqByReadID(int64_t readID, graphtype *graph)
{
	char readseq[1000], readseqRev[1000];
	int32_t i, readBlockID, rowNumInReadBlock, maxItemNumPerReadBlock, seqLen;
	int64_t readID_paired;

	readBlock_t *readBlockArr;
	read_t *pRead;
	readseqBlock_t *readseqBlockArr;
	uint64_t *pReadseq;

	if(graph->readSet)
	{
		readBlockArr = graph->readSet->readBlockArr;
		maxItemNumPerReadBlock = graph->readSet->maxItemNumPerReadBlock;
		readseqBlockArr = graph->readSet->readseqBlockArr;
	}else
	{
		printf("readSet is NULL, error!\n");
		return FAILED;
	}

	readBlockID = (readID - 1) / maxItemNumPerReadBlock;
	rowNumInReadBlock = (readID - 1) % maxItemNumPerReadBlock;
	pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;

	if(pRead->validFlag==YES)
	{
		pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
		seqLen = pRead->seqlen;

		getReadBaseFromPosByInt(readseq, pReadseq, seqLen, 0, seqLen);
		getReverseReadBaseFromPosByInt(readseqRev, pReadseq, seqLen, 0, seqLen);

		printf("read: %ld, len=%d\n", readID, seqLen);
		printf("\t%s\n\t%s\n", readseq, readseqRev);
	}else
	{
		printf("read %ld is invalid.\n", readID);
	}


	// paired read
	if(readID%2==1)
		readID_paired = readID + 1;
	else
		readID_paired = readID - 1;

	readBlockID = (readID_paired - 1) / maxItemNumPerReadBlock;
	rowNumInReadBlock = (readID_paired - 1) % maxItemNumPerReadBlock;
	pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;

	if(pRead->validFlag==YES)
	{
		pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
		seqLen = pRead->seqlen;

		getReadBaseFromPosByInt(readseq, pReadseq, seqLen, 0, seqLen);
		getReverseReadBaseFromPosByInt(readseqRev, pReadseq, seqLen, 0, seqLen);

		printf("paired read: %ld, len=%d\n", readID_paired, seqLen);
		printf("\t%s\n\t%s\n", readseq, readseqRev);
	}else
	{
		printf("read %ld is invalid.\n", readID_paired);
	}

	return SUCCESSFUL;
}

/**
 * Check the k-mer hash table including the read set and the kmers.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short checkGraph(graphtype *graph)
{
	printf("Begin checking k-mer hash table, please wait ...\n");

	uint64_t i, j, hashcode, *seqInt;
	kmertype *kmer;

	readBlock_t *pReadBlockArray;
	read_t *pRead, *pReadArray;
	int64_t rid;
	kmerHashBucket_t *pKmerBucket;


	readSet = graph->readSet;

	// check the read set
	rid = 0;
	for(i=0; i<readSet->blocksNumRead; i++)
	{
		pReadBlockArray = readSet->readBlockArr + i;
		pReadArray = pReadBlockArray->readArr;
		for(j=0; j<pReadBlockArray->itemNum; j++)
		{
			rid ++;
			pRead = pReadArray + j;
			if(pRead->successFlag==YES)
			{
				printf("line=%d, In %s(), error read %ld\n", __LINE__, __func__, rid);
				return FAILED;
			}
		}
	}

	// check the k-mer hash table
	for(i=0; i<hashTableSize; i++)
	{
		pKmerBucket = graph->kmerHashtable + i;
		if(pKmerBucket->kmerBlockID>0)
		{
			kmer = graph->kmerBlockArr[pKmerBucket->kmerBlockID-1].kmerArr + pKmerBucket->itemRowKmerBlock;
			while(kmer)
			{
				seqInt = graph->kmerSeqBlockArr[kmer->kmerseqBlockID-1].kmerSeqArr + kmer->itemRowKmerseqBlock * graph->entriesPerKmer;
				hashcode = kmerhashInt(seqInt);
				if(hashcode!=i)
				{
					printf("line=%d, In %s(), hashcode=%lu, i=%lu, error!\n", __LINE__, __func__, hashcode, i);
					return FAILED;
				}

				if(kmer->arraysize!=kmer->multiplicity)
				{
					printf("line=%d, In %s(), arraysize=%u != multiplicity=%u, error!\n", __LINE__, __func__, kmer->arraysize, kmer->multiplicity);
					return FAILED;
				}

				for(j=0; j<kmer->arraysize-1; j++)
				{
					if(kmer->ppos[j].rid>kmer->ppos[j+1].rid)
					{
						printf("line=%d, In %s(), rID Error.\n", __LINE__, __func__);
						return FAILED;
					}
				}

				// get next kmer
				if(kmer->nextKmerBlockID>0)
					kmer = graph->kmerBlockArr[kmer->nextKmerBlockID-1].kmerArr + kmer->nextItemRowKmerBlock;
				else
					kmer = NULL;
			}
		}
	}
	printf("End checking the graph, congratulations.\n");

	return SUCCESSFUL;
}

/**
 * Check the PEHash table.
 *  @return:
 *  	If succeeds, return SUCCESFUL; otherwise, return FAILED.
 */
short checkPEHashtable(PERead_t **PEHashArray)
{
	int i;

	for(i=0; i<TABLE_SIZE_HASH_PE; i++)
	{
		if(PEHashArray[i])
		{
			printf("line=%d, In %s(), the PEHash table is not empty, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}


short outputContigToTmpFile(char *fileOut, contigtype *contigArray, int64_t contigNodesNum, int outFileType)
{
	FILE *fpContigTmp;
	successRead_t *ridposorientation;
	int32_t i, j, num;
	char tmp_base, orient, fileName[256];

	strcpy(fileName, fileOut);
	if(outFileType==HANGING_READ_TYPE_CONTIG_FILE)
		strcat(fileName, ".hang");

	fpContigTmp = fopen(fileName, "w");
	if(fpContigTmp==NULL)
	{
		printf("line=%d, In %s(), can not open file [%s].\n", __LINE__, __func__, fileName);
		return FAILED;
	}

	if(outFileType==BASE_TYPE_FASTA_CONTIG_FILE)
	{
		fprintf(fpContigTmp, ">%d length: %ld\n", 1, contigNodesNum);
		for(i=0; i<contigNodesNum; i++)
		{
			switch(contigArray[i].base)
			{
				case 0: tmp_base = 'A'; break;
				case 1: tmp_base = 'C'; break;
				case 2: tmp_base = 'G'; break;
				case 3: tmp_base = 'T'; break;
				default: printf("line=%d, In %s(), error base int: %d\n", __LINE__, __func__, contigArray[i].base); return FAILED;
			}

			fprintf(fpContigTmp, "%c", tmp_base);
		}
		fprintf(fpContigTmp, "\n");
	}else
	{
		fprintf(fpContigTmp, ">%d length: %ld\n", 1, contigNodesNum);

		for(i=0; i<contigNodesNum; i++)
		{
			switch(contigArray[i].base)
			{
				case 0: tmp_base = 'A'; break;
				case 1: tmp_base = 'C'; break;
				case 2: tmp_base = 'G'; break;
				case 3: tmp_base = 'T'; break;
				default: printf("line=%d, In %s(), error base int: %d\n", __LINE__, __func__, contigArray[i].base); return FAILED;
			}

			fprintf(fpContigTmp, "%d\t%c\t%d", contigArray[i].index, tmp_base, contigArray[i].ridposnum);
			ridposorientation = contigArray[i].pridposorientation;
			num = contigArray[i].ridposnum;
			for(j=0; j<num; j++)
			{
				switch(ridposorientation[j].orientation)
				{
					case ORIENTATION_PLUS: orient = '+'; break;
					case ORIENTATION_MINUS: orient = '-'; break;
					default: printf("line=%d, In %s(), j=%d, invalid read orientation=%d, error!\n", __LINE__, __func__, j, (int32_t)ridposorientation[j].orientation);
				}

				fprintf(fpContigTmp, "\t(%lu,%d,%d,%c)", (uint64_t)ridposorientation[j].rid, (int32_t)ridposorientation[j].pos, (int32_t)ridposorientation[j].matchnum, orient);
			}
			fprintf(fpContigTmp, "\n");
		}
	}

	fclose(fpContigTmp);
	fpContigTmp = NULL;

	return SUCCESSFUL;
}


short outputPEHashArray(PERead_t **PEHashArray)
{
	int i, totalReadNum;
	PERead_t *tmpRead;

	totalReadNum = 0;
	for(i=0; i<TABLE_SIZE_HASH_PE; i++)
	{
		if(PEHashArray[i])
		{
			tmpRead = PEHashArray[i];
			while(tmpRead)
			{
				printf("[%d]: rid=%lu, cpos=%u, orient=%u, seqlen=%d\n", i, tmpRead->rid, tmpRead->cpos, tmpRead->orient, tmpRead->seqlen);
				totalReadNum ++;
				tmpRead = tmpRead->next;
			}
		}
	}

	printf("totalReadNum=%d\n", totalReadNum);

	return SUCCESSFUL;
}


short checkReadListArr(readList_t *readListArray, int64_t itemNumInReadListArray)
{
	int64_t i;

	for(i=0; i<itemNumInReadListArray-1; i++)
	{
		if(readListArray[i].readID >= readListArray[i+1].readID)
		{
			printf("line=%d, In %s(), the data was not correctly ordered, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}


short outputReadListToFile(char *readListFile, readList_t *readListArray, readPos_t *readPosArray, int64_t itemNumInReadListArray)
{
	uint64_t i, j, matchNum, readID;
	FILE *fpReadList;
	readPos_t *pReadPos;
	char orient;

	fpReadList = fopen(readListFile, "w");
	if(fpReadList==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readListFile);
		return FAILED;
	}

	for(i=0; i<itemNumInReadListArray; i++)
	{
		readID = readListArray[i].readID;
		matchNum = readListArray[i].matchNum;
		pReadPos = readPosArray + readListArray[i].firstRow;
		for(j=0; j<matchNum; j++)
		{
			switch(pReadPos[j].orientation)
			{
				case ORIENTATION_PLUS: orient = '+'; break;
				case ORIENTATION_MINUS: orient = '-'; break;
				default: printf("line=%d, In %s(), j=%lu, invalid read orientation=%u, error!\n", __LINE__, __func__, j, pReadPos[j].orientation);
			}

			fprintf(fpReadList, "[%lu]:\t%u\t%u\t%c\n", readID, pReadPos[j].contigPos, pReadPos[j].matchBaseNum, orient);
		}
	}

	fclose(fpReadList);
	fpReadList = NULL;

	return SUCCESSFUL;
}


short outputMatedReadsInReadListToFile(char *readListFile, readList_t *readListArray, readPos_t *readPosArray, int64_t itemNumInReadListArray)
{
	uint64_t i, j, matchNum, readID, seqLen;
	FILE *fpReadList;
	readPos_t *pReadPos;
	char orient;

	fpReadList = fopen(readListFile, "w");
	if(fpReadList==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readListFile);
		return FAILED;
	}

	i = 0;
	while(i<itemNumInReadListArray-1)
	{
		if(readListArray[i].readID % 2 == 1 && readListArray[i].readID+1==readListArray[i+1].readID)
		{
			if(readPosArray[readListArray[i].firstRow].matchBaseNum==readListArray[i].seqlen && readPosArray[readListArray[i+1].firstRow].matchBaseNum==readListArray[i+1].seqlen)
			{
				readID = readListArray[i].readID;
				seqLen = readListArray[i].seqlen;
				matchNum = readListArray[i].matchNum;
				pReadPos = readPosArray + readListArray[i].firstRow;
				for(j=0; j<matchNum; j++)
				{
					switch(pReadPos[j].orientation)
					{
						case ORIENTATION_PLUS: orient = '+'; break;
						case ORIENTATION_MINUS: orient = '-'; break;
						default: printf("line=%d, In %s(), j=%lu, invalid read orientation=%u, error!\n", __LINE__, __func__, j, pReadPos[j].orientation);
					}

					fprintf(fpReadList, "[%lu]:\t%lu\t%u\t%u\t%c\n", readID, seqLen, pReadPos[j].contigPos, pReadPos[j].matchBaseNum, orient);
				}

				readID = readListArray[i+1].readID;
				seqLen = readListArray[i+1].seqlen;
				matchNum = readListArray[i+1].matchNum;
				pReadPos = readPosArray + readListArray[i+1].firstRow;
				for(j=0; j<matchNum; j++)
				{
					switch(pReadPos[j].orientation)
					{
						case ORIENTATION_PLUS: orient = '+'; break;
						case ORIENTATION_MINUS: orient = '-'; break;
						default: printf("line=%d, In %s(), j=%lu, invalid read orientation=%u, error!\n", __LINE__, __func__, j, pReadPos[j].orientation);
					}

					fprintf(fpReadList, "[%lu]:\t%lu\t%u\t%u\t%c\n", readID, seqLen, pReadPos[j].contigPos, pReadPos[j].matchBaseNum, orient);
				}
			}

			i += 2;
		}else
		{
			i ++;
		}
	}

	fclose(fpReadList);
	fpReadList = NULL;

	return SUCCESSFUL;
}


short convertFragmentSizeFileInText(const char *fragmentSizeFile)
{
	int32_t fragmentSize;
	char fileName[256];
	FILE *fpFragSize_In, *fpFragSize_Out;

	strcpy(fileName, fragmentSizeFile);
	strcat(fileName, ".txt");

	fpFragSize_In = fopen(fragmentSizeFile, "rb");
	if(fpFragSize_In==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fragmentSizeFile);
		return FAILED;
	}

	fpFragSize_Out = fopen(fileName, "w");
	if(fpFragSize_Out==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fileName);
		return FAILED;
	}

	while(1)
	{
		if(fread(&fragmentSize, sizeof(int32_t), 1, fpFragSize_In)!=1)
		{
			if(feof(fpFragSize_In))
			{
				break;
			}else
			{
				printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		fprintf(fpFragSize_Out, "%d\n", fragmentSize);
	}

	fclose(fpFragSize_Out);
	fpFragSize_Out = NULL;
	fclose(fpFragSize_In);
	fpFragSize_Out = NULL;

	return SUCCESSFUL;
}

short outputReadPosInGraph(int64_t readID, graphtype *graph)
{
	uint64_t i, j, *seqInt;
	kmertype *kmer;
	ridpostype *ridpostable;
	int posNum;
	kmerHashBucket_t *pKmerBucket;


	for(i=0; i<hashTableSize; i++)
	{
		pKmerBucket = graph->kmerHashtable + i;
		if(pKmerBucket->kmerBlockID>0)
		{
			kmer = graph->kmerBlockArr[pKmerBucket->kmerBlockID-1].kmerArr + pKmerBucket->itemRowKmerBlock;
			while(kmer)
			{
				ridpostable = kmer->ppos;
				posNum = kmer->arraysize;
				for(j=0; j<posNum; j++)
				{
					if(ridpostable[j].rid==readID)
					{
						seqInt = graph->kmerSeqBlockArr[kmer->kmerseqBlockID-1].kmerSeqArr + kmer->itemRowKmerseqBlock * graph->entriesPerKmer;
						printf("[%lu]: rid=%lu, pos=%lu, kmerSeq=%s\n", j, (uint64_t)ridpostable[j].rid, (uint64_t)ridpostable[j].pos, getKmerBaseByInt(seqInt));
					}
				}

				if(kmer->nextKmerBlockID>0)
					kmer = graph->kmerBlockArr[kmer->nextKmerBlockID-1].kmerArr + kmer->nextItemRowKmerBlock;
				else
					kmer = NULL;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Output items in dtRowHashtable.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputDtRowHashtable(dtRowIndex_t **pDtRowHashtable)
{
	int32_t i;
	int64_t itemNum;
	dtRowIndex_t *pDtRow;

	itemNum = 0;
	for(i=0; i<TABLE_SIZE_HASH_DTROWINDEX; i++)
	{
		pDtRow = pDtRowHashtable[i];
		while(pDtRow)
		{
			printf("hash %d: rid=%ld, dtRow=%d\n", i, (int64_t)pDtRow->rid, (int32_t)pDtRow->dtRow);
			itemNum ++;

			pDtRow = pDtRow->next;
		}
	}

	printf("itemNum=%ld\n", itemNum);

	return SUCCESSFUL;
}

/**
 * Visit kmers in k-mer hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short visitKmers(graphtype *graph)
{
	int i, kmerNum;
	kmertype *kmer;
	kmerHashBucket_t *pKmerBucket;


	kmerNum = 0;
	for(i=0; i<hashTableSize; i++)
	{
		pKmerBucket = graph->kmerHashtable + i;
		if(pKmerBucket->kmerBlockID>0)
		{
			kmer = graph->kmerBlockArr[pKmerBucket->kmerBlockID-1].kmerArr + pKmerBucket->itemRowKmerBlock;
			while(kmer)
			{
				kmerNum ++;

				if(kmer->nextKmerBlockID>0)
					kmer = graph->kmerBlockArr[kmer->nextKmerBlockID-1].kmerArr + kmer->nextItemRowKmerBlock;
				else
					kmer = NULL;
			}
		}
	}

	printf("kmerNum=%d\n", kmerNum);

	return SUCCESSFUL;
}

/**
 * Output read sequences in read set.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputReadseqInReadset(char *outfile, readSet_t *readSet)
{
	int32_t i, j;
	readBlock_t *pReadBlockArray;
	read_t *pRead, *pReadArray;
	uint64_t *readseqInt;
	int64_t rid;
	FILE *fpOut;

	fpOut = fopen(outfile, "w");
	if(fpOut==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, outfile);
		return FAILED;
	}

	rid = 0;
	for(i=0; i<readSet->blocksNumRead; i++)
	{
		pReadBlockArray = readSet->readBlockArr + i;
		pReadArray = pReadBlockArray->readArr;
		for(j=0; j<pReadBlockArray->itemNum; j++)
		{
			rid ++;
			pRead = pReadArray + j;
			readseqInt = readSet->readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;

			fprintf(fpOut, "rid=%ld, seqlen=%d, nBaseNum=%d, validFlag=%d, successFlag=%d, seq=%s\n", rid, (int32_t)pRead->seqlen, (int32_t)pRead->nBaseNum, (int32_t)pRead->validFlag, (int32_t)pRead->successFlag, getReadBaseByInt(readseqInt, pRead->seqlen));
		}
	}

	fclose(fpOut);
	fpOut = NULL;

	return SUCCESSFUL;
}


/**
 * Output read sequences in read set.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputLinkedReads(int32_t contigID1, int32_t contigEnd1, int32_t contigID2, int32_t contigEnd2, contigGraph_t *contigGraph, readSet_t *readSet)
{
	uint64_t *readSeqInt, *readSeqIntRev;
	char readseqTmp[MAX_READ_LEN_IN_BUF+1];

	readBlock_t *readBlockArr;
	read_t *pRead;
	readseqBlock_t *readseqBlockArr;
	int32_t readBlockID, rowNumInReadBlock, maxItemNumPerReadBlock, seqLen; // block id starts from 0

	readMatchInfoBlock_t *readMatchInfoBlockArr;
	readMatchInfo_t *pReadMatchInfo;
	int32_t readMatchInfoBlockID, rowNumInReadMatchInfoBlock, maxItemNumPerReadMatchInfoBlock;
	contigRead_t *contigReadArray;
	int32_t i, j, contigReadNum, readOrient1, readOrient2, contigPos1, contigPos2;
	int64_t readID, readID_paired;

	readBlockArr = readSet->readBlockArr;
	maxItemNumPerReadBlock = readSet->maxItemNumPerReadBlock;
	readseqBlockArr = readSet->readseqBlockArr;

	readMatchInfoBlockArr = readSet->readMatchInfoBlockArr;
	maxItemNumPerReadMatchInfoBlock = readSet->maxItemNumPerReadMatchInfoBlock;

	if(contigEnd1==1)
	{
		if(contigGraph->contigItemArray[contigID1-1].alignRegSizeEnd3>0)
		{
			contigReadArray = contigGraph->contigItemArray[contigID1-1].contigReadArrayEnd3;
			contigReadNum = contigGraph->contigItemArray[contigID1-1].contigReadNumEnd3;
		}else
		{
			printf("short contigLen1: %d\n", contigGraph->contigItemArray[contigID1-1].contigLen);
			return FAILED;
		}
	}else
	{
		contigReadArray = contigGraph->contigItemArray[contigID1-1].contigReadArrayEnd5;
		contigReadNum = contigGraph->contigItemArray[contigID1-1].contigReadNumEnd5;
	}

	if(contigEnd2==1)
	{
		if(contigGraph->contigItemArray[contigID2-1].alignRegSizeEnd3<=0)
		{
			printf("short contigLen2: %d\n", contigGraph->contigItemArray[contigID1-1].contigLen);
			contigEnd2 = 2;
		}
	}

	for(i=0; i<contigReadNum; i++)
	{
		readID = contigReadArray[i].readID;
		readOrient1 = contigReadArray[i].orientation;
		contigPos1 = contigReadArray[i].contigPos;

		if(readID%2==1)
			readID_paired = readID + 1;
		else
			readID_paired = readID - 1;

		readMatchInfoBlockID = (readID_paired - 1) / maxItemNumPerReadMatchInfoBlock;
		rowNumInReadMatchInfoBlock = (readID_paired - 1) % maxItemNumPerReadMatchInfoBlock;
		pReadMatchInfo = readMatchInfoBlockArr[readMatchInfoBlockID].readMatchInfoArr + rowNumInReadMatchInfoBlock;

		readOrient2 = pReadMatchInfo->readOrientation;
		contigPos2 = pReadMatchInfo->contigPos;

		if(pReadMatchInfo->contigID==contigID2 && pReadMatchInfo->contigEnd==contigEnd2)
		{
			// output the read1
			readBlockID = (readID - 1) / maxItemNumPerReadBlock;
			rowNumInReadBlock = (readID - 1) % maxItemNumPerReadBlock;
			pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;
			readSeqInt = readSet->readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
			seqLen = pRead->seqlen;

			if(readOrient1==ORIENTATION_PLUS)
			{
				getReadBaseFromPosByInt(readseqTmp, readSeqInt, seqLen, 0, seqLen);
				printf("%ld\t%d\t%d\t%d\t+\t%d\t%s\n", readID, contigID1, contigEnd1, contigPos1, seqLen, readseqTmp);
			}else
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
				printf("%ld\t%d\t%d\t%d\t-\t%d\t%s\n", readID, contigID1, contigEnd1, contigPos1, seqLen, readseqTmp);

				free(readSeqIntRev);
			}

			// output the read2
			readBlockID = (readID_paired - 1) / maxItemNumPerReadBlock;
			rowNumInReadBlock = (readID_paired - 1) % maxItemNumPerReadBlock;
			pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;
			readSeqInt = readSet->readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
			seqLen = pRead->seqlen;

			if(readOrient2==ORIENTATION_PLUS)
			{
				getReadBaseFromPosByInt(readseqTmp, readSeqInt, seqLen, 0, seqLen);
				printf("%ld\t%d\t%d\t%d\t+\t%d\t%s\n", readID_paired, contigID2, contigEnd2, contigPos2, seqLen, readseqTmp);
			}else
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
				printf("%ld\t%d\t%d\t%d\t-\t%d\t%s\n", readID_paired, contigID2, contigEnd2, contigPos2, seqLen, readseqTmp);

				free(readSeqIntRev);
			}
		}
	}

	return SUCCESSFUL;
}
