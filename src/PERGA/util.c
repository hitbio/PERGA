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
	int i, j, num;
	char tmp_base;

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
			printf("\t(%lu,%d,%d,%c)", (int64_t)ridposorientation[j].rid, (int32_t)ridposorientation[j].pos, (int32_t)ridposorientation[j].matchnum, (int32_t)ridposorientation[j].orientation);
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
	int i, j, num;
	char tmp_base;

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
			printf("\t(%lu,%d,%d,%c)", (uint64_t)ridposorientation[i].rid, (int32_t)ridposorientation[i].pos, ridposorientation[i].matchnum, ridposorientation[i].orientation);
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
	int i, j, num, nodeNum;
	char tmp_base;

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
			printf("\t(%lu,%d,%d,%c)", (uint64_t)ridposorientation[i].rid, (int32_t)ridposorientation[i].pos, (int32_t)ridposorientation[i].matchnum, ridposorientation[i].orientation);
		printf("\n");

		nodeNum ++;
	}
	printf("There are %d nodes.\n", nodeNum);

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
	assemblingreadtype *this_assemblingRead = decisionTable;
	int i, lockedNum, appearNum, matedNum;

	lockedNum = 0;
	appearNum = 0;
	matedNum = 0;
	for(i=0; i<readsNum; i++)
	{
//		printf("rid=%lu, firstpos=%d, orientation=%c, status=%d, kmerappeartimes=%d, kmerunappearblocks=%d, kmerunappeartimes=%d, lastappearpos=%d, lastpos=%d, locked=%d, matedFlag=%d\n",
//				(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstpos, this_assemblingRead->orientation, this_assemblingRead->status,
//				this_assemblingRead->kmerappeartimes, this_assemblingRead->kmerunappearblocks, this_assemblingRead->kmerunappeartimes,
//				this_assemblingRead->lastappearpos,	this_assemblingRead->lastpos, this_assemblingRead->locked, this_assemblingRead->matedFlag);

		printf(	"rid=%lu, firstBasePos=%d, firstContigPos=%d, orientation=%c, status=%d, "
				"matchBaseNum=%d, unmatchBaseNum=%d, basePos=%d, lastMatchedBasePos=%d, "
				"seqlen=%d, entriesNumReadseq=%d, baseNumLastEentryReadseq=%d, "
				"kmerNumEnd5=%d, kmerNumEnd3=%d, "
				"locked=%d, matedFlag=%d\n",
				(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstBasePos, this_assemblingRead->firstContigPos, this_assemblingRead->orientation, this_assemblingRead->status,
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
	assemblingreadtype *this_assemblingRead = decisionTable;
	int i, lockedNum, appearNum, failedNum, matedNum;
	//char fileName[256];

	//sprintf(fileName, "readsInDT_%d_%d.txt", contigID, contigNodesNum);

	fpReads = fopen(outfile, "w");
	if(fpReads==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, outfile);
		return FAILED;
	}

	lockedNum = 0;
	appearNum = 0;
	matedNum = 0;
	for(i=0; i<readsNum; i++)
	{
//		fprintf(fpReads, "rid=%lu, firstpos=%d, orientation=%c, status=%d, kmerappeartimes=%d, kmerunappearblocks=%d, kmerunappeartimes=%d, lastappearpos=%d, lastpos=%d, locked=%d, matedFlag=%d, status=%d\n",
//				(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstpos, this_assemblingRead->orientation, this_assemblingRead->status,
//				this_assemblingRead->kmerappeartimes, this_assemblingRead->kmerunappearblocks, this_assemblingRead->kmerunappeartimes,
//				this_assemblingRead->lastappearpos,	this_assemblingRead->lastpos, this_assemblingRead->locked, this_assemblingRead->matedFlag, this_assemblingRead->status);

		fprintf(fpReads, "rid=%lu, firstBasePos=%d, firstContigPos=%d, orientation=%c, status=%d, "
				"matchBaseNum=%d, unmatchBaseNum=%d, basePos=%d, lastMatchedBasePos=%d, "
				"seqlen=%d, entriesNumReadseq=%d, baseNumLastEentryReadseq=%d, "
				"kmerNumEnd5=%d, kmerNumEnd3=%d, "
				"locked=%d, matedFlag=%d\n",
				(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstBasePos, this_assemblingRead->firstContigPos, this_assemblingRead->orientation, this_assemblingRead->status,
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

/*
	failedNum = 0;
	fprintf(fpReads, "\n\nfailed reads:\n");
	for(i=0; i<readsNum; i++)
	{
		if(this_assemblingRead->status==FAILED_STATUS)
		{
//			fprintf(fpReads, "rid=%lu, firstpos=%d, orientation=%c, status=%d, kmerappeartimes=%d, kmerunappearblocks=%d, kmerunappeartimes=%d, lastappearpos=%d, lastpos=%d, locked=%d, matedFlag=%d, status=%d\n",
//					(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstpos, this_assemblingRead->orientation, this_assemblingRead->status,
//					this_assemblingRead->kmerappeartimes, this_assemblingRead->kmerunappearblocks, this_assemblingRead->kmerunappeartimes,
//					this_assemblingRead->lastappearpos,	this_assemblingRead->lastpos, this_assemblingRead->locked, this_assemblingRead->matedFlag, this_assemblingRead->status);

			fprintf(fpReads, "rid=%lu, firstpos=%d, orientation=%c, status=%d, "
					"basePos=%d, lastMatchedBasePos=%d, locked=%d, matedFlag=%d, status=%d\n",
					(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstpos, this_assemblingRead->orientation, this_assemblingRead->status,
					this_assemblingRead->basePos, this_assemblingRead->lastMatchedBasePos, this_assemblingRead->locked, this_assemblingRead->matedFlag, this_assemblingRead->status);


			failedNum ++;
		}

		this_assemblingRead ++;
	}
	fprintf(fpReads, "failedNum=%d\n", failedNum);
*/

	fclose(fpReads);
	fpReads = NULL;

	return SUCCESSFUL;
}

short outputFailedReadsInDecisionTable(assemblingreadtype *decisionTable, int itemNumDecisionTable, int contigID, int contigNodesNum)
{
	assemblingreadtype *this_assemblingRead = decisionTable;
	int i, failedNum;

	failedNum = 0;
	for(i=0; i<itemNumDecisionTable; i++)
	{
		if(this_assemblingRead->status==FAILED_STATUS)
		{
			printf("\trid=%lu, firstBasePos=%d, firstContigPos=%d, orientation=%c, basePos=%d, lastMatchedBasePos=%d, locked=%d, matedFlag=%d, status=%d\n",
					(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstBasePos, this_assemblingRead->firstContigPos, this_assemblingRead->orientation,
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
	if(lockedReadsNum==0)
	{
		return SUCCESSFUL;
	}
	assemblingreadtype *this_assemblingRead = decisionTable;
	printf("Output locked reads:\n");
	int i, lockedNum = 0;
	for(i=0; i<readsNum; i++)
	{
		if(this_assemblingRead->locked==1)
		{
			printf("rid=%lu, firstBasePos=%d, firstContigPos=%d, orientation=%c, status=%d, basePos=%d, lastMatchedBasePos=%d, matedFlag=%d\n",
					(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstBasePos, this_assemblingRead->firstContigPos, this_assemblingRead->orientation, this_assemblingRead->status,
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
	assemblingreadtype *this_assemblingRead = decisionTable;
	printf("Output mated reads in decision table:\n");
	int i, matedNum = 0;
	for(i=0; i<readsNum; i++)
	{
		if(this_assemblingRead->matedFlag==YES)
		{
			printf("rid=%lu, firstBasePos=%d, firstContigPos=%d, orientation=%c, status=%d, basePos=%d, lastMatchedBasePos=%d, matedFlag=%d\n",
					(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstBasePos, this_assemblingRead->firstContigPos, this_assemblingRead->orientation, this_assemblingRead->status,
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


short outputRidpos(ridpostype *ridpos, int posNum)
{
	printf("arraysize=%d\n", posNum);
	int i = 0;
	for(; i<posNum; i++)
	{
		printf("\trid=%lu, pos=%u, delsign=%u, reserved=%u\n", (uint64_t)ridpos->rid, ridpos->pos, ridpos->delsign, ridpos->reserved);
		ridpos++;
	}
	return SUCCESSFUL;
}

/**
 * Output the successful reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
void outputSuccessReads(successRead_t *successReadArray, int successReadNum)
{
	int i = 0;
	for(; i<successReadNum; i++)
	{
		printf("Success Reads[%d]: rid=%lu, pos=%u, matchnum=%u, orientation=%c, seqlen=%d\n",
				i+1, (uint64_t)successReadArray[i].rid, successReadArray[i].pos, successReadArray[i].matchnum, successReadArray[i].orientation, successReadArray[i].seqlen);
	}
	printf("There are %d success reads\n", i);
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
	int i, j, num;
	char tmp_base, fileName[256];

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
				fprintf(fpContigTmp, "\t(%lu,%d,%d,%c)", (uint64_t)ridposorientation[j].rid, (int32_t)ridposorientation[j].pos, (int32_t)ridposorientation[j].matchnum, ridposorientation[j].orientation);
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
				printf("[%d]: rid=%lu, cpos=%u, orient=%u\n", i, tmpRead->rid, tmpRead->cpos, tmpRead->orient);
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
			fprintf(fpReadList, "[%lu]:\t%u\t%u\t%c\n", readID, pReadPos[j].contigPos, pReadPos[j].matchBaseNum, pReadPos[j].orientation);
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
					fprintf(fpReadList, "[%lu]:\t%lu\t%u\t%u\t%c\n", readID, seqLen, pReadPos[j].contigPos, pReadPos[j].matchBaseNum, pReadPos[j].orientation);
				}

				readID = readListArray[i+1].readID;
				seqLen = readListArray[i+1].seqlen;
				matchNum = readListArray[i+1].matchNum;
				pReadPos = readPosArray + readListArray[i+1].firstRow;
				for(j=0; j<matchNum; j++)
				{
					fprintf(fpReadList, "[%lu]:\t%lu\t%u\t%u\t%c\n", readID, seqLen, pReadPos[j].contigPos, pReadPos[j].matchBaseNum, pReadPos[j].orientation);
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
 * Output items in navigation occurrence queue.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputNaviOccQueue(double *naviOccQueuePara, int itemNumNaviOccQueuePara, int frontRowNaviOccQueuePara)
{
	int i, j;

	j = frontRowNaviOccQueuePara;
	for(i=0; i<itemNumNaviOccQueuePara; i++)
	{
		printf("naviOccQueue[%d]: %.2f\n", i, naviOccQueuePara[j]);
		j = (j+1) % maxItemNumNaviOccQueue;
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
 * test SVM.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short testSVM()
{
	struct sample_data_node
	{
		double data[5];
	};

	struct sample_group_node
	{
		int32_t labelID[3];
	};

	int32_t outClass, i, j, colsNumTmp, rowsNumTestData, dataCols[5] = {1, 1, 0, 0, 0};
	svmModel_t *svmModel;
	svmSampleVector_t svmSample;
	double tmpData[5];
	struct sample_data_node *testData;
	struct sample_group_node *testDataGroup;
	char ch;
	FILE *fpTestData;

	fpTestData = fopen("../model/tmp/svmTestResult.txt", "r");
	if(fpTestData==NULL)
	{
		printf("fopen file error!\n");
		return FAILED;
	}

	// get the number of rows of the test data
	rowsNumTestData = 0;
	ch = fgetc(fpTestData);
	while(ch!=EOF)
	{
		if(ch=='\n') rowsNumTestData ++;
		ch = fgetc(fpTestData);
	}

	testData = (struct sample_data_node*) calloc (rowsNumTestData, sizeof(struct sample_data_node));
	if(testData==NULL)
	{
		printf("error!\n");
		return FAILED;
	}
	testDataGroup = (struct sample_group_node*) calloc (rowsNumTestData, sizeof(struct sample_group_node));
	if(testDataGroup==NULL)
	{
		printf("error!\n");
		return FAILED;
	}

	// load the test data
	rewind(fpTestData);
	for(i=0; i<rowsNumTestData; i++)
		fscanf(fpTestData, "%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\n", testData[i].data, testData[i].data+1, testData[i].data+2, testData[i].data+3, testData[i].data+4, testDataGroup[i].labelID, testDataGroup[i].labelID+1);
	fclose(fpTestData);

	// outClass=0
//	tmpData[0] = 2;
//	tmpData[1] = 1;
//	tmpData[2] = 1.0391;

	// outClass=1
//	tmpData[0] = 4;
//	tmpData[1] = 1;
//	tmpData[2] = 0.6732;

	// outClass=1
	tmpData[0] = 20;
	tmpData[1] = 6;
//	tmpData[2] = 1.4627;

	printf("Begin to load svm ...\n");

	if(loadSvmModel(&svmModel, "../model/tmp/svmKF.txt", "../model/tmp/svmSV.txt", "../model/tmp/svmAlpha.txt", "../model/tmp/svmBias.txt", "../model/tmp/svmScaleData.txt")==FAILED)
	{
		printf("error!\n");
		return FAILED;
	}

	printf("Begin to classify ...\n");

	svmSample.vectorData = tmpData;
	svmSample.colsNum = svmModel->colsNumSupportVector;

	outClass = -1;
	mySvmClassifySE(&outClass, &svmSample, svmModel);
	printf("outClass=%d\n", outClass);

/*
	//outClass = -1;
	for(i=0; i<rowsNumTestData; i++)
	{
		colsNumTmp = 0;
		for(j=0; j<5; j++)
		{
			if(dataCols[j]==1)
				svmSample.vectorData[colsNumTmp++] = testData[i].data[j];
		}

		mySvmClassifySE(testDataGroup[i].labelID+2, &svmSample, svmModel);
	}

	//printf("outClass=%d\n", outClass);

	// check the classification
	for(i=0; i<rowsNumTestData; i++)
	{
		if(testDataGroup[i].labelID[1]!=testDataGroup[i].labelID[2])
		{
			printf("row=%d, classification error!\n", i);
		}
	}
*/
	freeSvmModel(&svmModel);

	free(testData);
	free(testDataGroup);

	return SUCCESSFUL;
}

/**
 * Output the classification result by SVM to file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputClassificationResult()
{
	int32_t i, rowsNumTestData;
	svmModel_t *svmModel;
	svmSampleVector_t svmSample;
	FILE *fpTestData, *fpClassificationResult;
	char ch;

	fpTestData = fopen("../model/testData.txt", "r");
	if(fpTestData==NULL)
	{
		printf("fopen file error!\n");
		return FAILED;
	}

	rowsNumTestData = 0;
	ch = fgetc(fpTestData);
	while(!feof(fpTestData))
	{
		if(ch=='\n')
			rowsNumTestData ++;
		ch = fgetc(fpTestData);
	}

	printf("rowsNumTestData=%d\n", rowsNumTestData);

	double tmpData[3], testData[rowsNumTestData][3];
	int32_t testDataGroup[rowsNumTestData][2];

	rewind(fpTestData);
	for(i=0; i<rowsNumTestData; i++)
		fscanf(fpTestData, "%lf\t%lf\t%lf\t%d\n", testData[i], testData[i]+1, testData[i]+2, testDataGroup[i]);
	fclose(fpTestData);

	printf("Begin load svm ...\n");

	if(loadSvmModel(&svmModel, "../model/tmp/svmKF.txt", "../model/tmp/svmSV.txt", "../model/tmp/svmAlpha.txt", "../model/tmp/svmBias.txt", "../model/tmp/svmScaleData.txt")==FAILED)
	{
		printf("error!\n");
		return FAILED;
	}

	printf("Begin to calssify ...\n");

	svmSample.vectorData = tmpData;
	svmSample.colsNum = svmModel->colsNumSupportVector;

	for(i=0; i<rowsNumTestData; i++)
	{
		svmSample.vectorData[0] = testData[i][0];
		svmSample.vectorData[1] = testData[i][1];
		svmSample.vectorData[2] = testData[i][2];

		mySvmClassifySE(testDataGroup[i]+1, &svmSample, svmModel);
	}

	// output the classification result to file
	fpClassificationResult = fopen("../model/classificationResult.txt", "w");
	if(fpClassificationResult==NULL)
	{
		printf("fopen error!\n");
		return FAILED;
	}

	for(i=0; i<rowsNumTestData; i++)
	{
		fprintf(fpClassificationResult, "%d\t%d\t%.4lf\t%d\t%d\n", (int32_t)testData[i][0], (int32_t)testData[i][1], testData[i][2], testDataGroup[i][0], testDataGroup[i][1]);
	}

	int32_t errClasNum;
	errClasNum = 0;
	for(i=0; i<rowsNumTestData; i++)
		if(testDataGroup[i][0]!=testDataGroup[i][1])
			errClasNum ++;

	printf("errRate=%.4f\n", (float)errClasNum/rowsNumTestData);

	fclose(fpClassificationResult);

	freeSvmModel(&svmModel);

	return SUCCESSFUL;
}

/**
 * Result statistics.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short resultStatistics()
{
	char ch, classFiles[6][256];
	int32_t i, itemNum;
	FILE *fpStatistic, *fpClass;

	strcpy(classFiles[0], "../model/classificationResult_1.txt");
	strcpy(classFiles[1], "../model/classificationResult_2.txt");
	strcpy(classFiles[2], "../model/classificationResult_3.txt");
	strcpy(classFiles[3], "../model/classificationResult_4.txt");
	strcpy(classFiles[4], "../model/classificationResult_5.txt");
	strcpy(classFiles[5], "../model/classificationResult_6.txt");

	fpClass = fopen(classFiles[0], "r");
	if(fpClass==NULL)
	{
		printf("error!\n");
		return FAILED;
	}

	itemNum = 0;
	ch = fgetc(fpClass);
	while(!feof(fpClass))
	{
		if(ch=='\n')
			itemNum ++;
		ch = fgetc(fpClass);
	}
	fclose(fpClass);


	int32_t errClassNumArray[itemNum], errNumArray[7] = {0};

	for(i=0; i<itemNum; i++)
		errClassNumArray[i] = 0;

	// calculate the errNums of each file
	for(i=0; i<6; i++)
	{
		// calculate the errorNum of single file
		if(calcErrNumSingleFile(errClassNumArray, classFiles[i])==FAILED)
		{
			printf("error!\n");
			return FAILED;
		}
	}

	// make statistics
	for(i=0; i<itemNum; i++)
		errNumArray[errClassNumArray[i]] ++;

	int32_t totalErrNum;
	totalErrNum = 0;
	for(i=0; i<itemNum; i++)
	{
		if(errClassNumArray[i]>0)
			totalErrNum ++;
	}
	printf("itemNum=%d, totalErrNum=%d, correctNum=%d\n", itemNum, totalErrNum, itemNum-totalErrNum);

	for(i=0; i<7; i++)
		printf("errNumArray[%d]=%d, ratio=%.4f\n", i, errNumArray[i], (float)errNumArray[i]/totalErrNum);


	// output the number of wrong prediction of each point
	fpStatistic = fopen("../model/errNums.txt", "w");
	if(fpStatistic==NULL)
	{
		printf("error!\n");
		return FAILED;
	}
	for(i=0; i<itemNum; i++)
		fprintf(fpStatistic, "%d\n", errClassNumArray[i]);

	fclose(fpStatistic);

	return SUCCESSFUL;
}

/**
 * Calculate the errorNum of single file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short calcErrNumSingleFile(int32_t *errClassNumArray, char *svmClassFile)
{
	int32_t i, group[2];
	double tmp[3];
	FILE *fpSvmClass;

	fpSvmClass = fopen(svmClassFile, "r");
	if(fpSvmClass==NULL)
	{
		printf("error!\n");
		return FAILED;
	}

	i = 0;
	while(1)
	{
//		if(feof(fpSvmClass))
//			break;

		if(fscanf(fpSvmClass, "%lf\t%lf\t%lf\t%d\t%d\n", tmp, tmp+1, tmp+2, group, group+1)!=5)
		{
//			printf("error!\n");
//			return FAILED;
			break;
		}

		if(group[0]!=group[1])
			errClassNumArray[i] ++;

		i ++;
	}

	fclose(fpSvmClass);

	return SUCCESSFUL;
}

/**
 * Calculate the errorNum of single file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputErrPoints()
{
	int32_t i, j, tmp, rowsNumTestData;
	FILE *fpTestData, *fpErrNum, *fpErrPoint;
	char ch, errPointFile[7][256];

	fpTestData = fopen("../model/testData.txt", "r");
	if(fpTestData==NULL)
	{
		printf("fopen file error!\n");
		return FAILED;
	}

	rowsNumTestData = 0;
	ch = fgetc(fpTestData);
	while(!feof(fpTestData))
	{
		if(ch=='\n')
			rowsNumTestData ++;
		ch = fgetc(fpTestData);
	}

	double testData[rowsNumTestData][3];
	int32_t testDataGroup[rowsNumTestData], errClassNumArray[rowsNumTestData];

	rewind(fpTestData);
	for(i=0; i<rowsNumTestData; i++)
		fscanf(fpTestData, "%lf\t%lf\t%lf\t%d\n", testData[i], testData[i]+1, testData[i]+2, testDataGroup+i);
	fclose(fpTestData);

	fpErrNum = fopen("../model/errNums.txt", "r");
	if(fpErrNum==NULL)
	{
		printf("fopen file error!\n");
		return FAILED;
	}

	i = 0;
	while(1)
	{
		if(fscanf(fpErrNum, "%d\n", &tmp)!=1)
			break;

		errClassNumArray[i++] = tmp;
	}
	fclose(fpErrNum);


	// output the wrongly predicted points
	strcpy(errPointFile[0], "../model/errPoints_0.txt");
	strcpy(errPointFile[1], "../model/errPoints_1.txt");
	strcpy(errPointFile[2], "../model/errPoints_2.txt");
	strcpy(errPointFile[3], "../model/errPoints_3.txt");
	strcpy(errPointFile[4], "../model/errPoints_4.txt");
	strcpy(errPointFile[5], "../model/errPoints_5.txt");
	strcpy(errPointFile[6], "../model/errPoints_6.txt");
	for(i=0; i<7; i++)
	{
		fpErrPoint = fopen(errPointFile[i], "w");
		if(fpErrPoint==NULL)
		{
			printf("error!\n");
			return FAILED;
		}

		for(j=0; j<rowsNumTestData; j++)
			if(errClassNumArray[j]==i)
				fprintf(fpErrPoint, "%d\t%d\t%.4f\t%d\t%d\n", (int32_t)testData[j][0], (int32_t)testData[j][1], testData[j][2], testDataGroup[j], errClassNumArray[j]);

		fclose(fpErrPoint);
	}

	return SUCCESSFUL;
}

/**
 * Test the multiple alignment.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short testMSAMalign()
{

}
