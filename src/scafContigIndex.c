/*
 * scafContigIndex.c
 *
 *  Created on: Jan 9, 2014
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Build the scafContigIndex.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short buildScafContigIndex(scafContigIndex_t **scafContigIndex, contigGraph_t *contigGraph)
{
	// initialize the scafContigIndex
	if(initScafContigIndex(scafContigIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the scafContigIndex, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// count the scafKmer occurrences
	if(countScafKmerOccs(*scafContigIndex, contigGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot count the kmers from read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// add the contigpos information
	if(addScafKmerContigpos(*scafContigIndex, contigGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot count the kmers from read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Initialize scafContigIndex.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initScafContigIndex(scafContigIndex_t **scafContigIndex)
{
	*scafContigIndex = (scafContigIndex_t *) calloc(1, sizeof(scafContigIndex_t));
	if((*scafContigIndex)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize k-mer block
	if(initScafKmerBlockInGraph(*scafContigIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the scafKmer blocks in k-mer hash table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize kmerseq block
	if(initScafKmerseqBlockInGraph(*scafContigIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the kmerseq blocks in k-mer hash table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Release the contig index.
 */
void releaseScafContigIndex(scafContigIndex_t **scafContigIndex)
{
	int32_t i;

	if((*scafContigIndex)==NULL)
		return;

	// free k-mer blocks
	if((*scafContigIndex)->kmerBlockArr)
	{
		for(i=0; i<(*scafContigIndex)->blocksNumKmer; i++)
			free((*scafContigIndex)->kmerBlockArr[i].kmerArr);
		(*scafContigIndex)->blocksNumKmer = 0;
		free((*scafContigIndex)->kmerBlockArr);
		(*scafContigIndex)->kmerBlockArr = NULL;
	}
	if((*scafContigIndex)->kmerHashtable)
	{
		free((*scafContigIndex)->kmerHashtable);
		(*scafContigIndex)->kmerHashtable = NULL;
	}

	// free kmerseq blocks
	if((*scafContigIndex)->kmerSeqBlockArr)
	{
		for(i=0; i<(*scafContigIndex)->blocksNumKmerSeq; i++)
			free((*scafContigIndex)->kmerSeqBlockArr[i].kmerSeqArr);
		(*scafContigIndex)->blocksNumKmerSeq = 0;
		free((*scafContigIndex)->kmerSeqBlockArr);
		(*scafContigIndex)->kmerSeqBlockArr = NULL;
	}

	// free contigpos blocks
	if((*scafContigIndex)->contigposBlockArr)
	{
		for(i=0; i<(*scafContigIndex)->blocksNumContigpos; i++)
			free((*scafContigIndex)->contigposBlockArr[i].contigposArr);
		(*scafContigIndex)->blocksNumContigpos = 0;
		free((*scafContigIndex)->contigposBlockArr);
		(*scafContigIndex)->contigposBlockArr = NULL;
	}

	// free index node
	free(*scafContigIndex);
	*scafContigIndex = NULL;
}

/**
 * Initialize scafKmer blocks.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initScafKmerBlockInGraph(scafContigIndex_t *scafContigIndex)
{
	scafContigIndex->bytesPerKmer = sizeof(scafKmer_t);
	scafContigIndex->totalItemNumKmer = 0;
	scafContigIndex->maxBlocksNumKmer = MAX_BLOCKS_NUM_KMER;
	scafContigIndex->maxItemNumPerKmerBlock = BLOCK_SIZE_PER_KMER / scafContigIndex->bytesPerKmer;

	scafContigIndex->hashTableSize = hashTableSize;

	scafContigIndex->kmerHashtable = (kmerHashBucket_t *) calloc (scafContigIndex->hashTableSize, sizeof(kmerHashBucket_t));
	if( scafContigIndex->kmerHashtable == NULL )
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate scafKmer blocks
	scafContigIndex->kmerBlockArr = (scafKmerBlock_t *) calloc(scafContigIndex->maxBlocksNumKmer, sizeof(scafKmerBlock_t));
	if( scafContigIndex->kmerBlockArr == NULL )
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	scafContigIndex->blocksNumKmer = 0;

	// add new scafKmer block
	if(addNewBlockScafKmer(scafContigIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Add new scafKmer block.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewBlockScafKmer(scafContigIndex_t *scafContigIndex)
{
	if(scafContigIndex->blocksNumKmer>=scafContigIndex->maxBlocksNumKmer)
	{
		scafContigIndex->kmerBlockArr = (scafKmerBlock_t *) realloc (scafContigIndex->kmerBlockArr, scafContigIndex->maxBlocksNumKmer*2*sizeof(scafKmerBlock_t));
		if(scafContigIndex->kmerBlockArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		scafContigIndex->maxBlocksNumKmer *= 2;
	}

	scafContigIndex->kmerBlockArr[scafContigIndex->blocksNumKmer].kmerArr = (scafKmer_t *) calloc (scafContigIndex->maxItemNumPerKmerBlock, scafContigIndex->bytesPerKmer);
	if(scafContigIndex->kmerBlockArr[scafContigIndex->blocksNumKmer].kmerArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	scafContigIndex->kmerBlockArr[scafContigIndex->blocksNumKmer].blockID = scafContigIndex->blocksNumKmer + 1;
	scafContigIndex->kmerBlockArr[scafContigIndex->blocksNumKmer].itemNum = 0;

	scafContigIndex->blocksNumKmer ++;

	return SUCCESSFUL;
}

/**
 * Initialize kmerseq blocks in contig index.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initScafKmerseqBlockInGraph(scafContigIndex_t *scafContigIndex)
{
	scafContigIndex->entriesPerKmer = entriesPerKmer;
	scafContigIndex->bytesPerKmerseq = entriesPerKmer * sizeof(uint64_t);
	scafContigIndex->totalItemNumKmerSeq = 0;
	scafContigIndex->maxBlocksNumKmerSeq = MAX_BLOCKS_NUM_KMER_SEQ;
	scafContigIndex->maxItemNumPerKmerSeqBlock = BLOCK_SIZE_PER_KMER_SEQ / scafContigIndex->bytesPerKmerseq;

	// allocate kmerseq blocks
	scafContigIndex->kmerSeqBlockArr = (kmerseqBlock_t *) malloc(scafContigIndex->maxBlocksNumKmerSeq * sizeof(kmerseqBlock_t));
	if( scafContigIndex->kmerSeqBlockArr == NULL )
	{
		printf("line=%d, In %s(), cannot allocate the kmerSeqBlockArr for the k-mer hash table, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	scafContigIndex->blocksNumKmerSeq = 0;

	// add new kmerseq block
	if(addNewBlockScafKmerSeq(scafContigIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Add new kmerseq block in contig index.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewBlockScafKmerSeq(scafContigIndex_t *scafContigIndex)
{
	if(scafContigIndex->blocksNumKmerSeq>=scafContigIndex->maxBlocksNumKmerSeq)
	{
		scafContigIndex->kmerSeqBlockArr = (kmerseqBlock_t *) realloc (scafContigIndex->kmerSeqBlockArr, scafContigIndex->maxBlocksNumKmerSeq*2*sizeof(kmerseqBlock_t));
		if(scafContigIndex->kmerSeqBlockArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		scafContigIndex->maxBlocksNumKmerSeq *= 2;
	}

	scafContigIndex->kmerSeqBlockArr[scafContigIndex->blocksNumKmerSeq].kmerSeqArr = (uint64_t *) calloc (scafContigIndex->maxItemNumPerKmerSeqBlock, scafContigIndex->bytesPerKmerseq);
	if(scafContigIndex->kmerSeqBlockArr[scafContigIndex->blocksNumKmerSeq].kmerSeqArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	scafContigIndex->kmerSeqBlockArr[scafContigIndex->blocksNumKmerSeq].blockID = scafContigIndex->blocksNumKmerSeq + 1;
	scafContigIndex->kmerSeqBlockArr[scafContigIndex->blocksNumKmerSeq].itemNum = 0;

	scafContigIndex->blocksNumKmerSeq ++;

	return SUCCESSFUL;
}

/**
 * Count the scafKmer occurrences.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short countScafKmerOccs(scafContigIndex_t *scafContigIndex, contigGraph_t *contigGraph)
{
	int32_t i, j, contigLen, startContigPos, endContigPos, basePos, baseInt, contigEndFlag;
	char *contigSeq;
	uint64_t *pKmerSeqIntDone, *pKmerSeqIntDoing, hashcode;
	kmerseqBlock_t *pKmerseqBlock;

	pKmerseqBlock = scafContigIndex->kmerSeqBlockArr + scafContigIndex->blocksNumKmerSeq - 1;
	scafContigIndex->pKmerSeqAdding = pKmerseqBlock->kmerSeqArr;

	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
		contigSeq = contigGraph->contigItemArray[i].contigSeq;
		contigLen = contigGraph->contigItemArray[i].contigLen;

		// 5' end
		startContigPos = 0;
		if(contigLen>2*contigAlignRegSize)
		{
			endContigPos = contigAlignRegSize - 1;
			contigEndFlag = 0;
			contigGraph->contigItemArray[i].alignRegSizeEnd5 = contigAlignRegSize;
		}
		else if(contigLen>contigAlignRegSize)
		{
			endContigPos = (contigLen - 1) / 2;
			contigEndFlag = 0;
			contigGraph->contigItemArray[i].alignRegSizeEnd5 = endContigPos + 1;
		}else
		{
			endContigPos = contigLen - 1;
			contigEndFlag = 2;	// only the 5' end bacause of the short (<=2000 bp) contig
			contigGraph->contigItemArray[i].alignRegSizeEnd5 = contigLen;
		}

		pKmerSeqIntDoing = scafContigIndex->pKmerSeqAdding;
		// generate the first kmerseqInt
		if(generateReadseqInt(pKmerSeqIntDoing, contigSeq, kmerSize, entriesPerKmer)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the integer sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}
		pKmerSeqIntDone = pKmerSeqIntDoing;

		// count the first scafKmer
		hashcode = kmerhashInt(pKmerSeqIntDoing);
		if(countScafKmer(hashcode, pKmerSeqIntDoing, scafContigIndex)==FAILED)
		{
			printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// process other scafKmers
		for(basePos=kmerSize; basePos<=endContigPos; basePos++)
		{
			switch(contigSeq[basePos])
			{
				case 'A':
				case 'a': baseInt = 0; break;
				case 'C':
				case 'c': baseInt = 1; break;
				case 'G':
				case 'g': baseInt = 2; break;
				case 'T':
				case 't': baseInt = 3; break;
				default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, contigSeq[basePos]); return FAILED;
			}

			pKmerSeqIntDoing = scafContigIndex->pKmerSeqAdding;

			// generate the kmer integer sequence
			if(entriesPerKmer>=2)
			{
				for(j=0; j<entriesPerKmer-2; j++)
				{
					pKmerSeqIntDoing[j] = (pKmerSeqIntDone[j] << 2) | (pKmerSeqIntDone[j+1] >> 62);
				}
				pKmerSeqIntDoing[entriesPerKmer-2] = (pKmerSeqIntDone[entriesPerKmer-2] << 2) | (pKmerSeqIntDone[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
			}
			pKmerSeqIntDoing[entriesPerKmer-1] = ((pKmerSeqIntDone[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

			pKmerSeqIntDone = pKmerSeqIntDoing;

			hashcode = kmerhashInt(pKmerSeqIntDoing);
			if(countScafKmer(hashcode, pKmerSeqIntDoing, scafContigIndex)==FAILED)
			{
				printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		// 3' end
		if(contigLen<=contigAlignRegSize)
		{ // do not the 3' end
			contigGraph->contigItemArray[i].alignRegSizeEnd3 = 0;
			continue;
		}else if(contigLen<=2*contigAlignRegSize)
		{ // handle the 3' end
			startContigPos = contigLen - (contigLen - 1) / 2 - 1;
			contigEndFlag = 1;  // the 3' end of the contig
			contigGraph->contigItemArray[i].alignRegSizeEnd3 = contigLen - startContigPos;
		}else
		{ // handle the 3' end
			startContigPos = contigLen - contigAlignRegSize;
			contigEndFlag = 1;  // the 3' end of the contig
			contigGraph->contigItemArray[i].alignRegSizeEnd3 = contigAlignRegSize;
		}
		endContigPos = contigLen - 1;

		pKmerSeqIntDoing = scafContigIndex->pKmerSeqAdding;
		// generate the first kmerseqInt
		if(generateReadseqInt(pKmerSeqIntDoing, contigSeq+startContigPos, kmerSize, entriesPerKmer)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the integer sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}
		pKmerSeqIntDone = pKmerSeqIntDoing;

		// count the first scafKmer
		hashcode = kmerhashInt(pKmerSeqIntDoing);
		if(countScafKmer(hashcode, pKmerSeqIntDoing, scafContigIndex)==FAILED)
		{
			printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// process other scafKmers
		for(basePos=startContigPos+kmerSize; basePos<=endContigPos; basePos++)
		{
			switch(contigSeq[basePos])
			{
				case 'A':
				case 'a': baseInt = 0; break;
				case 'C':
				case 'c': baseInt = 1; break;
				case 'G':
				case 'g': baseInt = 2; break;
				case 'T':
				case 't': baseInt = 3; break;
				default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, contigSeq[basePos]); return FAILED;
			}

			pKmerSeqIntDoing = scafContigIndex->pKmerSeqAdding;

			// generate the kmer integer sequence
			if(entriesPerKmer>=2)
			{
				for(j=0; j<entriesPerKmer-2; j++)
				{
					pKmerSeqIntDoing[j] = (pKmerSeqIntDone[j] << 2) | (pKmerSeqIntDone[j+1] >> 62);
				}
				pKmerSeqIntDoing[entriesPerKmer-2] = (pKmerSeqIntDone[entriesPerKmer-2] << 2) | (pKmerSeqIntDone[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
			}
			pKmerSeqIntDoing[entriesPerKmer-1] = ((pKmerSeqIntDone[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

			pKmerSeqIntDone = pKmerSeqIntDoing;

			hashcode = kmerhashInt(pKmerSeqIntDoing);
			if(countScafKmer(hashcode, pKmerSeqIntDoing, scafContigIndex)==FAILED)
			{
				printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	// the last k-mer block is empty, remove it
	if(scafContigIndex->kmerBlockArr[scafContigIndex->blocksNumKmer-1].itemNum==0)
	{
		free(scafContigIndex->kmerBlockArr[scafContigIndex->blocksNumKmer-1].kmerArr);
		scafContigIndex->kmerBlockArr[scafContigIndex->blocksNumKmer-1].kmerArr = NULL;
		scafContigIndex->blocksNumKmer --;
	}

	// the last kmerseq block is empty, remove it
	if(scafContigIndex->kmerSeqBlockArr[scafContigIndex->blocksNumKmerSeq-1].itemNum==0)
	{
		free(scafContigIndex->kmerSeqBlockArr[scafContigIndex->blocksNumKmerSeq-1].kmerSeqArr);
		scafContigIndex->kmerSeqBlockArr[scafContigIndex->blocksNumKmerSeq-1].kmerSeqArr = NULL;
		scafContigIndex->blocksNumKmerSeq --;
	}

	return SUCCESSFUL;
}

/**
 * Count the k-mer occurrences.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short countScafKmer(uint64_t hashcode, uint64_t *kmerSeqInt, scafContigIndex_t *scafContigIndex)
{
	scafKmer_t *kmer;
	scafKmerBlock_t *pKmerBlock;
	kmerseqBlock_t *pKmerseqBlock;

	kmer = getScafKmerByHash(hashcode, kmerSeqInt, scafContigIndex);
	if(!kmer)
	{
		// process k-mer blocks
		pKmerBlock = scafContigIndex->kmerBlockArr + scafContigIndex->blocksNumKmer - 1;
		kmer = pKmerBlock->kmerArr + pKmerBlock->itemNum;

		kmer->kmerseqBlockID = pKmerBlock->blockID;
		kmer->itemRowKmerseqBlock = pKmerBlock->itemNum;
		kmer->arraysize = 1;
		kmer->nextKmerBlockID = scafContigIndex->kmerHashtable[hashcode].kmerBlockID;
		kmer->nextItemRowKmerBlock = scafContigIndex->kmerHashtable[hashcode].itemRowKmerBlock;
		scafContigIndex->kmerHashtable[hashcode].kmerBlockID = kmer->kmerseqBlockID;
		scafContigIndex->kmerHashtable[hashcode].itemRowKmerBlock = kmer->itemRowKmerseqBlock;

		pKmerBlock->itemNum ++;
		scafContigIndex->totalItemNumKmer ++;

		if(pKmerBlock->itemNum >= scafContigIndex->maxItemNumPerKmerBlock)
		{
			// add new kmer block
			if(addNewBlockScafKmer(scafContigIndex)==FAILED)
			{
				printf("line=%d, In %s(), cannot add new k-mer block, Error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		// process kmerseq blocks
		pKmerseqBlock = scafContigIndex->kmerSeqBlockArr + scafContigIndex->blocksNumKmerSeq - 1;
		kmer->kmerseqBlockID = pKmerseqBlock->blockID;
		kmer->itemRowKmerseqBlock = pKmerseqBlock->itemNum;
		pKmerseqBlock->itemNum ++;
		scafContigIndex->totalItemNumKmerSeq ++;
		scafContigIndex->pKmerSeqAdding += scafContigIndex->entriesPerKmer;

		if(pKmerseqBlock->itemNum >= scafContigIndex->maxItemNumPerKmerSeqBlock)
		{
			// add new kmerseq block
			if(addNewBlockScafKmerSeq(scafContigIndex)==FAILED)
			{
				printf("line=%d, In %s(), cannot add new kmerseq block, Error!\n", __LINE__, __func__);
				return FAILED;
			}
			scafContigIndex->pKmerSeqAdding = scafContigIndex->kmerSeqBlockArr[scafContigIndex->blocksNumKmerSeq-1].kmerSeqArr;
		}
	}else
	{
		kmer->arraysize ++;
	}

	return SUCCESSFUL;
}

scafKmer_t *getScafKmerByHash(uint64_t hashvalue, uint64_t *kmerSeqInt, scafContigIndex_t *scafContigIndex)
{
	scafKmer_t *kmer;
	uint64_t *kmerseq;
	kmerHashBucket_t *pKmerBucket;

	pKmerBucket = scafContigIndex->kmerHashtable + hashvalue;
	if(pKmerBucket->kmerBlockID>0)
	{
		kmer = scafContigIndex->kmerBlockArr[pKmerBucket->kmerBlockID-1].kmerArr + pKmerBucket->itemRowKmerBlock;
		while(kmer)
		{
			kmerseq = scafContigIndex->kmerSeqBlockArr[kmer->kmerseqBlockID-1].kmerSeqArr + kmer->itemRowKmerseqBlock * scafContigIndex->entriesPerKmer;
			if(identicalKmerSeq(kmerSeqInt, kmerseq)==YES)
				break;

			if(kmer->nextKmerBlockID>0)
				kmer = scafContigIndex->kmerBlockArr[kmer->nextKmerBlockID-1].kmerArr + kmer->nextItemRowKmerBlock;
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
 * Add the k-mer contigpos information.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addScafKmerContigpos(scafContigIndex_t *scafContigIndex, contigGraph_t *contigGraph)
{
	int32_t i, j, contigID, contigLen, startContigPos, contigPos, endContigPos, basePos, baseInt, contigEndFlag;
	char *contigSeq;
	uint64_t pKmerSeqInt[entriesPerKmer], hashcode;

	// initialize contigpos blocks and the ridpos regions in that blocks
	if(initContigposBlocksInContigIndex(scafContigIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<contigGraph->itemNumContigItemArray; i++)
	{
		contigID = contigGraph->contigItemArray[i].contigID;
		contigSeq = contigGraph->contigItemArray[i].contigSeq;
		contigLen = contigGraph->contigItemArray[i].contigLen;

		// 5' end
		startContigPos = 0;
		contigPos = 1;
		endContigPos = contigGraph->contigItemArray[i].alignRegSizeEnd5 - 1;
		if(contigGraph->contigItemArray[i].alignRegSizeEnd3>0)
			contigEndFlag = 0;
		else
			contigEndFlag = 2;

		// generate the first kmerseqInt
		if(generateReadseqInt(pKmerSeqInt, contigSeq, kmerSize, entriesPerKmer)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the integer sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// count the first scafKmer
		hashcode = kmerhashInt(pKmerSeqInt);
		if(addScafKmer(hashcode, pKmerSeqInt, contigID, contigPos, contigEndFlag, scafContigIndex)==FAILED)
		{
			printf("line=%d, In %s(), cannot add kmer occurrence, error!\n", __LINE__, __func__);
			return FAILED;
		}
		contigPos ++;

		// process other scafKmers
		for(basePos=kmerSize; basePos<=endContigPos; basePos++, contigPos++)
		{
			switch(contigSeq[basePos])
			{
				case 'A':
				case 'a': baseInt = 0; break;
				case 'C':
				case 'c': baseInt = 1; break;
				case 'G':
				case 'g': baseInt = 2; break;
				case 'T':
				case 't': baseInt = 3; break;
				default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, contigSeq[basePos]); return FAILED;
			}

			// generate the kmer integer sequence
			if(entriesPerKmer>=2)
			{
				for(j=0; j<entriesPerKmer-2; j++)
				{
					pKmerSeqInt[j] = (pKmerSeqInt[j] << 2) | (pKmerSeqInt[j+1] >> 62);
				}
				pKmerSeqInt[entriesPerKmer-2] = (pKmerSeqInt[entriesPerKmer-2] << 2) | (pKmerSeqInt[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
			}
			pKmerSeqInt[entriesPerKmer-1] = ((pKmerSeqInt[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

			hashcode = kmerhashInt(pKmerSeqInt);
			if(addScafKmer(hashcode, pKmerSeqInt, contigID, contigPos, contigEndFlag, scafContigIndex)==FAILED)
			{
				printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		// 3' end
		if(contigGraph->contigItemArray[i].alignRegSizeEnd3<=0)
		{
			continue;
		}else
		{
			contigEndFlag = 1;  // the 3' end of the contig
		}

		//startContigPos = contigLen - contigAlignRegSize;
		startContigPos = contigLen - contigGraph->contigItemArray[i].alignRegSizeEnd3;
		endContigPos = contigLen - 1;
		contigPos = startContigPos + 1;

		// generate the first kmerseqInt
		if(generateReadseqInt(pKmerSeqInt, contigSeq+startContigPos, kmerSize, entriesPerKmer)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the integer sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// count the first scafKmer
		hashcode = kmerhashInt(pKmerSeqInt);
		if(addScafKmer(hashcode, pKmerSeqInt, contigID, contigPos, contigEndFlag, scafContigIndex)==FAILED)
		{
			printf("line=%d, In %s(), cannot add kmer occurrence, error!\n", __LINE__, __func__);
			return FAILED;
		}
		contigPos ++;

		// process other scafKmers
		for(basePos=startContigPos+kmerSize; basePos<=endContigPos; basePos++, contigPos++)
		{
			switch(contigSeq[basePos])
			{
				case 'A':
				case 'a': baseInt = 0; break;
				case 'C':
				case 'c': baseInt = 1; break;
				case 'G':
				case 'g': baseInt = 2; break;
				case 'T':
				case 't': baseInt = 3; break;
				default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, contigSeq[basePos]); return FAILED;
			}

			// generate the kmer integer sequence
			if(entriesPerKmer>=2)
			{
				for(j=0; j<entriesPerKmer-2; j++)
				{
					pKmerSeqInt[j] = (pKmerSeqInt[j] << 2) | (pKmerSeqInt[j+1] >> 62);
				}
				pKmerSeqInt[entriesPerKmer-2] = (pKmerSeqInt[entriesPerKmer-2] << 2) | (pKmerSeqInt[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
			}
			pKmerSeqInt[entriesPerKmer-1] = ((pKmerSeqInt[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

			hashcode = kmerhashInt(pKmerSeqInt);
			if(addScafKmer(hashcode, pKmerSeqInt, contigID, contigPos, contigEndFlag, scafContigIndex)==FAILED)
			{
				printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Initialize contigpos block in contig index.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initContigposBlocksInContigIndex(scafContigIndex_t *scafContigIndex)
{
	int64_t i, j, totalNum, sum;
	scafKmer_t *kmer;
	uint32_t blocksNumKmer, itemNumKmerBlock;

	scafContigIndex->bytesPerContigpos = sizeof(scafContigposBlock_t);
	scafContigIndex->blocksNumContigpos = 0;
	scafContigIndex->totalItemNumContigpos = 0;
	scafContigIndex->maxBlocksNumContigpos = MAX_BLOCKS_NUM_RIDPOS;
	scafContigIndex->maxItemNumPerContigposBlock = BLOCK_SIZE_PER_RIDPOS / scafContigIndex->bytesPerContigpos;

	// add contigpos block head nodes
	scafContigIndex->contigposBlockArr = (scafContigposBlock_t *) malloc (scafContigIndex->maxBlocksNumContigpos * sizeof(scafContigposBlock_t));
	if(scafContigIndex->contigposBlockArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// add new contigpos block
	if(addNewBlockContigpos(scafContigIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// set start ridpos regions and allocate ridpos blocks
	totalNum = sum = 0;
	blocksNumKmer = scafContigIndex->blocksNumKmer;
	for(i=0; i<blocksNumKmer; i++)
	{
		kmer = scafContigIndex->kmerBlockArr[i].kmerArr;
		itemNumKmerBlock = scafContigIndex->kmerBlockArr[i].itemNum;
		for(j=0; j<itemNumKmerBlock; j++)
		{
			kmer->ppos = scafContigIndex->contigposBlockArr[scafContigIndex->blocksNumContigpos-1].contigposArr + sum;

			sum += kmer->arraysize;
			if(sum >= scafContigIndex->maxItemNumPerContigposBlock)
			{
				if(sum > scafContigIndex->maxItemNumPerContigposBlock)
				{
					sum -= kmer->arraysize;
					scafContigIndex->contigposBlockArr[scafContigIndex->blocksNumContigpos-1].itemNum = sum;
					totalNum += sum;

					// add new ridpos block
					if(addNewBlockContigpos(scafContigIndex)==FAILED)
					{
						printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
						return FAILED;
					}

					kmer->ppos = scafContigIndex->contigposBlockArr[scafContigIndex->blocksNumContigpos-1].contigposArr;

					sum = kmer->arraysize;
				}else
				{
					scafContigIndex->contigposBlockArr[scafContigIndex->blocksNumContigpos-1].itemNum = sum;
					totalNum += sum;

					// add new contigpos block
					if(addNewBlockContigpos(scafContigIndex)==FAILED)
					{
						printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
						return FAILED;
					}
					sum = 0;
				}
			}
			kmer ++;
		}
	}

	// process the last contigpos block
	if(sum==0)
	{ // remove the last contigpos block if it is empty
		free(scafContigIndex->contigposBlockArr[scafContigIndex->blocksNumContigpos-1].contigposArr);
		scafContigIndex->contigposBlockArr[scafContigIndex->blocksNumContigpos-1].contigposArr = NULL;
		scafContigIndex->contigposBlockArr[scafContigIndex->blocksNumContigpos-1].blockID = 0;
		scafContigIndex->contigposBlockArr[scafContigIndex->blocksNumContigpos-1].itemNum = 0;

		scafContigIndex->blocksNumContigpos --;
	}else
	{
		scafContigIndex->contigposBlockArr[scafContigIndex->blocksNumContigpos-1].itemNum = sum;

		totalNum += sum;
	}

	scafContigIndex->totalItemNumContigpos = totalNum;

	return SUCCESSFUL;
}

/**
 * Add new contigpos block.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewBlockContigpos(scafContigIndex_t *scafContigIndex)
{
	if(scafContigIndex->blocksNumContigpos>=scafContigIndex->maxBlocksNumContigpos)
	{
		scafContigIndex->contigposBlockArr = (scafContigposBlock_t *) realloc (scafContigIndex->contigposBlockArr, scafContigIndex->maxBlocksNumContigpos*2*sizeof(scafContigposBlock_t));
		if(scafContigIndex->contigposBlockArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		scafContigIndex->maxBlocksNumContigpos *= 2;
	}

	scafContigIndex->contigposBlockArr[scafContigIndex->blocksNumContigpos].contigposArr = (scafContigpos_t *) calloc (scafContigIndex->maxItemNumPerContigposBlock, scafContigIndex->bytesPerContigpos);
	if(scafContigIndex->contigposBlockArr[scafContigIndex->blocksNumContigpos].contigposArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	scafContigIndex->contigposBlockArr[scafContigIndex->blocksNumContigpos].blockID = scafContigIndex->blocksNumContigpos + 1;
	scafContigIndex->contigposBlockArr[scafContigIndex->blocksNumContigpos].itemNum = 0;

	scafContigIndex->blocksNumContigpos ++;

	return SUCCESSFUL;
}

/**
 * Add a kmer to graph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addScafKmer(uint64_t hashcode, uint64_t *kmerSeqInt, int32_t contigID, int32_t contigPos, int32_t contigEndFlag, scafContigIndex_t *scafContigIndex)
{
	scafKmer_t *kmer;
	scafContigpos_t *contigpos;

	kmer = getScafKmerByHash(hashcode, kmerSeqInt, scafContigIndex);
	if(kmer)
	{
		contigpos = kmer->ppos + kmer->multiplicity;
		contigpos->contigID = contigID;
		contigpos->contigpos = contigPos;
		contigpos->contigEndFlag = contigEndFlag;

		kmer->multiplicity ++;

	}else
	{
		printf("line=%d, In %s(), can not get the k-mer %s, error!\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt));
		return FAILED;
	}

	return SUCCESSFUL;
}

