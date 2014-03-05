/*
 * kmerBlock.c
 *
 *  Created on: Dec 5, 2012
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"



/**
 * Initialize k-mer blocks.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initKmerBlockInGraph(graphtype *graph)
{
	graph->bytesPerKmer = sizeof(kmertype);
	graph->totalItemNumKmer = 0;
	graph->maxBlocksNumKmer = MAX_BLOCKS_NUM_KMER;
	graph->maxItemNumPerKmerBlock = BLOCK_SIZE_PER_KMER / graph->bytesPerKmer;
	graph->kmerSampleInterval = KMER_SAMPLE_INTERVAL;

	graph->hashTableSize = hashTableSize;

	graph->kmerHashtable = (kmerHashBucket_t *) calloc (graph->hashTableSize, sizeof(kmerHashBucket_t));
	if( graph->kmerHashtable == NULL )
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate k-mer blocks
	graph->kmerBlockArr = (kmerBlock_t *) calloc(graph->maxBlocksNumKmer, sizeof(kmerBlock_t));
	if( graph->kmerBlockArr == NULL )
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	graph->blocksNumKmer = 0;

	// add new k-mer block
	if(addNewBlockKmer(graph)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}


/**
 * Add new k-mer block.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewBlockKmer(graphtype *graph)
{
	if(graph->blocksNumKmer>=graph->maxBlocksNumKmer)
	{
		graph->kmerBlockArr = (kmerBlock_t *) realloc (graph->kmerBlockArr, graph->maxBlocksNumKmer*2*sizeof(kmerBlock_t));
		if(graph->kmerBlockArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		graph->maxBlocksNumKmer *= 2;
	}

	graph->kmerBlockArr[graph->blocksNumKmer].kmerArr = (kmertype *) calloc (graph->maxItemNumPerKmerBlock, graph->bytesPerKmer);
	if(graph->kmerBlockArr[graph->blocksNumKmer].kmerArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	graph->kmerBlockArr[graph->blocksNumKmer].blockID = graph->blocksNumKmer + 1;
	graph->kmerBlockArr[graph->blocksNumKmer].itemNum = 0;

	graph->blocksNumKmer ++;

	return SUCCESSFUL;
}
