/*
 * kmerseqBlock.c
 *
 *  Created on: Dec 3, 2012
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"



/**
 * Initialize kmerseq blocks.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initKmerseqBlockInGraph(graphtype *graph)
{
	graph->entriesPerKmer = entriesPerKmer;
	graph->bytesPerKmerseq = entriesPerKmer * sizeof(uint64_t);
	graph->totalItemNumKmerSeq = 0;
	graph->maxBlocksNumKmerSeq = MAX_BLOCKS_NUM_KMER_SEQ;
	graph->maxItemNumPerKmerSeqBlock = BLOCK_SIZE_PER_KMER_SEQ / graph->bytesPerKmerseq;

	// allocate kmerseq blocks
	graph->kmerSeqBlockArr = (kmerseqBlock_t *) malloc(graph->maxBlocksNumKmerSeq * sizeof(kmerseqBlock_t));
	if( graph->kmerSeqBlockArr == NULL )
	{
		printf("line=%d, In %s(), cannot allocate the kmerSeqBlockArr for the k-mer hash table, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	graph->blocksNumKmerSeq = 0;

	// add new kmerseq block
	if(addNewBlockKmerSeq(graph)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}


/**
 * Add new kmerseq block.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewBlockKmerSeq(graphtype *graph)
{
	if(graph->blocksNumKmerSeq>=graph->maxBlocksNumKmerSeq)
	{
		graph->kmerSeqBlockArr = (kmerseqBlock_t *) realloc (graph->kmerSeqBlockArr, graph->maxBlocksNumKmerSeq*2*sizeof(kmerseqBlock_t));
		if(graph->kmerSeqBlockArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		graph->maxBlocksNumKmerSeq *= 2;
	}

	graph->kmerSeqBlockArr[graph->blocksNumKmerSeq].kmerSeqArr = (uint64_t *) calloc (graph->maxItemNumPerKmerSeqBlock, graph->bytesPerKmerseq);
	if(graph->kmerSeqBlockArr[graph->blocksNumKmerSeq].kmerSeqArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	graph->kmerSeqBlockArr[graph->blocksNumKmerSeq].blockID = graph->blocksNumKmerSeq + 1;
	graph->kmerSeqBlockArr[graph->blocksNumKmerSeq].itemNum = 0;

	graph->blocksNumKmerSeq ++;

	return SUCCESSFUL;
}
