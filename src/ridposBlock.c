/*
 * ridposBlock.c
 *
 *  Created on: Dec 4, 2012
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"



/**
 * Initialize ridpos block in graph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initRidposBlocksInGraph(graphtype *graph)
{
	int64_t i, j, totalNum, sum;
	kmertype *kmer;
	uint32_t blocksNumKmer, itemNumKmerBlock;

	graph->bytesPerRidpos = sizeof(ridpostype);
	graph->blocksNumRidpos = 0;
	graph->totalItemNumRidpos = 0;
	graph->maxBlocksNumRidpos = MAX_BLOCKS_NUM_RIDPOS;
	graph->maxItemNumPerRidposBlock = BLOCK_SIZE_PER_RIDPOS / graph->bytesPerRidpos;


	// add ridpos block head nodes
	graph->ridposBlockArr = (ridposBlock_t *) malloc (graph->maxBlocksNumRidpos * sizeof(ridposBlock_t));
	if(graph->ridposBlockArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// add new ridpos block
	if(addNewBlockRidpos(graph)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// set start ridpos regions and allocate ridpos blocks
	totalNum = sum = 0;
	blocksNumKmer = graph->blocksNumKmer;
	for(i=0; i<blocksNumKmer; i++)
	{
		kmer = graph->kmerBlockArr[i].kmerArr;
		itemNumKmerBlock = graph->kmerBlockArr[i].itemNum;
		for(j=0; j<itemNumKmerBlock; j++)
		{
			kmer->ppos = graph->ridposBlockArr[graph->blocksNumRidpos-1].ridposArr + sum;

			sum += kmer->arraysize;
			if(sum >= graph->maxItemNumPerRidposBlock)
			{
				if(sum > graph->maxItemNumPerRidposBlock)
				{
					sum -= kmer->arraysize;
					graph->ridposBlockArr[graph->blocksNumRidpos-1].itemNum = sum;
					totalNum += sum;

					// add new ridpos block
					if(addNewBlockRidpos(graph)==FAILED)
					{
						printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
						return FAILED;
					}

					kmer->ppos = graph->ridposBlockArr[graph->blocksNumRidpos-1].ridposArr;

					sum = kmer->arraysize;
				}else
				{
					graph->ridposBlockArr[graph->blocksNumRidpos-1].itemNum = sum;
					totalNum += sum;

					// add new ridpos block
					if(addNewBlockRidpos(graph)==FAILED)
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

	// process the last ridpos block
	if(sum==0)
	{ // remove the last ridpos block if it is empty
		free(graph->ridposBlockArr[graph->blocksNumRidpos-1].ridposArr);
		graph->ridposBlockArr[graph->blocksNumRidpos-1].ridposArr = NULL;
		graph->ridposBlockArr[graph->blocksNumRidpos-1].blockID = 0;
		graph->ridposBlockArr[graph->blocksNumRidpos-1].itemNum = 0;

		graph->blocksNumRidpos --;
	}else
	{
		graph->ridposBlockArr[graph->blocksNumRidpos-1].itemNum = sum;

		totalNum += sum;
	}

	graph->totalItemNumRidpos = totalNum;

	return SUCCESSFUL;
}

/**
 * Add new ridpos block.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewBlockRidpos(graphtype *graph)
{
	if(graph->blocksNumRidpos>=graph->maxBlocksNumRidpos)
	{
		graph->ridposBlockArr = (ridposBlock_t *) realloc (graph->ridposBlockArr, graph->maxBlocksNumRidpos*2*sizeof(ridposBlock_t));
		if(graph->ridposBlockArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		graph->maxBlocksNumRidpos *= 2;
	}

	graph->ridposBlockArr[graph->blocksNumRidpos].ridposArr = (ridpostype *) calloc (graph->maxItemNumPerRidposBlock, graph->bytesPerRidpos);
	if(graph->ridposBlockArr[graph->blocksNumRidpos].ridposArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	graph->ridposBlockArr[graph->blocksNumRidpos].blockID = graph->blocksNumRidpos + 1;
	graph->ridposBlockArr[graph->blocksNumRidpos].itemNum = 0;

	graph->blocksNumRidpos ++;

	return SUCCESSFUL;
}

