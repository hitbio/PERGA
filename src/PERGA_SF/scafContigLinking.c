/*
 * scafContigLinking.c
 *
 *  Created on: Jun 18, 2011
 *      Author: zhuxiao
 */

#include "inc/scafStdinc.h"
#include "inc/scafExtvab.h"


/**
 * Contigs scaffolding.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short contigsLinking(const char *linkResultFile, const char *averLinkNumFile, const char *contigFileName, const char *sharedReadListFile, const char *contigListFile)
{
	printf("=========== Begin linking contigs, please wait ...\n");

	// initialize memory for linking of contigs
	if(initMemLinking(contigFileName, sharedReadListFile, contigListFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// construct the contigGraph
	if(constructContigGraph()==FAILED)
	{
		printf("line=%d, In %s(), cannot construct contigGraoh, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize the parameters for contigs linking
	if(initParaLinking(averLinkNumFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the parameters for contigs linking, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// begin linking contigs
	if(linkContigs(linkResultFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot link contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// free memory for linking of contigs
	freeMemLinking();

	printf("=========== End linked contigs.\n");

	return SUCCESSFUL;
}

/**
 * Initialize memory for the linking of contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */

short initMemLinking(const char *contigFileName, const char *sharedReadListFile, const char *contigListFile)
{
	// Allocate memory and load the data of shared read list (SRL)
	if(loadSingleReadList(sharedReadListFile, &sharedReadListArr, &readItemNumInSRL, &sharedReadPosArr, &matchItemNumInSRP)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the data of shared read list (SRL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Allocate memory and load the data of contig list (CL)
	if(loadContigList(contigListFile, &contigListArr, &contigItemNumInCL, &contigReadArr, &contigReadItemNumInCR)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the data of shared read list (SRL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize the contig information array
	if(initContigInfoArray(contigFileName)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the contig information array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize the contig link array
	if(initContigLinkArr(contigsNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the contig link array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize the contigGraph
	if(initMemContigGraph(contigsNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the contig graph, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Release memory for linking of contigs.
 */
void freeMemLinking()
{
	// reset the linking variables
	minLinksNumContigsThres = 0;
	maxRatioSecondFirstLinkNum = 0;
	secondLinkNumFactor = 0;

	// free shared read list (SRL)
	freeSingleReadList(&sharedReadListArr, &readItemNumInSRL, &sharedReadPosArr, &matchItemNumInSRP);

	// free contig list (CL)
	freeContigList(&contigListArr, &contigItemNumInCL, &contigReadArr, &contigReadItemNumInCR);

	// free contig information array
	freeMemContigInfo(&contigInfoArr, &contigsNum);

	// free contig link array
	freeMemContigLinkArr();

	// free ContigGraph
	freeMemContigGraph();
}

/**
 * Initialize the contig information array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initContigInfoArray(const char *contigFileName)
{
	//get the contigs number
	contigsNum = getContigsNum(contigFileName);
	if(contigsNum<=0)
	{
		printf("line=%d, In %s(), contigsNum=%d, error!\n", __LINE__, __func__, contigsNum);
		return FAILED;
	}

	contigInfoArr = (contigInfo *) calloc(contigsNum, sizeof(contigInfo));
	if(contigInfoArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}


	if(fillContigInfoArr(contigFileName, contigInfoArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigFileName);
		return FAILED;
	}


	return SUCCESSFUL;
}

/**
 * Release the memory of contig information array.
 */
void freeMemContigInfo(contigInfo **pContigInfoArr, int *contigsNum)
{
	int i;
	for(i=0; i<*contigsNum; i++)
	{
		free((*pContigInfoArr)[i].contigSeq);
		(*pContigInfoArr)[i].contigSeq = NULL;
		free((*pContigInfoArr)[i].headTitle);
		(*pContigInfoArr)[i].headTitle = NULL;
	}

	free((*pContigInfoArr));
	*pContigInfoArr = NULL;

	*contigsNum = 0;
}

/**
 * Initialize contig link array and its memory.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initContigLinkArr(const int contigsNum)
{
	itemNumContigLinkArr = 0;
	headRowContigLinkArr = 0;
	tailRowContigLinkArr = 0;

	contigLinkArr = (contigLink *) calloc(contigsNum, sizeof(contigLink));
	if(contigLinkArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Release memory contig link array.
 */
void freeMemContigLinkArr()
{
	itemNumContigLinkArr = 0;
	headRowContigLinkArr = 0;
	tailRowContigLinkArr = 0;

	free(contigLinkArr);
	contigLinkArr = NULL;
}

/**
 * Initialize memory for contigGraph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemContigGraph(const int contigsNum)
{
	itemNumContigGraphArr = 2 * contigsNum;
	contigGraphArr = (ContigGraph*) calloc (itemNumContigGraphArr, sizeof(ContigGraph));
	if(contigGraphArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory for contigGraphArr, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Free the memory of contigGraph and tripletArray.
 */
void freeMemContigGraph()
{
	int64_t i;

	for(i=0; i<itemNumContigGraphArr; i++)
	{
		if(contigGraphArr[i].arraySize>0)
		{
			free(contigGraphArr[i].pEdgeArray);
			contigGraphArr[i].pEdgeArray = NULL;
		}
	}

	free(contigGraphArr);
	contigGraphArr = NULL;
	itemNumContigGraphArr = 0;
}

/**
 * construct the contig graph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructContigGraph()
{
	// initialize the memory for arrLocArr
	if(initMemArrLocArr()==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize memory for arrLocArr, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the arrLocArr
	if(fillArrLocArr()==FAILED)
	{
		printf("line=%d, In %s(), cannot fill arrLocArr, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// sort the arrLocArr
	if(sortArrLocArr(arrLocArr, itemNumArrLocArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot sort arrLocArr, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// generate the contigGraph
	if(generateContigGraphEdge(contigGraphArr, itemNumContigGraphArr, arrLocArr, itemNumArrLocArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate contigGraph, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(FillContigGraphEdge(contigGraphArr, itemNumContigGraphArr, contigInfoArr, contigListArr, contigItemNumInCL, contigReadArr, sharedReadListArr, readItemNumInSRL, sharedReadPosArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot Fill contigGraphEdges, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(removeInvalidGraphEdge(contigGraphArr, itemNumContigGraphArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot remove invalid contigGraphEdges, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ###################### Debug information #######################
#if DEBUG_FLAG
//	if(checkArrLocArray(arrLocArr, itemNumArrLocArr)==NO)
//	{
//		printf("line=%d, In %s(), sort error!\n", __LINE__, __func__);
//		return FAILED;
//	}
#endif
	// ###################### Debug information #######################

	// fill the arrLocArr
	freeMemArrLocArr();

	return SUCCESSFUL;
}

/**
 * Initialize memory for arrLocArr.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemArrLocArr()
{
	maxItemNumArrLocArr = matchItemNumInSRP;
	arrLocArr = (arrLoc*) calloc (maxItemNumArrLocArr, sizeof(arrLoc));
	if(arrLocArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory for arrLocArr, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Free the memory of contigGraph and tripletArray.
 */
void freeMemArrLocArr()
{
	free(arrLocArr);
	arrLocArr = NULL;
	itemNumArrLocArr = 0;
}

/**
 * Fill the arrLocArr.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillArrLocArr()
{
	int64_t i;
	int firstRow_1, firstRow_2, endFlag_1, endFlag_2, contigID_1, contigID_2;
	int row, col;

	itemNumArrLocArr = 0;
	for(i=0; i<readItemNumInSRL; i+=2)
	{
		// Only the paired reads both having just one occurrence will be considered.
		if(sharedReadListArr[i].curNum==1 && sharedReadListArr[i+1].curNum==1)
		{ // only the unique reads are taken consideration
			firstRow_1 = sharedReadListArr[i].firstRow;
			firstRow_2 = sharedReadListArr[i+1].firstRow;

			contigID_1 = sharedReadPosArr[firstRow_1].contigID;
			contigID_2 = sharedReadPosArr[firstRow_2].contigID;
			endFlag_1 = sharedReadPosArr[firstRow_1].contigEnd;
			endFlag_2 = sharedReadPosArr[firstRow_2].contigEnd;

			if(contigID_1<=0 || contigID_2<=0)
				continue;

			// Only paired reads been matched to different contigs are
			// considered, except very short contigs (e.g., lengths < 300 bp)
			//if(contigID_1==contigID_2 && endFlag_1!=2)
			if(contigID_1==contigID_2)
				continue;

			// check their contig ends to get the row and column
			if(endFlag_1==1)
			{
				row = 2 * contigID_1 - 1;
			}else
			{
				row = 2 * contigID_1 - 2;
			}

			if(endFlag_2==1)
			{
				col = 2 * contigID_2 - 1;
			}else
			{
				col = 2 * contigID_2 - 2;
			}

			// save the arrLoc item
			arrLocArr[itemNumArrLocArr].row = row;
			arrLocArr[itemNumArrLocArr].col = col;
			arrLocArr[itemNumArrLocArr+1].row = col;
			arrLocArr[itemNumArrLocArr+1].col = row;
			itemNumArrLocArr += 2;
		}
	}

	if(itemNumArrLocArr>maxItemNumArrLocArr)
	{
		printf("line=%d, In %s(), itemNumArrLocArr=%ld > maxItemNumArrLocArr=%ld\n", __LINE__, __func__, itemNumArrLocArr, maxItemNumArrLocArr);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Sort the arrLocArr by ascending according to row.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short sortArrLocArr(arrLoc *pArrLocArray, int64_t itemNumArrLocArray)
{
	struct partNode
	{
		int curItemNum;
		int totalItemNum;
		int firstRow;
	};

	int i, step, total;
	arrLoc *data, *buf, *pArrLocArrayBuf;
	struct partNode *part;
	int partArrSize, stepBits, maxStepLen;
	unsigned int bitMask;
	unsigned int hashcode, firstRow, curItemNum;

	stepBits = 16;
	maxStepLen = 32;
	partArrSize = 1 << stepBits;
	bitMask = (1 << stepBits) - 1;

	part = (struct partNode *) malloc(partArrSize * sizeof(struct partNode));
	if(part==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	pArrLocArrayBuf = (arrLoc *) malloc(itemNumArrLocArray * sizeof(arrLoc));
	if(pArrLocArrayBuf==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Begin to sort
	step = 0;
	while(step!=maxStepLen)
	{
		// set the data and buf
		if(step==stepBits)
		{
			buf = pArrLocArray;
			data = pArrLocArrayBuf;
		}else
		{
			data = pArrLocArray;
			buf = pArrLocArrayBuf;
		}

		// count the number
		if(memset(part, 0, partArrSize * sizeof(struct partNode))==NULL)
		{
			printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
			free(part);
			free(pArrLocArrayBuf);
			return FAILED;
		}
		for(i=0; i<itemNumArrLocArray; i++)
			part[ (data[i].row >> step) & bitMask ].totalItemNum ++;
			//part[ bitMask - ((data[i].row >> step) & bitMask) ].totalItemNum ++;

		// initialize the part array
		total = 0;
		for(i=0; i<partArrSize; i++)
		{
			part[i].firstRow = total;
			total += part[i].totalItemNum;
		}

		// copy the data to the right place
		for(i=0; i<itemNumArrLocArray; i++)
		{
			//hashcode = bitMask - ((data[i].row >> step) & bitMask);
			hashcode = (data[i].row >> step) & bitMask;
			firstRow = part[hashcode].firstRow;
			curItemNum = part[hashcode].curItemNum;
			if(memcpy(buf+firstRow+curItemNum, data+i, sizeof(arrLoc))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				free(part);
				free(pArrLocArrayBuf);
				return FAILED;
			}
			part[hashcode].curItemNum ++;
		}

		step += stepBits;

		//######################## Debug information #######################
#if DEBUG_FLAG
		for(i=0; i<partArrSize; i++)
		{
			if(part[i].curItemNum!=part[i].totalItemNum)
			{
				printf("line=%d, In %s(), in part[%d], curItemNum=%d != totalItemNum=%d, error!\n", __LINE__, __func__, i, part[i].curItemNum, part[i].totalItemNum);
				free(part);
				free(pArrLocArrayBuf);
				return FAILED;
			}
		}
#endif
		//######################## Debug information #######################
	}


	free(part);
	part = NULL;
	free(pArrLocArrayBuf);
	pArrLocArrayBuf = NULL;

	return SUCCESSFUL;
}

/**
 * Generate the contigGraph from arrLocArray.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateContigGraphEdge(ContigGraph *pContigGraphArray, int64_t itemNumContigGraphArray, arrLoc *pArrLocArray, int64_t itemNumArrLocArray)
{
	int64_t i, j, rowValue, itemNum, tmpLinkedContigsNum, tmpContigItemNum;
	int32_t *arrLocBucketArray;
	ContigEdge *pEdgeArray;

	arrLocBucketArray = (int32_t *) malloc (itemNumContigGraphArray * sizeof(int32_t));
	if(arrLocBucketArray==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	i = 0;
	while(i<itemNumArrLocArray)
	{
		if(memset(arrLocBucketArray, 0, itemNumContigGraphArray * sizeof(int32_t))==NULL)
		{
			printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// get the itemNum of a row
		if(getItemNumSingleRow(&itemNum, i, pArrLocArray, itemNumArrLocArray)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the itemNum of a row, error!\n", __LINE__, __func__);
			return FAILED;
		}

		for(j=0; j<itemNum; j++)
		{
			arrLocBucketArray[ pArrLocArray[i+j].col ] ++;
		}

		// get the arraySize
		rowValue = pArrLocArray[i].row; // get the rowValue
		tmpLinkedContigsNum = 0;
		for(j=0; j<itemNumContigGraphArray; j++)
		{
			if(arrLocBucketArray[j] > 0)
				tmpLinkedContigsNum ++;
		}
		pContigGraphArray[rowValue].arraySize = tmpLinkedContigsNum;

		// allocate the edge array
		pEdgeArray = (ContigEdge *) calloc (tmpLinkedContigsNum, sizeof(ContigEdge));
		if(pEdgeArray==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		pContigGraphArray[rowValue].pEdgeArray = pEdgeArray;

		// allocate contigGraph node
		tmpContigItemNum = 0;
		for(j=0; j<itemNumContigGraphArray; j++)
		{
			if(arrLocBucketArray[j] > 0)
			{
				pEdgeArray[tmpContigItemNum].col = j;
				pEdgeArray[tmpContigItemNum].used = NO;
				pEdgeArray[tmpContigItemNum].pairedNum = arrLocBucketArray[j];

				//printf("pEdgeArray[%ld]: (%ld, %ld, %d)\n", tmpContigItemNum, rowValue, j, arrLocBucketArray[j]);

				tmpContigItemNum ++;
			}
		}
		pContigGraphArray[rowValue].contigsNum = tmpContigItemNum;

		// ######################## Debug information #####################
#if DEBUG_FLAG
		if(pContigGraphArray[rowValue].contigsNum!=pContigGraphArray[rowValue].arraySize)
		{
			printf("line=%d, In %s(), contigsNum=%d != arraySize=%d, error!\n",  __LINE__, __func__, pContigGraphArray[rowValue].contigsNum, pContigGraphArray[rowValue].arraySize);
			return FAILED;
		}
#endif
		// ######################## Debug information #####################

		i += itemNum;
	}


	free(arrLocBucketArray);
	arrLocBucketArray = NULL;

	return SUCCESSFUL;
}

/**
 * Get the itemNum of a single row.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getItemNumSingleRow(int64_t *itemNum, int64_t row, arrLoc *pArrLocArray, int64_t itemNumArrLocArray)
{
	int64_t i;

	*itemNum = 1;
	for(i=row; i<itemNumArrLocArray-1; i++)
	{
		if(pArrLocArray[i].row==pArrLocArray[i+1].row)
		{
			(*itemNum) ++;
		}else
		{
			break;
		}
	}

	return SUCCESSFUL;
}

/**
 * Fill the item values of contigGraphEdge.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short FillContigGraphEdge(ContigGraph *pContigGraphArray, int64_t itemNumContigGraphArray, contigInfo *contigInfoArray, ContigList *contigListArray, int64_t contigItemNumInCLArray, ContigRead *contigReadArray, ReadList *sharedReadListArray, int64_t readItemNumInSRLArray, ReadPos *sharedReadPosArray)
{
	int64_t j, row, col, arraySize;
	ContigEdge *pEdgeArray;
	uint64_t contigID1, contigID2, endFlag1, endFlag2;


	for(row=0; row<itemNumContigGraphArray; row++)
	{
		arraySize = pContigGraphArray[row].arraySize;
		if(arraySize>0)
		{ // there are some linked contigs

			contigID1 = row / 2 + 1;
			if((row&1)==0)
			{ // 5' contig end
				if(contigInfoArray[contigID1-1].onlyEnd5==NO)
					endFlag1 = 0;
				else
					endFlag1 = 2;
			}else
			{ // 3' contig end
				endFlag1 = 1;
			}

			pEdgeArray = pContigGraphArray[row].pEdgeArray;
			for(j=0; j<arraySize; j++)
			{
				col = pEdgeArray[j].col;
				contigID2 = col / 2 + 1;
				if((col&1)==0)
				{ // 5' contig end
					if(contigInfoArray[contigID2-1].onlyEnd5==NO)
						endFlag2 = 0;
					else
						endFlag2 = 2;
				}else
				{ // 3' contig end
					endFlag2 = 1;
				}

				// ######################## Debug information ########################
				//if(contigID1==119 && contigID2==916)
				//{
				//	printf("contigID1=%lu, endFlag1=%lu, contigID2=%lu, endFlag2=%lu\n", contigID1, endFlag1, contigID2, endFlag2);
				//}
				// ######################## Debug information ########################

				// fill the validNum[4] of contigEdge
				if(fillSituationArray(pEdgeArray+j, contigID1, contigID2, endFlag1, endFlag2, contigInfoArray, contigListArray, contigReadArray, sharedReadListArray, readItemNumInSRLArray, sharedReadPosArray)==FAILED)
				{
					printf("line=%d, In %s(), cannot fill situation array, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// compute the gapSize of the two contigs


			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Remove invalid graph edges.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeInvalidGraphEdge(ContigGraph *pContigGraphArray, int64_t itemNumContigGraphArray)
{
	int64_t i, top, bot;
	ContigEdge *pEdgeArray;

	for(i=0; i<itemNumContigGraphArray; i++)
	{
		if(pContigGraphArray[i].arraySize>0)
		{
			pEdgeArray = pContigGraphArray[i].pEdgeArray;
			top = 0;
			bot = pContigGraphArray[i].arraySize - 1;
			while(top<=bot)
			{
				if(pEdgeArray[top].totalSituationNum==0)
				{
					if(top<bot)
					{
						if(memcpy(pEdgeArray+top, pEdgeArray+bot, sizeof(ContigEdge))==NULL)
						{
							printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
					bot --;
					pContigGraphArray[i].contigsNum --;
				}else
				{
					top ++;
				}
			}

			// update the arraySize of the graphEdge node
			pContigGraphArray[i].arraySize = pContigGraphArray[i].contigsNum;

			// free the invalid memory of invalid graphEdge nodes
			if(pContigGraphArray[i].arraySize==0)
			{
				free(pContigGraphArray[i].pEdgeArray);
				pContigGraphArray[i].pEdgeArray = NULL;
			}else if(pContigGraphArray[i].arraySize<0)
			{
				printf("line=%d, In %s(), pContigGraphArray[%ld].arraySize=%d, error!\n", __LINE__, __func__, i, pContigGraphArray[i].arraySize);
				return FAILED;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Fill contig information array.
 *  Note:
 *  	The element [0] corresponds to the contig 1.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillContigInfoArr(const char *contigFileName, contigInfo *pContigInfoArr)
{
	FILE *fpContig;
	char *contigSeq, *contigHeadTitle, *tmpContigSeq, *tmpContigHeadTitle;
	int64_t contigID, contigLen, maxContigLen, tmpContigNum, headTitleLen, returnFlag;

	fpContig = fopen(contigFileName, "r");
	if(fpContig==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigFileName);
		return FAILED;
	}

	contigHeadTitle = (char *) malloc(1000 * sizeof(char));
	if(contigHeadTitle==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(getMaxContigLenFromFile(&maxContigLen, contigFileName)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the maximal contig length, error!\n", __LINE__, __func__);
		return FAILED;
	}

	contigSeq = (char *) malloc ((maxContigLen+1)*sizeof(char));
	if(contigSeq==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//get contigs in fasta format
	contigID = 0;
	tmpContigNum = 0;
	while((returnFlag=getSingleContigFastaFromFile(fpContig, contigHeadTitle, contigSeq, &contigLen))==SUCCESSFUL)
	{
		// fill contig head title
		headTitleLen = strlen(contigHeadTitle);
		tmpContigHeadTitle = (char *) malloc ((headTitleLen+1) * sizeof(char));
		if(tmpContigHeadTitle==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(memcpy(tmpContigHeadTitle, contigHeadTitle, (headTitleLen+1) * sizeof(char))==NULL)
		{
			printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// fill short sequences
		tmpContigSeq = (char *) malloc ((contigLen+1) * sizeof(char));
		if(tmpContigSeq==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(memcpy(tmpContigSeq, contigSeq, (contigLen+1) * sizeof(char))==NULL)
		{
			printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		pContigInfoArr[tmpContigNum].contigID = ++contigID;
		pContigInfoArr[tmpContigNum].contigLen = contigLen;
		pContigInfoArr[tmpContigNum].headTitle = tmpContigHeadTitle;
		pContigInfoArr[tmpContigNum].contigSeq = tmpContigSeq;
		pContigInfoArr[tmpContigNum].used = 0;
		pContigInfoArr[tmpContigNum].usedTimeEnd5 = 0;
		pContigInfoArr[tmpContigNum].usedTimeEnd3 = 0;
		pContigInfoArr[tmpContigNum].overlapKindEnd5 = 0;
		pContigInfoArr[tmpContigNum].maxOverlapLenEnd5 = 0;
		pContigInfoArr[tmpContigNum].overlapKindEnd3 = 0;
		pContigInfoArr[tmpContigNum].maxOverlapLenEnd3 = 0;

		if(contigLen<=contigEndLen)
			pContigInfoArr[tmpContigNum].onlyEnd5 = YES;
		else
			pContigInfoArr[tmpContigNum].onlyEnd5 = NO;

		//if(contigLen<=1.5*contigEndLen)
		if(contigLen<=MAX_SHORT_LEN_THRES)		// added 2012-11-21
			pContigInfoArr[tmpContigNum].shortFlag = YES;
		else
			pContigInfoArr[tmpContigNum].shortFlag = NO;

		tmpContigNum ++;
	}

	free(contigSeq);
	contigSeq = NULL;

	free(contigHeadTitle);
	contigHeadTitle = NULL;

	fclose(fpContig);
	fpContig = NULL;

	// handle the error situation
	if(returnFlag==ERROR)
	{
		printf("line=%d, In %s(), cannot get the contig: %ld, error!\n", __LINE__, __func__, contigID);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Initialize the parameters for contigs linking.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initParaLinking(const char *averLinkNumFile)
{
	FILE *fpAverLinkNum;

	fpAverLinkNum = fopen(averLinkNumFile, "w");
	if(fpAverLinkNum==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, averLinkNumFile);
		return FAILED;
	}

	// get the average linked pairs per contigEdge
	if(computeAverPairsEachContigEdge(&averLinkNum, contigGraphArr, itemNumContigGraphArr, minLinksNumContigsThres)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the average linked pairs per contigEdge, error!\n", __LINE__, __func__);
		return FAILED;
	}

	fprintf(fpAverLinkNum, "%lf\n", averLinkNum);

	fclose(fpAverLinkNum);
	fpAverLinkNum = NULL;

	maxLinksNumContigsThres = averLinkNum * MAX_FIRST_LINKNUM_FACTOR;
	minLinksNumContigsThres = averLinkNum * MIN_FIRST_LINKNUM_FACTOR;
	maxRatioSecondFirstLinkNum = SECOND_FIRST_SCORE_RATIO_LINKING;
	secondLinkNumFactor = SECOND_LINKNUM_FACTOR;

	if(minLinksNumContigsThres>MIN_FIRST_LINKNUM_THRES)
	{
		minLinksNumContigsThres = MIN_FIRST_LINKNUM_THRES;
	}
	//else if(minLinksNumContigsThres<MIN_FIRST_LINKNUM_THRES)
	//{
	//	minLinksNumContigsThres = MIN_FIRST_LINKNUM_THRES;
	//}

	maxSecondLinkNumThres = minLinksNumContigsThres * secondLinkNumFactor;
	if(maxSecondLinkNumThres > MAX_SECOND_LINKNUM_THRES)
	{
		maxSecondLinkNumThres = MAX_SECOND_LINKNUM_THRES;
	}

#if (DEBUG_FLAG==YES)
	printf("maxLinksNumContigsThres=%.2f\n", maxLinksNumContigsThres);
	printf("minLinksNumContigsThres=%.2f\n", minLinksNumContigsThres);
	printf("maxSecondLinkNumThres=%.2f\n", maxSecondLinkNumThres);
	printf("maxRatioSecondFirstLinkNum=%.2f\n", maxRatioSecondFirstLinkNum);
#endif

	return SUCCESSFUL;
}

/**
 * Compute average pairs each contigEdge.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeAverPairsEachContigEdge(double *averageLinkNum, ContigGraph *contigGraphArray, int64_t itemNumContigGraphArray, int minLinksNumContigsThreshold)
{
	int64_t i, j, k;
	ContigEdge *pEdgeArray;
	int64_t totalPairs, totalEdges;

	totalPairs = 0;
	totalEdges = 0;
	for(i=0; i<itemNumContigGraphArray; i++)
	{
		if(contigGraphArray[i].arraySize>0)
		{
			pEdgeArray = contigGraphArray[i].pEdgeArray;
			for(j=0; j<contigGraphArray[i].arraySize; j++)
			{
				for(k=0; k<pEdgeArray[j].totalSituationNum; k++)
				{
					//if(pEdgeArray[j].validNums[ pEdgeArray[j].maxIndexes[k] ] >= minLinksNumContigsThreshold)
					if(pEdgeArray[j].validNums[ pEdgeArray[j].maxIndexes[k] ] >= 1)
					{
						totalPairs += pEdgeArray[j].validNums[ pEdgeArray[j].maxIndexes[k] ];
					}
				}

				totalEdges ++;
			}
		}
	}

	totalPairs /= 2;
	totalEdges /= 2;

	*averageLinkNum = (double)totalPairs / totalEdges;

#if (DEBUG_FLAG==YES)
	printf("totalPairs=%ld, totalEdges=%ld, averageLinkNum=%.2f\n", totalPairs, totalEdges, *averageLinkNum);
#endif

	return SUCCESSFUL;
}


/**
 * Link contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short linkContigs(const char *linkResultFile)
{
	FILE *fpLinkResult;
	int firstContigID;
	int linkID, linkRound, linkStatus, returnCode;
	maxRowCol *maxRowColNode;

	maxRowColNode = (maxRowCol*) malloc (sizeof(maxRowCol));
	if(maxRowColNode==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	fpLinkResult = fopen(linkResultFile, "w");
	if(fpLinkResult==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, linkResultFile);
		return FAILED;
	}

	// begin linking
	firstContigID = 1;
	linkID = 1;
	while(firstContigID<=contigsNum)
	{
		//############################ Debug information ######################
#if DEBUG_FLAG
		printf("============ Begin linking scaffolds: %d ============\n", linkID);
		if(linkID==9)
		{
			printf("$$$$$$$$$$$$$$$$$$$$$ linkID=%d!\n", linkID);
		}
#endif
		//############################ Debug information ######################

		linkRound = FIRST_LINK_ROUND;

		// get first contigID
		if(getFirstLinkedContigs(&firstContigID, maxRowColNode, contigInfoArr, contigsNum, contigGraphArr)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the first linked contigs in single scaffold, error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(maxRowColNode->contigID1==-1 || maxRowColNode->contigID2==-1)
			break;

		// link the first two contigs
		if(addContigToContigLinkArr(&linkStatus, contigLinkArr, &itemNumContigLinkArr, &headRowContigLinkArr, &tailRowContigLinkArr, contigInfoArr, maxRowColNode, 2, linkRound)==FAILED)
		{
			printf("line=%d, In %s(), cannot add contigs to contigLinkArr, error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(linkStatus==FAILED)
		{
			continue;
		}

		// reset the corresponding rows and columns
		if(markContigGraphEdge(contigGraphArr, maxRowColNode)==FAILED)
		{
			printf("line=%d, In %s(), cannot mark contigGraph node, error!\n", __LINE__, __func__);
			return FAILED;
		}
		contigInfoArr[maxRowColNode->contigID1-1].used = YES;
		contigInfoArr[maxRowColNode->contigID2-1].used = YES;

		// link contigs given a contig
		while(linkRound<=SECOND_LINK_ROUND)
		{
			// given a contig with end

			if(linkRound==FIRST_LINK_ROUND)
			{ // the first link round

				if(changeMaxRowCol(maxRowColNode, contigLinkArr, headRowContigLinkArr, contigInfoArr, linkRound, NO)==FAILED)
				{
					printf("line=%d, In %s(), cannot change maxRowCol, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// get the column having the maximal value
				if(getMaxColsOfSingleRow(contigGraphArr, maxRowColNode, contigInfoArr)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the maximal values for single row %d, error!\n", __LINE__, __func__, maxRowColNode->maxRow);
					return FAILED;
				}


				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || maxRowColNode->secondMaxValue>minLinksNumContigsThres*secondLinkNumFactor)))
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres)))
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || (maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3))))
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || maxRowColNode->secondMaxValue>0.5*averLinkNum))))  // added 2012-11-19
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres)))))  // added 2012-11-22
				//if(maxRowColNode->maxValue==0 || (maxRowColNode->maxValue<minLinksNumContigsThres && maxRowColNode->secondMaxValue>0) || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres)))))  // added 2012-11-22, deleted 2012-11-24
				if(maxRowColNode->maxValue==0 || maxRowColNode->maxValue>maxLinksNumContigsThres || (maxRowColNode->maxValue<minLinksNumContigsThres && maxRowColNode->secondMaxValue>0) || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres)))))  // added 2012-11-24
				{
					if(changeMaxRowCol(maxRowColNode, contigLinkArr, headRowContigLinkArr, contigInfoArr, linkRound, YES)==FAILED)
					{
						printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
						return FAILED;
					}

					linkRound ++;
					continue;
					//break;
				//}else if(contigInfoArr[maxRowColNode->contigID2-1].onlyEnd5==NO)
				}else if(contigInfoArr[maxRowColNode->contigID2-1].shortFlag==NO)
				{
					if(contigInfoArr[maxRowColNode->contigID2-1].used==NO)
					{
						returnCode = isLinkSingleton(maxRowColNode, contigGraphArr, linkRound);
						if(returnCode==NO)
						{
							if(changeMaxRowCol(maxRowColNode, contigLinkArr, headRowContigLinkArr, contigInfoArr, linkRound, YES)==FAILED)
							{
								printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
								return FAILED;
							}

							linkRound ++;
							continue;
						}else if(returnCode==ERROR)
						{
							printf("line=%d, In %s(), cannot check link singleton, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}else
					{
						if(changeMaxRowCol(maxRowColNode, contigLinkArr, headRowContigLinkArr, contigInfoArr, linkRound, YES)==FAILED)
						{
							printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
							return FAILED;
						}

						linkRound ++;
						continue;
					}
				}else
				{
					//if(maxRowColNode->secondMaxValue>0.1*averLinkNum && contigInfoArr[maxRowColNode->contigID1-1].onlyEnd5==YES)
					if(maxRowColNode->secondMaxValue>0.1*averLinkNum && contigInfoArr[maxRowColNode->contigID1-1].shortFlag==YES)
					{
						if(changeMaxRowCol(maxRowColNode, contigLinkArr, headRowContigLinkArr, contigInfoArr, linkRound, YES)==FAILED)
						{
							printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
							return FAILED;
						}

						linkRound ++;
						continue;
					}else
					{
						returnCode = isLinkSingleton(maxRowColNode, contigGraphArr, linkRound);
						if(returnCode==NO)
						{
							if(changeMaxRowCol(maxRowColNode, contigLinkArr, headRowContigLinkArr, contigInfoArr, linkRound, YES)==FAILED)
							{
								printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
								return FAILED;
							}

							linkRound ++;
							continue;
						}else if(returnCode==ERROR)
						{
							printf("line=%d, In %s(), cannot check link singleton, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}



//					if(maxRowColNode->maxValue<15)
//					{
//						if(changeMaxRowCol(maxRowColNode, contigLinkArr, headRowContigLinkArr, contigInfoArr, linkRound, YES)==FAILED)
//						{
//							printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
//							return FAILED;
//						}
//
//						linkRound ++;
//						continue;
//					}
				}

				// orient the contigs of the first round
				if(addContigToContigLinkArr(&linkStatus, contigLinkArr, &itemNumContigLinkArr, &headRowContigLinkArr, &tailRowContigLinkArr, contigInfoArr, maxRowColNode, 1, linkRound)==FAILED)
				{
					printf("line=%d, In %s(), cannot add contigs to contigLinkArr, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(linkStatus==FAILED)
				{
					if(changeMaxRowCol(maxRowColNode, contigLinkArr, headRowContigLinkArr, contigInfoArr, linkRound, YES)==FAILED)
					{
						printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
						return FAILED;
					}

					linkRound ++;
					continue;
				}

				// reset the corresponding rows and columns
				if(markContigGraphEdge(contigGraphArr, maxRowColNode)==FAILED)
				{
					printf("line=%d, In %s(), cannot mark contigGraph node, error!\n", __LINE__, __func__);
					return FAILED;
				}
				contigInfoArr[maxRowColNode->contigID2-1].used = YES;

			}else
			{ // the second link round

				if(changeMaxRowCol(maxRowColNode, contigLinkArr, headRowContigLinkArr, contigInfoArr, linkRound, NO)==FAILED)
				{
					printf("line=%d, In %s(), cannot change maxRowCol, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// get the column having the maximal value
				if(getMaxRowsOfSingleCol(contigGraphArr, maxRowColNode, contigInfoArr)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the maximal rows for column %u, error!\n", __LINE__, __func__, maxRowColNode->maxCol);
					return FAILED;
				}


				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || maxRowColNode->secondMaxValue>minLinksNumContigsThres*secondLinkNumFactor)))
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres)))
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || (maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3))))
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || maxRowColNode->secondMaxValue>0.5*averLinkNum)))) // added 2012-11-19
				//if(maxRowColNode->maxValue<minLinksNumContigsThres || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres))))) // added 2012-11-22
				//if(maxRowColNode->maxValue==0 || (maxRowColNode->maxValue<minLinksNumContigsThres && maxRowColNode->secondMaxValue>0) || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres))))) // added 2012-11-22, deleted 2012-11-24
				if(maxRowColNode->maxValue==0 || maxRowColNode->maxValue>maxLinksNumContigsThres || (maxRowColNode->maxValue<minLinksNumContigsThres && maxRowColNode->secondMaxValue>0) || (maxRowColNode->secondMaxValue>0 && ((double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>maxRatioSecondFirstLinkNum || ((maxRowColNode->secondMaxValue>maxSecondLinkNumThres && (double)maxRowColNode->secondMaxValue/maxRowColNode->maxValue>0.3) || (maxRowColNode->secondMaxValue>0.5*averLinkNum || maxRowColNode->secondMaxValue>maxSecondLinkNumThres))))) // added 2012-11-24
				{
					linkRound ++;
					break;
				//}else if(contigInfoArr[maxRowColNode->contigID1-1].onlyEnd5==NO)
				}else if(contigInfoArr[maxRowColNode->contigID1-1].shortFlag==NO)
				{
					if(contigInfoArr[maxRowColNode->contigID1-1].used==NO)
					{
						returnCode = isLinkSingleton(maxRowColNode, contigGraphArr, linkRound);
						if(returnCode==NO)
						{
							linkRound ++;
							break;
						}else if(returnCode==ERROR)
						{
							printf("line=%d, In %s(), cannot check link singleton, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}else
					{
						linkRound ++;
						break;
					}
				}else
				{
					//if(maxRowColNode->secondMaxValue>0.1*averLinkNum && contigInfoArr[maxRowColNode->contigID2-1].onlyEnd5==YES)
					if(maxRowColNode->secondMaxValue>0.1*averLinkNum && contigInfoArr[maxRowColNode->contigID2-1].shortFlag==YES)
					{
						linkRound ++;
						break;
					}else
					{
						returnCode = isLinkSingleton(maxRowColNode, contigGraphArr, linkRound);
						if(returnCode==NO)
						{
							linkRound ++;
							break;
						}else if(returnCode==ERROR)
						{
							printf("line=%d, In %s(), cannot check link singleton, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}


//					if(maxRowColNode->maxValue<15)
//					{
//						if(changeMaxRowCol(maxRowColNode, contigLinkArr, headRowContigLinkArr, contigInfoArr, linkRound, YES)==FAILED)
//						{
//							printf("line=%d, In %s(), cannot change the first round to the second round, error!\n", __LINE__, __func__);
//							return FAILED;
//						}
//
//						linkRound ++;
//						break;
//					}
				}

				// orient the contigs of the second round
				if(addContigToContigLinkArr(&linkStatus, contigLinkArr, &itemNumContigLinkArr, &headRowContigLinkArr, &tailRowContigLinkArr, contigInfoArr, maxRowColNode, 1, linkRound)==FAILED)
				{
					printf("line=%d, In %s(), cannot add contigs to contigLinkArr, error!\n", __LINE__, __func__);
					return FAILED;
				}
				if(linkStatus==FAILED)
				{
					linkRound ++;
					break;
				}

				// reset the corresponding rows and columns
				if(markContigGraphEdge(contigGraphArr, maxRowColNode)==FAILED)
				{
					printf("line=%d, In %s(), cannot mark contigGraph node, error!\n", __LINE__, __func__);
					return FAILED;
				}
				contigInfoArr[maxRowColNode->contigID1-1].used = YES;

			}
		} // while(1)

		// save the linked contigs
		saveLinkResultToFile(fpLinkResult, linkID);
		linkID ++;
	}

	saveUnlinkedContigsToFile(fpLinkResult, &linkID);

	free(maxRowColNode);
	maxRowColNode = NULL;

	printf("The scaffolds number is: %d\n", linkID-1);

	fclose(fpLinkResult);
	fpLinkResult = NULL;

	return SUCCESSFUL;
}


/**
 * Get the maximal columns and its value of a single row.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxColsOfSingleRow(ContigGraph *pContigGraphArray, maxRowCol *pMaxRowColNode, contigInfo *contigInfoArray)
{
	int i, j, situationNum;
	ContigGraph *pContigGraphNode;
	ContigEdge *pEdgeArray;

	pMaxRowColNode->maxValue = 0;
	pMaxRowColNode->secondMaxValue = 0;
	pMaxRowColNode->maxCol = -1;
	pMaxRowColNode->secondMaxCol = -1;

	pContigGraphNode = pContigGraphArray + pMaxRowColNode->maxRow;
	pEdgeArray = pContigGraphNode->pEdgeArray;
	for(i=0; i<pContigGraphNode->arraySize; i++)
	{
		if(pEdgeArray[i].used==NO)
		{
			situationNum = pEdgeArray[i].totalSituationNum;
			for(j=0; j<situationNum; j++)
			{
				if(pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ]> pMaxRowColNode->maxValue)
				{
					pMaxRowColNode->secondMaxValue = pMaxRowColNode->maxValue;
					pMaxRowColNode->secondMaxArrIndex = pMaxRowColNode->maxArrIndex;
					pMaxRowColNode->secondMaxCol = pMaxRowColNode->maxCol;
					pMaxRowColNode->maxValue = pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ];
					pMaxRowColNode->maxArrIndex = pEdgeArray[i].maxIndexes[j];
					pMaxRowColNode->maxCol = pEdgeArray[i].col;
				}else if(pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ]> pMaxRowColNode->secondMaxValue)
				{
					pMaxRowColNode->secondMaxValue = pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ];
					pMaxRowColNode->secondMaxArrIndex = pEdgeArray[i].maxIndexes[j];
					pMaxRowColNode->secondMaxCol = pEdgeArray[i].col;
				}

				// #################### Debug information #####################
#if DEBUG_FLAG
				if(pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ] > 0)
				{
					printf("line=%d, In %s(), contigs(%d, %d), row=%d, col=%d, value=%d, situationID=%d\n", __LINE__, __func__, pMaxRowColNode->maxRow/2+1, pEdgeArray[i].col/2+1, pMaxRowColNode->maxRow, pEdgeArray[i].col, pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ], pEdgeArray[i].maxIndexes[j]);
				}
#endif
				// #################### Debug information #####################
			}
		}
	}

	if(pMaxRowColNode->maxCol>=0)
	{
		pMaxRowColNode->contigID2 = pMaxRowColNode->maxCol / 2 + 1;
		pMaxRowColNode->endFlag2 = pMaxRowColNode->maxCol % 2;
		if(contigInfoArray[pMaxRowColNode->contigID2-1].onlyEnd5==YES)
			pMaxRowColNode->endFlag2 = 2;

		// #################### Debug information #####################
#if DEBUG_OUT_FLAG
		printf("line=%d, In %s(), contigID1=%d, contigID2=%d, endFlag1=%d, endFlag2=%d, maxValue=%d, secondMaxValue=%d\n", __LINE__, __func__, pMaxRowColNode->contigID1, pMaxRowColNode->contigID2, pMaxRowColNode->endFlag1, pMaxRowColNode->endFlag2, pMaxRowColNode->maxValue, pMaxRowColNode->secondMaxValue);
#endif
		// #################### Debug information #####################
	}else
	{
		pMaxRowColNode->contigID2 = -1;
		pMaxRowColNode->endFlag2 = -1;
	}

	return SUCCESSFUL;
}


/**
 * Get the maximal columns and its value of a single row.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxRowsOfSingleCol(ContigGraph *pContigGraphArray, maxRowCol *pMaxRowColNode, contigInfo *contigInfoArray)
{
	int i, j, situationNum, newSituationID;
	ContigGraph *pContigGraphNode;
	ContigEdge *pEdgeArray;

	pMaxRowColNode->maxValue = 0;
	pMaxRowColNode->secondMaxValue = 0;
	pMaxRowColNode->maxRow = -1;
	pMaxRowColNode->secondMaxRow = -1;

	pContigGraphNode = pContigGraphArray + pMaxRowColNode->maxCol;
	pEdgeArray = pContigGraphNode->pEdgeArray;
	for(i=0; i<pContigGraphNode->arraySize; i++)
	{
		if(pEdgeArray[i].used==NO)
		{
			situationNum = pEdgeArray[i].totalSituationNum;
			for(j=0; j<situationNum; j++)
			{
				switch(pEdgeArray[i].maxIndexes[j])
				{
					case 0: newSituationID = 2; break;
					case 1: newSituationID = 1; break;
					case 2: newSituationID = 0; break;
					case 3: newSituationID = 3; break;
					default: printf("line=%d, In %s(), maxIndexes[%d]=%d, error!\n", __LINE__, __func__, j, pEdgeArray[i].maxIndexes[j]); return FAILED;
				}

				if(pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ]> pMaxRowColNode->maxValue)
				{
					pMaxRowColNode->secondMaxValue = pMaxRowColNode->maxValue;
					pMaxRowColNode->secondMaxArrIndex = pMaxRowColNode->maxArrIndex;
					pMaxRowColNode->secondMaxRow = pMaxRowColNode->maxRow;
					pMaxRowColNode->maxValue = pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ];
					pMaxRowColNode->maxArrIndex = newSituationID;
					pMaxRowColNode->maxRow = pEdgeArray[i].col;
				}else if(pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ]> pMaxRowColNode->secondMaxValue)
				{
					pMaxRowColNode->secondMaxValue = pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ];
					pMaxRowColNode->secondMaxArrIndex = newSituationID;
					pMaxRowColNode->secondMaxRow = pEdgeArray[i].col;
				}

				// #################### Debug information #####################
#if DEBUG_FLAG
				if(pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ] > 0)
				{
					printf("line=%d, In %s(), contigs(%d, %d), row=%d, col=%d, value=%d, situationID=%d\n", __LINE__, __func__, pEdgeArray[i].col/2+1, pMaxRowColNode->maxCol/2+1, pEdgeArray[i].col, pMaxRowColNode->maxCol, pEdgeArray[i].validNums[ pEdgeArray[i].maxIndexes[j] ], newSituationID);
				}
#endif
				// #################### Debug information #####################

			}
		}
	}

	if(pMaxRowColNode->maxRow>=0)
	{
		pMaxRowColNode->contigID1 = pMaxRowColNode->maxRow / 2 + 1;
		pMaxRowColNode->endFlag1 = pMaxRowColNode->maxRow % 2;
		if(contigInfoArray[pMaxRowColNode->contigID1-1].onlyEnd5==YES)
			pMaxRowColNode->endFlag1 = 2;

		// #################### Debug information #####################
#if DEBUG_OUT_FLAG
		printf("line=%d, In %s(), contigID1=%d, contigID2=%d, endFlag1=%d, endFlag2=%d, maxValue=%d, secondMaxValue=%d\n", __LINE__, __func__, pMaxRowColNode->contigID1, pMaxRowColNode->contigID2, pMaxRowColNode->endFlag1, pMaxRowColNode->endFlag2, pMaxRowColNode->maxValue, pMaxRowColNode->secondMaxValue);
#endif
		// #################### Debug information #####################
	}else
	{
		pMaxRowColNode->contigID1 = -1;
		pMaxRowColNode->endFlag1 = -1;
	}

	return SUCCESSFUL;
}

/**
 * Change the maximal row and columns.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short changeMaxRowCol(maxRowCol *pMaxRowColNode, contigLink *contigLinkArray, int headRowContigLinkArray, contigInfo *contigInfoArray, int linkRound, int turnRoundFlag)
{
	if(turnRoundFlag==NO)
	{
		if(linkRound==1)
		{
			if(pMaxRowColNode->endFlag2!=2)
			{
				if(pMaxRowColNode->maxCol%2==0)
				{ // 5' end -> 3 'end
					pMaxRowColNode->maxRow = pMaxRowColNode->maxCol + 1;
					pMaxRowColNode->endFlag1 = 1;
				}else
				{ // 3' end -> 5 'end
					pMaxRowColNode->maxRow = pMaxRowColNode->maxCol - 1;
					pMaxRowColNode->endFlag1 = 0;
				}
				pMaxRowColNode->contigID1 = pMaxRowColNode->contigID2;

			}else
			{
				pMaxRowColNode->maxRow = pMaxRowColNode->maxCol;
				pMaxRowColNode->endFlag1 = 2;
				pMaxRowColNode->contigID1 = pMaxRowColNode->contigID2;
			}
		}else
		{
			if(pMaxRowColNode->endFlag1!=2)
			{
				if(pMaxRowColNode->maxRow%2==0)
				{ // 5' end -> 3 'end
					pMaxRowColNode->maxCol = pMaxRowColNode->maxRow + 1;
					pMaxRowColNode->endFlag2 = 1;
				}else
				{
					pMaxRowColNode->maxCol = pMaxRowColNode->maxRow - 1;
					pMaxRowColNode->endFlag2 = 0;
				}
				pMaxRowColNode->contigID2 = pMaxRowColNode->contigID1;

			}else
			{
				pMaxRowColNode->maxCol = pMaxRowColNode->maxRow;
				pMaxRowColNode->endFlag2 = 2;
				pMaxRowColNode->contigID2 = pMaxRowColNode->contigID1;
			}
		}
	}else
	{ // when turn the first round to the second round
		pMaxRowColNode->contigID1 = contigLinkArray[headRowContigLinkArray].contigID;
		if(contigInfoArray[pMaxRowColNode->contigID1-1].onlyEnd5==YES)
		{
			printf("line=%d, In %s(), the first contig in the linked list should not be short contigs, error!\n", __LINE__, __func__);
			return FAILED;

			//pMaxRowColNode->maxRow = (pMaxRowColNode->contigID1 - 1) * 2;
			//pMaxRowColNode->endFlag1 = 2;
		}else
		{
			if(contigLinkArray[headRowContigLinkArray].orientation==ORIENTATION_PLUS)
			{ // plus orientation
				pMaxRowColNode->maxRow = (pMaxRowColNode->contigID1 - 1) * 2 + 1;
				pMaxRowColNode->endFlag1 = 1;
			}else
			{ // minus orienataion
				pMaxRowColNode->maxRow = (pMaxRowColNode->contigID1 - 1) * 2;
				pMaxRowColNode->endFlag1 = 0;
			}
		}
	}

	return SUCCESSFUL;
}


/**
 * Get the maximal valid links from situtationArray.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxLinksSituationArray(double *pArray, int *maxRowID, int *secondRowID, double *maxValue, double *secondValue)
{
	int i;

	*secondValue = 0;
	*maxValue = 0;
	*maxRowID = -1;
	*secondRowID = -1;

	for(i=0; i<4; i++)
	{
		if(pArray[i]> *maxValue)
		{
			*secondValue = *maxValue;
			*maxValue = pArray[i];
			*secondRowID = *maxRowID;
			*maxRowID = i;
		}else if(pArray[i]> *secondValue)
		{
			*secondValue = pArray[i];
			*secondRowID = i;
		}
	}

	return SUCCESSFUL;
}

/**
 * Get first two linked contigs given the first contigID.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getFirstLinkedContigs(int *firstContigID, maxRowCol *pMaxRowColNode, contigInfo *contigInfoArray, int tmpContigsNum, ContigGraph *pContigGraphArray)
{
	int i, satisfiedFlag, returnCode;

	// get first contigID
	satisfiedFlag = NO;
	while((*firstContigID) <= tmpContigsNum)
	{
		//if(contigInfoArray[(*firstContigID)-1].used==0 && contigInfoArray[(*firstContigID)-1].onlyEnd5==NO)
		//if(contigInfoArray[(*firstContigID)-1].used==0 && contigInfoArray[(*firstContigID)-1].shortFlag==NO)
		if(contigInfoArray[(*firstContigID)-1].used==0 && contigInfoArray[(*firstContigID)-1].shortFlag==NO && contigInfoArray[(*firstContigID)-1].onlyEnd5==NO)
		{ // the first contig is unused and not too short

			pMaxRowColNode->maxRow = (*firstContigID) * 2 - 1;
			pMaxRowColNode->contigID1 = pMaxRowColNode->maxRow / 2 + 1;
			pMaxRowColNode->endFlag1 = pMaxRowColNode->maxRow % 2;
			if(contigInfoArray[pMaxRowColNode->contigID1-1].onlyEnd5==YES)
				pMaxRowColNode->endFlag1 = 2;

			for(i=0; i<2; i++)
			{
				if(getMaxColsOfSingleRow(pContigGraphArray, pMaxRowColNode, contigInfoArray)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the maximal values for single row %d, error!\n", __LINE__, __func__, pMaxRowColNode->maxRow);
					return FAILED;
				}

				//if(maxValue>=minLinksNumContigsThres || (secondMaxValue>0 && (secondMaxValue/maxValue<maxRatioSecondFirstLinkNum || secondMaxValue<minLinksNumContigsThres*secondLinkNumFactor)))
				//if(maxValue>=minLinksNumContigsThres && (secondMaxValue/maxValue<maxRatioSecondFirstLinkNum && secondMaxValue<minLinksNumContigsThres*secondLinkNumFactor))
				//if(pMaxRowColNode->maxValue>=minLinksNumContigsThres && ((double)pMaxRowColNode->secondMaxValue/pMaxRowColNode->maxValue<maxRatioSecondFirstLinkNum && pMaxRowColNode->secondMaxValue<minLinksNumContigsThres*secondLinkNumFactor))
				//if(pMaxRowColNode->maxValue>=minLinksNumContigsThres && ((double)pMaxRowColNode->secondMaxValue/pMaxRowColNode->maxValue<maxRatioSecondFirstLinkNum && pMaxRowColNode->secondMaxValue<maxSecondLinkNumThres))
				//if(pMaxRowColNode->maxValue>=minLinksNumContigsThres && ((double)pMaxRowColNode->secondMaxValue/pMaxRowColNode->maxValue<maxRatioSecondFirstLinkNum && pMaxRowColNode->secondMaxValue==0))  // deleted 2012-11-24
				if(pMaxRowColNode->maxValue>=minLinksNumContigsThres && pMaxRowColNode->maxValue<=maxLinksNumContigsThres && ((double)pMaxRowColNode->secondMaxValue/pMaxRowColNode->maxValue<maxRatioSecondFirstLinkNum && pMaxRowColNode->secondMaxValue==0))  // added 2012-11-24
				{
					//if(contigInfoArray[pMaxRowColNode->contigID2-1].onlyEnd5==NO)
					if(contigInfoArray[pMaxRowColNode->contigID2-1].shortFlag==NO)
					{
						returnCode = isLinkSingleton(pMaxRowColNode, pContigGraphArray, FIRST_LINK_ROUND);
						if(returnCode==YES)
						{
							if((contigInfoArray[pMaxRowColNode->contigID1-1].contigLen<5*contigEndLen || contigInfoArray[pMaxRowColNode->contigID2-1].contigLen<5*contigEndLen) && pMaxRowColNode->maxValue>2*averLinkNum)
								satisfiedFlag = NO;
							else
								satisfiedFlag = YES;
						}else if(returnCode==NO)
							satisfiedFlag = NO;
						else
						{
							printf("line=%d, In %s(), cannot check link singleton, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}else
					{
						//satisfiedFlag = YES;
						satisfiedFlag = NO;
					}
				}

				if(satisfiedFlag==YES)
					break;

				pMaxRowColNode->maxRow --;
			}
		}

		(*firstContigID) ++;

		if(satisfiedFlag==YES)
			break;
	}

	if(satisfiedFlag==NO)
	{
		pMaxRowColNode->contigID1 = pMaxRowColNode->contigID2 = -1;
		pMaxRowColNode->endFlag1 = pMaxRowColNode->endFlag2 = -1;
	}

	return SUCCESSFUL;
}

/**
 * Check whether the link is singleton or not.
 *  @return:
 *   (1) If singleton, return YES; otherwise, return NO.
 *   (2) If errors occurred, return ERROR.
 */
short isLinkSingleton(maxRowCol *pMaxRowColNode, ContigGraph *pContigGraphArray, int linkRound)
{
	int64_t j, maxRow, maxCol;
	ContigGraph *pContigGraphNode;
	ContigEdge *pEdgeArray;
	int singletonFlag;

	singletonFlag = YES;

	if(linkRound==1)
	{ // the first link round
		maxCol = pMaxRowColNode->maxCol;
		pContigGraphNode = pContigGraphArray + maxCol;
		if(pContigGraphNode->arraySize==1)
		{
			singletonFlag = YES;
		}else if(pContigGraphNode->arraySize>1)
		{
			maxRow = pMaxRowColNode->maxRow;
			pEdgeArray = pContigGraphNode->pEdgeArray;
			for(j=0; j<pContigGraphNode->arraySize; j++)
			{

#if (DEBUG_FLAG==YES)
				//printf("num=%d\n", pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]]);
				printf("line=%d, In %s(), contigs(%d, %d), row=%d, col=%d, value=%d, situationID=%d\n", __LINE__, __func__, pMaxRowColNode->maxCol/2+1, pEdgeArray[j].col/2+1, pMaxRowColNode->maxCol, pEdgeArray[j].col, pEdgeArray[j].validNums[ pEdgeArray[j].maxIndexes[0] ], pEdgeArray[j].maxIndexes[0]);
#endif

				//if(pEdgeArray[j].col != maxRow && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] >= minLinksNumContigsThres*secondLinkNumFactor)
				//if(pEdgeArray[j].col != maxRow && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] >= maxSecondLinkNumThres)
				//if(pEdgeArray[j].col != maxRow && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] >= minLinksNumContigsThres)
				//if(singletonFlag==YES && pEdgeArray[j].col != maxRow && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] > 0)  // added 2012-11-21
				if(singletonFlag==YES && pEdgeArray[j].col != maxRow && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] > maxSecondLinkNumThres)  // added 2012-11-21
				{
					singletonFlag = NO;
				}
			}
		}else
		{
			printf("line=%d, In %s(), linkRound=%d, pContigGraphArray[%ld].arraySize=%d, error!\n", __LINE__, __func__, linkRound, maxCol, pContigGraphNode->arraySize);
			return ERROR;
		}
	}else
	{ // the second link round
		maxRow = pMaxRowColNode->maxRow;
		pContigGraphNode = pContigGraphArray + maxRow;
		if(pContigGraphNode->arraySize==1)
		{
			singletonFlag = YES;
		}else if(pContigGraphNode->arraySize>1)
		{
			maxCol = pMaxRowColNode->maxCol;
			pEdgeArray = pContigGraphNode->pEdgeArray;
			for(j=0; j<pContigGraphNode->arraySize; j++)
			{

#if (DEBUG_FLAG==YES)
				//printf("num=%d\n", pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]]);
				printf("line=%d, In %s(), contigs(%d, %d), row=%d, col=%d, value=%d, situationID=%d\n", __LINE__, __func__, pMaxRowColNode->maxRow/2+1, pEdgeArray[j].col/2+1, pMaxRowColNode->maxRow, pEdgeArray[j].col, pEdgeArray[j].validNums[ pEdgeArray[j].maxIndexes[0] ], pEdgeArray[j].maxIndexes[0]);
#endif

				//if(pEdgeArray[j].col != maxCol && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] >= minLinksNumContigsThres*secondLinkNumFactor)
				//if(pEdgeArray[j].col != maxCol && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] >= maxSecondLinkNumThres)
				//if(pEdgeArray[j].col != maxCol && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] >= minLinksNumContigsThres)
				//if(singletonFlag==YES && pEdgeArray[j].col != maxCol && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] > 0)  // added 2012-11-21
				if(singletonFlag==YES && pEdgeArray[j].col != maxCol && pEdgeArray[j].validNums[pEdgeArray[j].maxIndexes[0]] > maxSecondLinkNumThres)  // added 2012-11-21
				{
					singletonFlag = NO;
				}
			}
		}else
		{
			printf("line=%d, In %s(), linkRound=%d, pContigGraphArray[%ld].arraySize=%d, error!\n", __LINE__, __func__, linkRound, maxRow, pContigGraphNode->arraySize);
			return ERROR;
		}
	}

	return singletonFlag;
}

/**
 * Fill the situationArray.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillSituationArray(ContigEdge *pEdgeNode, int contigID1, int contigID2, int endFlag1, int endFlag2, contigInfo *contigInfoArray, ContigList *contigListArray, ContigRead *contigReadArray, ReadList *sharedReadListArray, int64_t itemNumInSRL, ReadPos *sharedReadPosArray)
{
	uint64_t i, j, readID1, readID2;
	ContigList *pContigList2;
	ContigRead *pContigRead2;
	int64_t hitRow;
	int rowNum2, orientation1, orientation2, endRead1, endRead2;
	ReadList *pReadList1, *pReadList2;
	ReadPos *pReadPos1;
	int tmp_contigID1, tmp_contigPos1;
	int32_t countArray[4];

	for(i=0; i<4; i++) pEdgeNode->validNums[i] = 0;  // reset the situation array

	// check each paired read between the two contig ends
	pContigList2 = contigListArray + contigID2 - 1;
	if(endFlag2==1)
	{
		pContigRead2 = contigReadArray + pContigList2->firstRow3;
		rowNum2 = pContigList2->EndNum3;
	}else
	{
		pContigRead2 = contigReadArray + pContigList2->firstRow5;
		rowNum2 = pContigList2->EndNum5;
	}

	for(i=0; i<rowNum2; i++)
	{
		if(pContigRead2->readID>0)
		{
			readID2 = pContigRead2->readID;
			orientation2 = pContigRead2->orientation;
			endRead2 = pContigRead2->contigEnd;

			if(endRead2==endFlag2)
			{
				hitRow = getReadRowFromReadList(readID2, sharedReadListArray, itemNumInSRL);
				if(hitRow>=0)
				{
					pReadList2 = sharedReadListArray + hitRow;
					if(pReadList2->curNum==1)
					{
						// get its paired end read
						if(readID2%2==1)
						{ // odd number
							readID1 = readID2 + 1;
						}else
						{ // even number
							readID1 = readID2 - 1;
						}

						hitRow = getReadRowFromReadList(readID1, sharedReadListArray, itemNumInSRL);
						if(hitRow>=0)
						{
							pReadList1 = sharedReadListArray + hitRow;
							if(pReadList1->curNum==1)
							{
								pReadPos1 = sharedReadPosArray + pReadList1->firstRow;

								orientation1 = pReadPos1->orientation;
								endRead1 = pReadPos1->contigEnd;

								tmp_contigPos1 = pReadPos1->contigPos;

								tmp_contigID1 = pReadPos1->contigID;
								if(tmp_contigID1==contigID1 && endRead1==endFlag1 && contigInfoArray[contigID1-1].used==0)
								{
									// fill the situations array
									if(endFlag1!=0 && orientation1==ORIENTATION_PLUS && endFlag2!=1 && orientation2==ORIENTATION_MINUS)
									{ // (3', +) , (5', -) ==> (+, +)
										pEdgeNode->validNums[0] ++;
									}else if(endFlag1!=0 && orientation1==ORIENTATION_PLUS && endFlag2!=0 && orientation2==ORIENTATION_PLUS)
									{ // (3', +) , (3', +) ==> (+, -)
										pEdgeNode->validNums[1] ++;
									}else if(endFlag1!=1 && orientation1==ORIENTATION_MINUS && endFlag2!=0 && orientation2==ORIENTATION_PLUS)
									{ // (5', -) , (3', +) ==> (-, -)
										pEdgeNode->validNums[2] ++;
									}else if(endFlag1!=1 && orientation1==ORIENTATION_MINUS && endFlag2!=1 && orientation2==ORIENTATION_MINUS)
									{ // (5', -) , (5', -) ==> (-, +)
										pEdgeNode->validNums[3] ++;
									}else
									{
#if DEBUG_OUT_FLAG
										printf("endRead1=%d, endRead2=%d, readID1=%lu, readID2=%lu, (endFlag1=%d, orientation1=%d, endFlag2=%d, orientation2=%d)\n", endRead1, endRead2, readID1, readID2, endFlag1, orientation1, endFlag2, orientation2);
#endif
									}
								}
							}
						}
					}
				}
			}
		}

		pContigRead2 ++;
	} // end for(i=0; i<rowNum2; i++)

	// compute the maxIndexes[4] by counting sort method
	for(i=0; i<4; i++) countArray[i] = 0;
	for(i=0; i<3; i++)
	{
		for(j=i+1; j<4; j++)
		{
			if(pEdgeNode->validNums[i] >= pEdgeNode->validNums[j])
				countArray[j] ++;
			else if(pEdgeNode->validNums[i] < pEdgeNode->validNums[j])
				countArray[i] ++;
		}
	}

	for(i=0; i<4; i++)
		pEdgeNode->maxIndexes[ countArray[i] ] = i;

	for(i=0; i<4; i++)
	{
		if(pEdgeNode->validNums[ pEdgeNode->maxIndexes[i] ] > 0)
		{
			pEdgeNode->totalSituationNum ++;
		}else
		{
			pEdgeNode->maxIndexes[i] = -1;
		}
	}

	return SUCCESSFUL;
}

/**
 * Add linked contig to contigLinkArray.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addContigToContigLinkArr(int *linkStatus, contigLink *contigLinkArray, int *itemNumContigLinkArray, int *headRowContigLinkArray, int *tailRowContigLinkArray, contigInfo *contigInfoArray, maxRowCol *pMaxRowColNode, int newContigNum, int linkRound)
{
	if(linkRound==1)
	{ // the first link round
		if(newContigNum==2)
		{ // orient the first two contigs

			*itemNumContigLinkArray = 0;
			*headRowContigLinkArray = -1;
			*tailRowContigLinkArray = -1;

			if(pMaxRowColNode->maxArrIndex==0)
			{
				contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID1;
				contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_PLUS;
				contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID1-1].contigLen;
				contigLinkArray[*itemNumContigLinkArray].previous = -1;
				contigLinkArray[*itemNumContigLinkArray].next = -1;
				*headRowContigLinkArray = *tailRowContigLinkArray = *itemNumContigLinkArray;
				(*itemNumContigLinkArray) ++;

				contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID2;
				contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_PLUS;
				contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID2-1].contigLen;
				contigLinkArray[*itemNumContigLinkArray].previous = -1;
				contigLinkArray[*itemNumContigLinkArray].next = -1;
				*tailRowContigLinkArray = *itemNumContigLinkArray;
				(*itemNumContigLinkArray) ++;

				contigLinkArray[*headRowContigLinkArray].next = *tailRowContigLinkArray;
				contigLinkArray[*tailRowContigLinkArray].previous = *headRowContigLinkArray;

				contigInfoArray[pMaxRowColNode->contigID1-1].used = 1;
				contigInfoArray[pMaxRowColNode->contigID2-1].used = 1;

			}else if(pMaxRowColNode->maxArrIndex==1)
			{
				contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID1;
				contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_PLUS;
				contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID1-1].contigLen;
				contigLinkArray[*itemNumContigLinkArray].previous = -1;
				contigLinkArray[*itemNumContigLinkArray].next = -1;
				*headRowContigLinkArray = *tailRowContigLinkArray = *itemNumContigLinkArray;
				(*itemNumContigLinkArray) ++;

				contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID2;
				contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_MINUS;
				contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID2-1].contigLen;
				contigLinkArray[*itemNumContigLinkArray].previous = -1;
				contigLinkArray[*itemNumContigLinkArray].next = -1;
				*tailRowContigLinkArray = *itemNumContigLinkArray;
				(*itemNumContigLinkArray) ++;

				contigLinkArray[*headRowContigLinkArray].next = *tailRowContigLinkArray;
				contigLinkArray[*tailRowContigLinkArray].previous = *headRowContigLinkArray;

				contigInfoArray[pMaxRowColNode->contigID1-1].used = 1;
				contigInfoArray[pMaxRowColNode->contigID2-1].used = 1;

			}else if(pMaxRowColNode->maxArrIndex==2)
			{
				contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID1;
				contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_MINUS;
				contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID1-1].contigLen;
				contigLinkArray[*itemNumContigLinkArray].previous = -1;
				contigLinkArray[*itemNumContigLinkArray].next = -1;
				*headRowContigLinkArray = *tailRowContigLinkArray = *itemNumContigLinkArray;
				(*itemNumContigLinkArray) ++;

				contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID2;
				contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_MINUS;
				contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID2-1].contigLen;
				contigLinkArray[*itemNumContigLinkArray].previous = -1;
				contigLinkArray[*itemNumContigLinkArray].next = -1;
				*tailRowContigLinkArray = *itemNumContigLinkArray;
				(*itemNumContigLinkArray) ++;

				contigLinkArray[*headRowContigLinkArray].next = *tailRowContigLinkArray;
				contigLinkArray[*tailRowContigLinkArray].previous = *headRowContigLinkArray;

				contigInfoArray[pMaxRowColNode->contigID1-1].used = 1;
				contigInfoArray[pMaxRowColNode->contigID2-1].used = 1;

			}else if(pMaxRowColNode->maxArrIndex==3)
			{
				contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID1;
				contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_MINUS;
				contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID1-1].contigLen;
				contigLinkArray[*itemNumContigLinkArray].previous = -1;
				contigLinkArray[*itemNumContigLinkArray].next = -1;
				*headRowContigLinkArray = *tailRowContigLinkArray = *itemNumContigLinkArray;
				(*itemNumContigLinkArray) ++;

				contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID2;
				contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_PLUS;
				contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID2-1].contigLen;
				contigLinkArray[*itemNumContigLinkArray].previous = -1;
				contigLinkArray[*itemNumContigLinkArray].next = -1;
				*tailRowContigLinkArray = *itemNumContigLinkArray;
				(*itemNumContigLinkArray) ++;

				contigLinkArray[*headRowContigLinkArray].next = *tailRowContigLinkArray;
				contigLinkArray[*tailRowContigLinkArray].previous = *headRowContigLinkArray;

				contigInfoArray[pMaxRowColNode->contigID1-1].used = 1;
				contigInfoArray[pMaxRowColNode->contigID2-1].used = 1;

			}else
			{
				printf("line=%d, In %s(), linkRound=%d, situationID=%d\n", __LINE__, __func__, linkRound, pMaxRowColNode->maxArrIndex);
				return FAILED;
			}

			*linkStatus = SUCCESSFUL;

		}else
		{ // link the second contig
			if(pMaxRowColNode->maxArrIndex==0)
			{
				if(contigLinkArray[*tailRowContigLinkArray].orientation==ORIENTATION_PLUS)
				{ // the orientation of linked contig is changed
					contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID2;
					contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_PLUS;
					contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID2-1].contigLen;
					contigLinkArray[*itemNumContigLinkArray].previous = *tailRowContigLinkArray;
					contigLinkArray[*itemNumContigLinkArray].next = -1;

					contigLinkArray[*tailRowContigLinkArray].next = *itemNumContigLinkArray;
					*tailRowContigLinkArray = *itemNumContigLinkArray;
					(*itemNumContigLinkArray) ++;

					contigInfoArray[pMaxRowColNode->contigID2-1].used = 1;

					*linkStatus = SUCCESSFUL;

				}else
				{
					*linkStatus = FAILED;
				}

			}else if(pMaxRowColNode->maxArrIndex==1)
			{
				if(contigLinkArray[*tailRowContigLinkArray].orientation==ORIENTATION_PLUS)
				{ // the orientation of linked contig is changed
					contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID2;
					contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_MINUS;
					contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID2-1].contigLen;
					contigLinkArray[*itemNumContigLinkArray].previous = *tailRowContigLinkArray;
					contigLinkArray[*itemNumContigLinkArray].next = -1;

					contigLinkArray[*tailRowContigLinkArray].next = *itemNumContigLinkArray;
					*tailRowContigLinkArray = *itemNumContigLinkArray;
					(*itemNumContigLinkArray) ++;

					contigInfoArray[pMaxRowColNode->contigID2-1].used = 1;

					*linkStatus = SUCCESSFUL;

				}else
				{
					*linkStatus = FAILED;
				}

			}else if(pMaxRowColNode->maxArrIndex==2)
			{
				if(contigLinkArray[*tailRowContigLinkArray].orientation==ORIENTATION_MINUS)
				{ // the orientation of linked contig is changed
					contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID2;
					contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_MINUS;
					contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID2-1].contigLen;
					contigLinkArray[*itemNumContigLinkArray].previous = *tailRowContigLinkArray;
					contigLinkArray[*itemNumContigLinkArray].next = -1;

					contigLinkArray[*tailRowContigLinkArray].next = *itemNumContigLinkArray;
					*tailRowContigLinkArray = *itemNumContigLinkArray;
					(*itemNumContigLinkArray) ++;

					contigInfoArray[pMaxRowColNode->contigID2-1].used = 1;

					*linkStatus = SUCCESSFUL;

				}else
				{
					*linkStatus = FAILED;
				}

			}else if(pMaxRowColNode->maxArrIndex==3)
			{
				if(contigLinkArray[*tailRowContigLinkArray].orientation==ORIENTATION_MINUS)
				{ // the orientation of linked contig is changed
					contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID2;
					contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_PLUS;
					contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID2-1].contigLen;
					contigLinkArray[*itemNumContigLinkArray].previous = *tailRowContigLinkArray;
					contigLinkArray[*itemNumContigLinkArray].next = -1;

					contigLinkArray[*tailRowContigLinkArray].next = *itemNumContigLinkArray;
					*tailRowContigLinkArray = *itemNumContigLinkArray;
					(*itemNumContigLinkArray) ++;

					contigInfoArray[pMaxRowColNode->contigID2-1].used = 1;

					*linkStatus = SUCCESSFUL;

				}else
				{
					*linkStatus = FAILED;
				}

			}else
			{
				printf("line=%d, In %s(), linkRound=%d, situationID=%d\n", __LINE__, __func__, linkRound, pMaxRowColNode->maxArrIndex);
				return FAILED;
			}
		}
	}else
	{ // the second link round

		if(pMaxRowColNode->maxArrIndex==0)
		{
			if(contigLinkArray[*headRowContigLinkArray].orientation==ORIENTATION_PLUS)
			{ // the orientation of linked contig is changed
				contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID1;
				contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_PLUS;
				contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID1-1].contigLen;
				contigLinkArray[*itemNumContigLinkArray].previous = -1;
				contigLinkArray[*itemNumContigLinkArray].next = *headRowContigLinkArray;

				contigLinkArray[*headRowContigLinkArray].previous = *itemNumContigLinkArray;
				*headRowContigLinkArray = *itemNumContigLinkArray;
				(*itemNumContigLinkArray) ++;

				contigInfoArray[pMaxRowColNode->contigID1-1].used = 1;

				*linkStatus = SUCCESSFUL;

			}else
			{
				*linkStatus = FAILED;
			}

		}else if(pMaxRowColNode->maxArrIndex==1)
		{
			if(contigLinkArray[*headRowContigLinkArray].orientation==ORIENTATION_MINUS)
			{ // the orientation of linked contig is changed
				contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID1;
				contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_PLUS;
				contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID1-1].contigLen;
				contigLinkArray[*itemNumContigLinkArray].previous = -1;
				contigLinkArray[*itemNumContigLinkArray].next = *headRowContigLinkArray;

				contigLinkArray[*headRowContigLinkArray].previous = *itemNumContigLinkArray;
				*headRowContigLinkArray = *itemNumContigLinkArray;
				(*itemNumContigLinkArray) ++;

				contigInfoArray[pMaxRowColNode->contigID1-1].used = 1;

				*linkStatus = SUCCESSFUL;

			}else
			{
				*linkStatus = FAILED;
			}

		}else if(pMaxRowColNode->maxArrIndex==2)
		{
			if(contigLinkArray[*headRowContigLinkArray].orientation==ORIENTATION_MINUS)
			{ // the orientation of linked contig is changed
				contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID1;
				contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_MINUS;
				contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID1-1].contigLen;
				contigLinkArray[*itemNumContigLinkArray].previous = -1;
				contigLinkArray[*itemNumContigLinkArray].next = *headRowContigLinkArray;

				contigLinkArray[*headRowContigLinkArray].previous = *itemNumContigLinkArray;
				*headRowContigLinkArray = *itemNumContigLinkArray;
				(*itemNumContigLinkArray) ++;

				contigInfoArray[pMaxRowColNode->contigID1-1].used = 1;

				*linkStatus = SUCCESSFUL;

			}else
			{
				*linkStatus = FAILED;
			}

		}else if(pMaxRowColNode->maxArrIndex==3)
		{
			if(contigLinkArray[*headRowContigLinkArray].orientation==ORIENTATION_PLUS)
			{ // the orientation of linked contig is changed
				contigLinkArray[*itemNumContigLinkArray].contigID = pMaxRowColNode->contigID1;
				contigLinkArray[*itemNumContigLinkArray].orientation = ORIENTATION_MINUS;
				contigLinkArray[*itemNumContigLinkArray].contigLen = contigInfoArray[pMaxRowColNode->contigID1-1].contigLen;
				contigLinkArray[*itemNumContigLinkArray].previous = -1;
				contigLinkArray[*itemNumContigLinkArray].next = *headRowContigLinkArray;

				contigLinkArray[*headRowContigLinkArray].previous = *itemNumContigLinkArray;
				*headRowContigLinkArray = *itemNumContigLinkArray;
				(*itemNumContigLinkArray) ++;

				contigInfoArray[pMaxRowColNode->contigID1-1].used = 1;

				*linkStatus = SUCCESSFUL;

			}else
			{
				*linkStatus = FAILED;
			}

		}else
		{
			printf("line=%d, In %s(), linkRound=%d, situationID=%d\n", __LINE__, __func__, linkRound, pMaxRowColNode->maxArrIndex);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Mark contigGraph node.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short markContigGraphEdge(ContigGraph *pContigGraphArray, maxRowCol *pMaxRowColNode)
{
	int i;
	ContigGraph *pGraphNode;
	ContigEdge *pEdgeArray;

	pGraphNode = pContigGraphArray + pMaxRowColNode->maxRow;
	pEdgeArray = pGraphNode->pEdgeArray;
	for(i=0; i<pGraphNode->arraySize; i++)
	{
		if(pEdgeArray[i].col==pMaxRowColNode->maxCol && pEdgeArray[i].used==NO)
		{
			pEdgeArray[i].used = YES;
		}
	}

	pGraphNode = pContigGraphArray + pMaxRowColNode->maxCol;
	pEdgeArray = pGraphNode->pEdgeArray;
	for(i=0; i<pGraphNode->arraySize; i++)
	{
		if(pEdgeArray[i].col==pMaxRowColNode->maxRow && pEdgeArray[i].used==NO)
		{
			pEdgeArray[i].used = YES;
		}
	}

	return SUCCESSFUL;
}

/**
 * Save the linked result to file.
 */
void saveLinkResultToFile(FILE *fpLinkResult, int linkID)
{
	int tmpRow;
	tmpRow = headRowContigLinkArr;

	fprintf(fpLinkResult, ">%d\t%d\n", linkID, itemNumContigLinkArr);
	while(1)
	{
		fprintf(fpLinkResult, "%d\t%d\t%d\n", contigLinkArr[tmpRow].contigID, contigLinkArr[tmpRow].orientation, contigLinkArr[tmpRow].contigLen);

		if(tmpRow==tailRowContigLinkArr)
			break;

		tmpRow = contigLinkArr[tmpRow].next;
	}

	//#######################################
	//fflush(fpLinkResult);
	//#######################################
}

/**
 * Save unlinked contigs to file.
 */
void saveUnlinkedContigsToFile(FILE *fpLinkResult, int *startLinkID)
{
	int i, tmpRow;
	tmpRow = headRowContigLinkArr;

	for(i=0; i<contigsNum; i++)
	{
		if(contigInfoArr[i].used==0)
		{
			fprintf(fpLinkResult, ">%d\t%d\n", (*startLinkID) ++, 1);
			fprintf(fpLinkResult, "%d\t%d\t%d\n", i+1, ORIENTATION_PLUS, contigInfoArr[i].contigLen);
		}
	}
	//#######################################
	//fflush(fpLinkResult);
	//#######################################
}


/**
 * Get the total number of scaffolds.
 *  @return:
 *   if succeed, return the number; otherwise, return ERROR.
 */
int getScaffoldsNum(const char *linkResultFile)
{
	FILE *fpLinkResult;
	int totalNum;
	char ch;

	fpLinkResult = fopen(linkResultFile, "r");
	if(fpLinkResult==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, linkResultFile);
		return ERROR;
	}

	// compute the total number
	totalNum = 0;
	ch = fgetc(fpLinkResult);
	while(ch!=EOF)
	{
		if(ch=='>')
			totalNum ++;

		ch = fgetc(fpLinkResult);
	}

	fclose(fpLinkResult);

	return totalNum;
}
