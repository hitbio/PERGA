/*
 * scafPerga.c
 *
 *  Created on: Dec 17, 2013
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Start the scaffolding.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short startScaffolding(char *scafSeqFile, char *readMatchInfoFile, contigGraph_t *contigGraph, readSet_t *readSet)
{
	printf("\n============= Begin scaffolding, please wait ... =============\n");

	struct timeval tp_start,tp_end;
	double time_used;
	gettimeofday(&tp_start,NULL);


	// initialize the global paras
	if(initMemScaffolding(&scaffoldSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory for scaffolding, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// build the scafContigIndex
	if(buildScafContigIndex(&scafContigIndex, contigGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot build contig index for scaffolding, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// map reads: map all the reads onto contig ends
	if(mapReads(contigGraph, readSet, scafContigIndex)==FAILED)
	{
		printf("line=%d, In %s(), cannot map reads to contigs for scaffolding, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the reads match information to contig ends
	if(fillReadMatchInfoContigEnds(contigGraph, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill read match information to read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ############################## Debug information ##########################
//	if(outputContigReadArrayInScaf(contigGraph)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot output the contig reads, error!\n", __LINE__, __func__);
//		return FAILED;
//	}
	// ############################## Debug information ##########################


	// construct scaffold graph and link contigs
	if(contigsLinking(scaffoldSet, contigGraph, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot link contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ############################## Debug information ##########################
//	if(outputScaffoldSetToFile("../tmpScafLinked_merge2.txt", scaffoldSet, contigGraph)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot output linked contigs, error!\n", __LINE__, __func__);
//		return FAILED;
//	}
	// ############################## Debug information ##########################

	// compute the overlaps
	if(overlapContigsInScaf(scaffoldSet, contigGraph, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot overlap contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ############################## Debug information ##########################
//	if(outputScaffoldSetToFile("../tmpScafOverlapped_merge2.txt", scaffoldSet, contigGraph)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot output linked contigs, error!\n", __LINE__, __func__);
//		return FAILED;
//	}
	// ############################## Debug information ##########################


	// gap filling
	if(gapFillFlag==YES)
	{
		if(gapFilling(scaffoldSet, contigGraph, readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot fill gaps, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// ############################## Debug information ##########################
//	if(outputScaffoldSetToFile("../tmpScafFilled_merge2.txt", scaffoldSet, contigGraph)==FAILED)
//	{
//		printf("line=%d, In %s(), cannot output linked contigs, error!\n", __LINE__, __func__);
//		return FAILED;
//	}
	// ############################## Debug information ##########################

	// generate the scaffold sequences
	if(generateScafSeq(scafSeqFile, scaffoldSet, contigGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the statistics of contig lengths
	if(contigsLenStatisticsScaf(scaffoldSet, minContigLen)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute contig length statistics, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// free the memory
	freeMemScaffolding(&scaffoldSet, &scafContigIndex);

	gettimeofday(&tp_end,NULL);
	time_used = tp_end.tv_sec-tp_start.tv_sec + (double)(tp_end.tv_usec-tp_start.tv_usec)/1000000;
	printf("Scaffolding Used Time: %.2f Seconds.\n", time_used);

	printf("============= End scaffolding. =============\n");

	return SUCCESSFUL;
}

/**
 * Initialize the parameters for contigs linking.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemScaffolding(scaffoldSet_t **scaffoldSet)
{
	*scaffoldSet = (scaffoldSet_t *) calloc (1, sizeof(scaffoldSet_t));
	if((*scaffoldSet)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Free the scaffolding memory.
 */
void freeMemScaffolding(scaffoldSet_t **scaffoldSet, scafContigIndex_t **scafContigIndex)
{
	releaseScaffoldSet(scaffoldSet);
	releaseScafContigIndex(scafContigIndex);
}

/**
 * Release the scaffold set.
 */
void releaseScaffoldSet(scaffoldSet_t **scaffoldSet)
{
	int32_t i, j;
	scaffoldItem_t *scaffoldItem, *scaffoldItemTmp;

	scaffoldItem = (*scaffoldSet)->scaffoldItemList;
	while(scaffoldItem)
	{
		scaffoldItemTmp = scaffoldItem->next;
		free(scaffoldItem->contigOverlapArray);
		free(scaffoldItem);
		scaffoldItem = scaffoldItemTmp;
	}

	free(*scaffoldSet);
	*scaffoldSet = NULL;
}

/**
 * Load the contig graph from the contigs file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short loadContigGraph(contigGraph_t **contigGraph, char *contigsFile)
{
	if(initContigGraph(contigGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the contig graph, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill contigs
	if(fillContigItems(*contigGraph, contigsFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the contig graph items, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Fill the contig graph items from the contigs file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillContigItems(contigGraph_t *contigGraph, char *contigsFile)
{
	FILE *fpContigs;
	int64_t contigID, contigLen, maxContigLen, returnFlag, contigHeadLen;
	char *contigSeq, *contigHeadTitle;
	contigGraphItem_t *contigItem;

	fpContigs = fopen(contigsFile, "r");
	if(fpContigs==NULL)
	{
		printf("line=%d, In %s(), cannot open the contig file [ %s ], error!\n", __LINE__, __func__, contigsFile);
		return FAILED;
	}

	contigHeadTitle = (char *) malloc(1000 * sizeof(char));
	if(contigHeadTitle==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(getMaxContigLenFromFile(&maxContigLen, contigsFile)==FAILED)
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
	while((returnFlag=getSingleContigFastaFromFile(fpContigs, contigHeadTitle, contigSeq, &contigLen))==SUCCESSFUL)
	{
		if(contigLen>=minContigLen)
		{
			contigID ++;

			if(contigGraph->itemNumContigItemArray>=contigGraph->maxItemNumContigItemArray)
			{
				contigGraph->maxItemNumContigItemArray *= 2;
				contigGraph->contigItemArray = realloc(contigGraph->contigItemArray, contigGraph->maxItemNumContigItemArray * sizeof(contigGraphItem_t));
				if(contigGraph->contigItemArray==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}

			contigItem = contigGraph->contigItemArray + contigGraph->itemNumContigItemArray;
			contigHeadLen = strlen(contigHeadTitle);
			contigItem->contigTitle = (char *) calloc (contigHeadLen+1, sizeof(char));
			if(contigItem->contigTitle==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			contigItem->contigSeq = (char *) calloc (contigLen+1, sizeof(char));
			if(contigItem->contigSeq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			strcpy(contigItem->contigTitle, contigHeadTitle);
			strcpy(contigItem->contigSeq, contigSeq);
			contigItem->contigID = contigID;
			contigItem->contigLen = contigLen;

			contigGraph->itemNumContigItemArray ++;
		}
	}

	free(contigSeq);
	contigSeq = NULL;

	free(contigHeadTitle);
	contigHeadTitle = NULL;

	fclose(fpContigs);

	return SUCCESSFUL;
}

/**
 * Get the maximal contig length from contig file.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxContigLenFromFile(int64_t *maxContigLen, const char *contigFileName)
{
	FILE *fpContigs;
	int64_t maxSeqLen, contigLen;
	char ch, headTitle[1000];

	fpContigs = fopen(contigFileName, "r");
	if(fpContigs==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigFileName);
		return FAILED;
	}

	*maxContigLen = 0;
	ch = fgetc(fpContigs);
	if(ch!='>')
	{
		printf("The contig file [ %s ] was not in fasta format, error!\n", contigFileName);
		return FAILED;
	}

	// get the head title of a contig
	fscanf(fpContigs, "%s\n", headTitle);

	maxSeqLen = 0;
	contigLen = 0;
	while((ch=fgetc(fpContigs))!=EOF)
	{
		if(ch!='>')
		{ // base sequence
			if(ch!='\n')
			{
				contigLen ++;
			}
		}else
		{ // ends of a contig
			if(contigLen>maxSeqLen)
			{
				maxSeqLen = contigLen;
			}
			contigLen = 0;
		}
	}

	// process the last contig
	if(contigLen>0)
	{
		if(contigLen>maxSeqLen)
		{
			maxSeqLen = contigLen;
		}
		contigLen = 0;
	}else
	{
		printf("line=%d, In %s(), contigLen=%ld, error!\n", __LINE__, __func__, contigLen);
		return FAILED;
	}

	fclose(fpContigs);
	fpContigs = NULL;

	*maxContigLen = maxSeqLen;
	if((*maxContigLen)<=0)
	{
		printf("line=%d, In %s(), maxContigLen=%ld, error!\n", __LINE__, __func__, *maxContigLen);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get a contig from the fasta file.
 *  @return:
 *  	if succeed, return SUCCESSFUL; if end of file, return FAILED; else, return ERROR.
 */
short getSingleContigFastaFromFile(FILE *fp_contig, char *contigHeadTitle, char *contigSeq, int64_t *contigLen)
{
	char ch, *pseq;

	*contigLen = 0;
	ch = fgetc(fp_contig);
	if(ch=='>')//start of a new contig
	{
		ch = fgetc(fp_contig);
		pseq = contigHeadTitle;
		while(ch!='\n' && ch!=EOF)
		{
			*pseq = ch;
			pseq ++;

			ch = fgetc(fp_contig);
		}
		*pseq = '\0';

	}else // end of file, return FAILED
	{
		return FAILED;
	}

	pseq = contigSeq;
	ch = fgetc(fp_contig);
	while(ch!='>' && ch!=EOF)
	{
		switch(ch)
		{
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case 'N':
			case 'n':
			case '.':
				*pseq = ch;
				pseq ++;
				(*contigLen) ++;
		}
		ch = fgetc(fp_contig);
	}
	*pseq = '\0';

	if(ch=='>')
	{
		fseek(fp_contig, -1, SEEK_CUR);
	}

	return SUCCESSFUL;
}

/**
 * Output contig overlap information to text file.
 *  Format:
 *  	(1) Header fields: >scaffoldID, linkedNum, rowsNum, which are separated by tab character;
 *  	(2)   Body fields: contigID1, orientation1, contigLen1, contigID2, orientation2, contigLen2, mergeFlag, overlapLen, gapSize, breakFlag, which are separated by tab character.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputScaffoldSetToFile(char *tmpScafFile, scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph)
{
	FILE *fpScaf;
	int32_t i, j;
	scaffoldItem_t *scaffoldItem;
	contigOverlap_t *pContigOverlapInfo;
	int32_t scaffoldID, linkedContigNum, rowsNum, startRow;
	int32_t subRowsNum, subStartRow, subEndRow, subLinkedContigNum;

	fpScaf = fopen(tmpScafFile, "w");
	if(fpScaf==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, tmpScafFile);
		return FAILED;
	}

	// break the erroneous links and output the overlap information
	scaffoldItem = scaffoldSet->scaffoldItemList;
	while(scaffoldItem)
	{
		scaffoldID = scaffoldItem->scaffoldID;
		rowsNum = scaffoldItem->itemNumContigOverlapArray;
		linkedContigNum = scaffoldItem->linkedContigsNum;
		pContigOverlapInfo = scaffoldItem->contigOverlapArray;

		//printf("scaffoldID=%d\n", scaffoldID);

		fprintf(fpScaf, ">%d\t%d\t%d\n", scaffoldID, linkedContigNum, rowsNum);

		if(linkedContigNum>=2)
		{
			for(i=0; i<rowsNum; i++)
			{
				fprintf(fpScaf, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
						pContigOverlapInfo[i].contigID1, pContigOverlapInfo[i].orientation1, contigGraph->contigItemArray[pContigOverlapInfo[i].contigID1-1].contigLen,
						pContigOverlapInfo[i].contigID2, pContigOverlapInfo[i].orientation2, contigGraph->contigItemArray[pContigOverlapInfo[i].contigID2-1].contigLen,
						pContigOverlapInfo[i].mergeFlag, pContigOverlapInfo[i].overlapLen, pContigOverlapInfo[i].gapSize, pContigOverlapInfo[i].breakFlag, pContigOverlapInfo[i].pairNum);
			}
		}else
		{
			fprintf(fpScaf, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
					pContigOverlapInfo->contigID1, ORIENTATION_PLUS, contigGraph->contigItemArray[pContigOverlapInfo->contigID1-1].contigLen,
					0, 0, 0, 0, 0, 0, 0, 0);
		}

		scaffoldItem = scaffoldItem->next;
	}

	fclose(fpScaf);

	return SUCCESSFUL;
}

