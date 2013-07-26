/*
 * scafContigIndex.c
 *
 *  Created on: Jun 11, 2011
 *      Author: zhuxiao
 */

#include "inc/scafStdinc.h"
#include "inc/scafExtvab.h"



//======= the following section is contig index (CI) building methods ========

/**
 * Build contig index from the contigs with fasta format.
 *  @return:
 *  	if succeed, return SUCCESSFUL; otherwise, return FAILED;
 */
short buildContigIndex(const char *contigFileName, const char *contigIndexFile)
{
	printf("=========== Begin building Contig Index (CI), please wait ...\n");

	//initialize Memory
	if(initMemory(contigFileName)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//extract sequences with length READ_LEN from both ends of contigs
	if(extractSequencesFromContigEnds(contigFileName)==FAILED)
	{
		printf("line=%d, In %s(), cannot extract sequences from both ends of contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//sort the base sequence in ascending order
	if(seqSort()==FAILED)
	{
		printf("line=%d, In %s(), cannot sort the sequences, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ####################### Debug information #################################
#if DEBUG_FLAG
	// check the sorting result
	if(checkSortResult()==FAILED)
	{
		printf("line=%d, In %s(), error sorting result!\n", __LINE__, __func__);
		return FAILED;
	}
#endif
	// ####################### Debug information #################################


	// allocate the memory
	if(initConvertingMem()==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the converting memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//convert the sequence indices to contig index (CI)
	if(convertToContigIndex()==FAILED)
	{
		printf("line=%d, In %s(), cannot convert the sequence indeices to contig indices (CI), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ####################### Debug information #################################
#if DEBUG_OUT_FLAG
	printf("itemNumSeqArr=%d, uniqueSeqNum=%d, uniqueRatio=%f\n", itemNumSeqArr, uniqueSeqNum, (double)uniqueSeqNum/itemNumSeqArr);
	printf("last item: startRow=%d, rowsNum=%d\n", seqRowArr[uniqueSeqNum-1].startRow, seqRowArr[uniqueSeqNum-1].rowsNum);
#endif
	// ####################### Debug information #################################

	//save contig index to file
	if(saveContigIndexToFile(contigIndexFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot save contig index (CI), error!\n", __LINE__, __func__);
		return FAILED;
	}

	//free the memory of base sequence index array
	freeMemAfterConveting();

	freeMemContigIndex(&uniqueSeqKmerArr, &uniqueSeqArr, &seqRowArr, &uniqueSeqNum, &contigMatchInfoArr, &itemNumSeqArr);


	printf("=========== End Built Contig Index (CI).\n");

	return SUCCESSFUL;
}

/**
 * Initialize the memory.
 *  @return:
 *  	if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemory(const char *contigFileName)
{
	//get the contigs number
	contigsNum = getContigsNum(contigFileName);
	if(contigsNum<=0)
	{
		printf("line=%d, In %s(), contigsNum=%d, error!\n", __LINE__, __func__, contigsNum);
		return FAILED;
	}

	//calculate the memory size of base sequence array and allocate it
	maxItemNumSeqArr = ((contigEndLen-readLen+1) * 2) * contigsNum;

	baseSeqIndexArr = (baseSeqIndex*) calloc(maxItemNumSeqArr, sizeof(baseSeqIndex));
	if(baseSeqIndexArr==NULL)
	{
		printf("line=%d, In %s(), cannot calloc the sequence index array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	seqArr = (uint64_t*) calloc(maxItemNumSeqArr, entriesPerRead * sizeof(uint64_t));
	if(seqArr==NULL)
	{
		printf("line=%d, In %s(), cannot calloc the sequence array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	itemNumSeqArr = 0; // set the initialized item number to zero

	return SUCCESSFUL;
}


/**
 * Get the number of contigs.
 * 	@return:
 * 		if succeed, return the number; otherwise, return ERROR.
 */
int getContigsNum(const char *contigFileName)
{
	int tmp_contigsNum = 0;
	FILE *fp_contig;
	char ch;
	fp_contig = fopen(contigFileName, "r");
	if(fp_contig==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigFileName);
		return ERROR;
	}

	ch = fgetc(fp_contig);
	while(ch!=EOF)
	{
		if(ch=='>')
			tmp_contigsNum ++;

		ch = fgetc(fp_contig);
	}

	fclose(fp_contig);

	return tmp_contigsNum;
}


/**
 * Extract sequences with length READ_LEN from both ends of contigs.
 * 	@return:
 * 		if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short extractSequencesFromContigEnds(char *contigFileName)
{
	FILE *fp_contig;
	int64_t contigID, contigLen, maxContigLen, returnFlag;
	char *contigSeq, *contigHeadTitle;

	fp_contig = fopen(contigFileName, "r");
	if(fp_contig==NULL)
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
	while((returnFlag=getSingleContigFastaFromFile(fp_contig, contigHeadTitle, contigSeq, &contigLen))==SUCCESSFUL)
	{
		contigID ++;

		// ############################ Debug information #######################
		//if(contigID==10355)
		//{
		//	printf("line=%d, In %s(), contigID=%ld, contigLen=%ld\n", __LINE__, __func__, contigID, contigLen);
		//}
		// ############################ Debug information #######################

		// fill short sequences
		if(fillBaseSeqIndex(contigID, contigSeq, contigLen)==FAILED)
		{
			printf("line=%d, In %s(), cannot fill the base sequence index array, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	free(contigSeq);
	contigSeq = NULL;

	free(contigHeadTitle);
	contigHeadTitle = NULL;

	fclose(fp_contig);
	fp_contig = NULL;

	// handle the error situation
	if(returnFlag==ERROR)
	{
		printf("line=%d, In %s(), cannot get the contig: %ld, error!\n", __LINE__, __func__, contigID);
		return FAILED;
	}

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
	if(maxContigLen<=0)
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
 * Fill the base sequences of a contig to the sequence array.
 *  @return:
 *  	if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillBaseSeqIndex(int contigID, char *contigSeq, int contigLen)
{
	uint64_t *pSeqArr, *p_newSeqArr;
	int i, startContigPos, endContigPos, basePos, tmp_baseInt, contigEndFlag;

	// handle 5' end sequences
	pSeqArr = seqArr + itemNumSeqArr * entriesPerRead;
	startContigPos = 1;
	if(contigLen>contigEndLen)
	{
		endContigPos = contigEndLen;
		contigEndFlag = 0;
	}else
	{
		endContigPos = contigLen;
		contigEndFlag = 2;	// only the 5' end bacause of the short (<=300 bp) contig
	}

	if(seqHash(pSeqArr, contigSeq, startContigPos)==FAILED) // hash the first sequence
	{
		printf("line=%d, In %s(), cannot hash the sequence of contig: %d\n", __LINE__, __func__, contigID);
		return FAILED;
	}
	baseSeqIndexArr[itemNumSeqArr].contigEnd = contigEndFlag;
	baseSeqIndexArr[itemNumSeqArr].contigID = contigID;
	baseSeqIndexArr[itemNumSeqArr].contigPos = startContigPos ++;
	baseSeqIndexArr[itemNumSeqArr].startRow = itemNumSeqArr * entriesPerRead;
	itemNumSeqArr ++;

	// begin hashing other sequences
	basePos = readLen;
	while(basePos<endContigPos)  // hash other sequences
	{
		switch(contigSeq[ basePos ])
		{
			case 'A':
			case 'a':
				tmp_baseInt = 0;
				break;
			case 'C':
			case 'c':
				tmp_baseInt = 1;
				break;
			case 'G':
			case 'g':
				tmp_baseInt = 2;
				break;
			case 'T':
			case 't':
				tmp_baseInt = 3;
				break;
			default:
				printf("line=%d, In %s(), unknown base %c, Error!\n", __LINE__, __func__, contigSeq[basePos]);
				return FAILED;
		}
		basePos ++;

		p_newSeqArr = pSeqArr + entriesPerRead;

		// update the base sequence
		for(i=0; i<entriesPerRead-2; i++)
		{
			p_newSeqArr[i] = (pSeqArr[i] << 2) | (pSeqArr[i+1] >> 62);
		}
		p_newSeqArr[entriesPerRead-2] = (pSeqArr[entriesPerRead-2] << 2) | (pSeqArr[entriesPerRead-1] >> ((lastEntryBaseNumRead-1) * 2));
		p_newSeqArr[entriesPerRead-1] = ((pSeqArr[entriesPerRead-1] << 2) | tmp_baseInt) & lastEntryMaskRead;

		//record the sequence
		baseSeqIndexArr[itemNumSeqArr].contigEnd = contigEndFlag;
		baseSeqIndexArr[itemNumSeqArr].contigID = contigID;
		baseSeqIndexArr[itemNumSeqArr].contigPos = startContigPos ++;
		baseSeqIndexArr[itemNumSeqArr].startRow = itemNumSeqArr * entriesPerRead;
		itemNumSeqArr ++;
		pSeqArr = p_newSeqArr;
	}

	// check the contig length to decide whether the 3' end should be handled
	if(contigLen<=contigEndLen)
	{ // do not the 3' end
		return SUCCESSFUL;
	}else
	{ // handle the 3' end
		contigEndFlag = 1;  // the 3' end of the contig
	}

	// handle 3' end sequences
	pSeqArr = seqArr + itemNumSeqArr * entriesPerRead;
	startContigPos = contigLen - contigEndLen + 1;
	endContigPos = contigLen;

	if(seqHash(pSeqArr, contigSeq, startContigPos)==FAILED) // hash the first sequence
	{
		printf("line=%d, In %s(), cannot hash the sequence of contig: %d\n", __LINE__, __func__, contigID);
		return FAILED;
	}
	baseSeqIndexArr[itemNumSeqArr].contigEnd = contigEndFlag;
	baseSeqIndexArr[itemNumSeqArr].contigID = contigID;
	baseSeqIndexArr[itemNumSeqArr].contigPos = startContigPos ++;
	baseSeqIndexArr[itemNumSeqArr].startRow = itemNumSeqArr * entriesPerRead;
	itemNumSeqArr ++;

	// begin hashing other sequences
	basePos = startContigPos + readLen - 2;
	while(basePos<endContigPos)  // hash other sequences
	{
		switch(contigSeq[ basePos ])
		{
			case 'A':
			case 'a':
				tmp_baseInt = 0;
				break;
			case 'C':
			case 'c':
				tmp_baseInt = 1;
				break;
			case 'G':
			case 'g':
				tmp_baseInt = 2;
				break;
			case 'T':
			case 't':
				tmp_baseInt = 3;
				break;
			default:
				printf("line=%d, In %s(), unknown base %c, Error!\n", __LINE__, __func__, contigSeq[basePos]);
				return FAILED;
		}
		basePos ++;

		p_newSeqArr = pSeqArr + entriesPerRead;

		// update the base sequence
		for(i=0; i<entriesPerRead-2; i++)
		{
			p_newSeqArr[i] = (pSeqArr[i] << 2) | (pSeqArr[i+1] >> 62);
		}
		p_newSeqArr[entriesPerRead-2] = (pSeqArr[entriesPerRead-2] << 2) | (pSeqArr[entriesPerRead-1] >> ((lastEntryBaseNumRead-1) * 2));
		p_newSeqArr[entriesPerRead-1] = ((pSeqArr[entriesPerRead-1] << 2) | tmp_baseInt) & lastEntryMaskRead;

		//record the sequence
		baseSeqIndexArr[itemNumSeqArr].contigEnd = contigEndFlag;
		baseSeqIndexArr[itemNumSeqArr].contigID = contigID;
		baseSeqIndexArr[itemNumSeqArr].contigPos = startContigPos ++;
		baseSeqIndexArr[itemNumSeqArr].startRow = itemNumSeqArr * entriesPerRead;
		itemNumSeqArr ++;
		pSeqArr = p_newSeqArr;
	}

	return SUCCESSFUL;
}


/**
 * Sequence hash function.
 *  @return:
 *  	if succeed, return SUCCESSFUL; otherwise, return FAILED.
 *
 *  Notes: startContigPos starts from 1.
 */
short seqHash(uint64_t *hashArray, char *contigSeq, int startContigPos)
{
	int i, entryPos, baseNum, tmp_startPos, tmp_endPos;

	entryPos = 0;
	baseNum = 0;
	tmp_startPos = startContigPos - 1;
	tmp_endPos = tmp_startPos + readLen - 1;

	// reset the hashArray
	for(i=0; i<entriesPerRead; i++)
		hashArray[i] = 0;

	// begin hashing
	for(i=tmp_startPos; i<=tmp_endPos; i++)
	{
		switch(contigSeq[ i ])
		{
			case 'A':
			case 'a':
				hashArray[entryPos] = hashArray[entryPos] << 2;
				break;
			case 'C':
			case 'c':
				hashArray[entryPos] = (hashArray[entryPos] << 2) | 1;
				break;
			case 'G':
			case 'g':
				hashArray[entryPos] = (hashArray[entryPos] << 2) | 2;
				break;
			case 'T':
			case 't':
				hashArray[entryPos] = (hashArray[entryPos] << 2) | 3;
				break;
			default:
				printf("In seqHash(), unknown base %c, Error!\n", contigSeq[i]);
				return FAILED;
		}

		baseNum ++;

		if(baseNum==32)
		{
			entryPos ++;
			baseNum = 0;
		}
	}

	return SUCCESSFUL;
}

/**
 * Sequence sort.
 *  @return:
 *  	if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short seqSort()
{
	//printf("Begin sorting the sequences ...\n");

	//initialize the sort memory
	if(initSortMem()==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the sort memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill k-mer data
	if(fillKmerData()==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the row data of each k-mer, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Begin sorting
	if(kmersSort()==FAILED)
	{
		printf("line=%d, In %s(), cannot sort k-mers, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// update the base sequence indices
	if(updateBaseSeqArr()==FAILED)
	{
		printf("line=%d, In %s(), cannot update base sequence index array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Free sort memory
	freeSortMem();

	//printf("End sorting the sequences.\n");

	return SUCCESSFUL;
}

/**
 * Initialize the sort memory.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short initSortMem()
{
	sortKmerArr = (sortkmertype *) calloc(SORT_KMER_ARR_SIZE, sizeof(sortkmertype));
	if(sortKmerArr==NULL)
	{
		printf("In initSortMem(), cannot allocate the memory, error!\n");
		return FAILED;
	}

	sortRowIndexArr = (int *) calloc(maxItemNumSeqArr, sizeof(int));
	if(sortRowIndexArr==NULL)
	{
		printf("In initSortMem(), cannot allocate the memory, error!\n");
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Free the sort memory.
 */
void freeSortMem()
{
	free(sortKmerArr);
	sortKmerArr = NULL;

	free(sortRowIndexArr);
	sortRowIndexArr = NULL;
}


/**
 * Fill the k-mers.
 *  @return:
 *  	if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillKmerData()
{
	int i;
	uint64_t *pSeqArr;
	sortkmertype *pSortKmer;
	int total;

	// initialize the right shift bits number
	if(initRightShiftBitsNum()==FAILED)
	{
		printf("In fillKmerData(), cannot initialize the right shift bits number, error!\n");
		return FAILED;
	}

	// compute the array size of each k-mer
	for(i=0; i<itemNumSeqArr; i++)
	{
		pSeqArr = seqArr + baseSeqIndexArr[i].startRow;
		sortKmerArr[ pSeqArr[0] >> rightShift_bitsNum ].arraysize ++;
	}

	// compute the start row of row index array for each k-mer
	total = 0;
	for(i=0; i<SORT_KMER_ARR_SIZE; i++)
	{
		if(sortKmerArr[i].arraysize>0)
		{
			sortKmerArr[i].startSortRow = total;
			total += sortKmerArr[i].arraysize;
		}
	}

	// fill data
	for(i=0; i<itemNumSeqArr; i++)
	{
		pSeqArr = seqArr + baseSeqIndexArr[i].startRow;

		pSortKmer = sortKmerArr + (pSeqArr[0] >> rightShift_bitsNum);
		sortRowIndexArr[ pSortKmer->startSortRow + pSortKmer->num ] = i;
		pSortKmer->num ++;
	}

	return SUCCESSFUL;
}

/**
 * Initialize the right shift bits number.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short initRightShiftBitsNum()
{
	// initialize the right shift bits number
	if(entriesPerRead>1)
	{
		rightShift_bitsNum = 64 - 2 * SORT_KMER_SIZE;
	}else if(lastEntryBaseNumRead>=SORT_KMER_SIZE)
	{
		rightShift_bitsNum = (lastEntryBaseNumRead - SORT_KMER_SIZE) * 2;
	}else
	{
		printf("line=%d, In %s(), too short sequence length: %d\n", __LINE__, __func__, lastEntryBaseNumRead);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Sort all k-mers.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short kmersSort()
{
	sortkmertype *pSortKmer;
	int i;

	for(i=0; i<SORT_KMER_ARR_SIZE; i++)
	{
		pSortKmer = sortKmerArr + i;
		if(pSortKmer->arraysize > 1)
		{ // sort the k-mer having array size larger than 1
			if(kmerSelectionSort(sortRowIndexArr+pSortKmer->startSortRow, pSortKmer->arraysize)==FAILED)
			{
				printf("line=%d, In %s(), cannot sort kmer: %d\n", __LINE__, __func__, i);
				return FAILED;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * K-mer selection sort.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short kmerSelectionSort(int *rowArr, int rowsNum)
{
	//printf("Begin k-mer selection sort ...\n");

	int i, j, k, tmp;
	int minRow;
	uint64_t *pSeqArr_i, *pSeqArr_j, *pSeqArr_min;

	for(i=0; i<rowsNum-1; i++)
	{
		pSeqArr_i = seqArr + baseSeqIndexArr[ rowArr[i] ].startRow;

		minRow = i;
		pSeqArr_min = pSeqArr_i;

		for(j=i+1; j<rowsNum; j++)
		{
			pSeqArr_j = seqArr + baseSeqIndexArr[ rowArr[j] ].startRow;
			for(k=0; k<entriesPerRead; k++)
			{
				if(pSeqArr_min[k] > pSeqArr_j[k])
				{
					pSeqArr_min = pSeqArr_j;
					minRow = j;

					break;
				}else if(pSeqArr_min[k] < pSeqArr_j[k])
				{
					break;
				}
			}
		}

		if(minRow!=i)
		{

			for(k=0; k<entriesPerRead; k++)
			{
				if(pSeqArr_min[k]<pSeqArr_i[k])
				{
					break;
				}else if(pSeqArr_min[k]>pSeqArr_i[k])
				{
					printf("In kmerSelectionSort(), error!\n");
					return FAILED;
				}
			}

			if(k<entriesPerRead)
			{ // Exchange the two items
				tmp = rowArr[i];
				rowArr[i] = rowArr[minRow];
				rowArr[minRow] = tmp;
			}else
			{
				printf("In kmerSelectionSort(), error!\n");
				return FAILED;
			}
		}
	}

	//printf("End k-mer selection sort.\n");

	return SUCCESSFUL;
}

/**
 * Update the base sequence index array.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateBaseSeqArr()
{
	int i;
	baseSeqIndex *tmp_baseSeqIndexArr;

	// allocate the auxiliary array
	tmp_baseSeqIndexArr = (baseSeqIndex *) calloc(maxItemNumSeqArr, sizeof(baseSeqIndex));
	if(tmp_baseSeqIndexArr==NULL)
	{
		printf("In updateBaseSeqArr(), cannot allocate the memory, error!\n");
		return FAILED;
	}

	// fill the sorted data
	for(i=0; i<itemNumSeqArr; i++)
	{
		if(memcpy(tmp_baseSeqIndexArr+i, baseSeqIndexArr+sortRowIndexArr[i], sizeof(baseSeqIndex))==NULL)
		{
			printf("line=%d, In updateBaseSeqArr(), cannot copy memory, error!\n", __LINE__);
			return FAILED;
		}
	}

	// copy the auxiliary array to the base sequence index array
	if(memcpy(baseSeqIndexArr, tmp_baseSeqIndexArr, itemNumSeqArr*sizeof(baseSeqIndex))==NULL)
	{
		printf("line=%d, In updateBaseSeqArr(), cannot copy memory, error!\n", __LINE__);
		return FAILED;
	}

	// release the auxiliary array
	free(tmp_baseSeqIndexArr);
	tmp_baseSeqIndexArr = NULL;


	//================= update the sequence array ==============
	// allocate the auxiliary array
	uint64_t *tmp_seqArr;
	tmp_seqArr = (uint64_t *) calloc(maxItemNumSeqArr, entriesPerRead*sizeof(uint64_t));
	if(tmp_seqArr==NULL)
	{
		printf("In updateBaseSeqArr(), cannot allocate the memory, error!\n");
		return FAILED;
	}

	// copy the sorted data to temporary sequence array
	int startRow = 0;
	for(i=0; i<itemNumSeqArr; i++)
	{
		if(memcpy(tmp_seqArr+startRow, seqArr+baseSeqIndexArr[i].startRow, entriesPerRead*sizeof(uint64_t))==NULL)
		{
			printf("line=%d, In updateBaseSeqArr(), cannot copy memory, error!\n", __LINE__);
			return FAILED;
		}

		baseSeqIndexArr[i].startRow = startRow;
		startRow += entriesPerRead;
	}

	// copy the auxiliary array to the base sequence array
	if(memcpy(seqArr, tmp_seqArr, itemNumSeqArr*entriesPerRead*sizeof(uint64_t))==NULL)
	{
		printf("line=%d, In updateBaseSeqArr(), cannot copy memory, error!\n", __LINE__);
		return FAILED;
	}

	// release the auxiliary array
	free(tmp_baseSeqIndexArr);
	tmp_baseSeqIndexArr = NULL;
	free(tmp_seqArr);
	tmp_seqArr = NULL;

	return SUCCESSFUL;
}

/**
 * Convert the sequence indices to contig index (CI).
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short convertToContigIndex()
{
	int i, j, k;
	int sameItemNum;
	int tmp_rowsNum, tmp_totalSeqNum;
	uint64_t *pSeqArr_i, *pSeqArr_j, *pUniqueSeq;

	// initialize the right shift bits number
	if(initRightShiftBitsNum()==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the right shift bits number, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// begin converting
	uniqueSeqNum = 0;
	i = 0;
	while(i<itemNumSeqArr)
	{
		sameItemNum = 1;
		pSeqArr_i = seqArr + baseSeqIndexArr[i].startRow;
		for(j=i+1; j<itemNumSeqArr; j++)
		{
			pSeqArr_j = seqArr + baseSeqIndexArr[j].startRow;
			for(k=0; k<entriesPerRead; k++)
			{
				if(pSeqArr_i[k] != pSeqArr_j[k])
				{
					break;
				}
			}

			if(k==entriesPerRead)
			{
				sameItemNum ++;
			}else
			{
				break;
			}
		}

		// ========= record the sequence data to contig index (CI) ==========
		if(memcpy(uniqueSeqArr+uniqueSeqNum*entriesPerRead, pSeqArr_i, entriesPerRead*sizeof(uint64_t))==NULL)
		{
			printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		tmp_rowsNum = sameItemNum;

		seqRowArr[uniqueSeqNum].startRow = i;
		seqRowArr[uniqueSeqNum].rowsNum = tmp_rowsNum;
		uniqueSeqNum ++;

		for(j=0; j<tmp_rowsNum; j++)
		{
			contigMatchInfoArr[i+j].contigID = baseSeqIndexArr[i+j].contigID;
			contigMatchInfoArr[i+j].contigPos = baseSeqIndexArr[i+j].contigPos;
			contigMatchInfoArr[i+j].contigEnd = baseSeqIndexArr[i+j].contigEnd;
		}

		i += tmp_rowsNum;
	}

	// compute the uniqueSeqKmerArr
	pUniqueSeq = uniqueSeqArr;
	for(i=0; i<uniqueSeqNum; i++)
	{
		uniqueSeqKmerArr[ pUniqueSeq[0] >> rightShift_bitsNum ].itemNum ++;
		pUniqueSeq += entriesPerRead;
	}

	tmp_totalSeqNum = 0;
	for(i=0; i<SORT_KMER_ARR_SIZE; i++)
	{
		if(uniqueSeqKmerArr[i].itemNum>0)
		{
			uniqueSeqKmerArr[i].itemRow = tmp_totalSeqNum;
			tmp_totalSeqNum += uniqueSeqKmerArr[i].itemNum;
		}
	}

	//####################### Debug information #####################
#if DEBUG_OUT_FLAG
	printf("uniqueSeqNum=%d, tmp_totalSeqNum=%d\n", uniqueSeqNum, tmp_totalSeqNum);
#endif
	//####################### Debug information #####################

	//####################### Debug information #####################
#if DEBUG_FLAG
	for(i=0; i<SORT_KMER_ARR_SIZE; i++)
	{
		if(uniqueSeqKmerArr[i].itemNum>0)
		{
			pUniqueSeq = uniqueSeqArr + uniqueSeqKmerArr[i].itemRow * entriesPerRead;
			if((pUniqueSeq[0]>>rightShift_bitsNum)!=i)
			{
				printf("line=%d, In %s(), error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}
	//printf("congratulations.\n");
#endif
	//####################### Debug information #####################

	return SUCCESSFUL;
}

/**
 * Initialize the converting memory.
 *  @return:
 *  	If succees, return SUCCESSFUL; otherwise, return FAILED.
 */
short initConvertingMem()
{
	uniqueSeqKmerArr = (uniqueSeqKmer*) calloc (SORT_KMER_ARR_SIZE, sizeof(uniqueSeqKmer));
	if(uniqueSeqKmerArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	uniqueSeqArr = (uint64_t*) calloc (itemNumSeqArr, entriesPerRead * sizeof(uint64_t));
	if(uniqueSeqArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	seqRowArr = (seqRow*) calloc (itemNumSeqArr, sizeof(seqRow));
	if(seqRowArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	contigMatchInfoArr = (contigMatchInfo*) calloc (itemNumSeqArr, sizeof(contigMatchInfo));
	if(contigMatchInfoArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Free the memory of temp sequence index array.
 */
void freeMemAfterConveting()
{
	free(seqArr);
	seqArr = NULL;

	free(baseSeqIndexArr);
	baseSeqIndexArr = NULL;
}



/**
 * Search the seq from contig index (CI).
 *  @return:
 *  	If hits, return its row number in sequence row array;
 *  	otherwise, return -1;
 */
int seqSearch(uint64_t *seq)
{
	int i, k, itemNum, uniqueSeqRow;
	uniqueSeqKmer *pUinqueSeqKmer;
	uint64_t *pUniqueSeq;

	pUinqueSeqKmer = uniqueSeqKmerArr +  (seq[0] >> rightShift_bitsNum);

	if(pUinqueSeqKmer->itemNum==0)
	{
		return -1;
	}

	itemNum = pUinqueSeqKmer->itemNum;
	uniqueSeqRow = pUinqueSeqKmer->itemRow;
	pUniqueSeq = uniqueSeqArr + uniqueSeqRow * entriesPerRead;

	// here we only allow the exact match, which also can be replaced by inexact match later.
	for(i=0; i<itemNum; i++)
	{
		for(k=0; k<entriesPerRead; k++)
		{
			if(seq[k]!=pUniqueSeq[k])
			{
				break;
				//return -1;
			}
		}

		if(k==entriesPerRead)
		{ // hits the target sequence, then return the row number in sequence row array
			return uniqueSeqRow + i;
		}

		pUniqueSeq += entriesPerRead;
	}

	return -1;
}


/**
 * Save the contig index to file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 *
 *  Format:
 *  	Four parts:
 *  	(1) three item numbers (int):
 *  		[0]- the valid k-mer number; and
 *  		[1]- the unique sequence number; and
 *  		[2]- the contig information item number.
 *  	(2) the unique sequence array data; and
 *  	(3) the unique unique sequence row array data; and
 *  	(4) the contig information item array data.
 */
short saveContigIndexToFile(const char *contigIndexFile)
{
	FILE *fpContigIndex;
	int i, validKmerNum;

	fpContigIndex = fopen(contigIndexFile, "wb");
	if(fpContigIndex==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigIndexFile);
		return FAILED;
	}

	// compute the valid k-mer number
	validKmerNum = 0;
	for(i=0; i<SORT_KMER_ARR_SIZE; i++)
	{
		if(uniqueSeqKmerArr[i].itemNum>0)
			validKmerNum ++;
	}

	// save the vaild k-mer number , the unique sequence number, and the total sequence item number in sequence number
	if(fwrite(&validKmerNum, sizeof(int), 1, fpContigIndex)!=1)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fwrite(&uniqueSeqNum, sizeof(int), 1, fpContigIndex)!=1)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fwrite(&itemNumSeqArr, sizeof(int), 1, fpContigIndex)!=1)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	// save the k-mer array
	// format:
	//	hash code (int), uniquence sequence k-mer node (item row in unique sequence array, and item number)
	for(i=0; i<SORT_KMER_ARR_SIZE; i++)
	{
		if(uniqueSeqKmerArr[i].itemNum>0)
		{
			if(fwrite(&i, sizeof(int), 1, fpContigIndex)!=1)
			{
				printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(fwrite(uniqueSeqKmerArr+i, sizeof(uniqueSeqKmer), 1, fpContigIndex)!=1)
			{
				printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	// save the unique sequence array and the sequence row array
	if(fwrite(uniqueSeqArr, entriesPerRead*sizeof(uint64_t), uniqueSeqNum, fpContigIndex)!=uniqueSeqNum)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fwrite(seqRowArr, sizeof(seqRow), uniqueSeqNum, fpContigIndex)!=uniqueSeqNum)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	// save the contig information array
	if(fwrite(contigMatchInfoArr, sizeof(contigMatchInfo), itemNumSeqArr, fpContigIndex)!=itemNumSeqArr)
	{
		printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
		return FAILED;
	}

	fclose(fpContigIndex);
	fpContigIndex = NULL;


	return SUCCESSFUL;
}


/**
 * Load contig index (CI).
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 *
 *  Format:
 *  	Four parts:
 *  	(1) three item numbers (int):
 *  		[0]- the valid k-mer number, and
 *  		[1]- the unique sequence number, and
 *  		[2]- the contig information item number;
 *  	(2) the unique sequence array data; and
 *  	(3) the unique unique sequence row array data; and
 *  	(4) the contig information item array data.
 */
short loadContigIndex(const char *contigIndexFile, uniqueSeqKmer **pUniqueSeqKmerArr, uint64_t **pUniqueSeqArr, seqRow **pSeqRowArr, int *uniqueSeqNum, contigMatchInfo **pContigMatchInfoArr, int *itemNumContigMatchInfo)
{
	FILE *fpContigIndex;
	int i, validKmerNum, hashcode;

	fpContigIndex = fopen(contigIndexFile, "rb");
	if(fpContigIndex==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigIndexFile);
		return FAILED;
	}


	// initialize the global variables
	if(initGlobalVariables()==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the right shift bits number, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the three numbers
	if(fread(&validKmerNum, sizeof(int), 1, fpContigIndex)!=1)
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fread(uniqueSeqNum, sizeof(int), 1, fpContigIndex)!=1)
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(fread(itemNumContigMatchInfo, sizeof(int), 1, fpContigIndex)!=1)
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}


	// allocate memory
	*pUniqueSeqKmerArr = (uniqueSeqKmer *) calloc (SORT_KMER_ARR_SIZE, sizeof(uniqueSeqKmer));
	if(*pUniqueSeqKmerArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	*pUniqueSeqArr = (uint64_t *) calloc (*uniqueSeqNum, entriesPerRead*sizeof(uint64_t));
	if(*pUniqueSeqArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	*pSeqRowArr = (seqRow *) calloc (*uniqueSeqNum, sizeof(seqRow));
	if(*pSeqRowArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	*pContigMatchInfoArr = (contigMatchInfo *) calloc (*itemNumContigMatchInfo, sizeof(contigMatchInfo));
	if(*pContigMatchInfoArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the unique sequence k-mer data
	for(i=0; i<validKmerNum; i++)
	{
		if(fread(&hashcode, sizeof(int), 1, fpContigIndex)!=1)
		{
			printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(fread((*pUniqueSeqKmerArr)+hashcode, sizeof(uniqueSeqKmer), 1, fpContigIndex)!=1)
		{
			printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// fill the sequence array
	if(fread(*pUniqueSeqArr, entriesPerRead*sizeof(uint64_t), *uniqueSeqNum, fpContigIndex)!=(*uniqueSeqNum))
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the sequence row array
	if(fread(*pSeqRowArr, sizeof(seqRow), *uniqueSeqNum, fpContigIndex)!=(*uniqueSeqNum))
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the contig information array
	if(fread(*pContigMatchInfoArr, sizeof(contigMatchInfo), *itemNumContigMatchInfo, fpContigIndex)!=(*itemNumContigMatchInfo))
	{
		printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
		return FAILED;
	}

	fclose(fpContigIndex);
	fpContigIndex = NULL;


	return SUCCESSFUL;
}

/**
 * Initialize the global variables.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initGlobalVariables()
{
	if(initRightShiftBitsNum()==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the right shift bits number, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Release memory of contig index (CI).
 */
void freeMemContigIndex(uniqueSeqKmer **pUniqueSeqKmerArr, uint64_t **pUniqueSeqArr, seqRow **pSeqRowArr, int *uniqueSeqNum, contigMatchInfo **pContigMatchInfoArr, int *itemNumContigInfo)
{
	*uniqueSeqNum = 0;
	*itemNumContigInfo = 0;

	free(*pUniqueSeqKmerArr);
	*pUniqueSeqKmerArr = NULL;

	free(*pUniqueSeqArr);
	*pUniqueSeqArr = NULL;

	free(*pSeqRowArr);
	*pSeqRowArr = NULL;

	free(*pContigMatchInfoArr);
	*pContigMatchInfoArr = NULL;
}
