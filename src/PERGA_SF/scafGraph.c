/*
 * scafGraph.c
 *
 *  Created on: Jul 26, 2011
 *      Author: zhuxiao
 */

#include "inc/scafStdinc.h"
#include "inc/scafExtvab.h"


/**
 * Build de Bruijn graph according to read lists(RL) and paired end files.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short buildGraph(const char *graphFile, const char *monoReadSeqFile, const char *monoReadListFile, const char *readListFile1, const char *readListFile2, const char *readSeqFile)
{
	printf("Begin building the graph in scaffolding ...\n");

	// generate the mono read list (MRL) to a binary file
	if(generateMonoReadList(monoReadListFile, readListFile1, readListFile2)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate the mono read list (MRL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the sequence of reads in mono read list (MRL)
	if(getMonoReadSequence(monoReadSeqFile, monoReadListFile, readSeqFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the mono read sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// construct the graph, and output it to a binary file
	if(constructScafGraph(graphFile, monoReadSeqFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot construct the scafGraphDeBruijn, error!\n", __LINE__, __func__);
		return FAILED;
	}

	printf("End building the graph in scaffolding.\n");

	return SUCCESSFUL;
}

/**
 * generate mono read list (MRL) and output it to a binary file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateMonoReadList(const char *monoReadListFile, const char *readListFile1, const char *readListFile2)
{
	printf("\tBegin generating the mono read list (MRL) in scaffolding ...\n");

	// initialize the memory for RL1, RL2 and MRL
	if(initMemGeMonoReadList(readListFile1, readListFile2)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory for generating the mono read list (MRL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// generating the mono read list (MRL)
	if(fillMonoReadList()==FAILED)
	{
		printf("line=%d, In %s(), cannot fill the data of mono read list (MRL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// output the mono read list (MRL) to file
	if(saveMonoReadListToFile(monoReadListFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot output the data of mono read list (MRL) to file [ %s ], error!\n", __LINE__, __func__, monoReadListFile);
		return FAILED;
	}

	freeMemGeMonoReadList();

	printf("\tEnd generating the mono read list (MRL) in scaffolding.\n");

	return SUCCESSFUL;
}

/**
 * Initialize the memory for generating the mono read list (MRL).
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemGeMonoReadList(const char *readListFile1, const char *readListFile2)
{
	// load readList1 and readList2
	if(loadSingleReadList(readListFile1, readListArr+1, readNumInRL+1, readPosArr+1, matchItemNumInRP+1)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the data of a read list (RL), error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(loadSingleReadList(readListFile2, readListArr+2, readNumInRL+2, readPosArr+2, matchItemNumInRP+2)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the data of a read list (RL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate the memory for mono read list (MRL)
	readNumInRL[0] = readNumInRL[1] + readNumInRL[2];
	matchItemNumInRP[0] = matchItemNumInRP[1] + matchItemNumInRP[2];

	readListArr[0] = (ReadList*) calloc(readNumInRL[0], sizeof(ReadList));
	if(readListArr[0]==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	readPosArr[0] = (ReadPos*) calloc(matchItemNumInRP[0], sizeof(ReadPos));
	if(readPosArr[0]==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Release the memory for generating the mono read list (MRL).
 */
void freeMemGeMonoReadList()
{
	// Release the memory of read list 1 and read list 2
	freeSingleReadList(readListArr+1, readNumInRL+1, readPosArr+1, matchItemNumInRP+1);
	freeSingleReadList(readListArr+2, readNumInRL+2, readPosArr+2, matchItemNumInRP+2);

	// free the mono read list (MRL)
	freeSingleReadList(readListArr, readNumInRL, readPosArr, matchItemNumInRP);
	monoReadListArr = NULL;
	monoReadPosArr = NULL;
}

/**
 * get mono read sequences and output them to a text file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMonoReadSequence(const char *monoReadSeqFile, const char *monoReadListFile, const char *readSeqFile)
{
	printf("\tBegin getting the sequences of mono reads ...\n");

	// Allocate memory and load the data of mono read list (MRL)
	if(loadSingleReadList(monoReadListFile, &monoReadListArr, &readItemNumInMRL, &monoReadPosArr, &matchItemNumInMRP)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the data of mono read list (MRL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	//Get the sequences of PE1 and PE2
	if(getMonoReadSeqFromReadSeqFile(monoReadSeqFile, readSeqFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the sequences of mono reads, error!\n", __LINE__, __func__);
		return FAILED;
	}


	// Release the memory of mono read list (MRL)
	freeSingleReadList(&monoReadListArr, &readItemNumInMRL, &monoReadPosArr, &matchItemNumInMRP);

	printf("\tEnd getting the sequences of mono reads.\n");

	return SUCCESSFUL;
}


/**
 * get mono read sequences from single PE file and output them to a text file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMonoReadSeqFromReadSeqFile(const char *monoReadSeqFile, const char *readSeqFile)
{
	FILE *fpMonoSeq, *fpReadSeq;
	uint64_t readID, paired_readID;
	int64_t hitRow;
	char *seq_data;

	seq_data = (char*) calloc((readLen+1), sizeof(char));
	if(seq_data==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	fpReadSeq = fopen(readSeqFile, "r");
	if(!fpReadSeq)
	{
		printf("line=%d, In %s(), Data File [%s] cannot open! Error Code:%d\n", __LINE__, __func__, readSeqFile, errno);
		return FAILED;
	}

	fpMonoSeq = fopen(monoReadSeqFile, "w");
	if(fpMonoSeq==NULL)
	{
		printf("line=%d, In %s(), cannot open the file [ %s ] , error!\n", __LINE__, __func__, monoReadSeqFile);
		return FAILED;
	}

	while(fscanf(fpReadSeq, "%lu\t%s\n", &readID, seq_data)!=EOF)
	{
		// generate the paired read ID
		if(readID%2==1) // odd number -> even number
		{
			paired_readID = readID + 1;
		}else	//even number -> odd number
		{
			paired_readID = readID - 1;
		}

		// search the read from mono read list (MRL)
		hitRow = getReadRowFromReadList(paired_readID, monoReadListArr, readItemNumInMRL);
		if(hitRow>=0)
		{ // find the paired read, then output itself
			fprintf(fpMonoSeq, "%lu\t%s\n", readID, seq_data);
		}
	}

	free(seq_data);
	seq_data = NULL;

	fclose(fpMonoSeq);
	fpMonoSeq = NULL;

	fclose(fpReadSeq);
	fpReadSeq = NULL;

	return SUCCESSFUL;
}

/**
 * Check whether the sequence contains unknown bases 'N' or '.'.
 *  @ return:
 *  	If it contains 'N' or '.', return YES; otherwise, return NO.
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
 *  Construct the graph in scaffolding and output it to a binary file.
 *  @ return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructScafGraph(const char *graphFile, const char *monoReadSeqFile)
{
	printf("\tBegin constructing the graph in scaffolding ...\n");

	//Initialize the De Bruijn graph in scaffolding.
	if(initGraphInScaf(&scafGrapDeBruijn) == FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the graph for gap filling, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//Count the element number of each k-mer in scaffolding.
	if(countKmersInScaf(monoReadSeqFile) == FAILED)
	{
		printf("line=%d, In %s(), cannot count k-mers from the file [%s], error!\n", __LINE__, __func__, monoReadSeqFile);
		return FAILED;
	}

	// add k-mers into graph
	if(addKmersInScaf(monoReadSeqFile) == FAILED)
	{
		printf("line=%d, In %s(), cannot add k-mers from the file [%s], error!\n", __LINE__, __func__, monoReadSeqFile);
		return FAILED;
	}


	//############################ debug information ######################################
#if DEBUG_FLAG
	if(checkScafKmersInGrapInScaf(monoReadSeqFile, scafGrapDeBruijn) == FAILED)
	{
		printf("line=%d, In %s(), there are errors in graph, error!\n", __LINE__, __func__);
		return FAILED;
	}
#endif
	//############################ debug information ######################################


	// save the graph to a binary file
	if(saveGraphInScaf(graphFile) == FAILED)
	{
		printf("line=%d, In %s(), cannot save the graph into the file [%s] in scaffolding, error!\n", __LINE__, __func__, graphFile);
		return FAILED;
	}


	// release the graph in scaffolding
	freeGraphInScaf(&scafGrapDeBruijn);

	printf("\tEnd constructing the graph in scaffolding.\n");

	return SUCCESSFUL;
}

/**
 * Get read length from readSeq file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadLenFromReadSeqFile(const char *monoReadSeqFile)
{
	FILE *fp;
	uint64_t tmp_readID;
	char tmp_readSeq[5000];

	fp = fopen(monoReadSeqFile, "r");
	if(fp==NULL)
	{
		printf("line=%d, In %s(), cannot open mono read sequence file [ %s ], error!\n", __LINE__, __func__, monoReadSeqFile);
		return FAILED;
	}

	if(fscanf(fp, "%lu\t%s\n", &tmp_readID, tmp_readSeq)!=2)
	{
		printf("line=%d, In %s(), cannot read mono read sequence file [ %s ], error!\n", __LINE__, __func__, monoReadSeqFile);
		return FAILED;
	}

	fclose(fp);
	fp = NULL;

	readLen = strlen(tmp_readSeq);
	if(readLen<kmerSize)
	{
		printf("line=%d, In %s(), readLen=%d, kmerSize=%d, error!\n", __LINE__, __func__, readLen, kmerSize);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 *  Count the element number of each k-mer in scaffolding.
 *  @ return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short countKmersInScaf(const char* pchFileName)
{
	FILE* srcfp = fopen(pchFileName, "r");
	if (!srcfp)
	{
		printf("line=%d, In %s(), Data File [%s] cannot open! Error Code:%d\n", __LINE__, __func__, pchFileName, errno);
		return FAILED;
	}

	char *readSeq = (char*) calloc(readLen+1, sizeof(char));
	if(readSeq==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	uint64_t readID;

	while( fscanf(srcfp, "%ld\t%s\n", &readID, readSeq) != EOF)
	{
		if(countKmersInReadInScaf(readSeq)==FAILED)
		{
			printf("line=%d, In %s(), cannot count the k-mer information for the read with sequence [ %s ], error!\n", __LINE__, __func__, readSeq);
			freeGraphInScaf(&scafGrapDeBruijn);
			return FAILED;
		}
	}

	free(readSeq);
	readSeq = NULL;

	fclose(srcfp);
	srcfp = NULL;

	return SUCCESSFUL;
}

/**
 *  Initialize the De Bruijn graph in scaffolding.
 *  @ return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initGraphInScaf(scafGraph **pScafGrapDeBruijn)
{
	*pScafGrapDeBruijn = (scafGraph *) malloc(sizeof(scafGraph));
	if(*pScafGrapDeBruijn==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory of graph, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	(*pScafGrapDeBruijn)->pkmers = ( scafKmer**)calloc(hashTableSize, sizeof( scafKmer* ));
	if( (*pScafGrapDeBruijn)->pkmers == NULL )
	{
		printf("line=%d, In %s(), cannot malloc the graph, Error!\n", __LINE__, __func__);
		free(*pScafGrapDeBruijn); *pScafGrapDeBruijn = NULL;
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 *  Count the element number of each k-mer in a single read in scaffolding.
 *  @ return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short countKmersInReadInScaf(char *seq)
{
	int i, j, baseInt;
	uint64_t hashcode;

	// generate the kmer integer sequence
	if(generateKmerSeqIntInScaf(kmerSeqInt, seq)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the hashcode and count the kmer occurrence
	hashcode = kmerhashInScaf(kmerSeqInt);
	if(countSingleKmerInScaf(hashcode, kmerSeqInt)==FAILED)
	{
		printf("line=%d, In %s(), cannot count single k-mer for the sequence [ %s ], error!\n", __LINE__, __func__, seq);
		return FAILED;
	}


	for(i=kmerSize; i<readLen; i++)
	{
		// get the baseInt
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

		// generate the kmer integer sequence
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				kmerSeqInt[j] = (kmerSeqInt[j] << 2) | (kmerSeqInt[j+1] >> 62);
			}

			kmerSeqInt[entriesPerKmer-2] = (kmerSeqInt[entriesPerKmer-2] << 2) | (kmerSeqInt[entriesPerKmer-1] >> (2*lastEntryBaseNumKmer-2));
		}
		kmerSeqInt[entriesPerKmer-1] = ((kmerSeqInt[entriesPerKmer-1] << 2) | baseInt) & lastEntryMaskKmer;

		// get the hashcode and count the kmer occurrence
		hashcode = kmerhashInScaf(kmerSeqInt);
		if(countSingleKmerInScaf(hashcode, kmerSeqInt)==FAILED)
		{
			printf("line=%d, In %s(), cannot count single k-mer for the sequence [ %s ], error!\n", __LINE__, __func__, seq);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Generate the kmer interger sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateKmerSeqIntInScaf(uint64_t *seqInt, char *seq)
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
 *  Get the hash code of a k-mer in scaffolding.
 *  @ return:
 *  	If succeeds, return the hash code; otherwise, return -1.
 */
uint64_t kmerhashInScaf(uint64_t *seqInt)
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

	for(j=0; j<lastEntryBaseNumKmer; j++)
	{
		hashcode += (hashcode << 5) | ((seqInt[entriesPerKmer-1] >> (2*(lastEntryBaseNumKmer-j-1))) & 3);
	}

	//return (hashcode & 0x7FFFFFFF) % hashTableSize;
	return hashcode % hashTableSize;
}

/**
 *  Count the element number of a single k-mer in scaffolding.
 *  @ return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short countSingleKmerInScaf(int hashcode, uint64_t *kmerSeqInt)
{
	uint64_t *tmpSeq;
	scafKmer *kmer;

	kmer = getKmerByHashInScaf(hashcode, kmerSeqInt, scafGrapDeBruijn);

	if(!kmer)
	{
		kmer = ( scafKmer* )malloc( sizeof( scafKmer ) );
		if( kmer == NULL )
		{
			printf("line=%d, In %s(), cannot allocate memory for the hash array, Error!\n", __LINE__, __func__);
			return FAILED;
		}

		kmer -> arraysize = 1;
		kmer -> multiplicity = 0;
		kmer -> ppos = NULL;

		tmpSeq = (uint64_t*) malloc(entriesPerKmer * sizeof(uint64_t));
		if(tmpSeq==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(memcpy(tmpSeq, kmerSeqInt, entriesPerKmer * sizeof(uint64_t))==NULL)
		{
			printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		kmer -> scafKmerseq = tmpSeq;

		kmer->next = scafGrapDeBruijn->pkmers[ hashcode ];
		scafGrapDeBruijn->pkmers[ hashcode ] = kmer;

	}else
	{
		kmer ->arraysize ++;
	}

	return SUCCESSFUL;
}

/**
 *  Add k-mers into graph in scaffolding.
 *  @ return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addKmersInScaf(const char* pchFileName)
{
	uint64_t readID;
	FILE* srcfp;

	srcfp = fopen(pchFileName, "r");
	if (!srcfp)
	{
		printf("line=%d, In %s(), Data File [%s] cannot open! Error Code:%d\n", __LINE__, __func__, pchFileName, errno);
		freeGraphInScaf(&scafGrapDeBruijn);
		return FAILED;;
	}

	char *readSeq = (char*) calloc(readLen+1, sizeof(char));
	if(readSeq==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	while(fscanf(srcfp, "%ld\t%s\n", &readID, readSeq) != EOF)
	{
		if(addReadInScaf(readSeq, readID)==FAILED)
		{
			printf("line=%d, In %s(), cannot add read with sequence [ %s ], error!\n", __LINE__, __func__, readSeq);
			freeGraphInScaf(&scafGrapDeBruijn);
			return FAILED;
		}
	}

	free(readSeq);
	readSeq = NULL;

	fclose(srcfp);
	srcfp = NULL;

	return SUCCESSFUL;
}

/**
 * Add read into the graph in scaffolding.
 *   @ return:
 *   	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addReadInScaf(char *seq, uint64_t rid)
{
	int i, j, baseInt, pos;
	uint64_t hashcode;

	// generate the kmer integer sequence
	if(generateKmerSeqIntInScaf(kmerSeqInt, seq)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the hashcode and count the kmer occurrence
	hashcode = kmerhashInScaf(kmerSeqInt);
	if(addKmerInScaf(hashcode, kmerSeqInt, rid, 1)==FAILED)
	{ // add kmer
		printf("line=%d, In %s(), addKmerInScaf() error!\n", __LINE__, __func__);
		return FAILED;
	}

	pos = 2;
	for(i=kmerSize; i<readLen; i++, pos++)
	{
		// get the baseInt
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

		// generate the kmer integer sequence
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				kmerSeqInt[j] = (kmerSeqInt[j] << 2) | (kmerSeqInt[j+1] >> 62);
			}
			kmerSeqInt[entriesPerKmer-2] = (kmerSeqInt[entriesPerKmer-2] << 2) | (kmerSeqInt[entriesPerKmer-1] >> (2*lastEntryBaseNumKmer-2));
		}
		kmerSeqInt[entriesPerKmer-1] = ((kmerSeqInt[entriesPerKmer-1] << 2) | baseInt) & lastEntryMaskKmer;

		// get the hashcode and add the kmer
		hashcode = kmerhashInScaf(kmerSeqInt);
		if(addKmerInScaf(hashcode, kmerSeqInt, rid, pos)==FAILED)
		{
			printf("line=%d, In %s(), error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * add a single k-mer to the de Bruijn graph in scaffolding。
 *   @ return:
 *     If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addKmerInScaf(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, uint16_t rpos)
{
	scafKmer *kmer;

	kmer = getKmerByHashInScaf(hashcode, kmerSeqInt, scafGrapDeBruijn);
	if(kmer)
	{
		if(kmer->ppos==NULL)
		{
			kmer->ppos = (scafRidpos*)calloc(kmer->arraysize, sizeof(scafRidpos));
			if(kmer->ppos==NULL)
			{
				printf("line=%d, In %s(), can not malloc the memory. Error and exit.\n", __LINE__, __func__);
				return FAILED;
			}
		}
		kmer->ppos[ kmer->multiplicity ].used = NO;
		kmer->ppos[ kmer->multiplicity ].pos = rpos;
		kmer->ppos[ kmer->multiplicity ].rid = rid;
		kmer->ppos[ kmer->multiplicity ].reserved = 0;
		kmer->multiplicity ++;
	}else
	{
		printf("line=%d, In %s(), can not get the kmer %s. Error and exit.\n", __LINE__, __func__, getKmerBaseByIntInScaf(kmerSeqInt));
		return FAILED;;
	}

	return SUCCESSFUL;
}

/**
 * Get the kmer bases from integer.
 */
char *getKmerBaseByIntInScaf(uint64_t *kmerSeqInt)
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
				case 0: kmerBaseSeq[k]='A'; break;
				case 1: kmerBaseSeq[k]='C'; break;
				case 2: kmerBaseSeq[k]='G'; break;
				case 3: kmerBaseSeq[k]='T'; break;
			}
			k ++;
		}
	}

	for(j=0; j<lastEntryBaseNumKmer; j++)
	{
		baseInt = (kmerSeqInt[entriesPerKmer-1] >> (2*(lastEntryBaseNumKmer-j-1))) & 3;
		switch(baseInt)
		{
			case 0: kmerBaseSeq[k]='A'; break;
			case 1: kmerBaseSeq[k]='C'; break;
			case 2: kmerBaseSeq[k]='G'; break;
			case 3: kmerBaseSeq[k]='T'; break;
		}
		k ++;
	}

	kmerBaseSeq[k] = '\0';

	return kmerBaseSeq;
}

scafKmer *getKmerInScaf(uint64_t *kmerSeqInt, scafGraph *graph)
{
	uint64_t hashcode;

	hashcode = kmerhashInScaf(kmerSeqInt);

	return getKmerByHashInScaf(hashcode, kmerSeqInt, graph);
}

/**
 * Get the scafKmer by hash code in local assembly.
 */
scafKmer *getKmerByHashInScaf(uint64_t hashvalue, uint64_t *kmerSeqInt, scafGraph *graph)
{
	scafKmer *kmer;

	kmer = graph->pkmers[hashvalue];
	while(kmer)
	{
		if(identicalKmerSeqInScaf(kmerSeqInt, kmer->scafKmerseq)==YES)
			break;

		kmer = kmer->next;
	}

	return kmer;
}

/**
 * Check whether the two sequence is identical.
 *  @return:
 *  	If identical, return YES; otherwise, return NO.
 */
short identicalKmerSeqInScaf(uint64_t *kmerSeqInt1, uint64_t *kmerSeqInt2)
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
 * 按照kmer碱基序列, 取得kmer.
 */
scafKmer *getKmerByBaseInScaf(char *str, scafGraph *graph)
{
	uint64_t hashcode;

	// generate the kmer integer sequence
	if(generateKmerSeqIntInScaf(kmerSeqInt, str)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
		return NULL;
	}

	hashcode = kmerhashInScaf(kmerSeqInt);

	return getKmerByHashInScaf(hashcode, kmerSeqInt, graph);
}

/**
 * Get the kmer integer sequence of a kmer by its kmer integer sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReverseKmerseqIntInScaf(uint64_t *kmerSeqIntRev, uint64_t *kmerSeqInt)
{
	int i, j, baseInt;
	int entryIndexRev, baseNum;

	for(i=0; i<entriesPerKmer; i++) kmerSeqIntRev[i] = 0;

	entryIndexRev = baseNum = 0;
	for(j=0; j<lastEntryBaseNumKmer; j++)
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
scafKmer *getReverseKmerInScaf(uint64_t *kmerSeqIntRev, uint64_t *kmerSeqInt, scafGraph *graph)
{
	uint64_t hashcode;

	if(getReverseKmerseqIntInScaf(kmerSeqIntRev, kmerSeqInt)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the reverse kmer, error!\n", __LINE__, __func__);
		exit(1);
	}

	hashcode = kmerhashInScaf(kmerSeqIntRev);

	return getKmerByHashInScaf(hashcode, kmerSeqIntRev, graph);
}


/**
 * Delete the scafKmer of read rid by marking its used flag.
 *   @ return:
 *     If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short delScafKmerByHashInScaf(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, uint32_t rpos, scafGraph *graph)
{
	scafRidpos *ridpostable = NULL;
	scafKmer *kmer;

	kmer = getKmerByHashInScaf(hashcode, kmerSeqInt, graph);
	if(!kmer)
	{ // the scafKmer does not exist, error
		//printf("line=%d, In %s(), can not delete the kmer (%lu,%u), it does not exist in the graph. Error!\n", __LINE__, __func__, rid, rpos);
		return FAILED;
	}

	ridpostable = kmer->ppos;
	if(!ridpostable)
	{ // the scafRidPos table does not exist, error
		printf("line=%d, In %s(), can not delete the kmer (%lu,%u), the position array does not exist in the graph. Error!\n", __LINE__, __func__, rid, rpos);
		return FAILED;
	}

	if(kmer->multiplicity==0)
	{ // all the scafKmers have been deleted, error
		//printf("line=%d, In %s(), can not delete the kmer (%lu, %u), all of its positions have been deleted.\n", __LINE__, __func__, rid, rpos);
		return FAILED;
	}

	int index = findDirectScafRidposIndexInScaf(rid, rpos, ridpostable, kmer->arraysize);//二分查找rid_pos_table
	if(index>=0)
	{  //find the scafRidpos

		if(ridpostable[index].used==0)
		{ // delete the scafKmer
			ridpostable[index].used = 1;
			kmer->multiplicity--;
		}else //the scafKmer has been deleted, error
		{
			printf("line=%d, In %s(), can not delete the kmer (%lu, %u), it has been deleted in the graph.\n", __LINE__, __func__, rid, rpos);
			return FAILED;
		}

	}else
	{  // the scafRidpos does not exist, error
		//printf("line=%d, In %s(), can not delete the kmer (%lu, %u), it does not exist in the graph.\n", __LINE__, __func__, rid, rpos);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Search the (rid, rpos) from scafRidpos table by binary search.
 *  @return:
 *  	If finds, return the entry index; otherwise, return -1;
 */
int32_t findDirectScafRidposIndexInScaf(uint64_t rid, uint32_t rpos, scafRidpos *ridpostable, uint32_t posNum)
{
	int32_t index = -1;
	int32_t startIndex = findStartScafRidposIndexInScaf(rid, ridpostable, posNum);//二分查找rid_pos_table
	if(startIndex>=0)
	{
		index = getExectScafRidposIndexInScaf(rid, rpos, startIndex, ridpostable, posNum);
	}
	return index;
}

/**
 * Search the start entry index of the read rid in scafRidpos table by binary search.
 *  @return:
 *  	If finds, return the start entry index; otherwise, return -1;
 */
int32_t findStartScafRidposIndexInScaf(uint64_t rid, scafRidpos *rid_pos_table, uint32_t posNum)
{
	int32_t left, right, middle, existFlag = 0;
	left = 0;
	right = posNum - 1;
	while(left<=right)
	{
		middle = (left+right)/2;
		if(rid==rid_pos_table[middle].rid)
		{
			existFlag = 1;
			break;
		}else if(rid>rid_pos_table[middle].rid)
			left = middle+1;
		else
			right = middle-1;
	}

	if(existFlag)
	{
		while(middle>0 && rid_pos_table[middle-1].rid==rid)  //find the first entry index of the read rid
			middle--;
		return middle; // return the entry index
	}
	return -1;
}

/**
 * Find the (rid, rpos) from the start entry index in scafRidpos table.
 *  @return:
 *  	If finds, return the start entry index; otherwise, return -1;
 */
int32_t getExectScafRidposIndexInScaf(uint64_t rid, uint32_t rpos, int32_t startIndex, scafRidpos *ridpostable, uint32_t posNum)
{
	int32_t i = startIndex, exectIndex = -1;
	while(i<posNum && ridpostable[i].rid==rid)
	{
		if(ridpostable[i].pos==rpos)
		{
			exectIndex = i;
			break;
		}
		//ridpostable[i].reserved = 1;
		i++;
	}

	return exectIndex;
}

/**
 *  Output the graph to file.
 *   File format:
 *   	(1) arraySize for kmer array, arraySize for ridposArray, and readLen, kmerSize, hashtableSize;
 *   	(2) kmerArray;
 *   	(3) ridposArray.
 *   	(4) kmerseqArray.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short saveGraphInScaf(const char* graphFile)
{
	struct kmerArrayNode{
		uint64_t hashCode;
		uint64_t arraysize;
	};
	struct ridposArrayNode{
		uint64_t rid:48;
		uint64_t rpos:16;
	};

	uint64_t kmerArraySize, ridposArraySize, tmp[5];
	uint64_t i, j, rowsNumKmerArray, rowsNumRidposArray, itemNumKmerSeq;
	struct kmerArrayNode kmerArray;
	struct ridposArrayNode ridposArray;
	scafKmer *kmer;
	scafRidpos *ridpos;
	FILE *fpGraph;

	fpGraph = fopen(graphFile, "wb");
	if(fpGraph==NULL)
	{
		printf("line=%d, In %s(), can not open the file [ %s ], Error!\n", __LINE__, __func__, graphFile);
		return FAILED;
	}

	kmerArraySize = ridposArraySize = 0;
	for(i=0; i<hashTableSize; i++)
	{
		kmer = scafGrapDeBruijn->pkmers[i];
		while(kmer)
		{
			kmerArraySize ++;
			ridposArraySize += kmer->arraysize;

			kmer = kmer->next;
		}
	}

	tmp[0] = kmerArraySize;
	tmp[1] = ridposArraySize;
	tmp[2] = readLen;
	tmp[3] = kmerSize;
	tmp[4] = hashTableSize;

	if(fwrite(tmp, sizeof(uint64_t), 5, fpGraph)!=5)
	{
		printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
		return FAILED;
	}

	rowsNumKmerArray = 0;
	for(i=0; i<hashTableSize; i++)
	{
		kmer = scafGrapDeBruijn->pkmers[i];
		while(kmer)
		{
			kmerArray.hashCode = i;
			kmerArray.arraysize = kmer->arraysize;
			rowsNumKmerArray ++;

			if(fwrite(&kmerArray, sizeof(struct kmerArrayNode), 1, fpGraph)!=1)
			{
				printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
				return FAILED;
			}

			kmer = kmer->next;
		}
	}

	rowsNumRidposArray = 0;
	for(i=0; i<hashTableSize; i++)
	{
		kmer = scafGrapDeBruijn->pkmers[i];
		while(kmer)
		{
			ridpos = kmer->ppos;
			for(j=0; j<kmer->arraysize; j++)
			{
				ridposArray.rid = ridpos[j].rid;
				ridposArray.rpos = ridpos[j].pos;
				rowsNumRidposArray ++;

				if(fwrite(&ridposArray, sizeof(struct ridposArrayNode), 1, fpGraph)!=1)
				{
					printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
					return FAILED;
				}
			}

			kmer = kmer->next;
		}
	}

	itemNumKmerSeq = 0;
	for(i=0; i<hashTableSize; i++)
	{
		kmer = scafGrapDeBruijn->pkmers[i];
		while(kmer)
		{
			if(fwrite(kmer->scafKmerseq, sizeof(uint64_t), entriesPerKmer, fpGraph)!=entriesPerKmer)
			{
				printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
				return FAILED;
			}
			itemNumKmerSeq ++;

			kmer = kmer->next;
		}
	}

	// ############################ Debug information ##############################
	if(kmerArraySize!=rowsNumKmerArray)
	{
		printf("line=%d, In %s(),kmerArraySize=%lu != rowsNumKmerArray=%lu. Error\n", __LINE__, __func__, kmerArraySize, rowsNumKmerArray);
		return FAILED;
	}
	if(itemNumKmerSeq!=rowsNumKmerArray)
	{
		printf("line=%d, In %s(),itemNumKmerSeq=%lu != rowsNumKmerArray=%lu. Error\n", __LINE__, __func__, itemNumKmerSeq, rowsNumKmerArray);
		return FAILED;
	}
	if(ridposArraySize!=rowsNumRidposArray)
	{
		printf("line=%d, In %s(),ridposArraySize=%lu != rowsNumRidposArray=%lu. Error\n", __LINE__, __func__, ridposArraySize, rowsNumRidposArray);
		return FAILED;
	}
	// ############################ Debug information ##############################

	fclose(fpGraph);
	fpGraph = NULL;

	return SUCCESSFUL;
}

/**
 *  Load the graph to memory.
 *   File format:
 *   	(1) arraySize for kmer array, arraySize for ridposArray, and readLen, kmerSize, hashtableSize;
 *   	(2) kmerArray;
 *   	(3) ridposArray.
 *   	(4) kmerseqArray.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short loadGraphInScaf(scafGraph **pScafGrapDeBruijn, const char* graphFile)
{
	struct kmerArrayNode{
		uint64_t hashCode;
		uint64_t arraysize;
	};
	struct ridposArrayNode{
		uint64_t rid:48;
		uint64_t rpos:16;
	};

	FILE *fpGraph;
	uint64_t kmerArraySize, ridposArraySize, tmp[5], hashcode, arraysize;
	uint64_t preHashCode;
	scafKmer *preKmer;
	uint64_t i, j;
	struct kmerArrayNode kmerArray;
	struct ridposArrayNode ridposArray;
	scafKmer *kmer;
	scafRidpos *ridpos;

	printf("Begin loading the graph ...\n");

	fpGraph = fopen(graphFile, "rb");
	if(fpGraph==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, graphFile);
		return FAILED;
	}

	// get the item number of the kmer array and ridpos array, respectively
	if(fread(tmp, sizeof(uint64_t), 5, fpGraph)!=5)
	{
		printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
		return FAILED;
	}

	kmerArraySize = tmp[0];
	ridposArraySize = tmp[1];
	readLen = tmp[2];
	kmerSize = tmp[3];
	hashTableSize = tmp[4];

	if(initGraphInScaf(pScafGrapDeBruijn) == FAILED)
	{
		printf("line=%d, In %s(), initialize graph error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the data
	preHashCode = -1;
	preKmer = NULL;
	for(i=0; i<kmerArraySize; i++)
	{
		if(fread(&kmerArray, sizeof(struct kmerArrayNode), 1, fpGraph)!=1)
		{
			printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
			return FAILED;
		}

		hashcode = kmerArray.hashCode;
		kmer = (scafKmer*) malloc(sizeof(scafKmer));
		if(!kmer)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		kmer->arraysize = kmer->multiplicity = kmerArray.arraysize;
		kmer->ppos = NULL;
		kmer->scafKmerseq = NULL;
		kmer->next = NULL;

		if(hashcode==preHashCode)
		{ // they are in the same list
			preKmer->next = kmer;
		}else
		{ // the first node in the list
			preHashCode = hashcode;
			(*pScafGrapDeBruijn)->pkmers[hashcode] = kmer;
		}
		//preHashCode = hashcode;
		preKmer = kmer;
	}

	for(i=0; i<hashTableSize; i++)
	{
		kmer = (*pScafGrapDeBruijn)->pkmers[i];
		while(kmer)
		{
			arraysize = kmer->arraysize;
			ridpos = (scafRidpos*) malloc(arraysize * sizeof(scafRidpos));
			if(!ridpos)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			for(j=0; j<arraysize; j++)
			{
				if(fread(&ridposArray, sizeof(struct ridposArrayNode), 1, fpGraph)!=1)
				{
					printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
					return FAILED;
				}

				ridpos[j].rid = ridposArray.rid;
				ridpos[j].pos = ridposArray.rpos;
				ridpos[j].used = NO;
				ridpos[j].reserved = 0;
			}
			kmer->ppos = ridpos;

			kmer = kmer->next;
		}
	}

	for(i=0; i<hashTableSize; i++)
	{
		kmer = (*pScafGrapDeBruijn)->pkmers[i];
		while(kmer)
		{
			kmer->scafKmerseq = (uint64_t*) malloc(entriesPerKmer*sizeof(uint64_t));
			if(!kmer->scafKmerseq)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(fread(kmer->scafKmerseq, sizeof(uint64_t), entriesPerKmer, fpGraph)!=entriesPerKmer)
			{
				printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
				return FAILED;
			}

			kmer = kmer->next;
		}
	}

	fclose(fpGraph);
	fpGraph = NULL;

	printf("End loading the graph.\n");

	return SUCCESSFUL;
}

/**
 * Release the memory of de Bruijn graph in scaffolding.
 */
void freeGraphInScaf(scafGraph **pScafGrapDeBruijn)
{
	uint64_t i;
	scafKmer *kmer, *head;

	for(i=0; i<hashTableSize; i++)
	{
		head = (*pScafGrapDeBruijn)->pkmers[i];  //取得kmer
		while(head)
		{
			kmer = head->next;
			free(head->ppos);
			free(head->scafKmerseq);
			free(head);
			head = kmer;
		}
		(*pScafGrapDeBruijn)->pkmers[i] = NULL;
	}

	free((*pScafGrapDeBruijn)->pkmers);
	(*pScafGrapDeBruijn)->pkmers = NULL;
	free(*pScafGrapDeBruijn);
	*pScafGrapDeBruijn = NULL;
}

