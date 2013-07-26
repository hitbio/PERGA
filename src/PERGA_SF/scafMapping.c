/*
 * scafMapping.c
 *
 *  Created on: Jun 16, 2011
 *      Author: zhuxiao
 */

#include "inc/scafStdinc.h"
#include "inc/scafExtvab.h"


/**
 * Map the paired reads in two files, and save results to their corresponding files.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short mapPEs(const char *readSeqFile, const char *readMatchFile1, const char *readMatchFile2, char **readFilesInput, int readFileNum, int readsFileFormatType, int pairedMode, const char *contigIndexFile)
{
	printf("=========== Begin mapping the paired end files, please wait ...\n");

	if(loadContigIndex(contigIndexFile, &uniqueSeqKmerArr, &uniqueSeqArr, &seqRowArr, &uniqueSeqNum, &contigMatchInfoArr, &itemNumSeqArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot build contig index (CI), error!\n", __LINE__, __func__);
		return FAILED;
	}

	//testSeqSearch();

	if(readsFileFormatType==FILE_FORMAT_FASTA)
	{
		if(pairedMode==1)
		{
			if(mapPEFastaSeparate(readSeqFile, readMatchFile1, readMatchFile2, readFilesInput, readFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot map paired end reads, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else if(pairedMode==2)
		{
			if(mapPEFastaInterleaved(readSeqFile, readMatchFile1, readMatchFile2, readFilesInput, readFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot map paired end reads, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else
		{
			printf("line=%d, In %s(), invalid paired mode: %d, error!\n", __LINE__, __func__, pairedMode);
			return FAILED;
		}

	}else if(readsFileFormatType==FILE_FORMAT_FASTQ)
	{
		if(pairedMode==1)
		{
			if(mapPEFastqSeparate(readSeqFile, readMatchFile1, readMatchFile2, readFilesInput, readFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot map paired end reads, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else if(pairedMode==2)
		{
			if(mapPEFastqInterleaved(readSeqFile, readMatchFile1, readMatchFile2, readFilesInput, readFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot map paired end reads, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else
		{
			printf("line=%d, In %s(), invalid paired mode: %d, error!\n", __LINE__, __func__, pairedMode);
			return FAILED;
		}
	}else
	{
		printf("Reads files should be in fastq format.\n");
		return FAILED;
	}

	freeMemContigIndex(&uniqueSeqKmerArr, &uniqueSeqArr, &seqRowArr, &uniqueSeqNum, &contigMatchInfoArr, &itemNumSeqArr);

	printf("=========== End mapped the paired end file.\n");

	return SUCCESSFUL;
}

/**
 * Map the paired reads in two separate files, and save results to file with a given name.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 *
 *  Format:
 *		Each matched information item is saved as a readMatchInfoTemp structure.
 */
short mapPEFastaSeparate(const char *readSeqFile, const char *readMatchFile1, const char *readMatchFile2, char **readFilesInput, int readFileNum)
{
	FILE *fpPE[2], *fpPE_Result[2], *fpReadSeq;
	uint64_t matchedReadNum[2], matchedContigNum[2];
	int matchFlag[2];
	uint64_t i, j, fileID, readID;		// readID starts from 1

	readBuf_t *pReadBuf1, *pReadBuf2;
	uint64_t readsNum1, readsNum2;

	short errorFlag[2];
	char *readSeq[2];

	struct tmpReadNode
	{
		uint64_t readID;
		char *seq;
	};

	struct tmpReadNode *tmpReadBuf;
	uint64_t tmpReadNum;

	fpReadSeq = fopen(readSeqFile, "w");
	if(fpReadSeq==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readSeqFile);
		return FAILED;
	}


	fpPE_Result[0] = fopen(readMatchFiles[0], "w");
	if(fpPE_Result[0]==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readMatchFiles[0]);
		return FAILED;
	}
	fpPE_Result[1] = fopen(readMatchFiles[1], "w");
	if(fpPE_Result[1]==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readMatchFiles[1]);
		return FAILED;
	}

	// initialize the memory of read buffers
	if(initMemReadsBuf(&pReadBuf1)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the read buffer, error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(initMemReadsBuf(&pReadBuf2)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the read buffer, error!\n", __LINE__, __func__);
		return FAILED;
	}

	tmpReadBuf = (struct tmpReadNode *) calloc (MAX_READ_BUF_SIZE, sizeof(struct tmpReadNode));
	if(tmpReadBuf==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory for temporary reads buffer, error!\n", __LINE__, __func__);
		return FAILED;
	}
	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		tmpReadBuf[i].seq = (char *) calloc (readLen+1, sizeof(char));
		if(tmpReadBuf[i].seq==NULL)
		{
			printf("line=%d, In %s(), cannot allocate the memory for temporary reads, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	matchedReadNum[0] = matchedContigNum[0] = 0;
	matchedReadNum[1] = matchedContigNum[1] = 0;
	readID = 1;
	tmpReadNum = 0;

	for(fileID=0; fileID<readFileNum; fileID+=2)
	{
		fpPE[0] = fopen(readFilesInput[fileID], "r");
		if(fpPE[0]==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readFilesInput[fileID]);
			return FAILED;
		}
		fpPE[1] = fopen(readFilesInput[fileID+1], "r");
		if(fpPE[1]==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readFilesInput[fileID+1]);
			return FAILED;
		}

		//================== begin Mapping =================
		//get single read
		while(1)
		{
			// check the end of files
			if(feof(fpPE[0]) && feof(fpPE[1]))
			{
				break;
			}else if(feof(fpPE[0]) || feof(fpPE[1]))
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// fille the reads to reads buffers
			if(fillReadsToBufFasta(fpPE[0], pReadBuf1, &readsNum1, readLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(fillReadsToBufFasta(fpPE[1], pReadBuf2, &readsNum2, readLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(readsNum1!=readsNum2)
			{
				printf("line=%d, In %s(), readID=%lu, readsNum1=%lu != readsNum2=%lu, error!\n", __LINE__, __func__, readID, readsNum1, readsNum2);
				return FAILED;
			}

			// map the reads in reads buffers
			for(i=0; i<readsNum1; i++)
			{
				readSeq[0] = pReadBuf1[i].seq;
				readSeq[1] = pReadBuf2[i].seq;

				// ####################### Debug information #######################
				//if(readID==6947)
				//{
				//	printf("readID=%lu\n", readID);
				//}
				// ####################### Debug information #######################

				errorFlag[0] = errorFlag[1] = NO;

				// check unknown bases
				if(pReadBuf1[i].len<readLen || containUnknownBase(readSeq[0])==YES)
				{
					errorFlag[0] = YES;
				}
				if(pReadBuf2[i].len<readLen || containUnknownBase(readSeq[1])==YES)
				{
					errorFlag[1] = YES;
				}

				if((errorFlag[0]==YES) && (errorFlag[1]==YES))
				{
					readID += 2;
					continue;
				}

				matchFlag[0] = matchFlag[1] = NO;	// set the match flag of the read to NO
				if(errorFlag[0]==NO)
				{
					if(mapSingleRead(readID, readSeq[0], matchFlag, fpPE_Result[0], matchedContigNum, matchedReadNum)==FAILED)
					{
						printf("line=%d, In %s(), cannot map single read, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

				if(errorFlag[1]==NO)
				{
					if(mapSingleRead(readID+1, readSeq[1], matchFlag+1, fpPE_Result[1], matchedContigNum+1, matchedReadNum+1)==FAILED)
					{
						printf("line=%d, In %s(), cannot map single read, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

				if(matchFlag[0]==YES || matchFlag[1]==YES)
				{
					if(errorFlag[0]==NO)
					{
						//fprintf(fpReadSeq, "%lu\t%s\n", readID, readSeq[0]);

						tmpReadBuf[tmpReadNum].readID = readID;
						strcpy(tmpReadBuf[tmpReadNum].seq, readSeq[0]);
						tmpReadNum ++;

						if(tmpReadNum==MAX_READ_BUF_SIZE)
						{
							for(j=0; j<tmpReadNum; j++)
								fprintf(fpReadSeq, "%lu\t%s\n", tmpReadBuf[j].readID, tmpReadBuf[j].seq);
							tmpReadNum = 0;
						}
					}

					if(errorFlag[1]==NO)
					{
						//fprintf(fpReadSeq, "%lu\t%s\n", readID+1, readSeq[1]);

						tmpReadBuf[tmpReadNum].readID = readID+1;
						strcpy(tmpReadBuf[tmpReadNum].seq, readSeq[1]);
						tmpReadNum ++;

						if(tmpReadNum==MAX_READ_BUF_SIZE)
						{
							for(j=0; j<tmpReadNum; j++)
								fprintf(fpReadSeq, "%lu\t%s\n", tmpReadBuf[j].readID, tmpReadBuf[j].seq);
							tmpReadNum = 0;
						}
					}
				}

				readID += 2;
			}
		}

		fclose(fpPE[0]);
		fpPE[0] = NULL;
		fclose(fpPE[1]);
		fpPE[1] = NULL;
	}

	// process the remained temporary reads
	if(tmpReadNum>0)
	{
		for(j=0; j<tmpReadNum; j++)
			fprintf(fpReadSeq, "%lu\t%s\n", tmpReadBuf[j].readID, tmpReadBuf[j].seq);
		tmpReadNum = 0;
	}

	// free memory of read buffers
	freeMemReadsBuf(&pReadBuf1);
	freeMemReadsBuf(&pReadBuf2);

	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		free(tmpReadBuf[i].seq);
		tmpReadBuf[i].seq = NULL;
	}
	free(tmpReadBuf);
	tmpReadBuf = NULL;


	fclose(fpReadSeq);
	fpReadSeq = NULL;
	fclose(fpPE_Result[0]);
	fpPE_Result[0] = NULL;
	fclose(fpPE_Result[1]);
	fpPE_Result[1] = NULL;


#if DEBUG_FLAG
	printf("matchedReadNum[0]=%lu, matchedContigNum[0]=%lu\n", matchedReadNum[0], matchedContigNum[0]);
	printf("matchedReadNum[1]=%lu, matchedContigNum[1]=%lu\n", matchedReadNum[1], matchedContigNum[1]);
	printf("total reads: %lu\n", readID-1);
#endif

	return SUCCESSFUL;
}

/**
 * Map the paired reads in two separate files, and save results to file with a given name.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 *
 *  Format:
 *		Each matched information item is saved as a readMatchInfoTemp structure.
 */
short mapPEFastqSeparate(const char *readSeqFile, const char *readMatchFile1, const char *readMatchFile2, char **readFilesInput, int readFileNum)
{
	FILE *fpPE[2], *fpPE_Result[2], *fpReadSeq;
	uint64_t matchedReadNum[2], matchedContigNum[2];
	int matchFlag[2];
	uint64_t i, j, fileID, readID;		// readID starts from 1

	readBuf_t *pReadBuf1, *pReadBuf2;
	uint64_t readsNum1, readsNum2;

	short errorFlag[2];
	char *readSeq[2];

	struct tmpReadNode
	{
		uint64_t readID;
		char *seq;
	};

	struct tmpReadNode *tmpReadBuf;
	uint64_t tmpReadNum;

	fpReadSeq = fopen(readSeqFile, "w");
	if(fpReadSeq==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readSeqFile);
		return FAILED;
	}


	fpPE_Result[0] = fopen(readMatchFiles[0], "w");
	if(fpPE_Result[0]==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readMatchFiles[0]);
		return FAILED;
	}
	fpPE_Result[1] = fopen(readMatchFiles[1], "w");
	if(fpPE_Result[1]==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readMatchFiles[1]);
		return FAILED;
	}

	// initialize the memory of read buffers
	if(initMemReadsBuf(&pReadBuf1)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the read buffer, error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(initMemReadsBuf(&pReadBuf2)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the read buffer, error!\n", __LINE__, __func__);
		return FAILED;
	}

	tmpReadBuf = (struct tmpReadNode *) calloc (MAX_READ_BUF_SIZE, sizeof(struct tmpReadNode));
	if(tmpReadBuf==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory for temporary reads buffer, error!\n", __LINE__, __func__);
		return FAILED;
	}
	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		tmpReadBuf[i].seq = (char *) calloc (readLen+1, sizeof(char));
		if(tmpReadBuf[i].seq==NULL)
		{
			printf("line=%d, In %s(), cannot allocate the memory for temporary reads, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	matchedReadNum[0] = matchedContigNum[0] = 0;
	matchedReadNum[1] = matchedContigNum[1] = 0;
	readID = 1;
	tmpReadNum = 0;

	for(fileID=0; fileID<readFileNum; fileID+=2)
	{
		fpPE[0] = fopen(readFilesInput[fileID], "r");
		if(fpPE[0]==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readFilesInput[fileID]);
			return FAILED;
		}
		fpPE[1] = fopen(readFilesInput[fileID+1], "r");
		if(fpPE[1]==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readFilesInput[fileID+1]);
			return FAILED;
		}

		//================== begin Mapping =================
		//get single read
		while(1)
		{
			// check the end of files
			if(feof(fpPE[0]) && feof(fpPE[1]))
			{
				break;
			}else if(feof(fpPE[0]) || feof(fpPE[1]))
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// fille the reads to reads buffers
			if(fillReadsToBufFastq(fpPE[0], pReadBuf1, &readsNum1, readLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(fillReadsToBufFastq(fpPE[1], pReadBuf2, &readsNum2, readLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(readsNum1!=readsNum2)
			{
				printf("line=%d, In %s(), readID=%lu, readsNum1=%lu != readsNum2=%lu, error!\n", __LINE__, __func__, readID, readsNum1, readsNum2);
				return FAILED;
			}

			// map the reads in reads buffers
			for(i=0; i<readsNum1; i++)
			{
				readSeq[0] = pReadBuf1[i].seq;
				readSeq[1] = pReadBuf2[i].seq;

				// ####################### Debug information #######################
				//if(readID==6947)
				//{
				//	printf("readID=%lu\n", readID);
				//}
				// ####################### Debug information #######################

				errorFlag[0] = errorFlag[1] = NO;

				// check unknown bases
				if(pReadBuf1[i].len<readLen || containUnknownBase(readSeq[0])==YES)
				{
					errorFlag[0] = YES;
				}
				if(pReadBuf2[i].len<readLen || containUnknownBase(readSeq[1])==YES)
				{
					errorFlag[1] = YES;
				}

				if((errorFlag[0]==YES) && (errorFlag[1]==YES))
				{
					readID += 2;
					continue;
				}

				matchFlag[0] = matchFlag[1] = NO;	// set the match flag of the read to NO
				if(errorFlag[0]==NO)
				{
					if(mapSingleRead(readID, readSeq[0], matchFlag, fpPE_Result[0], matchedContigNum, matchedReadNum)==FAILED)
					{
						printf("line=%d, In %s(), cannot map single read, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

				if(errorFlag[1]==NO)
				{
					if(mapSingleRead(readID+1, readSeq[1], matchFlag+1, fpPE_Result[1], matchedContigNum+1, matchedReadNum+1)==FAILED)
					{
						printf("line=%d, In %s(), cannot map single read, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

				if(matchFlag[0]==YES || matchFlag[1]==YES)
				{
					if(errorFlag[0]==NO)
					{
						//fprintf(fpReadSeq, "%lu\t%s\n", readID, readSeq[0]);

						tmpReadBuf[tmpReadNum].readID = readID;
						strcpy(tmpReadBuf[tmpReadNum].seq, readSeq[0]);
						tmpReadNum ++;

						if(tmpReadNum==MAX_READ_BUF_SIZE)
						{
							for(j=0; j<tmpReadNum; j++)
								fprintf(fpReadSeq, "%lu\t%s\n", tmpReadBuf[j].readID, tmpReadBuf[j].seq);
							tmpReadNum = 0;
						}
					}

					if(errorFlag[1]==NO)
					{
						//fprintf(fpReadSeq, "%lu\t%s\n", readID+1, readSeq[1]);

						tmpReadBuf[tmpReadNum].readID = readID+1;
						strcpy(tmpReadBuf[tmpReadNum].seq, readSeq[1]);
						tmpReadNum ++;

						if(tmpReadNum==MAX_READ_BUF_SIZE)
						{
							for(j=0; j<tmpReadNum; j++)
								fprintf(fpReadSeq, "%lu\t%s\n", tmpReadBuf[j].readID, tmpReadBuf[j].seq);
							tmpReadNum = 0;
						}
					}
				}

				readID += 2;
			}
		}

		fclose(fpPE[0]);
		fpPE[0] = NULL;
		fclose(fpPE[1]);
		fpPE[1] = NULL;
	}

	// process the remained temporary reads
	if(tmpReadNum>0)
	{
		for(j=0; j<tmpReadNum; j++)
			fprintf(fpReadSeq, "%lu\t%s\n", tmpReadBuf[j].readID, tmpReadBuf[j].seq);
		tmpReadNum = 0;
	}

	// free memory of read buffers
	freeMemReadsBuf(&pReadBuf1);
	freeMemReadsBuf(&pReadBuf2);

	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		free(tmpReadBuf[i].seq);
		tmpReadBuf[i].seq = NULL;
	}
	free(tmpReadBuf);
	tmpReadBuf = NULL;


	fclose(fpReadSeq);
	fpReadSeq = NULL;
	fclose(fpPE_Result[0]);
	fpPE_Result[0] = NULL;
	fclose(fpPE_Result[1]);
	fpPE_Result[1] = NULL;


#if DEBUG_FLAG
	printf("matchedReadNum[0]=%lu, matchedContigNum[0]=%lu\n", matchedReadNum[0], matchedContigNum[0]);
	printf("matchedReadNum[1]=%lu, matchedContigNum[1]=%lu\n", matchedReadNum[1], matchedContigNum[1]);
	printf("total reads: %lu\n", readID-1);
#endif

	return SUCCESSFUL;
}

/**
 * Map the paired reads in single files, and save results to file with a given name.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 *
 *  Format:
 *		Each matched information item is saved as a readMatchInfoTemp structure.
 */
short mapPEFastaInterleaved(const char *readSeqFile, const char *readMatchFile1, const char *readMatchFile2, char **readFilesInput, int readFileNum)
{
	FILE *fpPE, *fpPE_Result[2], *fpReadSeq;
	uint64_t matchedReadNum[2], matchedContigNum[2];
	int matchFlag[2];
	uint64_t i, j, fileID, readID;		// readID starts from 1

	readBuf_t *pReadBuf;
	uint64_t readsNum;

	short errorFlag[2];
	char *readSeq[2];

	struct tmpReadNode
	{
		uint64_t readID;
		char *seq;
	};

	struct tmpReadNode *tmpReadBuf;
	uint64_t tmpReadNum;

	fpReadSeq = fopen(readSeqFile, "w");
	if(fpReadSeq==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readSeqFile);
		return FAILED;
	}


	fpPE_Result[0] = fopen(readMatchFiles[0], "w");
	if(fpPE_Result[0]==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readMatchFiles[0]);
		return FAILED;
	}
	fpPE_Result[1] = fopen(readMatchFiles[1], "w");
	if(fpPE_Result[1]==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readMatchFiles[1]);
		return FAILED;
	}

	// initialize the memory of read buffers
	if(initMemReadsBuf(&pReadBuf)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the read buffer, error!\n", __LINE__, __func__);
		return FAILED;
	}


	tmpReadBuf = (struct tmpReadNode *) calloc (MAX_READ_BUF_SIZE, sizeof(struct tmpReadNode));
	if(tmpReadBuf==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory for temporary reads buffer, error!\n", __LINE__, __func__);
		return FAILED;
	}
	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		tmpReadBuf[i].seq = (char *) calloc (readLen+1, sizeof(char));
		if(tmpReadBuf[i].seq==NULL)
		{
			printf("line=%d, In %s(), cannot allocate the memory for temporary reads, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	matchedReadNum[0] = matchedContigNum[0] = 0;
	matchedReadNum[1] = matchedContigNum[1] = 0;
	readID = 1;
	tmpReadNum = 0;

	for(fileID=0; fileID<readFileNum; fileID++)
	{
		fpPE = fopen(readFilesInput[fileID], "r");
		if(fpPE==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readFilesInput[fileID]);
			return FAILED;
		}

		//================== begin Mapping =================
		//get single read
		while(1)
		{
			// check the end of files
			if(feof(fpPE))
			{
				break;
			}

			// fille the reads to reads buffers
			if(fillReadsToBufFasta(fpPE, pReadBuf, &readsNum, readLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// map the reads in reads buffers
			for(i=0; i<readsNum; i+=2)
			{
				readSeq[0] = pReadBuf[i].seq;
				readSeq[1] = pReadBuf[i+1].seq;

				// ####################### Debug information #######################
				//if(readID==6947)
				//{
				//	printf("readID=%lu\n", readID);
				//}
				// ####################### Debug information #######################

				errorFlag[0] = errorFlag[1] = NO;

				// check unknown bases
				if(pReadBuf[i].len<readLen || containUnknownBase(readSeq[0])==YES)
				{
					errorFlag[0] = YES;
				}
				if(pReadBuf[i+1].len<readLen || containUnknownBase(readSeq[1])==YES)
				{
					errorFlag[1] = YES;
				}

				if((errorFlag[0]==YES) && (errorFlag[1]==YES))
				{
					readID += 2;
					continue;
				}

				matchFlag[0] = matchFlag[1] = NO;	// set the match flag of the read to NO
				if(errorFlag[0]==NO)
				{
					if(mapSingleRead(readID, readSeq[0], matchFlag, fpPE_Result[0], matchedContigNum, matchedReadNum)==FAILED)
					{
						printf("line=%d, In %s(), cannot map single read, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

				if(errorFlag[1]==NO)
				{
					if(mapSingleRead(readID+1, readSeq[1], matchFlag+1, fpPE_Result[1], matchedContigNum+1, matchedReadNum+1)==FAILED)
					{
						printf("line=%d, In %s(), cannot map single read, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

				if(matchFlag[0]==YES || matchFlag[1]==YES)
				{
					if(errorFlag[0]==NO)
					{
						//fprintf(fpReadSeq, "%lu\t%s\n", readID, readSeq[0]);

						tmpReadBuf[tmpReadNum].readID = readID;
						strcpy(tmpReadBuf[tmpReadNum].seq, readSeq[0]);
						tmpReadNum ++;

						if(tmpReadNum==MAX_READ_BUF_SIZE)
						{
							for(j=0; j<tmpReadNum; j++)
								fprintf(fpReadSeq, "%lu\t%s\n", tmpReadBuf[j].readID, tmpReadBuf[j].seq);
							tmpReadNum = 0;
						}
					}

					if(errorFlag[1]==NO)
					{
						//fprintf(fpReadSeq, "%lu\t%s\n", readID+1, readSeq[1]);

						tmpReadBuf[tmpReadNum].readID = readID+1;
						strcpy(tmpReadBuf[tmpReadNum].seq, readSeq[1]);
						tmpReadNum ++;

						if(tmpReadNum==MAX_READ_BUF_SIZE)
						{
							for(j=0; j<tmpReadNum; j++)
								fprintf(fpReadSeq, "%lu\t%s\n", tmpReadBuf[j].readID, tmpReadBuf[j].seq);
							tmpReadNum = 0;
						}
					}
				}

				readID += 2;
			}
		}

		fclose(fpPE);
		fpPE = NULL;
	}

	// process the remained temporary reads
	if(tmpReadNum>0)
	{
		for(j=0; j<tmpReadNum; j++)
			fprintf(fpReadSeq, "%lu\t%s\n", tmpReadBuf[j].readID, tmpReadBuf[j].seq);
		tmpReadNum = 0;
	}

	// free memory of read buffers
	freeMemReadsBuf(&pReadBuf);

	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		free(tmpReadBuf[i].seq);
		tmpReadBuf[i].seq = NULL;
	}
	free(tmpReadBuf);
	tmpReadBuf = NULL;


	fclose(fpReadSeq);
	fpReadSeq = NULL;
	fclose(fpPE_Result[0]);
	fpPE_Result[0] = NULL;
	fclose(fpPE_Result[1]);
	fpPE_Result[1] = NULL;


#if DEBUG_FLAG
	printf("matchedReadNum[0]=%lu, matchedContigNum[0]=%lu\n", matchedReadNum[0], matchedContigNum[0]);
	printf("matchedReadNum[1]=%lu, matchedContigNum[1]=%lu\n", matchedReadNum[1], matchedContigNum[1]);
	printf("total reads: %lu\n", readID-1);
#endif

	return SUCCESSFUL;
}

/**
 * Map the paired reads in single files, and save results to file with a given name.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 *
 *  Format:
 *		Each matched information item is saved as a readMatchInfoTemp structure.
 */
short mapPEFastqInterleaved(const char *readSeqFile, const char *readMatchFile1, const char *readMatchFile2, char **readFilesInput, int readFileNum)
{
	FILE *fpPE, *fpPE_Result[2], *fpReadSeq;
	uint64_t matchedReadNum[2], matchedContigNum[2];
	int matchFlag[2];
	uint64_t i, j, fileID, readID;		// readID starts from 1

	readBuf_t *pReadBuf;
	uint64_t readsNum;

	short errorFlag[2];
	char *readSeq[2];

	struct tmpReadNode
	{
		uint64_t readID;
		char *seq;
	};

	struct tmpReadNode *tmpReadBuf;
	uint64_t tmpReadNum;

	fpReadSeq = fopen(readSeqFile, "w");
	if(fpReadSeq==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readSeqFile);
		return FAILED;
	}


	fpPE_Result[0] = fopen(readMatchFiles[0], "w");
	if(fpPE_Result[0]==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readMatchFiles[0]);
		return FAILED;
	}
	fpPE_Result[1] = fopen(readMatchFiles[1], "w");
	if(fpPE_Result[1]==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readMatchFiles[1]);
		return FAILED;
	}

	// initialize the memory of read buffers
	if(initMemReadsBuf(&pReadBuf)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the read buffer, error!\n", __LINE__, __func__);
		return FAILED;
	}


	tmpReadBuf = (struct tmpReadNode *) calloc (MAX_READ_BUF_SIZE, sizeof(struct tmpReadNode));
	if(tmpReadBuf==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory for temporary reads buffer, error!\n", __LINE__, __func__);
		return FAILED;
	}
	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		tmpReadBuf[i].seq = (char *) calloc (readLen+1, sizeof(char));
		if(tmpReadBuf[i].seq==NULL)
		{
			printf("line=%d, In %s(), cannot allocate the memory for temporary reads, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	matchedReadNum[0] = matchedContigNum[0] = 0;
	matchedReadNum[1] = matchedContigNum[1] = 0;
	readID = 1;
	tmpReadNum = 0;

	for(fileID=0; fileID<readFileNum; fileID++)
	{
		fpPE = fopen(readFilesInput[fileID], "r");
		if(fpPE==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readFilesInput[fileID]);
			return FAILED;
		}

		//================== begin Mapping =================
		//get single read
		while(1)
		{
			// check the end of files
			if(feof(fpPE))
			{
				break;
			}

			// fille the reads to reads buffers
			if(fillReadsToBufFastq(fpPE, pReadBuf, &readsNum, readLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// map the reads in reads buffers
			for(i=0; i<readsNum; i+=2)
			{
				readSeq[0] = pReadBuf[i].seq;
				readSeq[1] = pReadBuf[i+1].seq;

				// ####################### Debug information #######################
				//if(readID==6947)
				//{
				//	printf("readID=%lu\n", readID);
				//}
				// ####################### Debug information #######################

				errorFlag[0] = errorFlag[1] = NO;

				// check unknown bases
				if(pReadBuf[i].len<readLen || containUnknownBase(readSeq[0])==YES)
				{
					errorFlag[0] = YES;
				}
				if(pReadBuf[i+1].len<readLen || containUnknownBase(readSeq[1])==YES)
				{
					errorFlag[1] = YES;
				}

				if((errorFlag[0]==YES) && (errorFlag[1]==YES))
				{
					readID += 2;
					continue;
				}

				matchFlag[0] = matchFlag[1] = NO;	// set the match flag of the read to NO
				if(errorFlag[0]==NO)
				{
					if(mapSingleRead(readID, readSeq[0], matchFlag, fpPE_Result[0], matchedContigNum, matchedReadNum)==FAILED)
					{
						printf("line=%d, In %s(), cannot map single read, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

				if(errorFlag[1]==NO)
				{
					if(mapSingleRead(readID+1, readSeq[1], matchFlag+1, fpPE_Result[1], matchedContigNum+1, matchedReadNum+1)==FAILED)
					{
						printf("line=%d, In %s(), cannot map single read, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

				if(matchFlag[0]==YES || matchFlag[1]==YES)
				{
					if(errorFlag[0]==NO)
					{
						//fprintf(fpReadSeq, "%lu\t%s\n", readID, readSeq[0]);

						tmpReadBuf[tmpReadNum].readID = readID;
						strcpy(tmpReadBuf[tmpReadNum].seq, readSeq[0]);
						tmpReadNum ++;

						if(tmpReadNum==MAX_READ_BUF_SIZE)
						{
							for(j=0; j<tmpReadNum; j++)
								fprintf(fpReadSeq, "%lu\t%s\n", tmpReadBuf[j].readID, tmpReadBuf[j].seq);
							tmpReadNum = 0;
						}
					}

					if(errorFlag[1]==NO)
					{
						//fprintf(fpReadSeq, "%lu\t%s\n", readID+1, readSeq[1]);

						tmpReadBuf[tmpReadNum].readID = readID+1;
						strcpy(tmpReadBuf[tmpReadNum].seq, readSeq[1]);
						tmpReadNum ++;

						if(tmpReadNum==MAX_READ_BUF_SIZE)
						{
							for(j=0; j<tmpReadNum; j++)
								fprintf(fpReadSeq, "%lu\t%s\n", tmpReadBuf[j].readID, tmpReadBuf[j].seq);
							tmpReadNum = 0;
						}
					}
				}

				readID += 2;
			}
		}

		fclose(fpPE);
		fpPE = NULL;
	}

	// process the remained temporary reads
	if(tmpReadNum>0)
	{
		for(j=0; j<tmpReadNum; j++)
			fprintf(fpReadSeq, "%lu\t%s\n", tmpReadBuf[j].readID, tmpReadBuf[j].seq);
		tmpReadNum = 0;
	}

	// free memory of read buffers
	freeMemReadsBuf(&pReadBuf);

	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		free(tmpReadBuf[i].seq);
		tmpReadBuf[i].seq = NULL;
	}
	free(tmpReadBuf);
	tmpReadBuf = NULL;


	fclose(fpReadSeq);
	fpReadSeq = NULL;
	fclose(fpPE_Result[0]);
	fpPE_Result[0] = NULL;
	fclose(fpPE_Result[1]);
	fpPE_Result[1] = NULL;


#if DEBUG_FLAG
	printf("matchedReadNum[0]=%lu, matchedContigNum[0]=%lu\n", matchedReadNum[0], matchedContigNum[0]);
	printf("matchedReadNum[1]=%lu, matchedContigNum[1]=%lu\n", matchedReadNum[1], matchedContigNum[1]);
	printf("total reads: %lu\n", readID-1);
#endif

	return SUCCESSFUL;
}

/**
 * Initialize the memory of reads buffer.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemReadsBuf(readBuf_t **pBuf)
{
	uint64_t i;

	*pBuf = (readBuf_t*) calloc(MAX_READ_BUF_SIZE, sizeof(readBuf_t));
	if(*pBuf==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		(*pBuf)[i].seq = (char*) calloc(readLen+1, sizeof(char));
		if((*pBuf)[i].seq==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Release memory of reads buffer.
 */
void freeMemReadsBuf(readBuf_t **pBuf)
{
	uint64_t i;
	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		free((*pBuf)[i].seq);
		(*pBuf)[i].seq = NULL;
	}

	free(*pBuf);
	*pBuf = NULL;
}

/**
 * Fill reads to reads buffer.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillReadsToBufFasta(FILE *fpReads, readBuf_t *pBuf, uint64_t *readsNum, int readLenThreshold)
{
	uint64_t i;
	int returnCode;

	*readsNum = 0;
	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		returnCode = getSingleReadFasta(fpReads, pBuf+i, readLenThreshold);
		if(returnCode==FAILED)
			break;
		else if(returnCode==ERROR)
		{
			printf("line=%d, In %s(), cannot get single read, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	*readsNum = i;

	return SUCCESSFUL;
}

/**
 * Fill reads to reads buffer.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillReadsToBufFastq(FILE *fpReads, readBuf_t *pBuf, uint64_t *readsNum, int readLenThreshold)
{
	uint64_t i;
	int returnCode;

	*readsNum = 0;
	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		returnCode = getSingleReadFastq(fpReads, pBuf+i, readLenThreshold);
		if(returnCode==FAILED)
			break;
		else if(returnCode==ERROR)
		{
			printf("line=%d, In %s(), cannot get single read, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	*readsNum = i;

	return SUCCESSFUL;
}

/**
 * map single read.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED;
 */
short mapSingleRead(uint64_t readID, char *readSeq, int *matchFlag, FILE *fpPE_Result, uint64_t *matchedContigNum, uint64_t *matchedReadNum)
{
	uint64_t *readSeqInt, *reverse_readSeqInt;
	int hitRow, contigItemNum, j;
	contigMatchInfo *pContigMatchInfo;
	readMatchInfoTemp tmp_matchInfo;

	readSeqInt = (uint64_t *) calloc(entriesPerRead, sizeof(uint64_t));
	if(readSeqInt==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	reverse_readSeqInt = (uint64_t *) calloc(entriesPerRead, sizeof(uint64_t));
	if(reverse_readSeqInt==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	*matchFlag = NO;

	if(seqHash(readSeqInt, readSeq, 1)==FAILED)
	{
		printf("line=%d, In %s(), cannot hash the sequence [ %s ], error!\n", __LINE__, __func__, readSeq);
		return FAILED;
	}

	hitRow = seqSearch(readSeqInt);
	if(hitRow>=0)
	{
		*matchFlag = YES;

		pContigMatchInfo = contigMatchInfoArr + seqRowArr[hitRow].startRow;
		contigItemNum = seqRowArr[hitRow].rowsNum;

		for(j=0; j<contigItemNum; j++)
		{
			//fprintf(fpPE_Result, "%lu\t%u\t%u\t%u\t%c\n", readID, pContigInfo[j].contigID, pContigInfo[j].contigPos, pContigInfo[j].contigEnd, ORIENTATION_PLUS);

			tmp_matchInfo.readID = readID;
			tmp_matchInfo.contigID = pContigMatchInfo[j].contigID;
			tmp_matchInfo.contigPos = pContigMatchInfo[j].contigPos;
			tmp_matchInfo.orientation = ORIENTATION_PLUS;
			tmp_matchInfo.contigEnd = pContigMatchInfo[j].contigEnd;
			if(fwrite(&tmp_matchInfo, sizeof(readMatchInfoTemp), 1, fpPE_Result)!=1)
			{
				printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		(*matchedContigNum) += contigItemNum;
	}

	// map the reverse complements of the read
	if(getReverseComplements(reverse_readSeqInt, readSeqInt)==FAILED)
	{
		printf("line=%d, In %s(), cannot reverse the sequence [ %s ], error!\n", __LINE__, __func__, readSeq);
		return FAILED;
	}

	hitRow = seqSearch(reverse_readSeqInt);
	if(hitRow>=0)
	{
		*matchFlag = YES;

		pContigMatchInfo = contigMatchInfoArr + seqRowArr[hitRow].startRow;
		contigItemNum = seqRowArr[hitRow].rowsNum;

		for(j=0; j<contigItemNum; j++)
		{
			//fprintf(fpPE_Result, "%lu\t%u\t%u\t%u\t%c\n", readID, pContigInfo[j].contigID, pContigInfo[j].contigPos, pContigInfo[j].contigEnd, ORIENTATION_MINUS);

			tmp_matchInfo.readID = readID;
			tmp_matchInfo.contigID = pContigMatchInfo[j].contigID;
			tmp_matchInfo.contigPos = pContigMatchInfo[j].contigPos;
			tmp_matchInfo.orientation = ORIENTATION_MINUS;
			tmp_matchInfo.contigEnd = pContigMatchInfo[j].contigEnd;
			if(fwrite(&tmp_matchInfo, sizeof(readMatchInfoTemp), 1, fpPE_Result)!=1)
			{
				printf("line=%d, In %s(), fwrite error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		(*matchedContigNum) += contigItemNum;
	}
	if(*matchFlag==YES)
	{
		(*matchedReadNum) ++;
	}

	free(readSeqInt);
	readSeqInt = NULL;

	free(reverse_readSeqInt);
	reverse_readSeqInt = NULL;

	return SUCCESSFUL;
}

/**
 * get the single read from the given file.
 *  @return:
 *  	If succeed, return SUCCESSFUL;
 *  	else if the file end is reached, return FAILED;
 *  	otherwise, return ERROR.
 */
short getSingleReadFasta(FILE *fpPE, readBuf_t *pBuf, int readLenThreshold)
{
	int tmpLen;
	char ch, *readSeq;

	readSeq = pBuf->seq;

	ch = fgetc(fpPE);
	if(feof(fpPE))// the file end is reached.
	{
		return FAILED;
	}

	while(!feof(fpPE))
	{
		ch = fgetc(fpPE);
		while(ch!='\n') ch = fgetc(fpPE);

		tmpLen = 0;
		while(ch!='>' && ch!=-1)
		{
			if(ch!='\n')
			{
				if(tmpLen<readLenThreshold)
					readSeq[tmpLen++] = ch;
			}
			ch = fgetc(fpPE);
		}
		readSeq[tmpLen] = '\0';

		if(ch=='>')  //a read is read finished
		{
			pBuf->len = tmpLen;
			break;
		}
	}


	return SUCCESSFUL;
}

/**
 * get the single read from the given file.
 *  @return:
 *  	If succeed, return SUCCESSFUL;
 *  	else if the file end is reached, return FAILED;
 *  	otherwise, return ERROR.
 */
short getSingleReadFastq(FILE *fpPE, readBuf_t *pBuf, int readLenThreshold)
{
	unsigned int tmpLen, line_index;
	char *readSeq;

	readSeq = pBuf->seq;

	char ch = fgetc(fpPE);
	if(feof(fpPE))// the file end is reached.
	{
		return FAILED;
	}

	line_index = 0;
	while(!feof(fpPE))
	{
		if(line_index==0)  //the sequence name line
		{
			ch = fgetc(fpPE);
			while(ch!='\n' && ch!=-1)
				ch = fgetc(fpPE);
		}else if(line_index==1)  //the sequence line
		{
			tmpLen = 0;
			ch = fgetc(fpPE);
			while(ch!='\n')
			{
				if(tmpLen<readLenThreshold)
					readSeq[tmpLen++] = ch;
				ch = fgetc(fpPE);
			}
			readSeq[tmpLen] = '\0';
		}else if(line_index==2)  //the sequence name line
		{
			ch = fgetc(fpPE);
			while(ch!='\n')
				ch = fgetc(fpPE);
		}else
		{
			ch = fgetc(fpPE);
			while(ch!='\n'  && ch!=-1)
				ch = fgetc(fpPE);
		}
		line_index++;

		if(line_index==4)  //a read is read finished
		{
			pBuf->len = tmpLen;
			break;
		}
	}


	// check the return code
	if(line_index==4)
	{
		return SUCCESSFUL;
	}else
	{
		return ERROR;
	}
}

/**
 * Get the reverse complementary sequence integer of a read.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReverseComplements(uint64_t *reverse_readSeqInt, uint64_t *readSeqInt)
{
	int i, j, k, baseNum;

	for(k=0; k<entriesPerRead; k++)
	{
		reverse_readSeqInt[k] = 0;
	}

	// inverse the sequence
	i = 0;
	baseNum = 0;
	for(j=0; j<lastEntryBaseNumRead; j++)
	{
		reverse_readSeqInt[i] = (reverse_readSeqInt[i] << 2) | ((readSeqInt[entriesPerRead-1] >> (j<<1)) & 3);
		baseNum ++;
		if(baseNum==32)
		{
			i ++;
			baseNum = 0;
		}
	}
	for(k=entriesPerRead-2; k>=0; k--)
	{
		for(j=0; j<32; j++)
		{
			reverse_readSeqInt[i] = (reverse_readSeqInt[i] << 2) | ((readSeqInt[k] >> (j<<1)) & 3);
			baseNum ++;
			if(baseNum==32)
			{
				i ++;
				baseNum = 0;
			}
		}
	}

	// complement the sequence
	for(k=0; k<entriesPerRead-1; k++)
	{
		reverse_readSeqInt[k] = ~(reverse_readSeqInt[k]);
	}
	reverse_readSeqInt[entriesPerRead-1] = (~(reverse_readSeqInt[entriesPerRead-1])) & lastEntryMaskRead;

	return SUCCESSFUL;
}


/**
 * Get the match result file name.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short getResultFileName(char *resultFileName, const char *readFileName)
{
	int len, splitPos;
	len = strlen(readFileName);

	if(readFileName[len-1]=='/')
	{
		printf("line=%d, In %s(), error read file name [ %s ]!\n", __LINE__, __func__, readFileName);
		return FAILED;
	}

	for(splitPos=len-1; splitPos>=0; splitPos--)
	{
		if(readFileName[splitPos]=='/')
		{
			break;
		}
	}

	strcpy(resultFileName, outputPathStr);
	strcat(resultFileName, readFileName+splitPos+1);
	strcat(resultFileName, ".match.bin");

	return SUCCESSFUL;
}
