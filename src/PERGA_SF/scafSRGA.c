/*
 * scafSRGA.c
 *
 *  Created on: Mar 31, 2011
 *      Author: zhuxiao
 */

#include "inc/scafStdinc.h"
#include "inc/scafExtvab.h"


short startSRGA_SF(char *contigsFilePara, int alignRegSizePara, int minContigLenPara, int pairedModePara, char **readFilesPara, int readFileNumPara, double meanSizeInsertPara, double standardDevPara, char *gapFillFlagPara, char *outputPathPara, char *outputPrefixPara)
{
	struct timeval tp_start,tp_end;
	double time_used;
	unsigned long start_cpu, end_cpu; // CPU time

	gettimeofday(&tp_start,NULL);

	//set the start time: CPU time
	start_cpu = clock();


	if(setGlobalParas(outputPathPara, outputPrefixPara, contigsFilePara, minContigLenPara, alignRegSizePara, readFilesPara, readFileNumPara, pairedModePara, meanSizeInsertPara, standardDevPara, gapFillFlagPara)==FAILED)
	//if(setGlobalParas("../output/", "../input/contigs.fa", "../input/reads/SRR001665_1.fastq", "../input/reads/SRR001665_2.fastq")==FAILED)
	//if(setGlobalParas("/home/plzeng/longlong/test/SRGA/output/", "/home/plzeng/longlong/test/SRGA/contigs.fa", "/home/plzeng/longlong/data/P.sringae/ERR005143_1.fastq", "/home/plzeng/longlong/data/P.sringae/ERR005143_2.fastq")==FAILED)
	//if(setGlobalParas("/home/zx/study/scaffolding/P.sringae/output/", "/home/zx/study/scaffolding/P.sringae/input/contigs.fa", "/home/zx/study/scaffolding/P.sringae/input/reads/ERR005143_1.fastq", "/home/zx/study/scaffolding/P.sringae/input/reads/ERR005143_2.fastq")==FAILED)
	//if(setGlobalParas("../output/", "../contigs.fa", "/home/plzeng/longlong/data/P.sringae/ERR005143_1.fastq", "/home/plzeng/longlong/data/P.sringae/ERR005143_2.fastq")==FAILED)
	{
		printf("line=%d, In %s(), cannot set global variables, error!\n", __LINE__, __func__);
		return FAILED;
	}


	if(startScaffolding()==FAILED)
	{
		printf("line=%d, In %s(), cannot build scaffolds, error!\n", __LINE__, __func__);
		return FAILED;
	}


	freeGlobalParas();

	// set the end time: CPU time
	end_cpu = clock();

	// print the CPU time
	printf("Total CPU time used: %.2f seconds.\n", (double)(end_cpu-start_cpu)/CLOCKS_PER_SEC);

	// print the system time
	gettimeofday(&tp_end,NULL);
	time_used = tp_end.tv_sec-tp_start.tv_sec + (double)(tp_end.tv_usec-tp_start.tv_usec)/1000000;
	printf("Total Used Time By gettimeofday(): %.2f Seconds.\n", time_used);

	return 0;
}


/**
 * Set the global variables.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short setGlobalParas(const char *outputPathName, const char *outputPrefixName, const char *contigFileName, int minContigLenThreshold, int alignRegSizeThreshold, char **readFilesPara, int readFileNumPara, int pairedModePara, double meanSizeInsertPara, double standardDevPara, char *gapFillFlagPara)
{
	int i, prefixLen;

	if(setGlobalPath(outputPathName)==FAILED)
	{
		printf("line=%d, In %s(), cannot set global paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	strcpy(outputPrefix, outputPrefixName);

	if(strlen(contigFileName)>0)
	{
		if(setContigFileNames(outputPathStr, outputPrefix, contigFileName)==FAILED)
		{
			printf("line=%d, In %s(), cannot set contig file, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		printf("please specify contig file.\n");
		return FAILED;
	}

	readFileNum = readFileNumPara;
	for(i=0; i<readFileNum; i++)
	{
		readFilesInput[i] = (char *) calloc (256, sizeof(char));
		if(readFilesInput[i]==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		strcpy(readFilesInput[i], readFilesPara[i]);
	}


	if(strlen(gapFillFlagPara)==0)
	{
		gapFillFlag = YES;
	}else if(strcmp(gapFillFlagPara, "YES")==0 || strcmp(gapFillFlagPara, "yes")==0 || strcmp(gapFillFlagPara, "Y")==0 || strcmp(gapFillFlagPara, "y")==0)
	{
		gapFillFlag = YES;
	}else if(strcmp(gapFillFlagPara, "NO")==0 || strcmp(gapFillFlagPara, "no")==0 || strcmp(gapFillFlagPara, "N")==0 || strcmp(gapFillFlagPara, "n")==0)
		gapFillFlag = NO;
	else
	{
		printf("Exception: unknown gap filling flag: %s.\n", gapFillFlagPara);
		return FAILED;
	}

	if(pairedModePara==0)
		pairedMode = 1;
	else
		pairedMode = pairedModePara;

	if(pairedMode==1 && readFileNum%2==1)
	{
		printf("Exception: please specify the correct read files or correct paired mode.\n");
		return FAILED;
	}

	if(meanSizeInsertPara==0)
	{
		PEGivenType = NONE_PE_GIVEN_TYPE;
		meanSizeInsert = 0;
		stardardDeviationInsert = 0;
	}else if(meanSizeInsertPara!=0 && standardDevPara==0)
	{
		PEGivenType = INSERT_PE_GIVEN_TYPE;
		meanSizeInsert = meanSizeInsertPara;
		stardardDeviationInsert = meanSizeInsertPara * DRAFT_SDEV_FACTOR;
	}else if(meanSizeInsertPara!=0 && standardDevPara!=0)
	{
		PEGivenType = BOTH_PE_GIVEN_TYPE;
		meanSizeInsert = meanSizeInsertPara;
		stardardDeviationInsert = standardDevPara;
	}

	prefixLen = strlen(outputPrefix);

	// get the match result file name
	strcpy(readMatchFiles[0], outputPathStr);
	if(prefixLen>0)
	{
		strcat(readMatchFiles[0], outputPrefix);
		strcat(readMatchFiles[0], "_");
	}
	strcat(readMatchFiles[0], "reads_match_1.bin");
	strcpy(readMatchFiles[1], outputPathStr);
	if(prefixLen>0)
	{
		strcat(readMatchFiles[1], outputPrefix);
		strcat(readMatchFiles[1], "_");
	}
	strcat(readMatchFiles[1], "reads_match_2.bin");

	// readSeqFile
	strcpy(readSeqFile, outputPathStr);
	if(prefixLen>0)
	{
		strcat(readSeqFile, outputPrefix);
		strcat(readSeqFile, "_");
	}
	strcat(readSeqFile, "readSeq.txt");

	// monoReadSeqFile
	strcpy(monoReadSeqFile, outputPathStr);
	if(prefixLen>0)
	{
		strcat(monoReadSeqFile, outputPrefix);
		strcat(monoReadSeqFile, "_");
	}
	strcat(monoReadSeqFile, "monoReadSeq.txt");

	// contigIndexFile
	strcpy(contigIndexFile, outputPathStr);
	if(prefixLen>0)
	{
		strcat(contigIndexFile, outputPrefix);
		strcat(contigIndexFile, "_");
	}
	strcat(contigIndexFile, "contigIndex.bin");

	// readListFiles[4]: //[0]-- RL1, [1]-- RL2, [2]-- SRL, [3]-- MRL
	strcpy(readListFiles[0], outputPathStr);
	if(prefixLen>0)
	{
		strcat(readListFiles[0], outputPrefix);
		strcat(readListFiles[0], "_");
	}
	strcat(readListFiles[0], "readList1.bin");
	strcpy(readListFiles[1], outputPathStr);
	if(prefixLen>0)
	{
		strcat(readListFiles[1], outputPrefix);
		strcat(readListFiles[1], "_");
	}
	strcat(readListFiles[1], "readList2.bin");
	strcpy(readListFiles[2], outputPathStr);
	if(prefixLen>0)
	{
		strcat(readListFiles[2], outputPrefix);
		strcat(readListFiles[2], "_");
	}
	strcat(readListFiles[2], "shared_RL.bin");
	strcpy(readListFiles[3], outputPathStr);
	if(prefixLen>0)
	{
		strcat(readListFiles[3], outputPrefix);
		strcat(readListFiles[3], "_");
	}
	strcat(readListFiles[3], "monoReadList.bin");

	// readListFilesAfterOverlap[4]: //[0]-- RL1, [1]-- RL2, [2]-- SRL, [3]-- MRL
	strcpy(readListFilesAfterOverlap[0], readListFiles[0]);
	strcat(readListFilesAfterOverlap[0], ".new");
	strcpy(readListFilesAfterOverlap[1], readListFiles[1]);
	strcat(readListFilesAfterOverlap[1], ".new");
	strcpy(readListFilesAfterOverlap[2], readListFiles[2]);
	strcat(readListFilesAfterOverlap[2], ".new");
	strcpy(readListFilesAfterOverlap[3], readListFiles[3]);
	strcat(readListFilesAfterOverlap[3], ".new");

	// contigListFile
	strcpy(contigListFile, outputPathStr);
	if(prefixLen>0)
	{
		strcat(contigListFile, outputPrefix);
		strcat(contigListFile, "_");
	}
	strcat(contigListFile, "contigList.bin");

	// contigListFileAfterOverlap
	strcpy(contigListFileAfterOverlap, contigListFile);
	strcat(contigListFileAfterOverlap, ".new");

	// meanSdevFile
	strcpy(meanSdevFile, outputPathStr);
	if(prefixLen>0)
	{
		strcat(meanSdevFile, outputPrefix);
		strcat(meanSdevFile, "_");
	}
	strcat(meanSdevFile, "meanSdev");

	// contigOverlapFile
	strcpy(contigOverlapFile, outputPathStr);
	if(prefixLen>0)
	{
		strcat(contigOverlapFile, outputPrefix);
		strcat(contigOverlapFile, "_");
	}
	strcat(contigOverlapFile, "contigOverlap.txt");

	// contigOverlapFile after gap filling
	if(gapFillFlag==YES)
	{
		strcpy(contigOverlapFileAfterFapFilling, contigOverlapFile);
		strcat(contigOverlapFileAfterFapFilling, ".new");
	}

	// linkInfoFile
	strcpy(linkInfoFile, outputPathStr);
	if(prefixLen>0)
	{
		strcat(linkInfoFile, outputPrefix);
		strcat(linkInfoFile, "_");
	}
	strcat(linkInfoFile, "contigLinkInfo.txt");

	// averLinkNumFile
	strcpy(averLinkNumFile, outputPathStr);
	if(prefixLen>0)
	{
		strcat(averLinkNumFile, outputPrefix);
		strcat(averLinkNumFile, "_");
	}
	strcat(averLinkNumFile, "averLinkNum.txt");

	// scafGraphFile
	strcpy(scafGraphFile, outputPathStr);
	if(prefixLen>0)
	{
		strcat(scafGraphFile, outputPrefix);
		strcat(scafGraphFile, "_");
	}
	strcat(scafGraphFile, "scafGraph.bin");

	// scafSeqFile
	strcpy(scafSeqFile, outputPathStr);
	if(prefixLen>0)
	{
		strcat(scafSeqFile, outputPrefix);
		strcat(scafSeqFile, "_");
	}
	strcat(scafSeqFile, "scaffolds.fa");


	// get the reads file format
	if(getReadsFileFormatInScaf(&readsFileFormatType, readFilesInput, readFileNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot get read file format, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the readLen from fastq file
	if(getReadLenFromFilesInScaf(&readLen, &averReadLen, readFilesInput, readFileNum, readsFileFormatType)==FAILED)
	{
		printf("line=%d, In %s(), cannot get readLen, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(readLen>MAX_READ_LEN)
		readLen = MAX_READ_LEN;

	//initialize the variables
	entriesPerRead = ((readLen-1)/32) + 1;
	if(readLen%32==0)
	{
		lastEntryMaskRead = (uint64_t) -1;
		lastEntryBaseNumRead = 32;
	}else
	{
		lastEntryMaskRead = (1LLU << ((readLen%32)<<1)) - 1;
		lastEntryBaseNumRead = readLen % 32;
	}

	if(alignRegSizeThreshold>0)
	{
		contigEndLen = alignRegSizeThreshold;
	}else
	{
		contigEndLen = CONTIG_END_LEN;
		if(contigEndLen<readLen*CONTIG_END_LEN_FACTOR)
		{
			contigEndLen = readLen*CONTIG_END_LEN_FACTOR;
			printf("The minContigLenThres will be set to %d instead.\n", contigEndLen);
		}
	}

	if(minContigLenThreshold>0)
	{
		minContigLenThres = minContigLenThreshold;
	}else
	{
		minContigLenThres = MIN_CONTIG_LEN_THRES;
		if(minContigLenThres<readLen)
		{
			minContigLenThres = readLen;
			printf("The minContigLenThres will be set to %d instead.\n", minContigLenThres);
		}
	}


	// output the variables
	printf("\ncontigs file       : %s\n", contigsFile);
	printf("align region size  : %d\n", contigEndLen);
	printf("min contig size    : %d\n", minContigLenThres);
	printf("paired mode        : %d\n", pairedMode);
	for(i=0; i<readFileNumPara; i++)
		printf("read files[%d]      : %s\n", i, readFilesInput[i]);
	if(meanSizeInsert>0)
	{
		printf("insert size:       : %.2f\n", meanSizeInsert);
		printf("insert size sdev.  : %.2f\n", stardardDeviationInsert);
	}
	if(gapFillFlag==YES)
		printf("gap filling flag   : yes\n");
	else
		printf("gap filling flag   : no\n");
	printf("output directory   : %s\n\n", outputPathStr);
	//printf("output file prefix : %s\n\n", outputPrefixPara);


	//######################### Debug information #########################
	//printf("readLen=%d, entriesPerRead=%d, lastEntryBaseNumRead=%d, lastEntryMaskRead=%lX\n", readLen, entriesPerRead, lastEntryBaseNumRead, lastEntryMaskRead);
	//######################### Debug information #########################


	return SUCCESSFUL;
}

/**
 * Free the global parameters.
 */
void freeGlobalParas()
{
	int i;
	for(i=0; i<readFileNum; i++) { free(readFilesInput[i]); readFilesInput[i] = NULL; }

	if(DELETE_FILES)
	{
		// delete files
		remove(averLinkNumFile);
		remove(contigIndexFile);
		remove(linkInfoFile);
		remove(contigListFile);
		remove(contigListFileAfterOverlap);
		remove(contigOverlapFile);
		remove(contigOverlapFileAfterFapFilling);
		remove(meanSdevFile);
		remove(monoReadSeqFile);

		remove(contigsFileAfterOverlap);
		remove(contigsFileAfterGapFilling);

		for(i=0; i<4; i++) remove(readListFiles[i]);
		for(i=0; i<3; i++) remove(readListFilesAfterOverlap[i]);

		remove(readSeqFile);

		for(i=0; i<2; i++) remove(readMatchFiles[i]);

		remove(scafGraphFile);
	}

}

/**
 * Set the global input and output path directory.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short setGlobalPath(const char *outPathStr)
{
	int inPathLen, outPathLen;
	struct stat st;

	strcpy(outputPathStr, outPathStr);

	inPathLen = strlen(inputPathStr);
	if(inPathLen>=255)
	{
		printf("line=%d, In %s(), input path length=%d, error!\n", __LINE__, __func__, inPathLen);
		return FAILED;
	}
	outPathLen = strlen(outputPathStr);
	if(outPathLen>=255)
	{
		printf("line=%d, In %s(), output path length=%d, error!\n", __LINE__, __func__, outPathLen);
		return FAILED;
	}

	if(inputPathStr[inPathLen-1]!='/')
	{
		inputPathStr[inPathLen] = '/';
		inputPathStr[inPathLen+1] = '\0';
	}

	if(outPathLen>0)
	{
		if(outputPathStr[outPathLen-1]!='/')
		{
			outputPathStr[outPathLen] = '/';
			outputPathStr[outPathLen+1] = '\0';
		}

		if(stat(outPathStr, &st)==-1)
		{
			if(mkdir(outPathStr, 0755)==-1)
			{
				printf("line=%d, In %s(), cannot create directory [ %s ], error!\n", __LINE__, __func__, outPathStr);
				return FAILED;
			}
		}
	}else
	{
		strcpy(outputPathStr, "./");
	}

	return SUCCESSFUL;
}

/**
 * Set the new contig file names.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short setContigFileNames(const char *outputPathStr, const char *outputPrefixName, const char *contigFileName)
{
	int i, pathLen, prefixLen;

	prefixLen = strlen(outputPrefixName);

	pathLen = strlen(contigFileName);
	for(i=pathLen-1; i>=0; i--)
	{
		if(contigFileName[i]=='/')
			break;
	}
	i ++;

	strcpy(contigsFile, contigFileName);

	strcpy(contigsFileFiltered, outputPathStr);
	if(prefixLen>0)
	{
		strcat(contigsFileFiltered, outputPrefixName);
		strcat(contigsFileFiltered, "_");
	}
	strcat(contigsFileFiltered, contigFileName+i);
	strcat(contigsFileFiltered, ".filtered");

	strcpy(contigsFileUltraShort, outputPathStr);
	if(prefixLen>0)
	{
		strcat(contigsFileUltraShort, outputPrefixName);
		strcat(contigsFileUltraShort, "_");
	}
	strcat(contigsFileUltraShort, contigFileName+i);
	strcat(contigsFileUltraShort, ".ultrashort");

	strcpy(contigsFileAfterOverlap, contigsFileFiltered);
	strcat(contigsFileAfterOverlap, ".new");
	strcpy(contigsFileAfterGapFilling, contigsFileAfterOverlap);
	strcat(contigsFileAfterGapFilling, ".new");

	return SUCCESSFUL;
}

/**
 * Get reads file format.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadsFileFormatInScaf(int *readsFileFormatType, char **readFilesInput, int readFileNum)
{
	int i, readFormat;
	FILE *fpRead;
	char ch;

	readFormat = -1;
	for(i=0; i<readFileNum; i++)
	{
		fpRead = fopen(readFilesInput[i], "r");
		if(fpRead==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readFilesInput[i]);
			return FAILED;
		}

		ch = fgetc(fpRead);
		if(ch=='>')
		{
			if(readFormat!=-1 && readFormat!=FILE_FORMAT_FASTA)
			{
				printf("Reads files cannot be in multiple format comblined together, error!\n");
				return FAILED;
			}

			readFormat = FILE_FORMAT_FASTA;
		}else if(ch=='@')
		{
			if(readFormat!=-1 && readFormat!=FILE_FORMAT_FASTQ)
			{
				printf("Reads files cannot be in multiple format combined together, error!\n");
				return FAILED;
			}

			readFormat = FILE_FORMAT_FASTQ;
		}else
		{
			printf("Reads files must be in fasta or fastq format, error!\n");
			return FAILED;
		}

		fclose(fpRead);
	}

	if(readFormat!=-1)
		*readsFileFormatType = readFormat;
	else
	{
		*readsFileFormatType = -1;
		printf("Reads files must be in fasta or fastq format, error!\n");
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get the read length from read files.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadLenFromFilesInScaf(int *readLen, int *averReadLen, char **readFilesInput, int readFileNum, int readsFileFormatType)
{
	if(readsFileFormatType==FILE_FORMAT_FASTA)
	{
		if(getMinReadLenFromFastaFilesInScaf(readLen, averReadLen, readFilesInput, readFileNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot get readlenInFile, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else if(readsFileFormatType==FILE_FORMAT_FASTQ)
	{
		if(getMinReadLenFromFastqFilesInScaf(readLen, averReadLen, readFilesInput, readFileNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot get readlenInFile, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the read length from fasta file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMinReadLenFromFastaFilesInScaf(int *readLenInFile, int *averReadLen, char **readFilesInput, int readFileNum)
{
	int i, tmpReadLen, sumReadLenTmp, sumNumTmp;

	sumReadLenTmp = 0;
	sumNumTmp = 0;
	*readLenInFile = INT_MAX;
	for(i=0; i<readFileNum; i++)
	{
		if(getReadLenFromFastaInScaf(&tmpReadLen, readFilesInput[i])==FAILED)
		{
			printf("line=%d, In %s(), cannot get readlenInFile, error!\n", __LINE__, __func__);
			return FAILED;
		}

#if (DEBUG_FLAG==YES)
		printf("tmpReadLen=%d\n", tmpReadLen);
#endif

		if(tmpReadLen<*readLenInFile)
			*readLenInFile = tmpReadLen;

		if(tmpReadLen>0)
		{
			sumReadLenTmp += tmpReadLen;
			sumNumTmp ++;
		}
	}

	if(sumNumTmp>0)
		*averReadLen = (double)sumReadLenTmp / sumNumTmp;
	else
		*averReadLen = 0;

#if (DEBUG_FLAG==YES)
		printf("readLen=%d\n", *readLenInFile);
		printf("averReadLen=%d\n", *averReadLen);
#endif

	return SUCCESSFUL;
}

/**
 * Get the read length from fastq file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMinReadLenFromFastqFilesInScaf(int *readLenInFile, int *averReadLen, char **readFilesInput, int readFileNum)
{
	int i, tmpReadLen, sumReadLenTmp, sumNumTmp;

	sumReadLenTmp = 0;
	sumNumTmp = 0;
	*readLenInFile = INT_MAX;
	for(i=0; i<readFileNum; i++)
	{
		if(getReadLenFromFastqInScaf(&tmpReadLen, readFilesInput[i])==FAILED)
		{
			printf("line=%d, In %s(), cannot get readlenInFile, error!\n", __LINE__, __func__);
			return FAILED;
		}

#if (DEBUG_FLAG==YES)
		printf("tmpReadLen=%d\n", tmpReadLen);
#endif

		if(tmpReadLen<*readLenInFile)
			*readLenInFile = tmpReadLen;

		if(tmpReadLen>0)
		{
			sumReadLenTmp += tmpReadLen;
			sumNumTmp ++;
		}
	}

	if(sumNumTmp>0)
		*averReadLen = (double)sumReadLenTmp / sumNumTmp;
	else
		*averReadLen = 0;

#if (DEBUG_FLAG==YES)
		printf("readLen=%d\n", *readLenInFile);
		printf("averReadLen=%d\n", *averReadLen);
#endif

	return SUCCESSFUL;
}

/**
 * Get the read length from fastq file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadLenFromFastaInScaf(int *tmpReadLen, const char *fastaFile)
{
	FILE *fpFasta;
	char ch, seq_data[5000];
	int64_t tmpReadsNum, tmpLen, sumLen;

	fpFasta = fopen(fastaFile, "r");
	if(fpFasta==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fastaFile);
		return FAILED;
	}

	tmpReadsNum = 0;
	sumLen = 0;
	while(!feof(fpFasta))
	{
		ch = fgetc(fpFasta);
		while(ch!='\n') ch = fgetc(fpFasta);

		tmpLen = 0;
		while(ch!='>' && ch!=-1)
		{
			if(ch!='\n')
			{
				seq_data[tmpLen++] = ch;
			}
			ch = fgetc(fpFasta);
		}
		seq_data[tmpLen] = '\0';

		if(ch=='>')  //a read is read finished
		{
			sumLen += tmpLen;

			tmpReadsNum ++;
			if(tmpReadsNum>=MAX_READ_NUM_READ_LEN_SAMPLE)
				break;
		}
	}

	if(tmpReadsNum>0)
	{
		*tmpReadLen = sumLen / tmpReadsNum;
	}else
	{
		*tmpReadLen = 0;

		fclose(fpFasta);
		fpFasta = NULL;

		printf("Reads in data sets should be in equal size.\n");
		return FAILED;
	}

	fclose(fpFasta);
	fpFasta = NULL;

	return SUCCESSFUL;
}

/**
 * Get the read length from fastq file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadLenFromFastqInScaf(int *tmpReadLen, const char *fastqFile)
{
	FILE *fpFastq;
	char ch, seq_data[5000];
	int64_t i, line_index, tmpReadsNum, tmpLen, sumLen, validFlag;

	fpFastq = fopen(fastqFile, "r");
	if(fpFastq==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fastqFile);
		return FAILED;
	}

	sumLen = 0;
	tmpReadsNum = 0;
	line_index = 0;
	while(!feof(fpFastq))
	{
		if(line_index==0)  //the sequence name line
		{
			ch = fgetc(fpFastq);
			while(ch!='\n' && ch!=-1)
			{
				ch = fgetc(fpFastq);
			}
		}else if(line_index==1)  //the sequence line
		{
			tmpLen = 0;
			ch = fgetc(fpFastq);
			while(ch!='\n')
			{
				seq_data[tmpLen++] = ch;
				ch = fgetc(fpFastq);
			}
			seq_data[tmpLen] = '\0';
		}else if(line_index==2)  //the sequence name line
		{
			ch = fgetc(fpFastq);
			while(ch!='\n')
			{
				ch = fgetc(fpFastq);
			}
		}else
		{
			ch = fgetc(fpFastq);
			while(ch!='\n'  && ch!=-1) ch = fgetc(fpFastq);
		}
		line_index++;

		if(line_index==4)  //a read is read finished
		{
			sumLen += tmpLen;

			tmpReadsNum ++;
			if(tmpReadsNum>=MAX_READ_NUM_READ_LEN_SAMPLE)
				break;

			line_index = 0;
		}
	}

	if(tmpReadsNum>0)
	{
		*tmpReadLen = sumLen / tmpReadsNum;
	}else
	{
		*tmpReadLen = 0;

		fclose(fpFastq);
		fpFastq = NULL;

		printf("Reads in data sets should be in equal size.\n");
		return FAILED;
	}

	fclose(fpFastq);
	fpFastq = NULL;

	return SUCCESSFUL;
}

/**
 * Start scaffolding.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short startScaffolding()
{
	//To do ...
/*
	// Warning: this function has been depreciated
	if(generateContigsFasta("contigs.fa", "contig.txt")==FAILED)
	{
		printf("line=%d, In %s(), cannot get contigs in fasta, error!\n", __LINE__, __func__);
		return FAILED;
	}
*/


	if(filterShortContigs(contigsFileFiltered, contigsFileUltraShort, contigsFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot filter contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}




	if(buildContigIndex(contigsFileFiltered, contigIndexFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot build contig index (CI), error!\n", __LINE__, __func__);
		return FAILED;
	}



	if(mapPEs(readSeqFile, readMatchFiles[0], readMatchFiles[1], readFilesInput, readFileNum, readsFileFormatType, pairedMode, contigIndexFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot map PEs, error!\n", __LINE__, __func__);
		return FAILED;
	}



	if(buildReadLists(readListFiles[0], readListFiles[1], readListFiles[2], readMatchFiles[0], readMatchFiles[1])==FAILED)
	{
		printf("line=%d, In %s(), cannot build read lists, error!\n", __LINE__, __func__);
		return FAILED;
	}



	if(buildContigList(contigListFile, readListFiles[2], contigsFileFiltered)==FAILED)
	{
		printf("line=%d, In %s(), cannot build Contig List (CL), error!\n", __LINE__, __func__);
		return FAILED;
	}




	if(contigsLinking(linkInfoFile, averLinkNumFile, contigsFileFiltered, readListFiles[2], contigListFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot link contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}



	if(overlapContigsInScaffolds(contigOverlapFile, meanSdevFile, contigsFileAfterOverlap, readListFilesAfterOverlap[2], readListFilesAfterOverlap[0], readListFilesAfterOverlap[1], contigListFileAfterOverlap, linkInfoFile, contigsFileFiltered, readListFiles[2], readListFiles[0], readListFiles[1], contigListFile, averLinkNumFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot overlap contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}



	if(gapFillFlag==YES)
	{
		if(gapFilling(scafGraphFile, monoReadSeqFile, readListFiles[3], contigOverlapFileAfterFapFilling, contigsFileAfterGapFilling, readListFilesAfterOverlap[0], readListFilesAfterOverlap[1], readSeqFile, meanSdevFile, contigOverlapFile, contigsFileAfterOverlap)==FAILED)
		{
			printf("line=%d, In %s(), cannot fill gaps, error!\n", __LINE__, __func__);
			return FAILED;
		}



		if(generateScaffoldSequence(scafSeqFile, contigOverlapFileAfterFapFilling, contigsFileAfterGapFilling)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate scaffold sequences, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{

		if(generateScaffoldSequence(scafSeqFile, contigOverlapFile, contigsFileAfterOverlap)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate scaffold sequences, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}



	return SUCCESSFUL;
}
