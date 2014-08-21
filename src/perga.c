
#include "inc/stdinc.h"
#include "inc/extvab.h"

/**
 *  @return:
 *  	If succeeds, return SUCCESSFUL;
 *  	if the job is executed and failed, return ERROR;
 *  	otherwise, return FAILED.
 **/
short startPERGA(int operationModePara, int kmerSizePara, int readLenCutOffPara, int pairedModePara, char **readFilesPara, int readFileNumPara, int singleBaseQualThresPara, char *graphFilePara, char *contigsFilePara, char *readMatchInfoFilePara, double meanSizeInsertPara, double standardDevPara, char *outputPathPara, char *outputPrefixPara, int minContigLenPara, int contigAlignRegLenPara, char *gapFillFlagPara)
{
	struct timeval tp_start,tp_end;
	double time_used;
	gettimeofday(&tp_start,NULL);

	// initialize the global parameters
	if(initGlobalParas(operationModePara, outputPathPara, outputPrefixPara, readFilesPara, readFileNumPara, singleBaseQualThresPara, pairedModePara, kmerSizePara, readLenCutOffPara, graphFilePara, contigsFilePara, readMatchInfoFilePara, meanSizeInsertPara, standardDevPara, minContigLenPara, contigAlignRegLenPara, gapFillFlagPara)==FAILED)
	{
		//printf("line=%d, In %s(), cannot initialize the global parameters, error!\n", __LINE__, __func__);
		return FAILED;
	}


	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_HASHTABLE)
	{
		if(constructGraph(graphFile, readFilesInput, readFileNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot construct the graph, error!\n", __LINE__, __func__);
			return ERROR;
		}

		if(operationMode==OPERATION_MODE_HASHTABLE)
			releaseGraph(&deBruijnGraph);  //free k-mer hash table
	}


	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_ASSEMBLE || operationMode==OPERATION_MODE_ASSEM_SCAF)
	{
		if(operationMode==OPERATION_MODE_ASSEMBLE || operationMode==OPERATION_MODE_ASSEM_SCAF)
		{
			// load the graph to memory
			if(loadGraph(&deBruijnGraph, graphFile)==FAILED)
			{
				printf("line=%d, In %s(), cannot load graph to memory, error!\n", __LINE__, __func__);
				return ERROR;
			}

			// ############################ Debug information ##############################
			//if(checkGraph(deBruijnGraph)==FAILED)
			//{
			//	printf("line=%d, In %s(), checking graph error!\n", __LINE__, __func__);
			//	return FAILED;
			//}
			// ############################ Debug information ##############################
		}


		if(initContigGraph(&contigGraph)==FAILED)
		{
			printf("line=%d, In %s(), cannot initialize the contig graph, error!\n", __LINE__, __func__);
			return ERROR;
		}

		// build contigs
		if(buildContigs(contigsFileFasta, graphFile, readMatchInfoFile)==FAILED)
		{
			printf("line=%d, In %s(), cannot build contigs, error!\n", __LINE__, __func__);
			return ERROR;
		}

		if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_ASSEM_SCAF)
		{ // 'all' or 'assem_sf'

			readSet = deBruijnGraph->readSet;

			//clean k-mer hash table
			if(cleanKmerInfoInGraph(&deBruijnGraph)==FAILED)
			{
				printf("line=%d, In %s(), cannot clean k-mer hash table, error!\n", __LINE__, __func__);
				return ERROR;
			}
		}else if(operationMode==OPERATION_MODE_ASSEMBLE)
		{
			releaseGraph(&deBruijnGraph);  //free k-mer hash table
			releaseContigGraph(&contigGraph);
		}
	}

	// do the scaffolding
	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_SCAFFOLDING || operationMode==OPERATION_MODE_ASSEM_SCAF)
	{
		if(operationMode==OPERATION_MODE_SCAFFOLDING)
		{ // 'sf'
			// load the readSet from k-mer hash table
			if(loadReadSetFromGraphFile(&readSet, graphFile)==FAILED)
			{
				printf("line=%d, In %s(), cannot get the read set from k-mer hash table, error!\n", __LINE__, __func__);
				return ERROR;
			}

			// load the contig graph
			if(loadContigGraph(&contigGraph, contigsFileFasta)==FAILED)
			{
				printf("line=%d, In %s(), cannot the contig graph, error!\n", __LINE__, __func__);
				return ERROR;
			}
		}

		// start scaffolding
		if(startScaffolding(scafSeqFile, readMatchInfoFile, contigGraph, readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the scaffolds, error!\n", __LINE__, __func__);
			return ERROR;
		}

		releaseReadset(&readSet);
		releaseContigGraph(&contigGraph);

	}

	freeGlobalParas();

	gettimeofday(&tp_end,NULL);
	time_used = tp_end.tv_sec-tp_start.tv_sec + (double)(tp_end.tv_usec-tp_start.tv_usec)/1000000;

	printf("\nTotal Used Time: %f Seconds.\n", time_used);

    return SUCCESSFUL;
}


/**
 * Initialize the global parameters.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initGlobalParas(int operationModePara, char *outputPathName, char *prefix, char **readFilesPara, int readFileNumPara, int singleBaseQualThresPara, int pairedModePara, int kmerLen, int readLenCut, char *graphFilePara, char *contigsFilePara, char *readMatchInfoFilePara, double meanSizeInsertPara, double standardDevPara, int minContigLenPara, int contigAlignRegLenPara, char *gapFillFlagPara)
{
	int i, prefixLen, memInGB;
	char kmerSizeStr[20], readLenStr[20];
	struct stat st;

	//printf("\n============= Begin setting global parameters, please wait ... =============\n");

	prefixLen = strlen(prefix);


	operationMode = operationModePara;

	if(setGlobalPath(outputPathName)==FAILED)
	{
		printf("line=%d, In %s(), cannot set global paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_HASHTABLE)
	{ // 'all' or 'ht'
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

		// get reads file format type
		if(getReadsFileFormat(&readsFileFormatType, readFilesInput, readFileNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot get reads file formats, error!\n", __LINE__, __func__);
			return FAILED;
		}


		// get the read length from fastq file
		if(readsFileFormatType==FILE_FORMAT_FASTA)
		{
			if(getReadLenAndCountFromFilesFasta(&totalReadNumSample, &averReadLenInFileSample, &maxReadLenInFileSample, &minReadLenInFileSample, readFilesInput, readFileNum)==FAILED)
			{
				//printf("line=%d, In %s(), cannot get read length and count, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else if(readsFileFormatType==FILE_FORMAT_FASTQ)
		{
			if(getReadLenAndCountFromFilesFastq(&totalReadNumSample, &averReadLenInFileSample, &maxReadLenInFileSample, &minReadLenInFileSample, readFilesInput, readFileNum)==FAILED)
			{
				//printf("line=%d, In %s(), cannot get read length and count, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		if(readLenCut==0)
		{ // no read length cut
			readLenCutOff = MAX_READ_LEN_IN_BUF;
			readLen = averReadLenInFileSample;
		}else
		{
			if(readLenCut>MAX_READ_LEN_IN_BUF)
			{
				readLenCut = MAX_READ_LEN_IN_BUF;
				printf("The read length cut off threshold is set to be %d if you do not mind.\n", readLenCut);
			}

			readLenCutOff = readLenCut;
			if(readLenCutOff>averReadLenInFileSample)
				readLen = averReadLenInFileSample;
			else
				readLen = readLenCutOff;
		}


		if(kmerLen==0)
		{
			kmerSize = DEFAULT_KMER_SIZE;
		}else
		{
			kmerSize = kmerLen;
		}

		if(kmerSize%2==0)
			kmerSize --;

		if(kmerSize>readLen)
		{
			printf("Exception: invalid k-mer size, it is larger than the read length.\n");
			return FAILED;
		}


		hashTableSizeReadseq = HASH_TABLE_SIZE_SMALL;

		memInGB = getSysMemorySize();
		if(memInGB>100)
			hashTableSize = HASH_TABLE_SIZE_LARGE;
		else if(memInGB>10)
			hashTableSize = HASH_TABLE_SIZE_MEDIUM;
		else
			hashTableSize = HASH_TABLE_SIZE_SMALL;


	}else if(operationMode==OPERATION_MODE_ASSEMBLE || operationMode==OPERATION_MODE_ASSEM_SCAF || operationMode==OPERATION_MODE_SCAFFOLDING)
	{ // 'assem'
		strcpy(graphFile, graphFilePara);
		if(stat(graphFile, &st)==-1)
		{
			printf("Exception: please specify correct graph file.\n");
			return FAILED;
		}

		// load the graph to memory
		if(GlobalParasFromGraph(&readLen, &kmerSize, &hashTableSize, &hashTableSizeReadseq, &pairedMode, &PEGivenType, &meanSizeInsert, &standardDev, graphFile)==FAILED)
		{
			//printf("line=%d, In %s(), cannot load graph to memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(operationMode==OPERATION_MODE_SCAFFOLDING)
		{
			strcpy(contigsFileFasta, contigsFilePara);
			if(stat(contigsFileFasta, &st)==-1)
			{
				printf("Exception: please specify correct contigs file.\n");
				return FAILED;
			}

			// the reads match information file
			strcpy(readMatchInfoFile, readMatchInfoFilePara);
			//if(stat(readMatchInfoFile, &st)==-1)
			//{
			//	printf("Exception: please specify correct reads match information file.\n");
			//	return FAILED;
			//}

			// scafSeqFile
			strcpy(scafSeqFile, outputPathStr);
			if(prefixLen>0)
			{
				strcat(scafSeqFile, prefix);
				strcat(scafSeqFile, "_");
			}
			strcat(scafSeqFile, "scaffolds.fa");
		}

	}else
	{
		printf("Exception: please specify correct command, error!\n");
		return FAILED;
	}

	kmerRegLenRatioEnd5 = KMER_REG_LEN_RATIO_END5;
	kmerRegLenRatioEnd3 = KMER_REG_LEN_RATIO_END3;
	maxUnknownBaseNumPerRead = MAX_UNKNOWN_BASE_NUM;
	reserveHashItemBlocksFlag = RESERVE_HASH_ITEM_READ_SET;

	if(singleBaseQualThresPara>0)
		singleBaseQualThres = singleBaseQualThresPara;
	else
		singleBaseQualThres = SINGLE_QUAL_THRESHOLD;

	readSimilarityThres = READ_SIMILARITY_THRES;

	//sprintf(kmerSizeStr, "%d", kmerSize);
	//sprintf(readLenStr, "%d", readLen);

	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_ASSEMBLE || operationMode==OPERATION_MODE_ASSEM_SCAF)
	{
		strcpy(contigsFileFasta, outputPathStr);
		if(prefixLen>0)
		{
			strcat(contigsFileFasta, prefix);
			strcat(contigsFileFasta, "_");
		}
		strcat(contigsFileFasta, "contigs.fa");
		strcpy(contigsFileHanging, outputPathStr);
		if(prefixLen>0)
		{
			strcat(contigsFileHanging, prefix);
			strcat(contigsFileHanging, "_");
		}
		strcat(contigsFileHanging, "contigs_hanging.txt");

		strcpy(fragmentSizeFile, outputPathStr);
		if(prefixLen>0)
		{
			strcat(fragmentSizeFile, prefix);
			strcat(fragmentSizeFile, "_");
		}
		strcat(fragmentSizeFile, "fragmentSize.bin");

		strcpy(sampleContigsFile, outputPathStr);
		if(prefixLen>0)
		{
			strcat(sampleContigsFile, prefix);
			strcat(sampleContigsFile, "_");
		}
		strcat(sampleContigsFile, "sampleContigs.fa");

		strcpy(readMatchInfoFile, outputPathStr);
		if(prefixLen>0)
		{
			strcat(readMatchInfoFile, prefix);
			strcat(readMatchInfoFile, "_");
		}
		strcat(readMatchInfoFile, "readsMatchInfo.bin");

		// scafSeqFile
		strcpy(scafSeqFile, outputPathStr);
		if(prefixLen>0)
		{
			strcat(scafSeqFile, prefix);
			strcat(scafSeqFile, "_");
		}
		strcat(scafSeqFile, "scaffolds.fa");
	}

	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_HASHTABLE)
	{
		strcpy(graphFile, outputPathStr);
		if(prefixLen>0)
		{
			strcat(graphFile, prefix);
			strcat(graphFile, "_");
		}
		strcat(graphFile, "hashtable.bin");
	}

	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_SCAFFOLDING || operationMode==OPERATION_MODE_ASSEM_SCAF)
	{
		// the contig align region size
		if(contigAlignRegLenPara>0)
		{
			contigAlignRegSize = contigAlignRegLenPara;
		}else
		{
			contigAlignRegSize = CONTIG_ALIGN_REG_SIZE;
			if(contigAlignRegSize<readLen*CONTIG_ALIGN_REG_SIZE_FACTOR)
			{
				contigAlignRegSize = readLen*CONTIG_ALIGN_REG_SIZE_FACTOR;
				printf("The minContigLenThres will be set to %d instead.\n", contigAlignRegSize);
			}
		}

		// the gap filling flag
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
	}


#if (DEBUG_PARA_PRINT==YES)
	printf("\nhashTableSize=%lu\n", hashTableSize);
	printf("kmerRegLenRatioEnd5=%.2f\n", kmerRegLenRatioEnd5);
	printf("kmerRegLenRatioEnd3=%.2f\n", kmerRegLenRatioEnd3);
	printf("kmerSampleInterval=%d\n", KMER_SAMPLE_INTERVAL);
	printf("maxUnknownBaseNumPerRead=%d\n", maxUnknownBaseNumPerRead);
	printf("reserveHashItemBlocksFlag=%d\n\n", reserveHashItemBlocksFlag);
#endif


	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_HASHTABLE)
		pairedMode = pairedModePara;

	//if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_ASSEMBLE || operationMode==OPERATION_MODE_SCAFFOLDING || operationMode==OPERATION_MODE_ASSEM_SCAF)
	{
		if(pairedMode==0)
		{
			PEGivenType = SE_GIVEN_TYPE;
			meanSizeInsert = 0;
			standardDev = 0;
		}

		if(meanSizeInsertPara==0 && meanSizeInsert==0)
		{
			PEGivenType = NONE_PE_GIVEN_TYPE;
			meanSizeInsert = 0;
			standardDev = 0;
		}else if(meanSizeInsertPara!=0 && standardDevPara==0)
		{
			PEGivenType = INSERT_PE_GIVEN_TYPE;
			meanSizeInsert = meanSizeInsertPara;
			standardDev = meanSizeInsert * DRAFT_SDEV_FACTOR;
		}else if(meanSizeInsertPara!=0 && standardDevPara!=0)
		{
			PEGivenType = BOTH_PE_GIVEN_TYPE;
			meanSizeInsert = meanSizeInsertPara;
			standardDev = standardDevPara;
		}

		minContigLen = minContigLenPara;
		if(minContigLen==0)
		{
			minContigLen = CONTIG_LEN_THRESHOLD;
			if(minContigLen<CONTIG_LEN_FACTOR*readLen)
				minContigLen = CONTIG_LEN_FACTOR * readLen;
		}
	}


	entriesPerKmer = ((kmerSize-1) / 32) + 1;
	if(kmerSize%32==0)
	{
		lastEntryMask = (uint64_t) -1;
		lastEntryBaseNum = 32;
	}else
	{
		lastEntryMask = (1LLU << ((kmerSize%32)<<1)) - 1;
		lastEntryBaseNum = kmerSize % 32;
	}

	// allocate memory
	kmerSeqInt = (uint64_t*) malloc(entriesPerKmer*sizeof(uint64_t));
	if(kmerSeqInt==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	kmerSeqIntRev = (uint64_t*) malloc(entriesPerKmer*sizeof(uint64_t));
	if(kmerSeqIntRev==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//===================== SVM model begin ======================
#if (SVM_NAVI==YES)
	if(loadSvmModel(&svmModelPE, svmKF_PE, svmSV_PE, svmAlpha_PE, svmBias_PE, svmScaleData_PE, rowsNumSupportVector_PE, colsNumSupportVector_PE, svmUseScaleData_PE)==FAILED)
	{
		printf("Load SVM model error!\n");
		return ERROR;
	}
	if(loadSvmModel(&svmModelSE, svmKF_SE, svmSV_SE, svmAlpha_SE, svmBias_SE, svmScaleData_SE, rowsNumSupportVector_SE, colsNumSupportVector_SE, svmUseScaleData_SE)==FAILED)
	{
		printf("Load SVM model error!\n");
		return ERROR;
	}

	if(initSampleMemSvm(svmModelPE->colsNumSupportVector)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory for sample data when using SVM model, error!\n", __LINE__, __func__);
		return ERROR;
	}
#endif
	//===================== SVM model end ======================

	// output the variables
	if(operationMode==OPERATION_MODE_ALL)
	{
		printf("\nkmer size          : %d\n", kmerSize);
		//printf("read length cutoff : %d\n", readLenCutOff);
		printf("read length        : %d\n", readLen);
		printf("paired mode        : %d\n", pairedMode);
		for(i=0; i<readFileNumPara; i++)
		printf("read files[%d]      : %s\n", i, readFilesInput[i]);
		//printf("single qual thres  : %d\n", singleBaseQualThres);
		if(meanSizeInsert>0)
		{
			printf("insert size        : %.2f\n", meanSizeInsert);
			printf("insert size sdev.  : %.2f\n", standardDev);
		}
		printf("output directory   : %s\n", outputPathStr);
		printf("hash table file    : %s\n", graphFile);
		//printf("reads match file   : %s\n", readMatchInfoFile);
		printf("contig file        : %s\n", contigsFileFasta);
		printf("scaffolds file     : %s\n", scafSeqFile);
		printf("minimal contig size: %d\n", minContigLen);
		printf("align region size  : %d\n", contigAlignRegSize);
		if(gapFillFlag==YES)
			printf("gap filling flag   : yes\n\n");
		else
			printf("gap filling flag   : no\n\n");
	}else if(operationMode==OPERATION_MODE_HASHTABLE)
	{
		printf("\nkmer size          : %d\n", kmerSize);
		//printf("read length cutoff : %d\n", readLenCutOff);
		printf("read length        : %d\n", readLen);
		printf("paired mode        : %d\n", pairedMode);
		for(i=0; i<readFileNumPara; i++)
		printf("read files[%d]      : %s\n", i, readFilesInput[i]);
		//printf("single qual thres  : %d\n", singleBaseQualThres);
		printf("output directory   : %s\n", outputPathStr);
		printf("hash table file    : %s\n", graphFile);
	}else if(operationMode==OPERATION_MODE_ASSEMBLE)
	{
		printf("\nhash table file    : %s\n", graphFile);
		printf("kmer size          : %d\n", kmerSize);
		printf("read length        : %d\n", readLen);
		printf("paired mode        : %d\n", pairedMode);
		if(meanSizeInsert>0)
		{
			printf("insert size        : %.2f\n", meanSizeInsert);
			printf("insert size sdev.  : %.2f\n", standardDev);
		}
		printf("output directory   : %s\n", outputPathStr);
		//printf("reads match file   : %s\n", readMatchInfoFile);
		printf("contig file        : %s\n", contigsFileFasta);
		printf("minimal contig size: %d\n\n", minContigLen);
	}else if(operationMode==OPERATION_MODE_SCAFFOLDING || operationMode==OPERATION_MODE_ASSEM_SCAF)
	{
		printf("\noutput directory   : %s\n", outputPathStr);
		//printf("reads match file   : %s\n", readMatchInfoFile);
		printf("contig file        : %s\n", contigsFileFasta);
		printf("scaffolds file     : %s\n", scafSeqFile);
		printf("minimal contig size: %d\n", minContigLen);
		if(meanSizeInsert>0)
		{
			printf("insert size        : %.2f\n", meanSizeInsert);
			printf("insert size sdev.  : %.2f\n", standardDev);
		}
		printf("align region size  : %d\n", contigAlignRegSize);
		if(gapFillFlag==YES)
			printf("gap filling flag   : yes\n\n");
		else
			printf("gap filling flag   : no\n\n");
	}else
	{
		printf("Exception: invalid command %d\n", operationMode);
		return FAILED;
	}


	return SUCCESSFUL;
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
 * Free the global parameters.
 */
void freeGlobalParas()
{
	int i;
	if(operationMode==OPERATION_MODE_ALL || operationMode==OPERATION_MODE_HASHTABLE)
		for(i=0; i<readFileNum; i++) { free(readFilesInput[i]); readFilesInput[i] = NULL; }
	free(kmerSeqInt);
	kmerSeqInt = NULL;
	free(kmerSeqIntRev);
	kmerSeqIntRev = NULL;

	//===================== SVM model begin ======================
#if (SVM_NAVI==YES)
	freeSvmModel(&svmModelPE);
	freeSvmModel(&svmModelSE);
	freeSampleMemSvm();
#endif
	//===================== SVM model end ======================

}

/**
 *
 */
short getReadsFileFormat(int *readsFileFormatType, char **readFilesInput, int readFileNum)
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
 * Get the read length and total count from fastq files.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadLenAndCountFromFilesFasta(int64_t *tmpTotalReadNum, int *tmpAverReadLenInFile, int *tmpMaxReadLenInFile, int *tmpMinReadLenInFile, char **readsFileNames, int readsFileNum)
{
	int64_t i, tmpReadCount, sumReadLen, tmpSumReadLen, tmpMaxReadLen, tmpMinReadLen;

	(*tmpTotalReadNum) = 0;
	sumReadLen = 0;
	(*tmpMaxReadLenInFile) = 0;
	(*tmpMinReadLenInFile) = INT_MAX;

//	tmpReadCount = 0;
//	tmpSumReadLen = 0;
//	tmpMaxReadLen = 0;
//	tmpMinReadLen = INT_MAX;

	for(i=0; i<readFileNum; i++)
	{
		if(getReadLenAndCountFromSingleFasta(&tmpReadCount, &tmpSumReadLen, &tmpMaxReadLen, &tmpMinReadLen, readsFileNames[i])==FAILED)
		{
			//printf("line=%d, In %s(), cannot get readlenInFile, error!\n", __LINE__, __func__);
			return FAILED;
		}
		(*tmpTotalReadNum) += tmpReadCount;
		sumReadLen += tmpSumReadLen;

		if(tmpMaxReadLen>(*tmpMaxReadLenInFile)) (*tmpMaxReadLenInFile) = tmpMaxReadLen;
		if(tmpMinReadLen<(*tmpMinReadLenInFile)) (*tmpMinReadLenInFile) = tmpMinReadLen;
	}

	if((*tmpTotalReadNum)>0)
		(*tmpAverReadLenInFile) = round((double)sumReadLen / (*tmpTotalReadNum));
	else
	{
		*tmpAverReadLenInFile = 0;
		return FAILED;
	}

#if (DEBUG_PARA_PRINT==YES)
	printf("Total reads: %lu\n", (*tmpTotalReadNum));
	printf("averReadLen: %d\n", (*tmpAverReadLenInFile));
	printf("maxReadLen: %d\n", (*tmpMaxReadLenInFile));
	printf("minReadLen: %d\n", (*tmpMinReadLenInFile));
#endif

	return SUCCESSFUL;
}

/**
 * Get the read length and total count from fastq files.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadLenAndCountFromFilesFastq(int64_t *tmpTotalReadNum, int *tmpAverReadLenInFile, int *tmpMaxReadLenInFile, int *tmpMinReadLenInFile, char **readsFileNames, int readsFileNum)
{
	int64_t i, tmpReadCount, sumReadLen, tmpSumReadLen, tmpMaxReadLen, tmpMinReadLen;

	(*tmpTotalReadNum) = 0;
	sumReadLen = 0;
	(*tmpMaxReadLenInFile) = 0;
	(*tmpMinReadLenInFile) = INT_MAX;

	for(i=0; i<readFileNum; i++)
	{
		if(getReadLenAndCountFromSingleFastq(&tmpReadCount, &tmpSumReadLen, &tmpMaxReadLen, &tmpMinReadLen, readsFileNames[i])==FAILED)
		{
			//printf("line=%d, In %s(), cannot get readlenInFile, error!\n", __LINE__, __func__);
			return FAILED;
		}
		(*tmpTotalReadNum) += tmpReadCount;
		sumReadLen += tmpSumReadLen;

		if(tmpMaxReadLen>(*tmpMaxReadLenInFile)) (*tmpMaxReadLenInFile) = tmpMaxReadLen;
		if(tmpMinReadLen<(*tmpMinReadLenInFile)) (*tmpMinReadLenInFile) = tmpMinReadLen;
	}

	if((*tmpTotalReadNum)>0)
		(*tmpAverReadLenInFile) = ceil((double)sumReadLen / (*tmpTotalReadNum));
	else
	{
		*tmpAverReadLenInFile = 0;
		return FAILED;
	}

#if (DEBUG_PARA_PRINT==YES)
	printf("Total reads: %lu\n", (*tmpTotalReadNum));
	printf("averReadLen: %d\n", (*tmpAverReadLenInFile));
	printf("maxReadLen: %d\n", (*tmpMaxReadLenInFile));
	printf("minReadLen: %d\n", (*tmpMinReadLenInFile));
#endif

	return SUCCESSFUL;
}

/**
 * Get the read length and total count from single fasta file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadLenAndCountFromSingleFasta(int64_t *tmpReadCount, int64_t *tmpSumReadLen, int64_t *tmpMaxReadLen, int64_t *tmpMinReadLen, char *fastaFile)
{
	FILE *fpFasta;
	char ch;
	int tmpLen;

	fpFasta = fopen(fastaFile, "r");
	if(fpFasta==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fastaFile);
		return FAILED;
	}

	*tmpReadCount = 0;
	*tmpSumReadLen = 0;
	*tmpMaxReadLen = 0;
	*tmpMinReadLen = INT_MAX;

	ch = fgetc(fpFasta);
	while(!feof(fpFasta))
	{
		while(ch!='\n') ch = fgetc(fpFasta);

		tmpLen = 0;
		ch = fgetc(fpFasta);
		while(ch!='>' && ch!=-1)
		{
			if(ch!='\n')
				tmpLen ++;
			ch = fgetc(fpFasta);
		}

		(*tmpReadCount) ++;
		(*tmpSumReadLen) += tmpLen;

		if(tmpLen>(*tmpMaxReadLen)) *tmpMaxReadLen = tmpLen;
		if(tmpLen<(*tmpMinReadLen)) *tmpMinReadLen = tmpLen;

		if((*tmpReadCount)>=READS_NUM_PER_FILE_SAMPLE)
			break;
	}

	fclose(fpFasta);
	fpFasta = NULL;

	return SUCCESSFUL;
}

/**
 * Get the read length and total count from single fastq file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReadLenAndCountFromSingleFastq(int64_t *tmpReadCount, int64_t *tmpSumReadLen, int64_t *tmpMaxReadLen, int64_t *tmpMinReadLen, char *fastqFile)
{
	FILE *fpFastq;
	int line_index, tmpLen;
	char ch;

	fpFastq = fopen(fastqFile, "r");
	if(fpFastq==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fastqFile);
		return FAILED;
	}

	*tmpReadCount = 0;
	*tmpSumReadLen = 0;
	*tmpMaxReadLen = 0;
	*tmpMinReadLen = INT_MAX;

	line_index = 0;
	ch = fgetc(fpFastq);
	while(!feof(fpFastq))
	{
		if(line_index==0)  //the sequence name line
		{
			ch = fgetc(fpFastq);
			while(ch!='\n' && ch!=-1)
				ch = fgetc(fpFastq);
		}else if(line_index==1)  //the sequence line
		{
			tmpLen = 0;
			ch = fgetc(fpFastq);
			while(ch!='\n' && ch!=-1)
			{
				tmpLen ++;
				ch = fgetc(fpFastq);
			}
		}else if(line_index==2)  //the sequence name line
		{
			ch = fgetc(fpFastq);
			while(ch!='\n' && ch!=-1)
				ch = fgetc(fpFastq);
		}else
		{
			ch = fgetc(fpFastq);
			while(ch!='\n' && ch!=-1)
				ch = fgetc(fpFastq);
		}
		line_index++;

		if(line_index==4)  //the sequence is read finished, construct the read
		{
			(*tmpReadCount) ++;
			(*tmpSumReadLen) += tmpLen;

			if(tmpLen>(*tmpMaxReadLen)) *tmpMaxReadLen = tmpLen;
			if(tmpLen<(*tmpMinReadLen)) *tmpMinReadLen = tmpLen;

			line_index = 0;

			if((*tmpReadCount)>=READS_NUM_PER_FILE_SAMPLE)
				break;
		}
	}

	fclose(fpFastq);
	fpFastq = NULL;

	return SUCCESSFUL;
}

/**
 * Get system memory in GB.
 *  @return:
 *  	Return the memory in GB.
 */
int32_t getSysMemorySize()
{

	int64_t page_size;
    int64_t num_pages;
    int64_t free_pages;
    int64_t  mem_size;

    page_size = sysconf (_SC_PAGESIZE);
    num_pages = sysconf (_SC_PHYS_PAGES);
    free_pages = sysconf (_SC_AVPHYS_PAGES);

    mem_size = (int64_t)num_pages * (int64_t)page_size;
    mem_size >>= 30;

    return  mem_size;
}
