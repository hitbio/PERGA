
#include "inc/stdinc.h"
#include "inc/extvab.h"

/**
 *  @return:
 *  	If succeeds, return SUCCESSFUL;
 *  	if the job is executed and failed, return ERROR;
 *  	otherwise, return FAILED.
 **/
short startSRGA(int operationModePara, int kmerSizePara, int readLenCutOffPara, int pairedModePara, char **readFilesPara, int readFileNumPara, int singleBaseQualThresPara, int errorCorrectionPara, char *graphFilePara, double meanSizeInsertPara, double standardDevPara, char *outputPathPara, char *outputPrefixPara, int minContigLenPara)
{
	struct timeval tp_start,tp_end;
	double time_used;
	gettimeofday(&tp_start,NULL);


//	if(testSVM()==FAILED)
//	{
//		printf("test SVM error!\n");
//		return ERROR;
//	}

//	if(outputClassificationResult()==FAILED)
//	{
//		printf("output classification result error!\n");
//		return ERROR;
//	}

//	if(resultStatistics()==FAILED)
//	{
//		printf("result statistics error!\n");
//		return ERROR;
//	}

//	if(outputErrPoints()==FAILED)
//	{
//		printf("output error points error!\n");
//		return ERROR;
//	}



	// initialize the global parameters
	if(initGlobalParas(operationModePara, outputPathPara, outputPrefixPara, readFilesPara, readFileNumPara, singleBaseQualThresPara, errorCorrectionPara, pairedModePara, kmerSizePara, readLenCutOffPara, graphFilePara, meanSizeInsertPara, standardDevPara, minContigLenPara)==FAILED)
	{
		//printf("line=%d, In %s(), cannot initialize the global parameters, error!\n", __LINE__, __func__);
		return FAILED;
	}


	if(operationMode!=2)
	{
		if(constructGraph(graphFile, readFilesInput, readFileNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot construct the graph, error!\n", __LINE__, __func__);
			return ERROR;
		}
	}


	if(operationMode!=1 || errorCorrectionFlag==YES)
	{
		if(operationMode==2)
		{
			// load the graph to memory
			if(loadGraph(&deBruijnGraph, graphFile)==FAILED)
			{
				printf("line=%d, In %s(), cannot load graph to memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// ############################ Debug information ##############################
			//if(checkGraph(deBruijnGraph)==FAILED)
			//{
			//	printf("line=%d, In %s(), checking graph error!\n", __LINE__, __func__);
			//	return FAILED;
			//}
			// ############################ Debug information ##############################
		}


		// build contigs
		if(buildContigs(contigsFileFasta, graphFile, occPointFileArr1)==FAILED)
		{
			printf("line=%d, In %s(), cannot build contigs, error!\n", __LINE__, __func__);
			return ERROR;
		}
	}


	if(errorCorrectionFlag==YES)
	{
		cleanKmerInfoInGraph(deBruijnGraph);  //clean k-mer hash table

		errorCorrectionFlag = NO;

		// construct k-mer hash table by corrected reads
		if(constructGraphAfterCorrection(deBruijnGraph, graphFileCorrected, readCorrectedFile)==FAILED)
		{
			printf("line=%d, In %s(), cannot construct graph by corrected reads, error!\n", __LINE__, __func__);
			return ERROR;
		}

		if(operationMode!=1)
		{
			// build contigs
			if(buildContigs(contigsFileFastaCorrected, graphFileCorrected, occPointFileArr2)==FAILED)
			{
				printf("line=%d, In %s(), cannot build contigs by corrected reads, error!\n", __LINE__, __func__);
				return ERROR;
			}
		}
	}

	releaseGraph(deBruijnGraph);  //free k-mer hash table

	freeGlobalParas();

	if(operationMode==0 || operationMode==2)
	{
//		printf("totalReadNum=%ld, validReadNum=%ld, deledReadNum=%ld, validRatio=%.4f, "
//				"successReadNum=%ld, failedReadNum=%ld, failedRatio=%.4f\n",
//				totalReadNum, validReadNum, totalReadNum-validReadNum, (float)validReadNum/totalReadNum,
//				successReadNum, validReadNum-successReadNum, (float)(validReadNum-successReadNum)/validReadNum);

		//printf("Corrected read number: %d\n", number_of_corrected_reads);
		//printf("Overlaps less than %d: %d\n", MIN_OVERLAP_LEN, number_of_overlap_less_than_threshold);
		//printf("searchNum=%ld, totalSearchLen=%llu, fold=%lld\n", searchNum, totalSearchLen, totalSearchLen/searchNum);
	}


	gettimeofday(&tp_end,NULL);
	time_used = tp_end.tv_sec-tp_start.tv_sec + (double)(tp_end.tv_usec-tp_start.tv_usec)/1000000;

	//printf("Graph Used Mem : %f MB.\n", (double)countMemory/1024/1024);

	printf("\nTotal Used Time: %f Seconds.\n", time_used);


    return SUCCESSFUL;
}


/**
 * Initialize the global parameters.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initGlobalParas(int operationModePara, char *outputPathName, char *prefix, char **readFilesPara, int readFileNumPara, int singleBaseQualThresPara, int errorCorrectionPara, int pairedModePara, int kmerLen, int readLenCut, char *graphFilePara, double meanSizeInsertPara, double standardDevPara, int minContigLenPara)
{
	int i, prefixLen, memInGB;
	char kmerSizeStr[20], readLenStr[20];
	struct stat st;

	//printf("SRGA version : %s\n", VERSION_STR);
	//printf("Release date : %s\n", RELEASE_DATE_STR);

	//printf("\n============= Begin setting global parameters, please wait ... =============\n");

	operationMode = operationModePara;

	if(setGlobalPath(outputPathName)==FAILED)
	{
		printf("line=%d, In %s(), cannot set global paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(operationMode==0 || operationMode==1)
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


//		readLenCutOff = readLenCut;
//		if(readLenCutOff==0)
//			readLen = averReadLenInFileSample;
//		else if(readLenCutOff>averReadLenInFileSample)
//			readLen = averReadLenInFileSample;
//		else
//			readLen = readLenCutOff;

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

		//hashTableSize = HASH_TABLE_SIZE;
		//hashTableSize = 1LLU << (kmerSize << 1);

		memInGB = getSysMemorySize();
		if(memInGB>100)
			hashTableSize = HASH_TABLE_SIZE_LARGE;
		else if(memInGB>10)
			hashTableSize = HASH_TABLE_SIZE_MEDIUM;
		else
			hashTableSize = HASH_TABLE_SIZE_SMALL;


	}else if(operationMode==2)
	{ // 'assemble'
		strcpy(graphFile, graphFilePara);

		if(stat(graphFile, &st)==-1)
		{
			printf("Exception: please specify correct graph file.\n");
			return FAILED;
		}

		// load the graph to memory
		if(GlobalParasFromGraph(&readLen, &averReadLenInFileSample, &kmerSize, &hashTableSize, &hashTableSizeReadseq, &pairedMode, &PEGivenType, &meanSizeInsert, &standardDev, graphFile)==FAILED)
		{
			//printf("line=%d, In %s(), cannot load graph to memory, error!\n", __LINE__, __func__);
			return FAILED;
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

	trimReadLenFlag = TRIM_READ_LEN_FLAG;
	//qualityBaseNumEnd3 = ceil(readLen * QUAL_BASE_NUM_3End_FACTOR);
	//qualityBaseNumEnd5 = readLen - qualityBaseNumEnd3;
	//errorRegLenEnd3 = ceil(readLen * ERROR_REGION_LEN_3End_FACTOR);
	if(singleBaseQualThresPara>0)
		singleBaseQualThres = singleBaseQualThresPara;
	else
		singleBaseQualThres = SINGLE_QUAL_THRESHOLD;

	if(errorCorrectionPara!=0)
		errorCorrectionFlag = YES;
	else
		errorCorrectionFlag = NO;

	readSimilarityThres = READ_SIMILARITY_THRES;

	prefixLen = strlen(prefix);

	sprintf(kmerSizeStr, "%d", kmerSize);
	sprintf(readLenStr, "%d", readLen);

	if(operationMode==0 || operationMode==2)
	{
		strcpy(contigsFileFasta, outputPathStr);
		if(prefixLen>0)
		{
			strcat(contigsFileFasta, prefix);
			strcat(contigsFileFasta, "_");
		}
		strcat(contigsFileFasta, "contigs");
		strcat(contigsFileFasta, "_K");
		strcat(contigsFileFasta, kmerSizeStr);
		strcat(contigsFileFasta, "_R");
		strcat(contigsFileFasta, readLenStr);
		strcat(contigsFileFasta, ".fa");
		strcpy(contigsFileHanging, outputPathStr);
		if(prefixLen>0)
		{
			strcat(contigsFileHanging, prefix);
			strcat(contigsFileHanging, "_");
		}
		strcat(contigsFileHanging, "contigs");
		strcat(contigsFileHanging, "_K");
		strcat(contigsFileHanging, kmerSizeStr);
		strcat(contigsFileHanging, "_R");
		strcat(contigsFileHanging, readLenStr);
		strcat(contigsFileHanging, "_hanging.fa");

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

		strcpy(readCorrectedFile, outputPathStr);
		if(prefixLen>0)
		{
			strcat(readCorrectedFile, prefix);
			strcat(readCorrectedFile, "_");
		}
		strcat(readCorrectedFile, "readsCorrected");
	}

	if(operationMode==0 || operationMode==1)
	{
		strcpy(graphFile, outputPathStr);
		if(prefixLen>0)
		{
			strcat(graphFile, prefix);
			strcat(graphFile, "_");
		}
		strcat(graphFile, "hashtable");
		strcat(graphFile, "_K");
		strcat(graphFile, kmerSizeStr);
		strcat(graphFile, "_R");
		strcat(graphFile, readLenStr);
		strcat(graphFile, ".bin");
	}

	if(errorCorrectionFlag==YES)
	{
		strcpy(graphFileCorrected, graphFile);
		strcat(graphFileCorrected, "_corrected");

		strcpy(contigsFileFastaCorrected, contigsFileFasta);
		strcat(contigsFileFastaCorrected, "_corrected");
	}

#if (DEBUG_PARA_PRINT==YES)
//	printf("outputPathStr=%s\n", outputPathStr);
//	printf("contigsFileFasta=%s\n", contigsFileFasta);
//	printf("contigsFileHanging=%s\n", contigsFileHanging);
//	printf("fragmentSizeFile=%s\n", fragmentSizeFile);
//	printf("graphFile=%s\n", graphFile);
//	printf("sampleContigsFile=%s\n", sampleContigsFile);
//	printf("readCorrectedFile=%s\n", readCorrectedFile);

	//printf("memInGB=%d\n", memInGB);
	printf("hashTableSize=%lu\n", hashTableSize);
	printf("trimReadLenFlag=%d\n", trimReadLenFlag);
	//printf("qualityBaseNumEnd3=%d\n", qualityBaseNumEnd3);
	//printf("qualityBaseNumEnd5=%d\n", qualityBaseNumEnd5);
	//printf("errorRegLenEnd3=%d\n", errorRegLenEnd3);
	printf("kmerRegLenRatioEnd5=%.2f\n", kmerRegLenRatioEnd5);
	printf("kmerRegLenRatioEnd3=%.2f\n", kmerRegLenRatioEnd3);
	printf("maxUnknownBaseNumPerRead=%d\n", maxUnknownBaseNumPerRead);
	printf("reserveHashItemBlocksFlag=%d\n", reserveHashItemBlocksFlag);

//	printf("AVERAGE_QUAL_THRESHOLD_3End=%.4f\n", AVERAGE_QUAL_THRESHOLD_3End);
//	printf("AVERAGE_QUAL_THRESHOLD_5End=%.4f\n", AVERAGE_QUAL_THRESHOLD_5End);
//	printf("SINGLE_QUAL_THRESHOLD=%d\n", SINGLE_QUAL_THRESHOLD);
//	printf("ARTIFACTS_BASE_A_THRESHOLD=%.4f\n", ARTIFACTS_BASE_A_THRESHOLD);
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//PEGivenType = SE_GIVEN_TYPE;			// 0

	//PEGivenType = NONE_PE_GIVEN_TYPE;		// 1

	//PEGivenType = INSERT_PE_GIVEN_TYPE;		// 2
	//meanSizeInsert = 200;

	//PEGivenType = BOTH_PE_GIVEN_TYPE;		// 3
	//meanSizeInsert = 213.59;
	//standardDev = 14.53;

	//kmerSize = 12;
	//kmerSize = 13;
	//kmerSize = 15;
	//kmerSize = 17;
	//kmerSize = 19;
	//kmerSize = 21;
	//kmerSize = 23;
	//kmerSize = 25;
	//kmerSize = 27;
	//kmerSize = 29;
	//kmerSize = 31;

	if(operationMode==0 || operationMode==1)
		pairedMode = pairedModePara;

	if(operationMode==0 || operationMode==2)
	{
		if(pairedMode==0)
		{
			PEGivenType = SE_GIVEN_TYPE;
			meanSizeInsert = 0;
			standardDev = 0;
		}else if(meanSizeInsertPara==0)
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



	//printf("PEGivenType=%d\n", PEGivenType);
	//if(PEGivenType==INSERT_PE_GIVEN_TYPE)
	//{
		//printf("meanSizeInsert=%.2f\n", meanSizeInsert);
	//}else if(PEGivenType==BOTH_PE_GIVEN_TYPE)
	//{
		//printf("meanSizeInsert=%.2f\n", meanSizeInsert);
		//printf("standardDev=%.2f\n", standardDev);
	//}

	//printf("readLen=%d\n", readLen);
	//printf("kmerSize=%d\n", kmerSize);
	//printf("entriesPerKmer=%d\n", entriesPerKmer);
	//printf("lastEntryBaseNum=%d\n", lastEntryBaseNum);
	//printf("lastEntryMask=0x%lX\n", lastEntryMask);
	//printf("hashTableSize=%lu\n", hashTableSize);

	//bytesPerKmerseq = entriesPerKmer * sizeof(uint64_t);
	//itemNumIncrementKmerSeqArr = BLOCK_SIZE_INCREMENT_KMER_SEQ / bytesPerKmerseq;

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

	//===================== draw curve variable begin ======================
#if(DRAW_CURVE_FLAG==YES)
	strcpy(occPointFile, "../output/occPoints.txt");
	strcpy(occExtensionCorrectFile[0], "../output/occExtensionCorrect.txt");
	strcpy(occExtensionCorrectFile[1], "../output/occExtensionCorrect_PE.txt");
	strcpy(occExtensionCorrectFile[2], "../output/occExtensionCorrect_SPE.txt");
	strcpy(occExtensionCorrectFile[3], "../output/occExtensionCorrect_SE.txt");
	strcpy(occExtensionIncorrectFile[0], "../output/occExtensionIncorrect.txt");
	strcpy(occExtensionIncorrectFile[1], "../output/occExtensionIncorrect_PE.txt");
	strcpy(occExtensionIncorrectFile[2], "../output/occExtensionIncorrect_SPE.txt");
	strcpy(occExtensionIncorrectFile[3], "../output/occExtensionIncorrect_SE.txt");
	strcpy(occStopCorrectFile[0], "../output/occStopCorrect.txt");
	strcpy(occStopCorrectFile[1], "../output/occStopCorrect_PE.txt");
	strcpy(occStopCorrectFile[2], "../output/occStopCorrect_SPE.txt");
	strcpy(occStopCorrectFile[3], "../output/occStopCorrect_SE.txt");
	strcpy(occStopIncorrectFile[0], "../output/occStopIncorrect.txt");
	strcpy(occStopIncorrectFile[1], "../output/occStopIncorrect_PE.txt");
	strcpy(occStopIncorrectFile[2], "../output/occStopIncorrect_SPE.txt");
	strcpy(occStopIncorrectFile[3], "../output/occStopIncorrect_SE.txt");
	strcpy(refPosFile, "../output/refPos.txt");

	strcpy(occPointFile2, "../output/occPoints_cor.txt");
	strcpy(occExtensionCorrectFile2[0], "../output/occExtensionCorrect_cor.txt");
	strcpy(occExtensionCorrectFile2[1], "../output/occExtensionCorrect_cor_PE.txt");
	strcpy(occExtensionCorrectFile2[2], "../output/occExtensionCorrect_cor_SPE.txt");
	strcpy(occExtensionCorrectFile2[3], "../output/occExtensionCorrect_cor_SE.txt");
	strcpy(occExtensionIncorrectFile2[0], "../output/occExtensionIncorrect_cor.txt");
	strcpy(occExtensionIncorrectFile2[1], "../output/occExtensionIncorrect_cor_PE.txt");
	strcpy(occExtensionIncorrectFile2[2], "../output/occExtensionIncorrect_cor_SPE.txt");
	strcpy(occExtensionIncorrectFile2[3], "../output/occExtensionIncorrect_cor_SE.txt");
	strcpy(occStopCorrectFile2[0], "../output/occStopCorrect_cor.txt");
	strcpy(occStopCorrectFile2[1], "../output/occStopCorrect_cor_PE.txt");
	strcpy(occStopCorrectFile2[2], "../output/occStopCorrect_cor_SPE.txt");
	strcpy(occStopCorrectFile2[3], "../output/occStopCorrect_cor_SE.txt");
	strcpy(occStopIncorrectFile2[0], "../output/occStopIncorrect_cor.txt");
	strcpy(occStopIncorrectFile2[1], "../output/occStopIncorrect_cor_PE.txt");
	strcpy(occStopIncorrectFile2[2], "../output/occStopIncorrect_cor_SPE.txt");
	strcpy(occStopIncorrectFile2[3], "../output/occStopIncorrect_cor_SE.txt");
	strcpy(refPosFile2, "../output/refPos.txt");

	occPointFileArr1[0][0] = occPointFile;
	for(i=0; i<4; i++)
	{
		occPointFileArr1[i][1] = occExtensionCorrectFile[i];
		occPointFileArr1[i][2] = occExtensionIncorrectFile[i];
		occPointFileArr1[i][3] = occStopCorrectFile[i];
		occPointFileArr1[i][4] = occStopIncorrectFile[i];
	}
	occPointFileArr1[0][5] = refPosFile;

	occPointFileArr2[0][0] = occPointFile2;
	for(i=0; i<4; i++)
	{
		occPointFileArr2[i][1] = occExtensionCorrectFile2[i];
		occPointFileArr2[i][2] = occExtensionIncorrectFile2[i];
		occPointFileArr2[i][3] = occStopCorrectFile2[i];
		occPointFileArr2[i][4] = occStopIncorrectFile2[i];
	}
	occPointFileArr2[0][5] = refPosFile2;
#endif
	//===================== draw curve variable end ======================

	//===================== SVM model begin ======================
#if (SVM_NAVI==YES)
	strcpy(svmKernelFuncFile[0], "model/PE/svmKF.txt");
	strcpy(svmSupportVectorFile[0], "model/PE/svmSV.txt");
	strcpy(svmAlphaFile[0], "model/PE/svmAlpha.txt");
	strcpy(svmBiasFile[0], "model/PE/svmBias.txt");
	strcpy(svmScaleDataFile[0], "model/PE/svmScaleData.txt");

	strcpy(svmKernelFuncFile[1], "model/SE/svmKF.txt");
	strcpy(svmSupportVectorFile[1], "model/SE/svmSV.txt");
	strcpy(svmAlphaFile[1], "model/SE/svmAlpha.txt");
	strcpy(svmBiasFile[1], "model/SE/svmBias.txt");
	strcpy(svmScaleDataFile[1], "model/SE/svmScaleData.txt");

	if(loadSvmModel(&svmModelPE, svmKernelFuncFile[0], svmSupportVectorFile[0], svmAlphaFile[0], svmBiasFile[0], svmScaleDataFile[0])==FAILED)
	{
		printf("Load SVM model error!\n");
		return ERROR;
	}
	if(loadSvmModel(&svmModelSE, svmKernelFuncFile[1], svmSupportVectorFile[1], svmAlphaFile[1], svmBiasFile[1], svmScaleDataFile[1])==FAILED)
	{
		printf("Load SVM model error!\n");
		return ERROR;
	}

	if(initSampleMemSvm(svmModelPE)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory for sample data when using SVM model, error!\n", __LINE__, __func__);
		return ERROR;
	}
	if(initSampleMemSvm(svmModelSE)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory for sample data when using SVM model, error!\n", __LINE__, __func__);
		return ERROR;
	}

#endif
	//===================== SVM model end ======================

	// output the variables
	if(operationMode==0)
	{
		printf("\nkmer size          : %d\n", kmerSize);
		//printf("read length cutoff : %d\n", readLenCutOff);
		printf("read length        : %d\n", readLen);
		printf("paired mode        : %d\n", pairedMode);
		for(i=0; i<readFileNumPara; i++)
		printf("read files[%d]      : %s\n", i, readFilesInput[i]);
		//printf("single qual thres  : %d\n", singleBaseQualThres);
		if(errorCorrectionFlag==YES)
			printf("correction flag    : %d\n", errorCorrectionFlag);
		if(meanSizeInsert>0)
		{
			printf("insert size        : %.2f\n", meanSizeInsert);
			printf("insert size sdev.  : %.2f\n", standardDev);
		}
		printf("output directory   : %s\n", outputPathStr);
		printf("hash table file    : %s\n", graphFile);
		if(errorCorrectionFlag==YES)
			printf("new hash table file: %s\n", graphFileCorrected);
		printf("contig file        : %s\n", contigsFileFasta);
		printf("minimal contig size: %d\n\n", minContigLen);
	}else if(operationMode==1)
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
		if(errorCorrectionFlag==YES)
			printf("new hash table file: %s\n", graphFileCorrected);
	}else if(operationMode==2)
	{
		printf("\nhash table file    : %s\n", graphFile);
		if(errorCorrectionFlag==YES)
			printf("new hash table file: %s\n", graphFileCorrected);
		printf("kmer size          : %d\n", kmerSize);
		//printf("read length cutoff : %d\n", readLenCutOff);
		printf("read length        : %d\n", readLen);
		printf("paired mode        : %d\n", pairedMode);
		if(errorCorrectionFlag==YES)
			printf("correction flag    : %d\n", errorCorrectionFlag);
		if(meanSizeInsert>0)
		{
			printf("insert size        : %.2f\n", meanSizeInsert);
			printf("insert size sdev.  : %.2f\n", standardDev);
		}
		printf("output directory   : %s\n", outputPathStr);
		printf("contig file        : %s\n", contigsFileFasta);
		printf("minimal contig size: %d\n\n", minContigLen);
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
	if(operationMode==0 || operationMode==1)
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
