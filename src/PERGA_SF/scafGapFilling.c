/*
 * scafGapFilling.c
 *
 *  Created on: Jul 20, 2011
 *      Author: zhuxiao
 */

#include "inc/scafStdinc.h"
#include "inc/scafExtvab.h"

/**
 * Construct graphs and fill the gaps between contigs in ccaffolds.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short gapFilling(const char *graphFile, const char *monoReadSeqFile, const char *monoReadListFile, const char *newContigOverlapInfoFile, const char *newContigFile, const char *readListFile1, const char *readListFile2, const char *readSeqFile, const char *meanSdevFile, const char *contigOverlapInfoFile, const char *contigFile)
{
	printf("=========== Begin filling gaps between contigs, please wait ...\n");

	if(initGapFillingParas()==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the parameters for gap filling, error!\n", __LINE__, __func__);
		return 1;
	}

	if(buildGraph(graphFile, monoReadSeqFile, monoReadListFile, readListFile1, readListFile2, readSeqFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot build graph for gap filling, error!\n", __LINE__, __func__);
		return 1;
	}

	if(fillGaps(newContigOverlapInfoFile, newContigFile, meanSdevFile, contigOverlapInfoFile, contigFile, monoReadListFile, graphFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill gaps, error!\n", __LINE__, __func__);
		return 1;
	}

	freeGapFillingParas();

	printf("=========== End filled gaps between contigs.\n");

	return SUCCESSFUL;
}

/**
 * Fill the gaps between contigs in ccaffolds.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillGaps(const char *newContigOverlapInfoFile, const char *newContigFile, const char *meanSdevFile, const char *contigOverlapInfoFile, const char *contigFile, const char *monoReadListFile, const char *graphFile)
{
	// initialize the memory for gap filling
	if(initMemGapFilling(contigOverlapInfoFile, meanSdevFile, contigFile, monoReadListFile, graphFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory for gap filling, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ########################## Debug information #######################
#if DEBUG_FLAG
	if(checkScafKmersInGrapInScaf(monoReadSeqFile, scafGrapDeBruijn)==FAILED)
	{
		printf("line=%d, In %s(), there are some errors in scafGrap in local assembly, error!\n", __LINE__, __func__);
		return FAILED;
	}
#endif
	// ########################## Debug information #######################


	// print the total number of gap regions, and their total length
	printf("Before filling gaps:\n");
	if(outputGapSizeInScaf(contigOverlapIndexArr, scaffoldsNumInCOI, contigOverlapInfoArr, minBaseNumInGap)==FAILED)
	{
		printf("line=%d, In %s(), cannot output the gap size, error!\n", __LINE__, __func__);
		return FAILED;
	}


	// fill the gaps in scaffolds by local assembly
	if(localAssemblyInScaf()==FAILED)
	{
		printf("line=%d, In %s(), cannot fill gaps, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// print the total number of gap regions, and their total length
	printf("After filled gaps:\n");
	if(outputGapSizeInScaf(contigOverlapIndexArr, scaffoldsNumInCOI, contigOverlapInfoArr, minBaseNumInGap)==FAILED)
	{
		printf("line=%d, In %s(), cannot output the gap size, error!\n", __LINE__, __func__);
		return FAILED;
	}


	// output the contigs
	if(outputContigInfoArrToFile(newContigFile, contigInfoArr, contigsNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot output the contig information array to file [ %s ], error!\n", __LINE__, __func__, newContigFile);
		return FAILED;
	}

	// output updated contig overlap information to text file.
	//  Format:
	//  	(1) Header fields: >scaffoldID, linkedNum, rowsNum, which are separated by tab character;
	//  	(2)   Body fields: contigID1, orientation1, contigLen1, contigID2, orientation2, contigLen2, mergeFlag, overlapLen, gapSize, breakFlag, which are separated by tab character.
	if(outputContigOverlapInfoToFile(newContigOverlapInfoFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot output the contig overlap information to file [%s ], error!\n", __LINE__, __func__, newContigOverlapInfoFile);
		return FAILED;
	}


	freeMemGapFilling();

	return SUCCESSFUL;
}

/**
 * Initialize the global parameters for gap filling.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initGapFillingParas()
{
	//PEGivenType = SE_GIVEN_TYPE;			// 0

	PEGivenType = NONE_PE_GIVEN_TYPE;		// 1

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
	kmerSize = 21;
	//kmerSize = 23;
	//kmerSize = 25;
	//kmerSize = 27;
	//kmerSize = 29;
	//kmerSize = 31;

	errorRegLenEnd3 = ceil(readLen * ERROR_REGION_LEN_3End_FACTOR);
	if(kmerSize%2==0)
		kmerSize --;

	if(kmerSize>readLen-errorRegLenEnd3)
	{
		kmerSize = readLen - errorRegLenEnd3;
		if(kmerSize%2==0)
			kmerSize --;
	}

	entriesPerKmer = ((kmerSize-1) / 32) + 1;
	if(kmerSize%32==0)
	{
		lastEntryMaskKmer = (uint64_t) -1;
		lastEntryBaseNumKmer = 32;
	}else
	{
		lastEntryMaskKmer = (1LLU << ((kmerSize%32)<<1)) - 1;
		lastEntryBaseNumKmer = kmerSize % 32;
	}

	hashTableSize = 15728661;
	//hashTableSize = 1LLU << (kmerSize << 1);

	//printf("PEGivenType=%d\n", PEGivenType);
	if(PEGivenType==INSERT_PE_GIVEN_TYPE)
	{
		printf("meanSizeInsert=%.2f\n", meanSizeInsert);
	}else if(PEGivenType==BOTH_PE_GIVEN_TYPE)
	{
		printf("meanSizeInsert=%.2f\n", meanSizeInsert);
		printf("standardDev=%.2f\n", stardardDeviationInsert);
	}

	//printf("kmerSize=%d\n", kmerSize);
	//printf("entriesPerKmer=%d\n", entriesPerKmer);
	//printf("lastEntryBaseNumKmer=%d\n", lastEntryBaseNumKmer);
	//printf("lastEntryMaskKmer=0x%lX\n", lastEntryMaskKmer);
	//printf("hashTableSize=%lu\n", hashTableSize);
	//printf("errorRegLenEnd3=%d\n", errorRegLenEnd3);


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

	return SUCCESSFUL;
}

/**
 * Free the global parameters for gap filling.
 */
void freeGapFillingParas()
{
	free(kmerSeqInt);
	kmerSeqInt = NULL;
	free(kmerSeqIntRev);
	kmerSeqIntRev = NULL;
}

/**
 * Initialize the memory for gap filling.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemGapFilling(const char *contigOverlapInfoFile, const char *meanSdevFile, const char *contigFile, const char *monoReadListFile, const char *graphFile)
{
	int i;

	// initialize the contig information array
	if(initContigInfoArray(contigFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the contig information array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// load the contig overlap information
	if(loadContigOverlapInfo(contigOverlapInfoFile, &contigOverlapIndexArr, &scaffoldsNumInCOI, &contigOverlapInfoArr, &itemNumInContigOverlapArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the contig overlap information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Allocate memory and load the data of mono read list (MRL)
	if(loadSingleReadList(monoReadListFile, &monoReadListArr, &readItemNumInMRL, &monoReadPosArr, &matchItemNumInMRP)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the data of mono read list (MRL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// load the graph in scaffolding
	if(loadGraphInScaf(&scafGrapDeBruijn, graphFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the scafGraph, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// compute the averKmerOcc
	if(computeAverKmerOcc(&averKmerOcc, scafGrapDeBruijn)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute averKmerOcc, error!\n", __LINE__, __func__);
		return FAILED;
	}


	longKmerSize = ceil((readLen - 2*errorRegLenEnd3 - kmerSize) * LONG_KMER_SIZE_FACTOR) + kmerSize;
	if((longKmerSize & 1) == 0)
		longKmerSize --;
	longKmerStepSize = ceil((longKmerSize - kmerSize) / 3.0);
	if((longKmerStepSize & 1) == 1)
		longKmerStepSize ++;
	if(longKmerStepSize<1)
		longKmerStepSize = 2;

	//printf("longKmerSize=%d, longKmerStepSize=%d\n", longKmerSize, longKmerStepSize);


	if(PEGivenType>=NONE_PE_GIVEN_TYPE)
	{
		minKmerOccPE = ceil(averKmerOcc * MIN_KMER_OCC_FACTOR);
		minKmerOccSE = ceil(averKmerOcc * MIN_KMER_OCC_FACTOR);

		if(minKmerOccPE<MIN_KMER_OCC_THRES)
			minKmerOccPE = MIN_KMER_OCC_THRES;
		if(minKmerOccSE<MIN_KMER_OCC_THRES)
			minKmerOccSE = MIN_KMER_OCC_THRES;

		//maxSecondOcc = minKmerOccSE * OCCS_NUM_FACTOR;
		maxSecondOcc = ceil(minKmerOccSE * MAX_SECOND_OCC_FACTOR);
		maxFirstOcc = ceil(minKmerOccSE * MAX_FIRST_OCC_FACTOR);
		minLongKmerOcc = floor(minKmerOccSE * LONG_KMER_OCC_FACTOR);
		//minReadsNumPEHashThres = ceil(averKmerOcc * MIN_READ_NUM_PE_HASH_FACTOR);

		if(maxSecondOcc>MAX_SECOND_OCC_THRES)
		{
			maxSecondOcc = MAX_SECOND_OCC_THRES;
		}
		if(minLongKmerOcc>MIN_LONG_KMER_OCC_THRES)
		{
			minLongKmerOcc = MIN_LONG_KMER_OCC_THRES;
		}


		maxOccNumFaiedPE = ceil(OCCS_NUM_SE_FAILED_PE_FACTOR * averKmerOcc);
		if(maxOccNumFaiedPE>MAX_OCC_NUM_FAILED_PE_THRES)
			maxOccNumFaiedPE = MAX_OCC_NUM_FAILED_PE_THRES;
//		else
//			maxOccNumFaiedPE *= 0.8;
		//maxOccNumFaiedPE = MAX_OCC_NUM_FAILED_PE_THRES;
		maxNavigationNumSE = MAX_NAVI_NUM_SE_THRES;

		//printf("minKmerOccSE=%.2f, minKmerOccPE=%.2f\n", minKmerOccSE, minKmerOccPE);
		//printf("maxSecondOcc=%.2f\n", maxSecondOcc);
		//printf("maxFirstOcc=%.2f\n", maxFirstOcc);
		//printf("minLongKmerOcc=%.2f\n", minLongKmerOcc);
		//printf("minReadsNumPEHashThres=%.2f\n", minReadsNumPEHashThres);
	//	printf("maxOccNumFaiedPE=%.2f\n", maxOccNumFaiedPE);
		//printf("maxNavigationNumSE=%d\n", maxNavigationNumSE);
	}else
	{
		minKmerOccSE = ceil(averKmerOcc * MIN_KMER_OCC_FACTOR);

		if(minKmerOccSE<MIN_KMER_OCC_THRES)
			minKmerOccSE = MIN_KMER_OCC_THRES;

		//maxSecondOcc = minKmerOccSE * OCCS_NUM_FACTOR;
		maxSecondOcc = ceil(minKmerOccSE * MAX_SECOND_OCC_FACTOR);
		maxFirstOcc = ceil(minKmerOccSE * MAX_FIRST_OCC_FACTOR);
		minLongKmerOcc = floor(minKmerOccSE * LONG_KMER_OCC_FACTOR);

		if(maxSecondOcc>MAX_SECOND_OCC_THRES)
		{
			maxSecondOcc = MAX_SECOND_OCC_THRES;
		}
		if(minLongKmerOcc>MIN_LONG_KMER_OCC_THRES)
		{
			minLongKmerOcc = MIN_LONG_KMER_OCC_THRES;
		}

		//printf("minKmerOccSE=%.2f\n", minKmerOccSE);
		//printf("maxSecondOcc=%.2f\n", maxSecondOcc);
		//printf("maxFirstOcc=%.2f\n", maxFirstOcc);
		//printf("minLongKmerOcc=%.2f\n", minLongKmerOcc);
	}

	lockedReadsNumThres = averKmerOcc;
	//printf("lockedReadsNumThres=%.2f\n", lockedReadsNumThres);

	// initialize the variables
	maxScafAssemblingReadsNumInArr = TABLE_SIZE_SCAFASSEMBLINGREAD;
	maxScafSuccessReadsNumInArr = TABLE_SIZE_SUCCESSFUL_READ_ARRAY_SCAF;
	//maxComparisonSeqLenInScaf = 2 * readLen;
	//maxScafContigEndSeqLen = 2 * readLen;
	standardDeviationFactor = SDEV_FACTOR;

	maxOverlapSeqLen = 2 * readLen;
	minOverlapThres = MIN_OVERLAP_THRESHOLD;
	minExactOverlapThres = MIN_EXACT_OVERLAP_THRESHOLD;
	mismatchThres = MAX_MISMATCH_THRESHOLD;
	minBaseNumInGap = MIN_BASENUM_IN_GAP;
	gapSizeSdevFactorGapFilling = GAP_SIZE_SDEV_FACTOR_GAPFILLING;
	//exactOverlapSdevThres = EXACT_OVERLAP_SDEV_THRESHOLD;
	minAdjustGapSizeThres = MIN_ADJUST_GAP_SIZE_THRESHOLD;
	maxAdjustGapSizeThres = MAX_ADJUST_GAP_SIZE_THRESHOLD;
	matchScore = MATCH_SCORE;
	mismatchScore = MISMATCH_SCORE;
	gapScore = GAP_SCORE;

	// load mean size and standard deviation of fragments
	if(loadMeanSizeAndSDev(meanSdevFile, &meanSizeInsert, &stardardDeviationInsert)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the mean size and standard deviation of fragments, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate the memory for scafAssemblingReadArr
	scafAssemblingReadArr = (scafAssemblingRead *) calloc(maxScafAssemblingReadsNumInArr, sizeof(scafAssemblingRead));
	if(scafAssemblingReadArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate the memory for scafSuccessReadArr
	scafSuccessReadArr = (scafSuccessRead *) calloc(maxScafSuccessReadsNumInArr, sizeof(scafSuccessRead));
	if(scafSuccessReadArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate the memory for sacfContigSeqLastReadLen
	scafContigSeqLastReadLen = (char *) calloc(readLen+1, sizeof(char));
	if(scafContigSeqLastReadLen==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate the memory for *contigEndSeqInScaf[2]
	for(i=0; i<2; i++)
	{
		scafContigEndSeqArr[i] = (char *) calloc((maxOverlapSeqLen+1)*2, sizeof(char));
		if(scafContigEndSeqArr[i]==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// allocate the memory for comparisonSeqInScaf
	comparisonSeqInScaf = (char *) calloc((maxOverlapSeqLen+1)*2, sizeof(char));
	if(comparisonSeqInScaf==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate the memory for exact sequence comparison and alignment
	scoreArr = (int *) calloc((maxOverlapSeqLen+1)*(maxOverlapSeqLen+1), sizeof(int));
	if(scoreArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<3; i++)
	{
		alignResultArr[i] = (char *) calloc (2*(maxOverlapSeqLen+1), sizeof(char));
		if(alignResultArr[i]==NULL)
		{
			printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	//+++++++++++++++++++++++++++++++++++++++
	kmerSeqIntAssembly = (uint64_t*) malloc(entriesPerKmer*sizeof(uint64_t));
	if(kmerSeqIntAssembly==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	kmerSeqIntAssemblyRev = (uint64_t*) malloc(entriesPerKmer*sizeof(uint64_t));
	if(kmerSeqIntAssemblyRev==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	tmpKmerSeqIntAssembly = (uint64_t*) malloc(entriesPerKmer*sizeof(uint64_t));
	if(tmpKmerSeqIntAssembly==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<4; i++) tabooSeqInt[i] = 0;
	for(i=0; i<32; i++)
	{
		tabooSeqInt[1] = (tabooSeqInt[1] << 2) | 1;
		tabooSeqInt[2] = (tabooSeqInt[2] << 2) | 2;
		tabooSeqInt[3] = (tabooSeqInt[3] << 2) | 3;
	}

	return SUCCESSFUL;
}

/**
 * Free the memory for gap filling.
 */
void freeMemGapFilling()
{
	int i;

	// free contig information array
	freeMemContigInfo(&contigInfoArr, &contigsNum);

	// free the contig overlap information
	freeContigOverlapInfo(&contigOverlapIndexArr, &scaffoldsNumInCOI, &contigOverlapInfoArr, &itemNumInContigOverlapArr);

	// Release the memory of mono read list (MRL)
	freeSingleReadList(&monoReadListArr, &readItemNumInMRL, &monoReadPosArr, &matchItemNumInMRP);

	// Release the memory of scafGraph
	freeGraphInScaf(&scafGrapDeBruijn);

	// reset the variables
	maxScafAssemblingReadsNumInArr = 0;
	maxScafSuccessReadsNumInArr = 0;
	//maxScafContigEndSeqLen = 0;
	//maxComparisonSeqLenInScaf = 0;
	meanSizeInsert = 0;
	stardardDeviationInsert = 0;
	standardDeviationFactor = 0;

	maxOverlapSeqLen = 0;
	minOverlapThres = 0;
	minExactOverlapThres = 0;
	mismatchThres = 0;
	minBaseNumInGap = 0;
	gapSizeSdevFactorGapFilling = 0;
	//exactOverlapSdevThres = 0;
	minAdjustGapSizeThres = 0;
	maxAdjustGapSizeThres = 0;
	matchScore = 0;
	mismatchScore = 0;
	gapScore = 0;

	// free memory of scafAssemblingReadArr
	free(scafAssemblingReadArr);
	scafAssemblingReadArr = NULL;

	// free memory of scafSuccessReadArr
	free(scafSuccessReadArr);
	scafSuccessReadArr = NULL;

	// free memory of sacfContigSeqLastReadLen
	free(scafContigSeqLastReadLen);
	scafContigSeqLastReadLen = NULL;

	// free memory of *contigEndSeqInScaf[2]
	for(i=0; i<2; i++)
	{
		free(scafContigEndSeqArr[i]);
		scafContigEndSeqArr[i] = NULL;
	}

	// free memory of comparisonSeqInScaf
	free(comparisonSeqInScaf);
	comparisonSeqInScaf = NULL;

	// free memory of scoreArr
	free(scoreArr);
	scoreArr = NULL;

	// free memory of alignResultArr[3]
	for(i=0; i<3; i++)
	{
		free(alignResultArr[i]);
		alignResultArr[i] = NULL;
	}

	//+++++++++++++++++++++++++++++++++
	free(kmerSeqIntAssembly);
	kmerSeqIntAssembly = NULL;

	free(kmerSeqIntAssemblyRev);
	kmerSeqIntAssemblyRev = NULL;

	free(tmpKmerSeqIntAssembly);
	tmpKmerSeqIntAssembly = NULL;

}

/**
 * Compute the averKmerOcc.
 */
short computeAverKmerOcc(double *averKmerOcc, scafGraph *pScafGrapDeBruijn)
{
	int64_t i, filledBucketNum;
	uint64_t tmp_kmerSum, tmp_kmerNum;
	scafKmer *kmer;

	filledBucketNum = 0;
	tmp_kmerSum = tmp_kmerNum = 0;
	for(i=0; i<hashTableSize; i++)
	{
		kmer = pScafGrapDeBruijn->pkmers[i];
		if(kmer)
			filledBucketNum ++;

		while(kmer)
		{
			tmp_kmerSum += kmer->arraysize;
			tmp_kmerNum ++;

			kmer = kmer->next;
		}
	}

	if(tmp_kmerNum>0)
	{
		*averKmerOcc = (double)tmp_kmerSum / tmp_kmerNum;
		if(*averKmerOcc<1)
			*averKmerOcc = 1;
	}else
	{
		printf("line=%d, In %s(), kmerNum=%lu, error!\n", __LINE__, __func__, tmp_kmerNum);
		return FAILED;
	}

	//printf("filledBucketNum=%ld, filledRatio=%.2f, kmerNum=%lu\n", filledBucketNum, (double)filledBucketNum/hashTableSize, tmp_kmerNum);
	//printf("averKmerOcc=%.2f.\n", *averKmerOcc);

	return SUCCESSFUL;
}

/**
 * Fill gaps between contigs in scaffolds by local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short localAssemblyInScaf()
{
	contigOverlap *pContigOverlapInfo;
	int scaffoldID, contigID[2], contigOrient[2], contigLen[2], gapSize, newGapSize;
	int j, k, startRowInCOI, rowsNumInCOI, linkedContigsNum;
	int assemblyRound, maxAssemblyLen, minAssemblyLen;	// assemblyRound: 0 -- first round; 1 -- second round
	scafContig *localScafContigheadArr[2];
	int localScafContigNodesNumArr[2];
	int overlapLen, oldEndSeqLenArray[2];
	int breakFlag;


	// fill gaps between contigs in scaffolds by local assembly
	for(scaffoldID=1; scaffoldID<=scaffoldsNumInCOI; scaffoldID++)
	{
		// ####################### Debug information ####################
#if DEBUG_OUT_FLAG
		printf("================== scaffoldID=%d ==============\n", scaffoldID);
//		if(scaffoldID==197)
//		{
//			printf("scaffoldID=%d\n", scaffoldID);
//		}
#endif
		// ####################### Debug information ####################

		linkedContigsNum = contigOverlapIndexArr[scaffoldID-1].linkedNum;
		if(linkedContigsNum==1)
		{
			continue;
		}else if(linkedContigsNum<=0)
		{
			printf("line=%d, In %s(), scaffoldID=%d, linkedContigsNum=%d, error!\n", __LINE__, __func__, scaffoldID, linkedContigsNum);
			return FAILED;
		}

		rowsNumInCOI = contigOverlapIndexArr[scaffoldID-1].rowsNum;
		startRowInCOI = contigOverlapIndexArr[scaffoldID-1].startRow;
		pContigOverlapInfo = contigOverlapInfoArr + startRowInCOI;
		for(j=0; j<rowsNumInCOI; j++)
		{
			if(pContigOverlapInfo[j].mergeFlag==NO)
			{ // a gap between the two contigs

				// get the information of the two contigs
				contigID[0] = pContigOverlapInfo[j].contigID1;
				contigID[1] = pContigOverlapInfo[j].contigID2;
				contigOrient[0] = pContigOverlapInfo[j].orientation1;
				contigOrient[1] = pContigOverlapInfo[j].orientation2;
				contigLen[0] = contigInfoArr[contigID[0]-1].contigLen;
				contigLen[1] = contigInfoArr[contigID[1]-1].contigLen;
				gapSize = pContigOverlapInfo[j].gapSize;

				localScafContigheadArr[0] = localScafContigheadArr[1] = NULL;
				localScafContigNodesNumArr[0] = localScafContigNodesNumArr[1] = 0;

				//###################### Debug information #####################
#if DEBUG_FLAG
				if(contigID[0]==42 && contigID[1]==52)
				//if(contigID[0]==602 && contigID[1]==2445)
				{
					printf("line=%d, In %s(), contigID1=%d, contigOrient1=%d, contigLen1=%d, contigID2=%d, contigOrient2=%d, contigLen2=%d, gapSize=%d\n", __LINE__, __func__, contigID[0], contigOrient[0], contigLen[0], contigID[1], contigOrient[1], contigLen[1], gapSize);
				}
				printf("line=%d, In %s(), contigID1=%d, contigOrient1=%d, contigLen1=%d, contigID2=%d, contigOrient2=%d, contigLen2=%d, gapSize=%d\n", __LINE__, __func__, contigID[0], contigOrient[0], contigLen[0], contigID[1], contigOrient[1], contigLen[1], gapSize);
#endif
				//###################### Debug information #####################

				// get the two contig ends and the ends length
				if(getScafContigEndSeqs(pContigOverlapInfo+j)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the scafContig end sequences of contigs [ %d, %d ], error!\n", __LINE__, __func__, contigID[0], contigID[1]);
					return FAILED;
				}

				// assembly two contig ends
				for(assemblyRound=FIRST_ROUND_ASSEMBLY; assemblyRound<=SECOND_ROUND_ASSEMBLY; assemblyRound++)
				{
					// initialize the variables
					scafContighead = NULL;
					scafContigtail = NULL;
					scafContigLastReadLen = NULL;
					scafContigIndex = 0;
					scafContigSeqLastReadLen[0] = '\0';
					scafContigSeqLenLastReadLen = 0;
					scafAssemblingReadsNum = 0;
					successScafContig = NULL;
					lockedReadsNumInScaf = 0;
					successFilledFlag = NO;
					breakFlag = NO;

					// record the ends length before local assembly
					oldEndSeqLenArray[0] = scafContigEndSeqLenArr[0];
					oldEndSeqLenArray[1] = scafContigEndSeqLenArr[1];

					// get the comparison sequence
					if(getComparisonSeqInScaf(assemblyRound)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the comparison sequence, error!\n", __LINE__, __func__);
						return FAILED;
					}

					// prepare the first scafContig nodes with length of READ_LEN before local assembly
					if(prepareAssemblyInScaf(contigID[assemblyRound], contigLen[assemblyRound], contigOrient[assemblyRound], assemblyRound)==FAILED)
					{
						printf("line=%d, In %s(), cannot prepare the assembly of contigs [ %d, %d ], error!\n", __LINE__, __func__, contigID[0], contigID[1]);
						return FAILED;
					}

					if(getMaxMinAssemblyLen(&maxAssemblyLen, &minAssemblyLen, assemblyRound, localScafContigNodesNumArr, gapSize)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the maximal local assembly length, error!\n", __LINE__, __func__);
						return FAILED;
					}

					// local assembly
					if(scafKmers[0] || scafKmers[1])
					{
						while(scafKmers[0] || scafKmers[1])
						{
							// get next kmers in local assembly
							if(PEGivenType>NONE_PE_GIVEN_TYPE)
							{
								//取正反向kmer
								if(getNextKmerByMixInScaf(scafContigIndex, contigID[assemblyRound], contigLen[assemblyRound], contigOrient[assemblyRound], assemblyRound)==FAILED)
								{
									printf("line=%d, In %s(), cannot get the next kmer by mix, error!\n", __LINE__, __func__);
									return FAILED;
								}
							}else
							{
								//取正反向kmer
								if(getNextKmerBySEInScaf(scafContigIndex, contigID[assemblyRound], contigLen[assemblyRound], contigOrient[assemblyRound], assemblyRound)==FAILED)
								{
									printf("line=%d, In %s(), cannot get next kmer, error!\n", __LINE__, __func__);
									return FAILED;
								}
							}

							if(scafKmers[0]==NULL && scafKmers[1]==NULL)
							{
								break;
							}

							scafContigIndex ++;

							// append the base to scafContig chain
							if(appendContigBaseInScaf(kmerSeqIntAssembly[entriesPerKmer-1] & 3, scafContigIndex)==FAILED)
							{
								printf("line=%d, In %s(), cannot add base to scafContig chain, error!\n", __LINE__, __func__);
								return FAILED;
							}

							// update the decision table
							if(updateScafAssemblingReadsInScaf(contigID[assemblyRound], contigLen[assemblyRound], contigOrient[assemblyRound], scafContigIndex, assemblyRound)==FAILED)
							{
								printf("line=%d, In %s(), cannot update the decision table, error!\n", __LINE__, __func__);
								return FAILED;
							}

							// update reads status in decision table
							if(updateAssemblingreadsStatusInScaf()==FAILED)
							{
								printf("line=%d, In %s(), cannot update the status of reads in decision table, error!\n", __LINE__, __func__);
								return FAILED;
							}

							// update the locked reads and their total number
							updateLockedReadsInScaf();

							// update finished reads in decision table
							if(updateFinisedScafReadsInScaf()==FAILED)
							{
								printf("line=%d, In %s(), cannot update the finished reads in decision table, error!\n", __LINE__, __func__);
								return FAILED;
							}

							// update the sequence with a length of READ_LEN at the end of the contig nodes
							if(scafContigIndex<=scafContigSeqLenLastReadLen)
							{ //如果contig的节点数目<=36，则contig36指向contighead节点，并将当次拼接kmer的序列的最后一个碱基添加进lastseq36数组
								scafContigLastReadLen = scafContighead;
								switch(kmerSeqIntAssembly[entriesPerKmer-1] & 3)
								{
									case 0: scafContigSeqLastReadLen[scafContigIndex-1] = 'A'; break;
									case 1: scafContigSeqLastReadLen[scafContigIndex-1] = 'C'; break;
									case 2: scafContigSeqLastReadLen[scafContigIndex-1] = 'G'; break;
									case 3: scafContigSeqLastReadLen[scafContigIndex-1] = 'T'; break;
								}
								scafContigSeqLastReadLen[scafContigIndex] = '\0';

							}else
							{  //如果contig中节点数目>36，则contig36后移一个位置，并将lastseq36数组的碱基前移一位，并将当次拼接的kmer的序列添加进lastseq36数组
								scafContigLastReadLen = scafContigLastReadLen->next;
								for(k=1; k<scafContigSeqLenLastReadLen; k++) scafContigSeqLastReadLen[k-1] = scafContigSeqLastReadLen[k];
								switch(kmerSeqIntAssembly[entriesPerKmer-1] & 3)
								{
									case 0: scafContigSeqLastReadLen[scafContigSeqLenLastReadLen-1] = 'A'; break;
									case 1: scafContigSeqLastReadLen[scafContigSeqLenLastReadLen-1] = 'C'; break;
									case 2: scafContigSeqLastReadLen[scafContigSeqLenLastReadLen-1] = 'G'; break;
									case 3: scafContigSeqLastReadLen[scafContigSeqLenLastReadLen-1] = 'T'; break;
								}
								scafContigSeqLastReadLen[scafContigSeqLenLastReadLen] = '\0';
							}

							// update the scafGraph
							scafContig *tmp_successContig = successScafContig;
							if(scafSuccessReadsNum>0)
							{
								// delete the reads from scafGraph
								if(delReadsFromGraphInScaf()==FAILED)
								{
									printf("line=%d, In %s(), contigID=%d, scafContigIndex=%d, cannot delete the reads from graph, error!\n", __LINE__, __func__, contigID[assemblyRound], scafContigIndex);
									outputSuccessReads(scafSuccessReadArr, scafSuccessReadsNum);
									return FAILED;
								}

								// add the scafRidpos information to scafContig nodes
								if(addRidposToContigInScaf(scafContigIndex)==FAILED)
								{
									printf("line=%d, In %s(), contigID=%d, cannot add scafRidpos information to scafContig nodes, error!\n", __LINE__, __func__, contigID[assemblyRound]);
									return FAILED;
								}

								// get the successful scafContig node in local assembly
								tmp_successContig = getSuccessContigInScaf(scafContigIndex);

								// update the successful scafContig node pointer
								if(tmp_successContig==NULL)
								{
									//printf("line=%d, In %s(), contigsNum=%d, assemblyCycleNum=%d, contigIndex=%d, assemblyStatus=%d, previousAssemblyStatus=%d, the tmp_successContig==NULL!\n", __LINE__, __func__, contigsNum+1, assemblyCycleNum, contigIndex, assemblyStatus, previousAssemblyStatus);
									//outputContig(contig36);
									//return FAILED;
								}else
								{
									successScafContig = tmp_successContig;
								}
							}


							// check the assembly
							//处理死循环
							if(successScafContig==NULL)
							{ //还未有成功的reads, 并且拼接长度已经超过60, 则该contig拼接失败, 退出
								if(scafContigIndex>2*readLen-kmerSize)
								{
									//printf("line=%d, In %s(), contigID=%d, scafcontigIndex=%d, successContig==NULL!\n", __LINE__, __func__, contigID[assemblyRound], scafContigIndex);
									break;
								}
							}else if(tmp_successContig==NULL)
							{
								break;
							}else if(scafContigIndex-successScafContig->index > readLen-MIN_OVERLAP_LEN)
							{ //已经有成功的reads, 则根据拼接的情况, 确定是否需要继续拼接
								break;
							}
							//===============================================================
							else if(scafContigIndex>maxAssemblyLen)
							{ // get the maximal assembly length and terminate
								//break;
							}
							//===============================================================
						} // end while(scafKmers[0] || scafKmers[1])
					} // end if(scafKmers[0] || scafKmers[1])

					// after the local assembly, trim the uncovered bases at 3' ends
					if(successScafContig)
					{
						// 将contig节点退回到最近成功的contigIndex的位置, 并将之后的contig节点删掉.
						if(updateContigtailnodesInScaf(&scafContigIndex)==FAILED)
						{
							printf("line=%d, In %s(), contigID=%d, contigIndex=%d, cannot update Contigtail nodes, Error!\n", __LINE__, __func__, contigID[assemblyRound], scafContigIndex);
							return FAILED;
						}

						//######################### Debug information #####################
#if DEBUG_OUT_FLAG
						printf("Before updateScafContigEndSeqs(), assemblyRound=%d, contigEndSeq[%d]=%s, len=%d\n", assemblyRound, assemblyRound, scafContigEndSeqArr[assemblyRound], scafContigEndSeqLenArr[assemblyRound]);
#endif
						//######################### Debug information #####################

						// update the contig end sequences
						if(updateScafContigEndSeqs(assemblyRound, scafContigIndex)==FAILED)
						{
							printf("line=%d, In %s(), cannot update the scafContig end sequences, error!\n", __LINE__, __func__);
							return FAILED;
						}

						// adjust the scafContig nodes number and the oldEndSeqLenArr[0,1]
						localScafContigNodesNumArr[assemblyRound] = scafContigIndex;
						oldEndSeqLenArray[assemblyRound] = scafContigEndSeqLenArr[assemblyRound];

						//######################### Debug information #####################
#if DEBUG_OUT_FLAG
						printf("After updateScafContigEndSeqs(), assemblyRound=%d, contigEndSeq[%d]=%s, len=%d\n", assemblyRound, assemblyRound, scafContigEndSeqArr[assemblyRound], scafContigEndSeqLenArr[assemblyRound]);
#endif
						//######################### Debug information #####################

						if(scafContigIndex >= minAssemblyLen)
						//if(scafContigIndex >= readLen)  // 2012-11-19
						{ // the contig nodes number is more than the minimal assembly length
							// detect overlaps of contig ends in local assembly
							if(detectOverlapsInScaf(&successFilledFlag, &overlapLen, &newGapSize, &breakFlag, gapSize, localScafContigNodesNumArr, assemblyRound)==FAILED)
							{
								printf("line=%d, In %s(), cannot detect overlaps between contig ends, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}
					}

					//######################### Debug information #####################
#if DEBUG_OUT_FLAG
					printf("assemblyRound=%d, contigSeq=", assemblyRound);
					outputScafContigSequence(scafContighead);
#endif
					//######################### Debug information #####################

					localScafContigheadArr[assemblyRound] = scafContighead;
					scafContighead = NULL;
					localScafContigNodesNumArr[assemblyRound] = scafContigIndex;

					// control sentence
					if(successFilledFlag==YES)
						break;

				} // end for(assemblyRound=FIRST_ROUND_ASSEMBLY; assemblyRound<=SECOND_ROUND_ASSEMBLY; assemblyRound++)


				// check filled status, and update the information in local assembly
				if(successFilledFlag==YES)
				{ // successfully fill the gap, then update the overlap information, contig information, ...

					if(updateContigOverlapInfoInScaf(pContigOverlapInfo+j, successFilledFlag, overlapLen, newGapSize, breakFlag)==FAILED)
					{
						printf("line=%d, In %s(), cannot update contig overlap information, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(updateContigInfoInScaf(pContigOverlapInfo+j, localScafContigheadArr, localScafContigNodesNumArr, oldEndSeqLenArray, assemblyRound)==FAILED) ////////////////////////////////////////////////////////
					{
						printf("line=%d, In %s(), cannot update contig information, error!\n", __LINE__, __func__);
						return FAILED;
					}

					//===========================================================
					// further update RL1, RL2, SRL, MRL, CI, CL

					//===========================================================
				}else
				{ // failed fill the gap, update the overlap informatipon and the contig information
					// ######################### Debug information ########################
#if DEBUG_OUT_FLAG
					printf("scaffoldID=%d, contigID1=%d, contigID2=%d, contigOrient1=%d, contigOrient2=%d, assemblyRound=%d, gapSize=%d\n", scaffoldID, contigID[0], contigID[1], contigOrient[0], contigOrient[1], assemblyRound, gapSize);
					printf("contigNodesNum1=%d, contigNodesNum2=%d, prepareAssemblyLen1=%d, prepareAssemblyLen2=%d\n", localScafContigNodesNumArr[0], localScafContigNodesNumArr[1], prepareAssemblyLenArr[0], prepareAssemblyLenArr[1]);
#endif
					// ######################### Debug information ########################

					//if(scafContigIndex >= minAssemblyLen)
					//if(localScafContigNodesNumArr[0]+localScafContigNodesNumArr[1] > readLen)  // 2012-11-19
					{ // the contig nodes number is more than the minimal assembly length
						// detect overlaps of contig ends in local assembly
						if(detectOverlapsInScaf(&successFilledFlag, &overlapLen, &newGapSize, &breakFlag, gapSize, localScafContigNodesNumArr, assemblyRound)==FAILED)
						{
							printf("line=%d, In %s(), cannot detect overlaps between contig ends, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}

					// To do ...
					if(localScafContigNodesNumArr[0] > prepareAssemblyLenArr[0] || localScafContigNodesNumArr[1] > prepareAssemblyLenArr[1])
					{ // new bases are appended
//						if(scafContigIndex < minAssemblyLen)
//						{
//							// compute the new gap size
//							if(computeNewGapSizeInScaf(&newGapSize, gapSize, localScafContigNodesNumArr, assemblyRound)==FAILED)
//							{
//								printf("line=%d, In %s(), cannot compute new gap size, error!\n", __LINE__, __func__);
//								return FAILED;
//							}
//						}

						// update overlap information
						if(updateContigOverlapInfoInScaf(pContigOverlapInfo+j, successFilledFlag, overlapLen, newGapSize, breakFlag)==FAILED)
						{
							printf("line=%d, In %s(), cannot update contig overlap information, error!\n", __LINE__, __func__);
							return FAILED;
						}

						// update contig information
						if(updateContigInfoInScaf(pContigOverlapInfo+j, localScafContigheadArr, localScafContigNodesNumArr, oldEndSeqLenArray, assemblyRound)==FAILED) ////////////////////////////////////////////////////////
						{
							printf("line=%d, In %s(), cannot update contig information, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
				}

				// ######################## Debug information ##########################
#if DEBUG_OUT_FLAG
				if(localScafContigheadArr[0])
				{
					printf("contigSeq[0]: ");
					if(outputScafContigSequence(localScafContigheadArr[0])==FAILED)
					{
						printf("line=%d, In %s(), cannot get the scafContig sequences, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
				if(localScafContigheadArr[1])
				{
					printf("contigSeq[1]: ");
					if(outputScafContigSequence(localScafContigheadArr[1])==FAILED)
					{
						printf("line=%d, In %s(), cannot get the scafContig sequences, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
#endif
				// ######################## Debug information ##########################

				// release the scafContig nodes
				releaseScafContigInScaf(localScafContigheadArr[0]);
				releaseScafContigInScaf(localScafContigheadArr[1]);

			}
		}
	}


	return SUCCESSFUL;
}

/**
 * Get scafContig end sequences.
 *  Note:
 *  	(1) The first scafContig end sequence is from the plus strand, and
 *  	(2) The second scafContig end sequence is from the minus strand.
 *
 *  	As shown below, the double lines in the two parts indicate the target end sequences of the two scafContigs.
 *
 * 				scafContig 1						scafContig 2
 *  		---------------======>				-------------------->
 *  		<---------------------				<======--------------
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getScafContigEndSeqs(contigOverlap *pContigOverlapInfo)
{
	int contigID1, contigID2, contigOrient1, contigOrient2, contigLen1, contigLen2;

	contigID1 = pContigOverlapInfo->contigID1;
	contigID2 = pContigOverlapInfo->contigID2;
	contigOrient1 = pContigOverlapInfo->orientation1;
	contigOrient2 = pContigOverlapInfo->orientation2;
	contigLen1 = contigInfoArr[contigID1-1].contigLen;
	contigLen2 = contigInfoArr[contigID2-1].contigLen;

	// get the lengths of the sequence at two contig ends
	if(contigLen1>=maxOverlapSeqLen)
	{
		scafContigEndSeqLenArr[0] = maxOverlapSeqLen;
	}else
	{
		scafContigEndSeqLenArr[0] = contigLen1;
	}

	if(contigOrient1==ORIENTATION_PLUS)
	{
		strcpy(scafContigEndSeqArr[0], contigInfoArr[contigID1-1].contigSeq+contigLen1-scafContigEndSeqLenArr[0]);
	}else
	{
		strncpy(scafContigEndSeqArr[0], contigInfoArr[contigID1-1].contigSeq, scafContigEndSeqLenArr[0]);
		scafContigEndSeqArr[0][ scafContigEndSeqLenArr[0] ] = '\0';

		if(reverseSeq(scafContigEndSeqArr[0], scafContigEndSeqLenArr[0])==FAILED) // reverse the sequence
		{
			printf("line=%d, In %s(), cannot reverse the sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	if(contigLen2>=maxOverlapSeqLen)
	{
		scafContigEndSeqLenArr[1] = maxOverlapSeqLen;
	}else
	{
		scafContigEndSeqLenArr[1] = contigLen2;
	}

	if(contigOrient2==ORIENTATION_PLUS)
	{
		strncpy(scafContigEndSeqArr[1], contigInfoArr[contigID2-1].contigSeq, scafContigEndSeqLenArr[1]);
		scafContigEndSeqArr[1][ scafContigEndSeqLenArr[1] ] = '\0';

		if(reverseSeq(scafContigEndSeqArr[1], scafContigEndSeqLenArr[1])==FAILED) // reverse the sequence
		{
			printf("line=%d, In %s(), cannot reverse the sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		strcpy(scafContigEndSeqArr[1], contigInfoArr[contigID2-1].contigSeq+contigLen2-scafContigEndSeqLenArr[1]);
	}

	return SUCCESSFUL;
}

/**
 * Get the comparison sequence in gap filling.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getComparisonSeqInScaf(int assemblyRound)
{
	int i;
	char reverseBase;
	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{ // the first round
		if(scafContigEndSeqLenArr[1] >= maxOverlapSeqLen)
		{
			comparisonSeqLenInScaf = maxOverlapSeqLen;
		}else
		{
			comparisonSeqLenInScaf = scafContigEndSeqLenArr[1];
		}

		for(i=0; i<comparisonSeqLenInScaf; i++)
		{
			switch(scafContigEndSeqArr[1][ scafContigEndSeqLenArr[1]-i-1 ])
			{
				case 'A': reverseBase = 'T'; break;
				case 'C': reverseBase = 'G'; break;
				case 'G': reverseBase = 'C'; break;
				case 'T': reverseBase = 'A'; break;
				default: printf("line=%d, In %s(), unknown base [ %c ], error!\n", __LINE__, __func__, scafContigEndSeqArr[1][ scafContigEndSeqLenArr[1]-i-1 ]); return FAILED;
			}
			comparisonSeqInScaf[i] = reverseBase;
		}
		comparisonSeqInScaf[comparisonSeqLenInScaf] = '\0';
	}else
	{ // the second round
		if(scafContigEndSeqLenArr[0] >= maxOverlapSeqLen)
		{
			comparisonSeqLenInScaf = maxOverlapSeqLen;
		}else
		{
			comparisonSeqLenInScaf = scafContigEndSeqLenArr[0];
		}

		//strncpy(comparisonSeqInScaf, scafContigEndSeqArr[0], comparisonSeqLenInScaf);
		strncpy(comparisonSeqInScaf, scafContigEndSeqArr[0]+scafContigEndSeqLenArr[0]-comparisonSeqLenInScaf, comparisonSeqLenInScaf);
		comparisonSeqInScaf[comparisonSeqLenInScaf] = '\0';

		if(reverseSeq(comparisonSeqInScaf, comparisonSeqLenInScaf)==FAILED) // reverse the sequence
		{
			printf("line=%d, In %s(), cannot reverse the sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Prepare the first scafContig nodes with length of READ_LEN before local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short prepareAssemblyInScaf(int contigID, int contigLen, int contigOrient, int assemblyRound)
{
	int64_t i, j;
	uint64_t hashcode, baseInt;

	// get the contig end sequence to sacfContigSeqLastReadLen
	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{ // the first round
		if(scafContigEndSeqLenArr[0] >= readLen)
			scafContigSeqLenLastReadLen = readLen;
		else
			scafContigSeqLenLastReadLen = scafContigEndSeqLenArr[0];

		strcpy(scafContigSeqLastReadLen, scafContigEndSeqArr[0]+scafContigEndSeqLenArr[0]-scafContigSeqLenLastReadLen);

		prepareAssemblyLenArr[0] = scafContigSeqLenLastReadLen;

	}else
	{ // the second round
		if(scafContigEndSeqLenArr[1] >= readLen)
			scafContigSeqLenLastReadLen = readLen;
		else
			scafContigSeqLenLastReadLen = scafContigEndSeqLenArr[1];

		strcpy(scafContigSeqLastReadLen, scafContigEndSeqArr[1]+scafContigEndSeqLenArr[1]-scafContigSeqLenLastReadLen);

		prepareAssemblyLenArr[1] = scafContigSeqLenLastReadLen;

	}

	// get first scafKmers and their integer sequences
	if(generateKmerSeqIntInScaf(kmerSeqIntAssembly, scafContigSeqLastReadLen)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate kmer integer sequence, error!\n", __LINE__, __func__);
	}
	hashcode = kmerhashInScaf( kmerSeqIntAssembly );
	scafKmers[0] = getKmerByHashInScaf(hashcode, kmerSeqIntAssembly, scafGrapDeBruijn);
	scafKmers[1] = getReverseKmerInScaf(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, scafGrapDeBruijn);

	// initialize the scafContig nodes
	if(initScafContig(scafContigSeqLastReadLen, scafContigSeqLenLastReadLen)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the scafContig nodes, error!\n", __LINE__, __func__);
		return FAILED;
	}

	scafContigIndex = kmerSize;
	scafContigLastReadLen = scafContighead;

	// initialize the scafAssemblingReads
	if(initScafAssemblingReads(contigID, contigLen, contigOrient, scafContigIndex, assemblyRound)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the scafAssembling reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=kmerSize; i<scafContigSeqLenLastReadLen; i++)
	{
		// get the integer form sequence of next scafKmers
		switch(scafContigSeqLastReadLen[i])
		{
			case 'A':
			case 'a': baseInt = 0; break;
			case 'C':
			case 'c': baseInt = 1; break;
			case 'G':
			case 'g': baseInt = 2; break;
			case 'T':
			case 't': baseInt = 3; break;
			default: printf("line=%d, In %s(), unknown base: %c, error!\n", __LINE__, __func__, scafContigSeqLastReadLen[i]); return FAILED;
		}

		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
				kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
			kmerSeqIntAssembly[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNumKmer-2));
		}
		kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | baseInt) & lastEntryMaskKmer;

		// get the scafKmers
		hashcode = kmerhashInScaf( kmerSeqIntAssembly );
		scafKmers[0] = getKmerByHashInScaf(hashcode, kmerSeqIntAssembly, scafGrapDeBruijn);
		scafKmers[1] = getReverseKmerInScaf(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, scafGrapDeBruijn);

		// update decision table
		if(updateScafAssemblingReadsInScaf(contigID, contigLen, contigOrient, scafContigIndex, assemblyRound)==FAILED)
		{
			printf("line=%d, In %s(), cannot update the scafAssembling reads, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// update the reads status in decision table
		if(updateAssemblingreadsStatusInScaf()==FAILED)
		{
			printf("line=%d, In %s(), cannot update the status of scafAssembling reads, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// remove the finished reads from decision table
		if(updateFinisedScafReadsInScaf()==FAILED)
		{
			printf("line=%d, In %s(), cannot remove the reads from decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}

		scafContigIndex ++;
	}

	return SUCCESSFUL;
}

/**
 * Initialize the scafContig nodes in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initScafContig(char *seq, int seq_len)
{
	int i, baseInt;
	scafContig *contig = NULL, *pre_contig = NULL;
	for(i=0; i<seq_len; i++)
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
			default: printf("line=%d, In %s(), unknown base %c, error!\n", __LINE__, __func__, seq[i]); return FAILED;
		}

		contig = (scafContig *)malloc(sizeof(scafContig));
		if(contig==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		contig->index = i + 1;
		contig->base = baseInt;
		contig->next = NULL;
		contig->ridposnum = 0;
		contig->pridposorientation = NULL;

		if(i==0)
		{
			scafContighead = contig;
			pre_contig = contig;
		}else
		{
			pre_contig->next = contig;
			pre_contig = contig;
		}
	}
	scafContigtail = pre_contig;

	return SUCCESSFUL;
}

/**
 * Initialize the scafAssembling reads in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initScafAssemblingReads(int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound)
{
	scafRidpos *ridpostable = NULL;
	unsigned int i, rpos = 0, posNum = 0;
	int returnCode, matedFlag;

	// add the reads having plus orientation
	if(scafKmers[0])
	{
		ridpostable = scafKmers[0]->ppos;
		posNum = scafKmers[0]->arraysize;
		for(i=0; i<posNum; i++)
		{
			rpos = ridpostable[i].pos;
			if(rpos==1 && ridpostable[i].used==0)
			{ // the first position and unused

				returnCode = validReadPairInScaf(ridpostable[i].rid, ORIENTATION_PLUS, contigID, contigLen, contigOrient, contigNodesNum, assemblyRound);
				if(returnCode==YES)
				{
					matedFlag = YES;
				}else if(returnCode==NO)
				{
					matedFlag = NO;
				}else
				{
					printf("line=%d, In %s(), cannot valid read pair, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(addAssemblingReadsInScaf(ridpostable[i].rid, rpos, ORIENTATION_PLUS, matedFlag)==FAILED)
				{
					printf("line=%d, In %s(), cannot add reads into decision table, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}

	// add the reads having minus orientation
	if(scafKmers[1])
	{
		ridpostable = scafKmers[1]->ppos;
		posNum = scafKmers[1]->arraysize;
		for(i=0; i<posNum; i++)
		{
			rpos = ridpostable[i].pos;
			if(rpos>readLen-kmerSize+1-errorRegLenEnd3 && ridpostable[i].used==0)
			//if(rpos==READ_LEN-KMER_SIZE+1 && ridpostable[i].used==0)
			{ // the last several positions and unused

				returnCode = validReadPairInScaf(ridpostable[i].rid, ORIENTATION_MINUS, contigID, contigLen, contigOrient, contigNodesNum, assemblyRound);
				if(returnCode==YES)
				{
					matedFlag = YES;
				}else if(returnCode==NO)
				{
					matedFlag = NO;
				}else
				{
					printf("line=%d, In %s(), cannot valid read pair, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(addAssemblingReadsInScaf(ridpostable[i].rid, rpos, ORIENTATION_MINUS, matedFlag)==FAILED)
				{
					printf("line=%d, In %s(), cannot add reads into decision table, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Add reads to scafAssemblingReads table in local assembly.
 */
short addAssemblingReadsInScaf(int rid, int rpos, int orientation, int matedFlag)
{
	scafAssemblingRead *this_assemblingRead = scafAssemblingReadArr + scafAssemblingReadsNum;

	this_assemblingRead->rid = rid;
	this_assemblingRead->firstpos = rpos;
	this_assemblingRead->orientation = orientation;
	this_assemblingRead->status = ASSEMBLING_STATUS;
	this_assemblingRead->kmerappeartimes = 1;
	this_assemblingRead->kmerunappeartimes = 0;
	this_assemblingRead->latestpos = rpos;
	this_assemblingRead->lastpos = rpos;
	this_assemblingRead->kmerunappearblocks = 0;
	this_assemblingRead->delsign = 0;
	this_assemblingRead->reserved = 0;
	this_assemblingRead->locked = 0;
	this_assemblingRead->matedFlag = matedFlag;

	scafAssemblingReadsNum ++;

	// check the allowed number of reads in the array
	if(scafAssemblingReadsNum==maxScafAssemblingReadsNumInArr)
	{
		if(reallocateDecisionTableInScaf()==FAILED)
		{
			printf("line=%d, In %s(), cannot reallocate memory for the reads in decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Reallocate the memory of scafAssembling reads in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short reallocateDecisionTableInScaf()
{
	scafAssemblingRead *tmp_pScafAssemblingReadArray;
	tmp_pScafAssemblingReadArray = (scafAssemblingRead*) malloc (2*maxScafAssemblingReadsNumInArr * sizeof(scafAssemblingRead));
	if(tmp_pScafAssemblingReadArray==NULL)
	{
		printf("line=%d, In %s(), cannot reallocate memory for the reads in decision table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// copy memory
	if(memcpy(tmp_pScafAssemblingReadArray, scafAssemblingReadArr, scafAssemblingReadsNum * sizeof(scafAssemblingRead))==NULL)
	{
		printf("line=%d, In %s(), cannot copy memory for the reads in decision table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	free(scafAssemblingReadArr);
	scafAssemblingReadArr = tmp_pScafAssemblingReadArray;
	maxScafAssemblingReadsNumInArr *= 2;

	return SUCCESSFUL;
}

/**
 * Get the maximal assembly length in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxMinAssemblyLen(int *maxAssemblyLen, int *minAssemblyLen, int assemblyRound, int *localScafContigNodesNumArr, int gapSize)
{
	int endLen1, endLen2, newSeqLen;
	int newEndLen1;

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{ // the first round
		endLen1 = prepareAssemblyLenArr[0];

		if(scafContigEndSeqLenArr[1]>=readLen)
			endLen2 = readLen;
		else
			endLen2 = scafContigEndSeqLenArr[1];

		*maxAssemblyLen = endLen1 + endLen2 + gapSize;
		*minAssemblyLen = (*maxAssemblyLen) - endLen2;

	}else
	{ // the second round
		endLen1 = prepareAssemblyLenArr[0];
		endLen2 = prepareAssemblyLenArr[1];

		newSeqLen = localScafContigNodesNumArr[0] - prepareAssemblyLenArr[0];
		*maxAssemblyLen = endLen1 + endLen2 + gapSize - newSeqLen;

		if(scafContigEndSeqLenArr[0]>=readLen)
			newEndLen1 = readLen;
		else
			newEndLen1 = scafContigEndSeqLenArr[0];
		*minAssemblyLen = (*maxAssemblyLen) - newEndLen1;
	}

	// adjust the values
	if(*maxAssemblyLen < 0)
		*maxAssemblyLen = 0;

	if(*minAssemblyLen < 0)
		*minAssemblyLen = 0;

	return SUCCESSFUL;
}

/**
 * Get next scafKmer in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNextKmerBySEInScaf(int contigNodesNum, int contigID, int contigLen, int contigOrient, int assemblyRound)
{
	uint64_t tmp_kmerseq[entriesPerKmer], tmp_kmerseqRev[entriesPerKmer];
	scafKmer *tmp_kmers[4][2] = {{0,0},{0,0},{0,0},{0,0}};  //tmp_kmers[i][0]为正向的kmer, tmp_kmers[i][1]为反向的kmer
	short validKmerNum = 0, base_index = 0; //有效的kmer数目
	int i = 0, j;
	double maxOcc = 0, secondOcc = 0;
	int maxOccIndex = -1, secondOccIndex = -1;
	kmer_len = 0;

	if(scafAssemblingReadsNum>MAX_DECISION_TABLE_SIZE_HTRES)
	{
		scafKmers[0] = scafKmers[1] = NULL;
		return SUCCESSFUL;
	}

	//将8个正反向kmer添加进临时数组tmp_kmers
	for(i=0; i<4; i++)
	{
		occsNumSE[i] = 0;

		//开始计算
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				tmp_kmerseq[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
			}
			tmp_kmerseq[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*(lastEntryBaseNumKmer-1)));
		}
		tmp_kmerseq[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | i) & lastEntryMaskKmer;

		tmp_kmers[i][0] = getKmerInScaf(tmp_kmerseq, scafGrapDeBruijn); //取得kmer的指针放入数组
		tmp_kmers[i][1] = getReverseKmerInScaf(tmp_kmerseqRev, tmp_kmerseq, scafGrapDeBruijn); //取得反向互补的kmer的指针放入数组
		if(tmp_kmers[i][0] || tmp_kmers[i][1])
		{  //两个kmer中只要有一个存在, 就统计其数量
			validKmerNum ++;
			base_index = i;
		}
	}

	//检测kmer的数量
	if(validKmerNum==1)
	{ //只有一个kmer, 则将该kmer返回
		//计算该kmer得分,
		//score[base_index] = computeKmerScore(tmp_kmers[base_index], occsNum+base_index, assemblingreads, numassemblingreads);
		if(computeKmerOccNumUnlockedBySEInScaf(tmp_kmers[base_index], occsNumSE+base_index, contigID, contigLen, contigOrient, contigNodesNum, assemblyRound)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute kmer occsNums by unlocked reads in local assembly, error!\n", __LINE__, __func__);
			return FAILED;
		}

		//========================= condition 1 ==============================
		//if(score[base_index]>0)
		if(occsNumSE[base_index]>=minKmerOccSE)
		{
			if(entriesPerKmer>=2)
			{
				for(j=0; j<entriesPerKmer-2; j++)
				{
					kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
				}
				kmerSeqIntAssembly[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNumKmer-2));
			}
			kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | base_index) & lastEntryMaskKmer;

			scafKmers[0] = tmp_kmers[base_index][0];
			scafKmers[1] = tmp_kmers[base_index][1];
		}else
		{
			scafKmers[0] = scafKmers[1] = NULL;
		}

		return SUCCESSFUL;
	}else if(validKmerNum==0)
	{ //没有可以扩展的kmer, 则返回NULL
		scafKmers[0] = scafKmers[1] = NULL;
		return SUCCESSFUL;
	}


	//****************************************************************************************************

	//kmer_len = LONG_KMER_SIZE;//kmer_len = LONG_KMER_SIZE;

	if(contigNodesNum>=longKmerSize)
		kmer_len = longKmerSize;
	else
		kmer_len = contigNodesNum;

	//while(kmer_len>kmerSize)
	while(kmer_len>=kmerSize)
	//while(kmer_len>=MIN_KMER_SIZE)
	{
		//开始计算每个kmer得分
		validKmerNum = 0;
		maxOcc = 0, secondOcc = 0;
		maxOccIndex = -1, secondOccIndex = -1;
		for(i=0; i<4; i++)
		{
			occsNumSE[i] = 0;

			//kmer连接数目限制
			//######################### begin #############################//
			if(tmp_kmers[i][0]==NULL && tmp_kmers[i][1]==NULL)
			{  //如果该kmer不存在, 并且反向互补的kmer也不存在, 则检测下一个kmer
				continue;
			}

			//========================= condition 2 ==============================
			else
			{
				//修剪kmer连接数小于阈值的kmers
				if(tmp_kmers[i][0]!=NULL && tmp_kmers[i][1]!=NULL)
				{
					if(tmp_kmers[i][0]->multiplicity+tmp_kmers[i][1]->multiplicity<minKmerOccSE)
					{
						continue;
					}
				}else if(tmp_kmers[i][0]!=NULL)
				{
					if(tmp_kmers[i][0]->multiplicity<minKmerOccSE)
					{
						continue;
					}
				}else if(tmp_kmers[i][1]!=NULL)
				{
					if(tmp_kmers[i][1]->multiplicity<minKmerOccSE)
					{
						continue;
					}
				}
			}

			//########################## end ########################//
			//if(contigNodesNum>READ_LEN)
			if(contigNodesNum>kmer_len)
			//if(contigNodesNum>=kmer_len)  //--bad result
			{
				if(computeLongKmerOccNumBySEInScaf(tmp_kmers[i], occsNumSE+i, kmer_len, contigID, contigLen, contigOrient, contigNodesNum, assemblyRound)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute occsNum of long k-mers, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}else
			{
				if(computeKmerOccNumBySEInScaf(tmp_kmers[i], occsNumSE+i, contigID, contigLen, contigOrient, contigNodesNum, assemblyRound)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute k-mer occsNums, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}

			//score[i] = computeKmerScore(tmp_kmers[i], occsNum+i, assemblingreads, numassemblingreads);

			//========================= condition 3 ==============================
			//if(score[i]>=MIN_SCORE_THRESHOLD/**MIN_SCORE_FACTOR*/)
			//if(score[i]>=MIN_SCORE_THRESHOLD && occsNum[i]>=MIN_CONNECT_KMER_NUM)
			//if(occsNum[i]>=MIN_CONNECT_KMER_NUM)
			if(occsNumSE[i]>0)
			{ //该kmer得分>0
				validKmerNum ++;  //有效的kmer数目+1

//				if(score[i]>tmp_max)
//				{
//					tmp_max = score[i];
//					tmp_maxIndex = i;
//				}
			}
//			else
//			{
//				score[i] = 0;
//				occsNum[i] = 0;
//			}
		}

		if(validKmerNum>0)
		{
			maxOcc = 0, maxOccIndex = -1, secondOcc = 0, secondOccIndex = -1;
			for(j=0; j<4; j++)
			{
				if(maxOcc<occsNumSE[j])
				{
					secondOcc = maxOcc;
					secondOccIndex = maxOccIndex;
					maxOcc = occsNumSE[j];
					maxOccIndex = j;
				}else if(secondOcc<occsNumSE[j])
				{
					secondOcc = occsNumSE[j];
					secondOccIndex = j;
				}
			}
		}

/*
		//========================= condition 4 ==============================
		//if(validKmerNum>1 && maxIndex1!=maxOccIndex)
		//if(validKmerNum>1 && occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM)
		if(validKmerNum>0 && occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM) //--best result
		//if(occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM || (validKmerNum>1 && maxIndex1!=maxOccIndex))
		//if(occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM || occsNum[secondOccIndex]>=MIN_CONNECT_KMER_NUM)
		{
			validKmerNum = 0;
		}
*/
/*
		//========================= condition 5 ==============================
		//if(validKmerNum>0 && (float)occsNum[maxIndex1]/numassemblingreads < VALID_OCC_RATIO)
		//if(validKmerNum>0 && occsNum[maxIndex1] > MAX_OCC_NUM)
		//if((float)occsNum[thiskmerseq & 3]/numassemblingreads < VALID_OCC_RATIO || occsNum[thiskmerseq & 3] > MAX_OCC_NUM)
		//if((float)occsNum[thiskmerseq & 3]/numassemblingreads < VALID_OCC_RATIO && occsNum[thiskmerseq & 3] > MAX_OCC_NUM)
		//if(validKmerNum>1 && (occsNum[secondIndex1] > 7 || occsNum[maxIndex1] > 80)) //--best result
		//if(validKmerNum>1 && (occsNum[secondIndex1] > MAX_SECOND_OCC_NUM || occsNum[maxIndex1] > MAX_FIRST_OCC_NUM))
		//if(validKmerNum>1 && kmer_len>=MIN_KMER_SIZE && (occsNum[secondIndex1] > 10 || occsNum[maxIndex1] > MAX_OCC_NUM))
		//if(validKmerNum>1 && (occsNum[secondIndex1] > MAX_SECOND_OCC_NUM || occsNum[maxIndex1] > MAX_FIRST_OCC_NUM
		//		|| ((float)occsNum[secondIndex1]/occsNum[maxIndex1]>0.5 && occsNum[secondOccIndex]>=MIN_CONNECT_KMER_NUM))) //-- good result
		if(validKmerNum>1 && (occsNum[secondIndex1] > MAX_SECOND_OCC_NUM || (occsNum[maxIndex1] > MAX_FIRST_OCC_NUM && occsNum[secondIndex1]>=2*MIN_CONNECT_KMER_NUM)
				|| ((float)occsNum[secondIndex1]/occsNum[maxIndex1]>0.5 && occsNum[secondIndex1]>=MIN_CONNECT_KMER_NUM))) //-- best result
		{
			validKmerNum = 0;
		}
*/
/*
		//========================= condition 6 ==============================
		if(validKmerNum>1 && ((float)occsNum[maxIndex1]/scafAssemblingReadsNum < VALID_OCC_RATIO && occsNum[secondIndex1] >= MIN_CONNECT_KMER_NUM))
		{
			validKmerNum = 0;
		}


		//========================= condition 7 ==============================
		if(validKmerNum>1 && second/max > SECOND_FIRST_SECORE_RATIO)
		{
			validKmerNum = 0;
		}
*/


		//========================= condition 8 ==============================
		//=====these several lines have bad result, thus they are omitted. =======//
		//if(validKmerNum>1 && (secondOcc/maxOcc>SECOND_FIRST_OCC_RATIO || (kmer_len==25 && secondOcc>SECOND_OCC_THRESHOLD)))
		//if(validKmerNum>1 && kmer_len==25 && ((secondOcc/maxOcc>SECOND_FIRST_OCC_RATIO && secondOcc>MIN_CONNECT_KMER_NUM) || (secondOcc>SECOND_OCC_THRESHOLD)))
		if(validKmerNum>1 && ((secondOcc/maxOcc>SECOND_FIRST_OCC_RATIO && secondOcc>minKmerOccSE) || (secondOcc>maxSecondOcc)))
		{
			validKmerNum = 0;
		}


		//========================= condition 9 ==============================
		if(kmer_len > longKmerSize - longKmerStepSize && maxOcc < minLongKmerOcc)
		{
			validKmerNum = 0;
		}

		if(validKmerNum==1)
		{
			//if(occsNum[tmp_maxIndex]<MIN_CONNECT_KMER_NUM)
//			if(occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM)
//			{
//				kmers[0] = kmers[1] = NULL;
//			}else
//			{
				//更新kmerseq，并返回得分最大的kmer
				if(entriesPerKmer>=2)
				{
					for(j=0; j<entriesPerKmer-2; j++)
					{
						kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
					}
					kmerSeqIntAssembly[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNumKmer-2));
				}
				kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | maxOccIndex) & lastEntryMaskKmer;

				scafKmers[0] = tmp_kmers[maxOccIndex][0];
				scafKmers[1] = tmp_kmers[maxOccIndex][1];
//			}

			return SUCCESSFUL;

		}
		//***********************************************************************************
		else if(validKmerNum==0)
		{
//			if(length_k<=KMER_SIZE+KMER_SIZE_STEP)
//				length_k -= KMER_SIZE_STEP / 2;
//			else
//				length_k -= KMER_SIZE_STEP;

			kmer_len -= longKmerStepSize;
			continue;
		}
		else
		{
			//更新kmerseq，并返回得分最大的kmer
			if(entriesPerKmer>=2)
			{
				for(j=0; j<entriesPerKmer-2; j++)
				{
					kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
				}
				kmerSeqIntAssembly[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNumKmer-2));
			}
			kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | maxOccIndex) & lastEntryMaskKmer;

			scafKmers[0] = tmp_kmers[maxOccIndex][0];
			scafKmers[1] = tmp_kmers[maxOccIndex][1];

			return SUCCESSFUL;
		}
		//***********************************************************************************
	}

	scafKmers[0] = scafKmers[1] = NULL;
	return SUCCESSFUL;
}

/**
 * Check whether the read has valid read pair in local assembly.
 *  @return:
 *  	(1) If the read has valid read pair, return YES; if invalid, return NO.
 *  	(2) If errors occurred, return ERROR.
 */
short validReadPairInScaf(uint64_t readID, int readOrient, int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound)
{
	uint64_t readID_paired;
	int readOrient_paired, contigID_paired, contigPos_paired;
	int fragmentSize, validReadOrient, validReadOrient_paired, conditionID;
	int64_t hitRow;
	ReadPos *pReadPosArray;
	int i, matchNum, tmp_contigEneSeqLen, tmp_ContigIndex, endSearchScafContigIndex, tmp_successReadNum;
	scafContig *tmp_contig;
	scafSuccessRead *tmp_successReads;

	// generate the paired readID
	if(readID%2==1)
	{ // odd --> even
		readID_paired = readID + 1;
	}else
	{ // enve --> odd
		readID_paired = readID - 1;
	}

	// get the paired read from mono read list (MRL)
	hitRow = getReadRowFromReadList(readID_paired, monoReadListArr, readItemNumInMRL);
	if(hitRow>=0)
	{ // find the paired read in mono read list (MRL)

		// get the contig sequence length and valid reads orientations
		conditionID = 0;
		if(assemblyRound==FIRST_ROUND_ASSEMBLY)
		{ // the first round
			tmp_contigEneSeqLen = scafContigEndSeqLenArr[0];

			// get valid reads orientations
			if(contigOrient==ORIENTATION_PLUS)
			{ // condition 1: the minus read, and its plus pair
				validReadOrient = ORIENTATION_MINUS;
				validReadOrient_paired = ORIENTATION_PLUS;
				conditionID = 1;
			}else
			{ // condition 2: the minus read, and its minus pair
				validReadOrient = ORIENTATION_MINUS;
				validReadOrient_paired = ORIENTATION_MINUS;
				conditionID = 2;
			}
		}else
		{ // the second round
			tmp_contigEneSeqLen = scafContigEndSeqLenArr[1];

			// get valid reads orientations
			if(contigOrient==ORIENTATION_PLUS)
			{ // condition 3: the minus read, and its minus pair
				validReadOrient = ORIENTATION_MINUS;
				validReadOrient_paired = ORIENTATION_MINUS;
				conditionID = 3;
			}else
			{ // condition 4: the minus read, and its plus pair
				validReadOrient = ORIENTATION_MINUS;
				validReadOrient_paired = ORIENTATION_PLUS;
				conditionID = 4;
			}
		}

		matchNum = monoReadListArr[hitRow].matchNum;
		pReadPosArray = monoReadPosArr + monoReadListArr[hitRow].firstRow;

		for(i=0; i<matchNum; i++)
		{
			contigID_paired = pReadPosArray[i].contigID;
			if(contigID_paired==contigID)
			{
				readOrient_paired = pReadPosArray[i].orientation;
				contigPos_paired = pReadPosArray[i].contigPos;

				//###################### Debug information ####################
				//printf("line=%d, In %s(), contigID=%d, readID_paired=%lu, readOrient_paired=%d, contigPos_paired=%d\n", __LINE__, __func__, contigID, readID_paired, readOrient_paired, contigPos_paired);
				//###################### Debug information ####################

				// check the two reads
				if(conditionID==1 || conditionID==4)
				{ // condition 1 or 4
					if(readOrient==validReadOrient && readOrient_paired==validReadOrient_paired)
					{
						fragmentSize = contigLen - contigPos_paired + 1 + contigNodesNum - kmerSize - tmp_contigEneSeqLen + readLen;
						if(fragmentSize<=meanSizeInsert+standardDeviationFactor*stardardDeviationInsert && fragmentSize>=meanSizeInsert-standardDeviationFactor*stardardDeviationInsert)
						{ // m - 3*sdev <= fragSize <= m + 3*sdev
							return YES;
						}
					}
				}else if(conditionID==2 || conditionID==3)
				{ // condition 2 or 3
					if(readOrient==validReadOrient && readOrient_paired==validReadOrient_paired)
					{
						fragmentSize = contigPos_paired + readLen - 1 + contigNodesNum - kmerSize - tmp_contigEneSeqLen + readLen;
						if(fragmentSize<=meanSizeInsert+standardDeviationFactor*stardardDeviationInsert && fragmentSize>=meanSizeInsert-standardDeviationFactor*stardardDeviationInsert)
						{ // m - 3*sdev <= fragSize <= m + 3*sdev
							return YES;
						}
					}
				}else
				{ // erroneous conditions
					printf("line=%d, In %s(), unknown conditionID=%d, error!\n", __LINE__, __func__, conditionID);
					return ERROR;
				}
			}
		}
	}else
	{ // cannot find the paired read in mono read list (MRL)
		if(readOrient==ORIENTATION_MINUS)
		{ // only consider the minus reads
			endSearchScafContigIndex = contigNodesNum - meanSizeInsert - standardDeviationFactor*stardardDeviationInsert + 1 - readLen + 1;
			if(endSearchScafContigIndex>0)
			{
				tmp_ContigIndex = 1;
				tmp_contig = scafContighead;
				while(tmp_contig)
				{
					if(tmp_ContigIndex>endSearchScafContigIndex)
						break;

					tmp_successReadNum = tmp_contig->ridposnum;
					tmp_successReads = tmp_contig->pridposorientation;
					for(i=0; i<tmp_successReadNum; i++)
					{ // check each successful reads
						if(tmp_successReads[i].rid==readID_paired)
						{
							readOrient_paired = tmp_successReads[i].orientation;
							if(readOrient_paired==ORIENTATION_PLUS)
							{ // the paired read has plus orientation

								//###################### Debug information ####################
								//printf("line=%d, In %s(), contigID=%d, readID_paired=%lu, readOrient_paired=%d, tmp_ContigIndex=%d\n", __LINE__, __func__, contigID, readID_paired, readOrient_paired, tmp_ContigIndex);
								//###################### Debug information ####################

								fragmentSize = contigNodesNum - kmerSize - tmp_ContigIndex + 1 + readLen;
								if(fragmentSize<=meanSizeInsert+standardDeviationFactor*stardardDeviationInsert && fragmentSize>=meanSizeInsert-standardDeviationFactor*stardardDeviationInsert)
								{ // m - 3*sdev <= fragSize <= m + 3*sdev
									return YES;
								}
							}
						}
					}

					tmp_ContigIndex ++;
					tmp_contig = tmp_contig->next;
				}
			}
		}
	}

	return NO;
}

/**
 * Compute the sores according to long K-mers.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeLongKmerOccNumBySEInScaf(scafKmer *tmp_kmers[2], int *occNum, int length_k, int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound)
{
	scafAssemblingRead *this_assemblingRead = NULL;
	scafRidpos *rid_pos = NULL;
	scafRidpos *rid_pos_table = NULL;
	int posNum = 0;
	unsigned int this_rid = 0, this_pos = 0;
	int i = 0, j = 0, properIndex = -1, startIndex = -1, exectIndex = -1, limitLastpos = -1, exitFlag = NO;

	*occNum = 0;

	for(i=0; i<scafAssemblingReadsNum; i++)
	{ //顺序访问决策表中的每行
		//找合适的lastpos的read
		if((scafAssemblingReadArr[i].reserved==0 && scafAssemblingReadArr[i].lastpos>0)
				&& (scafAssemblingReadArr[i].kmerunappeartimes==0 && scafAssemblingReadArr[i].kmerunappearblocks==0)
				&& ((scafAssemblingReadArr[i].orientation==ORIENTATION_PLUS && scafAssemblingReadArr[i].lastpos>length_k-kmerSize)
						||(scafAssemblingReadArr[i].orientation==ORIENTATION_MINUS && scafAssemblingReadArr[i].lastpos<readLen-kmerSize-(length_k-kmerSize)+1)))
		{
			properIndex = getProperIndexInScaf(scafAssemblingReadArr+i);

			//********************* 调试信息 ******************
			if(properIndex<0)
			{
				printf("line=%d, In %s(), properIndex=%d, i=%d, Error!\n", __LINE__, __func__, properIndex, i);
				return FAILED;
			}
			//********************* 调试信息 ******************

			this_assemblingRead = scafAssemblingReadArr + properIndex;
			limitLastpos = this_assemblingRead->lastpos;

			//********************* 调试信息 *******************
			//if(properIndex!=i)
			//{
			//	printf("line=%d, In %s(), properIndex=%d, i=%d\n", __LINE__, __func__, properIndex, i);
			//}
			//********************* 调试信息 *******************

			exitFlag = NO;
			while(exitFlag==NO)
			{
				for(j=0; j<2; j++) //j==0为正向kmer, j==1为反向互补kmer
				{
					if(tmp_kmers[j] && (this_assemblingRead->kmerunappeartimes==0 && this_assemblingRead->kmerunappearblocks==0)
							&& ((j==0 && this_assemblingRead->orientation==ORIENTATION_PLUS && this_assemblingRead->lastpos>length_k-kmerSize)
									||(j==1 && this_assemblingRead->orientation==ORIENTATION_MINUS && this_assemblingRead->lastpos<readLen-kmerSize-(length_k-kmerSize)+1)))
					{ //计算正向和反向互补kmer的得分
						//判断ridpostable表中该read是否存在
						rid_pos_table = tmp_kmers[j]->ppos;
						posNum = tmp_kmers[j]->arraysize;
						startIndex = findStartScafRidposIndexInScaf(this_assemblingRead->rid, rid_pos_table, posNum);//二分查找rid_pos_table
						if(startIndex>=0)
						{  //存在, 继续查找精确位置
							if(this_assemblingRead->lastpos>0)
							{ //该read上次拼接出现
								if(j==0) //正向kmer
									exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->lastpos+1, startIndex, rid_pos_table, posNum);
								else //反向互补kmer
									exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->lastpos-1, startIndex, rid_pos_table, posNum);
							}else
							{  //该read上次拼接未出现
								if(j==0) //正向kmer
									exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->firstpos+this_assemblingRead->kmerappeartimes+this_assemblingRead->kmerunappeartimes, startIndex, rid_pos_table, posNum);
								else //反向互补kmer
									exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->firstpos-this_assemblingRead->kmerappeartimes-this_assemblingRead->kmerunappeartimes, startIndex, rid_pos_table, posNum);
							}

							if(exectIndex>=0)
							{  //该read位置合理
								rid_pos = rid_pos_table + exectIndex;
								this_rid = rid_pos->rid;  //取得read的rid
								this_pos = rid_pos->pos;  //取得pos
								if((this_assemblingRead->orientation==ORIENTATION_PLUS && this_pos>length_k-kmerSize)
									|| (this_assemblingRead->orientation==ORIENTATION_MINUS && this_pos<readLen-kmerSize-(length_k-kmerSize)+1))
								{
									if(rid_pos->used==0) //该read未被删除
									{
										if(this_assemblingRead->lastpos==0)
										{ //该read上次拼接未出现
											if(this_assemblingRead->kmerunappeartimes==kmerSize)
											{
												(*occNum) ++;
											}else  //连续未出现的次数不等于kmer_size，说明该read拼接不合理
											{
												continue;
											}
										}else if(this_assemblingRead->kmerunappeartimes==0)
										{
											(*occNum) ++;
										}
									}
								}
							}
							exitFlag = YES;
						}else
						{
							exitFlag = YES;
						}
					}
				} //end for(j)

				if(exitFlag==NO)
				{
					//该read位置不合理, 寻找合理的位置
					properIndex = getProperIndexLimitedInScaf(scafAssemblingReadArr+i, limitLastpos);
					if(properIndex<0)
					{ //没有了合适的reads, 退出while循环
						exitFlag = YES; //退出标记置为YES
						break;
					}
					this_assemblingRead = scafAssemblingReadArr+properIndex;
					limitLastpos = this_assemblingRead->lastpos;
				}

			} //end while(exitFlag)
		} //end if(reserved)
	}// end for(i)

#if 1
	int returnCode;
	if(tmp_kmers[0])
	{ //正向kmer不为空
		rid_pos_table = tmp_kmers[0]->ppos;
		posNum = tmp_kmers[0]->arraysize;
		//计算ridpostable表中未考虑的reads的得分, 也即是新的reads的得分: 只要可以拼接，则按照1分计算
		for(i=0; i<posNum; i++)//遍历ridpostable表
		{
			if(rid_pos_table[i].used==0 && rid_pos_table[i].pos==1)
			{  //未被删除, 未标记并且位置合理,则将该read可以拼接, 并计算该read的得分

				returnCode = validReadPairInScaf(rid_pos_table[i].rid, ORIENTATION_PLUS, contigID, contigLen, contigOrient, contigNodesNum, assemblyRound);
				if(returnCode==YES)
				{
					(*occNum) ++;
				}else if(returnCode==ERROR)
				{
					printf("line=%d, In %s(), cannot valid read pair, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
			//rid_pos_table[i].reserved = 0; //ridpostable表中保留位清零
		}
	}

	if(tmp_kmers[1])
	{ //反向互补的kmer不为空
		rid_pos_table = tmp_kmers[1]->ppos;
		posNum = tmp_kmers[1]->arraysize;
		//计算ridpostable表中未考虑的reads的得分, 也即是新的reads的得分: 只要可以拼接，则按照1分计算
		for(i=0; i<posNum; i++)//遍历ridpostable表
		{
			if(rid_pos_table[i].used==0 && rid_pos_table[i].pos>readLen-kmerSize+1-errorRegLenEnd3)
			{  //未被删除, 未标记并且位置合理,则将该read可以拼接, 并计算该read的得分
				returnCode = validReadPairInScaf(rid_pos_table[i].rid, ORIENTATION_MINUS, contigID, contigLen, contigOrient, contigNodesNum, assemblyRound);
				if(returnCode==YES)
				{
					(*occNum) ++;
				}else if(returnCode==ERROR)
				{
					printf("line=%d, In %s(), cannot valid read pair, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
			//rid_pos_table[i].reserved = 0; //ridpostable表中保留位清零
		}
	}

#endif

	for(i=0; i<scafAssemblingReadsNum; i++) scafAssemblingReadArr[i].reserved = 0; //决策表中保留位清零

	return SUCCESSFUL;
}

/**
 *   计算kmer得分.
 *   如果决策表中含有锁定的reads, 则只考虑锁定的reads.
 *   否则,考虑全部的reads.
 *
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
*/
short computeKmerOccNumBySEInScaf(scafKmer *tmp_kmers[2], int *occNum, int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound)
{
	if(lockedReadsNumInScaf>=lockedReadsNumThres)
	{ //如果决策表中含有锁定的reads, 则只考虑锁定的reads.
		if(computeKmerOccNumLockedBySEInScaf(tmp_kmers, occNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute kmer occNum in local assembly, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{ //否则,考虑全部的reads.
		//计算决策表中的reads的得分
		if(computeKmerOccNumUnlockedBySEInScaf(tmp_kmers, occNum, contigID, contigLen, contigOrient, contigNodesNum, assemblyRound)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute kmer occNum in local assembly, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * 按照一定的策略计算kmer得分.
 * 考虑全部的reads.
 *
 *   当前只考虑:
 *   	(1) 上次拼接出现的reads;
 *   	(2) 当次拼接出现, 上次拼接未出现并连续12次未出现的情况.
 *   	(3) ridpostable表中未考虑的reads的得分, 也即是新的reads的得分.
 *
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeKmerOccNumUnlockedBySEInScaf(scafKmer *tmp_kmers[2], int *occNum, int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound)
{
	scafAssemblingRead *this_assemblingRead = NULL;
	scafRidpos *rid_pos = NULL;
	scafRidpos *rid_pos_table = NULL;
	int posNum = 0;
	unsigned int this_rid = 0, this_pos = 0;
	int i = 0, j = 0, properIndex = -1, startIndex = -1, exectIndex = -1, limitLastpos = -1, exitFlag = NO;

	*occNum = 0;

	for(i=0; i<scafAssemblingReadsNum; i++)
	{ //顺序访问决策表中的每行
		//找合适的lastpos的read
		if(scafAssemblingReadArr[i].reserved==0 && scafAssemblingReadArr[i].lastpos>0)
		{
			properIndex = getProperIndexInScaf(scafAssemblingReadArr + i);

			//********************* 调试信息 ******************
			if(properIndex<0)
			{
				printf("line=%d, In %s(), properIndex=%d, i=%d, Error!\n", __LINE__, __func__, properIndex, i);
				return FAILED;
			}
			//********************* 调试信息 ******************

			this_assemblingRead = scafAssemblingReadArr + properIndex;
			limitLastpos = this_assemblingRead->lastpos;

			//********************* 调试信息 *******************
			//if(properIndex!=i)
			//{
			//	printf("In computeKmerScoreUnlocked(), properIndex=%d, i=%d\n", properIndex, i);
			//}
			//********************* 调试信息 *******************

			exitFlag = NO;
			while(exitFlag==NO)
			{
				for(j=0; j<2; j++) //j==0为正向kmer, j==1为反向互补kmer
				{
					if(tmp_kmers[j] && ((j==0 && this_assemblingRead->orientation==ORIENTATION_PLUS)||(j==1 && this_assemblingRead->orientation==ORIENTATION_MINUS)))
					{ //计算正向和反向互补kmer的得分
						//判断ridpostable表中该read是否存在
						rid_pos_table = tmp_kmers[j]->ppos;
						posNum = tmp_kmers[j]->arraysize;
						startIndex = findStartScafRidposIndexInScaf(this_assemblingRead->rid, rid_pos_table, posNum);//二分查找rid_pos_table
						if(startIndex>=0)
						{  //存在, 继续查找精确位置
							if(this_assemblingRead->lastpos>0)
							{ //该read上次拼接出现
								if(j==0) //正向kmer
									exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->lastpos+1, startIndex, rid_pos_table, posNum);
								else //反向互补kmer
									exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->lastpos-1, startIndex, rid_pos_table, posNum);
							}else
							{  //该read上次拼接未出现
								if(j==0) //正向kmer
									exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->firstpos+this_assemblingRead->kmerappeartimes+this_assemblingRead->kmerunappeartimes, startIndex, rid_pos_table, posNum);
								else //反向互补kmer
									exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->firstpos-this_assemblingRead->kmerappeartimes-this_assemblingRead->kmerunappeartimes, startIndex, rid_pos_table, posNum);
							}

							if(exectIndex>=0)
							{  //该read位置合理
								rid_pos = rid_pos_table + exectIndex;
								if(rid_pos->used==0) //该read未被删除
								{
									this_rid = rid_pos->rid;  //取得read的rid
									this_pos = rid_pos->pos;  //取得pos
									if(this_assemblingRead->lastpos==0)
									{ //该read上次拼接未出现
										if(this_assemblingRead->kmerunappeartimes==kmerSize)
										{
											(*occNum) ++;
										}else  //连续未出现的次数不等于kmer_size，说明该read拼接不合理
										{
											continue;
										}
									}else
									{
										(*occNum) ++;
									}
								}
							}
							exitFlag = YES;
						}else
						{
							exitFlag = YES;
						}
					}
				} //end for(j)


				if(exitFlag==NO)
				{
					//该read位置不合理, 寻找合理的位置
					properIndex = getProperIndexLimitedInScaf(scafAssemblingReadArr+i, limitLastpos);
					if(properIndex<0)
					{ //没有了合适的reads, 退出while循环
						exitFlag = YES; //退出标记置为YES
						break;
					}
					this_assemblingRead = scafAssemblingReadArr + properIndex;
					limitLastpos = this_assemblingRead->lastpos;
				}

			} //end while(exitFlag)
		} //end if(reserved)
	}// end for(i)

#if 1
	if(tmp_kmers[0])
	{ //正向kmer不为空
		rid_pos_table = tmp_kmers[0]->ppos;
		posNum = tmp_kmers[0]->arraysize;
		//计算ridpostable表中未考虑的reads的得分, 也即是新的reads的得分: 只要可以拼接，则按照1分计算
		for(i=0; i<posNum; i++)//遍历ridpostable表
		{
			if(rid_pos_table[i].used==0 && rid_pos_table[i].pos==1)
			{  //未被删除, 未标记并且位置合理,则将该read可以拼接, 并计算该read的得分
				(*occNum) ++;
			}
			//rid_pos_table[i].reserved = 0; //ridpostable表中保留位清零
		}
	}

	if(tmp_kmers[1])
	{ //反向互补的kmer不为空
		rid_pos_table = tmp_kmers[1]->ppos;
		posNum = tmp_kmers[1]->arraysize;
		//计算ridpostable表中未考虑的reads的得分, 也即是新的reads的得分: 只要可以拼接，则按照1分计算
		for(i=0; i<posNum; i++)//遍历ridpostable表
		{
			if(rid_pos_table[i].used==0 && rid_pos_table[i].pos>readLen-kmerSize+1-errorRegLenEnd3)
			{  //未被删除, 未标记并且位置合理,则将该read可以拼接, 并计算该read的得分
				(*occNum) ++;
			}
			//rid_pos_table[i].reserved = 0; //ridpostable表中保留位清零
		}
	}
#endif

	for(i=0;i<scafAssemblingReadsNum;i++) scafAssemblingReadArr[i].reserved = 0; //决策表中保留位清零

	return SUCCESSFUL;
}

/**
 * 按照一定的策略计算kmer得分.
 * 决策表中含有锁定的reads, 则只考虑锁定的reads.
 *
 *   当前只考虑:
 *   	(1) 上次拼接出现的reads;
 *   	(2) 当次拼接出现, 上次拼接未出现并连续12次未出现的情况.
 *   	(3) ridpostable表中未考虑的reads的得分, 也即是新的reads的得分.
 *
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */

short computeKmerOccNumLockedBySEInScaf(scafKmer *tmp_kmers[2], int *occNum)
{
	scafAssemblingRead *this_assemblingRead = NULL;
	scafRidpos *rid_pos = NULL;
	scafRidpos *rid_pos_table = NULL;
	int posNum = 0;
	unsigned int this_rid = 0;
	unsigned short this_pos = 0;
	int i, j, properIndex, startIndex, exectIndex, limitLastpos = -1, exitFlag = NO;

	*occNum = 0;

	for(i=0; i<scafAssemblingReadsNum; i++)
	{ //顺序访问决策表中的每行
		//找合适的lastpos的read, 并计算得分
		if(scafAssemblingReadArr[i].locked)
		{
			//找合适的lastpos的read
			if(scafAssemblingReadArr[i].reserved==0 && scafAssemblingReadArr[i].lastpos>0)
			{
				properIndex = getProperIndexInScaf(scafAssemblingReadArr+i);

				//********************* 调试信息 ******************
				if(properIndex<0)
				{
					printf("line=%d, In %s(), properIndex=%d, i=%d, Error!\n", __LINE__, __func__, properIndex, i);
					return FAILED;
				}
				//********************* 调试信息 ******************

				this_assemblingRead = scafAssemblingReadArr + properIndex;
				limitLastpos = this_assemblingRead->lastpos;

				//********************* 调试信息 *******************
				//if(properIndex!=i)
				//{
				//	printf("line=%d, In %s(), properIndex=%d, i=%d\n", __LINE__, __func__, properIndex, i);
				//}
				//********************* 调试信息 *******************

				exitFlag = NO;
				while(exitFlag==NO)
				{
					for(j=0; j<2; j++)
					{ //j==0为正向kmer, j==1为反向互补kmer
						if(tmp_kmers[j] && ((j==0 && this_assemblingRead->orientation==ORIENTATION_PLUS) || (j==1 && this_assemblingRead->orientation==ORIENTATION_MINUS)))
						{ //计算正向和反向互补kmer的得分
							//判断ridpostable表中该read是否存在
							rid_pos_table = tmp_kmers[j]->ppos;
							posNum = tmp_kmers[j]->arraysize;
							startIndex = findStartScafRidposIndexInScaf(this_assemblingRead->rid, rid_pos_table, posNum);//二分查找rid_pos_table
							if(startIndex>=0)
							{ //存在, 继续查找精确位置
								if(this_assemblingRead->lastpos>0)
								{ //该read上次拼接出现
									if(j==0) //正向kmer
										exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->lastpos+1, startIndex, rid_pos_table, posNum);
									else //反向互补kmer
										exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->lastpos-1, startIndex, rid_pos_table, posNum);
								}else
								{  //该read上次拼接未出现
									if(j==0) //正向kmer
										exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->firstpos+this_assemblingRead->kmerappeartimes+this_assemblingRead->kmerunappeartimes, startIndex, rid_pos_table, posNum);
									else //反向互补kmer
										exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->firstpos-this_assemblingRead->kmerappeartimes-this_assemblingRead->kmerunappeartimes, startIndex, rid_pos_table, posNum);
								}

								if(exectIndex>=0)
								{ //该read位置合理
									rid_pos = rid_pos_table + exectIndex;
									if(rid_pos->used==0) //该read未被删除
									{
										this_rid = rid_pos->rid;  //取得read的rid
										this_pos = rid_pos->pos;  //取得pos
										if(this_assemblingRead->lastpos==0)
										{ //该read上次拼接未出现
											if(this_assemblingRead->kmerunappeartimes==kmerSize)
											{
												(*occNum) ++;
											}else  //连续未出现的次数不等于kmer_size，说明该read拼接不合理
											{
												continue;
											}
										}else
										{
											(*occNum) ++;
										}
									}
									exitFlag = YES; //退出标记置为YES
								}
							}else
							{
								exitFlag = YES;
							}
						}

					} //end for(j)

					if(exitFlag==NO)
					{
						//该read位置不合理, 寻找合理的位置
						properIndex = getProperIndexLimitedInScaf(scafAssemblingReadArr+i, limitLastpos);
						if(properIndex<0)
						{ //没有了合适的reads, 退出while循环
							exitFlag = YES; //退出标记置为YES
							break;
						}
						this_assemblingRead = scafAssemblingReadArr + properIndex;
						limitLastpos = this_assemblingRead->lastpos;
					}
				} //end while(exitFlag)

			} //end if(reserved)
		} //end if (locked)
	}// end for(i)

#if 1
	if(tmp_kmers[0])
	{ //正向kmer不为空
		rid_pos_table = tmp_kmers[0]->ppos;
		posNum = tmp_kmers[0]->arraysize;
		//计算ridpostable表中未考虑的reads的得分, 也即是新的reads的得分: 只要可以拼接，则按照1分计算
		for(i=0; i<posNum; i++)//遍历ridpostable表
		{
			if(rid_pos_table[i].used==0 && rid_pos_table[i].pos==1)
			{  //未被删除, 未标记并且位置合理,则将该read可以拼接, 并计算该read的得分
				(*occNum) ++;
			}
			//rid_pos_table[i].reserved = 0; //ridpostable表中保留位清零
		}
	}

	if(tmp_kmers[1])
	{ //反向互补的kmer不为空
		rid_pos_table = tmp_kmers[1]->ppos;
		posNum = tmp_kmers[1]->arraysize;
		//计算ridpostable表中未考虑的reads的得分, 也即是新的reads的得分: 只要可以拼接，则按照1分计算
		for(i=0; i<posNum; i++)//遍历ridpostable表
		{
			if(rid_pos_table[i].used==0 && rid_pos_table[i].pos>readLen-kmerSize+1-errorRegLenEnd3)
			{  //未被删除, 未标记并且位置合理,则将该read可以拼接, 并计算该read的得分
				(*occNum) ++;
			}
			//rid_pos_table[i].reserved = 0; //ridpostable表中保留位清零
		}
	}
#endif

	for(i=0; i<scafAssemblingReadsNum; i++) scafAssemblingReadArr[i].reserved = 0; //决策表中保留位清零

	return SUCCESSFUL;
}


/**
 * 从决策表中相同rid的read中找lastpos合适的read，返回其下标，并标记rid的这些read.
 *	正向reads取lastpos最大的, 反向reads取lastpos最小的.
 *   @ return:
 *     成功, 返回lastpos合适的read在决策表中的下标; 失败, 返回-1.
 */
int getProperIndexInScaf(scafAssemblingRead *assemblingread)
{
	int proLastpos = 0, i, proIndex = -1, rid = assemblingread->rid;
	if(assemblingread->orientation==ORIENTATION_PLUS)
	{ //read正向
		for(i=0; i<scafAssemblingReadsNum; i++)
		{
			if(scafAssemblingReadArr[i].rid==rid && scafAssemblingReadArr[i].orientation==ORIENTATION_PLUS)
			{
				scafAssemblingReadArr[i].reserved = 1;
				if(scafAssemblingReadArr[i].lastpos>proLastpos)
				{
					proLastpos = scafAssemblingReadArr[i].lastpos;
					proIndex = i;
				}
			}
		}
	}else
	{ //read反向
		proLastpos = INT_MAX;
		for(i=0; i<scafAssemblingReadsNum; i++)
		{
			if(scafAssemblingReadArr[i].rid==rid && scafAssemblingReadArr[i].orientation==ORIENTATION_MINUS)
			{
				scafAssemblingReadArr[i].reserved = 1;
				if(scafAssemblingReadArr[i].lastpos>0 && scafAssemblingReadArr[i].lastpos<proLastpos)
				{
					proLastpos = scafAssemblingReadArr[i].lastpos;
					proIndex = i;
				}
			}
		}
	}
	return proIndex;
}

/**
 * 取得具有限制的合适的位置.
 * 正向reads取lastpos最大的, 反向reads取lastpos最小的.
 *   @ return:
 *     成功, 返回lastpos合适的read的下标; 失败, 返回-1.
 */
int getProperIndexLimitedInScaf(scafAssemblingRead *assemblingread, int limitLastpos)
{
	int pro = limitLastpos, i, proIndex = -1, rid = assemblingread->rid;
	if(assemblingread->orientation==ORIENTATION_PLUS)
	{ //read正向
		for(i=0;i<scafAssemblingReadsNum;i++)
		{
			if(scafAssemblingReadArr[i].rid==rid)
			{
				//scafAssemblingReadArr[i].reserved = 1;
				if(scafAssemblingReadArr[i].lastpos>pro)
				{
					pro = scafAssemblingReadArr[i].lastpos;
					proIndex = i;
				}
			}
		}
	}else
	{ //read反向
		for(i=0;i<scafAssemblingReadsNum;i++)
		{
			if(scafAssemblingReadArr[i].rid==rid)
			{
				//scafAssemblingReadArr[i].reserved = 1;
				if(scafAssemblingReadArr[i].lastpos>0 && scafAssemblingReadArr[i].lastpos<pro)
				{
					pro = scafAssemblingReadArr[i].lastpos;
					proIndex = i;
				}
			}
		}
	}
	return proIndex;
}

/**
 *   将碱基添加进contig中, 并返回末端contig节点.
 *     @ return:
 *       成功, 返回contig末端节点; 失败, 返回NULL.
*/
short appendContigBaseInScaf(unsigned char base, int contigIndex)
{
	scafContig *contignode = (scafContig*) malloc(sizeof(scafContig));
	if(contignode==NULL)
	{  //分配内存失败，返回NULL
		printf("line=%d, In %s(), can not malloc the memory. Error and Exit!\n", __LINE__, __func__);
		return FAILED;
	}

	contignode->index = contigIndex;
	contignode->base = base;
	contignode->ridposnum = 0;
	contignode->pridposorientation = NULL;
	contignode->next = NULL;
	contignode->reserved = 0;

	scafContigtail->next = contignode;
	scafContigtail = contignode;

	return SUCCESSFUL;
}

/**
 * 更新锁定的reads数量: lockedReadsNum.
 * 	@ return:
 * 		锁定的reads数量: lockedReadsNum.
 */
void updateLockedReadsInScaf()
{
	scafAssemblingRead *this_assemblingRead;
	int i = 0, appearNum = 0;

	if(lockedReadsNumInScaf>0)
	{
		//更新锁定的reads数量
		for(i=0; i<scafAssemblingReadsNum; i++)
		{
			if(scafAssemblingReadArr[i].locked==1)
			{ //针对锁定的reads
				if(scafAssemblingReadArr[i].status!=ASSEMBLING_STATUS)
				{ //状态不为拼接状态, 则锁定的数目-1
					lockedReadsNumInScaf --;
				}else if(scafAssemblingReadArr[i].lastpos>0)
				{ //本次拼接出现, 出现次数+1
					appearNum++;
				}
			}
		}
	}

	//更新锁定标记
	if(appearNum<minKmerOccSE*2 || lockedReadsNumInScaf<lockedReadsNumThres)
	{ //出现的reads数量为0, 或者锁定的reads数量小于阈值, 则重新更新reads的锁定标记

//		if(appearNum>=MIN_CONNECT_KMER_NUM*2)
//		{
//			printf("appearNum=%d, lockedReadsNumInScaf=%d\n", appearNum, lockedReadsNumInScaf);
//		}else
//		{
//			printf("++++++++++appearNum=%d, lockedReadsNumInScaf=%d\n", appearNum, lockedReadsNumInScaf);
//		}

		lockedReadsNumInScaf = 0; //重新设定初始锁定数量为0

		//更新决策表中的reads的锁定标记
		this_assemblingRead = scafAssemblingReadArr;
		for(i = 0; i<scafAssemblingReadsNum; i++)
		{ //遍历决策表, 寻找可以锁定的reads

			//******************* 调试信息 ******************
			//if(this_assemblingRead->rid==1655993)
			//{
			//	printf("line=%d, In %s(), rid=%d\n", __LINE__, __func__, this_assemblingRead->rid);
			//}
			//******************* 调试信息 ******************

			if(this_assemblingRead->lastpos==0 && this_assemblingRead->locked==1)
			{ //当次拼接未出现, 并且处于锁定状态, 则将该read解锁
				this_assemblingRead->locked = 0;
			}else if(this_assemblingRead->lastpos>0	&& this_assemblingRead->status==ASSEMBLING_STATUS)
			{ //该read当次拼接出现, 并且处于拼接中, 则需要锁定该read, 并更新锁定数量
				this_assemblingRead->locked = 1; //锁定该read
				lockedReadsNumInScaf ++; //锁定的reads数量+1
			}
			this_assemblingRead ++;
		}
	}
}

/**
 * 将成功的reads从de Bruijn图中删除。成功，返回SUCCESSFUL；失败，返回FAILED。
 */
short delReadsFromGraphInScaf()
{
	int i, j = 0, rpos = 0, len = scafContigSeqLenLastReadLen;
	char *reversed_lastseq;
	short Num, perfectMatchFlag;  //kmer数目, read完全与contig匹配的标记
	uint64_t hashcode = 0, rid = 0;
	uint64_t tmp_kmerseq[entriesPerKmer];
	int baseInt, first_kmer_index = 0;//记录read在lastseq36中第一个kmer的下标

	reversed_lastseq = (char *) calloc(readLen+1, sizeof(char));
	if(reversed_lastseq==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//将lastseq36复制到临时数组
	strcpy(reversed_lastseq, scafContigSeqLastReadLen);
	reversed_lastseq[len] = '\0';

	//将临时碱基反向互补
	if(reverseSeq(reversed_lastseq, len)==FAILED)
	{
		printf("line=%d, In %s(), cannot reverse the read sequence: %s", __LINE__, __func__, reversed_lastseq);
		return FAILED;
	}

	//取得正向的序列的第一个kmer的哈希值
	//strncpy(tmpseq, lastseq36, KMER_SIZE); //取得该kmer的碱基序列
	//tmpseq[KMER_SIZE] = '\0';
	//hashcode = kmerhash(tmpseq);
	//first_hashcode = kmerhash(lastseq36);

	//取得反向互补的序列的第一个kmer的哈希值
	//strncpy(tmpseq, reversed_lastseq36, KMER_SIZE); //取得该kmer的碱基序列
	//tmpseq[KMER_SIZE] = '\0';

	//reversed_first_hashcode = kmerhashInScaf(reversed_lastseq);

	i = 0;
	Num = len - kmerSize + 1;
	while(i<scafSuccessReadsNum)
	{
		perfectMatchFlag = YES;
		rid = scafSuccessReadArr[i].rid; //取得read ID

		/******************** 调试信息 ******************/
		//if(rid==12735474)
		//{
		//	printf("line=%d, In %s(), rid=%lu, startmatchpos=%d, matchnum=%d, orientation=%d\n", __LINE__, __func__, rid, scafSuccessReadArr[i].startmatchpos, scafSuccessReadArr[i].matchnum, scafSuccessReadArr[i].orientation);
		//}
		/******************** 调试信息 ******************/

		if(scafSuccessReadArr[i].orientation==ORIENTATION_PLUS)
		{ //正向read
			//开始删除该read中的kmers
			//hashcode = first_hashcode;
			if(scafSuccessReadArr[i].matchnum==readLen-kmerSize+1)
			{
				first_kmer_index = Num - scafSuccessReadArr[i].matchnum;
				if(first_kmer_index<0)
				{
					printf("line=%d, In %s(), first_kmer_index=%d, error!\n", __LINE__, __func__, first_kmer_index);
				}

			}else
			{
				first_kmer_index = Num - scafSuccessReadArr[i].matchnum - 1;
				if(first_kmer_index<0)
				{
					printf("line=%d, In %s(), first_kmer_index=%d, error!\n", __LINE__, __func__, first_kmer_index);
				}
			}

			// generate the kmer integer sequence
			if(generateKmerSeqIntInScaf(tmp_kmerseq, scafContigSeqLastReadLen+first_kmer_index)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			rpos = 1;
			while(rpos<=scafSuccessReadArr[i].matchnum)
			{
				// 删除该read的一个kmer, 也即是该read的一个ridpos位置
				hashcode = kmerhashInScaf(tmp_kmerseq);
				if(delScafKmerByHashInScaf(hashcode, tmp_kmerseq, rid, rpos, scafGrapDeBruijn)==FAILED)
				{
					perfectMatchFlag = NO; //该read与contig并不是完全匹配
					if(delRemainedKmersInScaf(scafContigSeqLastReadLen+first_kmer_index+rpos-1, tmp_kmerseq, rid, rpos, scafGrapDeBruijn)==FAILED)
					{
						printf("line=%d, In %s(), can not delete the read [%lu, %d] from position %d. Error!\n", __LINE__, __func__, rid, scafSuccessReadArr[i].orientation, rpos);
						return FAILED;
					}

					break; //该read的删除操作终止
				}

				if(rpos==scafSuccessReadArr[i].matchnum)
					break;

				//取下一个kmer
				switch(scafContigSeqLastReadLen[first_kmer_index+kmerSize+rpos-1])
				{
					case 'A':
					case 'a':
						baseInt = 0;
						break;
					case 'C':
					case 'c':
						baseInt = 1;
						break;
					case 'G':
					case 'g':
						baseInt = 2;
						break;
					case 'T':
					case 't':
						baseInt = 3;
						break;
					default:
						printf("line=%d, In %s(), unknown base, Error!\n", __LINE__, __func__);
						return FAILED;
				}

				if(entriesPerKmer>=2)
				{
					for(j=0; j<entriesPerKmer-2; j++)
						tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
					tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNumKmer-2));
				}
				tmp_kmerseq[entriesPerKmer-1] = ((tmp_kmerseq[entriesPerKmer-1] << 2) | baseInt) & lastEntryMaskKmer;

				rpos++;
			}

		}else
		{ //反向read

			// generate the kmer integer sequence
			if(generateKmerSeqIntInScaf(tmp_kmerseq, reversed_lastseq)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			//开始删除该read中的kmers
			rpos = 1;
			while(rpos<=Num)
			{
				/******************** 调试信息 ******************/
				//if(rid==19399227 && rpos==24)
				//{
				//	printf("rid=%d,rpos=%d\n", rid, rpos);
				//}
				/******************** 调试信息 ******************/

				//删除该read的一个kmer, 也即是该read的一个ridpos位置
				hashcode = kmerhashInScaf(tmp_kmerseq);
				if(delScafKmerByHashInScaf(hashcode, tmp_kmerseq, rid, rpos, scafGrapDeBruijn)==FAILED)
				{
					perfectMatchFlag = NO; //该read与contig并不是完全匹配
					if(delRemainedKmersInScaf(reversed_lastseq+rpos-1, tmp_kmerseq, rid, rpos, scafGrapDeBruijn)==FAILED)
					{
						printf("line=%d, In %s(), can not delete the read [%lu, %d] from position %d. Error!\n", __LINE__, __func__, rid, scafSuccessReadArr[i].orientation, rpos);
						return FAILED;
					}

					break; //该read的删除操作终止
				}

				if(rpos==Num)
					break;

				//取下一个kmer
				switch(reversed_lastseq[kmerSize+rpos-1])
				{
					case 'A':
					case 'a':
						baseInt = 0;
						break;
					case 'C':
					case 'c':
						baseInt = 1;
						break;
					case 'G':
					case 'g':
						baseInt = 2;
						break;
					case 'T':
					case 't':
						baseInt = 3;
						break;
					default:
						printf("line=%d, In %s(), unknown base, Error!\n", __LINE__, __func__);
						return FAILED;
				}

				if(entriesPerKmer>=2)
				{
					for(j=0; j<entriesPerKmer-2; j++)
						tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
					tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNumKmer-2));
				}
				tmp_kmerseq[entriesPerKmer-1] = ((tmp_kmerseq[entriesPerKmer-1] << 2) | baseInt) & lastEntryMaskKmer;

				rpos++;
			}
		}

		//在lastseq36长度小于36个碱基的情况下, 将read中的其他kmer删除
		if(rpos<readLen-kmerSize+1 && perfectMatchFlag==YES)
		{
			rpos++;
			while(rpos<=readLen-kmerSize+1)
			{
				//下一个kmer的碱基的整数表示
				if(entriesPerKmer>=2)
				{
					for(j=0; j<entriesPerKmer-2; j++)
						tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
					tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNumKmer-2));
				}
				tmp_kmerseq[entriesPerKmer-1] = (tmp_kmerseq[entriesPerKmer-1] << 2) & lastEntryMaskKmer;

				for(j=0; j<4; j++)
				{ //循环探测该kmer, 直到找到并将其删除为止
					hashcode = kmerhashInScaf(tmp_kmerseq);
					if(delScafKmerByHashInScaf(hashcode, tmp_kmerseq, rid, rpos, scafGrapDeBruijn)==SUCCESSFUL)
						break;
					tmp_kmerseq[entriesPerKmer-1] ++;
				}

				//******************** 调试信息 **********************
				if(j==4)
				{
					printf("line=%d, In %s(), cannot delete the kmer (%lu,%d), Error!\n", __LINE__, __func__, rid, rpos);
					return FAILED;
				}
				//******************** 调试信息 **********************

				rpos++;
			}
		}

		i++;
	}

	free(reversed_lastseq);
	reversed_lastseq = NULL;

	return SUCCESSFUL;
}

/**
 * 删除read中从rpos开始的剩余的kmers.
 * 	@ return:
 * 		成功, 返回成功标记; 失败, 返回失败标记.
 */
short delRemainedKmersInScaf(char *seq, uint64_t *tmp_kmerseq, uint64_t rid, uint32_t rpos, scafGraph *graph)
{
	int len = strlen(seq);
	int i, j, baseInt;
	uint64_t hash;

	//删除该read的从rpos开始的第1个kmer
	for(j=0; j<3; j++)
	{
		if((tmp_kmerseq[entriesPerKmer-1] & 3)==3)
			tmp_kmerseq[entriesPerKmer-1] -= 3;
		else
			tmp_kmerseq[entriesPerKmer-1] ++;

		hash = kmerhashInScaf(tmp_kmerseq);
		if(delScafKmerByHashInScaf(hash, tmp_kmerseq, rid, rpos, graph)==SUCCESSFUL)
			break;
	}

	//******************** 调试信息 **********************
	if(j==3)
	{ //该kmer不能删除, 打印出错
		printf("line=%d, In %s(), can not delete this kmer (%lu, %d) Error!\n", __LINE__, __func__, rid, rpos);
		return FAILED;
	}
	//******************** 调试信息 **********************

	//删除该read的从rpos开始的其他kmer
	i = 1;
	rpos++;
	while(i+kmerSize-1<len)
	{
		//取下一个kmer
		switch(seq[i+kmerSize-1])
		{
			case 'A':
			case 'a':
				baseInt = 0;
				break;
			case 'C':
			case 'c':
				baseInt = 1;
				break;
			case 'G':
			case 'g':
				baseInt = 2;
				break;
			case 'T':
			case 't':
				baseInt = 3;
				break;
			default:
				printf("line=%d, In %s(), unknown base, Error!\n", __LINE__, __func__);
				return FAILED;
		}

		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
			}
			tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNumKmer-2));
		}
		tmp_kmerseq[entriesPerKmer-1] = ((tmp_kmerseq[entriesPerKmer-1] << 2) | baseInt) & lastEntryMaskKmer;

		//删除该kmer
		hash = kmerhashInScaf(tmp_kmerseq);
		if(delScafKmerByHashInScaf(hash, tmp_kmerseq, rid, rpos, graph)==FAILED)
		{ //删除该kmer失败, 该kmer需要重新确定
			//重新确定该kmer
			for(j=0; j<3; j++)
			{
				if((tmp_kmerseq[entriesPerKmer-1] & 3)==3)
					tmp_kmerseq[entriesPerKmer-1] -= 3;
				else
					tmp_kmerseq[entriesPerKmer-1] ++;

				hash = kmerhashInScaf(tmp_kmerseq);
				if(delScafKmerByHashInScaf(hash, tmp_kmerseq, rid, rpos, graph)==SUCCESSFUL)
					break;
			}

			//******************** 调试信息 **********************
			if(j==3)
			{ //该kmer不能删除, 打印出错
				printf("line=%d, In %s(), can not delete this kmer (%lu, %d) Error!\n", __LINE__, __func__, rid, rpos);
				return FAILED;
			}
			//******************** 调试信息 **********************
		}
		i++;
		rpos++;
	}

	//如果存在剩余的kmers, 则将其删除
	while(rpos <= readLen-kmerSize+1)
	{
		//下一个kmer的碱基的整数表示
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
			}
			tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNumKmer-2));
		}
		tmp_kmerseq[entriesPerKmer-1] = (tmp_kmerseq[entriesPerKmer-1] << 2) & lastEntryMaskKmer;

		for(j=0; j<4; j++)
		{ //循环探测该kmer, 直到找到并将其删除为止
			hash = kmerhashInScaf(tmp_kmerseq);
			if(delScafKmerByHashInScaf(hash, tmp_kmerseq, rid, rpos, graph)==SUCCESSFUL)
				break;
			hash ++;
		}

		//******************** 调试信息 **********************
		if(j==4)
		{
			printf("line=%d, In %s(), cannot delete the kmer (%lu,%d), Error!\n", __LINE__, __func__, rid, rpos);
			return FAILED;
		}
		//******************** 调试信息 **********************

		rpos++;
	}

	return SUCCESSFUL;
}

/**
    在contig中，追加拼接成功结束的reads的信息.
		@return:
			成功, 返回成功标记; 失败, 返回失败标记.
*/
short addRidposToContigInScaf(int contigNodesNum)
{
	scafContig *contig = scafContigLastReadLen;
	scafSuccessRead *ridposorient = NULL;
	int i, j = 0, num, k, Num = readLen - kmerSize + 1, successNum = 0, sumNewReads = 0, oldReadsNum, //成功的reads数量, 当次添加的新的reads的总数量, 该节点上的原有的reads数量
			checkContigNum = 0; //限制检测contig的数量
	short *indicator;

	indicator = (short *) malloc (scafSuccessReadsNum * sizeof(short));
	if(indicator==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory for indicator array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//如果contig节点个数大于36, 则contig36所指向的节点个数为36
	if(contigNodesNum>readLen)
		contigNodesNum = readLen;

	i = 0;
	while(contig)
	{
		i++;

		checkContigNum ++;
		if(checkContigNum>errorRegLenEnd3+1)  //限制检测碱基数量为5
			break;

		if(memset(indicator, 0, scafSuccessReadsNum*sizeof(short))==NULL)
		{
			printf("line=%d, In %s(), cannot reset indicator, error!\n", __LINE__, __func__);
			return FAILED;
		}

		num = 0;
		j = 0;
		//newReadsNum = 0; //设该节点上的新的reads的数量为0
		oldReadsNum = 0; //设该节点上的原有的reads数量为0
		while(j<scafSuccessReadsNum)
		{ //查找该contig中的read个数

			//******************* 调试信息 ******************
			//if(ridposorientation[j].rid==942710)
			//{
			//	printf("line=%d, In %s(), rid=%d, startmatchpos=%d, matchnum=%d, orientation=%c, i=%d\n", __LINE__, __func__, scafSuccessReadArr[j].rid, scafSuccessReadArr[j].startmatchpos, scafSuccessReadArr[j].matchnum, scafSuccessReadArr[j].orientation, i);
			//}
			//******************* 调试信息 ******************

			//if((scafSuccessReadArr[j].orientation==ORIENTATION_PLUS && scafSuccessReadArr[j].pos==i)
					//||(scafSuccessReadArr[j].orientation==ORIENTATION_MINUS && contigNodesNum-KMER_SIZE+2-scafSuccessReadArr[j].pos==i))
			//if((i==1 && scafSuccessReadArr[j].orientation==ORIENTATION_PLUS)
					//|| (i==Num-scafSuccessReadArr[j].matchnum+1 && scafSuccessReadArr[j].orientation==ORIENTATION_MINUS)) --2010年8月22日23点注释掉
			if((i==1 && scafSuccessReadArr[j].matchnum==Num)
					|| (i==contigNodesNum-kmerSize+1-scafSuccessReadArr[j].matchnum && scafSuccessReadArr[j].orientation==ORIENTATION_PLUS)
					|| (i==contigNodesNum-kmerSize+1-scafSuccessReadArr[j].matchnum+1 && scafSuccessReadArr[j].orientation==ORIENTATION_MINUS))
			{
				num++;
				indicator[j] = 1;
			}
			j++;
		}

		successNum += oldReadsNum; //将原有的reads算作成功的数量

		//******************** 调试信息 *********************
		if(successNum>scafSuccessReadsNum)
		{
			printf("line=%d, In %s(), successNum(%d) > scafSuccessReadsNum(%d), Error!\n", __LINE__, __func__, successNum, scafSuccessReadsNum);
			return FAILED;
		}
		//******************** 调试信息 *********************

		if(num>0)
		{ //如果read个数大于0，则添加这些reads的信息
			//如果该contig节点已经有成功的reads, 则需要将这些已经成功的reads合并到新的reads集合中
			ridposorient = (scafSuccessRead*) malloc((contig->ridposnum + num)*sizeof(scafSuccessRead));
			if(ridposorient==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(contig->ridposnum>0)
			{ //该contig已经有成功的reads, 则将其复制到新开辟的数组中
				//复制原来的成功的reads
				if(memcpy(ridposorient, contig->pridposorientation, contig->ridposnum*sizeof(scafSuccessRead))==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				//释放掉原来的位置数组占用的内存
				free(contig->pridposorientation);
				contig->pridposorientation = NULL;
			}

			//将开辟出来的ridposorientation结构添加到contig中
			contig->pridposorientation = ridposorient;

			j = 0;
			k = contig->ridposnum;
			while(j<scafSuccessReadsNum) //添加新的成功的reads
			{ //循环遍历成功的reads表，将以该contig开始的reads添加进该contig
				if(indicator[j]==1)
				{ //该read需要添加到该contig节点

					if(scafSuccessReadArr[j].orientation==ORIENTATION_PLUS)
					{ //正向的reads, 位置直接赋值
						ridposorient[k].startmatchpos = scafSuccessReadArr[j].startmatchpos;
					}else
					{ //反向的reads, 位置需要变换
						ridposorient[k].startmatchpos = scafSuccessReadArr[j].startmatchpos + kmerSize - 1;
					}
					ridposorient[k].orientation = scafSuccessReadArr[j].orientation;
					ridposorient[k].rid = scafSuccessReadArr[j].rid;
					//ridposorient[k].matchnum = scafSuccessReadArr[j].matchnum; //按照read中的kmer位置计算
					ridposorient[k].matchnum = scafSuccessReadArr[j].matchnum + kmerSize - 1; //按照read中的碱基位置计算

					//统计新的reads数量
					sumNewReads ++;

					successNum ++;

					k++;
				}
				j++;
			}

			contig->ridposnum = k; //更新成功的reads数量

			if(successNum==scafSuccessReadsNum)
			{ //所有的reads都已经添加完毕, 则退出循环
				break;
			}

		} //end if(num>0)

		contig = contig->next; //指向下一个contig
	} //end while(contig)

	scafSuccessReadsNum = sumNewReads; //新有的reads数量

	free(indicator);
	indicator = NULL;

	return SUCCESSFUL;
}

/**
 * 取得successContig的地址.
 *
 *   @ return:
 *     找到, 返回successContig地址; 否则, 返回NULL.
 */
scafContig *getSuccessContigInScaf(int contigIndex)
{
	scafContig *contig, *new_successContig = NULL;
	scafSuccessRead *ridposorient = NULL;
	int i, numridposorient, tmp_contigIndex;

	int checkNum = 0;
	tmp_contigIndex = 0;
	contig = scafContigLastReadLen;
	while(contig)
	{ //遍历contig链表
		checkNum ++;
		if(checkNum>errorRegLenEnd3+1)  //限制检测5次
			break;

		if(contig->ridposnum > 0)
		{ //找到标记有成功reads的contig节点
			tmp_contigIndex = contig->index;
			break;
		}

		contig = contig->next;
	}

	if(successScafContig!=NULL && tmp_contigIndex>0 && tmp_contigIndex>successScafContig->index-MIN_OVERLAP_LEN+1)
	{
		return NULL;
	}


	int maxContigIndex = 0;
	int maxMatchNum = 0;
	checkNum = 0;
	tmp_contigIndex = 1;
	contig = scafContigLastReadLen;
	while(contig)
	{ //遍历contig链表
		checkNum ++;
		if(checkNum>errorRegLenEnd3+1)  //限制检测5次
			break;

		if(contig->ridposnum > 0)
		{ //找到标记有成功reads的contig节点
			maxMatchNum = 0;
			ridposorient = contig->pridposorientation;
			numridposorient = contig->ridposnum;
			for(i=0; i<numridposorient; i++)
			{
				if(maxMatchNum < ridposorient[i].matchnum)
				{
					maxMatchNum = ridposorient[i].matchnum;
				}
			}

			//找最大的successContigIndex
			if(tmp_contigIndex+maxMatchNum-1 > maxContigIndex)
				maxContigIndex = tmp_contigIndex + maxMatchNum - 1;
		}
		contig = contig->next;
		tmp_contigIndex ++;
	}

	//*********************** 调试信息 **********************
	if(maxContigIndex<=0)
	{
		printf("line=%d, In %s(), the maxContigIndex==%d\n", __LINE__, __func__, maxContigIndex);
		return NULL;
	}
	//*********************** 调试信息 **********************

	//找最大index的成功contig节点
	contig = scafContigLastReadLen;
	tmp_contigIndex = 1;
	while(contig)
	{
		if(tmp_contigIndex==maxContigIndex)
		{ //找到最大index的成功的reads的contig节点
			new_successContig = contig;
			break;
		}
		contig = contig->next;
		tmp_contigIndex ++;
	}

	//检测新的successContig是否合理
	if(successScafContig==NULL)
	{
		return new_successContig;
	}else if(successScafContig->index > new_successContig->index)
	{
		return successScafContig;
	}else
	{
		return new_successContig;
	}
}

/**
 * 拼接结束时在尾部更新contig节点.
 *   也就是根据成功reads间隔与contigIndex, 将contig尾部的节点删除.
 *
 *   @ return:
 *     成功, 返回成功标记; 失败, 返回失败标记.
 */
short updateContigtailnodesInScaf(int *scafContigIndex)
{
	scafContig *tmp_contig;
	int startIndex = successScafContig->index - readLen + 1;
	int i, j, divisionIndex;
	scafSuccessRead *ridposorient, *tmp_ridposorient;
	int numridposorient;

	if(startIndex<=0)
		startIndex = 1;

	tmp_contig = scafContighead;
	while(tmp_contig)
	{
		if(tmp_contig->index==startIndex)
			break;

		tmp_contig = tmp_contig->next;
	}

	divisionIndex = successScafContig->index;
	int delNum;
	while(tmp_contig)
	{
		if(tmp_contig->index>divisionIndex)
			break;

		if(tmp_contig->ridposnum>0)
		{
			delNum = 0;
			ridposorient = tmp_contig->pridposorientation;
			numridposorient = tmp_contig->ridposnum;
			short indicator[numridposorient];
			if(memset(indicator, 0, sizeof(short)*numridposorient)==NULL)
			{
				printf("line=%d, In %s(), cannot memset indicator to zero, error!\n", __LINE__, __func__);
				return FAILED;
			}

			for(i=0; i<numridposorient; i++)
			{
				if(tmp_contig->index+ridposorient->matchnum-1 > divisionIndex)
				{
					indicator[i] = 1;
					delNum ++;
				}
			}

			if(delNum==numridposorient)
			{
				free(tmp_contig->pridposorientation);
				tmp_contig->pridposorientation = NULL;
				tmp_contig->ridposnum = 0;

			}else if(delNum>0)
			{
				tmp_ridposorient = (scafSuccessRead*) malloc(sizeof(scafSuccessRead)*(numridposorient-delNum));
				if(tmp_ridposorient==NULL)
				{
					printf("line=%d, In %s(), cannot malloc tmp_ridposorient, error!\n", __LINE__, __func__);
					return FAILED;
				}
				tmp_contig->ridposnum = numridposorient - delNum;
				tmp_contig->pridposorientation = tmp_ridposorient;

				j = 0;
				for(i=0; i<numridposorient; i++)
				{
					if(indicator[i]==0)
					{
						if(memcpy(tmp_ridposorient+j, ridposorient+i, sizeof(scafSuccessRead))==NULL)
						{
							printf("line=%d, In %s(), cannot memcpy ridposorient, error!\n", __LINE__, __func__);
							return FAILED;
						}
						j ++;
					}
				}

				free(ridposorient);
			}
		}

		tmp_contig = tmp_contig->next;
	}

	int count = 0;
	while(successScafContig->next!=NULL)
	{
		tmp_contig = successScafContig->next;
		successScafContig->next = tmp_contig->next;
		if(tmp_contig->ridposnum>0)
			free(tmp_contig->pridposorientation);
		free(tmp_contig);
		count ++;
	}
	(*scafContigIndex) -= count;

	return SUCCESSFUL;
}

/**
 * Update contig end sequences in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateScafContigEndSeqs(int assemblyRound, int scafContigNodesNum)
{
	char *pEndSeq;
	int i, *pEndSeqLen;
	int prepareScafContigLen, startScafContigIndex, startEndSeqIndex;
	char base;
	scafContig *tmp_contig;

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{ // the first round
		pEndSeq = scafContigEndSeqArr[0];
		pEndSeqLen = scafContigEndSeqLenArr;
		prepareScafContigLen = prepareAssemblyLenArr[0];
	}else
	{ // the second round
		pEndSeq = scafContigEndSeqArr[1];
		pEndSeqLen = scafContigEndSeqLenArr + 1;
		prepareScafContigLen = prepareAssemblyLenArr[1];
	}


	if(scafContigNodesNum>=maxOverlapSeqLen)
	{
		startEndSeqIndex = -1; // -1 indicates that the index is invalid
		startScafContigIndex = scafContigNodesNum - maxOverlapSeqLen + 1;
	}else
	{
		if((*pEndSeqLen)+scafContigNodesNum-prepareScafContigLen >= maxOverlapSeqLen)
			startEndSeqIndex = (*pEndSeqLen) + scafContigNodesNum - prepareScafContigLen - maxOverlapSeqLen;
		else
			startEndSeqIndex = 0;

		startScafContigIndex = prepareScafContigLen + 1;
	}

	// update the scafContig end sequence
	if(startEndSeqIndex==-1)
	{
		i = 0;
	}else if(startEndSeqIndex==0)
	{
		i = (*pEndSeqLen);
	}else if(startEndSeqIndex>0)
	{
		for(i=startEndSeqIndex; i<(*pEndSeqLen); i++) pEndSeq[i-startEndSeqIndex] = pEndSeq[i];
		pEndSeq[(*pEndSeqLen)-startEndSeqIndex] = '\0';
		i = (*pEndSeqLen) - startEndSeqIndex;
	}else
	{
		printf("line=%d, In %s(), startEndSeqIndex=%d, error!\n", __LINE__, __func__, startEndSeqIndex);
		return FAILED;
	}

	tmp_contig = scafContighead;
	while(tmp_contig)
	{
		if(tmp_contig->index>=startScafContigIndex)
		{
			switch(tmp_contig->base)
			{
				case 0: base = 'A'; break;
				case 1: base = 'C'; break;
				case 2: base = 'G'; break;
				case 3: base = 'T'; break;
				default: printf("line=%d, In %s(), unknown base [ %d ], error!\n", __LINE__, __func__, tmp_contig->base); return FAILED;
			}

			pEndSeq[i++] = base;
		}

		tmp_contig = tmp_contig->next;
	}
	pEndSeq[i] = '\0';
	*pEndSeqLen = i;

	// ############################## Debug information ###########################
	if((*pEndSeqLen) != strlen(pEndSeq))
	{
		printf("line=%d, In %s(), (*pEndSeqLen)=%d != strlen(pEndSeq)=%d, error!\n", __LINE__, __func__, *pEndSeqLen, (int)strlen(pEndSeq));
		return FAILED;
	}
	// ############################## Debug information ###########################

	return SUCCESSFUL;
}

/**
 * Detect overlaps between contig ends in local zssembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short detectOverlapsInScaf(int *successFilledFlag, int *overlapLen, int *newGapSize, int *breakFlag, int gapSize, int *localScafContigNodesNumArray, int assemblyRound)
{
	char *contigEndSeq, *comparisonSeq, *otherContigEndSeq;
	int *contigEndSeqLen, comparisonSeqLen, *otherContigEndSeqLen;
	int i, j, overlapLenExact, overlapLenAlignment, mismatchNum, overlapLenAdjust, tmpGapSize;

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{ // the first round
		contigEndSeq = scafContigEndSeqArr[0];
		contigEndSeqLen = scafContigEndSeqLenArr;
		otherContigEndSeq = scafContigEndSeqArr[1];
		otherContigEndSeqLen = scafContigEndSeqLenArr + 1;
	}else
	{ // the second round
		contigEndSeq = scafContigEndSeqArr[1];
		contigEndSeqLen = scafContigEndSeqLenArr + 1;
		otherContigEndSeq = scafContigEndSeqArr[0];
		otherContigEndSeqLen = scafContigEndSeqLenArr;
	}
	comparisonSeq = comparisonSeqInScaf;
	comparisonSeqLen = comparisonSeqLenInScaf;

	// compute the new gap size
	if(computeNewGapSizeInScaf(&tmpGapSize, gapSize, localScafContigNodesNumArray, assemblyRound)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute new gap size, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if((gapSize>0.6*meanSizeInsert && gapSize>2*averReadLen) && tmpGapSize>averReadLen) // 2012-11-19
	{
#if (DEBUG_FLAG==YES)
		printf("line=%d, In %s(), broken.\n", __LINE__, __func__);
#endif
		*overlapLen = 0;
		*newGapSize = tmpGapSize;
		*successFilledFlag = NO;
		*breakFlag = YES;
	}
	else if(tmpGapSize<=stardardDeviationInsert)
	{
		// compute the overlap length by exact alignment
		if(computeSeqOverlapLenExact(&overlapLenExact, contigEndSeq, *contigEndSeqLen, comparisonSeq, comparisonSeqLen, scoreArr, tmpGapSize)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute sequence overlap length by sequence comparison, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// check the contig overlap
		//if(overlapLenExact>=minOverlapThres)
		//if(overlapLenExact>=minOverlapThres  && (overlapLenExact>=-tmpGapSize-stardardDeviationInsert && overlapLenExact<=-tmpGapSize+stardardDeviationInsert))
		if(overlapLenExact>=minOverlapThres && tmpGapSize>-gapSizeSdevFactorGapFilling*stardardDeviationInsert)
		{ // exact overlap is found by sequences comparison
			*overlapLen = overlapLenExact;
			*newGapSize = 0;
			*successFilledFlag = YES;
			*breakFlag = NO;
		}else
		{ // alignment the two sequences
			// pairwise alignment
			if(computeSeqOverlapLenByAlignment(contigEndSeq, *contigEndSeqLen, comparisonSeq, comparisonSeqLen, scoreArr, alignResultArr, &overlapLenAlignment, &mismatchNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot detect sequence overlap length by sequence alignment, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(overlapLenAlignment>=minOverlapThres && mismatchNum <= mismatchThres)
			{ // there is an overlap between the two contigs
				overlapLenAdjust = overlapLenAlignment;
				// adjust overlapped sequences
				if(adjustOverlapSeq(contigEndSeq, comparisonSeq, alignResultArr, &overlapLenAdjust)==FAILED)
				{
					printf("line=%d, In %s(), cannot adjust the overlapped sequences, error!\n", __LINE__, __func__);
					return FAILED;
				}

				// update the contig sequences
				if(overlapLenAdjust>0)
				{ // valid overlap
					*overlapLen = overlapLenAdjust;
					*newGapSize = 0;
					*successFilledFlag = YES;
					*breakFlag = NO;

					*contigEndSeqLen = strlen(contigEndSeq);
					comparisonSeqLenInScaf = strlen(comparisonSeq);

					// update the other contig end sequence
					for(i=comparisonSeqLenInScaf-1, j=0; i>=0; i--, j++)
					{
						switch(comparisonSeq[i])
						{
							case 'A': otherContigEndSeq[j] = 'T'; break;
							case 'C': otherContigEndSeq[j] = 'G'; break;
							case 'G': otherContigEndSeq[j] = 'C'; break;
							case 'T': otherContigEndSeq[j] = 'A'; break;
							default: printf("line=%d, In %s(), unknown base %c, error!\n", __LINE__, __func__, comparisonSeq[i]);
									return FAILED;
						}
					}
					otherContigEndSeq[comparisonSeqLenInScaf] = '\0';
					*otherContigEndSeqLen = comparisonSeqLenInScaf;
				}else
				{ // there is a gap or a short overlap between the two contigs
					// detect the short overlap of sequence comparison
					if(tmpGapSize<minAdjustGapSizeThres)
					{ // gapSize < -10

						*overlapLen = 0;
						*newGapSize = tmpGapSize;
						*successFilledFlag = NO;
						*breakFlag = NO;

					}else if(tmpGapSize<maxAdjustGapSizeThres)
					{ // -10 =< gapSize < 10
						// detect the short overlap or the short gap
						if(overlapLenAlignment>0 && mismatchNum==0)
						{
							*overlapLen = overlapLenAlignment;
							*newGapSize = 0;
							*successFilledFlag = YES;
							*breakFlag = NO;

						//}else if(overlapLenExact>0)
						}else if(overlapLenExact>=minExactOverlapThres)
						{
							*overlapLen = overlapLenExact;
							*newGapSize = 0;
							*successFilledFlag = YES;
							*breakFlag = NO;

						}else
						{
							*overlapLen = 0;
							*newGapSize = tmpGapSize;
							*successFilledFlag = NO;
							*breakFlag = NO;
						}

					}else
					{ // gapSize >= 10
						// there is a gap
						*overlapLen = 0;
						*newGapSize = tmpGapSize;
						*successFilledFlag = NO;
						*breakFlag = NO;
					}

				}
			}else // if(overlapLenAlignment<minOverlapThres || mismatchNum > mismatchThres)
			{ // there is a gap between the two contigs with high probability
				// detect the short overlap of sequence comparison
				if(tmpGapSize<minAdjustGapSizeThres)
				{ // gapSize < -10

					*overlapLen = 0;
					*newGapSize = tmpGapSize;
					*successFilledFlag = NO;

					if(tmpGapSize<-0.5*meanSizeInsert)
						*breakFlag = YES;
					else
						*breakFlag = NO;

				}else if(tmpGapSize<maxAdjustGapSizeThres)
				{ // -10 =< gapSize < 10
					// detect the short overlap or the short gap
					if(overlapLenAlignment>0 && mismatchNum==0)
					{
						*overlapLen = overlapLenAlignment;
						*newGapSize = 0;
						*successFilledFlag = YES;
						*breakFlag = NO;

					//}else if(overlapLenExact>0)
					}else if(overlapLenExact>=minExactOverlapThres)
					{
						*overlapLen = overlapLenExact;
						*newGapSize = 0;
						*successFilledFlag = YES;
						*breakFlag = NO;

					}else
					{
						*overlapLen = 0;
						*newGapSize = tmpGapSize;
						*successFilledFlag = NO;
						*breakFlag = NO;
					}
				}else
				{ // gapSize >= 10
					// there is a gap
					*overlapLen = 0;
					*newGapSize = tmpGapSize;
					*successFilledFlag = NO;
					*breakFlag = NO;
				}
			}
		}
	}else
	{
		// there is a gap
		*overlapLen = 0;
		*newGapSize = tmpGapSize;
		*successFilledFlag = NO;
		*breakFlag = NO;
	}

	return SUCCESSFUL;
}

/**
 * Compute new gap size in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeNewGapSizeInScaf(int *tmpGapSize, int gapSize, int *localScafContigNodesNumArray, int assemblyRound)
{
	int newSeqLen1, newSeqLen2;

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{ // the first round
		newSeqLen1 = localScafContigNodesNumArray[0] - prepareAssemblyLenArr[0];
		*tmpGapSize = gapSize - newSeqLen1;

	}else
	{ // the second round
		newSeqLen1 = localScafContigNodesNumArray[0] - prepareAssemblyLenArr[0];
		newSeqLen2 = localScafContigNodesNumArray[1] - prepareAssemblyLenArr[1];
		*tmpGapSize = gapSize - newSeqLen1- newSeqLen2;
	}

	return SUCCESSFUL;
}

/**
 * Update contig overlap information in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateContigOverlapInfoInScaf(contigOverlap *pContigOverlapInfo, int successFilledFlag, int overlapLen, int newGapSize, int breakFlag)
{
	if(pContigOverlapInfo==NULL)
	{
		printf("line=%d, In %s(), the contigOverlap pointer is NULL, error!\n", __LINE__, __func__);
		return FAILED;
	}else if(pContigOverlapInfo->mergeFlag==YES || pContigOverlapInfo->breakFlag==YES)
	{
		printf("line=%d, In %s(), the contig overlap information is invalid, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(successFilledFlag==YES)
	{
		pContigOverlapInfo->overlapLen = overlapLen;
		pContigOverlapInfo->mergeFlag = YES;
		pContigOverlapInfo->gapSize = 0;
		pContigOverlapInfo->breakFlag = NO;
	}else
	{
		if(breakFlag==NO)
		{
			pContigOverlapInfo->overlapLen = 0;
			pContigOverlapInfo->mergeFlag = NO;
			pContigOverlapInfo->gapSize = newGapSize;
			pContigOverlapInfo->breakFlag = NO;
		}else
		{
#if (DEBUG_FLAG==YES)
		printf("line=%d, In %s(), broken.\n", __LINE__, __func__);
#endif
			pContigOverlapInfo->overlapLen = 0;
			pContigOverlapInfo->mergeFlag = NO;
			pContigOverlapInfo->gapSize = 0;
			pContigOverlapInfo->breakFlag = YES;
		}
	}

	return SUCCESSFUL;
}

/**
 * Update contig information in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateContigInfoInScaf(contigOverlap *pContigOverlapInfo, scafContig **localScafContigheadArr, int *localScafContigNodesNumArr, int *oldEndSeqLenArr, int assemblyRound)
{
	int contigID1, contigID2, contigOrient1, contigOrient2;

	// get the overlap information
	contigID1 = pContigOverlapInfo->contigID1;
	contigID2 = pContigOverlapInfo->contigID2;
	contigOrient1 = pContigOverlapInfo->orientation1;
	contigOrient2 = pContigOverlapInfo->orientation2;

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{ // the first round
		if(localScafContigNodesNumArr[0] > prepareAssemblyLenArr[0])
		{ // new bases are appended
			if(updateSingleContigInfoInScaf(contigInfoArr+contigID1-1, contigOrient1, 0, localScafContigheadArr[0], localScafContigNodesNumArr[0], scafContigEndSeqArr[0], scafContigEndSeqLenArr[0], oldEndSeqLenArr[0], prepareAssemblyLenArr[0])==FAILED)
			{
				printf("line=%d, In %s(), cannot update single contig information in local assembly for contig %d, error!\n", __LINE__, __func__, contigID1);
				return FAILED;
			}

			// update the other contig information
			if(updateOtherContigInfoInScaf(contigInfoArr+contigID2-1, contigOrient2, scafContigEndSeqArr[1], scafContigEndSeqLenArr[1], oldEndSeqLenArr[1])==FAILED)
			{
				printf("line=%d, In %s(), cannot update the other contig information in local assembly for contig %d, error!\n", __LINE__, __func__, contigID2);
				return FAILED;
			}
		}else
		{
			printf("line=%d, In %s(), no bases are appended when the gap are successfully filled, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{ // the second round

		if(localScafContigNodesNumArr[0] <= prepareAssemblyLenArr[0] && localScafContigNodesNumArr[1] <= prepareAssemblyLenArr[1])
		{
			printf("line=%d, In %s(), assemblyRound=%d, no bases are appended, error!\n", __LINE__, __func__, assemblyRound);
			return FAILED;
		}

		// process the first contig
		if(updateSingleContigInfoInScaf(contigInfoArr+contigID1-1, contigOrient1, 0, localScafContigheadArr[0], localScafContigNodesNumArr[0], scafContigEndSeqArr[0], scafContigEndSeqLenArr[0], oldEndSeqLenArr[0], prepareAssemblyLenArr[0])==FAILED)
		{
			printf("line=%d, In %s(), cannot update single contig information in local assembly for contig %d, error!\n", __LINE__, __func__, contigID1);
			return FAILED;
		}

		// process the second contig
		if(updateSingleContigInfoInScaf(contigInfoArr+contigID2-1, contigOrient2, 1, localScafContigheadArr[1], localScafContigNodesNumArr[1], scafContigEndSeqArr[1], scafContigEndSeqLenArr[1], oldEndSeqLenArr[1], prepareAssemblyLenArr[1])==FAILED)
		{
			printf("line=%d, In %s(), cannot update single contig information in local assembly for contig %d, error!\n", __LINE__, __func__, contigID2);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Update single contig information in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateSingleContigInfoInScaf(contigInfo *pContigInfo, int contigOrient, int contigIndex, scafContig *scafContighead, int scafContigNodesNum, char *endSeq, int endSeqLen, int oldEndSeqLen, int prepareAssemblyLen)
{
	int i, j;
	int contigLen;
	char *contigSeq, *newContigSeq, *tmpSeq, tmpBase;
	int newEndSeqLen, tmpSeqLen;
	int newContigLen;

	int endScafContigIndex, startEndSeqIndex;
	scafContig *tmp_contig;

	contigLen = pContigInfo->contigLen;
	contigSeq = pContigInfo->contigSeq;

	// allocate the memory for tmpStr
	newEndSeqLen = endSeqLen;
	tmpSeqLen = scafContigNodesNum + newEndSeqLen - oldEndSeqLen;
	tmpSeq = (char *) malloc((tmpSeqLen+1)* sizeof(char));
	if(tmpSeq==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(scafContigNodesNum>=maxOverlapSeqLen)
	{
		endScafContigIndex = scafContigNodesNum - maxOverlapSeqLen;
		startEndSeqIndex = 0;
	}else
	{
		endScafContigIndex = -1;
		startEndSeqIndex = oldEndSeqLen - scafContigNodesNum;
	}

	// copy the new end sequence to tmpSeq
	if(endScafContigIndex<=0)
	{
		strcpy(tmpSeq, endSeq+startEndSeqIndex);
	}else
	{
		i = 0;
		tmp_contig = scafContighead;
		while(tmp_contig)
		{
			if(tmp_contig->index<=endScafContigIndex)
			{
				switch(tmp_contig->base)
				{
					case 0: tmpBase = 'A'; break;
					case 1: tmpBase = 'C'; break;
					case 2: tmpBase = 'G'; break;
					case 3: tmpBase = 'T'; break;
					default: printf("line=%d, In %s(), unknown base %d, error!\n", __LINE__, __func__, tmp_contig->base); return FAILED;
				}
				tmpSeq[ i++ ] = tmpBase;
			}else
			{
				break;
			}

			tmp_contig = tmp_contig->next;
		}
		tmpSeq[i] = '\0';

		strcpy(tmpSeq+i, endSeq);

		// ############################# Debug information ###########################
		if(strlen(tmpSeq)!=tmpSeqLen)
		{
			printf("line=%d, In %s(), strlen(tmpSeq)=%d != tmpSeqLen=%d, error!\n", __LINE__, __func__, (int)strlen(tmpSeq), tmpSeqLen);
			return FAILED;
		}
		// ############################# Debug information ###########################
	}

	// update the contig base sequence
	if(contigIndex==0)
	{ // update the first contig
		if(contigOrient==ORIENTATION_PLUS)
		{ // plus orientation
			newContigLen = contigLen + tmpSeqLen - prepareAssemblyLen;
			if(newContigLen==contigLen)
			{ // only update end sequence
				strcpy(contigSeq+contigLen-prepareAssemblyLen, tmpSeq);

			}else if(newContigLen<contigLen)
			{ // update the end sequence and remained sequence
				strcpy(contigSeq+contigLen-prepareAssemblyLen, tmpSeq);
				pContigInfo->contigLen = newContigLen;
			}else
			{ // allocate memory for new contig sequence, and update its sequence
				newContigSeq = (char *) malloc((newContigLen+1)* sizeof(char));
				if(newContigSeq==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}

				strcpy(newContigSeq, contigSeq);
				strcpy(newContigSeq+contigLen-prepareAssemblyLen, tmpSeq);

				// update the sequence pointer and contig length
				free(pContigInfo->contigSeq);
				pContigInfo->contigSeq = newContigSeq;
				pContigInfo->contigLen = newContigLen;
			}
		}else
		{ // minus orientation
			// get the reverse complements of new contig sequence
			if(reverseSeq(tmpSeq, tmpSeqLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot reverse new contig sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			newContigLen = contigLen + tmpSeqLen - prepareAssemblyLen;
			if(newContigLen==contigLen)
			{ // only update end sequence
				for(i=0; i<tmpSeqLen; i++) contigSeq[i] = tmpSeq[i];

			}else if(newContigLen<contigLen)
			{ // update the end sequence and remained sequence
				for(i=0; i<tmpSeqLen; i++) contigSeq[i] = tmpSeq[i];
				//strcpy(contigSeq+tmpSeqLen, contigSeq+prepareAssemblyLen);
				for(i=tmpSeqLen, j=prepareAssemblyLen; i<newContigLen; i++, j++) contigSeq[i] = contigSeq[j];
				contigSeq[newContigLen] = '\0';

				// update contig length
				pContigInfo->contigLen = newContigLen;

			}else
			{ // allocate memory for new contig sequence, and update its sequence
				newContigSeq = (char *) malloc((newContigLen+1)* sizeof(char));
				if(newContigSeq==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}

				strcpy(newContigSeq, tmpSeq);
				strcpy(newContigSeq+tmpSeqLen, contigSeq+prepareAssemblyLen);

				// update the sequence pointer and contig length
				free(pContigInfo->contigSeq);
				pContigInfo->contigSeq = newContigSeq;
				pContigInfo->contigLen = newContigLen;
			}
		}
	}else
	{ // update the second contig
		if(contigOrient==ORIENTATION_PLUS)
		{ // plus orientation

			// get the reverse complements of new contig sequence
			if(reverseSeq(tmpSeq, tmpSeqLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot reverse new contig sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			newContigLen = contigLen + tmpSeqLen - prepareAssemblyLen;
			if(newContigLen==contigLen)
			{ // only update end sequence
				for(i=0; i<tmpSeqLen; i++) contigSeq[i] = tmpSeq[i];

			}else if(newContigLen<contigLen)
			{ // update the end sequence and remained sequence
				for(i=0; i<tmpSeqLen; i++) contigSeq[i] = tmpSeq[i];
				//strcpy(contigSeq+tmpSeqLen, contigSeq+prepareAssemblyLen);
				for(i=tmpSeqLen, j=prepareAssemblyLen; i<newContigLen; i++, j++) contigSeq[i] = contigSeq[j];
				contigSeq[newContigLen] = '\0';

				// update contig length
				pContigInfo->contigLen = newContigLen;

			}else
			{ // allocate memory for new contig sequence, and update its sequence
				newContigSeq = (char *) malloc((newContigLen+1)* sizeof(char));
				if(newContigSeq==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}

				strcpy(newContigSeq, tmpSeq);
				strcpy(newContigSeq+tmpSeqLen, contigSeq+prepareAssemblyLen);

				// update the sequence pointer and contig length
				free(pContigInfo->contigSeq);
				pContigInfo->contigSeq = newContigSeq;
				pContigInfo->contigLen = newContigLen;
			}
		}else
		{ // minus orientation
			newContigLen = contigLen + tmpSeqLen - prepareAssemblyLen;
			if(newContigLen==contigLen)
			{ // only update end sequence
				strcpy(contigSeq+contigLen-prepareAssemblyLen, tmpSeq);

			}else if(newContigLen<contigLen)
			{ // update the end sequence and remained sequence
				strcpy(contigSeq+contigLen-prepareAssemblyLen, tmpSeq);
				pContigInfo->contigLen = newContigLen;

			}else
			{ // allocate memory for new contig sequence, and update its sequence
				newContigSeq = (char *) malloc((newContigLen+1)* sizeof(char));
				if(newContigSeq==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}

				strcpy(newContigSeq, contigSeq);
				strcpy(newContigSeq+contigLen-prepareAssemblyLen, tmpSeq);

				// update the sequence pointer and contig length
				free(pContigInfo->contigSeq);
				pContigInfo->contigSeq = newContigSeq;
				pContigInfo->contigLen = newContigLen;
			}
		}
	}

	// ############################# Debug information ###########################
	if(strlen(pContigInfo->contigSeq)!=pContigInfo->contigLen)
	{
		printf("line=%d, In %s(), strlen(contigSeq)=%d != contigLen=%d, error!\n", __LINE__, __func__, (int)strlen(pContigInfo->contigSeq), pContigInfo->contigLen);
		return FAILED;
	}
	// ############################# Debug information ###########################

	free(tmpSeq);
	tmpSeq = NULL;

	return SUCCESSFUL;
}

/**
 * Update other contig end sequence and contig length.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateOtherContigInfoInScaf(contigInfo *pContigInfo, int contigOrient, char *contigEndSeq, int endSeqLen, int oldEndSeqLen)
{
	int i, j, newEndSeqLen;
	int contigLen, newContigLen;
	char *contigSeq, *newContigSeq;

	contigLen = pContigInfo->contigLen;
	contigSeq = pContigInfo->contigSeq;

	newEndSeqLen = endSeqLen;

	if(contigOrient==ORIENTATION_PLUS)
	{ // plus orientation
		if(newEndSeqLen==oldEndSeqLen)
		{ // only update the end sequence
			for(i=0; i<newEndSeqLen; i++)
			{
				switch(contigEndSeq[i])
				{
					case 'A': contigSeq[newEndSeqLen-i-1] = 'T'; break;
					case 'C': contigSeq[newEndSeqLen-i-1] = 'G'; break;
					case 'G': contigSeq[newEndSeqLen-i-1] = 'C'; break;
					case 'T': contigSeq[newEndSeqLen-i-1] = 'A'; break;
					default: printf("line=%d, In %s(), unknown base %c, error!\n", __LINE__, __func__, contigEndSeq[i]);
							return FAILED;
				}
			}
		}else if(newEndSeqLen<oldEndSeqLen)
		{ // update the end sequence and remained sequence
			newContigLen = contigLen - oldEndSeqLen + newEndSeqLen;
			for(i=0; i<newEndSeqLen; i++)
			{
				switch(contigEndSeq[i])
				{
					case 'A': contigSeq[newEndSeqLen-i-1] = 'T'; break;
					case 'C': contigSeq[newEndSeqLen-i-1] = 'G'; break;
					case 'G': contigSeq[newEndSeqLen-i-1] = 'C'; break;
					case 'T': contigSeq[newEndSeqLen-i-1] = 'A'; break;
					default: printf("line=%d, In %s(), unknown base %c, error!\n", __LINE__, __func__, contigEndSeq[i]);
							return FAILED;
				}
			}

			//strcpy(contigSeq+newEndSeqLen, contigSeq+oldEndSeqLen);
			for(i=newEndSeqLen, j=oldEndSeqLen; i<newContigLen; i++, j++) contigSeq[i] = contigSeq[j];
			contigSeq[newContigLen] = '\0';

			pContigInfo->contigLen = newContigLen;

		}else
		{ // allocate new sequence memory, and update the sequence
			newContigLen = contigLen - oldEndSeqLen + newEndSeqLen;
			newContigSeq = (char *) malloc((newContigLen+1) * sizeof(char));
			if(newContigSeq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			for(i=0; i<newEndSeqLen; i++)
			{
				switch(contigEndSeq[i])
				{
					case 'A': newContigSeq[newEndSeqLen-i-1] = 'T'; break;
					case 'C': newContigSeq[newEndSeqLen-i-1] = 'G'; break;
					case 'G': newContigSeq[newEndSeqLen-i-1] = 'C'; break;
					case 'T': newContigSeq[newEndSeqLen-i-1] = 'A'; break;
					default: printf("line=%d, In %s(), unknown base %c, error!\n", __LINE__, __func__, contigEndSeq[i]);
							return FAILED;
				}
			}
			newContigSeq[newEndSeqLen] = '\0';

			strcpy(newContigSeq+newEndSeqLen, contigSeq+oldEndSeqLen);

			free(pContigInfo->contigSeq);
			pContigInfo->contigSeq = newContigSeq;
			pContigInfo->contigLen = newContigLen;
		}
	}else
	{ // minus orientation
		if(newEndSeqLen==oldEndSeqLen)
		{ // only update the end sequence
			strcpy(contigSeq+contigLen-oldEndSeqLen, contigEndSeq);

		}else if(newEndSeqLen<oldEndSeqLen)
		{ // update the end sequence and remained sequence
			strcpy(contigSeq+contigLen-oldEndSeqLen, contigEndSeq);

			pContigInfo->contigLen -= oldEndSeqLen - newEndSeqLen;

		}else
		{ // allocate new sequence memory, and update the sequence
			newContigLen = contigLen - oldEndSeqLen + newEndSeqLen;
			newContigSeq = (char *) malloc((newContigLen+1) * sizeof(char));
			if(newContigSeq==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			strcpy(newContigSeq, contigSeq);
			strcpy(newContigSeq+contigLen-oldEndSeqLen, contigEndSeq);

			free(pContigInfo->contigSeq);
			pContigInfo->contigSeq = newContigSeq;
			pContigInfo->contigLen = newContigLen;
		}
	}

	// ############################ Debug information ######################
	if(pContigInfo->contigLen!=strlen(pContigInfo->contigSeq))
	{
		printf("line=%d, In %s(), contigLen=%d != strlen(contigSeq)=%d, error!\n", __LINE__, __func__, pContigInfo->contigLen, (int)strlen(pContigInfo->contigSeq));
		return FAILED;
	}
	// ############################ Debug information ######################

	return SUCCESSFUL;
}

/**
 * Release the scafContig nodes in local assembly.
 */
void releaseScafContigInScaf(scafContig *contighead)
{
	scafContig *contig = contighead;
	while(contig)
	{
		contig = contighead->next;

		if(contighead->ridposnum>0)
			free(contighead->pridposorientation);
		free(contighead);

		contighead = contig;
	}
	//contighead = NULL;
}
