/*
 * scafGap.c
 *
 *  Created on: Dec 19, 2013
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Construct graphs and fill the gaps between contigs in ccaffolds.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short gapFilling(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph, readSet_t *readSet)
{
	printf("Begin filling gaps between contigs, please wait ...\n");

	if(constructScafGraphFromReadset(&deBruijnGraph, contigGraph, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot build graph for gap filling, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill gaps
	if(fillGaps(scaffoldSet, contigGraph, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill gaps, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// split broken scaffolds
	if(splitScaffolds(scaffoldSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate contig overlap information, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//clean k-mer hash table
	if(cleanKmerInfoInGraph(&deBruijnGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot clean k-mer hash table, error!\n", __LINE__, __func__);
		return ERROR;
	}

	printf("End filled gaps between contigs.\n");

	return SUCCESSFUL;
}


/**
 * Build de Bruijn graph in scaffolding.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructScafGraphFromReadset(graphtype **deBruijnGraph, contigGraph_t *contigGraph, readSet_t *readSet)
{
	double kmerOccTmp;
	int32_t endKmerNumTmp, kmerNumTmp;

	// initialize the k-mer hash table
	if(initgraph(deBruijnGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the graph, error!\n", __LINE__, __func__);
		return FAILED;
	}
	(*deBruijnGraph)->readSet = readSet;

	// count the k-mer occurrences
	if(countKmerOccsFromReadsetInScaf(*deBruijnGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot count the kmers from read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// add the k-mer ridpos information
	if(addKmerRidposFromReadsetInScaf(*deBruijnGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot count the kmers from read set, error!\n", __LINE__, __func__);
		return FAILED;
	}


#if (DEBUG_PARA_PRINT==YES)
	kmerOccTmp = (double)(*deBruijnGraph)->totalItemNumRidpos / (*deBruijnGraph)->totalItemNumKmer;
	kmerNumTmp = readLen - kmerSize + 1;
	endKmerNumTmp = ceil(readLen * kmerRegLenRatioEnd5) + ceil(readLen * kmerRegLenRatioEnd3);
	//endKmerNumTmp = (ceil(averReadLenInFileSample * kmerRegLenRatioEnd5)/deBruijnGraph->kmerSampleInterval + 1) + (ceil(averReadLenInFileSample * kmerRegLenRatioEnd3)/deBruijnGraph->kmerSampleInterval + 1);
	kmerOccTmp *= (double)kmerNumTmp / endKmerNumTmp;
	//kmerOccTmp *= (double)readLen / averReadLenInFileSample;
	//kmerOccTmp /= 1.4;
	kmerOccTmp /= 1.6;
	printf("kmerNum=%ld, ridposNum=%ld, averKmerOcc=%.2f, adjustKmerOcc=%.2f\n", (*deBruijnGraph)->totalItemNumKmer, (*deBruijnGraph)->totalItemNumRidpos, (double)(*deBruijnGraph)->totalItemNumRidpos/(*deBruijnGraph)->totalItemNumKmer, kmerOccTmp);
#endif

	return SUCCESSFUL;
}

/**
 * Count the k-mer occurrences from read set in scaffolding.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short countKmerOccsFromReadsetInScaf(graphtype *deBruijnGraph)
{
	int32_t i, j;
	readSet_t *readSet;
	readBlock_t *pReadBlockArray;
	read_t *pRead, *pReadArray;
	int64_t contigID1, contigID2, contigEnd1, contigEnd2, readsCount;
	readMatchInfoBlock_t *readMatchInfoBlockArray;
	readMatchInfo_t *readMatchInfoArray, *pReadMatchInfo;

	struct timeval tpstart, tpend;
	double timeused_counting;
	gettimeofday(&tpstart, NULL);


	printf("Filling k-mer information ...\n");

	readsCount = 0;
	readSet = deBruijnGraph->readSet;

	pKmerBlockTmp = deBruijnGraph->kmerBlockArr + deBruijnGraph->blocksNumKmer - 1;
	pKmerTmpDoing = pKmerBlockTmp->kmerArr;

	pKmerseqBlockTmp = deBruijnGraph->kmerSeqBlockArr + deBruijnGraph->blocksNumKmerSeq - 1;
	pKmerSeqTmpDoing = pKmerseqBlockTmp->kmerSeqArr;

	// get the mono reads from read set
	for(i=0; i<readSet->blocksNumRead; i++)
	{
		pReadBlockArray = readSet->readBlockArr + i;
		pReadArray = pReadBlockArray->readArr;
		readMatchInfoBlockArray = readSet->readMatchInfoBlockArr + i;
		readMatchInfoArray = readMatchInfoBlockArray->readMatchInfoArr;

		if(pReadBlockArray->itemNum%2!=0)
		{
			printf("line=%d, In %s(), itemNum=%d\n", __LINE__, __func__, pReadBlockArray->itemNum);
			return FAILED;
		}

		for(j=0; j<pReadBlockArray->itemNum; j+=2)
		{
			pRead = pReadArray + j;
			pReadMatchInfo = readMatchInfoArray + j;

			if(pRead->validFlag==YES && (pRead+1)->validFlag==YES)
			{
				contigID1 = pReadMatchInfo->contigID;
				contigEnd1 = pReadMatchInfo->contigEnd;
				contigID2 = (pReadMatchInfo+1)->contigID;
				contigEnd2 = (pReadMatchInfo+1)->contigEnd;

				//count the kmers
				if(contigID1==0 && ((contigID2>0 && contigEnd2!=-1) || contigID2==0))
				{
					if(addReadPreFromReadset(pRead, deBruijnGraph)==FAILED)
					{
						printf("line=%d, In %s(), cannot count the kmer, error!\n", __LINE__, __func__);
						return FAILED;
					}

					readsCount ++;
				}else if(((contigID1>0 && contigEnd1!=-1) || contigID1==0) && contigID2==0)
				{
					if(addReadPreFromReadset(pRead+1, deBruijnGraph)==FAILED)
					{
						printf("line=%d, In %s(), cannot count the kmer, error!\n", __LINE__, __func__);
						return FAILED;
					}

					readsCount ++;
				}
			}

			//fprintf(fpOut, "rid=%ld, seqlen=%d, nBaseNum=%d, validFlag=%d, successFlag=%d, seq=%s\n", rid, (int32_t)pRead->seqlen, (int32_t)pRead->nBaseNum, (int32_t)pRead->validFlag, (int32_t)pRead->successFlag, getReadBaseByInt(readseqInt, pRead->seqlen));
		}
	}

	// the last k-mer block is empty, remove it
	if(pKmerBlockTmp->itemNum==0)
	{
		free(pKmerBlockTmp->kmerArr);
		deBruijnGraph->kmerBlockArr[deBruijnGraph->blocksNumKmer-1].kmerArr = NULL;
		deBruijnGraph->blocksNumKmer --;
	}

	// the last kmerseq block is empty, remove it
	if(pKmerseqBlockTmp->itemNum==0)
	{
		free(pKmerseqBlockTmp->kmerSeqArr);
		deBruijnGraph->kmerSeqBlockArr[deBruijnGraph->blocksNumKmerSeq-1].kmerSeqArr = NULL;
		deBruijnGraph->blocksNumKmerSeq --;
	}


	gettimeofday(&tpend, NULL);
	timeused_counting = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Filling k-mers used time: %.2f seconds, readsCount=%ld.\n", timeused_counting, readsCount);

	return SUCCESSFUL;
}

/**
 * Add the k-mer ridpos information from read set.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addKmerRidposFromReadsetInScaf(graphtype *graph)
{
	int32_t i, j;
	readSet_t *readSet;
	readBlock_t *pReadBlockArray;
	read_t *pRead, *pReadArray;
	int64_t readID1, readID2, contigID1, contigID2, contigEnd1, contigEnd2, readsCount;
	readMatchInfoBlock_t *readMatchInfoBlockArray;
	readMatchInfo_t *readMatchInfoArray, *pReadMatchInfo;

	struct timeval tpstart, tpend;
	double timeused_adding;
	gettimeofday(&tpstart, NULL);


	printf("Filling reads information ...\n");

	// initialize ridpos blocks and the ridpos regions in that blocks
	if(initRidposBlocksInGraph(graph)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the mono reads from read set
	readID1 = 0;
	readsCount = 0;
	readSet = graph->readSet;
	for(i=0; i<readSet->blocksNumRead; i++)
	{
		pReadBlockArray = readSet->readBlockArr + i;
		pReadArray = pReadBlockArray->readArr;
		readMatchInfoBlockArray = readSet->readMatchInfoBlockArr + i;
		readMatchInfoArray = readMatchInfoBlockArray->readMatchInfoArr;

		if(pReadBlockArray->itemNum%2!=0)
		{
			printf("line=%d, In %s(), itemNum=%d\n", __LINE__, __func__, pReadBlockArray->itemNum);
			return FAILED;
		}

		for(j=0; j<pReadBlockArray->itemNum; j+=2)
		{
			readID1 ++;
			readID2 = readID1 + 1;
			pRead = pReadArray + j;
			pReadMatchInfo = readMatchInfoArray + j;

			if(pRead->validFlag==YES && (pRead+1)->validFlag==YES)
			{
				contigID1 = pReadMatchInfo->contigID;
				contigEnd1 = pReadMatchInfo->contigEnd;
				contigID2 = (pReadMatchInfo+1)->contigID;
				contigEnd2 = (pReadMatchInfo+1)->contigEnd;

				//add the kmers
				if(contigID1==0 && ((contigID2>0 && contigEnd2!=-1) || contigID2==0))
				{
					if(addReadFromReadset(readID1, pRead, graph)==FAILED)
					{
						printf("line=%d, In %s(), cannot count the kmer, error!\n", __LINE__, __func__);
						return FAILED;
					}
					readsCount ++;
				}else if(((contigID1>0 && contigEnd1!=-1) || contigID1==0) && contigID2==0)
				{
					if(addReadFromReadset(readID2, pRead+1, graph)==FAILED)
					{
						printf("line=%d, In %s(), cannot count the kmer, error!\n", __LINE__, __func__);
						return FAILED;
					}
					readsCount ++;
				}
			}
			//fprintf(fpOut, "rid=%ld, seqlen=%d, nBaseNum=%d, validFlag=%d, successFlag=%d, seq=%s\n", rid, (int32_t)pRead->seqlen, (int32_t)pRead->nBaseNum, (int32_t)pRead->validFlag, (int32_t)pRead->successFlag, getReadBaseByInt(readseqInt, pRead->seqlen));

			readID1 = readID2;
		}
	}

	gettimeofday(&tpend, NULL);
	timeused_adding = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Filling reads information used time: %.2f seconds, readsCount=%ld.\n", timeused_adding, readsCount);

	return SUCCESSFUL;
}

/**
 * Fill the gaps between contigs in ccaffolds.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillGaps(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph, readSet_t *readSet)
{
	if(initGapFillingParas()==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the parameters for gap filling, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// local assembly to fill gaps
	if(localAssemblyInScaf(scaffoldSet, contigGraph, readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill gaps, error!\n", __LINE__, __func__);
		return FAILED;
	}

	freeGapFillingParas();

	return SUCCESSFUL;
}


/**
 * Initialize the global parameters for gap filling.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initGapFillingParas()
{
	int32_t i;

	initFirstKmerThreshold();

	if(initMemory()==FAILED)
	{
		printf("line=%d, In %s(), cannot init the memory for filling gaps, error!\n", __LINE__, __func__);
		return FAILED;
	}

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

	return SUCCESSFUL;
}

/**
 * Free the global parameters for gap filling.
 */
void freeGapFillingParas()
{
	int32_t i;

	freeMemory();

	// free memory of *contigEndSeqInScaf[2]
	for(i=0; i<2; i++)
	{
		free(scafContigEndSeqArr[i]);
		scafContigEndSeqArr[i] = NULL;
	}

	// free memory of comparisonSeqInScaf
	free(comparisonSeqInScaf);
	comparisonSeqInScaf = NULL;

	free(scoreArr);
	scoreArr = NULL;

	for(i=0; i<3; i++)
	{
		free(alignResultArr[i]);
		alignResultArr[i] = NULL;
	}
}

/**
 * Fill gaps between contigs in scaffolds by local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short localAssemblyInScaf(scaffoldSet_t *scaffoldSet, contigGraph_t *contigGraph, readSet_t *readSet)
{
	scaffoldItem_t *scaffoldItem;
	contigOverlap_t *pContigOverlapInfo, *contigOverlapArray;
	int32_t scaffoldID, contigID[2], contigOrient[2], contigLen[2], gapSize, newGapSize;
	int32_t i, j, k, startRowInCOI, itemNumInCOA, linkedContigsNum;
	int32_t assemblyCycle, maxAssemblyLen, minAssemblyLen;	// assemblyRound: 0 -- first round; 1 -- second round
	contigtype *localScafContigheadArr[2];
	char *localScafContigSeqArr[2], tmpSeq1[MAX_READ_LEN_IN_BUF+1], tmpSeq2[MAX_READ_LEN_IN_BUF+1];
	int32_t localScafContigNodesNumArr[2];
	int32_t overlapLen, oldEndSeqLenArray[2];
	int32_t breakFlag, successFilledFlag, tmp_gapSize, validSuccessReadFlag;

	localScafContigSeqArr[0] = tmpSeq1;
	localScafContigSeqArr[1] = tmpSeq2;

	if(initPEHashParas()==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the PEHash table parameters, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<11; i++) errNumArr[i] = 0;


	// fill gaps between contigs in scaffolds by local assembly
	scaffoldItem = scaffoldSet->scaffoldItemList;
	localContigID = 0;
	while(scaffoldItem)
	{
		// ####################### Debug information ####################
#if(DEBUG_SCAF_FILLGAP_FLAG==YES)
		//printf("================== scaffoldID=%d ==============\n", scaffoldItem->scaffoldID);
		if(scaffoldItem->scaffoldID==101)
		{
			printf("scaffoldID=%d\n", scaffoldItem->scaffoldID);
		}
#endif
		// ####################### Debug information ####################

		scaffoldID = scaffoldItem->scaffoldID;
		linkedContigsNum = scaffoldItem->linkedContigsNum;

		if(linkedContigsNum>=2)
		{
			contigOverlapArray = scaffoldItem->contigOverlapArray;
			itemNumInCOA = scaffoldItem->itemNumContigOverlapArray;
			for(j=0; j<itemNumInCOA; j++)
			{
				pContigOverlapInfo = contigOverlapArray + j;
				if(pContigOverlapInfo->mergeFlag==NO)
				{ // a gap between the two contigs

					// get the information of the two contigs
					contigID[0] = pContigOverlapInfo->contigID1;
					contigID[1] = pContigOverlapInfo->contigID2;
					contigOrient[0] = pContigOverlapInfo->orientation1;
					contigOrient[1] = pContigOverlapInfo->orientation2;
					contigLen[0] = contigGraph->contigItemArray[contigID[0]-1].contigLen;
					contigLen[1] = contigGraph->contigItemArray[contigID[1]-1].contigLen;
					gapSize = pContigOverlapInfo->gapSize;

					//localScafContigheadArr[0] = localScafContigheadArr[1] = NULL;
					localScafContigNodesNumArr[0] = localScafContigNodesNumArr[1] = 0;

					//###################### Debug information #####################
#if(DEBUG_SCAF_FILLGAP_FLAG==YES)
					if(contigID[0]==4053 && contigID[1]==2842)
					{
						printf("line=%d, In %s(), contigID1=%d, contigOrient1=%d, contigLen1=%d, contigID2=%d, contigOrient2=%d, contigLen2=%d, gapSize=%d\n", __LINE__, __func__, contigID[0], contigOrient[0], contigLen[0], contigID[1], contigOrient[1], contigLen[1], gapSize);
					}
					//printf("line=%d, In %s(), contigID1=%d, contigOrient1=%d, contigLen1=%d, contigID2=%d, contigOrient2=%d, contigLen2=%d, gapSize=%d\n", __LINE__, __func__, contigID[0], contigOrient[0], contigLen[0], contigID[1], contigOrient[1], contigLen[1], gapSize);
#endif
					//###################### Debug information #####################

					// get the two contig ends and the ends length
					if(getScafContigEndSeqs(scafContigEndSeqArr, scafContigEndSeqLenArr, pContigOverlapInfo, contigGraph)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the scafContig end sequences of contigs [ %d, %d ], error!\n", __LINE__, __func__, contigID[0], contigID[1]);
						return FAILED;
					}

					// assembly two contig ends
					for(assemblyCycle=0; assemblyCycle<2; assemblyCycle++)
					{
						// initialize the variables
						successFilledFlag = NO;
						breakFlag = NO;

						assemblyRound = SECOND_ROUND_ASSEMBLY;
						itemNumContigArr = 0;
						itemNumDecisionTable = 0;
						successContigIndex = -1;
						lockedReadsNum = 0;
						readsNumInPEHashArr = 0;
						regLenPEHash = 0;
						localContigID ++;
						readsNumRatio = 1;
						naviTandFlag = NAVI_UNUSED;

						// record the ends length before local assembly
						oldEndSeqLenArray[0] = scafContigEndSeqLenArr[0];
						oldEndSeqLenArray[1] = scafContigEndSeqLenArr[1];

						// get the comparison sequence
						if(getComparisonSeqInScaf(comparisonSeqInScaf, &comparisonSeqLenInScaf, scafContigEndSeqArr, scafContigEndSeqLenArr, assemblyCycle)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the comparison sequence, error!\n", __LINE__, __func__);
							return FAILED;
						}

						// prepare the first scafContig nodes with length of READ_LEN before local assembly
						if(prepareAssemblyInScaf(pContigOverlapInfo, contigGraph, scafContigEndSeqArr, scafContigEndSeqLenArr, assemblyCycle)==FAILED)
						{
							printf("line=%d, In %s(), cannot prepare the assembly of contigs [ %d, %d ], error!\n", __LINE__, __func__, contigID[0], contigID[1]);
							return FAILED;
						}

						if(getMaxMinAssemblyLen(&maxAssemblyLen, &minAssemblyLen, assemblyCycle, localScafContigNodesNumArr, gapSize)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the maximal local assembly length, error!\n", __LINE__, __func__);
							return FAILED;
						}

						while(naviSuccessFlag==NAVI_SUCCESS)
						{

							// initialize or update the PE hash table
							if(PEGivenType>NONE_PE_GIVEN_TYPE)
							{
								if(updatePEHashTable(itemNumContigArr, assemblyRound)==FAILED)
								{
									printf("line=%d, In %s(), localContigID=%ld, contigID=%d, itemNumContigArr=%ld, cannot update the PE hash table, error!\n", __LINE__, __func__, localContigID, contigsNum+1, itemNumContigArr);
									return FAILED;
								}

								//if(readsNumInPEHashArr>=minReadsNumPEHashThres && regLenPEHash>=minRegLenUsingPE)
								if(readsNumInPEHashArr>0)
								{
									if(getNextKmerByMix(itemNumContigArr, assemblyRound)==FAILED)
									{
										printf("line=%d, In %s(), localContigID=%ld, cannot get the next kmer by mix, error!\n", __LINE__, __func__, localContigID);
										return FAILED;
									}

								}else
								{
									navigationFlag = NAVI_SE_FLAG;
									naviTandFlag = NAVI_UNUSED;
									if(getNextKmerBySE(itemNumContigArr)==FAILED)
									{
										printf("line=%d, In %s(), localContigID=%ld, cannot get next kmer, error!\n", __LINE__, __func__, localContigID);
										return FAILED;
									}

									for(i=0; i<4; i++) { occsNumPE[i] = 0; occsNumIndexPE[i] = -1; }

#if (SVM_NAVI==YES)
									//if((successContigIndex>0 && itemNumContigArr-successContigIndex>50) || readsNumRatio<0.3*minReadsNumRatioThres)
									//if((successContigIndex>0 && itemNumContigArr-successContigIndex>30) || readsNumRatio<0.3*minReadsNumRatioThres) // 2013-11-12
									if(contigPath->naviSuccessSize<2*readLen && ((successContigIndex>0 && itemNumContigArr-successContigIndex>30) || readsNumRatio<0.3*minReadsNumRatioThres)) // 2013-11-12
										naviSuccessFlag = NAVI_FAILED;
									else if(secondOccSE>0)								// deleted 2013-02-26
									//if(naviSuccessFlag==NAVI_SUCCESS && secondOccSE>0 && successContigIndex>0)		// added 2013-02-26
									{
				//							sumSecondOccSE = 0;
				//							for(i=0; i<4; i++) if(i!=occsNumIndexSE[0]) sumSecondOccSE += occsNumSE[i];

										svmFeatureArr[0] = occsNumSE[occsNumIndexSE[0]];
										svmFeatureArr[1] = occsNumSE[occsNumIndexSE[1]];
										//svmFeatureArr[1] = sumSecondOccSE;
										svmFeatureArr[2] = readsNumRatio;
										//svmFeatureArr[2] = svmFeatureArr[0] / svmFeatureArr[1];
										if(successContigIndex>0)
										{
											// compute the maximal gap size in contig tail region
											if(computeGapSizeInContig(&tmp_gapSize, contigArr, itemNumContigArr, assemblyRound)==FAILED)
											{
												printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
												return FAILED;
											}

											//svmFeatureArr[3] = itemNumContigArr - successContigIndex;
											svmFeatureArr[3] = tmp_gapSize;
										}else
											svmFeatureArr[3] = 0;

										if(fillSampleDataSVM(svmSampleSE, svmFeatureArr)==FAILED)
										{
											printf("line=%d, In %s(), cannot fill the sample data for SVM model, error!\n", __LINE__, __func__);
											return FAILED;
										}

										if(extensionDecisionBySvm(&naviSuccessFlag, svmSampleSE, svmModelSE)==FAILED)
										{
											printf("line=%d, In %s(), cannot decide the extension by SVM calssifier, error!\n", __LINE__, __func__);
											return FAILED;
										}


										//if(successContigIndex==-1 && (secondOccSE>=0.7*maxOccSE || maxOccSE-secondOccSE<2))		// added 2013-04-01
										//if(maxOccSE-secondOccSE<2 && (successContigIndex>0 &&  itemNumContigArr-successContigIndex>10))
										if(maxOccSE-secondOccSE<2)
											naviSuccessFlag = NAVI_FAILED;

										if(svmFeatureArr[1]>averKmerOcc || (svmFeatureArr[1]>=0.5*averKmerOcc && svmFeatureArr[1]/svmFeatureArr[0]>=0.6)) //2013-11-12
											naviSuccessFlag = NAVI_FAILED;

				//						if(naviSuccessFlag==NO && svmFeatureArr[1]/svmFeatureArr[0]<0.5)
				//							naviSuccessFlag = NAVI_SUCCESS;

				//						if(svmFeatureArr[0]==svmFeatureArr[1] || (svmFeatureArr[1]>=2 && svmFeatureArr[1]/svmFeatureArr[0]>=0.7))
				//							naviSuccessFlag = NAVI_FAILED;

										//if(localContigID==13 /*&& itemNumContigArr>=14805*/ && assemblyRound==FIRST_ROUND_ASSEMBLY)
										//if(naviSuccessFlag==NAVI_FAILED) //2013-11-12
										if(naviSuccessFlag==NAVI_FAILED && itemNumContigArr>readLen) //2013-11-12
										{
											if(decideByCandPathSE(&naviSuccessFlag, &maxBaseIndexAfterCandPathSE, &incorrectBaseNumCandPathSE, occsNumSE, occsNumIndexSE, NULL, NULL, decisionTable, itemNumDecisionTable, contigPath)==FAILED)
											{
												printf("line=%d, In %s(), cannot navigate by candidate paths using single ends, error!\n", __LINE__, __func__);
												return FAILED;
											}

											if(naviSuccessFlag==NAVI_FAILED && incorrectBaseNumCandPathSE>MAX_INCOR_BASE_NUM_CANDPATH)
											{
												// check the tandem repeats
												if(decideByCheckingTandemRepeatSE(&naviSuccessFlag, &maxBaseIndexAfterTandPathSE, &incorrectBaseNumTandPathSE, &newCandBaseNumAfterTandPathSE, occsNumSE, occsNumIndexSE, itemNumContigArr-1, decisionTable, &itemNumDecisionTable, contigPath, deBruijnGraph)==FAILED)
												{
													printf("line=%d, In %s(), cannot check the tandem repeats, error!\n", __LINE__, __func__);
													return FAILED;
												}
												naviTandFlag = naviSuccessFlag;
											}
										}
									}
#else
									//if(maxOccSE==secondOccSE)
									if((maxOccSE==secondOccSE) || (secondOccSE>=2 && secondOccSE>=SEC_MAX_OCC_RATIO_SVM*maxOccSE))
									{
										naviSuccessFlag = NAVI_FAILED;
										//kmers[0] = kmers[1] = NULL;
									}
#endif
								}
							}else
							{
								navigationFlag = NAVI_SE_FLAG;
								naviTandFlag = NAVI_UNUSED;
								if(getNextKmerBySE(itemNumContigArr)==FAILED)
								{
									printf("line=%d, In %s(), localContigID=%ld, cannot get next kmer, error!\n", __LINE__, __func__, localContigID);
									return FAILED;
								}
								for(i=0; i<4; i++) { occsNumPE[i] = 0; occsNumIndexPE[i] = -1; }


#if (SVM_NAVI==YES)
								//if((successContigIndex>0 && itemNumContigArr-successContigIndex>50) || readsNumRatio<0.3*minReadsNumRatioThres)
								//if((successContigIndex>0 && itemNumContigArr-successContigIndex>30) || readsNumRatio<0.3*minReadsNumRatioThres)
								if(contigPath->naviSuccessSize<2*readLen && ((successContigIndex>0 && itemNumContigArr-successContigIndex>30) || readsNumRatio<0.3*minReadsNumRatioThres)) // 2014-01-31
									naviSuccessFlag = NAVI_FAILED;
								else if(secondOccSE>0)								// deleted 2013-02-26
								//if(naviSuccessFlag==NAVI_SUCCESS && secondOccSE>0 && successContigIndex>0)		// added 2013-02-26
								{
				//					sumSecondOccSE = 0;
				//					for(i=0; i<4; i++) if(i!=occsNumIndexSE[0]) sumSecondOccSE += occsNumSE[i];

									svmFeatureArr[0] = occsNumSE[occsNumIndexSE[0]];
									svmFeatureArr[1] = occsNumSE[occsNumIndexSE[1]];
									//svmFeatureArr[1] = sumSecondOccSE;
									svmFeatureArr[2] = readsNumRatio;
									//svmFeatureArr[2] = svmFeatureArr[0] / svmFeatureArr[1];
									if(successContigIndex>0)
									{
										// compute the maximal gap size in contig tail region
				//						if(computeGapSizeInContig(&tmp_gapSize, contigArr, itemNumContigArr, assemblyRound)==FAILED)
				//						{
				//							printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
				//							return FAILED;
				//						}

										svmFeatureArr[3] = itemNumContigArr - successContigIndex;
										//svmFeatureArr[3] = tmp_gapSize;
									}else
										svmFeatureArr[3] = 0;

									if(fillSampleDataSVM(svmSampleSE, svmFeatureArr)==FAILED)
									{
										printf("line=%d, In %s(), cannot fill the sample data for SVM model, error!\n", __LINE__, __func__);
										return FAILED;
									}

									if(extensionDecisionBySvm(&naviSuccessFlag, svmSampleSE, svmModelSE)==FAILED)
									{
										printf("line=%d, In %s(), cannot decide the extension by SVM calssifier, error!\n", __LINE__, __func__);
										return FAILED;
									}


				//					if(naviSuccessFlag==NAVI_FAILED)
				//					{
				//						decideByMalignSE(&naviSuccessFlag, decisionTable, itemNumDecisionTable);
				//					}

									//if(successContigIndex==-1 && (secondOccSE>=0.7*maxOccSE || maxOccSE-secondOccSE<2))		// added 2013-04-01
									//if(maxOccSE-secondOccSE<=2 && (successContigIndex>0 &&  itemNumContigArr-successContigIndex>10))
									if(maxOccSE-secondOccSE<2)
										naviSuccessFlag = NAVI_FAILED;

									if(svmFeatureArr[1]>averKmerOcc || (svmFeatureArr[1]>=0.5*averKmerOcc && svmFeatureArr[1]/svmFeatureArr[0]>=0.6))  //2013-11-12
										naviSuccessFlag = NAVI_FAILED;

				//					if(naviSuccessFlag==NAVI_FAILED && svmFeatureArr[1]/svmFeatureArr[0]<0.5)
				//						naviSuccessFlag = NAVI_SUCCESS;

				//					if(svmFeatureArr[0]==svmFeatureArr[1] || (svmFeatureArr[1]>=2 && svmFeatureArr[1]/svmFeatureArr[0]>=0.7))
				//						naviSuccessFlag = NAVI_FAILED;

									if(naviSuccessFlag==NAVI_FAILED) //2013-11-12
									{
										if(decideByCandPathSE(&naviSuccessFlag, &maxBaseIndexAfterCandPathSE, &incorrectBaseNumCandPathSE, occsNumSE, occsNumIndexSE, NULL, NULL, decisionTable, itemNumDecisionTable, contigPath)==FAILED)
										{
											printf("line=%d, In %s(), cannot navigate by candidate paths using single ends, error!\n", __LINE__, __func__);
											return FAILED;
										}

										if(naviSuccessFlag==NAVI_FAILED && incorrectBaseNumCandPathSE>MAX_INCOR_BASE_NUM_CANDPATH)
										{
											// check the tandem repeats
											if(decideByCheckingTandemRepeatSE(&naviSuccessFlag, &maxBaseIndexAfterTandPathSE, &incorrectBaseNumTandPathSE, &newCandBaseNumAfterTandPathSE, occsNumSE, occsNumIndexSE, itemNumContigArr-1, decisionTable, &itemNumDecisionTable, contigPath, deBruijnGraph)==FAILED)
											{
												printf("line=%d, In %s(), cannot check the tandem repeats, error!\n", __LINE__, __func__);
												return FAILED;
											}
											naviTandFlag = naviSuccessFlag;
										}
									}
								}
#else
								//if(maxOccSE==secondOccSE)
								if((maxOccSE==secondOccSE) || (secondOccSE>=2 && secondOccSE>=SEC_MAX_OCC_RATIO_SVM*maxOccSE))
								{
									naviSuccessFlag = NAVI_FAILED;
									//kmers[0] = kmers[1] = NULL;
								}
#endif
							}

							// check the reads number in the sub read region
							if(naviSuccessFlag==NAVI_SUCCESS)
							{
								if(itemNumContigArr>=minContigLenCheckingReadsNum)
								{
									if(updateReadsNumReg(itemNumSuccessReadsArr, itemNumContigArr, assemblyRound)==FAILED)
									{
										printf("line=%d, In %s(), localContigID=%ld, itemNumContigArr=%ld, cannnot check the reads number in reads number region, error!\n", __LINE__, __func__, localContigID, itemNumContigArr);
										return FAILED;
									}
								}

#if(SVM_NAVI==YES)
								// added 2013-10-14
								// confirm the navigation of one single-end read
								if(confirmSingleEndNavi(&naviSuccessFlag, contigArr, itemNumContigArr)==FAILED)
								{
									printf("line=%d, In %s(), localContigID=%ld, itemNumContigArr=%ld, cannnot check the navigation of single-end read, error!\n", __LINE__, __func__, localContigID, itemNumContigArr);
									return FAILED;
								}
#endif
							}

							//if(kmers[0]==NULL && kmers[1]==NULL)
							if(naviSuccessFlag==NAVI_FAILED)
							{
								break;
							}

							itemNumContigArr ++;

							// Append a base to contig tail
							if(addContigBase(kmerSeqIntAssembly[entriesPerKmer-1] & 3)==FAILED)
							{
								printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, cannot add a contig base, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
								return FAILED;
							}

							// added 2013-10-14
							// add the read number
							if(addOccsToContig(contigArr, itemNumContigArr, naviTandFlag, newCandBaseNumAfterTandPathPE, contigPath)==FAILED)
							{
								printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, cannot add occ number to contig, error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
								return FAILED;
							}

							// update the decision table according to kmers
							if(updateDecisionTable(kmers, kmerSeqIntAssembly[entriesPerKmer-1] & 3)==FAILED)
							{
								printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, cannot update decision table, error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
								return FAILED;
							}

							// Update the reads status in decision table
							if(updateAssemblingreadsStatus()==FAILED)
							{
								printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, cannot update reads status in decision table, error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
								return FAILED;
							}

							// Update the locked reads and their total number
							if(updateLockedReads()==FAILED)
							{
								printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, cannot update locked reads, error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
								return FAILED;
							}

							// update the finished reads in decision table, and record the successful reads into successful reads array
							if(updateFinishedReadsInDecisionTable()==FAILED)
							{
								printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, cannot update finished reads in decision table, error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
								return FAILED;
							}

							// update the contig path
							if(updateContigPath(contigPath, navigationFlag, kmers, decisionTable, &itemNumDecisionTable, dtRowHashtable, contigArr, itemNumContigArr)==FAILED)
							{
								printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, cannot update contig path, error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
								return FAILED;
							}

							if(itemNumSuccessReadsArr>0)
							{

								// delete reads from De Bruijn graph
								if(delReadsFromGraph(successReadsArr, itemNumSuccessReadsArr)==FAILED)
								{
									printf("line=%d, In %s(), localContigID=%ld, contigID=%d, itemNumContigArr=%ld, assemblyRound=%d, cannot delete the reads from graph, error!\n", __LINE__, __func__, localContigID, contigsNum+1, itemNumContigArr, assemblyRound);
									outputSuccessReads(successReadsArr, itemNumSuccessReadsArr);
									return FAILED;
								}

								// add the successful reads information to contig chain
								if(addRidposToContig(successReadsArr, itemNumSuccessReadsArr, itemNumContigArr)==FAILED)
								{
									printf("line=%d, In %s(), localContigID=%ld, contigID=%d, cannot add a contig base, error!\n", __LINE__, __func__, localContigID, contigsNum+1);
									return FAILED;
								}

								// add the success reads to PEHashtable, 2014-01-26
								if(PEGivenType>NONE_PE_GIVEN_TYPE && itemNumContigArr>=minContigLenUsingPE)
								{
									if(addSuccessReadsToPEHashtable(successReadsArr, itemNumSuccessReadsArr, assemblyRound)==FAILED)
									{
										printf("line=%d, In %s(), localContigID=%ld, contigID=%d, cannot add success reads to PEHashtable, error!\n", __LINE__, __func__, localContigID, contigsNum+1);
										return FAILED;
									}
								}

								// get the valid success read flag
								if(getValidSuccessReadFlag(&validSuccessReadFlag, successReadsArr, itemNumSuccessReadsArr, minMatchNumSuccessRead)==FAILED)
								{
									printf("line=%d, In %s(), localContigID=%ld, contigID=%d, cannot get the success read flag, error!\n", __LINE__, __func__, localContigID, contigsNum+1);
									return FAILED;
								}
								if(validSuccessReadFlag==YES)
									successContigIndex = itemNumContigArr;

								this_successReadNum += itemNumSuccessReadsArr;
							}

							// process the dead loops
							if(successContigIndex<=0)
							{
								if(itemNumContigArr>2*readLen-kmerSize)
								{
									//printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, successContigIndex<=0!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
									break;
								}
							}else if(itemNumContigArr-successContigIndex > readLen-MIN_OVERLAP_LEN)
							{
								break;
							}
							//===============================================================
							else if(itemNumContigArr>maxAssemblyLen)
							{ // get the maximal assembly length and terminate
								break;
							}
							//===============================================================

						} // while(naviSuccessFlag==NAVI_SUCCESS)

						if(successContigIndex>0)
						{
							if(updateContigtailnodes(contigArr, successContigIndex, &itemNumContigArr, assemblyRound)==FAILED)  // deleted 2012-12-30
							{
								printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, cannot update Contigtail nodes, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
								return FAILED;
							}

							if(PEGivenType>NONE_PE_GIVEN_TYPE && readsNumInPEHashArr>0)
							{
								if(cleanReadsFromPEHashtable()==FAILED)
								{
									printf("line=%d, In %s(), localContigID=%ld, cannot clean PE hash table, error!\n", __LINE__, __func__, localContigID);
									return FAILED;
								}
							}

							// generate the new assembled contig sequence
							if(getNewContigSeq(localScafContigSeqArr[assemblyCycle], localScafContigNodesNumArr+assemblyCycle, contigArr, itemNumContigArr)==FAILED)
							{
								printf("line=%d, In %s(), cannot get the new assembled contig sequence, error!\n", __LINE__, __func__);
								return FAILED;
							}

							//######################### Debug information #####################
#if(DEBUG_SCAF_FILLGAP_FLAG==YES)
							printf("Before updateScafContigEndSeqs(), assemblyCycle=%d, contigEndSeq[%d]=%s, len=%d\n", assemblyCycle, assemblyCycle, scafContigEndSeqArr[assemblyCycle], scafContigEndSeqLenArr[assemblyCycle]);
#endif
							//######################### Debug information #####################

							// update the contig end sequences
							if(updateScafContigEndSeqs(localScafContigSeqArr[assemblyCycle], localScafContigNodesNumArr[assemblyCycle], assemblyCycle)==FAILED)
							{
								printf("line=%d, In %s(), cannot update the scafContig end sequences, error!\n", __LINE__, __func__);
								return FAILED;
							}

							// adjust the scafContig nodes number and the oldEndSeqLenArr[0,1]
							//localScafContigNodesNumArr[assemblyCycle] = itemNumContigArr;
							oldEndSeqLenArray[assemblyCycle] = scafContigEndSeqLenArr[assemblyCycle];

							//######################### Debug information #####################
#if(DEBUG_SCAF_FILLGAP_FLAG==YES)
							printf("After updateScafContigEndSeqs(), assemblyCycle=%d, contigEndSeq[%d]=%s, len=%d\n", assemblyCycle, assemblyCycle, scafContigEndSeqArr[assemblyCycle], scafContigEndSeqLenArr[assemblyCycle]);
#endif
							//######################### Debug information #####################

							if(itemNumContigArr >= minAssemblyLen)
							//if(scafContigIndex >= readLen)  // 2012-11-19
							{ // the contig nodes number is more than the minimal assembly length
								// detect overlaps of contig ends in local assembly
								if(detectOverlapsInScaf(&successFilledFlag, &overlapLen, &newGapSize, &breakFlag, gapSize, localScafContigNodesNumArr, assemblyCycle)==FAILED)
								{
									printf("line=%d, In %s(), cannot detect overlaps between contig ends, error!\n", __LINE__, __func__);
									return FAILED;
								}
							}
						}

						localScafContigNodesNumArr[assemblyCycle] = itemNumContigArr;

						// clean the PE hash table
						if(cleanReadsFromPEHashtable()==FAILED)
						{
							printf("line=%d, In %s, cannot clean PE hash table, error!\n", __LINE__, __func__);
							return FAILED;
						}

						// clean the DTRow hash table
						if(cleanDTRowIndexHashtable(dtRowHashtable)==FAILED)
						{
							printf("line=%d, In %s(), contigsNum=%d, cannot clean DTRow hash table, Error!\n", __LINE__, __func__, contigsNum+1);
							return ERROR;
						}

						// clean the contig path
						if(cleanContigPath(contigPath)==FAILED)
						{
							printf("line=%d, In %s(), cannot clean the contig path, error!\n", __LINE__, __func__);
							return FAILED;
						}

						// clean the contig array
						cleanContigArray(contigArr, &itemNumContigArr);

						// control sentence
						if(successFilledFlag==YES)
							break;

					} // end for(assemblyRound=FIRST_ROUND_ASSEMBLY; assemblyRound<=SECOND_ROUND_ASSEMBLY; assemblyRound++)

					// check filled status, and update the information in local assembly
					if(successFilledFlag==YES)
					{ // successfully fill the gap, then update the overlap information, contig information, ...

						if(updateContigOverlapInfoInScaf(pContigOverlapInfo, successFilledFlag, overlapLen, newGapSize, breakFlag)==FAILED)
						{
							printf("line=%d, In %s(), cannot update contig overlap information, error!\n", __LINE__, __func__);
							return FAILED;
						}

						if(updateContigInfoInScaf(pContigOverlapInfo, localScafContigSeqArr, localScafContigNodesNumArr, oldEndSeqLenArray, contigGraph, assemblyCycle)==FAILED) ////////////////////////////////////////////////////////
						{
							printf("line=%d, In %s(), cannot update contig information, error!\n", __LINE__, __func__);
							return FAILED;
						}

					}else
					{ // failed fill the gap, update the overlap informatipon and the contig information
						// ######################### Debug information ########################
#if(DEBUG_SCAF_FILLGAP_FLAG==YES)
						printf("scaffoldID=%d, contigID1=%d, contigID2=%d, contigOrient1=%d, contigOrient2=%d, assemblyRound=%d, gapSize=%d\n", scaffoldID, contigID[0], contigID[1], contigOrient[0], contigOrient[1], assemblyRound, gapSize);
						printf("contigNodesNum1=%d, contigNodesNum2=%d, prepareAssemblyLen1=%d, prepareAssemblyLen2=%d\n", localScafContigNodesNumArr[0], localScafContigNodesNumArr[1], prepareAssemblyLenArr[0], prepareAssemblyLenArr[1]);
#endif
						// ######################### Debug information ########################

						// detect overlaps of contig ends in local assembly
						if(detectOverlapsInScaf(&successFilledFlag, &overlapLen, &newGapSize, &breakFlag, gapSize, localScafContigNodesNumArr, assemblyCycle)==FAILED)
						{
							printf("line=%d, In %s(), cannot detect overlaps between contig ends, error!\n", __LINE__, __func__);
							return FAILED;
						}

						// To do ...
						if(localScafContigNodesNumArr[0] > prepareAssemblyLenArr[0] || localScafContigNodesNumArr[1] > prepareAssemblyLenArr[1])
						{ // new bases are appended

							// update overlap information
							if(updateContigOverlapInfoInScaf(pContigOverlapInfo, successFilledFlag, overlapLen, newGapSize, breakFlag)==FAILED)
							{
								printf("line=%d, In %s(), cannot update contig overlap information, error!\n", __LINE__, __func__);
								return FAILED;
							}

							// update contig information
							if(updateContigInfoInScaf(pContigOverlapInfo, localScafContigSeqArr, localScafContigNodesNumArr, oldEndSeqLenArray, contigGraph, assemblyCycle)==FAILED) ////////////////////////////////////////////////////////
							{
								printf("line=%d, In %s(), cannot update contig information, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}
					}
				}
			}
		}else if(linkedContigsNum<=0)
		{
			printf("line=%d, In %s(), scaffoldID=%d, linkedContigsNum=%d, error!\n", __LINE__, __func__, scaffoldID, linkedContigsNum);
			return FAILED;
		}

		scaffoldItem = scaffoldItem->next;
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
short getScafContigEndSeqs(char *endSeqArray[2], int32_t *endSeqLenArray, contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph)
{
	int contigID1, contigID2, contigOrient1, contigOrient2, contigLen1, contigLen2;

	contigID1 = pContigOverlapInfo->contigID1;
	contigID2 = pContigOverlapInfo->contigID2;
	contigOrient1 = pContigOverlapInfo->orientation1;
	contigOrient2 = pContigOverlapInfo->orientation2;
	contigLen1 = contigGraph->contigItemArray[contigID1-1].contigLen;
	contigLen2 = contigGraph->contigItemArray[contigID2-1].contigLen;

	// get the lengths of the sequence at two contig ends
	if(contigLen1>=maxOverlapSeqLen)
	{
		endSeqLenArray[0] = maxOverlapSeqLen;
	}else
	{
		endSeqLenArray[0] = contigLen1;
	}

	if(contigOrient1==ORIENTATION_PLUS)
	{
		strcpy(endSeqArray[0], contigGraph->contigItemArray[contigID1-1].contigSeq+contigLen1-endSeqLenArray[0]);
	}else
	{
		strncpy(endSeqArray[0], contigGraph->contigItemArray[contigID1-1].contigSeq, endSeqLenArray[0]);
		endSeqArray[0][ endSeqLenArray[0] ] = '\0';

		if(reverseSeq(endSeqArray[0], endSeqLenArray[0])==FAILED) // reverse the sequence
		{
			printf("line=%d, In %s(), cannot reverse the sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	if(contigLen2>=maxOverlapSeqLen)
	{
		endSeqLenArray[1] = maxOverlapSeqLen;
	}else
	{
		endSeqLenArray[1] = contigLen2;
	}

	if(contigOrient2==ORIENTATION_PLUS)
	{
		strncpy(endSeqArray[1], contigGraph->contigItemArray[contigID2-1].contigSeq, endSeqLenArray[1]);
		endSeqArray[1][ endSeqLenArray[1] ] = '\0';

		if(reverseSeq(endSeqArray[1], endSeqLenArray[1])==FAILED) // reverse the sequence
		{
			printf("line=%d, In %s(), cannot reverse the sequence, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		strcpy(endSeqArray[1], contigGraph->contigItemArray[contigID2-1].contigSeq+contigLen2-endSeqLenArray[1]);
	}

	return SUCCESSFUL;
}

/**
 * Get the comparison sequence in gap filling.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getComparisonSeqInScaf(char *comparisonSeqInScaf, int32_t *comparisonSeqLenInScaf, char *endSeqArray[2], int32_t *endSeqLenArray, int32_t assemblyCycle)
{
	int32_t i;
	char reverseBase;

	if(assemblyCycle==0)
	{ // the first round
		if(endSeqLenArray[1] >= maxOverlapSeqLen)
		{
			*comparisonSeqLenInScaf = maxOverlapSeqLen;
		}else
		{
			*comparisonSeqLenInScaf = endSeqLenArray[1];
		}

		for(i=0; i<(*comparisonSeqLenInScaf); i++)
		{
			switch(endSeqArray[1][ endSeqLenArray[1]-i-1 ])
			{
				case 'A': reverseBase = 'T'; break;
				case 'C': reverseBase = 'G'; break;
				case 'G': reverseBase = 'C'; break;
				case 'T': reverseBase = 'A'; break;
				default: printf("line=%d, In %s(), unknown base [ %c ], error!\n", __LINE__, __func__, endSeqArray[1][ endSeqLenArray[1]-i-1 ]); return FAILED;
			}
			comparisonSeqInScaf[i] = reverseBase;
		}
		comparisonSeqInScaf[*comparisonSeqLenInScaf] = '\0';
	}else
	{ // the second round
		if(endSeqLenArray[0] >= maxOverlapSeqLen)
		{
			*comparisonSeqLenInScaf = maxOverlapSeqLen;
		}else
		{
			*comparisonSeqLenInScaf = endSeqLenArray[0];
		}

		//strncpy(comparisonSeqInScaf, scafContigEndSeqArr[0], comparisonSeqLenInScaf);
		strncpy(comparisonSeqInScaf, endSeqArray[0]+endSeqLenArray[0]-(*comparisonSeqLenInScaf), *comparisonSeqLenInScaf);
		comparisonSeqInScaf[*comparisonSeqLenInScaf] = '\0';

		if(reverseSeq(comparisonSeqInScaf, *comparisonSeqLenInScaf)==FAILED) // reverse the sequence
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
short prepareAssemblyInScaf(contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph, char *endSeqArray[2], int32_t *endSeqLenArray, int32_t assemblyCycle)
{
	int32_t tmpSeqLen, startRow;
	char *pseq;

	// get the contig end sequence to sacfContigSeqLastReadLen
	if(assemblyCycle==0)
	{ // the first round
		if(endSeqLenArray[0] >= readLen)
			tmpSeqLen = readLen;
		else
			tmpSeqLen = endSeqLenArray[0];
		startRow = endSeqLenArray[0] - tmpSeqLen;
		pseq = endSeqArray[0] + startRow;
		prepareAssemblyLenArr[0] = tmpSeqLen;
	}else
	{ // the second round
		if(endSeqLenArray[1] >= readLen)
			tmpSeqLen = readLen;
		else
			tmpSeqLen = endSeqLenArray[1];
		startRow = endSeqLenArray[0] - tmpSeqLen;
		pseq = endSeqArray[1] + startRow;
		prepareAssemblyLenArr[1] = tmpSeqLen;
	}

	// initialize the scafContig nodes
	if(initScafContig(contigArr, &itemNumContigArr, pseq, tmpSeqLen)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the scafContig nodes, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//initialize PEhashTable
	if(initPEHashtableInScaf(pContigOverlapInfo, contigGraph, assemblyCycle)==FAILED)
	{
		if(cleanReadsFromPEHashtable()==FAILED)
		{
			printf("line=%d, In %s, cannot clean PE hash table, error!\n", __LINE__, __func__);
			return FAILED;
		}
		printf("line=%d, In %s, cannot initialize the PE hash table, error!\n", __LINE__, __func__);
		return FAILED;
	}


	// clean the DTRow hash table
	if(cleanDTRowIndexHashtable(dtRowHashtable)==FAILED)
	{
		printf("line=%d, In %s(), contigsNum=%d, cannot clean DTRow hash table, Error!\n", __LINE__, __func__, contigsNum+1);
		return FAILED;
	}

	// initialize the reads number region
	if(itemNumContigArr>=minContigLenCheckingReadsNum)
	{
		if(initReadsNumRegSecondAssembly(itemNumContigArr)==FAILED)
		{
			printf("line=%d, In %s(), cannot initialize the reads number region, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		readsNumRatio = 1;
	}

	// initialize the decision table
	if(initAssemblingTableSecondAssembly(contigArr, itemNumContigArr, deBruijnGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize assembly table when the second assembly, error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(itemNumDecisionTable<0)
	{
		printf("line=%d, In buildContigs(), cannot init assembly table when the second assembly, error!\n", __LINE__);
		return FAILED;
	}


	// set the PEHash table margin
	if(PEGivenType>NONE_PE_GIVEN_TYPE && itemNumContigArr>minContigLenUsingPE)
	{
		//initialize PEhashTable
		if(initPEHashtableSecondAssembly(contigArr, itemNumContigArr, NO)==FAILED)
		{
			if(cleanReadsFromPEHashtable()==FAILED)
			{
				printf("In %s, cannot clean PE hash table, error!\n", __func__);
				return FAILED;
			}
			printf("In %s, cannot initialize the PE hash table before second round assembly, error!\n", __func__);
			return FAILED;
		}
	}else
	{
		// clean the PE hash table
		if(cleanReadsFromPEHashtable()==FAILED)
		{
			printf("line=%d, In %s, cannot clean PE hash table, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Initialize the scafContig nodes in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initScafContig(contigtype *contigArray, int64_t *contigNodesNum, char *seq, int32_t seq_len)
{
	int32_t i, j, baseInt;

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

		contigArray[i].index = i + 1;
		contigArray[i].base = baseInt;
		contigArray[i].ridposnum = 0;
		contigArray[i].pridposorientation = NULL;

		// added 2013-10-14
		contigArray[i].naviFlag = NAVI_SE_FLAG;
		for(j=0; j<4; j++)
		{
			contigArray[i].occNumPE[j] = contigArr[i].occNumSE[j] = 0;
			contigArray[i].occIndexPE[j] = contigArr[i].occIndexSE[j] = -1;
		}
		contigArray[i].occNumSE[baseInt] = 2;
		contigArray[i].occIndexSE[0] = baseInt;

	}

	*contigNodesNum = seq_len;

	return SUCCESSFUL;
}

/**
 * Initialize the PE hash table for gap filling.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initPEHashtableInScaf(contigOverlap_t *pContigOverlapInfo, contigGraph_t *contigGraph, int32_t assemblyCycle)
{
	int32_t i, j, tmpLeftContigIndex, tmpRightContigIndex;
	contigRead_t *contigReadArray;
	int32_t contigID, contigLen, contigEndFlag, contigEndReadNum, startContigPos, endContigPos, tmpRegLen, tmpPos;
	int32_t ridposnum, validReadOrientPEHash;
	successRead_t tmpRidposorient;


	if(cleanReadsFromPEHashtable()==FAILED)
	{
		printf("line=%d, In %s, cannot clean PE hash table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(assemblyCycle==0)
	{ // the first contig
		contigID = pContigOverlapInfo->contigID1;
		contigLen = contigGraph->contigItemArray[contigID-1].contigLen;
		if(pContigOverlapInfo->orientation1==ORIENTATION_PLUS)
		{
			// 3' end
			contigReadArray = contigGraph->contigItemArray[contigID-1].contigReadArrayEnd3;
			contigEndReadNum = contigGraph->contigItemArray[contigID-1].contigReadNumEnd3;
			contigEndFlag = 1;
		}else
		{
			// 5' end
			contigReadArray = contigGraph->contigItemArray[contigID-1].contigReadArrayEnd5;
			contigEndReadNum = contigGraph->contigItemArray[contigID-1].contigReadNumEnd5;
			contigEndFlag = 0;
		}
	}else
	{ // the second contig
		contigID = pContigOverlapInfo->contigID2;
		contigLen = contigGraph->contigItemArray[contigID-1].contigLen;
		if(pContigOverlapInfo->orientation2==ORIENTATION_PLUS)
		{
			// 3' end
			contigReadArray = contigGraph->contigItemArray[contigID-1].contigReadArrayEnd5;
			contigEndReadNum = contigGraph->contigItemArray[contigID-1].contigReadNumEnd5;
			contigEndFlag = 0;
		}else
		{
			// 5' end
			contigReadArray = contigGraph->contigItemArray[contigID-1].contigReadArrayEnd3;
			contigEndReadNum = contigGraph->contigItemArray[contigID-1].contigReadNumEnd3;
			contigEndFlag = 1;
		}
	}


	if(maxRegLenPEHash<contigAlignRegSize)
		tmpRegLen = maxRegLenPEHash;
	else
		tmpRegLen = contigAlignRegSize;

	if(contigEndFlag==1)
	{ // 3' end

		startContigPos = contigLen - tmpRegLen + 1;
		endContigPos = contigLen;
		validReadOrientPEHash = ORIENTATION_PLUS;

		for(i=contigEndReadNum-1; i>=0; i--)
		{
			if(contigReadArray[i].contigPos>=startContigPos)
			{
				if(contigReadArray[i].orientation==validReadOrientPEHash)
				{
					tmpRidposorient.rid = contigReadArray[i].readID;
					tmpRidposorient.seqlen = contigReadArray[i].seqlen;
					if(addReadToPEHashtable(&tmpRidposorient, contigReadArray[i].contigPos, SECOND_ROUND_ASSEMBLY)==FAILED)
					{
						printf("line=%d, In %s(), cannot add read %ld to PE hash table, error!\n", __LINE__, __func__, (int64_t)tmpRidposorient.rid);
						return FAILED;
					}
				}
			}else
			{
				break;
			}
		}

	}else
	{ // 5' end
		startContigPos = 1;
		endContigPos = tmpRegLen;
		validReadOrientPEHash = ORIENTATION_MINUS;

		for(i=0; i<contigEndReadNum; i++)
		{
			tmpPos = contigReadArray[i].contigPos + contigReadArray[i].seqlen - 1;
			if(tmpPos<=endContigPos && contigReadArray[i].orientation==validReadOrientPEHash)
			{
				tmpRidposorient.rid = contigReadArray[i].readID;
				tmpRidposorient.seqlen = contigReadArray[i].seqlen;

				if(addReadToPEHashtable(&tmpRidposorient, tmpPos, SECOND_ROUND_ASSEMBLY)==FAILED)
				{
					printf("line=%d, In %s(), cannot add read %lu to PE hash table, error!\n", __LINE__, __func__, (uint64_t)tmpRidposorient.rid);
					return FAILED;
				}
			}else
			{
				break;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the maximal assembly length in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMaxMinAssemblyLen(int32_t *maxAssemblyLen, int32_t *minAssemblyLen, int32_t assemblyCycle, int32_t *localScafContigNodesNumArr, int32_t gapSize)
{
	int32_t endLen1, endLen2, newSeqLen;
	int32_t newEndLen1;

	if(assemblyCycle==0)
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
 * Get the new assembled contig sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNewContigSeq(char *localScafContigSeq, int32_t *localScafContigNodesNum, contigtype *contigArray, int32_t itemNumContigArray)
{
	int32_t i, j;
	char baseCh;

	*localScafContigNodesNum = 0;
	for(i=0; i<itemNumContigArray; i++)
	{
		switch(contigArray[i].base)
		{
			case 0: baseCh = 'A'; break;
			case 1: baseCh = 'C'; break;
			case 2: baseCh = 'G'; break;
			case 3: baseCh = 'T'; break;
			default: printf("line=%d, In %s(), invalid baseInt=%d, error!\n", __LINE__, __func__, contigArray[i].base); return FAILED;
		}

		localScafContigSeq[*localScafContigNodesNum] = baseCh;
		(*localScafContigNodesNum) ++;
	}

	return SUCCESSFUL;
}

/**
 * Update contig end sequences in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateScafContigEndSeqs(char *localScafContigSeq, int32_t localScafContigNodesNum, int32_t assemblyCycle)
{
	char *pEndSeq;
	int32_t i, j, startSeqRow, *pEndSeqLen;
	int32_t prepareScafContigLen, startScafContigRow, startEndSeqRow;

	if(assemblyCycle==0)
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


	if(localScafContigNodesNum>=maxOverlapSeqLen)
	{
		startEndSeqRow = -1; // -1 indicates that the index is invalid
		startScafContigRow = localScafContigNodesNum - maxOverlapSeqLen;
	}else
	{
		if((*pEndSeqLen)+localScafContigNodesNum-prepareScafContigLen >= maxOverlapSeqLen)
			startEndSeqRow = (*pEndSeqLen) + localScafContigNodesNum - prepareScafContigLen - maxOverlapSeqLen;
		else
			startEndSeqRow = 0;

		startScafContigRow = prepareScafContigLen;
	}

	// update the scafContig end sequence
	if(startEndSeqRow==-1)
	{
		startSeqRow = 0;
	}else if(startEndSeqRow==0)
	{
		startSeqRow = (*pEndSeqLen);
	}else if(startEndSeqRow>0)
	{
		for(i=startEndSeqRow; i<(*pEndSeqLen); i++) pEndSeq[i-startEndSeqRow] = pEndSeq[i];
		pEndSeq[(*pEndSeqLen)-startEndSeqRow] = '\0';
		startSeqRow = (*pEndSeqLen) - startEndSeqRow;
	}else
	{
		printf("line=%d, In %s(), startEndSeqRow=%d, error!\n", __LINE__, __func__, startEndSeqRow);
		return FAILED;
	}

	j = startSeqRow;
	for(i=startScafContigRow; i<localScafContigNodesNum; i++)
		pEndSeq[j++] = localScafContigSeq[i];
	pEndSeq[j] = '\0';
	*pEndSeqLen = j;

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
short detectOverlapsInScaf(int32_t *successFilledFlag, int32_t *overlapLen, int32_t *newGapSize, int32_t *breakFlag, int32_t gapSize, int32_t *localScafContigNodesNumArray, int32_t assemblyCycle)
{
	char *contigEndSeq, *comparisonSeq, *otherContigEndSeq;
	int32_t *contigEndSeqLen, comparisonSeqLen, *otherContigEndSeqLen;
	int32_t i, j, overlapLenExact, overlapLenAlignment, mismatchNum, overlapLenAdjust, tmpGapSize;

	if(assemblyCycle==0)
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

	//if((gapSize>0.6*meanSizeInsert && gapSize>2*readLen) && tmpGapSize>readLen) // 2012-11-19
	if((gapSize>meanSizeInsert && gapSize>2*readLen) && tmpGapSize>meanSizeInsert) // 2014-01-13
	{
#if (DEBUG_FLAG==YES)
		printf("line=%d, In %s(), broken.\n", __LINE__, __func__);
#endif
		*overlapLen = 0;
		*newGapSize = tmpGapSize;
		*successFilledFlag = NO;
		*breakFlag = YES;
	}
	else if(tmpGapSize<=standardDev)
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
		if(overlapLenExact>=minOverlapThres && tmpGapSize>-gapSizeSdevFactorGapFilling*standardDev)
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

					if(tmpGapSize<-meanSizeInsert)
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
short computeNewGapSizeInScaf(int32_t *tmpGapSize, int32_t gapSize, int32_t *localScafContigNodesNumArray, int32_t assemblyCycle)
{
	int32_t newSeqLen1, newSeqLen2;

	if(assemblyCycle==0)
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
short updateContigOverlapInfoInScaf(contigOverlap_t *pContigOverlapInfo, int32_t successFilledFlag, int32_t overlapLen, int32_t newGapSize, int32_t breakFlag)
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
short updateContigInfoInScaf(contigOverlap_t *pContigOverlapInfo, char **localScafContigSeqArr, int32_t *localScafContigNodesNumArr, int32_t *oldEndSeqLenArr, contigGraph_t *contigGraph, int32_t assemblyCycle)
{
	int32_t contigID1, contigID2, contigOrient1, contigOrient2;

	// get the overlap information
	contigID1 = pContigOverlapInfo->contigID1;
	contigID2 = pContigOverlapInfo->contigID2;
	contigOrient1 = pContigOverlapInfo->orientation1;
	contigOrient2 = pContigOverlapInfo->orientation2;

	if(assemblyCycle==0)
	{ // the first round
		if(localScafContigNodesNumArr[0] > prepareAssemblyLenArr[0])
		{ // new bases are appended
			if(updateSingleContigInfoInScaf(contigGraph->contigItemArray+contigID1-1, contigOrient1, 0, localScafContigSeqArr[0], localScafContigNodesNumArr[0], scafContigEndSeqArr[0], scafContigEndSeqLenArr[0], oldEndSeqLenArr[0], prepareAssemblyLenArr[0])==FAILED)
			{
				printf("line=%d, In %s(), cannot update single contig information in local assembly for contig %d, error!\n", __LINE__, __func__, contigID1);
				return FAILED;
			}

			// update the other contig information
			if(updateOtherContigInfoInScaf(contigGraph->contigItemArray+contigID2-1, contigOrient2, scafContigEndSeqArr[1], scafContigEndSeqLenArr[1], oldEndSeqLenArr[1])==FAILED)
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
		if(updateSingleContigInfoInScaf(contigGraph->contigItemArray+contigID1-1, contigOrient1, 0, localScafContigSeqArr[0], localScafContigNodesNumArr[0], scafContigEndSeqArr[0], scafContigEndSeqLenArr[0], oldEndSeqLenArr[0], prepareAssemblyLenArr[0])==FAILED)
		{
			printf("line=%d, In %s(), cannot update single contig information in local assembly for contig %d, error!\n", __LINE__, __func__, contigID1);
			return FAILED;
		}

		// process the second contig
		if(updateSingleContigInfoInScaf(contigGraph->contigItemArray+contigID2-1, contigOrient2, 1, localScafContigSeqArr[1], localScafContigNodesNumArr[1], scafContigEndSeqArr[1], scafContigEndSeqLenArr[1], oldEndSeqLenArr[1], prepareAssemblyLenArr[1])==FAILED)
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
short updateSingleContigInfoInScaf(contigGraphItem_t *contigItem, int32_t contigOrient, int32_t contigIndex, char *localScafContigSeq, int32_t scafContigNodesNum, char *endSeq, int32_t endSeqLen, int32_t oldEndSeqLen, int32_t prepareAssemblyLen)
{
	int32_t i, j, contigLen;
	char *contigSeq, *newContigSeq, *tmpSeq, tmpBase;
	int32_t newEndSeqLen, tmpSeqLen;
	int32_t newContigLen;
	int32_t endScafContigRow, startEndSeqRow;

	contigLen = contigItem->contigLen;
	contigSeq = contigItem->contigSeq;

	// allocate the memory for tmpStr
	newEndSeqLen = endSeqLen;
	tmpSeqLen = scafContigNodesNum + newEndSeqLen - oldEndSeqLen;
	tmpSeq = (char *) malloc((tmpSeqLen+1)* sizeof(char));
	if(tmpSeq==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(scafContigNodesNum>=oldEndSeqLen)
		endScafContigRow = scafContigNodesNum - oldEndSeqLen - 1;
	else
		endScafContigRow = -1;
	startEndSeqRow = oldEndSeqLen - scafContigNodesNum;

	// copy the new end sequence to tmpSeq
	if(startEndSeqRow>=0)
	{
		strcpy(tmpSeq, endSeq+startEndSeqRow);
	}else
	{
		for(i=0; i<=endScafContigRow; i++)
			tmpSeq[i] = localScafContigSeq[i];
		tmpSeq[i] = '\0';

		strcpy(tmpSeq+i, endSeq);
	}

	// ############################# Debug information ###########################
	if(strlen(tmpSeq)!=tmpSeqLen)
	{
		printf("line=%d, In %s(), strlen(tmpSeq)=%d != tmpSeqLen=%d, error!\n", __LINE__, __func__, (int)strlen(tmpSeq), tmpSeqLen);
		return FAILED;
	}
	// ############################# Debug information ###########################


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
				contigItem->contigLen = newContigLen;
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
				free(contigItem->contigSeq);
				contigItem->contigSeq = newContigSeq;
				contigItem->contigLen = newContigLen;
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
				contigItem->contigLen = newContigLen;

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
				free(contigItem->contigSeq);
				contigItem->contigSeq = newContigSeq;
				contigItem->contigLen = newContigLen;
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
				contigItem->contigLen = newContigLen;

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
				free(contigItem->contigSeq);
				contigItem->contigSeq = newContigSeq;
				contigItem->contigLen = newContigLen;
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
				contigItem->contigLen = newContigLen;

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
				free(contigItem->contigSeq);
				contigItem->contigSeq = newContigSeq;
				contigItem->contigLen = newContigLen;
			}
		}
	}

	// ############################# Debug information ###########################
	if(strlen(contigItem->contigSeq)!=contigItem->contigLen)
	{
		printf("line=%d, In %s(), strlen(contigSeq)=%d != contigLen=%d, error!\n", __LINE__, __func__, (int32_t)strlen(contigItem->contigSeq), contigItem->contigLen);
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
short updateOtherContigInfoInScaf(contigGraphItem_t *contigItem, int32_t contigOrient, char *contigEndSeq, int32_t endSeqLen, int32_t oldEndSeqLen)
{
	int32_t i, j, newEndSeqLen;
	int32_t contigLen, newContigLen;
	char *contigSeq, *newContigSeq;

	contigLen = contigItem->contigLen;
	contigSeq = contigItem->contigSeq;

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

			contigItem->contigLen = newContigLen;

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

			free(contigItem->contigSeq);
			contigItem->contigSeq = newContigSeq;
			contigItem->contigLen = newContigLen;
		}
	}else
	{ // minus orientation
		if(newEndSeqLen==oldEndSeqLen)
		{ // only update the end sequence
			strcpy(contigSeq+contigLen-oldEndSeqLen, contigEndSeq);

		}else if(newEndSeqLen<oldEndSeqLen)
		{ // update the end sequence and remained sequence
			strcpy(contigSeq+contigLen-oldEndSeqLen, contigEndSeq);

			contigItem->contigLen -= oldEndSeqLen - newEndSeqLen;

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

			free(contigItem->contigSeq);
			contigItem->contigSeq = newContigSeq;
			contigItem->contigLen = newContigLen;
		}
	}

	// ############################ Debug information ######################
	if(contigItem->contigLen!=strlen(contigItem->contigSeq))
	{
		printf("line=%d, In %s(), contigLen=%d != strlen(contigSeq)=%d, error!\n", __LINE__, __func__, contigItem->contigLen, (int32_t)strlen(contigItem->contigSeq));
		return FAILED;
	}
	// ############################ Debug information ######################

	return SUCCESSFUL;
}
