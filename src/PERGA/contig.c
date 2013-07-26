
#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Build contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short buildContigs(char *contigFile, char *graphFileName, char *(*occPointFileArray)[6])
{
	printf("\n============= Begin building contigs, please wait ... =============\n");

	struct timeval tp_start,tp_end;
	double time_used;
	gettimeofday(&tp_start,NULL);

	char hangingFile[256];
	int32_t i, turnContigIndex;
	int64_t validReadNum;

	int32_t percent, percentNum, tmp_gapSize;

	initFirstKmerThreshold();
	if(initMemory(occPointFileArray)==FAILED)
	{
		printf("line=%d, In %s(), cannot init the memory of two tables, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//++++++++++++++++++++++++++++++++++++
	if(PEGivenType>=NONE_PE_GIVEN_TYPE)
	{
		if(PEGivenType!=BOTH_PE_GIVEN_TYPE)
		{
			// estimate the insert size and standard deviation of fragment library
			if(estimateInsertSizeAndSdev(graphFileName)==FAILED)
			{
				printf("line=%d, In %s(), cannot estimate the insert size and standard deviation of fragment library, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		if(initPEHashParas()==FAILED)
		{
			printf("line=%d, In %s(), cannot initialize the PEHash table parameters, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	fpContigsBase = fopen(contigFile, "w");
	if(fpContigsBase==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigFile);
		return FAILED;
	}
	if(hangingContigOutFlag==YES)
	{
		strcpy(hangingFile, contigFile);
		strcat(hangingFile, ".hang");
		fpContigsHanging = fopen(hangingFile, "w");
		if(fpContigsHanging==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, hangingFile);
			return FAILED;
		}
	}

	if(errorCorrectionFlag==YES)
	{
		fpReadCorrected = fopen(readCorrectedFile, "w");
		if(fpReadCorrected==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readCorrectedFile);
			return FAILED;
		}
	}

	basesNum = 0;
	localContigID = 0;
	percent = 0;
	percentNum = 0;
	successReadNum = 0;
	contigsNum = 0;
	validReadNum = deBruijnGraph->readSet->totalValidItemNumRead;
	totalAlignedSuccessReadNum = 0;
	totalErrReadNum = 0;
	totalReadsNumCorreted = 0;
	itemNumContigsLenArr = 0;


	if(errorCorrectionFlag==YES)
		correctionAllowed = YES;
	else
		correctionAllowed = NO;

	for(i=0; i<11; i++) errNumArr[i] = 0;

	for(assemblyCycle=1; assemblyCycle<=2; assemblyCycle++)
	{
		if(initFirstKmerBounder(&lowerBoundFirstKmer, &upperBoundFirstKmer, assemblyCycle, averKmerOcc)==FAILED)
		{
			printf("line=%d, In %s(), cannot initialize the bounder for first k-mers, error!\n", __LINE__, __func__);
			return FAILED;
		}

		kmerIndex = 0;
		while(kmerIndex < hashTableSize)
		{
			itemNumContigArr = 0;
			itemNumDecisionTable = 0;
			validHeadRowContigArr = 0;
			validTailRowContigArr = 0;
			successContigIndex = -1;
			assemblyRound = FIRST_ROUND_ASSEMBLY;  // first round assembly
			lockedReadsNum = 0;
			this_successReadNum = 0;
			readsNumInPEHashArr = 0;
			regLenPEHash = 0;
			turnContigIndex = 0;
			allowedUpdatePEHashArrFlag = YES;
			localContigID ++;
			refStrandContig = -1;
			refPosSoildFlag = NO;
			refPosContig = -1;
			refPosMatchFlag = NO;
			refPosFirstBase[0] = NO;
			refPosFirstBase[1] = -1;
			refPosFirstBase[2] = -1;
			readsNumRatio = 1;

			// get the first kmers and its sequences for assembly of a contig
			if(getFirstKmers(&kmerIndex, &firstKmer)==FAILED)
			{
				printf("line=%d, In %s(), cannot get first kmers, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(kmerIndex>=hashTableSize)
			{ // reached the bottom of the k-mer hash table
				break;
			}

			// initialize the contig chain
			if(initContig()==FAILED)
			{
				printf("line=%d, In %s(), cannot initialize the contig nodes, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// initialize the decision table
			if(addFirstKmerToDecisionTable(kmers)==FAILED)
			{
				printf("line=%d, In %s(), cannot initialize the decision table, error!\n", __LINE__, __func__);
				return FAILED;
			}


			if(setEmptyNaviOccQueue(naviOccQueue, &itemNumNaviOccQueue, &frontRowNaviOccQueue, &rearRowNaviOccQueue)==FAILED)
			{
				printf("line=%d, In %s(), cannot initialize the empty navigation occurrence queue, error!\n", __LINE__, __func__);
				return FAILED;
			}

			//while(kmers[0]||kmers[1])
			while(naviSuccessFlag==NAVI_SUCCESS)
			{
#if (DEBUG_CONTIG_CHECK==YES)
				// ############################ Debug information ##############################
				if(localContigID==8 && itemNumContigArr>=105500 && assemblyRound!=FIRST_ROUND_ASSEMBLY)
				{
					printf("localContigID=%ld, contigID=%d, itemNumContigArr=%ld, assemblyRound=%d\n", localContigID, contigsNum+1, itemNumContigArr, assemblyRound);
				}
				// ############################ Debug information ##############################
#endif


				// initialize or update the PE hash table
				if(PEGivenType>NONE_PE_GIVEN_TYPE && itemNumContigArr>=minContigLenUsingPE)
				{
					if(updatePEHashTable(itemNumContigArr, assemblyRound)==FAILED)
					{
						printf("line=%d, In %s(), localContigID=%ld, contigID=%d, itemNumContigArr=%ld, cannot update the PE hash table, error!\n", __LINE__, __func__, localContigID, contigsNum+1, itemNumContigArr);
						return FAILED;
					}

					// ########################### Debug information ########################
					//if(assemblyRound==SECOND_ROUND_ASSEMBLY && itemNumContigArr-hashRegRightContig->index<minContigLenUsingPE-1)
					//{
					//	printf("line=%d, In %s(), assemblyRound=%d, itemNumContigArr=%d, hashRegRightContig->index=%d, error!\n", __LINE__, __func__, assemblyRound, itemNumContigArr, hashRegRightContig->index);
					//	return FAILED;
					//}
					// ########################### Debug information ########################


					//if(readsNumInPEHashArr>=minReadsNumPEHashThres && regLenPEHash>=minRegLenUsingPE)
					if(readsNumInPEHashArr>0 && regLenPEHash>=minRegLenUsingPE)
					{
						if(getNextKmerByMix(itemNumContigArr, assemblyRound)==FAILED)
						{
							printf("line=%d, In %s(), localContigID=%ld, cannot get the next kmer by mix, error!\n", __LINE__, __func__, localContigID);
							return FAILED;
						}

					}else
					{
						navigationFlag = NAVI_SE_FLAG;
						if(getNextKmerBySE(itemNumContigArr)==FAILED)
						{
							printf("line=%d, In %s(), localContigID=%ld, cannot get next kmer, error!\n", __LINE__, __func__, localContigID);
							return FAILED;
						}


#if (SVM_NAVI==YES)
						if((successContigIndex>0 && itemNumContigArr-successContigIndex>50) || readsNumRatio<0.3*minReadsNumRatioThres)
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
							//if(maxOccSE-secondOccSE<2 && (successContigIndex>0 &&  itemNumContigArr-successContigIndex>15))
//							if(maxOccSE-secondOccSE<2)
//								naviSuccessFlag = NAVI_FAILED;

//							if(naviSuccessFlag==NO && svmFeatureArr[1]/svmFeatureArr[0]<0.5)
//								naviSuccessFlag = NAVI_SUCCESS;

//							if(svmFeatureArr[0]==svmFeatureArr[1] || (svmFeatureArr[1]>=2 && svmFeatureArr[1]/svmFeatureArr[0]>=0.7))
//								naviSuccessFlag = NAVI_FAILED;

							//if(localContigID==13 /*&& itemNumContigArr>=14805*/ && assemblyRound==FIRST_ROUND_ASSEMBLY)
//							if(naviSuccessFlag==NAVI_FAILED)
//							{
//								decideByMalignSE(&naviSuccessFlag, decisionTable, itemNumDecisionTable);
//							}
						}
#else
						//if(maxOccSE==secondOccSE)
						if((maxOccSE==secondOccSE) || (secondOccSE>=2 && secondOccSE>=SEC_MAX_OCC_RATIO_SVM*maxOccSE))
						{
							naviSuccessFlag = NAVI_FAILED;
							//kmers[0] = kmers[1] = NULL;
						}
#endif

#if(USING_RULES==YES)
						if(successContigIndex>0 &&  itemNumContigArr-successContigIndex > 0.6*readLen)  // added 2012-11-08
						{
							naviSuccessFlag = NAVI_FAILED;
							kmers[0] = kmers[1] = NULL;
						}
						//else if(itemNumContigArr<3*readLen && (successContigIndex>0 &&  (itemNumContigArr-successContigIndex > 0.3*readLen || itemNumContigArr-successContigIndex > 15)))  // added 2012-11-11, deleted 2012-11-23
						else if(averKmerOcc>15 && itemNumContigArr<3*readLen && (successContigIndex>0 &&  (itemNumContigArr-successContigIndex > 0.3*readLen || itemNumContigArr-successContigIndex > 15)))  // added 2012-11-23, deleted 2012-12-29
						{
							naviSuccessFlag = NAVI_FAILED;
							kmers[0] = kmers[1] = NULL;
						}
						//else if(itemNumContigArr<3*readLen && secondOccSE>15*minKmerOccSE && secondOccSE/maxOccSE>0.6)  // added 2012-11-13
						//else if(itemNumContigArr<3*readLen && secondOccSE>5*minKmerOccSE && secondOccSE/maxOccSE>0.75 && readsNumRatio<0.4)  // added 2012-11-16, deleted 2012-11-28
						else if(itemNumContigArr<3*readLen && secondOccSE>3*minKmerOccSE && secondOccSE/maxOccSE>0.6)  // added 2012-11-28
						{
							naviSuccessFlag = NAVI_FAILED;
							kmers[0] = kmers[1] = NULL;
						}
//						else if(itemNumContigArr>=3*readLen && maxOccSE>2*maxOccNumFaiedPE && readsNumRatio<2.5*minReadsNumRatioThres)  // added 2012-11-14
//						{
//							// compute the maximal gap size in contig tail region
//							if(computeGapSizeInContig(&tmp_gapSize, contigArr, itemNumContigArr, assemblyRound)==FAILED)
//							{
//								printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
//								return FAILED;
//							}
//
//							if(tmp_gapSize>10)
//							{
//								printf("==line=%d, itemNumContigArr=%ld, tmp_gapSize=%d\n", __LINE__, itemNumContigArr, tmp_gapSize);
//								naviSuccessFlag = NAVI_FAILED;
//								kmers[0] = kmers[1] = NULL;
//							}
//						}

						else if(itemNumContigArr>=3*readLen && secondOccSE>0 && readsNumRatio<3*minReadsNumRatioThres)  // added 2013-01-20
						{
							// compute the maximal gap size in contig tail region
							if(computeGapSizeInContig(&tmp_gapSize, contigArr, itemNumContigArr, assemblyRound)==FAILED)
							{
								printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
								return FAILED;
							}

							if(tmp_gapSize>0.3*readLen)		// added 2013-01-20
							{
#if (DEBUG_OUTPUT==YES)
								printf("==line=%d, itemNumContigArr=%ld, tmp_gapSize=%d\n", __LINE__, itemNumContigArr, tmp_gapSize);
#endif
								naviSuccessFlag = NAVI_FAILED;
								kmers[0] = kmers[1] = NULL;
							}
						}
#endif
					}
				}else
				{
					navigationFlag = NAVI_SE_FLAG;
					if(getNextKmerBySE(itemNumContigArr)==FAILED)
					{
						printf("line=%d, In %s(), localContigID=%ld, cannot get next kmer, error!\n", __LINE__, __func__, localContigID);
						return FAILED;
					}


#if (SVM_NAVI==YES)
					if((successContigIndex>0 && itemNumContigArr-successContigIndex>50) || readsNumRatio<0.3*minReadsNumRatioThres)
						naviSuccessFlag = NAVI_FAILED;
					else if(secondOccSE>0)								// deleted 2013-02-26
					//if(naviSuccessFlag==NAVI_SUCCESS && secondOccSE>0 && successContigIndex>0)		// added 2013-02-26
					{
//						sumSecondOccSE = 0;
//						for(i=0; i<4; i++) if(i!=occsNumIndexSE[0]) sumSecondOccSE += occsNumSE[i];

						svmFeatureArr[0] = occsNumSE[occsNumIndexSE[0]];
						svmFeatureArr[1] = occsNumSE[occsNumIndexSE[1]];
						//svmFeatureArr[1] = sumSecondOccSE;
						svmFeatureArr[2] = readsNumRatio;
						//svmFeatureArr[2] = svmFeatureArr[0] / svmFeatureArr[1];
						if(successContigIndex>0)
						{
							// compute the maximal gap size in contig tail region
//							if(computeGapSizeInContig(&tmp_gapSize, contigArr, itemNumContigArr, assemblyRound)==FAILED)
//							{
//								printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
//								return FAILED;
//							}

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

//						if(naviSuccessFlag==NAVI_FAILED)
//						{
//							decideByMalignSE(&naviSuccessFlag, decisionTable, itemNumDecisionTable);
//						}

						//if(successContigIndex==-1 && (secondOccSE>=0.7*maxOccSE || maxOccSE-secondOccSE<2))		// added 2013-04-01
						//if(maxOccSE-secondOccSE<=2 && (successContigIndex>0 &&  itemNumContigArr-successContigIndex>=8))
//						if(maxOccSE-secondOccSE<2)
//							naviSuccessFlag = NAVI_FAILED;

//						if(naviSuccessFlag==NAVI_FAILED && svmFeatureArr[1]/svmFeatureArr[0]<0.5)
//							naviSuccessFlag = NAVI_SUCCESS;

//						if(svmFeatureArr[0]==svmFeatureArr[1] || (svmFeatureArr[1]>=2 && svmFeatureArr[1]/svmFeatureArr[0]>=0.7))
//							naviSuccessFlag = NAVI_FAILED;
					}
#else
					//if(maxOccSE==secondOccSE)
					if((maxOccSE==secondOccSE) || (secondOccSE>=2 && secondOccSE>=SEC_MAX_OCC_RATIO_SVM*maxOccSE))
					{
						naviSuccessFlag = NAVI_FAILED;
						//kmers[0] = kmers[1] = NULL;
					}
#endif

#if(USING_RULES==YES)
					if(successContigIndex>0 &&  itemNumContigArr-successContigIndex > 0.6*readLen)  // added 2012-11-08
					{
						naviSuccessFlag = NAVI_FAILED;
						kmers[0] = kmers[1] = NULL;
					}
					//else if(itemNumContigArr<3*readLen && (successContigIndex>0 && (itemNumContigArr-successContigIndex > 0.3*readLen || itemNumContigArr-successContigIndex > 15)))  // added 2012-11-11, deleted 2012-11-23
					else if(averKmerOcc>15 && itemNumContigArr<3*readLen && (successContigIndex>0 && (itemNumContigArr-successContigIndex > 0.3*readLen || itemNumContigArr-successContigIndex > 15)))  // added 2012-11-23, deleted 2012-12-29
					{
						naviSuccessFlag = NAVI_FAILED;
						kmers[0] = kmers[1] = NULL;
					}
					//else if(itemNumContigArr<3*readLen && secondOccSE>15*minKmerOccSE && secondOccSE/maxOccSE>0.6)  // added 2012-11-13
					//else if(itemNumContigArr<3*readLen && secondOccSE>5*minKmerOccSE && secondOccSE/maxOccSE>0.75 && readsNumRatio<0.4)  // added 2012-11-16, deleted 2012-11-28
					else if(itemNumContigArr<3*readLen && secondOccSE>3*minKmerOccSE && secondOccSE/maxOccSE>0.6)  // added 2012-11-28
					{
						naviSuccessFlag = NAVI_FAILED;
						kmers[0] = kmers[1] = NULL;
					}
//					else if(itemNumContigArr>=100*readLen && maxOccSE>2*maxOccNumFaiedPE && readsNumRatio<2*minReadsNumRatioThres)  // added 2012-11-14
//					{
//						// compute the maximal gap size in contig tail region
//						if(computeGapSizeInContig(&tmp_gapSize, contigArr, itemNumContigArr, assemblyRound)==FAILED)
//						{
//							printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
//							return FAILED;
//						}
//
//						if(tmp_gapSize>10)
//						{
//							printf("==line=%d, itemNumContigArr=%ld, tmp_gapSize=%d\n", __LINE__, itemNumContigArr, tmp_gapSize);
//							naviSuccessFlag = NAVI_FAILED;
//							kmers[0] = kmers[1] = NULL;
//						}
//					}

					else if(itemNumContigArr>=3*readLen && secondOccSE>0 && readsNumRatio<3*minReadsNumRatioThres)  // added 2013-01-20
					{
						// compute the maximal gap size in contig tail region
						if(computeGapSizeInContig(&tmp_gapSize, contigArr, itemNumContigArr, assemblyRound)==FAILED)
						{
							printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
							return FAILED;
						}

						if(tmp_gapSize>0.3*readLen)		// added 2013-01-20
						{
#if (DEBUG_OUTPUT==YES)
							printf("==line=%d, itemNumContigArr=%ld, tmp_gapSize=%d\n", __LINE__, itemNumContigArr, tmp_gapSize);
#endif
							naviSuccessFlag = NAVI_FAILED;
							kmers[0] = kmers[1] = NULL;
						}
					}
#endif
				}

				// check the reads number in the sub read region
				//if(kmers[0] || kmers[1])
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
				}

				// draw points
#if(DRAW_CURVE_FLAG==YES)
				if(refPosContig>0)
				{
					if(updateRefPosContig(successReadsArr, itemNumSuccessReadsArr, refArr)==FAILED)
					{
						printf("line=%d, In %s(), cannot update ref position of contig, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

				// draw the points for maxOcc and secondOcc when there are branches
				//if(successContigIndex>0 && ((navigationFlag==NAVI_PE_FLAG && secondOccPE>0) || (navigationFlag!=NAVI_PE_FLAG && secondOccSE>0))) // deleted 2013-02-27
				//if((navigationFlag==NAVI_PE_FLAG && secondOccPE>0) || (navigationFlag!=NAVI_PE_FLAG && secondOccSE>0))		// added 2013-02-27
				if((navigationFlag==NAVI_PE_FLAG && secondOccPE>0) || (navigationFlag!=NAVI_PE_FLAG && secondOccSE>0))		// added 2013-02-27
				{
					if(drawOccPoints(navigationFlag, naviSuccessFlag, refArr)==FAILED)
					{
						printf("line=%d, In %s(), cannot draw points for maxOcc and secondOcc, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}
#endif


				//if(kmers[0]==NULL && kmers[1]==NULL)
				if(naviSuccessFlag==NAVI_FAILED)
				{
					if(successContigIndex<0)
					{ // no successful reads, then assembly failed
						//printf("line=%d, In %s(), localContigID=%ld, contigsNum=%d, assemblyRound=%d, itemNumContigArr=%ld, the successContigIndex<=0!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
						break;
					}

#if (DEBUG_OUTPUT==YES)
					printf("localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, itemNumDecisionTable=%d\n", localContigID, contigsNum+1, assemblyRound, itemNumContigArr, itemNumDecisionTable);
					if(navigationFlag==NAVI_PE_FLAG)
						printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
					else if(navigationFlag==NAVI_SE_FLAG)
						printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
					else
					{
						printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
						printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
					}
					if(successContigIndex>0)
						printf("\tdistance=%ld, readsNumRatio=%.2f\n", itemNumContigArr-successContigIndex, readsNumRatio);
					if(DRAW_CURVE_FLAG==YES)
						printf("refPosSoildFlag=%d, refStrandContig=%d, refPosContig=%d, refPosMatchFlag=%d\n", refPosSoildFlag, refStrandContig, refPosContig, refPosMatchFlag);
#endif

					// prepare for next round assembly
					if(assemblyRound==FIRST_ROUND_ASSEMBLY)
					{ // finished first round, then start the second round

						assemblyRound ++;
						turnContigIndex = itemNumContigArr;

						int returnCode = initSecondAssembly();
						if(returnCode==FAILED)
						{
							break;
						}else if(returnCode==ERROR)
						{
							return FAILED;
						}

						continue;

					}else
					{ // finish the second round assembly, then the whole assembly finished

						// ############################ Debug information ##############################
						//if(successContigIndex>0 && itemNumContigArr>=CONTIG_LEN_THRESHOLD)
						//{
						//	printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%d, successContigIndex>0, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
						//	return FAILED;
						//}
						// ############################ Debug information ##############################

						break;
					}
				}

				if(navigationFlag==NAVI_PE_FLAG)
				{
					if(updateNaviOccQueue(naviOccQueue, &itemNumNaviOccQueue, &frontRowNaviOccQueue, &rearRowNaviOccQueue, maxOccPE)==FAILED)
					{
						printf("line=%d, In %s(), localContigID=%ld, cannot update the navigation occurrence queue, error!\n", __LINE__, __func__, localContigID);
						return FAILED;
					}
				}else
				{
					if(updateNaviOccQueue(naviOccQueue, &itemNumNaviOccQueue, &frontRowNaviOccQueue, &rearRowNaviOccQueue, maxOccSE)==FAILED)
					{
						printf("line=%d, In %s(), localContigID=%ld, cannot update the navigation occurrence queue, error!\n", __LINE__, __func__, localContigID);
						return FAILED;
					}
				}

#if (DEBUG_CONTIG_CHECK==YES)
				// ############################ Debug information ##############################
				if(localContigID==8 && itemNumContigArr>=105500 && assemblyRound!=FIRST_ROUND_ASSEMBLY)
				{
					printf("localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, itemNumDecisionTable=%d\n", localContigID, contigsNum+1, assemblyRound, itemNumContigArr, itemNumDecisionTable);
					if(navigationFlag==NAVI_PE_FLAG)
						printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
					else if(navigationFlag==NAVI_SE_FLAG)
						printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
					else
					{
						printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
						printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
					}
					if(successContigIndex)
						printf("\tdistance=%ld, readsNumRatio=%.2f\n", itemNumContigArr-successContigIndex, readsNumRatio);
					if(DRAW_CURVE_FLAG==YES)
						printf("refPosSoildFlag=%d, refStrandContig=%d, refPosContig=%d, refPosMatchFlag=%d\n", refPosSoildFlag, refStrandContig, refPosContig, refPosMatchFlag);
				}
				// ############################ Debug information ##############################
#endif

				itemNumContigArr ++;

				// Append a base to contig tail
				if(addContigBase(kmerSeqIntAssembly[entriesPerKmer-1] & 3)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, cannot add a contig base, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
					return FAILED;
				}

				// update the decision table according to kmers
				if(updateDecisionTable(kmers, kmerSeqIntAssembly[entriesPerKmer-1] & 3)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, cannot update decision table, error!\n", __LINE__, __func__, localContigID);
					return FAILED;
				}

				// Update the reads status in decision table
				if(updateAssemblingreadsStatus()==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, cannot update reads status in decision table, error!\n", __LINE__, __func__, localContigID);
					return FAILED;
				}

				// Update the locked reads and their total number
				if(updateLockedReads()==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, cannot update locked reads, error!\n", __LINE__, __func__, localContigID);
					return FAILED;
				}

				// update the finished reads in decision table, and record the successful reads into successful reads array
				if(updateFinishedReadsInDecisionTable()==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, cannot update finished reads in decision table, error!\n", __LINE__, __func__, localContigID);
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

					successContigIndex = itemNumContigArr;
					this_successReadNum += itemNumSuccessReadsArr;
				}

#if(DRAW_CURVE_FLAG==YES)
				//if(refPosContig==-1)		// deleted 2013-03-01
				if(refPosSoildFlag==NO)		// added 2013-03-01
				{
					if(itemNumSuccessReadsArr>0)
					{
						// compute the ref position of contig
						if(computeRefPosContig(successReadsArr, itemNumSuccessReadsArr, refArr)==FAILED)
						{
							printf("line=%d, In %s(), cannot compute ref position of contig, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
				}
#endif

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

					number_of_overlap_less_than_threshold ++;

#if (DEBUG_OUTPUT==YES)
					printf("===localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, itemNumDecisionTable=%d\n", localContigID, contigsNum+1, assemblyRound, itemNumContigArr, itemNumDecisionTable);
					if(navigationFlag==NAVI_PE_FLAG)
						printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
					else if(navigationFlag==NAVI_SE_FLAG)
						printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
					else
					{
						printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
						printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
					}
					if(successContigIndex>0)
						printf("\tdistance=%ld, readsNumRatio=%.2f\n", itemNumContigArr-successContigIndex, readsNumRatio);
					if(DRAW_CURVE_FLAG==YES)
						printf("refPosSoildFlag=%d, refStrandContig=%d, refPosContig=%d, refPosMatchFlag=%d\n", refPosSoildFlag, refStrandContig, refPosContig, refPosMatchFlag);
#endif


					if(assemblyRound==SECOND_ROUND_ASSEMBLY && itemNumContigArr<CONTIG_LEN_THRESHOLD)
					{
						break;
					}

					// prepare for the second round assembly
					if(assemblyRound==FIRST_ROUND_ASSEMBLY)
					{

						assemblyRound ++;
						turnContigIndex = itemNumContigArr;

						int returnCode = initSecondAssembly();
						if(returnCode==FAILED)
						{
							break;
						}else if(returnCode==ERROR)
						{
							return FAILED;
						}

					}else
					{

						// ############################ Debug information ##############################
						//if(successContigIndex>0 && itemNumContigArr>=CONTIG_LEN_THRESHOLD)
						//{
						//	printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, successContigIndex<=0, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
						//	return FAILED;
						//}
						// ############################ Debug information ##############################

						break;
					}
				}
/*
				else if(itemNumContigArr-successContigIndex >= 7)
				{

					if(calcAverOccNaviOccQueue(&averOccNumNaviOccQueue, naviOccQueue, itemNumNaviOccQueue)==FAILED)
					{
						printf("line=%d, In %s(), localContigID=%ld, contigID=%d, cannot compute the average occurrence in navigation occurrence queue, error!\n", __LINE__, __func__, localContigID, contigsNum+1);
						return FAILED;
					}

					// get the length of low occurrence region
					if(getLowOccLenNaviOccQueue(&lowOccNum, naviOccQueue, itemNumNaviOccQueue, frontRowNaviOccQueue)==FAILED)
					{
						printf("line=%d, In %s(), localContigID=%ld, contigID=%d, cannot get the low occurrence number in navigation occurrence queue, error!\n", __LINE__, __func__, localContigID, contigsNum+1);
						return FAILED;
					}
					printf("#### itemNumContigArr=%d, distance=%d, readsNumRatio=%.2f, averOccNumNaviOccQueue=%.2f, lowOccNum=%d\n", itemNumContigArr, itemNumContigArr-successContigIndex, readsNumRatio, averOccNumNaviOccQueue, lowOccNum);

					//if((lowOccNum>2) || ((navigationFlag==NAVI_PE_FLAG && averOccNumNaviOccQueue<2*minKmerOccPE) || (navigationFlag==NAVI_SE_FLAG && averOccNumNaviOccQueue<2*minKmerOccSE)) || (readsNumRatio>2 || readsNumRatio<0.3) || ((navigationFlag==NAVI_PE_FLAG && secondOccPE>=minKmerOccPE && secondOccPE/maxOccPE>=0.2) || (navigationFlag==NAVI_SE_FLAG && secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>=0.15)))
					//if((lowOccNum>2) || (readsNumRatio>2 || readsNumRatio<0.3) || ((navigationFlag==NAVI_PE_FLAG && secondOccPE>=minKmerOccPE && secondOccPE/maxOccPE>=0.2) || (navigationFlag==NAVI_SE_FLAG && secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>=0.15)))
					//if((averOccNumNaviOccQueue<2.5*minKmerOccSE) || (readsNumRatio>2 || readsNumRatio<0.3) || ((navigationFlag==NAVI_PE_FLAG && secondOccPE>=minKmerOccPE && secondOccPE/maxOccPE>=0.2) || (navigationFlag==NAVI_SE_FLAG && secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>=0.15)))
					if((readsNumRatio>3 || readsNumRatio<0.3) || ((navigationFlag==NAVI_PE_FLAG && secondOccPE>=minKmerOccPE && secondOccPE/maxOccPE>=0.2) || (navigationFlag==NAVI_SE_FLAG && secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>=0.15)))
					{
						number_of_overlap_less_than_threshold ++;

						printf("===localContigID=%ld, assemblyRound=%d, itemNumContigArr=%d, itemNumDecisionTable=%d\n", localContigID, assemblyRound, itemNumContigArr, itemNumDecisionTable);
						if(navigationFlag==NAVI_PE_FLAG)
							printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
						else if(navigationFlag==NAVI_SE_FLAG)
						{
							printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
							printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
						}
						printf("\tdistance=%d, readsNumRatio=%.2f, averOccNumNaviOccQueue=%.2f, lowOccNum=%d\n", itemNumContigArr-successContigIndex, readsNumRatio, averOccNumNaviOccQueue, lowOccNum);


						if(assemblyRound==SECOND_ROUND_ASSEMBLY && itemNumContigArr<CONTIG_LEN_THRESHOLD)
						{
							break;
						}

						// prepare for the second round assembly
						if(assemblyRound==FIRST_ROUND_ASSEMBLY)
						{

							assemblyRound ++;
							turnContigIndex = itemNumContigArr;

							int returnCode = initSecondAssembly();
							if(returnCode==FAILED)
							{
								break;
							}else if(returnCode==ERROR)
							{
								return FAILED;
							}

						}else
						{

							// ############################ Debug information ##############################
							//if(successContigIndex>0 && itemNumContigArr>=CONTIG_LEN_THRESHOLD)
							//{
							//	printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, successContigIndex<=0, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
							//	return FAILED;
							//}
							// ############################ Debug information ##############################

							break;
						}
					}
				}
*/
			}//end while(kmer)

			if(successContigIndex>0)
			{
				// ############################ Debug information ##############################
//				if(localContigID==98)
//				{
//					outputContigToTmpFile(contighead, HANGING_READ_TYPE_CONTIG_FILE);
//				}
				// ############################ Debug information ##############################

				//====================================================
				// trim a read length of contig nodes at tail
				//if(averKmerOcc>10 && itemNumContigArr>=3*readLen)  // deleted 2012-11-28
				if(trimReadLenFlag==YES)						// added 2012-11-28
				{
					if(trimContigTail(&successContigIndex, &itemNumContigArr, readLen, SECOND_ROUND_ASSEMBLY)==FAILED)
					{
						printf("line=%d, In %s(), localContigID=%ld, cannot trim contig nodes at contig tail by a read length, error!\n", __LINE__, __func__, localContigID);
						return ERROR;
					}
				}
				//====================================================


				if(updateContigtailnodes(contigArr, successContigIndex, &itemNumContigArr)==FAILED)  // deleted 2012-12-30
				//if(updateContigNodes(contigArr, validHeadRowContigArr, validTailRowContigArr, &successContigIndex, &itemNumContigArr)==FAILED)  // deleted 2012-12-30
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


				if(itemNumContigArr>=minContigLen)
				{
					contigsNum ++;

					// ############################ Debug information ##############################
					//printf("contigID=%d, contigLen=%d, turnContigIndex=%d.\n", contigsNum, itemNumContigArr, turnContigIndex);
					// ############################ Debug information ##############################

					successReadNum += this_successReadNum;

					// output contig nodes to file
					if(outputContigToFile(fpContigsBase, BASE_TYPE_FASTA_CONTIG_FILE, contigsNum, contigArr, itemNumContigArr)==FAILED)
					{
						printf("line=%d, In %s(), localContigID=%ld, cannot output contig nodes to file, error!\n", __LINE__, __func__, localContigID);
						return FAILED;
					}
					if(hangingContigOutFlag==YES)
					{
						if(outputContigToFile(fpContigsHanging, HANGING_READ_TYPE_CONTIG_FILE, contigsNum, contigArr, itemNumContigArr)==FAILED)
						{
							printf("line=%d, In %s(), localContigID=%ld, cannot output contig nodes to file, error!\n", __LINE__, __func__, localContigID);
							return FAILED;
						}
					}

					basesNum += itemNumContigArr;

					contigsLenArr[itemNumContigsLenArr] = itemNumContigArr;
					itemNumContigsLenArr ++;
					if(itemNumContigsLenArr>=maxItemNumContigsLenArr)
					{
						contigsLenArr = (int32_t *) realloc(contigsLenArr, 2*maxItemNumContigsLenArr*sizeof(int32_t));
						if(contigsLenArr==NULL)
						{
							printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}
				}else
				{
					successReadNum -= this_successReadNum;
				}
			}else
			{
				// ############################ Debug information ##############################
				//if(itemNumContigArr>=CONTIG_LEN_THRESHOLD)
				//{
				//	printf("line=%d, In %s(), localContigID=%ld, contigsNum=%d, assemblyRound=%d, itemNumContigArr=%ld, successContigIndex<=0, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
				//	return FAILED;
				//}
				// ############################ Debug information ##############################
			}

			// clean the PE hash table
			if(PEGivenType>NONE_PE_GIVEN_TYPE && readsNumInPEHashArr>0)
			{
				if(cleanReadsFromPEHashtable()==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, cannot clean reads from PE hash table, error!\n", __LINE__, __func__, localContigID);
					return FAILED;
				}
			}

			// clean the DTRow hash table
			if(cleanDTRowIndexHashtable(dtRowHashtable)==FAILED)
			{
				printf("line=%d, In %s(), contigsNum=%d, cannot clean DTRow hash table, Error!\n", __LINE__, __func__, contigsNum+1);
				return ERROR;
			}

			// clean the contig array
			cleanContigArray(contigArr, &itemNumContigArr);

			// update the assembly percentage
			if((int)((double) successReadNum / validReadNum * 100) > percent)
			{
				percentNum ++;
				percent = (int)((double) successReadNum / validReadNum * 100);

#if(DEBUG_OUTPUT==YES)
			printf("%d%%\n", percent);
#else
			printf("%d%%", percent);
			if(percentNum%10==0)
				printf("\n");
			else
				printf("\t");
#endif

				fflush(stdout);
			}

		} //end while(kmerIndex < TABLE_SIZE_DE_BRUIJN)
	}

	if(percent!=100)
	{
		percent = 100;
		printf("%d%%\n", percent);
		fflush(stdout);
		percentNum ++;
	}


	fclose(fpContigsBase);
	fpContigsBase = NULL;
	if(hangingContigOutFlag==YES)
	{
		fclose(fpContigsHanging);
		fpContigsHanging = NULL;
	}

	if(errorCorrectionFlag==YES)
	{
		fclose(fpReadCorrected);
		fpReadCorrected = NULL;
	}

	// get the statistics of contig lengths
	if(contigsLenStatistics(contigsLenArr, itemNumContigsLenArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute contig length statistics, error!\n", __LINE__, __func__);
		return FAILED;
	}

	freeMemory();

	printf("contigsNum=%d, basesNum=%ld\n", contigsNum, basesNum);

#if (DEBUG_PARA_PRINT==YES)
	if(errorCorrectionFlag==YES)
		printf("totalReadsNumCorreted=%ld\n", totalReadsNumCorreted);
	printf("successReadNum=%ld\n", successReadNum);
	printf("totalAlignedSuccessReadNum=%d\n", totalAlignedSuccessReadNum);
	printf("totalErrReadNum=%d\n", totalErrReadNum);
	for(i=0; i<11; i++) printf("errNumArr[%d]=%d\n", i, errNumArr[i]);
#endif

	gettimeofday(&tp_end,NULL);
	time_used = tp_end.tv_sec-tp_start.tv_sec + (double)(tp_end.tv_usec-tp_start.tv_usec)/1000000;
	printf("Assembly Used Time: %.2f Seconds.\n", time_used);

	printf("============= End build contigs. =============\n");

	return SUCCESSFUL;
}


/**
 * Initialize the memory for assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemory(char *(*occPointFileArray)[6])
{
	int32_t i;

	//longKmerSize = ceil((readLen - 2*errorRegLenEnd3 - kmerSize) * LONG_KMER_SIZE_FACTOR) + kmerSize;
	//longKmerSize = ceil((readLen - 1.5*errorRegLenEnd3 - kmerSize) * LONG_KMER_SIZE_FACTOR) + kmerSize; //=====================
	//longKmerSize = ceil((readLen - errorRegLenEnd3 - kmerSize) * LONG_KMER_SIZE_FACTOR) + kmerSize;
	//longKmerSize = ceil((readLen - errorRegLenEnd3) * 0.9);

	longKmerSize = ceil(readLen * LONG_KMER_SIZE_FACTOR);

	if((longKmerSize & 1) == 0)
		longKmerSize --;
	longKmerStepSize = floor((longKmerSize - kmerSize) / 5.0);
	if((longKmerStepSize & 1) == 1)
		longKmerStepSize ++;
	if(longKmerStepSize<1)
		longKmerStepSize = 2;

	if(PEGivenType>=NONE_PE_GIVEN_TYPE)
	{
		minKmerOccPE = ceil(averKmerOcc * MIN_KMER_OCC_FACTOR);
		minKmerOccSE = ceil(averKmerOcc * MIN_KMER_OCC_FACTOR);

		if(minKmerOccPE<MIN_KMER_OCC_THRES)
			minKmerOccPE = MIN_KMER_OCC_THRES;
		//else if(minKmerOccPE>MAX_KMER_OCC_THRES)
		//	minKmerOccPE = MAX_KMER_OCC_THRES;
		if(minKmerOccSE<MIN_KMER_OCC_THRES)
			minKmerOccSE = MIN_KMER_OCC_THRES;
		//else if(minKmerOccSE>MAX_KMER_OCC_THRES)
		//	minKmerOccSE = MAX_KMER_OCC_THRES;

		//maxSecondOcc = minKmerOccSE * OCCS_NUM_FACTOR;
		//maxSecondOcc = ceil(minKmerOccSE * MAX_SECOND_OCC_FACTOR);
		maxSecondOcc = averKmerOcc;
		//maxFirstOcc = ceil(minKmerOccSE * MAX_FIRST_OCC_FACTOR);
		//minLongKmerOcc = floor(minKmerOccSE * LONG_KMER_OCC_FACTOR);
		minLongKmerOcc = floor(averKmerOcc * 0.5);
		minReadsNumPEHashThres = ceil(averKmerOcc * MIN_READ_NUM_PE_HASH_FACTOR);

		//if(maxSecondOcc>MAX_SECOND_OCC_THRES)
		//{
		//	maxSecondOcc = MAX_SECOND_OCC_THRES;
		//}

		if(minLongKmerOcc<minKmerOccSE)
			minLongKmerOcc = minKmerOccSE;
		if(minLongKmerOcc>MIN_LONG_KMER_OCC_THRES)
		{
			minLongKmerOcc = MIN_LONG_KMER_OCC_THRES;
		}


		maxOccNumFaiedPE = ceil(OCCS_NUM_SE_FAILED_PE_FACTOR * averKmerOcc);
		if(maxOccNumFaiedPE>MAX_OCC_NUM_FAILED_PE_THRES)
			maxOccNumFaiedPE = MAX_OCC_NUM_FAILED_PE_THRES;
//		else
//			maxOccNumFaiedPE *= 0.8;
		//maxNavigationNumSE = MAX_NAVI_NUM_SE_THRES;

#if(DEBUG_PARA_PRINT==YES)
		printf("averKmerOcc=%.2f\n", averKmerOcc);
		printf("longKmerSize=%d, longKmerStepSize=%d\n", longKmerSize, longKmerStepSize);
		printf("minKmerOccSE=%.2f, minKmerOccPE=%.2f\n", minKmerOccSE, minKmerOccPE);
		printf("maxSecondOcc=%.2f\n", maxSecondOcc);
		//printf("maxFirstOcc=%.2f\n", maxFirstOcc);
		printf("minLongKmerOcc=%.2f\n", minLongKmerOcc);
		printf("minReadsNumPEHashThres=%.2f\n", minReadsNumPEHashThres);
		printf("maxOccNumFaiedPE=%.2f\n", maxOccNumFaiedPE);
		//printf("maxNavigationNumSE=%d\n", maxNavigationNumSE);
#endif
	}else
	{
		minKmerOccSE = ceil(averKmerOcc * MIN_KMER_OCC_FACTOR);

		if(minKmerOccSE<MIN_KMER_OCC_THRES)
			minKmerOccSE = MIN_KMER_OCC_THRES;

		//maxSecondOcc = minKmerOccSE * OCCS_NUM_FACTOR;
		//maxSecondOcc = ceil(minKmerOccSE * MAX_SECOND_OCC_FACTOR);
		maxSecondOcc = averKmerOcc;
		//maxFirstOcc = ceil(minKmerOccSE * MAX_FIRST_OCC_FACTOR);
		//minLongKmerOcc = floor(minKmerOccSE * LONG_KMER_OCC_FACTOR);
		minLongKmerOcc = floor(averKmerOcc * 0.5);


		//if(maxSecondOcc>MAX_SECOND_OCC_THRES)
		//{
		//	maxSecondOcc = MAX_SECOND_OCC_THRES;
		//}

		if(minLongKmerOcc<minKmerOccSE)
			minLongKmerOcc = minKmerOccSE;
		if(minLongKmerOcc>MIN_LONG_KMER_OCC_THRES)
		{
			minLongKmerOcc = MIN_LONG_KMER_OCC_THRES;
		}

#if(DEBUG_PARA_PRINT==YES)
		printf("averKmerOcc=%.2f\n", averKmerOcc);
		printf("longKmerSize=%d, longKmerStepSize=%d\n", longKmerSize, longKmerStepSize);
		printf("minKmerOccSE=%.2f\n", minKmerOccSE);
		printf("maxSecondOcc=%.2f\n", maxSecondOcc);
		//printf("maxFirstOcc=%.2f\n", maxFirstOcc);
		printf("minLongKmerOcc=%.2f\n", minLongKmerOcc);
#endif
	}

	lockedReadsNumThres = averKmerOcc;

	// the global variables of reads number region
	maxRegLenReadsNumReg = ceil(readLen * REG_LEN_READS_NUM_REG_FACTOR);
	minContigLenCheckingReadsNum = readLen + maxRegLenReadsNumReg;
	maxReadsNumRatioThres = MAX_READS_NUM_RATIO_THRES;
	minReadsNumRatioThres = MIN_READS_NUM_RATIO_THRES;
	solvedRepeatsNum = 0;

	maxUnmatchBaseNumPerRead = MAX_UNMATCH_BASE_NUM_PER_READ;
	minUnmatchBaseNumAlign = MIN_UNMATCH_BASE_NUM_ALIGN_FACTOR * maxUnmatchBaseNumPerRead;
	maxUnmatchBaseNumAfterAlign = MAX_UNMATCH_BASE_NUM_AFTER_ALIGN_FACTOR * maxUnmatchBaseNumPerRead;
	maxErrBaseNumInCorrection = MAX_ERR_NUM_IN_CORRECTION;
	minSuccessiveAppearedBaseNum = MIN_SUCCESSIVE_APPEARED_BASE_NUM;
	maxSuccessiveAppearedBaseNum = kmerSize;

	maxSeqLenAlign = deBruijnGraph->readSet->maxReadLen;
	matchScore = MATCH_SCORE;
	mismatchScore = MISMATCH_SCORE;
	gapScore = GAP_SCORE;
	refPosContig = -1;

	lowOccThresNaviOccQueue = 0.4 * averKmerOcc;

#if(DEBUG_PARA_PRINT==YES)
	printf("lockedReadsNumThres=%.2f\n", lockedReadsNumThres);
	printf("maxRegLenReadsNumReg=%d\n", maxRegLenReadsNumReg);
	printf("minContigLenCheckingReadsNum=%d\n", minContigLenCheckingReadsNum);
	printf("maxReadsNumRatioThres=%.2f\n", maxReadsNumRatioThres);
	printf("minReadsNumRatioThres=%.2f\n", minReadsNumRatioThres);
	printf("maxUnmatchBaseNumPerRead=%d\n", maxUnmatchBaseNumPerRead);
	printf("minUnmatchBaseNumAlign=%d\n", minUnmatchBaseNumAlign);
	printf("maxUnmatchBaseNumAfterAlign=%d\n", maxUnmatchBaseNumAfterAlign);
	printf("maxSuccessiveAppearedBaseNum=%d\n", maxSuccessiveAppearedBaseNum);
	printf("minSuccessiveAppearedBaseNum=%d\n", minSuccessiveAppearedBaseNum);
	//printf("lowOccThresNaviOccQueue=%.2f\n", lowOccThresNaviOccQueue);
#endif

	hangingContigOutFlag = HANGING_CONTIG_OUT_FLAG;
	maxItemNumDecisionTable = TABLE_SIZE_ASSEMBLINGREAD;
	// decision table
	decisionTable = (assemblingreadtype*) malloc(maxItemNumDecisionTable*sizeof(assemblingreadtype));
	if(decisionTable==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the decision table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	itemNumContigArr = 0;
	maxItemNumContigArr = ITEM_NUM_CONTIG_ARR;
	contigArr = (contigtype *) malloc (maxItemNumContigArr * sizeof(contigtype));
	if(contigArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	dtRowHashtable = (dtRowIndex_t **) calloc(TABLE_SIZE_HASH_DTROWINDEX, sizeof(dtRowIndex_t*));
	if(dtRowHashtable==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the decision table row index hash table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	maxItemNumSuccessReadsArr = TABLE_SIZE_RIDPOSORIENTATION;
	successReadsArr = (successRead_t*) malloc(maxItemNumSuccessReadsArr*sizeof(successRead_t));
	if(successReadsArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	countingHangingBucketArr = (int32_t *) calloc (ITEM_NUM_COUNTHANGING_BUCKETS, sizeof(int32_t));
	if(countingHangingBucketArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//++++++++++++++++++++++++++++++++++++
	PEHashArr = (PERead_t **) calloc(TABLE_SIZE_HASH_PE, sizeof(PERead_t **));
	if(PEHashArr==NULL)
	{
		printf("In %s(), cannot allocate memory, error!\n", __func__);
		return FAILED;
	}

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

	//maxItemNumNaviOccQueue = MAX_ITEM_NUM_NAVI_OCC_QUEUE;
	maxItemNumNaviOccQueue = readLen * 0.5;
	naviOccQueue = (double*) malloc(maxItemNumNaviOccQueue*sizeof(double));
	if(naviOccQueue==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// alignments
	readSeqAlign = (char *) malloc (2*(maxSeqLenAlign+1) *sizeof(char));
	if(readSeqAlign==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	contigSeqAlign = (char *) malloc (2*(maxSeqLenAlign+1) *sizeof(char));
	if(contigSeqAlign==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<3; i++)
	{
		alignResultSeqArr[i] = (char *) malloc (2*(maxSeqLenAlign+1) *sizeof(char));
		if(alignResultSeqArr[i]==NULL)
		{
			printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	matchFlagArr = (char *) malloc (2*(maxSeqLenAlign+1) *sizeof(char));
	if(matchFlagArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	alignScoreArr = (int32_t *) malloc (2*(maxSeqLenAlign+1)*2*(maxSeqLenAlign+1) *sizeof(int32_t));
	if(alignScoreArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	readseqCorrected = (uint64_t *) malloc (((MAX_READ_LEN_IN_BUF>>5)+1) *sizeof(uint64_t));
	if(readseqCorrected==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	itemNumContigsLenArr = 0;
	maxItemNumContigsLenArr = MAX_ITEM_NUM_CONTIGS_LEN_ARR;
	contigsLenArr = (int32_t *) malloc (maxItemNumContigsLenArr*sizeof(int32_t));
	if(contigsLenArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//==============================================
	// multi-align variables
	maxItemNumReadsMalignArr = MAX_READS_MALIGN;
	readseqMalignArr = (char **) calloc(maxItemNumReadsMalignArr, sizeof(char*));
	if(readseqMalignArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	alignResultsMalign = (char **) calloc(maxItemNumReadsMalignArr, sizeof(char*));
	if(alignResultsMalign==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	for(i=0; i<maxItemNumReadsMalignArr; i++)
	{
		readseqMalignArr[i] = (char *) calloc(MAX_SEQ_LEN_MALIGN+1, sizeof(char));
		if(readseqMalignArr[i]==NULL)
		{
			printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		alignResultsMalign[i] = (char *) calloc(MAX_SEQ_LEN_MALIGN+1, sizeof(char));
		if(alignResultsMalign[i]==NULL)
		{
			printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}
	readIDMalignArr = (int64_t *) calloc(maxItemNumReadsMalignArr, sizeof(int64_t));
	if(readIDMalignArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}


	for(i=0; i<4; i++) tabooSeqInt[i] = 0;
	for(i=0; i<32; i++)
	{
		tabooSeqInt[1] = (tabooSeqInt[1] << 2) | 1;
		tabooSeqInt[2] = (tabooSeqInt[2] << 2) | 2;
		tabooSeqInt[3] = (tabooSeqInt[3] <<  2) | 3;
	}


	//===================== draw curve variable begin ======================
#if(DRAW_CURVE_FLAG==YES)
	// load the reference
	strcpy(refFile, "../model/NC_000913.fasta");
	if(loadReferenceFromFile(refFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot load reference, error!\n", __LINE__, __func__);
		return FAILED;
	}

	fpOccPoint = fopen(occPointFileArray[0][0], "w");
	if(fpOccPoint==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, occPointFileArray[0][0]);
		return FAILED;
	}
	for(i=0; i<4; i++)
	{
		fpOccExtensionCorrect[i] = fopen(occPointFileArray[i][1], "w");
		if(fpOccExtensionCorrect==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, occPointFileArray[i][1]);
			return FAILED;
		}
		fpOccExtensionIncorrect[i] = fopen(occPointFileArray[i][2], "w");
		if(fpOccExtensionIncorrect==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, occPointFileArray[i][2]);
			return FAILED;
		}
		//=====
		fpOccStopCorrect[i] = fopen(occPointFileArray[i][3], "w");
		if(fpOccStopCorrect==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, occPointFileArray[i][3]);
			return FAILED;
		}
		fpOccStopIncorrect[i] = fopen(occPointFileArray[i][4], "w");
		if(fpOccStopIncorrect==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, occPointFileArray[i][4]);
			return FAILED;
		}
	}
	fpRefPos = fopen(occPointFileArray[0][5], "w");
	if(fpRefPos==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, occPointFileArray[0][5]);
		return FAILED;
	}

#endif
	//===================== draw curve variable end ======================

	return SUCCESSFUL;
}

void freeMemory()
{
	int32_t i;

	//maxItemNumDecisionTable = 0;
	//maxItemNumSuccessReadsArr = 0;

	free(decisionTable);
	decisionTable = NULL;

	free(contigArr);
	contigArr = NULL;

	free(dtRowHashtable);
	dtRowHashtable = NULL;

	free(successReadsArr);
	successReadsArr = NULL;

	free(countingHangingBucketArr);
	countingHangingBucketArr = NULL;

	free(PEHashArr);
	PEHashArr = NULL;

	free(kmerSeqIntAssembly);
	kmerSeqIntAssembly = NULL;

	free(kmerSeqIntAssemblyRev);
	kmerSeqIntAssemblyRev = NULL;

	free(tmpKmerSeqIntAssembly);
	tmpKmerSeqIntAssembly = NULL;

	free(naviOccQueue);
	naviOccQueue = NULL;

	free(readSeqAlign);
	readSeqAlign = NULL;

	free(contigSeqAlign);
	contigSeqAlign = NULL;

	for(i=0; i<3; i++)
	{
		free(alignResultSeqArr[i]);
		alignResultSeqArr[i] = NULL;
	}

	free(matchFlagArr);
	matchFlagArr = NULL;

	free(alignScoreArr);
	alignScoreArr = NULL;

	free(readseqCorrected);
	readseqCorrected = NULL;

	free(contigsLenArr);
	contigsLenArr = NULL;

	for(i=0; i<maxItemNumReadsMalignArr; i++)
	{
		free(readseqMalignArr[i]);
		free(alignResultsMalign[i]);
	}
	free(readseqMalignArr);
	free(alignResultsMalign);
	free(readIDMalignArr);


#if(DRAW_CURVE_FLAG==YES)
	freeReference();
	fclose(fpOccPoint);
	for(i=0; i<4; i++)
	{
		fclose(fpOccExtensionCorrect[i]);
		fclose(fpOccExtensionIncorrect[i]);
		fclose(fpOccStopCorrect[i]);
		fclose(fpOccStopIncorrect[i]);
	}
	fclose(fpRefPos);
#endif

}

/**
 * Initialize the bounder for first k-mers.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initFirstKmerBounder(double *lowerBoundFirstKmer, double *upperBoundFirstKmer, short assemblyCycle, double averKmerOccNum)
{
	if(assemblyCycle==1)
	{ // the first assembly cycle
		*lowerBoundFirstKmer = averKmerOccNum * LOWER_BOUND_FACTOR_CYCLE1;
		*upperBoundFirstKmer = averKmerOccNum * UPPER_BOUND_FACTOR_CYCLE1;

		if((*lowerBoundFirstKmer)>MIN_LOWER_BOUND)
			*lowerBoundFirstKmer = MIN_LOWER_BOUND;
		else if((*lowerBoundFirstKmer)<MIN_LOWER_BOUND_RESCUE)
			*lowerBoundFirstKmer = MIN_LOWER_BOUND_RESCUE;
	}
	else if(assemblyCycle==2)
	{ // the second assembly cycle
		*lowerBoundFirstKmer = averKmerOccNum * UPPER_BOUND_FACTOR_CYCLE1;
		*upperBoundFirstKmer = averKmerOccNum * UPPER_BOUND_FACTOR_CYCLE2;
	}
//	else if(assemblyCycle==3)
//	{
//		*lowerBoundFirstKmer = LOWER_BOUND_CYCLE3;
//		*upperBoundFirstKmer = averKmerOccNum * LOWER_BOUND_FACTOR_CYCLE1;
//	}
	else
	{
		printf("line=%d, In %s(), assemblyCycle=%d, error!\n", __LINE__, __func__, assemblyCycle);
		return FAILED;
	}

#if(DEBUG_PARA_PRINT==YES)
	printf("assemblyCycle=%d: lowerBoundFirstKmer=%.2f, upperBoundFirstKmer=%.2f, averKmerOcc=%.2f\n", assemblyCycle, *lowerBoundFirstKmer, *upperBoundFirstKmer, averKmerOccNum);
#endif

	return SUCCESSFUL;
}

/**
 * Initialize the contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initContig()
{
	int32_t i, maxNumPlus, ridMaxNumPlus, maxNumMinus, ridMaxNumMinus;
	int32_t posNum, rpos, seqLen, basePos, entriesNumTmp, baseNumLastEntryTmp, baseInt, entryRow, entryPos;
	ridpostype *ridpostable;

	int32_t readBlockID, rowNumInReadBlock, maxItemNumPerReadBlock, errorRegLenEnd3; // block id starts from 0

	readBlock_t *readBlockArr;
	read_t *pRead;
	readseqBlock_t *readseqBlockArr;
	uint64_t *pReadseq;

	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;
	readseqBlockArr = deBruijnGraph->readSet->readseqBlockArr;

	maxNumPlus = 0;
	ridMaxNumPlus = -1;
	maxNumMinus = 0;
	ridMaxNumMinus = -1;

	// get the end omitted bases
	if(kmers[0])
	{
		ridpostable = kmers[0]->ppos;
		posNum = kmers[0]->arraysize;
		for(i=0; i<posNum; i++)
		{
			rpos = ridpostable[i].pos;
			if(rpos==1 && ridpostable[i].delsign==0)
			{
				if(rpos-1>maxNumPlus)
				{
					maxNumPlus = rpos - 1;
					ridMaxNumPlus = ridpostable[i].rid;
				}
			}
		}
	}

	if(kmers[1])
	{
		ridpostable = kmers[1]->ppos;
		posNum = kmers[1]->arraysize;
		for(i=0; i<posNum; i++)
		{
			rpos = ridpostable[i].pos;

			readBlockID = (ridpostable[i].rid - 1) / maxItemNumPerReadBlock;
			rowNumInReadBlock = (ridpostable[i].rid - 1) % maxItemNumPerReadBlock;
			pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;
			seqLen = pRead->seqlen;
			errorRegLenEnd3 = ceil(seqLen * ERROR_REGION_LEN_3End_FACTOR);

			//printf("rid=%ld, rpos=%d\n", (int64_t)ridpostable[i].rid, rpos);

			if(rpos>seqLen-kmerSize+1-errorRegLenEnd3 && ridpostable[i].delsign==0)
			{
				if(seqLen-(rpos+kmerSize-1) > maxNumMinus)
				{
					maxNumMinus = seqLen - (rpos + kmerSize - 1);
					ridMaxNumMinus = ridpostable[i].rid;
				}
			}
		}
	}

	if(maxNumPlus>0 || maxNumMinus>0)
	{
		if(maxNumPlus>=maxNumMinus)
		{
			readBlockID = (ridMaxNumPlus - 1) / maxItemNumPerReadBlock;
			rowNumInReadBlock = (ridMaxNumPlus - 1) % maxItemNumPerReadBlock;
			pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;
			pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
			seqLen = pRead->seqlen;

			entriesNumTmp = ((seqLen - 1) >> 5) + 1;
			baseNumLastEntryTmp = ((seqLen - 1) % 32) + 1;

			for(i=0; i<maxNumPlus; i++)
			{
				entryRow = i >> 5;
				entryPos = i % 32;

				contigArr[i].index = i + 1;
				if(entryRow==entriesNumTmp-1)
					contigArr[i].base = (pReadseq[entryRow] >> (2*(baseNumLastEntryTmp-entryPos-1))) & 3;
				else
					contigArr[i].base = (pReadseq[entryRow] >> (64 - 2*(entryPos+1))) & 3;
				contigArr[i].ridposnum = 0;
				contigArr[i].pridposorientation = NULL;
			}

			itemNumContigArr += maxNumPlus;
		}else
		{
			readBlockID = (ridMaxNumMinus - 1) / maxItemNumPerReadBlock;
			rowNumInReadBlock = (ridMaxNumMinus - 1) % maxItemNumPerReadBlock;
			pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;
			pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
			seqLen = pRead->seqlen;

			entriesNumTmp = ((seqLen - 1) >> 5) + 1;
			baseNumLastEntryTmp = ((seqLen - 1) % 32) + 1;

			basePos = seqLen - 1;
			for(i=0; i<maxNumMinus; i++, basePos--)
			{
				entryRow = basePos >> 5;
				entryPos = basePos % 32;

				if(entryRow==entriesNumTmp-1)
					baseInt = (pReadseq[entryRow] >> (2*(baseNumLastEntryTmp-entryPos-1))) & 3;
				else
					baseInt = (pReadseq[entryRow] >> (64 - 2*(entryPos+1))) & 3;

				contigArr[i].index = i + 1;
				contigArr[i].base = (~baseInt) & 3;
				contigArr[i].ridposnum = 0;
				contigArr[i].pridposorientation = NULL;
			}

			itemNumContigArr += maxNumMinus;
		}
	}


	// add the k-mer bases
	for(i=0; i<kmerSize; i++, itemNumContigArr++)
	{
		entryRow = i >> 5;
		entryPos = i % 32;

		contigArr[itemNumContigArr].index = itemNumContigArr + 1;
		if(entryRow==entriesPerKmer-1)
			contigArr[itemNumContigArr].base = (kmerSeqIntAssembly[entryRow] >> (2*(lastEntryBaseNum-entryPos-1))) & 3;
		else
			contigArr[itemNumContigArr].base = (kmerSeqIntAssembly[entryRow] >> (64 - 2*(entryPos+1))) & 3;
		contigArr[itemNumContigArr].ridposnum = 0;
		contigArr[itemNumContigArr].pridposorientation = NULL;
	}

	return SUCCESSFUL;
}

/**
 * set first kmer threshold.
 */
short initFirstKmerThreshold()
{
	double kmerOccTmp;
	int32_t endKmerNumTmp, kmerNumTmp;

	kmerOccTmp = (double)deBruijnGraph->totalItemNumRidpos / deBruijnGraph->totalItemNumKmer;

	kmerNumTmp = averReadLenInFileSample - kmerSize + 1;
	endKmerNumTmp = ceil(averReadLenInFileSample * kmerRegLenRatioEnd5) + ceil (averReadLenInFileSample * kmerRegLenRatioEnd3);

	averKmerOcc = kmerOccTmp;
	averKmerOcc *= (double)kmerNumTmp / endKmerNumTmp;
	averKmerOcc *= (double)readLen / averReadLenInFileSample;
	averKmerOcc /= 1.4;

	//averKmerOcc = (double)deBruijnGraph->totalItemNumRidpos / deBruijnGraph->totalItemNumKmer;
	if(averKmerOcc<1)
		averKmerOcc = 1;

	firstKmerThres = averKmerOcc;

#if(DEBUG_PARA_PRINT==YES)
	printf("old averKmerOcc=%.2f, new averKmerOcc=%.2f\n", kmerOccTmp, averKmerOcc);
#endif

	return SUCCESSFUL;
}


/**
 * Initialize the first kmers.
 */
short getFirstKmers(uint64_t *kmerIndex, kmertype **firstKmer)
{
	kmertype *kmer;
	uint64_t i, *seqInt, kmerBlockID, itemRowKmerBlock;


	i = *kmerIndex;
	kmer = *firstKmer;
	naviSuccessFlag = NAVI_FAILED;
	kmers[0] = kmers[1] = NULL;
	while(i<hashTableSize)
	{
		if(kmer==NULL)
		{
			kmerBlockID = deBruijnGraph->kmerHashtable[i].kmerBlockID;
			itemRowKmerBlock = deBruijnGraph->kmerHashtable[i].itemRowKmerBlock;
			if(kmerBlockID>0)
				kmer = deBruijnGraph->kmerBlockArr[kmerBlockID-1].kmerArr + itemRowKmerBlock;
			else
				kmer = NULL;
		}
		else
		{
			kmerBlockID = kmer->nextKmerBlockID;
			itemRowKmerBlock = kmer->nextItemRowKmerBlock;
			if(kmerBlockID>0)
				kmer = deBruijnGraph->kmerBlockArr[kmerBlockID-1].kmerArr + itemRowKmerBlock;
			else
				kmer = NULL;
		}

		while(kmer)
		{
			//########################## Debug information #############################
//			if(kmer->multiplicity>=firstKmerThres && kmer->multiplicity<=FIRSTKMER_FACTOR*firstKmerThres)
//				printf("########### i=%lu, multiplicity=%u\n", i, kmer->multiplicity);
			//########################## Debug information #############################

			seqInt = deBruijnGraph->kmerSeqBlockArr[kmer->kmerseqBlockID-1].kmerSeqArr + kmer->itemRowKmerseqBlock * deBruijnGraph->entriesPerKmer;
			if(isValidFirstKmer(kmer, seqInt)==YES)
			{
				//########################## Debug information #############################
//				printf("########### i=%lu\n",i);
				//########################## Debug information #############################

				//if(kmer->multiplicity>=firstKmerThres && kmer->multiplicity>=FIRSTKMER_SUBTRACT_THRESHOLD*kmer->arraysize)
				//if(kmer->multiplicity>=firstKmerThres && kmer->multiplicity<=FIRSTKMER_FACTOR*firstKmerThres && kmer->multiplicity>=FIRSTKMER_SUBTRACT_THRESHOLD*kmer->arraysize)  // deleted 2012-11-28
				if(kmer->multiplicity>=lowerBoundFirstKmer && kmer->multiplicity<=upperBoundFirstKmer && kmer->multiplicity>=FIRSTKMER_SUBTRACT_THRESHOLD*kmer->arraysize)  // added 2012-11-28
				{
					if(containFirstPos(kmer)==YES)
					{
						kmers[0] = kmer;

						if(memcpy(kmerSeqIntAssembly, seqInt, deBruijnGraph->bytesPerKmerseq)==NULL)
						{
							printf("line=%d, In %s(), can not memcpy the kmerseq. error!\n", __LINE__, __func__);
							return FAILED;
						}
						kmers[1] = getReverseKmer(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, deBruijnGraph);

						*kmerIndex = i;
						*firstKmer = kmer;

						naviSuccessFlag = NAVI_SUCCESS;

						return SUCCESSFUL;
					}
					else if(containLastPos(kmer)==YES)
					{
						kmers[1] = kmer;

						if(memcpy(kmerSeqIntAssemblyRev, seqInt, deBruijnGraph->bytesPerKmerseq)==NULL)
						{
							printf("line=%d, In %s(), can not memcpy the kmerseq. error!\n", __LINE__, __func__);
							return FAILED;
						}
						kmers[0] = getReverseKmer(kmerSeqIntAssembly, kmerSeqIntAssemblyRev, deBruijnGraph);

						*kmerIndex = i;
						*firstKmer = kmer;

						naviSuccessFlag = NAVI_SUCCESS;

						return SUCCESSFUL;
					}
				}
			}

			kmerBlockID = kmer->nextKmerBlockID;
			itemRowKmerBlock = kmer->nextItemRowKmerBlock;
			if(kmerBlockID>0)
				kmer = deBruijnGraph->kmerBlockArr[kmerBlockID-1].kmerArr + itemRowKmerBlock;
			else
				kmer = NULL;
		}
		i++;
	}
	if(i>=hashTableSize)
	{
		*kmerIndex = i;
		*firstKmer = NULL;
	}

	return SUCCESSFUL;
}

/**
 * Whether a kmer is a valid first kmer.
 *  @return:
 *  	If valid, return YES; otherwise return NO.
 *  	If errors, return ERROR.
 */
short isValidFirstKmer(kmertype *kmer, uint64_t *seqInt)
{
	int i, j, tabooFlag;
	uint64_t tmpTabooSeq;

	if(kmer==NULL)
		return NO;

	if(seqInt==NULL)
	{
		printf("line=%d, In %s(), seqInt==NULL. Error.\n", __LINE__, __func__);
		return ERROR;
	}

	tabooFlag = NO;
	if(entriesPerKmer >= 2)
	{
		for(j=0; j<4; j++)
		{
			if(seqInt[0] == tabooSeqInt[j])
			{
				tabooFlag = YES;
				tmpTabooSeq = seqInt[0];
				break;
			}
		}

		if(tabooFlag==YES)
		{
			for(i=1; i<entriesPerKmer-1; i++)
			{
				if(seqInt[i]!=tmpTabooSeq)
					return YES;
			}

			if(seqInt[entriesPerKmer-1] == (tmpTabooSeq & lastEntryMask))
				return NO;
		}else
			return YES;
	}else
	{
		for(j=0; j<4; j++)
		{
			if(seqInt[0] == (tabooSeqInt[j] & lastEntryMask))
				return NO;
		}
	}

	return YES;
}

/**
 * Check whether the k-mer contain the first read position.
 *  @return:
 *  	If it contains the first position, return YES; otherwise, return NO.
 */
short containFirstPos(kmertype *kmer)
{
	int i, posNum;
	ridpostype *ridpostable;

	if(kmer->multiplicity==0)
		return NO;

	posNum = kmer->arraysize;
	ridpostable = kmer->ppos;
	for(i=0; i<posNum; i++)
	{
		if(ridpostable[i].delsign==0 && ridpostable[i].pos==1)
			return YES;
	}
	return NO;
}

/**
 * Check whether the k-mer contain the last read position.
 *  @return:
 *  	If it contains the last position, return YES; otherwise, return NO.
 */
short containLastPos(kmertype *kmer)
{
	int i, pos, posNum;
	ridpostype *ridpostable;
	readBlock_t *readBlockArr;
	read_t *pRead;

	int32_t rid, maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3;


	if(kmer->multiplicity==0)
		return NO;

	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;
	posNum = kmer->arraysize;
	ridpostable = kmer->ppos;
	for(i=0; i<posNum; i++)
	{
		pos = ridpostable[i].pos;

		// get seqLen and errorRegLenEnd3
		rid = ridpostable[i].rid;
		blockRowRead = (rid - 1) / maxItemNumPerReadBlock;
		itemRowBlockRead = (rid - 1) % maxItemNumPerReadBlock;
		pRead = readBlockArr[blockRowRead].readArr + itemRowBlockRead;
		seqLen = pRead->seqlen;
		errorRegLenEnd3 = ceil(seqLen * ERROR_REGION_LEN_3End_FACTOR);

		if(ridpostable[i].delsign==0 && ridpostable[i].pos==seqLen-kmerSize+1)
			return YES;
	}
	return NO;
}

/**
 * Add first kmer of a contig to decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addFirstKmerToDecisionTable(kmertype **kmers)
{
	ridpostype *ridpostable;
	uint32_t i, j, rpos, posNum, entriesNumTmp, baseNumLastEntryTmp, basePos, entryRow, entryPos, baseInt, seqLen;

	int32_t readBlockID, rowNumInReadBlock, maxItemNumPerReadBlock, errorRegLenEnd3; // block id starts from 0
	int32_t contigPosTmp;  // starts from 1

	readBlock_t *readBlockArr;
	read_t *pRead;
	readseqBlock_t *readseqBlockArr;
	uint64_t *pReadseq;

	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;
	readseqBlockArr = deBruijnGraph->readSet->readseqBlockArr;


	if(kmers[0])
	{
		ridpostable = kmers[0]->ppos;
		posNum = kmers[0]->arraysize;
		for(i=0; i<posNum; i++)
		{
			rpos = ridpostable[i].pos;
			if(rpos==1 && ridpostable[i].delsign==0)
			{
				// ############################ Debug information ##############################
				//if(ridpostable[i].rid==664454)
				//{
				//	printf("line=%d, In %s(), rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
				//}
				// ############################ Debug information ##############################

				readBlockID = (ridpostable[i].rid - 1) / maxItemNumPerReadBlock;
				rowNumInReadBlock = (ridpostable[i].rid - 1) % maxItemNumPerReadBlock;
				pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;
				pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
				seqLen = pRead->seqlen;


				if(addReadToDecisionTable(ridpostable[i].rid, rpos, ORIENTATION_PLUS, NO, seqLen, pReadseq)==FAILED)
				{
					printf("line=%d, In %s(), cannot add read %lu, to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
					return FAILED;
				}

				if(rpos>1)
				{
					entriesNumTmp = ((seqLen - 1) >> 5) + 1;
					baseNumLastEntryTmp = ((seqLen - 1) % 32) + 1;

					// get the adjusted rpos
					basePos = rpos - 2;
					//contigPosTmp = itemNumContigArr - readLen;
					contigPosTmp = itemNumContigArr - kmerSize - 1;
					if(contigPosTmp<0)
					{
						//contigPosTmp = 0;
						printf("line=%d, In %s(), contigPosTmp=%d, error!\n", __LINE__, __func__, contigPosTmp);
						return FAILED;
					}

					//for(j=0; basePos>=0 && itemNumContigArr-kmerSize-1-j>=0; j++, basePos--)
					for(; basePos>=0 && contigPosTmp>=0; contigPosTmp--, basePos--)
					{
						// match to contig
						entryRow = basePos >> 5;
						entryPos = basePos % 32;

						if(entryRow==entriesNumTmp-1)
							baseInt = (pReadseq[entryRow] >> (2*(baseNumLastEntryTmp-entryPos-1))) & 3;
						else
							baseInt = (pReadseq[entryRow] >> (62 - 2*entryPos)) & 3;

						//if(baseInt==contigArr[itemNumContigArr-kmerSize-1-j].base)
						if(baseInt==contigArr[contigPosTmp].base)
						{
							decisionTable[itemNumDecisionTable-1].firstBasePos --;
							decisionTable[itemNumDecisionTable-1].firstContigPos --;
							decisionTable[itemNumDecisionTable-1].matchBaseNum ++;
						}else
						{
							break;
						}
					}
				}
			}
		}
	}


	if(kmers[1])
	{
		ridpostable = kmers[1]->ppos;
		posNum = kmers[1]->arraysize;
		for(i=0; i<posNum; i++)
		{
			rpos = ridpostable[i].pos;
			readBlockID = (ridpostable[i].rid - 1) / maxItemNumPerReadBlock;
			rowNumInReadBlock = (ridpostable[i].rid - 1) % maxItemNumPerReadBlock;
			pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;
			pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
			seqLen = pRead->seqlen;
			errorRegLenEnd3 = ceil(seqLen * ERROR_REGION_LEN_3End_FACTOR);

			if(rpos>seqLen-kmerSize+1-errorRegLenEnd3 && ridpostable[i].delsign==0)
			{
				// ############################ Debug information ##############################
				//if(ridpostable[i].rid==664454)
				//{
				//	printf("line=%d, In %s(), rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
				//}
				// ############################ Debug information ##############################

				if(addReadToDecisionTable(ridpostable[i].rid, rpos, ORIENTATION_MINUS, NO, seqLen, pReadseq)==FAILED)
				{
					printf("line=%d, In %s(), cannot add read %lu, to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
					return FAILED;
				}

				if(rpos<seqLen-kmerSize+1)
				{
					entriesNumTmp = ((seqLen - 1) >> 5) + 1;
					baseNumLastEntryTmp = ((seqLen - 1) % 32) + 1;

					basePos = rpos + kmerSize - 1;
					//contigPosTmp = itemNumContigArr - readLen;
					contigPosTmp = itemNumContigArr - kmerSize - 1;
					if(contigPosTmp<0)
					{
						//contigPosTmp = 0;
						printf("line=%d, In %s(), contigPosTmp=%d, error!\n", __LINE__, __func__, contigPosTmp);
						return FAILED;
					}

					//for(j=0; basePos<seqLen && itemNumContigArr-kmerSize-1-j>=0; j++, basePos++)
					for(; basePos<seqLen && contigPosTmp>=0; contigPosTmp--, basePos++)
					{
						// match to contig
						entryRow = basePos >> 5;
						entryPos = basePos % 32;

						if(entryRow==entriesNumTmp-1)
							baseInt = (pReadseq[entryRow] >> (2*(baseNumLastEntryTmp-entryPos-1))) & 3;
						else
							baseInt = (pReadseq[entryRow] >> (62 - 2*entryPos)) & 3;

						baseInt = (~baseInt) & 3;

						if(baseInt==contigArr[contigPosTmp].base)
						{
							decisionTable[itemNumDecisionTable-1].firstBasePos ++;
							decisionTable[itemNumDecisionTable-1].firstContigPos --;
							decisionTable[itemNumDecisionTable-1].matchBaseNum ++;
						}else
						{
							break;
						}
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Add a read into decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addReadToDecisionTable(uint64_t rid, int32_t rpos, int32_t orientation, int32_t matedFlag, int32_t seqLen, uint64_t *readseq)
{
	readBlock_t *readBlockArr;
	readseqBlock_t *readseqBlockArr;
	uint64_t maxItemNumPerReadBlock;
	assemblingreadtype *this_assemblingRead;

	this_assemblingRead = decisionTable + itemNumDecisionTable;

	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;
	readseqBlockArr = deBruijnGraph->readSet->readseqBlockArr;

	// ####################### Debug information ##########################
	//if(rid==544250)
	//{
	//	printf("line=%d, In %s(), rid=%lu, rpos=%d, orient=%d, matedFlag=%d\n", __LINE__, __func__, rid, rpos, orientation, matedFlag);
	//}
	// ####################### Debug information ##########################

	this_assemblingRead->rid = rid;
	if(orientation==ORIENTATION_PLUS)
		this_assemblingRead->firstBasePos = rpos - 1;
	else
		this_assemblingRead->firstBasePos = rpos + kmerSize - 2;
	this_assemblingRead->firstContigPos = itemNumContigArr - kmerSize;
	this_assemblingRead->orientation = orientation;
	this_assemblingRead->status = ASSEMBLING_STATUS;
	//this_assemblingRead->kmerappeartimes = 1;
	//this_assemblingRead->kmerunappeartimes = 0;
	//this_assemblingRead->lastappearpos = rpos;
	//this_assemblingRead->lastpos = rpos;
	//this_assemblingRead->kmerunappearblocks = 0;
	this_assemblingRead->delsign = 0;
	this_assemblingRead->reserved = 0;
	this_assemblingRead->locked = 0;
	this_assemblingRead->matedFlag = matedFlag;

	this_assemblingRead->successiveAppearBases = kmerSize;
	this_assemblingRead->successiveUnappearBases = 0;
	this_assemblingRead->unappearBlocksNum = 0;

	this_assemblingRead->matchBaseNum = kmerSize;
	this_assemblingRead->unmatchBaseNum = 0;
	this_assemblingRead->alignNum = 0;
	if(orientation==ORIENTATION_PLUS)
		this_assemblingRead->basePos = rpos + kmerSize - 2;
	else
		this_assemblingRead->basePos = rpos - 1;
	this_assemblingRead->lastMatchedBasePos = this_assemblingRead->basePos;
	this_assemblingRead->seqlen = seqLen;
	this_assemblingRead->readseq = readseq;
	this_assemblingRead->entriesNumReadseq = ((seqLen - 1) >> 5) + 1;
	this_assemblingRead->baseNumLastEentryReadseq = ((seqLen - 1) % 32) + 1;
	this_assemblingRead->kmerNumEnd5 = ceil(seqLen * kmerRegLenRatioEnd5);
	this_assemblingRead->kmerNumEnd3 = ceil(seqLen * kmerRegLenRatioEnd3);

	this_assemblingRead->errBaseNum = 0;
	this_assemblingRead->newErrNumAfterCorrection = 0;
	this_assemblingRead->alignSuccessTimes = 0;
	this_assemblingRead->alignSuccessTimesOld = 0;

	// add the read to DTRow hash table
	if(addReadToDTRowHashtable(rid, itemNumDecisionTable, dtRowHashtable)==FAILED)
	{
		printf("line=%d, In %s(), cannot add read rid=%lu to dtRowHashtable, error!\n", __LINE__, __func__, rid);
		return FAILED;
	}

	itemNumDecisionTable ++;

	// check the allowed number of reads in the array
	if(itemNumDecisionTable==maxItemNumDecisionTable)
	{
		if(reallocateDecisionTable()==FAILED)
		{
			printf("line=%d, In %s(), cannot reallocate memory for the reads in decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

    return SUCCESSFUL;
}

/**
 * Reallocate the memory of decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short reallocateDecisionTable()
{
	assemblingreadtype *tmp_pDecisionTable;
	tmp_pDecisionTable = (assemblingreadtype*) malloc (2*maxItemNumDecisionTable * sizeof(assemblingreadtype));
	if(tmp_pDecisionTable==NULL)
	{
		printf("In %s(), cannot reallocate memory for the reads in decision table, error!\n", __func__);
		return FAILED;
	}

	// copy memory
	if(memcpy(tmp_pDecisionTable, decisionTable, maxItemNumDecisionTable * sizeof(assemblingreadtype))==NULL)
	{
		printf("In %s(), cannot copy memory for the reads in decision table, error!\n", __func__);
		return FAILED;
	}

	free(decisionTable);
	decisionTable = tmp_pDecisionTable;
	maxItemNumDecisionTable *= 2;

	return SUCCESSFUL;
}

/**
 * Get the next k-mers for navigation combining the k-mer hash table and the decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNextKmerBySE(int contigNodesNum)
{
	uint64_t tmp_kmerseq[entriesPerKmer], tmp_kmerseqRev[entriesPerKmer];
	kmertype *tmp_kmers[4][2];
	int32_t i, j, validKmerNum, iterateID, successiveAppearBaseNum, ignoreEndBaseFlag, base_index;

	kmer_len = 0;

//	if(itemNumDecisionTable>MAX_DECISION_TABLE_SIZE_HTRES)
//	{
//		naviSuccessFlag = NAVI_FAILED;
//		kmers[0] = kmers[1] = NULL;
//		return SUCCESSFUL;
//	}

	maxOccIndexSE = -1;
	maxOccSE = 0;
	secondOccIndexSE = -1;
	secondOccSE = 0;


//	validKmerNum = 0;
	for(i=0; i<4; i++)
	{
		//occsNumSE[i] = 0;

		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				tmp_kmerseq[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
			}
			tmp_kmerseq[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*(lastEntryBaseNum-1)));
		}
		tmp_kmerseq[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | i) & lastEntryMask;

		tmp_kmers[i][0] = getKmer(tmp_kmerseq, deBruijnGraph);
		tmp_kmers[i][1] = getReverseKmer(tmp_kmerseqRev, tmp_kmerseq, deBruijnGraph);

//		if(tmp_kmers[i][0] || tmp_kmers[i][1])
//		{
//			validKmerNum ++;
//			base_index = i;
//		}
	}
/*
	if(validKmerNum==1)
	{
		//score[base_index] = computeKmerScore(tmp_kmers[base_index], occsNum+base_index, assemblingreads, numassemblingreads);
		//score[base_index] = computeKmerScoreUnlocked(tmp_kmers[base_index], occsNum+base_index, assemblingreads, numassemblingreads);
		if(computeKmerOccNumUnlocked(occsNumSE+base_index, tmp_kmers[base_index], base_index, minSuccessiveAppearedBaseNum, NO)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute kmer occNum by all reads, error!\n", __LINE__, __func__);
			return FAILED;
		}
		//========================= condition 1 ==============================
		if(occsNumSE[base_index]>0)
		//if(occsNumSE[base_index]>=minKmerOccSE)
		{
			if(entriesPerKmer>=2)
			{
				for(j=0; j<entriesPerKmer-2; j++)
				{
					kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
				}
				kmerSeqIntAssembly[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
			}
			kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | base_index) & lastEntryMask;

			kmers[0] = tmp_kmers[base_index][0];
			kmers[1] = tmp_kmers[base_index][1];

			maxOccIndexSE = base_index;
			maxOccSE = occsNumSE[maxOccIndexSE];
			secondOccIndexSE = -1;
			secondOccSE = 0;
			naviSuccessFlag = NAVI_SUCCESS;
		}else
		{
			naviSuccessFlag = NAVI_FAILED;
			kmers[0] = kmers[1] = NULL;
		}

		return SUCCESSFUL;

	}else if(validKmerNum==0)
	{
		naviSuccessFlag = NAVI_FAILED;
		kmers[0] = kmers[1] = NULL;
		return SUCCESSFUL;
	}
*/

	//****************************************************************************************************

	//kmer_len = longKmerSize;

	for(iterateID=0; iterateID<2; iterateID++)
	{
		if(iterateID==0)
		{
			successiveAppearBaseNum = maxSuccessiveAppearedBaseNum;
			ignoreEndBaseFlag = YES;
			//ignoreEndBaseFlag = NO;
		}else
		{
			successiveAppearBaseNum = minSuccessiveAppearedBaseNum;
			ignoreEndBaseFlag = NO;
		}

		if(contigNodesNum>=longKmerSize)
			kmer_len = longKmerSize;
		else
			kmer_len = contigNodesNum;

		//while(kmer_len>kmerSize)
		while(kmer_len>=kmerSize)
		//while(kmer_len>=MIN_KMER_SIZE)
		{
			validKmerNum = 0;
			maxOccSE = 0, secondOccSE = 0;
			maxOccIndexSE = -1, secondOccIndexSE = -1;
			for(i=0; i<4; i++)
			{
				occsNumSE[i] = 0;
	/*
				//######################### begin #############################//
				if(tmp_kmers[i][0]==NULL && tmp_kmers[i][1]==NULL)
				{
					continue;
				}

				//========================= condition 2 ==============================
				else
				{
					if(tmp_kmers[i][0]!=NULL && tmp_kmers[i][1]!=NULL)
					{
						//if(tmp_kmers[i][0]->multiplicity+tmp_kmers[i][1]->multiplicity<minKmerOccSE)
						if(tmp_kmers[i][0]->multiplicity+tmp_kmers[i][1]->multiplicity<=0)
						{
							continue;
						}
					}else if(tmp_kmers[i][0]!=NULL)
					{
						//if(tmp_kmers[i][0]->multiplicity<minKmerOccSE)
						if(tmp_kmers[i][0]->multiplicity<=0)
						{
							continue;
						}
					}else if(tmp_kmers[i][1]!=NULL)
					{
						//if(tmp_kmers[i][1]->multiplicity<minKmerOccSE)
						if(tmp_kmers[i][1]->multiplicity<=0)
						{
							continue;
						}
					}
				}
	*/
				//########################## end ########################//

				//if(contigNodesNum>READ_LEN)
				if(contigNodesNum>kmer_len)
				//if(contigNodesNum>=kmer_len)  //--bad result
				{
					if(computeLongKmerOccNum(occsNumSE+i, tmp_kmers[i], i, kmer_len, successiveAppearBaseNum, ignoreEndBaseFlag)==FAILED)
					{
						printf("line=%d, In %s(), cannot compute kmer occNum by long kmer, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}else
				{
					if(computeKmerOccNum(occsNumSE+i, tmp_kmers[i], i, successiveAppearBaseNum, ignoreEndBaseFlag)==FAILED)
					{
						printf("line=%d, In %s(), cannot compute kmer occNum, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}

				//score[i] = computeKmerScore(tmp_kmers[i], occsNum+i, assemblingreads, numassemblingreads);

				//========================= condition 3 ==============================
				//if(score[i]>=MIN_SCORE_THRESHOLD/**MIN_SCORE_FACTOR*/)
				//if(score[i]>=MIN_SCORE_THRESHOLD && occsNum[i]>=MIN_CONNECT_KMER_NUM)
				//if(occsNumSE[i]>=MIN_CONNECT_KMER_NUM)
				if(occsNumSE[i]>0)
				{
					validKmerNum ++;

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
				maxOccSE = 0, maxOccIndexSE = -1, secondOccSE = 0, secondOccIndexSE = -1;
				for(j=0; j<4; j++)
				{
					if(maxOccSE<occsNumSE[j])
					{
						secondOccSE = maxOccSE;
						secondOccIndexSE = maxOccIndexSE;
						maxOccSE = occsNumSE[j];
						maxOccIndexSE = j;
					}else if(secondOccSE<occsNumSE[j])
					{
						secondOccSE = occsNumSE[j];
						secondOccIndexSE = j;
					}
				}

				occsNumIndexSE[0] = maxOccIndexSE;
				occsNumIndexSE[1] = secondOccIndexSE;
			}

	/*
			//========================= condition 4 ==============================
			//if(validKmerNum>1 && maxIndex1!=maxOccIndex)
			//if(validKmerNum>1 && occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM)
			//if(validKmerNum>0 && occsNumSE[maxIndex1]<minKmerOccSE) //--best result
			//if(occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM || (validKmerNum>1 && maxIndex1!=maxOccIndex))
			//if(occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM || occsNum[secondOccIndex]>=MIN_CONNECT_KMER_NUM)
			if(validKmerNum>0 && maxOcc<minKmerOccSE)  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
			if(validKmerNum>1 && (occsNumSE[secondIndex1] > maxSecondOcc || (occsNumSE[maxIndex1] > maxFirstOcc && occsNumSE[secondIndex1]>=maxSecondOcc)
					|| ((float)occsNumSE[secondIndex1]/occsNumSE[maxIndex1]>SECOND_FIRST_OCC_RATIO && occsNumSE[secondIndex1]>=maxSecondOcc))) //-- best result
			{
				validKmerNum = 0;
			}
	*/

	/*
			//========================= condition 6 ==============================
			//if(validKmerNum>1 && ((float)occsNumSE[maxIndex1]/itemNumDecisionTable < VALID_OCC_RATIO && occsNumSE[secondIndex1] >= maxSecondOcc))
			if(validKmerNum>1 && ((float)maxOcc/itemNumDecisionTable < VALID_OCC_RATIO && secondOcc > maxSecondOcc))
			{
				validKmerNum = 0;
			}


			//========================= condition 7 ==============================
			if(validKmerNum>1 && second/max > SECOND_FIRST_SECORE_RATIO)
			{
				validKmerNum = 0;
			}
	*/

#if(USING_RULES==YES)
			//========================= condition 8 ==============================
			//=====these several lines have bad result, thus they are omitted. =======//
			//if(validKmerNum>1 && (secondOcc/maxOcc>SECOND_FIRST_OCC_RATIO || (kmer_len==25 && secondOcc>SECOND_OCC_THRESHOLD)))
			//if(validKmerNum>1 && kmer_len==25 && ((secondOcc/maxOcc>SECOND_FIRST_OCC_RATIO && secondOcc>MIN_CONNECT_KMER_NUM) || (secondOcc>SECOND_OCC_THRESHOLD)))
			//if(validKmerNum>1 && ((secondOccSE/maxOccSE>SECOND_FIRST_OCC_RATIO && secondOccSE>minKmerOccSE) || (secondOccSE>maxSecondOcc)))
			if(maxOccSE<minKmerOccSE || (validKmerNum>1 && ((secondOccSE/maxOccSE>SECOND_FIRST_OCC_RATIO && secondOccSE>minKmerOccSE) || (secondOccSE>maxSecondOcc))))
			//if(maxOccSE<minKmerOccSE || (validKmerNum>1 && (secondOccSE/maxOccSE>SECOND_FIRST_OCC_RATIO || (secondOccSE>maxSecondOcc))))
			{
				if(secondOccSE/maxOccSE>SECOND_FIRST_OCC_FAILED_RATIO)
					validKmerNum = 0;
				else if(secondOccSE>0 && maxOccSE-secondOccSE<2)
					validKmerNum = 0;
			}


			//========================= condition 9 ==============================
			//if(kmer_len > longKmerSize - longKmerStepSize && occsNumSE[maxIndex1] <= minLongKmerOcc)
			//if(kmer_len > longKmerSize - longKmerStepSize && maxOccSE < minLongKmerOcc)
			if(kmer_len > kmerSize && maxOccSE < minLongKmerOcc)
			{
				validKmerNum = 0;
			}
#endif


			if(validKmerNum==1)
			{
				//if(occsNum[tmp_maxIndex]<MIN_CONNECT_KMER_NUM)
	//			if(occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM)
	//			{
	//				kmers[0] = kmers[1] = NULL;
	//			}else
	//			{
					if(entriesPerKmer>=2)
					{
						for(j=0; j<entriesPerKmer-2; j++)
						{
							kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
						}
						kmerSeqIntAssembly[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
					}
					kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | maxOccIndexSE) & lastEntryMask;

					kmers[0] = tmp_kmers[maxOccIndexSE][0];
					kmers[1] = tmp_kmers[maxOccIndexSE][1];

					naviSuccessFlag = NAVI_SUCCESS;
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

#if(USING_RULES==YES)
				if(secondOccSE>4*minKmerOccSE && secondOccSE/maxOccSE>SECOND_FIRST_OCC_FAILED_RATIO)
				{
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;

					return SUCCESSFUL;
				}else
				{
					kmer_len -= longKmerStepSize;
					if(kmer_len<kmerSize && kmer_len>kmerSize-longKmerStepSize)
						kmer_len = kmerSize;

					continue;
				}
#endif

				kmer_len -= longKmerStepSize;
				if(kmer_len<kmerSize && kmer_len>kmerSize-longKmerStepSize)
					kmer_len = kmerSize;

				continue;
			}
			else
			{
				if(entriesPerKmer>=2)
				{
					for(j=0; j<entriesPerKmer-2; j++)
					{
						kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
					}
					kmerSeqIntAssembly[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
				}
				kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | maxOccIndexSE) & lastEntryMask;

				kmers[0] = tmp_kmers[maxOccIndexSE][0];
				kmers[1] = tmp_kmers[maxOccIndexSE][1];

				naviSuccessFlag = NAVI_SUCCESS;

				return SUCCESSFUL;
			}
			//***********************************************************************************
		}
	}

	naviSuccessFlag = NAVI_FAILED;
	kmers[0] = kmers[1] = NULL;

	return SUCCESSFUL;
}


/**
 * Compute the occNum according to long K-mers.
 * 	@return:
 * 		If succeeds, return SUCCSEEFUL; otherwise, return FAILED.
 */
short computeLongKmerOccNum(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t length_k, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag)
{
	assemblingreadtype *this_assemblingRead;
	ridpostype *rid_pos_table;
	readBlock_t *readBlockArr;
	read_t *pRead;

	int32_t i, properRow, posNum, basePos, baseInt_read;
	int32_t entryIDReadseq, entryPosReadseq, baseNumEntryReadseq;
	int32_t rid, maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3;


	*occNum = 0;
	for(i=0; i<itemNumDecisionTable; i++)
	{
		if((decisionTable[i].reserved==0)
				//&& (decisionTable[i].unmatchBaseNum==0)
				&& (decisionTable[i].unmatchBaseNum<=maxUnmatchBaseNumPerRead)
				&& (decisionTable[i].successiveAppearBases>=successiveAppearBaseNum)	// added 2013-01-12
				&& decisionTable[i].basePos==decisionTable[i].lastMatchedBasePos
				&& ((decisionTable[i].orientation==ORIENTATION_PLUS && decisionTable[i].basePos>=length_k-2)
					||(decisionTable[i].orientation==ORIENTATION_MINUS && decisionTable[i].basePos<=decisionTable[i].seqlen-length_k+1)))
		{
			properRow = getProperDtRow(decisionTable+i, decisionTable, dtRowHashtable);

			// ############################ Debug information ##############################
			if(properRow<0 || properRow>=itemNumDecisionTable || decisionTable[properRow].rid!=decisionTable[i].rid)
			{
				printf("line=%d, In %s(), properRow=%d, i=%d, Error!\n", __LINE__, __func__, properRow, i);
				continue;
			}
			// ############################ Debug information ##############################

			this_assemblingRead = decisionTable + properRow;

			// ignore the end bases to avoid erroneous bases
			if(ignoreEndBaseFlag==YES)
			{
				if((this_assemblingRead->orientation==ORIENTATION_PLUS && this_assemblingRead->basePos>this_assemblingRead->seqlen-this_assemblingRead->kmerNumEnd3)
					|| (this_assemblingRead->orientation==ORIENTATION_MINUS && this_assemblingRead->basePos<this_assemblingRead->kmerNumEnd5))
					continue;
			}

			// get the base
			if(this_assemblingRead->orientation==ORIENTATION_PLUS)
				basePos = this_assemblingRead->basePos + 1;
			else
				basePos = this_assemblingRead->basePos - 1;

			entryIDReadseq = basePos >> 5;
			entryPosReadseq = basePos % 32;

			if(entryIDReadseq==this_assemblingRead->entriesNumReadseq-1)
				baseNumEntryReadseq = this_assemblingRead->baseNumLastEentryReadseq;
			else
				baseNumEntryReadseq = 32;

			baseInt_read = (this_assemblingRead->readseq[entryIDReadseq] >> (2*(baseNumEntryReadseq-entryPosReadseq-1))) & 3;

			if(this_assemblingRead->orientation==ORIENTATION_PLUS)
			{ // plus orientation
				if(baseInt_read==baseInt_kmer)
				{ // base match
					(*occNum) ++;
				}
			}else
			{ // minus orientation
				if(((~baseInt_read)&3)==baseInt_kmer)
				{ // base match
					(*occNum) ++;
				}
			}
		} //end if(reserved)
	}// end for(i)

#if 1
	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;

	if(length_k<kmerSize+longKmerStepSize)
	//if(length_k==kmerSize)
	{
		if(tmp_kmers[0])
		{
			rid_pos_table = tmp_kmers[0]->ppos;
			posNum = tmp_kmers[0]->arraysize;
			for(i=0; i<posNum; i++)
			{
				if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
				{
					(*occNum) ++;
				}
				//rid_pos_table[i].reserved = 0;
			}
		}

		if(tmp_kmers[1])
		{
			rid_pos_table = tmp_kmers[1]->ppos;
			posNum = tmp_kmers[1]->arraysize;
			for(i=0; i<posNum; i++)
			{
				// get seqLen and errorRegLenEnd3
				rid = rid_pos_table[i].rid;
				blockRowRead = (rid - 1) / maxItemNumPerReadBlock;
				itemRowBlockRead = (rid - 1) % maxItemNumPerReadBlock;
				pRead = readBlockArr[blockRowRead].readArr + itemRowBlockRead;
				seqLen = pRead->seqlen;
				errorRegLenEnd3 = ceil(seqLen * ERROR_REGION_LEN_3End_FACTOR);

				if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos>seqLen-kmerSize+1-errorRegLenEnd3)
				{
					if(existReadInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
						(*occNum) ++;
				}
				//rid_pos_table[i].reserved = 0;
			}
		}
	}
#endif

	for(i=0; i<itemNumDecisionTable; i++) decisionTable[i].reserved = 0;

	return SUCCESSFUL;
}


/**
 * Compute kmer occurrence number.
 *  If there are some locked reads, the occNum will be computed by the locked reads;
 *  otherwise, consider all reads.
 *
 *  @return:
 *   	If succeeds, return SUCCSEEFUL; otherwise, return FAILED.
*/
short computeKmerOccNum(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag)
{
	//if(lockedReadsNum>=KOCKED_READS_NUM_THRESHOLD)
	if(lockedReadsNum>=lockedReadsNumThres)
	{
		if(computeKmerOccNumLocked(occNum, tmp_kmers, baseInt_kmer, successiveAppearBaseNum, ignoreEndBaseFlag)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute kmer occNum by locked reads, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		if(computeKmerOccNumUnlocked(occNum, tmp_kmers, baseInt_kmer, successiveAppearBaseNum, ignoreEndBaseFlag)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute kmer occNum by all reads, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Compute kmer occNum by all reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeKmerOccNumUnlocked(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag)
{
	assemblingreadtype *this_assemblingRead;
	ridpostype *rid_pos_table;
	readBlock_t *readBlockArr;
	read_t *pRead;

	int32_t i, properRow, posNum, basePos, baseInt_read;
	int32_t entryIDReadseq, entryPosReadseq, baseNumEntryReadseq;
	int32_t rid, maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3;


	*occNum = 0;
	for(i=0; i<itemNumDecisionTable; i++)
	{
		//if(decisionTable[i].reserved==0)
		//if(decisionTable[i].reserved==0 && decisionTable[i].unmatchBaseNum<=3)

		if(decisionTable[i].reserved==0
			&& decisionTable[i].unmatchBaseNum<=maxUnmatchBaseNumPerRead
			&& (decisionTable[i].successiveAppearBases>=successiveAppearBaseNum)	// added 2013-01-12
			&& decisionTable[i].basePos==decisionTable[i].lastMatchedBasePos)
		{
			properRow = getProperDtRow(decisionTable+i, decisionTable, dtRowHashtable);

			// ############################ Debug information ##############################
			if(properRow<0 || properRow>=itemNumDecisionTable || decisionTable[properRow].rid!=decisionTable[i].rid)
			{
				printf("line=%d, In %s(), properRow=%d, i=%d, Error!\n", __LINE__, __func__, properRow, i);
				//continue;
				return FAILED;
			}
			// ############################ Debug information ##############################

			this_assemblingRead = decisionTable + properRow;

			// ignore the end bases to avoid erroneous bases
			if(ignoreEndBaseFlag==YES)
			{
				if((this_assemblingRead->orientation==ORIENTATION_PLUS && this_assemblingRead->basePos>this_assemblingRead->seqlen-this_assemblingRead->kmerNumEnd3)
					|| (this_assemblingRead->orientation==ORIENTATION_MINUS && this_assemblingRead->basePos<this_assemblingRead->kmerNumEnd5))
					continue;
			}

			// get the base
			if(this_assemblingRead->orientation==ORIENTATION_PLUS)
				basePos = this_assemblingRead->basePos + 1;
			else
				basePos = this_assemblingRead->basePos - 1;

			entryIDReadseq = basePos >> 5;
			entryPosReadseq = basePos % 32;

			if(entryIDReadseq==this_assemblingRead->entriesNumReadseq-1)
				baseNumEntryReadseq = this_assemblingRead->baseNumLastEentryReadseq;
			else
				baseNumEntryReadseq = 32;

			baseInt_read = (this_assemblingRead->readseq[entryIDReadseq] >> (2*(baseNumEntryReadseq-entryPosReadseq-1))) & 3;

			if(this_assemblingRead->orientation==ORIENTATION_PLUS)
			{ // plus orientation
				if(baseInt_read==baseInt_kmer)
				{ // base match
					(*occNum) ++;
				}
			}else
			{ // minus orientation
				if(((~baseInt_read)&3)==baseInt_kmer)
				{ // base match
					(*occNum) ++;
				}
			}
		} //end if(reserved)
	}// end for(i)


#if 1
	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;

	if(tmp_kmers[0])
	{
		rid_pos_table = tmp_kmers[0]->ppos;
		posNum = tmp_kmers[0]->arraysize;
		for(i=0; i<posNum; i++)
		{
			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
			{
				(*occNum) ++;
			}
			//rid_pos_table[i].reserved = 0;
		}
	}

	if(tmp_kmers[1])
	{
		rid_pos_table = tmp_kmers[1]->ppos;
		posNum = tmp_kmers[1]->arraysize;
		for(i=0; i<posNum; i++)
		{
			// get seqLen and errorRegLenEnd3
			rid = rid_pos_table[i].rid;
			blockRowRead = (rid - 1) / maxItemNumPerReadBlock;
			itemRowBlockRead = (rid - 1) % maxItemNumPerReadBlock;
			pRead = readBlockArr[blockRowRead].readArr + itemRowBlockRead;
			seqLen = pRead->seqlen;
			errorRegLenEnd3 = ceil(seqLen * ERROR_REGION_LEN_3End_FACTOR);

			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos>seqLen-kmerSize+1-errorRegLenEnd3)
			{
				if(existReadInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
					(*occNum) ++;
			}
			//rid_pos_table[i].reserved = 0;
		}
	}
#endif

	for(i=0; i<itemNumDecisionTable; i++) decisionTable[i].reserved = 0;

	return SUCCESSFUL;
}


/**
 * Compute the kmer occurrence number by locked reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeKmerOccNumLocked(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag)
{
	assemblingreadtype *this_assemblingRead;
	ridpostype *rid_pos_table;
	readBlock_t *readBlockArr;
	read_t *pRead;

	int32_t i, properRow, posNum, basePos, baseInt_read;
	int32_t entryIDReadseq, entryPosReadseq, baseNumEntryReadseq;
	int32_t rid, maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3;


	*occNum = 0;

	for(i=0; i<itemNumDecisionTable; i++)
	{
		if(decisionTable[i].locked)
		{
			//if(decisionTable[i].reserved==0)
			//if(decisionTable[i].reserved==0 && decisionTable[i].unmatchBaseNum<=3)
			if(decisionTable[i].reserved==0
				&& decisionTable[i].unmatchBaseNum<=maxUnmatchBaseNumPerRead
				&& (decisionTable[i].successiveAppearBases>=successiveAppearBaseNum)	// added 2013-01-12
				&& decisionTable[i].basePos==decisionTable[i].lastMatchedBasePos)
			{
				properRow = getProperDtRow(decisionTable+i, decisionTable, dtRowHashtable);

				// ############################ Debug information ##############################
				if(properRow<0 || properRow>=itemNumDecisionTable || decisionTable[properRow].rid!=decisionTable[i].rid)
				{
					printf("line=%d, In %s(), properRow=%d, i=%d, Error!\n", __LINE__, __func__, properRow, i);
					//continue;
					return FAILED;
				}
				// ############################ Debug information ##############################

				this_assemblingRead = decisionTable + properRow;

				// ignore the end bases to avoid erroneous bases
				if(ignoreEndBaseFlag==YES)
				{
					if((this_assemblingRead->orientation==ORIENTATION_PLUS && this_assemblingRead->basePos>this_assemblingRead->seqlen-this_assemblingRead->kmerNumEnd3)
						|| (this_assemblingRead->orientation==ORIENTATION_MINUS && this_assemblingRead->basePos<this_assemblingRead->kmerNumEnd5))
						continue;
				}

				// get the base
				if(this_assemblingRead->orientation==ORIENTATION_PLUS)
					basePos = this_assemblingRead->basePos + 1;
				else
					basePos = this_assemblingRead->basePos - 1;

				entryIDReadseq = basePos >> 5;
				entryPosReadseq = basePos % 32;

				if(entryIDReadseq==this_assemblingRead->entriesNumReadseq-1)
					baseNumEntryReadseq = this_assemblingRead->baseNumLastEentryReadseq;
				else
					baseNumEntryReadseq = 32;

				baseInt_read = (this_assemblingRead->readseq[entryIDReadseq] >> (2*(baseNumEntryReadseq-entryPosReadseq-1))) & 3;

				if(this_assemblingRead->orientation==ORIENTATION_PLUS)
				{ // plus orientation
					if(baseInt_read==baseInt_kmer)
					{ // base match
						(*occNum) ++;
					}
				}else
				{ // minus orientation
					if(((~baseInt_read)&3)==baseInt_kmer)
					{ // base match
						(*occNum) ++;
					}
				}
			} //end if(reserved)
		} //end if (locked)
	}// end for(i)

#if 1
	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;

	if(tmp_kmers[0])
	{
		rid_pos_table = tmp_kmers[0]->ppos;
		posNum = tmp_kmers[0]->arraysize;
		for(i=0; i<posNum; i++)
		{
			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
			{
				(*occNum) ++;
			}
			//rid_pos_table[i].reserved = 0;
		}
	}

	if(tmp_kmers[1])
	{
		rid_pos_table = tmp_kmers[1]->ppos;
		posNum = tmp_kmers[1]->arraysize;
		for(i=0; i<posNum; i++)
		{
			// get seqLen and errorRegLenEnd3
			rid = rid_pos_table[i].rid;
			blockRowRead = (rid - 1) / maxItemNumPerReadBlock;
			itemRowBlockRead = (rid - 1) % maxItemNumPerReadBlock;
			pRead = readBlockArr[blockRowRead].readArr + itemRowBlockRead;
			seqLen = pRead->seqlen;
			errorRegLenEnd3 = ceil(seqLen * ERROR_REGION_LEN_3End_FACTOR);

			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos>seqLen-kmerSize+1-errorRegLenEnd3)
			{
				if(existReadInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
					(*occNum) ++;
			}
			//rid_pos_table[i].reserved = 0;
		}
	}
#endif

	for(i=0; i<itemNumDecisionTable; i++) decisionTable[i].reserved = 0;

	return SUCCESSFUL;
}

/**
 * Append a base to contig tail.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
*/
short addContigBase(uint32_t baseInt)
{
	if(itemNumContigArr>maxItemNumContigArr)
	{
		contigArr = (contigtype *) realloc(contigArr, 2*maxItemNumContigArr*sizeof(contigtype));
		if(contigArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		maxItemNumContigArr = 2 * maxItemNumContigArr;
	}

	contigArr[itemNumContigArr-1].index = itemNumContigArr;
	contigArr[itemNumContigArr-1].base = baseInt;
	contigArr[itemNumContigArr-1].ridposnum = 0;
	contigArr[itemNumContigArr-1].pridposorientation = NULL;
	contigArr[itemNumContigArr-1].reserved = 0;

	return SUCCESSFUL;
}

/**
 * Add ridpos to contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
*/
short addRidposToContig(successRead_t *successReadArray, int successReadNum, int contigNodesNum)
{
	successRead_t *ridposorient;
	int32_t i, j, num, k;
	int64_t hangingIndex, minHangingIndex, maxHangingIndex;

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{ // first round assembly
		// allocate memory
		ridposorient = (successRead_t*) malloc(successReadNum*sizeof(successRead_t));
		if(ridposorient==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		contigArr[contigNodesNum-1].pridposorientation = ridposorient;
		contigArr[contigNodesNum-1].ridposnum = successReadNum;

		// fill the reads
		for(i=0; i<successReadNum; i++)
		{
			if(successReadArray[i].orientation==ORIENTATION_PLUS)
				ridposorient[i].orientation = ORIENTATION_MINUS;
			else
				ridposorient[i].orientation = ORIENTATION_PLUS;
			ridposorient[i].rid = successReadArray[i].rid;
			ridposorient[i].pos = successReadArray[i].seqlen - successReadArray[i].pos + 1;
			ridposorient[i].matchnum = successReadArray[i].matchnum;
			ridposorient[i].matchlen = successReadArray[i].matchlen;
			ridposorient[i].seqlen = successReadArray[i].seqlen;
			ridposorient[i].hangingIndex = contigNodesNum;
			ridposorient[i].readseq = successReadArray[i].readseq;
			ridposorient[i].entriesNumReadseq = successReadArray[i].entriesNumReadseq;
			ridposorient[i].baseNumLastEentryReadseq = successReadArray[i].baseNumLastEentryReadseq;
			ridposorient[i].kmerNumEnd5 = successReadArray[i].kmerNumEnd5;
			ridposorient[i].kmerNumEnd3 = successReadArray[i].kmerNumEnd3;
		}
	}
	else
	{ // second round assembly
		if(successReadNum==1)
		{ // only one read
			hangingIndex = contigNodesNum - successReadArray[0].matchlen + 1;
			if(hangingIndex<=0)
				hangingIndex = 1;
			successReadArray[0].hangingIndex = hangingIndex;

			ridposorient = (successRead_t*) malloc((contigArr[hangingIndex-1].ridposnum + 1)*sizeof(successRead_t));
			if(ridposorient==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(contigArr[hangingIndex-1].ridposnum>0)
			{
				if(memcpy(ridposorient, contigArr[hangingIndex-1].pridposorientation, contigArr[hangingIndex-1].ridposnum*sizeof(successRead_t))==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}

				free(contigArr[hangingIndex-1].pridposorientation);
				contigArr[hangingIndex-1].pridposorientation = NULL;
			}

			contigArr[hangingIndex-1].pridposorientation = ridposorient;
			if(memcpy(ridposorient+contigArr[hangingIndex-1].ridposnum, successReadArray, sizeof(successRead_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			contigArr[hangingIndex-1].ridposnum ++;

		}else
		{ // more than one reads
			maxHangingIndex = 0;
			minHangingIndex = contigNodesNum;
			for(i=0; i<successReadNum; i++)
			{
				// ########################## Debug information ########################
				//if(successReadArray[i].rid==23356022)
				//{
				//	printf("line=%d, In %s(), rid=%ld\n", __LINE__, __func__, (int64_t)successReadArray[i].rid);
				//}
				// ########################## Debug information ########################

				hangingIndex = contigNodesNum - successReadArray[i].matchlen + 1;
				if(hangingIndex<=0)
					hangingIndex = 1;
				successReadArray[i].hangingIndex = hangingIndex;

				if(maxHangingIndex<hangingIndex)
					maxHangingIndex = hangingIndex;
				if(minHangingIndex>hangingIndex)
					minHangingIndex = hangingIndex;
			}

			for(i=0; i<successReadNum; i++)
				countingHangingBucketArr[successReadArray[i].hangingIndex-minHangingIndex] ++;

			for(i=minHangingIndex; i<=maxHangingIndex; i++)
			{
				if(countingHangingBucketArr[i-minHangingIndex]>0)
				{
					num = countingHangingBucketArr[i-minHangingIndex];

					ridposorient = (successRead_t*) malloc((contigArr[i-1].ridposnum + num)*sizeof(successRead_t));
					if(ridposorient==NULL)
					{
						printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(contigArr[i-1].ridposnum>0)
					{
						if(memcpy(ridposorient, contigArr[i-1].pridposorientation, contigArr[i-1].ridposnum*sizeof(successRead_t))==NULL)
						{
							printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
							return FAILED;
						}

						free(contigArr[i-1].pridposorientation);
						contigArr[i-1].pridposorientation = NULL;
					}

					contigArr[i-1].pridposorientation = ridposorient;

					k = contigArr[i-1].ridposnum;
					for(j=0; j<successReadNum; j++)
					{
						if(successReadArray[j].hangingIndex==i)
						{
							if(memcpy(ridposorient+k, successReadArray+j, sizeof(successRead_t))==NULL)
							{
								printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
								return FAILED;
							}
							k ++;
						}
					}

					contigArr[i-1].ridposnum = k;
					countingHangingBucketArr[i-minHangingIndex] = 0;
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Clean contig nodes.
 */
void cleanContigArray(contigtype *contigArr, int64_t *contigNodesNum)
{
	int64_t i;

	for(i=0; i<*contigNodesNum; i++)
	{
		if(contigArr[i].ridposnum>0)
			free(contigArr[i].pridposorientation);
	}

	*contigNodesNum = 0;
}

/**
 * Output the contig nodes to file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputContigToFile(FILE *fpContig, int outFileType, int contigID, contigtype *contigArr, int64_t contigNodeNum)
{
	successRead_t *ridposorientation;
	int i, j, num, indexTmp;
	char tmp_base;

	if(fpContig)
	{
		if(outFileType==BASE_TYPE_FASTA_CONTIG_FILE)
		{
			fprintf(fpContig, ">%d length: %ld\n", contigID, contigNodeNum);
			for(i=0; i<contigNodeNum; i++)
			{
				switch(contigArr[i].base)
				{
					case 0: tmp_base = 'A'; break;
					case 1: tmp_base = 'C'; break;
					case 2: tmp_base = 'G'; break;
					case 3: tmp_base = 'T'; break;
					default: printf("line=%d, In %s(), error baseInt: %d\n", __LINE__, __func__, contigArr[i].base); return FAILED;
				}

				fprintf(fpContig, "%c", tmp_base);
			}
			fprintf(fpContig, "\n");
			//fflush(pFileResult);
		}else if(outFileType==HANGING_READ_TYPE_CONTIG_FILE)
		{ // output the contigs in hanging reads format
			fprintf(fpContig, ">%d length: %ld\n", contigID, contigNodeNum);

			indexTmp = 1;
			for(i=0; i<contigNodeNum; i++)
			{
				switch(contigArr[i].base)
				{
					case 0: tmp_base = 'A'; break;
					case 1: tmp_base = 'C'; break;
					case 2: tmp_base = 'G'; break;
					case 3: tmp_base = 'T'; break;
					default: printf("line=%d, In %s(), error base int: %d\n", __LINE__, __func__, contigArr[i].base); return FAILED;
				}

				fprintf(fpContig, "%d\t%c\t%d", indexTmp++, tmp_base, contigArr[i].ridposnum);
				ridposorientation = contigArr[i].pridposorientation;
				num = contigArr[i].ridposnum;
				for(j=0; j<num; j++)
					fprintf(fpContig, "\t(%lu,%d,%d,%c)", (uint64_t)ridposorientation[j].rid, ridposorientation[j].pos, ridposorientation[j].matchnum, ridposorientation[j].orientation);
				fprintf(fpContig, "\n");

			}
			//fflush(pFileResult);
		}else
		{
			printf("line=%d, In %s(), unknown output contig file type %d, error!\n", __LINE__, __func__, outFileType);
			return FAILED;
		}
	}else
	{
		printf("line=%d, In %s(), file pointer is NULL, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Delete the successful reads from de Bruijn graph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short delReadsFromGraph(successRead_t *successReadArray, int successReadNum)
{
	uint64_t tmp_kmerseq[deBruijnGraph->entriesPerKmer];
	int32_t i, j, k, rpos, basePos, startKmerPos, startBasePos, endBasePos, endKmerNum;
	int32_t baseInt, seqLen, entriesNum, baseNumLastEntry, entryRow, entryPos;
	uint64_t hashcode, rid, *readseqInt, maxItemNumPerReadBlock;
	readBlock_t *readBlockArr;
	read_t *pRead;

	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;

	for(i=0; i<successReadNum; i++)
	{
		rid = successReadArray[i].rid;
		readseqInt = successReadArray[i].readseq;
		seqLen = successReadArray[i].seqlen;
		entriesNum = successReadArray[i].entriesNumReadseq;
		baseNumLastEntry = successReadArray[i].baseNumLastEentryReadseq;

		//if(rid==239089)
		//{
		//	printf("line=%d, In %s(), rid=%lu\n", __LINE__, __func__, rid);
		//}

		pRead = readBlockArr[(rid - 1) / maxItemNumPerReadBlock].readArr + ((rid - 1) % maxItemNumPerReadBlock);
		pRead->successFlag = YES;

		if(ceil(seqLen * kmerRegLenRatioEnd5) + ceil(seqLen * kmerRegLenRatioEnd3) >= seqLen - kmerSize + 1)
		{ // the 5' + 3' >= kmerNum
			// generate the kmer integer sequence
			if(generateKmerSeqIntFromReadset(tmp_kmerseq, successReadArray[i].readseq, 0, entriesNum, baseNumLastEntry)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			rpos = 1;
			hashcode = kmerhashInt(tmp_kmerseq);
			if(delKmerByHash(hashcode, tmp_kmerseq, rid, rpos, deBruijnGraph)==FAILED)
			{
				printf("line=%d, In %s(), can not delete the read [%lu, %c]. Error!\n", __LINE__, __func__, rid, successReadArray[i].orientation);
				return FAILED;
			}

			rpos ++;
			entryRow = kmerSize >> 5;
			entryPos = kmerSize % 32;
			for(basePos=kmerSize; basePos<seqLen; basePos++, rpos++)
			{
				// get the baseInt
				if(entryRow<entriesNum-1)
					baseInt = (readseqInt[entryRow] >> (62-2*entryPos)) & 3;
				else
					baseInt = (readseqInt[entryRow] >> (2*(baseNumLastEntry-1)-2*entryPos)) & 3;

				entryPos ++;
				if(entryPos==32)
				{
					entryRow ++;
					entryPos = 0;
				}

				// generate the kmer integer sequence
				if(entriesPerKmer>=2)
				{
					for(j=0; j<entriesPerKmer-2; j++)
					{
						tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
					}

					tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
				}
				tmp_kmerseq[entriesPerKmer-1] = ((tmp_kmerseq[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

				hashcode = kmerhashInt(tmp_kmerseq);
				if(delKmerByHash(hashcode, tmp_kmerseq, rid, rpos, deBruijnGraph)==FAILED)
				{
					printf("line=%d, In %s(), can not delete the read [%lu, %c]. Error!\n", __LINE__, __func__, rid, successReadArray[i].orientation);
					return FAILED;
				}
			}

		}else
		{
			for(k=0; k<2; k++)
			{
				if(k==0)
				{
					endKmerNum = ceil(seqLen * kmerRegLenRatioEnd5);
					startKmerPos = 0;
					startBasePos = kmerSize;
					endBasePos = startBasePos + endKmerNum - 2;
				}
				else
				{
					endKmerNum = ceil(seqLen * kmerRegLenRatioEnd3);
					startKmerPos = seqLen - endKmerNum - kmerSize + 1;
					startBasePos = seqLen - endKmerNum + 1;
					endBasePos = seqLen - 1;
				}

				// generate the kmer integer sequence
				if(generateKmerSeqIntFromReadset(tmp_kmerseq, readseqInt, startKmerPos, entriesNum, baseNumLastEntry)==FAILED)
				{
					printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
					return FAILED;
				}

				rpos = startKmerPos + 1;
				hashcode = kmerhashInt(tmp_kmerseq);
				if(delKmerByHash(hashcode, tmp_kmerseq, rid, rpos, deBruijnGraph)==FAILED)
				{
					printf("line=%d, In %s(), can not delete the read [%lu, %c]. Error!\n", __LINE__, __func__, rid, successReadArray[i].orientation);
					return FAILED;
				}

				rpos ++;
				entryRow = startBasePos >> 5;
				entryPos = startBasePos % 32;
				for(basePos=startBasePos; basePos<=endBasePos; basePos++, rpos++)
				{
					// get the baseInt
					if(entryRow<entriesNum-1)
						baseInt = (readseqInt[entryRow] >> (62-2*entryPos)) & 3;
					else
						baseInt = (readseqInt[entryRow] >> (2*(baseNumLastEntry-1)-2*entryPos)) & 3;

					entryPos ++;
					if(entryPos==32)
					{
						entryRow ++;
						entryPos = 0;
					}

					// generate the kmer integer sequence
					if(entriesPerKmer>=2)
					{
						for(j=0; j<entriesPerKmer-2; j++)
						{
							tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
						}

						tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
					}
					tmp_kmerseq[entriesPerKmer-1] = ((tmp_kmerseq[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

					hashcode = kmerhashInt(tmp_kmerseq);
					if(delKmerByHash(hashcode, tmp_kmerseq, rid, rpos, deBruijnGraph)==FAILED)
					{
						printf("line=%d, In %s(), can not delete the read [%lu, %c]. Error!\n", __LINE__, __func__, rid, successReadArray[i].orientation);
						return FAILED;
					}
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Reverse the read sequence.
 *  @return:
 *  	If succeeds, return SUCCESFUL; otherwise, return FAILED.
 */
short reverseReadseq(char *str)
{
	if(!str) { printf("line=%d, In %s(): the tuple is NULL! Error!\n", __LINE__, __func__); return FAILED; }
	int i = 0, len = strlen(str);
	char tmp[len+1];
	strcpy(tmp, str);

	while(i<len)
	{
		switch(tmp[i])
		{
			case 'A': str[len-i-1] = 'T'; break;
			case 'C': str[len-i-1] = 'G'; break;
			case 'G': str[len-i-1] = 'C'; break;
			case 'T': str[len-i-1] = 'A'; break;
			case 'a': str[len-i-1] = 't'; break;
			case 'c': str[len-i-1] = 'g'; break;
			case 'g': str[len-i-1] = 'c'; break;
			case 't': str[len-i-1] = 'a'; break;
			default:  str[len-i-1] = tmp[i];  return FAILED;
		}
		i++;
	}
	return SUCCESSFUL;
}

/**
 * Update the locked reads and their total number.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateLockedReads()
{
	assemblingreadtype *this_assemblingRead;

	int i = 0, appearNum = 0;
	if(lockedReadsNum>0)
	{
		this_assemblingRead = decisionTable;
		for(i=0; i<itemNumDecisionTable; i++)
		{
			if(this_assemblingRead->locked==1)
			{
				if(this_assemblingRead->status!=ASSEMBLING_STATUS)
				{
					lockedReadsNum --;
				}else if(this_assemblingRead->basePos==this_assemblingRead->lastMatchedBasePos)
				{
					appearNum++;
				}
			}
			this_assemblingRead ++;
		}
	}

	//if(appearNum<minKmerOccSE*2 || lockedReadsNum<KOCKED_READS_NUM_THRESHOLD)
	if(appearNum<minKmerOccSE*2 || lockedReadsNum<lockedReadsNumThres)
	{

//		if(appearNum>=MIN_CONNECT_KMER_NUM*2)
//		{
//			printf("appearNum=%d, lockedReadsNum=%d\n", appearNum, lockedReadsNum);
//		}else
//		{
//			printf("++++++++++appearNum=%d, lockedReadsNum=%d\n", appearNum, lockedReadsNum);
//		}

		lockedReadsNum = 0;

		this_assemblingRead = decisionTable;
		for(i=0; i<itemNumDecisionTable; i++)
		{

			// ############################ Debug information ##############################
			//if(this_assemblingRead->rid==1655993)
			//{
			//	printf("line=%d, In %s(), rid=%d\n", __LINE__, __func__, this_assemblingRead->rid);
			//}
			// ############################ Debug information ##############################

			if(this_assemblingRead->basePos!=this_assemblingRead->lastMatchedBasePos && this_assemblingRead->locked==1)
			{
				this_assemblingRead->locked = 0;
			}else if(this_assemblingRead->basePos==this_assemblingRead->lastMatchedBasePos && this_assemblingRead->status==ASSEMBLING_STATUS)
			{
				this_assemblingRead->locked = 1;
				lockedReadsNum ++;
			}
			this_assemblingRead ++;
		}
	}

	return SUCCESSFUL;
}

/**
 * Trim contig tail nodes.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateContigtailnodes(contigtype *contigArr, int64_t successContigIndex, int64_t *contigNodesNum)
{
	int i;

	for(i=successContigIndex-1; i<(*contigNodesNum); i++)
	{
		if(contigArr[i].ridposnum)
		{
			free(contigArr[i].pridposorientation);
			contigArr[i].pridposorientation = NULL;

			printf("line=%d, In %s, ridpos>0, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	*contigNodesNum = successContigIndex;

	return SUCCESSFUL;
}


/**
 * Trim contig head and tail nodes.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateContigNodes(contigtype *contigArr, int64_t *validHeadRowContigArray, int64_t *validTailRowContigArray, int64_t *successContigIndex, int64_t *contigNodesNum)
{
	int32_t i, newContigIndex;

	// update head nodes
	for(i=0; i<*contigNodesNum; i++)
	{
		if(contigArr[i].ridposnum>0)
		{
			*validHeadRowContigArray = i;
			break;
		}
	}

	// update contig nodes
	if((*validHeadRowContigArray)>0)
	{
		newContigIndex = 1;
		for(i=*validHeadRowContigArray; i<*contigNodesNum; i++, newContigIndex++)
			contigArr[i].index = newContigIndex;
	}

	// update tail nodes
	for(i=(*successContigIndex)-1; i<(*contigNodesNum); i++)
	{
		if(contigArr[i].ridposnum>0)
		{
			free(contigArr[i].pridposorientation);
			contigArr[i].pridposorientation = NULL;

			printf("line=%d, In %s, ridpos>0, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	*contigNodesNum = (*successContigIndex) - (*validHeadRowContigArray);
	*validTailRowContigArray = (*successContigIndex) - 1;
	*successContigIndex -= *validHeadRowContigArray;

	return SUCCESSFUL;
}
/**
 * Trim contig nodes before second round assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short trimContigBeforeCycle2(contigtype *contigArr, int64_t *successContigIndex, int64_t *contigNodeNum)
{
	int64_t i, row, newIndex, newContigHeadIndex;

	if(getNewHeadContigIndex(&newContigHeadIndex, contigArr, *contigNodeNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the new head contig node, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//if(newContigHeadIndex<=0)
	//{
	//	printf("line=%d, In %s(), cannot get the new head contig node, error!\n", __LINE__, __func__);
	//	return FAILED;
	//}

	newIndex = 1;
	for(i=newContigHeadIndex-1; i<*contigNodeNum; i++, newIndex++)
		contigArr[i].index = newIndex;

	if(newContigHeadIndex>1)
	{
		row = 0;
		for(i=newContigHeadIndex-1; i<*contigNodeNum; i++, row++)
		{
			if(memcpy(contigArr+row, contigArr+i, sizeof(contigtype))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	*contigNodeNum = (*successContigIndex) - newContigHeadIndex + 1;
	*successContigIndex -= newContigHeadIndex - 1;

	return SUCCESSFUL;
}

/**
 * Get the first successful contig nodes.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNewHeadContigIndex(int64_t *contigHeadIndexTmp, contigtype *contigArr, int64_t contigNodesNum)
{
	successRead_t *ridposorient = NULL;
	short successFlag;
	int i = 0, j, numridposorient = 0, this_successContigIndex = INT_MAX;
	int32_t maxMatchLen, endIndex, minSuccessContigIndex;

	if(contigNodesNum>MAX_READ_LEN_IN_BUF)
		endIndex = MAX_READ_LEN_IN_BUF;
	else
		endIndex = contigNodesNum;

	successFlag = NO;
	minSuccessContigIndex = INT_MAX;
	for(i=0; i<endIndex; i++)
	{
		if(contigArr[i].ridposnum > 0)
		{

			successFlag = YES;

			maxMatchLen = 0;
			ridposorient = contigArr[i].pridposorientation;
			numridposorient = contigArr[i].ridposnum;
			for(j=0; j<numridposorient; j++)
			{
				if(maxMatchLen < ridposorient[j].matchlen)
				{
					maxMatchLen = ridposorient[j].matchlen;
				}
			}

			this_successContigIndex = i+1 - maxMatchLen + 1;

			// ############################ Debug information ##############################
			if(this_successContigIndex<=0)
			{
				//printf("line=%d, In %s(), this_successContigIndex=%d <=0, error!\n", __LINE__, __func__, this_successContigIndex);
				//return FAILED;
				this_successContigIndex = 1;
			}
			// ############################ Debug information ##############################

			if(this_successContigIndex < minSuccessContigIndex)
			{
				minSuccessContigIndex = this_successContigIndex;
			}
		}
	}

	if(successFlag==YES)
		*contigHeadIndexTmp = minSuccessContigIndex;
	else
		*contigHeadIndexTmp = -1;

	return SUCCESSFUL;
}

/**
 * Get the k-mers for second assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getSecondAssemblyFirstKmers(contigtype *contigArr, int64_t contigNodesNum, graphtype *graph)
{
	int i, j, startRow;

	for(i=0; i<graph->entriesPerKmer; i++)	kmerSeqIntAssembly[i] = 0;

	j = 0;
	startRow = contigNodesNum - readLen;
	if(startRow<0)
		startRow = 0;
	for(i=1; i<=kmerSize; i++, startRow++)
	{
		kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | contigArr[startRow].base;
		if(i%32==0)
			j++;
	}

	kmers[0] = getKmer(kmerSeqIntAssembly, graph);
	kmers[1] = getReverseKmer(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, graph);

	naviSuccessFlag = NAVI_SUCCESS;

	//if(kmers[0]==NULL && kmers[1]==NULL)
	//	return FAILED;

	return SUCCESSFUL;
}

/**
 * Reverse the contig nodes.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short reverseContig(contigtype *contigArr, int64_t contigNodesNum)
{
	contigtype tmpContig;
	int i, startRow, endRow;

	startRow = 0;
	endRow = contigNodesNum - 1;
	while(startRow<endRow)
	{
		if(memcpy(&tmpContig, contigArr+startRow, sizeof(contigtype))==NULL)
		{
			printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(memcpy(contigArr+startRow, contigArr+endRow, sizeof(contigtype))==NULL)
		{
			printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(memcpy(contigArr+endRow, &tmpContig, sizeof(contigtype))==NULL)
		{
			printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		startRow ++;
		endRow --;
	}

	// reset the index
	for(i=0; i<contigNodesNum; i++)
	{
		contigArr[i].index = i + 1;
		contigArr[i].base = (~(contigArr[i].base)) & 3;
	}

	return SUCCESSFUL;
}

/**
 * Initialize the decision table for second round assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initAssemblingTableSecondAssembly(contigtype *contigArr, int64_t contigNodeNum, graphtype *graph)
{
	int32_t i, j, row, tmpNodesNum;

	// get the k-mers for second round assembly
	if(getSecondAssemblyFirstKmers(contigArr, contigNodeNum, graph)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the second round assembly first kmer, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize the decision table
	itemNumDecisionTable = 0;
	if(addFirstKmerToDecisionTable(kmers)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the decision table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	row = contigNodeNum - readLen + kmerSize;
	if(row<kmerSize)
		row = kmerSize;

	lockedReadsNum = 0;
	for(tmpNodesNum=row; tmpNodesNum<contigNodeNum; tmpNodesNum++)
	{
/*
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
			}
			kmerSeqIntAssembly[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
		}
		kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | contigArr[i].base) & lastEntryMask;

		// get k-mers
		kmers[0] = getKmer(kmerSeqIntAssembly, graph);
		kmers[1] = getReverseKmer(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, graph);
*/

		// initialize or update the PE hash table
		if(PEGivenType>NONE_PE_GIVEN_TYPE && tmpNodesNum>=minContigLenUsingPE)
		{
			//if(readsNumInPEHashArr>=minReadsNumPEHashThres && regLenPEHash>=minRegLenUsingPE)
			if(readsNumInPEHashArr>=0 && regLenPEHash>=minRegLenUsingPE)
			{
				if(getNextKmerByMix(tmpNodesNum, assemblyRound)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the next kmer by mix, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}else
			{
				navigationFlag = NAVI_SE_FLAG;
				if(getNextKmerBySE(tmpNodesNum)==FAILED)
				{
					printf("line=%d, In %s(), cannot get next kmer, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}else
		{
			navigationFlag = NAVI_SE_FLAG;
			if(getNextKmerBySE(tmpNodesNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot get next kmer, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		if(naviSuccessFlag==NAVI_FAILED)
		{
			if(entriesPerKmer>=2)
			{
				for(j=0; j<entriesPerKmer-2; j++)
				{
					kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
				}
				kmerSeqIntAssembly[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
			}
			kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | contigArr[tmpNodesNum].base) & lastEntryMask;

			// get k-mers
			kmers[0] = getKmer(kmerSeqIntAssembly, graph);
			kmers[1] = getReverseKmer(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, graph);
		}else if((kmerSeqIntAssembly[entriesPerKmer-1] & 3) != contigArr[tmpNodesNum].base)
		{
			contigArr[tmpNodesNum].base = kmerSeqIntAssembly[entriesPerKmer-1] & 3;
		}

		// update the decision table according to k-mers
		if(updateDecisionTable(kmers, kmerSeqIntAssembly[entriesPerKmer-1] & 3)==FAILED)
		{
			printf("line=%d, In %s(), cannot update decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// Update the reads status in decision table
		if(updateAssemblingreadsStatus()==FAILED)
		{
			printf("line=%d, In %s(), cannot update reads status in decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// Update the locked reads and their total number
		if(updateLockedReads()==FAILED)
		{
			printf("line=%d, In %s(), cannot update locked reads, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// delete the reads from decision table
		if(removeFinishedReadsFromDecisionTable()==FAILED)
		{
			printf("line=%d, In %s(), cannot resort decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Get the reversed sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReversedSeq(char *reversed_seq, char *seq, int seq_len)
{
	short i = 0;
	for(; i<seq_len; i++)
	{
		switch(seq[i])
		{
			case 'A':
			case 'a':
				reversed_seq[seq_len-i-1] = 'T';
				break;
			case 'C':
			case 'c':
				reversed_seq[seq_len-i-1] = 'G';
				break;
			case 'G':
			case 'g':
				reversed_seq[seq_len-i-1] = 'C';
				break;
			case 'T':
			case 't':
				reversed_seq[seq_len-i-1] = 'A';
				break;
			default:
				printf("line=%d, In %s(), unknown base: %c, Error!\n", __LINE__, __func__, seq[i]);
				return FAILED;
		}
	}

	return SUCCESSFUL;
}


/**
 * Initialize the second round assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; if it cannot do second round assembly, return FAILED;
 *  	otherwise, return ERROR.
 */
short initSecondAssembly()
{
	if(successContigIndex>0)
	{
		//====================================================
		// trim a read length of contig nodes at tail
		//if(averKmerOcc>10 && itemNumContigArr>=3*readLen)  // deleted 2012-11-28
		if(trimReadLenFlag==YES)						// added 2012-11-28
		{
			if(trimContigTail(&successContigIndex, &itemNumContigArr, readLen, FIRST_ROUND_ASSEMBLY)==FAILED)
			{
				printf("line=%d, In %s(), cannot trim contig nodes at contig tail by a read length, error!\n", __LINE__, __func__);
				return ERROR;
			}
		}
		//====================================================

		if(trimContigBeforeCycle2(contigArr, &successContigIndex, &itemNumContigArr)==FAILED)  // deleted 2012-12-30
		{
			printf("line=%d, In %s(), contigsNum==%d, itemNumContigArr=%ld, cannot trim the Contig before 2nd round assembly, Error!\n", __LINE__, __func__, contigsNum+1, itemNumContigArr);
			return ERROR;
		}
	}else
	{
		return FAILED;
	}

	//reverse contig nodes
	if(reverseContig(contigArr, itemNumContigArr)==FAILED)
	{
		printf("line=%d, In %s(), contigsNum=%d, cannot reverse Contig, Error!\n", __LINE__, __func__, contigsNum+1);
		return ERROR;
	}


	if(PEGivenType>NONE_PE_GIVEN_TYPE && itemNumContigArr>=minContigLenUsingPE)
	{
		//initialize PEhashTable
		//if(initPEHashtableSecondAssembly(contigArr, itemNumContigArr-readLen+kmerSize)==FAILED)
		if(initPEHashtableSecondAssembly(contigArr, itemNumContigArr)==FAILED)
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

	// clean the DTRow hash table
	if(cleanDTRowIndexHashtable(dtRowHashtable)==FAILED)
	{
		printf("line=%d, In %s(), contigsNum=%d, cannot clean DTRow hash table, Error!\n", __LINE__, __func__, contigsNum+1);
		return ERROR;
	}

	if(initAssemblingTableSecondAssembly(contigArr, itemNumContigArr, deBruijnGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize assembly table when the second assembly, error!\n", __LINE__, __func__);
		return ERROR;
	}
	if(itemNumDecisionTable<0)
	{
		printf("line=%d, In buildContigs(), cannot init assembly table when the second assembly, error!\n", __LINE__);
		return ERROR;
	}

	// initialize the reads number region
//	if(itemNumContigArr>=minContigLenCheckingReadsNum)
//	{
		if(initReadsNumRegSecondAssembly(itemNumContigArr)==FAILED)
		{
			printf("line=%d, In %s(), cannot initialize the reads number region, error!\n", __LINE__, __func__);
			return FAILED;
		}
//	}

	if(setEmptyNaviOccQueue(naviOccQueue, &itemNumNaviOccQueue, &frontRowNaviOccQueue, &rearRowNaviOccQueue)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the empty navigation occurrence queue, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//lockedReadsNum = 0;

	if(PEGivenType>NONE_PE_GIVEN_TYPE && itemNumContigArr>=minContigLenUsingPE)
		allowedUpdatePEHashArrFlag = NO;

#if (DRAW_CURVE_FLAG==YES)
	refPosSoildFlag = refPosFirstBase[0];
	if(refPosFirstBase[1]==STRAND_PLUS)
		refStrandContig = STRAND_MINUS;
	else if(refPosFirstBase[1]==STRAND_MINUS)
		refStrandContig = STRAND_PLUS;
	refPosContig = refPosFirstBase[2];
#endif

	return SUCCESSFUL;
}

/**
 * Check the reads number in the reads number region.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateReadsNumReg(int itemNumSuccessReadsArr, int contigNodesNum, int assemblyRound)
{
	int i, contigIndexLeft, contigIndexRight;
	//double averOccNumNaviOccQueue;

	if(contigNodesNum>minContigLenCheckingReadsNum)
	{ // only slide the region

		// ############################ Debug information ########################
		if(regLenReadsNumReg<maxRegLenReadsNumReg)
		{
			printf("line=%d, In %s(), regLenReadsNumReg=%d < maxRegLenReadsNumReg=%d, error!\n", __LINE__, __func__, regLenReadsNumReg, maxRegLenReadsNumReg);
			return FAILED;
		}
		// ############################ Debug information ########################

		// get the number of reads leaving the region, and subtract the number
		//readsNumReadsNumReg -= leftContigReadsNumReg->ridposnum;
		//leftContigReadsNumReg ++;
		readsNumReadsNumReg -= contigArr[leftContigRowReadsNumReg].ridposnum;
		leftContigRowReadsNumReg ++;

		// get the number of reads entering the region, and add the number
		readsNumReadsNumReg += itemNumSuccessReadsArr;
		//rightContigReadsNumReg ++;
		rightContigRowReadsNumReg ++;

		readsNumTotal += itemNumSuccessReadsArr;


	}else if(contigNodesNum==minContigLenCheckingReadsNum)
	{ // initialize the region with the maximal size
		readsNumTotal = readsNumReadsNumReg = 0;
		regLenReadsNumReg = 0;
		if(assemblyRound==FIRST_ROUND_ASSEMBLY)
		{ // the first round
			contigIndexLeft = contigNodesNum - maxRegLenReadsNumReg + 1;
			contigIndexRight = contigNodesNum;
		}else
		{ // the second round
			contigIndexLeft = contigNodesNum - maxRegLenReadsNumReg + 1 - readLen + 1;
			contigIndexRight = contigNodesNum - readLen + 1;
		}

		if(contigIndexLeft<=0)
			contigIndexLeft = 1;
		if(contigIndexRight<contigIndexLeft)
			contigIndexRight = contigIndexRight;

		regLenReadsNumReg = contigIndexRight - contigIndexLeft + 1;

		// ############################ Debug information ########################
		if(regLenReadsNumReg!=maxRegLenReadsNumReg)
		{
			printf("line=%d, In %s(), regLenReadsNumReg=%d != maxRegLenReadsNumReg=%d, error!\n", __LINE__, __func__, regLenReadsNumReg, maxRegLenReadsNumReg);
			return FAILED;
		}
		// ############################ Debug information ########################

		//leftContigReadsNumReg = contigArr + contigIndexLeft - 1;
		//rightContigReadsNumReg = contigArr + contigIndexRight - 1;
		leftContigRowReadsNumReg = contigIndexLeft - 1;
		rightContigRowReadsNumReg = contigIndexRight - 1;
		readsNumReadsNumReg = 0;
		for(i=contigIndexLeft-1; i<contigIndexRight; i++)
			if(contigArr[i].ridposnum>0)
				readsNumReadsNumReg += contigArr[i].ridposnum;

		readsNumTotal = 0;
		for(i=0; i<contigNodesNum; i++)
			if(contigArr[i].ridposnum>0)
				readsNumTotal += contigArr[i].ridposnum;
	}

	// compute the ratio of read numbers, and decide whether the assembly should be terminated
	readsNumRatio = (double)readsNumReadsNumReg * contigNodesNum / (readsNumTotal * regLenReadsNumReg);

#if(USING_RULES==YES)
	//if(readsNumRatio>maxReadsNumRatioThres || readsNumRatio<minReadsNumRatioThres)
	if(readsNumRatio>maxReadsNumRatioThres || readsNumRatio<minReadsNumRatioThres)
	{
		//if((navigationFlag==NAVI_PE_FLAG && secondOccIndexPE>minKmerOccPE) || (navigationFlag==NAVI_SE_FLAG && secondOccIndexSE>minKmerOccSE))
		if((navigationFlag==NAVI_PE_FLAG && secondOccPE/maxOccPE>0.5) || (navigationFlag!=NAVI_PE_FLAG && secondOccSE/maxOccSE>0.5))
		//if((navigationFlag==NAVI_PE_FLAG && secondOccPE/maxOccPE>SECOND_FIRST_OCC_FAILED_RATIO) || (navigationFlag!=NAVI_PE_FLAG && secondOccSE/maxOccSE>SECOND_FIRST_OCC_FAILED_RATIO))
		{
			if(navigationFlag==NAVI_MIX_FLAG)
			{
				if(maxOccIndexPE!=maxOccIndexSE)
				{
					solvedRepeatsNum ++;
					readsNumTotal = readsNumReadsNumReg = 0;
					leftContigRowReadsNumReg = rightContigRowReadsNumReg = -1;
					kmers[0] = kmers[1] = NULL;
					naviSuccessFlag = NAVI_FAILED;
				}
			}else
			{
				solvedRepeatsNum ++;
				readsNumTotal = readsNumReadsNumReg = 0;
				leftContigRowReadsNumReg = rightContigRowReadsNumReg = -1;
				kmers[0] = kmers[1] = NULL;
				naviSuccessFlag = NAVI_FAILED;
			}
		}
		else if(successContigIndex>0 && contigNodesNum-successContigIndex > 0.4*readLen)
		{
/*
			if(calcAverOccNaviOccQueue(&averOccNumNaviOccQueue, naviOccQueue, itemNumNaviOccQueue)==FAILED)
			{
				printf("line=%d, In %s(), cannot compute the average occurrence in navigation occurrence queue, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(averOccNumNaviOccQueue < 2 && readsNumRatio < 0.2*minReadsNumRatioThres)
			{
				solvedRepeatsNum ++;
				readsNumTotal = readsNumReadsNumReg = 0;
				leftContigReadsNumReg = rightContigReadsNumReg = NULL;
				kmers[0] = kmers[1] = NULL;
				naviSuccessFlag = NAVI_FAILED;
			}
*/
/*
			else if(readsNumRatio<0.2*minReadsNumRatioThres)
			{
				// compute the maximal gap size in contig tail region
				if(computeGapSizeInContig(&tmp_gapSize, contighead, contig36, contigNodesNum, assemblyRound)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(tmp_gapSize>0.6*readLen)
				{
					kmers[0] = kmers[1] = NULL;
					naviSuccessFlag = NAVI_FAILED;
				}
			}
*/

/*
			if(navigationFlag==NAVI_MIX_FLAG)
			{
				if(((maxOccIndexPE!=-1 && maxOccPE < minKmerOccPE) || (maxOccIndexSE!=-1 && maxOccSE < minKmerOccSE)) && readsNumRatio < 0.2*minReadsNumRatioThres && averOccNumNaviOccQueue<2)
				//if(((maxOccIndexPE!=-1 && maxOccPE < minKmerOccPE) || (maxOccIndexSE!=-1 && maxOccSE < minKmerOccSE)) && readsNumRatio < 0.2*minReadsNumRatioThres)
				{
					solvedRepeatsNum ++;
					readsNumTotal = readsNumReadsNumReg = 0;
					leftContigReadsNumReg = rightContigReadsNumReg = NULL;
					kmers[0] = kmers[1] = NULL;
					naviSuccessFlag = NAVI_FAILED;
				}
			}else if(navigationFlag==NAVI_PE_FLAG)
			{
				if((maxOccIndexPE!=-1 && maxOccPE < minKmerOccPE) && readsNumRatio < 0.2*minReadsNumRatioThres && averOccNumNaviOccQueue<2)
				//if((maxOccIndexPE!=-1 && maxOccPE < minKmerOccPE) && readsNumRatio < 0.2*minReadsNumRatioThres)
				{
					solvedRepeatsNum ++;
					readsNumTotal = readsNumReadsNumReg = 0;
					leftContigReadsNumReg = rightContigReadsNumReg = NULL;
					kmers[0] = kmers[1] = NULL;
					naviSuccessFlag = NAVI_FAILED;
				}
			}else if(navigationFlag==NAVI_SE_FLAG)
			{
				if((maxOccIndexSE!=-1 && maxOccSE < minKmerOccSE) && readsNumRatio < 0.2*minReadsNumRatioThres && averOccNumNaviOccQueue<2)
				//if((maxOccIndexSE!=-1 && maxOccSE < minKmerOccSE) && readsNumRatio < 0.2*minReadsNumRatioThres)
				{
					solvedRepeatsNum ++;
					readsNumTotal = readsNumReadsNumReg = 0;
					leftContigReadsNumReg = rightContigReadsNumReg = NULL;
					kmers[0] = kmers[1] = NULL;
					naviSuccessFlag = NAVI_FAILED;
				}
			}
*/

		}
	}
#endif

	return SUCCESSFUL;
}


/**
 * Initialize the reads number region when second assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initReadsNumRegSecondAssembly(int contigNodesNum)
{
	int i, contigIndexLeft, contigIndexRight;

	if(contigNodesNum>=minContigLenCheckingReadsNum)
	{ // initialize the region with the maximal size

		readsNumTotal = readsNumReadsNumReg = 0;
		regLenReadsNumReg = 0;
		if(assemblyRound==FIRST_ROUND_ASSEMBLY)
		{ // the first round
			contigIndexLeft = contigNodesNum - maxRegLenReadsNumReg + 1;
			contigIndexRight = contigNodesNum;
		}else
		{ // the second round
			contigIndexLeft = contigNodesNum - maxRegLenReadsNumReg + 1 - readLen + 1;
			contigIndexRight = contigNodesNum - readLen + 1;
		}

		if(contigIndexLeft<=0)
			contigIndexLeft = 1;
		if(contigIndexRight<contigIndexLeft)
			contigIndexRight = contigIndexRight;

		regLenReadsNumReg = contigIndexRight - contigIndexLeft + 1;

		// ############################ Debug information ########################
		if(regLenReadsNumReg!=maxRegLenReadsNumReg)
		{
			printf("line=%d, In %s(), regLenReadsNumReg=%d != maxRegLenReadsNumReg=%d, error!\n", __LINE__, __func__, regLenReadsNumReg, maxRegLenReadsNumReg);
			return FAILED;
		}
		// ############################ Debug information ########################

		//leftContigReadsNumReg = contigArr + contigIndexLeft - 1;
		//rightContigReadsNumReg = contigArr + contigIndexRight - 1;
		leftContigRowReadsNumReg = contigIndexLeft - 1;
		rightContigRowReadsNumReg = contigIndexRight - 1;
		readsNumReadsNumReg = 0;
		for(i=contigIndexLeft-1; i<contigIndexRight; i++)
			if(contigArr[i].ridposnum>0)
				readsNumReadsNumReg += contigArr[i].ridposnum;

		readsNumTotal = 0;
		for(i=0; i<contigNodesNum; i++)
			if(contigArr[i].ridposnum>0)
				readsNumTotal += contigArr[i].ridposnum;

		// compute the ratio of read numbers, and decide whether the assembly should be terminated
		readsNumRatio = (double)readsNumReadsNumReg * contigNodesNum / (readsNumTotal * regLenReadsNumReg);
		if(readsNumRatio>maxReadsNumRatioThres || readsNumRatio<minReadsNumRatioThres)
		//if(readsNumRatio>maxReadsNumRatioThres)
		{
			readsNumTotal = readsNumReadsNumReg = 0;
			//leftContigReadsNumReg = rightContigReadsNumReg = NULL;
			leftContigRowReadsNumReg = rightContigRowReadsNumReg = -1;
			//kmers[0] = kmers[1] = NULL;
			naviSuccessFlag = NAVI_FAILED;

			//====================================================
			// trim a read length of contig nodes at tail
			if(trimReadLenFlag==YES)
			{
				if(trimContigTail(&successContigIndex, &itemNumContigArr, readLen, SECOND_ROUND_ASSEMBLY)==FAILED)
				{
					printf("line=%d, In %s(), cannot trim contig nodes at contig tail by a read length, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
			//====================================================
		}
	}else
	{
		readsNumRatio = 1;
		leftContigRowReadsNumReg = rightContigRowReadsNumReg = -1;
	}

	return SUCCESSFUL;
}


