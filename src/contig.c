
#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Build contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short buildContigs(char *contigFile, char *graphFileName, char *readMatchInfoFile)
{
	printf("\n============= Begin building contigs, please wait ... =============\n");

	struct timeval tp_start,tp_end;
	double time_used;
	gettimeofday(&tp_start,NULL);

	int32_t i, percent, percentNum, tmp_gapSize, validSuccessReadFlag;
	int64_t validReadNum;

	initFirstKmerThreshold();
	if(initMemory()==FAILED)
	{
		printf("line=%d, In %s(), cannot init the memory of two tables, error!\n", __LINE__, __func__);
		return FAILED;
	}

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

	// record the reads match information
	fpReadMatchInfo = fopen(readMatchInfoFile, "w");
	if(fpReadMatchInfo==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readMatchInfoFile);
		return FAILED;
	}

#if(DEBUG_OUTPUT==YES)
	fpContigHang = fopen(contigsFileHanging, "w");
	if(fpContigHang==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigsFileHanging);
		return FAILED;
	}
#endif

	basesNum = 0;
	localContigID = 0;
	percent = 0;
	percentNum = 0;
	successReadNum = 0;
	contigsNum = 0;
	validReadNum = deBruijnGraph->readSet->totalValidItemNumRead;
	totalErrReadNum = 0;

	correctNumCandPath = wrongNumCandPath = 0;
	naviBeforeCandPathPE = naviBeforeCandPathSE = NAVI_UNUSED;


	for(i=0; i<11; i++) errNumArr[i] = 0;

	if(initFirstKmerBounder(&lowerBoundFirstKmer, &upperBoundFirstKmer, averKmerOcc)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the bounder for first k-mers, error!\n", __LINE__, __func__);
		return FAILED;
	}

	kmerIndex = 0;
	while(kmerIndex < hashTableSize)
	{
		itemNumContigArr = 0;
		itemNumDecisionTable = 0;
		successContigIndex = -1;
		assemblyRound = FIRST_ROUND_ASSEMBLY;  // first round assembly
		lockedReadsNum = 0;
		this_successReadNum = 0;
		readsNumInPEHashArr = 0;
		regLenPEHash = 0;
		validReadOrientPEHash = -1;
		turnContigIndex = 0;
		localContigID ++;
		readsNumRatio = 1;
		naviTandFlag = NAVI_UNUSED;

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

		if(cleanReadsFromPEHashtable()==FAILED)
		{
			printf("line=%d, In %s(), cannot clean PE hash table, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// clean the contig path
		if(cleanContigPath(contigPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot clean the contig path, error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(updateContigPath(contigPath, NAVI_SE_FLAG, kmers, decisionTable, &itemNumDecisionTable, dtRowHashtable, contigArr, itemNumContigArr)==FAILED)
		{
			printf("line=%d, In %s(), cannot clean the contig path, error!\n", __LINE__, __func__);
			return FAILED;
		}

		while(naviSuccessFlag==NAVI_SUCCESS)
		{
#if (DEBUG_CONTIG_CHECK==YES)
			// ############################ Debug information ##############################
			if(localContigID==2 && itemNumContigArr>=166 && assemblyRound==FIRST_ROUND_ASSEMBLY)
			//if(localContigID==1 && itemNumContigArr>=100 && assemblyRound==FIRST_ROUND_ASSEMBLY)
			//if(localContigID==451 && itemNumContigArr>=2400 && assemblyRound==FIRST_ROUND_ASSEMBLY)
			//if(localContigID==532 && itemNumContigArr>=27700 && assemblyRound==FIRST_ROUND_ASSEMBLY)
			//if(localContigID==221 && itemNumContigArr>=2700 && assemblyRound==FIRST_ROUND_ASSEMBLY)
			//if(localContigID==1 && itemNumContigArr>=18900 && assemblyRound!=FIRST_ROUND_ASSEMBLY)
			//if(localContigID==2 && itemNumContigArr>=26150 && assemblyRound!=FIRST_ROUND_ASSEMBLY)
			//if(localContigID==11 && itemNumContigArr>=1600 && assemblyRound!=FIRST_ROUND_ASSEMBLY)
			{
				printf("localContigID=%ld, contigID=%d, itemNumContigArr=%ld, assemblyRound=%d\n", localContigID, contigsNum+1, itemNumContigArr, assemblyRound);
				outputContigPath(contigPath, YES);
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

				//if(readsNumInPEHashArr>=minReadsNumPEHashThres && regLenPEHash>=minRegLenUsingPE)
				if(readsNumInPEHashArr>0 && regLenPEHash>=minRegLenUsingPE)
				{
					if(getNextKmerByMix(itemNumContigArr, assemblyRound)==FAILED)
					{
						printf("line=%d, In %s(), localContigID=%ld, contigID=%d, itemNumContigArr=%ld, assemblyRound=%d, cannot get the next kmer by mix, error!\n", __LINE__, __func__, localContigID, contigsNum+1, itemNumContigArr, assemblyRound);
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
					//if(contigPath->naviSuccessSize<2*readLen && ((successContigIndex>0 && itemNumContigArr-successContigIndex>30) || readsNumRatio<0.3*minReadsNumRatioThres)) // 2014-01-31
					if(contigPath->naviSuccessSize<2*readLen && ((successContigIndex>0 && itemNumContigArr-successContigIndex>0.3*readLen) || readsNumRatio<0.3*minReadsNumRatioThres)) // 2014-12-23
						naviSuccessFlag = NAVI_FAILED;
					else if(itemNumContigArr<0.8*readLen && contigPath->itemNumPathItemList>=2 && contigPath->maxPathItem->supportReadsNum<3*contigPath->secPathItem->supportReadsNum) // 2014-01-19
					{
						naviSuccessFlag = NAVI_FAILED;
					}
					else
						if(secondOccSE>0)								// deleted 2013-02-26
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
				naviBeforeCandPathPE = naviAfterCandPathPE = naviBeforeCandPathSE = naviAfterCandPathSE = NAVI_UNUSED;
				naviBeforeTandPathPE = naviAfterTandPathPE = naviBeforeTandPathSE = naviAfterTandPathSE = NAVI_UNUSED;

				if(getNextKmerBySE(itemNumContigArr)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, cannot get next kmer, error!\n", __LINE__, __func__, localContigID);
					return FAILED;
				}
				for(i=0; i<4; i++) { occsNumPE[i] = 0; occsNumIndexPE[i] = -1; }


#if (SVM_NAVI==YES)
				//if((successContigIndex>0 && itemNumContigArr-successContigIndex>50) || readsNumRatio<0.3*minReadsNumRatioThres)
				//if((successContigIndex>0 && itemNumContigArr-successContigIndex>30) || readsNumRatio<0.3*minReadsNumRatioThres)
				//if(contigPath->naviSuccessSize<2*readLen && ((successContigIndex>0 && itemNumContigArr-successContigIndex>30) || readsNumRatio<0.3*minReadsNumRatioThres)) // 2014-01-31
				if(contigPath->naviSuccessSize<2*readLen && ((successContigIndex>0 && itemNumContigArr-successContigIndex>0.3*readLen) || readsNumRatio<0.3*minReadsNumRatioThres)) // 2014-12-23
					naviSuccessFlag = NAVI_FAILED;
				else if(itemNumContigArr<0.8*readLen && contigPath->itemNumPathItemList>=2 && contigPath->maxPathItem->supportReadsNum<3*contigPath->secPathItem->supportReadsNum) // 2014-01-19
				{
					naviSuccessFlag = NAVI_FAILED;
				}
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
				if(successContigIndex<0)
				{ // no successful reads, then assembly failed
					//printf("line=%d, In %s(), localContigID=%ld, contigsNum=%d, assemblyRound=%d, itemNumContigArr=%ld, the successContigIndex<=0!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
					break;
				}

#if (DEBUG_OUTPUT==YES)
				printf("localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, itemNumDecisionTable=%d, naviTandFlag=%d, validCandTandPathFlagPE=%d\n", localContigID, contigsNum+1, assemblyRound, itemNumContigArr, itemNumDecisionTable, naviTandFlag, contigPath->validCandPathTandPathFlagPE);
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
					break;
				}
			}

#if (DEBUG_CONTIG_CHECK==YES)
			// ############################ Debug information ##############################
			if(localContigID==62 && itemNumContigArr>=43550 && assemblyRound==FIRST_ROUND_ASSEMBLY)
			//if(localContigID==219 && itemNumContigArr>=16000 && assemblyRound!=FIRST_ROUND_ASSEMBLY)
			//if(localContigID==451 && itemNumContigArr>=2400 && assemblyRound==FIRST_ROUND_ASSEMBLY)
			//if(localContigID==532 && itemNumContigArr>=27700 && assemblyRound==FIRST_ROUND_ASSEMBLY)
			//if(localContigID==221 && itemNumContigArr>=2700 && assemblyRound==FIRST_ROUND_ASSEMBLY)
			//if(localContigID==1 && itemNumContigArr>=18900 && assemblyRound!=FIRST_ROUND_ASSEMBLY)
			//if(localContigID==2 && itemNumContigArr>=26150 && assemblyRound!=FIRST_ROUND_ASSEMBLY)
			//if(localContigID==11 && itemNumContigArr>=1600 && assemblyRound!=FIRST_ROUND_ASSEMBLY)
			{
				printf("localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, itemNumDecisionTable=%d, naviTandFlag=%d\n", localContigID, contigsNum+1, assemblyRound, itemNumContigArr, itemNumDecisionTable, naviTandFlag);
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

				number_of_overlap_less_than_threshold ++;

#if (DEBUG_OUTPUT==YES)
				printf("===localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, itemNumDecisionTable=%d, naviTandFlag=%d, validCandTandPathFlagPE=%d\n", localContigID, contigsNum+1, assemblyRound, itemNumContigArr, itemNumDecisionTable, naviTandFlag, contigPath->validCandPathTandPathFlagPE);
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
					break;
				}
			}
		}//end while(naviSuccessFlag==NAVI_SUCCESS)

		if(successContigIndex>0)
		{
			// ############################ Debug information ##############################
//			if(localContigID==98)
//			{
//				outputContigToTmpFile(contighead, HANGING_READ_TYPE_CONTIG_FILE);
//			}
			// ############################ Debug information ##############################

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


			if(itemNumContigArr>=minContigLen)
			{
				contigsNum ++;
				successReadNum += this_successReadNum;
				basesNum += itemNumContigArr;

				// add the contig item to contig graph, added 2013-10-28
				if(addContigItemToContigGraph(contigGraph, contigsNum, localContigID, contigArr, itemNumContigArr)==FAILED)
				{
					printf("line=%d, In %s(), cannot add contig item to contig graph, contigID=%d, localContigID=%ld, error!\n", __LINE__, __func__, contigsNum, localContigID);
					return FAILED;
				}

				// save the reads match information
				if(saveReadsMatchInfo(fpReadMatchInfo, contigsNum, contigArr, itemNumContigArr, contigAlignRegSize)==FAILED)
				{
					printf("line=%d, In %s(), cannot record the reads match information, contigID=%d, localContigID=%ld, error!\n", __LINE__, __func__, contigsNum, localContigID);
					return FAILED;
				}

#if(DEBUG_OUTPUT==YES)
				if(outputContigToFile(fpContigHang, HANGING_READ_TYPE_CONTIG_FILE, contigsNum, contigArr, itemNumContigArr)==FAILED)
				{
					printf("line=%d, In %s(), cannot output contig item to file, contigID=%d, localContigID=%ld, error!\n", __LINE__, __func__, contigsNum, localContigID);
					return FAILED;
				}
#endif
			}else
			{
				successReadNum -= this_successReadNum;
			}
		}

		// clean the contig array
		cleanContigArray(contigArr, &itemNumContigArr);

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

		// clean the contig path
		if(cleanContigPath(contigPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot clean the contig path, error!\n", __LINE__, __func__);
			return FAILED;
		}

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

	if(percent!=100)
	{
		percent = 100;
		printf("%d%%\n", percent);
		fflush(stdout);
		percentNum ++;
	}

	freeMemory();

	// close the file for the reads match information
	fclose(fpReadMatchInfo);
	fpReadMatchInfo = NULL;

#if(DEBUG_OUTPUT==YES)
	fclose(fpContigHang);
	fpContigHang = NULL;
#endif

	//printf("contigsNum=%d, basesNum=%ld\n", contigsNum, basesNum);


	// merge overlapped contigs
	if(mergeOverlappedContigs(readMatchInfoFile, contigGraph, deBruijnGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot merge overlapped contigs, error!\n", __LINE__, __func__);
		return FAILED;
	}


	// output the contigs from contig graph
	if(outputContigFromContigGraph(contigsFileFasta, contigGraph, minContigLen)==FAILED)
	{
		printf("line=%d, In %s(), cannot add contig item to contig graph, contigID=%d, localContigID=%ld, error!\n", __LINE__, __func__, contigsNum, localContigID);
		return FAILED;
	}

#if (DEBUG_PARA_PRINT==YES)
	printf("maxItemNumContigPath=%d, maxItemNumContigPathAdjusted=%d\n", maxItemNumContigPath, maxItemNumContigPathAdjusted);
	printf("successReadNum=%ld\n", successReadNum);
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
short initMemory()
{
	int32_t i;

	//longKmerSize = ceil(readLen * LONG_KMER_SIZE_FACTOR);
	//if(longKmerSize<=kmerSize)
		longKmerSize = kmerSize + ceil((readLen - kmerSize) * LONG_KMER_SIZE_FACTOR);

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

		minMatchNumSuccessRead = 0.8 * readLen;

#if(DEBUG_PARA_PRINT==YES)
		printf("averKmerOcc=%.2f\n", averKmerOcc);
		printf("longKmerSize=%d, longKmerStepSize=%d\n", longKmerSize, longKmerStepSize);
		printf("minKmerOccSE=%.2f, minKmerOccPE=%.2f\n", minKmerOccSE, minKmerOccPE);
		printf("maxSecondOcc=%.2f\n", maxSecondOcc);
		//printf("maxFirstOcc=%.2f\n", maxFirstOcc);
		printf("minLongKmerOcc=%.2f\n", minLongKmerOcc);
		printf("minReadsNumPEHashThres=%.2f\n", minReadsNumPEHashThres);
		printf("minMatchNumSuccessRead=%d\n", minMatchNumSuccessRead);
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
	minSuccessiveAppearedBaseNum = MIN_SUCCESSIVE_APPEARED_BASE_NUM;
	maxSuccessiveAppearedBaseNum = kmerSize;

	//lowOccThresNaviOccQueue = 0.4 * averKmerOcc;

#if(DEBUG_PARA_PRINT==YES)
	printf("lockedReadsNumThres=%.2f\n", lockedReadsNumThres);
	printf("maxRegLenReadsNumReg=%d\n", maxRegLenReadsNumReg);
	printf("minContigLenCheckingReadsNum=%d\n", minContigLenCheckingReadsNum);
	printf("maxReadsNumRatioThres=%.2f\n", maxReadsNumRatioThres);
	printf("minReadsNumRatioThres=%.2f\n", minReadsNumRatioThres);
	printf("maxUnmatchBaseNumPerRead=%d\n", maxUnmatchBaseNumPerRead);
	printf("maxSuccessiveAppearedBaseNum=%d\n", maxSuccessiveAppearedBaseNum);
	printf("minSuccessiveAppearedBaseNum=%d\n", minSuccessiveAppearedBaseNum);
	//printf("lowOccThresNaviOccQueue=%.2f\n", lowOccThresNaviOccQueue);
#endif

	hangingContigOutFlag = HANGING_CONTIG_OUT_FLAG;
	maxItemNumDecisionTable = TABLE_SIZE_ASSEMBLINGREAD;
	// decision table
	decisionTable = (assemblingreadtype*) calloc(maxItemNumDecisionTable, sizeof(assemblingreadtype));
	if(decisionTable==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the decision table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	itemNumContigArr = 0;
	maxItemNumContigArr = ITEM_NUM_CONTIG_ARR;
	contigArr = (contigtype *) calloc (maxItemNumContigArr, sizeof(contigtype));
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
	successReadsArr = (successRead_t*) calloc(maxItemNumSuccessReadsArr, sizeof(successRead_t));
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


	if(initMemCandPath(&candPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the candPath, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(initMemContigPath(&contigPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the contig path, error!\n", __LINE__, __func__);
		return FAILED;
	}


	for(i=0; i<4; i++) tabooSeqInt[i] = 0;
	for(i=0; i<32; i++)
	{
		tabooSeqInt[1] = (tabooSeqInt[1] << 2) | 1;
		tabooSeqInt[2] = (tabooSeqInt[2] << 2) | 2;
		tabooSeqInt[3] = (tabooSeqInt[3] <<  2) | 3;
	}

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

	releaseMemCandPath(&candPath);

	releaseMemContigPath(&contigPath);
}

/**
 * Initialize the bounder for first k-mers.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initFirstKmerBounder(double *lowerBoundFirstKmer, double *upperBoundFirstKmer, double averKmerOccNum)
{
		*lowerBoundFirstKmer = averKmerOccNum * KMER_OCC_LOWER_BOUND_FACTOR;
		*upperBoundFirstKmer = averKmerOccNum * KMER_OCC_UPPER_BOUND_FACTOR;

		if((*lowerBoundFirstKmer)>MIN_LOWER_BOUND)
			*lowerBoundFirstKmer = MIN_LOWER_BOUND;
		else if((*lowerBoundFirstKmer)<MIN_LOWER_BOUND_RESCUE)
			*lowerBoundFirstKmer = MIN_LOWER_BOUND_RESCUE;

#if(DEBUG_PARA_PRINT==YES)
	printf("lowerBoundFirstKmer=%.2f, upperBoundFirstKmer=%.2f, averKmerOcc=%.2f\n", *lowerBoundFirstKmer, *upperBoundFirstKmer, averKmerOccNum);
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
	int32_t i, j, maxNumPlus, ridMaxNumPlus, maxNumMinus, ridMaxNumMinus, kmerNumEnd5;
	int32_t posNum, rpos, seqLen, basePos, entriesNumTmp, baseNumLastEntryTmp, baseInt, entryRow, entryPos;
	ridpostype *ridpostable;

	int32_t readBlockID, rowNumInReadBlock, maxItemNumPerReadBlock, errorRegLenEnd3; // block id starts from 0
	int32_t readsNumFristKmer;

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
	readsNumFristKmer = 0;

	// get the end omitted bases
	if(kmers[0])
	{
		ridpostable = kmers[0]->ppos;
		posNum = kmers[0]->arraysize;
		for(i=0; i<posNum; i++)
		{
			rpos = ridpostable[i].pos;

			readBlockID = (ridpostable[i].rid - 1) / maxItemNumPerReadBlock;
			rowNumInReadBlock = (ridpostable[i].rid - 1) % maxItemNumPerReadBlock;
			pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;
			seqLen = pRead->seqlen;
			kmerNumEnd5 = ceil(seqLen * kmerRegLenRatioEnd5);

			//if(rpos==1 && ridpostable[i].delsign==0)
			if(rpos<=kmerNumEnd5 && ridpostable[i].delsign==0)
			{
				if(rpos-1>maxNumPlus)
				{
					maxNumPlus = rpos - 1;
					ridMaxNumPlus = ridpostable[i].rid;
				}

				readsNumFristKmer ++;
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

				readsNumFristKmer ++;
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

				if(entryRow==entriesNumTmp-1)
					baseInt = (pReadseq[entryRow] >> (2*(baseNumLastEntryTmp-entryPos-1))) & 3;
				else
					baseInt = (pReadseq[entryRow] >> (64 - 2*(entryPos+1))) & 3;

				contigArr[i].index = i + 1;
				contigArr[i].base = baseInt;
				contigArr[i].ridposnum = 0;
				contigArr[i].pridposorientation = NULL;

				// added 2013-10-14
				contigArr[i].naviFlag = NAVI_SE_FLAG;
				contigArr[i].naviTandFlag = NAVI_UNUSED;
				for(j=0; j<4; j++)
				{
					contigArr[i].occNumPE[j] = contigArr[i].occNumSE[j] = 0;
					contigArr[i].occIndexPE[j] = contigArr[i].occIndexSE[j] = -1;
				}
				contigArr[i].occNumSE[baseInt] = readsNumFristKmer;
				contigArr[i].occIndexSE[0] = baseInt;
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

				// added 2013-10-14
				contigArr[i].naviFlag = NAVI_SE_FLAG;
				contigArr[i].naviTandFlag = NAVI_UNUSED;
				for(j=0; j<4; j++)
				{
					contigArr[i].occNumPE[j] = contigArr[i].occNumSE[j] = 0;
					contigArr[i].occIndexPE[j] = contigArr[i].occIndexSE[j] = -1;
				}
				contigArr[i].occNumSE[(~baseInt)&3] = readsNumFristKmer;
				contigArr[i].occIndexSE[0] = (~baseInt) & 3;
			}

			itemNumContigArr += maxNumMinus;
		}
	}


	// add the k-mer bases
	for(i=0; i<kmerSize; i++, itemNumContigArr++)
	{
		entryRow = i >> 5;
		entryPos = i % 32;

		if(entryRow==entriesPerKmer-1)
			baseInt = (kmerSeqIntAssembly[entryRow] >> (2*(lastEntryBaseNum-entryPos-1))) & 3;
		else
			baseInt = (kmerSeqIntAssembly[entryRow] >> (64 - 2*(entryPos+1))) & 3;

		contigArr[itemNumContigArr].index = itemNumContigArr + 1;
		contigArr[itemNumContigArr].base = baseInt;
		contigArr[itemNumContigArr].ridposnum = 0;
		contigArr[itemNumContigArr].pridposorientation = NULL;

		// added 2013-10-14
		contigArr[itemNumContigArr].naviFlag = NAVI_SE_FLAG;
		contigArr[itemNumContigArr].naviTandFlag = NAVI_UNUSED;
		for(j=0; j<4; j++)
		{
			contigArr[itemNumContigArr].occNumPE[j] = contigArr[itemNumContigArr].occNumSE[j] = 0;
			contigArr[itemNumContigArr].occIndexPE[j] = contigArr[itemNumContigArr].occIndexSE[j] = -1;
		}
		contigArr[itemNumContigArr].occNumSE[baseInt] = readsNumFristKmer;
		contigArr[itemNumContigArr].occIndexSE[0] = baseInt;
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

	kmerNumTmp = readLen - kmerSize + 1;
	endKmerNumTmp = ceil(readLen * kmerRegLenRatioEnd5) + ceil (readLen * kmerRegLenRatioEnd3);
	//endKmerNumTmp = (ceil(readLen * kmerRegLenRatioEnd5)/deBruijnGraph->kmerSampleInterval + 1) + (ceil(readLen * kmerRegLenRatioEnd3)/deBruijnGraph->kmerSampleInterval + 1);

	averKmerOcc = kmerOccTmp;
	averKmerOcc *= (double)kmerNumTmp / endKmerNumTmp;
	//averKmerOcc /= 1.4;
	averKmerOcc /= 1.6;

/*
	kmerNumTmp = averReadLenInFileSample - kmerSize + 1;
	endKmerNumTmp = ceil(averReadLenInFileSample * kmerRegLenRatioEnd5) + ceil (averReadLenInFileSample * kmerRegLenRatioEnd3);
	//endKmerNumTmp = (ceil(averReadLenInFileSample * kmerRegLenRatioEnd5)/deBruijnGraph->kmerSampleInterval + 1) + (ceil(averReadLenInFileSample * kmerRegLenRatioEnd3)/deBruijnGraph->kmerSampleInterval + 1);

	averKmerOcc = kmerOccTmp;
	averKmerOcc *= (double)kmerNumTmp / endKmerNumTmp;
	averKmerOcc *= (double)readLen / averReadLenInFileSample;
	//averKmerOcc /= 1.4;
	averKmerOcc /= 1.6;
*/

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
	int32_t multiCopyFlag;


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

						// check the multiCopy flag
						if(getMultiReadCopyFlag(&multiCopyFlag, kmers)==FAILED)
						{
							printf("line=%d, In %s(), cannot update the contig path using paired-ends, error!\n", __LINE__, __func__);
							return FAILED;
						}

						if(multiCopyFlag==NO)
						{
							*kmerIndex = i;
							*firstKmer = kmer;
							naviSuccessFlag = NAVI_SUCCESS;
							return SUCCESSFUL;
						}
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

						// check the multiCopy flag
						if(getMultiReadCopyFlag(&multiCopyFlag, kmers)==FAILED)
						{
							printf("line=%d, In %s(), cannot update the contig path using paired-ends, error!\n", __LINE__, __func__);
							return FAILED;
						}

						if(multiCopyFlag==NO)
						{
							*kmerIndex = i;
							*firstKmer = kmer;
							naviSuccessFlag = NAVI_SUCCESS;
							return SUCCESSFUL;
						}
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
	uint32_t i, rpos, posNum, entriesNumTmp, baseNumLastEntryTmp, basePos, entryRow, entryPos, baseInt, seqLen;

	int32_t readBlockID, rowNumInReadBlock, maxItemNumPerReadBlock, errorRegLenEnd3; // block id starts from 0
	int32_t contigPosTmp;  // starts from 1

	readBlock_t *readBlockArr;
	read_t *pRead;
	readseqBlock_t *readseqBlockArr;
	uint64_t *pReadseq;
	int32_t returnCode, matedFlag, dtPathExtractedFlag, copyNum, bestRowInDT, mismatchNum;
	assemblingreadtype *dtReadPaired;

	mismatchNum = 1;
	dtPathExtractedFlag = NO;
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


				if(addReadToDecisionTable(ridpostable[i].rid, rpos, ORIENTATION_PLUS, NO, seqLen, pReadseq, mismatchNum, itemNumContigArr)==FAILED)
				{
					printf("line=%d, In %s(), cannot add read %lu, to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
					return FAILED;
				}

				// if the new read is has multiple copies in decision table, it will compute the correct one and remove other wrong ones
				if(getCopyNumOfReadInDecisionTable(&copyNum, ridpostable[i].rid, dtRowHashtable)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the copy number of read %lu from decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
					return FAILED;
				}

				if(copyNum>=2)
				{  // check the correct match information
					if(dtPathExtractedFlag==NO)
					{
						// get the candidate paths from single-copy reads in decision table
						if(getCandPathFromSingleCopyReadInDT(candPath, decisionTable, itemNumDecisionTable, dtRowHashtable)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the candidate paths from unique reads in decision table, error!\n", __LINE__, __func__);
							return FAILED;
						}
						dtPathExtractedFlag = YES;
					}

					// remove the incorrect row in decision table
					if(rmIncorrectRowInDecisionTable(&bestRowInDT, ridpostable[i].rid, seqLen, candPath, decisionTable, &itemNumDecisionTable, dtRowHashtable)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the correct row of a multiple-copied read in decision table, error!\n", __LINE__, __func__);
						return FAILED;
					}
				}else
				{
					bestRowInDT = itemNumDecisionTable - 1;
				}

				if(bestRowInDT==itemNumDecisionTable-1)
				{
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
							continue;

							//printf("line=%d, In %s(), contigPosTmp=%d, error!\n", __LINE__, __func__, contigPosTmp);
							//return FAILED;
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

				if(existReadWithPosInDecisionTable(ridpostable[i].rid, rpos-1, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
				{
					returnCode = validReadPair(&dtReadPaired, ridpostable[i].rid, kmerSize, seqLen, itemNumContigArr, decisionTable, dtRowHashtable);
					if(returnCode==YES)
					{
						if(dtReadPaired)
							dtReadPaired->matedFlag = YES;
						matedFlag = YES;
					}else if(returnCode==NO)
					{
						matedFlag = NO;
					}else
					{
						printf("line=%d, In %s(), cannot valid read pair, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(addReadToDecisionTable(ridpostable[i].rid, rpos, ORIENTATION_MINUS, NO, seqLen, pReadseq, mismatchNum, itemNumContigArr)==FAILED)
					{
						printf("line=%d, In %s(), cannot add read %lu, to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
						return FAILED;
					}

					// if the new read is has multiple copies in decision table, it will compute the correct one and remove other wrong ones
					if(getCopyNumOfReadInDecisionTable(&copyNum, ridpostable[i].rid, dtRowHashtable)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the copy number of read %lu from decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
						return FAILED;
					}

					if(copyNum>=2)
					{  // check the correct match information
						if(dtPathExtractedFlag==NO)
						{
							// get the candidate paths from single-copy reads in decision table
							if(getCandPathFromSingleCopyReadInDT(candPath, decisionTable, itemNumDecisionTable, dtRowHashtable)==FAILED)
							{
								printf("line=%d, In %s(), cannot get the candidate paths from unique reads in decision table, error!\n", __LINE__, __func__);
								return FAILED;
							}
							dtPathExtractedFlag = YES;
						}

						// remove the incorrect row in decision table
						if(rmIncorrectRowInDecisionTable(&bestRowInDT, ridpostable[i].rid, seqLen, candPath, decisionTable, &itemNumDecisionTable, dtRowHashtable)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the correct row of a multiple-copied read in decision table, error!\n", __LINE__, __func__);
							return FAILED;
						}
					}else
					{
						bestRowInDT = itemNumDecisionTable - 1;
					}

					if(bestRowInDT==itemNumDecisionTable-1)
					{
						if(rpos<seqLen-kmerSize+1)
						{
							entriesNumTmp = ((seqLen - 1) >> 5) + 1;
							baseNumLastEntryTmp = ((seqLen - 1) % 32) + 1;

							basePos = rpos + kmerSize - 1;
							//contigPosTmp = itemNumContigArr - readLen;
							contigPosTmp = itemNumContigArr - kmerSize - 1;
							if(contigPosTmp<0)
							{
								continue;

								//printf("line=%d, In %s(), contigPosTmp=%d, error!\n", __LINE__, __func__, contigPosTmp);
								//return FAILED;
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
		}
	}

	return SUCCESSFUL;
}

/**
 * Add a read into decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addReadToDecisionTable(uint64_t rid, int32_t rpos, int32_t orientation, int32_t matedFlag, int32_t seqLen, uint64_t *readseq, int32_t mismatchNum, int64_t itemNumContigArray)
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
	this_assemblingRead->firstContigPos = itemNumContigArray - kmerSize;
	this_assemblingRead->orientation = orientation;
	this_assemblingRead->status = ASSEMBLING_STATUS;
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
	this_assemblingRead->contigPathItem = NULL;
	this_assemblingRead->contigPathItemRead = NULL;
	this_assemblingRead->entriesNumReadseq = ((seqLen - 1) >> 5) + 1;
	this_assemblingRead->baseNumLastEentryReadseq = ((seqLen - 1) % 32) + 1;
	this_assemblingRead->kmerNumEnd5 = ceil(seqLen * kmerRegLenRatioEnd5);
	this_assemblingRead->kmerNumEnd3 = ceil(seqLen * kmerRegLenRatioEnd3);
	this_assemblingRead->mismatchNumWithContigPath = mismatchNum;

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
 * Replace a read into decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short replaceReadInDecisionTable(uint64_t rid, int32_t rpos, int32_t orientation, int32_t matedFlag, int32_t mismatchNum, assemblingreadtype *dtReadOld, int64_t itemNumContigArray, contigPath_t *contigPath)
{
	if(dtReadOld->rid!=rid)
	{
		printf("line=%d, In %s(), invalid old rid=%ld, new rid=%lu, when place the read, error!\n", __LINE__, __func__, (int64_t)dtReadOld->rid, rid);
		return FAILED;
	}

	//=======================
	if(dtReadOld->locked==YES)
		lockedReadsNum --;

	// delete read from the contig path item
	if(delReadFromContigPath(dtReadOld, contigPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot delete read from contig path item, error!\n", __LINE__, __func__);
		return FAILED;
	}

	dtReadOld->contigPathItem = NULL;
	dtReadOld->contigPathItemRead = NULL;

	if(orientation==ORIENTATION_PLUS)
		dtReadOld->firstBasePos = rpos - 1;
	else
		dtReadOld->firstBasePos = rpos + kmerSize - 2;
	dtReadOld->firstContigPos = itemNumContigArray - kmerSize;
	dtReadOld->orientation = orientation;
	dtReadOld->status = ASSEMBLING_STATUS;
	dtReadOld->delsign = 0;
	dtReadOld->reserved = 0;
	dtReadOld->locked = 0;
	dtReadOld->matedFlag = matedFlag;

	dtReadOld->successiveAppearBases = kmerSize;
	dtReadOld->successiveUnappearBases = 0;
	dtReadOld->unappearBlocksNum = 0;

	dtReadOld->matchBaseNum = kmerSize;
	dtReadOld->unmatchBaseNum = 0;
	dtReadOld->alignNum = 0;
	if(orientation==ORIENTATION_PLUS)
		dtReadOld->basePos = rpos + kmerSize - 2;
	else
		dtReadOld->basePos = rpos - 1;
	dtReadOld->lastMatchedBasePos = dtReadOld->basePos;

	dtReadOld->mismatchNumWithContigPath = mismatchNum;

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
	int32_t tmpOccIndexArray[2], itemNumTmpOccIndexArray;

	kmer_len = 0;
	maxOccIndexSE = -1;
	maxOccSE = 0;
	secondOccIndexSE = -1;
	secondOccSE = 0;

	for(i=0; i<4; i++) {occsNumSE[i] = 0; occsNumIndexSE[i] = -1;}

	for(i=0; i<4; i++)
	{
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
	}

	for(iterateID=0; iterateID<2; iterateID++)
	{
		if(iterateID==0)
		{
			successiveAppearBaseNum = maxSuccessiveAppearedBaseNum;
			//ignoreEndBaseFlag = YES;
			ignoreEndBaseFlag = NO;
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
				occsNumIndexSE[i] = -1;

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

				if(secondOccSE>0)
				{
					itemNumTmpOccIndexArray = 0;
					for(j=0; j<4; j++)
					{
						if(j!=maxOccIndexSE && j!=secondOccIndexSE && occsNumSE[j]>0)
							tmpOccIndexArray[itemNumTmpOccIndexArray++] = j;
					}

					if(itemNumTmpOccIndexArray==1)
					{
						occsNumIndexSE[2] = occsNumSE[ tmpOccIndexArray[itemNumTmpOccIndexArray-1] ];
					}else if(itemNumTmpOccIndexArray==2)
					{
						if(tmpOccIndexArray[0]>=tmpOccIndexArray[1])
						{
							occsNumIndexSE[2] = occsNumSE[ tmpOccIndexArray[0] ];
							occsNumIndexSE[3] = occsNumSE[ tmpOccIndexArray[1] ];
						}else
						{
							occsNumIndexSE[2] = occsNumSE[ tmpOccIndexArray[1] ];
							occsNumIndexSE[3] = occsNumSE[ tmpOccIndexArray[0] ];
						}
					}else if(itemNumTmpOccIndexArray>2 || itemNumTmpOccIndexArray<0)
					{
						printf("line=%d, In %s(), invalid itemNumTmpOccIndexArray=%d, error!\n", __LINE__, __func__, itemNumTmpOccIndexArray);
						return FAILED;
					}
				}
			}

			if(validKmerNum==1)
			{
				if(maxOccSE==1 && kmer_len > kmerSize)
				{
					kmer_len -= longKmerStepSize;
					if(kmer_len<kmerSize && kmer_len>kmerSize-longKmerStepSize)
						kmer_len = kmerSize;

					continue;
				}

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
				kmer_len -= longKmerStepSize;
				if(kmer_len<kmerSize && kmer_len>kmerSize-longKmerStepSize)
					kmer_len = kmerSize;

				continue;
			}
			else
			{
				if(maxOccSE<2 && kmer_len > kmerSize)
				{
					kmer_len -= longKmerStepSize;
					if(kmer_len<kmerSize && kmer_len>kmerSize-longKmerStepSize)
						kmer_len = kmerSize;

					continue;
				}else
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
			}
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
	readseqBlock_t *readseqBlockArr;
	uint64_t *pReadseq;

	int32_t i, properRow, posNum, basePos, baseInt_read;
	int32_t entryIDReadseq, entryPosReadseq, baseNumEntryReadseq;
	int32_t rid, maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3, mismatchNum;

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

#if 0
	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;
	readseqBlockArr = deBruijnGraph->readSet->readseqBlockArr;

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
					rid = rid_pos_table[i].rid;
					blockRowRead = (rid - 1) / maxItemNumPerReadBlock;
					itemRowBlockRead = (rid - 1) % maxItemNumPerReadBlock;
					pRead = readBlockArr[blockRowRead].readArr + itemRowBlockRead;
					pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;

					if(getMismatchNumWithContigPath(&mismatchNum, pReadseq, pRead->seqlen, rid_pos_table[i].pos, ORIENTATION_PLUS, itemNumDecisionTable, contigPath)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)rid_pos_table[i].rid);
						return FAILED;
					}

					if(mismatchNum<=contigPath->maxMismatchNumThres)
					{
						(*occNum) ++;
					}
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
				pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
				seqLen = pRead->seqlen;
				errorRegLenEnd3 = ceil(seqLen * ERROR_REGION_LEN_3End_FACTOR);

				if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos>seqLen-kmerSize+1-errorRegLenEnd3)
				{
					if(existReadWithPosInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
					{
						if(getMismatchNumWithContigPath(&mismatchNum, pReadseq, seqLen, rid_pos_table[i].pos, ORIENTATION_MINUS, itemNumDecisionTable, contigPath)==FAILED)
						{
							printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)rid_pos_table[i].rid);
							return FAILED;
						}

						if(mismatchNum<=contigPath->maxMismatchNumThres)
						{
							(*occNum) ++;
						}
					}
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
	readseqBlock_t *readseqBlockArr;
	uint64_t *pReadseq, rid;

	int32_t i, properRow, posNum, basePos, baseInt_read;
	int32_t entryIDReadseq, entryPosReadseq, baseNumEntryReadseq;
	int32_t maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3, mismatchNum;


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


#if 0
	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;
	readseqBlockArr = deBruijnGraph->readSet->readseqBlockArr;

	if(tmp_kmers[0])
	{
		rid_pos_table = tmp_kmers[0]->ppos;
		posNum = tmp_kmers[0]->arraysize;
		for(i=0; i<posNum; i++)
		{
			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
			{
				(*occNum) ++;

				rid = rid_pos_table[i].rid;
				blockRowRead = (rid - 1) / maxItemNumPerReadBlock;
				itemRowBlockRead = (rid - 1) % maxItemNumPerReadBlock;
				pRead = readBlockArr[blockRowRead].readArr + itemRowBlockRead;
				pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;

				if(getMismatchNumWithContigPath(&mismatchNum, pReadseq, pRead->seqlen, rid_pos_table[i].pos, ORIENTATION_PLUS, itemNumDecisionTable, contigPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)rid);
					return FAILED;
				}

				if(mismatchNum<=contigPath->maxMismatchNumThres)
				{
					(*occNum) ++;
				}
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
			pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
			seqLen = pRead->seqlen;
			errorRegLenEnd3 = ceil(seqLen * ERROR_REGION_LEN_3End_FACTOR);

			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos>seqLen-kmerSize+1-errorRegLenEnd3)
			{
				if(existReadWithPosInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
				{
					if(getMismatchNumWithContigPath(&mismatchNum, pReadseq, seqLen, rid_pos_table[i].pos, ORIENTATION_MINUS, itemNumDecisionTable, contigPath)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)rid);
						return FAILED;
					}

					if(mismatchNum<=contigPath->maxMismatchNumThres)
					{
						(*occNum) ++;
					}
				}
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
	readseqBlock_t *readseqBlockArr;
	uint64_t *pReadseq, rid;

	int32_t i, properRow, posNum, basePos, baseInt_read;
	int32_t entryIDReadseq, entryPosReadseq, baseNumEntryReadseq;
	int32_t maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3, mismatchNum;


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

#if 0
	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;
	readseqBlockArr = deBruijnGraph->readSet->readseqBlockArr;

	if(tmp_kmers[0])
	{
		rid_pos_table = tmp_kmers[0]->ppos;
		posNum = tmp_kmers[0]->arraysize;
		for(i=0; i<posNum; i++)
		{
			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
			{
				rid = rid_pos_table[i].rid;
				blockRowRead = (rid - 1) / maxItemNumPerReadBlock;
				itemRowBlockRead = (rid - 1) % maxItemNumPerReadBlock;
				pRead = readBlockArr[blockRowRead].readArr + itemRowBlockRead;
				pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;

				if(getMismatchNumWithContigPath(&mismatchNum, pReadseq, pRead->seqlen, rid_pos_table[i].pos, ORIENTATION_PLUS, itemNumDecisionTable, contigPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)rid);
					return FAILED;
				}

				if(mismatchNum<=contigPath->maxMismatchNumThres)
				{
					(*occNum) ++;
				}
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
			pReadseq = readseqBlockArr[pRead->readseqBlockID-1].readseqArr + pRead->rowReadseqInBlock;
			seqLen = pRead->seqlen;
			errorRegLenEnd3 = ceil(seqLen * ERROR_REGION_LEN_3End_FACTOR);

			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos>seqLen-kmerSize+1-errorRegLenEnd3)
			{
				if(existReadWithPosInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
				{
					if(getMismatchNumWithContigPath(&mismatchNum, pReadseq, seqLen, rid_pos_table[i].pos, ORIENTATION_MINUS, itemNumDecisionTable, contigPath)==FAILED)
					{
						printf("line=%d, In %s(), cannot get the mismatch number of read %lu with contig path, error!\n", __LINE__, __func__, (uint64_t)rid);
						return FAILED;
					}

					if(mismatchNum<=contigPath->maxMismatchNumThres)
					{
						(*occNum) ++;
					}
				}
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
short addRidposToContig(successRead_t *successReadArray, int32_t successReadNum, int32_t contigNodesNum)
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
			ridposorient[i].hangingIndex = successReadArray[i].hangingIndex = contigNodesNum;
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

			// ########################## Debug information ########################
//			if(successReadArray[0].rid==16798181)
//			{
//				printf("line=%d, In %s(), rid=%ld\n", __LINE__, __func__, (int64_t)successReadArray[0].rid);
//			}
			// ########################## Debug information ########################

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
//				if(successReadArray[i].rid==16798181)
//				{
//					printf("line=%d, In %s(), rid=%ld\n", __LINE__, __func__, (int64_t)successReadArray[i].rid);
//				}
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

	for(i=0; i<(*contigNodesNum); i++)
	{
		if(contigArr[i].ridposnum>0)
		{
			free(contigArr[i].pridposorientation);
			contigArr[i].pridposorientation = NULL;
			contigArr[i].ridposnum = 0;
		}
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
	int32_t i, j, num, indexTmp, orient;
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
				{
					switch(ridposorientation[j].orientation)
					{
						case ORIENTATION_PLUS: orient = '+'; break;
						case ORIENTATION_MINUS: orient = '-'; break;
						default: printf("line=%d, In %s(), invalid orientation %d, error!\n", __LINE__, __func__, ridposorientation[j].orientation); return FAILED;
					}
					fprintf(fpContig, "\t(%lu,%d,%d,%c)", (uint64_t)ridposorientation[j].rid, (int32_t)ridposorientation[j].pos, ridposorientation[j].matchnum, orient);
				}
				fprintf(fpContig, "\n");

			}
			//fflush(fpContig);
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
	int32_t baseInt, seqLen, entriesNum, baseNumLastEntry, entryRow, entryPos, kmerIntervalTmp;
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

			kmerIntervalTmp = 0;
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

				kmerIntervalTmp ++;
				if(kmerIntervalTmp==deBruijnGraph->kmerSampleInterval || basePos==seqLen-1)
				{
					hashcode = kmerhashInt(tmp_kmerseq);
					if(delKmerByHash(hashcode, tmp_kmerseq, rid, rpos, deBruijnGraph)==FAILED)
					{
						printf("line=%d, In %s(), can not delete the read [%lu, %c]. Error!\n", __LINE__, __func__, rid, successReadArray[i].orientation);
						return FAILED;
					}
					kmerIntervalTmp = 0;
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
					//startKmerPos = seqLen - endKmerNum - kmerSize + 1;
					//startBasePos = seqLen - endKmerNum + 1;
					startKmerPos = seqLen - kmerSize - 1 - ((endKmerNum-1)/deBruijnGraph->kmerSampleInterval)*deBruijnGraph->kmerSampleInterval + 1;
					startBasePos = startKmerPos + kmerSize;
					endBasePos = seqLen - 1;
				}

				// generate the kmer integer sequence
				if(generateKmerSeqIntFromReadset(tmp_kmerseq, readseqInt, startKmerPos, entriesNum, baseNumLastEntry)==FAILED)
				{
					printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
					return FAILED;
				}

				kmerIntervalTmp = 0;
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

					kmerIntervalTmp ++;
					if(kmerIntervalTmp==deBruijnGraph->kmerSampleInterval)
					{
						hashcode = kmerhashInt(tmp_kmerseq);
						if(delKmerByHash(hashcode, tmp_kmerseq, rid, rpos, deBruijnGraph)==FAILED)
						{
							printf("line=%d, In %s(), can not delete the read [%lu, %c]. Error!\n", __LINE__, __func__, rid, successReadArray[i].orientation);
							return FAILED;
						}
						kmerIntervalTmp = 0;
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
short updateContigtailnodes(contigtype *contigArr, int64_t successContigIndex, int64_t *contigNodesNum, int32_t assemblyRound)
{
	int64_t i, j, startContigRow, newSuccessContigIndex, readNum;
	successRead_t *successReadArray;

	if(successContigIndex!=(*contigNodesNum))
	{
		startContigRow = (*contigNodesNum) - 2 * readLen;
		if(startContigRow<0)
			startContigRow = 0;

		newSuccessContigIndex = successContigIndex;
		for(i=startContigRow; i<(*contigNodesNum); i++)
		{
			if(contigArr[i].ridposnum>0)
			{
				successReadArray = contigArr[i].pridposorientation;
				readNum = contigArr[i].ridposnum;
				for(j=0; j<readNum; j++)
					if(i+successReadArray[j].matchlen>newSuccessContigIndex)
						newSuccessContigIndex = i + successReadArray[j].matchlen;
			}
		}

		for(i=newSuccessContigIndex-1; i<(*contigNodesNum); i++)
		{
			if(contigArr[i].ridposnum>0)
			{
				free(contigArr[i].pridposorientation);
				contigArr[i].pridposorientation = NULL;
				contigArr[i].ridposnum = 0;

				//printf("line=%d, In %s(), ridpos>0, error!\n", __LINE__, __func__);
				//return FAILED;
			}
		}
		*contigNodesNum = newSuccessContigIndex;
	}

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
short trimContigBeforeCycle2(contigtype *contigArray, int64_t *successContigIndex, int64_t *contigNodeNum)
{
	int64_t i, row, newIndex, newContigHeadIndex;

	if((*successContigIndex)!=(*contigNodeNum))
	{
		if(getSuccessContigIndex(successContigIndex, contigArray, *contigNodeNum, FIRST_ROUND_ASSEMBLY)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the successContigIndex, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	if(getNewHeadContigIndex(&newContigHeadIndex, contigArray, *contigNodeNum)==FAILED)
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
		contigArray[i].index = newIndex;

	if(newContigHeadIndex>1)
	{
		row = 0;
		for(i=newContigHeadIndex-1; i<*contigNodeNum; i++, row++)
		{
			if(memcpy(contigArray+row, contigArray+i, sizeof(contigtype))==NULL)
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
 * Get the successful contig nodes.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getSuccessContigIndex(int64_t *successContigIndex, contigtype *contigArray, int64_t contigNodeNum, int32_t assemblyRound)
{
	int64_t i, j, startContigRow, endContigRow, posNum, newSuccIndex;
	successRead_t *successReadArray;

	newSuccIndex = *successContigIndex;

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{
		startContigRow = contigNodeNum - 1;
		endContigRow = contigNodeNum - 2 *readLen;
		if(endContigRow<0)
			endContigRow = 0;
		for(i=startContigRow; i>=endContigRow; i--)
		{
			if(contigArray[i].ridposnum>0)
			{
				newSuccIndex = contigArray[i].index;
				break;
			}
		}
	}else
	{
		startContigRow = contigNodeNum - 2 * readLen;
		if(startContigRow<0)
			startContigRow = 0;

		for(i=startContigRow; i<contigNodeNum; i++)
		{
			if(contigArray[i].ridposnum>0)
			{
				successReadArray = contigArray[i].pridposorientation;
				posNum = contigArray[i].ridposnum;
				for(j=0; j<posNum; j++)
					if(i+successReadArray[j].matchlen>newSuccIndex)
						newSuccIndex = i + successReadArray[j].matchlen;
			}
		}
	}

	*successContigIndex = newSuccIndex;

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
	int32_t i, j, startRow, endRow;
	uint32_t tmp;

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

		// added 2013-10-14
		if(contigArr[i].naviFlag==NAVI_PE_FLAG)
		{
			// exchange occNum of A and T
			tmp = contigArr[i].occNumPE[0];
			contigArr[i].occNumPE[0] = contigArr[i].occNumPE[3];
			contigArr[i].occNumPE[3] = tmp;
			// exchange occNum of C and G
			tmp = contigArr[i].occNumPE[1];
			contigArr[i].occNumPE[1] = contigArr[i].occNumPE[2];
			contigArr[i].occNumPE[2] = tmp;

			// exchange occIndex
			for(j=0; j<4; j++)
				if(contigArr[i].occIndexPE[j]!=-1)
					contigArr[i].occIndexPE[j] = (~(contigArr[i].occIndexPE[j])) & 3;
		}else if(contigArr[i].naviFlag==NAVI_SE_FLAG)
		{
			// exchange occNum of A and T
			tmp = contigArr[i].occNumSE[0];
			contigArr[i].occNumSE[0] = contigArr[i].occNumSE[3];
			contigArr[i].occNumSE[3] = tmp;
			// exchange occNum of C and G
			tmp = contigArr[i].occNumSE[1];
			contigArr[i].occNumSE[1] = contigArr[i].occNumSE[2];
			contigArr[i].occNumSE[2] = tmp;

			// exchange occIndex
			for(j=0; j<4; j++)
				if(contigArr[i].occIndexSE[j]!=-1)
					contigArr[i].occIndexSE[j] = (~(contigArr[i].occIndexSE[j])) & 3;
		}else if(contigArr[i].naviFlag==NAVI_MIX_FLAG)
		{
			if(contigArr[i].occIndexPE[0]!=-1)
			{
				// exchange occNum of A and T
				tmp = contigArr[i].occNumPE[0];
				contigArr[i].occNumPE[0] = contigArr[i].occNumPE[3];
				contigArr[i].occNumPE[3] = tmp;
				// exchange occNum of C and G
				tmp = contigArr[i].occNumPE[1];
				contigArr[i].occNumPE[1] = contigArr[i].occNumPE[2];
				contigArr[i].occNumPE[2] = tmp;

				// exchange occIndex
				for(j=0; j<4; j++)
					if(contigArr[i].occIndexPE[j]!=-1)
						contigArr[i].occIndexPE[j] = (~(contigArr[i].occIndexPE[j])) & 3;
			}

			if(contigArr[i].occIndexSE[0]!=-1)
			{
				// exchange occNum of A and T
				tmp = contigArr[i].occNumSE[0];
				contigArr[i].occNumSE[0] = contigArr[i].occNumSE[3];
				contigArr[i].occNumSE[3] = tmp;
				// exchange occNum of C and G
				tmp = contigArr[i].occNumSE[1];
				contigArr[i].occNumSE[1] = contigArr[i].occNumSE[2];
				contigArr[i].occNumSE[2] = tmp;

				// exchange occIndex
				for(j=0; j<4; j++)
					if(contigArr[i].occIndexSE[j]!=-1)
						contigArr[i].occIndexSE[j] = (~(contigArr[i].occIndexSE[j])) & 3;
			}else
			{
				printf("line=%d, In %s(), invalid occIndexSE[0]=%d, error!\n", __LINE__, __func__, contigArr[i].occIndexSE[0]);
				return FAILED;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Initialize the decision table for second round assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initAssemblingTableSecondAssembly(contigtype *contigArray, int64_t contigNodeNum, graphtype *graph)
{
	int32_t j, row, disagreeNum;
	int64_t tmpNodesNum;
	int32_t newContigBaseIntArray[2*readLen], newBaseNum;

	// get the k-mers for second round assembly
	if(getSecondAssemblyFirstKmers(contigArray, contigNodeNum, graph)==FAILED)
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

	// clean the contig path and update it again
	if(cleanContigPath(contigPath)==FAILED)
	{
		printf("line=%d, In %s(), cannot clean the contig path, error!\n", __LINE__, __func__);
		return FAILED;
	}
	if(updateContigPath(contigPath, NAVI_SE_FLAG, kmers, decisionTable, &itemNumDecisionTable, dtRowHashtable, contigArray, row)==FAILED)
	{
		printf("line=%d, In %s(), cannot update the contig path, error!\n", __LINE__, __func__);
		return FAILED;
	}

	disagreeNum = 0;
	lockedReadsNum = 0;
	newBaseNum = 0;
	for(tmpNodesNum=row; tmpNodesNum<contigNodeNum; tmpNodesNum++)
	{
		// initialize or update the PE hash table
		//if(PEGivenType>NONE_PE_GIVEN_TYPE && tmpNodesNum>=minContigLenUsingPE)  // deleted 2013-12-20
		if(PEGivenType>NONE_PE_GIVEN_TYPE)   // added 2013-12-20
		{
			//if(readsNumInPEHashArr>=minReadsNumPEHashThres && regLenPEHash>=minRegLenUsingPE)
			//if(readsNumInPEHashArr>=0 && regLenPEHash>=minRegLenUsingPE)  // deleted 2013-12-20
			if(readsNumInPEHashArr>=0)   // added 2013-12-20
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
			if(tmpNodesNum>=row+kmerSize)
			{
				if(navigationFlag==NAVI_SE_FLAG)
				{
					if(maxOccIndexSE==-1)
						break;
				}else if(navigationFlag==NAVI_MIX_FLAG)
				{
					if(maxOccIndexPE==-1 && maxOccIndexSE==-1)
						break;
				}
			}


			naviSuccessFlag = NAVI_SUCCESS;

			if(entriesPerKmer>=2)
			{
				for(j=0; j<entriesPerKmer-2; j++)
				{
					kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
				}
				kmerSeqIntAssembly[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
			}
			kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | contigArray[tmpNodesNum].base) & lastEntryMask;

			// get k-mers
			kmers[0] = getKmer(kmerSeqIntAssembly, graph);
			kmers[1] = getReverseKmer(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, graph);

			newContigBaseIntArray[newBaseNum++] = contigArray[tmpNodesNum].base;

		}else if((kmerSeqIntAssembly[entriesPerKmer-1] & 3) != contigArray[tmpNodesNum].base)
		{
			//printf("line=%d, In %s(), old base=%d, new base=%d\n", __LINE__, __func__, contigArray[tmpNodesNum].base, (int32_t)(kmerSeqIntAssembly[entriesPerKmer-1] & 3));
			//contigArray[tmpNodesNum].base = kmerSeqIntAssembly[entriesPerKmer-1] & 3;
			newContigBaseIntArray[newBaseNum++] = kmerSeqIntAssembly[entriesPerKmer-1] & 3;
			disagreeNum ++;
		}else
		{
			newContigBaseIntArray[newBaseNum++] = contigArray[tmpNodesNum].base;
		}

/*
#if (DEBUG_CONTIG_CHECK==YES)
		if(localContigID==189)
		{
			printf("localContigID=%ld, contigID=%d, assemblyRound=%d, tmpNodesNum=%ld, itemNumDecisionTable=%d\n", localContigID, contigsNum+1, assemblyRound, tmpNodesNum, itemNumDecisionTable);
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
				printf("\tdistance=%ld, readsNumRatio=%.2f\n", tmpNodesNum-successContigIndex, readsNumRatio);

			//outputContigPath(contigPath);
		}
#endif
*/
		if(contigNodeNum<minContigLen && contigPath->itemNumPathItemList>=2)
		{
			if(navigationFlag==NAVI_PE_FLAG && (double)occsNumPE[occsNumIndexPE[1]]/occsNumPE[occsNumIndexPE[0]]>0.6)
			{
				naviSuccessFlag = NAVI_FAILED;
				break;
			}else if(navigationFlag==NAVI_SE_FLAG && (double)occsNumSE[occsNumIndexSE[1]]/occsNumSE[occsNumIndexSE[0]]>0.6)
			{
				naviSuccessFlag = NAVI_FAILED;
				break;
			}
			else if((occsNumSE[occsNumIndexSE[0]]>0 && (double)occsNumSE[occsNumIndexSE[1]]/occsNumSE[occsNumIndexSE[0]]>0.6) || (occsNumPE[occsNumIndexPE[0]]>0 && (double)occsNumPE[occsNumIndexPE[1]]/occsNumPE[occsNumIndexPE[0]]>0.6))
			{
				naviSuccessFlag = NAVI_FAILED;
				break;
			}
		}


		// update the decision table according to k-mers
		if(updateDecisionTable(kmers, kmerSeqIntAssembly[entriesPerKmer-1] & 3)==FAILED)
		{
			printf("line=%d, In %s(), cannot update decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// add the read number
		//if(addOccsToContig(contigArray, tmpNodesNum)==FAILED)
		//{
		//	printf("line=%d, In %s(), localContigID=%ld, cannot add occ number to contig, error!\n", __LINE__, __func__, localContigID);
		//	return FAILED;
		//}

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
		if(removeFinishedReadsFromDecisionTable(decisionTable, &itemNumDecisionTable)==FAILED)
		{
			printf("line=%d, In %s(), cannot resort decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// update the contig path
		if(updateContigPath(contigPath, navigationFlag, kmers, decisionTable, &itemNumDecisionTable, dtRowHashtable, contigArray, tmpNodesNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot update the contig path, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	if(disagreeNum>=3)
	{
		naviSuccessFlag = NAVI_FAILED;
	}else
	{
		if(naviSuccessFlag==NAVI_SUCCESS)
		{
			for(j=0; j<newBaseNum; j++)
				contigArray[row+j].base = newContigBaseIntArray[j];
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
		if(initPEHashtableSecondAssembly(contigArr, itemNumContigArr, YES)==FAILED)
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

	return SUCCESSFUL;
}

/**
 * Check the reads number in the reads number region.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateReadsNumReg(int32_t itemNumSuccessReadsArr, int32_t contigNodesNum, int32_t assemblyRound)
{
	int32_t i, contigIndexLeft, contigIndexRight;
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

	return SUCCESSFUL;
}


/**
 * Initialize the reads number region when second assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initReadsNumRegSecondAssembly(int32_t contigNodesNum)
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
	}else
	{
		readsNumRatio = 1;
		leftContigRowReadsNumReg = rightContigRowReadsNumReg = -1;
	}

	return SUCCESSFUL;
}

/**
 * Add the occ number of reads to contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addOccsToContig(contigtype *contigArray, int32_t contigNodesNum, int32_t naviTandFlag, int32_t newCandBaseNumAfterTandPathPE, contigPath_t *contigPath)
{
	int32_t i, maxPathLen, secPathLen, pathLen, baseIndex, bothNullUnequalFlag;
	char *maxPathseq, *secPathseq, *pathseq;

	if(contigPath->itemNumPathItemList==1)
		contigPath->sameBaseMaxSecContigPathItem = YES;
	else if(contigPath->itemNumPathItemList==2)
	{
		if(contigPath->maxPathItem->contigPathLen-contigPath->startRowNewBase>0 && contigPath->secPathItem->contigPathLen-contigPath->startRowNewBase>0
			&& contigPath->maxPathItem->contigPathStr[contigPath->startRowNewBase]==contigPath->secPathItem->contigPathStr[contigPath->startRowNewBase])
		{
			contigPath->sameBaseMaxSecContigPathItem = YES;
		}
	}else
	{
		contigPath->sameBaseMaxSecContigPathItem = NO;
	}

	bothNullUnequalFlag = NO;
	if(contigPath->itemNumPathItemList>=2 && contigPath->maxPathItem->supportReadsNum==0 && contigPath->maxPathItem->supportReadsNum==0)
	{
		maxPathseq = contigPath->maxPathItem->contigPathStr + contigPath->startRowNewBase;
		maxPathLen = contigPath->maxPathItem->contigPathLen - contigPath->startRowNewBase;
		secPathseq = contigPath->secPathItem->contigPathStr + contigPath->startRowNewBase;
		secPathLen = contigPath->secPathItem->contigPathLen - contigPath->startRowNewBase;
		if(maxPathLen>0 && secPathLen>0 && maxPathseq[0]==secPathseq[0])
		{
			switch(maxPathseq[0])
			{
				case 'A': baseIndex = 0; break;
				case 'C': baseIndex = 1; break;
				case 'G': baseIndex = 2; break;
				case 'T': baseIndex = 3; break;
				default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, maxPathseq[0]); return FAILED;
			}

			if(baseIndex!=contigArray[contigNodesNum-1].base)
				bothNullUnequalFlag = YES;
		}
	}
	contigArray[contigNodesNum-1].bothNullUnequalFlag = bothNullUnequalFlag;

	contigArray[contigNodesNum-1].itemNumContigPath = contigPath->itemNumPathItemList;


	// add the read number, added 2013-10-14
	if(navigationFlag==NAVI_PE_FLAG)
	{
		contigArray[contigNodesNum-1].naviFlag = NAVI_PE_FLAG;
		contigArray[contigNodesNum-1].naviTandFlag = naviTandFlag;
		contigArray[contigNodesNum-1].newCandBaseNumAfterTandPathPE = newCandBaseNumAfterTandPathPE;
		contigArray[contigNodesNum-1].naviSuccessSize = contigPath->naviSuccessSize;
		contigArray[contigNodesNum-1].sameBaseMaxSecContigPathItem = contigPath->sameBaseMaxSecContigPathItem;
		for(i=0; i<4; i++)
		{
			contigArray[contigNodesNum-1].occNumPE[i] = occsNumPE[i];
			contigArray[contigNodesNum-1].occIndexPE[i] = occsNumIndexPE[i];
			contigArray[contigNodesNum-1].occNumSE[i] = 0;
			contigArray[contigNodesNum-1].occIndexSE[i] = -1;
		}
	}else if(navigationFlag==NAVI_SE_FLAG)
	{
		contigArray[contigNodesNum].naviFlag = NAVI_SE_FLAG;
		contigArray[contigNodesNum-1].naviTandFlag = naviTandFlag;
		contigArray[contigNodesNum-1].newCandBaseNumAfterTandPathPE = -1;
		contigArray[contigNodesNum-1].naviSuccessSize = contigPath->naviSuccessSize;
		contigArray[contigNodesNum-1].sameBaseMaxSecContigPathItem = contigPath->sameBaseMaxSecContigPathItem;
		for(i=0; i<4; i++)
		{
			contigArray[contigNodesNum-1].occNumPE[i] = 0;
			contigArray[contigNodesNum-1].occIndexPE[i] = -1;
			contigArray[contigNodesNum-1].occNumSE[i] = occsNumSE[i];
			contigArray[contigNodesNum-1].occIndexSE[i] = occsNumIndexSE[i];
		}
	}else if(navigationFlag==NAVI_MIX_FLAG)
	{
		contigArray[contigNodesNum-1].naviFlag = NAVI_MIX_FLAG;
		contigArray[contigNodesNum-1].naviTandFlag = naviTandFlag;
		contigArray[contigNodesNum-1].newCandBaseNumAfterTandPathPE = newCandBaseNumAfterTandPathPE;
		contigArray[contigNodesNum-1].naviSuccessSize = contigPath->naviSuccessSize;
		contigArray[contigNodesNum-1].sameBaseMaxSecContigPathItem = contigPath->sameBaseMaxSecContigPathItem;
		for(i=0; i<4; i++)
		{
			contigArray[contigNodesNum-1].occNumPE[i] = occsNumPE[i];
			contigArray[contigNodesNum-1].occIndexPE[i] = occsNumIndexPE[i];
			contigArray[contigNodesNum-1].occNumSE[i] = occsNumSE[i];
			contigArray[contigNodesNum-1].occIndexSE[i] = occsNumIndexSE[i];
		}
	}else
	{
		printf("line=%d, In %s(), unknown navigationFlag=%d, error!\n", __LINE__, __func__, navigationFlag);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get the valid success read flag.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getValidSuccessReadFlag(int32_t *validSuccessReadFlag, successRead_t *successReadsArray, int32_t itemNumSuccessReadsArray, int32_t minMatchNum)
{
	int32_t i;

	*validSuccessReadFlag = NO;
	for(i=0; i<itemNumSuccessReadsArray; i++)
	{
		if(successReadsArray[i].matchnum>=minMatchNum)
		{
			*validSuccessReadFlag = YES;
			break;
		}
	}

	return SUCCESSFUL;
}

/**
 * Confirm navigation by one single-end read.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short confirmSingleEndNavi(int32_t *thisNaviSuccessFlag, contigtype *contigArray, int32_t itemNumContigArray)
{
	// added 2013-10-14
	int32_t i, startPos, endPos, allOneSingleEnd, allWithoutPE, lastPEPos, multiExtFlag;
	int32_t lowRatioNum, lowRatioNum2, maxValue, secValue, maxIndex, secIndex, largeGapNum, gapLen;
	double averOccBeforeSE, averOccSE, baseNumAT, baseRatioAT;
	int32_t baseInt, pathLen, candPathLen, maxPathLen, secPathLen, shareLen, baseIntMax, baseIntSec, shiftSize, shiftType, difNumLocal, difNumLocal2;
	char *pathseq, *candPathseq, *maxPathseq, *secPathseq;
	int32_t baseNumArray[5], baseNumIndexArray[4], kmerBaseInt;
	double maxBaseRatio;

	// allow one single-end to extend less than MAX_NAVI_LEN_SINGLE_READ
	startPos = itemNumContigArray - MAX_NAVI_LEN_SINGLE_READ;
	if(startPos<0)
		startPos = 0;

	allOneSingleEnd = YES;
	for(i=startPos; i<itemNumContigArray; i++)
	{
		if(contigArray[i].naviFlag==NAVI_PE_FLAG)
		{
			allOneSingleEnd = NO;
			break;
		}else
		{
			if(contigArray[i].occNumSE[ contigArray[i].occIndexSE[0] ]>1)
			{
				allOneSingleEnd = NO;
				break;
			}
		}
	}

	if(allOneSingleEnd==YES)
	{
		*thisNaviSuccessFlag = NAVI_FAILED;
		return SUCCESSFUL;
	}

	// allow single-end to extend contig without paired-end less than MAX_NAVI_LEN_WITHOUT_PE
	//if((*thisNaviSuccessFlag)==NAVI_SUCCESS && (successContigIndex>0 && itemNumContigArr-successContigIndex>15) && itemNumContigArray>5*minContigLenUsingPE)
	//if(PEGivenType>NONE_PE_GIVEN_TYPE && (*thisNaviSuccessFlag)==NAVI_SUCCESS && (successContigIndex>0 && itemNumContigArr-successContigIndex>15) && (itemNumContigArray>5*minContigLenUsingPE+turnContigIndex))
	if(PEGivenType>NONE_PE_GIVEN_TYPE && (*thisNaviSuccessFlag)==NAVI_SUCCESS && (itemNumContigArray>5*minContigLenUsingPE+turnContigIndex))
	{
		endPos = itemNumContigArray - 1 - MAX_NAVI_LEN_WITHOUT_PE;
		if(endPos<0)
			endPos = 0;

		allWithoutPE = YES;
		for(i=itemNumContigArray-1; i>=endPos; i--)
		{
			if(contigArray[i].naviFlag==NAVI_PE_FLAG)
			{
				allWithoutPE = NO;
				break;
			}else
			{
				if(contigArray[i].occIndexPE[0]!=-1)
				{
					allWithoutPE = NO;
					break;
				}
			}
		}

		if(allWithoutPE==YES)
		{
			// get lastPEPos
			lastPEPos = 0;
			for(i=endPos-1; i>=0; i--)
			{
				if(contigArray[i].naviFlag==NAVI_PE_FLAG)
				{
					lastPEPos = i;
					break;
				}else if(contigArray[i].occIndexPE[0]!=-1)
				{
					lastPEPos = i;
					break;
				}
			}

			// get averOcc1
			startPos = lastPEPos - 2 * readLen;
			if(startPos<0)
				startPos = 0;
			averOccBeforeSE = 0;
			for(i=0; i<readLen; i++)
			{
				if(contigArray[i+startPos].naviFlag==NAVI_PE_FLAG)
				{
					averOccBeforeSE += contigArray[i+startPos].occNumPE[contigArray[i+startPos].occIndexPE[0]];
				}else if(contigArray[i].occIndexSE[0]!=-1)
				{
					averOccBeforeSE += contigArray[i+startPos].occNumSE[contigArray[i+startPos].occIndexSE[0]];
				}
			}
			averOccBeforeSE /= readLen;

			// get averOcc2
			averOccSE = 0;
			startPos = lastPEPos + 1;
			if(startPos<0)
				startPos = 0;
			endPos = startPos + 10;
			for(i=startPos; i<=endPos; i++)
				averOccSE += contigArray[i].occNumSE[contigArray[i].occIndexSE[0]];
			averOccSE /= 10;

			if(averOccBeforeSE/averOccSE>4)
			{
				*thisNaviSuccessFlag = NAVI_FAILED;
				//printf(" ******** line=%d, In %s(), localContigID=%ld, itemNumContigArr=%ld, turnContigIndex=%d!\n", __LINE__, __func__, localContigID, itemNumContigArr, turnContigIndex);
				return SUCCESSFUL;
			}


			// if there are some bases that secOcc>=3 && secOcc/maxOcc>0.2, stop the extension
			startPos = itemNumContigArray - 1 - MAX_NAVI_LEN_WITHOUT_PE;
			if(startPos<0)
				startPos = 0;
			lowRatioNum = 0;
			for(i=startPos; i<itemNumContigArray; i++)
			{
				if(contigArray[i].occNumSE[contigArray[i].occIndexSE[0]]>=3 && (double)contigArray[i].occNumSE[contigArray[i].occIndexSE[1]]/contigArray[i].occNumSE[contigArray[i].occIndexSE[0]]>0.2)
					lowRatioNum ++;
				else if(contigArray[i].occIndexSE[1]!=-1 && (double)contigArray[i].occNumSE[contigArray[i].occIndexSE[1]]/contigArray[i].occNumSE[contigArray[i].occIndexSE[0]]>0.3)
					lowRatioNum ++;
			}

			if(lowRatioNum>=1)
			{
				*thisNaviSuccessFlag = NAVI_FAILED;
				//printf(" ******** line=%d, In %s(), localContigID=%ld, itemNumContigArr=%ld, turnContigIndex=%d!\n", __LINE__, __func__, localContigID, itemNumContigArr, turnContigIndex);
				return SUCCESSFUL;
			}

		}
	}

	//if(PEGivenType>NONE_PE_GIVEN_TYPE && (*thisNaviSuccessFlag)==NAVI_SUCCESS && (itemNumContigArray>5*minContigLenUsingPE+turnContigIndex))
	//if(PEGivenType>NONE_PE_GIVEN_TYPE && (*thisNaviSuccessFlag)==NAVI_SUCCESS && (itemNumContigArray>3*minContigLenUsingPE+turnContigIndex)) // 2014-01-22, contig 800, 893
	if(PEGivenType>NONE_PE_GIVEN_TYPE && (*thisNaviSuccessFlag)==NAVI_SUCCESS && (itemNumContigArray>minContigLenCheckGap+turnContigIndex)) // 2014-01-24
	{
		//endPos = itemNumContigArray - 1 - MAX_NAVI_LEN_WITHOUT_PE * 2;
		endPos = itemNumContigArray - 1 - MAX_NAVI_LEN_WITHOUT_PE * 1.3;
		if(endPos<0)
			endPos = 0;

		allWithoutPE = YES;
		for(i=itemNumContigArray-1; i>=endPos; i--)
		{
			if(contigArray[i].naviFlag==NAVI_PE_FLAG)
			{
				allWithoutPE = NO;
				break;
			}else
			{
				if(contigArray[i].occIndexPE[0]!=-1)
				{
					allWithoutPE = NO;
					break;
				}
			}
		}

		if(allWithoutPE==YES)
		{
			*thisNaviSuccessFlag = NAVI_FAILED;
			//printf(" ******** line=%d, In %s(), localContigID=%ld, itemNumContigArr=%ld, turnContigIndex=%d!\n", __LINE__, __func__, localContigID, itemNumContigArr, turnContigIndex);
			return SUCCESSFUL;
		}
	}

/*
	// limit the number of multi-extensions in a region for paired-ends
	if((*thisNaviSuccessFlag)==NAVI_SUCCESS && itemNumContigArray>2*minContigLenUsingPE)
	{
		startPos = itemNumContigArray - 30;
		if(startPos<0)
			startPos = 0;
		multiBaseNnum = 0;
		for(i=startPos; i<itemNumContigArray; i++)
		{
			if(contigArray[i].occIndexPE[1]!=-1)
			{
				maxIndex = contigArray[i].occIndexPE[0];
				secIndex = contigArray[i].occIndexPE[1];
				maxValue = contigArray[i].occNumPE[maxIndex];
				secValue = contigArray[i].occNumPE[secIndex];

				if(secValue>=2)
					multiBaseNnum ++;
			}else if(contigArray[i].occIndexSE[1]!=-1)
			{
				maxIndex = contigArray[i].occIndexSE[0];
				secIndex = contigArray[i].occIndexSE[1];
				maxValue = contigArray[i].occNumSE[maxIndex];
				secValue = contigArray[i].occNumSE[secIndex];

				if(secValue>=2)
					multiBaseNnum ++;
			}
		}

		if(multiBaseNnum>=5)
		{
			*thisNaviSuccessFlag = NAVI_FAILED;
			return SUCCESSFUL;
		}
	}
*/

	// limit low ratio of paired ends
	if(PEGivenType>NONE_PE_GIVEN_TYPE && (*thisNaviSuccessFlag)==NAVI_SUCCESS && itemNumContigArray>2*minContigLenUsingPE)
	{
		startPos = itemNumContigArray - readLen;
		if(startPos<0)
			startPos = 0;
		lowRatioNum = 0;
		for(i=startPos; i<itemNumContigArray; i++)
		{
			if(contigArray[i].occIndexPE[1]!=-1)
			{
				maxIndex = contigArray[i].occIndexPE[0];
				secIndex = contigArray[i].occIndexPE[1];
				maxValue = contigArray[i].occNumPE[maxIndex];
				secValue = contigArray[i].occNumPE[secIndex];

				if((double)secValue/maxValue>=0.6)
					lowRatioNum ++;
			}
		}

		if(lowRatioNum>=3)
		{
			*thisNaviSuccessFlag = NAVI_FAILED;
			return SUCCESSFUL;
		}
/*
		// check the low ratio and large gap of paired-ends, 2014-01-06
		if(readsNumRatio<0.3 && (successContigIndex>0 && itemNumContigArr-successContigIndex>5))
		{
			// check the low ratio
			startPos = itemNumContigArray - 1.5*readLen;
			if(startPos<0)
				startPos = 0;
			lowRatioNum = 0;
			for(i=startPos; i<itemNumContigArray; i++)
			{
				if(contigArray[i].occIndexPE[1]!=-1)
				{
					maxIndex = contigArray[i].occIndexPE[0];
					secIndex = contigArray[i].occIndexPE[1];
					maxValue = contigArray[i].occNumPE[maxIndex];
					secValue = contigArray[i].occNumPE[secIndex];

					if((double)secValue/maxValue>=0.5 && maxValue-secValue<2)
						lowRatioNum ++;
				}
			}

			// check the gap region numbers, 2014-01-06
			startPos = itemNumContigArray - 1.5*readLen;
			if(startPos<readLen)
				startPos = readLen;
			largeGapNum = 0;

			if(contigArray[startPos].ridposnum>0)
				gapLen = 0;
			else
				gapLen = 1;
			for(i=startPos+1; i<itemNumContigArray; i++)
			{
				if(contigArray[i].ridposnum==0)
				{
					gapLen ++;
				}else
				{
					if(gapLen>25)
						largeGapNum ++;
					gapLen = 0;
				}
			}

			if(gapLen>25)
				largeGapNum ++;

			if(lowRatioNum>=1 && largeGapNum>=1)
			{
				*thisNaviSuccessFlag = NAVI_FAILED;
				return SUCCESSFUL;
			}
		}
*/
	}

	// check the long gapLen region, added 2013-11-09
	if(PEGivenType>NONE_PE_GIVEN_TYPE && (*thisNaviSuccessFlag)==NAVI_SUCCESS && (successContigIndex>0 && itemNumContigArr-successContigIndex>40) && readsNumRatio<0.6)
	{
		startPos = itemNumContigArray - readLen;
		if(startPos<0)
			startPos = 0;
		lowRatioNum = 0;
		for(i=startPos; i<itemNumContigArray; i++)
		{
			if(contigArray[i].occIndexPE[1]!=-1)
			{
				maxIndex = contigArray[i].occIndexPE[0];
				secIndex = contigArray[i].occIndexPE[1];
				maxValue = contigArray[i].occNumPE[maxIndex];
				secValue = contigArray[i].occNumPE[secIndex];

				if((double)secValue/maxValue>=0.5 && maxValue-secValue<2)
					lowRatioNum ++;
			}

			if(lowRatioNum>=1)
			{
				*thisNaviSuccessFlag = NAVI_FAILED;
				return SUCCESSFUL;
			}
		}

		// check the gap region numbers
		startPos = itemNumContigArray - 2*readLen;
		if(startPos<readLen)
			startPos = readLen;
		largeGapNum = 0;

		if(contigArray[startPos].ridposnum>0)
			gapLen = 0;
		else
			gapLen = 1;
		for(i=startPos+1; i<itemNumContigArray; i++)
		{
			if(contigArray[i].ridposnum==0)
			{
				gapLen ++;
			}else
			{
				if(gapLen>=25)
					largeGapNum ++;
				gapLen = 0;
			}
		}

		if(gapLen>=25)
			largeGapNum ++;

		if(largeGapNum>=3)
		//if(lowRatioNum>=1 && largeGapNum>=3) // 2014-03-12
		{
			*thisNaviSuccessFlag = NAVI_FAILED;
			return SUCCESSFUL;
		}
	}

	// 2013-12-30
	if(PEGivenType>NONE_PE_GIVEN_TYPE && (*thisNaviSuccessFlag)==NAVI_SUCCESS && (successContigIndex>0 && itemNumContigArr-successContigIndex>40) && (readsNumRatio>2 || readsNumRatio<0.3))
	{
		startPos = itemNumContigArray - 2*readLen;
		if(startPos<0)
			startPos = 0;
		lowRatioNum = 0;
		for(i=startPos; i<itemNumContigArray; i++)
		{
			if(contigArray[i].occIndexPE[1]!=-1)
			{
				maxIndex = contigArray[i].occIndexPE[0];
				secIndex = contigArray[i].occIndexPE[1];
				maxValue = contigArray[i].occNumPE[maxIndex];
				secValue = contigArray[i].occNumPE[secIndex];

				if((double)secValue/maxValue>=0.5)
					lowRatioNum ++;
			}
		}

		// check the gap region numbers
		startPos = itemNumContigArray - 2*readLen;
		if(startPos<readLen)
			startPos = readLen;
		largeGapNum = 0;

		if(contigArray[startPos].ridposnum>0)
			gapLen = 0;
		else
			gapLen = 1;
		for(i=startPos+1; i<itemNumContigArray; i++)
		{
			if(contigArray[i].ridposnum==0)
			{
				gapLen ++;
			}else
			{
				if(gapLen>=25)
					largeGapNum ++;
				gapLen = 0;
			}
		}

		if(gapLen>=25)
			largeGapNum ++;

		if(lowRatioNum>=2 && largeGapNum>=2)
		{
			*thisNaviSuccessFlag = NAVI_FAILED;
			return SUCCESSFUL;
		}
	}

	// check the large gap, 2014-01-15
	if(PEGivenType>NONE_PE_GIVEN_TYPE && (*thisNaviSuccessFlag)==NAVI_SUCCESS && (successContigIndex>0 && itemNumContigArr-successContigIndex>50) && (readsNumRatio<0.7))
	{
		startPos = itemNumContigArray - readLen;
		if(startPos<0)
			startPos = 0;
		lowRatioNum = 0;
		for(i=startPos; i<itemNumContigArray; i++)
		{
			if(contigArray[i].occIndexPE[1]!=-1)
			{
				maxIndex = contigArray[i].occIndexPE[0];
				secIndex = contigArray[i].occIndexPE[1];
				maxValue = contigArray[i].occNumPE[maxIndex];
				secValue = contigArray[i].occNumPE[secIndex];

				if((double)secValue/maxValue>=0.25)
					lowRatioNum ++;
			}
		}

		if(lowRatioNum>=1)
		{
			*thisNaviSuccessFlag = NAVI_FAILED;
			return SUCCESSFUL;
		}
	}


	// check the contig path
	if(PEGivenType>NONE_PE_GIVEN_TYPE && (*thisNaviSuccessFlag)==NAVI_SUCCESS && itemNumContigArray>10*minContigLenUsingPE)
	{
		if(contigPath->naviSuccessSize==0 && contigPath->itemNumPathItemList>=2
			&& contigPath->maxPathItem->supportReadsNum<0.5*contigPath->secPathItem->supportReadsNum
			&& (contigPath->maxPathItem->contigPathLen-contigPath->startRowNewBase<10)
			&& (successContigIndex>0 && itemNumContigArr-successContigIndex>20))
		{
			*thisNaviSuccessFlag = NAVI_FAILED;
			return SUCCESSFUL;
		}


		// check the naviSuccessSize, 2014-01-20
		if((contigPath->preNaviSuccessSize>0 && contigPath->naviSuccessSize<0.5*readLen) && contigPath->itemNumPathItemList>=2 && naviTandFlag!=NAVI_SUCCESS)
		{

			// contig 5
			if(navigationFlag==NAVI_PE_FLAG && occsNumIndexPE[1]!=-1)
			{
				maxIndex = occsNumIndexPE[0];
				secIndex = occsNumIndexPE[1];
				maxValue = occsNumPE[maxIndex];
				secValue = occsNumPE[secIndex];
			}
			//else if(navigationFlag!=NAVI_PE_FLAG && occsNumIndexSE[1]!=-1)
			//{
			//	maxIndex = occsNumIndexSE[0];
			//	secIndex = occsNumIndexSE[1];
			//	maxValue = occsNumSE[maxIndex];
			//	secValue = occsNumSE[secIndex];
			//}
			else
			{
				maxValue = secValue = 0;
			}

			if((double)secValue/maxValue>=0.3 && secValue>3)
			{
				startPos = itemNumContigArray - readLen;
				if(startPos<0)
					startPos = 0;
				lowRatioNum = 0;
				lowRatioNum2 = 0;
				if((double)secValue/maxValue>=0.8)
				{
					lowRatioNum ++;
					lowRatioNum2 ++;
				}else if(secValue>3 && (double)secValue/maxValue>0.3)
				{
					lowRatioNum2 ++;
				}

				for(i=startPos; i<itemNumContigArray; i++)
				{
					if(contigArray[i].naviFlag==NAVI_PE_FLAG && contigArray[i].occIndexPE[1]!=-1 && contigArray[i].naviTandFlag!=NAVI_SUCCESS)
					{
						maxIndex = contigArray[i].occIndexPE[0];
						secIndex = contigArray[i].occIndexPE[1];
						maxValue = contigArray[i].occNumPE[maxIndex];
						secValue = contigArray[i].occNumPE[secIndex];

						if(secValue>3 && (double)secValue/maxValue>0.8)
						{
							lowRatioNum ++;
							lowRatioNum2 ++;
						}else if(secValue>3 && (double)secValue/maxValue>0.3)
						{
							lowRatioNum2 ++;
						}
					}
					//else if(contigArray[i].naviFlag!=NAVI_PE_FLAG && contigArray[i].occIndexSE[1]!=-1 && contigArray[i].naviTandFlag!=NAVI_SUCCESS)
					//{
					//	maxIndex = contigArray[i].occIndexSE[0];
					//	secIndex = contigArray[i].occIndexSE[1];
					//	maxValue = contigArray[i].occNumPE[maxIndex];
					//	secValue = contigArray[i].occNumPE[secIndex];

					//	if((double)secValue/maxValue>0.8)
					//		lowRatioNum ++;
					//}
				}

				if(lowRatioNum>=1 && lowRatioNum2>=2)
				{
					*thisNaviSuccessFlag = NAVI_FAILED;
					return SUCCESSFUL;
				}
			}


			// contig 46
			// check the two naviSuccessSize
			if(contigPath->naviSuccessSize==0 && (contigPath->preNaviSuccessSize>0 && contigPath->preNaviSuccessSize<readLen))
			{
				if(navigationFlag==NAVI_PE_FLAG && occsNumIndexPE[1]!=-1)
				{
					maxIndex = occsNumIndexPE[0];
					secIndex = occsNumIndexPE[1];
					maxValue = occsNumPE[maxIndex];
					secValue = occsNumPE[secIndex];
				}
				//else if(navigationFlag!=NAVI_PE_FLAG && occsNumIndexSE[1]!=-1)
				//{
				//	maxIndex = occsNumIndexSE[0];
				//	secIndex = occsNumIndexSE[1];
				//	maxValue = occsNumSE[maxIndex];
				//	secValue = occsNumSE[secIndex];
				//}
				else
				{
					maxValue = secValue = 0;
				}

				startPos = itemNumContigArray - readLen;
				if(startPos<0)
					startPos = 0;
				lowRatioNum = 0;
				if(secValue>3 && (double)secValue/maxValue>=0.7)
					lowRatioNum ++;
				for(i=startPos; i<itemNumContigArray; i++)
				{
					if(contigArray[i].naviFlag==NAVI_PE_FLAG && contigArray[i].occIndexPE[1]!=-1 && contigArray[i].naviTandFlag!=NAVI_SUCCESS)
					{
						maxIndex = contigArray[i].occIndexPE[0];
						secIndex = contigArray[i].occIndexPE[1];
						maxValue = contigArray[i].occNumPE[maxIndex];
						secValue = contigArray[i].occNumPE[secIndex];

						if(secValue>3 && (double)secValue/maxValue>0.7)
							lowRatioNum ++;
					}
					//else if(contigArray[i].naviFlag!=NAVI_PE_FLAG && contigArray[i].occIndexSE[1]!=-1 && contigArray[i].naviTandFlag!=NAVI_SUCCESS)
					//{
					//	maxIndex = contigArray[i].occIndexSE[0];
					//	secIndex = contigArray[i].occIndexSE[1];
					//	maxValue = contigArray[i].occNumPE[maxIndex];
					//	secValue = contigArray[i].occNumPE[secIndex];

					//	if((double)secValue/maxValue>0.7)
					//		lowRatioNum ++;
					//}
				}

				if(lowRatioNum>=1)
				{
					*thisNaviSuccessFlag = NAVI_FAILED;
					return SUCCESSFUL;
				}
			}
		}


		// contig 800
		if(contigPath->naviSuccessSize<0.5*readLen && contigPath->preNaviSuccessSize<0.5*readLen && contigPath->itemNumPathItemList>=2 && (successContigIndex>0 && itemNumContigArr-successContigIndex>0.4*readLen) && (readsNumRatio>1.5))
		{

			startPos = itemNumContigArray - 2*readLen;
			if(startPos<0)
				startPos = 0;
			lowRatioNum = 0;
			for(i=startPos; i<itemNumContigArray; i++)
			{
				if(contigArray[i].naviFlag==NAVI_PE_FLAG && contigArray[i].occIndexPE[1]!=-1 && contigArray[i].naviSuccessSize<readLen)
				{
					maxIndex = contigArray[i].occIndexPE[0];
					secIndex = contigArray[i].occIndexPE[1];
					maxValue = contigArray[i].occNumPE[maxIndex];
					secValue = contigArray[i].occNumPE[secIndex];

					if(secValue>3 && (double)secValue/maxValue>0.3)
						lowRatioNum ++;
				}
				//else if(contigArray[i].naviFlag!=NAVI_PE_FLAG && contigArray[i].occIndexSE[1]!=-1 && contigArray[i].naviSuccessSize<readLen)
				//{
				//	maxIndex = contigArray[i].occIndexSE[0];
				//	secIndex = contigArray[i].occIndexSE[1];
				//	maxValue = contigArray[i].occNumPE[maxIndex];
				//	secValue = contigArray[i].occNumPE[secIndex];

				//	if((double)secValue/maxValue>0.3)
				//		lowRatioNum ++;
				//}
			}

			if(lowRatioNum>=2)
			{
				*thisNaviSuccessFlag = NAVI_FAILED;
				return SUCCESSFUL;
			}
		}

		// contig 714, check the tandPath operation
		if((contigPath->itemNumPathItemList>=2 && contigPath->secPathItem->supportReadsNum>3*averKmerOcc && contigPath->maxPathItem->supportReadsNum>2*contigPath->secPathItem->supportReadsNum && readsNumRatio>2)
			&& (contigPath->naviPathItem && contigPath->naviPathItem==contigPath->maxPathItem))
		{
			if(navigationFlag==NAVI_PE_FLAG && naviTandFlag==NAVI_SUCCESS && occsNumIndexPE[1]!=-1)
			{
				maxIndex = occsNumIndexPE[0];
				secIndex = occsNumIndexPE[1];
				maxValue = occsNumPE[maxIndex];
				secValue = occsNumPE[secIndex];
			}
			//else if(navigationFlag!=NAVI_PE_FLAG && occsNumIndexSE[1]!=-1)
			//{
			//	maxIndex = occsNumIndexSE[0];
			//	secIndex = occsNumIndexSE[1];
			//	maxValue = occsNumSE[maxIndex];
			//	secValue = occsNumSE[secIndex];
			//}
			else
			{
				maxIndex = secIndex = -1;
				maxValue = secValue = 0;
			}

			if(secValue>averKmerOcc && (double)secValue/maxValue>=0.4)
			{
				pathseq = contigPath->maxPathItem->contigPathStr + contigPath->startRowNewBase;
				pathLen = contigPath->maxPathItem->contigPathLen - contigPath->startRowNewBase;
				if(pathLen>0)
				{
					switch(pathseq[0])
					{
						case 'A': baseInt = 0; break;
						case 'C': baseInt = 1; break;
						case 'G': baseInt = 2; break;
						case 'T': baseInt = 3; break;
						default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, pathseq[0]); return FAILED;
					}

					if(baseInt!=maxIndex)
					{
						*thisNaviSuccessFlag = NAVI_FAILED;
						return SUCCESSFUL;
					}
				}
			}
		}


		// contig 799
		if(contigPath->naviSuccessSize<0.5*readLen && contigPath->preNaviSuccessSize<0.5*readLen && (successContigIndex>0 && itemNumContigArr-successContigIndex>0.3*readLen))
		{
			startPos = itemNumContigArray - 4*readLen;
			if(startPos<0)
				startPos = 0;
			lowRatioNum = 0;
			for(i=startPos; i<itemNumContigArray; i++)
			{
				if(contigArray[i].occIndexPE[1]!=-1 && contigArray[i].naviSuccessSize<2*readLen)
				{
					maxIndex = contigArray[i].occIndexPE[0];
					secIndex = contigArray[i].occIndexPE[1];
					maxValue = contigArray[i].occNumPE[maxIndex];
					secValue = contigArray[i].occNumPE[secIndex];

					if(secValue>=3 && (double)secValue/maxValue>=0.2)
						lowRatioNum ++;
				}
			}

			// check the gap region numbers
			startPos = itemNumContigArray - 4*readLen;
			if(startPos<readLen)
				startPos = readLen;
			largeGapNum = 0;

			if(contigArray[startPos].ridposnum>0)
				gapLen = 0;
			else
				gapLen = 1;
			for(i=startPos+1; i<itemNumContigArray; i++)
			{
				if(contigArray[i].ridposnum==0)
				{
					gapLen ++;
				}else
				{
					if(gapLen>=25)
						largeGapNum ++;
					gapLen = 0;
				}
			}

			if(gapLen>=25)
				largeGapNum ++;

			if(lowRatioNum>=5 && largeGapNum>=2)
			{
				*thisNaviSuccessFlag = NAVI_FAILED;
				return SUCCESSFUL;
			}
		}


		// contig 1427
		if(contigPath->naviSuccessSize<0.5*readLen && contigPath->itemNumPathItemList>=3)
		{

			startPos = itemNumContigArray - 3*readLen;
			if(startPos<0)
				startPos = 0;
			lowRatioNum = 0;
			for(i=startPos; i<itemNumContigArray; i++)
			{
				if(contigArray[i].naviFlag==NAVI_PE_FLAG && contigArray[i].occIndexPE[1]!=-1)
				{
					maxIndex = contigArray[i].occIndexPE[0];
					secIndex = contigArray[i].occIndexPE[1];
					maxValue = contigArray[i].occNumPE[maxIndex];
					secValue = contigArray[i].occNumPE[secIndex];

					if(secValue>=3 && (double)secValue/maxValue>0.2)
						lowRatioNum ++;
				}
				//else if(contigArray[i].naviFlag!=NAVI_PE_FLAG && contigArray[i].occIndexSE[1]!=-1 && contigArray[i].naviSuccessSize<readLen)
				//{
				//	maxIndex = contigArray[i].occIndexSE[0];
				//	secIndex = contigArray[i].occIndexSE[1];
				//	maxValue = contigArray[i].occNumPE[maxIndex];
				//	secValue = contigArray[i].occNumPE[secIndex];

				//	if((double)secValue/maxValue>0.3)
				//		lowRatioNum ++;
				//}
			}


			startPos = itemNumContigArray - 0.5*readLen;
			if(startPos<0)
				startPos = 0;
			lowRatioNum2 = 0;
			for(i=startPos; i<itemNumContigArray; i++)
			{
				if(contigArray[i].naviFlag==NAVI_PE_FLAG && contigArray[i].occIndexPE[1]!=-1 && contigArray[i].naviSuccessSize<0.5*readLen)
				{
					maxIndex = contigArray[i].occIndexPE[0];
					secIndex = contigArray[i].occIndexPE[1];
					maxValue = contigArray[i].occNumPE[maxIndex];
					secValue = contigArray[i].occNumPE[secIndex];

					if(secValue>=5 && (double)secValue/maxValue>0.25)
						lowRatioNum2 ++;
				}
				//else if(contigArray[i].naviFlag!=NAVI_PE_FLAG && contigArray[i].occIndexSE[1]!=-1 && contigArray[i].naviSuccessSize<0.5*readLen)
				//{
				//	maxIndex = contigArray[i].occIndexSE[0];
				//	secIndex = contigArray[i].occIndexSE[1];
				//	maxValue = contigArray[i].occNumPE[maxIndex];
				//	secValue = contigArray[i].occNumPE[secIndex];

				//	if((double)secValue/maxValue>0.3)
				//		lowRatioNum ++;
				//}
			}

			if(lowRatioNum>=5 && lowRatioNum2>=2)
			{
				*thisNaviSuccessFlag = NAVI_FAILED;
				return SUCCESSFUL;
			}
		}


		if(contigPath->naviSuccessSize<3*readLen && contigPath->preNaviSuccessSize<readLen && contigPath->itemNumPathItemList>=2)
		{
			if(navigationFlag==NAVI_PE_FLAG && occsNumIndexPE[1]!=-1)
			{
				maxIndex = occsNumIndexPE[0];
				secIndex = occsNumIndexPE[1];
				maxValue = occsNumPE[maxIndex];
				secValue = occsNumPE[secIndex];
			}
			//else if(navigationFlag!=NAVI_PE_FLAG && occsNumIndexSE[1]!=-1)
			//{
			//	maxIndex = occsNumIndexSE[0];
			//	secIndex = occsNumIndexSE[1];
			//	maxValue = occsNumSE[maxIndex];
			//	secValue = occsNumSE[secIndex];
			//}
			else
			{
				maxValue = secValue = 0;
			}

			startPos = itemNumContigArray - readLen;
			if(startPos<0)
				startPos = 0;
			lowRatioNum = 0;
			if(secValue>=3 && (double)secValue/maxValue>0.3)
				lowRatioNum ++;
			for(i=startPos; i<itemNumContigArray; i++)
			{
				if(contigArray[i].naviFlag==NAVI_PE_FLAG && contigArray[i].occIndexPE[1]!=-1 && contigArray[i].naviTandFlag!=NAVI_SUCCESS)
				{
					maxIndex = contigArray[i].occIndexPE[0];
					secIndex = contigArray[i].occIndexPE[1];
					maxValue = contigArray[i].occNumPE[maxIndex];
					secValue = contigArray[i].occNumPE[secIndex];

					if(secValue>=3 && (double)secValue/maxValue>0.3)
						lowRatioNum ++;
				}
				//else if(contigArray[i].naviFlag!=NAVI_PE_FLAG && contigArray[i].occIndexSE[1]!=-1 && contigArray[i].naviTandFlag!=NAVI_SUCCESS)
				//{
				//	maxIndex = contigArray[i].occIndexSE[0];
				//	secIndex = contigArray[i].occIndexSE[1];
				//	maxValue = contigArray[i].occNumPE[maxIndex];
				//	secValue = contigArray[i].occNumPE[secIndex];

				//	if((double)secValue/maxValue>0.7)
				//		lowRatioNum ++;
				//}
			}

			// check the gap region numbers
			startPos = itemNumContigArray - readLen;
			if(startPos<readLen)
				startPos = readLen;
			largeGapNum = 0;

			if(contigArray[startPos].ridposnum>0)
				gapLen = 0;
			else
				gapLen = 1;
			for(i=startPos+1; i<itemNumContigArray; i++)
			{
				if(contigArray[i].ridposnum==0)
				{
					gapLen ++;
				}else
				{
					if(gapLen>30)
						largeGapNum ++;
					gapLen = 0;
				}
			}

			if(gapLen>30)
				largeGapNum ++;

			if(lowRatioNum>=2 && largeGapNum>=1)
			{
				*thisNaviSuccessFlag = NAVI_FAILED;
				return SUCCESSFUL;
			}
		}

		if(contigPath->naviSuccessSize<readLen && contigPath->itemNumPathItemList>=2)
		{
			startPos = itemNumContigArray - readLen;
			if(startPos<0)
				startPos = 0;
			lowRatioNum = 0;
			lowRatioNum2 = 0;
			for(i=startPos; i<itemNumContigArray; i++)
			{
				if(contigArray[i].naviFlag==NAVI_PE_FLAG && contigArray[i].occIndexPE[1]!=-1 && contigArray[i].naviTandFlag!=NAVI_SUCCESS)
				{
					maxIndex = contigArray[i].occIndexPE[0];
					secIndex = contigArray[i].occIndexPE[1];
					maxValue = contigArray[i].occNumPE[maxIndex];
					secValue = contigArray[i].occNumPE[secIndex];

					if(secValue>3 && (double)secValue/maxValue>0.75)
					{
						lowRatioNum ++;
						lowRatioNum2 ++;
					}else if(secValue>3 && (double)secValue/maxValue>0.4)
						lowRatioNum2 ++;
				}
				//else if(contigArray[i].naviFlag!=NAVI_PE_FLAG && contigArray[i].occIndexSE[1]!=-1 && contigArray[i].naviTandFlag!=NAVI_SUCCESS)
				//{
				//	maxIndex = contigArray[i].occIndexSE[0];
				//	secIndex = contigArray[i].occIndexSE[1];
				//	maxValue = contigArray[i].occNumPE[maxIndex];
				//	secValue = contigArray[i].occNumPE[secIndex];

				//	if((double)secValue/maxValue>0.7)
				//		lowRatioNum ++;
				//}
			}

			if(lowRatioNum>=1 && lowRatioNum2>=2)
			{
				*thisNaviSuccessFlag = NAVI_FAILED;
				return SUCCESSFUL;
			}

			// 2014-01-25
			if(contigPath->secPathItem->supportReadsNum>8*averKmerOcc && contigPath->maxPathItem->supportReadsNum>2*contigPath->secPathItem->supportReadsNum && (successContigIndex>0 && itemNumContigArr-successContigIndex>0.3*readLen))
			{
				*thisNaviSuccessFlag = NAVI_FAILED;
				return SUCCESSFUL;
			}
		}
	}


	if((*thisNaviSuccessFlag)==NAVI_SUCCESS && contigArray[itemNumContigArray-1].naviFlag==NAVI_PE_FLAG && contigArray[itemNumContigArray-1].occIndexPE[1]==-1 && contigArray[itemNumContigArray-1].occIndexPE[0]!=-1)
	{
		// fill the ratio array
		if(getBaseNumContigPath(baseNumArray, baseNumIndexArray, contigPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the base ratio by contigPath, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// check the navigation by the ratio of the base_PE
		kmerBaseInt = kmerSeqIntAssembly[entriesPerKmer-1] & 3;

		maxIndex = contigArray[itemNumContigArray-1].occIndexPE[0];
		maxValue = contigArray[itemNumContigArray-1].occNumPE[maxIndex];
		if(maxValue==1 && contigPath->naviPathItem==NULL && contigPath->itemNumPathItemList>2 && contigPath->naviSuccessSize<readLen && baseNumArray[4]>0)
		{
			maxBaseRatio = (double)baseNumArray[baseNumIndexArray[0]] / baseNumArray[4];
			if(kmerBaseInt!=baseNumIndexArray[0] && maxBaseRatio<0.6)
			{
				//printf("=*==*==*==*==N localContigID=%ld, assemblyRound=%d, itemNumContigArr=%ld, kmer_len=%d, occsNumPE:(%d,%d,%d,%d), maxBaseRatio=%.4f\n", localContigID, assemblyRound, itemNumContigArr, kmer_len, contigArray[itemNumContigArray-1].occNumPE[0], contigArray[itemNumContigArray-1].occNumPE[1], contigArray[itemNumContigArray-1].occNumPE[2], contigArray[itemNumContigArray-1].occNumPE[3], maxBaseRatio);
				//outputContigPath(contigPath, YES);


				*thisNaviSuccessFlag = NAVI_FAILED;
			}
		}
	}


	return SUCCESSFUL;
}
