/*
 * PEAssembly.c
 *
 *  Created on: Dec 8, 2011
 *      Author: xiao
 */


#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Build contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short buildEstContigs(char *contigFile)
{
	int32_t i, turnContigIndex, validFlag, tmp_gapSize, validSuccessReadFlag;
	FILE *fpFragmentSize, *fpContigBaseEst, *fpContigHangingEst;;
	char hangingFile[256];

	fpFragmentSize = fopen(fragmentSizeFile, "wb");
	if(fpFragmentSize==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fragmentSizeFile);
		return FAILED;
	}

#if	(RESERVE_EST_FILES==YES)
	// ############################ Debug information ##############################
	strcpy(hangingFile, contigFile);
	strcat(hangingFile, ".hang");
	fpContigBaseEst = fopen(contigFile, "w");
	if(fpContigBaseEst==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigFile);
		return FAILED;
	}
	fpContigHangingEst = fopen(hangingFile, "w");
	if(fpContigHangingEst==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, hangingFile);
		return FAILED;
	}
	// ############################ Debug information ##############################
#endif

	contigNumEstContigArr = 0;
	basesNum = 0;

	localContigID = 0;
	contigsNum = 0;

	validFlag = NO;

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

#if (DEBUG_EST_CONTIG_CHECK==YES)
			// ############################ Debug information ##############################
			if(localContigID==2 && itemNumContigArr>=165 && assemblyRound==FIRST_ROUND_ASSEMBLY)
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
					printf("line=%d, In %s(), localContigID=%ld, cannot update the PE hash table, error!\n", __LINE__, __func__, localContigID);
					return FAILED;
				}

				//if(readsNumInPEHashArr>=minReadsNumPEHashThres && regLenPEHash>=minRegLenUsingPE)
				if(readsNumInPEHashArr>=0 && regLenPEHash>=minRegLenUsingPE)
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
						printf("line=%d, In %s(), localContigID=%ld cannot get next kmer, error!\n", __LINE__, __func__, localContigID);
						return FAILED;
					}
					for(i=0; i<4; i++) { occsNumPE[i] = 0; occsNumIndexPE[i] = -1; }

#if (SVM_NAVI==YES)
					//if((successContigIndex>0 && itemNumContigArr-successContigIndex>50) || readsNumRatio<0.3*minReadsNumRatioThres)
					//if((successContigIndex>0 && itemNumContigArr-successContigIndex>30) || readsNumRatio<0.3*minReadsNumRatioThres) // 2013-11-12
					//if(contigPath->naviSuccessSize<2*readLen && ((successContigIndex>0 && itemNumContigArr-successContigIndex>30) || readsNumRatio<0.3*minReadsNumRatioThres)) // 2013-11-12
					if(contigPath->naviSuccessSize<2*readLen && ((successContigIndex>0 && itemNumContigArr-successContigIndex>0.3*readLen) || readsNumRatio<0.3*minReadsNumRatioThres)) // 2014-12-23
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
				//if(contigPath->naviSuccessSize<2*readLen && ((successContigIndex>0 && itemNumContigArr-successContigIndex>30) || readsNumRatio<0.3*minReadsNumRatioThres)) // 2014-01-31
				if(contigPath->naviSuccessSize<2*readLen && ((successContigIndex>0 && itemNumContigArr-successContigIndex>0.3*readLen) || readsNumRatio<0.3*minReadsNumRatioThres)) // 2014-01-31
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
			//if(kmers[0] || kmers[1])
			if(naviSuccessFlag==NAVI_SUCCESS)
			{
				if(itemNumContigArr>=minContigLenCheckingReadsNum)
				{
					if(updateReadsNumReg(itemNumSuccessReadsArr, itemNumContigArr, assemblyRound)==FAILED)
					{
						printf("line=%d, In %s(), localContigID=%ld, cannnot check the reads number in reads number region, error!\n", __LINE__, __func__, localContigID);
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

#if (DEBUG_EST_CONTIG_CHECK==YES)
			// ############################ Debug information ##############################
			if(localContigID==1 && itemNumContigArr>=21453 && assemblyRound==FIRST_ROUND_ASSEMBLY)
			{
				printf("localContigID=%ld, assemblyRound=%d, itemNumContigArr=%ld, itemNumDecisionTable=%d\n", localContigID, assemblyRound, itemNumContigArr, itemNumDecisionTable);
				printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
				if(successContigIndex>0)
					printf("\tdistance=%ld, readsNumRatio=%.2f\n", itemNumContigArr-successContigIndex, readsNumRatio);
			}
			// ############################ Debug information ##############################
#endif

			//if(kmers[0]==NULL && kmers[1]==NULL)
			if(naviSuccessFlag==NAVI_FAILED)
			{
				if(successContigIndex<=0)
				{
					//printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, the successContig<=0!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
					break;
				}

#if (DEBUG_EST_OUTPUT==YES)
				printf("localContigID=%ld, assemblyRound=%d, itemNumContigArr=%ld, itemNumDecisionTable=%d\n", localContigID, assemblyRound, itemNumContigArr, itemNumDecisionTable);
				printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
				if(successContigIndex>0)
					printf("\tdistance=%ld, readsNumRatio=%.2f\n", itemNumContigArr-successContigIndex, readsNumRatio);
#endif

				// prepare for next round assembly
				if(assemblyRound==FIRST_ROUND_ASSEMBLY)
				{  // finished first round, then start the second round

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
					//	printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, successContigIndex<=0, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
					//	return FAILED;
					//}
					// ############################ Debug information ##############################

					break;
				}
			}

			itemNumContigArr ++;

			// Append a base to contig tail
			if(addContigBase(kmerSeqIntAssembly[entriesPerKmer-1] & 3)==FAILED)
			{
				printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, cannot add a contig base, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
				return FAILED;
			}

			// added 2013-11-27
			// add the read number
			if(addOccsToContig(contigArr, itemNumContigArr, naviTandFlag, newCandBaseNumAfterTandPathPE, contigPath)==FAILED)
			{
				printf("line=%d, In %s(), localContigID=%ld, cannot add occ number to contig, error!\n", __LINE__, __func__, localContigID);
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
				printf("line=%d, In %s(), localContigID=%ld, itemNumContigArr=%ld, cannot update finished reads in decision table, error!\n", __LINE__, __func__, localContigID, itemNumContigArr);
				return FAILED;
			}

			// update the contig path
			if(updateContigPath(contigPath, navigationFlag, kmers, decisionTable, &itemNumDecisionTable, dtRowHashtable, contigArr, itemNumContigArr)==FAILED)
			{
				printf("line=%d, In %s(), localContigID=%ld, cannot update contig path, error!\n", __LINE__, __func__, localContigID);
				return FAILED;
			}

			if(itemNumSuccessReadsArr>0)
			{

				// delete reads from De Bruijn graph
				if(delReadsFromGraph(successReadsArr, itemNumSuccessReadsArr)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, contigID=%d, itemNumContigArr=%ld, cannot delete the reads from graph, error!\n", __LINE__, __func__, localContigID, contigsNum+1, itemNumContigArr);
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

				//successContigIndex = itemNumContigArr;
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
			}else if(itemNumContigArr-successContigIndex > readLen-MIN_OVERLAP_LEN )
			{

				number_of_overlap_less_than_threshold ++;

#if (DEBUG_EST_OUTPUT==YES)
				printf("===localContigID=%ld, assemblyRound=%d, itemNumContigArr=%ld, itemNumDecisionTable=%d\n", localContigID, assemblyRound, itemNumContigArr, itemNumDecisionTable);
				printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
				if(successContigIndex>0)
					printf("\tdistance=%ld, readsNumRatio=%.2f\n", itemNumContigArr-successContigIndex, readsNumRatio);
#endif


				if(assemblyRound==SECOND_ROUND_ASSEMBLY && itemNumContigArr<CONTIG_LEN_THRESHOLD)
				{
					printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld < %d!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr, CONTIG_LEN_THRESHOLD);
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
		}//end while(kmer)

		if(successContigIndex>0)
		{
			if(updateContigtailnodes(contigArr, successContigIndex, &itemNumContigArr, assemblyRound)==FAILED)
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

			if(itemNumContigArr>=minContigLenEst)
			{ // if the contig length is larger than minContigLenEst, then save it to the contig array
				contigsNum ++;

#if	(RESERVE_EST_FILES==YES)
				// ############################ Debug information ##############################
				// output contig nodes to file
				if(outputContigToFile(fpContigBaseEst, BASE_TYPE_FASTA_CONTIG_FILE, contigsNum, contigArr, itemNumContigArr)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, cannot output contig nodes to file, error!\n", __LINE__, __func__, localContigID);
					return FAILED;
				}
				if(outputContigToFile(fpContigHangingEst, HANGING_READ_TYPE_CONTIG_FILE, contigsNum, contigArr, itemNumContigArr)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, cannot output contig nodes to file, error!\n", __LINE__, __func__, localContigID);
					return FAILED;
				}
				// ############################ Debug information ##############################
#endif

				// process single contig
				if(getPairedEndsFromSingleContig(fpFragmentSize, contigArr, itemNumContigArr)==FAILED)
				{
					printf("line=%d, In %s(), cannot get paired reads from contig: %d, error!\n", __LINE__, __func__, contigsNum);
					return FAILED;
				}

				successReadNum += this_successReadNum;
				basesNum += itemNumContigArr;
				contigNumEstContigArr ++;

				if(contigNumEstContigArr>=MAX_NUM_EST_CONTIG || basesNum>=TOTAL_CONTIG_LEN_EST_THRES || successReadNum>=TOTAL_READS_NUM_EST_THRES)
				{
					validFlag = YES;
				}
			}else
			{
				successReadNum -= this_successReadNum;
			}
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

		// clean the contig path
		if(cleanContigPath(contigPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot clean the contig path, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// clean the contig array
		cleanContigArray(contigArr, &itemNumContigArr);

		if(validFlag==YES)
			break;

	} //end while(kmerIndex < TABLE_SIZE_DE_BRUIJN)


	// ############################ Debug information ##############################
	if(contigNumEstContigArr!=contigsNum)
	{
		printf("line=%d, In %s(), contigNumEstContigArr=%d != contigsNum=%d, error!\n", __LINE__, __func__, contigNumEstContigArr, contigsNum);
		return FAILED;
	}
	// ############################ Debug information ##############################

	fclose(fpFragmentSize);
	fpFragmentSize = NULL;

#if	(RESERVE_EST_FILES==YES)
	// ############################ Debug information ##############################
	fclose(fpContigBaseEst);
	fpContigBaseEst = NULL;
	fclose(fpContigHangingEst);
	fpContigHangingEst = NULL;
	// ############################ Debug information ##############################
#endif

	if(contigNumEstContigArr>0)
		return SUCCESSFUL;
	else
		return FAILED;
}

/**
 * Get the next kmer by mixture of PE and SE.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNextKmerByMix(int32_t contigNodesNum, int32_t assemblyRound)
{
	kmertype *tmp_kmers[2];
	int32_t tmp_gapSize, naviPE, naviSE;
	int32_t useCandPathPE, useCandPathSE, returnFlag;
	int32_t naviCandPathPE, naviCandPathSE, naviTandPathPE, naviTandPathSE, naviContigPathPE, naviContigPathSE;

	useCandPathPE = useCandPathSE = NAVI_UNUSED;
	naviTandFlag = NAVI_UNUSED;
	newCandBaseNumAfterTandPathPE = -1;

	tmp_gapSize = -1;
	naviPE = naviSE = NAVI_UNUSED;
	naviCandPathPE = naviCandPathSE = NAVI_UNUSED;
	navigationFlag = NAVI_PE_FLAG;
	tmp_kmers[0] = kmers[0];
	tmp_kmers[1] = kmers[1];
	if(memcpy(tmpKmerSeqIntAssembly, kmerSeqIntAssembly, deBruijnGraph->bytesPerKmerseq)==NULL)
	{
		printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the score and occsNum of PE
	if(getNextKmerByPE(contigNodesNum)==FAILED)
	{
		printf("line=%d, In %s(), cannot get next kmer, error!\n", __LINE__, __func__);
		return FAILED;
	}

#if (SVM_NAVI==YES)
	useCandPathPE = NO;
	if(readsNumRatio<0.3*minReadsNumRatioThres && contigPath->itemNumPathItemList>=2)		// added 2013-04-01
	{
		naviSuccessFlag = NAVI_FAILED;
		return SUCCESSFUL;
	}
	//else if((successContigIndex>0 && contigNodesNum-successContigIndex>50) && readsNumRatio<0.4 && maxOccPE<=5)	// added 2013-11-08
	else if(contigPath->itemNumPathItemList>=2 && (successContigIndex>0 && contigNodesNum-successContigIndex>50) && readsNumRatio<0.4 && maxOccPE<=5)	// added 2014-02-01
	{
		naviSuccessFlag = NAVI_FAILED;
		//return SUCCESSFUL;
	}
//	else if(secondOccPE==0 && (successContigIndex>0 && contigNodesNum-successContigIndex>30) && readsNumRatio<0.6 && maxOccPE<2)	// added 2013-11-12
//	{
//		naviSuccessFlag = NAVI_FAILED;
//		//return SUCCESSFUL;
//	}
	else if(secondOccPE==0 && readsNumRatio<0.5 && maxOccPE<3 && maxOccPE<0.1*averKmerOcc) // 2014-01-15
	{
		if(contigPath->secPathItem && contigPath->secPathItem->supportReadsNum>0 && contigPath->maxPathItem->supportReadsNum>20*contigPath->secPathItem->supportReadsNum
			&& contigPath->maxPathItem->contigPathStr[contigPath->startRowNewBase]!=contigPath->secPathItem->contigPathStr[contigPath->startRowNewBase])
			naviSuccessFlag = NAVI_FAILED;
	}
	else if(secondOccPE>0)
	{
//		sumSecondOccPE = 0;
//		for(i=0; i<4; i++) if(i!=occsNumIndexPE[0]) sumSecondOccPE += occsNumPE[i];

		svmFeatureArr[0] = occsNumPE[occsNumIndexPE[0]];
		svmFeatureArr[1] = occsNumPE[occsNumIndexPE[1]];
		//svmFeatureArr[1] = sumSecondOccPE;
		svmFeatureArr[2] = readsNumRatio;
		//svmFeatureArr[2] = svmFeatureArr[0] / svmFeatureArr[1];
		if(successContigIndex>0)
		{
			// compute the maximal gap size in contig tail region
			if(computeGapSizeInContig(&tmp_gapSize, contigArr, contigNodesNum, assemblyRound)==FAILED)
			{
				printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
				return FAILED;
			}

			//svmFeatureArr[3] = contigNodesNum - successContigIndex;
			svmFeatureArr[3] = tmp_gapSize;
		}else
			svmFeatureArr[3] = 0;

		if(fillSampleDataSVM(svmSamplePE, svmFeatureArr)==FAILED)
		{
			printf("line=%d, In %s(), cannot fill the sample data for SVM model, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(extensionDecisionBySvm(&naviPE, svmSamplePE, svmModelPE)==FAILED)
		{
			printf("line=%d, In %s(), cannot decide the extension by SVM calssifier, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(naviPE==NAVI_FAILED && svmFeatureArr[0]/svmFeatureArr[1]>=3)
			naviPE = NAVI_SUCCESS;

		if(secondOccPE/maxOccPE>=0.25 && maxOccPE<0.3*averKmerOcc)
		{
			naviPE = NAVI_FAILED;
		}

		//if(secondOccPE>0.5*averKmerOcc && secondOccPE/maxOccPE>=0.4)  // added 2013-11-12
		if(secondOccPE/maxOccPE>=0.7)   // added 2013-11-14
		{
			useCandPathPE = YES;
		}

		// malign conditions
		//if(naviPE==NAVI_FAILED || (secondOccPE>10 && readsNumRatio>5))
		//if(naviPE==NAVI_FAILED || (secondOccPE>3 && readsNumRatio>1.5))  // 2013-10-18, reused on 2013-11-16
		//if(naviPE==NAVI_FAILED || (secondOccPE>3 && readsNumRatio>1.5) || (secondOccPE/maxOccPE>=0.5 && maxOccPE-secondOccPE<=3))  // 2013-10-18, deleted 2013-11-16
		//if(naviPE==NAVI_FAILED || (secondOccPE>3 && readsNumRatio>1.2) || (secondOccPE/maxOccPE>=0.5 && maxOccPE-secondOccPE<=3))  // 2013-12-13
		if(naviPE==NAVI_FAILED || (secondOccPE>3) || (secondOccPE/maxOccPE>=0.5 && maxOccPE-secondOccPE<=3))  // 2013-12-26
		{
			useCandPathPE = YES;
		}

		// use look ahead approach to check the extension
		if(useCandPathPE==YES)
		{
			// determine the extension by candidate paths, added 2012-11-12
			if(decideByCandPathPE(&naviCandPathPE, &maxBaseIndexAfterCandPathPE, &incorrectBaseNumCandPathPE, occsNumPE, occsNumIndexPE, decisionTable, itemNumDecisionTable, contigPath)==FAILED)
			{
				printf("line=%d, In %s(), cannot navigate by candidate paths using paired ends, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(naviCandPathPE==NAVI_FAILED && incorrectBaseNumCandPathPE>0)  // added 2013-12-02
			{
				// check the tandem repeats
				if(decideByCheckingTandemRepeatPE(&naviTandPathPE, &maxBaseIndexAfterTandPathPE, &incorrectBaseNumTandPathPE, &newCandBaseNumAfterTandPathPE, occsNumPE, occsNumIndexPE, itemNumContigArr-1, decisionTable, &itemNumDecisionTable, dtRowHashtable, contigPath, deBruijnGraph)==FAILED)
				{
					printf("line=%d, In %s(), cannot check the tandem repeats, error!\n", __LINE__, __func__);
					return FAILED;
				}

				naviTandFlag = naviTandPathPE;

				naviPE = naviTandPathPE;
			}else
			{
				naviPE = naviCandPathPE;
			}

			if(naviPE==NAVI_FAILED)  // 2014-01-04
			{
				// navigate using contigPath
				if(decideByContigPath(&naviContigPathPE, contigPath, occsNumPE, occsNumIndexPE, deBruijnGraph, OCC_RATIO_THRES_NAVI_PATH_PE)==FAILED)
				{
					printf("line=%d, In %s(), cannot navigate by the contig path, error!\n", __LINE__, __func__);
					return FAILED;
				}

				naviTandFlag = naviContigPathPE;
				naviPE = naviContigPathPE;
			}

		}

		naviSuccessFlag = naviPE;

	}

	// 2014-03-14
	if(occsNumIndexPE[0]!=-1)
	{
		if(confirmNaviPE(&naviSuccessFlag, occsNumPE, occsNumIndexPE, contigPath, deBruijnGraph)==FAILED)
		{
			printf("line=%d, In %s(), cannot confirm paired-end navigation, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

#else
	//if(maxOccPE==secondOccPE)
	if((maxOccPE==secondOccPE) || (secondOccPE>=2 && secondOccPE>=SEC_MAX_OCC_RATIO_SVM*maxOccPE))
	{
		naviSuccessFlag = NAVI_FAILED;
		//kmers[0] = kmers[1] = NULL;
	}
#endif


	if(naviSuccessFlag==NAVI_FAILED)
	{
		// get the score and occsNum of SE
		navigationFlag = NAVI_MIX_FLAG;
		kmers[0] = tmp_kmers[0];
		kmers[1] = tmp_kmers[1];
		if(memcpy(kmerSeqIntAssembly, tmpKmerSeqIntAssembly, entriesPerKmer*sizeof(uint64_t))==NULL)
		{
			printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(getNextKmerBySE(contigNodesNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot get next kmer, error!\n", __LINE__, __func__);
			return FAILED;
		}


#if (SVM_NAVI==YES)
		useCandPathSE = NO;
		//if((successContigIndex>0 && itemNumContigArr-successContigIndex>60))
		if((successContigIndex>0 && itemNumContigArr-successContigIndex>50))   // added 2013-10-14
			naviSuccessFlag = NAVI_FAILED;
		else if(maxOccPE==0 && (successContigIndex>0 && itemNumContigArr-successContigIndex>40) && readsNumRatio<0.5) // added 2013-11-08
		//else if(contigPath->itemNumPathItemList>=2 && maxOccPE==0 && (successContigIndex>0 && itemNumContigArr-successContigIndex>40) && readsNumRatio<0.5) // added 2014-03-12
		{
			naviSuccessFlag = NAVI_FAILED;
		}
		else if(naviCandPathPE==NAVI_FAILED && (successContigIndex>0 && itemNumContigArr-successContigIndex>25))  // 2014-01-04
		{
			naviSuccessFlag = NAVI_FAILED;
		}
//		else if(successContigIndex>0 && (readsNumRatio<0.3 && maxOccPE==0)) // added 2013-04-23
//		{
//			// compute the maximal gap size in contig tail region
//			if(computeGapSizeInContig(&tmp_gapSize,contigArr, contigNodesNum, assemblyRound)==FAILED)
//			{
//				printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
//				return FAILED;
//			}
//
//			if(tmp_gapSize>15)
//			{
//				naviSuccessFlag = NAVI_FAILED;
//			}
//		}
		else if(secondOccSE==0) // added 2013-11-08
		{
			if(secondOccPE>0 && (secondOccPE/maxOccPE>=0.5 && maxOccPE-secondOccPE<=3) && maxOccIndexPE!=maxOccIndexSE)
			{
				//naviSuccessFlag = NAVI_FAILED;

				// compute the maximal gap size in contig tail region
				if(computeGapSizeInContig(&tmp_gapSize,contigArr, contigNodesNum, assemblyRound)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(decideByCandPathSE(&naviCandPathSE, &maxBaseIndexAfterCandPathSE, &incorrectBaseNumCandPathSE, occsNumSE, occsNumIndexSE, occsNumPE, occsNumIndexPE, decisionTable, itemNumDecisionTable, contigPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot navigate by candidate paths using single ends, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(naviCandPathSE==NAVI_FAILED && incorrectBaseNumCandPathSE>0) // added 2013-12-02
				{
					// check the tandem repeats
					if(decideByCheckingTandemRepeatSE(&naviTandPathSE, &maxBaseIndexAfterTandPathSE, &incorrectBaseNumTandPathSE, &newCandBaseNumAfterTandPathSE, occsNumSE, occsNumIndexSE, itemNumContigArr-1, decisionTable, &itemNumDecisionTable, contigPath, deBruijnGraph)==FAILED)
					{
						printf("line=%d, In %s(), cannot check the tandem repeats, error!\n", __LINE__, __func__);
						return FAILED;
					}

					naviTandFlag = naviTandPathSE;

					naviSE = naviTandPathSE;
				}else
				{
					naviSE = naviCandPathSE;
				}

//				if(naviSE==NAVI_FAILED)  // 2014-01-04
//				{
//					// navigate using contigPath
//					if(decideByContigPath(&naviContigPathSE, contigPath, occsNumSE, occsNumIndexSE, deBruijnGraph)==FAILED)
//					{
//						printf("line=%d, In %s(), cannot navigate by the contig path, error!\n", __LINE__, __func__);
//						return FAILED;
//					}
//
//					naviSE = naviContigPathSE;
//				}


				naviSuccessFlag = naviSE;
			}
		}
		else if(secondOccSE>0)								// added 2013-02-26
		{
//			sumSecondOccSE = 0;
//			for(i=0; i<4; i++) if(i!=occsNumIndexSE[0]) sumSecondOccSE += occsNumSE[i];

			svmFeatureArr[0] = occsNumSE[occsNumIndexSE[0]];
			svmFeatureArr[1] = occsNumSE[occsNumIndexSE[1]];
			//svmFeatureArr[1] = sumSecondOccSE;
			svmFeatureArr[2] = readsNumRatio;
			//svmFeatureArr[2] = svmFeatureArr[0] / svmFeatureArr[1];
			tmp_gapSize = -1;
			if(successContigIndex>0)
			{
				// compute the maximal gap size in contig tail region
				if(computeGapSizeInContig(&tmp_gapSize,contigArr, contigNodesNum, assemblyRound)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
					return FAILED;
				}

				//svmFeatureArr[3] = contigNodesNum - successContigIndex;
				svmFeatureArr[3] = tmp_gapSize;
			}else
				svmFeatureArr[3] = 0;

			if(fillSampleDataSVM(svmSampleSE, svmFeatureArr)==FAILED)
			{
				printf("line=%d, In %s(), cannot fill the sample data for SVM model, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(extensionDecisionBySvm(&naviSE, svmSampleSE, svmModelSE)==FAILED)
			{
				printf("line=%d, In %s(), cannot decide the extension by SVM calssifier, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(secondOccSE>0 && secondOccSE<3 && maxOccSE>40*secondOccSE) // 2014-01-29
			{
				naviSE = NAVI_SUCCESS;
			}
			else if(maxOccIndexPE!=-1 && maxOccIndexPE!=maxOccIndexSE)
			{ // failed
				naviSE = NAVI_FAILED;
			}
			//else if(tmp_gapSize>50 && svmFeatureArr[1]/svmFeatureArr[0]>=0.5)
			//else if(tmp_gapSize>=30)
			else if(tmp_gapSize>=30 && secondOccSE/maxOccSE>=0.3)
			{ // failed
				naviSE = NAVI_FAILED;
			}
			else if(maxOccPE==0 && secondOccSE/maxOccSE>=0.3 && maxOccSE<0.3*averKmerOcc)
			{
				naviSE = NAVI_FAILED;
			}
			//else if(maxOccSE-secondOccSE<2 && readsNumRatio>2)
			else if(maxOccSE-secondOccSE<=3 && readsNumRatio>1.5) // 2014-02-05
			{
				naviSE = NAVI_FAILED;
			}
			//else if((secondOccPE>10 && readsNumRatio>6 && tmp_gapSize>10) || (secondOccPE/maxOccPE>=0.7))
			else if((naviCandPathPE==NAVI_FAILED) || (secondOccPE>0 && secondOccPE>=maxOccPE*0.7))
			{
				//if(kmer_len<longKmerSize-longKmerStepSize || (secondOccSE/maxOccSE>0.1 || secondOccSE>3))
				//if(maxOccPE<10*secondOccPE || maxOccSE<50*secondOccSE)   // added 2014-01-15
					naviSE = NAVI_FAILED;
			}
			else if(maxOccPE==0 && secondOccSE>3*averKmerOcc && secondOccSE/maxOccSE>=0.3 && readsNumRatio>5) // 2013-12-28
			{
				naviSE = NAVI_FAILED;
			}

//			else if(maxOccPE==0 && secondOccSE>3 && secondOccSE/maxOccSE>=0.3)  // 2014-03-25
//			{
//				printf("^^^^^^^^^^^^^^^ contigNodesNum=%d, maxOccSE=%d, secondOccSE=%d\n", contigNodesNum, (int32_t)maxOccSE, (int32_t)secondOccSE);
//				useCandPathSE = YES;
//			}

			//else if(secondOccSE>averKmerOcc && secondOccSE/maxOccSE>=0.4 && readsNumRatio<0.5)  // 2014-02-14
			//{
			//	naviSE = NAVI_FAILED;
			//}


			//if(naviSE==NAVI_FAILED && tmp_gapSize<30)	// best
			//if(naviSE==NAVI_FAILED)
			if((naviSE==NAVI_FAILED || useCandPathSE==YES) && tmp_gapSize<30)
			{
				if(decideByCandPathSE(&naviCandPathSE, &maxBaseIndexAfterCandPathSE, &incorrectBaseNumCandPathSE, occsNumSE, occsNumIndexSE, occsNumPE, occsNumIndexPE, decisionTable, itemNumDecisionTable, contigPath)==FAILED)
				{
					printf("line=%d, In %s(), cannot navigate by candidate paths using single ends, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(naviCandPathSE==NAVI_FAILED && incorrectBaseNumCandPathSE>0)  // added 2013-12-02
				{
					// check the tandem repeats
					if(decideByCheckingTandemRepeatSE(&naviTandPathSE, &maxBaseIndexAfterTandPathSE, &incorrectBaseNumTandPathSE, &newCandBaseNumAfterTandPathSE, occsNumSE, occsNumIndexSE, itemNumContigArr-1, decisionTable, &itemNumDecisionTable, contigPath, deBruijnGraph)==FAILED)
					{
						printf("line=%d, In %s(), cannot check the tandem repeats, error!\n", __LINE__, __func__);
						return FAILED;
					}

					naviTandFlag = naviTandPathSE;

					naviSE = naviTandPathSE;
				}else
				{
					naviSE = naviCandPathSE;
				}

				if(naviSE==NAVI_FAILED)  // 2014-01-16
				{
					// navigate using contigPath
					if(decideByContigPath(&naviContigPathSE, contigPath, occsNumSE, occsNumIndexSE, deBruijnGraph, OCC_RATIO_THRES_NAVI_PATH_SE)==FAILED)
					{
						printf("line=%d, In %s(), cannot navigate by the contig path, error!\n", __LINE__, __func__);
						return FAILED;
					}

					naviTandFlag = naviContigPathSE;
					naviSE = naviContigPathSE;
				}

			}

			naviSuccessFlag = naviSE;
		}
#else
		//if(maxOccSE==secondOccSE)
		if((maxOccSE==secondOccSE) || (secondOccSE>=2 && secondOccSE>=SEC_MAX_OCC_RATIO_SVM*maxOccSE))
		{
			naviSuccessFlag = NAVI_FAILED;
			//kmers[0] = kmers[1] = NULL;
		}
//		else if(maxOccIndexPE!=-1 && maxOccIndexPE!=maxOccIndexSE)
//		{
//			naviSuccessFlag = NAVI_FAILED;
//		}
#endif
	}

	return SUCCESSFUL;
}

/**
 * Get the next k-mers for navigation combining the k-mer hash table and the decision table using paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNextKmerByPE(int32_t contigNodesNum)
{
	if(shortInsertFlag==NO)
	{
		if(getNextKmerByPEEqualOverlap(contigNodesNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot get next kmer, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		if(getNextKmerByPEVariableOverlap(contigNodesNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot get next kmer, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}


/**
 * Get the next k-mers for navigation combining the k-mer hash table and the decision table using paired ends with equal overlap size.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNextKmerByPEEqualOverlap(int32_t contigNodesNum)
{
	uint64_t tmp_kmerseq[entriesPerKmer], tmp_kmerseqRev[entriesPerKmer];
	kmertype *tmp_kmers[4][2];
	int32_t i, j, validKmerNum, iterateID, successiveAppearBaseNum, ignoreEndBaseFlag, base_index;
	int32_t tmpOccIndexArray[2], itemNumTmpOccIndexArray;

	kmer_len = 0;
	maxOccIndexPE = -1;
	maxOccPE = 0;
	secondOccIndexPE = -1;
	secondOccPE = 0;

	for(i=0; i<4; i++) {occsNumPE[i] = 0; occsNumIndexPE[i] = -1;}
	for(i=0; i<4; i++)
	{
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				tmp_kmerseq[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
			}
			tmp_kmerseq[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
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
			ignoreEndBaseFlag = NO;
		}else
		{
			successiveAppearBaseNum = minSuccessiveAppearedBaseNum;
			ignoreEndBaseFlag = NO;
		}

//		if(contigNodesNum>=longKmerSize)
//			kmer_len = longKmerSize;
//		else
//			kmer_len = contigNodesNum;
//		while(kmer_len>=kmerSize)
//		{
			validKmerNum = 0;
			maxOccPE = 0, secondOccPE = 0;
			maxOccIndexPE = -1, secondOccIndexPE = -1;
//			occsNumIndexPE[0] = -1;
//			occsNumIndexPE[1] = -1;
			for(i=0; i<4; i++)
			{
				occsNumPE[i] = 0;
				occsNumIndexPE[i] = -1;

				//if(contigNodesNum>READ_LEN)
//				if(contigNodesNum>kmer_len)
//				//if(contigNodesNum>=kmer_len)  //--bad result
//				{
//					if(computeLongKmerOccNumByPE(occsNumPE+i, tmp_kmers[i],  i, kmer_len, successiveAppearBaseNum, ignoreEndBaseFlag)==FAILED)
//					{
//						printf("line=%d, In %s(), cannot compute score and occsNum of long k-mers, error!\n", __LINE__, __func__);
//						return FAILED;
//					}
//				}else
//				{

					//if(computeKmerScoreByPE(tmp_kmers[i], score+i, occsNum+i, contigNodesNum, assemblyRound)==FAILED)		// 100_1
					if(computeKmerOccNumUnlockedByPE(occsNumPE+i, tmp_kmers[i], i, successiveAppearBaseNum, ignoreEndBaseFlag)==FAILED)	//100_2
					//if(computeLongKmerOccNumByPE(occsNumPE+i, tmp_kmers[i],  i, kmer_len)==FAILED)
					{
						printf("line=%d, In %s(), cannot compute occsNum of k-mers, error!\n", __LINE__, __func__);
						return FAILED;
					}
//				}

				//score[i] = computeKmerScore(tmp_kmers[i], occsNum+i, assemblingreads, numassemblingreads);

				//========================= condition 3 ==============================
				//if(score[i]>=MIN_SCORE_THRESHOLD/**MIN_SCORE_FACTOR*/)
				//if(score[i]>=MIN_SCORE_THRESHOLD && occsNum[i]>=MIN_CONNECT_KMER_NUM)
				//if(occsNum[i]>=MIN_CONNECT_KMER_NUM)
				if(occsNumPE[i]>0)
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
				maxOccPE = 0, maxOccIndexPE = -1, secondOccPE = 0, secondOccIndexPE = -1;
				for(j=0; j<4; j++)
				{
					if(maxOccPE<occsNumPE[j])
					{
						secondOccPE = maxOccPE;
						secondOccIndexPE = maxOccIndexPE;
						maxOccPE = occsNumPE[j];
						maxOccIndexPE = j;
					}else if(secondOccPE<occsNumPE[j])
					{
						secondOccPE = occsNumPE[j];
						secondOccIndexPE = j;
					}
				}

				occsNumIndexPE[0] = maxOccIndexPE;
				occsNumIndexPE[1] = secondOccIndexPE;

				if(secondOccPE>0)
				{
					itemNumTmpOccIndexArray = 0;
					for(j=0; j<4; j++)
					{
						if(j!=maxOccIndexPE && j!=secondOccIndexPE && occsNumPE[j]>0)
							tmpOccIndexArray[itemNumTmpOccIndexArray++] = j;
					}

					if(itemNumTmpOccIndexArray==1)
					{
						occsNumIndexPE[2] = occsNumPE[ tmpOccIndexArray[itemNumTmpOccIndexArray-1] ];
					}else if(itemNumTmpOccIndexArray==2)
					{
						if(tmpOccIndexArray[0]>=tmpOccIndexArray[1])
						{
							occsNumIndexPE[2] = occsNumPE[ tmpOccIndexArray[0] ];
							occsNumIndexPE[3] = occsNumPE[ tmpOccIndexArray[1] ];
						}else
						{
							occsNumIndexPE[2] = occsNumPE[ tmpOccIndexArray[1] ];
							occsNumIndexPE[3] = occsNumPE[ tmpOccIndexArray[0] ];
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
				//if(occsNum[tmp_maxIndex]<MIN_CONNECT_KMER_NUM)
		//		if(occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM)
		//		{
		//			kmers[0] = kmers[1] = NULL;
		//		}else
		//		{

					if(entriesPerKmer>=2)
					{
						for(j=0; j<entriesPerKmer-2; j++)
						{
							kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
						}
						kmerSeqIntAssembly[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
					}
					kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | maxOccIndexPE) & lastEntryMask;

					kmers[0] = tmp_kmers[maxOccIndexPE][0];
					kmers[1] = tmp_kmers[maxOccIndexPE][1];

					naviSuccessFlag = NAVI_SUCCESS;
		//		}

				return SUCCESSFUL;

			}
			//***********************************************************************************
			else if(validKmerNum==0)
			{
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
				kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | maxOccIndexPE) & lastEntryMask;

				kmers[0] = tmp_kmers[maxOccIndexPE][0];
				kmers[1] = tmp_kmers[maxOccIndexPE][1];

				naviSuccessFlag = NAVI_SUCCESS;

				return SUCCESSFUL;
			}
//		}
	}

	naviSuccessFlag = NAVI_FAILED;
	kmers[0] = kmers[1] = NULL;

	return SUCCESSFUL;
}

/**
 * Get the next k-mers for navigation combining the k-mer hash table and the decision table using paired ends with variable overlap size.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNextKmerByPEVariableOverlap(int32_t contigNodesNum)
{
	uint64_t tmp_kmerseq[entriesPerKmer], tmp_kmerseqRev[entriesPerKmer];
	kmertype *tmp_kmers[4][2];
	int32_t i, j, validKmerNum, iterateID, successiveAppearBaseNum, ignoreEndBaseFlag, base_index;
	int32_t tmpOccIndexArray[2], itemNumTmpOccIndexArray;

	kmer_len = 0;
	maxOccIndexPE = -1;
	maxOccPE = 0;
	secondOccIndexPE = -1;
	secondOccPE = 0;

	for(i=0; i<4; i++) {occsNumPE[i] = 0; occsNumIndexPE[i] = -1;}
	for(i=0; i<4; i++)
	{
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				tmp_kmerseq[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
			}
			tmp_kmerseq[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
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
		while(kmer_len>=kmerSize)
		{
			validKmerNum = 0;
			maxOccPE = 0, secondOccPE = 0;
			maxOccIndexPE = -1, secondOccIndexPE = -1;
			occsNumIndexPE[0] = -1;
			occsNumIndexPE[1] = -1;
			for(i=0; i<4; i++)
			{
				occsNumPE[i] = 0;
				occsNumIndexPE[i] = -1;

				//if(contigNodesNum>READ_LEN)
//				if(contigNodesNum>kmer_len)
//				//if(contigNodesNum>=kmer_len)  //--bad result
//				{
					if(computeLongKmerOccNumByPE(occsNumPE+i, tmp_kmers[i],  i, kmer_len, successiveAppearBaseNum, ignoreEndBaseFlag)==FAILED)
					{
						printf("line=%d, In %s(), cannot compute score and occsNum of long k-mers, error!\n", __LINE__, __func__);
						return FAILED;
					}
//				}else
//				{

					//if(computeKmerOccNumUnlockedByPE(occsNumPE+i, tmp_kmers[i], i, successiveAppearBaseNum, ignoreEndBaseFlag)==FAILED)	//100_2
					//{
					//	printf("line=%d, In %s(), cannot compute occsNum of k-mers, error!\n", __LINE__, __func__);
					//	return FAILED;
					//}
//				}


				//========================= condition 3 ==============================
				//if(score[i]>=MIN_SCORE_THRESHOLD/**MIN_SCORE_FACTOR*/)
				//if(score[i]>=MIN_SCORE_THRESHOLD && occsNum[i]>=MIN_CONNECT_KMER_NUM)
				//if(occsNum[i]>=MIN_CONNECT_KMER_NUM)
				if(occsNumPE[i]>0)
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
				maxOccPE = 0, maxOccIndexPE = -1, secondOccPE = 0, secondOccIndexPE = -1;
				for(j=0; j<4; j++)
				{
					if(maxOccPE<occsNumPE[j])
					{
						secondOccPE = maxOccPE;
						secondOccIndexPE = maxOccIndexPE;
						maxOccPE = occsNumPE[j];
						maxOccIndexPE = j;
					}else if(secondOccPE<occsNumPE[j])
					{
						secondOccPE = occsNumPE[j];
						secondOccIndexPE = j;
					}
				}

				occsNumIndexPE[0] = maxOccIndexPE;
				occsNumIndexPE[1] = secondOccIndexPE;

				if(secondOccPE>0)
				{
					itemNumTmpOccIndexArray = 0;
					for(j=0; j<4; j++)
					{
						if(j!=maxOccIndexPE && j!=secondOccIndexPE && occsNumPE[j]>0)
							tmpOccIndexArray[itemNumTmpOccIndexArray++] = j;
					}

					if(itemNumTmpOccIndexArray==1)
					{
						occsNumIndexPE[2] = occsNumPE[ tmpOccIndexArray[itemNumTmpOccIndexArray-1] ];
					}else if(itemNumTmpOccIndexArray==2)
					{
						if(tmpOccIndexArray[0]>=tmpOccIndexArray[1])
						{
							occsNumIndexPE[2] = occsNumPE[ tmpOccIndexArray[0] ];
							occsNumIndexPE[3] = occsNumPE[ tmpOccIndexArray[1] ];
						}else
						{
							occsNumIndexPE[2] = occsNumPE[ tmpOccIndexArray[1] ];
							occsNumIndexPE[3] = occsNumPE[ tmpOccIndexArray[0] ];
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
				if(maxOccPE==1 && kmer_len > kmerSize)
				{
					kmer_len -= longKmerStepSize;
					if(kmer_len<kmerSize && kmer_len>kmerSize-longKmerStepSize)
						kmer_len = kmerSize;

					continue;
				}

				//if(occsNum[tmp_maxIndex]<MIN_CONNECT_KMER_NUM)
		//		if(occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM)
		//		{
		//			kmers[0] = kmers[1] = NULL;
		//		}else
		//		{

					if(entriesPerKmer>=2)
					{
						for(j=0; j<entriesPerKmer-2; j++)
						{
							kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
						}
						kmerSeqIntAssembly[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
					}
					kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | maxOccIndexPE) & lastEntryMask;

					kmers[0] = tmp_kmers[maxOccIndexPE][0];
					kmers[1] = tmp_kmers[maxOccIndexPE][1];

					naviSuccessFlag = NAVI_SUCCESS;
		//		}

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
				if(maxOccPE<2 && kmer_len > kmerSize)
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
					kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | maxOccIndexPE) & lastEntryMask;

					kmers[0] = tmp_kmers[maxOccIndexPE][0];
					kmers[1] = tmp_kmers[maxOccIndexPE][1];

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
 * Compute the kmer occNum.
 *  	If there some locked reads, then the occNum will be computed by considering the locked reads;
 *  	otherwise, consider all the reads in decision table.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
*/
short computeKmerOccNumByPE(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag)
{
	if(lockedReadsNum>=lockedReadsNumThres)
	{
		if(computeKmerOccNumLockedByPE(occNum, tmp_kmers, baseInt_kmer, successiveAppearBaseNum, ignoreEndBaseFlag)==FAILED)
		{
			printf("In %s(), cannot compute kmer occNum in local assembly, error!\n", __func__);
			return FAILED;
		}
	}else
	{
		if(computeKmerOccNumUnlockedByPE(occNum, tmp_kmers, baseInt_kmer, successiveAppearBaseNum, ignoreEndBaseFlag)==FAILED)
		{
			printf("In %s(), cannot compute kmer occNum in local assembly, error!\n", __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Compute the occurrence number for unlocked paired reads in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeKmerOccNumUnlockedByPE(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag)
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
	int32_t seqLen, errorRegLenEnd3, validFlag, mismatchNum;
	assemblingreadtype *dtReadPaired;


	*occNum = 0;
	for(i=0; i<itemNumDecisionTable; i++)
	{
		//if(decisionTable[i].matedFlag==YES && decisionTable[i].reserved==0)
		//if(decisionTable[i].matedFlag==YES && decisionTable[i].reserved==0 && decisionTable[i].unmatchBaseNum<=3)
		if(decisionTable[i].matedFlag==YES && decisionTable[i].reserved==0
			&& decisionTable[i].unmatchBaseNum<=maxUnmatchBaseNumPerRead
			&& (decisionTable[i].successiveAppearBases>=successiveAppearBaseNum)	// added 2013-01-12
			&& decisionTable[i].basePos==decisionTable[i].lastMatchedBasePos)
		{
			properRow = getProperDtRow(decisionTable+i, decisionTable, dtRowHashtable);

			// ############################ Debug information ##############################
			if(properRow<0 || properRow>=itemNumDecisionTable || decisionTable[properRow].rid!=decisionTable[i].rid)
			{
				printf("line=%d, In %s(), properRow=%d, i=%d, error!\n", __LINE__, __func__, properRow, i);
				//continue;
				return FAILED;
			}
			// ############################ Debug information ##############################

			this_assemblingRead = decisionTable + properRow;

			// ######################### Debug information ###########################
			//if(this_assemblingRead->matchBaseNum+this_assemblingRead->unmatchBaseNum>=60
			//	&& (double)this_assemblingRead->unmatchBaseNum/this_assemblingRead->matchBaseNum>=0.5)
			//{
			//	printf("line=%d, In %s(), rid=%lu, matchBaseNum=%d, unmatchBaseNum=%d\n", __LINE__, __func__, (uint64_t)this_assemblingRead->rid, this_assemblingRead->matchBaseNum, this_assemblingRead->unmatchBaseNum);
			//}
			// ######################### Debug information ###########################

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
	int returnCode;
	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;
	readseqBlockArr = deBruijnGraph->readSet->readseqBlockArr;

//	if(tmp_kmers[0])
//	{
//		rid_pos_table = tmp_kmers[0]->ppos;
//		posNum = tmp_kmers[0]->arraysize;
//		for(i=0; i<posNum; i++)
//		{
//			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
//			{
//				returnCode = validReadPair(&dtReadPaired, rid_pos_table[i].rid, decisionTable, dtRowHashtable);
//				if(returnCode==YES)
//				{
//					(*occNum) ++;
//					printf("kkkkkkkkkkkkkkkkk!\n");
//				}else if(returnCode==ERROR)
//				{
//					printf("In %s(), cannot valid read pair, error!\n", __func__);
//					return FAILED;
//				}
//			}
//			//rid_pos_table[i].reserved = 0;
//		}
//	}

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
					returnCode = validReadPair(&dtReadPaired, rid_pos_table[i].rid, kmerSize, seqLen, itemNumContigArr, decisionTable, dtRowHashtable);
					if(returnCode==YES)
					{
						validFlag = YES;
					}else if(returnCode==NO)
					{
						validFlag = NO;
					}else if(returnCode==ERROR)
					{
						printf("In %s(), cannot valid read pair, error!\n", __func__);
						return FAILED;
					}

					if(validFlag==YES)
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
			}
			//rid_pos_table[i].reserved = 0;
		}
	}
#endif

	for(i=0; i<itemNumDecisionTable; i++) decisionTable[i].reserved = 0;

	return SUCCESSFUL;
}

/**
 * Compute the occurrence number for locked paired reads in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeKmerOccNumLockedByPE(int32_t *occNum, kmertype *tmp_kmers[2], int32_t baseInt_kmer, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag)
{
	assemblingreadtype *this_assemblingRead;
	ridpostype *rid_pos_table;
	readBlock_t *readBlockArr;
	read_t *pRead;

	int32_t i, properRow, posNum, basePos, baseInt_read;
	int32_t entryIDReadseq, entryPosReadseq, baseNumEntryReadseq;
	int32_t rid, maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3;
	assemblingreadtype *dtReadPaired;

	*occNum = 0;
	for(i=0; i<itemNumDecisionTable; i++)
	{
		if(decisionTable[i].locked)
		{
			//if(decisionTable[i].matedFlag==YES && decisionTable[i].reserved==0)
			//if(decisionTable[i].matedFlag==YES && decisionTable[i].reserved==0 && decisionTable[i].unmatchBaseNum<=3)
			if(decisionTable[i].matedFlag==YES && decisionTable[i].reserved==0
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
	int returnCode;
	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;

//	if(tmp_kmers[0])
//	{
//		rid_pos_table = tmp_kmers[0]->ppos;
//		posNum = tmp_kmers[0]->arraysize;
//		for(i=0; i<posNum; i++)
//		{
//			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
//			{
//				returnCode = validReadPair(&dtReadPaired, rid_pos_table[i].rid, decisionTable, dtRowHashtable);
//				if(returnCode==YES)
//				{
//					(*occNum) ++;
//					printf("kkkkkkkkkkkkkkkkk!\n");
//				}else if(returnCode==ERROR)
//				{
//					printf("In %s(), cannot valid read pair, error!\n", __func__);
//					return FAILED;
//				}
//			}
//			//rid_pos_table[i].reserved = 0;
//		}
//	}

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
				if(existReadWithPosInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
				{
					returnCode = validReadPair(&dtReadPaired, rid_pos_table[i].rid, kmerSize, seqLen, itemNumContigArr, decisionTable, dtRowHashtable);
					if(returnCode==YES)
					{
						(*occNum) ++;
					}else if(returnCode==ERROR)
					{
						printf("In %s(), cannot valid read pair, error!\n", __func__);
						return FAILED;
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
 * Compute the occNum according to long K-mers.
 * 	@return:
 * 		If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeLongKmerOccNumByPE(int32_t *occNum, kmertype *tmp_kmers[2],  int32_t baseInt_kmer, int32_t length_k, int32_t successiveAppearBaseNum, int32_t ignoreEndBaseFlag)
{
	assemblingreadtype *this_assemblingRead = NULL;
	ridpostype *rid_pos_table;
	readBlock_t *readBlockArr;
	read_t *pRead;

	int32_t i, properRow, posNum, basePos, baseInt_read;
	int32_t entryIDReadseq, entryPosReadseq, baseNumEntryReadseq;
	int32_t rid, maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3;
	assemblingreadtype *dtReadPaired;

	*occNum = 0;

	for(i=0; i<itemNumDecisionTable; i++)
	{
		if((decisionTable[i].matedFlag==YES && decisionTable[i].reserved==0)
				//&& (decisionTable[i].unmatchBaseNum<=maxUnmatchBaseNumPerRead)
				&& (decisionTable[i].unmatchBaseNum<=1)
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
	int returnCode;
	readBlockArr = deBruijnGraph->readSet->readBlockArr;
	maxItemNumPerReadBlock = deBruijnGraph->readSet->maxItemNumPerReadBlock;

	//if(length_k<kmerSize+longKmerStepSize)
	if(length_k==kmerSize)
	{
//		if(tmp_kmers[0])
//		{
//			rid_pos_table = tmp_kmers[0]->ppos;
//			posNum = tmp_kmers[0]->arraysize;
//			for(i=0; i<posNum; i++)
//			{
//				if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
//				{
//					returnCode = validReadPair(&dtReadPaired, rid_pos_table[i].rid, decisionTable, dtRowHashtable);
//					if(returnCode==YES)
//					{
//						(*occNum) ++;
//					}else if(returnCode==ERROR)
//					{
//						printf("In %s(), cannot valid read pair, error!\n", __func__);
//						return FAILED;
//					}
//				}
//				//rid_pos_table[i].reserved = 0;
//			}
//		}

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
					if(existReadWithPosInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
					{
						returnCode = validReadPair(&dtReadPaired, rid_pos_table[i].rid, kmerSize, seqLen, itemNumContigArr, decisionTable, dtRowHashtable);
						if(returnCode==YES)
						{
							(*occNum) ++;
						}else if(returnCode==ERROR)
						{
							printf("In %s(), cannot valid read pair, error!\n", __func__);
							return FAILED;
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
 * Check whether the read has valid read pair.
 *  @return:
 *  	(1) If the read has valid read pair, return YES; if invalid, return NO.
 *  	(2) If errors occurred, return ERROR.
 */
short validReadPair(assemblingreadtype **dtReadPaired, uint64_t readID, int32_t overlapSize, int32_t seqLen, int64_t contigNodesNum, assemblingreadtype *decisionTable, dtRowIndex_t **dtRowHashtable)
{
	uint64_t readID_paired;
	int32_t validReadOrient, validReadOrient_paired;
	PERead_t *pPERead;
	double fragSize, fragSizeDif;

	*dtReadPaired = NULL;

	// generate the paired readID
	if(readID%2==1)
	{ // odd --> even
		readID_paired = readID + 1;
	}else
	{ // even --> odd
		readID_paired = readID - 1;
	}

	if(getReadFromPEHashtable(&pPERead, readID_paired)==FAILED)
	{
		printf("In %s(), cannot get the paired read %lu, error!\n", __func__, readID_paired);
		return ERROR;
	}

	if(pPERead)
	{
//		fragSize = contigNodesNum - overlapSize - pPERead->cpos + seqLen;
//		fragSizeDif = fragSize - meanSizeInsert;
//		if(fabs(fragSizeDif)<3*standardDev && (pPERead->seqlen>0.5*readLen && seqLen>40))
//		//if(fabs(fragSizeDif)<3*standardDev)
//		{
//			if(localContigID==52)
//			{
//				printf("localContigID=%ld, contigNodesNum=%ld, itemNumDT=%d, rid=%lu, cpos=%d, fragSize=%.2f, fragSizeDif=%.2f, seqlen_Paired=%d, seqLen=%d\n", localContigID, itemNumContigArr, itemNumDecisionTable, readID, pPERead->cpos, fragSize, fragSizeDif, pPERead->seqlen, seqLen);
//			}

//			if(itemNumDecisionTable>5000 && pPERead->seqlen<0.6*readLen)
//			{
//				printf("---------- localContigID=%ld, contigNodesNum=%ld, itemNumDT=%d, rid=%lu, cpos=%d, fragSize=%.2f, fragSizeDif=%.2f, seqlen_Paired=%d, seqLen=%d\n", localContigID, itemNumContigArr, itemNumDecisionTable, readID, pPERead->cpos, fragSize, fragSizeDif, pPERead->seqlen, seqLen);
//				return NO;
//			}

			return YES;
//		}else
//		{
//			//printf("============ localContigID=%ld, contigNodesNum=%ld, itemNumDT=%d, rid=%lu, cpos=%d, fragSize=%.2f, fragSizeDif=%.2f, seqlen_Paired=%d, seqLen=%d\n", localContigID, itemNumContigArr, itemNumDecisionTable, readID, pPERead->cpos, fragSize, fragSizeDif, pPERead->seqlen, seqLen);
//			return NO;
//		}
	}else if(shortInsertFlag==YES)
	{
		// get the exist read
		if(getExistReadInDT(dtReadPaired, readID_paired, decisionTable, dtRowHashtable)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the exist read %lu in decision table, error!\n", __LINE__, __func__, readID_paired);
			return FAILED;
		}

		//if((*dtReadPaired) && (*dtReadPaired)->orientation==ORIENTATION_PLUS && (*dtReadPaired)->unmatchBaseNum==0)
		if((*dtReadPaired) && (*dtReadPaired)->orientation==ORIENTATION_PLUS && (*dtReadPaired)->unmatchBaseNum==0)
		{
			fragSize = contigNodesNum - overlapSize - (*dtReadPaired)->firstContigPos + seqLen;
			fragSizeDif = fragSize - meanSizeInsert;

			//if(fragSizeDif>=-2*standardDev)
			if(fragSizeDif>=-2*standardDev || fragSizeDif>-0.1*meanSizeInsert)
				return YES;
			else
				return NO;
		}else
			return NO;
	}

	return NO;
}

/**
 * Trim contig nodes at tail.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short trimContigTail(int64_t *successContigIndex, int64_t *contigNodesNum, int32_t trimLen, int assemblyRound)
{
	int64_t i, j, k, divContigIndex, startCheckContigIndex, ridposnum, newSuccessContigIndex;
	successRead_t *ridposorient;
	char seq[MAX_READ_LEN_IN_BUF+1] = {0}, reversed_seq[MAX_READ_LEN_IN_BUF+1] = {0};
	int seq_len;

	// ######################### Debug information ##############################
	//printf("Before trim tail: contigNodesNum=%d, assemblyRound=%d\n", *contigNodesNum, assemblyRound);
	// ######################### Debug information ##############################

	divContigIndex = (*contigNodesNum) - trimLen + (contigArr[0].index - 1);
	startCheckContigIndex = (*contigNodesNum) - MAX_READ_LEN_IN_BUF + 1 + (contigArr[0].index - 1);
	if(startCheckContigIndex<=contigArr[0].index)
		startCheckContigIndex = contigArr[0].index;

	seq_len = 0;
	for(i=startCheckContigIndex-1; i<(*contigNodesNum); i++, seq_len++)
	{
		switch(contigArr[i].base)
		{
			case 0: seq[seq_len] = 'A'; break;
			case 1: seq[seq_len] = 'C'; break;
			case 2: seq[seq_len] = 'G'; break;
			case 3: seq[seq_len] = 'T'; break;
			default:
				printf("line=%d, In %s(), error base int: %d\n", __LINE__, __func__, contigArr[i].base);
				return FAILED;
		}
	}
	seq[seq_len] = '\0';

	if(getReversedSeq(reversed_seq, seq, seq_len)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the reverse contig sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{ // the first assembly round
		for(i=startCheckContigIndex-1; i<(*contigNodesNum); i++)
		{
			if(contigArr[i].index>divContigIndex && contigArr[i].ridposnum>0)
			{ // remove the reads from the contig array
				contigArr[i].ridposnum = 0;
				free(contigArr[i].pridposorientation);
				contigArr[i].pridposorientation = NULL;
			}
		}
	}else
	{ // the second assembly round
		for(i=startCheckContigIndex-1; i<(*contigNodesNum); i++)
		{
			if(contigArr[i].ridposnum>0)
			{
				ridposnum = contigArr[i].ridposnum;
				ridposorient = contigArr[i].pridposorientation;
				for(j=0; j<ridposnum; j++)
				{
					if(contigArr[i].index+ridposorient[j].matchlen-1>divContigIndex)
					{
						// remove the read from the contig node list
						for(k=j+1; k<ridposnum; k++)
						{
							if(memcpy(ridposorient+k-1, ridposorient+k, sizeof(successRead_t))==NULL)
							{
								printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}

						ridposnum --;
						continue;
					}
				}

				contigArr[i].ridposnum = ridposnum;
				if(contigArr[i].ridposnum==0)
				{ // free the memory
					free(contigArr[i].pridposorientation);
					contigArr[i].pridposorientation = NULL;
				}else if(contigArr[i].ridposnum<0)
				{ // error
					printf("line=%d, In %s(), ridposnum=%d, error!\n", __LINE__, __func__, contigArr[i].ridposnum);
					return FAILED;
				}
			}
		}
	}

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{ // the first assembly round
		// get newSuccessContig
		newSuccessContigIndex = -1;
		for(i=startCheckContigIndex-1; i<(*contigNodesNum); i++)
			if(contigArr[i].ridposnum>0 && contigArr[i].index>newSuccessContigIndex)
				newSuccessContigIndex = contigArr[i].index;

		// ##################### Debug information #######################
		if(newSuccessContigIndex>divContigIndex)
		{
			printf("line=%d, In %s(), newSuccessContigIndex=%ld > divContigIndex=%ld, error!\n", __LINE__, __func__, newSuccessContigIndex, divContigIndex);
			return FAILED;
		}
		// ##################### Debug information #######################
	}else
	{ // the second assembly round
		// get newSuccessContig
		newSuccessContigIndex = -1;
		for(i=startCheckContigIndex-1; i<(*contigNodesNum); i++)
		{
			if(contigArr[i].ridposnum>0)
			{
				ridposnum = contigArr[i].ridposnum;
				ridposorient = contigArr[i].pridposorientation;
				for(j=0; j<ridposnum; j++)
				{
					if(contigArr[i].index+ridposorient[j].matchlen-1>newSuccessContigIndex)
					{
						newSuccessContigIndex = contigArr[i].index + ridposorient[j].matchlen - 1;
					}
				}
			}
		}

		// ##################### Debug information #######################
		if(newSuccessContigIndex>divContigIndex)
		{
			printf("line=%d, In %s(), newSuccessContigIndex=%ld > divContigIndex=%ld, error!\n", __LINE__, __func__, newSuccessContigIndex, divContigIndex);
			return FAILED;
		}
		// ##################### Debug information #######################
	}

	// free the memory of tail contig nodes
	if(newSuccessContigIndex>0)
	{
		// delete the contig nodes at tail
		free(contigArr+newSuccessContigIndex+1);
		*contigNodesNum = newSuccessContigIndex - (contigArr[0].index - 1);
		*successContigIndex = newSuccessContigIndex;
	}else
	{
		*successContigIndex = -1;
		*contigNodesNum = 0;
	}

	// ######################### Debug information ##############################
	//printf("After trim tail: contigNodesNum=%d, assemblyRound=%d\n", *contigNodesNum, assemblyRound);
	// ######################### Debug information ##############################

	return SUCCESSFUL;
}

/**
 * Compute the maximal gap size in contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeGapSizeInContig(int *gapSize, contigtype *contigArray, int64_t contigNodesNum, int assemblyRound)
{
	int64_t i, tmp_GapLen, maxLen, startIndex, endIndex;

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{
		startIndex = contigNodesNum - 1.5*readLen + 1;
		if(startIndex<=0)
			startIndex = 1;
		else
		{
			for(i=startIndex-1; i>0; i--)
				if(contigArray[i].ridposnum>0)
					break;
			if(i==0)
				startIndex = 1;
			else
				startIndex = i + 1;
		}
		endIndex = contigNodesNum;
	}else
	{
		endIndex = contigNodesNum - 0.85 * readLen;
		if(endIndex<=0)
			endIndex = 1;
		startIndex = endIndex - 1.5*readLen + 1;
		if(startIndex<=0)
			startIndex = 1;
		else
		{
			for(i=startIndex-1; i>0; i--)
				if(contigArray[i].ridposnum>0)
					break;
			if(i==0)
				startIndex = 1;
			else
				startIndex = i + 1;
		}
	}

	tmp_GapLen = 0;
	maxLen = 0;
	for(i=startIndex-1; i<endIndex; i++)
	{
		if(contigArr[i].ridposnum==0)
		{
			tmp_GapLen ++;
		}else
		{
			if(tmp_GapLen>maxLen)
				maxLen = tmp_GapLen;
			tmp_GapLen = 0;
		}
	}

	if(tmp_GapLen>maxLen)
		maxLen = tmp_GapLen;

	*gapSize = maxLen;

	return SUCCESSFUL;
}


/**
 * Confirm the paired-end navigation.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short confirmNaviPE(int32_t *this_naviSuccessFlag, int32_t *occsArrayPE, int32_t *occsIndexArrayPE, contigPath_t *contigPath, graphtype *graph)
{
	int32_t i, j, baseNumArray[5], baseNumIndexArray[4], kmerBaseInt, naviBaseInt, naviPathLen;
	int32_t maxBaseInt, secBaseInt, maxPathLen, secPathLen;
	double maxBaseRatio;
	char *naviPathseq, *maxPathseq, *secPathseq;

	if(occsIndexArrayPE[1]==-1 && occsArrayPE[occsIndexArrayPE[0]]==2)
	{
		// fill the ratio array
		if(getBaseNumContigPath(baseNumArray, baseNumIndexArray, contigPath)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute the base ratio by contigPath, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// check the navigation by the ratio of the base_PE
		kmerBaseInt = kmerSeqIntAssembly[entriesPerKmer-1] & 3;
		if(baseNumIndexArray[0]!=-1 && kmerBaseInt!=baseNumIndexArray[0] && contigPath->naviPathItem && contigPath->itemNumPathItemList>=2 && kmer_len<0.3*readLen)
		{
			naviPathseq = contigPath->naviPathItem->contigPathStr + contigPath->startRowNewBase;
			naviPathLen = contigPath->naviPathItem->contigPathLen - contigPath->startRowNewBase;
			if(naviPathLen>0)
			{
				switch(naviPathseq[0])
				{
					case 'A': naviBaseInt = 0; break;
					case 'C': naviBaseInt = 1; break;
					case 'G': naviBaseInt = 2; break;
					case 'T': naviBaseInt = 3; break;
					default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, naviPathseq[0]); return FAILED;
				}

				if(naviBaseInt!=kmerBaseInt)
				{
					maxBaseRatio = (double)baseNumArray[baseNumIndexArray[0]] / baseNumArray[4];
					if(maxBaseRatio>=0.9)
					{
						//printf("############## localContigID=%ld, assemblyRound=%d, itemNumContigArr=%ld, kmer_len=%d, occsArrayPE:(%d,%d,%d,%d), maxBaseRatio=%.4f\n", localContigID, assemblyRound, itemNumContigArr, kmer_len, occsArrayPE[0], occsArrayPE[1], occsArrayPE[2], occsArrayPE[3], maxBaseRatio);
						//outputContigPath(contigPath, YES);

						kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] >> 2) << 2) | baseNumIndexArray[0];

						kmers[0] = getKmer(kmerSeqIntAssembly, graph);
						kmers[1] = getReverseKmer(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, graph);
					}
				}
			}
		}

		else if(contigPath->naviPathItem==NULL && contigPath->naviSuccessSize<readLen && contigPath->itemNumPathItemList>2 && contigPath->maxPathItem->supportReadsNum>0 && contigPath->secPathItem->supportReadsNum>0 && baseNumArray[4]>0)
		{
			maxPathseq = contigPath->maxPathItem->contigPathStr + contigPath->startRowNewBase;
			secPathseq = contigPath->secPathItem->contigPathStr + contigPath->startRowNewBase;
			maxPathLen = contigPath->maxPathItem->contigPathLen - contigPath->startRowNewBase;
			secPathLen = contigPath->secPathItem->contigPathLen - contigPath->startRowNewBase;
			if(maxPathLen>0 && secPathLen>0 && maxPathseq[0]==secPathseq[0])
			{
				switch(maxPathseq[0])
				{
					case 'A': maxBaseInt = 0; break;
					case 'C': maxBaseInt = 1; break;
					case 'G': maxBaseInt = 2; break;
					case 'T': maxBaseInt = 3; break;
					default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, maxPathseq[0]); return FAILED;
				}
				switch(secPathseq[0])
				{
					case 'A': secBaseInt = 0; break;
					case 'C': secBaseInt = 1; break;
					case 'G': secBaseInt = 2; break;
					case 'T': secBaseInt = 3; break;
					default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, secPathseq[0]); return FAILED;
				}


				maxBaseRatio = (double)baseNumArray[baseNumIndexArray[0]] / baseNumArray[4];
				//if(kmerBaseInt!=baseNumIndexArray[0] && maxBaseRatio>0.8 && kmer_len<0.3*readLen)
				if(kmerBaseInt!=baseNumIndexArray[0] && maxBaseInt!=kmerBaseInt && maxBaseRatio>0.8)
				{

					//printf("#*##*##*##*## localContigID=%ld, assemblyRound=%d, itemNumContigArr=%ld, kmer_len=%d, occsArrayPE:(%d,%d,%d,%d), maxBaseRatio=%.4f\n", localContigID, assemblyRound, itemNumContigArr, kmer_len, occsArrayPE[0], occsArrayPE[1], occsArrayPE[2], occsArrayPE[3], maxBaseRatio);
					//outputContigPath(contigPath, YES);


	//				kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] >> 2) << 2) | baseNumIndexArray[0];
	//
	//				kmers[0] = getKmer(kmerSeqIntAssembly, graph);
	//				kmers[1] = getReverseKmer(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, graph);

					*this_naviSuccessFlag = NAVI_FAILED;
				}
			}
		}

		else if(contigPath->naviPathItem && contigPath->itemNumPathItemList>2 && baseNumArray[4]>0)
		{
			if(baseNumIndexArray[1]!=-1 && (double)baseNumArray[baseNumIndexArray[1]]/baseNumArray[baseNumIndexArray[0]]>0.8)
			{
				naviPathseq = contigPath->naviPathItem->contigPathStr + contigPath->startRowNewBase;
				naviPathLen = contigPath->naviPathItem->contigPathLen - contigPath->startRowNewBase;
				if(naviPathLen>0)
				{
					switch(naviPathseq[0])
					{
						case 'A': naviBaseInt = 0; break;
						case 'C': naviBaseInt = 1; break;
						case 'G': naviBaseInt = 2; break;
						case 'T': naviBaseInt = 3; break;
						default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, naviPathseq[0]); return FAILED;
					}

					if(naviBaseInt!=kmerBaseInt)
					{

						//printf("-*--*--*--*-- localContigID=%ld, assemblyRound=%d, itemNumContigArr=%ld, kmer_len=%d, occsArrayPE:(%d,%d,%d,%d), maxBaseRatio=%.4f\n", localContigID, assemblyRound, itemNumContigArr, kmer_len, occsArrayPE[0], occsArrayPE[1], occsArrayPE[2], occsArrayPE[3], maxBaseRatio);
						//outputContigPath(contigPath, YES);



						*this_naviSuccessFlag = NAVI_FAILED;
					}

				}
			}
		}

	}

	return SUCCESSFUL;
}

/**
 * Compute the base ratio by contigPath.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getBaseNumContigPath(int32_t *baseNumArray, int32_t *baseNumIndexArray, contigPath_t *contigPath)
{
	int32_t i, j, checkNum, pathLen, baseInt, maxValue, secValue, maxBaseIndex, secBaseIndex;
	int32_t tmpOccIndexArray[2], itemNumTmpOccIndexArray;
	char *pathseq;
	contigPathItem_t *pathItem;

	if(contigPath->itemNumPathItemList>10)
		checkNum = 10;
	else
		checkNum = contigPath->itemNumPathItemList;

	for(i=0; i<4; i++) { baseNumArray[i] = 0; baseNumIndexArray[i] = -1; }

	i = 0;
	baseNumArray[4] = 0;
	pathItem = contigPath->contigPathItemList;
	while(pathItem)
	{
		pathseq = pathItem->contigPathStr + contigPath->startRowNewBase;
		pathLen = pathItem->contigPathLen - contigPath->startRowNewBase;
		if(pathLen>0 && pathItem->supportReadsNum>0)
		{
			switch(pathseq[0])
			{
				case 'A': baseInt = 0; break;
				case 'C': baseInt = 1; break;
				case 'G': baseInt = 2; break;
				case 'T': baseInt = 3; break;
				default: printf("line=%d, In %s(), invalid base %c, error!\n", __LINE__, __func__, pathseq[0]); return FAILED;
			}
			baseNumArray[baseInt] += pathItem->supportReadsNum;
			baseNumArray[4] += pathItem->supportReadsNum;
		}

		i ++;
		if(i>checkNum)
			break;

		pathItem = pathItem->nextPathItem;
	}

	if(baseNumArray[4]>0)
	{
		maxValue = 0, maxBaseIndex = -1, secValue = 0, secBaseIndex = -1;
		for(i=0; i<4; i++)
		{
			if(maxValue<baseNumArray[i])
			{
				secValue = maxValue;
				secBaseIndex = maxBaseIndex;
				maxValue = baseNumArray[i];
				maxBaseIndex = i;
			}else if(secValue<baseNumArray[i])
			{
				secValue = baseNumArray[i];
				secBaseIndex = i;
			}
		}
		baseNumIndexArray[0] = maxBaseIndex;
		baseNumIndexArray[1] = secBaseIndex;

		if(secondOccSE>0)
		{
			itemNumTmpOccIndexArray = 0;
			for(i=0; i<4; i++)
			{
				if(i!=maxBaseIndex && i!=secBaseIndex && baseNumArray[i]>0)
					tmpOccIndexArray[itemNumTmpOccIndexArray++] = i;
			}

			if(itemNumTmpOccIndexArray==1)
			{
				baseNumIndexArray[2] = baseNumArray[ tmpOccIndexArray[itemNumTmpOccIndexArray-1] ];
			}else if(itemNumTmpOccIndexArray==2)
			{
				if(tmpOccIndexArray[0]>=tmpOccIndexArray[1])
				{
					baseNumIndexArray[2] = baseNumArray[ tmpOccIndexArray[0] ];
					baseNumIndexArray[3] = baseNumArray[ tmpOccIndexArray[1] ];
				}else
				{
					baseNumIndexArray[2] = baseNumArray[ tmpOccIndexArray[1] ];
					baseNumIndexArray[3] = baseNumArray[ tmpOccIndexArray[0] ];
				}
			}else if(itemNumTmpOccIndexArray>2 || itemNumTmpOccIndexArray<0)
			{
				printf("line=%d, In %s(), invalid itemNumTmpOccIndexArray=%d, error!\n", __LINE__, __func__, itemNumTmpOccIndexArray);
				return FAILED;
			}
		}
	}
	else
	{
		//printf("line=%d, In %s(), localContigID=%ld, itemNumContigArr=%ld, totalNum=%d, occsNumPE(%d,%d,%d,%d)!\n", __LINE__, __func__, localContigID, itemNumContigArr, baseNumArray[4], occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
		//return FAILED;
	}

	return SUCCESSFUL;
}
