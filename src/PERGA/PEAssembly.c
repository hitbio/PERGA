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
	int i, turnContigIndex, validFlag, tmp_gapSize;
	//int64_t localContigID;
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
	correctionAllowed = NO;

	validFlag = NO;

	for(assemblyCycle=1; assemblyCycle<=2 && validFlag==NO; assemblyCycle++)
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
			successContigIndex = -1;
			assemblyRound = FIRST_ROUND_ASSEMBLY;  // first round assembly
			lockedReadsNum = 0;
			this_successReadNum = 0;
			readsNumInPEHashArr = 0;
			regLenPEHash = 0;
			turnContigIndex = 0;
			allowedUpdatePEHashArrFlag = YES;
			localContigID ++;


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

#if (DEBUG_EST_CONTIG_CHECK==YES)
				// ############################ Debug information ##############################
				if(localContigID==1 && itemNumContigArr>=15003 && assemblyRound==FIRST_ROUND_ASSEMBLY)
				{
					printf("localContigID=%ld, contigID=%d, itemNumContigArr=%ld\n", localContigID, contigsNum+1, itemNumContigArr);
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
						if(getNextKmerBySE(itemNumContigArr)==FAILED)
						{
							printf("line=%d, In %s(), localContigID=%ld cannot get next kmer, error!\n", __LINE__, __func__, localContigID);
							return FAILED;
						}

						if(successContigIndex>0 && itemNumContigArr-successContigIndex > 0.5*readLen)  // added 2012-11-08
						{
							naviSuccessFlag = NAVI_FAILED;
							kmers[0] = kmers[1] = NULL;
						}
						//else if(itemNumContigArr<3*readLen && (successContigIndex>0 && (itemNumContigArr-successContigIndex > 0.3*readLen || itemNumContigArr-successContigIndex > 15)))  // added 2012-11-11, deleted 2012-11-23
						//else if(averKmerOcc>15 && itemNumContigArr<3*readLen && (successContigIndex>0 && (itemNumContigArr-successContigIndex > 0.3*readLen || itemNumContigArr-successContigIndex > 15)))  // added 2012-11-23, deleted 2012-12-29
						//{
						//	naviSuccessFlag = NAVI_FAILED;
						//	kmers[0] = kmers[1] = NULL;
						//}
						//else if(itemNumContigArr<3*readLen && secondOccSE>15*minKmerOccSE && secondOccSE/maxOccSE>0.6)  // added 2012-11-13
						//else if(itemNumContigArr<3*readLen && secondOccSE>5*minKmerOccSE && secondOccSE/maxOccSE>0.75 && readsNumRatio<0.4)  // added 2012-11-16, deleted 2012-11-28
						else if(itemNumContigArr<3*readLen && secondOccSE>3*minKmerOccSE && secondOccSE/maxOccSE>0.6)  // added 2012-11-28
						{
							naviSuccessFlag = NAVI_FAILED;
							kmers[0] = kmers[1] = NULL;
						}

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
					}
				}else
				{
					navigationFlag = NAVI_SE_FLAG;
					if(getNextKmerBySE(itemNumContigArr)==FAILED)
					{
						printf("line=%d, In %s(), localContigID=%ld, cannot get next kmer, error!\n", __LINE__, __func__, localContigID);
						return FAILED;
					}

					if(successContigIndex>0 &&  itemNumContigArr-successContigIndex > 0.5*readLen)  // added 2012-11-08
					{
						naviSuccessFlag = NAVI_FAILED;
						kmers[0] = kmers[1] = NULL;
					}
					//else if(itemNumContigArr<3*readLen && (successContigIndex>0 && (itemNumContigArr-successContigIndex > 0.3*readLen || itemNumContigArr-successContigIndex > 15)))  // added 2012-11-11, deleted 2012-11-23
					//else if(averKmerOcc>15 && itemNumContigArr<3*readLen && (successContigIndex>0 && (itemNumContigArr-successContigIndex > 0.3*readLen || itemNumContigArr-successContigIndex > 15)))  // added 2012-11-23, deleted 2012-12-29
					//{
					//	naviSuccessFlag = NAVI_FAILED;
					//	kmers[0] = kmers[1] = NULL;
					//}
					//else if(itemNumContigArr<3*readLen && secondOccSE>15*minKmerOccSE && secondOccSE/maxOccSE>0.6)  // added 2012-11-13
					//else if(itemNumContigArr<3*readLen && secondOccSE>5*minKmerOccSE && secondOccSE/maxOccSE>0.75 && readsNumRatio<0.4)  // added 2012-11-16, deleted 2012-11-28
					else if(itemNumContigArr<3*readLen && secondOccSE>3*minKmerOccSE && secondOccSE/maxOccSE>0.6)  // added 2012-11-28
					{
						naviSuccessFlag = NAVI_FAILED;
						kmers[0] = kmers[1] = NULL;
					}

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
				}

#if (DEBUG_EST_CONTIG_CHECK==YES)
				// ############################ Debug information ##############################
				if(localContigID==1 && itemNumContigArr>=15003 && assemblyRound==FIRST_ROUND_ASSEMBLY)
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

				//if(localContigID==2)
				//	outputReadsInDecisionTableToFile(decisionTable, itemNumDecisionTable, (int)localContigID, itemNumContigArr);
				//outputFailedReadsInDecisionTable(decisionTable, itemNumDecisionTable, (int)localContigID, itemNumContigArr);

				// update the finished reads in decision table, and record the successful reads into successful reads array
				if(updateFinishedReadsInDecisionTable()==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, itemNumContigArr=%ld, cannot update finished reads in decision table, error!\n", __LINE__, __func__, localContigID, itemNumContigArr);
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

					if((averOccNumNaviOccQueue<2.5*minKmerOccSE) || (readsNumRatio>2 || readsNumRatio<0.3) || ((navigationFlag==NAVI_PE_FLAG && secondOccPE>=minKmerOccPE && secondOccPE/maxOccPE>=0.2) || (navigationFlag==NAVI_SE_FLAG && secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>=0.15)))
					{
						number_of_overlap_less_than_threshold ++;

						printf("===localContigID=%ld, assemblyRound=%d, itemNumContigArr=%d, itemNumDecisionTable=%d\n", localContigID, assemblyRound, itemNumContigArr, itemNumDecisionTable);
						printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
						printf("\tdistance=%d, readsNumRatio=%.2f, averOccNumNaviOccQueue=%.2f\n", itemNumContigArr-successContigIndex, readsNumRatio, averOccNumNaviOccQueue);


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

				if(updateContigtailnodes(contigArr, successContigIndex, &itemNumContigArr)==FAILED)
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
						break;
					}
				}else
				{
					successReadNum -= this_successReadNum;

					// clean the contig array
					cleanContigArray(contigArr, &itemNumContigArr);
				}
			}else
			{
				// ############################ Debug information ##############################
				//if(itemNumContigArr>=CONTIG_LEN_THRESHOLD)
				//{
				//	printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, itemNumContigArr=%ld, successContigIndex<=0, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, itemNumContigArr);
				//	return FAILED;
				//}
				// ############################ Debug information ##############################

				// clean the contig array
				cleanContigArray(contigArr, &itemNumContigArr);
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

		} //end while(kmerIndex < TABLE_SIZE_DE_BRUIJN)

	}


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
short getNextKmerByMix(int contigNodesNum, int assemblyRound)
{
	unsigned int tmp_kmerseqInt, tmp_kmerseqPE, tmp_kmerseqSE;
	kmertype *tmp_kmers[2];
	int32_t i, tmp_gapSize, naviPE, naviSE, naviMalignPE, naviMalignSE, incorrectNumPE, incorrectNumSE;
	double averOccNumNaviOccQueue;

//	if(itemNumDecisionTable>MAX_DECISION_TABLE_SIZE_HTRES)
//	{
//		naviSuccessFlag = NAVI_FAILED;
//		kmers[0] = kmers[1] = NULL;
//		return SUCCESSFUL;
//	}

	tmp_gapSize = -1;
	naviPE = naviSE = NAVI_UNUSED;
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
	naviMalignPE = NAVI_UNUSED;
	incorrectNumPE = -1;
	if(readsNumRatio<0.3*minReadsNumRatioThres)		// added 2013-04-01
	{
		naviSuccessFlag = NAVI_FAILED;
		return SUCCESSFUL;
	}
	else if(secondOccPE>0)								// deleted 2013-02-26
	//if(naviSuccessFlag==NAVI_SUCCESS && secondOccPE>0 && successContigIndex>0)		// added 2013-02-26
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

//		if(secondOccPE>10 && readsNumRatio>6 && tmp_gapSize>10)
//			naviPE = NAVI_FAILED;

		//else if(tmp_gapSize>40)
		//	naviPE = NAVI_FAILED;

		if(secondOccSE/maxOccSE>=0.25 && maxOccSE<0.3*averKmerOcc)
		{
			naviSE = NAVI_FAILED;
		}

		//if((secondOccPE/maxOccPE>0.3 /*&& secondOccPE/maxOccPE<=0.7*/) || secondOccPE>10)
		if(naviPE==NAVI_FAILED || (secondOccPE>10 && readsNumRatio>5))  // best
		{
#if(MALIGN_OUTPUT==YES)
			printf("naviPE=%d\n", naviPE);
#endif
			if(decideByMalignPE(&naviMalignPE, &incorrectNumPE, occsNumIndexPE[0], decisionTable, itemNumDecisionTable, deBruijnGraph)==FAILED)
			{
				printf("line=%d, In %s(), cannot navigate by multiple alignment using paired ends, error!\n", __LINE__, __func__);
				return FAILED;
			}
			naviPE = naviMalignPE;
		}

//		if(naviPE==NAVI_FAILED && svmFeatureArr[1]/svmFeatureArr[0]<0.5)
//			naviPE = NAVI_SUCCESS;

//		if(svmFeatureArr[0]==svmFeatureArr[1] || (svmFeatureArr[1]>=2 && svmFeatureArr[1]/svmFeatureArr[0]>=0.7))
//			naviPE = NAVI_FAILED;
//		else if(tmp_gapSize>50 && svmFeatureArr[1]/svmFeatureArr[0]>=0.5)
//		{
//			naviPE = NAVI_FAILED;
//		}

		naviSuccessFlag = naviPE;

	}
#else
	//if(maxOccPE==secondOccPE)
	if((maxOccPE==secondOccPE) || (secondOccPE>=2 && secondOccPE>=SEC_MAX_OCC_RATIO_SVM*maxOccPE))
	{
		naviSuccessFlag = NAVI_FAILED;
		//kmers[0] = kmers[1] = NULL;
	}
#endif

#if(USING_RULES==YES)
	sumSecondOccPE = 0;
	if(secondOccPE>0)
		for(i=0; i<4; i++) if(i!=maxOccIndexPE) sumSecondOccPE += occsNumPE[i];

	//if(kmers[0] || kmers[1])
	if(naviSuccessFlag==NAVI_SUCCESS)
	{
		//if((successContigIndex>0 && contigNodesNum-successContigIndex >= readLen-kmerSize) && readsNumRatio<0.4) // new added 2012-11-08
		//if((successContigIndex>0 && contigNodesNum-successContigIndex >= readLen-kmerSize)) // new added 2012-11-08
		//if((successContigIndex>0 && contigNodesNum-successContigIndex >= 0.5*readLen)) // new added 2012-11-10
		//if((successContigIndex>0 && contigNodesNum-successContigIndex >= 0.6*readLen)) // new added 2012-11-13, deleted 2012-11-25
		if((averKmerOcc<20 && (successContigIndex>0 && contigNodesNum-successContigIndex >= 0.6*readLen)) || (averKmerOcc>20 && (successContigIndex>0 && contigNodesNum-successContigIndex >= 0.7*readLen))) // new added 2012-11-25
		{
#if (DEBUG_OUTPUT==YES)
			printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
			naviSuccessFlag = NAVI_FAILED;
			kmers[0] = kmers[1] = NULL;
		}

		else if(averKmerOcc<20 && (successContigIndex>0 && contigNodesNum-successContigIndex >= 0.5*readLen) && readsNumRatio<0.3)  // added 2012-11-22
		{
#if (DEBUG_OUTPUT==YES)
			printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
			naviSuccessFlag = NAVI_FAILED;
			kmers[0] = kmers[1] = NULL;
		}

		//if((secondOccPE>0 && secondOccPE/maxOccPE>0.2) || (maxOccPE<minKmerOccPE && secondOccPE>0))
		else if((secondOccPE>0 && secondOccPE/maxOccPE>SECOND_FIRST_OCC_FAILED_RATIO) || (maxOccPE<minKmerOccPE && secondOccPE>0))
		{
			if((readsNumRatio<0.2 || readsNumRatio>3) && (successContigIndex>0 && contigNodesNum-successContigIndex >= 15))
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}else if((maxOccPE>0 && maxOccPE<=minKmerOccPE) && (readsNumRatio<0.3 || readsNumRatio>3))
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}else if(readsNumRatio>2.5 && (successContigIndex>0 && contigNodesNum-successContigIndex >= 20))
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				kmers[0] = kmers[1] = NULL;
			}else if(readsNumRatio>6 && (maxOccPE>maxOccNumFaiedPE && secondOccPE/maxOccPE>0.2))
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}
			else if(readsNumRatio>5)
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}else if(maxOccPE==1 && readsNumRatio < 0.7)
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}
			//else if(secondOccPE>10*minKmerOccPE && secondOccPE/maxOccPE>SECOND_FIRST_OCC_RATIO) // new added 2012-11-10
			else if(secondOccPE>3*minKmerOccPE && secondOccPE/maxOccPE>SECOND_FIRST_OCC_RATIO) // new added 2012-11-16
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}
			else if((maxOccPE>0 && maxOccPE<minKmerOccPE && secondOccPE>0) && maxOccPE-secondOccPE<2) // new added 2012-11-14
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}
			else if((maxOccPE<minKmerOccPE && secondOccPE<0.3*minKmerOccPE && secondOccPE/maxOccPE>0.3))  // added 2012-11-16
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}

		}else
		{
			if(sumSecondOccPE>maxOccPE && secondOccPE>maxOccNumFaiedPE)
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum,readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}
			//else if(sumSecondOccPE/maxOccPE>0.5 && secondOccPE>maxOccNumFaiedPE && readsNumRatio>5)
			else if(sumSecondOccPE/maxOccPE>SECOND_FIRST_OCC_FAILED_RATIO && secondOccPE>maxOccNumFaiedPE && readsNumRatio>7)
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}

			//else if(readsNumRatio>2.5 && secondOccPE>20 && sumSecondOccPE/maxOccPE>0.4) // added 2012-11-08
			else if(readsNumRatio>2.5 && secondOccPE>4*minKmerOccPE && sumSecondOccPE/maxOccPE>0.4) // added 2012-11-08
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}

			//else if(readsNumRatio<0.4 && secondOccPE>10 && sumSecondOccPE/maxOccPE>0.5 && (successContigIndex>0 &&  contigNodesNum-successContigIndex >= 15)) // added 2012-11-10
			//else if(readsNumRatio<0.4 && secondOccPE>3*minKmerOccPE && sumSecondOccPE/maxOccPE>0.5) // added 2012-11-10, deleted 2013-01-19
			else if(secondOccPE>=2*minKmerOccPE && sumSecondOccPE/maxOccPE>0.5) // added 2013-01-19
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}

			//else if((maxOccPE<3*minKmerOccPE && secondOccPE<minKmerOccPE && secondOccPE/maxOccPE>0.3))  // added 2012-11-15
			else if((maxOccPE<minKmerOccPE && secondOccPE<0.3*minKmerOccPE && secondOccPE/maxOccPE>0.3))  // added 2012-11-16
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}
			else if(secondOccPE>20*minKmerOccPE && secondOccPE/maxOccPE>0.2 && readsNumRatio>7)  // added 2012-11-16
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}

			//else if(readsNumRatio>7 && (secondOccPE>maxOccNumFaiedPE && secondOccPE/maxOccPE>0.25))  // added 2012-11-17
			//else if(readsNumRatio>8.5 && (secondOccPE>maxOccNumFaiedPE && secondOccPE/maxOccPE>0.25))  // added 2012-11-17
			//else if(readsNumRatio>8.5 && (secondOccPE>maxOccNumFaiedPE && secondOccPE/maxOccPE>0.15))  // added 2012-11-17
			//else if(readsNumRatio>7.5 && (secondOccPE>maxOccNumFaiedPE && secondOccPE/maxOccPE>0.15))  // added 2012-11-17
			else if(readsNumRatio>7 && (secondOccPE>maxOccNumFaiedPE && secondOccPE/maxOccPE>0.15))  // added 2012-11-18
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}
			//else if(readsNumRatio>7 && (secondOccPE>4*minKmerOccPE && secondOccPE/maxOccPE>0.1))  // added 2012-11-22, deleted 2012-11-23
			//else if(readsNumRatio>6 && (secondOccPE>4*minKmerOccPE && secondOccPE/maxOccPE>0.1))  // added 2012-11-23, deleted 2013-01-20
			else if(readsNumRatio>6 && (secondOccPE>3*minKmerOccPE && secondOccPE/maxOccPE>0.1))  // added 2013-01-20
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}

			else if(readsNumRatio>6.5 && (maxOccPE>4*maxOccNumFaiedPE && secondOccPE>0.5*maxOccNumFaiedPE && secondOccPE/maxOccPE>0.1))  // added 2012-11-20
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}

			else if(secondOccPE>2*minKmerOccPE && secondOccPE/maxOccPE>=0.5)  // added 2012-11-22
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}

			else if(secondOccPE>4*minKmerOccPE && secondOccPE/maxOccPE>0.3)  // added 2013-01-21
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}

			else if(secondOccPE>minKmerOccPE && secondOccPE/maxOccPE>0.15)  // added 2012-11-20
			{
				// compute the maximal gap size in contig tail region
				if(computeGapSizeInContig(&tmp_gapSize, contigArr, contigNodesNum, assemblyRound)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
					return FAILED;
				}

				//if(tmp_gapSize>10)  // added 2012-11-20
				if(tmp_gapSize>=15)  // added 2012-11-20
				{
#if (DEBUG_OUTPUT==YES)
					printf("==line=%d, contigNodesNum=%d, tmp_gapSize=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, tmp_gapSize, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}
			}

			//else if(readsNumRatio<0.2*minReadsNumRatioThres)
			//else if(readsNumRatio<0.2)
			else if(readsNumRatio<minReadsNumRatioThres)
			{
				// compute the maximal gap size in contig tail region
				if(computeGapSizeInContig(&tmp_gapSize, contigArr, contigNodesNum, assemblyRound)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(tmp_gapSize>0.7*readLen)
				//if(tmp_gapSize>0.6*readLen)  // added 2012-11-16
				{
#if (DEBUG_OUTPUT==YES)
					printf("==line=%d, contigNodesNum=%d, tmp_gapSize=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, tmp_gapSize, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}
				else if(readsNumRatio<0.3*minReadsNumRatioThres && tmp_gapSize>0.4*readLen)   // added 2012-11-14
				{
#if (DEBUG_OUTPUT==YES)
					printf("==line=%d, contigNodesNum=%d, tmp_gapSize=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, tmp_gapSize, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}
				else if(tmp_gapSize>0.3*readLen && maxOccPE>5*minKmerOccPE && (successContigIndex>0 &&  contigNodesNum-successContigIndex >= 10))
				{
#if (DEBUG_OUTPUT==YES)
					printf("==line=%d, contigNodesNum=%d, tmp_gapSize=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, tmp_gapSize, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}
			}
		}
	}
#endif

	//if(kmers[0]==NULL && kmers[1]==NULL)
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
		if((successContigIndex>0 && itemNumContigArr-successContigIndex>60))
			naviSuccessFlag = NAVI_FAILED;
		else if(successContigIndex>0 && (readsNumRatio<0.3 && maxOccPE==0)) // added 2013-04-23
		{
			// compute the maximal gap size in contig tail region
			if(computeGapSizeInContig(&tmp_gapSize,contigArr, contigNodesNum, assemblyRound)==FAILED)
			{
				printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(tmp_gapSize>15)
			{
				naviSuccessFlag = NAVI_FAILED;
			}
		}
		else if(secondOccSE>0)								// deleted 2013-02-26
		//if(naviSuccessFlag==NAVI_SUCCESS && secondOccSE>0 && successContigIndex>0)		// added 2013-02-26
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

			if(maxOccIndexPE!=-1 && maxOccIndexPE!=maxOccIndexSE)
			{
				//if(secondOccPE!=maxOccPE)
				{
//						printf("occsNumPE: (%d,%d,%d,%d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
//						printf("occsNumSE: (%d,%d,%d,%d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
					naviSE = NAVI_FAILED;
				}
			}
			//else if(tmp_gapSize>50 && svmFeatureArr[1]/svmFeatureArr[0]>=0.5)
			//else if(tmp_gapSize>=30)
			else if(tmp_gapSize>=30 && secondOccSE/maxOccSE>=0.3)
			{
				naviSE = NAVI_FAILED;
			}
			else if(maxOccPE==0 && secondOccSE/maxOccSE>=0.3 && maxOccSE<0.3*averKmerOcc)
			{
				naviSE = NAVI_FAILED;
			}
			else if(maxOccSE-secondOccSE<2 && readsNumRatio>2)
			{
				naviSE = NAVI_FAILED;
			}
			//else if((secondOccPE>10 && readsNumRatio>6 && tmp_gapSize>10) || (secondOccPE/maxOccPE>=0.7))
			else if((naviMalignPE==NAVI_FAILED) || (secondOccPE>0 && secondOccPE>=maxOccPE*0.7))
			{
				naviSE = NAVI_FAILED;
			}

//			if((secondOccPE/maxOccPE>=0.7) && naviSE!=NAVI_FAILED)
//			{
				if(naviSE==NAVI_FAILED && tmp_gapSize<30)	// best
				//if(naviSE==NAVI_FAILED)
				{
#if(MALIGN_OUTPUT==YES)
					printf("naviSE=%d\n", naviSE);
#endif
					if(decideByMalignSE(&naviMalignSE, incorrectNumPE, maxOccIndexPE, decisionTable, itemNumDecisionTable, deBruijnGraph)==FAILED)
					{
						printf("line=%d, In %s(), cannot navigate by multiple alignment using single ends, error!\n", __LINE__, __func__);
						return FAILED;
					}
					naviSE = naviMalignSE;
				}
//			}


			naviSuccessFlag = naviSE;

//			if((secondOccPE/maxOccPE>0.3 /*&& secondOccPE/maxOccPE<=0.7*/) || secondOccPE>10)
//				decideByMalignPE(decisionTable, itemNumDecisionTable);
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

#if(USING_RULES==YES)
		//if(kmers[0] || kmers[1])
		if(naviSuccessFlag==NAVI_SUCCESS)
		{
			//if(successContigIndex>0 &&  contigNodesNum-successContigIndex > readLen-kmerSize)  // added 2012-11-08
			//if(successContigIndex>0 &&  contigNodesNum-successContigIndex >= 0.6*readLen)  // added 2012-11-13, deleted 2012-11-25
			if((averKmerOcc<20 && (successContigIndex>0 &&  contigNodesNum-successContigIndex >= 0.6*readLen)) || (averKmerOcc>=20 && (successContigIndex>0 &&  contigNodesNum-successContigIndex >= 0.7*readLen)))  // added 2012-11-25
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}

			else if(averKmerOcc<20 && (successContigIndex>0 && contigNodesNum-successContigIndex >= 0.5*readLen) && readsNumRatio<0.3)  // added 2012-11-22
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}

			else if(averKmerOcc<15 && (successContigIndex>0 && contigNodesNum-successContigIndex > 0.4*readLen) && readsNumRatio<0.5)  // added 2012-11-23
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}

			else if(secondOccSE>0)
			{
				sumSecondOccSE = 0;
				for(i=0; i<4; i++) if(i!=maxOccIndexSE) sumSecondOccSE += occsNumSE[i];

				// check the scores and occsNums of PE and SE
				//if(maxOccIndexPE==-1 && maxOccSE>OCCS_NUM_SE_FAILED_PE_FACTOR*averKmerOcc)  //==================================
				//if(maxOccIndexPE==-1 && maxOccSE>maxFirstOcc)
				//if(maxOccIndexPE==-1 && maxOccSE>maxOccNumFaiedPE)
				//if((maxOccIndexPE==-1 && (secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>SECOND_FIRST_OCC_FAILED_RATIO && maxOccSE>maxOccNumFaiedPE)) && navigationID==0)
				//if((maxOccIndexPE==-1 && ((secondOccSE>=minKmerOccSE && maxOccSE>maxOccNumFaiedPE) || maxOccSE<minKmerOccSE)) && navigationID==0)
				//if((maxOccIndexPE==-1 && ((secondOccSE>=minKmerOccSE && (maxOccSE>maxOccNumFaiedPE || secondOccSE/maxOccSE>SECOND_FIRST_OCC_FAILED_RATIO)) || maxOccSE<minKmerOccSE)) && navigationID==0)
				//if((maxOccIndexPE==-1 && ((secondOccSE>=minKmerOccSE && (maxOccSE>maxOccNumFaiedPE || secondOccSE/maxOccSE>0.5)) || maxOccSE<minKmerOccSE)) && navigationID==0)
				//if(((secondOccSE>=minKmerOccSE && (maxOccSE>maxOccNumFaiedPE || secondOccSE/maxOccSE>0.5)) || maxOccSE<minKmerOccSE) && navigationID==0)
				//if((secondOccSE>=minKmerOccSE && (maxOccSE>maxOccNumFaiedPE || secondOccSE/maxOccSE>0.5)) || maxOccSE<minKmerOccSE)
				//if((secondOccSE>=minKmerOccSE && (maxOccSE>maxOccNumFaiedPE || secondOccSE/maxOccSE>SECOND_FIRST_OCC_FAILED_RATIO)) || (secondOccSE>0 && maxOccSE<minKmerOccSE)) // deleted 2012-11-23
				if((secondOccSE>=minKmerOccSE && (secondOccSE/maxOccSE>SECOND_FIRST_OCC_FAILED_RATIO)) || (secondOccSE>0 && maxOccSE<minKmerOccSE))  // added 2012-11-23
				{
					if(maxOccIndexPE==-1)
					{
						if(secondOccSE/maxOccSE>SECOND_FIRST_OCC_FAILED_RATIO)
						{
#if (DEBUG_OUTPUT==YES)
							printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
							naviSuccessFlag = NAVI_FAILED;
							kmers[0] = kmers[1] = NULL;
						}
						else if(maxOccSE==1 && (successContigIndex>0 &&  contigNodesNum-successContigIndex >= 3))
						{
							if(maxOccIndexPE!=maxOccIndexSE)
							{
#if (DEBUG_OUTPUT==YES)
							printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
							naviSuccessFlag = NAVI_FAILED;
							kmers[0] = kmers[1] = NULL;
							}
						}
					}else
					{
						if(maxOccPE>=minKmerOccPE && maxOccSE>=minKmerOccSE && maxOccIndexPE==maxOccIndexSE && secondOccSE/maxOccSE<=SECOND_FIRST_OCC_FAILED_RATIO)
						{
							//if(readsNumRatio>6 && secondOccSE/maxOccSE>0.15)
							//if((maxOccIndexPE!=maxOccIndexSE || (sumSecondOccPE>maxOccPE && sumSecondOccSE>maxOccSE && secondOccPE>maxOccNumFaiedPE && secondOccSE>maxOccNumFaiedPE)) && readsNumRatio>6 && secondOccSE/maxOccSE>0.15)
							if((sumSecondOccPE>maxOccPE && sumSecondOccSE>maxOccSE) && readsNumRatio>5 && secondOccSE/maxOccSE>0.15)
							{
#if (DEBUG_OUTPUT==YES)
								printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
								naviSuccessFlag = NAVI_FAILED;
								kmers[0] = kmers[1] = NULL;
							}
							else if(sumSecondOccPE>0.5*maxOccPE && secondOccPE>maxOccNumFaiedPE && readsNumRatio>5)
							{
#if (DEBUG_OUTPUT==YES)
								printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
								naviSuccessFlag = NAVI_FAILED;
								kmers[0] = kmers[1] = NULL;
							}
							//else if(secondOccSE>10*minKmerOccPE && secondOccPE>15*minKmerOccPE && secondOccPE/maxOccPE>0.2 && readsNumRatio>7)  // added 2012-11-16
							//else if((readsNumRatio>8.5 && secondOccPE>maxOccNumFaiedPE) && secondOccSE>0.4*maxOccNumFaiedPE)  // added 2012-11-17
							//else if((readsNumRatio>7.5 && secondOccPE>maxOccNumFaiedPE) && secondOccSE>0.8*maxOccNumFaiedPE && secondOccPE/maxOccPE>0.15)  // added 2012-11-17
							else if((readsNumRatio>7 && secondOccPE>maxOccNumFaiedPE) && secondOccSE>0.8*maxOccNumFaiedPE && secondOccPE/maxOccPE>0.15)  // added 2012-11-18
							{
#if (DEBUG_OUTPUT==YES)
								printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
								naviSuccessFlag = NAVI_FAILED;
								kmers[0] = kmers[1] = NULL;
							}
							else if(readsNumRatio>6.5 && (maxOccPE>4*maxOccNumFaiedPE && secondOccPE>0.5*maxOccNumFaiedPE && secondOccPE/maxOccPE>0.1) && (maxOccSE>2*maxOccNumFaiedPE && secondOccSE>0.3*maxOccNumFaiedPE && secondOccSE/maxOccSE>0.15))  // added 2012-11-20
							{
#if (DEBUG_OUTPUT==YES)
								printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
								naviSuccessFlag = NAVI_FAILED;
								kmers[0] = kmers[1] = NULL;
							}
						}else
						{
							if(secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>0.5 && maxOccIndexPE!=maxOccIndexSE)
							{
#if (DEBUG_OUTPUT==YES)
								printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
								naviSuccessFlag = NAVI_FAILED;
								kmers[0] = kmers[1] = NULL;
							}else if(secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>SECOND_FIRST_OCC_FAILED_RATIO)
							{
#if (DEBUG_OUTPUT==YES)
								printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
								naviSuccessFlag = NAVI_FAILED;
								kmers[0] = kmers[1] = NULL;
							}
							else if((maxOccIndexPE!=-1 && maxOccIndexPE!=maxOccIndexSE) && (secondOccPE/maxOccPE>SECOND_FIRST_OCC_FAILED_RATIO && secondOccPE>maxOccNumFaiedPE && readsNumRatio>5))
							//else if((maxOccIndexPE!=-1 && maxOccIndexPE!=maxOccIndexSE) && (secondOccPE/maxOccPE>SECOND_FIRST_OCC_FAILED_RATIO))
							{
#if (DEBUG_OUTPUT==YES)
								printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
								naviSuccessFlag = NAVI_FAILED;
								kmers[0] = kmers[1] = NULL;
							}
							else if((maxOccIndexPE!=-1 && maxOccIndexPE!=maxOccIndexSE) && secondOccPE>minKmerOccPE && secondOccSE>minKmerOccSE)
							{
#if (DEBUG_OUTPUT==YES)
								printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
								naviSuccessFlag = NAVI_FAILED;
								kmers[0] = kmers[1] = NULL;
							}
							else if(maxOccIndexPE!=-1 && maxOccPE<minKmerOccPE && maxOccSE<minKmerOccSE && secondOccPE/maxOccPE>0.9)
							{
#if (DEBUG_OUTPUT==YES)
								printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
								naviSuccessFlag = NAVI_FAILED;
								kmers[0] = kmers[1] = NULL;
							}
						}
					}
				}

				else if(secondOccPE>3*minKmerOccPE && secondOccSE>=3*minKmerOccSE && sumSecondOccPE>maxOccPE && sumSecondOccSE>maxOccSE)  // added 2012-11-21
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				//else if((maxOccPE>0 && maxOccPE<minKmerOccPE && secondOccPE>0 && maxOccPE-secondOccPE<2) && secondOccSE>minKmerOccSE)   // added 2012-11-14
				//else if(maxOccIndexPE==-1 && secondOccSE>0 && maxOccSE-secondOccSE<2)   // added 2012-11-16, deleted 2012-11-23
				//else if(maxOccIndexPE==-1 && maxOccSE<minKmerOccSE && maxOccSE-secondOccSE<2)   // added 2012-11-23, deleted 2012-11-24
				else if(maxOccIndexPE==-1 && maxOccSE-secondOccSE<2)   // added 2012-11-24
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				else if(maxOccIndexPE==-1 && secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>0.1)  // added 2012-11-24
				{
					// compute the maximal gap size in contig tail region
					if(computeGapSizeInContig(&tmp_gapSize,contigArr, contigNodesNum, assemblyRound)==FAILED)
					{
						printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
						return FAILED;
					}

					//if(tmp_gapSize>10)  // added 2012-11-24, deleted 2012-11-26
					if(secondOccSE/maxOccSE>0.4 && tmp_gapSize>=15)  // added 2012-11-26
					//if((averKmerOcc<20 && tmp_gapSize>=13) || (averKmerOcc>=20 && secondOccSE/maxOccSE>0.5 && tmp_gapSize>=15))  // added 2012-11-26, deleted 2012-11-26
					{
#if (DEBUG_OUTPUT==YES)
						printf("==line=%d, contigNodesNum=%d, tmp_gapSize=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, tmp_gapSize, readsNumRatio);
#endif
						naviSuccessFlag = NAVI_FAILED;
						kmers[0] = kmers[1] = NULL;
					}
				}

				else if((maxOccPE>0 && maxOccPE<minKmerOccPE && secondOccPE>0 && maxOccPE-secondOccPE<2) && secondOccSE/maxOccSE>0.35)   // added 2012-11-16
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				//else if(maxOccIndexPE==-1 && secondOccSE>minKmerOccSE && sumSecondOccSE/maxOccSE>0.6)   // added 2012-11-21
				//else if(maxOccIndexPE==-1 && secondOccSE>minKmerOccSE  && sumSecondOccSE/maxOccSE>0.6 && readsNumRatio<0.7)   // added 2012-11-21, deleted 2013-01-20
				else if(secondOccSE>=minKmerOccSE && sumSecondOccSE/maxOccSE>0.7)   // added 2013-01-20
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				else if(maxOccIndexPE!=-1 && maxOccIndexPE!=maxOccIndexSE && secondOccSE>=minKmerOccSE)   // added 2013-01-20
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				else if(secondOccPE>3*minKmerOccSE && sumSecondOccPE/maxOccPE>0.15 && secondOccSE>2*minKmerOccSE && sumSecondOccSE/maxOccSE>0.15)   // added 2013-01-21
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				else if(maxOccIndexPE==-1 && secondOccSE<1.5*minKmerOccSE && secondOccSE/maxOccSE>0.5 && (successContigIndex>0 &&  contigNodesNum-successContigIndex >= 15))   // added 2012-11-21
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				//else if(maxOccIndexPE==-1 && readsNumRatio<0.3)
				else if(maxOccIndexPE==-1 && readsNumRatio<minReadsNumRatioThres)
				{
					// compute the maximal gap size in contig tail region
					if(computeGapSizeInContig(&tmp_gapSize, contigArr, contigNodesNum, assemblyRound)==FAILED)
					{
						printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(tmp_gapSize>0.7*readLen)
					//if(tmp_gapSize>0.6*readLen)   // added 2012-11-16
					{
#if (DEBUG_OUTPUT==YES)
						printf("==line=%d, contigNodesNum=%d, tmp_gapSize=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, tmp_gapSize, readsNumRatio);
#endif
						naviSuccessFlag = NAVI_FAILED;
						kmers[0] = kmers[1] = NULL;
					}
					else if(maxOccPE==0 && tmp_gapSize>0.3*readLen)   // added 2012-11-14
					{
#if (DEBUG_OUTPUT==YES)
						printf("==line=%d, contigNodesNum=%d, tmp_gapSize=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, tmp_gapSize, readsNumRatio);
#endif
						naviSuccessFlag = NAVI_FAILED;
						kmers[0] = kmers[1] = NULL;
					}
					else if(readsNumRatio<0.3*minReadsNumRatioThres && tmp_gapSize>0.4*readLen)   // added 2012-11-14
					{
#if (DEBUG_OUTPUT==YES)
						printf("==line=%d, contigNodesNum=%d, tmp_gapSize=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, tmp_gapSize, readsNumRatio);
#endif
						naviSuccessFlag = NAVI_FAILED;
						kmers[0] = kmers[1] = NULL;
					}
				}
				else if(maxOccIndexPE!=-1 && maxOccIndexPE!=maxOccIndexSE && secondOccPE>minKmerOccPE && secondOccPE/maxOccPE>0.5 && secondOccSE/maxOccSE>0.5)
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}
//				//else if(readsNumRatio>1.9 && secondOccPE>5*minKmerOccPE && secondOccSE>5*minKmerOccSE && secondOccPE/maxOccPE>0.5 && secondOccSE/maxOccSE>0.5)  // added 2012-11-09
//				else if(maxOccIndexPE!=-1 && readsNumRatio>1.9 && secondOccPE>10 && secondOccSE>15 && secondOccPE/maxOccPE>0.55 && secondOccSE/maxOccSE>0.5)  // added 2012-11-09
//				{
//#if (DEBUG_OUTPUT==YES)
//					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
//#endif
//					kmers[0] = kmers[1] = NULL;
//				}
//				else if(maxOccIndexPE!=-1 && readsNumRatio>2.5 && secondOccPE>25 && secondOccSE>10 && secondOccPE/maxOccPE>0.6 && secondOccSE/maxOccSE>0.6)  // added 2012-11-09
//				{
//#if (DEBUG_OUTPUT==YES)
//					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
//#endif
//					kmers[0] = kmers[1] = NULL;
//				}
				//else if(readsNumRatio>1.5 && secondOccPE>6*minKmerOccPE && secondOccSE>5*minKmerOccSE && secondOccPE/maxOccPE>SECOND_FIRST_OCC_RATIO)  // added 2012-11-09
				else if(maxOccIndexPE!=-1 && readsNumRatio>1.5 && secondOccPE>6*minKmerOccPE && secondOccSE>5*minKmerOccSE && secondOccPE/maxOccPE>0.9)  // added 2012-11-09
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}
				else if(maxOccIndexPE!=-1 && secondOccSE>minKmerOccSE && secondOccPE/maxOccPE>SECOND_FIRST_OCC_RATIO)  // added 2012-11-10
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				else if(secondOccPE>4*minKmerOccPE && secondOccPE/maxOccPE>0.9 && maxOccIndexPE!=maxOccIndexSE)  // added 2012-11-23
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				else if(maxOccIndexPE!=maxOccIndexSE && secondOccPE>2*minKmerOccPE && secondOccSE>2*minKmerOccPE && secondOccSE/maxOccSE>0.6) // added 2012-11-26
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				//else if((secondOccPE>3*minKmerOccPE && secondOccPE/maxOccPE>0.6) && (secondOccSE>minKmerOccSE && secondOccSE/maxOccSE>0.2))  // added 2012-11-21, deleted 2012-11-25
				else if((secondOccPE>3*minKmerOccPE && secondOccPE/maxOccPE>0.7) && (secondOccSE>minKmerOccSE && secondOccSE/maxOccSE>0.2))  // added 2012-11-25
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				else if((secondOccPE>3*minKmerOccPE && secondOccPE/maxOccPE>0.6) && (secondOccSE>minKmerOccSE && secondOccSE/maxOccSE>0.4))  // added 2012-11-26
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				else if((secondOccPE>2*minKmerOccPE && secondOccPE/maxOccPE>=0.5) && (secondOccSE>minKmerOccSE && secondOccSE/maxOccSE>=0.5))  // added 2012-11-22
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				else if((secondOccPE>4*minKmerOccPE && secondOccPE/maxOccPE>0.4) && (secondOccSE>4*minKmerOccSE && secondOccSE/maxOccSE>=0.5) && readsNumRatio>3)  // added 2012-11-22
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				//else if(readsNumRatio>7 && (secondOccPE>4*minKmerOccPE && secondOccPE/maxOccPE>0.15) && secondOccSE>minKmerOccSE)   // added 2012-11-22, deleted 2012-11-23
				//else if(readsNumRatio>6 && (secondOccPE>4*minKmerOccPE && secondOccPE/maxOccPE>0.1) && secondOccSE>minKmerOccSE)   // added 2012-11-23, deleted 2013-01-20
				else if(readsNumRatio>6 && (secondOccPE>3*minKmerOccPE && secondOccPE/maxOccPE>0.1) && (secondOccSE>3*minKmerOccSE && secondOccSE/maxOccSE>0.2))   // added 2013-01-20
				//else if(readsNumRatio>6 && (secondOccPE>4*minKmerOccPE && secondOccPE/maxOccPE>0.1) && (secondOccSE>minKmerOccSE && secondOccSE/maxOccSE>0.2))   // added 2012-11-25, deleted 2012-11-25
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				//else if(maxOccIndexPE!=-1 && maxOccPE<minKmerOccPE && secondOccPE/maxOccPE>0.9 && ((successContigIndex>0 &&  contigNodesNum-successContigIndex >= 10)))  // added 2012-11-10 ===========================
				//else if(maxOccPE>minKmerOccPE && secondOccPE/maxOccPE>0.9 && secondOccPE/maxOccPE>0.3 && ((successContigIndex>0 &&  contigNodesNum-successContigIndex >= 10)))  // added 2012-11-15, deleted 2012-11-23
				else if(maxOccPE>minKmerOccPE && secondOccPE/maxOccPE>0.9 && ((successContigIndex>0 &&  contigNodesNum-successContigIndex >= 10)))  // added 2012-11-23
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}
				else if(maxOccIndexPE!=-1 && secondOccPE>5*minKmerOccPE && secondOccPE/maxOccPE>SECOND_FIRST_OCC_FAILED_RATIO && ((successContigIndex>0 &&  contigNodesNum-successContigIndex >= 20)))  // added 2012-11-10 ===========================
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}

				navigationID = 0;
				navigationNumSE ++;
			}

			//else if(secondOccPE<minKmerOccPE && secondOccPE/maxOccPE>SECOND_FIRST_OCC_RATIO && (successContigIndex>0 &&  contigNodesNum-successContigIndex >= 20))  // added 2012-11-09
			else if(secondOccPE<minKmerOccPE && secondOccPE/maxOccPE>SECOND_FIRST_OCC_RATIO && (successContigIndex>0 &&  contigNodesNum-successContigIndex >= 30))  // added 2012-11-16
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}
			else if(readsNumRatio<0.4 && secondOccPE>3*minKmerOccPE && sumSecondOccPE/maxOccPE>0.5) // added 2012-11-10
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}
			//else if(secondOccPE>0 && secondOccPE/maxOccPE>0.9) // added 2012-11-14
			//else if(secondOccPE>minKmerOccPE && secondOccPE/maxOccPE>0.9) // added 2012-11-16, deleted 2012-11-23
			else if(maxOccSE<minKmerOccSE && secondOccPE>minKmerOccPE && secondOccPE/maxOccPE>0.9) // added 2012-11-23
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}
			else if(readsNumRatio<minReadsNumRatioThres && maxOccSE>maxOccNumFaiedPE && (successContigIndex>0 &&  contigNodesNum-successContigIndex >= 0.4*readLen)) // added 2012-11-14
			{
#if (DEBUG_OUTPUT==YES)
				printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
				naviSuccessFlag = NAVI_FAILED;
				kmers[0] = kmers[1] = NULL;
			}

			//else if(maxOccSE>8*minKmerOccSE && (maxOccPE<3*minKmerOccPE && secondOccPE<minKmerOccPE && secondOccPE/maxOccPE>0.3))  // added 2012-11-15
			//else if(maxOccSE<4*minKmerOccSE && (maxOccPE<minKmerOccPE && secondOccPE<0.3*minKmerOccPE && secondOccPE/maxOccPE>0.3))  // added 2012-11-16, deleted 2012-11-25
			//else if(averKmerOcc>80 && maxOccSE<4*minKmerOccSE && (maxOccPE<minKmerOccPE && secondOccPE<0.3*minKmerOccPE && secondOccPE/maxOccPE>0.3))  // added 2012-11-25, deleted 2012-11-26
			else if(maxOccSE<2*minKmerOccSE && (maxOccPE>0 && maxOccPE<minKmerOccPE && secondOccPE<0.3*minKmerOccPE && secondOccPE/maxOccPE>0.3))  // added 2012-11-26
			{
				if((maxOccIndexPE!=maxOccIndexSE && maxOccPE-secondOccPE>1) || secondOccSE>0)
				{
#if (DEBUG_OUTPUT==YES)
					printf("####line=%d, contigNodesNum=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}
			}

			//else if(readsNumRatio<0.2*minReadsNumRatioThres)
			//else if(readsNumRatio<0.2)
			else if(readsNumRatio<minReadsNumRatioThres)
			{
				// compute the maximal gap size in contig tail region
				if(computeGapSizeInContig(&tmp_gapSize, contigArr, contigNodesNum, assemblyRound)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute the gap size, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(tmp_gapSize>0.7*readLen)
				//if(tmp_gapSize>0.6*readLen)   // added 2012-11-16
				{
#if (DEBUG_OUTPUT==YES)
					printf("==line=%d, contigNodesNum=%d, tmp_gapSize=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, tmp_gapSize, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}
				else if(maxOccPE==0 && tmp_gapSize>0.3*readLen)   // added 2012-11-14
				{
#if (DEBUG_OUTPUT==YES)
					printf("==line=%d, contigNodesNum=%d, tmp_gapSize=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, tmp_gapSize, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}
				else if(readsNumRatio<0.3*minReadsNumRatioThres && tmp_gapSize>0.4*readLen)   // added 2012-11-14
				{
#if (DEBUG_OUTPUT==YES)
					printf("==line=%d, contigNodesNum=%d, tmp_gapSize=%d, readsNumRatio=%.2f\n", __LINE__, contigNodesNum, tmp_gapSize, readsNumRatio);
#endif
					naviSuccessFlag = NAVI_FAILED;
					kmers[0] = kmers[1] = NULL;
				}
			}
		}
#endif
	}else
	{
		navigationID = 1;
		navigationNumSE = 0;
	}

	return SUCCESSFUL;
}

/**
 * Get the next k-mers for navigation combining the k-mer hash table and the decision table using paired ends.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNextKmerByPE(int contigNodesNum)
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

	maxOccIndexPE = -1;
	maxOccPE = 0;
	secondOccIndexPE = -1;
	secondOccPE = 0;
	occsNumIndexPE[0] = -1;
	occsNumIndexPE[1] = -1;


//	validKmerNum = 0;
	for(i=0; i<4; i++)
	{
		//occsNumPE[i] = 0;

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
//		if(tmp_kmers[i][0] || tmp_kmers[i][1])
//		{
//			validKmerNum ++;
//			base_index = i;
//		}

	}
/*
	if(validKmerNum==1)
	{
		//score[base_index] = computeKmerScoreByPE(tmp_kmers[base_index], occsNum+base_index, assemblingreads, numassemblingreads);
		if(computeKmerOccNumUnlockedByPE(occsNumPE+base_index, tmp_kmers[base_index], base_index, minSuccessiveAppearedBaseNum, NO)==FAILED)
		{
			printf("In %s(), cannot compute the occsNum of a kmer, error!\n", __func__);
			return FAILED;
		}

		//========================= condition 1 ==============================
		if(occsNumPE[base_index]>0)
		//if(occsNumPE[base_index]>=minKmerOccPE)
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

			maxOccIndexPE = base_index;
			maxOccPE = occsNumPE[maxOccIndexPE];
			secondOccIndexPE = -1;
			secondOccPE = 0;
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
	for(iterateID=0; iterateID<2; iterateID++)
	{
		if(iterateID==0)
		{
			successiveAppearBaseNum = maxSuccessiveAppearedBaseNum;
			ignoreEndBaseFlag = NO;
		}else
		{
			successiveAppearBaseNum = minSuccessiveAppearedBaseNum;
			//if(maxOccPE<minKmerOccPE)
			//	ignoreEndBaseFlag = YES;
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
						//if(tmp_kmers[i][0]->multiplicity+tmp_kmers[i][1]->multiplicity<minKmerOccPE)
						if(tmp_kmers[i][0]->multiplicity+tmp_kmers[i][1]->multiplicity<=0)
						{
							continue;
						}
					}else if(tmp_kmers[i][0]!=NULL)
					{
						//if(tmp_kmers[i][0]->multiplicity<minKmerOccPE)
						if(tmp_kmers[i][0]->multiplicity<=0)
						{
							continue;
						}
					}else if(tmp_kmers[i][1]!=NULL)
					{
						//if(tmp_kmers[i][1]->multiplicity<minKmerOccPE)
						if(tmp_kmers[i][1]->multiplicity<=0)
						{
							continue;
						}
					}
				}
		*/
				//########################## end ########################//


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
			}

		//	if(validKmerNum>1)
		//	{
		//		printf("validKmerNum=%d, contigNodexNum=%d\n", validKmerNum, contigNodesNum);
		//	}

		/*
			//========================= condition 4 ==============================
			//if(validKmerNum>1 && maxIndex1!=maxOccIndex)
			//if(validKmerNum>1 && occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM)
			//if(validKmerNum>0 && occsNumPE[maxIndex1]<MIN_CONNECT_KMER_NUM) //--best result
			if(validKmerNum>0 && occsNumPE[maxIndex1]<minKmerOccPE)
			//if(occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM || (validKmerNum>1 && maxIndex1!=maxOccIndex))
			//if(occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM || occsNum[secondOccIndex]>=MIN_CONNECT_KMER_NUM)
			{
				validKmerNum = 0;
			}


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
			if(validKmerNum>1 && (occsNumPE[secondIndex1] > maxSecondOcc || (occsNumPE[maxIndex1] > maxFirstOcc && occsNumPE[secondIndex1]>=maxSecondOcc)
					|| ((float)occsNumPE[secondIndex1]/occsNumPE[maxIndex1]>SECOND_FIRST_OCC_RATIO && occsNumPE[secondIndex1]>=maxSecondOcc))) //-- best result
			{
				validKmerNum = 0;
			}


			//========================= condition 6 ==============================
			if(validKmerNum>1 && ((float)occsNumPE[maxIndex1]/itemNumDecisionTable < VALID_OCC_RATIO && occsNumPE[secondIndex1] >= maxSecondOcc))
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
			//if(validKmerNum>1 && ((secondOccPE/maxOccPE>SECOND_FIRST_OCC_RATIO/2 && secondOccPE>minKmerOccPE) || (secondOccPE>maxSecondOcc)))
			//if(maxOccPE<minKmerOccPE || (validKmerNum>1 && ((secondOccPE/maxOccPE>SECOND_FIRST_OCC_RATIO/2 && secondOccPE>minKmerOccPE) || (secondOccPE>maxSecondOcc))))
			if(maxOccPE<minKmerOccPE || (validKmerNum>1 && (secondOccPE/maxOccPE>SECOND_FIRST_OCC_RATIO/2 || (secondOccPE>maxSecondOcc))))
			{
				if(secondOccPE/maxOccPE>SECOND_FIRST_OCC_FAILED_RATIO)
					validKmerNum = 0;
				else if(secondOccPE>0 && maxOccPE-secondOccPE<2)
					validKmerNum = 0;
				//else if(maxOccPE<minKmerOccPE && itemNumDecisionTable>3*minKmerOccPE) // added 2013-01-19, deleted 2013-01-20
				else if(maxOccPE==1 && itemNumDecisionTable>4*minKmerOccPE) // added 2013-01-19
					validKmerNum = 0;
			}
#endif
	/*
			//========================= condition 9 ==============================
			//if(kmer_len > longKmerSize - longKmerStepSize && occsNumSE[maxIndex1] <= minLongKmerOcc)
			if(kmer_len > longKmerSize - longKmerStepSize && maxOccPE < minLongKmerOcc)
			{
				validKmerNum = 0;
			}
	*/

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

				//naviSuccessFlag = NAVI_FAILED;
				//kmers[0] = kmers[1] = NULL;
				//return SUCCESSFUL;

//				kmer_len -= longKmerStepSize;
//				if(kmer_len<kmerSize && kmer_len>kmerSize-longKmerStepSize)
//					kmer_len = kmerSize;
//				continue;
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
	//***********************************************************************************

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

	int32_t i, properRow, posNum, basePos, baseInt_read;
	int32_t entryIDReadseq, entryPosReadseq, baseNumEntryReadseq;
	int32_t rid, maxItemNumPerReadBlock, blockRowRead, itemRowBlockRead;
	int32_t seqLen, errorRegLenEnd3;


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

#if 1
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
//				returnCode = validReadPair(rid_pos_table[i].rid);
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
				if(existReadInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
				{
					returnCode = validReadPair(rid_pos_table[i].rid);
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

#if 1
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
//				returnCode = validReadPair(rid_pos_table[i].rid);
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
				if(existReadInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
				{
					returnCode = validReadPair(rid_pos_table[i].rid);
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

	*occNum = 0;

	for(i=0; i<itemNumDecisionTable; i++)
	{
		if((decisionTable[i].matedFlag==YES && decisionTable[i].reserved==0)
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
//					returnCode = validReadPair(rid_pos_table[i].rid);
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
					if(existReadInDecisionTable(rid_pos_table[i].rid, rid_pos_table[i].pos, ORIENTATION_MINUS, decisionTable, dtRowHashtable)==NO)
					{
						returnCode = validReadPair(rid_pos_table[i].rid);
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
short validReadPair(uint64_t readID)
{
	uint64_t readID_paired;
	int readOrient_paired, contigPos_paired;
	int fragmentSize, validReadOrient, validReadOrient_paired;
	PERead_t *pPERead;

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
		return YES;
	else
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
 * Set the navigation occurrence queue to empty.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short setEmptyNaviOccQueue(double *naviOccQueuePara, int *itemNumNaviOccQueuePara, int *frontRowNaviOccQueuePara, int *rearRowNaviOccQueuePara)
{
	int i;

	for(i=0; i<*itemNumNaviOccQueuePara; i++) naviOccQueuePara[i] = 0;
	*itemNumNaviOccQueuePara = 0;
	*frontRowNaviOccQueuePara = 0;
	*rearRowNaviOccQueuePara = 0;

	return SUCCESSFUL;
}

/**
 * Update the navigation occurrence queue.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateNaviOccQueue(double *naviOccQueuePara, int *itemNumNaviOccQueuePara, int *frontRowNaviOccQueuePara, int *rearRowNaviOccQueuePara, double maxOccNum)
{
	if(*itemNumNaviOccQueuePara>=maxItemNumNaviOccQueue)
	{ // full queue, then add the first item
		(*rearRowNaviOccQueuePara) = ((*rearRowNaviOccQueuePara) + 1) % maxItemNumNaviOccQueue;
		(*frontRowNaviOccQueuePara) = ((*frontRowNaviOccQueuePara) + 1) % maxItemNumNaviOccQueue;
		naviOccQueuePara[*rearRowNaviOccQueuePara] = maxOccNum;
	}else if(*itemNumNaviOccQueuePara==0)
	{ // empty queue, then add the first item
		*frontRowNaviOccQueuePara = *rearRowNaviOccQueuePara = 0;
		naviOccQueuePara[*rearRowNaviOccQueuePara] = maxOccNum;
		*itemNumNaviOccQueuePara = 1;
	}else if(*itemNumNaviOccQueuePara<maxItemNumNaviOccQueue)
	{ // non-full queue, add new item
		(*rearRowNaviOccQueuePara) ++;
		naviOccQueuePara[*rearRowNaviOccQueuePara] = maxOccNum;
		(*itemNumNaviOccQueuePara) ++;
	}else
	{ // Exception errors
		printf("line=%d, In %s(), itemNumNaviOccQueue=%d, error!\n", __LINE__, __func__, *itemNumNaviOccQueuePara);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Compute the average number in navigation occurrence queue.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short calcAverOccNaviOccQueue(double *averOccNum, double *naviOccQueuePara, int itemNumNaviOccQueuePara)
{
	int i;

	*averOccNum = 0;
	for(i=0; i<itemNumNaviOccQueuePara; i++)
		*averOccNum += naviOccQueuePara[i];

	if(itemNumNaviOccQueuePara>0)
		*averOccNum /= itemNumNaviOccQueuePara;
	else
		*averOccNum = 0;

	return SUCCESSFUL;
}

/**
 * Get the length of low occurrence region in navigation occurrence queue.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getLowOccLenNaviOccQueue(int *lowLen, double *naviOccQueuePara, int itemNumNaviOccQueuePara, int frontRowNaviOccQueuePara)
{
	int i, maxLen, len;

	*lowLen = 0;
	for(i=0; i<itemNumNaviOccQueuePara; i++)
	{
		if(naviOccQueuePara[i]<=lowOccThresNaviOccQueue)
		{
			(*lowLen) ++;
		}
	}

/*
	len = 0;
	maxLen = 0;
	j = frontRowNaviOccQueuePara;
	for(i=0; i<itemNumNaviOccQueuePara; i++)
	{
		if(naviOccQueuePara[j]<=lowOccThresNaviOccQueue)
		{
			len ++;
		}else
		{
			if(len>maxLen)
				maxLen = len;
			len = 0;
		}
		j = (j+1) % maxItemNumNaviOccQueue;
	}

	*lowLen = maxLen;
*/

	return SUCCESSFUL;
}


/**
 * Compute the maximal gap size in contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeGapSizeInContig(int *gapSize, contigtype *contigArr, int64_t contigNodesNum, int assemblyRound)
{
	int64_t i, tmp_GapLen, maxLen, startIndex, endIndex;

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{
		startIndex = contigNodesNum - 1.5*readLen + 1;
		if(startIndex<=0)
			startIndex = 1;
		endIndex = contigNodesNum;
	}else
	{
		endIndex = contigNodesNum - 0.85 * readLen;
		if(endIndex<=0)
			endIndex = 1;
		startIndex = endIndex - 1.5*readLen + 1;
		if(startIndex<=0)
			startIndex = 1;
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
