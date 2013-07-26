/*
 * scafUpdate.c
 *
 *  Created on: Aug 2, 2011
 *      Author: zhuxiao
 */

#include "inc/scafStdinc.h"
#include "inc/scafExtvab.h"


/**
 * Update the scafAssembling reads in decision table in local assembly.
 *  Note:
 *  	For the successful reads, their assembly status will be updated.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateScafAssemblingReadsInScaf(int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound)
{
	int i, j, startIndex, exectIndex, posNum;
	scafRidpos *ridpostable;
	scafRidpos *rid_pos;
	scafAssemblingRead *this_assemblingRead = scafAssemblingReadArr;
	int returnCode, matedFlag;

	// traverse decision table, update its reads information
	for(i=0; i<scafAssemblingReadsNum; i++)
	{
		for(j=0; j<2; j++)
		{ //j==0 for scafKmers[0], j==1 for scafKmers[1]
			if(scafKmers[j] && ((j==0 && this_assemblingRead->orientation==ORIENTATION_PLUS) || (j==1 && this_assemblingRead->orientation==ORIENTATION_MINUS)))
			{ // plus scafKmer exists, and the read has plus orientation; or minus scafKmer exists, and the read has minus orientation
				// check whether the read exist in scafRispos table
				ridpostable = scafKmers[j]->ppos;
				posNum = scafKmers[j]->arraysize;
				startIndex = findStartScafRidposIndexInScaf(this_assemblingRead->rid, ridpostable, posNum); //binary search the start entry index of the read in scafRidpos table
				if(startIndex>=0)
				{ // find the start entry index, then get the exact entry index
					if(this_assemblingRead->lastpos>0)
					{ // the read take part in local assembly last time
						if(j==0) // plus scafKmer
							exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->lastpos+1, startIndex, ridpostable, posNum);
						else // minus scafKmer
							exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->lastpos-1, startIndex, ridpostable, posNum);
					}else
					{ // the read does not take part in local assembly last time
						if(j==0) // plus scafKmer
							exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->firstpos+this_assemblingRead->kmerappeartimes+this_assemblingRead->kmerunappeartimes, startIndex, ridpostable, posNum);
						else // minus scafKmer
							exectIndex = getExectScafRidposIndexInScaf(this_assemblingRead->rid, this_assemblingRead->firstpos-this_assemblingRead->kmerappeartimes-this_assemblingRead->kmerunappeartimes, startIndex, ridpostable, posNum);
					}

					if(exectIndex>=0)
					{  // valid read position
						rid_pos = ridpostable + exectIndex;
						if(rid_pos->used==0)
						{ // unused read
							if(this_assemblingRead->lastpos!=0)
							{ // the read takes part in assembly last time, then it take part the assembly this time
								this_assemblingRead->kmerappeartimes++;
								this_assemblingRead->lastpos = rid_pos->pos;
								this_assemblingRead->latestpos = rid_pos->pos;
							}else
							{ // the read does not take part in assembly last time, then the read take part in assembly this time
								//if(assemblingreadtable[i].kmerunapperblocks>0)
								//{ // disappeared block, then the new disappeared block is open
									//assemblingreadtable[i].kmerappeartimes += assemblingreadtable[i].kmerunappeartimes + 1;
									//assemblingreadtable[i].kmerunappeartimes = 0;

									this_assemblingRead->kmerappeartimes++;
									this_assemblingRead->lastpos = rid_pos->pos;
									this_assemblingRead->latestpos = rid_pos->pos;
								//}
							}
							rid_pos->reserved = 1;
						}else
						{ // the used read
							this_assemblingRead->delsign = 1;
						}

					}else
					{ // invalid read position, then the read does not take part in local assembly
						if(this_assemblingRead->lastpos!=0)
						{ // new disappeared block is open, and the appear time is recomputed
							this_assemblingRead->kmerappeartimes += this_assemblingRead->kmerunappeartimes;
							this_assemblingRead->kmerunappeartimes = 1;
							this_assemblingRead->lastpos = 0;
							this_assemblingRead->kmerunappearblocks++;
						}else
						{ // the read disappears continuous two times
							this_assemblingRead->kmerunappeartimes++;
						}

					} // end if(exectIndex>=0)

				}else
				{ // does not find the read, then the read does not take part in the local assembly this time
					if(this_assemblingRead->lastpos!=0)
					{ // new disappeared block is open, and the appear time is recomputed
						this_assemblingRead->kmerappeartimes += this_assemblingRead->kmerunappeartimes;
						this_assemblingRead->kmerunappeartimes = 1;
						this_assemblingRead->lastpos = 0;
						this_assemblingRead->kmerunappearblocks++;
					}else
					{
						this_assemblingRead->kmerunappeartimes++;
					}
				}

			}else if(scafKmers[j]==NULL && ((j==0 && this_assemblingRead->orientation==ORIENTATION_PLUS) || (j==1 && this_assemblingRead->orientation==ORIENTATION_MINUS)))
			{ // plus scafKmer does not exist, but the read has plus orientation; or minus scafKmer does not exist, but the read has minus orientation
				if(this_assemblingRead->lastpos!=0)
				{ // new disappeared block is open, and the appear time is recomputed
					this_assemblingRead->kmerappeartimes += this_assemblingRead->kmerunappeartimes;
					this_assemblingRead->kmerunappeartimes = 1;
					this_assemblingRead->lastpos = 0;
					this_assemblingRead->kmerunappearblocks++;
				}else
				{
					this_assemblingRead->kmerunappeartimes++;
				}
			}
		}

		this_assemblingRead ++;
	}

	// add new reads to decision table
	int pos = 0;
	if(scafKmers[0])
	{
		ridpostable = scafKmers[0]->ppos;
		posNum = scafKmers[0]->arraysize;
		for(i=0; i<posNum; i++)
		{
			pos = ridpostable[i].pos;
			if(ridpostable[i].used==0 && ridpostable[i].reserved==0 && pos==1)
			{ // unused read, unmarked, and pos==1, then add to decision table

				//######################## Debug information #######################
				//if(ridpostable[i].rid==5323225)
				//{
				//	printf("line=%d, In %s(), rid=%d, pos=%d, used=%d, reserved=%d, assemblyStatus=%d\n", __LINE__, __func__, ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].used, ridpostable[i].reserved, assemblyStatus);
				//}
				//######################## Debug information #######################

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

				if(addAssemblingReadsInScaf(ridpostable[i].rid, pos, ORIENTATION_PLUS, matedFlag)==FAILED)
				{
					printf("line=%d, In %s(), cannot add read to decision table, error!\n", __LINE__, __func__);
					return FAILED;
				}

			}//end if
			ridpostable[i].reserved = 0;
		}//end for
	}

	if(scafKmers[1])
	{
		ridpostable = scafKmers[1]->ppos;
		posNum = scafKmers[1]->arraysize;
		for(i=0; i<posNum; i++)
		{
			pos = ridpostable[i].pos;
			//this_assemblingRead = &assemblingreadtable[assemblingtable_size];
			if(ridpostable[i].used==0 && ridpostable[i].reserved==0 && pos>=readLen-kmerSize+1-errorRegLenEnd3)
			{ // unused read, unmarked, and pos at 3' end, then add to decision table

				//######################## Debug information #######################
				//if(ridpostable[i].rid==5323225)
				//{
				//	printf("line=%d, In %s(), rid=%d, pos=%d, delsign=%d, reserved=%d, assemblyStatus=%d\n", __LINE__, __func__, ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved, assemblyStatus);
				//}
				//######################## Debug information #######################

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

				if(addAssemblingReadsInScaf(ridpostable[i].rid, pos, ORIENTATION_MINUS, matedFlag)==FAILED)
				{
					printf("line=%d, In %s(), cannot add read to decision table, error!\n", __LINE__, __func__);
					return FAILED;
				}

			}//end if
			ridpostable[i].reserved = 0;
		}//end for
	}

	return SUCCESSFUL;
}

/**
 * Update the status of reads in decision table in local assembly.
 *  Three statuses:
 *  	(1) ASSEMBLING_STATUS for the reads being taking part in the assembly;
 *  	(2) SUCCESSFUL_STATUS for the reads successfully assembled;
 *  	(3) FAILED_STATUS for the reads failed in the assembly.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateAssemblingreadsStatusInScaf()
{
	scafAssemblingRead *this_assemblingRead = scafAssemblingReadArr;
	int i;

	//if(this_assemblingRead->rid==3674534)
	//{
	//	printf("rid=%u\n", this_assemblingRead->rid);
	//}

	for(i=0; i<scafAssemblingReadsNum; i++)
	{
		if(this_assemblingRead->orientation==ORIENTATION_PLUS)
		{ // plus orientation read
			if(this_assemblingRead->delsign==1)
			{ // the read has been used, then assembly finished and failed
				this_assemblingRead->status = FAILED_STATUS;

			}else if(this_assemblingRead->kmerunappearblocks>1)
			{ // the disappeared blocks number more than 1, then assembly finished and failed
				this_assemblingRead->status = FAILED_STATUS;

			}else if(this_assemblingRead->kmerunappearblocks==1)
			{ // only one disappeared block
				if(this_assemblingRead->lastpos==0 && this_assemblingRead->latestpos>=readLen-kmerSize+1-errorRegLenEnd3)
				{ // 2 bases near to the 3' end disappeared, then assembly finished and succeeded
					this_assemblingRead->status = SUCCESSFUL_STATUS;
					updateSameReadStatusToFailedInScaf(this_assemblingRead->rid, i);

				}else if(this_assemblingRead->kmerappeartimes+this_assemblingRead->kmerunappeartimes==readLen-kmerSize+1)
				{ // get the maximal assembly times, then assembly finished
					if(this_assemblingRead->lastpos==readLen-kmerSize+1)
					{ // the read appeared this time and it is the last scafKmer of the read, then succeeded
						this_assemblingRead->status = SUCCESSFUL_STATUS;
						updateSameReadStatusToFailedInScaf(this_assemblingRead->rid, i);
					}
					else
					{ // otherwise, failed
						this_assemblingRead->status = FAILED_STATUS;
					}
				}else
				{ // get the maximal assembly times of a read
					if(this_assemblingRead->lastpos>0 && this_assemblingRead->kmerunappeartimes!=kmerSize)
					{ // the read appeared this time and disappeared times is not KMER_SIZE, then finished and failed
						this_assemblingRead->status = FAILED_STATUS;

					}else if(this_assemblingRead->lastpos==0 && this_assemblingRead->kmerunappeartimes>kmerSize)
					{ // the read disappeared this time and disappeared times is more than KMER_SIZE, then finished and failed
						this_assemblingRead->status = FAILED_STATUS;
					}
				}
			}else
			{ // no disappeared blocks
				if(this_assemblingRead->kmerappeartimes==readLen-kmerSize+1)
				{ // the maximal assembly times of a read, then finished and succeeded
					this_assemblingRead->status = SUCCESSFUL_STATUS;
					updateSameReadStatusToFailedInScaf(this_assemblingRead->rid, i);
				}
			}
		}else
		{ // minus orientation read
			if(this_assemblingRead->delsign==1)
			{ // the read had been used, then its assembly finished and failed
				this_assemblingRead->status = FAILED_STATUS;

			}else if(this_assemblingRead->kmerunappearblocks>1)
			{ // more than one disappeared blocks, then assembly finished and failed
				this_assemblingRead->status = FAILED_STATUS;

			}else if(this_assemblingRead->kmerunappearblocks==1)
			{ // only one disappeared block
				if(this_assemblingRead->firstpos!=readLen-kmerSize+1)
				{ // the first scafKmer position was not the last one in the read, then assembly finised and failed
					this_assemblingRead->status = FAILED_STATUS;

				}else
				{ // the first scafKmer position is the last one in the read
					if(this_assemblingRead->kmerappeartimes+this_assemblingRead->kmerunappeartimes==readLen-kmerSize+1)
					{ // get the maximal assembly times, then finished
						if(this_assemblingRead->lastpos==1)
						{ // the read appeared this time, and its position is the first one in the read, then successful
							this_assemblingRead->status = SUCCESSFUL_STATUS;
							updateSameReadStatusToFailedInScaf(this_assemblingRead->rid, i);
						}
						else
						{ // otherwise, failed
							this_assemblingRead->status = FAILED_STATUS;
						}
					}else
					{ // had not get the maximal assembly times of the read
						if(this_assemblingRead->lastpos>0 && this_assemblingRead->kmerunappeartimes!=kmerSize)
						{ // the read appeared this time, and disappeared times is not KMER_SIZE, then finished and failed
							this_assemblingRead->status = FAILED_STATUS;

						}else if(this_assemblingRead->lastpos==0 && this_assemblingRead->kmerunappeartimes>kmerSize)
						{ // the read disappeared this time, and disappeared times was more than KMER_SIZE, then finished and failed
							this_assemblingRead->status = FAILED_STATUS;
						}
					}
				}
			}else
			{  // no disappeared block
				if(this_assemblingRead->firstpos!=readLen-kmerSize+1 && this_assemblingRead->lastpos==1)
				{ // the firstpos was not the last position in a read, and the position this time is the first one in a read, then finished and succeeded
					this_assemblingRead->status = SUCCESSFUL_STATUS;
					updateSameReadStatusToFailedInScaf(this_assemblingRead->rid, i);

				}else if(this_assemblingRead->firstpos==readLen-kmerSize+1 && this_assemblingRead->kmerappeartimes==readLen-kmerSize+1)
				{ // the firstpos was the last position in a read, and all the scafKmers were assembled onto contig, then finished and succeeded
					this_assemblingRead->status = SUCCESSFUL_STATUS;
					updateSameReadStatusToFailedInScaf(this_assemblingRead->rid, i);
				}
			}
		}

		this_assemblingRead ++;
	}

	return SUCCESSFUL;
}

/**
 * Change status of the read with the same read ID rid besides the given entry index to failed. *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateSameReadStatusToFailedInScaf(uint64_t rid, int assemblingreadIndex)
{
	scafAssemblingRead *this_assemblingRead = scafAssemblingReadArr;
	int i;
	for(i=0; i<scafAssemblingReadsNum; i++)
	{
		if(this_assemblingRead->rid==rid && i!=assemblingreadIndex)
		{
			this_assemblingRead->status = FAILED_STATUS;
		}
		this_assemblingRead ++;
	}
	return SUCCESSFUL;
}

/**
 * Update finished reads in decision table in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateFinisedScafReadsInScaf()
{
	// get the successful reads from decision table
	if(getScafSuccessReadsInScaf()==FAILED)
	{
		printf("line=%d, In %s(), cannot get the successful reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// remove the finished reads from decision table
	if(removeFinisedReadsInScaf()==FAILED)
	{
		printf("line=%d, In %s(), cannot remove the finished reads from decision table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get successful reads in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getScafSuccessReadsInScaf()
{
	int i;

	scafSuccessReadsNum = 0;
	for(i=0; i<scafAssemblingReadsNum; i++)
	{
		if(scafAssemblingReadArr[ i ].status==SUCCESSFUL_STATUS)
		{  //该read成功结束
			scafSuccessReadArr[ scafSuccessReadsNum ].rid = scafAssemblingReadArr[i].rid;
			scafSuccessReadArr[ scafSuccessReadsNum ].startmatchpos = scafAssemblingReadArr[i].firstpos;
			if(scafAssemblingReadArr[i].orientation==ORIENTATION_PLUS)
			{ //正向reads的kmer匹配数量, 为出现的两端的kmer相减+1
				scafSuccessReadArr[ scafSuccessReadsNum ].matchnum = scafAssemblingReadArr[i].latestpos - scafAssemblingReadArr[i].firstpos + 1;
			}else
			{
				scafSuccessReadArr[ scafSuccessReadsNum ].matchnum = scafAssemblingReadArr[i].firstpos - scafAssemblingReadArr[i].latestpos + 1;
			}
			scafSuccessReadArr[ scafSuccessReadsNum ].orientation = scafAssemblingReadArr[i].orientation;

			scafSuccessReadsNum ++;

			// check the allowed maximal number of reads
			if(scafSuccessReadsNum==maxScafSuccessReadsNumInArr)
			{ // reallocate the memory for scafSuccessReadArr
				if(reallocateScafSuccessReadsArrInScaf()==FAILED)
				{
					printf("line=%d, In %s(), cannot reallocate memory for successful reads, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}
	}

	return SUCCESSFUL;
}


/**
 * Reallocate the memory for the successful reads in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short reallocateScafSuccessReadsArrInScaf()
{
	scafSuccessRead *tmp_scafSuccessReadArray;

	tmp_scafSuccessReadArray = (scafSuccessRead *) malloc (2*maxScafSuccessReadsNumInArr * sizeof(scafSuccessRead));
	if(tmp_scafSuccessReadArray==NULL)
	{
		printf("line=%d, In %s(), cannot reallocate memory for successful reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// copy memory
	if(memcpy(tmp_scafSuccessReadArray, scafSuccessReadArr, scafSuccessReadsNum * sizeof(scafSuccessRead))==NULL)
	{
		printf("line=%d, In %s(), cannot copy memory for successful reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	free(scafSuccessReadArr);
	scafSuccessReadArr = tmp_scafSuccessReadArray;
	maxScafSuccessReadsNumInArr *= 2;

	return SUCCESSFUL;
}

/**
 * Remove the finished reads from decision table in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short removeFinisedReadsInScaf()
{
	int i, j;
	i = 0;
	j= scafAssemblingReadsNum - 1;
	while(i <= j)
	{
		if(scafAssemblingReadArr[i].status !=  ASSEMBLING_STATUS && scafAssemblingReadArr[j].status == ASSEMBLING_STATUS)
		{ // [i] is finished, [j] is unfinished, then [i] is replaced by [j], and head move down, tail move up
			if(memcpy(scafAssemblingReadArr+i, scafAssemblingReadArr+j, sizeof(scafAssemblingRead))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory for scafAssembling reads, error!\n", __LINE__, __func__);
				return FAILED;
			}
			scafAssemblingReadArr[i].status = ASSEMBLING_STATUS;
			scafAssemblingReadsNum --;
			j--;
			i++;
		}

		if(scafAssemblingReadArr[i].status == ASSEMBLING_STATUS)
			i++;

		if(scafAssemblingReadArr[j].status != ASSEMBLING_STATUS)
		{
			j--;
			scafAssemblingReadsNum --;
		}
	}

	return SUCCESSFUL;
}

