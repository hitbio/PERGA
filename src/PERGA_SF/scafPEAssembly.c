/*
 * scafPEAssembly.c
 *
 *  Created on: Feb 13, 2012
 *      Author: xiao
 */

#include "inc/scafStdinc.h"
#include "inc/scafExtvab.h"


/**
 * Get the next kmer by mixture of PE and SE.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNextKmerByMixInScaf(int contigNodesNum, int contigID, int contigLen, int contigOrient, int assemblyRound)
{
	//unsigned int tmp_kmerseqInt, tmp_kmerseqPE, tmp_kmerseqSE;
	scafKmer *tmp_kmers[2]/*, *tmp_kmersPE[2], *tmp_kmersSE[2]*/;
	int i, maxOccPE, maxOccSE, maxOccIndexPE, maxOccIndexSE;

	if(scafAssemblingReadsNum>MAX_DECISION_TABLE_SIZE_HTRES)
	{
		scafKmers[0] = scafKmers[1] = NULL;
		return SUCCESSFUL;
	}

	tmp_kmers[0] = scafKmers[0];
	tmp_kmers[1] = scafKmers[1];
	if(memcpy(tmpKmerSeqIntAssembly, kmerSeqIntAssembly, entriesPerKmer*sizeof(uint64_t))==NULL)
	{
		printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the score and occsNum of PE
	if(getNextKmerByPEInScaf(contigNodesNum, contigID, contigLen, contigOrient, assemblyRound)==FAILED)
	{
		printf("line=%d, In %s(), cannot get next kmer, error!\n", __LINE__, __func__);
		return FAILED;
	}
//	tmp_kmersPE[0] = kmers[0];
//	tmp_kmersPE[1] = kmers[1];
//	tmp_kmerseqPE = thiskmerseq;

	if(scafKmers[0]==NULL && scafKmers[1]==NULL)
	{
		// get the score and occsNum of SE
		scafKmers[0] = tmp_kmers[0];
		scafKmers[1] = tmp_kmers[1];
		if(memcpy(kmerSeqIntAssembly, tmpKmerSeqIntAssembly, entriesPerKmer*sizeof(uint64_t))==NULL)
		{
			printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(getNextKmerBySEInScaf(contigNodesNum, contigID, contigLen, contigOrient, assemblyRound)==FAILED)
		{
			printf("line=%d, In %s(), cannot get next kmer, error!\n", __LINE__, __func__);
			return FAILED;
		}
//		tmp_kmersSE[0] = kmers[0];
//		tmp_kmersSE[1] = kmers[1];
//		tmp_kmerseqSE = thiskmerseq;

		maxOccIndexSE = -1;
		maxOccSE = 0;
		for(i=0; i<4; i++)
		{
			if(occsNumSE[i]>maxOccSE)
			{
				maxOccSE = occsNumSE[i];
				maxOccIndexSE = i;
			}
		}

		maxOccIndexPE = -1;
		maxOccPE = 0;
		for(i=0; i<4; i++)
		{
			if(occsNumPE[i]>maxOccPE)
			{
				maxOccPE = occsNumPE[i];
				maxOccIndexPE = i;
			}
		}

		// check the scores and occsNums of PE and SE
		//if(maxOccIndexPE==-1 && maxOccSE>OCCS_NUM_SE_FAILED_PE_FACTOR*averKmerOcc)  //==================================
		//if(maxOccIndexPE==-1 && maxOccSE>maxFirstOcc)
		//if(maxOccIndexPE==-1 && maxOccSE>maxOccNumFaiedPE)
		if((maxOccIndexPE==-1 && maxOccSE>maxOccNumFaiedPE) && navigationID==0)
		//if((maxOccIndexPE==-1 && maxOccSE>maxOccNumFaiedPE) && (navigationNumSE>maxNavigationNumSE))
		{
			scafKmers[0] = scafKmers[1] = NULL;
		}

		navigationID = 0;
		navigationNumSE ++;
	}else
	{
		navigationID = 1;
		navigationNumSE = 0;
	}

/*
	maxOccIndexSE = -1;
	maxOccSE = 0;
	for(i=0; i<4; i++)
	{
		if(occsNumSE[i]>maxOccSE)
		{
			maxOccSE = occsNumSE[i];
			maxOccIndexSE = i;
		}
	}

	maxOccIndexPE = -1;
	maxOccPE = 0;
	for(i=0; i<4; i++)
	{
		if(occsNumPE[i]>maxOccPE)
		{
			maxOccPE = occsNumPE[i];
			maxOccIndexPE = i;
		}
	}

	// check the scores and occsNums of PE and SE
	if(maxOccIndexPE==-1 && maxOccIndexSE==-1)
	{
		kmers[0] = kmers[1] = NULL;
		thiskmerseq = 0;
	}else if(maxOccIndexPE!=-1 && maxOccIndexSE==-1)
	{
		kmers[0] = tmp_kmersPE[0];
		kmers[1] = tmp_kmersPE[1];
		thiskmerseq = tmp_kmerseqPE;
	}else if(maxOccIndexPE==-1 && maxOccIndexSE!=-1)
	{
		kmers[0] = tmp_kmersSE[0];
		kmers[1] = tmp_kmersSE[1];
		thiskmerseq = tmp_kmerseqSE;
	}else // if(maxOccIndexPE!=-1 && maxOccIndexSE!=-1)
	{
		if(maxOccIndexPE!=maxOccIndexSE)
		{
			if(maxOccPE>=MIN_CONNECT_KMER_NUM)
			{
				kmers[0] = tmp_kmersPE[0];
				kmers[1] = tmp_kmersPE[1];
				thiskmerseq = tmp_kmerseqPE;
			}else if(maxOccSE>=MIN_CONNECT_KMER_NUM)
			{
				kmers[0] = tmp_kmersSE[0];
				kmers[1] = tmp_kmersSE[1];
				thiskmerseq = tmp_kmerseqSE;
			}else
			{
				kmers[0] = kmers[1] = NULL;
				thiskmerseq = 0;
			}
		}else
		{
			kmers[0] = tmp_kmersPE[0];
			kmers[1] = tmp_kmersPE[1];
			thiskmerseq = tmp_kmerseqPE;
		}
	}
*/

	return SUCCESSFUL;
}


/**
 * Get next scafKmer in local assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNextKmerByPEInScaf(int contigNodesNum, int contigID, int contigLen, int contigOrient, int assemblyRound)
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
		occsNumPE[i] = 0;

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
		if(computeKmerOccNumUnlockedByPEInScaf(tmp_kmers[base_index], occsNumPE+base_index, contigID, contigLen, contigOrient, contigNodesNum, assemblyRound)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute kmer occsNums by unlocked reads in local assembly, error!\n", __LINE__, __func__);
			return FAILED;
		}

		//========================= condition 1 ==============================
		//if(score[base_index]>0)
		if(occsNumPE[base_index]>=minKmerOccPE)
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

	//开始计算每个kmer得分
	validKmerNum = 0;
	maxOcc = 0, secondOcc = 0;
	maxOccIndex = -1, secondOccIndex = -1;
	for(i=0; i<4; i++)
	{
		occsNumPE[i] = 0;

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
				if(tmp_kmers[i][0]->multiplicity+tmp_kmers[i][1]->multiplicity<minKmerOccPE)
				{
					continue;
				}
			}else if(tmp_kmers[i][0]!=NULL)
			{
				if(tmp_kmers[i][0]->multiplicity<minKmerOccPE)
				{
					continue;
				}
			}else if(tmp_kmers[i][1]!=NULL)
			{
				if(tmp_kmers[i][1]->multiplicity<minKmerOccPE)
				{
					continue;
				}
			}
		}

		//########################## end ########################//
		//if(contigNodesNum>READ_LEN)
//		if(contigNodesNum>kmer_len)
//		//if(contigNodesNum>=kmer_len)  //--bad result
//		{
//			if(computeLongKmerScoreByPEInScaf(tmp_kmers[i], scorePE+i, occsNumPE+i, kmer_len, contigID, contigLen, contigOrient, contigNodesNum, assemblyRound)==FAILED)
//			{
//				printf("line=%d, In %s(), cannot compute occsNum of long k-mers, error!\n", __LINE__, __func__);
//				return FAILED;
//			}
//		}else
//		{
			if(computeKmerOccNumByPEInScaf(tmp_kmers[i], occsNumPE+i, contigID, contigLen, contigOrient, contigNodesNum, assemblyRound)==FAILED)
			{
				printf("line=%d, In %s(), cannot compute k-mer occsNums, error!\n", __LINE__, __func__);
				return FAILED;
			}
//		}

		//score[i] = computeKmerScore(tmp_kmers[i], occsNum+i, assemblingreads, numassemblingreads);

		//========================= condition 3 ==============================
		//if(score[i]>=MIN_SCORE_THRESHOLD/**MIN_SCORE_FACTOR*/)
		//if(score[i]>=MIN_SCORE_THRESHOLD && occsNum[i]>=MIN_CONNECT_KMER_NUM)
		//if(occsNum[i]>=MIN_CONNECT_KMER_NUM)
		if(occsNumPE[i]>0)
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
			if(maxOcc<occsNumPE[j])
			{
				secondOcc = maxOcc;
				secondOccIndex = maxOccIndex;
				maxOcc = occsNumPE[j];
				maxOccIndex = j;
			}else if(secondOcc<occsNumPE[j])
			{
				secondOcc = occsNumPE[j];
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
	if(validKmerNum>1 && ((secondOcc/maxOcc>SECOND_FIRST_OCC_RATIO && secondOcc>minKmerOccPE) || (secondOcc>maxSecondOcc)))
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
		scafKmers[0] = scafKmers[1] = NULL;
		return SUCCESSFUL;
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


/**
 *   计算kmer得分.
 *   如果决策表中含有锁定的reads, 则只考虑锁定的reads.
 *   否则,考虑全部的reads.
 *
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
*/
short computeKmerOccNumByPEInScaf(scafKmer *tmp_kmers[2], int *occNum, int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound)
{
	if(lockedReadsNumInScaf>=lockedReadsNumThres)
	{ //如果决策表中含有锁定的reads, 则只考虑锁定的reads.
		if(computeKmerOccNumLockedByPEInScaf(tmp_kmers, occNum, contigID, contigLen, contigOrient, contigNodesNum, assemblyRound)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute kmer occNum in local assembly, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{ //否则,考虑全部的reads.
		//计算决策表中的reads的得分
		if(computeKmerOccNumUnlockedByPEInScaf(tmp_kmers, occNum, contigID, contigLen, contigOrient, contigNodesNum, assemblyRound)==FAILED)
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
short computeKmerOccNumUnlockedByPEInScaf(scafKmer *tmp_kmers[2], int *occNum, int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound)
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
		if(scafAssemblingReadArr[i].matedFlag==YES && scafAssemblingReadArr[i].reserved==0 && scafAssemblingReadArr[i].lastpos>0)
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

short computeKmerOccNumLockedByPEInScaf(scafKmer *tmp_kmers[2], int *occNum, int contigID, int contigLen, int contigOrient, int contigNodesNum, int assemblyRound)
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
			if(scafAssemblingReadArr[i].matedFlag==YES && scafAssemblingReadArr[i].reserved==0 && scafAssemblingReadArr[i].lastpos>0)
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
