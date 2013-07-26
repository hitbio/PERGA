/*
 * correction.c
 *
 *  Created on: Dec 28, 2012
 *      Author: zhuxiao
 */


#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Correct read sequence and output the new bases to file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short correctRead(assemblingreadtype *assemblingRead, contigtype *contigArray, graphtype *graph, FILE *fpReadCorrected)
{
	// ################################## Debug information ###############################
	//if(assemblingRead->rid==1371453)
	//{
	//	printf("line=%d, In %s(), rid=%ld\n", __LINE__, __func__, (int64_t)assemblingRead->rid);
	//}
	// ################################## Debug information ###############################

	if((assemblingRead->orientation==ORIENTATION_PLUS && assemblingRead->firstBasePos!=0)
		|| (assemblingRead->orientation==ORIENTATION_MINUS && assemblingRead->firstBasePos!=assemblingRead->seqlen-1))
	{
		// get omitted read bases
		if(getOmittedReadBases(readSeqAlign, &readSeqLenAlign, assemblingRead)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the omitted read bases for correction, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// get omitted contig bases
		if(getOmittedContigBases(contigSeqAlign, &contigSeqLenAlign, assemblingRead, contigArray)==FAILED)
		{
			printf("line=%d, In %s(), cannot get the omitted contig bases for correction, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// alignment
		if(generateAlignment(alignResultSeqArr, matchFlagArr, &alignResultSeqLen, readSeqAlign, readSeqLenAlign, contigSeqAlign, contigSeqLenAlign)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the alignment between read bases and contig bases, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// fill erroneous bases and sort them
		if(fillErrBasesSingleRead(assemblingRead, matchFlagArr, alignResultSeqLen, contigSeqAlign, contigSeqLenAlign)==FAILED)
		{
			printf("line=%d, In %s(), cannot generate the alignment between read bases and contig bases, error!\n", __LINE__, __func__);
			return FAILED;
		}

	}

	// generate new read sequence
	if(generateNewReadseq(readseqCorrected, &readseqLenCorrected, assemblingRead)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate the new read bases, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// update the read in read set and output the new read sequence to file
	if(updateAndSaveReadToFile(readseqCorrected, readseqLenCorrected, assemblingRead, graph, fpReadCorrected)==FAILED)
	{
		printf("line=%d, In %s(), cannot update the read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Get omitted bases of a read.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getOmittedReadBases(char *readSeqAlign, int32_t *readSeqLenAlign, assemblingreadtype *assemblingRead)
{
	int32_t i, firstBasePos, baseInt, entriesNum, baseNumLastEntry, entryRow, entryPos, seqPos;
	uint64_t *readseqInt;

	readseqInt = assemblingRead->readseq;
	entriesNum = assemblingRead->entriesNumReadseq;
	baseNumLastEntry = assemblingRead->baseNumLastEentryReadseq;
	firstBasePos = assemblingRead->firstBasePos;

	if(assemblingRead->orientation==ORIENTATION_PLUS)
	{ // plus orientation
		// get the read bases from positions [0, ..., firstBasePos-1]
		entryRow = entryPos = 0;
		for(i=0; i<firstBasePos; i++)
		{
			if(entryRow==entriesNum-1)
				baseInt = (readseqInt[entryRow] >> (2*(baseNumLastEntry-1)-2*entryPos)) & 3;
			else
				baseInt = (readseqInt[entryRow] >> (62-2*entryPos)) & 3;

			entryPos ++;
			if(entryPos==32)
			{
				entryRow ++;
				entryPos = 0;
			}
			switch(baseInt)
			{
				case 0: readSeqAlign[i] = 'A'; break;
				case 1: readSeqAlign[i] = 'C'; break;
				case 2: readSeqAlign[i] = 'G'; break;
				case 3: readSeqAlign[i] = 'T'; break;
			}
		}

		readSeqAlign[i] = '\0';
		*readSeqLenAlign = i;
	}else
	{ // minus orientation
		// get the read bases from positions [seqLen-1 ... firstBasePos-1]
		seqPos = 0;
		entryRow = (assemblingRead->seqlen - 1) >> 5;
		entryPos = (assemblingRead->seqlen - 1) % 32;
		for(i=assemblingRead->seqlen-1; i>firstBasePos; i--, seqPos++)
		{
			if(entryRow==entriesNum-1)
				baseInt = (readseqInt[entryRow] >> (2*(baseNumLastEntry-1)-2*entryPos)) & 3;
			else
				baseInt = (readseqInt[entryRow] >> (62-2*entryPos)) & 3;

			entryPos --;
			if(entryPos<0)
			{
				entryRow --;
				entryPos = 31;
			}

			switch(baseInt)
			{
				case 0: readSeqAlign[seqPos] = 'T'; break;
				case 1: readSeqAlign[seqPos] = 'G'; break;
				case 2: readSeqAlign[seqPos] = 'C'; break;
				case 3: readSeqAlign[seqPos] = 'A'; break;
			}
		}

		readSeqAlign[seqPos] = '\0';
		*readSeqLenAlign = seqPos;
	}

	return SUCCESSFUL;
}

/**
 * Get omitted bases of contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getOmittedContigBases(char *contigSeqAlign, int32_t *contigSeqLenAlign, assemblingreadtype *assemblingRead, contigtype *contigArray)
{
	int32_t i, firstContigPos, seqPos, startContigPos;

	firstContigPos = assemblingRead->firstContigPos;
	if(assemblingRead->orientation==ORIENTATION_PLUS)
		startContigPos = firstContigPos - assemblingRead->firstBasePos;
	else
		startContigPos = firstContigPos - (assemblingRead->seqlen - 1 - assemblingRead->firstBasePos);

	seqPos = 0;
	for(i=startContigPos; i<firstContigPos; i++, seqPos++)
	{
		switch(contigArray[i].base)
		{
			case 0: contigSeqAlign[seqPos] = 'A'; break;
			case 1: contigSeqAlign[seqPos] = 'C'; break;
			case 2: contigSeqAlign[seqPos] = 'G'; break;
			case 3: contigSeqAlign[seqPos] = 'T'; break;
		}
	}
	contigSeqAlign[seqPos] = '\0';
	*contigSeqLenAlign = seqPos;

	return SUCCESSFUL;
}

/**
 * Generate new read sequence and output it to file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateNewReadseq(uint64_t *readseqCorrected, int32_t *newSeqLen, assemblingreadtype *assemblingRead)
{
	int32_t i, basePos, newBasePos, seqLen, baseInt, itemNumErrBaseArray;
	int32_t entriesNum, baseNumLastEntry, entryRow, entryPos, newEntryRow, newEntryPos;
	int32_t rowErrBaseArr, errBasePos, newBase, errKind, endBasePos, errExistFlag;
	uint64_t *readseqInt;
	errBase_t *pErrBaseArray;

	seqLen = assemblingRead->seqlen;
	pErrBaseArray = assemblingRead->errBaseArr;
	itemNumErrBaseArray = assemblingRead->errBaseNum;

	readseqInt = assemblingRead->readseq;
	entriesNum = assemblingRead->entriesNumReadseq;
	baseNumLastEntry = assemblingRead->baseNumLastEentryReadseq;

	// generate new readseq and stored into readseqCorrected
	if(assemblingRead->orientation==ORIENTATION_PLUS)
		rowErrBaseArr = 0;
	else
		rowErrBaseArr = itemNumErrBaseArray - 1;

	entryRow = newEntryRow = 0;
	entryPos = newEntryPos = 0;
	basePos = newBasePos = 0;
	readseqCorrected[newEntryRow] = 0;
	while(basePos<seqLen)
	{
		if(assemblingRead->orientation==ORIENTATION_PLUS)
		{ // plus orientation
			if(rowErrBaseArr<itemNumErrBaseArray)
			{
				errBasePos = pErrBaseArray[rowErrBaseArr].basePos;
				errKind = pErrBaseArray[rowErrBaseArr].errKind;
				newBase = pErrBaseArray[rowErrBaseArr].newBase;
				rowErrBaseArr ++;

				endBasePos = errBasePos - 1;
				errExistFlag = YES;
			}else
			{
				endBasePos = seqLen - 1;
				errExistFlag = NO;
			}
		}else
		{ // minus orientation
			if(rowErrBaseArr>=0)
			{
				errBasePos = pErrBaseArray[rowErrBaseArr].basePos;
				errKind = pErrBaseArray[rowErrBaseArr].errKind;
				newBase = pErrBaseArray[rowErrBaseArr].newBase;
				rowErrBaseArr --;

				endBasePos = errBasePos - 1;
				errExistFlag = YES;
			}else
			{
				endBasePos = seqLen - 1;
				errExistFlag = NO;
			}
		}

		// copy the correct read bases
		for(i=basePos; i<=endBasePos; i++)
		{
			if(entryRow==entriesNum-1)
				baseInt = (readseqInt[entryRow] >> (2*(baseNumLastEntry-1)-2*entryPos)) & 3;
			else
				baseInt = (readseqInt[entryRow] >> (62-2*entryPos)) & 3;

			entryPos ++;
			if(entryPos==32)
			{
				entryRow ++;
				entryPos = 0;
			}

			readseqCorrected[newEntryRow] = (readseqCorrected[newEntryRow] << 2) | baseInt;
			newEntryPos ++;
			if(newEntryPos==32)
			{
				newEntryRow ++;
				newEntryPos = 0;
				readseqCorrected[newEntryRow] = 0;
			}

			basePos ++;
			newBasePos ++;
		}

		if(errExistFlag==YES)
		{
			// copy the new base
			if(errKind==1)
			{
				readseqCorrected[newEntryRow] = (readseqCorrected[newEntryRow] << 2) | newBase;

				entryPos ++;
				if(entryPos==32)
				{
					entryRow ++;
					entryPos = 0;
				}

				newEntryPos ++;
				if(newEntryPos==32)
				{
					newEntryRow ++;
					newEntryPos = 0;
					readseqCorrected[newEntryRow] = 0;
				}

				basePos ++;
				newBasePos ++;
			}else if(errKind==2)
			{
				readseqCorrected[newEntryRow] = (readseqCorrected[newEntryRow] << 2) | newBase;

				newEntryPos ++;
				if(newEntryPos==32)
				{
					newEntryRow ++;
					newEntryPos = 0;
					readseqCorrected[newEntryRow] = 0;
				}

				newBasePos ++;
			}else if(errKind==3)
			{
				entryPos ++;
				if(entryPos==32)
				{
					entryRow ++;
					entryPos = 0;
				}

				basePos ++;
			}else
			{
				printf("line=%d, In %s(), invalid errKind=%d, rid=%ld, error!\n", __LINE__, __func__, errKind, (int64_t)assemblingRead->rid);
				return FAILED;
			}
		}
	}

	*newSeqLen = newBasePos;

	return SUCCESSFUL;
}

/**
 * Update the read in read set and output it to file.
 *  File format:
 *  	(1) read ID;
 *  	(2) seqLen;
 *  	(3) sequence.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateAndSaveReadToFile(uint64_t *readseqCorrected, int32_t newSeqLen, assemblingreadtype *assemblingRead, graphtype *graph, FILE *fpReadCorrected)
{
	int32_t readBlockID, rowNumInReadBlock, maxItemNumPerReadBlock, entriesNum;
	readBlock_t *readBlockArr;
	read_t *pRead;
	correctedRead_t correctedRead;

	readBlockArr = graph->readSet->readBlockArr;
	maxItemNumPerReadBlock = graph->readSet->maxItemNumPerReadBlock;

	readBlockID = (assemblingRead->rid - 1) / maxItemNumPerReadBlock;
	rowNumInReadBlock = (assemblingRead->rid - 1) % maxItemNumPerReadBlock;
	pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;

	pRead->modified = YES;
	//pRead->seqlen = newSeqLen;

	correctedRead.rid = assemblingRead->rid;
	correctedRead.seqlen = newSeqLen;
	//correctedRead.successFlag = YES;
	//correctedRead.validFlag = YES;


	if(fwrite(&correctedRead, sizeof(correctedRead_t), 1, fpReadCorrected)!=1)
	{
		printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
		return FAILED;
	}

	entriesNum = ((newSeqLen - 1) >> 5) + 1;
	if(fwrite(readseqCorrected, sizeof(uint64_t), entriesNum, fpReadCorrected)!=entriesNum)
	{
		printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
		return FAILED;
	}

/*
	fprintf(fpReadCorrected, "%ld\t%d\n", (int64_t)correctedRead.rid, correctedRead.seqlen);
	fprintf(fpReadCorrected, "\t%s\n", getReadBaseByInt(assemblingRead->readseq, assemblingRead->seqlen));
	fprintf(fpReadCorrected, "\t%s\n", getReadBaseByInt(readseqCorrected, correctedRead.seqlen));

	fflush(fpReadCorrected);
*/

	return SUCCESSFUL;
}

/**
 * Fill erroneous bases in a read.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillErrBasesSingleRead(assemblingreadtype *assemblingRead, char *matchFlagArr, int32_t alignResultSeqLen, char *contigSeqAlign, int32_t contigSeqLenAlign)
{
	int32_t i, j, errBaseNum, newErrBaseNum, newTotalErrBaseNum, baseInt, readPosTmp, contigPosTmp;
	errBase_t *errBaseArr;

	errBaseArr = assemblingRead->errBaseArr;
	errBaseNum = assemblingRead->errBaseNum;

	newErrBaseNum = 0;
	for(i=0; i<alignResultSeqLen; i++)
		if(matchFlagArr[i] != '0')
			newErrBaseNum ++;

	if(newErrBaseNum>0)
	{
		newTotalErrBaseNum = errBaseNum + newErrBaseNum;
		if(newTotalErrBaseNum>2*maxErrBaseNumInCorrection+1)
		{
			printf("line=%d, In %s(), newTotalErrBaseNum=%d > %d, error!\n", __LINE__, __func__, newTotalErrBaseNum, 2*maxErrBaseNumInCorrection+1);
			return FAILED;
		}
		assemblingRead->errBaseNum = newTotalErrBaseNum;

		// copy the erroneous bases to new locations
		for(i=errBaseNum-1; i>=0; i--)
		{
			if(memcpy(errBaseArr+i+newErrBaseNum, errBaseArr+i, sizeof(errBase_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		// fill erroneous bases information
		if(assemblingRead->orientation==ORIENTATION_PLUS)
		{ // plus orientation
			j = 0;
			readPosTmp = 0;
			contigPosTmp = 0;
			for(i=0; i<alignResultSeqLen; i++)
			{
				if(matchFlagArr[i]=='0')
				{
					readPosTmp ++;
					contigPosTmp ++;
				}else if(matchFlagArr[i]=='1')
				{ // mismatch
					switch(contigSeqAlign[contigPosTmp])
					{
						case 'A': baseInt = 0; break;
						case 'C': baseInt = 1; break;
						case 'G': baseInt = 2; break;
						case 'T': baseInt = 3; break;
						default: printf("line=%d, In %s(), error base %c\n", __LINE__, __func__, contigSeqAlign[contigPosTmp]); return FAILED;
					}

					errBaseArr[j].basePos = readPosTmp;
					errBaseArr[j].newBase = baseInt;
					errBaseArr[j].errKind = 1;
					j ++;

					readPosTmp ++;
					contigPosTmp ++;
				}else if(matchFlagArr[i]=='2')
				{ // deletion in read, gap in read
					switch(contigSeqAlign[contigPosTmp])
					{
						case 'A': baseInt = 0; break;
						case 'C': baseInt = 1; break;
						case 'G': baseInt = 2; break;
						case 'T': baseInt = 3; break;
						default: printf("line=%d, In %s(), error base %c\n", __LINE__, __func__, contigSeqAlign[contigPosTmp]); return FAILED;
					}

					errBaseArr[j].basePos = readPosTmp;
					errBaseArr[j].newBase = baseInt;
					errBaseArr[j].errKind = 2;
					j ++;

					contigPosTmp ++;
				}else if(matchFlagArr[i]=='3')
				{
					errBaseArr[j].basePos = readPosTmp;
					errBaseArr[j].newBase = '-';
					errBaseArr[j].errKind = 3;
					j ++;

					readPosTmp ++;
				}else
				{
					printf("line=%d, In %s(), error matchFlag %c\n", __LINE__, __func__, matchFlagArr[i]);
					return FAILED;
				}
			}

			// ########################## Debug information ########################
			if(j!=newErrBaseNum)
			{
				printf("line=%d, In %s(), j=%d != newErrBaseNum=%d, error!\n", __LINE__, __func__, j, newErrBaseNum);
				return FAILED;
			}
			// ########################## Debug information ########################
		}else
		{ // minus orientation

			readPosTmp = assemblingRead->seqlen - 1;
			contigPosTmp = 0;
			j = 0;
			for(i=0; i<alignResultSeqLen; i++)
			{
				if(matchFlagArr[i]=='0')
				{
					readPosTmp --;
					contigPosTmp ++;
				}else if(matchFlagArr[i]=='1')
				{ // mismatch
					switch(contigSeqAlign[contigPosTmp])
					{
						case 'A': baseInt = 3; break;
						case 'C': baseInt = 2; break;
						case 'G': baseInt = 1; break;
						case 'T': baseInt = 0; break;
						default: printf("line=%d, In %s(), error base %c\n", __LINE__, __func__, contigSeqAlign[contigPosTmp]); return FAILED;
					}

					errBaseArr[j].basePos = readPosTmp;
					errBaseArr[j].newBase = baseInt;
					errBaseArr[j].errKind = 1;
					j ++;

					readPosTmp --;
					contigPosTmp ++;
				}else if(matchFlagArr[i]=='2')
				{
					switch(contigSeqAlign[contigPosTmp])
					{
						case 'A': baseInt = 3; break;
						case 'C': baseInt = 2; break;
						case 'G': baseInt = 1; break;
						case 'T': baseInt = 0; break;
						default: printf("line=%d, In %s(), error base %c\n", __LINE__, __func__, contigSeqAlign[contigPosTmp]); return FAILED;
					}

					errBaseArr[j].basePos = readPosTmp;
					errBaseArr[j].newBase = baseInt;
					errBaseArr[j].errKind = 2;
					j ++;

					contigPosTmp ++;

				}else if(matchFlagArr[i]=='3')
				{
					errBaseArr[j].basePos = readPosTmp;
					errBaseArr[j].newBase = '-';
					errBaseArr[j].errKind = 3;
					j ++;

					readPosTmp --;
				}else
				{
					printf("line=%d, In %s(), error matchFlag %c\n", __LINE__, __func__, matchFlagArr[i]);
					return FAILED;
				}
			}

			// ########################## Debug information ########################
			if(j!=newErrBaseNum)
			{
				printf("line=%d, In %s(), j=%d != newErrBaseNum=%d, error!\n", __LINE__, __func__, j, newErrBaseNum);
				return FAILED;
			}
			// ########################## Debug information ########################
		}
	}

	return SUCCESSFUL;
}


//========================================================================================
//========================================================================================
//========================================================================================

/**
 * Adjust read set by corrected reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short AdjustReadsetByCorrectedReads(graphtype *graph, char *correctedReadsFile)
{
	FILE *fpCorrectedReadsTmp;
	correctedRead_t correctedReadTmp;
	int64_t rid, hashcode;
	int32_t seqLen, entriesNum, maxItemNumPerReadBlock, maxEntryNumReadseqBlock, readBlockID, rowNumInReadBlock;
	readSet_t *readSet;
	readBlock_t *readBlockArr;
	readseqBlock_t *readseqBlockArr;
	readseqHashItem_t *pReadseqHashItem;
	read_t *pRead;

	fpCorrectedReadsTmp = fopen(correctedReadsFile, "r");
	if(fpCorrectedReadsTmp==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, correctedReadsFile);
		return FAILED;
	}

	readSet = graph->readSet;
	readBlockArr = readSet->readBlockArr;
	readseqBlockArr = readSet->readseqBlockArr;
	maxItemNumPerReadBlock = readSet->maxItemNumPerReadBlock;
	maxEntryNumReadseqBlock = readSet->maxEntryNumReadseqBlock;


	// prepare the tmp new read sequence buffer and tmp readseqHashItem buffer
	pReadseqBlockTmp = readseqBlockArr + readSet->blocksNumReadseq - 1;
	pReadseqTmpDoing = pReadseqBlockTmp->readseqArr + pReadseqBlockTmp->rowsNum;

	pReadseqHashItemBlockTmp = readSet->readseqHashItemBlockArr + readSet->blocksNumReadseqHashItem - 1;
	pReadseqHashItemTmpDoing = pReadseqHashItemBlockTmp->readseqHashItemArr + pReadseqHashItemBlockTmp->itemNum;


	while(1)
	{
		if((fread(&correctedReadTmp, sizeof(correctedRead_t), 1, fpCorrectedReadsTmp)!=1)) // get the chrPos and index
		{
			if(feof(fpCorrectedReadsTmp))
			{
				break;
			}
			else
			{
				printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		rid = correctedReadTmp.rid;
		seqLen = correctedReadTmp.seqlen;
		entriesNum = ((seqLen - 1) >> 5) + 1;

		readBlockID = (rid - 1) / maxItemNumPerReadBlock;
		rowNumInReadBlock = (rid - 1) % maxItemNumPerReadBlock;
		pRead = readBlockArr[readBlockID].readArr + rowNumInReadBlock;

		// set the valid pointer
		if(pReadseqBlockTmp->rowsNum + entriesNum > maxEntryNumReadseqBlock)
		{
			// add new readseq block
			if(addNewBlockReadseq(readSet)==FAILED)
			{
				printf("line=%d, In %s(), cannot add new readseq block, Error!\n", __LINE__, __func__);
				return FAILED;
			}
			pReadseqBlockTmp = readSet->readseqBlockArr + readSet->blocksNumReadseq - 1;
			pReadseqTmpDoing = pReadseqBlockTmp->readseqArr;
		}

		// get the read sequence from file
		if((fread(pReadseqTmpDoing, sizeof(uint64_t), entriesNum, fpCorrectedReadsTmp)!=entriesNum)) // get the chrPos and index
		{
			printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
			return FAILED;
		}

		// get the hashcode and add the readseq
		hashcode = readseqHashInt(pReadseqTmpDoing, seqLen, entriesNum);
		pReadseqHashItem = getReadseqHashItemByHash(hashcode, pReadseqTmpDoing, seqLen, entriesNum, readSet);
		if(pReadseqHashItem)
		{ // set the readseq to the readseq blocks
			pRead->readseqBlockID = pReadseqHashItem->readseqBlockID;
			pRead->rowReadseqInBlock = pReadseqHashItem->rowReadseqInBlock;
		}
		else
		{ // add the new readseq item into readseq blocks and readseq hash item blocks
			// read block operations
			pRead->readseqBlockID = pReadseqBlockTmp->blockID;
			pRead->rowReadseqInBlock = pReadseqBlockTmp->rowsNum;

			// readseq hash item block operations
			// add new item in readseq hash item block
			pReadseqHashItemTmpDoing->readseqBlockID = pReadseqBlockTmp->blockID;
			pReadseqHashItemTmpDoing->rowReadseqInBlock = pReadseqBlockTmp->rowsNum;
			pReadseqHashItemTmpDoing->seqlen = seqLen;
			pReadseqHashItemTmpDoing->nextHashItemBlockID = readSet->readseqHashtable[hashcode].hashItemBlockID;
			pReadseqHashItemTmpDoing->nextItemRowHashItemBlock = readSet->readseqHashtable[hashcode].itemRowHashItemBlock;
			readSet->readseqHashtable[hashcode].hashItemBlockID = pReadseqHashItemBlockTmp->blockID;
			readSet->readseqHashtable[hashcode].itemRowHashItemBlock = pReadseqHashItemBlockTmp->itemNum;

			readSet->totalItemNumReadseqHashItem ++;
			pReadseqHashItemBlockTmp->itemNum ++;
			pReadseqHashItemTmpDoing ++;
			if(pReadseqHashItemBlockTmp->itemNum >= readSet->maxItemNumPerReadseqHashItemBlock)
			{
				// add new readseq hash item block
				if(addNewBlockReadseqHashItem(readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot add new readseq hash item block, Error!\n", __LINE__, __func__);
					return FAILED;
				}
				pReadseqHashItemBlockTmp = readSet->readseqHashItemBlockArr + readSet->blocksNumReadseqHashItem - 1;
				pReadseqHashItemTmpDoing = pReadseqHashItemBlockTmp->readseqHashItemArr;
			}

			// readseq block operations
			pReadseqTmpDoing += entriesNum;
			readSet->totalItemNumReadseq ++;
			pReadseqBlockTmp->rowsNum += entriesNum;

			// set the valid pointer
			if(pReadseqBlockTmp->rowsNum >= readSet->maxEntryNumReadseqBlock)
			{
				// add new readseq block
				if(addNewBlockReadseq(readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot add new readseq block, Error!\n", __LINE__, __func__);
					return FAILED;
				}
				pReadseqBlockTmp = readSet->readseqBlockArr + readSet->blocksNumReadseq - 1;
				pReadseqTmpDoing = pReadseqBlockTmp->readseqArr;
			}
		}

		pRead->seqlen = seqLen;
	}

	// the last readseq block is empty, remove it
	if(pReadseqBlockTmp->rowsNum==0)
	{
		free(pReadseqBlockTmp->readseqArr);
		graph->readSet->readseqBlockArr[graph->readSet->blocksNumReadseq-1].readseqArr = NULL;
		graph->readSet->blocksNumReadseq --;
	}

	// the last readseq hash item block is empty, remove it
	if(pReadseqHashItemBlockTmp->itemNum==0)
	{
		free(pReadseqHashItemBlockTmp->readseqHashItemArr);
		graph->readSet->readseqHashItemBlockArr[graph->readSet->blocksNumReadseqHashItem-1].readseqHashItemArr = NULL;
		graph->readSet->blocksNumReadseqHashItem --;
	}


	fclose(fpCorrectedReadsTmp);
	fpCorrectedReadsTmp = NULL;

	return SUCCESSFUL;
}





