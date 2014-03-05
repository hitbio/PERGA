/*
 * dtRowHash.c
 *
 *  Created on: Dec 2, 2012
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Add a read to DTRow hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addReadToDTRowHashtable(int64_t rid, int32_t dtRow, dtRowIndex_t **pDtRowHashtable)
{
	int32_t hashcode;
	dtRowIndex_t *pDTRow, *head;

	hashcode = rid & RID_LOW_BITS_MASK;
	head = pDtRowHashtable[hashcode];

	pDTRow = (dtRowIndex_t *) malloc(sizeof(dtRowIndex_t));
	if(pDTRow==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	pDTRow->rid = rid;
	pDTRow->dtRow = dtRow;
	pDTRow->next = head;
	pDtRowHashtable[hashcode] = pDTRow;

	return SUCCESSFUL;
}

/**
 * Delete a read from DTRow hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short delReadFromDTRowHashtable(int64_t rid, int32_t dtRow, dtRowIndex_t **pDtRowHashtable)
{
	int32_t hashcode, successFlag;
	dtRowIndex_t *pDTRow, *pNextDTRow;

	successFlag = NO;
	hashcode = rid & RID_LOW_BITS_MASK;
	pDTRow = pDtRowHashtable[hashcode];

	if(pDTRow)
	{
		pNextDTRow = pDTRow->next;
		if(pDTRow->rid==rid && pDTRow->dtRow==dtRow)
		{
			pDtRowHashtable[hashcode] = pNextDTRow;
			free(pDTRow);
			successFlag = YES;
		}else
		{
			while(pNextDTRow)
			{
				if(pNextDTRow->rid==rid && pNextDTRow->dtRow==dtRow)
				{
					pDTRow->next = pNextDTRow->next;
					free(pNextDTRow);
					successFlag = YES;
					break;
				}

				pDTRow = pNextDTRow;
				pNextDTRow = pNextDTRow->next;
			}
		}
	}else
	{
		printf("line=%d, In %s(), the read rid=%ld having dtRow=%d does not exist in dtRowHashtable, error!\n", __LINE__, __func__, rid, dtRow);
		return FAILED;
	}

	if(successFlag==YES)
		return SUCCESSFUL;
	else
	{
		printf("line=%d, In %s(), the read rid=%ld having dtRow=%d does not exist in dtRowHashtable, error!\n", __LINE__, __func__, rid, dtRow);
		return FAILED;
	}
}

/**
 * Update the read in DTRowIndex hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateReadInDTRowHashtable(int64_t rid, int32_t dtRowOld, int32_t dtRowNew, dtRowIndex_t **pDtRowHashtable)
{
	int32_t hashcode, successFlag;
	dtRowIndex_t *pDTRow;

	successFlag = NO;
	hashcode = rid & RID_LOW_BITS_MASK;
	pDTRow = pDtRowHashtable[hashcode];

	while(pDTRow)
	{
		if(pDTRow->rid==rid && pDTRow->dtRow==dtRowOld)
		{
			pDTRow->dtRow = dtRowNew;
			successFlag = YES;
			break;
		}
		pDTRow = pDTRow->next;
	}

	if(successFlag==YES)
		return SUCCESSFUL;
	else
		return FAILED;
}


/**
 * Clean the DTRowIndex hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short cleanDTRowIndexHashtable(dtRowIndex_t **pDtRowHashtable)
{
	int i;
	dtRowIndex_t *pDTRow, *head;

	for(i=0; i<TABLE_SIZE_HASH_DTROWINDEX; i++)
	{
		head = pDtRowHashtable[i];
		while(head)
		{
			pDTRow = head->next;
			free(head);
			head = pDTRow;
		}
		pDtRowHashtable[i] = NULL;
	}

	return SUCCESSFUL;
}

/**
 * Get the proper dtRow in the DTRowIndex hash table.
 *  @return:
 *  	If succeeds, return proper dtRow; otherwise, return -1.
 */
int32_t getProperDtRow(assemblingreadtype *assemblingread, assemblingreadtype *assemblingreads, dtRowIndex_t **pDtRowHashtable)
{
	int32_t proBasepos, proDtRow, rid, hashcode;
	dtRowIndex_t *pDtRow;
	assemblingreadtype *pReadDt;

	proDtRow = -1;
	rid = assemblingread->rid;
	hashcode = rid & RID_LOW_BITS_MASK;
	pDtRow = pDtRowHashtable[hashcode];

	if(assemblingread->orientation==ORIENTATION_PLUS)
	{
		proBasepos = 0;
		while(pDtRow)
		{
			if(pDtRow->rid==rid)
			{
				pReadDt = assemblingreads + pDtRow->dtRow;
				if(pReadDt->orientation==ORIENTATION_PLUS)
				{
					pReadDt->reserved = 1;
					if(pReadDt->basePos>proBasepos)
					{
						proBasepos = pReadDt->basePos;
						proDtRow = pDtRow->dtRow;
					}
				}
			}
			pDtRow = pDtRow->next;
		}
	}else
	{
		proBasepos = INT_MAX;
		while(pDtRow)
		{
			if(pDtRow->rid==rid)
			{
				pReadDt = assemblingreads + pDtRow->dtRow;
				if(pReadDt->orientation==ORIENTATION_MINUS)
				{
					pReadDt->reserved = 1;
					//if(pReadDt->lastpos>0 && pReadDt->lastpos<proBasepos)
					if(pReadDt->basePos<proBasepos)
					{
						proBasepos = pReadDt->basePos;
						proDtRow = pDtRow->dtRow;
					}
				}
			}
			pDtRow = pDtRow->next;
		}
	}

	return proDtRow;
}

/**
 * Get the limited proper dtRow in the DTRowIndex hash table.
 *  @return:
 *  	If succeeds, return proper dtRow; otherwise, return -1.
 */
int32_t getProperDtRowLimited(assemblingreadtype *assemblingread, assemblingreadtype *assemblingreads, dtRowIndex_t **pDtRowHashtable, int limitLastpos)
{
	int32_t proBasePos, proDtRow, rid, hashcode;
	dtRowIndex_t *pDtRow;
	assemblingreadtype *pReadDt;

	proBasePos = limitLastpos;
	proDtRow = -1;
	rid = assemblingread->rid;
	hashcode = rid & RID_LOW_BITS_MASK;
	pDtRow = pDtRowHashtable[hashcode];

	if(assemblingread->orientation==ORIENTATION_PLUS)
	{
		while(pDtRow)
		{
			if(pDtRow->rid==rid)
			{
				pReadDt = assemblingreads + pDtRow->dtRow;
				if(pReadDt->orientation==ORIENTATION_PLUS)
				{
					if(pReadDt->basePos>proBasePos)
					{
						proBasePos = pReadDt->basePos;
						proDtRow = pDtRow->dtRow;
					}
				}
			}
			pDtRow = pDtRow->next;
		}
	}else
	{
		while(pDtRow)
		{
			if(pDtRow->rid==rid)
			{
				pReadDt = assemblingreads + pDtRow->dtRow;
				if(pReadDt->orientation==ORIENTATION_MINUS)
				{
					if(pReadDt->basePos==pReadDt->lastMatchedBasePos && pReadDt->basePos<proBasePos)
					{
						proBasePos = pReadDt->basePos;
						proDtRow = pDtRow->dtRow;
					}
				}
			}
			pDtRow = pDtRow->next;
		}
	}

	return proDtRow;
}

/**
 * Check whether the read already exist in decision table.
 *  @return:
 *  	If exists, return YES; otherwise, return NO.
 */
int32_t existReadInDecisionTable(int64_t rid, dtRowIndex_t **pDtRowHashtable)
{
	int32_t hashcode, existFlag;
	dtRowIndex_t *pDtRow;

	hashcode = rid & RID_LOW_BITS_MASK;
	pDtRow = pDtRowHashtable[hashcode];

	existFlag = NO;
	while(pDtRow)
	{
		if(pDtRow->rid==rid)
		{
			existFlag = YES;
			break;
		}

		pDtRow = pDtRow->next;
	}

	return existFlag;
}

/**
 * Get the the read already exist in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getExistReadInDT(assemblingreadtype **pReadDt, int64_t rid, assemblingreadtype *decisionTable, dtRowIndex_t **pDtRowHashtable)
{
	int32_t hashcode;
	dtRowIndex_t *pDtRow;

	hashcode = rid & RID_LOW_BITS_MASK;
	pDtRow = pDtRowHashtable[hashcode];

	*pReadDt = NULL;
	while(pDtRow)
	{
		if(pDtRow->rid==rid)
		{
			*pReadDt = decisionTable + pDtRow->dtRow;
			break;
		}

		pDtRow = pDtRow->next;
	}

	return SUCCESSFUL;
}

/**
 * Check whether the read with matched position already exist in decision table.
 *  @return:
 *  	If exists, return YES; otherwise, return NO.
 */
int32_t existReadWithPosInDecisionTable(int64_t rid, int32_t basePos, int32_t orientation, assemblingreadtype *decisionTable, dtRowIndex_t **pDtRowHashtable)
{
	int32_t hashcode, existFlag;
	dtRowIndex_t *pDtRow;
	assemblingreadtype *pReadDt;

	hashcode = rid & RID_LOW_BITS_MASK;
	pDtRow = pDtRowHashtable[hashcode];

	existFlag = NO;
	while(pDtRow)
	{
		if(pDtRow->rid==rid)
		{
			pReadDt = decisionTable + pDtRow->dtRow;
			if(pReadDt->orientation==orientation && pReadDt->basePos==basePos)
			{
				existFlag = YES;
				break;
			}
		}

		pDtRow = pDtRow->next;
	}

	return existFlag;
}

/**
 * Get the copy number of a read in decision table.
 *  @return:
 *  	If exists, return YES; otherwise, return NO.
 */
int32_t getCopyNumOfReadInDecisionTable(int32_t *copyNum, int64_t rid, dtRowIndex_t **pDtRowHashtable)
{
	int32_t hashcode;
	dtRowIndex_t *pDtRow;

	hashcode = rid & RID_LOW_BITS_MASK;
	pDtRow = pDtRowHashtable[hashcode];

	*copyNum = 0;
	while(pDtRow)
	{
		if(pDtRow->rid==rid)
			(*copyNum) ++;
		pDtRow = pDtRow->next;
	}

	return SUCCESSFUL;
}

