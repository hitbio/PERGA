/*
 * scaffoldingContig.c
 *
 *  Created on: Apr 1, 2011
 *      Author: zhuxiao
 */

#include "inc/scafStdinc.h"
#include "inc/scafExtvab.h"


#define LEN 100

/**
 * 从文件中取得contig.
 */
FILE *contigFp;
char *contig;
int contigID;
char tmp2[LEN];
unsigned long length = 0;

FILE *fp_out;

/**
 * 从文件中取得contig.
 */
short generateSingleContigFasta()
{

	char ch, *pseq;

	ch = fgetc(contigFp);
	if(ch=='>')//start point of a new contig
	{
		fscanf(contigFp,"%d %s %ld", &contigID, tmp2, &length);
	}
	else
	{
		return 1;
	}

	pseq = contig = malloc(length+1);
	ch = fgetc(contigFp);
	while(ch!='>' && ch!=EOF)
	{
		switch(ch)
		{
			case 'A':
			case 'C':
			case 'G':
			case 'T':
				*pseq = ch;
				pseq++;
		}
		ch = fgetc(contigFp);
	}
	*pseq = '\0';

	if(ch=='>')
	{
		fseek(contigFp, -1, SEEK_CUR);
	}
	return 0;
}

/**
 * get contigs and output them in fasta format.
 *  @return:
 *  	if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateContigsFasta(char *outContigFileName, char *srcContigFileName)
{
	printf("begin generating contigs in fasta ...\n");

	char fullSrcContigFileName[256];
	char fullOutContigFileName[256];

	strcpy(fullSrcContigFileName, inputPathStr);
	strcat(fullSrcContigFileName, srcContigFileName);
	strcpy(fullOutContigFileName, outputPathStr);
	strcat(fullOutContigFileName, outContigFileName);


	contigFp = fopen(fullSrcContigFileName, "r");
	if(contigFp==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fullSrcContigFileName);
		return FAILED;
	}

	fp_out = fopen(fullOutContigFileName, "w");
	if(fp_out==NULL)
	{
		printf("line=%d, In %s(), cannot create [ %s ], error!\n", __LINE__, __func__, fullOutContigFileName);
		return FAILED;
	}

	while(!generateSingleContigFasta())
	{
		fprintf(fp_out, ">%d\t%s\t%lu\n%s\n", contigID, tmp2, length, contig);
		free(contig);
	}
	fclose(fp_out);

	fclose(contigFp);

	printf("generating contigs end.\n");

	return SUCCESSFUL;
}


/**
 * Filter the short contigs.
 *  @return:
 *  	if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short filterShortContigs(char *contigsFileFiltered, char *contigsFileUltraShort, char *contigsFile)
{
	printf("=========== Begin filtering ultra-short contigs, please wait ...\n");

	FILE *fpContigs, *fpContigsFiltered, *fpContigsUltraShort;
	int64_t contigID, contigLen, maxContigLen, returnFlag;
	char *contigSeq, *contigHeadTitle;

	fpContigs = fopen(contigsFile, "r");
	if(fpContigs==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigsFile);
		return FAILED;
	}

	fpContigsFiltered = fopen(contigsFileFiltered, "w");
	if(fpContigsFiltered==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigsFileFiltered);
		return FAILED;
	}

	fpContigsUltraShort = fopen(contigsFileUltraShort, "w");
	if(fpContigsUltraShort==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigsFileUltraShort);
		return FAILED;
	}

	contigHeadTitle = (char *) malloc(1000 * sizeof(char));
	if(contigHeadTitle==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(getMaxContigLenFromFile(&maxContigLen, contigsFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the maximal contig length, error!\n", __LINE__, __func__);
		return FAILED;
	}

	contigSeq = (char *) malloc ((maxContigLen+1)*sizeof(char));
	if(contigSeq==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//get contigs in fasta format
	contigID = 0;
	while((returnFlag=getSingleContigFastaFromFile(fpContigs, contigHeadTitle, contigSeq, &contigLen))==SUCCESSFUL)
	{
		contigID ++;

		if(contigLen>=minContigLenThres)
		{
			// output the filtered contigs to file
			if(outputSingleContigToFile(fpContigsFiltered, contigID, contigHeadTitle, contigSeq, contigLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot output the filtered contig to file, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else
		{
			// output the ultra-short contigs to file
			if(outputSingleContigToFile(fpContigsUltraShort, contigID, contigHeadTitle, contigSeq, contigLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot output the filtered contig to file, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}

	free(contigSeq);
	contigSeq = NULL;

	free(contigHeadTitle);
	contigHeadTitle = NULL;

	fclose(fpContigsFiltered);
	fpContigsFiltered = NULL;
	fclose(fpContigsUltraShort);
	fpContigsUltraShort = NULL;

	fclose(fpContigs);
	fpContigs = NULL;

	// handle the error situation
	if(returnFlag==ERROR)
	{
		printf("line=%d, In %s(), cannot get the contig: %ld, error!\n", __LINE__, __func__, contigID);
		return FAILED;
	}

	printf("=========== End filtered ultra-short contigs.\n");

	return SUCCESSFUL;
}

/**
 * Output single contig to file.
 *  Format:
 *  	(1) head line: >head title;
 *  	(2) base sequence.
 *
 *  @return:
 *  	if succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputSingleContigToFile(FILE *fpContigs, int64_t contigID, char *contigHeadTitle, char *contigSeq, int64_t contigLen)
{
	fprintf(fpContigs, ">%s\n", contigHeadTitle);
	fprintf(fpContigs, "%s\n", contigSeq);

	return SUCCESSFUL;
}
