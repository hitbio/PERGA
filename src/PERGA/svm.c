/*
 * svm.c
 *
 *  Created on: Feb 6, 2013
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Extension decision by SVM model.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short extensionDecisionBySvm(int32_t *naviExtensionFlag, svmSampleVector_t *svmSample, svmModel_t *svmModel)
{
	int32_t outClass;

	// classify the data
	if(mySvmClassifySE(&outClass, svmSample, svmModel)==FAILED)
	{
		printf("line=%d, In %s(), cannot classify data by SVM model, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(outClass==YES)
		*naviExtensionFlag = YES;
	else
		*naviExtensionFlag = NO;

	return SUCCESSFUL;
}

/**
 * Load the SVM model.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short loadSvmModel(svmModel_t **svmModel, char *kernelFuncName, char *supportVectorFile, char *alphaFile, char *biasFile, char *scaleDataFile)
{
	// initialize the memory for the SVM classifier model
	if(initMemSvmModel(svmModel, supportVectorFile, scaleDataFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the memory for SVM model, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// load kernel function name
	if(loadKernelFuncName(*svmModel, kernelFuncName)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the support vectors, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// load the support vectors
	if(loadSupportVectorData(*svmModel, supportVectorFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the support vectors, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// load the Alpha vector data
	if(loadAlphaVectorData(*svmModel, alphaFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the alpha vectors, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// load the bias
	if(loadBiasData(*svmModel, biasFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the bias data, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// load the scale data
	if(loadScaleData(*svmModel, scaleDataFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot load the scale data, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Initialize the memory for sample data before using SVM classifier.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemSvmModel(svmModel_t **svmModel, char *supportVectorFile, char *scaleDataFile)
{
	int32_t i, rowsNum, colsNum, tmp;
	FILE *fpTmp;

	// allocate the SVM parameter node
	*svmModel = (svmModel_t *) calloc (1, sizeof(svmModel_t));
	if((*svmModel)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get rows number and columns number
	if(getRowsColsNumSupportVectors(&((*svmModel)->rowsNumSupportVector), &((*svmModel)->colsNumSupportVector), supportVectorFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the number of columns for support vectors, error!", __LINE__, __func__);
		return FAILED;
	}
	rowsNum = (*svmModel)->rowsNumSupportVector;
	colsNum = (*svmModel)->colsNumSupportVector;

	// allocate memory for support vectors
	(*svmModel)->supportVectors = (svmSupportVector_t*) calloc(rowsNum, sizeof(svmSupportVector_t));
	if((*svmModel)->supportVectors==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=0; i<rowsNum; i++)
	{
		(*svmModel)->supportVectors[i].vectorData = (double*) malloc(sizeof(double)*colsNum);
		if((*svmModel)->supportVectors[i].vectorData==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	// allocate memory for alpha vector
	(*svmModel)->rowsNumAlphaVector = (*svmModel)->rowsNumSupportVector;
	(*svmModel)->alphaVector = (double*) calloc((*svmModel)->rowsNumAlphaVector, sizeof(double));
	if((*svmModel)->alphaVector==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate memory for scaleData
	fpTmp = fopen(scaleDataFile, "r");
	if(fpTmp==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, scaleDataFile);
		return FAILED;
	}

	fscanf(fpTmp, "%d\n", &tmp);

	if(tmp==1)
	{ // initialize the memory for scale data if necessary
		(*svmModel)->scaleData = (struct scaleDataNode *) malloc (sizeof(struct scaleDataNode));
		if((*svmModel)->scaleData==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		(*svmModel)->scaleData->colsNum = (*svmModel)->colsNumSupportVector;

		(*svmModel)->scaleData->scaleFactorVector = (double *) malloc ((*svmModel)->scaleData->colsNum * sizeof(double));
		if((*svmModel)->scaleData->scaleFactorVector==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		(*svmModel)->scaleData->shiftVector = (double *) malloc ((*svmModel)->scaleData->colsNum * sizeof(double));
		if((*svmModel)->scaleData->shiftVector==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		(*svmModel)->scaleData = NULL;
	}

	fclose(fpTmp);

	return SUCCESSFUL;
}

/**
 * initialize the memory for sample data before using SVM classifier.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initSampleMemSvm(svmModel_t *svmModel)
{
	int32_t colsNum;

	colsNum = svmModel->colsNumSupportVector;

	// allocate svmSample
	svmSample = (svmSampleVector_t*) malloc (sizeof(svmSampleVector_t));
	if(svmSample==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	svmSample->colsNum = colsNum;
	svmSample->vectorData = (double *) malloc(colsNum*sizeof(double));
	if(svmSample->vectorData==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate svmSampleSE
	svmSampleSE = (svmSampleVector_t*) malloc (sizeof(svmSampleVector_t));
	if(svmSampleSE==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	svmSampleSE->colsNum = colsNum;
	svmSampleSE->vectorData = (double *) malloc(colsNum*sizeof(double));
	if(svmSampleSE->vectorData==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// allocate svmSamplePE
	svmSamplePE = (svmSampleVector_t*) malloc (sizeof(svmSampleVector_t));
	if(svmSamplePE==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	svmSamplePE->colsNum = colsNum;
	svmSamplePE->vectorData = (double *) malloc(colsNum*sizeof(double));
	if(svmSamplePE->vectorData==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Free sample memory after using SVM model.
 */
void freeSampleMemSvm()
{
	free(svmSample);
	svmSample = NULL;

	free(svmSampleSE);
	svmSampleSE = NULL;

	free(svmSamplePE);
	svmSamplePE = NULL;
}

/**
 * Free SVM parameters.
 */
short freeSvmModel(svmModel_t **svmModel)
{
	int32_t i;

	for(i=0; i<(*svmModel)->rowsNumSupportVector; i++)
		free((*svmModel)->supportVectors[i].vectorData);
	free((*svmModel)->supportVectors);
	(*svmModel)->rowsNumSupportVector = 0;
	(*svmModel)->colsNumSupportVector = 0;

	free((*svmModel)->alphaVector);
	(*svmModel)->rowsNumAlphaVector = 0;

	(*svmModel)->bias = 0;

	if((*svmModel)->scaleData)
	{
		free((*svmModel)->scaleData->scaleFactorVector);
		free((*svmModel)->scaleData->shiftVector);
		(*svmModel)->scaleData->colsNum = 0;
		free((*svmModel)->scaleData);
		(*svmModel)->scaleData = NULL;
	}

	free(*svmModel);
	*svmModel = NULL;

	return SUCCESSFUL;
}

/**
 * Get the number of rows and columns for support vectors before using SVM classifier.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getRowsColsNumSupportVectors(int32_t *rowsNumSupportVector, int32_t *colsNumSupportVector, char *supportVectorFile)
{
	int32_t tabNum, rowsNum;
	FILE *fpSupportVector;
	char ch;

	fpSupportVector = fopen(supportVectorFile, "r");
	if(fpSupportVector==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, supportVectorFile);
		return FAILED;
	}

	// get one line and compute the columns number
	tabNum = 0;
	ch = fgetc(fpSupportVector);
	while(ch!=EOF && ch!='\n')
	{
		if(ch=='\t') tabNum ++;
		ch = fgetc(fpSupportVector);
	}

	*colsNumSupportVector = tabNum + 1;

	// compute the rows number
	rewind(fpSupportVector);
	rowsNum = 0;
	ch = fgetc(fpSupportVector);
	while(ch!=EOF)
	{
		if(ch=='\n') rowsNum ++;
		ch = fgetc(fpSupportVector);
	}

	*rowsNumSupportVector = rowsNum;

	fclose(fpSupportVector);

	return SUCCESSFUL;
}

/**
 * Load the kernel function name.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short loadKernelFuncName(svmModel_t *svmModel, char *kernelFuncName)
{
	int32_t i, rowsNum, colsNum, lineLen, colID;
	FILE *fpKernelFunc;
	char ch, *pch;

	fpKernelFunc = fopen(kernelFuncName, "r");
	if(fpKernelFunc==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, kernelFuncName);
		return FAILED;
	}

	i = 0;
	pch = svmModel->kernelFunc;
	ch = fgetc(fpKernelFunc);
	while(ch!=EOF && ch!='\n')
	{
		pch[i] = ch;
		i++;
		ch = fgetc(fpKernelFunc);
	}
	pch[i] = '\0';

	fclose(fpKernelFunc);

	return SUCCESSFUL;
}

/**
 * Load the support vectors data.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short loadSupportVectorData(svmModel_t *svmModel, char *supportVectorFile)
{
	int32_t i, rowsNum, colsNum, lineLen, colID;
	FILE *fpSupportVector;
	char *pch, lineStr[5000];

	rowsNum = svmModel->rowsNumSupportVector;
	colsNum = svmModel->colsNumSupportVector;

	fpSupportVector = fopen(supportVectorFile, "r");
	if(fpSupportVector==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, supportVectorFile);
		return FAILED;
	}

	for(i=0; i<rowsNum; i++)
	{
		if(readLine(lineStr, &lineLen, fpSupportVector)==FAILED)
		{
			printf("line=%d, In %s(), cannot read one line, error!\n", __LINE__, __func__);
			return FAILED;
		}

		colID = 0;
		pch = strtok(lineStr, "\t");
		while(pch)
		{
			svmModel->supportVectors[i].vectorData[colID++] = atof(pch);
			if(colID>colsNum)
			{
				printf("line=%d, In %s(), invalid column number, error!\n", __LINE__, __func__);
				return FAILED;
			}

			pch = strtok(NULL, "\t");
		}
	}

	fclose(fpSupportVector);

	return SUCCESSFUL;
}

/**
 * Load the alpha vectors data.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short loadAlphaVectorData(svmModel_t *svmModel, char *alphaVectorFile)
{
	int32_t i, rowsNum;
	FILE *fpAlphaVector;

	rowsNum = svmModel->rowsNumAlphaVector;

	fpAlphaVector = fopen(alphaVectorFile, "r");
	if(fpAlphaVector==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, alphaVectorFile);
		return FAILED;
	}

	for(i=0; i<rowsNum; i++)
		fscanf(fpAlphaVector, "%lf\n", svmModel->alphaVector+i);

	fclose(fpAlphaVector);

	return SUCCESSFUL;
}

/**
 * Load the bias data.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short loadBiasData(svmModel_t *svmModel, char *biasFile)
{
	FILE *fpBias;

	fpBias = fopen(biasFile, "r");
	if(fpBias==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, biasFile);
		return FAILED;
	}

	fscanf(fpBias, "%lf\n", &(svmModel->bias));

	fclose(fpBias);

	return SUCCESSFUL;
}

/**
 * Load the scale data.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short loadScaleData(svmModel_t *svmModel, char *scaleDataFile)
{
	int32_t i, colsNum, colID, lineLen;
	FILE *fpScaleData;
	char ch, *pch, lineStr[5000];
	double *pVector;

	fpScaleData = fopen(scaleDataFile, "r");
	if(fpScaleData==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, scaleDataFile);
		return FAILED;
	}

	// skip one line
	ch = fgetc(fpScaleData);
	while(ch!=EOF && ch!='\n') ch = fgetc(fpScaleData);

	// load the scale data
	if(svmModel->scaleData)
	{
		colsNum = svmModel->scaleData->colsNum;

		for(i=0; i<2; i++)
		{
			if(readLine(lineStr, &lineLen, fpScaleData)==FAILED)
			{
				printf("line=%d, In %s(), cannot read one line, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(i==0)
				pVector = svmModel->scaleData->scaleFactorVector;
			else
				pVector = svmModel->scaleData->shiftVector;

			colID = 0;
			pch = strtok(lineStr, "\t");
			while(pch)
			{
				pVector[colID++] = atof(pch);
				if(colID>colsNum)
				{
					printf("line=%d, In %s(), invalid column number, error!\n", __LINE__, __func__);
					return FAILED;
				}

				pch = strtok(NULL, "\t");
			}
		}
	}else
	{
		printf("line=%d, In %s(), scale data is empty, error!\n", __LINE__, __func__);
		return FAILED;
	}

	fclose(fpScaleData);

	return SUCCESSFUL;
}

/**
 * Read one line of a file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short readLine(char *lineStr, int32_t *lineLen, FILE *fp)
{
	char ch;
	int32_t len;

	if(fp==NULL)
	{
		printf("line=%d, In %s(), invalid file operation, error!\n", __LINE__, __func__);
		return FAILED;
	}

	len = 0;
	ch = fgetc(fp);
	while(ch!=EOF && ch!='\n')
	{
		lineStr[len++] = ch;
		ch = fgetc(fp);
	}
	lineStr[len] = '\0';

	*lineLen = len;

	return SUCCESSFUL;
}

/**
 * Use SVM to classify the a single end sample data.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short mySvmClassifySE(int32_t *outClass, svmSampleVector_t *svmSampleSE, svmModel_t *svmModelSE)
{
	int32_t i, supportVectorNum, dimNum;
	double cvalue;

	// use the scale data to adjust the sample data if necessary
	dimNum = svmSampleSE->colsNum;
	if(svmModelSE->scaleData)
		for(i=0; i<dimNum; i++)
			svmSampleSE->vectorData[i] = svmModelSE->scaleData->scaleFactorVector[i] * (svmSampleSE->vectorData[i] + svmModelSE->scaleData->shiftVector[i]);

	cvalue = 0;
	supportVectorNum = svmModelSE->rowsNumSupportVector;
	if(strcmp(svmModelSE->kernelFunc, "linear")==0)
		for(i=0; i<supportVectorNum; i++)
			cvalue += myLinearKernelFuntion(svmSampleSE->vectorData, svmModelSE->supportVectors[i].vectorData, dimNum) * svmModelSE->alphaVector[i];
	else if(strcmp(svmModelSE->kernelFunc, "rbf")==0)
		for(i=0; i<supportVectorNum; i++)
			cvalue += myRBFKernelFuntion(svmSampleSE->vectorData, svmModelSE->supportVectors[i].vectorData, dimNum) * svmModelSE->alphaVector[i];
	else if(strcmp(svmModelSE->kernelFunc, "polynomial")==0)
		for(i=0; i<supportVectorNum; i++)
			cvalue += myPolyKernelFuntion(svmSampleSE->vectorData, svmModelSE->supportVectors[i].vectorData, dimNum) * svmModelSE->alphaVector[i];
	else
	{
		printf("line=%d, In %s(), invalid kernel function: %s, error!\n", __LINE__, __func__, svmModelSE->kernelFunc);
		return FAILED;
	}
	cvalue += svmModelSE->bias;

	if(cvalue>=0)
		*outClass = 0;
	else
		*outClass = 1;

	return SUCCESSFUL;
}

/**
 * Linear kernel function K(x,y) for dimNum dimensional vectors.
 *  K(x,y) = <x,y>, where <x,y> is the dot product of x, y.
 *
 *  @return:
 *  	The function value should >= 0.
 *  	If succeeds, return the function value; otherwise, return -1.
 */
double myLinearKernelFuntion(double *x, double *y, int32_t dimNum)
{
	int32_t i;
	double dotPro;

	// compute the Norm-2 value
	dotPro = 0;
	for(i=0; i<dimNum; i++)
		dotPro += x[i] * y[i];

	return dotPro;
}

/**
 * RBF kernel function K(x,y) for dimNum dimensional vectors.
 *  K(x,y) = exp(-0.5*||x-y||^2), where ||x-y|| is the Norm-2 of x-y.
 *
 *  @return:
 *  	The function value should >= 0.
 *  	If succeeds, return the function value; otherwise, return -1.
 */
double myRBFKernelFuntion(double *x, double *y, int32_t dimNum)
{
	int32_t i;
	double fvalue, normSqureValue, dif;

	// compute the Norm-2 value
	normSqureValue = 0;
	for(i=0; i<dimNum; i++)
	{
		dif = x[i] - y[i];
		normSqureValue += dif * dif;
	}

	fvalue = exp(-0.5*normSqureValue);

	if(fvalue>=0)
		return fvalue;
	else
		return -1;
}

/**
 * Polynomials kernel function K(x,y) for dimNum dimensional vectors.
 *  K(x,y) = <x,y>*(1+<x,y>)^2, where <x,y> is the dot product of of x,y.
 *
 *  @return:
 *  	The function value.
 */
double myPolyKernelFuntion(double *x, double *y, int32_t dimNum)
{
	int32_t i;
	double dotProduct, K;

	// compute the product value
	dotProduct = 0;
	for(i=0; i<dimNum; i++)
		dotProduct += x[i] * y[i];

	K = dotProduct;
	for(i=2; i<=3; i++)
		K = K * (1 + dotProduct);

	//K = dotProduct * pow(1+dotProduct, 2);

	return K;
}

/**
 * Fill the sample data.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillSampleDataSVM(svmSampleVector_t *svmSample, double *svmFeatureArray)
{
	int32_t i;

	for(i=0; i<svmSample->colsNum; i++)
		svmSample->vectorData[i] = svmFeatureArray[i];

	return SUCCESSFUL;
}
