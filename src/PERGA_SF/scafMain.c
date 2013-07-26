/*
 * scafMain.c
 *
 *  Created on: Jul 17, 2012
 *      Author: zhuxiao
 */

#include "inc/scafStdinc.h"
#include "inc/scafGlobal.h"


int main(int argc, char** argv)
{
	if(parseCommandParasAndExe(argc, argv)==ERROR)
	{
		printf("Try '-help' for more information.\n");
		return FAILED;
	}

	return 0;
}

/**
 * Parse the input parameters.
 *  @parameters:
 *  	(1) -c <FILE>
 *  		Input contigs file in fasta format.
 *  	(2) -s <INT>
 *  		Region size of align region at contg ends.
 *  	(3) -m <INT>
 *  		Minimum size of contigs for scaffolding. Default is 100 bp.
 *  	(4) -p <INT>
 *  		Paired end mode for read files.
 *  			0 - do not treat reads as paired (default);
 *  			1 - reads are paired with the first read in the first file, and the second read in the second file;
 *  			2 - reads are paired and are interleaved within a single file.
 *  	(5)-f <FILES>
 *  		Read files in fastq format.
 *  	(6) -ins_len <FLOAT>
 *  		Insert size for paired end library.
 *  	(7) -ins_sdev <FLOAT>
 *  		Standard deviation of insert size for paired end library. If not specified and the insert size is specified,
 *  		this value is set to 0.1*insertSize.
 *  	(8) -o <STR>
 *  		Prefix of the output files.
 *  	(9) -d <STR>
 *  		Output directory for the output files. Default is "./".
 *  	(10) -h (-help)
 *  		Show help information.
 *
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL;
 *  	if the job is executed and failed, return ERROR;
 *  	otherwise, return FAILED.
 */
short parseCommandParasAndExe(int argc, char **argv)
{
	int i;

	char contigsFilePara[256];
	char outputPathPara[256];
	char outputPrefixPara[256];
	char *readFilesPara[256];
	char readFilesBuf[256][256];
	int readFileNumPara;
	int pairedModePara;
	int alignRegSizePara;
	double meanSizeInsertPara, standardDevPara;
	int minContigLenPara;
	char gapFillFlagPara[256];


	if(argc==1)
	{
		if(showUsageInfo()==FAILED)
		{
			printf("line=%d, In %s(), cannot show the usage information for SRGA, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		// reset the parameters
		contigsFilePara[0] = '\0';
		outputPathPara[0] = '\0';
		outputPrefixPara[0] = '\0';
		readFileNumPara = 0;
		pairedModePara = -1;
		alignRegSizePara = 0;
		meanSizeInsertPara = standardDevPara = 0;
		minContigLenPara = 0;
		gapFillFlagPara[0] = '\0';

		i = 1;
		while(i<argc)
		{

			if(strcmp(argv[i], "-c")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						strcpy(contigsFilePara, argv[i+1]);
					}else
					{
						printf("Exception: please specify the contigs file.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the contigs file.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-s")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						alignRegSizePara = atol(argv[i+1]);
						if(alignRegSizePara<0)
						{
							printf("Exception: please specify the correct align region size.\n");
							return FAILED;
						}
					}else
					{
						printf("Exception: please specify the align region size.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the align region size.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-m")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						minContigLenPara = atol(argv[i+1]);
						if(minContigLenPara<=0)
						{
							printf("Exception: please specify the correct minimal contig length.\n");
							return FAILED;
						}
					}else
					{
						printf("Exception: please specify the correct minimal contig length.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the minimal contig length.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-p")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						pairedModePara = atol(argv[i+1]);
						if(pairedModePara<0 || pairedModePara>2)
						{
							printf("Exception: please specify the correct paired mode.\n");
							return FAILED;
						}
					}else
					{
						printf("Exception: please specify the paired mode.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the paired mode.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-f")==0)
			{ // read files
				readFileNumPara = 0;

				i++;
				while(i<argc)
				{
					if(argv[i][0]!='-')
					{
						strcpy(readFilesBuf[readFileNumPara], argv[i]);
						readFilesPara[readFileNumPara] = readFilesBuf[readFileNumPara];
						readFileNumPara ++;
					}else
					{ // next option
//						printf("Exception: please specify the correct paired mode.\n");
//						return FAILED;
						break;
					}

					i ++;
				}

				if(readFileNumPara==0)
				{
					printf("Exception: please specify the read files.\n");
					return FAILED;
				}
			}else if(strcmp(argv[i], "-ins_len")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						meanSizeInsertPara = atof(argv[i+1]);
						if(meanSizeInsertPara<0)
						{
							printf("Exception: please specify the correct insert size.\n");
							return FAILED;
						}
					}else
					{
						printf("Exception: please specify the insert size.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the insert size.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-ins_sdev")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						standardDevPara = atof(argv[i+1]);
						if(standardDevPara<0)
						{
							printf("Exception: please specify the correct standard deviation of insert size.\n");
							return FAILED;
						}
					}else
					{
						printf("Exception: please specify the standard deviation of insert size.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the standard deviation of insert size.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-d")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						strcpy(outputPathPara, argv[i+1]);
					}else
					{
						printf("Exception: please specify the output directory.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the output directory.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-o")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						strcpy(outputPrefixPara, argv[i+1]);
					}else
					{
						printf("Exception: please specify the prefix for output files.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the prefix for output files.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-g")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						strcpy(gapFillFlagPara, argv[i+1]);
					}else
					{
						printf("Exception: please specify the gap filling flag.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the gap filling flag.\n");
					return FAILED;
				}
				i += 2;
			}
			else if(strcmp(argv[i], "-h")==0 || strcmp(argv[i], "-help")==0)
			{
				if(showUsageInfo()==FAILED)
				{
					printf("line=%d, In %s(), cannot show the usage information for globalMetrics, error!\n", __LINE__, __func__);
					return FAILED;
				}

				return FAILED;
			}else
			{
				if(argv[i][0]=='-')
				{
					printf("%s : unknown argument\n", argv[i]);
					printf("Please use -h or -help for more information.\n");
				}else
				{
					printf("%s : invalid argument\n", argv[i]);
					printf("Please use -h or -help for more information.\n");
				}

				return FAILED;
			}
		}

		if(strlen(contigsFilePara)==0)
		{
			printf("Exception: please specify the contig file.\n");
			return FAILED;
		}

		if(readFileNumPara==0)
		{
			printf("Exception: please specify the read files.\n");
			return FAILED;
		}

		// check whether the parameter is valid
		if(pairedModePara==-1)
		{
			pairedModePara = 2;
		}else if(pairedModePara<=0 || pairedModePara>2)
		{
			printf("Exception: please specify the correct paired end mode.\n");
			return FAILED;
		}

		if(meanSizeInsertPara==0)
			standardDevPara = 0;

		if(pairedModePara==1 && readFileNumPara%2==1)
		{
			printf("Please specify the correct read files or paired mode.\n");
			return FAILED;
		}

		if((strlen(gapFillFlagPara)!=0) && (strcmp(gapFillFlagPara, "YES")!=0 && strcmp(gapFillFlagPara, "yes")!=0 && strcmp(gapFillFlagPara, "Y")!=0 && strcmp(gapFillFlagPara, "y")!=0
				&& strcmp(gapFillFlagPara, "NO")!=0 && strcmp(gapFillFlagPara, "no")!=0 && strcmp(gapFillFlagPara, "N")!=0 && strcmp(gapFillFlagPara, "n")!=0))
		{
			printf("Exception: please specify the correct gap filling flag.\n");
			return FAILED;
		}


		// begin to do the job
		if(startSRGA_SF(contigsFilePara, alignRegSizePara, minContigLenPara, pairedModePara, readFilesPara, readFileNumPara, meanSizeInsertPara, standardDevPara, gapFillFlagPara, outputPathPara, outputPrefixPara)==FAILED)
		{
			printf("line=%d, In %s(), cannot build scaffolds, error!\n", __LINE__, __func__);
			return ERROR;
		}
	}

	return SUCCESSFUL;
}

/**
 * Show the usage information.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise return FAILED.
 */
short showUsageInfo()
{
	printf("\nUsage: SRGA_SF [option]\n");

	printf("\nPROGRAM OPTIONS:\n");
	printf("    -c <FILE>          Input contigs file in fasta format. It is necessary to be\n"
		   "                       specified.\n");
	printf("    -f <FILES>         Read files. It is necessary to be specified.\n");
	printf("    -p <INT>           Paired end mode for read files.\n"
		   "                       1 - reads are paired with the first read in the first\n"
		   "                       file, and the second read in the second file;\n"
		   "                       2 - reads are paired and interleaved within a single file\n"
		   "                       (default).\n");
	printf("    -s <INT>           Region size of align region at contg ends. Default is\n"
		   "                       1000.\n");
	printf("    -m <INT>           Minimum size of contigs for scaffolding. Default is 100.\n");
	printf("    -ins_len <FLOAT>   Insert size for paired end library.\n");
	printf("    -ins_sdev <FLOAT>  Standard deviation of insert size for paired end library.\n"
		   "                       If not specified and the insert size is specified, this\n"
		   "                       value is set to 0.1*insertSize.\n");
	printf("    -g <STR>           Gap filling flag that indicates whether intra-scaffold\n"
		   "                       gaps should be filled in scaffolding.\n"
		   "                       Y(es) - fill gaps (default);\n"
		   "                       N(o) - do not fill gaps.\n");
	printf("    -o <STR>           Prefix of the output files.\n");
	printf("    -d <STR>           Output directory for the output files. Default is \"./\"\n");
	printf("    -h\n");
	printf("    -help              Show help information.\n");

	printf("\nREPORT BUGS to ydwang@hit.edu.cn, zhuxiao.hit@gmail.com\n\n");

	return SUCCESSFUL;
}
