/*
 * pergaMain.c
 *
 *  Created on: Jul 14, 2012
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/global.h"

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
 *  	(1) -k <INT>
 *  		K-mer size for assembly.
 *  	(2) -r <INT>
 *  		Read length cutoff. Reads longer than INT will be cut to INT to avoid base errors at 3' end. Default is the read length.
 *  	(3) -p <INT>
 *  		Paired end mode for read files.
 *  			0 - do not treat reads as paired (default);
 *  			1 - reads are paired with the first read in the first file, and the second read in the second file;
 *  			2 - reads are paired and are interleaved within a single file.
 *  	(4) -f <FILES>
 *  		Read files.
 *  	(5) -q <INT>
 *  		Single base quality threshold. Reads with single bases Phred quality < INT will be discarded. Default is 2.
 *  	(6) -ec <INT>
 *  		Whether operate the correction of reads.
 *  			0 - do not correct reads;
 *  			1 - correct reads (default).
 *  	(7) -ins_len <FLOAT>
 *  		Insert size for paired end library.
 *  	(8) -ins_sdev <FLOAT>
 *  		Standard deviation of insert size for paired end library. If not specified and the insert size is specified,
 *  		this value is set to 0.1*insertSize.
 *  	(9) -g <FILE>
 *  			Hash table file for assembly. It is required for the command 'assemble'.
 *  	(10) -o <STR>
 *  		Prefix of the output files.
 *  	(11) -d <STR>
 *  		Output directory for the output files. Default is "./".
 *  	(12) -m <INT>
 *  		Minimum size of the contigs to output.
 *  	(13) -h (-help)
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
	int i, returnCode;

	char outputPathPara[256];
	char outputPrefixPara[256];
	char *readFilesPara[256];
	char readFilesBuf[256][256];
	char graphFilePara[256];
	char contigsFilePara[256];
	char gapFillFlagPara[256];
	char readMatchInfoFilePara[256];
	int readFileNumPara;
	int operationModepara;
	int kmerSizePara;
	int readLenCutOffPara;
	int singleBaseQualThresPara;
	int pairedModePara;
	double meanSizeInsertPara, standardDevPara;
	int minContigLenPara;
	int contigAlignRegLenPara;

	// version information
	printf("\nPERGA  : A Paired-End Reads Guided De Novo Assembler\n");
	printf("Version: %s\n", VERSION_STR);
	printf("Release: %s\n\n", RELEASE_DATE_STR);

	if(argc==1)
	{
		if(showUsageInfo()==FAILED)
		{
			printf("line=%d, In %s(), cannot show the usage information for PERGA, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else if(argc==2)
	{
		i = 1;
		while(i<argc)
		{
			if(strcmp(argv[i], "all")==0 || strcmp(argv[i], "ht")==0 || strcmp(argv[i], "assem")==0 || strcmp(argv[i], "sf")==0 || strcmp(argv[i], "assem_sf")==0)
			{
				if(showUsageInfo()==FAILED)
				{
					printf("line=%d, In %s(), cannot show the usage information for PERGA, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}else if(strcmp(argv[i], "-h")==0 || strcmp(argv[i], "-help")==0)
			{
				if(showUsageInfo()==FAILED)
				{
					printf("line=%d, In %s(), cannot show the usage information for PERGA, error!\n", __LINE__, __func__);
					return FAILED;
				}
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

			i ++;
		}

		return FAILED;
	}
	else
	{
		// reset the parameters
		operationModepara = -1;
		outputPathPara[0] = '\0';
		outputPrefixPara[0] = '\0';
		graphFilePara[0] = '\0';
		contigsFilePara[0] = '\0';
		gapFillFlagPara[0] = '\0';
		readMatchInfoFilePara[0] = '\0';
		readFileNumPara = 0;
		kmerSizePara = 0;
		readLenCutOffPara = 0;
		singleBaseQualThresPara = 0;
		pairedModePara = 0;
		meanSizeInsertPara = standardDevPara = 0;
		minContigLenPara = 0;
		contigAlignRegLenPara = 0;


		if(strcmp(argv[1], "all")==0)
		{
			operationModepara = OPERATION_MODE_ALL;
		}else if(strcmp(argv[1], "ht")==0)
		{
			operationModepara = OPERATION_MODE_HASHTABLE;
		}else if(strcmp(argv[1], "assem")==0)
		{
			operationModepara = OPERATION_MODE_ASSEMBLE;
		}else if(strcmp(argv[1], "sf")==0)
		{
			operationModepara = OPERATION_MODE_SCAFFOLDING;
		}else if(strcmp(argv[1], "assem_sf")==0)
		{
			operationModepara = OPERATION_MODE_ASSEM_SCAF;
		}else
		{
			if(strcmp(argv[1], "-h")==0 || strcmp(argv[1], "-help")==0)
			{
				if(showUsageInfo()==FAILED)
				{
					printf("line=%d, In %s(), cannot show the usage information for globalMetrics, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}else
			{
				if(argv[1][0]=='-')
				{
					printf("%s : unknown argument\n", argv[1]);
					printf("Please use -h or -help for more information.\n");
				}else
				{
					printf("%s : invalid argument\n", argv[1]);
					printf("Please use -h or -help for more information.\n");
				}
			}

			return FAILED;
		}


		i = 2;
		while(i<argc)
		{

			if(strcmp(argv[i], "-k")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						kmerSizePara = atol(argv[i+1]);
						if(kmerSizePara<=0)
						{
							printf("Exception: please specify the correct k-mer size.\n");
							return FAILED;
						}
					}else
					{
						printf("Exception: please specify the k-mer size.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the k-mer size.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-r")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						readLenCutOffPara = atol(argv[i+1]);
						if(readLenCutOffPara<=0)
						{
							printf("Exception: please specify the correct read length cutoff.\n");
							return FAILED;
						}
					}else
					{
						printf("Exception: please specify the read length cutoff.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the read length cutoff.\n");
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
			}else if(strcmp(argv[i], "-q")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						singleBaseQualThresPara = atoi(argv[i+1]);
						if(singleBaseQualThresPara<0)
						{
							printf("Exception: please specify the correct single base quality threshold.\n");
							return FAILED;
						}
					}else
					{
						printf("Exception: please specify the single base quality threshold.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the single base quality threshold.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-cf")==0)
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
			}else if(strcmp(argv[i], "-mf")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						strcpy(readMatchInfoFilePara, argv[i+1]);
					}else
					{
						//printf("Exception: please specify the reads match information file.\n");
						//return FAILED;
					}
				}else
				{
					printf("Exception: please specify the reads match information file.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-g")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						strcpy(graphFilePara, argv[i+1]);
					}else
					{
						printf("Exception: please specify the hash table file.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the hash table file.\n");
					return FAILED;
				}
				i += 2;
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
			}
			else if(strcmp(argv[i], "-m")==0)
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
			}else if(strcmp(argv[i], "-gf")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						strcpy(gapFillFlagPara, argv[i+1]);
					}else
					{
						printf("Exception: please specify the correct contig align region size.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the contig align region size.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-s")==0)
			{
				if(i+1<argc)
				{
					if(argv[i+1][0]!='-')
					{
						contigAlignRegLenPara = atoi(argv[i+1]);
						if(contigAlignRegLenPara<=0)
						{
							printf("Exception: please specify the correct contig align region size.\n");
							return FAILED;
						}
					}else
					{
						printf("Exception: please specify the correct contig align region size.\n");
						return FAILED;
					}
				}else
				{
					printf("Exception: please specify the contig align region size.\n");
					return FAILED;
				}
				i += 2;
			}else if(strcmp(argv[i], "-h")==0 || strcmp(argv[i], "-help")==0)
			{
				if(showUsageInfo()==FAILED)
				{
					printf("line=%d, In %s(), cannot show the usage information for PERGA, error!\n", __LINE__, __func__);
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


		if(kmerSizePara>0 && readLenCutOffPara>0 && kmerSizePara>readLenCutOffPara)
		{
			printf("Exception: invalid k-mer size or read length cutoff.\n");
			return FAILED;
		}

		// check whether the parameter is valid
		if(pairedModePara==1 && readFileNumPara%2==1)
		{
			printf("Exception: please specify the correct read files for paired end library.\n");
			return FAILED;
		}

		if((operationModepara==OPERATION_MODE_ALL || operationModepara==OPERATION_MODE_HASHTABLE) && readFileNumPara==0)
		{
			printf("Exception: please specify the read files.\n");
			return FAILED;
		}

		if(operationModepara==OPERATION_MODE_ASSEMBLE && strlen(graphFilePara)==0)
		{
			printf("Exception: please specify the hash table file.\n");
			return FAILED;
		}

		if((strlen(gapFillFlagPara)!=0) && (strcmp(gapFillFlagPara, "YES")!=0 && strcmp(gapFillFlagPara, "yes")!=0 && strcmp(gapFillFlagPara, "Y")!=0 && strcmp(gapFillFlagPara, "y")!=0
				&& strcmp(gapFillFlagPara, "NO")!=0 && strcmp(gapFillFlagPara, "no")!=0 && strcmp(gapFillFlagPara, "N")!=0 && strcmp(gapFillFlagPara, "n")!=0))
		{
			printf("Exception: please specify the correct gap filling flag.\n");
			return FAILED;
		}


		// begin to do the job
		returnCode = startPERGA(operationModepara, kmerSizePara, readLenCutOffPara, pairedModePara, readFilesPara, readFileNumPara, singleBaseQualThresPara, graphFilePara, contigsFilePara, readMatchInfoFilePara, meanSizeInsertPara, standardDevPara, outputPathPara, outputPrefixPara, minContigLenPara, contigAlignRegLenPara, gapFillFlagPara);
		if(returnCode==FAILED)
		{
			printf("Please use -h or -help for more information.\n");
			return FAILED;
		}else if(returnCode==ERROR)
		{
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
	printf("\nUsage: perga [command] [options] ...\n");
	printf("\nPROGRAM COMMANDS:\n");
	printf("    ht                 construct hash table prior to assembly\n");
	printf("    assem              assembly (assemble reads to contigs)\n");
	printf("    sf                 scaffolding (link contigs to scaffolds)\n");
	printf("    all                do all the above in turn\n");
	printf("    assem_sf           do the assembly and scaffolding in turn\n");

	printf("\nPROGRAM OPTIONS:\n");
	printf("  1) ht -- hash table construction:\n");
	printf("    -k <INT>           K-mer size for assembly. Default is 21.\n");
	printf("    -r <INT>           Read length cutoff. Reads longer than INT will be cut to\n"
		   "                       INT to eliminate base errors at 3' end. Default is the\n"
		   "                       read length.\n");
	printf("    -p <INT>           Paired end mode for read files.\n"
		   "                       0 - do not treat reads as paired (default);\n"
		   "                       1 - reads are paired with the first read in the first\n"
		   "                       file, and the second read in the second file;\n"
		   "                       2 - reads are paired and interleaved within a single\n"
		   "                       file.\n");
	printf("    -f <FILES>         Read files. It is necessary to be specified for commands\n"
		   "                       'all' and 'ht'.\n");
	printf("    -q <INT>           Single base quality threshold. Reads with single base\n"
		   "                       Phred quality < INT will be discarded. Default is 2.\n");
	printf("    -o <STR>           Prefix of the output files.\n");
	printf("    -d <STR>           Output directory for the output files. Default is \"./\"\n");
	printf("    -h\n");
	printf("    -help              Show help information.\n");

	printf("  2) assem -- assembly (assemble reads to contigs):\n");
	printf("    -g <FILE>          Hash table file for assembly. It is required for the\n"
		   "                       command 'assem'.\n");
	printf("    -ins_len <FLOAT>   Insert size for paired end library.\n");
	printf("    -ins_sdev <FLOAT>  Standard deviation of insert size for paired end library.\n"
		   "                       If it is not specified but the insert size is specified,\n"
		   "                       this value will be set to 0.15*insertSize automatically.\n");
	printf("    -o <STR>           Prefix of the output files.\n");
	printf("    -d <STR>           Output directory for the output files. Default is \"./\"\n");
	printf("    -m <STR>           Minimum size of contigs to output. Default is 100 bp.\n");
	printf("    -h\n");
	printf("    -help              Show help information.\n");

	printf("  3) sf -- scaffolding (link contigs to scaffolds):\n");
	printf("    -g <FILE>          Hash table file for assembly. It is required for the\n"
		   "                       command 'sf'.\n");
	printf("    -cf <FILE>         Contigs file in fasta format. It is required for the\n"
		   "                       command 'sf'.\n");
	printf("    -mf <FILE>         Read match information file.  It is required for the\n"
		   "                       command 'sf'.\n");
	printf("    -s <INT>           Region size of align region at contg ends. Default is\n"
		   "                       2000.\n");
	printf("    -ins_len <FLOAT>   Insert size for paired end library.\n");
	printf("    -ins_sdev <FLOAT>  Standard deviation of insert size for paired end library.\n"
		   "                       If it is not specified and the insert size is specified,\n"
		   "                       this value will be set to 0.15*insertSize automatically.\n");
	printf("    -gf <STR>          Gap filling flag that indicates whether intra-scaffold\n"
		   "                       gaps should be filled in scaffolding.\n"
		   "                       y(es) - fill gaps (default);\n"
		   "                       n(o) - do not fill gaps.\n");
	printf("    -o <STR>           Prefix of the output files.\n");
	printf("    -d <STR>           Output directory for the output files. Default is \"./\"\n");
	printf("    -h\n");
	printf("    -help              Show help information.\n");

	printf("\nREPORT BUGS to ydwang@hit.edu.cn\n");

	return SUCCESSFUL;
}


