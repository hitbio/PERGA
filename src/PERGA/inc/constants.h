#ifndef CONSTANTS_H_INCLUDED
#define CONSTANTS_H_INCLUDED 1

#define __USE_FILE_OFFSET64
#define __USE_LARGEFILE64
//#define _FILE_OFFSET_BITS 64


#define VERSION_STR						("v0.4.01.03")
#define RELEASE_DATE_STR				("May 6, 2013")

#ifndef NULL
#define NULL ((void *)0)
#endif

#define SUCCESSFUL	0
#define FAILED		1
#define ERROR		-1

#define YES		1
#define NO		0

// Debug output configurations
#define DEBUG_PARA_PRINT				(NO)
#define DEBUG_OUTPUT					(NO)
#define DEBUG_CONTIG_CHECK				(NO)
#define DEBUG_EST_OUTPUT 				(NO)
#define DEBUG_EST_CONTIG_CHECK			(NO)
#define RESERVE_EST_FILES				(NO)
#define HANGING_CONTIG_OUT_FLAG			(NO)
#define DEBUG_ALIGNMENT					(NO)
#define MALIGN_OUTPUT					(NO)

#define TRIM_READ_LEN_FLAG				(NO)
#define ERROR_CORRECTION_FLAG_DEFAULT	(NO)

#define DRAW_CURVE_FLAG					(NO)
#define USING_RULES						(NO)
#define SVM_NAVI						(YES)


/* the Strand flag */
#define ORIENTATION_PLUS				'+'
#define ORIENTATION_MINUS				'-'

#define STRAND_PLUS						0
#define STRAND_MINUS					1

#define DEFAULT_KMER_SIZE				21
#define LONG_KMER_SIZE_FACTOR			0.8f  // ========================
//#define LONG_KMER_SIZE_FACTOR			0.7f

#define HASH_TABLE_SIZE_LARGE					(1572866113LLU)  // not prime
#define HASH_TABLE_SIZE_MEDIUM					(157286611LLU)
#define HASH_TABLE_SIZE_SMALL					(15728681LLU)

#define BLOCK_SIZE_PER_READ				(1 << 26)	// 64 MB
#define MAX_BLOCKS_NUM_READ				1024

#define BLOCK_SIZE_PER_READ_SEQ			(1 << 26)	// 64 MB
#define MAX_BLOCKS_NUM_READ_SEQ			1024

#define BLOCK_SIZE_PER_READ_SEQ_HASH	(1 << 26)	// 64 MB
#define MAX_BLOCKS_NUM_READ_SEQ_HASH	1024

#define BLOCK_SIZE_PER_KMER_SEQ			(1 << 26)	// 64 MB
#define MAX_BLOCKS_NUM_KMER_SEQ			1024

#define BLOCK_SIZE_PER_KMER				(1 << 26)	// 64 MB
#define MAX_BLOCKS_NUM_KMER				1024

#define BLOCK_SIZE_PER_RIDPOS			(1 << 27)	// 128 MB
#define MAX_BLOCKS_NUM_RIDPOS			2048

#define READS_NUM_PER_FILE_SAMPLE		10000
#define UNKNOWN_BASE_REPLACE_CHAR		'C'
#define MAX_UNKNOWN_BASE_NUM			0
#define RESERVE_HASH_ITEM_READ_SET		(YES)


#define KMER_REG_LEN_RATIO_END5			0.03f
#define KMER_REG_LEN_RATIO_END3			0.08f

#define NAVI_SUCCESS					1
#define NAVI_FAILED						0
#define NAVI_UNUSED						-1


//=========================================================
#define TABLE_SIZE_ASSEMBLINGREAD		1000000
#define TABLE_SIZE_RIDPOSORIENTATION 	1000000
#define MAX_DECISION_TABLE_SIZE_HTRES	1000000
#define ITEM_NUM_CONTIG_ARR				500000

#define MAX_ITEM_NUM_CONTIGS_LEN_ARR	10000

#define ITEM_NUM_COUNTHANGING_BUCKETS	5000

#define FIRST_ROUND_ASSEMBLY			1		// the first assembly round
#define SECOND_ROUND_ASSEMBLY			2		// the second assembly round

#define FIRSTKMER_FACTOR 				5
//#define FIRSTKMER_FACTOR 				15

#define LOWER_BOUND_FACTOR_CYCLE1		0.5f
#define UPPER_BOUND_FACTOR_CYCLE1		30
#define UPPER_BOUND_FACTOR_CYCLE2		200
//#define LOWER_BOUND_FACTOR_CYCLE3		0.2
#define MIN_LOWER_BOUND					6
#define MIN_LOWER_BOUND_RESCUE			3

#define READ_SIMILARITY_THRES			0.9f

#define FIRSTKMER_SUBTRACT_THRESHOLD 	0.9f

//#define AVERAGE_QUAL_THRESHOLD_3End 	2.0f
//#define AVERAGE_QUAL_THRESHOLD_3End 	5.0f	//==========================
#define AVERAGE_QUAL_THRESHOLD_3End 	3.0f
//#define AVERAGE_QUAL_THRESHOLD_5End 	20.0f	//==========================
#define AVERAGE_QUAL_THRESHOLD_5End 	10.0f
#define SINGLE_QUAL_THRESHOLD			2
//#define SINGLE_QUAL_THRESHOLD			3

#define ARTIFACTS_BASE_A_THRESHOLD 		0.9f

#define QUAL_BASE_NUM_3End_FACTOR		0.2f
#define ERROR_REGION_LEN_3End_FACTOR	0.1f
//#define ERROR_REGION_LEN_3End_FACTOR	0.07f
//#define ERROR_REGION_LEN_3End_FACTOR	0.15f

//#define VALID_OCC_RATIO					0.05   //==================================
//#define MAX_SECOND_OCC_FACTOR			2  //==================================
#define MAX_SECOND_OCC_FACTOR			4.0f
//#define MAX_SECOND_OCC_THRES			8  //==================================
//#define MAX_SECOND_OCC_THRES			6
#define MAX_SECOND_OCC_THRES			16
//#define MAX_FIRST_OCC_FACTOR			30  //==================================
#define MAX_FIRST_OCC_FACTOR			20

//#define SECOND_FIRST_SECORE_RATIO		0.5f   //==================================
//#define SECOND_FIRST_SECORE_RATIO		0.7f
//#define SECOND_FIRST_OCC_RATIO			0.5f   //==================================
//#define SECOND_FIRST_OCC_RATIO			0.7f // best
#define SECOND_FIRST_OCC_RATIO			0.8f
//#define SECOND_FIRST_OCC_FAILED_RATIO	0.3f
//#define SECOND_FIRST_OCC_FAILED_RATIO	0.5f
#define SECOND_FIRST_OCC_FAILED_RATIO	0.7f

#define CONTIG_LEN_THRESHOLD			100
#define CONTIG_LEN_FACTOR				2


//#define MIN_OVERLAP_LEN					23   //==================================
//#define MIN_OVERLAP_LEN					11
#define MIN_OVERLAP_LEN					7

//#define MIN_KMER_OCC_FACTOR			0.2f  // best : 0.17f
//#define MIN_KMER_OCC_FACTOR			0.17f    //==================================
//#define MIN_KMER_OCC_FACTOR				0.15f
#define MIN_KMER_OCC_FACTOR				0.1f
//#define MIN_KMER_OCC_THRES				3    //==================================
#define MIN_KMER_OCC_THRES				2
#define MAX_KMER_OCC_THRES				5

//#define MIN_LONG_KMER_OCC_THRES			8  //==================================
#define MIN_LONG_KMER_OCC_THRES			6
//#define MIN_LONG_KMER_OCC_THRES			4
//#define LONG_KMER_OCC_FACTOR			2   //==================================
//#define LONG_KMER_OCC_FACTOR			2.5f
#define LONG_KMER_OCC_FACTOR			2.5f

//决策表中reads的3种状态
#define ASSEMBLING_STATUS				1
#define SUCCESSFUL_STATUS				2
#define FAILED_STATUS					3

#define NAVI_PE_FLAG					1
#define NAVI_SE_FLAG					2
#define NAVI_MIX_FLAG					3

//#define MAX_ITEM_NUM_NAVI_OCC_QUEUE		8

#define BASE_TYPE_FASTA_CONTIG_FILE				0
#define HANGING_READ_TYPE_CONTIG_FILE			1

//++++++++++++++++++++++++++++++++++++
#define MAX_READ_BUF_SIZE				10000
#define MAX_READ_LEN_IN_BUF				5000


#define FILE_FORMAT_FASTA				1
#define FILE_FORMAT_FASTQ				2

// PE given types: 0--singeEnd, 1--none, 2--only insert, 3-both insert and sdev.
#define SE_GIVEN_TYPE					0
#define NONE_PE_GIVEN_TYPE				1
#define INSERT_PE_GIVEN_TYPE			2
#define BOTH_PE_GIVEN_TYPE				3

#define RID_LOW_BITS_NUM				20
#define TABLE_SIZE_HASH_PE				(1 << RID_LOW_BITS_NUM)
#define RID_LOW_BITS_MASK				((1 << RID_LOW_BITS_NUM) - 1)

#define TABLE_SIZE_HASH_DTROWINDEX		(1 << RID_LOW_BITS_NUM)

#define MIN_READ_NUM_PE_HASH_FACTOR		0.5f

#define MAX_NUM_EST_CONTIG				100
//#define TOTAL_CONTIG_LEN_EST_THRES		500000
#define TOTAL_CONTIG_LEN_EST_THRES		200000
#define MIN_CONTIG_LEN_EST				1000
#define MIN_CONTIG_LEN_EST_FACTOR		5.0f
#define TOTAL_READS_NUM_EST_THRES		100000

#define MAX_INSERT_SIZE_THRES			5000

#define MAX_INSERT_SIZE_FACTOR			5.0f
#define SDEV_FACTOR						3.0f  //===================================
//#define SDEV_FACTOR						5.0f
#define DRAFT_SDEV_FACTOR				0.15f
#define REG_LEN_PE_HASH_FACTOR			0.6f


//++++++++++++++++++++++++++++++++++++
//#define REG_LEN_READS_NUM_REG_FACTOR			1.3f
//#define REG_LEN_READS_NUM_REG_FACTOR			2.5f
//#define REG_LEN_READS_NUM_REG_FACTOR			1.5f // =============================
//#define REG_LEN_READS_NUM_REG_FACTOR			0.5f
#define REG_LEN_READS_NUM_REG_FACTOR			1.0f
//#define REG_LEN_READS_NUM_REG_FACTOR			0.9f
//#define REG_LEN_READS_NUM_REG_FACTOR			0.8f
//#define MAX_READS_NUM_RATIO_THRES				2.0f   //===================================
//#define MAX_READS_NUM_RATIO_THRES				5.0f
#define MAX_READS_NUM_RATIO_THRES				20.0f
#define MIN_READS_NUM_RATIO_THRES				0.1f   //===================================
//#define MIN_READS_NUM_RATIO_THRES				0.2f
//#define MIN_READS_NUM_RATIO_THRES				0.05f

//#define OCCS_NUM_SE_FAILED_PE_FACTOR			1.0f
//#define OCCS_NUM_SE_FAILED_PE_FACTOR			5.0f
#define OCCS_NUM_SE_FAILED_PE_FACTOR			6.0f  //===================================
//#define OCCS_NUM_SE_FAILED_PE_FACTOR			3.0f  //??
//#define OCCS_NUM_SE_FAILED_PE_FACTOR			4.5f
//#define MAX_OCC_NUM_FAILED_PE_THRES				51.0f
//#define MAX_OCC_NUM_FAILED_PE_THRES				60.0f
#define MAX_OCC_NUM_FAILED_PE_THRES				100.0f  //============================
//#define MAX_OCC_NUM_FAILED_PE_THRES				200.0f
//#define MAX_OCC_NUM_FAILED_PE_THRES				300.0f

#define MAX_NAVI_NUM_SE_THRES					2  // === not used ===


//=============== correction.c ====================
//#define MAX_UNMATCH_BASE_NUM_PER_READ				4
#define MAX_UNMATCH_BASE_NUM_PER_READ				3  // best

#define MIN_UNMATCH_BASE_NUM_ALIGN_FACTOR			3
#define MAX_UNMATCH_BASE_NUM_AFTER_ALIGN_FACTOR		1

//#define MAX_ERR_NUM_IN_CORRECTION					10
#define MAX_ERR_NUM_IN_CORRECTION					5  // best
//#define MAX_ERR_NUM_IN_CORRECTION					4

#define MIN_SUCCESSIVE_APPEARED_BASE_NUM			8
//#define MIN_SUCCESSIVE_APPEARED_BASE_NUM			10

#define MATCH_SCORE			1
#define MISMATCH_SCORE		-1
#define GAP_SCORE			-1


//================= ref.c =======================
#define MAX_REF_POS_NEAR_THRES						1000
#define MAX_READS_MALIGN							10000


//================= malign.c =======================
#define SEC_MAX_OCC_RATIO_SVM						0.7f

#define MAX_INCOR_BASE_NUM_MALIGN					3
#define MAX_ALIGN_LEN_MALIGN						30
#define MXA_OCC_RATIO_MALIGN						0.6f
#define INCOR_RATIO_MALIGN							0.9f
//#define INCOR_RATIO_MALIGN							0.8f

#define MAX_SEQ_LEN_MALIGN							(MAX_ALIGN_LEN_MALIGN*10)


#endif // CONSTANTS_H_INCLUDED
