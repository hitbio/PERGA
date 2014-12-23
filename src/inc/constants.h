#ifndef CONSTANTS_H_INCLUDED
#define CONSTANTS_H_INCLUDED 1

#define __USE_FILE_OFFSET64
#define __USE_LARGEFILE64
//#define _FILE_OFFSET_BITS 64


#define VERSION_STR						("0.5.03.01")
#define RELEASE_DATE_STR				("Dec 23, 2014")

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
#define CANDPATH_OUTPUT					(NO)
#define TANDPATH_OUTPUT					(NO)

#define DEBUG_SCAF_FLAG					(NO)
#define DEBUG_SCAF_OUT_FLAG				(NO)
#define DEBUG_SCAF_OVERLAP_FLAG			(NO)
#define DEBUG_SCAF_FILLGAP_FLAG			(NO)

#define SVM_NAVI						(YES)


#define OPERATION_MODE_ALL				0
#define OPERATION_MODE_HASHTABLE		1
#define OPERATION_MODE_ASSEMBLE			2
#define OPERATION_MODE_SCAFFOLDING		3
#define OPERATION_MODE_ASSEM_SCAF		4


/* the Strand flag */
#define ORIENTATION_PLUS				0
#define ORIENTATION_MINUS				1

//#define DEFAULT_KMER_SIZE				21  // old
#define DEFAULT_KMER_SIZE				25
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
#define RESERVE_HASH_ITEM_READ_SET		(NO)


#define KMER_REG_LEN_RATIO_END5			0.03f
#define KMER_REG_LEN_RATIO_END3			0.08f
#define KMER_SAMPLE_INTERVAL			2

#define NAVI_SUCCESS					1
#define NAVI_FAILED						0
#define NAVI_UNUSED						-1


//=========================================================
#define rowsNumSupportVector_PE			16
#define colsNumSupportVector_PE			4
#define rowsNumSupportVector_SE			171
#define colsNumSupportVector_SE			4

#define svmKF_PE					("polynomial")
#define svmKF_SE					("polynomial")

#define svmUseScaleData_PE				YES
#define svmUseScaleData_SE				YES


//=========================================================
#define TABLE_SIZE_ASSEMBLINGREAD		1000000
#define TABLE_SIZE_RIDPOSORIENTATION 	1000000
#define ITEM_NUM_CONTIG_ARR				500000

#define MAX_ITEM_NUM_CONTIGS_LEN_ARR	10000

#define ITEM_NUM_COUNTHANGING_BUCKETS	5000

#define FIRST_ROUND_ASSEMBLY			1		// the first assembly round
#define SECOND_ROUND_ASSEMBLY			2		// the second assembly round

#define FIRSTKMER_FACTOR 				5
//#define FIRSTKMER_FACTOR 				15

#define KMER_OCC_LOWER_BOUND_FACTOR		0.5f
//#define KMER_OCC_UPPER_BOUND_FACTOR		30
#define KMER_OCC_UPPER_BOUND_FACTOR		10
#define MIN_LOWER_BOUND					6
#define MIN_LOWER_BOUND_RESCUE			3

#define READ_SIMILARITY_THRES			0.9f
//#define READ_SIMILARITY_THRES			0.95f

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

#define CONTIG_LEN_THRESHOLD			200
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

// Three status of reads in decision table
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
#define TOTAL_CONTIG_LEN_EST_THRES		20000
#define MIN_CONTIG_LEN_EST				1000
#define MIN_CONTIG_LEN_EST_FACTOR		5.0f
#define TOTAL_READS_NUM_EST_THRES		10000

#define MAX_INSERT_SIZE_THRES			5000

#define MAX_INSERT_SIZE_FACTOR			5.0f
#define SDEV_FACTOR						3.0f  //===================================
//#define SDEV_FACTOR						4.0f
#define DRAFT_SDEV_FACTOR				0.15f
#define REG_LEN_PE_HASH_FACTOR			0.6f
#define MAX_REG_LEN_PE_HASH_THRES		80


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

//#define MAX_UNMATCH_BASE_NUM_PER_READ				4
#define MAX_UNMATCH_BASE_NUM_PER_READ				3  // best

#define MIN_SUCCESSIVE_APPEARED_BASE_NUM			8
//#define MIN_SUCCESSIVE_APPEARED_BASE_NUM			10

#define SEC_MAX_OCC_RATIO_SVM						0.7f

//================= candPath.c =======================
// single-end navi length thresholds
#define MAX_NAVI_LEN_SINGLE_READ					30
//#define MAX_NAVI_LEN_WITHOUT_PE						60
#define MAX_NAVI_LEN_WITHOUT_PE						50
//#define MAX_NAVI_LEN_WITHOUT_PE						150

#define MAX_CANDIDATE_PATH_NUM						500
#define MAX_CANDIDATE_PATH_LEN						2000

//#define MAX_MISMATCH_NUM_CANDIDATE_PATH				2
#define MAX_MISMATCH_NUM_CANDIDATE_PATH				0

#define MAX_OCC_RATIO_CANDPATH						0.6f
//#define MAX_INCOR_BASE_NUM_CANDPATH					5
#define MAX_INCOR_BASE_NUM_CANDPATH					3
#define INCOR_RATIO_CANDPATH						0.9f

#define MIN_OVERLAP_SIZE_TANDPATH					10
#define MAX_INCOR_BASE_NUM_TANDPATH					3
#define MAX_OCC_RATIO_TANDPATH						0.7f
#define BEST_ITEMNUM_TANDPATH						20
#define MIN_MATCH_BASENUM_TANDPATH					40

#define MAX_UPDATE_INTERVAL_FACTOR_CONTIGPATH		0.4f
//#define MAX_UPDATE_INTERVAL_FACTOR_CONTIGPATH		0.2f
#define MAX_MISMATCHNUM_FACTOR_CONTIGPATH			0.1f
//#define MAX_MISMATCHNUM_FACTOR_CONTIGPATH			0.2f
//#define MAX_MISMATCHNUM_FACTOR_CONTIGPATH			0.05f
#define BEST_ITEM_NUM_CONTIGPATH					20
#define MAX_ITEM_NUM_CONTIGPATH						100
#define OVERLAP_SIZE_WITH_CONTIG_FACTOR				0.5f

//================= contiggraph.c =======================
#define MAX_ITEMNUM_CONTIGGRAPH_DEFAULT				5000
#define OCC_RATIO_THRES_NAVI_PATH_PE				0.8f
//#define OCC_RATIO_THRES_NAVI_PATH_PE				0.7f
#define OCC_RATIO_THRES_NAVI_PATH_SE				0.3f




// ================= scaffolding =====================
#define CONTIG_ALIGN_REG_SIZE				2000
#define CONTIG_ALIGN_REG_SIZE_FACTOR		10

//#define MAX_SHORT_LEN_THRES					500
#define MAX_SHORT_LEN_THRES					1000  // 2014-01-16

#define FIRST_LINK_ROUND			1
#define SECOND_LINK_ROUND			2

#define MIN_FIRST_LINKNUM_THRES				3
//#define SECOND_FIRST_SCORE_RATIO_LINKING	0.7f
#define SECOND_FIRST_SCORE_RATIO_LINKING	0.4f
#define SECOND_LINKNUM_FACTOR				3  //=================================
#define BREAK_LINKNUM_FACTOR				5

//#define MIN_FIRST_LINKNUM_FACTOR			0.2f
#define MIN_FIRST_LINKNUM_FACTOR			0.1f
//#define MIN_FIRST_LINKNUM_THRES				40
#define MAX_FIRST_LINKNUM_FACTOR			13
//#define MAX_SECOND_LINKNUM_THRES			15
//#define MAX_SECOND_LINKNUM_THRES			10  //=================================
//#define MAX_SECOND_LINKNUM_THRES			5  // 2012-11-19
//#define MAX_SECOND_LINKNUM_THRES			2  // 2012-11-22
#define MAX_SECOND_LINKNUM_THRES			1  // 2014-01-15
//#define MIN_BREAK_LINKNUM_THRES				15
//#define MIN_BREAK_LINKNUM_THRES				10
#define MIN_BREAK_LINKNUM_THRES				5		// 2012-11-19


#define MIN_OVERLAP_THRESHOLD				3
#define MIN_EXACT_OVERLAP_THRESHOLD			2
#define MAX_MISMATCH_THRESHOLD				3
#define MAX_MISMATCH_RATIO_THRESHOLD		0.1f
#define SUB_MISMATCH_FACTOR					2
#define MIN_ADJUST_GAP_SIZE_THRESHOLD		-10
#define MAX_ADJUST_GAP_SIZE_THRESHOLD		10
#define MIN_BASENUM_IN_GAP					2
//#define GAP_SIZE_SDEV_FACTOR_OVERLAP		0.8f
#define GAP_SIZE_SDEV_FACTOR_OVERLAP		1.0f
#define GAP_SIZE_SDEV_FACTOR_GAPFILLING		4.0f
//#define EXACT_OVERLAP_SDEV_THRESHOLD		0.5f

#define MATCH_SCORE			1
#define MISMATCH_SCORE		-2
#define GAP_SCORE			-4


#endif // CONSTANTS_H_INCLUDED
