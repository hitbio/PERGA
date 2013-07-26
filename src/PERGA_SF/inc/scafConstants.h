/*
 * scafConstants.h
 *
 *  Created on: Mar 31, 2011
 *      Author: zhuxiao
 */

#ifndef SCAF_CONSTANTS_H_
#define SCAF_CONSTANTS_H_ 1

#define __USE_FILE_OFFSET64
#define __USE_LARGEFILE64

#ifndef NULL
#define NULL ((void *)0)
#endif


#define SUCCESSFUL					0
#define FAILED						1
#define ERROR						-1

#define YES							1
#define NO							0

#define DEBUG_FLAG					(NO)		// YES / NO
#define DEBUG_OUT_FLAG				(NO)		// YES / NO

#define DELETE_FILES				(YES)		// YES / NO: delete the files when process finished?

/* the Strand flag */
#define ORIENTATION_PLUS			0
#define ORIENTATION_MINUS			1

#define FIRST_ROUND_ASSEMBLY			0		// the first assembly round
#define SECOND_ROUND_ASSEMBLY			1		// the second assembly round

#define CONTIG_END_LEN				1000  //===================================
#define CONTIG_END_LEN_FACTOR		10
#define MIN_CONTIG_LEN_THRES		100

#define SORT_KMER_SIZE				11
#define SORT_KMER_ARR_SIZE			(1 << (SORT_KMER_SIZE<<1))
//#define KMER_MASK					((1 << (KMER_SIZE<<1))-1)

#define MAX_READ_BUF_SIZE			1000000
#define MAX_READ_LEN				36
//#define MAX_READ_LEN				50
#define MAX_READ_NUM_READ_LEN_SAMPLE	10000

#define FIRST_LINK_ROUND			1
#define SECOND_LINK_ROUND			2

#define MAX_SHORT_LEN_THRES				500  // added 2012-11-21

#define MIN_FIRST_LINKNUM_THRES				3
#define SECOND_FIRST_SCORE_RATIO_LINKING	0.7f
#define SECOND_LINKNUM_FACTOR				3  //=================================
#define BREAK_LINKNUM_FACTOR				5

//#define MIN_FIRST_LINKNUM_FACTOR			0.2f
#define MIN_FIRST_LINKNUM_FACTOR			0.1f
//#define MIN_FIRST_LINKNUM_THRES				40
#define MAX_FIRST_LINKNUM_FACTOR			13
//#define MAX_SECOND_LINKNUM_THRES			15
//#define MAX_SECOND_LINKNUM_THRES			10  //=================================
//#define MAX_SECOND_LINKNUM_THRES			5  // 2012-11-19
#define MAX_SECOND_LINKNUM_THRES			2  // 2012-11-22
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
#define GAP_SIZE_SDEV_FACTOR_OVERLAP		0.8f
#define GAP_SIZE_SDEV_FACTOR_GAPFILLING		4.0f
//#define EXACT_OVERLAP_SDEV_THRESHOLD		0.5f

#define MATCH_SCORE			1
#define MISMATCH_SCORE		-1
#define GAP_SCORE			-1

//======================= gap filling =======================
#define TABLE_SIZE_SCAFASSEMBLINGREAD			500000
#define TABLE_SIZE_SUCCESSFUL_READ_ARRAY_SCAF	500000
#define MAX_DECISION_TABLE_SIZE_HTRES			500000

//#define ERROR_REGION_LEN_3End					2	//erroneous region length of 3' end of a read
//#define ERROR_REGION_LEN_3End					3
#define ERROR_REGION_LEN_3End_FACTOR			0.1f

#define MIN_OVERLAP_LEN							7

// 3 statuses of a read in decision table
#define ASSEMBLING_STATUS						1	// assembling
#define SUCCESSFUL_STATUS						2	// successfully assembled
#define FAILED_STATUS							3	// failed assembled


//++++++++++++++++++++++++++++++++++++
#define LONG_KMER_SIZE_FACTOR			0.7f

#define VALID_OCC_RATIO	0.05

#define MIN_KMER_OCC_FACTOR				0.15f
#define MIN_KMER_OCC_THRES				2

//#define SECOND_FIRST_SECORE_RATIO		0.7
#define SECOND_FIRST_OCC_RATIO			0.8

#define MAX_SECOND_OCC_FACTOR			4.0f
#define MAX_SECOND_OCC_THRES			6

#define MAX_FIRST_OCC_FACTOR			20

#define LONG_KMER_OCC_FACTOR			1.5f
#define MIN_LONG_KMER_OCC_THRES			6

//#define OCCS_NUM_SE_FAILED_PE_FACTOR			5.0f
#define OCCS_NUM_SE_FAILED_PE_FACTOR			6.0f
//#define MAX_OCC_NUM_FAILED_PE_THRES				60.0f
#define MAX_OCC_NUM_FAILED_PE_THRES				100.0f

#define MAX_NAVI_NUM_SE_THRES					2  // === not used ===


#define FILE_FORMAT_FASTA				1
#define FILE_FORMAT_FASTQ				2

// PE given types: 0--singeEnd, 1--none, 2--only insert, 3-both insert and sdev.
#define SE_GIVEN_TYPE					0
#define NONE_PE_GIVEN_TYPE				1
#define INSERT_PE_GIVEN_TYPE			2
#define BOTH_PE_GIVEN_TYPE				3

#define SDEV_FACTOR						3.0f  //===================================
//#define SDEV_FACTOR						5.0f

#define DRAFT_SDEV_FACTOR				0.15f

#endif /* SCAF_CONSTANTS_H_ */
