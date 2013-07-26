/*
 * scafStructure.h
 *
 *  Created on: Mar 31, 2011
 *      Author: zhuxiao
 */

#ifndef SCAFFOLDING_STRUCTURE_H_
#define SCAFFOLDING_STRUCTURE_H_ 1

#include <stdint.h>


//=================== structures for contig index (CI) ===================
// k-mer node for contig index (CI)
typedef struct uniqueSeqKmerNode
{
	unsigned int itemRow;
	unsigned int itemNum;
}uniqueSeqKmer;


// base sequence index
typedef struct baseSeqNode
{
	uint32_t contigID;			// starts from 1
	uint32_t contigPos: 30;		// starts from 1
	uint32_t contigEnd: 2;		// 0 -- 5' end; 1 -- 3' end, 2 -- only 5' end because of short contigs (<=300 bp)
	uint32_t startRow;			// the start row in sequence array
}baseSeqIndex;

// node for sequence row array
typedef struct seqRowNode
{
	uint32_t rowsNum;		// the number of rows
	uint32_t startRow;		// the start row in the contigInfo array
}seqRow;

// node for contigMatchInfo array
typedef struct contigmatchInfoNode
{
	uint32_t contigID;			// starts from 1
	uint32_t contigPos: 30;		// starts from 1
	uint32_t contigEnd: 2;		// 0 -- 5' end; 1 -- 3' end, 2 -- only 5' end because of short contigs (<=300 bp)
}contigMatchInfo;

// sort k-mer node for sorting
typedef struct sortKmernode
{
	unsigned int startSortRow;
	unsigned int num;
	unsigned int arraysize;
}sortkmertype;



//====================structure for matching =================
typedef struct readBufNode{
	char *seq;
	//char *qual;
	int len;
}readBuf_t;


// temporary match information of read
typedef struct readMatchInfoTempNode
{
	uint64_t readID;
	uint32_t contigID: 30;
	uint32_t orientation: 2;
	uint32_t contigPos: 30;
	uint32_t contigEnd: 2;		// 0 -- 5' end; 1 -- 3' end, 2 -- only 5' end because of short contigs (<=300 bp)
}readMatchInfoTemp;



//======================== structures for read list (RL) ===================
//read list Index
typedef struct ReadListNode
{
	uint64_t readID;			// read ID
	uint32_t matchNum;			// the match number of the read
	uint32_t curNum;			// the current match number of the read
	uint32_t firstRow;			// the firstRow of the read in the read position array
}ReadList;

//read position array
typedef struct ReadPosNode
{
	uint32_t contigID: 30;		// contig ID
	uint32_t orientation: 2;	// the read orientation
	uint32_t contigPos: 30;		// the contig position
	uint32_t contigEnd: 2;		// 0 -- 5' end; 1 -- 3' end, 2 -- only 5' end because of short contigs (<=300 bp)
}ReadPos;


//====================== structures for contig list (CL) =====================
/* contig list index */
typedef struct ContigListNode
{
	uint32_t contigID;
	uint32_t EndNum5;
	uint32_t curNum5;
	uint32_t firstRow5;
	uint32_t EndNum3;
	uint32_t curNum3;
	uint32_t firstRow3;
}ContigList;

/* read position node */
typedef struct ContigReadNode
{
	uint64_t readID;			// read ID
	uint32_t contigPos;			// the contig position
	uint16_t orientation;		// the read orientation
	uint16_t contigEnd;			// 0 -- 5' end; 1 -- 3' end, 2 -- only 5' end because of short contigs (<=300 bp)
}ContigRead;



//====================== structures for linking =====================
// arrLoc node
typedef struct arrLocNode
{
	int32_t row;
	int32_t col;
	//int32_t value;
}arrLoc;

// ContigEdge node
typedef struct contigEdgeNode
{
	//int32_t row;
	int32_t col: 30;
	int32_t used: 2;
	//int32_t value;						// need to be updated

	//========================================
//	uint32_t contigID1: 29;
//	uint32_t endFlag1: 2;
//	uint32_t orient1: 1;
//	uint32_t contigID2: 29;
//	uint32_t endFlag2: 2;
//	uint32_t orient2: 1;
	int32_t pairedNum;
	int32_t validNums[4];
	int32_t maxIndexes[4];
	int32_t totalSituationNum;
	double gapSize;

	//struct contigGraphNode *next;
}ContigEdge;

// ContigGraph node
typedef struct contigGraphNode
{
	int32_t contigsNum;
	int32_t arraySize;
	struct contigEdgeNode *pEdgeArray;
}ContigGraph;

/* general contig information*/
typedef struct contigInfoNode
{
	uint32_t contigID;
	uint32_t contigLen: 29;		// the length of contig
	uint32_t onlyEnd5: 1;		// whether only 5' end because of short length (<=300 bp):  YES (1); NO (0).
	uint32_t shortFlag: 1;		// short means : len < 2 * alignRegLen
	uint32_t used: 1;			// the contig is used: 0-- unused, 1-- used.
	uint32_t usedTimeEnd5;			// the used time of 5' end
	uint32_t usedTimeEnd3;			// the used time of 3' end
	uint16_t maxOverlapLenEnd5: 14;
	uint16_t overlapKindEnd5: 2;	// 0-- unused; 1-- overlaped with other contigs; 2-- gapped with other contigs.
	uint16_t maxOverlapLenEnd3: 14;
	uint16_t overlapKindEnd3: 2;	// 0-- unused; 1-- overlaped with other contigs; 2-- gapped with other contigs.
	char *headTitle;
	char *contigSeq;
}contigInfo;

// maxRowCol Node
typedef struct maxRowColNode
{
	int32_t maxRow;
	int32_t maxCol;
	int32_t maxValue;
	int32_t maxArrIndex;				// the index of maximal value of validNums[4]
	int32_t secondMaxRow;
	int32_t secondMaxCol;
	int32_t secondMaxValue;
	int32_t secondMaxArrIndex;			// the index of second maximal value of validNums[4]

	int32_t contigID1;
	int32_t contigID2;
	int8_t endFlag1;
	int8_t endFlag2;
	int8_t orient1;
	int8_t orient2;
}maxRowCol;


/* contig link information */
typedef struct contigLinkNode
{
	uint32_t contigID: 30;
	uint32_t orientation: 2;	// orientation of the contig
	int32_t contigLen;
	int32_t previous;
	int32_t next;
}contigLink;

// linkCandidateNode
typedef struct linkCandidateNode
{
	uint32_t contigID1: 29;
	uint32_t endFlag1: 2;
	uint32_t orient1: 1;
	uint32_t contigID2: 29;
	uint32_t endFlag2: 2;
	uint32_t orient2: 1;
	int32_t pairedNum;
	int32_t validNums[4];
	int32_t maxIndexes[4];
	int32_t totalSituationNum;
	double gapSize;
}linkCandidate;


//====================== structures for contig overlap =====================
/* contigOverlap index */
typedef struct contigOverlapIndexNode
{
	uint32_t scaffoldID;		// the ID of the scaffold
	uint32_t linkedNum;			// the number of linked contigs in the scaffold
	uint32_t startRow;			// the start row in contigOverlap array of the scaffold
	uint32_t rowsNum;			// the number of rows
}contigOverlapIndex;

/* contigOverlap information */
typedef struct contigOverlapNode
{
	uint32_t contigID1: 30;
	uint32_t orientation1: 2;
	uint32_t contigID2: 29;		// 0 indicates the contig 2 does not exist.
	uint32_t orientation2: 1;
	uint32_t mergeFlag: 1;		// YES or NO
	uint32_t breakFlag: 1;		// YES or NO
	int32_t overlapLen;			// the overlap length > 0 if mergeFlag is YES.
	int32_t gapSize;			// the gap size >= 0 if mergeFlag is NO
}contigOverlap;


//====================== structures for gap filling =====================
/*scafKmer in scaffolding */
typedef struct scafKmerNode
{
	struct scafRidposNode *ppos;		// ָpointer to scafRidpos table
	uint32_t multiplicity;				// scafKmer multiplicity
	uint32_t arraysize;					// the item number of the scafRidpos table
	uint64_t *scafKmerseq;				// ========= new added ========
	struct scafKmerNode *next;			// ========= new added ========
}scafKmer;

/*de Bruijn graph in scaffolding */
typedef struct scafGraphNode
{
	struct scafKmerNode **pkmers;		// pointer to the scafKmer array
}scafGraph;


/* scafRidpos in scaffolding */
typedef struct scafRidposNode
{
	uint64_t rid:46;					// read ID
	uint64_t pos:16;					// the scafKmer position in a read
	uint64_t used:1;					// the used flag of a scafKmer in a read
	uint64_t reserved:1;				// the reserved bit
}scafRidpos;


/* scafContigNode in scaffolding */
typedef struct scafContigNode
{
	uint32_t index;										// contig ID
	uint32_t base: 8;									// contig base: 0, 1, 2, 3
	uint32_t ridposnum: 23;								// item number of its scafRidposOrient array
	uint32_t reserved: 1;								// reserved bit
	struct scafSuccessReadNode *pridposorientation;		// pointer to its scafRidposOrient array
	struct scafContigNode *next;						// pointer to next scafContigNode
}scafContig;

/* scafSuccessReadNode in scaffolding */
typedef struct scafSuccessReadNode
{
	uint64_t rid: 62;				// read ID
	uint64_t orientation: 2;		// read orientation: ORIENTATION_PLUS, ORIENTATION_MINUS
	uint32_t startmatchpos;			// start match position of a read
	uint32_t matchnum;				// number of matched bases of a read
}scafSuccessRead;

/* scafAssemblingReadNode in scaffolding */
typedef struct scafAssemblingReadNode
{
	uint64_t rid: 49;					// read ID

	uint64_t status: 2;					// Status: ASSEMBLING_STATUS, SUCCESSFUL_STATUS, FAILED_STATUS
	uint64_t kmerunappearblocks: 8;		// scafKmer continuous disappeared blocks number
	uint64_t delsign: 1;				// delete flag
	uint64_t locked: 1;					// lock flag
	uint64_t matedFlag: 1;				// the mated read
	uint64_t reserved: 2;				// reserved bit

	uint32_t firstpos: 24;				// the first scafKmer position of a read participating the assembly
	uint32_t orientation: 8;			// read orientation: ORIENTATION_PLUS, ORIENTATION_MINUS
	uint32_t kmerappeartimes: 24;		// appeared times of a scafKmer
	uint32_t kmerunappeartimes: 8;		// disappeared times of scafKmer
	uint32_t latestpos;					// the latest appeared scafKmer position of a read
	uint32_t lastpos;					// the last scafKmer position

}scafAssemblingRead;



#endif /* SCAFFOLDING_STRUCTURE_H_ */
