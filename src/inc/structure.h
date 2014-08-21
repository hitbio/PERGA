#ifndef STRUCTURE_H_INCLUDED
#define STRUCTURE_H_INCLUDED 1


//read
typedef struct readNode
{
	uint64_t rowReadseqInBlock: 26;		// point to the row of readseqBlock.readseqArr
	uint64_t readseqBlockID: 16;
	uint64_t nBaseNum: 9;		// unknown base count: it will be changed to 'C' automatically
	uint64_t seqlen: 10;
	uint64_t validFlag: 1;
	uint64_t successFlag: 1;
	uint64_t modified: 1;
}read_t;

// read block
typedef struct readBlockNode
{
	uint32_t blockID;
	uint32_t itemNum;
	read_t *readArr;
}readBlock_t;

// readseq block
typedef struct readseqBlockNode
{
	uint32_t blockID;
	uint32_t rowsNum;
	uint64_t *readseqArr;
}readseqBlock_t;


//readseq hash node
typedef struct readseqHashItemNode
{
	//uint64_t *readseq;
	uint32_t rowReadseqInBlock;		// point to the row of readseqBlock.readseqArr
	uint16_t readseqBlockID;
	uint16_t seqlen;
	//struct readseqHashItemNode *next;
	uint32_t nextHashItemBlockID;
	uint32_t nextItemRowHashItemBlock;
}readseqHashItem_t;

// readseq hash block
typedef struct readseqHashItemBlockNode
{
	uint32_t blockID;
	uint32_t itemNum;
	readseqHashItem_t *readseqHashItemArr;
}readseqHashItemBlock_t;

//readseq hash node
typedef struct readseqHashBucketNode
{
	uint32_t hashItemBlockID;
	uint32_t itemRowHashItemBlock;
}readseqHashBucket_t;

//readSet
typedef struct readSetNode
{
	// read blocks
	readBlock_t *readBlockArr;		// point to kmerSeqBlock array
	int16_t blocksNumRead;
	int16_t maxBlocksNumRead;
	int16_t bytesPerRead;
	int16_t maxReadLen;
	int64_t totalItemNumRead;
	int64_t totalValidItemNumRead;
	int64_t maxItemNumPerReadBlock;

	// readseq blocks
	readseqBlock_t *readseqBlockArr;		// point to kmerSeqBlock array
	int32_t blocksNumReadseq;
	int16_t maxBlocksNumReadseq;
	int16_t bytesPerEntryReadseq;
	int64_t totalItemNumReadseq;
	int64_t maxEntryNumReadseqBlock;

	// readseq hash table
	readseqHashBucket_t *readseqHashtable;
	int32_t hashTableSizeReadseq;
	//int16_t baseNumHashingReadseq;					// the base number for generating hash code  ================???????????????

	// readseq hash item blocks
	readseqHashItemBlock_t *readseqHashItemBlockArr;		// point to readseqHashBlock array
	int32_t blocksNumReadseqHashItem;
	int16_t maxBlocksNumReadseqHashItem;
	int16_t bytesPerReadseqHashItem;
	int64_t totalItemNumReadseqHashItem;
	int64_t maxItemNumPerReadseqHashItemBlock;

	// read blocks
	struct readMatchInfoBlockNode *readMatchInfoBlockArr;		// point to kmerSeqBlock array
	int16_t blocksNumReadMatchInfo;
	//int16_t maxBlocksNumReadMatchInfo;
	int16_t bytesPerReadMatchInfo;
	int64_t totalItemNumReadMatchInfo;
	int64_t totalValidItemNumReadMatchInfo;
	int64_t maxItemNumPerReadMatchInfoBlock;

}readSet_t;

//corrected read
typedef struct correctedReadNode
{
	uint64_t rid: 48;
	uint64_t seqlen: 16;
	//uint64_t validFlag: 1;
	//uint64_t successFlag: 1;
}correctedRead_t;

//==============================================
// decision function: c = sum(alpha_i*K(s_i, x)) + b
// where, K(s_i, x) = <s_i, x> * (1 + <s_i, x>)^2
//
// If c >= 0, it belongs to class 1;
// otherwise, it belongs to class 2.

// SVM sample vector node
typedef struct svmSampleVectorNode
{
	double *vectorData;
	int32_t colsNum;
}svmSampleVector_t;

// SVM support vector node
typedef struct svmSupportVectorNode
{
	double *vectorData;
}svmSupportVector_t;

//SVM scale Data
struct scaleDataNode
{
	double *scaleFactorVector;
	double *shiftVector;
	int32_t colsNum;
};

//SVM paras node
typedef struct svmModelNode
{
	svmSupportVector_t *supportVectors;
	int32_t rowsNumSupportVector, colsNumSupportVector;
	double *alphaVector;
	int32_t rowsNumAlphaVector;
	double bias;
	struct scaleDataNode *scaleData;
	char kernelFunc[20];
}svmModel_t;


//=========================================
//ridpos
typedef struct ridposnode
{
	uint64_t delsign:1;
	uint64_t pos:16;
	uint64_t reserved:1;
	uint64_t rid:46;
}ridpostype;

//kmer
typedef struct kmerNode
{
	ridpostype *ppos;
	uint32_t multiplicity;
	uint32_t arraysize;
	uint32_t kmerseqBlockID;
	uint32_t itemRowKmerseqBlock;
	//struct kmerNode *next;
	uint32_t nextKmerBlockID;
	uint32_t nextItemRowKmerBlock;
}kmertype;

// k-mer block
typedef struct kmerBlockNode
{
	uint32_t blockID;
	uint32_t itemNum;
	kmertype *kmerArr;
}kmerBlock_t;

//kmer hash node
typedef struct kmerHashBucketNode
{
	uint32_t kmerBlockID;
	uint32_t itemRowKmerBlock;
}kmerHashBucket_t;

// kmerseq block
typedef struct kmerseqBlockNode
{
	uint32_t blockID;
	uint32_t itemNum;
	uint64_t *kmerSeqArr;
}kmerseqBlock_t;

// ridpos block
typedef struct ridposBlockNode
{
	uint32_t blockID;
	uint32_t itemNum;			// unused
	ridpostype *ridposArr;
}ridposBlock_t;

//k-mer hash table
typedef struct graphNode
{
	// read set
	readSet_t *readSet;

	// k-mer blocks
	kmerHashBucket_t *kmerHashtable;
	int32_t hashTableSize;

	kmerBlock_t *kmerBlockArr;				// point to kmerBlock array
	int32_t blocksNumKmer;
	int16_t maxBlocksNumKmer;
	int16_t bytesPerKmer;
	int64_t totalItemNumKmer;
	int32_t maxItemNumPerKmerBlock;
	int32_t kmerSampleInterval;

	// kmerseq blocks
	kmerseqBlock_t *kmerSeqBlockArr;		// point to kmerSeqBlock array
	int32_t blocksNumKmerSeq;
	int16_t maxBlocksNumKmerSeq;
	int8_t bytesPerKmerseq;
	int8_t entriesPerKmer;
	int64_t totalItemNumKmerSeq;
	int64_t maxItemNumPerKmerSeqBlock;

	// ridpos blocks
	ridposBlock_t *ridposBlockArr;		// point to ridposBlock array
	int16_t blocksNumRidpos;
	int16_t maxBlocksNumRidpos;
	int16_t bytesPerRidpos;
	int64_t totalItemNumRidpos;
	int64_t maxItemNumPerRidposBlock;

}graphtype;


// contigGraph node, added 2013-10-28
typedef struct contigGraphItemNode
{
	int32_t contigID, localContigID, contigLen, totalAlignedReadsNum;
	char *contigSeq, *contigTitle;

	uint32_t onlyEnd5: 1;		// whether only 5' end because of short length (<=300 bp):  YES (1); NO (0).
	uint32_t shortFlag: 1;		// short means : len < 2 * alignRegLen
	uint32_t used: 1;			// the contig is used: 0-- unused, 1-- used.

	uint32_t usedTimeEnd5;			// the used time of 5' end
	uint32_t usedTimeEnd3;			// the used time of 3' end
	uint16_t maxOverlapLenEnd5: 14;
	uint16_t overlapKindEnd5: 2;	// 0-- unused; 1-- overlaped with other contigs; 2-- gapped with other contigs.
	uint16_t maxOverlapLenEnd3: 14;
	uint16_t overlapKindEnd3: 2;	// 0-- unused; 1-- overlaped with other contigs; 2-- gapped with other contigs.

	// for scaffolding, linking
	int32_t validFlag, alignRegSizeEnd5, alignRegSizeEnd3;
	int8_t delReadFlagEnd5, delReadFlagEnd3;
	double averCovNumEnd5, averCovNumEnd3;
	struct contigReadNode *contigReadArrayEnd5, *contigReadArrayEnd3;
	int32_t contigReadNumEnd5, contigReadNumEnd3;
	struct contigEdgeNode *contigEdgeArrayEnd5, *contigEdgeArrayEnd3;
	int32_t itemNumContigEdgeArrayEnd5, contigEdgeArraySizeEnd5, itemNumContigEdgeArrayEnd3, contigEdgeArraySizeEnd3;

	// the overlaps

}contigGraphItem_t;

// contigGraph structure, added 2013-10-28
typedef struct contigGraphNode
{
	contigGraphItem_t *contigItemArray;
	int32_t itemNumContigItemArray, maxItemNumContigItemArray;
	int32_t averLinkNum;
	double averCovNum;
}contigGraph_t;

// contigEdge node
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
	uint32_t totalSituationNum: 4;
	uint32_t pairedNum_gapEstimate: 26;
	uint32_t sideFlag: 2;	 // 0: unused; 1: left; 2: right.
	double gapSize;

}contigEdge_t;


//contig node
typedef struct contignode
{
	int32_t index;
	uint8_t base;
	uint16_t ridposnum;
	uint8_t reserved;
	struct successReadNode *pridposorientation;

	// added 2013-10-14
	int8_t naviFlag, naviTandFlag;    // naviFlag: 1-PE, 2-SE, 3-Mix; naviTandFlag: 1-SUCCESS, 0-FAILED, -1-UNUSED
	int8_t newCandBaseNumAfterTandPathPE, newCandBaseNumAfterTandPathSE;
	int8_t bothNullUnequalFlag; // whether the contigPath Base equal to kmer when the maxContigPathItem and secContigPathItem has no reads
	int8_t sameBaseMaxSecContigPathItem, itemNumContigPath;
	int32_t naviSuccessSize;
	int32_t occNumPE[4], occNumSE[4];
	int8_t occIndexPE[4], occIndexSE[4];
}contigtype;

//successRead Node
typedef struct successReadNode
{
	uint64_t rid: 44;
	uint64_t pos: 14;
	uint64_t orientation: 6;
	uint16_t matchnum;
	uint16_t seqlen;
	uint32_t hangingIndex;

	uint64_t *readseq;
	uint16_t matchlen;
	uint8_t entriesNumReadseq;
	uint8_t baseNumLastEentryReadseq;
	uint16_t kmerNumEnd5;
	uint16_t kmerNumEnd3;
}successRead_t;


// readMatchInfo
typedef struct readMatchInfoNode
{
	int64_t contigID: 23;
	int64_t contigPos: 25;			// the contig position
	int64_t matchlen: 16;
	int64_t readID: 42;
	int64_t seqlen: 16;
	int64_t readOrientation: 2;
	int64_t contigEnd: 4;
}readMatchInfo_t;

// read block
typedef struct readMatchInfoBlockNode
{
	uint32_t blockID;
	uint32_t itemNum;
	struct readMatchInfoNode *readMatchInfoArr;
}readMatchInfoBlock_t;

// contigReadNode
typedef struct contigReadNode
{
	int64_t readID;				// read ID
	int32_t contigPos;			// the contig position
	int16_t seqlen;
	int8_t orientation;			// the read orientation
	int8_t contigEnd;			// -1 -- unaligned; 0 -- 5' end; 1 -- 3' end, 2 -- only 5' end because of short contigs (<=300 bp)
}contigRead_t;


//assemblingread
typedef struct assemblingreadnode
{
	uint64_t rid:48;
	uint64_t firstBasePos:16;
	int32_t firstContigPos;
	uint16_t orientation:8;
	uint16_t status: 2;
	uint16_t delsign: 1;
	uint16_t locked: 1;
	uint16_t reserved: 1;
	uint16_t matedFlag: 1;

	uint16_t successiveAppearBases, successiveUnappearBases, unappearBlocksNum;

	uint16_t matchBaseNum;
	uint8_t unmatchBaseNum, alignNum;
	uint16_t basePos, lastMatchedBasePos;		// from 0: [0, ..., seqlen-1]
	uint64_t *readseq;
	struct contigPathItemNode *contigPathItem;
	struct contigPathItemReadNode *contigPathItemRead;
	int16_t mismatchNumWithContigPath;

	uint16_t seqlen;
	uint16_t entriesNumReadseq;
	uint16_t baseNumLastEentryReadseq;
	uint8_t kmerNumEnd5, kmerNumEnd3;
}assemblingreadtype;

typedef struct multiCopyReadErrNumNode
{
	int64_t readID;
	char *hangingSeq;
	int16_t seqLen;
	int16_t errBaseNum;
	int32_t rowInDecisionTable;
}multiCopyReadErrNum_t;


//+++++++++++++++++++++++++++++++++  // added 2012-12-02
typedef struct dtRowIndexNode
{
	int64_t rid:44;
	int64_t dtRow:20;
	struct dtRowIndexNode *next;
}dtRowIndex_t;


//++++++++++++++++++++++++++++++++++++
typedef struct PEReadNode
{
	int64_t rid;
	int32_t cpos;
	int16_t orient;
	int16_t seqlen;
	struct PEReadNode *next;
}PERead_t;

typedef struct readBufNode{
	char *seq;
	char *qual;
	int len;
}readBuf_t;

// temporary match information of read
typedef struct readPosTempNode
{
	uint64_t readID: 48;
	uint64_t seqlen: 16;
	uint32_t contigPos;
	uint16_t orientation;
	uint16_t matchBaseNum;
}readPosTemp_t;

//======================== structures for read list (RL) ===================
//read list Index
typedef struct readListNode
{
	uint64_t readID: 48;
	uint64_t seqlen: 16;
	uint64_t matchNum: 16;
	uint64_t firstRow: 48;
}readList_t;

//read position array
typedef struct readPosNode
{
	uint32_t contigPos;			// the contig position
	uint16_t orientation;		// the read orientation
	uint16_t matchBaseNum;		// the matched base number
}readPos_t;


// candPathNode
typedef struct candPathNode
{
	struct candPathItemNode *candPathItemArray;
	int32_t maxPathLen, itemNumCandPathItemArray, maxItemNumCandPathItemArray, maxCandPathItemLen;
}candPath_t;

// candPathItemNode
typedef struct candPathItemNode
{
	char *candPathStr;
	int32_t pathLen, supportReadsNum;

}candPathItem_t;


// dtReadInTandemPathNode
typedef struct dtReadInTandemPathNode
{
	int32_t readseqLen;
	int32_t dtRow;
	int32_t contigtailReadPos;		// the position of the contigtail node in the read base sequence
	int32_t startRowInTandemPath;
	int32_t startContigPos, endContigPos;
	//int32_t newFirstReadPos;
	struct dtReadInTandemPathNode *next;
}dtReadInTandemPath_t;

// tandemPathNode
typedef struct tandemPathItemNode
{
	char tandemPathStr[2*MAX_READ_LEN_IN_BUF+1];
	int32_t tandemPathLen;
	int32_t contigtailPathPos;		// the position of the contigtail node in the path base sequence
	int32_t startContigPos, endContigPos, contigtailPos;
	int8_t shortFragFlag, validFlag, shiftMergeFlag;
	double fragSize, fragSizeDif;
	int32_t matchWithContigFlag, mismatchNum;
	int32_t itemNumDtReadsList;
	dtReadInTandemPath_t *dtReadsList;
	struct tandemPathItemNode *next;
}tandemPathItem_t;

//// tandemPathNode
//typedef struct tandemPathNode
//{
//	tandemPathItem_t *tandPahtItemList;
//	int32_t itemNumTandPahItemList, bestItemNumTandPahtItemList;
//
//}tandemPath_t;


// contigPathNode
typedef struct contigPathNode
{
	int32_t maxContigPathLen, maxContigPathLenThres;
	int32_t startRowNewBase;		// the start position of the contigtail node in the path base sequence
	int32_t updateInterval, updateIntervalThres;
	float mismatchFactor;
	int32_t maxMismatchNumThres, overlapWithContigThres, naviSuccessSize, preNaviSuccessSize, preNaviOverlapSize;
	struct contigPathItemNode *contigPathItemList, *tailPathItem, *maxPathItem, *secPathItem, *naviPathItem;
	int32_t itemNumPathItemList, bestItemNumPathItemList;
	int32_t sameBaseMaxSecContigPathItem;

	char candPathseqTandPathPE[MAX_READ_LEN_IN_BUF+1], candPathseqTandPathSE[MAX_READ_LEN_IN_BUF+1];
	int32_t candPathLenTandPathPE, candPathLenTandPathSE;
	int32_t startRowCandPathTandPathPE, startRowCandPathTandPathSE, appearTimesCandTandPathPE, appearTimesCandTandPathSE;
	int32_t validCandPathTandPathFlagPE, validCandPathTandPathFlagSE;

	// variables for contig tail sequence
	char *contigtailSeq;
	int32_t contigtailSeqLen, maxContigtailSeqLen;
	int32_t defaultInitContigtailSeqLen;		// default: readLen
}contigPath_t;

// contigPathItemNode
typedef struct contigPathItemNode
{
	char contigPathStr[MAX_READ_LEN_IN_BUF+1];
	int32_t contigPathLen, supportReadsNum, maxOverlapWithContig;
	struct contigPathItemReadNode *pathItemReadList, *tailPathItemRead;
	struct contigPathItemNode *prePathItem, *nextPathItem;
}contigPathItem_t;

// contigPathItemReadNode
typedef struct contigPathItemReadNode
{
	int64_t rid;
	struct contigPathItemReadNode *prePathItemRead, *nextPathItemRead;
}contigPathItemRead_t;

// contigPathItemSortNode
typedef struct contigPathItemSortNode
{
	struct contigPathItemNode *pathItem;
	int32_t supportReadsNum;
}contigPathItemSort_t;


//====================== structures for scaffolds =====================
//scafContigposNode
typedef struct scafContigposNode
{
	int32_t contigID;
	int32_t contigpos;
	int32_t contigEndFlag;
}scafContigpos_t;

//scafKmer
typedef struct scafKmerNode
{
	scafContigpos_t *ppos;
	uint32_t multiplicity;
	uint32_t arraysize;
	uint32_t kmerseqBlockID;
	uint32_t itemRowKmerseqBlock;
	uint32_t nextKmerBlockID;
	uint32_t nextItemRowKmerBlock;
}scafKmer_t;

// scafKmer block
typedef struct scafKmerBlockNode
{
	uint32_t blockID;
	uint32_t itemNum;
	scafKmer_t *kmerArr;
}scafKmerBlock_t;

// contigpos block
typedef struct scafContigposBlockNode
{
	uint32_t blockID;
	uint32_t itemNum;			// unused
	scafContigpos_t *contigposArr;
}scafContigposBlock_t;

// scafContigIndex
typedef struct scafContigIndexNode
{
	// k-mer blocks
	kmerHashBucket_t *kmerHashtable;
	int32_t hashTableSize;

	scafKmerBlock_t *kmerBlockArr;				// point to scafKmerBlock array
	int32_t blocksNumKmer;
	int16_t maxBlocksNumKmer;
	int16_t bytesPerKmer;
	int64_t totalItemNumKmer;
	int32_t maxItemNumPerKmerBlock;

	// kmerseq blocks
	kmerseqBlock_t *kmerSeqBlockArr;		// point to kmerSeqBlock array
	int32_t blocksNumKmerSeq;
	int16_t maxBlocksNumKmerSeq;
	int8_t bytesPerKmerseq;
	int8_t entriesPerKmer;
	int64_t totalItemNumKmerSeq;
	int64_t maxItemNumPerKmerSeqBlock;

	// ridpos blocks
	scafContigposBlock_t *contigposBlockArr;		// point to contigposBlock array
	int16_t blocksNumContigpos;
	int16_t maxBlocksNumContigpos;
	int16_t bytesPerContigpos;
	int64_t totalItemNumContigpos;
	int64_t maxItemNumPerContigposBlock;


	//scafKmer_t *pKmerAdding;
	uint64_t *pKmerSeqAdding;
	//scafKmerBlock_t *pKmerBlockAdding;

}scafContigIndex_t;



// maxRowCol Node
typedef struct maxRowColNode
{
	int32_t maxRow, maxCol, secondMaxRow, secondMaxCol;
	int32_t maxValue, secondMaxValue;
	int16_t maxArrIndex, secondMaxArrIndex;

	int32_t contigID1, contigID2;
	int8_t endFlag1, endFlag2, orient1, orient2;

	double gapSize;
}maxRowCol_t;

// contig link information
typedef struct contigLinkNode
{
	struct contigLinkItemNode *contigLinkItemArray;		// the contig link array
	int32_t itemNumContigLinkItemArray;			// item number of contig link array
	int32_t maxItemNumContigLinkItemArray;		// maximal item number of contig link array
	int32_t headRowContigLinkItemArray;			// head row of contig link array
	int32_t tailRowContigLinkItemArray;			// tail row of contig link array
}contigLink_t;

// contig link item information
typedef struct contigLinkItemNode
{
	uint32_t contigID: 30;
	uint32_t orientation: 2;	// orientation of the contig
	int32_t contigLen;
	int16_t pairNumWithPreLink, pairNumWithNextLink;
	int32_t previous, next;
}contigLinkItem_t;


/* contigOverlap information */
typedef struct contigOverlapNode
{
	uint32_t contigID1: 30;
	uint32_t orientation1: 2;
	uint32_t contigID2: 29;		// 0 indicates the contig 2 does not exist.
	uint32_t orientation2: 1;
	uint32_t mergeFlag: 1;		// YES or NO
	uint32_t breakFlag: 1;		// YES or NO
	int16_t overlapLen;			// the overlap length > 0 if mergeFlag is YES.
	int16_t gapSize;			// the gap size >= 0 if mergeFlag is NO
	int32_t pairNum;
}contigOverlap_t;

// scaffoldItemNode
typedef struct scaffoldItemNode
{
	int32_t scaffoldID, scaffoldLen;		// the ID of the scaffold
	int32_t linkedContigsNum;			// the number of linked contigs in the scaffold
	contigOverlap_t *contigOverlapArray;
	int32_t itemNumContigOverlapArray;
	struct scaffoldItemNode *previous, *next;
}scaffoldItem_t;

// scaffoldSetNode
typedef struct scaffoldSetNode
{
	scaffoldItem_t *scaffoldItemList, *tailScaffoldItem;
	int32_t scaffoldNum;
}scaffoldSet_t;


#endif // STRUCTURE_H_INCLUDED
