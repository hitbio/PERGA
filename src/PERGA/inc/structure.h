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
#if(DRAW_CURVE_FLAG==YES)
	//uint32_t chrSegID: 30;
	uint32_t refStrand;
	uint32_t refPos;
#endif
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
}readSet_t;

//corrected read
typedef struct correctedReadNode
{
	uint64_t rid: 48;
	uint64_t seqlen: 16;
	//uint64_t validFlag: 1;
	//uint64_t successFlag: 1;
}correctedRead_t;

//reference node
typedef struct refNode
{
	int64_t reflen: 62;
	int64_t circularFlag: 2;
	char *refseq;
}ref_t;

//==============================================
// decision function: c = sum(alpha_i*K(s_i, x)) + b
// where, K(s_i, x) = exp(-0.5 * norm(s_i, x))
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
	struct ridposnode *ppos;
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
	int64_t maxItemNumPerKmerBlock;

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


//contig node
typedef struct contignode
{
	int index;
	unsigned char base;
	unsigned short ridposnum;
	unsigned char reserved;
	struct successReadNode *pridposorientation;
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

//successRead Node
typedef struct errBaseNode
{
	uint16_t basePos;
	uint8_t newBase;
	uint8_t errKind;		// 0-- match; 1-- mismatch; 2-- deletion in read; 3-- insertion in read
}errBase_t;

//assemblingread
typedef struct assemblingreadnode
{
	uint64_t rid:48;
	uint64_t firstBasePos:16;
	int32_t firstContigPos;
	//uint16_t kmerappeartimes;
	//uint16_t kmerunappeartimes;
	//uint16_t lastappearpos;
	//uint16_t lastpos;
	uint16_t orientation:8;
	uint16_t status: 2;
	//uint16_t kmerunappearblocks: 2;
	uint16_t delsign: 1;
	uint16_t locked: 1;
	uint16_t reserved: 1;
	uint16_t matedFlag: 1;

	uint16_t successiveAppearBases, successiveUnappearBases, unappearBlocksNum;

	uint16_t matchBaseNum;
	uint8_t unmatchBaseNum, alignNum;
	uint16_t basePos, lastMatchedBasePos;		// from 0: [0, ..., seqlen-1]
	uint64_t *readseq;

	uint16_t seqlen;
	uint16_t entriesNumReadseq;
	uint16_t baseNumLastEentryReadseq;
	uint8_t kmerNumEnd5, kmerNumEnd3;

	errBase_t errBaseArr[(MIN_UNMATCH_BASE_NUM_ALIGN_FACTOR+1)*MAX_UNMATCH_BASE_NUM_PER_READ+1];
	uint8_t errBaseNum, newErrNumAfterCorrection;
	uint8_t alignSuccessTimes, alignSuccessTimesOld;			// the times of align successfully

}assemblingreadtype;


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
	int32_t orient;
	struct PEReadNode *next;
}PERead_t;

typedef struct readBufNode{
	char *seq;
	char *qual;
	int len;
#if(DRAW_CURVE_FLAG==YES)
	int headNameLen;
	char *headname;
#endif
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





#endif // STRUCTURE_H_INCLUDED
