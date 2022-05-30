#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include "uthash.h"
#define _LOC_LEN_ 100
#ifndef START
	#define START 4
#endif

#ifndef MAX_READ_LENGTH
        #define MAX_READ_LENGTH 400
#endif

#ifndef MAX_READ_LENGTH_
        #define MAX_READ_LENGTH_ MAX_READ_LENGTH+2
#endif

#ifndef FASLE
	#define FASLE 0
#endif
#ifndef TRUE
	#define TRUE 1
#endif
#ifndef NUM_GERM
	#define NUM_GERM 10000
#endif

typedef int STATUS;

typedef struct _NODE{
	int data;
	char seque[2*MAX_READ_LENGTH];
//	char *seque;
	struct _NODE* next;
}NODE;

typedef struct _HASH_TABLE{
	int count;
	NODE *value[NUM_GERM];
}HASH_TABLE;

typedef unsigned long long bit64_t;

struct Merged_II{

	char rEAd[MAX_READ_LENGTH*2];
	char qUANl[MAX_READ_LENGTH*2];
};

typedef struct TwoTermination{
	
	bit64_t seg1;
	bit64_t seg2;
}twoterm;

typedef struct LOCLEN{
		
	short count;
	unsigned long loc[_LOC_LEN_];
	unsigned len[_LOC_LEN_];
}loclen;

typedef struct _NODE_{
	twoterm sKey;
	loclen nValue;
	UT_hash_handle hh;
}HashNode;

typedef struct BE_LOCLEN{
	
	unsigned long loc;
	unsigned len;
}best_loclen;

typedef struct ScoreLen{
	
	float score;
	short intlen1;
	short intlen2;
}Score;

struct VJmem{
	float   Match_Rate;	
	int	length;
};

struct emp_freq{

	double prA;		/* freqA / total */
	double prC;		/* freqC / total */
	double prG;		/* freqG / total */
	double prT;		/* freqT / total */

	double q;		/* prA*prA + prC*prC + prG*prG +prT*prT */
};




struct DeterPara{

	int beloc[2];		/*  beloc[0] denotes the highest number in loc_num, beloc[1] denotes the highest number overlap's location. */

	int MatchNum;		/* the number of matched bases */
	int OvlapLength;	/* the best overlap's length */
	float MatchRate;	/* for the voerlap, calculate the match_rate whose formula is 1 - match_num / ovlap_length */
	double AS;	/* as_score calculates the overlap score and the score can be used to determine whether the overlap is selected. */
	int OvStart;	// Note the start position for assemble.
	char MergedRead[MAX_READ_LENGTH*2]; // merged bases
	char MergedQual[MAX_READ_LENGTH*2]; // merged qualities
};

struct OVLPstart{

	int Rostart;	//The overlap position in Read1 to start.
	int AssStart;	//For assemle start pisition.
	int R1Len;
	int R2Len;
	int OvlapLen;
	int OvlStaInR1;
};



/* this struct type is used to record the loction of alignment */
struct locat{

	int	number;
	int	posit;
	int	over;
};
