#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Istruct.h"
#include "uthash.h"

extern int TWO_LENGTH_;
void HashInit(void);
void Hash(twoterm *sKey, long id_num, int length);
void substr(char *pre_do, int start, int length, char *dest);

bit64_t StringToULL(const char * sz)
{//	printf("r:%s\n",sz);
	bit64_t u64Result = 1;
	while(*sz != '\0')
	{
		u64Result <<= 2;
		switch(*sz)
		{
			case 'A': u64Result |= 0;       break;
			case 'T': u64Result |= 1;       break;
			case 'C': u64Result |= 2;       break;
			case 'G': u64Result |= 3;       break;
		}
		sz++;
	}
	return u64Result;
}

twoterm GetTwoTerm(char *seq)
{
	int K = strlen(seq);
	twoterm key;
	char *frount, *behind;
	frount = (char *)malloc((TWO_LENGTH_+1)*2*sizeof(char));
	behind = (char *)malloc((TWO_LENGTH_+1)*2*sizeof(char));
	substr(seq, START, TWO_LENGTH_, frount);
	substr(seq, -START, TWO_LENGTH_, behind);
	frount[TWO_LENGTH_] = '\0';
	behind[TWO_LENGTH_] = '\0';
	key.seg1 = StringToULL(frount);
	key.seg2 = StringToULL(behind);
	free(frount);frount=NULL;
	free(behind);behind=NULL;
	return key;

}

void PRODUCE_HASH(FILE *file)
{
	int K=0,i=1;
	long L=0,l=-1;
	twoterm *key;
	key = (twoterm *)malloc(sizeof(twoterm));
	char sequence[MAX_READ_LENGTH*2+1];
	HashInit();
	while(feof(file) != 1)
	{
		L = ftell(file);
		fgets(sequence,MAX_READ_LENGTH*2,file);
		sequence[strlen(sequence)-1] = '\0';
		
		//if(strlen(sequence) <= 2*(TWO_LENGTH_+START))
		//	continue;

		if(i%4 == 2)
		{
			*key = GetTwoTerm(sequence);
			K = strlen(sequence);
//			L = ftell(file);
			Hash(key,L,K); 
		}	
		if(L-l == 0)
			break;
		l = L;
		i++;
	}
	free(key);
	key = NULL;
	return;
}
