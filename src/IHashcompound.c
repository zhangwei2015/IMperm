#include <stdlib.h>
#include <stdio.h>
#include "uthash.h"
#include "Istruct.h"

HashNode *records = NULL;
HashNode *p, *r, *tmp;

void HashInit(void)
{
//	HashNode l, *p, *r, *tmp;
	r = (HashNode*)malloc(sizeof(HashNode));
	memset(r,0,sizeof(HashNode));
	r->sKey.seg1 = 1000;
	r->sKey.seg2 = 2000;
	HASH_ADD(hh,records,sKey,sizeof(twoterm),r);
	return;
}

/*
	memset(&l,0,sizeof(HashNode));
	l.sKey.seg1 = 1000;
	l.sKey.seg2 = 2000;
*/

void Hash(twoterm *sKey, long id_num, int length)
{
	HASH_FIND(hh,records,sKey,sizeof(twoterm),p);
	if(p)
	{
		if(p->nValue.count<_LOC_LEN_)
		{
			p->nValue.loc[p->nValue.count] = id_num;
			p->nValue.len[p->nValue.count] = length;
			p->nValue.count += 1;
		}
//		printf("found %u %u %u\n",p->sKey.seg1, p->sKey.seg2,p->nValue.count);
	}
	else
	{
		r = (HashNode*)malloc(sizeof(HashNode));
		memset(r,0,sizeof(HashNode));
		r->sKey = *sKey;
		HASH_ADD(hh,records,sKey,sizeof(twoterm),r);
		HASH_FIND(hh,records,sKey,sizeof(twoterm),p);
		p->nValue.loc[p->nValue.count] = id_num;
		p->nValue.len[p->nValue.count] = length;
		p->nValue.count += 1;
	}
	p = NULL;
	return;
}


void HashClear(void)
{
	HASH_ITER(hh,records,p,tmp){
		HASH_DEL(records,p);
		free(p);
	}
	return;
}

