#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include "Istruct.h"
/*
#define FASLE 0
#define TRUE 1
#define NUM_GERM 100

typedef int STATUS;

typedef struct _NODE{
	
	int data;
//	char seque[2*MAX_READ_LENGTH];
	char *seque;
	struct _NODE* next;
}NODE;

typedef struct _HASH_TABLE{
	int count;
	NODE *value[NUM_GERM];
}HASH_TABLE;
*/
//#define NUM_GERM 100000
HASH_TABLE* pHashTbl;

HASH_TABLE* HInitialiZe(void)
{
//	extern FILE *ERR;
	pHashTbl = (HASH_TABLE*)malloc(sizeof(HASH_TABLE));
	if(pHashTbl == NULL)
	{
//		fprintf(ERR,"failed to allocate space for 'HASH_TABLE' in HInitialiZe()\n");
		exit(0);
	}
	memset(pHashTbl,0,sizeof(HASH_TABLE));
//	pHashTbl->count = 0;
	return pHashTbl;
}

// Find the exacting Germline sequence by the ID (num, from 1..138 for TRBV)
NODE* HFind(HASH_TABLE* pHashTbl, int data)
{
	NODE* pNode;

	if(pHashTbl == NULL)
		return NULL;
	if((pNode = pHashTbl->value[data]) == NULL)
		return NULL;
	
	while(pNode)
	{
		if(data == pNode->data)
			return pNode;
		pNode = pNode->next;
	}
	return NULL;
}

STATUS HInsert(HASH_TABLE* pHashTbl, int data,char *str)
{
	NODE* pNode = NULL;
	if(pHashTbl == NULL)
		return FASLE;

	if(pHashTbl->value[data] == NULL)
	{
		if(pHashTbl->count >= NUM_GERM)
		return FASLE;

		pHashTbl->count++;
		pNode = (NODE*)malloc(sizeof(NODE));
		if(pNode == NULL)
			return FASLE;
		memset(pNode,0,sizeof(NODE));
		pNode->data = data;
//		pNode->seque = (char *)malloc(sizeof(char)*(2*MAX_READ_LENGTH));
			memcpy(pNode->seque,str,strlen(str));
			pNode->seque[strlen(str)-1] = '\0';

		pHashTbl->value[data] = pNode;
		return TRUE;
	}
	
	if(NULL != HFind(pHashTbl,data))
		return FASLE;
	pNode = pHashTbl->value[data];
	while(NULL != pNode->next)
		pNode = pNode->next;

	pNode->next = (NODE*)malloc(sizeof(NODE));
	pNode->next->data = data;
	return TRUE;
}

STATUS HDelete(HASH_TABLE* pHashTbl, int data)
{
	NODE* pHead;
	NODE* pNode;
	if(NULL == pHashTbl || NULL == pHashTbl->value[data])
		return FASLE;

	if(NULL == (pNode = HFind(pHashTbl,data)))
		return FASLE;
	if(pNode == pHashTbl->value[data])
	{
		pHashTbl->value[data] = pNode->next;
		free(pNode);
		return TRUE;
	}

	pHead = pHashTbl->value[data];
	while(pNode != pHead->next)
		pHead = pHead->next;
	pHead->next = pNode->next;

}

STATUS HFree(HASH_TABLE* pHashTbl)
{
	if(pHashTbl == NULL)
		return FASLE;
	NODE *p = NULL, *q = NULL;
	int i;
	for(i=0;i<pHashTbl->count;i++)
	{
		p = pHashTbl->value[i];
		while(p != NULL)
		{
			p->data = 0;
			q = p->next;
		/*	if(p->seque)
			{
				free(p->seque);
				p->seque = NULL;
			}
		*/	
			if(p)
			{
				free(p);
				p = NULL;
			}
			p = q;
		}
	}
	if(pHashTbl)
	{
		free(pHashTbl);
		pHashTbl = NULL;
	}
	return TRUE;
}
/*
void main(void)
{
	HASH_TABLE* hashtable = HInitialiZe();

	char *p="Jason";
	HInsert(hashtable,1,p);
	HInsert(hashtable,4,p);
	HInsert(hashtable,11,p);
	HInsert(hashtable,21,p);
	NODE* node1 = HFind(hashtable,11);
	NODE* node2 = HFind(hashtable,21);
	printf("hashtable 1: %d ----%s\n",hashtable->value[1]->data,hashtable->value[1]->seque);
	if(hashtable->value[2] == NULL)
		printf("hashtable 2 is null\n");
	printf("hashtable 1: %d \n",node1->data);
	printf("hashtable 1: %d \n",node2->data);
	HDelete(hashtable,21);
	NODE* node3 = HFind(hashtable,21);
	if(node3==NULL)
		printf("21 is cancel\n");
	else
		printf("hashtable 1: %d",node3->data);
}
*/
