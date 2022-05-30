#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "uthash.h"
#include "Istruct.h"
#include <zlib.h>

char QuanlCeil(char La);
char QuanlFloor(char La);
//extern int AuxiMergedNum;
extern gzFile OutFile;
extern char Maximal_quality_c;
extern int TWO_LENGTH_;
extern int OL_threshold;
extern double MR_threshold_step2_other;
bit64_t StringToULL(const char * sz);
double MatchedScore(char Read,const char QuANL1, char QuANL2);
double uMatchedScore(const char Read1, const char Read2, const char QuANL1, const char QuANL2);

double OverlapScore(char *read1, char *quanl1,char *read2, char *quanl2)
{
	double AS = 0;
	int i, overlen = strlen(read1);
	for(i=0;i<overlen;i++)
	{
		if(read1[i] == read2[i])
			AS += MatchedScore(read1[i],quanl1[i],quanl2[i]);
		else{
			AS += uMatchedScore(read1[i],read2[i],quanl1[i],quanl2[i]);
		}
	}
	AS = (float)AS/overlen;
	return AS;
}


// function: extract the sub-sequence from *pre_do, by the start position and length
void substr(char *pre_do, int start, int length, char *dest)
{
	int i,str_len=strlen(pre_do);
	if(start>0)
	{
		if(str_len>=length+start){
			dest = strncpy(dest,pre_do+start-1,length);
		}else{
			dest = strncpy(dest,pre_do+start-1,str_len-start+1);
		}
	}
	else if(start<0)
	{
		if(str_len>=length-start){
			int new_start =str_len-length+start ;
			dest = strncpy(dest,pre_do+new_start+1,length);
		}else{
			dest = strncpy(dest,pre_do,str_len+start);
		}
	}
}

int InterMatctRate(best_loclen *best_c, char *read1, char *read2,FILE *file, int left_flag, int right_flag)
{
	fseek(file,best_c->loc,0);
	char tmp[MAX_READ_LENGTH*2+1];
	fgets(tmp,MAX_READ_LENGTH*2,file);

	if(left_flag == 1)
	{
		int left_1 = START+TWO_LENGTH_-1;
		int left_2 = best_c->len-strlen(read2);
		if(left_2 <= left_1+1)// no need to compare
			return 1;
		
		float max_mismatch = (left_2-left_1)*MR_threshold_step2_other;
		int match_num=0;
		
		for(int i=left_1;i<left_2;i++)
		{
			if(tmp[i] == read1[i])
				match_num += 1;
			if(i-left_1+1-match_num > max_mismatch)
				return -1;
		}
	}
	
	if(right_flag == 1)
        {
                int right_1 = START+TWO_LENGTH_-1;
                int right_2 = best_c->len-strlen(read1);
                if(right_1+1 >= right_2)
			return 1;

                float max_mismatch = (right_2-right_1)*MR_threshold_step2_other;
                int match_num=0;

                for(int i=right_1;i<right_2;i++)
                {
                        if(tmp[best_c->len-i-1] == read2[strlen(read2)-i-1])
                                match_num += 1;
                        if(i-right_1+1-match_num > max_mismatch)
                                return -1;
                }
        }

	return 1;
}

void AuxiConnect(char *read1,char *quanl1,char *read2,char *quanl2, int read_len, char *Name,char *Symbol)
{
	int i;
	int OV_len = 0, read1_len = strlen(read1), read2_len = strlen(read2);
	int rEAd_1_term,rEAd_2_term;

	struct Merged_II con;

	if(read_len<=read1_len && read_len<=read2_len)
	{
		OV_len = read_len;
		rEAd_1_term = 0;
		rEAd_2_term = read2_len - OV_len;
		for(i=0;i<OV_len;i++)
		{
			if(read1[i+rEAd_1_term] == read2[i+rEAd_2_term])
			{
				con.rEAd[i] = read1[i];
				con.qUANl[i] = Maximal_quality_c;
			}
			else
			{
				if(read1[i+rEAd_1_term] == 'N')
				{
					con.rEAd[i] = read2[i+rEAd_2_term];
					con.qUANl[i] = quanl2[i+rEAd_2_term];
				}else if(read2[i+rEAd_2_term] == 'N'){
					con.rEAd[i] = read1[i+rEAd_1_term];
					con.qUANl[i] = quanl1[i+rEAd_1_term];
				}
				else if(quanl1[i+rEAd_1_term]>quanl2[i+rEAd_2_term])
				{
					con.rEAd[i] = read1[i+rEAd_1_term];
					con.qUANl[i] = quanl1[i+rEAd_1_term];
				}
				else
				{
					con.rEAd[i] = read2[i+rEAd_2_term];
					con.qUANl[i] = quanl2[i+rEAd_2_term];
				}
			}
		}
	}
	else if(read_len >= read1_len && read_len >= read2_len)//type 2
	{
		OV_len = read1_len + read2_len - read_len;
		rEAd_1_term = read_len - read2_len;

		for(i=0;i<rEAd_1_term; i++)
		{
			con.rEAd[i] = read1[i];
			con.qUANl[i] = quanl1[i];
		}

		for(i=rEAd_1_term;i<read1_len;i++)
		{
			if(read1[i] == read2[read1_len-i])
			{
				con.rEAd[i] = read1[i];
				con.qUANl[i] = Maximal_quality_c;
			}
			else
			{
				if(read1[i] == 'N'){
					con.rEAd[i] = read2[read1_len-i];
					con.qUANl[i] = quanl2[read1_len-i];
				}else if(read2[read1_len-i] == 'N'){
					con.rEAd[i] = read1[i];
					con.qUANl[i] = quanl1[i];
				}else if(quanl1[i]>quanl2[read1_len-i])
				{
					con.rEAd[i] = read1[i];
					con.qUANl[i] = quanl1[i];
				}
				else
				{
					con.rEAd[i] = read2[read1_len-i];
					con.qUANl[i] = quanl2[read1_len-i];
				}
			}
		}

		for(i=0;i<read2_len-OV_len;i++)
		{
			con.qUANl[read1_len+i] = quanl2[OV_len+i];
			con.rEAd[read1_len+i] = read2[OV_len+i];
		}
	}

///
        else if(read_len >= read1_len && read_len <= read2_len)// type 3
        {
                OV_len = read1_len;
                rEAd_1_term = read2_len-read_len;

                for(i=0;i<read1_len;i++)
                {
                        if(read1[i] == read2[rEAd_1_term+i])
                        {
                                con.rEAd[i] = read1[i];
                                con.qUANl[i] = Maximal_quality_c;
                        }
                        else
                        {
                                if(read1[i] == 'N'){
                                        con.rEAd[i] = read2[rEAd_1_term+i];
                                        con.qUANl[i] = quanl2[rEAd_1_term+i];
                                }else if(read2[rEAd_1_term+i] == 'N'){
                                        con.rEAd[i] = read1[i];
                                        con.qUANl[i] = quanl1[i];
                                }else if(quanl1[i]>quanl2[rEAd_1_term+i])
                                {
                                        con.rEAd[i] = read1[i];
                                        con.qUANl[i] = quanl1[i];
                                }
                                else
                                {
                                        con.rEAd[i] = read2[rEAd_1_term+i];
                                        con.qUANl[i] = quanl2[rEAd_1_term+i];
                                }
                        }
                }

                for(i=0;i<read2_len-OV_len-rEAd_1_term;i++)
                {
                        con.qUANl[read1_len+i] = quanl2[OV_len+rEAd_1_term+i];
                        con.rEAd[read1_len+i] = read2[OV_len+rEAd_1_term+i];
                }
        }
	
//
	else if(read_len <= read1_len && read_len >= read2_len)// type 4
        {
                OV_len = read2_len;
                rEAd_1_term = read_len-read2_len;

		for(i=0; i<rEAd_1_term;i++)
		{
			con.rEAd[i] = read1[i];
			con.qUANl[i] = quanl1[i];
		}
                for(i=rEAd_1_term;i<read_len;i++)
                {
                        if(read1[i] == read2[i-rEAd_1_term])
                        {
                                con.rEAd[i] = read1[i];
                                con.qUANl[i] = Maximal_quality_c;
                        }
                        else
                        {
                                if(read1[i] == 'N'){
                                        con.rEAd[i] = read2[i-rEAd_1_term];
                                        con.qUANl[i] = quanl2[i-rEAd_1_term];
                                }else if(read2[i-rEAd_1_term] == 'N'){
                                        con.rEAd[i] = read1[i];
                                        con.qUANl[i] = quanl1[i];
                                }else if(quanl1[i]>quanl2[i-rEAd_1_term])
                                {
                                        con.rEAd[i] = read1[i];
                                        con.qUANl[i] = quanl1[i];
                                }
                                else
                                {
                                        con.rEAd[i] = read2[i-rEAd_1_term];
                                        con.qUANl[i] = quanl2[i-rEAd_1_term];
                                }
                        }
                }
        }
	
	con.qUANl[read_len] = '\0';
	con.rEAd[read_len]  = '\0';
	
	char *out_name,out_name2[MAX_READ_LENGTH];
       	if((out_name = strstr(Name,"/1")) != NULL){
        	strncpy(out_name2,Name,strlen(Name)-strlen(out_name));
		out_name2[strlen(Name)-strlen(out_name)] = '\0';
		gzputs(OutFile,out_name2);gzputs(OutFile,"\n");gzputs(OutFile,con.rEAd);gzputs(OutFile,"\n");
	}else{
		gzputs(OutFile,Name);gzputs(OutFile,con.rEAd);gzputs(OutFile,"\n");
	}
	gzputs(OutFile,Symbol);gzputs(OutFile,con.qUANl);gzputs(OutFile,"\n");
}


int Condition_Set(char *read1,char *quanl1,char *read2,char *quanl2,FILE *file,char *Name,char *Symbol)
{
	if(strlen(read1) <= 2*(TWO_LENGTH_+START) || strlen(read2) <= 2*(TWO_LENGTH_+START))
		return 0;

	char *subread1, *subread2;
	subread1 = (char *)calloc((TWO_LENGTH_+1)*2,sizeof(char));
	subread2 = (char *)calloc((TWO_LENGTH_+1)*2,sizeof(char));
	substr(read1,START,TWO_LENGTH_,subread1); // extract the sub-sequence(10bp) from position START
	substr(read2,-START,TWO_LENGTH_,subread2); // extract the last sub-sequence(10bp) from -STARt
	twoterm *key;
	best_loclen *best_c;
	//conne *con;
	//con = NULL;
	key = (twoterm *)malloc(sizeof(twoterm));
	key->seg1 = StringToULL(subread1);
	key->seg2 = StringToULL(subread2);
	extern HashNode *records;
	extern double  AS_threshold_overlap2;
	HashNode *p;
	HASH_FIND(hh,records,key,sizeof(twoterm),p);
	if(p)
	{
		float xAS = 0;
		int i, best = -1;

		// consider overlapped region
		for(i=0;i<p->nValue.count;i++)
		{
			// ovelap length between two reads
			int Overlap_L = 0;
                        // extract the overlapped sequence
                        char *Deed1,*Deed2, *Qdde1, *Qdde2;
                        Deed1 = (char *)calloc((p->nValue.len[i]+1)*2,sizeof(char));
                        Deed2 = (char *)calloc((p->nValue.len[i]+1)*2,sizeof(char));
                        Qdde1 = (char *)calloc((p->nValue.len[i]+1)*2,sizeof(char));
                        Qdde2 = (char *)calloc((p->nValue.len[i]+1)*2,sizeof(char));
			
			// four types of overlap
			if(p->nValue.len[i] >= strlen(read1) && p->nValue.len[i] >= strlen(read2))// normal
			{
				if(strlen(read1)+strlen(read2) < p->nValue.len[i]+OL_threshold){
					free(Deed1);Deed1=NULL;
	        	                free(Deed2);Deed2=NULL;
        	        	        free(Qdde1);Qdde1=NULL;
        		                free(Qdde2);Qdde2=NULL;
					continue;
				}else
					Overlap_L = strlen(read1)+strlen(read2)-p->nValue.len[i];

				substr(read1,-1,Overlap_L,Deed1);substr(quanl1,-1,Overlap_L,Qdde1);
				substr(read2,1,Overlap_L,Deed2);substr(quanl2,1,Overlap_L,Qdde2);
			}
			else if(p->nValue.len[i] >= strlen(read1) && p->nValue.len[i] <= strlen(read2)){
				Overlap_L = strlen(read1);
				strcpy(Deed1,read1);strcpy(Qdde1,quanl1);
				strncpy(Deed2,read2+(strlen(read2)-p->nValue.len[i]),Overlap_L);strncpy(Qdde2,quanl2+(strlen(read2)-p->nValue.len[i]),Overlap_L);
			}else if(p->nValue.len[i] <= strlen(read1) && p->nValue.len[i] >= strlen(read2)){
				Overlap_L = strlen(read2);
				strncpy(Deed1,read1+(p->nValue.len[i]-strlen(read2)),Overlap_L);strncpy(Qdde1,quanl1+(p->nValue.len[i]-strlen(read2)),Overlap_L);
				strcpy(Deed2,read2);strcpy(Qdde2,quanl2);
			}else if(p->nValue.len[i] < strlen(read1) && p->nValue.len[i] < strlen(read2)){// reads contaring adapter sequence 
				Overlap_L = p->nValue.len[i];
				substr(read1,1,Overlap_L,Deed1);substr(quanl1,1,Overlap_L,Qdde1);
				substr(read2,-1,Overlap_L,Deed2);substr(quanl2,-1,Overlap_L,Qdde2);
			}
			
			//calcualte the match soore for overlapped region
			float xAS_cp = OverlapScore(Deed1,Qdde1,Deed2,Qdde2);
			if(xAS_cp > xAS)
			{
				xAS = xAS_cp;
				best = i;
			}
			free(Deed1);Deed1=NULL;
			free(Deed2);Deed2=NULL;
			free(Qdde1);Qdde1=NULL;
			free(Qdde2);Qdde2=NULL;
		}
		// consider other region
		if(best != -1 && xAS > AS_threshold_overlap2)
		{
			best_c = (best_loclen*)malloc(sizeof(best_loclen));
			best_c->loc = p->nValue.loc[best];
			best_c->len = p->nValue.len[best];

			int other_flag = -1;

			if(p->nValue.len[best] <= strlen(read1) &&  p->nValue.len[best] <= strlen(read2))// un-normal
			{
				other_flag = 1;
			}
			else if(p->nValue.len[best] >=  strlen(read1) &&  p->nValue.len[best] >=  strlen(read2))
			{
				other_flag  = InterMatctRate(best_c,read1,read2,file,1,1);
			}
			else if(p->nValue.len[best] >=  strlen(read1) &&  p->nValue.len[best] <=  strlen(read2))
			{
				other_flag  = InterMatctRate(best_c,read1,read2,file,0,1);
			}
			else if(p->nValue.len[best] <=  strlen(read1) &&  p->nValue.len[best] >=  strlen(read2))
			{
				other_flag  = InterMatctRate(best_c,read1,read2,file,1,0);
			}
			
			if(other_flag  == 1)
			{
				AuxiConnect(read1,quanl1,read2,quanl2,p->nValue.len[best],Name,Symbol);// output
				free(best_c);best_c=NULL;
                                free(subread1);subread1=NULL;
                                free(subread2);subread2=NULL;
                                free(key);key=NULL;
                                return 1;
			}else{
				free(best_c);best_c=NULL;
                                free(subread1);subread1=NULL;
                                free(subread2);subread2=NULL;
                                free(key);key=NULL;
                                return 0;
			}
		}
		else{
			free(subread1);subread1=NULL;
                        free(subread2);subread2=NULL;
                        free(key);key=NULL;
			return 0;
		}
	}
}

