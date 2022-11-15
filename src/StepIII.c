#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Istruct.h"
#include <zlib.h> 

#define MIS_NUM_JL 1


int MAX_READ_LENGTH_germline = 500;
void Germline_seedsearch(char *read, char *seque, struct VJmem *Vvalue);
int Remove_len = 10; // the last 10bp of germline sequence is removed
int re_num =1;

struct locat Position2[500+2*MAX_READ_LENGTH];
extern FILE *Gfile;
extern gzFile OutFile;
extern double GERMLINE_MATCH_RATE;
//extern int AuxiGermlNum;
HASH_TABLE *pHashTbl;
extern void Upper_check(char *read);
extern int MinResidualLen;
extern char Maximal_quality_c;
extern int Out_gap_N_flag;
extern int Out_gap_R_flag;
extern int QuanlityUsed;

HASH_TABLE* HInitialiZe(void);
NODE* HFind(HASH_TABLE* pHashTbl, int data);
STATUS HInsert(HASH_TABLE* pHashTbl, int data,char *str);
STATUS HFree(HASH_TABLE* pHashTbl);

int Connect(char *Name, char *read1,char *quanl1,char *read2, char *quanl2, char * Symbol,int Vstart, int Jstart, char *seque);
//VJinfo_mem *Vfield_read2(char *read,char *seque, int position);
//VJinfo_mem *Vfield_read1(char *read,char *seque);
void substr(char *pre_do, int start, int length, char *dest);

void Germline_fa_Adjust(char *read, char *read_reverse);

// reverse the fasta sequence
void Germline_fa_Adjust(char *read, char *read_reverse)
{
        int i,sLen = strlen(read);
        char tmp;

        for(i=0; i<sLen; i++)
                switch(read[i])
                {
                        case 'A':       read_reverse[i] = 'T';
                                                break;
                        case 'T':       read_reverse[i] = 'A';
                                                break;
                        case 'G':       read_reverse[i] = 'C';
                                                break;
                        case 'C':       read_reverse[i] = 'G';
                                                break;
                        default:        read_reverse[i] = read[i];
						break;
                } 
        for(i=0; i<(sLen)/2; i++)
        {
                tmp = read_reverse[i];
                read_reverse[i] = read_reverse[sLen-i-1];
                read_reverse[sLen-i-1] = tmp;
        }
	read_reverse[sLen] = '\0';
	//free(read_reverse); read_reverse = NULL;
}


// read the germline sequences
// importantly, both forward and minus strands are required to store
void getGermlineRead(void)
{
	pHashTbl = HInitialiZe();
	int re = 1;
	while(feof(Gfile) != 1)
	{
		char *read_r;
		read_r = (char *)calloc(MAX_READ_LENGTH_germline,sizeof(char));
		
		if(re%2 == 0)
		{
			char *read_reverse, *read;
			read_reverse = (char *)calloc(MAX_READ_LENGTH_germline,sizeof(char));
			read = (char *)calloc(MAX_READ_LENGTH_germline,sizeof(char));

			fgets(read_r,MAX_READ_LENGTH_germline,Gfile);
			read_r[strlen(read_r)-1] = '\0';
			strncpy(read, read_r,strlen(read_r)-Remove_len);
			Upper_check(read);// change the lowercase to uppercase if has
			HInsert(pHashTbl,re_num,read);// store
			re_num++;
			// get the reverse of Germline sequence
			Germline_fa_Adjust(read, read_reverse);// get the reverse of Germline sequence
			HInsert(pHashTbl,re_num,read_reverse);
			re_num++;
			free(read_reverse);read_reverse = NULL;
			free(read);read = NULL;
		}
		else
			fgets(read_r,MAX_READ_LENGTH_germline,Gfile);
		
		free(read_r);read_r = NULL;
		re++;
	}
}


//  Merging the paired reads assisting by Germline V allele sequences
int Judge_Connect(char *Name, char *read1,char *quanl1,char * Symbol, char *read2, char *quanl2)
{
//	getGermlineRead();
	if(strlen(read1) < 15 || strlen(read2) <15)
		return -1;

	struct VJmem *Vvalue, *Jvalue;
	Vvalue = (struct VJmem *)malloc(sizeof(struct VJmem));
	Jvalue = (struct VJmem *)malloc(sizeof(struct VJmem));
	memset(Vvalue,0,sizeof(struct VJmem));
	memset(Jvalue,0,sizeof(struct VJmem));
	int i,connect_signal;
	NODE* tmp=NULL;	
	for(i=1;i<=re_num;i++)// re_num: the number of germline sequences. to iterating all germline sequences by ID(num)
	{
		if((tmp = HFind(pHashTbl, i))==NULL)// hash searching
		{
			continue;
		}
		Germline_seedsearch(read1,tmp->seque,Vvalue);// find the position of read1 on the reference
		if(Vvalue->Match_Rate >= GERMLINE_MATCH_RATE)
			Germline_seedsearch(read2,tmp->seque,Jvalue);
		if(Vvalue->Match_Rate >= GERMLINE_MATCH_RATE && Jvalue->Match_Rate >= GERMLINE_MATCH_RATE)
		{
			connect_signal = Connect(Name,read1,quanl1,read2,quanl2,Symbol,Vvalue->length,Jvalue->length,tmp->seque);// step3: further merging
			if(connect_signal == 1){
				free(Vvalue); Vvalue = NULL;
				free(Jvalue); Jvalue = NULL;
				return 1;
			}
		}
	}
//	HFree(pHashTbl);
	free(Vvalue); Vvalue = NULL;
	free(Jvalue); Jvalue = NULL;
	return -1;

}


// after finding the anchors of reads on Germline sequence, then this function is used to further merge and judge whether it can success
int Connect(char *Name, char *read1,char *quanl1,char *read2, char *quanl2, char * Symbol,int Vstart, int Jstart, char *seque)
{
	int distanCE = Jstart-Vstart,R1_OV_start, OV_LEN, i;
	int read1_len = strlen(read1), read2_len = strlen(read2);
	int left_loc = Jstart>Vstart?Jstart:Vstart;
	int right_loc = (Vstart+read1_len)<(Jstart+read2_len)?(Vstart+read1_len):(Jstart+read2_len);
	int overlap_len = right_loc-left_loc;
	extern int GERMLINE_OVERLAP_LEN;
	extern double GERMLINE_MAT_RATE_2;
	
//	char *Read,*Quanl;
//	Read  = (char *)calloc(distanCE+read2_len+1,sizeof(char));
//	Quanl = (char *)calloc(distanCE+read2_len+1,sizeof(char));
	if(distanCE+read2_len+1 <= 0)
		return -1;
	char Read[distanCE+read2_len+1],Quanl[distanCE+read2_len+1];
	if(distanCE >= read1_len && GERMLINE_OVERLAP_LEN == 0)// there is a gap between the paired reads because the PCR amplication is too long
	{
		if(Out_gap_N_flag || Out_gap_R_flag)// output the sequence, gap will be filled with "N"
		{
	                R1_OV_start = read1_len;
        	        OV_LEN = distanCE - read1_len;
                	substr(read1,1,read1_len,Read);
	                substr(quanl1,1,read1_len,Quanl);
			if(Out_gap_N_flag){
        	        	for(i=0;i<OV_LEN;i++)
                		{
	                        	Read[read1_len+i] = 'n';
					if(QuanlityUsed == 64){
						Quanl[read1_len+i] = 'E';
					}else{
	                        		Quanl[read1_len+i] = '&';
					}
        	        	}
			}
			else{
				for(i=0;i<OV_LEN;i++)
				{
					Read[read1_len+i] = seque[Vstart+read1_len+i-MAX_READ_LENGTH_germline]+32;
					if(QuanlityUsed == 64){
						Quanl[read1_len+i] = 'T';
					}else{
						Quanl[read1_len+i] = '5';
					}
				}
			}
                	for(i=0;i<read2_len;i++)
	                {
        	                Read[distanCE+i] = read2[i];
                	        Quanl[distanCE+i] = quanl2[i];
	                }
        	        Read[distanCE+i] = '\0';
	                Quanl[distanCE+i] = '\0';

		}else{
			return -1;
		}
	}
	else if(distanCE >= read1_len && GERMLINE_OVERLAP_LEN != 0)// gap is not allowed
	{
		return -1;
	}
	else if(distanCE < 0) // the amplication is too short, both paired reads contain the adaptive sequence
	{
		if(overlap_len < MinResidualLen)// the length of merged sequence is too short
		{
			return -1;
		}else{
			int R2_OV_start = -distanCE;
			OV_LEN = overlap_len;			
			float match_rate_1;
			int match_num_1 = 0;
			for(i=0;i<OV_LEN;i++)// overlap region
			{
				int k = i+R2_OV_start;
				if(read1[i] == read2[k]){
					match_num_1++;
					Read[i] = read2[k];
					Quanl[i] = Maximal_quality_c;
				}else{
					if(quanl1[i] == 'N'){
						Read[i] = read2[k];
						Quanl[i] = quanl2[k];
					}else if(quanl2[k] == 'N'){
						Read[i] = read1[i];
						Quanl[i] = quanl1[i];
					}else if(quanl1[i] > quanl2[k]){
						Read[i] = read1[i];
						Quanl[i] = quanl1[i];
					}else{
						Read[i] = read2[k];
						Quanl[i] = quanl2[k];
					}
				}
			}
                        match_rate_1 = (float)match_num_1/(float)OV_LEN;
			if(match_rate_1 < GERMLINE_MAT_RATE_2)// the match rate of overlapped region less than expect, then the reapd cannot merge
			{
				return -1;
			}
			if(Jstart+strlen(read2) > Vstart+strlen(read1))// at the end of 3', read2 length > read1 length
			{
				int lastpart = strlen(read2)-strlen(read1)-R2_OV_start;
				for(i=0; i<lastpart; i++){
					Read[OV_LEN+i] = read2[R2_OV_start+strlen(read1)+i];
					Quanl[OV_LEN+i] = quanl2[R2_OV_start+strlen(read1)+i];
				}
				Read[OV_LEN+lastpart] = '\0';
                                Quanl[OV_LEN+lastpart] = '\0';

			}else{
				Read[OV_LEN] = '\0';
				Quanl[OV_LEN] = '\0';
			}
		}
	}
	else// with overlap
	{
		if(overlap_len < GERMLINE_OVERLAP_LEN)// overlap length is too short
		{
			return -1;
		}else if(overlap_len > GERMLINE_OVERLAP_LEN)
		{
			R1_OV_start = distanCE;
			OV_LEN = overlap_len;
			substr(read1,1,R1_OV_start,Read);// 1. before overlap region
			substr(quanl1,1,R1_OV_start,Quanl);
			float match_rate_1;
			int match_num_1 = 0;
				
			for(i=0;i<OV_LEN;i++)// 2. overlap region
			{
				if(read1[i+R1_OV_start] == read2[i])
				{
					match_num_1++;
					Read[i+R1_OV_start] = read2[i];
					Quanl[i+R1_OV_start] = Maximal_quality_c;
				}
				else
				{
					if(read1[i+R1_OV_start] == 'N'){
						Read[i+R1_OV_start] = read2[i];
						Quanl[i+R1_OV_start] = quanl2[i];
					}else if(read2[i] == 'N'){
						Read[i+R1_OV_start] = read1[i+R1_OV_start];
						Quanl[i+R1_OV_start] = quanl1[i+R1_OV_start];
					}else if(quanl1[i+R1_OV_start] > quanl2[i])
					{
						Read[i+R1_OV_start] = read1[i+R1_OV_start];
						Quanl[i+R1_OV_start] = quanl1[i+R1_OV_start];
					}
					else
					{
						Read[i+R1_OV_start] = read2[i];
						Quanl[i+R1_OV_start] = quanl2[i];
					}
				}
			}
			match_rate_1 = (float)match_num_1/(float)OV_LEN;
			if(match_rate_1 < GERMLINE_MAT_RATE_2)// the match rate of overlapped region less than expect, then the reapd cannot merge
			{
				return -1;
			}
			if(Vstart+strlen(read1) > Jstart+strlen(read2)){
				Read[R1_OV_start+OV_LEN] = '\0';
				Quanl[R1_OV_start+OV_LEN] = '\0';
			}else{
				for(i=OV_LEN;i<read2_len;i++)// 3. after overlap region
				{
					Read[R1_OV_start+i] = read2[i];
					Quanl[R1_OV_start+i] = quanl2[i];
				}
				Read[R1_OV_start+i] = '\0';
				Quanl[R1_OV_start+i] = '\0';
			}
		}
	}
//	fputs(Name,OutFile);	fputs(Read,OutFile);	fputs(Symbol,OutFile);	fputs(Quanl,OutFile);
	if(strlen(Read) < MinResidualLen)
		return -1;
	char *out_name, out_name2[MAX_READ_LENGTH]; 
	if((out_name = strstr(Name,"/1")) != NULL){
        	strncpy(out_name2,Name,strlen(Name)-strlen(out_name));
		out_name2[strlen(Name)-strlen(out_name)] = '\0';
		gzputs(OutFile,out_name2);gzputs(OutFile,"\n");    gzputs(OutFile,Read);gzputs(OutFile,"\n");    gzputs(OutFile,Symbol);  gzputs(OutFile,Quanl);gzputs(OutFile,"\n");
	}else{
		gzputs(OutFile,Name); gzputs(OutFile,Read);gzputs(OutFile,"\n");    gzputs(OutFile,Symbol);  gzputs(OutFile,Quanl);gzputs(OutFile,"\n");
	}
	return 1;
}

/*
//	find the germline sequence for read1, return the position
VJinfo_mem *Vfield_read1(char *read,char *seque)
{
	VJinfo_mem *info_o;
	info_o = (VJinfo_mem *)malloc(sizeof(VJinfo_mem));
	memset(info_o,0,sizeof(VJinfo_mem));
	int mat_num, len=0;
	int read_len = strlen(read);
	int seque_len = strlen(seque);
	int GERMLINE_MIN_LEN = 15;
	int i,j;

	// Setp1: the left(5') part of read1 is leader sequence and so only 3' of read1 is V region
	for(i=read_len-GERMLINE_MIN_LEN; i>0 ; i--)
	{
		mat_num = 0;
		len = 0;
		for(j=i ; j<read_len ; j++)
		{
			if(seque[j-i] == read[j])
				mat_num++;
			len++;
			if((len-mat_num) > (read_len-j)*(1-GERMLINE_MATCH_RATE))
				break;
		}
		if(j == read_len){
			info_o->Match_Rate = (float)mat_num/(float)len;
			info_o->length = -i;
			return info_o;
		}
	}
	// Step2: the 5' part of read1 is in V region
	for(i=0; i <= seque_len-GERMLINE_MIN_LEN; i++)
	{
		len=0;
		mat_num = 0;
		int max_j = read_len>seque_len-i?seque_len-i:read_len;
		for(j=0; j < max_j; j++)
		{
			if(seque[i+j] == read[j])
				mat_num++;
			len++;
			if((len-mat_num) > max_j*(1-GERMLINE_MATCH_RATE))
				break;
		}
		if(j == max_j)
		{
			info_o->Match_Rate = (float)mat_num/(float)len;
			info_o->length = i;
			return info_o;
		}
	}
	return info_o;
}
*/

void Germline_seedsearch(char *read, char *seque, struct VJmem *info_o)
{
        //VJinfo_mem *info_o;
        //info_o = (VJinfo_mem *)malloc(sizeof(VJinfo_mem));
        //memset(info_o,0,sizeof(VJinfo_mem));
	int GERMLINE_MIN_LEN = 15;
	int SeedLength_2 = 5;
                int i,j,lammda;
                int R1len,Rlen;
                int Over_len = 0;
                char R_seed[SeedLength_2+2];
                memset(R_seed,'\0',sizeof(R_seed));
                //char *left_strstr, *Read1_cp;
                char *left_strstr, Read1_cp[MAX_READ_LENGTH_];
		int j_R1=0;
                //Read1_cp = (char *)calloc(MAX_READ_LENGTH_germline,sizeof(char));
                //left_strstr = (char *)calloc(MAX_READ_LENGTH_germline,sizeof(char));


                R1len = strlen(seque);
                Rlen = strlen(read);
// initialize the struct typeof array.
                for(i=0; i< 2*MAX_READ_LENGTH_germline; i++)
                {
                        Position2[i].posit = i;
                        Position2[i].number = 0;
                        Position2[i].over = 0;
                }

//identify the aligment position.
                for(i=0; i< Rlen-SeedLength_2; i++)// read2: get k-mer from each positioin
                {
                        strncpy(R_seed,read+i,SeedLength_2);// read2: get k-mer from start position i
                        strcpy(Read1_cp, seque); // Read1_cp store the sequence of read1 for further process
                        while(1)// to find the k-mer in read1 including all possible possible, and then break
                        {
                                if((left_strstr = strstr(Read1_cp,R_seed)) != NULL)// searching whether the k-mer in read1
                                {
                                        int ls_len = strlen(left_strstr);
                                        j_R1 = R1len - ls_len;
                                        Position2[MAX_READ_LENGTH_germline+(j_R1-i)].number++;
                                        if(ls_len > SeedLength_2+1)// cut one nucleotide and the left sequence for further searching the k-mer
                                        {
                                                memset(Read1_cp,'\0',sizeof(MAX_READ_LENGTH_germline));
						strcpy(Read1_cp,left_strstr+1);
                                        }else{
                                                break;
                                        }
                                }else{
                                        break;
                                }
                        }
                }
		int R1_max_over = -1;
		int R1_max_num = 0;
		int R1_max_posit = 0;
		int R1_match_n = 0;
		
                for(i=0; i< 2*MAX_READ_LENGTH_germline; i++)// get the best position
                {
                        if(Position2[i].number != 0){
                                if(Position2[i].posit>=MAX_READ_LENGTH_germline){
					if(Position2[i].posit+Rlen < MAX_READ_LENGTH_germline+R1len)
						Position2[i].over = Rlen;
					else
						Position2[i].over = MAX_READ_LENGTH_germline + R1len - Position2[i].posit;
				}else{
					if(Position2[i].posit+Rlen < MAX_READ_LENGTH_germline+R1len)
                                        	Position2[i].over = Position2[i].posit+Rlen-MAX_READ_LENGTH_germline;
					else
						Position2[i].over = R1len;
				}
				if(Position2[i].over != 0 && (float)Position2[i].number/Position2[i].over > (float)R1_max_num/R1_max_over){
					R1_max_over = Position2[i].over;
					R1_max_posit = Position2[i].posit;
					R1_max_num = Position2[i].number;
				}
			}
                }
		
		if(R1_max_over >= GERMLINE_MIN_LEN)// get the match rate of overlap
		{
			if(R1_max_posit >= MAX_READ_LENGTH_germline)
			{
				for(int i = 0; i<R1_max_over; i++)
				{
					if(seque[R1_max_posit-MAX_READ_LENGTH_germline+i] == read[i])
						R1_match_n++;
					if(i-R1_match_n+1 > R1_max_over*(1-GERMLINE_MATCH_RATE))
						break;
				}
			}else{
				for(int i= 0; i<R1_max_over; i++)
				{
					if(seque[i] == read[MAX_READ_LENGTH_germline-R1_max_posit+i])
                                                R1_match_n++;
                                        if(i-R1_match_n+1 > R1_max_over*(1-GERMLINE_MATCH_RATE))
                                                break;

				}
			}
			info_o->Match_Rate = (float)R1_match_n/R1_max_over;
			info_o->length = R1_max_posit;
			//free(Read1_cp);Read1_cp = NULL;
			//free(left_strstr);left_strstr=NULL;
			//return info_o;
		}
//		free(Read1_cp);Read1_cp = NULL;
//                free(left_strstr);left_strstr=NULL;

		//return info_o;
}

