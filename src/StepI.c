#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Istruct.h"
#include <zlib.h>
#include <ctype.h>

#ifndef MAX_READ_LENGTH
	#define MAX_READ_LENGTH 400
#endif
#ifndef MAX_READ_LENGTH_
        #define MAX_READ_LENGTH_ MAX_READ_LENGTH+2
#endif

void Upper_check(char read[MAX_READ_LENGTH_]);

struct locat Position[2*MAX_READ_LENGTH];
extern	int     BestNumber;
extern  int     SeedLength;
extern  int     QuanlityUsed;
extern double AS_threshold;
extern float   RewardVal;
extern float   PenaltyVal;
extern  int     OL_threshold;
extern gzFile	OutFile;
extern struct emp_freq LK;
extern char  Maximal_quality_c;
extern int MinResidualLen;
extern int Step2_flag,Step3_flag;
extern int R_end_Q_f;

extern double prAGC2,  prATC2, prTGC2, prTAG2,
	   pr2AGC, pr2ATC, pr2TGC, pr2TAG,
	   prAGC,  prATC,  prTGC,  prTAG,
	   prAG2,  prAT2,  prAC2,  prGT2,  prGC2,  prTC2,
	   pr2AG,  pr2AT,  pr2AC,  pr2GT,  pr2GC,  pr2TC;

extern	char IdName1[MAX_READ_LENGTH],  IdName2[MAX_READ_LENGTH],
		        Symbol1[MAX_READ_LENGTH],      Symbol2[MAX_READ_LENGTH],
		        Read1[MAX_READ_LENGTH_],       Read2[MAX_READ_LENGTH_],
		        Read[MAX_READ_LENGTH_],        QuanlVal[MAX_READ_LENGTH_],
				QuanlVal1[MAX_READ_LENGTH_],	QuanlVal2[MAX_READ_LENGTH_];



// The function separator is used to print the separator for detecting whether the corresponding module works.
// For instance, separator("JAON.DEN",30), the result will be like this: ********** JAON.DEN **********

void separator(char string[],int Lolen)
{
	int i,slen;

	slen = Lolen/2 - strlen(string)/2;

	for(i=0; i<slen-1; i++)
		printf("*");

	printf(" ");
	printf("%s",string);
	printf(" ");

	for(i=slen+strlen(string); i<Lolen; i++)
		printf("*");

	printf("\n");
}





//convert the read to merge and find seeds.

void ReadQuanlityAdjust(char *read, char *quanl)
{
//	printf("before Reverse:%s",read);
	int i,sLen = strlen(read);
	char tmp;
	for(i=0; i<sLen; i++)
		switch(read[i])
		{
			case 'A':	read[i] = 'T';
						break;
			case 'T':	read[i] = 'A';
						break;
			case 'G':	read[i] = 'C';
						break;
			case 'C':	read[i] = 'G';
						break;
			default:	break;
		}
//printf("after Reverse :%s",read);
	for(i=0; i<(sLen)/2; i++)
	{
		tmp = read[i];
		read[i] = read[sLen-i-1];
		read[sLen-i-1] = tmp;

		tmp = quanl[i];
		quanl[i] = quanl[sLen-i-1];
		quanl[sLen-i-1] = tmp;
	}
//	printf("after Converse:%s",read);
}



// The function seed_opt is created for calculating the aligment position.
// The aligment patameters should be opted as needed.

void SeedSelect(void)
{
		int i,j,lammda;
		int R1len,Rlen;
		int Over_len = 0;
		char R_seed[SeedLength+2];
		memset(R_seed,'\0',sizeof(R_seed));
		char *left_strstr, Read1_cp[MAX_READ_LENGTH_];
		int j_R1=0;
//		Read1_cp = (char *)calloc(MAX_READ_LENGTH_,sizeof(char));
//		left_strstr = (char *)calloc(MAX_READ_LENGTH_,sizeof(char));

		
		R1len = strlen(Read1);
		Rlen = strlen(Read);
// initialize the struct typeof array.
		for(i=0; i< 2*MAX_READ_LENGTH; i++)
		{
			Position[i].posit = i;
			Position[i].number = 0;
			Position[i].over = 0;
		}
		
//identify the aligment position.
		for(i=0; i< Rlen-SeedLength; i++)// read2: get k-mer from each positioin 
		{
			strncpy(R_seed,Read+i,SeedLength);// read2: get k-mer from start position i
			strcpy(Read1_cp, Read1); // Read1_cp store the sequence of read1 for further process 
			//Read1_cp = Read1;
			while(1)// to find the k-mer in read1 including all possible possible, and then break
			{
				if((left_strstr = strstr(Read1_cp,R_seed)) != NULL)// searching whether the k-mer in read1
				{
					int ls_len = strlen(left_strstr);
					j_R1 = R1len - ls_len;
					Position[MAX_READ_LENGTH+(j_R1-i)].number++;
					if(ls_len > SeedLength+1)// cut one nucleotide and the left sequence for further searching the k-mer
					{	
						memset(Read1_cp,'\0',sizeof(MAX_READ_LENGTH_));
						strcpy(Read1_cp,left_strstr+1);
					}else{
						break;
					}
				}else{
					break;
				}
			}
		}

//sort the highest overlap position in descending order. The Bubble ranking algorithm is selected.
/*		for(i=0; i<2*MAX_READ_LENGTH-1; i++)
		{
			for(j=0; j< 2*MAX_READ_LENGTH-1-i; j++)
			{
				if(Position[j].number < Position[j+1].number)
				{
					lammda = Position[j+1].number;
					Position[j+1].number = Position[j].number;
					Position[j].number   = lammda;

					lammda = Position[j+1].posit;
					Position[j+1].posit = Position[j].posit;
					Position[j].posit   = lammda;
				}
			}
		}
*/
		for(i=0; i< 2*MAX_READ_LENGTH; i++)
		{
			if(Position[i].number != 0){
				if(Position[i].posit>=MAX_READ_LENGTH){
					if(Position[i].posit+Rlen < MAX_READ_LENGTH+R1len)
						Position[i].over = Rlen;
					else
						Position[i].over = MAX_READ_LENGTH + R1len - Position[i].posit;
				}else{
					if(Position[i].posit+Rlen < MAX_READ_LENGTH+R1len)
						Position[i].over = Position[i].posit+Rlen-MAX_READ_LENGTH;
					else
						Position[i].over = R1len;
				}
			}
		}

		for(i=0;i<BestNumber;i++)
		{
			for(j=2*MAX_READ_LENGTH-2;j>=i;j--)
			{
				if((Position[j+1].over != 0 && Position[j].over == 0) || (Position[j+1].over != 0 && Position[j].over != 0 && (float)Position[j+1].number/Position[j+1].over > (float)Position[j].number/Position[j].over))
				{
					lammda = Position[j].number;
					Position[j].number = Position[j+1].number;
					Position[j+1].number   = lammda;

					lammda = Position[j].posit;
					Position[j].posit = Position[j+1].posit;
					Position[j+1].posit   = lammda;

					lammda = Position[j].over;
					Position[j].over = Position[j+1].over;
					Position[j+1].over   = lammda;
				}
			}
		}
//		free(Read1_cp);Read1_cp = NULL;
//		free(left_strstr);left_strstr=NULL;
}

// change the lower base into upper base
void Upper_check(char read[MAX_READ_LENGTH_])
{
	char Bases[5] = {'a','c','g','t','n'};
	char *loc;
//	loc = (char *)calloc(MAX_READ_LENGTH_,sizeof(char));
	int i;
	for(i=0; i<=4; i++){
		while(1)
		{
			if((loc = strchr(read,Bases[i])) != NULL)// find the lower base
				read[strlen(read)-strlen(loc)] = toupper(Bases[i]);// change to upper base
			else
				break;
		
		}	
	}
//	free(loc);loc = NULL;

}


// Read the block with four lines.
void BlockRead(FILE *File1, FILE *File2)
{
	int i;
//	printf("DOing Block reading\n");
	
//	memset(Read,'\0',sizeof(Read));
//	memset(QuanlVal,'\0',sizeof(QuanlVal));
	for(i=0; i<4; i++)
		switch (i)
		{
			case 0: fgets(IdName1,MAX_READ_LENGTH,File1);
					fgets(IdName2,MAX_READ_LENGTH,File2);
					break;
			case 1: fgets(Read1, MAX_READ_LENGTH_, File1);
					fgets(Read2, MAX_READ_LENGTH_, File2);
					Read1[strlen(Read1)-1] = '\0';
					Read2[strlen(Read2)-1] = '\0';
					Upper_check(Read1);// check whether there is some lower base, if has, then chang it into upper base
					Upper_check(Read2);
//					printf("%s%s",Read1,Read2);
					break;
			case 2: fgets(Symbol1,MAX_READ_LENGTH,File1);
					fgets(Symbol2,MAX_READ_LENGTH,File2);
//					printf("%s%s",Symbol1,Symbol2);
					break;
			case 3: fgets(QuanlVal1, MAX_READ_LENGTH_, File1);
					fgets(QuanlVal2, MAX_READ_LENGTH_, File2);
					QuanlVal1[strlen(QuanlVal1)-1] = '\0';
					QuanlVal2[strlen(QuanlVal2)-1] = '\0';
//					printf("%s%s",QuanlVal1,QuanlVal2);
					break;
		}
	if(R_end_Q_f > 0)// the low quality bases at the end of reads were trimmed
	{
		int L_raw = strlen(Read1);
                int R_raw = strlen(Read2);
		for(int i=L_raw-1; i>=L_raw/2 ; i--){
			if(QuanlVal1[i] - QuanlityUsed < R_end_Q_f){
				Read1[i] = '\0';
				QuanlVal1[i] = '\0';
			}else
				break;
		}
                for(int i=R_raw-1; i>=R_raw/2 ; i--){
                        if(QuanlVal2[i] - QuanlityUsed < R_end_Q_f){
                                Read2[i] = '\0';
                                QuanlVal2[i] = '\0';
                        }else
                                break;
                }
	}
	
}

// read the *.gz files
void BlockRead_gz(gzFile gzfp1, gzFile gzfp2)
{	
	int i;
	for(i=0; i<4; i++)
	{
		switch (i)
		{
			case 0: gzgets(gzfp1,IdName1,MAX_READ_LENGTH);
				gzgets(gzfp2,IdName2,MAX_READ_LENGTH);
				break;
			case 1: gzgets(gzfp1,Read1,MAX_READ_LENGTH);
				gzgets(gzfp2,Read2,MAX_READ_LENGTH);
				Read1[strlen(Read1)-1] = '\0';
				Read2[strlen(Read2)-1] = '\0';
				Upper_check(Read1);// check whether there is some lower base, if has, then chang it into upper base
				Upper_check(Read2);
				break;
			case 2: gzgets(gzfp1,Symbol1,MAX_READ_LENGTH);
				gzgets(gzfp2,Symbol2,MAX_READ_LENGTH);
				break;
			case 3: gzgets(gzfp1,QuanlVal1,MAX_READ_LENGTH);
			        gzgets(gzfp2,QuanlVal2,MAX_READ_LENGTH);
				QuanlVal1[strlen(QuanlVal1)-1] = '\0';
				QuanlVal2[strlen(QuanlVal2)-1] = '\0';
			        break;
		}
	}
        
	if(R_end_Q_f > 0)// the low quality bases at the end of reads were trimmed
	{
		int L_raw = strlen(Read1);
		int R_raw = strlen(Read2);
                for(int i=L_raw-1; i>=L_raw/2 ; i--){
                        if(QuanlVal1[i] - QuanlityUsed < R_end_Q_f){
                                Read1[i] = '\0';
                                QuanlVal1[i] = '\0';
                        }else
                                break;
                }
                for(int i=R_raw-1; i>=R_raw/2 ; i--){
                        if(QuanlVal2[i] - QuanlityUsed < R_end_Q_f){
                                Read2[i] = '\0';
                                QuanlVal2[i] = '\0';
                        }else
                                break;
                }
        }
}

void BlockRead_step2(FILE *File1, FILE *File2)
{
        int i;
        for(i=0; i<4; i++)
                switch (i)
                {
                        case 0: fgets(IdName1,MAX_READ_LENGTH,File1);
                                        fgets(IdName2,MAX_READ_LENGTH,File2);
                                        break;
                        case 1: fgets(Read1, MAX_READ_LENGTH_, File1);
                                        fgets(Read2, MAX_READ_LENGTH_, File2);
                                        Read1[strlen(Read1)-1] = '\0';
                                        Read2[strlen(Read2)-1] = '\0';
                                        Upper_check(Read1);// check whether there is some lower base, if has, then chang it into upper base
                                        Upper_check(Read2);
//                                      printf("%s%s",Read1,Read2);
                                        break;
                        case 2: fgets(Symbol1,MAX_READ_LENGTH,File1);
                                        fgets(Symbol2,MAX_READ_LENGTH,File2);
//                                      printf("%s%s",Symbol1,Symbol2);
                                        break;
                        case 3: fgets(QuanlVal1, MAX_READ_LENGTH_, File1);
                                        fgets(QuanlVal2, MAX_READ_LENGTH_, File2);
                                        QuanlVal1[strlen(QuanlVal1)-1] = '\0';
                                        QuanlVal2[strlen(QuanlVal2)-1] = '\0';
//                                      printf("%s%s",QuanlVal1,QuanlVal2);
                                        break;
                }
}

float Phred(char QuANL)
{
	float k=0;
	k = pow(10,(float)(QuanlityUsed - QuANL)/10);
	return k;
}


float uPhTIMEuPh(char QuANL1, char QuANL2)
{
	float k=0;
	k = (1-Phred(QuANL1))*(1-Phred(QuANL2));
	return k;
}


float uPhTIMEPh(char QuANL1, char QuANL2)
{
	float k=0;
	k = (1-Phred(QuANL1))*Phred(QuANL2);
	return k;
}


float PhTIMEPh(char QuANL1, char QuANL2)
{
	float k=0;
	k = Phred(QuANL1)*Phred(QuANL2);
	return k;
}


double MatchedScore(char Read,const char QuANL1, char QuANL2)
{
	double k=0;
	switch(Read)
	{
		case 'A':	k = ( uPhTIMEuPh(QuANL1,QuANL2)+ PhTIMEPh(QuANL1,QuANL2)*prTGC2/pr2TGC);
					break;
		case 'G':	k = ( uPhTIMEuPh(QuANL1,QuANL2)+ PhTIMEPh(QuANL1,QuANL2)*prATC2/pr2ATC);
					break;
		case 'T':	k = ( uPhTIMEuPh(QuANL1,QuANL2)+ PhTIMEPh(QuANL1,QuANL2)*prAGC2/pr2AGC);
					break;
		case 'C':	k = ( uPhTIMEuPh(QuANL1,QuANL2)+ PhTIMEPh(QuANL1,QuANL2)*prTAG2/pr2TAG);
					break;
		case 'N':	k = LK.q;
					break;
		default :	break;
	}
//printf("QuANL1: %c/%c  QuANL2:%c/%c MatchedScore:%f\n",Read,QuANL1,Read,QuANL2,k);
	return k;
}



double uMatchedScore(const char Read1, const char Read2, const char QuANL1, const char QuANL2)
{
	double k=0;
	switch(Read1)
	{
	case 'A':switch(Read2)
			 {
			case 'G':k = uPhTIMEPh(QuANL1,QuANL2)*LK.prA/prTGC + uPhTIMEPh(QuANL2,QuANL1)*LK.prG/prATC + PhTIMEPh(QuANL1,QuANL2)*prTC2/(prTGC*prATC);
					break;
			case 'T':k = uPhTIMEPh(QuANL1,QuANL2)*LK.prA/prTGC + uPhTIMEPh(QuANL2,QuANL1)*LK.prT/prAGC + PhTIMEPh(QuANL1,QuANL2)*prGC2/(prTGC*prAGC);
					break;
			case 'C':k = uPhTIMEPh(QuANL1,QuANL2)*LK.prA/prTGC + uPhTIMEPh(QuANL2,QuANL1)*LK.prC/prTAG + PhTIMEPh(QuANL1,QuANL2)*prGT2/(prTGC*prTAG);
					break;
			case 'N':k = LK.prA;
					 break;
			default :break;
			 }
			break;

	case 'G':switch(Read2)
			{
			case 'A':k = uPhTIMEPh(QuANL1,QuANL2)*LK.prG/prATC + uPhTIMEPh(QuANL2,QuANL1)*LK.prA/prTGC + PhTIMEPh(QuANL1,QuANL2)*prTC2/(prATC*prTGC);
					break;
			case 'T':k = uPhTIMEPh(QuANL1,QuANL2)*LK.prG/prATC + uPhTIMEPh(QuANL2,QuANL1)*LK.prT/prAGC + PhTIMEPh(QuANL1,QuANL2)*prAC2/(prATC*prAGC);
					break;
			case 'C':k = uPhTIMEPh(QuANL1,QuANL2)*LK.prG/prATC + uPhTIMEPh(QuANL2,QuANL1)*LK.prC/prTAG + PhTIMEPh(QuANL1,QuANL2)*prAT2/(prATC*prTAG);
					break;
			case 'N':k = LK.prG;
					 break;
			default :break;

			}
			break;

	case 'T':switch(Read2)
			{
			case 'A':k = uPhTIMEPh(QuANL1,QuANL2)*LK.prT/prAGC + uPhTIMEPh(QuANL2,QuANL1)*LK.prA/prTGC + PhTIMEPh(QuANL1,QuANL2)*prGC2/(prAGC*prTGC);
					break;
			case 'G':k = uPhTIMEPh(QuANL1,QuANL2)*LK.prT/prAGC + uPhTIMEPh(QuANL2,QuANL1)*LK.prG/prATC + PhTIMEPh(QuANL1,QuANL2)*prAC2/(prAGC*prATC);
					break;
			case 'C':k = uPhTIMEPh(QuANL1,QuANL2)*LK.prT/prAGC + uPhTIMEPh(QuANL2,QuANL1)*LK.prC/prTAG + PhTIMEPh(QuANL1,QuANL2)*prAG2/(prAGC*prTAG);
					break;
			case 'N':k = LK.prT;
					 break;
			default :break;
			}
			break;

	case 'C':switch(Read2)
			{
			case 'A':k = uPhTIMEPh(QuANL1,QuANL2)*LK.prC/prTAG + uPhTIMEPh(QuANL2,QuANL1)*LK.prA/prTGC + PhTIMEPh(QuANL1,QuANL2)*prGT2/(prTAG*prTGC);
					break;
			case 'G':k = uPhTIMEPh(QuANL1,QuANL2)*LK.prC/prTAG + uPhTIMEPh(QuANL2,QuANL1)*LK.prG/prATC + PhTIMEPh(QuANL1,QuANL2)*prAT2/(prTAG*prATC);
					break;
			case 'T':k = uPhTIMEPh(QuANL1,QuANL2)*LK.prC/prTAG + uPhTIMEPh(QuANL2,QuANL1)*LK.prT/prAGC + PhTIMEPh(QuANL1,QuANL2)*prAG2/(prTAG*prAGC);
					break;
			case 'N':k = LK.prC;
					 break;
			default :break;
			}
			break;
	case 'N':switch(Read2)
		 	{
			case 'A':k = LK.prA;
			case 'G':k = LK.prG;
			case 'C':k = LK.prC;
			case 'T':k = LK.prT;
		 	}
			 break;
	default:break;
	}
//printf("LK.prA: %f LK.prG:%f LK.prT: %f LK.prC: %f  QuANL1: %c/%c  QuANL2: %c/%c uMatchedScore:%f\n",LK.prA,LK.prG,LK.prT,LK.prC,Read1,QuANL1,Read2,QuANL2,k);
	return k;
}


// The function(best_ovlat) is used to select the best overlap.
void BestOverlap(struct DeterPara *para)
{
	int i,j,ostart;
	int Read1Len,ReadLen;

	memset(Read,'\0',sizeof(Read));
	memset(QuanlVal,'\0',sizeof(QuanlVal));
	strcpy(Read,Read2);
	strcpy(QuanlVal,QuanlVal2);
	Read1Len = strlen(Read1);
	ReadLen = strlen(Read);
	ReadQuanlityAdjust(Read,QuanlVal);
	SeedSelect();// get the overlapped position by k-mers
//	for(i=0; i< BestNumber; i++)
//	{
//		para[i].beloc[0] = Position[i].number;
//		para[i].beloc[1] = Position[i].posit;
//	}
	
	if(Position[0].number == 0) // no seed found in read1 and return
	{
		return;
	}
	//determine the overlap's length, matched rate, para[j].AS.
	for(j=0; j< BestNumber; j++)
	{
		if(Position[j].number == 0) //no seed found in read1 and break(because Position.number is sorted before)
			break;
		para[j].beloc[0] = Position[j].number;
		para[j].beloc[1] = Position[j].posit;
		ostart = 0;
		double AS = 0;

		if(para[j].beloc[1]>=MAX_READ_LENGTH)// 1. insert-size >= single read
		{
			if(MAX_READ_LENGTH + Read1Len <= para[j].beloc[1]+ReadLen){
				para[j].OvlapLength = MAX_READ_LENGTH + Read1Len - para[j].beloc[1];
			}else{
				para[j].OvlapLength = ReadLen;
			}
			if(para[j].OvlapLength < OL_threshold){
				continue;
			}

			char R[para[j].OvlapLength+1],Q[para[j].OvlapLength+1];
			memset(R,'\0',sizeof(R));memset(Q,'\0',sizeof(Q));
			int flag = 0;
			for(i=0; i<para[j].OvlapLength; i++)// calculate the total score for overalpped region
			{	
				if(Read1[para[j].beloc[1]-MAX_READ_LENGTH+i] == Read[i])
				{
					para[j].MatchNum++;
					AS += MatchedScore(Read[i],QuanlVal1[para[j].beloc[1]-MAX_READ_LENGTH+i],QuanlVal[i]);
					R[i] = Read[i];
					Q[i] = Maximal_quality_c;
				}
				else
				{
					//to improve the running speed
					if(i+1-para[j].MatchNum > para[j].OvlapLength*(1-AS_threshold)*1.2){
						flag = 1;
						break;
					}
					AS += uMatchedScore(Read1[para[j].beloc[1]-MAX_READ_LENGTH+i],Read[i],QuanlVal1[para[j].beloc[1]-MAX_READ_LENGTH+i],QuanlVal[i]);
					if(Read1[para[j].beloc[1]-MAX_READ_LENGTH+i] == 'N'){
						R[i] = Read[i];
						Q[i] = QuanlVal[i];
					}else if(Read[i] == 'N'){
						R[i] = Read1[para[j].beloc[1]-MAX_READ_LENGTH+i];
						Q[i] = QuanlVal1[para[j].beloc[1]-MAX_READ_LENGTH+i];
					}else if(QuanlVal1[para[j].beloc[1]-MAX_READ_LENGTH+i] > QuanlVal[i]){
						R[i] = Read1[para[j].beloc[1]-MAX_READ_LENGTH+i];
						Q[i] = QuanlVal1[para[j].beloc[1]-MAX_READ_LENGTH+i];
					}else{
						R[i] = Read[i];
						Q[i] = QuanlVal[i];
					}
				}
			}
			if(flag == 0) //good overlap
			{
				para[j].AS = (float)AS/para[j].OvlapLength;//calculate match score
				para[j].OvStart = ostart;
				//1) before overlapped region
				//strncpy(para[j].MergedRead,Read1,Read1Len-para[j].OvlapLength);
				//strncpy(para[j].MergedQual,QuanlVal1,Read1Len-para[j].OvlapLength);
				strncpy(para[j].MergedRead,Read1,para[j].beloc[1]-MAX_READ_LENGTH);
                                strncpy(para[j].MergedQual,QuanlVal1,para[j].beloc[1]-MAX_READ_LENGTH);
				//2) overlapped region
				strcat(para[j].MergedRead,R);
				strcat(para[j].MergedQual,Q);
				//3) after overlapped region
				strcat(para[j].MergedRead,Read+para[j].OvlapLength);para[j].MergedRead[strlen(para[j].MergedRead)] = '\0';
				strcat(para[j].MergedQual,QuanlVal+para[j].OvlapLength);para[j].MergedQual[strlen(para[j].MergedQual)] = '\0';
			}
		}
		else // 2. insert-size < single read
		{	
			if(MAX_READ_LENGTH + Read1Len >= para[j].beloc[1]+ReadLen){
				para[j].OvlapLength = para[j].beloc[1]+ReadLen-MAX_READ_LENGTH;
			}else{
				para[j].OvlapLength = Read1Len;
			}
			if(para[j].OvlapLength < OL_threshold){
				continue;
			}
			
			char R[para[j].OvlapLength+1],Q[para[j].OvlapLength+1];
			memset(R,'\0',sizeof(R));memset(Q,'\0',sizeof(Q));
			int flag = 0;
			for(i=0; i<para[j].OvlapLength; i++)
			{
				if(Read1[i] == Read[MAX_READ_LENGTH-para[j].beloc[1]+i])
				{
					para[j].MatchNum++;
					AS += MatchedScore(Read1[i],QuanlVal1[i],QuanlVal[MAX_READ_LENGTH-para[j].beloc[1]+i]);
					R[i] = Read1[i];
					Q[i] = Maximal_quality_c;
				}
				else{
					// to improve running speed
					if(i+1-para[j].MatchNum > para[j].OvlapLength*(1-AS_threshold)*1.2){
						flag = 1;
						break;
					}
					AS += uMatchedScore(Read1[i],Read[MAX_READ_LENGTH-para[j].beloc[1]+i],QuanlVal1[i],QuanlVal[MAX_READ_LENGTH-para[j].beloc[1]+i]);
					if(Read1[i] == 'N'){
						R[i] = Read[MAX_READ_LENGTH-para[j].beloc[1]+i];
						Q[i] = QuanlVal[MAX_READ_LENGTH-para[j].beloc[1]+i];
					}else if(Read[MAX_READ_LENGTH-para[j].beloc[1]+i] == 'N'){
						R[i] = Read1[i];
						Q[i] = QuanlVal1[i];
					}else if(QuanlVal1[i] >QuanlVal[MAX_READ_LENGTH-para[j].beloc[1]+i]){
						R[i] = Read1[i];
						Q[i] = QuanlVal1[i];
					}else{
						R[i] = Read[ReadLen - para[j].OvlapLength + i];
						Q[i] = QuanlVal[ReadLen - para[j].OvlapLength + i];
					}
				}
			}
			if(flag == 0)// calculate match score and get the overlapped read and quality
			{
				para[j].AS = (float)AS/para[j].OvlapLength;
				para[j].OvStart = 0;
				strcpy(para[j].MergedRead,R);
				strcpy(para[j].MergedQual,Q);
				if(MAX_READ_LENGTH+Read1Len < para[j].beloc[1]+ReadLen){
					int L_a = Read1Len+MAX_READ_LENGTH-para[j].beloc[1];
					strcat(para[j].MergedRead,Read+L_a);
					strcat(para[j].MergedQual,QuanlVal+L_a);
				}
				para[j].MergedRead[strlen(para[j].MergedRead)] = '\0';
				para[j].MergedQual[strlen(para[j].MergedQual)] = '\0';
			}
		}

	}
}




// This function aims to assemble the pair-end reads conditions satisfied.
// if TWO bases from different reads are identical then the quanlity of assembled reads is seted to be maximal quanlity of two plus the distance between two quanlity.
// if TWO bases are different then the base with higher quanlity,
/*
int Assemble(struct DeterPara BeOvlap, FILE *tmpMergeFile)
{
	char R[2*MAX_READ_LENGTH],Q[2*MAX_READ_LENGTH];
	int i,j, StartPo, Read1Len, ReadLen;

	Read1Len = strlen(Read1);
	ReadLen = strlen(Read);
	memset(R,'\0',sizeof(R));
	memset(Q,'\0',sizeof(Q));
	if(BeOvlap.beloc[1] >= MAX_READ_LENGTH)
	{
		j = 0;
		StartPo = Read1Len-BeOvlap.OvlapLength-1;
		for(i=0; i<StartPo+BeOvlap.OvStart; i++,j++)
		{
			R[j] = Read1[i];
			Q[j] = QuanlVal1[i];
		}
		for(i=StartPo+BeOvlap.OvStart; i<Read1Len-1; i++,j++)
		{
			if(Read1[i] == Read[j+BeOvlap.OvStart])
			{
				R[j] = Read1[i];
				Q[j] = Maximal_quality_c;
			}
			else
			{
				if(QuanlVal1[i] > QuanlVal[j+BeOvlap.OvStart])
				{
					R[j] = Read1[i];
					Q[j] = QuanlVal1[i];
				}
				else
				{
					R[j] = Read[j+BeOvlap.OvStart];
					Q[j] = QuanlVal[j+BeOvlap.OvStart];
				}
			}
		}
		for(i=BeOvlap.OvlapLength; i<ReadLen-1; i++,j++)
		{
			R[j] = Read[i];
			Q[j] = QuanlVal[i];
		}
		R[j] = '\0';	Q[j] = '\0';
	}
	else
	{
		StartPo = ReadLen-1-BeOvlap.OvlapLength;
		for(i=0; i<BeOvlap.OvlapLength; i++)
		{
			if(Read1[i] == Read[i+StartPo])
			{
				R[i] = Read1[i];
				Q[i] = Maximal_quality_c; 
			}
			else
			{
				if(QuanlVal[i+StartPo] > QuanlVal1[i])
				{
					R[i] = Read[i+StartPo];
					Q[i] = QuanlVal[i+StartPo];
				}
				else
				{
					R[i] = Read1[i];
					Q[i] = QuanlVal1[i];
				}
			}
		}
	}

	if(strlen(R) < MinResidualLen){
		return -1;
	}
	// put out the results assembled.
	gzputs(OutFile,IdName1);gzputs(OutFile,R);gzprintf(OutFile,"\n");gzputs(OutFile,Symbol1);gzputs(OutFile,Q);gzprintf(OutFile,"\n");
	if(Step2_flag == 1 || Step3_flag == 1){	
		fputs(IdName1,tmpMergeFile);
		fputs(R,tmpMergeFile);	fprintf(tmpMergeFile,"\n");
		fputs(Symbol1,tmpMergeFile);
		fputs(Q,tmpMergeFile);	fprintf(tmpMergeFile,"\n");
	}
	return 1;
}
*/
