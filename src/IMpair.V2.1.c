
/********************************************************/
/*The main function for merger*/
/********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "Istruct.h"
#include "uthash.h"
#include <zlib.h>

#ifndef MAX_READ_LENGTH
	#define MAX_READ_LENGTH 400
#endif

#ifndef MAX_READ_LENGTH_
	#define MAX_READ_LENGTH_ MAX_READ_LENGTH+2
#endif

float	RewardVal, PenaltyVal;
float	MR_threshold;
int	OL_threshold;
double	AS_threshold;
double	AS_threshold_overlap2;
double	MR_threshold_step2_other;
double	GERMLINE_MATCH_RATE;
double GERMLINE_MAT_RATE_2;
int	GERMLINE_OVERLAP_LEN;
int	SeedLength,AuxiMergedNum,AuxiGermlNum;
int	QuanlityUsed;
int	AlignLength, SimiL_threshold;
int	BestNumber;
int BestLocation;
int	MinResidualLen;
int	Maximal_quality;
char	Maximal_quality_c;
int	TWO_LENGTH_;
int Out_gap_N_flag;
int Out_gap_R_flag;
int Step2_flag;
int Step3_flag;
int R_end_Q_f;

struct	emp_freq LK;

struct locat Position[2*MAX_READ_LENGTH];
char IdName1[MAX_READ_LENGTH],	IdName2[MAX_READ_LENGTH],
	 Read1[MAX_READ_LENGTH_],	Read2[MAX_READ_LENGTH_],
	 Symbol1[MAX_READ_LENGTH],	Symbol2[MAX_READ_LENGTH],
	 Read[MAX_READ_LENGTH_],	QuanlVal[MAX_READ_LENGTH_],
	 QuanlVal1[MAX_READ_LENGTH_],QuanlVal2[MAX_READ_LENGTH_];

double	prAGC2,  prATC2, prTGC2, prTAG2,
	pr2AGC,  pr2ATC, pr2TGC, pr2TAG,
	prAGC,   prATC,  prTGC,  prTAG,
	prAG2,   prAT2,  prAC2,  prGT2,  prGC2,  prTC2,
	pr2AG,   pr2AT,  pr2AC,  pr2GT,  pr2GC,  pr2TC;

FILE *File1, *File2, *Gfile;
gzFile DisFile1, DisFile2, OutFile;
gzFile gzfp1, gzfp2;

extern char *optarg;
extern int optind;
extern int opterr;
extern int optopt;
int if_gz_file;

// the below subfunctions are followed in file InitialVal.c
//void FileStaticATGC(const char *File1, const char *File2);
void ComputBasicVal(void);
void Usage(void);
void optInit(int argc, char **argv);
void stro4f(char *str);

//produce the hash storaging the coresponding location and value
void PRODUCE_HASH(FILE *file);
void getGermlineRead(void);
int Judge_Connect(char *Name, char *read1,char *quanl1,char *Symbol, char *read2, char *quanl2);

// the below subfunctions are followed in file IOverlap.c
int Condition_Set(char *read1,char *quanl1,char *read2,char *quanl2,FILE *file,char *Name,char *Symbol);
// the below subfunctions are followed in file IAuxiAss.c
void AuxiAss(char *Read1, char *Read2, char *Quanl1, char *Quanl2,int R1Len, int R2Len, int Len, char *Read, char *QuanL);
float SimilariTy(int _Start, const char *MdRead, const char *Target, int Len);

// the below subfunctions are followed in file Iread.c
void separator(char string[],int Lolen);
void ReadQualityAdjust(char *read, char *quanl);
void SeedSelect(void);
void BlockRead(FILE *File1, FILE *File2);
void BlockRead_gz(gzFile gzfp1, gzFile gzfp2);
void BlockRead_step2(FILE *File1, FILE *File2);
float Phred(char QuANL);
double MatchedScore(char Read,const char QuANL1, char QuANL2);
double uMatchedScore(const char Read1, const char Read2, const char QuANL1, const char QuANL2);
void BestOverlap(struct DeterPara *para);
int  Assemble(struct DeterPara BeOvlap, FILE *tmpMergeFile);
void ReadQuanlityAdjust(char *read, char *quanl);

// the below subfunctions are followed in file IAuxilMg.c
void AuxilOfMgd(char *Pattern1, char *Pattern2, int *Next1, int *Next2, FILE *tmFile3, FILE *tmFile4, FILE *ObjFile);
float MRComputE(const char *Sequence, const char *read, int loc, int Len, int sign, int A);
int KMP(const char *Read1,const char *Pattern, int *next, int sign);
void getNext(char*pattern,int *next);
int getPattern(const char *Read,char *pattern, int sign);

// the below subfunctions are followed in file IAuxilGerm.c
void AuxiloMgdWGerml(char *Pattern1, char *Pattern2, int *Next1, int *Next2, FILE *ObjFile);


int main(int argc, char **argv)
{
	// initialize variables	
	OL_threshold = 10;	MR_threshold = 0.75;		AS_threshold = 0.85;	AS_threshold_overlap2 = 0.6;
	RewardVal = 1,		PenaltyVal = -1;
	SeedLength = 3,		AuxiMergedNum =0, AuxiGermlNum=0;
	QuanlityUsed = 64;	
	AlignLength = 10;	SimiL_threshold = 0;	Maximal_quality = 40;
	float mean_AS;
	MR_threshold_step2_other = 0.9;
	BestNumber = 3;	MinResidualLen = 50;	BestLocation = 3;
	GERMLINE_MATCH_RATE = 0.85;	GERMLINE_MAT_RATE_2 = 0.5;	GERMLINE_OVERLAP_LEN = 0;
	LK.prA = 0.25,		LK.prC = 0.25,		 LK.prT = 0.25,		LK.prG = 0.25;
	TWO_LENGTH_ = 10;	Out_gap_N_flag = 1;	Out_gap_R_flag = 0;	Step2_flag = 0;		Step3_flag = 0;
	R_end_Q_f = 0;
	
	optInit(argc,argv);
	
	Maximal_quality_c = Maximal_quality+QuanlityUsed;// get the quality character
	
	FILE *tmFile1, *tmFile2, *tmpMergeFile;
	FILE *tmFile3, *tmFile4;
	int i,ConnectedNum = 0,TotalNum = 0;
	time_t start, middle, final;
	struct DeterPara MaxOverlap, para[BestNumber];
	char PatterNin1[AlignLength],PatterNin2[AlignLength];
	int Next1[AlignLength], Next2[AlignLength];
//	long Offset1,Offset2;
	memset(Next1,0,sizeof(int)*AlignLength);
	memset(Next2,0,sizeof(int)*AlignLength);
//	struct OVLPstart BestOverlap;
	
	tmFile1 = tmpfile();	tmFile2 = tmpfile();
	tmFile3 = tmpfile();	tmFile4 = tmpfile();
	tmpMergeFile = tmpfile();
	if(tmFile1 == NULL && tmFile2 == NULL && tmpMergeFile != NULL)
	{
		perror("tmpfile");
		printf("Failed to creat tmpfile for step2");
		exit(1);
	}

	start = time(NULL);
	printf("Step I: Processing and connecting paired reads\n");
	int flag = 0;	
	while(1)
	{
		// initialize the struct typeof array.
		for(i=0; i<BestNumber; i++)
		{
			para[i].beloc[0]	= 0;
			para[i].beloc[1]	= 0;
			para[i].MatchNum	= 0;
			para[i].OvlapLength = 0;
			para[i].MatchRate	= 0;
			para[i].AS		= 0;
			memset(para[i].MergedRead,'\0',sizeof(para[i].MergedRead));
			memset(para[i].MergedQual,'\0',sizeof(para[i].MergedQual));
		}
		// input file: general fastq file			
		if(if_gz_file == 0){
			BlockRead(File1,File2);
			if(feof(File1) != 0 || feof(File2) != 0)
				break;
		}
		else if(if_gz_file == 1)// *.gz file
		{
			BlockRead_gz(gzfp1,gzfp2);
			if(gzeof(gzfp1) != 0 || gzeof(gzfp2) != 0)
				break;
		}
		TotalNum++;
		flag++;
			
//	Step I:  merge paired reads		
		BestOverlap(para);// 1. identify the position of reads by the kmers; 2. get the overlapped region and merged sequences;
		MaxOverlap = para[0];
		for(i=0; i<BestNumber; i++){
			if(para[i].AS > MaxOverlap.AS)
				MaxOverlap = para[i];
		}
		// AS_threshold: the average score of overlapped region;
		if(MaxOverlap.AS >= AS_threshold && strlen(MaxOverlap.MergedRead)>= MinResidualLen) // 1. final merged sequences for Step I
		{
			char *out_name,out_name2[MAX_READ_LENGTH];
		        if((out_name = strstr(IdName1,"/1")) != NULL){
				strncpy(out_name2,IdName1,strlen(IdName1)-strlen(out_name));
				out_name2[strlen(IdName1)-strlen(out_name)] = '\0';
				gzputs(OutFile,out_name2);gzputs(OutFile,"\n");gzputs(OutFile,MaxOverlap.MergedRead);gzputs(OutFile,"\n");gzputs(OutFile,Symbol1);gzputs(OutFile,MaxOverlap.MergedQual);gzputs(OutFile,"\n");
			}else{
				gzputs(OutFile,IdName1);gzputs(OutFile,MaxOverlap.MergedRead);gzputs(OutFile,"\n");gzputs(OutFile,Symbol1);gzputs(OutFile,MaxOverlap.MergedQual);gzputs(OutFile,"\n");
			}
			ConnectedNum++; 
			if(Step2_flag == 1){
				fputs(IdName1,tmpMergeFile);
				fputs(MaxOverlap.MergedRead,tmpMergeFile);fputs("\n",tmpMergeFile);
				fputs(Symbol1,tmpMergeFile);
				fputs(MaxOverlap.MergedQual,tmpMergeFile);fputs("\n",tmpMergeFile);
			}
		}else if(MaxOverlap.AS >= AS_threshold)// 2. the length of merged sequence is too short, output
		{
//			int len_flag = Assemble(MaxOverlap,tmpMergeFile);// merge the paired reads
			gzputs(DisFile1,IdName1);gzputs(DisFile1,Read1);gzputs(DisFile1,"\n");gzputs(DisFile1,Symbol1);gzputs(DisFile1,QuanlVal1);gzputs(DisFile1,"\n");
			gzputs(DisFile2,IdName2);gzputs(DisFile2,Read2);gzputs(DisFile2,"\n");gzputs(DisFile2,Symbol2);gzputs(DisFile2,QuanlVal2);gzputs(DisFile2,"\n");
		}
		else // 3. fialed merged sequences
		{
			// whether it will run StepII or StepIII ?
			if(Step2_flag == 1 || Step3_flag == 1)
			{
				fputs(IdName1,tmFile1);fputs(Read1,tmFile1);fputs("\n",tmFile1);fputs(Symbol1,tmFile1);fputs(QuanlVal1,tmFile1);fputs("\n",tmFile1);
				fputs(IdName2,tmFile2);fputs(Read2,tmFile2);fputs("\n",tmFile2);fputs(Symbol2,tmFile2);fputs(QuanlVal2,tmFile2);fputs("\n",tmFile2);
			}else{
				gzputs(DisFile1,IdName1);gzputs(DisFile1,Read1);gzputs(DisFile1,"\n");gzputs(DisFile1,Symbol1);gzputs(DisFile1,QuanlVal1);gzputs(DisFile1,"\n");
				gzputs(DisFile2,IdName2);gzputs(DisFile2,Read2);gzputs(DisFile2,"\n");gzputs(DisFile2,Symbol2);gzputs(DisFile2,QuanlVal2);gzputs(DisFile2,"\n");
			}
		}
	}
	middle = time(NULL);		printf("Running time I(s):\t%lf\n",difftime(middle,start));
	fseek(tmpMergeFile,0,0);	fseek(tmFile1,0,0);	fseek(tmFile2,0,0);
	
	
	printf("StepII and StepIII: Doing the auxilary connection by successed merged sequences and Germline V sequences\n");
	if(Step2_flag == 1)
		PRODUCE_HASH(tmpMergeFile);// store the merged sequences in hash
	
	//conne *result;
	//result = (conne *)malloc(sizeof(conne));
	if(Step3_flag == 1)
		getGermlineRead();
	
	while(1)
	{
		BlockRead_step2(tmFile1,tmFile2);// read one pair of unmerged reads
		if(feof(tmFile1) != 0 || feof(tmFile2) != 0)
			break;
		int flag_2 = 0;
		char Read2_inver[MAX_READ_LENGTH_],QuanlVal2_inver[MAX_READ_LENGTH_];
		strcpy(Read2_inver,Read2);
		strcpy(QuanlVal2_inver,QuanlVal2);
		ReadQuanlityAdjust(Read2_inver,QuanlVal2_inver); // inverse the read2 information, including read and quanlity
// Step II: assist unmberged reads by above merged sequences
		if(Step2_flag == 1)
		{
			int con_flag = 0;
			con_flag = Condition_Set(Read1,QuanlVal1,Read2_inver,QuanlVal2_inver,tmpMergeFile,IdName1,Symbol1);
			
			if(con_flag == 1)
			{
				flag_2 = 1;
				AuxiMergedNum++;
				//gzputs(OutFile,IdName1);gzputs(OutFile,con.rEAd);gzputs(OutFile,"\n");
				//gzputs(OutFile,Symbol1);gzputs(OutFile,con.qUANl);gzputs(OutFile,"\n");
			}else if(Step3_flag == 0)// this is the last step and so the failed reads required to output
			{
				gzputs(DisFile1,IdName1);gzputs(DisFile1,Read1);gzputs(DisFile1,"\n");gzputs(DisFile1,Symbol1);gzputs(DisFile1,QuanlVal1);gzputs(DisFile1,"\n");
				gzputs(DisFile2,IdName2);gzputs(DisFile2,Read2);gzputs(DisFile2,"\n");gzputs(DisFile2,Symbol2);gzputs(DisFile2,QuanlVal2);gzputs(DisFile2,"\n");
			}
		}
// Setp III: assist unmberged reads by germline sequences
		if(Step3_flag == 1 && flag_2 ==0)
		{
			if(Judge_Connect(IdName1,Read1,QuanlVal1,Symbol1, Read2_inver, QuanlVal2_inver) < 0)
			{
				gzputs(DisFile1,IdName1);gzputs(DisFile1,Read1);gzputs(DisFile1,"\n");gzputs(DisFile1,Symbol1);gzputs(DisFile1,QuanlVal1);gzputs(DisFile1,"\n");
				gzputs(DisFile2,IdName2);gzputs(DisFile2,Read2);gzputs(DisFile2,"\n");gzputs(DisFile2,Symbol2);gzputs(DisFile2,QuanlVal2);gzputs(DisFile2,"\n");
			}else{
				AuxiGermlNum++;
			}
		}
		
	}
	
	//free(result);	result = NULL;
	fclose(tmFile1);	fclose(tmFile2);	fclose(tmpMergeFile);

	final = time(NULL);
	printf("Running time II and III(s):\t%lf\n",difftime(final,middle));
	
	printf("Running Total Time(s): %lf\n",difftime(final,start));
	printf("Total_paired_reads(#)\tMerged_sequences(#)\tMerged_Rate(%%)\tMerged_by_StepI(#)\tMerged_by_StepII(#)\tMerged_by_StepIII(#)\n%-16d%-16d%-16.2f%-16d%-16d%-16d\n",TotalNum,ConnectedNum+AuxiGermlNum+AuxiMergedNum,(float)(ConnectedNum+AuxiGermlNum+AuxiMergedNum)/TotalNum*100,ConnectedNum,AuxiMergedNum,AuxiGermlNum);


	fclose(tmFile3); fclose(tmFile4); 
	if(Step3_flag == 1){
		fclose(Gfile);
		Gfile = NULL;
	}
	if(if_gz_file==0){
		fclose(File1); fclose(File2);
		File1 = NULL;  File2 = NULL;
	}else if(if_gz_file==1){
		gzclose(gzfp1);gzclose(gzfp2);
		gzfp1 = NULL; gzfp2 = NULL;
	}
	gzclose(OutFile);gzclose(DisFile1);gzclose(DisFile2);
	tmFile3 = NULL;  tmFile4 = NULL;  DisFile1 = NULL;  DisFile2 = NULL;  OutFile = NULL;
	return 0;
}
