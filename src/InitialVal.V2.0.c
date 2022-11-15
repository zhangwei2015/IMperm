#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>
#include <dirent.h>
#include<zlib.h>

#include <time.h>
#include <math.h>
#include "Istruct.h"

#ifndef MAX_READ_LENGTH
        #define MAX_READ_LENGTH 400
#endif

extern struct  emp_freq LK;
extern double  prAGC2,  prATC2, prTGC2, prTAG2,
	   pr2AGC,  pr2ATC, pr2TGC, pr2TAG,
	   prAGC,   prATC,  prTGC,  prTAG,
	   prAG2,   prAT2,  prAC2,  prGT2,  prGC2,  prTC2,
	   pr2AG,   pr2AT,  pr2AC,  pr2GT,  pr2GC,  pr2TC;

void stro4f(char *str);
void time_dir(char *s);
extern double AS_threshold_overlap2,R_threshold_step2_other,GERMLINE_MATCH_RATE,GERMLINE_MAT_RATE_2,MR_threshold_step2_other;
extern int TWO_LENGTH_, Step2_flag, Step3_flag, GERMLINE_OVERLAP_LEN,Out_gap_N_flag,Out_gap_R_flag;

// The necessary values used in calculating the as_score are computed in this module.

void ComputBasicVal(void)
{
//separator("ComputBasicVal",150);

	LK.q = LK.prA*LK.prA + LK.prC*LK.prC + LK.prT*LK.prT + LK.prG*LK.prG;

	prAGC2 = LK.prA*LK.prA + LK.prG*LK.prG + LK.prC*LK.prC;	prATC2 = LK.prA*LK.prA + LK.prT*LK.prT + LK.prC*LK.prC;
	prTGC2 = LK.prT*LK.prT + LK.prG*LK.prG + LK.prC*LK.prC;	prTAG2 = LK.prT*LK.prT + LK.prA*LK.prA + LK.prG*LK.prG;

	pr2AGC = (LK.prA + LK.prG + LK.prC)*(LK.prA + LK.prG + LK.prC);	 pr2ATC = (LK.prA + LK.prT + LK.prC)*(LK.prA + LK.prT + LK.prC);
	pr2TGC = (LK.prG + LK.prT + LK.prC)*(LK.prG + LK.prT + LK.prC);	pr2TAG = (LK.prG + LK.prT + LK.prA)*(LK.prG + LK.prT + LK.prA);

	prAGC = (LK.prA + LK.prG + LK.prC);	prATC = (LK.prA + LK.prT + LK.prC);
	prTGC = (LK.prG + LK.prT + LK.prC);	prTAG = (LK.prG + LK.prT + LK.prA);

	prAG2 = LK.prA*LK.prA + LK.prG*LK.prG;	prAT2 = LK.prA*LK.prA + LK.prT*LK.prT;	prAC2 = LK.prA*LK.prA + LK.prC*LK.prC;
	prGT2 = LK.prG*LK.prG + LK.prT*LK.prT;	prGC2 = LK.prG*LK.prG + LK.prC*LK.prC;	prTC2 = LK.prT*LK.prT + LK.prC*LK.prC;

	pr2AG = (LK.prA + LK.prG)*(LK.prA + LK.prG);	pr2AT = (LK.prA + LK.prT)*(LK.prA + LK.prT);	pr2AC = (LK.prA + LK.prC)*(LK.prA + LK.prC);
	pr2GT = (LK.prG + LK.prT)*(LK.prG + LK.prT);	pr2GC = (LK.prG + LK.prC)*(LK.prG + LK.prC);	pr2TC = (LK.prT + LK.prC)*(LK.prT + LK.prC);

//	printf("FileStaticATGC:   prAGC2: %.5lf\n",prAGC2);

//separator("END",150);

}


void Usage(void)
{
	fprintf(stdout,"\nProgram:        IMperm: Paired-end reads merger for immune repertoire data\n");
	fprintf(stdout,"Version		V1.0.0\n\n");
	fprintf(stdout, "Usage:\t IMpair -a -b -o -1 -2 [options]\n\n");
	fprintf(stdout, "\tCompulsory Parametors for Input/Output:\n");
	fprintf(stdout,"\t-a\t<str>\t Query read1 file with FASTQ format(*.fq or *.fq.gz)\n");
	fprintf(stdout,"\t-b\t<str>\t Query read2 file with FASTQ format(*.fq or *.fq.gz)\n");
	fprintf(stdout,"\t-o\t<str>\t Output the merged file (*.fq.gz)\n");
	fprintf(stdout,"\t-1\t<str>\t Output the failed read1 (*.fq.gz)\n");
	fprintf(stdout,"\t-2\t<str>\t Output the failed read2 (*.fq.gz)\n\n");

	fprintf(stdout, "\tGeneral Optional Parameters:\n");
	fprintf(stdout,"\t-Q\t<int>\t The value is used to decode the sequencing quality score for each nucleotide. Generally, 33 and 64 are common used values [64]\n");
	fprintf(stdout,"\t-L\t<int>\t The mimimum length of merged sequence [50]\n");//
	fprintf(stdout,"\t-T\t<int>\t cutff of base quality: The bases at the end of read will be trimmed if the base quality is less than -T [2]\n\n");
	
	fprintf(stdout, "\tStep I Optional Parameters:\n");
	fprintf(stdout,"\t-k\t<int>\t The length of seed(k-mer). One read is splitted into many k-mers that are used to identify the overlapped positions of paired reads [3]\n");
	fprintf(stdout,"\t-n\t<int>\t The number of top [-n] overlapped position candidates that are ranked by the number of mapped k-mers to another read [3]\n");
	fprintf(stdout,"\t-f\t<float>\t The piror probability(frequency) of four nucleotides:A/T/G/C that will be used for score calculation of overlapped region. The default value is: [0.25/0.25/0.25/0.25]\n");
	fprintf(stdout,"\t-m\t<float>\t The threshold of match score for overlapped region [0.85]\n");
	fprintf(stdout,"\t-l\t<int>\t The threshold of length(bp) for overlapped region [10]\n\n");

	fprintf(stdout, "\tStep II Optional Parameters:\n");
	fprintf(stdout, "\t-C\t use step II to merge the remaining reads from above step\n");
	fprintf(stdout,"\t-S\t<int>\t The length of subsequence extracted from the begining of reads.The subsequence is used to find the merged sequences [10]\n");//TWO_LENGTH_
	fprintf(stdout,"\t-M\t<float>\t The threshold of match score for overlapped region [0.6]\n");//AS_threshold_overlap2
	fprintf(stdout,"\t-X\t<float>\t The threshold of match rate for other region(except overlap and subsequences) [0.9]\n\n");//MR_threshold_step2_other

        fprintf(stdout, "\tStep III Optional Parameters:\n");
	fprintf(stdout, "\t-D\t use step III to merge the remaining reads from above steps\n");
	fprintf(stdout,"\t-F\t<str>\t The file contains Germline V reference sequences(*.fa: each sequence has two lines: ID(first)+sequence(second)).The default file is human Germline data, including TCRs and BCRs\n");
	fprintf(stdout,"\t-G\t<float>\t The threshold of match rate between read and Germline sequence [0.85]\n");//GERMLINE_MATCH_RATE
	fprintf(stdout,"\t-W\t<int>\t The threshold of length(bp) for overlapped region. if -W=0, it means no overlap is required [0]\n\n");//GERMLINE_OVERLAP_LEN
	fprintf(stdout,"\t-Y\t<float>\t The threshold of match rate for overlapped region [0.5]\n");//GERMLINE_MAT_RATE_2
	fprintf(stdout,"\t-N or -R\t if -W=0, there may be a gap between the paired reads. The gap is filled with -N:(base = 'n' && quality=5) or -R:(base(lowercase) from Germline sequence && quality=20). if -W>0, -N and -R are invalid.[-N]\n\n");
	
	fprintf(stdout,"\t-h\t For help\n\n");

	fprintf(stdout,"--------------------- Note ---------------------\n");
	fprintf(stdout,"Step I: mgering paired reads by its overlapped region\n");
	fprintf(stdout,"Step II: for the reads failed to merge in step I, we use the successful merged sequences(step I) to aid for connection\n");
	fprintf(stdout,"Step III: for the reads failed to merged in step I or II, we use Germline sequences of V genes to aid for connection\n\n");
	fprintf(stdout,"1. only use the Step I ot merge PE reads\n");
	fprintf(stdout,"\tIMpair -a -b -o -1 -2 [options]\n");
	fprintf(stdout,"2. use Step I and  Step II to merge PE reads\n");
	fprintf(stdout,"\tIMpair -a -b -o -1 -2 -C [options]; here -C is compulsory\n");
	fprintf(stdout,"3. use Step I and Step III to merge PE reads\n");
	fprintf(stdout,"\tIMpair -a -b -o -1 -2 -D [options]; here -D is compulsory\n");
	fprintf(stdout,"4. use Step I, Step II and Step III to merge PE reads\n");
	fprintf(stdout,"\tIMpair -a -b -o -1 -2 -C -D [options]; here -C and -D are compulsory\n");
	fprintf(stdout,"if we use Step III to merge PE reads and -W=0, the paired reads maybe no overlap, then:\n");
	fprintf(stdout,"\t-D -w=0 -N(optional): to output merged sequence, the gap(between the paired reads) will be filled with 'n'\n"); 
	fprintf(stdout,"\t-D -w=0 -R: to output merged sequences,the gap(between the paired reads) will be filled with matched Germline V sequences\n\n");
	fprintf(stdout,"Contact:        Wei Zhang: wzhang287-c@my.cityu.edu.hk\n\n");
}



/*Initialize Parameters*/
void optInit(int argc,char **argv)
{
	int result;
	int Readfile = 0;
	DIR *dir;
	extern char *optarg;
	extern int optind;
	extern int opterr;
	extern int optopt;
	extern	FILE *File1, *File2, *Gfile;
	extern gzFile gzfp1, gzfp2;
	extern	gzFile DisFile1, DisFile2, OutFile;
	extern	float	RewardVal, PenaltyVal;
	extern	float MR_threshold;
	extern	int  OL_threshold;
	extern	double	AS_threshold;
	extern	int SeedLength;
	extern	int QuanlityUsed;
	extern	int BestNumber;
	extern	int	BestLocation;
	extern	int	MinResidualLen;
	extern int if_gz_file;
	extern int R_end_Q_f;

	opterr = 1;
	int o_signal = 0, d_signal = 0, e_signal = 0, g_signal = 0, m_signal = 0, F_signal = 0;
	char *o_default = "Merged_file.fastq", *d_default = "Discard1.fastq", *e_default = "Discard2.fastq";
	char *HumanGermline  = "V_Germline_Ref/Human_TCR_BCR_Germline_V_reference.fa";

	if(argc == 1){Usage();exit(1);}

	int IdirSignal = 0;
//	time_dir(Idir);
	char  *oFile, *FileName, *dFile, *eFile, *file_flag;
	file_flag=(char*)calloc(4,sizeof(char));

	while((result=getopt(argc,argv,"a:b:o:1:2:Q:L:k:n:T:f:m:l:c:S:M:F:X:G:W:Y:NCDR"))!=-1)
	{
		switch(result)
		{
			//essential arguments
			case 'a':	strncpy(file_flag,optarg+strlen(optarg)-3,3);
					if(!strcmp(file_flag,".gz")){
						if_gz_file = 1;
						if((gzfp1 = gzopen(optarg,"r")) == NULL ){
							printf("Failed to open %s\n",optarg);	
							exit(1);}	
					}else{
						if_gz_file = 0;
						if((File1 = fopen(optarg,"r")) == NULL ){
							printf("Failed to open %s\n",optarg);
							exit(1);}
					}
					break;
			case 'b':	strncpy(file_flag,optarg+strlen(optarg)-3,3);
					if(!strcmp(file_flag,".gz")){
						if((gzfp2 = gzopen(optarg,"r")) == NULL ){
							printf("Failed to open %s\n",optarg);	
							exit(1);}
					}else{
						if((File2 = fopen(optarg,"r")) == NULL ){
							printf("Failed to open %s\n",optarg);
							exit(1);}
					}
					break;
			case 'F':	if((Gfile = fopen(optarg,"r")) == NULL){
						printf("Failed to open %s\n",optarg);	
						exit(1);
					}else{	
						g_signal = 1;
					}	
					break;
			case 'o':	if((oFile = (char *)malloc(sizeof(char)*(strlen(optarg)+1))) != NULL)
						{
							memset(oFile,'\0',strlen(optarg)+1);
							strcat(oFile,optarg);
							o_signal = 1;
						}
						else
						{
							printf("Failed to create the outfile of %s\n",optarg);
							exit(1);
						}
						break;
			case '1':	if((dFile = (char *)malloc(sizeof(char)*(strlen(optarg)+1))) != NULL)
						{
							memset(dFile,'\0',strlen(optarg)+1);
							strcat(dFile,optarg);
							d_signal = 1;
						}
						else
						{
							printf("Failed to create the outfile of %s\n",optarg);
							exit(1);
						}
						break;
			case '2':	e_signal = 1;
						if((eFile = (char *)malloc(sizeof(char)*(strlen(optarg)+1))) != NULL)
						{
							memset(eFile,'\0',strlen(optarg)+1);
							strcat(eFile,optarg);
						}
						else
						{
							printf("Failed to open %s\n",optarg);
							exit(1);
						}
						break;
			case 'T':	R_end_Q_f = atoi(optarg);	break;
			case 'f':	stro4f(optarg);	break;
			case 'm':	AS_threshold = atof(optarg);	break;
			case 'l':	OL_threshold = atoi(optarg);	break;
			case 'k':	SeedLength = atoi(optarg);	break;
			case 'L':   MinResidualLen = atoi(optarg);	break;
			case 'Q':	QuanlityUsed = atoi(optarg);	break;
			case 'n':	BestNumber = atoi(optarg);	break;
			case 'S':       TWO_LENGTH_ = atoi(optarg);    break;
			case 'M':       AS_threshold_overlap2 = atof(optarg);    break;
			case 'X':       MR_threshold_step2_other = atof(optarg);    break;
			case 'G':       GERMLINE_MATCH_RATE = atof(optarg);    break;
			case 'W':       GERMLINE_OVERLAP_LEN = atoi(optarg);    break;
			case 'Y':       GERMLINE_MAT_RATE_2 = atof(optarg);    break;
			case 'N':       Out_gap_N_flag = 1;    break;
			case 'R':	Out_gap_R_flag = 1;Out_gap_N_flag = 0;	break;
			case 'C':       Step2_flag = 1;    break;
			case 'D':       Step3_flag = 1;    break;

			case 'h':
			case '?': Usage();
					  exit(1);
		}
	}
	free(file_flag);file_flag=NULL;
// if no germline file is inputted and the default germline of Human and Monkey will be used.
		if(Step3_flag == 1 && g_signal != 1)
		{		
				char *dir_HumanGermline;
//				dir_HumanGermline=(char *)calloc((strlen(argv[0])+strlen(HumanGermline)+4),sizeof(char));
				dir_HumanGermline=(char *)malloc(sizeof(char)*(strlen(argv[0])+strlen(HumanGermline)+4));
				memset(dir_HumanGermline,'\0',(strlen(argv[0])+strlen(HumanGermline)+4));
				char *getdir;
//				(char *)calloc((strlen(argv[0])+4),sizeof(char));
				getdir=strrchr(argv[0],'/');
				strncpy(dir_HumanGermline,argv[0],strlen(argv[0])-strlen(getdir)+1);
				strcat(dir_HumanGermline,HumanGermline);
			if(access(dir_HumanGermline,0) == 0 && (Gfile = fopen(dir_HumanGermline,"r")) == NULL)
			{
				g_signal = 1;
				printf("Failed to open %s\n",optarg);
				exit(1);
			}
			free(dir_HumanGermline);dir_HumanGermline = NULL;
		
		}
		if((FileName = (char *)malloc(sizeof(char)*(strlen(oFile)+2))) != NULL)
		{
			memset(FileName,'\0',strlen(oFile)+2);
			strcat(FileName,oFile);
			free(oFile); oFile = NULL;
		}else{
			printf("Failed to create the outfile of %s\n",FileName);
			exit(1);
		}
		if((OutFile = gzopen(FileName,"w")) == NULL){printf("Failed to open %s\n",FileName);	exit(1);}
//		printf("Succeed to open the file for ERR\n");
		free(FileName);	FileName = NULL;

			if((FileName = (char *)malloc(sizeof(char)*(strlen(dFile)+2))) != NULL)
			{
				memset(FileName,'\0',strlen(dFile)+2);
				strcat(FileName,dFile);
				free(dFile);	dFile = NULL;
			}
			else
			{
				printf("Failed to create the outfile of %s\n",FileName);
				exit(1);
			}
			
		if((DisFile1 = gzopen(FileName,"w")) == NULL){printf("Failed to open %s\n",FileName);exit(1);}
		free(FileName); FileName = NULL;

			if((FileName = (char *)malloc(sizeof(char)*(strlen(eFile)+2))) != NULL)
			{
				memset(FileName,'\0',strlen(eFile)+2);
				strcat(FileName,eFile);
				free(eFile);	eFile = NULL;
			}
			else
			{
				printf("Failed to create the outfile of %s\n",FileName);
				exit(1);
			}
		if((DisFile2 = gzopen(FileName,"w")) == NULL){printf("Failed to open %s\n",FileName);exit(1);}
		free(FileName);	FileName = NULL;
//	 if(Readfile == 1 && if_gz_file ==0)
//		FileStaticATGC(File1,File2);
//	 else if(Readfile == 1 && if_gz_file == 1)
//		 FileStaticATGC(gzfp1,gzfp2);
	 ComputBasicVal();
}


void stro4f(char *str)
{
	char a[10];
	int i,j=0;
	int num = 0;
	float Freq[4];

	memset(a,'\0',10);

//	printf("in stro4f before for!\n");
	for(i=0; i<strlen(str); i++)
	{
		if(str[i] != '/')
		{
			a[j] = str[i];
			j++;
		}
		else
		{
//			printf("part:%s\n",a);
//			printf("float :%f\n",atof(a));
			Freq[num] = atof(a);
			j = 0;
			num++;
			memset(a,'\0',10);
		}
	}

//	printf("part:%s\n",a);
//	printf("float :%f\n",atof(a));
	Freq[num] = atof(a);

	if(num == 3)
	{
		LK.prA = Freq[0];
		LK.prT = Freq[1];
		LK.prG = Freq[2];
		LK.prC = Freq[3];
	}
	else
		printf("Four values are needed on command -f\n");
}
