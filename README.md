# IMperm
Program:        IMpair: Merging Paired-end Reads of Immune Repertoire Data
Version         V1.0.0

Usage:   IMpair -a -b -o -2 -3 [options]

        Compulsory Parametors for Input/Output:
        -a      <str>    Query read1 file with FASTQ format(*.fq or *.fq.gz)
        -b      <str>    Query read2 file with FASTQ format(*.fq or *.fq.gz)
        -o      <str>    Output the merged file (*.fq.gz)
        -2      <str>    Output the failed read1 (*.fq.gz)
        -3      <str>    Output the failed read2 (*.fq.gz)

        General Optional Parameters:
        -Q      <int>    The value is used to decode the sequencing quality score for each nucleotide. Generally, 33 and 64 are common used values [64]
        -L      <int>    The mimimum length of merged sequence [50]
        -T      <int>    cutff of base quality: The bases at the end of read will be trimmed if the base quality is less than -T [0]

        Step I Optional Parameters:
        -k      <int>    The length of seed(k-mer). One read is splitted into many k-mers that are used to identify the overlapped positions of paired reads [3]
        -n      <int>    The number of top [-n] overlapped position candidates that are ranked by the number of mapped k-mers to another read [3]
        -f      <float>  The piror probability(frequency) of four nucleotides:A/T/G/C that will be used for score calculation of overlapped region. The default value is: [0.25/0.25/0.25/0.25]
        -m      <float>  The threshold of match score for overlapped region [0.85]
        -l      <int>    The threshold of length(bp) for overlapped region [10]

        Step II Optional Parameters:
        -C       use step II to merge the remaining reads from above step
        -S      <int>    The length of subsequence extracted from the begining of reads.The subsequence is used to find the merged sequences [10]
        -M      <float>  The threshold of match score for overlapped region [0.6]
        -X      <float>  The threshold of match rate for other region(except overlap and subsequences) [0.9]

        Step III Optional Parameters:
        -D       use step III to merge the remaining reads from above steps
        -F      <str>    The file contains Germline V reference sequences(*.fa: each sequence has two lines: ID(first)+sequence(second)).The default file is human Germline data, including TCRs and BCRs
        -G      <float>  The threshold of match rate between read and Germline sequence [0.85]
        -W      <int>    The threshold of length(bp) for overlapped region. if -W=0, it means no overlap is required [0]

        -Y      <float>  The threshold of match rate for overlapped region [0.5]
        -N or -R         if -W=0, there may be a gap between the paired reads. The gap is filled with -N:(base = 'n' && quality=5) or -R:(base(lowercase) from Germline sequence && quality=20). if -W>0, -N and -R are invalid.[-N]

        -h       For help


--------------------- Note ---------------------


Step I: mgering paired reads by its overlapped region
Step II: for the reads failed to merge in step I, we use the successful merged sequences(step I) to aid for connection
Step III: for the reads failed to merged in step I or II, we use Germline sequences of V genes to aid for connection

1. only use the Step I ot merge PE reads
        IMpair -a -b -o -2 -3 [options]
2. use Step I and  Step II to merge PE reads
        IMpair -a -b -o -2 -3 -C [options]; here -C is compulsory
3. use Step I and Step III to merge PE reads
        IMpair -a -b -o -2 -3 -D [options]; here -D is compulsory
4. use Step I, Step II and Step III to merge PE reads
        IMpair -a -b -o -2 -3 -C -D [options]; here -C and -D are compulsory
if we use Step III to merge PE reads and -W=0, the paired reads maybe no overlap, then:
        -D -w=0 -N(optional): to output merged sequence, the gap(between the paired reads) will be filled with 'n'
        -D -w=0 -R: to output merged sequences,the gap(between the paired reads) will be filled with matched Germline V sequences
        
There is an test data(Test/\*.fq.gz) and example for running(Test/test.sh)

Contact:        Wei Zhang: wzhang287-c@my.cityu.edu.hk
