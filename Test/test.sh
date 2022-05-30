../IMperm -a test_1.fq.gz -b test_2.fq.gz -2 test.left.I.1.gz -3 test.left.I.3.gz -o test.out.I.gz -Q 33  > test.S1.I.log
../IMperm -a test_1.fq.gz -b test_2.fq.gz -2 test.left.I.II.1.gz -3 test.left.I.II.3.gz -o test.out.I.II.gz -Q 33 -C > test.S1.I.II.log
../IMperm -a test_1.fq.gz -b test_2.fq.gz -2 test.left.I.III.1.gz -3 test.left.I.III.3.gz -o test.out.I.III.gz -Q 33 -D -R -F ../V_Germline_Ref/Human_TRB_Germline_V_reference.fa > test.S1.I.III.log
../IMperm -a test_1.fq.gz -b test_2.fq.gz -2 test.left.I.II.III.1.gz -3 test.left.I.II.III.3.gz -o test.out.I.II.III.gz -Q 33 -C -D -R -F ../V_Germline_Ref/Human_TRB_Germline_V_reference.fa > test.S1.I.II.III.log

