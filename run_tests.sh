#!/bin/bash
cd examples
if [[ ! -d out.test ]]; then
mkdir out.test
fi
arrparms=( "annotation.gff -g genome.fa -w out.test/transcripts.fa" \
"annotation.gff -T -o out.test/annotation_T.gtf" \
'annotation.gff --table @id,@chr,@start,@end,@strand,@exons,Name,gene,product -o out.test/annotation.tbl' )
arrout=( "transcripts.fa" "annotation_T.gtf" "annotation.tbl" )

for i in ${!arrparms[@]}; do
 fout="out.test/${arrout[$i]}"
 /bin/rm -f $fout
 fcmp="out.exp/${arrout[$i]}"
 if [ ! -f $fcmp ]; then
   echo "Error: test file $fcmp does not exist!"
   exit 1
 fi
 n=$i
 ((n++))
 echo "Test ${n}: ${arrparms[$i]}"
 ../gffread ${arrparms[$i]}
 if [ ! -f $fout ]; then
   echo "Error: file $fout not created! Test failed."
   exit 1
 fi
 if diff -q -I '^#' $fout $fcmp &>/dev/null; then
    echo "  OK."
 else
   echo "Error: test failed, output $fout different than $fcmp!"
   #exit 1
 fi
done
