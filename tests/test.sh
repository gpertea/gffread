#!/bin/env bash
files=(hg38test.fa.fai hg38_t_sel.w.fa hg38_t_sel.x.fa hg38_t_sel.y.fa hg38_t_sel.gtf)
/bin/rm -f ${files[*]}
set -e
#memcheck (valgrind --leak-check=full --show-reachable=yes)
fa=${files[0]}
fa=${fa%.*}
#restore the gff3 name from output gtf
gff=${files[4]}
gff=${gff%.*}
gff=$gff.gff3
../gffread -g $fa -w ${files[1]} -x ${files[2]} -y ${files[3]} -T -o ${files[4]} $gff
for f in  ${files[@]}; do
 echo -n "checkink $f .."
 diff $f out_check/$f
 echo -e "\tOK"
done
