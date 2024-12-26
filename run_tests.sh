#!/usr/bin/env bash

function err_exit {
 echo -e "Error: $1"
 exit 1
}

prog="../gffread" ## actual program to run as gffread
prog_name=$(basename "$prog")
if [[ ! -x $prog_name ]]; then
 make -j4 || err_exit "Build failed."
fi

cd examples || err_exit "'examples' must exist in current directory'"


### test commands to be executed
# gffread -E -o ann_simple.gff annotation.gff 
# gffread -T -o annotation.gtf annotation.gff 
# gffread -w transcripts.fa -g genome.fa annotation.gff
# gffread -W -x transcripts_CDS.fa -g genome.fa annotation.gff
# gffread -y transcripts_prot.fa -g genome.fa annotation.gff
# gffread -w transcripts.fa -y transcripts_prot.fa -g genome.fa annotation.gff
# gffread --table @id,@chr,@start,@end,@strand,@exons,Name,gene,product -o annotation.tbl annotation.gff
### output filename are following -o, -w, -x or -y
### spread parameters into these arrays, which must have the same length


## output parameters
arr_outs=( "-o ann_simple.gff" "-o annotation.gtf" "-w transcripts.fa" \
           "-x transcripts_CDS.fa" "-y transcripts_prot.fa" \
           "-w transcripts.fa -y transcripts_prot.fa" "-o annotation.tbl" )

num_tests=${#arr_outs[@]}

arr_opts=( "-E" "-T" "-g genome.fa" "-W -g genome.fa" "-g genome.fa" "-g genome.fa" \
           "--table @id,@chr,@start,@end,@strand,@exons,Name,gene,product" )
## input parameters
# We need to repeat "annotation.gff" num_tests times
arr_ins=($(printf "%s " $(seq 1 $num_tests | xargs -I{} echo "annotation.gff")))

tfailed=0
for (( i=0; i<${#arr_ins[@]}; i++ )); do
  tnum=$((i+1))
  echo ">>> Running test $tnum"
  cmd="$prog ${arr_opts[$i]} ${arr_outs[$i]} ${arr_ins[$i]}"
  echo "   $cmd"
  eval "$cmd" || err_exit "Command failed: $cmd"

  # Extract output file names
  IFS=' ' read -ra out_params <<< "${arr_outs[$i]}"
  for (( j=0; j<${#out_params[@]}; j+=2 )); do
    if [[ "${out_params[$j]}" == "-o" || "${out_params[$j]}" == "-w" || \
          "${out_params[$j]}" == "-x" || "${out_params[$j]}" == "-y" ]]; then
      out_file="${out_params[$j+1]}"
      # echo "Comparing $out_file with exp_out/$out_file"

      # Filter out version lines from the first 3 lines if present
      #filter_cmd='{ head -n 3 | grep -v "^# .*'"$prog_name"' "; head -n 3 | grep -v "^#"; tail -n +4; }'
      # Use process substitution to apply the filter to both files
      #if diff -q <(eval "$filter_cmd < \"$out_file\"") <(eval "$filter_cmd < \"exp_out/$out_file\"") ; then
      if diff -q -I '^#' $out_file exp_out/$out_file &>/dev/null; then
        echo " OK."
      else
        echo " ERROR: test failed, output '$out_file' different from expected!"
        ((tfailed++))
      fi
    fi
  done
  echo "---------------------------------"
done
if ((tfailed > 0)); then
  echo "Error: $tfailed tests failed!"
else
  echo "All tests passed successfully!"
fi
