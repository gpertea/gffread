#!/usr/bin/env bash

err_exit() {
 echo -e "Error: $1" >&2
 exit 1
}

prog="./gffread"
gclibtest="./gclib/gclib-test"
shift_bases=5368709120
srcdir="examples"
expdir="examples/exp_out"
wrkdir="examples/wrk"
generated_files=(
 "$wrkdir/genome_large.fa.gz"
 "$wrkdir/genome_large.fa.gz.gzi"
 "$wrkdir/genome_large.fa.gz.fai"
 "$wrkdir/annotation_large.gff"
 "$wrkdir/annotation_large.gtf"
 "$wrkdir/annotation_large.gtf.expected"
 "$wrkdir/annotation_large.tbl"
 "$wrkdir/annotation_large.tbl.expected"
 "$wrkdir/transcripts_large.fa"
 "$wrkdir/transcripts_large2.fa"
 "$wrkdir/transcripts_CDS_large.fa"
 "$wrkdir/transcripts_CDS_large.expected.fa"
 "$wrkdir/transcripts_prot_large.fa"
 "$wrkdir/transcripts_prot_large2.fa"
)

if [[ ! -x $prog ]]; then
 make -j4 gffread || err_exit "Build failed."
fi
if [[ ! -x $gclibtest ]]; then
 make -C gclib -j4 gclib-test || err_exit "Build of gclib-test (bgzip helper) failed."
fi

[[ -d "$srcdir" ]] || err_exit "'$srcdir' must exist in current directory"
[[ -d "$expdir" ]] || err_exit "'$expdir' must exist in current directory"

cleanup_generated() {
 rm -f "${generated_files[@]}"
}

mkdir -p "$wrkdir" || err_exit "cannot create $wrkdir"
cleanup_generated
trap 'status=$?; if [[ $status -eq 0 ]]; then cleanup_generated; fi; exit $status' EXIT

shift_feature_file() {
 local infile="$1"
 local outfile="$2"
 SHIFT_BASES="$shift_bases" perl -F'\t' -lane '
  BEGIN { $shift=$ENV{SHIFT_BASES}+0; }
  if ($_ ne "" && $F[0] !~ /^#/ && @F >= 5) {
    $F[3]+=$shift;
    $F[4]+=$shift;
    $_=join("\t", @F);
  }
  print $_;
 ' "$infile" > "$outfile" || err_exit "could not create $outfile"
}

shift_table_expected() {
 local infile="$1"
 local outfile="$2"
 SHIFT_BASES="$shift_bases" perl -F'\t' -lane '
  BEGIN {
    $shift=$ENV{SHIFT_BASES}+0;
    sub shift_ranges {
      my ($s, $shift)=@_;
      return $s if !defined($s) || $s eq "." || $s eq "";
      my @ranges=split(/,/, $s);
      for my $r (@ranges) {
        my ($start, $end)=split(/-/, $r, 2);
        if (defined($start) && defined($end) && $start =~ /^\d+$/ && $end =~ /^\d+$/) {
          $r=($start+$shift)."-".($end+$shift);
        }
      }
      return join(",", @ranges);
    }
  }
  if (@F >= 4 && $F[2] =~ /^\d+$/ && $F[3] =~ /^\d+$/) {
    $F[2]+=$shift;
    $F[3]+=$shift;
  }
  if (@F >= 6) {
    $F[5]=shift_ranges($F[5], $shift);
  }
  print join("\t", @F);
 ' "$infile" > "$outfile" || err_exit "could not create $outfile"
}

shift_cds_expected() {
 local infile="$1"
 local outfile="$2"
 SHIFT_BASES="$shift_bases" perl -lane '
  BEGIN { $shift=$ENV{SHIFT_BASES}+0; }
  if (/^>/) {
    s/(loc:\S+\))(\d+)-(\d+)/$1.($2+$shift)."-".($3+$shift)/e;
  }
  print;
 ' "$infile" > "$outfile" || err_exit "could not create $outfile"
}

# Build a large-offset genome and store it BGZF-compressed (never materializing
# the multi-GB plain FASTA): the mostly-'N' sequence is streamed straight into
# our bgzip helper, which also writes the companion .gzi (block index) and .fai
# (record index with uncompressed offsets). This exercises 64-bit coordinates
# through BGZF random access while keeping the on-disk fixture only a few MB.
generate_large_fasta() {
 local infile="$1"
 local outgz="$2"
 local seqname=""
 local line_len=""
 local full_lines=0
 local remainder=0
 local nline=""
 seqname=$(sed -n '1{s/^>//;p;q;}' "$infile")
 line_len=$(awk 'NR==2 { print length($0); exit }' "$infile")
 [[ -n $seqname ]] || err_exit "could not read FASTA header from $infile"
 [[ -n $line_len && $line_len -gt 0 ]] || err_exit "could not determine FASTA line length from $infile"
 full_lines=$((shift_bases / line_len))
 remainder=$((shift_bases % line_len))
 nline=$(printf '%*s' "$line_len" '' | tr ' ' 'N')
 {
  printf '>%s\n' "$seqname"
  {
    yes "$nline" | head -n "$full_lines" | tr -d '\n'
    if (( remainder > 0 )); then
      printf '%*s' "$remainder" '' | tr ' ' 'N'
    fi
    tail -n +2 "$infile" | tr -d '\n'
  } | fold -w "$line_len"
 } | "$gclibtest" bgzip - "$outgz" || err_exit "could not create $outgz"
}

run_cmd() {
 local label="$1"
 shift
 echo ">>> Running $label"
 echo "   $prog $*"
 "$prog" "$@" || err_exit "command failed for $label"
}

compare_files() {
 local label="$1"
 local outfile="$2"
 local expfile="$3"
 if diff -q -I '^#' "$outfile" "$expfile" &>/dev/null; then
  echo " OK."
 else
  echo " ERROR: output mismatch for $label"
  diff -u -I '^#' "$expfile" "$outfile" | sed -n '1,120p'
  err_exit "output mismatch for $label"
 fi
}

echo "=== Building large-offset fixtures from committed example data"
generate_large_fasta "$srcdir/genome.fa" "$wrkdir/genome_large.fa.gz"
shift_feature_file "$srcdir/annotation.gff" "$wrkdir/annotation_large.gff"
shift_feature_file "$expdir/annotation.gtf" "$wrkdir/annotation_large.gtf.expected"
shift_table_expected "$expdir/annotation.tbl" "$wrkdir/annotation_large.tbl.expected"
shift_cds_expected "$expdir/transcripts_CDS.fa" "$wrkdir/transcripts_CDS_large.expected.fa"

run_cmd "large test 2" -T -o "$wrkdir/annotation_large.gtf" "$wrkdir/annotation_large.gff"
compare_files "large test 2" "$wrkdir/annotation_large.gtf" "$wrkdir/annotation_large.gtf.expected"
echo "---------------------------------"

run_cmd "large test 3" -g "$wrkdir/genome_large.fa.gz" -w "$wrkdir/transcripts_large.fa" "$wrkdir/annotation_large.gff"
compare_files "large test 3" "$wrkdir/transcripts_large.fa" "$expdir/transcripts.fa"
echo "---------------------------------"

run_cmd "large test 4" -W -g "$wrkdir/genome_large.fa.gz" -x "$wrkdir/transcripts_CDS_large.fa" "$wrkdir/annotation_large.gff"
compare_files "large test 4" "$wrkdir/transcripts_CDS_large.fa" "$wrkdir/transcripts_CDS_large.expected.fa"
echo "---------------------------------"

run_cmd "large test 5" -g "$wrkdir/genome_large.fa.gz" -y "$wrkdir/transcripts_prot_large.fa" "$wrkdir/annotation_large.gff"
compare_files "large test 5" "$wrkdir/transcripts_prot_large.fa" "$expdir/transcripts_prot.fa"
echo "---------------------------------"

run_cmd "large test 6" -g "$wrkdir/genome_large.fa.gz" -w "$wrkdir/transcripts_large2.fa" -y "$wrkdir/transcripts_prot_large2.fa" "$wrkdir/annotation_large.gff"
compare_files "large test 6 -w" "$wrkdir/transcripts_large2.fa" "$expdir/transcripts.fa"
compare_files "large test 6 -y" "$wrkdir/transcripts_prot_large2.fa" "$expdir/transcripts_prot.fa"
echo "---------------------------------"

run_cmd "large test 7" --table @id,@chr,@start,@end,@strand,@exons,Name,gene,product -o "$wrkdir/annotation_large.tbl" "$wrkdir/annotation_large.gff"
compare_files "large test 7" "$wrkdir/annotation_large.tbl" "$wrkdir/annotation_large.tbl.expected"
echo "---------------------------------"

echo "All large-offset tests passed successfully!"
