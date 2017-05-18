#!/bin/sh
ver=$(fgrep '#define VERSION ' gffread.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
pack=gffread-$ver
linpack=$pack.Linux_x86_64
echo "preparing $linpack.tar.gz"
echo "-------------------"
/bin/rm -rf $linpack
/bin/rm -f $linpack.tar.gz
mkdir $linpack
make clean
make release
cp gffread $linpack/
tar cvfz $linpack.tar.gz $linpack
ls -l $linpack.tar.gz
echo "On the salzX servers the web files can be updated like this:"
echo "cp $linpack.tar.gz $pack.tar.gz  ~/html/software/stringtie/dl/"
echo "perl -i -pe 's/gffread\-\d\.\d+\.\d+\w?\./gffread-$ver./g' ~/html/software/stringtie/gff.shtml"
