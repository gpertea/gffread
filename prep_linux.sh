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
echo "If you're on igmN machines you can also update the web files:"
echo "cp $linpack.tar.gz $pack.tar.gz  ~/html/software/stringtie/dl/"
echo "perl -i -pe 's/gffread\-\d\.\d+\.\d+\w?\./gffread-$ver./g' ~/html/software/stringtie/gff.shtml"
