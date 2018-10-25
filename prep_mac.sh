#!/bin/sh
ver=$(fgrep '#define VERSION ' gffread.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
pack=gffread-$ver
macpack=$pack.OSX_x86_64
echo "preparing $macpack.tar.gz"
echo "-------------------"
/bin/rm -rf $macpack
/bin/rm -f $macpack.tar.gz
mkdir $macpack
make clean
make release
cp gffread $macpack/
tar cvfz $macpack.tar.gz $macpack
ls -l $macpack.tar.gz
#echo "If you're on igmN machines you can also update the web files:"
echo "scp $macpack.tar.gz salz:~/html/software/stringtie/dl/"
#echo "perl -i -pe 's/gffread\-\d\.\d+\.\d+\./gffread-$ver./g' ~/html/software/gffutils/home.shtml"
