#!/bin/bash -e
ver=$(fgrep '#define VERSION ' gffread.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
srcpack=gffread-$ver
source prep_source.sh
linpack=$pack.Linux_x86_64
echo "preparing $linpack.tar.gz"
echo "-------------------"
/bin/rm -rf $linpack
/bin/rm -f $linpack.tar.gz
mkdir $linpack
cd $srcpack
make clean
make release
cp LICENSE README.md gffread ../$linpack/
cd ..
tar cvfz $linpack.tar.gz $linpack
ls -l $srcpack.tar.gz $linpack.tar.gz
echo "scp $linpack.tar.gz $pack.tar.gz  salz:~/html/software/stringtie/dl/"
echo ".. then on the server: "
echo "perl -i -pe 's/gffread\-\d\.\d+\.\d+\w?\./gffread-$ver./g' ~/html/software/stringtie/gff*.shtml"
