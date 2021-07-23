#!/bin/sh
ver=$(fgrep '#define VERSION ' gffread.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
pack=gffread-$ver
echo " preparing souce $pack.tar.gz"
echo "----------------------------------"
/bin/rm -rf $pack
/bin/rm -f $pack.tar.gz
mkdir $pack
mkdir $pack/gclib
libdir=$pack/gclib/

cp LICENSE README.md gffread.cpp gff_utils.{h,cpp} $pack/
sed 's|\.\./gclib|./gclib|' Makefile > $pack/Makefile

cp ../gclib/{GVec,GList,GHashMap,khashl}.hh ../gclib/xxhash.h ../gclib/wyhash.h ../gclib/GBitVec.h $libdir
cp ../gclib/{GArgs,GBase,gdna,GStr,gff,codons,GFaSeqGet,GFastaIndex}.{h,cpp} $libdir
tar cvfz $pack.tar.gz $pack
ls -l $pack.tar.gz

