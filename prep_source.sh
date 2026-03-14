#!/usr/bin/env bash
ver=$(fgrep '#define VERSION ' gffread.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
pack=gffread-$ver
echo " preparing source $pack.tar.gz"
echo "----------------------------------"
/bin/rm -rf $pack $pack.tar.gz
mkdir -p $pack/gclib

cp Makefile LICENSE README.md gffread.cpp gff_utils.{h,cpp} $pack/

GCL=./gclib
if [ ! -f $GCL/GBase.h ]; then
  if [ -d .git ]; then
    git submodule update --init gclib
  else
    echo "Error: $GCL/GBase.h not found"
    exit 1
  fi
fi

cp -p $GCL/{GVec,GList,GHashMap,khashl}.hh $GCL/{xxhash,wyhash,GBitVec}.h $pack/gclib/
cp -p $GCL/{GArgs,GBase,gdna,GStr,gff,codons,GFaSeqGet,GFastaIndex}.{h,cpp} $pack/gclib/
tar cvfz $pack.tar.gz $pack
ls -l $pack.tar.gz
