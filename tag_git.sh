#!/bin/bash -e
git checkout master
ver=$(fgrep '#define VERSION ' gffread.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
git tag -a "v$ver" -m "release $ver"
git push --tags
