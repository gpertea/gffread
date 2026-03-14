## GffRead 

GFF/GTF utility providing format conversions, filtering, FASTA sequence 
extraction and more.

More details and usage examples can be found in the paper [DOI: 10.12688/f1000research.23297.1](http://dx.doi.org/10.12688/f1000research.23297.1) which can be also used to cite this software.

The official webpage with download packages for this utility can be found online here: 
 http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread

Use `gffread -h` to see the command line usage options.

## Installation
Building this program from source requires the `gclib` core sources included as a git submodule.

```
  cd /some/build/dir
  git clone --recurse-submodules https://github.com/gpertea/gffread
  cd gffread
  make release
```
Plain `git clone ...` does not fetch required submodules by default.
If you already cloned without `--recurse-submodules`, run:
`git submodule update --init gclib`

Alternatively, running `make` in a normal git checkout will auto-initialize the
default `./gclib` submodule when needed.

If you downloaded the standalone source package `gffread-*.tar.gz`, just unpack
it and run `make release` in the unpacked directory.

This should create the **gffread** binary in the current directory.

