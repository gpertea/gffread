## GffRead 

GFF/GTF utility providing format conversions, filtering, FASTA sequence 
extraction and more.

More details and usage examples can be found in the online paper: [DOI: 10.12688/f1000research.23297.1](http://dx.doi.org/10.12688/f1000research.23297.1).

The official webpage and download packages for this utility can be found online here: 
 http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread

Use `gffread -h` to see the command line usage options.

## Installation
Building this program from source requires the [GCLib](../../../gclib) code 
library. 

```
  cd /some/build/dir
  git clone https://github.com/gpertea/gffread
  cd gffread
  make release
```
This should create the **gffread** binary in the current directory.


