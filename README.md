# gffread
GFF/GTF parsing utility providing format conversions, region filtering, 
FASTA sequence extraction and more.


Use gffread -h to check the usage options.

Compiling this program from source requires my other code 
library, [GCLib](../../../gclib). It can be done like this:

```
  cd /some/build/dir
  git clone https://github.com/gpertea/gclib
  git clone https://github.com/gpertea/gffread
  cd gffread
  make release
```
This should build the **gffread** binary in the current directory.
