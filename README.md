## GffRead 

GFF/GTF utility providing format conversions, filtering, FASTA sequence 
extraction and more.

More details and usage examples can be found in the paper [DOI: 10.12688/f1000research.23297.1](http://dx.doi.org/10.12688/f1000research.23297.1) which can be also used to cite this software.

The official webpage with download packages for this utility can be found online here: 
 http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread

Use `gffread -h` to see the command line usage options.

The genomic sequence given with `-g` may be a plain FASTA file or a
**BGZF-compressed** FASTA (`-g genome.fa.gz`), as long as the matching `.fai`
and `.gzi` index files are present next to it (as produced by `samtools faidx`,
`bgzip -r`, or the bundled `gclib/gclib-test bgzip` helper). BGZF support is
self-contained (built on zlib); it does not require htslib or samtools.

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

## Updating the bundled gclib core

The `gclib/` directory is a git submodule tracking the shared
[`gclib-core`](https://github.com/gpertea/gclib-core) repository (branch
`master`). **Do not edit the core files under `gclib/` directly** — shared-core
changes must be made in `gclib-core` first, then pulled in here as a submodule
pointer bump:

```
# after the change has been committed and pushed in gclib-core:
cd gffread
git submodule update --remote --checkout gclib   # advance gclib/ to latest master
git add gclib
git commit -m "Bump gclib submodule to latest core"
git push
```

To simply build against the currently pinned core commit, no submodule update is
needed — `make` auto-initializes `./gclib` when required.

