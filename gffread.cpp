#include "GArgs.h"
#include "gff_utils.h"
#include <ctype.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#define VERSION "0.12.7"

#define USAGE "gffread v" VERSION ". Usage:\n\
gffread [-g <genomic_seqs_fasta> | <dir>] [-s <seq_info.fsize>] \n\
 [-o <outfile>] [-t <trackname>] [-r [<strand>]<chr>:<start>-<end> [-R]]\n\
 [--jmatch <chr>:<start>-<end>] [--no-pseudo] \n\
 [-CTVNJMKQAFPGUBHZWTOLE] [-w <exons.fa>] [-x <cds.fa>] [-y <tr_cds.fa>]\n\
 [-j ][--ids <IDs.lst> | --nids <IDs.lst>] [--attrs <attr-list>] [-i <maxintron>]\n\
 [--stream] [--bed | --gtf | --tlf] [--table <attrlist>] [--sort-by <ref.lst>]\n\
 [<input_gff>] \n\n\
 Filter, convert or cluster GFF/GTF/BED records, extract the sequence of\n\
 transcripts (exon or CDS) and more.\n\
 By default (i.e. without -O) only transcripts are processed, discarding any\n\
 other non-transcript features. Default output is a simplified GFF3 with only\n\
 the basic attributes.\n\
 \n\
Options:\n\
 --ids discard records/transcripts if their IDs are not listed in <IDs.lst>\n\
 --nids discard records/transcripts if their IDs are listed in <IDs.lst>\n\
 -i   discard transcripts having an intron larger than <maxintron>\n\
 -l   discard transcripts shorter than <minlen> bases\n\
 -r   only show transcripts overlapping coordinate range <start>..<end>\n\
      (on chromosome/contig <chr>, strand <strand> if provided)\n\
 -R   for -r option, discard all transcripts that are not fully \n\
      contained within the given range\n\
 --jmatch only output transcripts matching the given junction\n\
 -U   discard single-exon transcripts\n\
 -C   coding only: discard mRNAs that have no CDS features\n\
 --nc non-coding only: discard mRNAs that have CDS features\n\
 --ignore-locus : discard locus features and attributes found in the input\n\
 -A   use the description field from <seq_info.fsize> and add it\n\
      as the value for a 'descr' attribute to the GFF record\n\
 -s   <seq_info.fsize> is a tab-delimited file providing this info\n\
      for each of the mapped sequences:\n\
      <seq-name> <seq-length> <seq-description>\n\
      (useful for -A option with mRNA/EST/protein mappings)\n\
Sorting: (by default, chromosomes are kept in the order they were found)\n\
 --sort-alpha : chromosomes (reference sequences) are sorted alphabetically\n\
 --sort-by : sort the reference sequences by the order in which their\n\
      names are given in the <refseq.lst> file\n\
Misc options: \n\
 -F   keep all GFF attributes (for non-exon features)\n\
 --keep-exon-attrs : for -F option, do not attempt to reduce redundant\n\
      exon/CDS attributes\n\
 -G   do not keep exon attributes, move them to the transcript feature\n\
      (for GFF3 output)\n\
 --attrs <attr-list> only output the GTF/GFF attributes listed in <attr-list>\n\
    which is a comma delimited list of attribute names to\n\
 --keep-genes : in transcript-only mode (default), also preserve gene records\n\
 --keep-comments: for GFF3 input/output, try to preserve comments\n\
 -O   process other non-transcript GFF records (by default non-transcript\n\
      records are ignored)\n\
 -V   discard any mRNAs with CDS having in-frame stop codons (requires -g)\n\
 -H   for -V option, check and adjust the starting CDS phase\n\
      if the original phase leads to a translation with an \n\
      in-frame stop codon\n\
 -B   for -V option, single-exon transcripts are also checked on the\n\
      opposite strand (requires -g)\n\
 -P   add transcript level GFF attributes about the coding status of each\n\
      transcript, including partialness or in-frame stop codons (requires -g)\n\
 --add-hasCDS : add a \"hasCDS\" attribute with value \"true\" for transcripts\n\
      that have CDS features\n\
 --adj-stop stop codon adjustment: enables -P and performs automatic\n\
      adjustment of the CDS stop coordinate if premature or downstream\n\
 -N   discard multi-exon mRNAs that have any intron with a non-canonical\n\
      splice site consensus (i.e. not GT-AG, GC-AG or AT-AC)\n\
 -J   discard any mRNAs that either lack initial START codon\n\
      or the terminal STOP codon, or have an in-frame stop codon\n\
      (i.e. only print mRNAs with a complete CDS)\n\
 --no-pseudo: filter out records matching the 'pseudo' keyword\n\
 --in-bed: input should be parsed as BED format (automatic if the input\n\
           filename ends with .bed*)\n\
 --in-tlf: input GFF-like one-line-per-transcript format without exon/CDS\n\
           features (see --tlf option below); automatic if the input\n\
           filename ends with .tlf)\n\
 --stream: fast processing of input GFF/BED transcripts as they are received\n\
           ((no sorting, exons must be grouped by transcript in the input data)\n\
Clustering:\n\
 -M/--merge : cluster the input transcripts into loci, discarding\n\
      \"redundant\" transcripts (those with the same exact introns\n\
      and fully contained or equal boundaries)\n\
 -d <dupinfo> : for -M option, write duplication info to file <dupinfo>\n\
 --cluster-only: same as -M/--merge but without discarding any of the\n\
      \"duplicate\" transcripts, only create \"locus\" features\n\
 -K   for -M option: also discard as redundant the shorter, fully contained\n\
       transcripts (intron chains matching a part of the container)\n\
 -Q   for -M option, no longer require boundary containment when assessing\n\
      redundancy (can be combined with -K); only introns have to match for\n\
      multi-exon transcripts, and >=80% overlap for single-exon transcripts\n\
 -Y   for -M option, enforce -Q but also discard overlapping single-exon \n\
      transcripts, even on the opposite strand (can be combined with -K)\n\
Output options:\n\
 --force-exons: make sure that the lowest level GFF features are considered\n\
       \"exon\" features\n\
 --gene2exon: for single-line genes not parenting any transcripts, add an\n\
       exon feature spanning the entire gene (treat it as a transcript)\n\
 --t-adopt:  try to find a parent gene overlapping/containing a transcript\n\
       that does not have any explicit gene Parent\n\
 -D    decode url encoded characters within attributes\n\
 -Z    merge very close exons into a single exon (when intron size<4)\n\
 -g   full path to a multi-fasta file with the genomic sequences\n\
      for all input mappings, OR a directory with single-fasta files\n\
      (one per genomic sequence, with file names matching sequence names)\n\
 -j    output the junctions and the corresponding transcripts\n\
 -w    write a fasta file with spliced exons for each transcript\n\
 --w-add <N> for the -w option, extract additional <N> bases\n\
       both upstream and downstream of the transcript boundaries\n\
 --w-nocds for -w, disable the output of CDS info in the FASTA file\n\
 -x    write a fasta file with spliced CDS for each GFF transcript\n\
 -y    write a protein fasta file with the translation of CDS for each record\n\
 -W    for -w, -x and -y options, write in the FASTA defline all the exon\n\
       coordinates projected onto the spliced sequence;\n\
 -S    for -y option, use '*' instead of '.' as stop codon translation\n\
 -L    Ensembl GTF to GFF3 conversion, adds version to IDs\n\
 -m    <chr_replace> is a name mapping table for converting reference \n\
       sequence names, having this 2-column format:\n\
       <original_ref_ID> <new_ref_ID>\n\
 -t    use <trackname> in the 2nd column of each GFF/GTF output line\n\
 -o    write the output records into <outfile> instead of stdout\n\
 -T    main output will be GTF instead of GFF3\n\
 --bed output records in BED format instead of default GFF3\n\
 --tlf output \"transcript line format\" which is like GFF\n\
       but with exons and CDS related features stored as GFF \n\
       attributes in the transcript feature line, like this:\n\
         exoncount=N;exons=<exons>;CDSphase=<N>;CDS=<CDScoords> \n\
       <exons> is a comma-delimited list of exon_start-exon_end coordinates;\n\
       <CDScoords> is CDS_start:CDS_end coordinates or a list like <exons>\n\
 --table output a simple tab delimited format instead of GFF, with columns\n\
       having the values of GFF attributes given in <attrlist>; special\n\
       pseudo-attributes (prefixed by @) are recognized:\n\
       @id, @geneid, @chr, @start, @end, @strand, @numexons, @exons, \n\
       @cds, @covlen, @cdslen\n\
       If any of -w/-y/-x FASTA output files are enabled, the same fields\n\
       (excluding @id) are appended to the definition line of corresponding\n\
       FASTA records\n\
 -v,-E expose (warn about) duplicate transcript IDs and other potential\n\
       problems with the given GFF/GTF records\n\
"

GStr sortBy; //file name with chromosomes listed in the desired order

bool BEDinput=false;
bool TLFinput=false;

//bool protmap=false;
//int maxintron=999000000;
//bool mergeCloseExons=false;
//range filter:

GffLoader gffloader;

GList<GenomicSeqData> g_data(true,true,true); //list of GFF records by genomic seq

void loadIDlist(FILE* f, GStrSet<> & idhash) {
  GLineReader fr(f);
  while (!fr.isEof()) {
      char* line=fr.getLine();
      if (line==NULL) break;
      if (line[0]=='#') continue; //skip comments
      GDynArray<char*> ids;
      strsplit(line, ids);
      for (uint i=0;i<ids.Count();i++) {
    	  if (strlen(ids[i])>0)
    		  idhash.Add(ids[i]);
      }
  }
}

void loadSeqInfo(FILE* f, GHash<SeqInfo*> &si) {
  GLineReader fr(f);
  while (!fr.isEof()) {
      char* line=fr.getLine();
      if (line==NULL) break;
      char* id=line;
      char* lenstr=NULL;
      char* text=NULL;
      char* p=line;
      while (*p!=0 && !isspace(*p)) p++;
      if (*p==0) continue;
      *p=0;p++;
      while (*p==' ' || *p=='\t') p++;
      if (*p==0) continue;
      lenstr=p;
      while (*p!=0 && !isspace(*p)) p++;
      if (*p!=0) { *p=0;p++; }
      while (*p==' ' || *p=='\t') p++;
      if (*p!=0) text=p; //else text remains NULL
      int len=0;
      if (!parseInt(lenstr,len)) {
         GMessage("Warning: could not parse sequence length: %s %s\n",
                  id, lenstr);
         continue;
         }
      // --- here we have finished parsing the line
      si.Add(id, new SeqInfo(len,text));
      } //while lines
}

void getAttrList(GStr& s) {
	 if (s.is_empty()) return;
	 s.startTokenize(",;:", tkCharSet);
	 GStr w;
	 while (s.nextToken(w)) {
		  if (w.length()>0)
    	    attrList.Add(w.chars());
	 }
}

void setTableFormat(GStr& s) {
	 if (s.is_empty()) return;
	 GHash<ETableFieldType> specialFields;
	 specialFields.Add("chr", ctfGFF_chr);
	 specialFields.Add("id", ctfGFF_ID);
	 specialFields.Add("geneid", ctfGFF_geneID);
	 specialFields.Add("genename", ctfGFF_geneName);
	 specialFields.Add("parent", ctfGFF_Parent);
	 specialFields.Add("feature", ctfGFF_feature);
	 specialFields.Add("start", ctfGFF_start);
	 specialFields.Add("end", ctfGFF_end);
	 specialFields.Add("strand", ctfGFF_strand);
	 specialFields.Add("numexons", ctfGFF_numexons);
	 specialFields.Add("exons", ctfGFF_exons);
	 specialFields.Add("cds", ctfGFF_cds);
	 specialFields.Add("covlen", ctfGFF_covlen);
	 specialFields.Add("cdslen", ctfGFF_cdslen);

	 s.startTokenize(" ,;.:", tkCharSet);
	 GStr w;
	 while (s.nextToken(w)) {
      if (w[0]=='@') {
    	  w=w.substr(1);
    	  w.lower();
    	  ETableFieldType* v=specialFields.Find(w.chars());
    	  if (v!=NULL) {
    		  CTableField tcol(*v);
    		  tableCols.Add(tcol);
    	  }
    	  else GMessage("Warning: table field '@%s' not recognized!\n",w.chars());
    	  continue;
      }
      if (w=="ID" || w=="transcript_id") {
    	  CTableField tcol(ctfGFF_ID);
    	  tableCols.Add(tcol);
    	  continue;
      }
      if (w=="geneID" || w=="gene_id") {
    	  CTableField tcol(ctfGFF_geneID);
    	  tableCols.Add(tcol);
    	  continue;
      }
      if (w=="Parent") {
    	  CTableField tcol(ctfGFF_Parent);
    	  tableCols.Add(tcol);
    	  continue;
      }
      CTableField col(w);
      tableCols.Add(col);
	 }
}

void loadRefTable(FILE* f, GHash<RefTran*>& rt) {
  GLineReader fr(f);
  char* line=NULL;
  while ((line=fr.getLine())) {
      char* orig_id=line;
      char* p=line;
      while (*p!=0 && !isspace(*p)) p++;
      if (*p==0) continue;
      *p=0;p++;//split the line here
      while (*p==' ' || *p=='\t') p++;
      if (*p==0) continue;
      rt.Add(orig_id, new RefTran(p));
      } //while lines
}

void openfw(FILE* &f, GArgs& args, char opt) {
	GStr s=args.getOpt(opt);
	if (!s.is_empty()) {
		if (s=='-')
			f=stdout;
		else {
			f=fopen(s,"w");
			if (f==NULL) GError("Error creating file: %s\n", s.chars());
		}
	}
}

#define FWCLOSE(fh) if (fh!=NULL && fh!=stdout) fclose(fh)

void printGff3Header(FILE* f, GArgs& args) {
  if (gffloader.keepGff3Comments) {
	for (int i=0;i<gffloader.headerLines.Count();i++) {
		fprintf(f, "%s\n", gffloader.headerLines[i]);
	}
  } else {
    fprintf(f, "##gff-version 3\n");
    fprintf(f, "# gffread v" VERSION "\n");
    fprintf(f, "# ");args.printCmdLine(f);
  }
}

void printGSeqHeader(FILE* f, GenomicSeqData* gdata) {
if (f && gffloader.keepGff3Comments && gdata->seqreg_start>0 && gdata->seqreg_end>0)
	 fprintf(f, "##sequence-region %s %d %d\n", gdata->gseq_name,
			 gdata->seqreg_start, gdata->seqreg_end);

}

void processGffComment(const char* cmline, GfList* gflst) {
 if (cmline[0]!='#') return;
 const char* p=cmline;
 while (*p=='#') p++;
 GStr s(p);
 //this can be called only after gffloader initialization
 // so we can use gffloader.names->gseqs.addName()
 s.startTokenize("\t ", tkCharSet);
 GStr w;
  if (s.nextToken(w) && w=="sequence-region") {
	 GStr chr, wend;
	 if (s.nextToken(chr) && s.nextToken(w) && s.nextToken(wend)) {
		 int gseq_id=gffloader.names->gseqs.addName(chr.chars());
		 if (gseq_id>=0) {
			 GenomicSeqData* gseqdata=getGSeqData(g_data, gseq_id);
			 gseqdata->seqreg_start=w.asInt();
			 gseqdata->seqreg_end=wend.asInt();
		 }
		 else GError("Error adding ref seq ID %s\n", chr.chars());
	 }
	 return;
 }
 if (gflst->Count()==0) {
    //initial Gff3 header, store it
	char* hl=Gstrdup(cmline);
    gffloader.headerLines.Add(hl);
 }
}

void printGffObj(FILE* f, GffObj* gfo, GStr& locname, GffPrintMode exonPrinting, int& out_counter) {
    GffObj& t=*gfo;
    GTData* tdata=(GTData*)(t.uptr);
    if (tdata->replaced_by!=NULL || !T_PRINTABLE(t.udata)) return;
    //if (t.exons.Count()==0 && t.children.Count()==0 && forceExons)
    //  t.addExonSegment(t.start,t.end);
    T_NO_PRINT(t.udata);
    if (!fmtGFF3 && !gfo->isTranscript())
    	return; //only GFF3 prints non-transcript records (incl. parent genes)
    t.addAttr("locus", locname.chars());
    out_counter++;
    if (fmtGFF3) {
         //print the parent first, if any and if not printed already
         if (t.parent!=NULL && T_PRINTABLE(t.parent->udata)) {
             GTData* pdata=(GTData*)(t.parent->uptr);
             if (pdata && pdata->geneinfo!=NULL)
                  pdata->geneinfo->finalize();
             t.parent->addAttr("locus", locname.chars());
             t.parent->printGxf(f, exonPrinting, tracklabel, NULL, decodeChars);
             T_NO_PRINT(t.parent->udata);
         }
    }
    t.printGxf(f, exonPrinting, tracklabel, NULL, decodeChars);
}

void printAsTable(FILE* f, GffObj* gfo, int* out_counter=NULL) {
    GffObj& t=*gfo;
    GTData* tdata=(GTData*)(t.uptr);
    if (tdata->replaced_by!=NULL || !T_PRINTABLE(t.udata)) return;
    T_NO_PRINT(t.udata);
    if (out_counter!=NULL) (*out_counter)++;
	 //print the parent first, if any and if not printed already
	 if (t.parent!=NULL && T_PRINTABLE(t.parent->udata)) {
		 GTData* pdata=(GTData*)(t.parent->uptr);
		 if (pdata && pdata->geneinfo!=NULL)
			  pdata->geneinfo->finalize();
		 //t.parent->addAttr("locus", locname.chars());
		 //(*out_counter)++; ?
		 printTableData(f, *t.parent);
		 T_NO_PRINT(t.parent->udata);
	 }
    printTableData(f, *gfo);
}

void shutDown() {
	seqinfo.Clear();
	//if (faseq!=NULL) delete faseq;
	//if (gcdb!=NULL) delete gcdb;
	delete fltRange;
	delete fltJunction;
	FWCLOSE(f_out);
	FWCLOSE(f_w);
	FWCLOSE(f_x);
	FWCLOSE(f_y);
	FWCLOSE(f_j);
}

int main(int argc, char* argv[]) {
 GArgs args(argc, argv,
   "version;debug;merge;stream;adj-stop;bed;in-bed;tlf;in-tlf;cluster-only;nc;cov-info;help;"
    "sort-alpha;keep-genes;w-nocds;attrs=;w-add=;ids=;nids=;jmatch=;gtf;keep-comments;keep-exon-attrs;force-exons;t-adopt;gene2exon;"
    "ignore-locus;no-pseudo;table=sort-by=hvOUNHPWCVJMKQYTDARSZFGLEBm:g:i:r:s:l:t:o:w:x:y:j:d:");
 args.printError(USAGE, true);
 int numfiles = args.startNonOpt();
 if (args.getOpt("version")) {
    printf(VERSION"\n");
    exit(0);
 }
 if (args.getOpt('h') || args.getOpt("help") || ( numfiles==0 && !haveStdInput())) {
    GMessage("%s",USAGE);
    exit(1);
 }

 debugMode=(args.getOpt("debug")!=NULL);
 decodeChars=(args.getOpt('D')!=NULL);
 gffloader.forceExons=(args.getOpt("force-exons")!=NULL);
 gffloader.streamIn=(args.getOpt("stream")!=NULL);
 gffloader.noPseudo=(args.getOpt("no-pseudo")!=NULL);
 gffloader.ignoreLocus=(args.getOpt("ignore-locus")!=NULL);
 gffloader.transcriptsOnly=(args.getOpt('O')==NULL);
 //sortByLoc=(args.getOpt('S')!=NULL);
 addDescr=(args.getOpt('A')!=NULL);
 verbose=(args.getOpt('v')!=NULL || args.getOpt('E')!=NULL);
 wCDSonly=(args.getOpt('C')!=NULL);
 wNConly=(args.getOpt("nc")!=NULL);
 addCDSattrs=(args.getOpt('P')!=NULL);
 add_hasCDS=(args.getOpt("add-hasCDS")!=NULL);
 adjustStop=(args.getOpt("adj-stop")!=NULL);
 if (adjustStop) addCDSattrs=true;
 validCDSonly=(args.getOpt('V')!=NULL);
 altPhases=(args.getOpt('H')!=NULL);
 fmtGTF=(args.getOpt('T')!=NULL || args.getOpt("gtf")!=NULL); //switch output format to GTF
 fmtBED=(args.getOpt("bed")!=NULL); //BED output
 fmtTLF=(args.getOpt("tlf")!=NULL); //TLF output
 if (fmtGTF || fmtBED || fmtTLF) {
	 if (!gffloader.transcriptsOnly) {
		 GMessage("Error: option -O is only supported with GFF3 output");
		 exit(1);
	 }
	 fmtGFF3=false;
 }

 BEDinput=(args.getOpt("in-bed")!=NULL);
 TLFinput=(args.getOpt("in-tlf")!=NULL);
 bothStrands=(args.getOpt('B')!=NULL);
 fullCDSonly=(args.getOpt('J')!=NULL);
 spliceCheck=(args.getOpt('N')!=NULL);
 StarStop=(args.getOpt('S')!=NULL);

 gffloader.keepGenes=(args.getOpt("keep-genes")!=NULL);
 gffloader.trAdoption=(args.getOpt("t-adopt")!=NULL);
 gffloader.keepGff3Comments=(args.getOpt("keep-comments")!=NULL);
 gffloader.sortRefsAlpha=(args.getOpt("sort-alpha")!=NULL);
 if (args.getOpt("sort-by")!=NULL) {
	  if (gffloader.sortRefsAlpha)
		  GError("Error: options --sort-by and --sort-alpha are mutually exclusive!\n");
	 sortBy=args.getOpt("sort-by");
 }
 if (!sortBy.is_empty())
	   gffloader.loadRefNames(sortBy);
 gffloader.gene2exon=(args.getOpt("gene2exon")!=NULL);
 gffloader.matchAllIntrons=(args.getOpt('K')==NULL);
 gffloader.fuzzSpan=(args.getOpt('Q')!=NULL);
 gffloader.dOvlSET=(args.getOpt('Y')!=NULL);
 if (args.getOpt('M') || args.getOpt("merge")) {
	 gffloader.doCluster=true;
	 gffloader.collapseRedundant=true;
  } else {
    if (!gffloader.matchAllIntrons || gffloader.fuzzSpan || gffloader.dOvlSET) {
      GMessage("%s",USAGE);
      GMessage("Error: options -K,-Q,-Y require -M/--merge option!\n");
      exit(1);
    }
 }
 if (args.getOpt("cluster-only")) {
	 gffloader.doCluster=true;
	 gffloader.collapseRedundant=false;
    if (!gffloader.matchAllIntrons || gffloader.fuzzSpan || gffloader.dOvlSET) {
      GMessage("%s",USAGE);
      GMessage("Error: option -K,-Q,-Y have no effect with --cluster-only.\n");
      exit(1);
    }
 }
 if (gffloader.dOvlSET)
	 gffloader.fuzzSpan=true; //-Q enforced by -Y
 covInfo=(args.getOpt("cov-info"));
 if (covInfo) gffloader.doCluster=true; //need to collapse overlapping exons
 if (fullCDSonly) validCDSonly=true;
 if (verbose) {
     fprintf(stderr, "Command line was:\n");
     args.printCmdLine(stderr);
     }
 gffloader.fullAttributes=(args.getOpt('F')!=NULL);
 gffloader.keep_AllExonAttrs=(args.getOpt("keep-exon-attrs")!=NULL);
 if (gffloader.keep_AllExonAttrs && !gffloader.fullAttributes) {
	 GMessage("Error: option --keep-exon-attrs requires option -F !\n");
	 exit(0);
 }
 if (args.getOpt('G')==NULL)
	 gffloader.gatherExonAttrs=!gffloader.fullAttributes;
 else {
	   gffloader.gatherExonAttrs=true;
	   gffloader.fullAttributes=true;
 }
 if (gffloader.noPseudo && !gffloader.fullAttributes) {
	 gffloader.gatherExonAttrs=true;
	 gffloader.fullAttributes=true;
 }
 gffloader.ensemblProc=(args.getOpt('L')!=NULL);
 if (gffloader.ensemblProc) {
    gffloader.fullAttributes=true;
    gffloader.gatherExonAttrs=false;
    //sortByLoc=true;
 }

 tableFormat=args.getOpt("table");
 if (!tableFormat.is_empty()) {
	 setTableFormat(tableFormat);
	 fmtTable=true;
	 fmtGFF3=false;
	 gffloader.fullAttributes=true;
 }

 gffloader.mergeCloseExons=(args.getOpt('Z')!=NULL);

 multiExon=(args.getOpt('U')!=NULL);
 writeExonSegs=(args.getOpt('W')!=NULL);
 tracklabel=args.getOpt('t');

 if (args.getOpt('g'))
	   gfasta.init(args.getOpt('g'));
 //if (gfasta.fastaPath!=NULL)
 //    sortByLoc=true; //enforce sorting by chromosome/contig
 GStr s=args.getOpt('i');
 if (!s.is_empty()) maxintron=s.asInt();
 s=args.getOpt('l');
 if (!s.is_empty()) minLen=s.asInt();
 TFilters=(multiExon || wCDSonly || wNConly); //TODO: all transcript filters should be included here through validateGffRec()
 FILE* f_repl=NULL; //duplicate/collapsing info output file
 s=args.getOpt('d');
 if (!s.is_empty()) {
   if (s=="-") f_repl=stdout;
   else {
       f_repl=fopen(s.chars(), "w");
       if (f_repl==NULL) GError("Error creating file %s\n", s.chars());
   }
 }

 s=args.getOpt("attrs");
 if (!s.is_empty()) {
	 getAttrList(s);
	 gffloader.attrsFilter=(attrList.Count()>1);
	 gffloader.fullAttributes=true;
 }
 rfltWithin=(args.getOpt('R')!=NULL);
 char* sz=args.getOpt('r');
 if (sz) {
	fltRange=new GRangeParser(sz);
 	if (fltRange->end==0) //end coordinate not given
 		fltRange->end=UINT_MAX;
 } else {
   if (rfltWithin)
     GError("Error: option -R requires -r!\n");
 }
 sz=args.getOpt("jmatch");
 if (sz) {
	//TODO: check if this is a file?
	fltJunction=new GRangeParser(sz);
	if (fltJunction->strand=='.') fltJunction->strand=0;
 } //gseq/range filtering

 s=args.getOpt('m');
 if (!s.is_empty()) {
   FILE* ft=fopen(s,"r");
   if (ft==NULL) GError("Error opening reference table: %s\n",s.chars());
   loadRefTable(ft, reftbl);
   fclose(ft);
   }
 s=args.getOpt('s');
 if (!s.is_empty()) {
   FILE* fsize=fopen(s,"r");
   if (fsize==NULL) GError("Error opening info file: %s\n",s.chars());
   loadSeqInfo(fsize, seqinfo);
   fclose(fsize);
 }
 s=args.getOpt("ids");
 if (s.is_empty()) {
	 s=args.getOpt("nids");
	 if (!s.is_empty())
		 IDflt=idFlt_Exclude;
 } else {
	 IDflt=idFlt_Only;
 }
 if (!s.is_empty()) {
   FILE* f=fopen(s,"r");
   if (f==NULL) GError("Error opening ID list file: %s\n",s.chars());
   loadIDlist(f, fltIDs);
   if (fltIDs.Count()==0) {
	   GMessage("Warning: no IDs were loaded from file %s\n", s.chars());
	   IDflt=idFlt_None;
   }
   fclose(f);
 }
 openfw(f_out, args, 'o');
 //if (f_out==NULL) f_out=stdout;
 if (gfasta.fastaPath==NULL && (validCDSonly || spliceCheck || args.getOpt('w')!=NULL || args.getOpt('x')!=NULL || args.getOpt('y')!=NULL))
  GError("Error: -g option is required for options -w, -x, -y, -V, -N, -M !\n");
 openfw(f_w, args, 'w');
 openfw(f_x, args, 'x');
 openfw(f_y, args, 'y');
 openfw(f_j, args, 'j');
 s=args.getOpt("w-add");
 if (!s.is_empty()) {
	 if (f_w==NULL) GError("Error: --w-add option requires -w option!\n");
	 wPadding=s.asInt();
 }

 if (f_w!=NULL && args.getOpt("w-nocds"))
	 wfaNoCDS=true;

 if (f_out==NULL && f_w==NULL && f_x==NULL && f_y==NULL && f_j==NULL && !covInfo)
	 f_out=stdout;

 //if (f_y!=NULL || f_x!=NULL) wCDSonly=true;
 //useBadCDS=useBadCDS || (fgtfok==NULL && fgtfbad==NULL && f_y==NULL && f_x==NULL);

 //GList<GffObj> gfkept(false,true); //unsorted, free items on delete
 int out_counter=0; //number of records printed

 if (fmtGTF)
	 exonPrinting = gffloader.forceExons ? pgtfBoth : pgtfAny;
 else if (fmtBED)
	 exonPrinting=pgffBED;
 else if (fmtTLF)
	exonPrinting=pgffTLF;
 else { //printing regular GFF3
	exonPrinting = gffloader.forceExons ? pgffBoth : pgffAny;
 }

 while (true) {
   GStr infile;
   if (numfiles) {
      infile=args.nextNonOpt();
      if (infile.is_empty()) break;
      if (infile=="-") { f_in=stdin; infile="stdin"; }
           else  if ((f_in=fopen(infile, "r"))==NULL)
                    GError("Error: cannot open input file %s!\n",infile.chars());
                 else fclose(f_in);
      numfiles--;
   }
   else infile="-";

   const char* fext=getFileExt(infile.chars());
   if (BEDinput || (Gstricmp(fext, "bed")==0))
	   gffloader.BEDinput=true;
   if (TLFinput || (Gstricmp(fext, "tlf")==0))
	   gffloader.TLFinput=true;
   gffloader.openFile(infile);
   if (gffloader.streamIn) { //streaming in - disable all bulk load features
  	 gffloader.transcriptsOnly=true;
  	 gffloader.doCluster=false;
  	 covInfo=false;
   }

   gffloader.load(g_data, &processGffComment);
   if (gffloader.streamIn) {
	   //we're done, GffLoader::load() took care of everything
	   shutDown();
	   return 0;
   }


   // will also place the transcripts in loci, if doCluster is enabled
   if (gffloader.doCluster)
     collectLocusData(g_data, covInfo);
   if (numfiles==0) break;
 }
 if (covInfo) {
	 //report coverage info at STDOUT
	 uint64 f_bases=0;
	 uint64 r_bases=0;
	 uint64 u_bases=0;
	 for (int g=0;g<g_data.Count();g++) {
		 f_bases+=g_data[g]->f_bases;
		 r_bases+=g_data[g]->r_bases;
		 u_bases+=g_data[g]->u_bases;
	 }
	 fprintf(stdout, "Total bases covered by transcripts:\n");
	 if (f_bases>0) fprintf(stdout, "\t%" PRIu64 " on + strand\n", f_bases);
	 if (r_bases>0) fprintf(stdout, "\t%" PRIu64 " on - strand\n", r_bases);
	 if (u_bases>0) fprintf(stdout, "\t%" PRIu64 " on . strand\n", u_bases);
 }
 GStr loctrack("gffcl");
 if (tracklabel) loctrack=tracklabel;
 if (gffloader.sortRefsAlpha)
    g_data.setSorted(&gseqCmpName);
 bool firstGff3Print=fmtGFF3;
 if (gffloader.doCluster) {
   //grouped in loci
   for (int g=0;g<g_data.Count();g++) {
     GenomicSeqData* gdata=g_data[g];
     bool firstGSeqHeader=fmtGFF3;
     if (f_out && fmtGFF3 && gffloader.keepGff3Comments && gdata->seqreg_start>0)
    	 fprintf(f_out, "##sequence-region %s %d %d\n", gdata->gseq_name,
    			 gdata->seqreg_start, gdata->seqreg_end);
     for (int l=0;l<gdata->loci.Count();l++) {
       bool firstLocusPrint=true;
       GffLocus& loc=*(gdata->loci[l]);
       //check all non-replaced transcripts in this locus:
       //int numvalid=0;
       //int idxfirstvalid=-1;
       for (int i=0;i<loc.rnas.Count();i++) {
         GffObj& t=*(loc.rnas[i]);
         GTData* tdata=(GTData*)(t.uptr);
         if (tdata->replaced_by!=NULL) {
            if (f_repl && T_DUPSHOWABLE(t.udata)) {
               fprintf(f_repl, "%s", t.getID());
               GTData* rby=tdata;
               while (rby->replaced_by!=NULL) {
                  fprintf(f_repl," => %s", rby->replaced_by->getID());
                  T_NO_DUPSHOW(rby->rna->udata);
                  rby=(GTData*)(rby->replaced_by->uptr);
               }
               fprintf(f_repl, "\n");
            }
            T_NO_PRINT(t.udata);
            if (verbose) {
            	GMessage("Info: %s discarded: superseded by %s\n",
            			t.getID(), tdata->replaced_by->getID());
            }
            continue;
         }
         //restore strand for dOvlSET
         char orig_strand=T_OSTRAND(t.udata);
         if (orig_strand!=0) t.strand=orig_strand;
         /* -- transcripts are filtered upon loading
         if (process_transcript(gfasta, t)) {
             numvalid++;
             if (idxfirstvalid<0) idxfirstvalid=i;
         }
         */
       } //for each transcript

       int rnas_i=0;
       //if (idxfirstvalid>=0) rnas_i=idxfirstvalid;
       int gfs_i=0;
       if (f_out) {
           GStr locname("RLOC_");
           locname.appendfmt("%08d",loc.locus_num);
           //GMessage("Locus: %s (%d-%d), %d rnas, %d gfs\n", locname.chars(), loc.start, loc.end,
           //	   loc.rnas.Count(), loc.gfs.Count());
		   while (gfs_i<loc.gfs.Count() || rnas_i<loc.rnas.Count()) {
			   if (gfs_i<loc.gfs.Count() && (rnas_i>=loc.rnas.Count() ||
					     loc.gfs[gfs_i]->start<=loc.rnas[rnas_i]->start) ) {
				   //print the gene object first
				   if (fmtGFF3) { //BED, TLF and GTF: only show transcripts
					   if (firstGff3Print) { printGff3Header(f_out, args);firstGff3Print=false; }
					   if (firstGSeqHeader) { printGSeqHeader(f_out, gdata); firstGSeqHeader=false; }
					   if (firstLocusPrint) {
						   //loc.print(f_out, idxfirstvalid, locname, loctrack);
						   loc.print(f_out, 0, locname, loctrack);
						   firstLocusPrint=false;
					   }
					   printGffObj(f_out, loc.gfs[gfs_i], locname, exonPrinting, out_counter);
				   }
				   ++gfs_i;
				   continue;
			   }
			   if (rnas_i<loc.rnas.Count()) {
				       //loc.rnas[rnas_i]->printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
				       if (fmtGFF3) {
					     if (firstGff3Print) { printGff3Header(f_out, args); firstGff3Print=false; }
					     if (firstGSeqHeader) { printGSeqHeader(f_out, gdata); firstGSeqHeader=false; }
					     if (firstLocusPrint) {
					    	 //loc.print(f_out, idxfirstvalid, locname, loctrack);
					    	 loc.print(f_out, 0, locname, loctrack);
					    	 firstLocusPrint=false;
					     }
				       }
				       if (fmtTable)  printAsTable(f_out, loc.rnas[rnas_i], &out_counter);
				       else printGffObj(f_out, loc.rnas[rnas_i], locname, exonPrinting, out_counter);
					   ++rnas_i;
			   }
		   }
       }
     }//for each locus
    } //for each genomic sequence
   } //if Clustering enabled
  else { //no clustering
   //not grouped into loci, print the rnas with their parents, if any
   //int numvalid=0;
   for (int g=0;g<g_data.Count();g++) {
     GenomicSeqData* gdata=g_data[g];
     bool firstGSeqHeader=fmtGFF3;
     int gfs_i=0;
     for (int m=0;m<gdata->rnas.Count();m++) {
        GffObj& t=*(gdata->rnas[m]);
        if (f_out && (fmtGFF3 || fmtTable)) {
         //print other non-transcript (gene?) feature that might be there before t
           while (gfs_i<gdata->gfs.Count() && gdata->gfs[gfs_i]->start<=t.start) {
             GffObj& gfst=*(gdata->gfs[gfs_i]);
             if (TFilters && gfst.isGene() && gfst.children.Count()==0) // gene with no children left, skip it if filters were applied
            	 { ++gfs_i; continue; }
             if T_PRINTABLE(gfst.udata) { //never printed
               T_NO_PRINT(gfst.udata);
               if (fmtGFF3) {
                 if (firstGff3Print) { printGff3Header(f_out, args);firstGff3Print=false; }
                 if (firstGSeqHeader) { printGSeqHeader(f_out, gdata); firstGSeqHeader=false; }
                 gfst.printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
               }
               else printTableData(f_out, gfst);
             }
             ++gfs_i;
           }
        }
        GTData* tdata=(GTData*)(t.uptr);
        if (tdata->replaced_by!=NULL) continue;
        //if (process_transcript(gfasta, t)) {
        //   numvalid++;
           if (f_out && T_PRINTABLE(t.udata) ) {
             T_NO_PRINT(t.udata);
             if (fmtGFF3 || fmtTable || t.isTranscript()) {
				 if (tdata->geneinfo)
					 tdata->geneinfo->finalize();
				 out_counter++;
				 if (fmtGFF3) {
				   if (firstGff3Print) { printGff3Header(f_out, args);firstGff3Print=false; }
				   if (firstGSeqHeader) { printGSeqHeader(f_out, gdata); firstGSeqHeader=false; }
				 }
				 //for GFF3 && table output, print the parent first, if any
				 if ((fmtGFF3 || fmtTable) && t.parent!=NULL && T_PRINTABLE(t.parent->udata)) {
					 //GTData* pdata=(GTData*)(t.parent->uptr);
					 //if (pdata && pdata->geneinfo!=NULL)
					 //  pdata->geneinfo->finalize();
					 if (fmtTable)
						 printTableData(f_out, *(t.parent));
					 else { //GFF3 output
						 t.parent->printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
					 }
					 T_NO_PRINT(t.parent->udata);
				 }
				 if (fmtTable)
					 printTableData(f_out, t);
				 else
					 t.printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
             }
           }//GFF/GTF output requested
        //} //valid transcript
     } //for each rna
     //print the rest of the isolated pseudo/gene/region features not printed yet
     if (f_out && (fmtGFF3 || fmtTable)) {
      while (gfs_i<gdata->gfs.Count()) {
         GffObj& gfst=*(gdata->gfs[gfs_i]);
         if (TFilters && gfst.isGene() && gfst.children.Count()==0) // gene with no children left, skip it if filters were applied
        	 { ++gfs_i; continue; }
         if T_PRINTABLE(gfst.udata) { //never printed
           T_NO_PRINT(gfst.udata);
           if (fmtGFF3) {
              if (firstGff3Print) { printGff3Header(f_out, args); firstGff3Print=false; }
              if (firstGSeqHeader) { printGSeqHeader(f_out, gdata); firstGSeqHeader=false; }
              gfst.printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
           } else
              printTableData(f_out, gfst);
         }
         ++gfs_i;
      }
     }
    } //for each genomic seq
   } //no clustering
 if (f_repl && f_repl!=stdout) fclose(f_repl);
 shutDown();
 }


