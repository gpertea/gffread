#include "GArgs.h"
#include "gff_utils.h"
#include <ctype.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#define VERSION "0.11.6"

#define USAGE "gffread v" VERSION ". Usage:\n\
gffread <input_gff> [-g <genomic_seqs_fasta> | <dir>][-s <seq_info.fsize>] \n\
 [-o <outfile>] [-t <trackname>] [-r [[<strand>]<chr>:]<start>..<end> [-R]]\n\
 [-CTVNJMKQAFPGUBHZWTOLE] [-w <exons.fa>] [-x <cds.fa>] [-y <tr_cds.fa>]\n\
 [-i <maxintron>] [--bed] [--table <attrlist>] [--sort-by <refseq_list.txt>]\n\
 \n\
 Filter, convert or cluster GFF/GTF/BED records, extract the sequence of\n\
 transcripts (exon or CDS) and more.\n\
 By default (i.e. without -O) only transcripts are processed, discarding any\n\
 other non-transcript features. Default output is a simplified GFF3 with only\n\
 the basic attributes.\n\
 \n\
 <input_gff> is a GFF file, use '-' for stdin\n\
 \n\
Options:\n\
 -i   discard transcripts having an intron larger than <maxintron>\n\
 -l   discard transcripts shorter than <minlen> bases\n\
 -r   only show transcripts overlapping coordinate range <start>..<end>\n\
      (on chromosome/contig <chr>, strand <strand> if provided)\n\
 -R   for -r option, discard all transcripts that are not fully \n\
      contained within the given range\n\
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
 -F   preserve all GFF attributes (for non-exon features)\n\
 --keep-exon-attrs : for -F option, do not attempt to reduce redundant\n\
      exon/CDS attributes\n\
 -G   do not keep exon attributes, move them to the transcript feature\n\
      (for GFF3 output)\n\
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
Clustering:\n\
 -M/--merge : cluster the input transcripts into loci, discarding\n\
      \"duplicated\" transcripts (those with the same exact introns\n\
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
 -w    write a fasta file with spliced exons for each transcript\n\
 --w-add <N> for the -w option, extract additional <N> bases\n\
       both upstream and downstream of the transcript boundaries\n\
 -x    write a fasta file with spliced CDS for each GFF transcript\n\
 -y    write a protein fasta file with the translation of CDS for each record\n\
 -W    for -w and -x options, write in the FASTA defline the exon\n\
       coordinates projected onto the spliced sequence;\n\
       for -y option, write transcript attributes in the FASTA defline\n\
 -S    for -y option, use '*' instead of '.' as stop codon translation\n\
 -L    Ensembl GTF to GFF3 conversion (implies -F; should be used with -m)\n\
 -m    <chr_replace> is a name mapping table for converting reference \n\
       sequence names, having this 2-column format:\n\
       <original_ref_ID> <new_ref_ID>\n\
 -t    use <trackname> in the 2nd column of each GFF/GTF output line\n\
 -o    write the records into <outfile> instead of stdout\n\
 -T    main output will be GTF instead of GFF3\n\
 --bed output records in BED format instead of default GFF3\n\
 --tlf output \"transcript line format\" which is like GFF\n\
       but exons, CDS features and related data are stored as GFF \n\
       attributes in the transcript feature line, like this:\n\
         exoncount=N;exons=<exons>;CDSphase=<N>;CDS=<CDScoords> \n\
       <exons> is a comma-delimited list of exon_start-exon_end coordinates;\n\
       <CDScoords> is CDS_start:CDS_end coordinates or a list like <exons>\n\
 --table output a simple tab delimited format instead of GFF, with columns\n\
       having the values of GFF attributes given in <attrlist>; special\n\
       pseudo-attributes (prefixed by @) are recognized:\n\
       @chr, @start, @end, @strand, @numexons, @exons, @cds, @covlen, @cdslen\n\
 -v,-E expose (warn about) duplicate transcript IDs and other potential\n\
       problems with the given GFF/GTF records\n\
"

class SeqInfo { //populated from the -s option of gffread
 public:
  int len;
  char* descr;
  SeqInfo( int l, char* s): len(l), descr(NULL) {
    if (s!=NULL)
      descr=Gstrdup(s);
  }
  ~SeqInfo() {
    GFREE(descr);
  }
};

class RefTran {
 public:
   char* new_name;
   RefTran(char *ns) {
      new_name=NULL;
      if (ns!=NULL)
         new_name=Gstrdup(ns);
      }
   ~RefTran() {
      GFREE(new_name);
      }
};

enum ETableFieldType {
  ctfGFF_Attr=0, // attribute name as is
  ctfGFF_ID, //ID or @id
  ctfGFF_Parent, //Parent or @parent
  ctfGFF_chr, //@chr
  ctfGFF_feature, //@feature
  ctfGFF_start, //@start
  ctfGFF_end, //@end
  ctfGFF_strand, //@strand
  ctfGFF_numexons, //@numexons
  ctfGFF_exons, //@exons
  ctfGFF_cds, //@cds
  ctfGFF_covlen, //@covlen
  ctfGFF_cdslen//@cdslen
};

class CTableField {
 public:
   ETableFieldType type;
   GStr name; //only for type ctfGFF_Attr
   CTableField(ETableFieldType atype=ctfGFF_Attr):type(atype) { }
   CTableField(GStr& attrname):type(ctfGFF_Attr),name(attrname) { }
};

FILE* ffasta=NULL;
FILE* f_in=NULL;
FILE* f_out=NULL;
FILE* f_w=NULL; //writing fasta with spliced exons (transcripts)
int wPadding = 0; //padding for -w option
FILE* f_x=NULL; //writing fasta with spliced CDS
FILE* f_y=NULL; //wrting fasta with translated CDS

bool wCDSonly=false;
bool wNConly=false;
int minLen=0; //minimum transcript length
bool validCDSonly=false; // translation with no in-frame STOP
bool bothStrands=false; //for single-exon mRNA validation, check the other strand too
bool altPhases=false; //if original phase fails translation validation,
                     //try the other 2 phases until one makes it
bool addCDSattrs=false;
bool add_hasCDS=false;
bool adjustStop=false; //automatic adjust the CDS stop coordinate
bool covInfo=false; // --cov-info : only report genome coverage
//bool transcriptsOnly=true;
//bool keepGenes=false; //for transcriptsOnly

//bool sortAlpha=false;
GStr sortBy; //file name with chromosomes listed in the desired order
//bool keepRefOrder=false; //sort within chromosomes, but follow the input chromosome order -- default!
GStr tableFormat; //list of "attributes" to print in tab delimited format
//bool NoPseudo=false;
bool spliceCheck=false; //only known splice-sites
bool decodeChars=false; //decode url-encoded chars in attrs (-D)
bool StarStop=false; //use * instead of . for stop codon translation

bool fullCDSonly=false; // starts with START, ends with STOP codon
//bool fullattr=false; //-F
//bool gatherExonAttrs=false; //-G

//bool sortByLoc=false; // if the GFF output should be sorted by location
bool ensembl_convert=false; //-L, assist in converting Ensembl GTF to GFF3
bool BEDinput=false;
bool TLFinput=false;

bool fmtGFF3=true; //default output: GFF3
//other formats only make sens in transcriptOnly mode
bool fmtGTF=false;
bool fmtBED=false;
bool fmtTLF=false;
bool fmtTable=false;
bool addDescr=false;
//bool protmap=false;
bool multiExon=false;
bool writeExonSegs=false;
char* tracklabel=NULL;
int maxintron=999000000;
//bool mergeCloseExons=false;
//range filter:
char* rfltGSeq=NULL;
char rfltStrand=0;
uint rfltStart=0;
uint rfltEnd=MAX_UINT;
bool rfltWithin=false; //check for full containment within given range

GffLoader gffloader;

GList<GenomicSeqData> g_data(true,true,true); //list of GFF records by genomic seq

//hash with sequence info
GHash<SeqInfo> seqinfo;
GHash<int> isoCounter; //counts the valid isoforms
GHash<RefTran> reftbl;
GHash<GeneInfo> gene_ids;
  //min-max gene span associated to chr|gene_id (mostly for Ensembl conversion)

bool debugMode=false;
//bool verbose=false;


GVec<CTableField> tableCols; //table output format fields

void loadSeqInfo(FILE* f, GHash<SeqInfo> &si) {
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

void setTableFormat(GStr& s) {
	 if (s.is_empty()) return;
	 GHash<ETableFieldType> specialFields;
	 specialFields.Add("chr", new ETableFieldType(ctfGFF_chr));
	 specialFields.Add("id", new ETableFieldType(ctfGFF_ID));
	 specialFields.Add("parent", new ETableFieldType(ctfGFF_Parent));
	 specialFields.Add("feature", new ETableFieldType(ctfGFF_feature));
	 specialFields.Add("start", new ETableFieldType(ctfGFF_start));
	 specialFields.Add("end", new ETableFieldType(ctfGFF_end));
	 specialFields.Add("strand", new ETableFieldType(ctfGFF_strand));
	 specialFields.Add("numexons", new ETableFieldType(ctfGFF_numexons));
	 specialFields.Add("exons", new ETableFieldType(ctfGFF_exons));
	 specialFields.Add("cds", new ETableFieldType(ctfGFF_cds));
	 specialFields.Add("covlen", new ETableFieldType(ctfGFF_covlen));
	 specialFields.Add("cdslen", new ETableFieldType(ctfGFF_cdslen));

	 s.startTokenize(" ,;.:", tkCharSet);
	 GStr w;
	 while (s.nextToken(w)) {
      if (w[0]=='@') {
    	  w=w.substr(1);
    	  ETableFieldType* v=specialFields.Find(w.chars());
    	  if (v!=NULL) {
    		  CTableField tcol(*v);
    		  tableCols.Add(tcol);
    	  }
    	  else GMessage("Warning: table field '@%s' not recognized!\n",w.chars());
    	  continue;
      }
      if (w=="ID") {
    	  CTableField tcol(ctfGFF_ID);
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

void loadRefTable(FILE* f, GHash<RefTran>& rt) {
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

char* getSeqDescr(char* seqid) {
 static char charbuf[128];
 if (seqinfo.Count()==0) return NULL;
 char* suf=rstrchr(seqid, '.');
 if (suf!=NULL) *suf=0;
 SeqInfo* seqd=seqinfo.Find(seqid);
 if (suf!=NULL) *suf='.';
 if (seqd!=NULL) {
  GStr s(seqd->descr);
  //cleanup some Uniref gunk
  if (s[0]=='[') {
    int r=s.index(']');
    if (r>=0 && r<8 && isdigit(s[1]))
       s.remove(0,r+1);
    }
  if (s.length()>80) {
    int r=s.index(';');
    if (r>5) s.cut(r);
    }
  if (s.length()>127) {
   s.cut(127);
   int r=s.rindex(' ');
   if (r>0) s.cut(r);
   }
  strcpy(charbuf, s.chars());
  return charbuf;
  }
 else return NULL;
}

char* getSeqName(char* seqid) {
  static char charbuf[128];
  char* suf=rstrchr(seqid, '.');
  if (suf!=NULL) *suf=0;
  strcpy(charbuf, seqid);
  if (suf!=NULL) *suf='.';
  return charbuf;
}


int adjust_stopcodon(GffObj& gffrec, int adj, GList<GSeg>* seglst=NULL) {
  //adj>0, extend CDS to include a potential stop codon
  //when CDS is expanded, the terminal exon might have to be adjusted too
  int realadj=0;
  if (gffrec.strand=='-') {
       if ((int)gffrec.CDstart>adj) {
           gffrec.CDstart-=adj;
           realadj=adj;
           if (gffrec.exons.First()->start>gffrec.CDstart) {
                 gffrec.covlen+=gffrec.exons.First()->start - gffrec.CDstart;
                 gffrec.exons.First()->start=gffrec.CDstart;
                 gffrec.start=gffrec.CDstart;
                 }
             }
          }
        else { // forward strand
         //expand beyond
         realadj=adj;
         gffrec.CDend+=adj;
         if (adj<0) {//restore
           if (gffrec.exons.Last()->end==gffrec.CDend-adj) {
                        gffrec.exons.Last()->end+=adj;
                        gffrec.end=gffrec.exons.Last()->end;
                        gffrec.covlen+=adj;
                        }
         }
         else if (gffrec.exons.Last()->end<gffrec.CDend) {
             gffrec.covlen+=gffrec.CDend-gffrec.exons.Last()->end;
             gffrec.exons.Last()->end=gffrec.CDend;
             gffrec.end=gffrec.CDend;
             }
         }
  if (seglst!=NULL) seglst->Last()->end+=realadj;
  return realadj;
 }

bool process_transcript(GFastaDb& gfasta, GffObj& gffrec) {
 if (!gffrec.isTranscript()) return false; //shouldn't call this function unless it's a transcript
 //returns true if the transcript passed the filter
 char* gname=gffrec.getGeneName();
 if (gname==NULL) gname=gffrec.getGeneID();
 if (ensembl_convert && startsWith(gffrec.getID(), "ENS")) {
      const char* biotype=gffrec.getAttr("gene_biotype");
      if (biotype) {
         gffrec.addAttr("type", biotype);
         gffrec.removeAttr("gene_biotype");
         }
       else { //old Ensembl files lacking gene_biotype
         gffrec.addAttr("type", gffrec.getTrackName());
         }

      //bool is_gene=false;
      bool is_pseudo=false;
      if (strcmp(biotype, "protein_coding")==0 || gffrec.hasCDS())
                gffrec.setFeatureName("mRNA");
       else {
          if (strcmp(biotype, "processed_transcript")==0)
              gffrec.setFeatureName("proc_RNA");
            else {
              //is_gene=endsWith(biotype, "gene");
              is_pseudo=strifind(biotype, "pseudo");
              if (is_pseudo) {
                   gffrec.setFeatureName("pseudo_RNA");
                   }
                else if (endsWith(biotype, "RNA")) {
                   gffrec.setFeatureName(biotype);
                   } else gffrec.setFeatureName("misc_RNA");
              }
          }
      }
 if (gname && strcmp(gname, gffrec.getID())!=0) {
   int* isonum=isoCounter.Find(gname);
   if  (isonum==NULL) {
       isonum=new int(1);
       isoCounter.Add(gname,isonum);
       }
      else (*isonum)++;
   //defline.appendfmt(" gene=%s", gname);
   }
  int seqlen=0;

  const char* tlabel=tracklabel;
  if (tlabel==NULL) tlabel=gffrec.getTrackName();
  //defline.appendfmt(" track:%s",tlabel);
  char* cdsnt = NULL;
  char* cdsaa = NULL;
  int aalen=0;
  for (int i=1;i<gffrec.exons.Count();i++) {
     int ilen=gffrec.exons[i]->start-gffrec.exons[i-1]->end-1;
     if (verbose && ilen>4000000)
            GMessage("Warning: very large intron (%d) for transcript %s\n",
                           ilen, gffrec.getID());
     if (ilen>maxintron) {
         return false;
     }
  }
  GMapSegments seglst(gffrec.strand);
  GFaSeqGet* faseq=NULL;
  if (f_x!=NULL || f_y!=NULL || f_w!=NULL || spliceCheck || validCDSonly || addCDSattrs) {
	  faseq=fastaSeqGet(gfasta, gffrec.getGSeqName());
      if (faseq==NULL)
	    	GError("Error: no genomic sequence available (check -g option!).\n");
  }
  if (spliceCheck && gffrec.exons.Count()>1) {
    //check introns for splice site consensi ( GT-AG, GC-AG or AT-AC )
    int glen=gffrec.end-gffrec.start+1;
    const char* gseq=faseq->subseq(gffrec.start, glen);
    bool revcompl=(gffrec.strand=='-');
    bool ssValid=true;
    for (int e=1;e<gffrec.exons.Count();e++) {
      const char* intron=gseq+gffrec.exons[e-1]->end+1-gffrec.start;
      int intronlen=gffrec.exons[e]->start-gffrec.exons[e-1]->end-1;
      GSpliceSite acceptorSite(intron,intronlen,true, revcompl);
      GSpliceSite    donorSite(intron,intronlen, false, revcompl);
      //GMessage("%c intron %d-%d : %s .. %s\n",
      //           gffrec.strand, istart, iend, donorSite.nt, acceptorSite.nt);
      if (acceptorSite=="AG") { // GT-AG or GC-AG
         if (!donorSite.canonicalDonor()) {
            ssValid=false;break;
            }
         }
      else if (acceptorSite=="AC") { //AT-AC also accepted
         if (donorSite!="AT") { ssValid=false; break; }
         }
      else { ssValid=false; break; }
      }
    if (!ssValid) {
      if (verbose)
         GMessage("Unrecognized splice sites found for '%s'\n",gffrec.getID());
      return false; //don't print this one!
    }
  }
  bool trprint=true;
  bool inframeStop=false;
  //int stopCodonAdjust=0;
  int mCDphase=0;
  bool fullCDS=false;
  bool endStop=false;
  bool stopAdjusted=false;
  if (add_hasCDS && gffrec.hasCDS()) gffrec.addAttr("hasCDS", "true");
  if (gffrec.CDphase=='1' || gffrec.CDphase=='2')
      mCDphase = gffrec.CDphase-'0';
  //CDS partialness only added when -y -x -V options are given
  if (gffrec.hasCDS() && (f_y!=NULL || f_x!=NULL || validCDSonly || addCDSattrs)) {
    int strandNum=0;
    int phaseNum=0;
  CDS_CHECK:
    uint cds_olen=0;
    cdsnt=gffrec.getSpliced(faseq, true, &seqlen, NULL, &cds_olen, &seglst, adjustStop);
    //if adjustStop, seqlen has the CDS+3'UTR length, but cds_olen still has the original CDS length
    if (cdsnt!=NULL && cdsnt[0]!='\0') { //has CDS
         cdsaa=translateDNA(cdsnt, aalen, seqlen);
         char* p=strchr(cdsaa,'.');
         int cds_aalen=aalen;
         if (adjustStop)
        	 cds_aalen=cds_olen/3; //originally stated CDS length
         endStop=false;
         if (p!=NULL) { //stop codon found
        	 if (p-cdsaa==cds_aalen-1) { //stop found as the stated last CDS codon
                  *p='\0';//remove it
                  endStop=true;
                  if (adjustStop) {
                	  seqlen=cds_aalen*3;
                	  aalen=cds_aalen;
                  }
                  cds_aalen--;
                  aalen--;
                  //no need to adjust stop codon
              }
              else {//stop found in a different position than the last codon
            	  if (p-cdsaa<cds_aalen-1 && !adjustStop) {
            		  inframeStop=true;
            	  }
            	  if (adjustStop) {
            		  *p='\0';
            		  cds_aalen=p-cdsaa+1; //adjusted CDS length
            		  seqlen=cds_aalen*3;
            		  aalen=cds_aalen;
            		  uint gc=seglst.gmap(seqlen);
            		  if (gffrec.strand=='-') gffrec.CDstart=gc;
            		  else gffrec.CDend=gc;
            		  endStop=true;
            		  stopAdjusted=true;
            	  }
              }
         }//stop codon found
         //if (trprint==false) { //failed CDS validity check
         if (inframeStop) {
           //in-frame stop codon found
           if (altPhases && phaseNum<3) {
              phaseNum++; //try a different phase
              gffrec.CDphase = '0'+((mCDphase+phaseNum)%3);
              GFREE(cdsaa);
              goto CDS_CHECK;
           }
           if (gffrec.exons.Count()==1 && bothStrands) {
              strandNum++;
              phaseNum=0;
              if (strandNum<2) {
                 GFREE(cdsaa);
                 gffrec.strand = (gffrec.strand=='-') ? '+':'-';
                 goto CDS_CHECK; //repeat the CDS check for a different frame
              }
           }
           if (verbose) GMessage("Warning: In-frame STOP found for '%s'\n",gffrec.getID());
           if (addCDSattrs) gffrec.addAttr("InFrameStop", "true");
         } //has in-frame STOP
         if (stopAdjusted) {
      	   if (addCDSattrs) gffrec.addAttr("CDStopAdjusted", "true");
      	   inframeStop=false; //pretend it's OK now that we've adjusted it
         }
         if (!inframeStop) {
			 bool hasStart=(cdsaa[0]=='M'); //for the regular eukaryotic translation table
			 fullCDS=(endStop && hasStart);
			 if (!fullCDS) {
				 const char* partialness=NULL;
				 if (hasStart) partialness="3";
				 else {
					partialness = endStop ? "5" : "5_3";
				 }
				 if (addCDSattrs) gffrec.addAttr("partialness", partialness);
			 }
         }
         if (trprint && ((fullCDSonly && !fullCDS) || (validCDSonly && inframeStop)) )
        	 trprint=false;
         //} // Valid CDS only requested?
      } //has CDS
  } //translation or codon check was requested
  if (!trprint) {
    GFREE(cdsnt);
    GFREE(cdsaa);
    //if (adjstop!=NULL) delete adjstop;
    return false;
  }
  /*
  if (validCDSonly) {
     int stopCodonAdjust=adjstop->restore();
     if (stopCodonAdjust!=0 && !endStop) {
        //restore stop codon location
        //adjust_stopcodon(gffrec, -stopCodonAdjust, &seglst);
	    if (seglst.Count()>0) seglst.Last()->end-=stopCodonAdjust;
        if (cdsnt!=NULL && seqlen>0) {
           seqlen-=stopCodonAdjust;
           cdsnt[seqlen]=0;
        }
        if (cdsaa!=NULL) aalen--;
     }
  }
  if (adjstop!=NULL) delete adjstop;
  */
  if (cdsnt!=NULL) { // && !inframeStop) {
	  if (f_y!=NULL) { //CDS translation fasta output requested
			 if (cdsaa==NULL) { //translate now if not done before
			   cdsaa=translateDNA(cdsnt, aalen, seqlen);
			 }
			 GStr defline(gffrec.getID());
			 if (gffrec.attrs!=NULL) {
				 //append all attributes found for each transcripts
				for (int i=0;i<gffrec.attrs->Count();i++) {
				  defline.append(" ");
				  defline.append(gffrec.getAttrName(i));
				  defline.append("=");
				  char* s=gffrec.getAttrValue(i);
				  if (s[0]=='"') defline.append(s);
				  else defline.appendQuoted(s, '{', true);
				}
			 }
			 if (aalen>0) {
			   if (cdsaa[aalen-1]=='.' || cdsaa[aalen-1]=='\0') --aalen; //avoid printing the stop codon
			   printFasta(f_y, defline, cdsaa, aalen, StarStop);
			 }
	  }
	  if (f_x!=NULL) { //CDS only
			 GStr defline(gffrec.getID());
			 if (writeExonSegs) {
				  defline.append(" loc:");
				  defline.append(gffrec.getGSeqName());
				  defline.appendfmt("(%c)",gffrec.strand);
				  //warning: not CDS coordinates are written here, but the exon ones
				  defline+=(int)gffrec.start;
				  defline+=(char)'-';
				  defline+=(int)gffrec.end;
				  // -- here these are CDS substring coordinates on the spliced sequence:
				  defline.append(" segs:");
				  for (int i=0;i<seglst.Count();i++) {
					  if (i>0) defline.append(",");
					  defline+=(int)seglst[i].start;
					  defline.append("-");
					  defline+=(int)seglst[i].end;
					  }
				  }
			 if (gffrec.attrs!=NULL) {
				 //append all attributes found for each transcript
				for (int i=0;i<gffrec.attrs->Count();i++) {
					defline.append(" ");
					defline.append(gffrec.getAttrName(i));
					defline.append("=");
					char* s=gffrec.getAttrValue(i);
					if (s[0]=='"') defline.append(s);
					else defline.appendQuoted(s, '{', true);
				}
			}
			printFasta(f_x, defline, cdsnt, seqlen);
	  }
	  GFREE(cdsnt);
	  GFREE(cdsaa);
  } //writing CDS or its translation
  if (f_w!=NULL) { //write spliced exons
	  uint cds_start=0;
	  uint cds_end=0;
	  seglst.Clear();
	  //TODO: ? if wPadding is set, *temporarily* change first and last exon coordinates ?!?
	  // or perhaps getSpliced() should take an additional padding parameter ?!?
	  int padLeft=0;
	  int padRight=0;
	  if (wPadding>0) {
		padLeft= (gffrec.start>(uint)wPadding) ? wPadding : gffrec.start - 1;
		int ediff=faseq->getseqlen()-gffrec.end;
	    padRight=(wPadding>ediff) ?  ediff : wPadding;
   	    gffrec.addPadding(padLeft, padRight);
	  }
	  char* exont=gffrec.getSpliced(faseq, false, &seqlen, &cds_start, &cds_end, &seglst);
	  //restore exons to normal (remove padding)
	  if (wPadding>0)
		  gffrec.removePadding(padLeft, padRight);

	  GStr defline(gffrec.getID());
	  if (exont!=NULL) {
		  if (gffrec.CDstart>0) {
			  defline.appendfmt(" CDS=%d-%d", cds_start, cds_end);
		  }
		  if (writeExonSegs) {
			  defline.append(" loc:");
			  defline.append(gffrec.getGSeqName());
			  defline+=(char)'|';
			  defline+=(int)gffrec.start;
			  defline+=(char)'-';
			  defline+=(int)gffrec.end;
			  defline+=(char)'|';
			  defline+=(char)gffrec.strand;
			  defline.append(" exons:");
			  for (int i=0;i<gffrec.exons.Count();i++) {
				  if (i>0) defline.append(",");
				  defline+=(int)gffrec.exons[i]->start;
				  defline.append("-");
				  defline+=(int)gffrec.exons[i]->end;
			  }
			if (wPadding>0) {
				defline.append(" padding:");
				defline.append(padLeft);
				defline+=(char)'|';
				defline.append(padRight);
			}

			defline.append(" segs:");
			for (int i=0;i<seglst.Count();i++) {
				if (i>0) defline.append(",");
				defline+=(int)seglst[i].start;
				defline.append("-");
				defline+=(int)seglst[i].end;
				}
		  }

		  if (gffrec.attrs!=NULL) {
			  //append all attributes found for each transcripts
			  for (int i=0;i<gffrec.attrs->Count();i++) {
				  defline.append(" ");
				  defline.append(gffrec.getAttrName(i));
				  defline.append("=");
				  char* s=gffrec.getAttrValue(i);
				  if (s[0]=='"') defline.append(s);
				  else defline.appendQuoted(s, '{', true);
			  }
		  }
		  printFasta(f_w, defline, exont, seqlen);
		  GFREE(exont);
	  }
  } //writing f_w (spliced exons)
  return true;
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
    fprintf(f, "# ");
    args.printCmdLine(f);
    fprintf(f, "# gffread v" VERSION "\n");
    fprintf(f, "##gff-version 3\n");
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

bool validateGffRec(GffObj* gffrec, GList<GffObj>* gfnew) {
	if (reftbl.Count()>0) { //check if we need to reject by ref seq filter
		GStr refname(gffrec->getRefName());
		RefTran* rt=reftbl.Find(refname.chars());
		if (rt==NULL && refname.length()>2 && refname[-2]=='.' && isdigit(refname[-1])) {
			//try removing the version suffix
			refname.cut(-2);
			//GMessage("[DEBUG] Trying ref name '%s'...\n", refname.chars());
			rt=reftbl.Find(refname.chars());
		}
		if (rt) {
			gffrec->setRefName(rt->new_name);
		}
		/* //no, do not discard non-matching entries, let them pass through!
		else {
			if (verbose)
				GMessage("Info: %s discarded due to reference %s not being mapped\n",
						gffrec->getID(), refname.chars());
			return false; //discard, ref seq not in the given translation table
		}*/
	}
	if (gffloader.transcriptsOnly && gffrec->isDiscarded()) {
		//discard generic "locus" features with no other detailed subfeatures
		//GMessage("Warning: discarding %s GFF generic gene/locus container %s\n",gffrec->getID());
		return false;
	}
	if (minLen>0 && gffrec->covlen<minLen) {
		if (verbose)
			GMessage("Info: %s discarded due to minimum length threshold %d\n",
					gffrec->getID(), minLen);
    	return false;
	}
	if (rfltGSeq!=NULL) { //filter by gseqName
		if (strcmp(gffrec->getGSeqName(),rfltGSeq)!=0) {
			return false;
		}
	}
	if (rfltStrand>0 && gffrec->strand !=rfltStrand) {
		return false;
	}
	//check coordinates
	if (rfltStart!=0 || rfltEnd!=MAX_UINT) {
		if (rfltWithin) {
			if (gffrec->start<rfltStart || gffrec->end>rfltEnd) {
				return false; //not within query range
			}
		}
		else {
			if (gffrec->start>rfltEnd || gffrec->end<rfltStart) {
				return false;
			}
		}
	}
	if (multiExon && gffrec->exons.Count()<=1) {
		return false;
	}
	if (wCDSonly && gffrec->CDstart==0) {
		return false;
	}
	if (wNConly && gffrec->hasCDS()) return false;
	if (ensembl_convert && startsWith(gffrec->getID(), "ENS")) {
		//keep track of chr|gene_id data -- coordinate range
		char* geneid=gffrec->getGeneID();
		if (geneid!=NULL) {
			GeneInfo* ginfo=gene_ids.Find(geneid);
			if (ginfo==NULL) {//first time seeing this gene ID
				GeneInfo* geneinfo=new GeneInfo(gffrec, ensembl_convert);
				gene_ids.Add(geneid, geneinfo);
				if (gfnew!=NULL)
					gfnew->Add(geneinfo->gf); //do we really need this?
			}
			else ginfo->update(gffrec);
		}
	}
	return true;
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


void printGxfTab(FILE* f, GffObj& g) {
 //using attribute list in tableCols
	char* av=NULL;
	for(int i=0;i<tableCols.Count();i++) {
		if (i>0) fprintf(f,"\t");
		switch(tableCols[i].type) {
		case ctfGFF_Attr:
			av=g.getAttr(tableCols[i].name.chars());
			if (av!=NULL) fprintf(f,"%s",av);
			else fprintf(f, ".");
			break;
		case ctfGFF_chr:
			fprintf(f,"%s",g.getGSeqName());
			break;
		case ctfGFF_ID:
			fprintf(f,"%s",g.getID());
			break;
		case ctfGFF_Parent:
			if (g.parent!=NULL) fprintf(f,"%s",g.parent->getID());
			else fprintf(f, ".");
			break;
		case ctfGFF_feature:
			fprintf(f,"%s",g.getFeatureName());
			break;
		case ctfGFF_start:
			fprintf(f,"%d",g.start);
			break;
		case ctfGFF_end:
			fprintf(f,"%d",g.end);
			break;
		case ctfGFF_strand:
			fprintf(f,"%c",g.strand);
			break;
		case ctfGFF_numexons:
			fprintf(f,"%d",g.exons.Count());
			break;
		case ctfGFF_exons:
			if (g.exons.Count()>0) {
				for (int x=0;x<g.exons.Count();x++) {
					if (x>0) fprintf(f,",");
					fprintf(f,"%d-%d",g.exons[x]->start, g.exons[x]->end);
				}
			} else fprintf(f,".");
			break;
		case ctfGFF_cds:
			if (g.hasCDS()) {
				GVec<GffExon> cds;
				g.getCDSegs(cds);
				for (int x=0;x<cds.Count();x++) {
					if (x>0) fprintf(f,",");
				    fprintf(f,"%d-%d",cds[x].start, cds[x].end);
				}
			}
			else fprintf(f,".");
			break;
		case ctfGFF_covlen:
			fprintf(f, "%d", g.covlen);
			break;
		case ctfGFF_cdslen:
			if (g.hasCDS()) {
				GVec<GffExon> cds;
				g.getCDSegs(cds);
				int clen=0;
				for (int x=0;x<cds.Count();x++)
				    clen+=cds[x].end-cds[x].start+1;
				fprintf(f, "%d", clen);
			}
			else fprintf(f, "0");
			break;
		}
	}
	fprintf(f,"\n");
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
		 printGxfTab(f, *t.parent);
		 T_NO_PRINT(t.parent->udata);
	 }
    printGxfTab(f, *gfo);
}

int main(int argc, char* argv[]) {
 GArgs args(argc, argv,
   "version;debug;merge;adj-stop;bed;in-bed;tlf;in-tlf;cluster-only;nc;cov-info;help;"
    "sort-alpha;keep-genes;w-add=;keep-comments;keep-exon-attrs;force-exons;t-adopt;gene2exon;"
    "ignore-locus;no-pseudo;table=sort-by=hvOUNHPWCVJMKQYTDARSZFGLEBm:g:i:r:s:l:t:o:w:x:y:d:");
 args.printError(USAGE, true);
 if (args.getOpt('h') || args.getOpt("help")) {
    GMessage("%s",USAGE);
    exit(1);
 }

 debugMode=(args.getOpt("debug")!=NULL);
 decodeChars=(args.getOpt('D')!=NULL);
 gffloader.forceExons=(args.getOpt("force-exons")!=NULL);
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
 fmtGTF=(args.getOpt('T')!=NULL); //switch output format to GTF
 fmtBED=(args.getOpt("bed")!=NULL);
 fmtTLF=(args.getOpt("tlf")!=NULL);
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
 if (args.getOpt("version")) {
  GMessage(VERSION"\n");
  exit(0);
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
 ensembl_convert=(args.getOpt('L')!=NULL);
 if (ensembl_convert) {
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
 GFastaDb gfasta(args.getOpt('g'));
 //if (gfasta.fastaPath!=NULL)
 //    sortByLoc=true; //enforce sorting by chromosome/contig
 GStr s=args.getOpt('i');
 if (!s.is_empty()) maxintron=s.asInt();
 s=args.getOpt('l');
 if (!s.is_empty()) minLen=s.asInt();

 FILE* f_repl=NULL;
 s=args.getOpt('d');
 if (!s.is_empty()) {
   if (s=="-") f_repl=stdout;
   else {
       f_repl=fopen(s.chars(), "w");
       if (f_repl==NULL) GError("Error creating file %s\n", s.chars());
   }
 }

 rfltWithin=(args.getOpt('R')!=NULL);
 s=args.getOpt('r');
 if (!s.is_empty()) {
   s.trim();
   if (s[0]=='+' || s[0]=='-') {
     rfltStrand=s[0];
     s.cut(0,1);
     }
   int isep=s.index(':');
   if (isep>0) { //gseq name given
      if (rfltStrand==0 && (s[isep-1]=='+' || s[isep-1]=='-')) {
        isep--;
        rfltStrand=s[isep];
        s.cut(isep,1);
        }
      if (isep>0)
          rfltGSeq=Gstrdup((s.substr(0,isep)).chars());
      s.cut(0,isep+1);
      }
   GStr gsend;
   char slast=s[s.length()-1];
   if (rfltStrand==0 && (slast=='+' || slast=='-')) {
      s.chomp(slast);
      rfltStrand=slast;
      }
   if (s.index("..")>=0) gsend=s.split("..");
                    else gsend=s.split('-');
   if (!s.is_empty()) rfltStart=(uint)s.asInt();
   if (!gsend.is_empty()) {
      rfltEnd=(uint)gsend.asInt();
      if (rfltEnd==0) rfltEnd=MAX_UINT;
      }
   } //gseq/range filtering
 else {
   if (rfltWithin)
     GError("Error: option -R requires -r!\n");
   //if (rfltWholeTranscript)
   //  GError("Error: option -P requires -r!\n");
   }
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

 openfw(f_out, args, 'o');
 //if (f_out==NULL) f_out=stdout;
 if (gfasta.fastaPath==NULL && (validCDSonly || spliceCheck || args.getOpt('w')!=NULL || args.getOpt('x')!=NULL || args.getOpt('y')!=NULL))
  GError("Error: -g option is required for options -w, -x, -y, -V, -N, -M !\n");
 openfw(f_w, args, 'w');
 openfw(f_x, args, 'x');
 openfw(f_y, args, 'y');
 s=args.getOpt("w-add");
 if (!s.is_empty()) {
	 if (f_w==NULL) GError("Error: --w-add option requires -w option!\n");
	 wPadding=s.asInt();
 }

 if (f_out==NULL && f_w==NULL && f_x==NULL && f_y==NULL && !covInfo)
	 f_out=stdout;

 //if (f_y!=NULL || f_x!=NULL) wCDSonly=true;
 //useBadCDS=useBadCDS || (fgtfok==NULL && fgtfbad==NULL && f_y==NULL && f_x==NULL);

 int numfiles = args.startNonOpt();
 //GList<GffObj> gfkept(false,true); //unsorted, free items on delete
 int out_counter=0; //number of records printed
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
   gffloader.load(g_data, &validateGffRec, &processGffComment);
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
 GffPrintMode exonPrinting;
 if (fmtGTF)
	 exonPrinting = pgtfAny;
 else if (fmtBED)
	 exonPrinting=pgffBED;
 else if (fmtTLF)
	exonPrinting=pgffTLF;
 else { //printing regular GFF3
	exonPrinting = gffloader.forceExons ? pgffBoth : pgffAny;
 }
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
       int numvalid=0;
       int idxfirstvalid=-1;
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

         if (process_transcript(gfasta, t)) {
             numvalid++;
             if (idxfirstvalid<0) idxfirstvalid=i;
         }
       } //for each transcript

       int rnas_i=0;
       if (idxfirstvalid>=0) rnas_i=idxfirstvalid;
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
						   loc.print(f_out, idxfirstvalid, locname, loctrack);
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
					    	 loc.print(f_out, idxfirstvalid, locname, loctrack);
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
   int numvalid=0;
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
             if T_PRINTABLE(gfst.udata) { //never printed
               T_NO_PRINT(gfst.udata);
               if (fmtGFF3) {
                 if (firstGff3Print) { printGff3Header(f_out, args);firstGff3Print=false; }
                 if (firstGSeqHeader) { printGSeqHeader(f_out, gdata); firstGSeqHeader=false; }
                 gfst.printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
               }
               else printGxfTab(f_out, gfst);
             }
             ++gfs_i;
           }
        }
        GTData* tdata=(GTData*)(t.uptr);
        if (tdata->replaced_by!=NULL) continue;
        if (process_transcript(gfasta, t)) {
           numvalid++;
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
					 GTData* pdata=(GTData*)(t.parent->uptr);
					 if (pdata && pdata->geneinfo!=NULL)
						  pdata->geneinfo->finalize();
					 if (fmtTable)
						 printGxfTab(f_out, *(t.parent));
					 else
						 t.parent->printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
					 T_NO_PRINT(t.parent->udata);
				 }
				 if (fmtTable)
					 printGxfTab(f_out, t);
				 else
					 t.printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
             }
           }//GFF/GTF output requested
        } //valid transcript
     } //for each rna
     //print the rest of the isolated pseudo/gene/region features not printed yet
     if (f_out && (fmtGFF3 || fmtTable)) {
      while (gfs_i<gdata->gfs.Count()) {
         GffObj& gfst=*(gdata->gfs[gfs_i]);
         if T_PRINTABLE(gfst.udata) { //never printed
           T_NO_PRINT(gfst.udata);
           if (fmtGFF3) {
              if (firstGff3Print) { printGff3Header(f_out, args); firstGff3Print=false; }
              if (firstGSeqHeader) { printGSeqHeader(f_out, gdata); firstGSeqHeader=false; }
              gfst.printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
           } else
              printGxfTab(f_out, gfst);
         }
         ++gfs_i;
      }
     }
    } //for each genomic seq
   } //no clustering
 if (f_repl && f_repl!=stdout) fclose(f_repl);
 seqinfo.Clear();
 //if (faseq!=NULL) delete faseq;
 //if (gcdb!=NULL) delete gcdb;
 GFREE(rfltGSeq);
 FWCLOSE(f_out);
 FWCLOSE(f_w);
 FWCLOSE(f_x);
 FWCLOSE(f_y);
 }


