#include "GArgs.h"
#include "gff_utils.h"
#include <ctype.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#define VERSION "0.10.3"

#define USAGE "gffread v" VERSION ". Usage:\n\
gffread <input_gff> [-g <genomic_seqs_fasta> | <dir>][-s <seq_info.fsize>] \n\
 [-o <outfile.gff>] [-t <tname>] [-r [[<strand>]<chr>:]<start>..<end> [-R]]\n\
 [-CTVNJMKQAFPGUBHZWTOLE] [-w <exons.fa>] [-x <cds.fa>] [-y <tr_cds.fa>]\n\
 [-i <maxintron>] \n\
 Filters and/or converts GFF3/GTF2 records.\n\
 <input_gff> is a GFF file, use '-' if the GFF records will be given at stdin\n\
 \n\
 Options:\n\
 -g   full path to a multi-fasta file with the genomic sequences\n\
      for all input mappings, OR a directory with single-fasta files\n\
      (one per genomic sequence, with file names matching sequence names)\n\
 -s   <seq_info.fsize> is a tab-delimited file providing this info\n\
      for each of the mapped sequences:\n\
      <seq-name> <seq-length> <seq-description>\n\
      (useful for -A option with mRNA/EST/protein mappings)\n\
 -i   discard transcripts having an intron larger than <maxintron>\n\
 -r   only show transcripts overlapping coordinate range <start>..<end>\n\
      (on chromosome/contig <chr>, strand <strand> if provided)\n\
 -R   for -r option, discard all transcripts that are not fully \n\
      contained within the given range\n\
 -U   discard single-exon transcripts\n\
 -C   coding only: discard mRNAs that have no CDS features\n\
 --nc non-coding only: discard mRNAs that have CDS features\n\
 -F   full GFF attribute preservation (all attributes are shown)\n\
 -G   only parse additional exon attributes from the first exon\n\
      and move them to the mRNA level (useful for GTF input)\n\
 -A   use the description field from <seq_info.fsize> and add it\n\
      as the value for a 'descr' attribute to the GFF record\n\
 \n\
 -O   process also non-transcript GFF records (by default non-transcript\n\
      records are ignored)\n\
 -V   discard any mRNAs with CDS having in-frame stop codons (requires -g)\n\
 -H   for -V option, check and adjust the starting CDS phase\n\
      if the original phase leads to a translation with an \n\
      in-frame stop codon\n\
 -B   for -V option, single-exon transcripts are also checked on the\n\
      opposite strand (requires -g)\n\
 -P   add transcript level GFF attributes about the coding status of each\n\
      transcript, including partialness or in-frame stop codons (requires -g)\n\
 -N   discard multi-exon mRNAs that have any intron with a non-canonical\n\
      splice site consensus (i.e. not GT-AG, GC-AG or AT-AC)\n\
 -J   discard any mRNAs that either lack initial START codon\n\
      or the terminal STOP codon, or have an in-frame stop codon\n\
      (i.e. only print mRNAs with a complete CDS)\n\
 --no-pseudo: filter out records matching the 'pseudo' keyword\n\
 \n\
 -M/--merge : cluster the input transcripts into loci, collapsing matching\n\
       transcripts (those with the same exact introns and fully contained)\n\
 -d <dupinfo> : for -M option, write collapsing info to file <dupinfo>\n\
 --cluster-only: same as --merge but without collapsing matching transcripts\n\
 -K    for -M option: also collapse shorter, fully contained transcripts\n\
       with fewer introns than the container\n\
 -Q    for -M option, remove the containment restriction:\n\
       (multi-exon transcripts will be collapsed if just their introns match,\n\
       while single-exon transcripts can partially overlap (80%))\n\
 \n\
 --force-exons: make sure that the lowest level GFF features are printed as \n\
       \"exon\" features\n\
 --gene2exon: for single-line genes not parenting any transcripts, add an\n\
       exon feature spanning the entire gene (treat as transcript)\n\
 -E    expose (warn about) duplicate transcript IDs and other potential \n\
       problems with the given GFF/GTF records\n\
 -D    decode url encoded characters within attributes\n\
 -Z    merge close exons into a single exon (for intron size<4)\n\
 -w    write a fasta file with spliced exons for each GFF transcript\n\
 -x    write a fasta file with spliced CDS for each GFF transcript\n\
 -W    for -w and -x options, write in the FASTA defline the exon\n\
       coordinates projected onto the spliced sequence;\n\
       for -y option, write transcript attributes in the FASTA defline\n\
 -y    write a protein fasta file with the translation of CDS for each record\n\
 -S    for -y option, use '*' instead of '.' as stop codon translation\n\
 -L    Ensembl GTF to GFF3 conversion (implies -F; should be used with -m)\n\
 -m    <chr_replace> is a reference (genomic) sequence replacement table with\n\
       this format:\n\
       <original_ref_ID> <new_ref_ID>\n\
       GFF records on reference sequences that are not found among the\n\
       <original_ref_ID> entries in this file will be filtered out\n\
 -o    the \"filtered\" GFF records will be written to <outfile.gff>\n\
        (use -o- to enable printing to stdout)\n\
 -t    use <trackname> in the 2nd column of each GFF/GTF output line\n\
 -T    output GTF instead of GFF3 (for -o) \n\
 --bed output BED format instead of GFF3 (for -o)\n\
 --tab for -o option, output transcripts in one-line tab delimited format:\n\
      tID  chr  +/-  t_start  t_end  exonCount  exons  CDScoords  GFF_attrs ...\n\
      (<exons> are shown as a comma-delimited list of start-end coordinates;\n\
      <CDScoords> is '.' if no CDS is present, or CDS_start:CDS_end otherwise)\n\
	  The other GFF attributes, if any, are appended as name=value columns\n\
"


class SeqInfo { //populated from the -s option of gffread
 public:
  int len;
  char* descr;
  SeqInfo( int l, char* s) {
   len=l;
   if (s==NULL) {
     descr=NULL;
     }   else {
     descr=Gstrdup(s);
     }
   }
  ~SeqInfo() {
   GFREE(descr);
   }
};

struct CStopAdjData {
	int gseqlen;
	GffObj* t;
	int CDS_shift;
	int exon_shift;
	CStopAdjData(uint glen=0, GffObj* gfo=NULL):gseqlen(glen), t(gfo),
		CDS_shift(0), exon_shift(0) {
	}
	int apply(int cdshift, bool reset=false) {
		if (cdshift==0) {
			if (reset) {
				CDS_shift=0;
				exon_shift=0;
			}
			return 0;
		}
		if (t->strand=='-') {
		  if (cdshift>0 && (int)t->CDstart>cdshift) {
			  CDS_shift+=cdshift;
			  t->CDstart-=cdshift;
			  if (t->exons.First()->start>t->CDstart) {
				  int eshift=t->exons.First()->start-t->CDstart;
				  exon_shift+=eshift;
				  t->exons.First()->start-=eshift;
				  t->start-=eshift;
				  t->covlen+=eshift;
			  }
		  }
		  else if (cdshift<0) { //shrinking CDS
			   //only used in order to undo a previous expansion
               CDS_shift+=cdshift;
               t->CDstart-=cdshift;
               if (exon_shift>0) { //undo previous expansion
            	   int eshift=-exon_shift;
            	   t->exons.First()->start-=eshift;
            	   t->start-=eshift;
            	   t->covlen+=eshift;
            	   exon_shift=0;
               }
		  }
		}
		else { //forward strand
		  if (cdshift>0 && (int)t->CDend+cdshift<=gseqlen) {
			  CDS_shift+=cdshift;
			  t->CDend+=cdshift;
			  if (t->exons.Last()->end<t->CDend) {
				  int eshift=t->CDend-t->exons.Last()->end;
				  exon_shift+=eshift;
				  t->exons.Last()->end+=eshift;
				  t->end+=eshift;
				  t->covlen+=eshift;
			  }

		  }
		  else if (cdshift<0) { //shrinking CDS
              CDS_shift+=cdshift;
              t->CDend+=cdshift;
              if (exon_shift>0) { //undo previous expansion
           	   int eshift=-exon_shift;
           	   t->exons.Last()->end+=eshift;
           	   t->end+=eshift;
           	   t->covlen+=eshift;
           	   exon_shift=0;
              }
		  }
		}

		if (reset) {
			CDS_shift=0;
			exon_shift=0;
		}
		return CDS_shift;
	}

	int restore() {
      if (CDS_shift==0) return 0;
      int r=CDS_shift;
	  if (t->strand=='-') {
		  t->CDstart+=CDS_shift;
		  if (exon_shift!=0) {
			t->exons.First()->start+=exon_shift;
			t->covlen-=exon_shift;
			t->start+=exon_shift;
		  }
	  } else {
		  t->CDend-=CDS_shift;
		  if (exon_shift!=0) {
			t->exons.Last()->end-=exon_shift;
			t->covlen-=exon_shift;
			t->end+=exon_shift;
		  }
	   }
	  CDS_shift=0;
	  exon_shift=0;
      return r;
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

FILE* ffasta=NULL;
FILE* f_in=NULL;
FILE* f_out=NULL;
FILE* f_w=NULL; //fasta with spliced exons (transcripts)
FILE* f_x=NULL; //fasta with spliced CDS
FILE* f_y=NULL; //fasta with translated CDS
bool wCDSonly=false;
bool wNConly=false;

bool validCDSonly=false; // translation with no in-frame STOP
bool bothStrands=false; //for single-exon mRNA validation, check the other strand too
bool altPhases=false; //if original phase fails translation validation,
                     //try the other 2 phases until one makes it
bool addCDSattrs=false;
bool covInfo=false; // --cov-info option -- report genome coverage per strand
bool mRNAOnly=true;
bool NoPseudo=false;
bool forceExons=false;
bool spliceCheck=false; //only known splice-sites
bool decodeChars=false; //decode url-encoded chars in attrs (-D)
bool StarStop=false; //use * instead of . for stop codon translation

bool fullCDSonly=false; // starts with START, ends with STOP codon
bool fullattr=false;
//bool sortByLoc=false; // if the GFF output should be sorted by location
bool ensembl_convert=false; //-L, assist in converting Ensembl GTF to GFF3


//GStr gseqpath;
//GStr gcdbfa;
//bool multiGSeq=false; //if a directory or a .cidx file was given to -g option
//GFaSeqGet* faseq=NULL;
//GCdbYank* gcdb=NULL;
//int gseq_id=-1; //current genome sequence ID -- the current GffObj::gseq_id
bool fmtGTF=false;
bool fmtBED=false;
bool fmtTab=false;
bool addDescr=false;
//bool protmap=false;
bool multiExon=false;
bool writeExonSegs=false;
char* tracklabel=NULL;
int maxintron=999000000;
bool mergeCloseExons=false;
//range filter:
char* rfltGSeq=NULL;
char rfltStrand=0;
uint rfltStart=0;
uint rfltEnd=MAX_UINT;
bool rfltWithin=false; //check for full containment within given range
bool noExonAttr=false;

bool doCluster=false;
bool doCollapseRedundant=false;

GList<GenomicSeqData> g_data(true,true,true); //list of GFF records by genomic seq

//hash with sequence info
GHash<SeqInfo> seqinfo;
GHash<int> isoCounter; //counts the valid isoforms
GHash<RefTran> reftbl;
GHash<GeneInfo> gene_ids;
  //min-max gene span associated to chr|gene_id (mostly for Ensembl conversion)

bool debugMode=false;
bool verbose=false;

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
 //returns true if the transcript passed the filter
 char* gname=gffrec.getGeneName();
 if (gname==NULL) gname=gffrec.getGeneID();
 GStr defline(gffrec.getID());
 if (f_out && !fmtGTF && !fmtBED && !fmtTab) {
     const char* tname=NULL;
     if ((tname=gffrec.getAttr("transcript_name"))!=NULL) {
        gffrec.addAttr("Name", tname);
        gffrec.removeAttr("transcript_name");
        }
     }
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
  GList<GSeg> seglst(false,true);
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
  //int stopCodonAdjust=0;
  int mCDphase=0;
  bool fullCDS=false;
  bool endStop=false;
  if (addCDSattrs && gffrec.hasCDS()) gffrec.addAttr("hasCDS", "true");
  if (gffrec.CDphase=='1' || gffrec.CDphase=='2')
      mCDphase = gffrec.CDphase-'0';
  //CDS partialness only added when -y -x -V options are given
  if (f_y!=NULL || f_x!=NULL || validCDSonly || addCDSattrs) {
    int strandNum=0;
    int phaseNum=0;
  CDS_CHECK:
    cdsnt=gffrec.getSpliced(faseq, true, &seqlen, NULL, NULL, &seglst);
    //if (cdsnt==NULL) trprint=false;
    if (cdsnt!=NULL) { //has CDS
      //if (validCDSonly) {
         cdsaa=translateDNA(cdsnt, aalen, seqlen);
         char* p=strchr(cdsaa,'.');
         endStop=false;
         if (p!=NULL) {
              //if (p-cdsaa>=aalen-2) { //stop found as the last OR prev-to-last codon
        	 if (p-cdsaa==aalen-1) { //stop found as the last codon
                  *p='\0';//remove it
                  endStop=true;
                  aalen--;
                  /*
                  if (p-cdsaa==aalen-2) {
                    //previous to last codon is the stop codon
                    //so correct the CDS stop accordingly
                    //adjstop->apply(-3, true);
                    if (seglst.Count()>0) seglst.Last()->end-=3;
                    //stopCodonAdjust=0; //clear artificial stop adjustment
                    seqlen-=3;
                    cdsnt[seqlen]=0;
                  }
                  else {
                	 //last codon is a stop codon
                	 adjstop->apply(0, true);
                  }
                  aalen=p-cdsaa;
                  */
              }
              else {//stop found before the last codon - not valid
                  trprint=false;
              }
         }//stop codon found
         if (trprint==false) { //failed CDS validity check
           //in-frame stop codon found
           if (altPhases && phaseNum<3) {
              phaseNum++;
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
           if (verbose) GMessage("In-frame STOP found for '%s'\n",gffrec.getID());
           gffrec.addAttr("InFrameStop", "true");
         } //has in-frame STOP
         bool hasStart=(cdsaa[0]=='M'); //could be a
         fullCDS=(endStop && hasStart);
         if (!fullCDS) {
        	 const char* partialness=NULL;
        	 if (hasStart) partialness="3";
        	 else {
        		partialness = endStop ? "5" : "5_3";
        	 }
        	 gffrec.addAttr("partialness", partialness);
         }
         if (trprint && fullCDSonly && !fullCDS)
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
  if (cdsnt!=NULL) {
	  if (f_y!=NULL) { //CDS translation fasta output requested
			 if (cdsaa==NULL) { //translate now if not done before
			   cdsaa=translateDNA(cdsnt, aalen, seqlen);
			   }
			 if (fullattr && gffrec.attrs!=NULL) {
				 //append all attributes found for each transcripts
				for (int i=0;i<gffrec.attrs->Count();i++) {
				  defline.append(" ");
				  defline.append(gffrec.getAttrName(i));
				  defline.append("=");
				  defline.append(gffrec.getAttrValue(i));
				}
			 }
			 if (validCDSonly && !fullattr && !fullCDS) {
				defline.append(" [partial]");
			 }
			 if (aalen>0) {
			   if (cdsaa[aalen-1]=='.') --aalen; //avoid printing the final stop codon
			   printFasta(f_y, defline, cdsaa, aalen, StarStop);
			 }
	  }
	  if (f_x!=NULL) { //CDS only
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
					  defline+=(int)seglst[i]->start;
					  defline.append("-");
					  defline+=(int)seglst[i]->end;
					  }
				  }
			 if (fullattr && gffrec.attrs!=NULL) {
				 //append all attributes found for each transcript
				  for (int i=0;i<gffrec.attrs->Count();i++) {
					defline.append(" ");
					defline.append(gffrec.getAttrName(i));
					defline.append("=");
					defline.append(gffrec.getAttrValue(i));
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
      char* exont=gffrec.getSpliced(faseq, false, &seqlen, &cds_start, &cds_end, &seglst);
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
        defline.append(" segs:");
        for (int i=0;i<seglst.Count();i++) {
            if (i>0) defline.append(",");
            defline+=(int)seglst[i]->start;
            defline.append("-");
            defline+=(int)seglst[i]->end;
            }
      }
      if (fullattr && gffrec.attrs!=NULL) {
       //append all attributes found for each transcripts
        for (int i=0;i<gffrec.attrs->Count();i++) {
          defline.append(" ");
          defline.append(gffrec.getAttrName(i));
          defline.append("=");
          defline.append(gffrec.getAttrValue(i));
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
  fprintf(f, "# ");
  args.printCmdLine(f);
  fprintf(f, "# gffread v" VERSION "\n");
  fprintf(f, "##gff-version 3\n");
  //for (int i=0;i<gseqdata.Count();i++) {
  //
  //}
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
		else return false; //discard, ref seq not in the given translation table
	}
	if (mRNAOnly && gffrec->isDiscarded()) {
		//discard generic "locus" features with no other detailed subfeatures
		//GMessage("Warning: discarding %s GFF generic gene/locus container %s\n",gffrec->getID());
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
					gfnew->Add(geneinfo->gf); //FIXME: do I really need this?
			}
			else ginfo->update(gffrec);
		}
	}
	return true;
}

void printGffObj(FILE* f, GffObj* gfo, GStr& locname, GffPrintMode exonPrinting, int& out_counter) {
    GffObj& t=*gfo;
    GTData* tdata=(GTData*)(t.uptr);
    if (tdata->replaced_by!=NULL || ((t.udata & 4)!=0)) return;
    if (t.exons.Count()==0 && t.children.Count()==0 && forceExons)
      t.addExon(t.start,t.end);
    t.udata|=4;
    t.addAttr("locus", locname.chars());
    out_counter++;
    if (fmtGTF) t.printGxf(f, exonPrinting, tracklabel, NULL, decodeChars);
        else {
         //print the parent first, if any, and if not printed already
         if (t.parent!=NULL && ((t.parent->udata & 4)==0)) {
             GTData* pdata=(GTData*)(t.parent->uptr);
             if (pdata && pdata->geneinfo!=NULL)
                  pdata->geneinfo->finalize();
             t.parent->addAttr("locus", locname.chars());
             t.parent->printGxf(f, exonPrinting, tracklabel, NULL, decodeChars);
             t.parent->udata|=4;
         }
         t.printGxf(f, exonPrinting, tracklabel, NULL, decodeChars);
    }

}


int main(int argc, char* argv[]) {
 GArgs args(argc, argv,
   "version;debug;merge;bed;tab;cluster-only;nc;cov-info;help;force-exons;gene2exon;no-pseudo;MINCOV=MINPID=hvOUNHPWCVJMKQTDARSZFGLEBm:g:i:r:s:t:o:w:x:y:d:");
 args.printError(USAGE, true);
 if (args.getOpt('h') || args.getOpt("help")) {
    GMessage("%s",USAGE);
    exit(1);
    }
 debugMode=(args.getOpt("debug")!=NULL);
 decodeChars=(args.getOpt('D')!=NULL);
 forceExons=(args.getOpt("force-exons")!=NULL);
 NoPseudo=(args.getOpt("no-pseudo")!=NULL);
 mRNAOnly=(args.getOpt('O')==NULL);
 //sortByLoc=(args.getOpt('S')!=NULL);
 addDescr=(args.getOpt('A')!=NULL);
 verbose=(args.getOpt('v')!=NULL);
 wCDSonly=(args.getOpt('C')!=NULL);
 wNConly=(args.getOpt("nc")!=NULL);
 addCDSattrs=(args.getOpt('P')!=NULL);
 validCDSonly=(args.getOpt('V')!=NULL);
 altPhases=(args.getOpt('H')!=NULL);
 fmtGTF=(args.getOpt('T')!=NULL); //switch output format to GTF
 fmtBED=(args.getOpt("bed")!=NULL);
 fmtTab=(args.getOpt("tab")!=NULL);
 bothStrands=(args.getOpt('B')!=NULL);
 fullCDSonly=(args.getOpt('J')!=NULL);
 spliceCheck=(args.getOpt('N')!=NULL);
 StarStop=(args.getOpt('S')!=NULL);
 bool gene2exon=(args.getOpt("gene2exon")!=NULL);
 bool matchAllIntrons=(args.getOpt('K')==NULL);
 bool fuzzSpan=(args.getOpt('Q')!=NULL);
 if (args.getOpt('M') || args.getOpt("merge")) {
    doCluster=true;
    doCollapseRedundant=true;
    }
   else {
    if (!matchAllIntrons || fuzzSpan) {
      GMessage("%s",USAGE);
      GMessage("Error: -K or -Q options require -M/--merge option!\n");
      exit(1);
      }
    }
 if (args.getOpt("cluster-only")) {
    doCluster=true;
    doCollapseRedundant=false;
    if (!matchAllIntrons || fuzzSpan) {
      GMessage("%s",USAGE);
      GMessage("Error: -K or -Q options have no effect with --cluster-only.\n");
      exit(1);
      }
 }
 covInfo=(args.getOpt("cov-info"));
 if (covInfo) doCluster=true;
 if (fullCDSonly) validCDSonly=true;
 if (verbose) {
     fprintf(stderr, "Command line was:\n");
     args.printCmdLine(stderr);
     }
 if (args.getOpt("version")) {
  GMessage(VERSION"\n");
  exit(0);
 }
 fullattr=(args.getOpt('F')!=NULL);
 if (args.getOpt('G')==NULL)
    noExonAttr=!fullattr;
   else {
     noExonAttr=true;
     fullattr=true;
     }
 if (NoPseudo && !fullattr) {
	 noExonAttr=true;
	 fullattr=true;
 }
 ensembl_convert=(args.getOpt('L')!=NULL);
 if (ensembl_convert) {
    fullattr=true;
    noExonAttr=false;
    //sortByLoc=true;
    }

 mergeCloseExons=(args.getOpt('Z')!=NULL);
 multiExon=(args.getOpt('U')!=NULL);
 writeExonSegs=(args.getOpt('W')!=NULL);
 tracklabel=args.getOpt('t');
 GFastaDb gfasta(args.getOpt('g'));
 //if (gfasta.fastaPath!=NULL)
 //    sortByLoc=true; //enforce sorting by chromosome/contig
 GStr s=args.getOpt('i');
 if (!s.is_empty()) maxintron=s.asInt();

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
               else
                 if ((f_in=fopen(infile, "r"))==NULL)
                    GError("Error: cannot open input file %s!\n",infile.chars());
                 else fclose(f_in);
          }
        else infile="-";
   GffLoader gffloader(infile.chars());
   gffloader.transcriptsOnly=mRNAOnly;
   gffloader.gene2exon=gene2exon;
   gffloader.fullAttributes=fullattr;
   gffloader.noExonAttrs=noExonAttr;
   gffloader.mergeCloseExons=mergeCloseExons;
   gffloader.showWarnings=(args.getOpt('E')!=NULL);
   gffloader.noPseudo=NoPseudo;
   gffloader.load(g_data, &validateGffRec, doCluster, doCollapseRedundant,
                             matchAllIntrons, fuzzSpan, forceExons);
   if (doCluster)
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
	 fprintf(stdout, "\t%" PRIu64 " on + strand\n", f_bases);
	 fprintf(stdout, "\t%" PRIu64 " on - strand\n", r_bases);
	 fprintf(stdout, "\t%" PRIu64 " on . strand\n", u_bases);
 }
 GStr loctrack("gffcl");
 if (tracklabel) loctrack=tracklabel;
 g_data.setSorted(&gseqCmpName);
 GffPrintMode exonPrinting;
 if (fmtGTF || fmtBED || fmtTab) {
	 exonPrinting = pgtfAny;
 } else {
	 exonPrinting = forceExons ? pgffBoth : pgffAny;
 }
 bool firstGff3Print=! (fmtGTF || fmtBED || fmtTab);
 if (doCluster) {
   //grouped in loci
   for (int g=0;g<g_data.Count();g++) {
     GenomicSeqData* gdata=g_data[g];
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
            if (f_repl && (t.udata & 8)==0) {
               //t.udata|=8;
               fprintf(f_repl, "%s", t.getID());
               GTData* rby=tdata;
               while (rby->replaced_by!=NULL) {
                  fprintf(f_repl," => %s", rby->replaced_by->getID());
                  rby->rna->udata|=8;
                  rby=(GTData*)(rby->replaced_by->uptr);
                  }
               fprintf(f_repl, "\n");
               }
            t.udata|=4; //not going to print this
            continue;
            }
         if (process_transcript(gfasta, t)) {
             //t.udata|=4; //tag it as valid
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
				   if (!(fmtBED || fmtGTF || fmtTab)) { //sorry, BED only shows transcripts
					   if (firstGff3Print) { printGff3Header(f_out, args);firstGff3Print=false; }
					   if (firstLocusPrint) { loc.print(f_out, idxfirstvalid, locname, loctrack, fmtGTF);firstLocusPrint=false; }
					   printGffObj(f_out, loc.gfs[gfs_i], locname, exonPrinting, out_counter);
				   }
				   ++gfs_i;
				   continue;
			   }
			   if (rnas_i<loc.rnas.Count()) {
				       if (fmtBED)
				    	   loc.rnas[rnas_i]->printBED(f_out);
				       else if (fmtTab)
				    	   printTabFormat(f_out, loc.rnas[rnas_i]);
				       else if (fmtGTF)
				    	   loc.rnas[rnas_i]->printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
				       else {
					     if (firstGff3Print) { printGff3Header(f_out, args);firstGff3Print=false; }
					     if (firstLocusPrint) { loc.print(f_out, idxfirstvalid, locname, loctrack, fmtGTF);firstLocusPrint=false; }
					     printGffObj(f_out, loc.rnas[rnas_i], locname, exonPrinting, out_counter);
				       }
					   ++rnas_i;
			   }
		   }
       }

     }//for each locus
    } //for each genomic sequence
   } //if Clustering enabled
  else {
   //not grouped into loci, print the rnas with their parents, if any
   int numvalid=0;
   for (int g=0;g<g_data.Count();g++) {
     GenomicSeqData* gdata=g_data[g];
     int gfs_i=0;
     for (int m=0;m<gdata->rnas.Count();m++) {
        GffObj& t=*(gdata->rnas[m]);
        if (f_out) {
         while (gfs_i<gdata->gfs.Count() && gdata->gfs[gfs_i]->start<=t.start) {
            GffObj& gfst=*(gdata->gfs[gfs_i]);
            if ((gfst.udata&4)==0) { //never printed
              gfst.udata|=4;
              if (firstGff3Print) { printGff3Header(f_out, args);firstGff3Print=false; }
              if (gfst.exons.Count()==0 && gfst.children.Count()==0 && forceExons)
               gfst.addExon(gfst.start,gfst.end);
              gfst.printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
            }
            ++gfs_i;
         }
        }
        GTData* tdata=(GTData*)(t.uptr);
        if (tdata->replaced_by!=NULL) continue;
        if (process_transcript(gfasta, t)) {
           numvalid++;
           if (f_out && (t.udata & 4) ==0 ) {
             if (tdata->geneinfo) tdata->geneinfo->finalize();
             out_counter++;
             t.udata |=4;
             if (fmtBED) t.printBED(f_out);
             else if (fmtTab) printTabFormat(f_out, &t);
             else if (fmtGTF) t.printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
             else {
                if (firstGff3Print) { printGff3Header(f_out, args);firstGff3Print=false; }
                //print the parent first, if any
                if (t.parent!=NULL && ((t.parent->udata & 4)==0)) {
                    GTData* pdata=(GTData*)(t.parent->uptr);
                    if (pdata && pdata->geneinfo!=NULL)
                         pdata->geneinfo->finalize();
                    t.parent->printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
                    t.parent->udata|=4;
                    }
                t.printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
                }
             }//GFF/GTF output requested
           } //valid transcript
        } //for each rna
     //print the rest of the isolated pseudo/gene/region features not printed yet
     if (f_out) {
      while (gfs_i<gdata->gfs.Count()) {
         GffObj& gfst=*(gdata->gfs[gfs_i]);
         if ((gfst.udata&4)==0) { //never printed
           gfst.udata|=4;
           if (firstGff3Print) { printGff3Header(f_out, args);firstGff3Print=false; }
           if (gfst.exons.Count()==0 && gfst.children.Count()==0 && forceExons)
            gfst.addExon(gfst.start,gfst.end);
           gfst.printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
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


