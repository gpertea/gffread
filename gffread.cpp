#include "GArgs.h"
#include "gff_utils.h"
#include <ctype.h>

#define VERSION "0.9.8"

#define USAGE "gffread v"VERSION". Usage:\n\
gffread <input_gff> [-g <genomic_seqs_fasta> | <dir>][-s <seq_info.fsize>] \n\
 [-o <outfile.gff>] [-t <tname>] [-r [[<strand>]<chr>:]<start>..<end> [-R]]\n\
 [-CTVNJMKQAFGUBHZWTOLE] [-w <exons.fa>] [-x <cds.fa>] [-y <tr_cds.fa>]\n\
 [-i <maxintron>] \n\
 Filters and/or converts GFF3/GTF2 records.\n\
 <input_gff> is a GFF file, use '-' if the GFF records will be given at stdin\n\
 \n\
 Options:\n\
  -g  full path to a multi-fasta file with the genomic sequences\n\
      for all input mappings, OR a directory with single-fasta files\n\
      (one per genomic sequence, with file names matching sequence names)\n\
  -s  <seq_info.fsize> is a tab-delimited file providing this info\n\
      for each of the mapped sequences:\n\
      <seq-name> <seq-length> <seq-description>\n\
      (useful for -A option with mRNA/EST/protein mappings)\n\
  -i  discard transcripts having an intron larger than <maxintron>\n\
  -r  only show transcripts overlapping coordinate range <start>..<end>\n\
      (on chromosome/contig <chr>, strand <strand> if provided)\n\
  -R  for -r option, discard all transcripts that are not fully \n\
      contained within the given range\n\
  -U  discard single-exon transcripts\n\
  -C  coding only: discard mRNAs that have no CDS feature\n\
  -F  full GFF attribute preservation (all attributes are shown)\n\
  -G  only parse additional exon attributes from the first exon\n\
      and move them to the mRNA level (useful for GTF input)\n\
  -A  use the description field from <seq_info.fsize> and add it\n\
      as the value for a 'descr' attribute to the GFF record\n\
  \n\
  -O  process also non-transcript GFF records (by default non-transcript\n\
      records are ignored)\n\
  -V  discard any mRNAs with CDS having in-frame stop codons\n\
  -H  for -V option, check and adjust the starting CDS phase\n\
      if the original phase leads to a translation with an \n\
      in-frame stop codon\n\
  -B  for -V option, single-exon transcripts are also checked on the\n\
      opposite strand\n\
  -N  discard multi-exon mRNAs that have any intron with a non-canonical\n\
      splice site consensus (i.e. not GT-AG, GC-AG or AT-AC)\n\
  -J  discard any mRNAs that either lack initial START codon\n\
      or the terminal STOP codon, or have an in-frame stop codon\n\
      (only print mRNAs with a fulll, valid CDS)\n\
  --no-pseudo: filter out records matching the 'pseudo' keyword\n\
 \n\
  -M/--merge : cluster the input transcripts into loci, collapsing matching\n\
       transcripts (those with the same exact introns and fully contained)\n\
  -d <dupinfo> : for -M option, write collapsing info to file <dupinfo>\n\
  --cluster-only: same as --merge but without collapsing matching transcripts\n\
  -K  for -M option: also collapse shorter, fully contained transcripts\n\
      with fewer introns than the container\n\
  -Q  for -M option, remove the containment restriction:\n\
      (multi-exon transcripts will be collapsed if just their introns match,\n\
      while single-exon transcripts can partially overlap (80%))\n\
 \n\
  --force-exons: make sure that the lowest level GFF features are printed as \n\
      \"exon\" features\n\
  -E  expose (warn about) duplicate transcript IDs and other potential \n\
      problems with the given GFF/GTF records\n\
  -D  decode url encoded characters within attributes\n\
  -Z  merge close exons into a single exon (for intron size<4)\n\
  -w  write a fasta file with spliced exons for each GFF transcript\n\
  -x  write a fasta file with spliced CDS for each GFF transcript\n\
  -W  for -w and -x options, also write for each fasta record the exon\n\
      coordinates projected onto the spliced sequence\n\
  -y  write a protein fasta file with the translation of CDS for each record\n\
  -L  Ensembl GTF to GFF3 conversion (implies -F; should be used with -m)\n\
  -m  <chr_replace> is a reference (genomic) sequence replacement table with\n\
      this format:\n\
      <original_ref_ID> <new_ref_ID>\n\
      GFF records on reference sequences that are not found among the\n\
      <original_ref_ID> entries in this file will be filtered out\n\
  -o  the \"filtered\" GFF records will be written to <outfile.gff>\n\
      (use -o- for printing to stdout)\n\
  -t  use <trackname> in the second column of each GFF output line\n\
  -T  -o option will output GTF format instead of GFF3\n\
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

bool validCDSonly=false; // translation with no in-frame STOP
bool bothStrands=false; //for single-exon mRNA validation, check the other strand too
bool altPhases=false; //if original phase fails translation validation,
                     //try the other 2 phases until one makes it
bool mRNAOnly=true;
bool NoPseudo=false;
bool forceExons=false;
bool spliceCheck=false; //only known splice-sites
bool decodeChars=false; //decode url-encoded chars in attrs (-D)
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

GFaSeqGet* fastaSeqGet(GFastaDb& gfasta, GffObj& gffrec) {
  if (gfasta.fastaPath==NULL) return NULL;
  return gfasta.fetch(gffrec.gseq_id);
}


int adjust_stopcodon(GffObj& gffrec, int adj, GList<GSeg>* seglst=NULL) {
 //adj>0 => extedn CDS,  adj<0 => shrink CDS
 //when CDS is expanded, exons have to be checked too and 
 // expanded accordingly if they had the same boundary
  int realadj=0;
  if (gffrec.strand=='-') {
       if ((int)gffrec.CDstart>adj) {

           gffrec.CDstart-=adj;
           realadj=adj;
           if (adj<0) { //restore
              if (gffrec.exons.First()->start==gffrec.CDstart+adj) {
                 gffrec.exons.First()->start-=adj;
                 gffrec.start=gffrec.exons.First()->start;
                 gffrec.covlen+=adj;
                 }
              }
           else if (gffrec.exons.First()->start>=gffrec.CDstart) {
                 gffrec.exons.First()->start-=adj;
                 gffrec.start=gffrec.exons.First()->start;
                 gffrec.covlen+=adj;
                 }
             }
          }
        else {
         realadj=adj;
         gffrec.CDend+=adj;
         if (adj<0) {//restore
           if (gffrec.exons.Last()->end==gffrec.CDend-adj) {
                        gffrec.exons.Last()->end+=adj;
                        gffrec.end=gffrec.exons.Last()->end;
                        gffrec.covlen+=adj;
                        }
          }
         else if (gffrec.exons.Last()->end<=gffrec.CDend) {
             gffrec.exons.Last()->end+=adj;
             gffrec.end=gffrec.exons.Last()->end;
             gffrec.covlen+=adj;
             }
         }
  if (seglst!=NULL) seglst->Last()->end+=adj;
  return realadj;
 }

bool process_transcript(GFastaDb& gfasta, GffObj& gffrec) {
 //returns true if the transcript passed the filter
 char* gname=gffrec.getGeneName();
 if (gname==NULL) gname=gffrec.getGeneID();
 GStr defline(gffrec.getID());
 if (f_out && !fmtGTF) {
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
   defline.appendfmt(" gene=%s", gname);
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
  GFaSeqGet* faseq=fastaSeqGet(gfasta, gffrec);
  if (spliceCheck && gffrec.exons.Count()>1) {
    //check introns for splice site consensi ( GT-AG, GC-AG or AT-AC )
    if (faseq==NULL) GError("Error: no genomic sequence available!\n");
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
      else if (acceptorSite=="AC") { //
         if (donorSite!="AT") { ssValid=false; break; }
         }
      else { ssValid=false; break; }
      }
    //GFREE(gseq);
    if (!ssValid) {
      if (verbose)
         GMessage("Invalid splice sites found for '%s'\n",gffrec.getID());
      return false; //don't print this one!
      }
    }

  bool trprint=true;
  int stopCodonAdjust=0;
  int mCDphase=0;
  bool hasStop=false;
  if (gffrec.CDphase=='1' || gffrec.CDphase=='2')
      mCDphase = gffrec.CDphase-'0';
  if (f_y!=NULL || f_x!=NULL || validCDSonly) {
    if (faseq==NULL) GError("Error: no genomic sequence provided!\n");
    //if (protmap && fullCDSonly) {
    //if (protmap && (fullCDSonly ||  (gffrec.qlen>0 && gffrec.qend==gffrec.qlen))) {
    
    if (validCDSonly) { //make sure the stop codon is always included 
      //adjust_stopcodon(gffrec,3);
      stopCodonAdjust=adjust_stopcodon(gffrec,3);
      }
    int strandNum=0;
    int phaseNum=0;
  CDS_CHECK:
    cdsnt=gffrec.getSpliced(faseq, true, &seqlen, NULL, NULL, &seglst);
    if (cdsnt==NULL) trprint=false;
    else { //has CDS
      if (validCDSonly) {
         cdsaa=translateDNA(cdsnt, aalen, seqlen);
         char* p=strchr(cdsaa,'.');
         hasStop=false;
         if (p!=NULL) {
              if (p-cdsaa>=aalen-2) { //stop found as the last codon
                      *p='0';//remove it
                      hasStop=true;
                      if (aalen-2==p-cdsaa) {
                        //previous to last codon is the stop codon
                        //so correct the CDS stop accordingly
                        adjust_stopcodon(gffrec,-3, &seglst);
                        stopCodonAdjust=0; //clear artificial stop adjustment
                        seqlen-=3;
                        cdsnt[seqlen]=0;
                        }
                      aalen=p-cdsaa;
                      }
                   else {//stop found before the last codon
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
           } //has in-frame STOP
         if (fullCDSonly) {
             if (!hasStop || cdsaa[0]!='M') trprint=false;
             }
         } // CDS check requested
      } //has CDS
    } //translation or codon check/output was requested
  if (!trprint) {
    GFREE(cdsnt);
    GFREE(cdsaa);
    return false;
    }
  if (stopCodonAdjust>0 && !hasStop) {
          //restore stop codon location
          adjust_stopcodon(gffrec, -stopCodonAdjust, &seglst);
          if (cdsnt!=NULL && seqlen>0) {
             seqlen-=stopCodonAdjust;
             cdsnt[seqlen]=0;
             }
          if (cdsaa!=NULL) aalen--;
          }

  if (f_y!=NULL) { //CDS translation fasta output requested
         //char* 
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
         printFasta(f_y, defline, cdsaa, aalen);
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
  fprintf(f, "# gffread v"VERSION"\n");
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


int main(int argc, char * const argv[]) {
 GArgs args(argc, argv, 
   "version;debug;merge;cluster-only;help;force-exons;no-pseudo;MINCOV=MINPID=hvOUNHWCVJMKQNSXTDAPRZFGLEm:g:i:r:s:t:a:b:o:w:x:y:d:");
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
 validCDSonly=(args.getOpt('V')!=NULL);
 altPhases=(args.getOpt('H')!=NULL);
 fmtGTF=(args.getOpt('T')!=NULL); //switch output format to GTF
 bothStrands=(args.getOpt('B')!=NULL);
 fullCDSonly=(args.getOpt('J')!=NULL);
 spliceCheck=(args.getOpt('N')!=NULL);
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
 if (f_y!=NULL || f_x!=NULL) wCDSonly=true;
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
   gffloader.fullAttributes=fullattr;
   gffloader.noExonAttrs=noExonAttr;
   gffloader.mergeCloseExons=mergeCloseExons;
   gffloader.showWarnings=(args.getOpt('E')!=NULL);
   gffloader.noPseudo=NoPseudo;
   gffloader.load(g_data, &validateGffRec, doCluster, doCollapseRedundant, 
                             matchAllIntrons, fuzzSpan, forceExons);
   if (doCluster) 
     collectLocusData(g_data);
   if (numfiles==0) break;
   }
   
 GStr loctrack("gffcl");
 if (tracklabel) loctrack=tracklabel;
 g_data.setSorted(&gseqCmpName);
 GffPrintMode exonPrinting;
 if (fmtGTF) {
	 exonPrinting = pgtfAny;
 } else {
	 exonPrinting = forceExons ? pgffBoth : pgffAny;
 }
 bool firstGff3Print=!fmtGTF;
 if (doCluster) {
   //grouped in loci
   for (int g=0;g<g_data.Count();g++) {
     GenomicSeqData* gdata=g_data[g];
     int gfs_i=0;
     for (int l=0;l<gdata->loci.Count();l++) {
       GffLocus& loc=*(gdata->loci[l]);
       //check all non-replaced transcripts in this locus:
       int numvalid=0;
       int idxfirstvalid=-1;
       for (int i=0;i<loc.rnas.Count();i++) {
         GffObj& t=*(loc.rnas[i]);
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
            continue;
            }
         if (process_transcript(gfasta, t)) {
             t.udata|=4; //tag it as valid
             numvalid++;
             if (idxfirstvalid<0) idxfirstvalid=i;
             }
         }
       if (f_out && numvalid>0) {
         GStr locname("RLOC_");
         locname.appendfmt("%08d",loc.locus_num);
         if (!fmtGTF) {
           if (firstGff3Print) { printGff3Header(f_out, args);firstGff3Print=false; }
           fprintf(f_out,"%s\t%s\tlocus\t%d\t%d\t.\t%c\t.\tID=%s;locus=%s",
                    loc.rnas[0]->getGSeqName(), loctrack.chars(), loc.start, loc.end, loc.strand,
                     locname.chars(), locname.chars());
           //const char* loc_gname=loc.getGeneName();
           if (loc.gene_names.Count()>0) { //print all gene names associated to this locus
              fprintf(f_out, ";genes=%s",loc.gene_names.First()->name.chars());
              for (int i=1;i<loc.gene_names.Count();i++) {
                fprintf(f_out, ",%s",loc.gene_names[i]->name.chars());
                }
              }
           if (loc.gene_ids.Count()>0) { //print all GeneIDs names associated to this locus
              fprintf(f_out, ";geneIDs=%s",loc.gene_ids.First()->name.chars());
              for (int i=1;i<loc.gene_ids.Count();i++) {
                fprintf(f_out, ",%s",loc.gene_ids[i]->name.chars());
                }
              }
           fprintf(f_out, ";transcripts=%s",loc.rnas[idxfirstvalid]->getID());
           for (int i=idxfirstvalid+1;i<loc.rnas.Count();i++) {
              fprintf(f_out, ",%s",loc.rnas[i]->getID());
              }
           fprintf(f_out, "\n");
           }
         //now print all valid, non-replaced transcripts in this locus:
         for (int i=0;i<loc.rnas.Count();i++) {
           GffObj& t=*(loc.rnas[i]);
           GTData* tdata=(GTData*)(t.uptr);
           if (tdata->replaced_by!=NULL || ((t.udata & 4)==0)) continue;
           t.addAttr("locus", locname.chars());
           out_counter++;
           if (fmtGTF) t.printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
               else {
                if (firstGff3Print) { printGff3Header(f_out, args);firstGff3Print=false; }
                //print the parent first, if any
                if (t.parent!=NULL && ((t.parent->udata & 4)==0)) {
                    GTData* pdata=(GTData*)(t.parent->uptr);
                    if (pdata && pdata->geneinfo!=NULL) 
                         pdata->geneinfo->finalize();
                    t.parent->addAttr("locus", locname.chars());
                    t.parent->printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
                    t.parent->udata|=4;
                    }
                t.printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
                }
            }
          } //have valid transcripts to print
       }//for each locus
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
    } //for each genomic sequence
   }
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
           t.udata|=4; //tag it as valid
           numvalid++;
           if (f_out) {
             if (tdata->geneinfo) tdata->geneinfo->finalize();
             out_counter++;
             if (fmtGTF) t.printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
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
   } //not clustered
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


