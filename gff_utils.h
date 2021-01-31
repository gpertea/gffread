#ifndef GFF_UTILS_H
#define GFF_UTILS_H
#include "gff.h"
#include "GStr.h"
#include "GVec.hh"
#include "GFaSeqGet.h"

extern bool verbose;
extern bool debugMode;
extern bool ensembl_convert;

extern FILE* ffasta;
extern FILE* f_in;
extern FILE* f_out;
extern FILE* f_w; //writing fasta with spliced exons (transcripts)
extern int wPadding; //padding for -w option
extern FILE* f_x; //writing fasta with spliced CDS
extern FILE* f_y; //wrting fasta with translated CDS

extern FILE* f_j; //wrting junctions (introns)

extern bool TFilters;

extern bool wfaNoCDS;

extern int maxintron;

extern bool wCDSonly;
extern bool wNConly;

enum ID_Flt_Type {
	idFlt_None=0,
	idFlt_Only,
	idFlt_Exclude
};

extern ID_Flt_Type IDflt;

extern int minLen; //minimum transcript length
extern bool validCDSonly; // translation with no in-frame STOP
extern bool bothStrands; //for single-exon mRNA validation, check the other strand too
extern bool altPhases; //if original phase fails translation validation,
                     //try the other 2 phases until one makes it
extern bool addCDSattrs;
extern bool add_hasCDS;

extern bool adjustStop; //automatic adjust the CDS stop coordinate
extern bool covInfo; // --cov-info : only report genome coverage
extern GStr tableFormat; //list of "attributes" to print in tab delimited format
extern bool spliceCheck; //only known splice-sites
extern bool decodeChars; //decode url-encoded chars in attrs (-D)
extern bool StarStop; //use * instead of . for stop codon translation
extern bool fullCDSonly; // starts with START, ends with STOP codon

extern bool multiExon;
extern bool writeExonSegs;
extern char* tracklabel;
extern bool rfltWithin; //check for full containment within given range
extern bool addDescr;


extern bool fmtGFF3; //output: GFF3
//other formats only make sense in transcriptOnly mode
extern bool fmtGTF;
extern bool fmtBED;
extern bool fmtTLF;
extern bool fmtTable;


extern GffPrintMode exonPrinting;

//typedef bool GFValidateFunc(GffObj* gf, GList<GffObj>* gfadd);
typedef bool GFValidateFunc(GffObj* gf);

//test if a transcript should be printed (and not printed yet)
#define T_PRINTABLE(d) (((d) & 0x100)==0)

//set a transcript to not be printed
#define T_NO_PRINT(d) d |= 0x100

//test if a duplicate transcript should be shown in the duplicate info file
#define T_DUPSHOWABLE(d) (((d) & 0x200)==0)

//set a duplicate transcript to not be shown in the duplicate info file
#define T_NO_DUPSHOW(d) d |= 0x200

//check original/old strand:
#define T_OSTRAND(d) (d & 0xFF)
//keep/set original/old strand
#define T_SET_OSTRAND(d, s) d |= s

extern GRangeParser* fltRange;

extern GRangeParser* fltJunction;

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

extern GFastaDb gfasta;
extern GHash<SeqInfo*> seqinfo;
extern GHash<int> isoCounter; //counts the valid isoforms
extern GHash<RefTran*> reftbl;

char* getSeqDescr(char* seqid);
char* getSeqName(char* seqid);
int adjust_stopcodon(GffObj& gffrec, int adj, GList<GSeg>* seglst=NULL);
void printTableData(FILE* f, GffObj& g, bool inFasta=false);

enum ETableFieldType {
  ctfGFF_Attr=0, // attribute name as is
  ctfGFF_ID, //ID or @id or transcript_id
  ctfGFF_geneID, //geneID or @gene_id or @geneid
  ctfGFF_geneName, //geneName or @gene_name or @genename
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


extern GVec<CTableField> tableCols; //table output format fields
extern GStrSet<> attrList;
extern GStrSet<> fltIDs;

class GffLocus;
class GenomicSeqData;
class GeneInfo;

struct CIntronData:public GSeg {
	char strand;//'.' < '-' < '+' (reverse ASCII order)
	GVec<GStr> ts; //list of transcript IDs sharing this intron
	CIntronData(uint istart, uint iend, char tstrand, const char* t_id=NULL):GSeg(istart, iend),
			strand(tstrand) {
		if (t_id!=NULL) {
			GStr tid(t_id);
		    ts.Add(tid);
		}
	}
	void add(const char* t_id) {
		 GStr tid(t_id);
		 ts.Add(tid);
	}
	bool operator==(CIntronData& d){
	  return (start==d.start && end==d.end && strand==d.strand);
	}
	bool operator<(CIntronData& d){
	 if (start==d.start) {
		 if (end==d.end) return strand>d.strand;
			else return (end<d.end);
	  }
	  else return (start<d.start);
	}

};

struct CIntronList {
	int gseq_id;
	uint last_t_start; //just to check if input is sorted properly!
	GList<CIntronData> jlst;
	CIntronList():gseq_id(-1),last_t_start(0), jlst(true, true) {}
	void add(GffObj& t) { //add all introns of t to jlst
		if (t.exons.Count()<2) return; //nothing to do
		if (gseq_id>=0 && gseq_id!=t.gseq_id)
			GError("Error: CIntronList::add(%s) on different ref seq!\n", t.getID());
		gseq_id=t.gseq_id;
		for (int i=1;i<t.exons.Count();++i) {
		  CIntronData* nintr = new CIntronData(t.exons[i-1]->end+1,
				  t.exons[i]->start-1, t.strand, t.getID());
		  int fidx=-1;
		  CIntronData* xintr=jlst.AddIfNew(nintr, true, &fidx);
		  if (xintr!=nintr) {
			  //nintr already exists,it was deallocated
			  xintr->add(t.getID());
		  }
		  last_t_start=t.start;
		} //for each intron
	}
	void clear() {
		gseq_id=-1;
		last_t_start=0;
		jlst.Clear();
	}
	void print(FILE* f) {
		//simple tab delimited format: chr, start, end, strand, transcriptIDs comma delimited
		const char* gseqname=GffObj::names->gseqs.getName(gseq_id);
		for (int i=0;i<jlst.Count();++i) {
			CIntronData& idata=*(jlst[i]);
			fprintf(f,"%s\t%d\t%d\t%c\t",gseqname, idata.start, idata.end, idata.strand);
			if (idata.ts.Count()>1)
				idata.ts.Sort();
			for (int t=0;t<idata.ts.Count();t++) {
				if (t) fprintf(f, ",%s", idata.ts[t].chars());
				  else fprintf(f,  "%s", idata.ts[t].chars());
			}
			fprintf(f, "\n");
		}
	}
};

class GTData { // transcript associated data
 public:
   GffObj* rna;
   GenomicSeqData* gdata;
   GffLocus* locus;
   GffObj* replaced_by;
   GeneInfo* geneinfo;
   GTData(GffObj* t=NULL, GenomicSeqData* gd=NULL);
   bool operator<(GTData& b) { return (rna < b.rna); }
   bool operator==(GTData& b) { return (rna==b.rna); }
};

class GeneInfo {
 public:
   int flag;
   GffObj* gf;
   GList<GStr> gene_names;
   GList<GStr> transcripts; //list of transcript IDs
   GeneInfo():gene_names(true, true, true), transcripts(true,true,true) {
     gf=NULL;
     flag=0;
   }

   GeneInfo(GffObj* gfrec, GenomicSeqData* gdata, bool ensembl_convert=false):flag(0), gf(NULL), gene_names(true, true, true),
                    transcripts(true,true,true) {
     if (gfrec->getGeneName())
        gene_names.Add(new GStr(gfrec->getGeneName()));
     transcripts.Add(new GStr(gfrec->getID()));
     create_gf(gfrec, gdata ,ensembl_convert);
   }

   ~GeneInfo() {
      delete gf;
   }

   void create_gf(GffObj* gfrec, GenomicSeqData* gdata, bool ensembl_convert) {
     gf=new GffObj(gfrec->getGeneID());
     GTData* gfdata=new GTData(gf, gdata);
     gfdata->geneinfo=this;
     gf->gseq_id=gfrec->gseq_id;
     gf->track_id=gfrec->track_id;
     gf->start=gfrec->start;
     gf->end=gfrec->end;
     gf->strand=gfrec->strand;
     gf->setFeatureName("gene");
     gf->isGene(true);
     gf->isUsed(true);
     //gf->uptr=gfdata; //for these new gene objects
     gfrec->incLevel();
     gfrec->parent=gf;
     gf->children.Add(gfrec);
     const char* s=NULL;
     if ((s=gfrec->getGeneName())) {
    	 gf->addAttr("Name", s);
    	 gf->copyAttrs(gfrec);
     }
     if (ensembl_convert) {
       //gf->addAttr("type", gf->getTrackName());
       const char* biotype=gfrec->getAttr("type");
       if (biotype) gf->addAttr("type", biotype);
       }
       // gf->children.Add(gfrec);
   }

   void update(GffObj* gfrec) {
     if (transcripts.AddedIfNew(new GStr(gfrec->getID()))<0)
       return;
     gene_names.AddedIfNew(new GStr(gfrec->getGeneName()));
     if (gf==NULL) {
        GError("GeneInfo::update() called on uninitialized gf!\n");
     }
     gfrec->parent=gf;
     gf->children.Add(gfrec);
     gfrec->incLevel();
     if (gf->start>gfrec->start)
           gf->start=gfrec->start;
     if (gf->end<gfrec->end)
           gf->end=gfrec->end;
     }

   void finalize() {
     //prepare attributes for printing
     //must be called right before printing
     if (gf==NULL || transcripts.Count()==0) return;
     if (gene_names.Count()>0) {
       gf->addAttr("Name", gene_names[0]->chars());
     } //has gene names
     GStr t(transcripts[0]->chars());
     for (int i=1;i<transcripts.Count();i++) {
          t.append(",");
          t.append(transcripts[i]->chars());
     }
     gf->addAttr("transcripts", t.chars());
   }
};
class GenomicSeqData {
  int gseq_id;
 public:
  const char* gseq_name;
  int seqreg_start; //if given by ##sequence-region comment
  int seqreg_end;
  GList<GffObj> gfs; //all non-transcript features -> usually gene features
  GList<GffObj> rnas; //all transcripts on this genomic sequence
  GList<GffLocus> loci; //all loci clusters
  GList<GTData> tdata; //transcript data (uptr holder for all rnas loaded here)
  uint64 f_bases;//base coverage on forward strand
  uint64 r_bases;//base coverage on reverse strand
  uint64 u_bases;//base coverage on undetermined strand
  //GenomicSeqData(int gid=-1):rnas(true,true,false),loci(true,true,true),
  GenomicSeqData(int gid=-1):gseq_id(gid), gseq_name(NULL), seqreg_start(0), seqreg_end(0),
		  gfs(true, true, false),rnas((GCompareProc*)gfo_cmpByLoc),loci(true,true,false),
		  tdata(false,true,false),  f_bases(0), r_bases(0), u_bases(0) {
  if (gseq_id>=0)
    gseq_name=GffObj::names->gseqs.getName(gseq_id);
  }

  bool operator==(GenomicSeqData& d){
    return gseq_id==d.gseq_id;
  }
  bool operator<(GenomicSeqData& d){
    return (gseq_id<d.gseq_id);
  }
};

class CGeneSym {
 public:
  GStr name;
  int freq;
  CGeneSym(const char* n=NULL, int f=0):name(n), freq(f) { }
  bool operator<(CGeneSym& b) {
     return (freq==b.freq)? ( (name.length()==b.name.length()) ? (name<b.name) :
         (name.length()<b.name.length()) ) : ( freq>b.freq );
     }
  bool operator==(CGeneSym& b) { return name==b.name; }
};

const char* getGeneDescr(const char* gsym);

void printLocus(GffLocus* loc, const char* pre=NULL);

class GffLocus:public GSeg {
public:
    int gseq_id; //id of underlying genomic sequence
    int locus_num;
    bool is_mrna;
    char strand;
    GffObj* t_maxcov;  //transcript with maximum coverage (for main "ref" transcript)
    GList<GffObj> gfs; //list of non-transcripts (genes) in this locus
    GList<GffObj> rnas; //list of transcripts (isoforms) for this locus
    GArray<GSeg> mexons; //list of merged exons in this region
    GList<CGeneSym> gene_names;
    GList<CGeneSym> gene_ids;
    int v; //user flag/data
   /*
   bool operator==(GffLocus& d){
       return (gseq_id==d.gseq_id && strand==d.strand && start==d.start && end==d.end);
       }
   bool operator<(GffLocus& d){
     if (gseq_id!=d.gseq_id) return (gseq_id<d.gseq_id);
     if (start==d.start) {
        if (end==d.end) return strand<d.strand;
                     else return end<d.end;
        } else return (start<d.start);
     }
    */
    const char* getGeneName() {
         if (gene_names.Count()==0) return NULL;
         return gene_names.First()->name.chars();
         }
    const char* get_tmax_id() {
         return t_maxcov->getID();
         }
    const char* get_descr() {
       if (gene_names.Count()>0) {
          for (int i=0;i<gene_names.Count();i++) {
            const char* gn=getGeneDescr(gene_names.First()->name.chars());
            if (gn!=NULL) return gn;
            }
          }
       char* s=t_maxcov->getAttr("product");
       if (s!=NULL) return s;
       s=t_maxcov->getAttr("descr");
       if (s!=NULL) return s;
       s=t_maxcov->getAttr("description");
       if (s!=NULL) return s;
       s=t_maxcov->getAttr("info");
       if (s!=NULL) return s;
       return NULL;
       }

    GffLocus(GffObj* t=NULL):gfs(true,false,false), rnas(true,false,false),mexons(true,true),
           gene_names(true,true,false), gene_ids(true,true,false) {
        //this will NOT free rnas!
        t_maxcov=NULL;
        gseq_id=-1;
        v=0;
        locus_num=0;
        start=0;
        end=0;
        strand=0;
        is_mrna=false;

        if (t!=NULL) {
           GSeg seg;
           bool is_t=(t->exons.Count()>0);
           if (is_t) {
             start=t->exons.First()->start;
             end=t->exons.Last()->end;
             for (int i=0;i<t->exons.Count();i++) {
               seg.start=t->exons[i]->start;
               seg.end=t->exons[i]->end;
               mexons.Add(seg);
             }
             rnas.Add(t);
           }
           else {
        	   start=t->start;
        	   end=t->end;
        	   seg.start=start;
        	   seg.end=end;
        	   mexons.Add(seg);
        	   gfs.Add(t);
           }
           gseq_id=t->gseq_id;
           ((GTData*)(t->uptr))->locus=this;
           t_maxcov=t;
           strand=t->strand;
           //if (t->ftype_id==gff_fid_mRNA) {
           if (t->isTranscript())
              is_mrna=true;
        }
    }

    void print(FILE *f, int idxfirstvalid, GStr& locname, GStr& loctrack) {
        const char* gseqname=NULL;
        if (rnas.Count()>0) gseqname=rnas[0]->getGSeqName();
        else gseqname=gfs[0]->getGSeqName();
        fprintf(f,"%s\t%s\tlocus\t%d\t%d\t.\t%c\t.\tID=%s",
                   gseqname, loctrack.chars(), this->start, this->end, this->strand,
                    locname.chars());
        //const char* loc_gname=loc.getGeneName();
        if (this->gene_names.Count()>0) { //print all gene names associated to this locus
             fprintf(f, ";genes=%s",this->gene_names.First()->name.chars());
             for (int i=1;i<this->gene_names.Count();i++) {
               fprintf(f, ",%s",this->gene_names[i]->name.chars());
             }
        }
        if (this->gene_ids.Count()>0) { //print all GeneIDs names associated to this locus
             fprintf(f, ";geneIDs=%s",this->gene_ids.First()->name.chars());
             for (int i=1;i<this->gene_ids.Count();i++) {
               fprintf(f, ",%s",this->gene_ids[i]->name.chars());
             }
        }
        if (idxfirstvalid>=0) {
        	GVec<int> tidx; //set of printable (non-discarded) rnas indexes
        	for (int i=idxfirstvalid;i<this->rnas.Count();i++)
        		if (((GTData*)this->rnas[i]->uptr)->replaced_by==NULL)
        			tidx.Add(i);
        	if (tidx.Count()>0) {
               fprintf(f, ";transcripts=%s",this->rnas[tidx[0]]->getID());
               for (int i=1;i<tidx.Count();i++)
                 fprintf(f, ",%s",this->rnas[tidx[i]]->getID());
        	}
        }
        fprintf(f, "\n");
    }

   void addMerge(GffLocus& locus, GffObj* lnkrna) {
     //add all the elements of the other locus (merging)
     //-- merge mexons
     GArray<int> ovlexons(true,true); //list of locus.mexons indexes overlapping existing mexons
     int i=0; //index of first mexons with a merge
     int j=0; //index current mrna exon
     while (i<mexons.Count() && j<locus.mexons.Count()) {
    	 uint istart=mexons[i].start;
    	 uint iend=mexons[i].end;
    	 uint jstart=locus.mexons[j].start;
    	 uint jend=locus.mexons[j].end;
    	 if (iend<jstart) { i++; continue; }
    	 if (jend<istart) { j++; continue; }
    	 ovlexons.Add(j);
    	 //extend mexons[i] as needed
		 if (jstart<istart) mexons[i].start=jstart;
		 if (jend>iend) { //mexons[i] end extend
			 mexons[i].end=jend;
			 //now this could overlap the next mexon(s), so we have to merge them all
			 while (i<mexons.Count()-1 && mexons[i].end>mexons[i+1].start) {
				 uint nextend=mexons[i+1].end;
				 mexons.Delete(i+1);
				 if (nextend>mexons[i].end) {
					 mexons[i].end=nextend;
					 break; //no need to check next mexons
				 }
			 } //while next mexons merge
		 } // mexons[i] end extend
		 j++; //check the next locus.mexon
     }
     //-- add the rest of the non-overlapping mexons:
     GSeg seg;
     for (int i=0;i<locus.mexons.Count();i++) {
            seg.start=locus.mexons[i].start;
            seg.end=locus.mexons[i].end;
            if (!ovlexons.Exists(i)) mexons.Add(seg);
     }
     // -- add locus.rnas
     for (int i=0;i<locus.rnas.Count();i++) {
          ((GTData*)(locus.rnas[i]->uptr))->locus=this;
          if (locus.rnas[i]!=lnkrna) rnas.Add(locus.rnas[i]);
     }
     for (int i=0;i<locus.gfs.Count();i++) {
          ((GTData*)(locus.gfs[i]->uptr))->locus=this;
          if (locus.gfs[i]!=lnkrna) gfs.Add(locus.gfs[i]);
     }

     // -- adjust start/end as needed
     if (start>locus.start) start=locus.start;
     if (end<locus.end) end=locus.end;
     if (locus.is_mrna) is_mrna=true;
     if (t_maxcov->covlen<locus.t_maxcov->covlen)
            t_maxcov=locus.t_maxcov;
    }

    bool add_gfobj(GffObj* t, bool adj) {
        //if (rnas.Count()==0) return true; //? should never be called on an empty locus
    	uint t_start=t->start;
    	uint t_end=t->end;
    	if (adj) {
    		t_start--;
    		t_end++;
    	}
        if (t->gseq_id!=gseq_id || /* t->strand!=strand || */ t_start>end || start>t_end)
              return false; //rna must be on the same genomic seq
        //check for exon overlap with existing mexons
        //also update mexons accordingly if t is to be added
        bool hasovl=false;
        if (t->exons.Count()>0) { //transcript-like entity
        	if (adj) {
        		t->exons.First()->start--;
        		t->exons.Last()->end++;
        	}
			int i=0; //index of first mexons with a merge
			int j=0; //index current t exon
			GArray<int> ovlexons(true,true); //list of mrna exon indexes overlapping mexons
			while (i<mexons.Count() && j<t->exons.Count()) {
				uint istart=mexons[i].start;
				uint iend=mexons[i].end;
				uint jstart=t->exons[j]->start;
				uint jend=t->exons[j]->end;
				if (iend<jstart) { i++; continue; }
				if (jend<istart) { j++; continue; }
				//exon overlap found if we're here:
				ovlexons.Add(j);
				hasovl=true;
				//extend mexons[i] as needed
				if (jstart<istart) mexons[i].start=jstart;
				if (jend>iend) { //mexon stretch up
					mexons[i].end=jend;
					//now this could overlap the next mexon(s), so we have to merge them all
					while (i<mexons.Count()-1 && mexons[i].end>mexons[i+1].start) {
						uint nextend=mexons[i+1].end;
						mexons.Delete(i+1);
						if (nextend>mexons[i].end) {
							mexons[i].end=nextend;
							break; //no need to check next mexons
						}
					} //while next mexons merge
				} //possible mexons merge

				j++; //check the next t exon
			}//all vs all exon check loop
        	if (adj) {
        		t->exons.First()->start++;
        		t->exons.Last()->end--;
        	}
	        if (hasovl) {
	            GSeg seg;
	             //add the rest of the non-overlapping exons
	            for (int i=0;i<t->exons.Count();i++) {
	                seg.start=t->exons[i]->start;
	                seg.end=t->exons[i]->end;
	                if (!ovlexons.Exists(i)) mexons.Add(seg);
	            }
	            t_add(t);
	            // add to rnas
	            ((GTData*)t->uptr)->locus=this;
	            gseq_id=t->gseq_id;
	       }
        } else {
        	//gene overlap check
			uint jstart=t->start;
			uint jend=t->end;
        	for (int i=0;i<mexons.Count();++i) {
				uint istart=mexons[i].start;
				uint iend=mexons[i].end;
				if (iend<jstart) continue;
				if (istart>jend) break;
				//exon overlap found:
				hasovl=true;
				//extend mexons[i] as needed
				if (jstart<istart) mexons[i].start=jstart;
				if (jend>iend) { //mexon stretch up
					mexons[i].end=jend;
					//now this could overlap the next mexon(s), so we have to merge them all
					while (i<mexons.Count()-1 && mexons[i].end>mexons[i+1].start) {
						uint nextend=mexons[i+1].end;
						mexons.Delete(i+1);
						if (nextend>mexons[i].end) {
							mexons[i].end=nextend;
							break; //no need to check next mexons
						}
					} //while next mexons merge
				} //possible mexons merge
        	}
        	if (hasovl) {
        		t_add(t); // add to locus rnas or gfs
        		((GTData*)t->uptr)->locus=this;
        		gseq_id=t->gseq_id;
        	}
        }
        return hasovl;
    }

    //basic adding of a GffObj to a locus
    void t_add(GffObj* t) {
      if (t->exons.Count()>0) rnas.Add(t);
      else gfs.Add(t);
      // adjust start/end
      //if (start==0 || start>t->start) start=t->start;
      if (start==0) start=t->start;
      else if (start>t->start) {
          start=t->start;
      }
      if (end<t->end) end=t->end;
      if (t_maxcov->covlen<t->covlen) t_maxcov=t;
      if (strand==0 || (strand=='.' && t->strand!='.')) strand=t->strand;
      //if (t->ftype_id==gff_fid_mRNA) is_mrna=true;
      if (t->isTranscript()) is_mrna=true;
    }
};


int gseqCmpName(const pointer p1, const pointer p2);

GenomicSeqData* getGSeqData(GList<GenomicSeqData>& seqdata, int gseq_id);

class GSpliceSite {
 public:
  char nt[3];
  GSpliceSite(const char* c, bool revc=false) {
    nt[2]=0;
    if (c==NULL) {
      nt[0]=0;
      nt[1]=0;
      return;
      }
    if (revc) {
      nt[0]=toupper(ntComplement(c[1]));
      nt[1]=toupper(ntComplement(c[0]));
      }
    else {
      nt[0]=toupper(c[0]);
      nt[1]=toupper(c[1]);
      }
    }

  GSpliceSite(const char* intron, int intronlen, bool getAcceptor, bool revc=false) {
    nt[2]=0;
    if (intron==NULL || intronlen==0)
       GError("Error: invalid intron or intron len for GSpliceSite()!\n");
    const char* c=intron;
    if (revc) {
      if (!getAcceptor) c+=intronlen-2;
      nt[0]=toupper(ntComplement(c[1]));
      nt[1]=toupper(ntComplement(c[0]));
      }
    else { //on forward strand
      if (getAcceptor) c+=intronlen-2;
      nt[0]=toupper(c[0]);
      nt[1]=toupper(c[1]);
      }//forward strand
    }

  GSpliceSite(const char n1, const char n2) {
    nt[2]=0;
    nt[0]=toupper(n1);
    nt[1]=toupper(n2);
    }
  bool canonicalDonor() {
    return (nt[0]=='G' && (nt[1]=='C' || nt[1]=='T'));
    }
  bool operator==(GSpliceSite& c) {
    return (c.nt[0]==nt[0] && c.nt[1]==nt[1]);
    }
  bool operator==(GSpliceSite* c) {
    return (c->nt[0]==nt[0] && c->nt[1]==nt[1]);
    }
  bool operator==(const char* c) {
    //return (nt[0]==toupper(c[0]) && nt[1]==toupper(c[1]));
    //assumes given const nucleotides are uppercase already!
    return (nt[0]==c[0] && nt[1]==c[1]);
    }
  bool operator!=(const char* c) {
    //assumes given const nucleotides are uppercase already!
    return (nt[0]!=c[0] || nt[1]!=c[1]);
    }
};

class GffLoader {
 public:
  GVec<char*> headerLines; //for GFF3 we keep the first few header lines (not the sequence-region one)
  GStr fname;
  FILE* f;
  GffNames* names;
  CIntronList intronList; // collect introns for -j output
  union {
	  unsigned int options;
	  struct {
		bool transcriptsOnly:1;
		bool gene2exon:1;
		bool fullAttributes:1;
		bool keep_AllExonAttrs:1;
		bool gatherExonAttrs:1;
		bool mergeCloseExons:1;
		bool ignoreLocus:1;
		bool noPseudo:1;
		bool BEDinput:1;
		bool TLFinput:1;
		bool keepGenes:1;
		bool trAdoption:1; //orphan transcript adoption by the container gene
		bool keepGff3Comments:1;
		bool sortRefsAlpha:1;
		bool doCluster:1;
		bool collapseRedundant:1; //discard "redundant" transcripts (-M/--merge activated)
		bool matchAllIntrons:1; //if true, contained transcripts are NOT discarded
		bool fuzzSpan:1; //matching/contained redundancy relaxed to disregard full boundary containment
		bool dOvlSET:1; //discard overlapping Single Exon Transcripts on any strand
		bool forceExons:1;
		bool streamIn:1;
		bool ensemblProc:1;
		bool attrsFilter:1;
	  };
  };

  GffLoader():fname(),f(NULL), names(NULL), intronList(), options(0) {
      transcriptsOnly=true;
      gffnames_ref(GffObj::names);
      names=GffObj::names;
  }

  void loadRefNames(GStr& flst);

  void openFile(GStr& file_name) {
	  //if (f!=NULL) closeFile();
	  fname=file_name;
      if (fname=="-" || fname=="stdin") {
         f=stdin;
         fname="stdin";
      }
      else {
         if ((f=fopen(fname.chars(), "r"))==NULL) {
           GError("Error: cannot open GFF file %s!\n",fname.chars());
         }
      }
  }

  bool validateGffRec(GffObj* gffrec);
  bool process_transcript(GFastaDb& gfasta, GffObj& gffrec);

  bool checkFilters(GffObj* gffrec);

  void collectIntrons(GffObj& t); //for -j output

  void load(GList<GenomicSeqData>&seqdata, GFFCommentParser* gf_parsecomment=NULL);

  bool placeGf(GffObj* t, GenomicSeqData* gdata);


  bool unsplContained(GffObj& ti, GffObj&  tj);
  GffObj* redundantTranscripts(GffObj& ti, GffObj&  tj);


  void terminate() {
	  //if (f!=NULL) closeFile(); GffReader is going to close the file
	  gffnames_unref(GffObj::names);
	  names=NULL;
  }
  void clearHeaderLines() {
	  if (headerLines.Count()>0) {
		  for (int i=0;i<headerLines.Count();i++) {
			  GFREE(headerLines[i]);
			  headerLines[i]=NULL;
		  }
	  }
  }
  ~GffLoader() {
	  this->terminate();
	  clearHeaderLines();
  }

};

void printFasta(FILE* f, GStr* defline, char* seq, int seqlen=-1, bool useStar=false);

//void printTabFormat(FILE* f, GffObj* t);

//"position" a given coordinate x within a list of transcripts sorted by their start (lowest)
//coordinate, using quick-search; the returned int is the list index of the closest *higher*
//GffObj - i.e. starting right *ABOVE* the given coordinate
//Convention: returns -1 if there is no such GffObj (i.e. last GffObj starts below x)
int qsearch_rnas(uint x, GList<GffObj>& rnas);
int qsearch_gloci(uint x, GList<GffLocus>& loci);

GffObj* redundantTranscripts(GffObj& ti, GffObj&  tj, bool matchAllIntrons=true, bool fuzzSpan=false);

//void loadGFF(FILE* f, GList<GenomicSeqData>& seqdata, const char* fname);

void collectLocusData(GList<GenomicSeqData>& ref_data, bool covInfo=false);

#endif
