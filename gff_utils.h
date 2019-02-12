#ifndef GFF_UTILS_H
#define GFF_UTILS_H
#include "gff.h"
#include "GStr.h"
#include "GFaSeqGet.h"

extern bool debugMode;

typedef bool GFValidateFunc(GffObj* gf, GList<GffObj>* gfadd);

class GeneInfo { //for Ensembl GTF conversion
 public:
   int flag;
   GffObj* gf;
   GList<GStr> gene_names;
   GList<GStr> transcripts; //list of transcript IDs
   GeneInfo():gene_names(true, true, true), transcripts(true,true,true) {
     gf=NULL;
     flag=0;
     }
   GeneInfo(GffObj* gfrec, bool ensembl_convert=false):gene_names(true, true, true),
                    transcripts(true,true,true) {
     flag=0;
     if (gfrec->getGeneName())
        gene_names.Add(new GStr(gfrec->getGeneName()));
     transcripts.Add(new GStr(gfrec->getID()));
     create_gf(gfrec, ensembl_convert);
     }

   void create_gf(GffObj* gfrec, bool ensembl_convert) {
     gf=new GffObj(gfrec->getGeneID());
     gf->gseq_id=gfrec->gseq_id;
     gf->track_id=gfrec->track_id;
     gf->start=gfrec->start;
     gf->end=gfrec->end;
     gf->strand=gfrec->strand;
     gf->setFeatureName("gene");
     gf->isGene(true);
     gf->isUsed(true);
     gf->uptr=this;
     gfrec->incLevel();
     gfrec->parent=gf;
     gf->children.Add(gfrec);
     if (ensembl_convert) {
       //gf->addAttr("type", gf->getTrackName());
       const char* biotype=gfrec->getAttr("type");
       if (biotype) gf->addAttr("type", biotype);
       }
     //gf->children.Add(gfrec);
     }
   //~GeneInfo() {
   //  }
   void update(GffObj* gfrec) {
     if (transcripts.AddedIfNew(new GStr(gfrec->getID()))<0)
       return;
     gene_names.AddedIfNew(new GStr(gfrec->getGeneName()));
     if (gf==NULL) {
        GError("GeneInfo::update() called on uninitialized gf!\n");
        //create_gf(gfrec);
        //return;
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
       /*
       GStr s(gene_names[0]->chars());
       for (int i=1;i<gene_names.Count();i++) {
          s.append(",");
          s.append(gene_names[i]->chars());
          }
       gf->addAttr("genes", s.chars());
       */
       } //has gene names
       GStr t(transcripts[0]->chars());
       for (int i=1;i<transcripts.Count();i++) {
          t.append(",");
          t.append(transcripts[i]->chars());
          }
       gf->addAttr("transcripts", t.chars());
     }
};


class GffLocus;

class GTData { //transcript associated data
 public:
    GffObj* rna;
    GffLocus* locus;
    GffObj* replaced_by;
    GeneInfo* geneinfo;
    int flag;
    GTData(GffObj* t=NULL) {
        rna=t;
        flag=0;
        locus=NULL;
        replaced_by=NULL;
        geneinfo=NULL;
        if (rna!=NULL) {
            geneinfo=(GeneInfo*)rna->uptr; //take over geneinfo, if there
            rna->uptr=this;
            }
        }
   bool operator<(GTData& b) { return (rna < b.rna); }
   bool operator==(GTData& b) { return (rna==b.rna); }
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
           if (t->ftype_id==gff_fid_mRNA) {
              is_mrna=true;
           }
        }
    }

    void print(FILE *f, int idxfirstvalid, GStr& locname, GStr& loctrack) {
        const char* gseqname=NULL;
        if (rnas.Count()>0) gseqname=rnas[0]->getGSeqName();
        else gseqname=gfs[0]->getGSeqName();
        fprintf(f,"%s\t%s\tlocus\t%d\t%d\t.\t%c\t.\tID=%s;locus=%s",
                   gseqname, loctrack.chars(), this->start, this->end, this->strand,
                    locname.chars(), locname.chars());
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
            fprintf(f, ";transcripts=%s",this->rnas[idxfirstvalid]->getID());
            for (int i=idxfirstvalid+1;i<this->rnas.Count();i++) {
              fprintf(f, ",%s",this->rnas[i]->getID());
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

    bool exonOverlap(GffLocus& loc) {
        //check if any mexons overlap!
        if (strand!=loc.strand || loc.start>end || start>loc.end) return false;
        int i=0;
        int j=0;
        while (i<mexons.Count() && j<loc.mexons.Count()) {
            uint istart=mexons[i].start;
            uint iend=mexons[i].end;
            uint jstart=loc.mexons[j].start;
            uint jend=loc.mexons[j].end;
            if (iend<jstart) { i++; continue; }
            if (jend<istart) { j++; continue; }
            //exon overlap found if we're here:
            return true;
        }
        return false;
    }

    bool add_gfobj(GffObj* t) {
        //if (rnas.Count()==0) return true; //? should never be called on an empty locus
        if (t->gseq_id!=gseq_id || t->strand!=strand || t->start>end || start>t->end)
              return false; //rna must be on the same genomic seq
        //check for exon overlap with existing mexons
        //also update mexons accordingly if t is to be added
        bool hasovl=false;
        if (t->exons.Count()>0) {
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
      if (strand==0) strand=t->strand;
      if (t->ftype_id==gff_fid_mRNA) is_mrna=true;
      }
};

class GenomicSeqData {
  int gseq_id;
 public:
  const char* gseq_name;
  GList<GffObj> gfs; //all non-transcript features -> usually gene features
  GList<GffObj> rnas; //all transcripts on this genomic sequence
  GList<GffLocus> loci; //all loci clusters
  GList<GTData> tdata; //transcript data (uptr holder for all rnas loaded here)
  uint64 f_bases;//base coverage on forward strand
  uint64 r_bases;//base coverage on reverse strand
  uint64 u_bases;//base coverage on undetermined strand
  //GenomicSeqData(int gid=-1):rnas(true,true,false),loci(true,true,true),
  GenomicSeqData(int gid=-1):gfs(true, true, false),rnas((GCompareProc*)gfo_cmpByLoc),loci(true,true,false),
       tdata(false,true,false),  f_bases(0), r_bases(0), u_bases(0) {
  gseq_id=gid;
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

int gseqCmpName(const pointer p1, const pointer p2);

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

struct GffLoader {
  GStr fname;
  FILE* f;
  bool transcriptsOnly;
  bool gene2exon;
  bool fullAttributes;
  bool gatherExonAttrs;
  bool mergeCloseExons;
  bool showWarnings;
  bool noPseudo;
  bool BEDinput;
  bool TLFinput;
  bool placeGf(GffObj* t, GenomicSeqData* gdata, bool doCluster=true, bool collapseRedundant=true,
                                    bool matchAllIntrons=true, bool fuzzSpan=false);
  void load(GList<GenomicSeqData>&seqdata, GFValidateFunc* gf_validate=NULL,
                      bool doCluster=true, bool doCollapseRedundant=true,
                      bool matchAllIntrons=true, bool fuzzSpan=false, bool forceExons=false);
  GffLoader(const char* filename):fname(filename) {
      f=NULL;
      transcriptsOnly=true;
      gene2exon=false;
      fullAttributes=false;
      gatherExonAttrs=false;
      mergeCloseExons=false;
      showWarnings=false;
      noPseudo=false;
      BEDinput=false;
      TLFinput=false;
      if (fname=="-" || fname=="stdin") {
         f=stdin;
         fname="stdin";
         }
        else {
          if ((f=fopen(fname.chars(), "r"))==NULL) {
            GError("Error: cannot open gff file %s!\n",fname.chars());
            }
          }
      }
};

void printFasta(FILE* f, GStr& defline, char* seq, int seqlen=-1, bool useStar=false);

void printTabFormat(FILE* f, GffObj* t);

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
