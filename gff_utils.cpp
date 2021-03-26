#include "gff_utils.h"

GHash<GeneInfo*> gene_ids;

bool verbose=false; //same with GffReader::showWarnings and GffLoader::beVserbose
bool debugMode=false;
bool ensembl_convert=false; //-L, assist in converting Ensembl GTF to GFF3

FILE* ffasta=NULL;
FILE* f_in=NULL;
FILE* f_out=NULL;
FILE* f_w=NULL; //writing fasta with spliced exons (transcripts)
int wPadding = 0; //padding for -w option
FILE* f_x=NULL; //writing fasta with spliced CDS
FILE* f_y=NULL; //wrting fasta with translated CDS

FILE* f_j=NULL; //wrting junctions (intron coordinates)

int maxintron=999000000;

ID_Flt_Type IDflt=idFlt_None;

bool TFilters=false;
bool wCDSonly=false;
bool wNConly=false;
int minLen=0; //minimum transcript length
bool validCDSonly=false; // translation with no in-frame STOP
bool bothStrands=false; //for single-exon mRNA validation, check the other strand too
bool altPhases=false; //if original phase fails translation validation,
                     //try the other 2 phases until one makes it
bool addCDSattrs=false;
bool add_hasCDS=false;
//bool streamIn=false; // --stream option
bool adjustStop=false; //automatic adjust the CDS stop coordinate
bool covInfo=false; // --cov-info : only report genome coverage
GStr tableFormat; //list of "attributes" to print in tab delimited format
bool spliceCheck=false; //only known splice-sites
bool decodeChars=false; //decode url-encoded chars in attrs (-D)
bool StarStop=false; //use * instead of . for stop codon translation
bool fullCDSonly=false; // starts with START, ends with STOP codon

bool multiExon=false;
bool writeExonSegs=false;
char* tracklabel=NULL;
/*
char* rfltGSeq=NULL;
char rfltStrand=0;
uint rfltStart=0;
uint rfltEnd=MAX_UINT;*/
GRangeParser* fltRange=NULL;

GRangeParser* fltJunction=NULL;

bool rfltWithin=false; //check for full containment within given range
bool addDescr=false;

bool wfaNoCDS=false;

bool fmtGFF3=true; //default output: GFF3
//other formats only make sense in transcriptOnly mode
bool fmtGTF=false;
bool fmtBED=false;
bool fmtTLF=false;
bool fmtTable=false;

GffPrintMode exonPrinting=pgffAny;

GFastaDb gfasta;

GHash<SeqInfo*> seqinfo;
GVec<CTableField> tableCols;
GHash<RefTran*> reftbl;
GStrSet<> fltIDs;

GStrSet<> attrList;

GHash<int> isoCounter; //counts the valid isoforms

void printFasta(FILE* f, GStr* defline, char* seq, int seqlen, bool useStar) {
 if (seq==NULL) return;
 int len=(seqlen>0)?seqlen:strlen(seq);
 if (len<=0) return;
 if (defline!=NULL)
     fprintf(f, ">%s\n",defline->chars());
 int ilen=0;
 for (int i=0; i < len; i++, ilen++) {
   if (ilen == 70) {
     fputc('\n', f);
     ilen = 0;
     }
   if (useStar && seq[i]=='.')
        putc('*', f);
   else putc(seq[i], f);
   } //for
 fputc('\n', f);
}

int qsearch_gloci(uint x, GList<GffLocus>& loci) {
  //binary search
  //do the simplest tests first:
  if (loci[0]->start>x) return 0;
  if (loci.Last()->start<x) return -1;
  uint istart=0;
  int i=0;
  int idx=-1;
  int maxh=loci.Count()-1;
  int l=0;
  int h = maxh;
  while (l <= h) {
     i = (l+h)>>1;
     istart=loci[i]->start;
     if (istart < x)  l = i + 1;
          else {
             if (istart == x) { //found matching coordinate here
                  idx=i;
                  while (idx<=maxh && loci[idx]->start==x) {
                     idx++;
                     }
                  return (idx>maxh) ? -1 : idx;
                  }
             h = i - 1;
             }
     } //while
 idx = l;
 while (idx<=maxh && loci[idx]->start<=x) {
    idx++;
    }
 return (idx>maxh) ? -1 : idx;
}

int qsearch_rnas(uint x, GList<GffObj>& rnas) {
  //binary search
  //do the simplest tests first:
  if (rnas[0]->start>x) return 0;
  if (rnas.Last()->start<x) return -1;
  uint istart=0;
  int i=0;
  int idx=-1;
  int maxh=rnas.Count()-1;
  int l=0;
  int h = maxh;
  while (l <= h) {
     i = (l+h)>>1;
     istart=rnas[i]->start;
     if (istart < x)  l = i + 1;
          else {
             if (istart == x) { //found matching coordinate here
                  idx=i;
                  while (idx<=maxh && rnas[idx]->start==x) {
                     idx++;
                     }
                  return (idx>maxh) ? -1 : idx;
                  }
             h = i - 1;
             }
     } //while
 idx = l;
 while (idx<=maxh && rnas[idx]->start<=x) {
    idx++;
    }
 return (idx>maxh) ? -1 : idx;
}

int cmpRedundant(GffObj& a, GffObj& b) {
  if (a.exons.Count()==b.exons.Count()) {
     if (a.covlen==b.covlen) {
       return strcmp(a.getID(), b.getID());
       }
     else return (a.covlen>b.covlen)? 1 : -1;
     }
   else return (a.exons.Count()>b.exons.Count())? 1: -1;
}

bool tMatch(GffObj& a, GffObj& b) {
  //strict intron chain match, or single-exon perfect match
  int imax=a.exons.Count()-1;
  int jmax=b.exons.Count()-1;
  int ovlen=0;
  if (imax!=jmax) return false; //different number of introns

  if (imax==0) { //single-exon mRNAs
    //if (equnspl) {
      //fuzz match for single-exon transfrags:
      // it's a match if they overlap at least 80% of max len
      ovlen=a.exons[0]->overlapLen(b.exons[0]);
      int maxlen=GMAX(a.covlen,b.covlen);
      return (ovlen>=maxlen*0.8);
    /*}
    else {
      //only exact match
      ovlen=a.covlen;
      return (a.exons[0]->start==b.exons[0]->start &&
          a.exons[0]->end==b.exons[0]->end);

       }*/
     }
  //check intron overlaps
  ovlen=a.exons[0]->end-(GMAX(a.start,b.start))+1;
  ovlen+=(GMIN(a.end,b.end))-a.exons.Last()->start;
  for (int i=1;i<=imax;i++) {
    if (i<imax) ovlen+=a.exons[i]->len();
    if ((a.exons[i-1]->end!=b.exons[i-1]->end) ||
      (a.exons[i]->start!=b.exons[i]->start)) {
            return false; //intron mismatch
    }
  }
  return true;
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


int adjust_stopcodon(GffObj& gffrec, int adj, GList<GSeg>* seglst) {
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

void printTableData(FILE* f, GffObj& g, bool inFasta) {
 //using attribute list in tableCols
	const int DBUF_LEN=1024; //there should not be attribute values larger than 1K!
	char dbuf[DBUF_LEN];
	char* av=NULL;
	for(int i=0;i<tableCols.Count();i++) {
		if (i>0 || inFasta) {
     	   if (!inFasta || tableCols[i].type!=ctfGFF_ID)
     		   fprintf(f,"\t");
		}
		switch(tableCols[i].type) {
		case ctfGFF_Attr:
			av=g.getAttr(tableCols[i].name.chars());
			if (av) {
				if (decodeChars) {
					GffObj::decodeHexChars(dbuf, av, DBUF_LEN-1);
					fprintf(f,"%s", dbuf);
				} else fprintf(f,"%s",av);
			} else fprintf(f,".");
			break;
		case ctfGFF_chr:
			fprintf(f,"%s",g.getGSeqName());
			break;
		case ctfGFF_ID:
			if (!inFasta)
			  fprintf(f,"%s",g.getID());
			break;
		case ctfGFF_geneID:
			fprintf(f,"%s",g.getGeneID()!=NULL ? g.getGeneID() : ".");
			break;
		case ctfGFF_geneName:
			fprintf(f,"%s",g.getGeneName()!=NULL ? g.getGeneName() : ".");
			break;
		case ctfGFF_Parent:
			fprintf(f,"%s",g.parent!=NULL ? g.parent->getID() : ".");
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
		} //switch
	}
	fprintf(f,"\n");
}

bool GffLoader::validateGffRec(GffObj* gffrec) {
	if (!checkFilters(gffrec)) {
		if (gffrec->isTranscript()) {
			TFilters=true;
			if (gffrec->parent!=NULL && keepGenes) {
			   GPVec<GffObj>& pchildren=gffrec->parent->children;
			   for (int c=0;c<pchildren.Count();c++) {
				   if (pchildren[c]==gffrec) {
					   pchildren.Delete(c);
					   break;
				   }
			   }
			}
		return false;
		}
		if (gffrec->isGene() && keepGenes) return true;
		return false;
	} //transcript rejected
	return true;
}

bool GffLoader::checkFilters(GffObj* gffrec) {
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
	if (transcriptsOnly && gffrec->isDiscarded()) {
		//discard generic "locus" features with no other detailed subfeatures
		//GMessage("Warning: discarding %s GFF generic gene/locus container %s\n",gffrec->getID());
		return false;
	}
	if (IDflt) {
		if (fltIDs.hasKey(gffrec->getID())) {
			if (IDflt==idFlt_Exclude) return false;
		}
		else if (IDflt==idFlt_Only) return false;
	}
	if (minLen>0 && gffrec->covlen<minLen) {
		if (verbose)
			GMessage("Info: %s discarded due to minimum length threshold %d\n",
					gffrec->getID(), minLen);
    	return false;
	}
	if (fltRange!=NULL) { //filter by gseqName
		if (fltRange->refName!=NULL && strcmp(gffrec->getGSeqName(),fltRange->refName)!=0) {
			return false;
		}
		if (fltRange->strand>0 && gffrec->strand!=fltRange->strand) {
			return false;
		}
		//check coordinates
		if (fltRange->start || fltRange->end<UINT_MAX) {
			if (rfltWithin) {
				if (gffrec->start<fltRange->start || gffrec->end>fltRange->end) {
					return false; //not within query range
				}
			}
			else {
				if (gffrec->start>fltRange->end || gffrec->end<fltRange->start) {
					return false;
				}
			}
		}
	}
	if (this->attrsFilter) { //mostly relevant for transcripts and gene records
		//remove attributes that are not in attrList
		gffrec->removeAttrs(attrList);
	}
    if (gffrec->isTranscript()) {    // && TFilters) ?
    	//these filters only apply to transcripts
		if (multiExon && gffrec->exons.Count()<=1) {
			return false;
		}
		if (wCDSonly && gffrec->CDstart==0) {
			return false;
		}
		if (wNConly && gffrec->hasCDS()) return false;
		if (fltJunction!=NULL) {
			if (gffrec->exons.Count()<=1) return false;
			if (fltJunction->refName!=NULL && strcmp(gffrec->getGSeqName(),fltJunction->refName)!=0) {
				return false;
			}
			if (fltJunction->strand && gffrec->strand!=fltJunction->strand) {
				return false;
			}
			//check coordinates
			uint jstart=fltJunction->start;
			uint jend=fltJunction->end;
			if (jstart==0) jstart=jend;
			if (jend==0)  jend=jstart;
			if (gffrec->start>=jstart || gffrec->end<=jend) {
						return false;
	        }

			bool noJMatch=true;
			for (int i=0;i<gffrec->exons.Count()-1;++i) {
				if (fltJunction->start && fltJunction->end) {
					if (gffrec->exons[i]->end+1==fltJunction->start &&
							gffrec->exons[i+1]->start-1==fltJunction->end)
						{ noJMatch=false; break; }
				} else if (fltJunction->start) { //end match not required
					if (gffrec->exons[i]->end+1==fltJunction->start)
						{ noJMatch=false; break; }
				} else { //only end match required:
					if (gffrec->exons[i+1]->start-1==fltJunction->end)
						{ noJMatch=false; break; }
				}
			}
			if (noJMatch) return false;
		}

		return process_transcript(gfasta, *gffrec);
    } //transcript filters check
	return true;
}

bool GffLoader::process_transcript(GFastaDb& gfasta, GffObj& gffrec) {
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
		   //isonum=new int(1);
		   isoCounter.Add(gname,1);
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
    if (gseq==NULL) {
    	GMessage("Error at GFF ID %s : could not retrieve subsequence %s:%d-%d !\n",
    			  gffrec.getID(), gffrec.getRefName(), gffrec.start, gffrec.end);
    	return false;
    }
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
    inframeStop=false;
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
	  GStr defline(gffrec.getID(), 94);
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
	  if (f_y!=NULL) { //CDS translation fasta output requested
			 if (cdsaa==NULL) { //translate now if not done before
			   cdsaa=translateDNA(cdsnt, aalen, seqlen);
			 }
			 if (aalen>0) {
			   if (cdsaa[aalen-1]=='.' || cdsaa[aalen-1]=='\0') --aalen; //avoid printing the stop codon
 			   fprintf(f_y, ">%s", defline.chars());
 			   if (fmtTable) printTableData(f_y, gffrec, true);
 			   else {
 				  if (gffrec.attrs!=NULL && gffrec.attrs->Count()>0) fprintf(f_y," ");
 				  gffrec.printAttrs(f_y, ";", false, decodeChars, false);
 				  fprintf(f_y, "\n");
 			   }
			   printFasta(f_y, NULL, cdsaa, aalen, StarStop);
			 }
	  }
	  if (f_x!=NULL) { //CDS only
			 fprintf(f_x, ">%s", defline.chars());
			 if (fmtTable) printTableData(f_x, gffrec, true);
			 else {
				 if (gffrec.attrs!=NULL && gffrec.attrs->Count()>0) fprintf(f_x," ");
				 gffrec.printAttrs(f_x, ";", false, decodeChars, false);
				 fprintf(f_x, "\n");
			 }
			 printFasta(f_x, NULL, cdsnt, seqlen);
	  }
	  GFREE(cdsnt);
	  GFREE(cdsaa);
  } //writing CDS or its translation
  if (f_w!=NULL) { //write spliced exons
	  uint cds_start=0;
	  uint cds_end=0;
	  seglst.Clear();
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
		  if (!wfaNoCDS && gffrec.CDstart>0) {
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

		  fprintf(f_w, ">%s", defline.chars());
		  if (fmtTable) printTableData(f_w, gffrec, true);
		    else {
		    	if (gffrec.attrs!=NULL && gffrec.attrs->Count()>0) fprintf(f_w," ");
		    	gffrec.printAttrs(f_w, ";", false, decodeChars, false);
		    	fprintf(f_w, "\n");
		    }
		  printFasta(f_w, NULL, exont, seqlen);
		  GFREE(exont);
	  }
  } //writing f_w (spliced exons)
  return true;
}


GTData::GTData(GffObj* t, GenomicSeqData* gd):rna(t),gdata(gd), locus(NULL), replaced_by(NULL), geneinfo(NULL) {
    if (rna!=NULL) {
        //geneinfo=(GeneInfo*)rna->uptr; //take over geneinfo, if there
        rna->uptr=this;
    }
    if (gdata!=NULL)
 	   gdata->tdata.Add(this);
}

bool GffLoader::unsplContained(GffObj& ti, GffObj&  tj) {
 //returns true only if ti (which MUST be single-exon) is "almost" contained in any of tj's exons
 //but it does not cross any intron-exon boundary of tj
  int imax=ti.exons.Count()-1;
  int jmax=tj.exons.Count()-1;
  if (imax>0) GError("Error: bad unsplContained() call, 1st parameter must be single-exon transcript!\n");
  if (fuzzSpan) {
    int maxIntronOvl=dOvlSET ? 25 : 0;
    //int minovl = dOvlSET ? 5 : (int)(0.8 * ti.len()); //minimum overlap to declare "redundancy"
    for (int j=0;j<=jmax;j++) {
       bool exonOverlap=false;
       if (dOvlSET) {
    	   exonOverlap= (tj.exons[j]->overlapLen(ti.start-1, ti.end+1) > 0);
       } else {
    	   exonOverlap=(ti.overlapLen(tj.exons[j])>=0.8 * ti.len());
       }
       if (exonOverlap) {
          //must not overlap the introns
          if ((j>0 && ti.start+maxIntronOvl<tj.exons[j]->start)
             || (j<jmax && ti.end>tj.exons[j]->end+maxIntronOvl))
             return false;
          return true;
       }
    } //for each exon
  } else { // not fuzzSpan, strict containment required
    for (int j=0;j<=jmax;j++) {
        if (ti.end<=tj.exons[j]->end && ti.start>=tj.exons[j]->start)
          return true;
    }
 }
 return false;
}

GffObj* GffLoader::redundantTranscripts(GffObj& ti, GffObj&  tj) {
  // matchAllIntrons==true:  transcripts are considered "redundant" only if
  //                   they have the exact same number of introns and same splice sites (or none)
  //                   (single-exon transcripts should be also fully contained to be considered matching)
  // matchAllIntrons==false: an intron chain could be a subset of a "container" chain,
  //                   as long as no intron-exon boundaries are violated; also, a single-exon
  //                   transcript will be collapsed if it's contained in one of the exons of the another transcript
  // fuzzSpan==false: the genomic span of one transcript MUST BE contained in or equal with the genomic
  //                  span of the other
  //
  // fuzzSpan==true: then genomic spans of transcripts are no longer required to be fully contained
  //                 (i.e. they may extend each-other in opposite directions)

  //if redundancy is detected, the "bigger" transcript is returned (otherwise NULL is returned)
 int adj=dOvlSET ? 1 : 0;
 if (ti.start>tj.end+adj || tj.start>ti.end+adj ||
		 (tj.strand!='.' && ti.strand!='.' && tj.strand!=ti.strand)) return NULL; //no span overlap
 int imax=ti.exons.Count()-1;
 int jmax=tj.exons.Count()-1;
 GffObj* bigger=NULL;
 GffObj* smaller=NULL;
 if (matchAllIntrons) { //full intron chain match expected, or full containment for SET
   if (imax!=jmax) return NULL; //must have the same number of exons!
   if (ti.covlen>tj.covlen) {
      bigger=&ti;
      if (!fuzzSpan && (ti.start>tj.start || ti.end<tj.end))
        return NULL; //no containment
   }
   else { //ti.covlen<=tj.covlen
      bigger=&tj;
      if (!fuzzSpan && (tj.start>ti.start || tj.end<ti.end))
         return NULL; //no containment
   }
   //check that all introns really match
   for (int i=0;i<imax;i++) {
     if (ti.exons[i]->end!=tj.exons[i]->end ||
         ti.exons[i+1]->start!=tj.exons[i+1]->start) return NULL;
     }
   return bigger;
 }
 //--- matchAllIntrons==false: intron-chain containment is also considered redundancy
 int minlen=0;
 if (ti.covlen>tj.covlen) {
      if (tj.exons.Count()>ti.exons.Count()) {
          //exon count override
          bigger=&tj;
          smaller=&ti;
      } else {
          bigger=&ti;
          smaller=&tj;
      }
      //maxlen=ti.covlen;
      minlen=tj.covlen;
 } else { //tj has more bases covered
      if (ti.exons.Count()>tj.exons.Count()) {
          //exon count override
          bigger=&ti;
          smaller=&tj;
      } else {
          bigger=&tj;
          smaller=&ti;
      }
      //maxlen=tj.covlen;
      minlen=ti.covlen;
 }
 if (imax==0 && jmax==0) {
     //single-exon transcripts: if fuzzSpan, at least 80% of the shortest one must be overlapped by the other
     if (fuzzSpan) {
       if (dOvlSET) {
           return (ti.exons[0]->overlapLen(tj.exons[0]->start-1, tj.exons[0]->end+1)>0) ? bigger : NULL;
       } else {
          return (ti.exons[0]->overlapLen(tj.exons[0])>=minlen*0.8) ? bigger : NULL;
       }
     } else { //boundary containment required
       return (smaller->start>=bigger->start && smaller->end<=bigger->end) ? bigger : NULL;
     }
 }
 //containment is also considered redundancy
 if (smaller->exons.Count()==1) {
   //check if this single exon is contained in any of tj exons
   //without violating any intron-exon boundaries
   return (unsplContained(*smaller, *bigger) ? bigger : NULL);
 }

 //--- from here on: both are multi-exon transcripts: imax>0 && jmax>0
  if (ti.exons[imax]->start<tj.exons[0]->end ||
     tj.exons[jmax]->start<ti.exons[0]->end )
         return NULL; //intron chains do not overlap at all
 //checking full intron chain containment
 uint eistart=0, eiend=0, ejstart=0, ejend=0; //exon boundaries
 int i=1; //exon idx to the right of the current intron of ti
 int j=1; //exon idx to the right of the current intron of tj
 //find the first intron overlap:
 while (i<=imax && j<=jmax) {
    eistart=ti.exons[i-1]->end;
    eiend=ti.exons[i]->start;
    ejstart=tj.exons[j-1]->end;
    ejend=tj.exons[j]->start;
    if (ejend<eistart) { j++; continue; }
    if (eiend<ejstart) { i++; continue; }
    //we found an intron overlap
    break;
 }
 if (!fuzzSpan && (bigger->start>smaller->start || bigger->end < smaller->end)) return NULL;
 if ((i>1 && j>1) || i>imax || j>jmax) {
     return NULL; //either no intron overlaps found at all
                  //or it's not the first intron for at least one of the transcripts
 }
 if (eistart!=ejstart || eiend!=ejend) return NULL; //not an exact intron match
 int maxIntronOvl=dOvlSET ? 25 : 0;
 if (j>i) {
   //i==1, ti's start must not conflict with the previous intron of tj
   if (ti.start+maxIntronOvl<tj.exons[j-1]->start) return NULL;
   //comment out the line above if you just want "intron compatibility" (i.e. extension of intron chains )
   //so i's first intron starts AFTER j's first intron
   // then j must contain i, so i's last intron must end with or before j's last intron
   if (ti.exons[imax]->start>tj.exons[jmax]->start) return NULL;
 }
 else if (i>j) {
   //j==1, tj's start must not conflict with the previous intron of ti
   if (tj.start+maxIntronOvl<ti.exons[i-1]->start) return NULL;
   //comment out the line above for just "intronCompatible()" check (allowing extension of intron chain)
   //so j's intron chain starts AFTER i's
   // then i must contain j, so j's last intron must end with or before j's last intron
   if (tj.exons[jmax]->start>ti.exons[imax]->start) return NULL;
 }
 //now check if the rest of the introns overlap, in the same sequence
 i++;
 j++;
 while (i<=imax && j<=jmax) {
   if (ti.exons[i-1]->end!=tj.exons[j-1]->end ||
      ti.exons[i]->start!=tj.exons[j]->start) return NULL;
   i++;
   j++;
 }
 i--;
 j--;
 if (i==imax && j<jmax) {
   // tj has more introns to the right, check if ti's end doesn't conflict with the current tj exon boundary
   if (ti.end>tj.exons[j]->end+maxIntronOvl) return NULL;
   }
 else if (j==jmax && i<imax) {
   if (tj.end>ti.exons[i]->end+maxIntronOvl) return NULL;
   }
 return bigger;
}


int gseqCmpName(const pointer p1, const pointer p2) {
 return strcmp(((GenomicSeqData*)p1)->gseq_name, ((GenomicSeqData*)p2)->gseq_name);
}


void printLocus(GffLocus* loc, const char* pre) {
  if (pre!=NULL) fprintf(stderr, "%s", pre);
  GMessage(" [%d-%d] : ", loc->start, loc->end);
  GMessage("%s",loc->rnas[0]->getID());
  for (int i=1;i<loc->rnas.Count();i++) {
    GMessage(",%s",loc->rnas[i]->getID());
    }
  GMessage("\n");
}

void preserveContainedCDS(GffObj* tcontainer, GffObj* t) {
 //transfer contained CDS info to the container if t has a CDS but container does not
 if (!t->hasCDS()) return;
  if (!tcontainer->hasCDS())//no CDS info on container, just copy it from the contained
	 tcontainer->setCDS(t);
}

bool exonOverlap2Gene(GffObj* t, GffObj& g) {
	if (t->exons.Count()>0) {
		return t->exonOverlap(g.start, g.end);
	}
	else return g.overlap(*t);
}
bool GffLoader::placeGf(GffObj* t, GenomicSeqData* gdata) {
  bool keep=false;
  GTData* tdata=NULL;
  //int tidx=-1;
  /*
  if (debug) {
     GMessage(">>Placing transcript %s\n", t->getID());
     debugState=true;
     }
    else debugState=false;
   */
  //dumb TRNA case for RefSeq: gene parent link missing
  //try to restore it here; BUT this only works if gene feature comes first
  ////DEBUG ONLY:
  //if (strcmp(t->getID(),"id24448")==0) { //&& t->start==309180) {
  //	 GMessage("placeGf %s (%d, %d) (%d exons)\n", t->getID(),t->start, t->end, t->exons.Count());
  //}
  //GMessage("DBG>>Placing transcript %s(%d-%d, %d exons)\n", t->getID(), t->start, t->end, t->exons.Count());

  if (t->parent==NULL && t->isTranscript() && trAdoption) {
  	int gidx=gdata->gfs.Count()-1;
  	while (gidx>=0 && gdata->gfs[gidx]->end>=t->start) {
  		GffObj& g = *(gdata->gfs[gidx]);
  		//try to find a container gene object for this transcript
  		//if (g.isGene() && t->strand==g.strand && exonOverlap2Gene(t, g)) {
  		if (g.isGene() && (t->strand=='.' || t->strand==g.strand) && g.exons.Count()==0
  				  && t->start>=g.start && t->end<=g.end) {
  			if (g.children.IndexOf(t)<0)
  				g.children.Add(t);
  			keep=true;
  			if (tdata==NULL) {
  		       tdata=new GTData(t, gdata); //additional transcript data
  			}
  			t->parent=&g;
  			//disable printing of gene if transcriptsOnly and --keep-genes wasn't given
  			if (transcriptsOnly && !keepGenes) {
  				T_NO_PRINT(g.udata); //tag it as non-printable
  				//keep gene ID and Name into transcript, when we don't print genes
  	  			const char* geneName=g.getAttr("Name");
  	  			if (t->getAttr("Name")==NULL && geneName) {
  	  				t->addAttr("Name", geneName);
  	  				if (t->getAttr("gene_name")==NULL)
  	  					t->addAttr("gene_name", geneName);
  	  			}
  	  			t->addAttr("geneID", g.getID());
  			}
  			break;
  		}
  		--gidx;
  	}
  }
  bool noexon_gfs=false;
  if (t->exons.Count()>0) { //treating this entry as a transcript
	gdata->rnas.Add(t); //added it in sorted order
	if (tdata==NULL) {
	   tdata=new GTData(t, gdata); //additional transcript data
	   //gdata->tdata.Add(tdata);
	}
	keep=true;
  }
   else {
    if (t->isGene() || !this->transcriptsOnly) {
	   gdata->gfs.Add(t);
	   keep=true;
	   if (tdata==NULL) {
		   tdata=new GTData(t, gdata); //additional transcript data
		   //gdata->tdata.Add(tdata);
	   }
	   noexon_gfs=true; //gene-like record, no exons defined
	   keep=true;
    } else {
       return false; //nothing to do with these non-transcript objects
    }
  }
  //keeping track of genes in special cases
	char* geneid=t->getGeneID();
	bool trackGenes=!t->isGene() && ( (keepGenes && t->parent==NULL) ||
			    (ensembl_convert && startsWith(t->getID(), "ENS") ) ) ;
	if (trackGenes) {
		GTData* tdata=(GTData*)(t->uptr);
		//keep track of chr|gene_id data and coordinate range
		if (geneid!=NULL) {
			GeneInfo* ginfo=gene_ids.Find(geneid);
			if (ginfo==NULL) {//first time seeing this gene ID
				GeneInfo* geneinfo=new GeneInfo(t, tdata->gdata, ensembl_convert);
				gene_ids.Add(geneid, geneinfo);
				//if (gfnew!=NULL) //new gene features
				//  gfnew->Add(geneinfo->gf);
			}
			else ginfo->update(t);
		}
	}



  if (!doCluster) return keep;

  if (!keep) return false;

  //---- place into a locus
  if (dOvlSET && t->exons.Count()==1) {
	  //for single exon transcripts temporarily set the strand to '.'
	  //so we can check both strands for overlap/locus
      T_SET_OSTRAND(t->udata, t->strand);
      t->strand='.';
  }
  if (gdata->loci.Count()==0) {
       gdata->loci.Add(new GffLocus(t));
       return true; //new locus on this ref seq
  }
  //--- look for any existing loci overlapping t
  uint t_end=t->end;
  uint t_start=t->start;
  if (dOvlSET) {
	  t_end++;
	  t_start--;
  }
  int nidx=qsearch_gloci(t_end, gdata->loci); //get index of nearest locus starting just ABOVE t->end
  //GMessage("\tlooking up end coord %d in gdata->loci.. (qsearch got nidx=%d)\n", t->end, nidx);
  if (nidx==0) {
     //cannot have any overlapping loci
     //if (debug) GMessage("  <<no ovls possible, create locus %d-%d \n",t->start, t->end);
     gdata->loci.Add(new GffLocus(t));
     return true;
  }
  if (nidx==-1) nidx=gdata->loci.Count();//all loci start below t->end
  int lfound=0; //count of parent loci
  GArray<int> mrgloci(false);
  GList<GffLocus> tloci(true); //candidate parent loci to adopt this
  //if (debug) GMessage("\tchecking all loci from %d to 0\n",nidx-1);
  for (int l=nidx-1;l>=0;l--) {
      GffLocus& loc=*(gdata->loci[l]);
      if ((loc.strand=='+' || loc.strand=='-') && t->strand!='.'&& loc.strand!=t->strand) continue;
      if (t_start>loc.end) {
           if (t->start-loc.start>GFF_MAX_LOCUS) break; //give up already
           continue;
      }
      if (loc.start>t_end) {
               //this should never be the case if nidx was found correctly
               GMessage("Warning: qsearch_gloci found loc.start>t.end!(t=%s)\n", t->getID());
               continue;
      }

      if (loc.add_gfobj(t, dOvlSET)) {
         //will add this transcript to loc
         lfound++;
         mrgloci.Add(l);
         if (collapseRedundant && !noexon_gfs) {
           //compare to every single transcript in this locus
           for (int ti=0;ti<loc.rnas.Count();ti++) {
                 if (loc.rnas[ti]==t) continue;
                 GTData* odata=(GTData*)(loc.rnas[ti]->uptr);
                 //GMessage("  ..redundant check vs overlapping transcript %s\n",loc.rnas[ti]->getID());
                 GffObj* container=NULL;
                 if (odata->replaced_by==NULL &&
                      (container=redundantTranscripts(*t, *(loc.rnas[ti])))!=NULL) {
                     if (container==t) {
                        odata->replaced_by=t;
                        preserveContainedCDS(t, loc.rnas[ti]);
                     }
                     else {// t is being replaced by previously defined transcript
                        tdata->replaced_by=loc.rnas[ti];
                        preserveContainedCDS(loc.rnas[ti], t);
                     }
                 }
           }//for each transcript in the exon-overlapping locus
         } //if doCollapseRedundant
      } //overlapping locus
  } //for each existing locus
  if (lfound==0) {
      //overlapping loci not found, create a locus with only this mRNA
      int addidx=gdata->loci.Add(new GffLocus(t));
      if (addidx<0) {
         //should never be the case!
         GMessage("  WARNING: new GffLocus(%s:%d-%d) not added!\n",t->getID(), t->start, t->end);
      }
   }
   else { //found at least one overlapping locus
     lfound--;
     int locidx=mrgloci[lfound];
     GffLocus& loc=*(gdata->loci[locidx]);
     //last locus index found is also the smallest index
     if (lfound>0) {
       //more than one loci found parenting this mRNA, merge loci
       /* if (debug)
          GMessage(" merging %d loci \n",lfound);
       */
       for (int l=0;l<lfound;l++) {
          int mlidx=mrgloci[l];
          loc.addMerge(*(gdata->loci[mlidx]), t);
          gdata->loci.Delete(mlidx); //highest indices first, so it's safe to remove
       }
     }
     int i=locidx;
     while (i>0 && loc<*(gdata->loci[i-1])) {
       //bubble down until it's in the proper order
       i--;
       gdata->loci.Swap(i,i+1);
     }
  }//found at least one overlapping locus
  return true;
}

void collectLocusData(GList<GenomicSeqData>& ref_data, bool covInfo) {
	int locus_num=0;
	for (int g=0;g<ref_data.Count();g++) {
		GenomicSeqData* gdata=ref_data[g];
		for (int l=0;l<gdata->loci.Count();l++) {
			GffLocus& loc=*(gdata->loci[l]);
			GHash<int> gnames; //gene names in this locus
			//GHash<int> geneids(true); //Entrez GeneID: numbers
			GHash<int> geneids;
			int fstrand=0,rstrand=0,ustrand=0;
			for (int i=0;i<loc.rnas.Count();i++) {
				GffObj& t=*(loc.rnas[i]);
				char tstrand=(char) T_OSTRAND(t.udata);
				if (tstrand==0) tstrand=t.strand;
				if (tstrand=='+') fstrand++;
				 else if (tstrand=='-') rstrand++;
				   else ustrand++;
				GStr gname(t.getGeneName());
				if (!gname.is_empty()) {
					gname.upper();
					int* prevg=gnames.Find(gname.chars());
					if (prevg!=NULL) (*prevg)++;
					else gnames.Add(gname.chars(), 1);
				}
				GStr geneid(t.getGeneID());
				if (!geneid.is_empty()) {
					int* prevg=gnames.Find(geneid.chars());
					if (prevg!=NULL) (*prevg)++;
					geneids.Add(geneid.chars(), 1);
				}
				//parse GeneID xrefs, if any (RefSeq):
				/*
				GStr xrefs(t.getAttr("xrefs"));
				if (!xrefs.is_empty()) {
					xrefs.startTokenize(",");
					GStr token;
					while (xrefs.nextToken(token)) {
						token.upper();
						if (token.startsWith("GENEID:")) {
							token.cut(0,token.index(':')+1);
							int* prevg=geneids.Find(token.chars());
							if (prevg!=NULL) (*prevg)++;
							else geneids.Add(token, new int(1));
						}
					} //for each xref
				} //xrefs parsing
				*/
			}//for each transcript
            if ((fstrand>0 && rstrand>0) ||
            		 (fstrand==0 && rstrand==0)) loc.strand='.';
            else if (fstrand==0 && rstrand>0) loc.strand='-';
            else loc.strand='+';
			for (int i=0;i<loc.gfs.Count();i++) {
				GffObj& nt=*(loc.gfs[i]);
				if (nt.isGene()) {
					GStr gname(nt.getGeneName());
					if (!gname.is_empty()) {
						gname.upper();
						int* prevg=gnames.Find(gname.chars());
						if (prevg!=NULL) (*prevg)++;
						else gnames.Add(gname, 1);
					}
					GStr geneid(nt.getID());
					if (!geneid.is_empty()) {
						int* prevg=gnames.Find(geneid.chars());
						if (prevg!=NULL) (*prevg)++;
						geneids.Add(geneid.chars(),1);
					}
				}
				//parse GeneID xrefs, if any (RefSeq):
				/*
				GStr xrefs(nt.getAttr("xrefs"));
				if (!xrefs.is_empty()) {
					xrefs.startTokenize(",");
					GStr token;
					while (xrefs.nextToken(token)) {
						token.upper();
						if (token.startsWith("GENEID:")) {
							token.cut(0,token.index(':')+1);
							int* prevg=geneids.Find(token.chars());
							if (prevg!=NULL) (*prevg)++;
							else geneids.Add(token, new int(1));
						}
					} //for each xref
				} //xrefs parsing
				*/
			}//for each non-transcript (genes?)
			if (covInfo) {
				for (int m=0;m<loc.mexons.Count();m++) {
					if (loc.strand=='+')
						gdata->f_bases+=loc.mexons[m].len();
					else if (loc.strand=='-')
						gdata->r_bases+=loc.mexons[m].len();
					else gdata->u_bases+=loc.mexons[m].len();
				}
			}
			locus_num++;
			loc.locus_num=locus_num;
			if (gnames.Count()>0) { //collect all gene names associated to this locus
				gnames.startIterate();
				int gfreq=0;
				const char* key=NULL;
				while ((key=gnames.Next(gfreq))!=NULL) {
					loc.gene_names.AddIfNew(new CGeneSym(key, gfreq));
				}
			} //added collected gene_names
			if (geneids.Count()>0) { //collect all GeneIDs names associated to this locus
				geneids.startIterate();
				int gfreq=0;
				const char* key=NULL;
				while ((key=geneids.Next(gfreq))!=NULL) {
					loc.gene_ids.AddIfNew(new CGeneSym(key, gfreq));
				}
			}
		} //for each locus
	}//for each genomic sequence
}

void GffLoader::loadRefNames(GStr& flst) {
 //load the whole file and split by (' \t\n\r,'
	int64_t fsize=fileSize(flst.chars());
	if (fsize<0) GError("Error: could not get file size for %s !\n",
			flst.chars());
	GStr slurp("", fsize+1);
	//sanity check for file size?
	FILE* f=fopen(flst.chars(), "r");
	if (f==NULL)
		GError("Error: could not open file %s !\n", flst.chars());
	slurp.read(f, NULL);
	fclose(f);
	slurp.startTokenize(" ,;\t\r\n", tkCharSet);
	GStr refname;
	while (slurp.nextToken(refname)) {
		if (refname.is_empty()) continue;
		names->gseqs.addName(refname.chars());
	}
}

GenomicSeqData* getGSeqData(GList<GenomicSeqData>& seqdata, int gseq_id) {
	int i=-1;
	GenomicSeqData f(gseq_id);
	GenomicSeqData* gdata=NULL;
	if (seqdata.Found(&f,i)) gdata=seqdata[i];
	else { //entry not created yet for this genomic seq
		gdata=new GenomicSeqData(gseq_id);
		seqdata.Add(gdata);
	}
	return gdata;
}

void warnPseudo(GffObj& m) {
	GMessage("Info: pseudo gene/transcript record with ID=%s discarded.\n",m.getID());
}
void GffLoader::collectIntrons(GffObj& t) {
	// assume t are coming grouped by chromosome and sorted by start coordinate!
	if (intronList.jlst.Count()>0) {
         if (t.gseq_id!=intronList.gseq_id ||
        		 t.start>=intronList.jlst.Last()->start) {
        	intronList.print(f_j);
        	intronList.clear();
         }
         else if (t.start<intronList.last_t_start)
        	 GError("Error collectIntrons(%s) called when last_t_start was %d\n",
        			 t.getID(), intronList.last_t_start );
         //add this transcript's introns
	}
    intronList.add(t);
}

void GffLoader::load(GList<GenomicSeqData>& seqdata, GFFCommentParser* gf_parsecomment) {
	if (f==NULL) GError("Error: GffLoader::load() cannot be called before ::openFile()!\n");
	GffReader* gffr=new GffReader(f, this->transcriptsOnly, true); //not only mRNA features, sorted
	clearHeaderLines();
	gffr->showWarnings(verbose);
	//           keepAttrs   mergeCloseExons  noExonAttr
	gffr->gene2Exon(gene2exon);
	if (BEDinput) gffr->isBED(true);
	//if (TLFinput) gffr->isTLF(true);
	gffr->mergeCloseExons(mergeCloseExons);
	gffr->keepAttrs(fullAttributes, gatherExonAttrs, keep_AllExonAttrs);
	gffr->keepGenes(keepGenes);
	gffr->setIgnoreLocus(ignoreLocus);
	gffr->setRefAlphaSorted(this->sortRefsAlpha);
	gffr->procEnsemblID(this->ensemblProc);
	if (keepGff3Comments && gf_parsecomment!=NULL) gffr->setCommentParser(gf_parsecomment);
    int outcounter=0;
	if (streamIn) { //this will ignore any clustering options
		GffObj* t=NULL;
		while ((t=gffr->readNext())!=NULL) {
			if (!validateGffRec(t)) {
				delete t;
				continue;
			}

			if (f_j!=NULL && t->isTranscript() && t->exons.Count()>1)
				collectIntrons(*t);

			outcounter++;
			if (f_out) {
			  if (fmtTable)
					printTableData(f_out, *t);
			  else //GFF3, GTF, BED, TLF
				t->printGxf(f_out, exonPrinting, tracklabel, NULL, decodeChars);
			}
			delete t;
		}
		if (f_j && intronList.jlst.Count()>0) {
        	intronList.print(f_j);
        	intronList.clear();
		}
		delete gffr;
		return;
	}

	gffr->readAll();
	GVec<int> pseudoFeatureIds; //feature type: pseudo*
	GVec<int> pseudoAttrIds;  // attribute: [is]pseudo*=true/yes/1
	GVec<int> pseudoTypeAttrIds;  // attribute: *_type=pseudo*

	if (this->noPseudo) {
		GffNameList& fnames = GffObj::names->feats; //gffr->names->feats;
		for (int i=0;i<fnames.Count();i++) {
			char* n=fnames[i]->name;
			if (startsWith(n, "pseudo")) {
				pseudoFeatureIds.Add(fnames[i]->idx);
			}
		}
		GffNameList& attrnames = GffObj::names->attrs;//gffr->names->attrs;
		for (int i=0;i<attrnames.Count();i++) {
			char* n=attrnames[i]->name;
			if (endsiWith(n, "type")) {
				pseudoTypeAttrIds.Add(attrnames[i]->idx);
			}// else {
			char* p=strifind(n, "pseudo");
			if (p==n || (p==n+2 && tolower(n[0])=='i' && tolower(n[1])=='s') ||
					(p==n+3 && startsiWith(n, "is_")) ) {
				pseudoAttrIds.Add(attrnames[i]->idx);
			}
			//}
		}
	}

	if (verbose) GMessage("   .. loaded %d genomic features from %s\n", gffr->gflst.Count(), fname.chars());
	//add to GenomicSeqData, adding to existing loci and identifying intron-chain duplicates
	for (int k=0;k<gffr->gflst.Count();k++) {
		GffObj* m=gffr->gflst[k];
		if (strcmp(m->getFeatureName(), "locus")==0 &&
				m->getAttr("transcripts")!=NULL) {
			continue; //discard locus meta-features
		}
		if (this->noPseudo) {
			bool is_pseudo=false;
			for (int i=0;i<pseudoFeatureIds.Count();++i) {
				if (pseudoFeatureIds[i]==m->ftype_id) {
					is_pseudo=true;
					break;
				}
			}
			if (is_pseudo) {
				if (verbose) warnPseudo(*m);
				continue;
			}
			for (int i=0;i<pseudoAttrIds.Count();++i) {
				char* attrv=NULL;
				if (m->attrs!=NULL) attrv=m->attrs->getAttr(pseudoAttrIds[i]);
				if (attrv!=NULL) {
					char fc=tolower(attrv[0]);
					if (fc=='t' || fc=='y' || fc=='1') {
						is_pseudo=true;
						break;
					}
				}
			}
			if (is_pseudo) {
				if (verbose) warnPseudo(*m);
				continue;
			}
			//  *type=*_pseudogene
            //find all attributes ending with _type and have value like: *_pseudogene
			for (int i=0;i<pseudoTypeAttrIds.Count();++i) {
				char* attrv=NULL;
				if (m->attrs!=NULL) attrv=m->attrs->getAttr(pseudoTypeAttrIds[i]);
				if (attrv!=NULL &&
						(startsWith(attrv, "pseudogene") || endsWith(attrv, "_pseudogene")) ) {
					is_pseudo=true;
					break;
				}
			}
			if (is_pseudo) {
				if (verbose) warnPseudo(*m);
				continue;
			}
		} //pseudogene detection requested
		char* rloc=m->getAttr("locus");
		if (rloc!=NULL && startsWith(rloc, "RLOC_")) {
			m->removeAttr("locus", rloc);
		}
		if (forceExons) {
			m->subftype_id=gff_fid_exon;
		}
		//GList<GffObj> gfadd(false,false); -- for gf_validate()?
		if (!validateGffRec(m)) { //this will also apply process_transcript() CDS filters etc.
			continue;
		}
		if (f_j!=NULL && m->isTranscript() && m->exons.Count()>1)
			collectIntrons(*m);

		m->isUsed(true); //so the gffreader won't destroy it
		GenomicSeqData* gdata=getGSeqData(seqdata, m->gseq_id);
		bool keep=placeGf(m, gdata);
		if (!keep) {
			m->isUsed(false);
			//DEBUG
			//GMessage("Feature %s(%d-%d) is going to be discarded..\n",m->getID(), m->start, m->end);
		}
	} //for each read gffObj
	if (f_j && intronList.jlst.Count()>0) {
    	intronList.print(f_j);
    	intronList.clear();
	}
	//if (verbose) GMessage("  .. %d records from %s clustered into loci.\n", gffr->gflst.Count(), fname.chars());
	//if (f && f!=stdin) { fclose(f); f=NULL; }
	delete gffr;
}
