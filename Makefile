GCLDIR := ../gclib
SEARCHDIRS := -I. -I${GCLDIR}

SYSTYPE :=     $(shell uname)

MACHTYPE :=     $(shell uname -m)
ifeq ($(MACHTYPE), i686)
    MARCH = -march=i686
else
    MARCH = 
endif    

CC      := g++

BASEFLAGS  := -Wall -Wextra ${SEARCHDIRS} $(MARCH) -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -D_REENTRANT -fno-strict-aliasing -fno-exceptions -fno-rtti

#add the link-time optimization flag if gcc version > 4.5
GCC_VERSION:=$(subst ., ,$(shell gcc -dumpversion))
GCC_MAJOR:=$(word 1,$(GCC_VERSION))
GCC_MINOR:=$(word 2,$(GCC_VERSION))
#GCC_SUB:=$(word 3,$(GCC_VERSION))
GCC_SUB:=x

GCC45OPTS :=
GCC45OPTMAIN :=

ifeq ($(findstring release,$(MAKECMDGOALS)),release)
  CFLAGS := -O2 -DNDEBUG $(BASEFLAGS)
  LDFLAGS :=
else
  CFLAGS := -g -DDEBUG $(BASEFLAGS)
  LDFLAGS := -g
endif

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

# C/C++ linker

LINKER  := g++
LIBS := 
OBJS := ${GCLDIR}/GBase.o ${GCLDIR}/GArgs.o ${GCLDIR}/GFaSeqGet.o \
 ${GCLDIR}/gdna.o ${GCLDIR}/codons.o ${GCLDIR}/gff.o ${GCLDIR}/GStr.o \
 ${GCLDIR}/GFastaIndex.o gff_utils.o
 
.PHONY : all
all:    gffread

version: ; @echo "GCC Version is: "$(GCC_MAJOR)":"$(GCC_MINOR)":"$(GCC_SUB)
	@echo "> GCC Opt. string is: "$(GCC45OPTS)
debug:  gffread
nodebug: release
release: gffread

$(OBJS) : $(GCLDIR)/GBase.h $(GCLDIR)/gff.h
gffread.o : gff_utils.h $(GCLDIR)/GBase.h $(GCLDIR)/gff.h
gff_utils.o : gff_utils.h $(GCLDIR)/gff.h
${GCLDIR}/gff.o : ${GCLDIR}/gff.h ${GCLDIR}/GFaSeqGet.h ${GCLDIR}/GList.hh ${GCLDIR}/GHash.hh
${GCLDIR}/GFaSeqGet.o : ${GCLDIR}/GFaSeqGet.h
gffread: $(OBJS) gffread.o
	${LINKER} ${LDFLAGS} $(GCC45OPTS) $(GCC45OPTMAIN) -o $@ ${filter-out %.a %.so, $^} ${LIBS}

# target for removing all object files

.PHONY : clean
clean:: 
	@${RM} gffread gffread.o* gffread.exe $(OBJS)
	@${RM} core.*


