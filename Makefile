GCLDIR := ../gclib
SEARCHDIRS := -I. -I${GCLDIR}

SYSTYPE :=     $(shell uname)

MACHTYPE :=     $(shell uname -m)
ifeq ($(MACHTYPE), i686)
    MARCH = -march=i686
else
    MARCH = 
endif    

CXX   := $(if $(CXX),$(CXX),g++)
LINKER  := $(if $(LINKER),$(LINKER),g++)

LDFLAGS := $(if $(LDFLAGS),$(LDFLAGS),)

#CC      := g++

BASEFLAGS  := -Wall -Wextra ${SEARCHDIRS} $(MARCH) -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -D_REENTRANT -fno-strict-aliasing -fno-exceptions -fno-rtti

ifneq (,$(filter %release %static, $(MAKECMDGOALS)))
  # -- release build
  #RELEASE_BUILD=1
  CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS), -O3)
  CXXFLAGS += -DNDEBUG $(BASEFLAGS)
else
  ifneq (,$(filter %memcheck %memdebug, $(MAKECMDGOALS)))
     #use sanitizer in gcc 4.9+
     MEMCHECK_BUILD=1
     GCCVER49 := $(shell expr `g++ -dumpversion | cut -f1,2 -d.` \>= 4.9)
     ifeq "$(GCCVER49)" "0"
       $(error gcc version 4.9 or greater is required for this build target)
     endif
     CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O0)
     CXXFLAGS += -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address $(BASEFLAGS)
     GCCVER5 := $(shell expr `g++ -dumpversion | cut -f1 -d.` \>= 5)
     ifeq "$(GCCVER5)" "1"
       CXXFLAGS += -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=vptr
       CXXFLAGS += -fsanitize=float-cast-overflow -fsanitize=object-size
       #CXXFLAGS += -fcheck-pointer-bounds -mmpx
     endif
     CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG -fno-common -fstack-protector
     LIBS := -lasan -lubsan -ldl $(LIBS)
  else
     #just plain debug build
     DEBUG_BUILD=1
     CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O0)
     CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG $(BASEFLAGS)
  endif
endif

#ifneq (,$(filter %memtrace %memusage %memuse, $(MAKECMDGOALS)))
#    CXXFLAGS += -DGMEMTRACE
#    OBJS += ${GDIR}/proc_mem.o
#endif

#ifdef DEBUG_BUILD
#  #$(warning Building DEBUG version of stringtie.. )
#  DBG_WARN=@echo
#  DBG_WARN+='WARNING: built DEBUG version [much slower], use "make clean release" for a faster, optimized version of the program.'
#endif

%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

# C/C++ linker

LINKER  := g++
LIBS := 
OBJS := ${GCLDIR}/GBase.o ${GCLDIR}/GArgs.o ${GCLDIR}/GFaSeqGet.o \
 ${GCLDIR}/gdna.o ${GCLDIR}/codons.o ${GCLDIR}/gff.o ${GCLDIR}/GStr.o \
 ${GCLDIR}/GFastaIndex.o gff_utils.o
 
.PHONY : all

nodebug: release
all release debug memcheck memdebug: gffread

$(OBJS) : $(GCLDIR)/GBase.h $(GCLDIR)/gff.h
gffread.o : gff_utils.h $(GCLDIR)/GBase.h $(GCLDIR)/gff.h
gff_utils.o : gff_utils.h $(GCLDIR)/gff.h
${GCLDIR}/gff.o : ${GCLDIR}/gff.h ${GCLDIR}/GFaSeqGet.h ${GCLDIR}/GList.hh ${GCLDIR}/GHash.hh
${GCLDIR}/GFaSeqGet.o : ${GCLDIR}/GFaSeqGet.h
gffread: $(OBJS) gffread.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
#	@echo
#	${DBG_WARN}


# target for removing all object files

.PHONY : clean
clean:
	@${RM} gffread gffread.o* gffread.exe $(OBJS)
	@${RM} core.*


