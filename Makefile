GCLDIR := $(if $(GCLDIR),$(GCLDIR),../gclib)

SEARCHDIRS := -I. -I${GCLDIR}

SYSTYPE :=     $(shell uname)

CXX   := $(if $(CXX),$(CXX),g++)
LINKER  := $(if $(LINKER),$(LINKER),g++)

LDFLAGS := $(if $(LDFLAGS),$(LDFLAGS),-g)

BASEFLAGS  := -Wall -Wextra -std=c++11 ${SEARCHDIRS} -D_FILE_OFFSET_BITS=64 \
 -D_LARGEFILE_SOURCE -D_REENTRANT -fno-strict-aliasing \
 -fno-exceptions -fno-rtti

GCCV8 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 8)
ifeq "$(GCCV8)" "1"
 BASEFLAGS += -Wno-class-memaccess
endif

CXXFLAGS := $(if $(CXXFLAGS),$(BASEFLAGS) $(CXXFLAGS),$(BASEFLAGS))

ifneq (,$(filter %release %static, $(MAKECMDGOALS)))
  # -- release build
  LIBS := 
  ifneq (,$(findstring static,$(MAKECMDGOALS)))
    LDFLAGS += -static-libstdc++ -static-libgcc
  endif
  CXXFLAGS := -O3 -DNDEBUG $(CXXFLAGS)
else #debug builds
  ifneq (,$(filter %profile %gprof %prof, $(MAKECMDGOALS)))
    CXXFLAGS += -pg -O0 -DNDEBUG
    LDFLAGS += -pg
  else
    #CXXFLAGS += -g -O0 -DNDEBUG
    CXXFLAGS += -g -O0 -DDEBUG -D_DEBUG -DGDEBUG
  endif
  ifneq (,$(filter %memcheck %memdebug, $(MAKECMDGOALS)))
     #use sanitizer in gcc 4.9+
     MEMCHECK_BUILD := 1
     GCCVER49 := $(shell expr `${CXX} -dumpversion | cut -f1,2 -d.` \>= 4.9)
     ifeq "$(GCCVER49)" "0"
       $(error gcc version 4.9 or greater is required for this build target)
     endif
     CXXFLAGS += -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address
     GCCVER5 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 5)
     ifeq "$(GCCVER5)" "1"
       CXXFLAGS += -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=vptr
       CXXFLAGS += -fsanitize=float-cast-overflow -fsanitize=object-size
       #CXXFLAGS += -fcheck-pointer-bounds -mmpx
     endif
     CXXFLAGS += -fno-common -fstack-protector
     LIBS := -lasan -lubsan -ldl $(LIBS)
  else
     #just plain debug build
     DEBUG_BUILD := 1
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

OBJS := ${GCLDIR}/GBase.o ${GCLDIR}/GArgs.o ${GCLDIR}/GFaSeqGet.o \
 ${GCLDIR}/gdna.o ${GCLDIR}/codons.o ${GCLDIR}/gff.o ${GCLDIR}/GStr.o \
 ${GCLDIR}/GFastaIndex.o gff_utils.o
 
.PHONY : all

all static release debug memcheck memdebug profile gprof prof: ../gclib gffread

../gclib:
	git clone https://github.com/gpertea/gclib.git ../gclib

$(OBJS) : $(GCLDIR)/GBase.h $(GCLDIR)/gff.h
gffread.o : gff_utils.h $(GCLDIR)/GBase.h $(GCLDIR)/gff.h
gff_utils.o : gff_utils.h $(GCLDIR)/gff.h
${GCLDIR}/gff.o : ${GCLDIR}/gff.h ${GCLDIR}/GFaSeqGet.h ${GCLDIR}/GList.hh
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


