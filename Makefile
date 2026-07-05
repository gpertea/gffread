GCLDIR := $(if $(GCLDIR),$(GCLDIR),./gclib)

SEARCHDIRS := -I.
ifdef STRICT_COORDS
SEARCHDIRS += -isystem ${GCLDIR}
else
SEARCHDIRS += -I${GCLDIR}
endif

SYSTYPE :=     $(shell uname)

CXX   := $(if $(CXX),$(CXX),g++)
LINKER  := $(if $(LINKER),$(LINKER),g++)

LDFLAGS := $(if $(LDFLAGS),$(LDFLAGS),-g)
LIBS := -lz

BASEFLAGS  := -Wall -Wextra -std=c++11 ${SEARCHDIRS} -D_FILE_OFFSET_BITS=64 \
 -D_LARGEFILE_SOURCE -D_REENTRANT -fno-strict-aliasing \
 -fno-exceptions -fno-rtti
STRICT_CHECK_FLAGS := -Wconversion -Wsign-conversion
STRICT_COORD_CXXFLAGS :=

ifdef STRICT_COORDS
STRICT_COORD_CXXFLAGS += $(STRICT_CHECK_FLAGS)
endif

GCCV8 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 8)
ifeq "$(GCCV8)" "1"
 BASEFLAGS += -Wno-class-memaccess
endif

CXXFLAGS := $(if $(CXXFLAGS),$(BASEFLAGS) $(CXXFLAGS),$(BASEFLAGS))

ifneq (,$(filter %release %static, $(MAKECMDGOALS)))
  # -- release build
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
     LIBS += -lasan -lubsan -ldl
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
 ${GCLDIR}/GFastaIndex.o ${GCLDIR}/GBgzf.o gff_utils.o
 
.PHONY : all gclib-init strict-coords large-tests

all static release debug memcheck memdebug profile gprof prof: gclib-init gffread

gclib-init:
	@if [ ! -f "${GCLDIR}/GBase.h" ]; then \
	  if [ "${GCLDIR}" = "./gclib" ] && [ -d .git ]; then \
	    git submodule sync -- gclib; \
	    git submodule update --init --checkout gclib; \
	    test -f "${GCLDIR}/GBase.h" || { \
	      echo "Error: gclib submodule init failed"; \
	      exit 1; \
	    }; \
	  else \
	    echo "Error: ${GCLDIR}/GBase.h not found"; \
	    echo "Hint: clone with --recurse-submodules or run: git submodule update --init gclib"; \
	    exit 1; \
	  fi; \
	fi

$(GCLDIR)/GBase.h $(GCLDIR)/gff.h:
	@$(MAKE) --no-print-directory gclib-init

$(OBJS) : $(GCLDIR)/GBase.h $(GCLDIR)/gff.h
gffread.o : gff_utils.h $(GCLDIR)/GBase.h $(GCLDIR)/gff.h
gff_utils.o : gff_utils.h $(GCLDIR)/gff.h
gff_utils.o gffread.o : CXXFLAGS += $(STRICT_COORD_CXXFLAGS)
${GCLDIR}/gff.o : ${GCLDIR}/gff.h ${GCLDIR}/GFaSeqGet.h ${GCLDIR}/GList.hh
${GCLDIR}/GFaSeqGet.o : ${GCLDIR}/GFaSeqGet.h ${GCLDIR}/GBgzf.h
${GCLDIR}/GBgzf.o : ${GCLDIR}/GBgzf.h
gffread: gclib-init $(OBJS) gffread.o
	${LINKER} ${LDFLAGS} -o $@ $(OBJS) gffread.o ${LIBS}
#	@echo
#	${DBG_WARN}

test tests: gffread
	@./run_tests.sh

large-tests: gffread
	@./run_large_tests.sh

strict-coords: gclib-init
	@$(MAKE) --no-print-directory clean debug STRICT_COORDS=1
	@tmp=$$(mktemp); \
	filt=$$(mktemp); \
	for src in gff_utils.cpp gffread.cpp; do \
	  echo "Checking $$src with $(STRICT_CHECK_FLAGS)"; \
	  ${CXX} -I. -isystem ${GCLDIR} -Wall -Wextra -std=c++11 -D_FILE_OFFSET_BITS=64 \
	    -D_LARGEFILE_SOURCE -D_REENTRANT -fno-strict-aliasing -fno-exceptions -fno-rtti \
	    -g -O0 -DDEBUG -D_DEBUG -DGDEBUG $(STRICT_CHECK_FLAGS) -fsyntax-only $$src >>$$tmp 2>&1 || true; \
	done; \
	rg -n "(gff_utils\\.h|gff_utils\\.cpp|gffread\\.cpp):[0-9]+:[0-9]+: (warning|error):.*(int64_t|%d)" $$tmp >$$filt || true; \
	if [ -s $$filt ]; then \
	  cat $$filt; \
	  rm -f $$tmp $$filt; \
	  echo "strict-coords: relevant coordinate/container conversion diagnostics found."; \
	  exit 1; \
	fi; \
	rm -f $$tmp $$filt; \
	echo "strict-coords: no relevant coordinate/container conversion diagnostics found."

# target for removing all object files

.PHONY : clean
clean:
	@${RM} gffread gffread.o* gffread.exe $(OBJS)
	@${RM} core.*
