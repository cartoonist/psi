## Directories:
# Binary files directory.
BINDIR       := ../bin
# Release binary directory (relative to `$BINDIR`).
RELDIR       := release
# Debug binary directory (relative to `$BINDIR`).
DEBDIR       := debug
# Build type
BUILD        ?= release  # Should be one of the binary directories: `debug` or `release`
# Make the binary files directory.
$(shell mkdir -p ${BINDIR}/${BUILD} >/dev/null)
# Dependencies Makefiles' directory.
DPMDIR       := .d
# Make the dependencies Makefiles' directory.
$(shell mkdir -p ${DPMDIR} >/dev/null)

##  Build variables:
# Compiler.
CXX          ?= g++
# Last commit date.
UPDATE_DATE  := $(shell git log -n 1 --pretty='format:%cd' --date='format:%b %d %Y')
# Git version tag.
GIT_VERSION  := $(shell git describe)
# Compiler flags.
ifeq ($(BUILD),debug)
MACROS       ?= -DGREM_DEBUG=1
endif
MACROS       ?= -DGREM_DEBUG=0
MACROS       += -DGREM_UPDATE_DATE="\"${UPDATE_DATE}\"" -DGREM_GIT_VERSION="\"${GIT_VERSION}\""

CPPFLAGS      = -Wall -Werror -pedantic -fopenmp -std=c++14 ${MACROS} ${LIBS_CFLAGS}
CXXFLAGS     ?= -g -O3
# FIXME: Manage the dependencies in a better way!
VGDIR        ?= /home/ghaffari/src/vg
LIBS_CFLAGS  := -I. \
                -I${VGDIR}/include \
                -I${VGDIR}/src \
                -I${VGDIR} \
                -I${VGDIR}/deps/sdsl-lite/build/include \
                $(shell pkg-config --cflags jansson) \
                $(shell pkg-config --cflags seqan-2) \
                $(shell pkg-config --cflags protobuf)
LDFLAG       := -fopenmp \
                -lpthread \
                -lm \
                -lz \
                -llzma \
                -lbz2 \
                -lrt \
                -latomic \
                -L${VGDIR}/lib -lvg \
                -L${VGDIR}/deps/vcflib/lib -lvcflib \
                -L${VGDIR}/deps/vcflib/tabixpp/htslib -lhts \
                -L${VGDIR}/deps/sdsl-lite/build/lib -lsdsl -ldivsufsort -ldivsufsort64 \
                $(shell pkg-config --libs jansson) \
                $(shell pkg-config --libs protobuf)

##  Specifying source files to automatically generate dependencies:
# Consider all *.cc files in the src directory as 'source' file.
SOURCES      := $(wildcard *.cc)

##  Targets:
TARGET       := ${BINDIR}/${BUILD}/grem

##  Recipes:
COMPILE.cxx   = ${CXX} -c ${CXXFLAGS} ${CPPFLAGS} -o $@
LINK          = ${CXX} $^ ${LDFLAG} -o $@
MAKE_DMAKES   = ${CXX} ${CPPFLAGS} -MM -MG -MF $@

# Specifying phony targets.
.PHONY: all clean distclean

# building target.
all: ${TARGET}

${BINDIR}/release/grem: grem.o
	${LINK}

${BINDIR}/debug/grem: grem.g.o
	${LINK}

# general rule for compiling cc files.
%.g.o: %.cc
	${COMPILE.cxx} $<
%.o: %.cc
	${COMPILE.cxx} $<

clean:
	rm -f *.o *~

distclean: clean
	rm -f -r ${BINDIR}; \
	rm -f -r ${DPMDIR}

# generating sources' prerequisites automatically by remaking *.d Makefiles.
${DPMDIR}/%.d: %.cc
	${MAKE_DMAKES} -MT $*.o -MT $@ $<

# generating sources' prerequisites automatically by remaking *.g.d Makefiles (debug mode).
${DPMDIR}/%.g.d: %.cc
	${MAKE_DMAKES} -MT $*.g.o -MT $@ $<

# include *.d Makefiles.
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(BUILD),debug)
-include $(addprefix ${DPMDIR}/, ${SOURCES:%.cc=%.d})
else
-include $(addprefix ${DPMDIR}/, ${SOURCES:%.cc=%.g.d})
endif
endif
endif
