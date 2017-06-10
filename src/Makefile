## Directories:
# External libraries directory.
EXTDIR       := ../ext
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
# Compiler flags.
ifeq ($(BUILD),debug)
MACROS       ?= -DGREM_DEBUG=1
CXXFLAGS     := -g
endif
MACROS       ?= -DGREM_DEBUG=0

CPPFLAGS      = -Wall -Werror -pedantic -fopenmp -std=c++14 ${MACROS} ${LIBS_CFLAGS}
CXXFLAGS     ?= -O3
LIBS_CFLAGS  := -I ${EXTDIR}/ \
	        $(shell pkg-config --cflags seqan-2) \
	        $(shell pkg-config --cflags protobuf)
LDFLAG       := -fopenmp \
	        -lz \
		-lbz2 \
	        -lrt \
	        $(shell pkg-config --libs protobuf)

##  Specifying source files to automatically generate dependencies:
# Consider all *.cc files in the src directory as 'source' file.
SOURCES      := $(wildcard *.cc)

##  Targets:
TARGET       := ${BINDIR}/${BUILD}/grem

##  Recipes:
COMPILE.cxx   = ${CXX} -c ${CXXFLAGS} ${CPPFLAGS} -o $@
LINK          = ${CXX} ${LDFLAG} -o $@
MAKE_DMAKES   = ${CXX} ${CPPFLAGS} -MM -MG -MF $@

# Specifying phony targets.
.PHONY: all clean distclean

# building target.
all: ${TARGET}

${BINDIR}/release/grem: grem.o vargraph.o vg.pb.o
	${LINK} $^

${BINDIR}/debug/grem: grem.g.o vargraph.g.o vg.pb.g.o
	${LINK} $^

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

# include *.d Makefiles.
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
-include $(addprefix ${DPMDIR}/, ${SOURCES:%.cc=%.d})
endif
endif
