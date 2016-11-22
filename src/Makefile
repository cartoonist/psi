## Directories:
# External libraries directory.
EXTDIR       := ../ext
# Binary files directory.
BINDIR       := ../bin
# Make the binary files directory.
$(shell mkdir -p ${BINDIR} >/dev/null)
# Dependencies Makefiles' directory.
DPMDIR       := .d
# Make the dependencies Makefiles' directory.
$(shell mkdir -p ${DPMDIR} >/dev/null)

##  Build variables:
# Compiler.
CXX          ?= g++
# Compiler flags.
CPPFLAGS      = -Wall -Werror -pedantic -fopenmp -std=c++14 ${LIBS_CFLAGS}
CXXFLAGS      = -O3
LIBS_CFLAGS  := -I ${EXTDIR}/ \
	        $(shell pkg-config --cflags seqan-2) \
	        $(shell pkg-config --cflags protobuf)
LDFLAG       := -lrt \
	        $(shell pkg-config --libs protobuf)

##  Specifying source files to automatically generate dependencies:
# Consider all *.cc files in the src directory as 'source' file.
SOURCES      := $(wildcard *.cc)

##  Targets:
TARGET       := ${BINDIR}/GREM

##  Recipes:
COMPILE.cxx   = ${CXX} -c ${CXXFLAGS} ${CPPFLAGS} -o $@
LINK          = ${CXX} ${LDFLAG} -o $@
MAKE_DMAKES   = ${CXX} ${CPPFLAGS} -MM -MF $@

# Specifying phony targets.
.PHONY: all clean dist-clean

# building target.
all: ${TARGET}

${BINDIR}/GREM: grem.o vargraph.o traverser.o vg.pb.o
	${LINK} $^

# general rule for compiling cc files.
%.o: %.cc
	${COMPILE.cxx} $<

clean:
	rm -f *.o *~

dist-clean: clean
	rm -f -r ${BINDIR}; \
	rm -f -r ${DPMDIR}

# generating sources' prerequisites automatically by remaking *.d Makefiles.
${DPMDIR}/%.d: %.cc
	${MAKE_DMAKES} -MT $*.o -MT $@ $<

# include *.d Makefiles.
-include $(addprefix ${DPMDIR}/, ${SOURCES:%.cc=%.d})
