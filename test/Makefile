## Directories:
# External libraries directory.
EXTDIR       := ../ext
# Package source files directory.
SRCDIR       := ../src
# Binary files directory.
BINDIR       := ./bin
# Make the binary files directory.
$(shell mkdir -p ${BINDIR} >/dev/null)
# Dependencies Makefiles' directory.
DPMDIR       := .d
# Make the dependencies Makefiles' directory.
$(shell mkdir -p ${DPMDIR} >/dev/null)

#Makefile absolute path.
TEST_PATH := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

##  Build variables:
# Compiler.
CXX          ?= g++
# Compiler flags.
MACROS       ?= -DGREM_DEBUG=1 -DTESTDIR=\"${TEST_PATH}\"
CPPFLAGS      = -Wall -Werror -pedantic -fopenmp -std=c++14 ${MACROS} ${LIBS_CFLAGS}
CXXFLAGS     ?= -O3
LIBS_CFLAGS  := -I ${SRCDIR}/ \
	        -I ${EXTDIR}/ \
	        $(shell pkg-config --cflags seqan-2) \
	        $(shell pkg-config --cflags xg) \
	        $(shell pkg-config --cflags protobuf)
LDFLAG       := -fopenmp \
	        -lz \
		-lbz2 \
	        -lrt \
	        $(shell pkg-config --libs xg) \
	        $(shell pkg-config --libs protobuf)

##  Specifying source files to automatically generate dependencies:
# Consider all *.cc files in the src directory as 'source' file.
SOURCES      := $(wildcard *.cc)

##  Targets:
TARGETS       := ${BINDIR}/tests_utils ${BINDIR}/tests_indexiter ${BINDIR}/tests_traverser ${BINDIR}/tests_vargraph ${BINDIR}/tests_bfs ${BINDIR}/test_linear

##  Recipes:
COMPILE.cxx   = ${CXX} -c ${CXXFLAGS} ${CPPFLAGS} -o $@
LINK          = ${CXX} $^ ${LDFLAG} -o $@
MAKE_DMAKES   = ${CXX} ${CPPFLAGS} -MM -MG -MF $@

# Specifying phony targets.
.PHONY: all clean distclean

# building target.
all: ${TARGETS}

debug: CXXFLAGS := -g
debug: $(addsuffix _d, ${TARGETS})

${BINDIR}/test_linear: linear.o
	${LINK}

${BINDIR}/test_linear_d: linear.g.o
	${LINK}

${BINDIR}/tests_bfs: tests_bfs.o tests_main.o ${SRCDIR}/vargraph.o
	${LINK}

${BINDIR}/tests_bfs_d: tests_bfs.g.o tests_main.g.o ${SRCDIR}/vargraph.g.o
	${LINK}

${BINDIR}/tests_vargraph: tests_vargraph.o tests_main.o ${SRCDIR}/vargraph.o
	${LINK}

${BINDIR}/tests_vargraph_d: tests_vargraph.g.o tests_main.g.o ${SRCDIR}/vargraph.g.o
	${LINK}

${BINDIR}/tests_traverser: tests_traverser.o tests_main.o ${SRCDIR}/vargraph.o
	${LINK}

${BINDIR}/tests_traverser_d: tests_traverser.g.o tests_main.g.o ${SRCDIR}/vargraph.g.o
	${LINK}

${BINDIR}/tests_indexiter: tests_indexiter.o tests_main.o
	${LINK}

${BINDIR}/tests_indexiter_d: tests_indexiter.g.o tests_main.g.o
	${LINK}

${BINDIR}/tests_utils: tests_utils.o tests_main.o
	${LINK}

${BINDIR}/tests_utils_d: tests_utils.g.o tests_main.g.o
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

# include *.d Makefiles.
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
-include $(addprefix ${DPMDIR}/, ${SOURCES:%.cc=%.d})
endif
endif
