##  Directories:
# Source files directory.
SRCDIR       := src
# Proto files directory.
PROTODIR     := proto
# Tests directory.
TESTDIR      := test
# Binary files directory.
BINDIR       := ./bin
# The location where the binary files should be installed.
PREFIX       ?= ~/.local

##  Sources:
# Proto files.
PROTOS       := $(wildcard ${PROTODIR}/*.proto)
# Sources to be generated from proto files.
PROTO_SRCS   := $(subst ${PROTODIR}, ${SRCDIR}, ${PROTOS:%.proto=%.pb.cc})
# All source files (including header files).
ALL_SOURCES  := $(wildcard ${SRCDIR}/*.cc)
ALL_SOURCES  += $(wildcard ${SRCDIR}/*.h)
ALL_SOURCES  += $(wildcard ${TESTDIR}/*.cc)
ALL_SOURCES  += $(wildcard ${TESTDIR}/*.h)

## Recipes:
COMPILE.proto = protoc -I=${PROTODIR}/ --cpp_out=${SRCDIR}/

# Specifying phony targets.
.PHONY: all debug test doc tags clean distclean
# Specifying precious targets.
.PRECIOUS: ${SRCDIR}/%.pb.cc ${SRCDIR}/%.pb.h

all: ${PROTO_SRCS}
	make -C ${SRCDIR}

debug: ${PROTO_SRCS}
	make -C ${SRCDIR} MACROS=-DGREM_DEBUG=1

${SRCDIR}/%.pb.cc ${SRCDIR}/%.pb.h:: ${PROTODIR}/%.proto
	${COMPILE.proto} $<

test: all
	@echo
	@echo "Building tests..."
	@echo
	@make -C ${TESTDIR}
	@echo
	@echo "Running tests..."
	@echo
	@find ${TESTDIR}/bin -maxdepth 1 -type f -name "tests_*" | xargs -I{} sh -c {}

doc:
	doxygen

tags:
	ctags ${ALL_SOURCES}

install:
	install -v ${BINDIR}/grem ${PREFIX}/bin

clean:
	make -C ${SRCDIR} $@
	make -C ${TESTDIR} $@

distclean: clean
	make -C ${SRCDIR} $@
	rm -f ${SRCDIR}/*.pb.cc ${SRCDIR}/*.pb.h
	make -C ${TESTDIR} $@
