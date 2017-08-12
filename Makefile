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
.PHONY: all release benchmark debug test doc tags clean distclean
# Specifying precious targets.
.PRECIOUS: ${SRCDIR}/%.pb.cc ${SRCDIR}/%.pb.h

## Functions:
define echotitle
	@echo
	@tput setaf 2
	@echo "$1"
	@tput sgr0
	@echo
endef

all: release

release: ${PROTO_SRCS}
	$(call echotitle,"Building sources...")
	@make -C ${SRCDIR} BUILD=release

benchmark: ${PROTO_SRCS}
	$(call echotitle,"Building sources for benchmark...")
	@make -C ${SRCDIR} BUILD=release MACROS=-DGREM_DEBUG=1

debug: ${PROTO_SRCS}
	$(call echotitle,"Building sources for debug...")
	@make -C ${SRCDIR} BUILD=debug

${SRCDIR}/%.pb.cc ${SRCDIR}/%.pb.h:: ${PROTODIR}/%.proto
	$(call echotitle,"Compiling protocol buffers...")
	${COMPILE.proto} $<

test: benchmark
	$(call echotitle,"Building tests...")
	@make -C ${TESTDIR}
	$(call echotitle,"Running tests...")
	@find ${TESTDIR}/bin -maxdepth 1 -type f -name "tests_*" | xargs -I{} sh -c {}

test-debug: debug
	$(call echotitle,"Building tests for debug...")
	@make -C ${TESTDIR} debug
	$(call echotitle,"Running tests...")
	@find ${TESTDIR}/bin -maxdepth 1 -type f -name "tests_*_d" | xargs -I{} sh -c {}

doc:
	$(call echotitle,"Generating source code documentation...")
	doxygen

tags:
	$(call echotitle,"Updating ctags...")
	ctags ${ALL_SOURCES}

install: all
	$(call echotitle,"Installing binaries...")
	@install -v ${BINDIR}/release/grem ${PREFIX}/bin

clean:
	@make -C ${SRCDIR} $@
	@make -C ${TESTDIR} $@

distclean:
	@make -C ${SRCDIR} $@
	rm -f ${SRCDIR}/*.pb.cc ${SRCDIR}/*.pb.h
	@make -C ${TESTDIR} $@
