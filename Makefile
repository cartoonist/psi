##  Directories:
# Source files directory.
SRCDIR       := src
# Proto files directory.
PROTODIR     := proto
# Tests directory.
TESTDIR      := test

##  Sources:
# Proto files.
PROTOS       := $(wildcard ${PROTODIR}/*.proto)
# Sources to be generated from proto files.
PROTO_SRCS   := $(subst ${PROTODIR}, ${SRCDIR}, ${PROTOS:%.proto=%.pb.cc})

## Recipes:
COMPILE.proto = protoc -I=${PROTODIR}/ --cpp_out=${SRCDIR}/

# Specifying phony targets.
.PHONY: all test clean dist-clean
# Specifying precious targets.
.PRECIOUS: ${SRCDIR}/%.pb.cc ${SRCDIR}/%.pb.h

all: ${PROTO_SRCS}
	make -C ${SRCDIR}

${SRCDIR}/%.pb.cc ${SRCDIR}/%.pb.h:: ${PROTODIR}/%.proto
	${COMPILE.proto} $<

test: all
	make -C ${TESTDIR}

clean:
	make -C ${SRCDIR} $@
	make -C ${TESTDIR} $@

dist-clean: clean
	make -C ${SRCDIR} $@
	rm -f ${SRCDIR}/*.pb.cc ${SRCDIR}/*.pb.h
	make -C ${TESTDIR} $@
