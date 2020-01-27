##  Directories:
# Source files directory.
SRCDIR       := src
# Tests directory.
TESTDIR      := test
# Binary files directory.
BINDIR       := ./bin
# The location where the binary files should be installed.
PREFIX       ?= ~/.local

##  Sources:
SOURCES  := $(wildcard ${SRCDIR}/*.cc)
SOURCES  += $(wildcard ${SRCDIR}/*.h)
SOURCES  += $(wildcard ${TESTDIR}/*.cc)
SOURCES  += $(wildcard ${TESTDIR}/*.h)

# Specifying phony targets.
.PHONY: all release debug test test-debug doc tags install install-debug clean distclean

## Functions:
define echotitle
	@echo
	@tput setaf 2
	@echo "$1"
	@tput sgr0
	@echo
endef

define make_test
	$(call echotitle,"Building tests...")
	@make -C ${TESTDIR}
endef

define run_test
	$(call echotitle,"Running tests...")
	@find ${TESTDIR}/bin -maxdepth 1 -type f -name "tests_*" | xargs -I{} sh -c {}
endef

all: release

release:
	$(call echotitle,"Building sources...")
	@make -C ${SRCDIR} BUILD=release

debug:
	$(call echotitle,"Building sources for debug...")
	@make -C ${SRCDIR} BUILD=debug CXXFLAGS=-g

test: release
	$(call make_test)
	$(call run_test)

test-debug: debug
	$(call make_test)
	$(call run_test)

doc:
	$(call echotitle,"Generating source code documentation...")
	doxygen

tags:
	$(call echotitle,"Updating ctags...")
	ctags ${SOURCES}

install: all
	$(call echotitle,"Installing binaries...")
	@install -v ${BINDIR}/release/grem ${PREFIX}/bin

install-debug: debug
	$(call echotitle,"Installing debug binaries...")
	@install -v ${BINDIR}/debug/grem ${PREFIX}/bin

clean:
	@make -C ${SRCDIR} $@
	@make -C ${TESTDIR} $@

distclean:
	@make -C ${SRCDIR} $@
	@make -C ${TESTDIR} $@
