# Top level Makefile

UTIL_MODULE = util
CHECKOUT_SCRIPT = checkout.sh
DIFF_SCRIPT = diff.sh
COMPILE_SCRIPT = compile.sh
CLEAN_SCRIPT = clean.sh
TEST_SCRIPT = test.sh
CLEAN_TEST_SCRIPT = clean_test.sh
DOC_SCRIPT = doc.sh
CLEAN_DOC_SCRIPT = clean_doc.sh
EXPORT_SCRIPT = export.sh

all: compile

build: checkout compile

clean: clean_build clean_test clean_doc

checkout:
	@echo 
	@echo -------- Getting version Latest of module util --------
	@cvs co $(UTIL_MODULE)
	@sh -c 'cd ./$(UTIL_MODULE); ./$(CHECKOUT_SCRIPT)'

diff:
	@sh -c 'cd ./$(UTIL_MODULE); ./$(DIFF_SCRIPT)'

compile:
	@sh -c 'cd ./$(UTIL_MODULE); ./$(COMPILE_SCRIPT)'

debug:
	@sh -c 'cd ./$(UTIL_MODULE); ./$(COMPILE_SCRIPT) debug'

clean_build: 
	@sh -c 'cd ./$(UTIL_MODULE); ./$(CLEAN_SCRIPT)'

test: compile
	@sh -c 'cd ./$(UTIL_MODULE); ./$(TEST_SCRIPT)'

clean_test: 
	@sh -c 'cd ./$(UTIL_MODULE); ./$(CLEAN_TEST_SCRIPT)'

doc:
	@sh -c 'cd ./$(UTIL_MODULE); ./$(DOC_SCRIPT)'

clean_doc:
	@sh -c 'cd ./$(UTIL_MODULE); ./$(CLEAN_DOC_SCRIPT)'

export: clean
	@sh -c 'cd ./$(UTIL_MODULE); ./$(EXPORT_SCRIPT)'

