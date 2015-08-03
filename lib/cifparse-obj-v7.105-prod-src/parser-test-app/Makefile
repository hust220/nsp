#
#        CIFPARSE-TEST module makefile
#
#----------------------------------------------------------------------------
# Project specific path defintions.
#----------------------------------------------------------------------------
M_INCL_DIR = ../include
M_LIB_DIR  = ../lib
M_BIN_DIR  = ../bin

PROJ_DIR   = .

L_INCL_DIR = $(PROJ_DIR)/include
SRC_DIR    = $(PROJ_DIR)/src
OBJ_DIR    = $(PROJ_DIR)/obj
L_BIN_DIR  = $(PROJ_DIR)/bin
TEST_DIR   = $(PROJ_DIR)/test

VPATH = $(L_BIN_DIR)

#----------------------------------------------------------------------------
# LINCLUDES and LDEFINES are appended to CFLAGS and C++FLAGS
#----------------------------------------------------------------------------
LDEFINES  =
LINCLUDES = -I$(L_INCL_DIR) -I$(M_INCL_DIR) 

#----------------------------------------------------------------------------
# Include the appropriate compiler/platform definitions ...
#----------------------------------------------------------------------------
include ../etc/Makefile.platform

# Dependent libraries
CIF_FILE_LIB      = $(M_LIB_DIR)/cif-file.a
CIFPARSE_LIB      = $(M_LIB_DIR)/cifparse-obj.a
TABLES_LIB        = $(M_LIB_DIR)/tables.a
COMMON_LIB        = $(M_LIB_DIR)/common.a
REGEX_LIB         = $(M_LIB_DIR)/regex.a

ALL_DEP_LIBS = $(CIF_FILE_LIB) $(CIFPARSE_LIB) $(TABLES_LIB) \
  $(COMMON_LIB) $(REGEX_LIB)

# Base file names. Must have ".ext" at the end of the file.
BASE_MAIN_FILES = FOtest1.ext \
                  FOtest2.ext \
                  FOtest3.ext \
                  FOtest4.ext \
                  FOtest5.ext \
                  FOtest6.ext \
                  FOtest7.ext \
                  FOtest8.ext \
                  FOtest9.ext \
                  FOtest10.ext \
                  FOtest11.ext \
                  FOtest12.ext \
                  FOtest13.ext \
                  FOtest14.ext \
                  DataChecking.ext \
                  checkdictionary.ext \
                  selective-reading.ext \
                  BigCifTest.ext \
                  CifReader.ext \
                  SdbReader.ext \
                  CifDiff.ext \
                  CifInfo.ext \
                  DictInfo.ext \
                  Search.ext

# Main source files. Replace ".ext" with ".C"
SRC_MAIN_FILES = ${BASE_MAIN_FILES:.ext=.C}

# All source files.
SRC_FILES = $(SRC_MAIN_FILES)

# All object files. Replace ".C" with ".o"
OBJ_FILES = ${SRC_FILES:.C=.o}

# Executables. Remove ".ext"
TARGETS = ${BASE_MAIN_FILES:.ext=}

# Test related files
TEST_FILES = $(TEST_DIR)/ADH041.cif \
             $(TEST_DIR)/ddl-mm \
             $(TEST_DIR)/dictionary.dic \
             $(TEST_DIR)/ndb_ddl \
             $(TEST_DIR)/ndb_dictionary.dic \
             $(TEST_DIR)/infile01.db \
             $(TEST_DIR)/infile02.db \
             $(TEST_DIR)/Test1.cif \
             $(TEST_DIR)/Test2.cif \
             $(TEST_DIR)/Test3.cif \
             $(TEST_DIR)/Test14.cif \
             $(TEST_DIR)/view1.cif \
             $(TEST_DIR)/view2.cif \
             $(TEST_DIR)/test.sh

.PHONY: ../etc/Makefile.platform all install test export clean clean_build clean_test 

.PRECIOUS: $(OBJ_DIR)/%.o

# All
all: install


# Installation
install: $(TARGETS)


# Test
test: all
	@sh -c 'cd $(TEST_DIR); ./test.sh'


export:
	mkdir -p $(EXPORT_DIR)
	@cp Makefile $(EXPORT_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(SRC_DIR)
	@cd $(SRC_DIR); $(EXPORT) $(EXPORT_LIST) $(SRC_FILES) ../$(EXPORT_DIR)/$(SRC_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(OBJ_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(L_BIN_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(TEST_DIR)
	@cp $(TEST_FILES) $(EXPORT_DIR)/$(TEST_DIR)


clean: clean_build clean_test


# Rule for making executables
%: $(OBJ_DIR)/%.o $(ALL_DEP_LIBS)
	$(CCC) $(LDFLAGS) $< $(ALL_DEP_LIBS) $(MALLOCLIB) -lm -o $(L_BIN_DIR)/$@
	@cp -f $(L_BIN_DIR)/$@ $(M_BIN_DIR)/$@


# Rule for build cleaning
clean_build:
	@rm -f $(OBJ_DIR)/*.o
	@rm -rf $(OBJ_DIR)/ii_files
	@cd $(L_BIN_DIR); rm -f $(TARGETS)
	@cd $(M_BIN_DIR); rm -f $(TARGETS)


# Rule for test results cleaning
clean_test:
	@rm -f $(TEST_DIR)/*-parser.log 
	@rm -f $(TEST_DIR)/*-diag.log
	@rm -f $(TEST_DIR)/test*.sdb $(TEST_DIR)/*.ocif $(TEST_DIR)/*.out
	@rm -f $(TEST_DIR)/Test15.dic $(TEST_DIR)/dictionary.sdb
	@rm -f $(TEST_DIR)/view13.dic
	@rm -f $(TEST_DIR)/exectime.txt


# Rule for making object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.C
	$(CCC) $(C++FLAGS) -c $< -o $@

