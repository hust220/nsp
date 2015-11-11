#
#        REGEX module makefile
#
#----------------------------------------------------------------------------
# Project specific path defintions.
#----------------------------------------------------------------------------
M_INCL_DIR = ../include
M_LIB_DIR  = ../lib

PROJ_DIR   = .

L_INCL_DIR = $(PROJ_DIR)/include
SRC_DIR    = $(PROJ_DIR)/src
OBJ_DIR    = $(PROJ_DIR)/obj
L_LIB_DIR  = $(PROJ_DIR)/lib
L_BIN_DIR  = $(PROJ_DIR)/bin
TEST_DIR    = $(PROJ_DIR)/test-regex

VPATH = $(OBJ_DIR)

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
ALL_DEP_LIBS =

# Module libraries
MOD_LIB = regex.a

# Agregate library
AGR_LIB = all.a

# Temporary library. Used to obtain the agregate library.
TMP_LIB = tmp.a

L_MOD_LIB = $(L_LIB_DIR)/$(MOD_LIB)
M_MOD_LIB = $(M_LIB_DIR)/$(MOD_LIB)
L_AGR_LIB = $(L_LIB_DIR)/$(AGR_LIB)
M_AGR_LIB = $(M_LIB_DIR)/$(AGR_LIB)


SRC_FILES = regcomp.c \
	    regexec.c \
	    regerror.c \
	    regfree.c

# Object files. Replace ".c" with ".o"
OBJ_FILES = ${SRC_FILES:.c=.o}

HEADER_FILES = regex.h

H = cclass.h \
    cname.h \
    regex2.h \
    utils.h

IH = regcomp.ih \
     engine.ih \
     regerror.ih \
     debug.ih \
     main.ih

SRC_INCLUDE_FILES = engine.c

TOPFILES = Makefile COPYRIGHT README WHATSNEW regex.3 regex.7

ALL_OBJ_FILES = *.o

.PHONY: ../etc/Makefile.platform all install export clean clean_build


all: install


install: $(M_MOD_LIB)


export:
	mkdir -p $(EXPORT_DIR)
	@cp $(TOPFILES) $(EXPORT_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(L_INCL_DIR)
	@cd $(L_INCL_DIR); $(EXPORT) $(EXPORT_LIST) $(HEADER_FILES) ../$(EXPORT_DIR)/$(L_INCL_DIR)
	@cd $(L_INCL_DIR); $(EXPORT) $(EXPORT_LIST) $(H) ../$(EXPORT_DIR)/$(L_INCL_DIR)
	@cd $(L_INCL_DIR); $(EXPORT) $(EXPORT_LIST) $(IH) ../$(EXPORT_DIR)/$(L_INCL_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(SRC_DIR)
	@cd $(SRC_DIR); $(EXPORT) $(EXPORT_LIST) $(SRC_FILES) ../$(EXPORT_DIR)/$(SRC_DIR)
	@cd $(SRC_DIR); $(EXPORT) $(EXPORT_LIST) $(SRC_INCLUDE_FILES) ../$(EXPORT_DIR)/$(SRC_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(OBJ_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(L_LIB_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(L_BIN_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(TEST_DIR)
	@cp $(TESTFILES) $(EXPORT_DIR)/$(TEST_DIR)


clean: clean_build clean_test


$(M_MOD_LIB): $(L_MOD_LIB)
#       Install header files
	@cd $(L_INCL_DIR); \
          ../$(INSTALL) $(INSTALLOPTS) $(HEADER_FILES) ../$(M_INCL_DIR)

#       Install module library
	$(INSTALL) $(INSTALLOPTS) $(L_MOD_LIB) $(M_LIB_DIR)

#       Create agregate library

	@cd $(M_LIB_DIR); ../etc/initlib.sh $(MOD_LIB)

	@cd $(L_LIB_DIR); cp ../$(M_AGR_LIB) $(TMP_LIB)
	@cd $(L_LIB_DIR); $(AR) $(AR_GETFLAGS) $(TMP_LIB)
	@cd $(L_LIB_DIR); rm -f $(TMP_LIB)

	@cd $(L_LIB_DIR); cp $(MOD_LIB) $(TMP_LIB)
	@cd $(L_LIB_DIR); $(AR) $(AR_GETFLAGS) $(TMP_LIB)
	@cd $(L_LIB_DIR); rm -f $(TMP_LIB)

	@cd $(L_LIB_DIR); $(AR) $(AR_PUTFLAGS) $(AGR_LIB) $(ALL_OBJ_FILES)
	@cd $(L_LIB_DIR); rm -f $(ALL_OBJ_FILES)

	$(INSTALL) $(INSTALLOPTS) $(L_AGR_LIB) $(M_LIB_DIR)
	@rm -f $(L_AGR_LIB)


clean_build:	
	@cd $(M_INCL_DIR); rm -f $(HEADER_FILES)
	@rm -f $(OBJ_DIR)/*.o
	@rm -rf $(OBJ_DIR)/ii_files
	@rm -f $(L_MOD_LIB)
	@rm -f $(M_MOD_LIB)
	@rm -f $(M_AGR_LIB)

clean_test:


$(L_MOD_LIB): $(OBJ_FILES)
#       Create module library
	@cd $(OBJ_DIR); $(AR) $(AR_PUTFLAGS) ../$@ $(OBJ_FILES)
	$(RANLIB) $@
	@echo $@ " is up to date."


# Rule for making object files
%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS_NONANSI) -DPOSIX_MISTAKE -c $< -o $(OBJ_DIR)/$@


# Extra, non-used in PDB build

ALLOBJS = $(OBJ_FILES) \
	$(OBJ_DIR)/split.o \
	$(OBJ_DIR)/debug.o \
	$(OBJ_DIR)/main.o

TESTFILES = $(TEST_DIR)/tests


REGSRC = $(SRC_DIR)/regcomp.c \
	$(SRC_DIR)/regerror.c \
	$(SRC_DIR)/regexec.c \
	$(SRC_DIR)/regfree.c

ALLSRC = $(REGSRC) \
	$(SRC_DIR)/engine.c \
	$(SRC_DIR)/debug.c \
	$(SRC_DIR)/main.c \
	$(SRC_DIR)/split.c

TARGET   = $(L_BIN_DIR)/re


$(OBJ_DIR)/main.o: $(H) $(HEADER_FILES) $(SRC_DIR)/main.c $(L_INCL_DIR)/main.ih
	$(CC) $(CFLAGS_NONANSI) -DPOSIX_MISTAKE -c $(SRC_DIR)/main.c -o $@

$(OBJ_DIR)/split.o: $(H) $(HEADER_FILES) $(SRC_DIR)/split.c 
	$(CC) $(CFLAGS_NONANSI) -DPOSIX_MISTAKE -c $(SRC_DIR)/split.c -o $@

debug.o: $(H) $(HEADER_FILES) $(SRC_DIR)/debug.c $(L_INCL_DIR)/debug.ih
	$(CC) $(CFLAGS_NONANSI) -DPOSIX_MISTAKE -c $(SRC_DIR)/debug.c -o $(OBJ_DIR)/$@


# tester
$(TARGET): $(ALLOBJS)
	$(CC) $(CFLAGS) -DPOSIX_MISTAKE $(LDFLAGS) $(ALLOBJS) $(LIBS) -o $@


# regression test
test:	$(TARGET) $(TEST_DIR)/tests
	$(TARGET) <$(TEST_DIR)/tests
	$(TARGET) -el <$(TEST_DIR)/tests
	$(TARGET) -er <$(TEST_DIR)/tests


# 57 variants, and other stuff, for development use -- not useful to you
ra:	$(TARGET) $(TEST_DIR)/tests
	-$(TARGET) <$(TEST_DIR)/tests
	-$(TARGET) -el <$(TEST_DIR)/tests
	-$(TARGET) -er <$(TEST_DIR)/tests


rx:	$(TARGET) $(TEST_DIR)/tests
	$(TARGET) -x <$(TEST_DIR)/tests
	$(TARGET) -x -el <$(TEST_DIR)/tests
	$(TARGET) -x -er <$(TEST_DIR)/tests


t:	$(TARGET) $(TEST_DIR)/tests
	-time $(TARGET) <$(TEST_DIR)/tests
	-time $(TARGET) -cs <$(TEST_DIR)/tests
	-time $(TARGET) -el <$(TEST_DIR)/tests
	-time $(TARGET) -cs -el <$(TEST_DIR)/tests


fullprint:
	ti README WHATSNEW notes todo | list
	ti *.h | list
	list *.c
	list regex.3 regex.7


print:
	ti README WHATSNEW notes todo | list
	ti *.h | list
	list reg*.c engine.c


