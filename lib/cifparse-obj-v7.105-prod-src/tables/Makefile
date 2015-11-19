#
#        CIF-TABLE-OBJ module makefile
#
#----------------------------------------------------------------------------
# Project specific path defintions.
#----------------------------------------------------------------------------
M_INCL_DIR  = ../include
M_LIB_DIR   = ../lib

PROJ_DIR    = .

L_INCL_DIR  = $(PROJ_DIR)/include
SRC_DIR     = $(PROJ_DIR)/src
OBJ_DIR     = $(PROJ_DIR)/obj
L_LIB_DIR   = $(PROJ_DIR)/lib

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
COMMON_LIB = $(M_LIB_DIR)/common.a

ALL_DEP_LIBS = $(COMMON_LIB)

# Module libraries
MOD_LIB = tables.a 

# Agregate library
AGR_LIB = all.a

# Temporary library. Used to obtain the agregate library.
TMP_LIB = tmp.a

L_MOD_LIB = $(L_LIB_DIR)/$(MOD_LIB)
M_MOD_LIB = $(M_LIB_DIR)/$(MOD_LIB)
L_AGR_LIB = $(L_LIB_DIR)/$(AGR_LIB)
M_AGR_LIB = $(M_LIB_DIR)/$(AGR_LIB)


# Base file names. Must have ".ext" at the end of the file.
BASE_FILES = TTable.ext \
             ITTable.ext \
             ISTable.ext \
             TableFile.ext

# Source files. Replace ".ext" with ".C"
SRC_FILES = ${BASE_FILES:.ext=.C}

# Object files. Replace ".ext" with ".o"
OBJ_FILES = ${BASE_FILES:.ext=.o}

# Base header files. Replace ".ext" with ".h"
BASE_HEADER_FILES = ${BASE_FILES:.ext=.h}

EXTRA_HEADER_FILES = TableError.h

HEADER_FILES = $(BASE_HEADER_FILES) $(EXTRA_HEADER_FILES)

ALL_OBJ_FILES = *.o

.PHONY: ../etc/Makefile.platform all install export clean clean_build


all: install


install: $(M_MOD_LIB)


export:
	mkdir -p $(EXPORT_DIR)
	@cp Makefile $(EXPORT_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(L_INCL_DIR)
	@cd $(L_INCL_DIR); $(EXPORT) $(EXPORT_LIST) $(HEADER_FILES) ../$(EXPORT_DIR)/$(L_INCL_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(SRC_DIR)
	@cd $(SRC_DIR); $(EXPORT) $(EXPORT_LIST) $(SRC_FILES) ../$(EXPORT_DIR)/$(SRC_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(OBJ_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(L_LIB_DIR)


clean: clean_build


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


$(L_MOD_LIB): $(OBJ_FILES)
#       Create module library
	@cd $(OBJ_DIR); $(AR) $(AR_PUTFLAGS) ../$@ $(OBJ_FILES)
	$(RANLIB) $@
	@echo $@ " is up to date."


# Rule for making object files
%.o: $(SRC_DIR)/%.C
	$(CCC) $(C++FLAGS) -c $< -o $(OBJ_DIR)/$@

