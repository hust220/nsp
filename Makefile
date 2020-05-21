##### Please set the BIN_DIR variable #####
#BIN_DIR   := /nas/longleaf/home/jianopt/work/programs/medusa/bin
PREFIX    := $(shell cat .prefix)
###########################################

BIN_DIR   := $(PREFIX)/bin
LIB_DIR   := $(PREFIX)/lib
INC_DIR   := $(PREFIX)/include

FLAGS     := -std=c++14 -pthread -lm -Isrc -MMD -lfftw3
TAIL_FLAGS := -lz
CC        := g++

ifeq (true,${MPI})
	FLAGS        := $(FLAGS) -DUSEMPI
endif

ifeq (true,${INFO})
	FLAGS        := $(FLAGS) -DJN_SHOW_INFO
endif

BUILD_TOP_PREFIX := build

ifeq (true,${DEBUG})
	BUILD_PREFIX := $(BUILD_TOP_PREFIX)/debug
	FLAGS        := $(FLAGS) -g -gdwarf-2
else
	BUILD_PREFIX := $(BUILD_TOP_PREFIX)/release
	FLAGS        := $(FLAGS) -DNDEBUG -O3
endif

SRC_DIR     := src
APPS_DIR    := apps

HEAD_DIR     := $(shell find $(SRC_DIR)/jnc -type d)
$(info $(HEAD_DIR))
FLAGS := $(FLAGS) $(addprefix -I, $(HEAD_DIR))
#vpath %.h   $(SRC_DIR) $(APPS_DIR)
#vpath %.hpp $(SRC_DIR) $(APPS_DIR)
#vpath %.c   $(SRC_DIR) $(APPS_DIR)
#vpath %.cpp $(SRC_DIR) $(APPS_DIR)

BUILD_DIRS  := $(addprefix $(BUILD_PREFIX)/, $(shell find $(SRC_DIR) $(APPS_DIR) -type d))

SRC_CPP     := $(shell find $(SRC_DIR)  -name "*.cpp" -or -name "*.c" -or -name "*.cc" -or -name "*.cxx")
APPS_CPP    := $(shell find $(APPS_DIR) -name "*.cpp" -or -name "*.c" -or -name "*.cc" -or -name "*.cxx")

HEAD_SRC    := $(shell find $(SRC_DIR) -type f ! -name "*.cpp" ! -name "*.c" ! -name "*.cc" ! -name "*.cxx")
HEAD_TARGET := $(patsubst $(SRC_DIR)/%, $(INC_DIR)/%, $(HEAD_SRC))
#$(info $$HEAD_TARGET is [${HEAD_TARGET}])

SRC_OBJ     := $(patsubst %, $(BUILD_PREFIX)/%.o, $(SRC_CPP))
APPS_OBJ    := $(patsubst %, $(BUILD_PREFIX)/%.o, $(APPS_CPP))
DEPS        := $(patsubst %, $(BUILD_PREFIX)/%.d, $(APPS_CPP) $(SRC_CPP))

#APPS        := $(patsubst %, $(basename $(notdir %)), $(APPS_CPP))
#APPS      := $(notdir $(APPS))
APPS      := $(foreach cpp, $(APPS_CPP), $(basename $(notdir $(cpp))))

all: checkdirs $(SRC_OBJ) $(APPS_OBJ)

install: checkdirs _lib _include $(APPS)

_lib: checkdirs $(LIB_DIR)/liball.a

_include: checkdirs $(HEAD_TARGET)
#_include: checkdirs

#$(info $(SRC_OBJ))

$(LIB_DIR)/liball.a: $(SRC_OBJ)
	ar rcs $@ $(SRC_OBJ)
#	ar ru $@ $(SRC_OBJ)
#	ranlib $@

define make-apps

app := $(notdir $(basename $(basename $(1))))

$$(app): checkdirs $(BIN_DIR)/$$(app)

$(BIN_DIR)/$$(app): $1 $(SRC_OBJ)
	$(CC) $(FLAGS) $$^ -o $$@ $(TAIL_FLAGS)

endef

$(foreach obj, $(APPS_OBJ), $(eval $(call make-apps,$(obj))))

$(BUILD_PREFIX)/%.o: %
	$(CC) $(FLAGS) -c $< -o $@ $(TAIL_FLAGS)

$(INC_DIR)/%: $(SRC_DIR)/%
	mkdir -p $(dir $@)
	cp $< $@

checkdirs: $(BIN_DIR) $(LIB_DIR) $(INC_DIR) $(BUILD_PREFIX) $(BUILD_DIRS)

$(BIN_DIR) $(LIB_DIR) $(INC_DIR) $(BUILD_PREFIX) $(BUILD_DIRS):
	@mkdir -p $@

clean:
	@rm -rf $(LIB_DIR) $(BIN_DIR) $(INC_DIR) $(BUILD_TOP_PREFIX)

.PHONY: all install _lib _include checkdirs clean $(APPS)

-include $(DEPS)
