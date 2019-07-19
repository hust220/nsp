##### Please set the BIN_DIR variable #####
BIN_DIR   := /home/share/apps/nsp/1.8.0
###########################################

FLAGS     := -std=c++14 -pthread -lm -Isrc -MMD -lfftw3
CC        := g++

ifeq (true,${MPI})
	FLAGS        := $(FLAGS) -DUSEMPI
endif

ifeq (true,${DEBUG})
	BUILD_PREFIX := build/debug
	FLAGS        := $(FLAGS) -g -gdwarf-2
else
	BUILD_PREFIX := build/release
	FLAGS        := $(FLAGS) -DNDEBUG -O3
endif

SRC_DIR   := src
APPS_DIR  := apps

vpath %.cpp $(SRC_DIR) $(APPS_DIR)

BUILD_DIR := $(addprefix $(BUILD_PREFIX)/, $(SRC_DIR) $(APPS_DIR))
#FLAGS     := $(FLAGS) $(foreach sdir, $(SRC_DIR) $(APPS_DIR), -I$(sdir))

SRC_CPP   := $(foreach sdir,$(SRC_DIR), $(notdir $(wildcard $(sdir)/*.cpp)))
APPS_CPP  := $(foreach sdir,$(APPS_DIR), $(notdir $(wildcard $(sdir)/*.cpp)))

SRC_OBJ   := $(patsubst %.cpp, $(BUILD_PREFIX)/%.o, $(SRC_CPP))
APPS_OBJ  := $(patsubst %.cpp, $(BUILD_PREFIX)/%.o, $(APPS_CPP))

APPS      := $(patsubst %.cpp, $(notdir %), $(APPS_CPP))

#$(BIN_DIR)/$1: $(BUILD_PREFIX)/$(notdir $1).o $(SRC_OBJ)

#$(CC) $$^ -o $$@ $(FLAGS) -L/hpc/home/juw1179/programs/fftw/3.3.8/lib

define make-apps

$1: checkdirs $(BIN_DIR)/$1

$(BIN_DIR)/$1: $(APPS_DIR)/$(notdir $1).cpp $(SRC_OBJ)
	$(CC) $$^ -o $$@ $(FLAGS)

endef

all: checkdirs $(SRC_OBJ) $(APPS_OBJ)

install: checkdirs $(APPS)

$(foreach app, $(APPS), $(eval $(call make-apps,$(app))))

$(BUILD_PREFIX)/%.o: %.cpp
	$(CC) $(FLAGS) -nostartfiles -c $< -o $@

.PHONY: all install checkdirs clean $(APPS)

checkdirs: $(BUILD_DIR) $(BIN_DIR) $(BUILD_PREFIX)

$(BUILD_DIR) $(BIN_DIR) $(BUILD_PREFIX):
	@mkdir -p $@

clean:
	@rm -rf bin build

-include $(BUILD_PREFIX)/*.d

