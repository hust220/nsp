##################################################################################################

PROJ="serial_debug"

SD=$(cd $(dirname ${0}); pwd) # source directory
BD=$SD/.build # building directory
NTHREADS=16 # number of threads
CXX_FLAGS="-std=c++14 -pthread -lm"
CXX_INCLUDES="-I$SD/ext -I$SD/ext/lua-5.3.3/include"
LD_LIB_PATHS="${SD}/ext/lua-5.3.3/lib"
LD_LIBS="lua dl fftw3"
TARGET=$BD/nsp

proj_serial_release() {
    PD="$BD/serial-release" # project directory

    BUILD_MODE="Serial" # Serial or Para
    BUILD_TYPE="Release" # Release or Debug

    CXX=g++

    #LD_FLAGS="-L${SD}/lib/lua-5.3.3/lib -Wl,-rpath,${SD}/lib/lua-5.3.3/lib: -rdynamic -llua -ldl"
}

proj_serial_debug() {
    PD="$BD/serial-debug" # project directory

    BUILD_MODE="Serial" # Serial or Para
    BUILD_TYPE="Debug" # Release or Debug

    CXX=g++
}

proj_para_release() {
    PD="$BD/para-release" # project directory

    BUILD_MODE="Para" # Serial or Para
    BUILD_TYPE="Release" # Release or Debug

    CXX=mpicxx
}

###################################################################################################

set_ld_lib_paths() {
    local libpaths="$@"
    local i
    for i in ${libpaths}; do
        echo -n -L${i}" "
    done
    echo -n -Wl,-rpath,
    for i in ${libpaths}; do
        echo -n ${i}:
    done
}

set_ld_libs() {
    local libs="$@"
    local i
    for i in ${libs}; do
        echo -n -l${i}" "
    done
}

set_proj() {
    proj_${1}
    if [[ $? -ne 0 ]]; then return 1; fi

    LD_FLAGS="$(set_ld_lib_paths ${LD_LIB_PATHS}) $(set_ld_libs ${LD_LIBS})"

    if [[ ${BUILD_TYPE} = "Release" ]]; then
        CXX_FLAGS="${CXX_FLAGS} -O3 -DNDEBUG"
    elif [[ ${BUILD_TYPE} = "Debug" ]]; then
        CXX_FLAGS="${CXX_FLAGS} -g -gdwarf-2"
    fi

    if [[ ${BUILD_MODE} = "Para" ]]; then
        CXX_FLAGS="${CXX_FLAGS} -DJN_PARA"
    fi

}

# path of directory
pdir() {
    echo $PD/objs${1}
}

# path of object
pobj() {
    echo $PD/objs${1}.o
}

check_dir() {
    local pdir=$(pdir $1)
    if [ ! -d ${pdir} ]; then
        mkdir -p ${pdir}
    fi
}

find_cpps() {
    local dir="";
    for dir in $@; do
        dir=$(cd ${dir}; pwd) # absolute path
        check_dir ${dir}
        for i in $(ls ${dir}); do
            local name=${dir}/${i}
            if [ -d ${name} ]; then
                find_cpps ${name}
            else
                local filename=${name%.*}
                local ext=${name##*.}
                if [[ ${ext} = "cpp" || ${ext} = "c" ]]; then
                    echo ${name}
                fi
            fi
        done
    done
}

build_lua() {
    cd $SD/ext/lua-5.3.3
    make linux
    make install
    cd $SD
}

build() {
    if [ ! -d ${PD} ]; then mkdir -p ${PD}; fi
    cd ${PD}
    rm -rf makefile

    local cpps=$(find_cpps ${SD}/src)
    local objs=$(for cpp in ${cpps}; do echo $(pobj ${cpp}); done)

    echo all: $objs >>makefile
    echo -e "\t"${CXX} ${CXX_FLAGS} ${CXX_INCLUDES} $objs -o ${TARGET} ${LD_FLAGS} >>makefile

    for cpp in ${cpps}; do
        local obj=$(pobj ${cpp})
        echo ${obj}: ${cpp} >>makefile
        echo -e "\t"${CXX} ${CXX_FLAGS} ${CXX_INCLUDES} -c ${cpp} -o ${obj} >>makefile
        echo >>makefile
    done

    make -j${NTHREADS}
    cd ${SD}
}

if [[ $# -eq 0 ]]; then
    set_proj ${PROJ}
elif [[ $# -eq 1 ]]; then
    set_proj ${1}
    if [[ $? -ne 0 ]]; then
        echo Error: unknown project \'$1\'!
        exit
    fi
fi
build_lua
build


