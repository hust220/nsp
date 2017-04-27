##################################################################################################

PROJ="serial_release"

PD=$(cd $(dirname ${0}); pwd) # project directory
TARGET=nsp
NTHREADS=16 # number of threads
CXX_FLAGS="-std=c++14 -pthread -lm"
CXX_INCLUDES="-I$PD -I$PD/lib -I$PD/lib/lua-5.3.3/include"
LD_LIB_PATHS="${PD}/lib/lua-5.3.3/lib"
LD_LIBS="lua dl"

proj_serial_release() {
    BD="$PD/.install/serial-release" # building directory

    BUILD_MODE="Serial" # Serial or Para
    BUILD_TYPE="Release" # Release or Debug

    CXX=g++

    #LD_FLAGS="-L${PD}/lib/lua-5.3.3/lib -Wl,-rpath,${PD}/lib/lua-5.3.3/lib: -rdynamic -llua -ldl"
}

proj_serial_debug() {
    BD="$PD/.install/serial-debug" # building directory

    BUILD_MODE="Serial" # Serial or Para
    BUILD_TYPE="Debug" # Release or Debug

    CXX=g++
}

proj_para_release() {
    BD="$PD/.install/para-release" # building directory

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
    echo $BD/objs${1}
}

# path of object
pobj() {
    echo $BD/objs${1}.o
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
    cd $PD/lib/lua-5.3.3
    make linux
    make install
    cd $PD
}

build() {
    if [ ! -d ${BD} ]; then mkdir -p ${BD}; fi
    cd ${BD}
    rm -rf makefile

    local cpps=$(find_cpps ${PD}/lib/jian ${PD}/nsp ${PD}/src)
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
    cd ${PD}
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

