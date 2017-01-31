## Installation requirements

*   g++ version is greater than 4.8
*   cmake version is greater than 2.8.7
*   git

## Compile and install

1.  upgrade g++
    
    If the g++ version is less than 4.8, you need to upgrade the g++ firstly.Here take gcc-4.9.3 as an example:
        
    *   root user
 
            wget http://mirror.hust.edu.cn/gnu/gcc/gcc-4.9.3/gcc-4.9.3.tar.gz
            tar xvzf gcc-4.9.3.tar.gz
            cd gcc-4.9.3
            ./contrib/download_prerequisites
            mkdir build
            cd build
            ../configure --enable-checking=release --enable-languages=c,c++ --disable-multilib
            make -j4
            sudo make install

    *   ordinary user

        in the step of the configure plus `--prefix=<PATH/TO/INSTALL/GCC>`

2.  download nsp

        git clone https://github.com/hust220/nsp.git

3.  compile and install nsp

    *   root user
 
            cd nsp
            python jnpack.py install
            mkdir build
            cd build
            cmake ..
            make install

    *   ordinary user

            cd nsp
            python jnpack.py install
            mkdir build
            cd build
            cmake -DP=<PATH> ..
            make install

4.  download and uncompress templates library

        wget http://biophy.hust.edu.cn/download/nsp-lib.tar.gz
        tar xvzf nsp-lib.tar.gz

5.  Set the enviroment
    
        export NSP=<PATH/TO/NSP/LIBRARY>

## Our laboratory

Members of our lab can enter the following commands directly on the server.

    source $HOME/../wangjian/wangjian.sh
 
