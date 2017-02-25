## <b id="requirements">Installation requirements</b>

*   g++ version is greater than 4.8
*   cmake version is greater than 2.8.7
*   git

## <b id="install">Compile and install</b>

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

            wget http://mirror.hust.edu.cn/gnu/gcc/gcc-4.9.3/gcc-4.9.3.tar.gz
            tar xvzf gcc-4.9.3.tar.gz
            cd gcc-4.9.3
            ./contrib/download_prerequisites
            mkdir build
            cd build
            ../configure --enable-checking=release --enable-languages=c,c++ --disable-multilib --prefix=<PATH/TO/INSTALL/GCC>
            make -j4
            sudo make install

2.  download nsp

    Browse [this page](https://github.com/hust220/nsp/releases) to get the latest version of NSP

    Example:

        wget https://codeload.github.com/hust220/nsp/tar.gz/1.2
        tar xvzf nsp-1.2.tar.gz

3.  compile and install nsp

            cd nsp-1.2/nsp-1.2
            python jnpack.py update
            python jnpack.py install

    The binary target named 'nsp' would be located in the '.jnp/bin' directory.

4.  download and uncompress templates library

        wget http://biophy.hust.edu.cn/download/nsp-lib.tar.gz
        tar xvzf nsp-lib.tar.gz

5.  Set the enviroment
    
        export NSP=<PATH/TO/NSP/LIBRARY>

## Our laboratory

Members of our lab can enter the following commands directly on the server.

    source $HOME/../wangjian/wangjian.sh
 
