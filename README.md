核酸分子二级结构及三级结构预测

Nucleic acid 2D and 3D Structure Prediction (NSP)

项目：<http://git.oschina.net/hust220/nsp>

1. <a href="#安装">安装</a>
 - [安装要求](#安装要求)
 - [编译安装](#编译安装)
2. [使用](#使用)
 - [二级结构预测](#二级结构预测)
 - [三级结构预测](#三级结构预测)
 - [聚类](#聚类)
 - [打分](#打分)
 - [RMSD](# RMSD)
 - [获取序列](#获取序列)

<h2 id="安装">安装</h2>
####安装要求
- g++版本大于4.8
- cmake版本大于2.8.7
- 安装有boost库

####编译安装
1. 如果g++的版本小于4.8，需要先升级g++

首先下载g++的大于4.8的版本，这里以gcc-4.9.3为例

 - root用户

```
wget http://mirror.hust.edu.cn/gnu/gcc/gcc-4.9.3/gcc-4.9.3.tar.gz
tar xvzf gcc-4.9.3.tar.gz
cd gcc-4.9.3
./contrib/download_prerequisites
mkdir build
cd build
../configure --enable-checking=release --enable-languages=c,c++ --disable-multilib
make -j4
sudo make install
```

 - 普通用户

在configure这一步的时候加上--prefix=<PATH/TO/INSTALL/GCC>

2. 下载nsp

`git clone https://git.oschina.net/hust220/nsp.git`

3. 编译安装

 - root用户：

```
    cd nsp
    mkdir build
    cd build
    cmake ..
    make install
```

 - 普通用户：

```
    cd nsp
    mkdir build
    cd build
    cmake -D CMAKE_INSTALL_PREFIX=<PATH> ..
    make install
```

####本实验室
本实验室可以直接在集群上里输入如下命令：

`source $HOME/../wangjian/wangjian.sh`

<h2 id='用法'>用法</h2>
使用之前需要设置环境变量NSP为模板库所在的文件夹
`export NSP=<PATH/OF/TEMPLATES/LIBRARY>`

<h4 id='二级结构预测'>二级结构预测</h4>
######用自由能最小方法预测
`nsp ss_pred -seq <SEQUENCE>`

######结合自由能最小方法以及DCA预测的DI值进行二级结构预测
`nsp ss_dca -seq <SEQUENCE> -di <DI_FILE> [-k <K>]`

k值是用来设置读取前k*L个DI值，如果k=1，就代表读取前L个，如果k=0.5，就代表读取前L/2个。

######计算MCC
`nsp mcc -nat "<SECONDARY_STRUCTURE_OF_NATIVE>" -pred "<SECONDARY_STRUCTURE_OF_PREDICTION>"`

######计算STY
`nsp sty -nat "<SECONDARY_STRUCTURE_OF_NATIVE>" -pred "<SECONDARY_STRUCTURE_OF_PREDICTION>"`

######计算PPV
`nsp ppv -nat "<SECONDARY_STRUCTURE_OF_NATIVE>" -pred "<SECONDARY_STRUCTURE_OF_PREDICTION>"`

<h4 id='三级结构预测'>三级结构预测</h4>
######组装
1. 组装
`nsp assemble -name <JOB_NAME> -seq <SEQUENCE> -ss "<SECONDARY_STRUCTURE>"`

2. 组装+采样
`nsp assemble -name <JOB_NAME> -seq <SEQUENCE> -ss "<SECONDARY_STRUCTURE>" -sample -num <NUMBER_OF_PREDICTIONS>`

######优化
`nsp mcpsb -name <JOB_NAME> -seq <SEQUENCE> -ss "<SECONDARY_STRUCTURE>" -pdb <INITIAL_PDB_FILE> 
-out <FINAL_STRUCTURE> [-traj <TRAJECTORY_FILE>] [-<constraints|c> <CONSTRAINTS_FILE>] [-seed <SEED>]`

用-name设置名字，-seq设置序列，用-ss设置二级结构，用-seed设置种子

用-out设置用于存放优化后的结构的文件

用-traj设置轨道文件

用-pdb设置起始结构，起始结构可以就用组装后得到的结构，也可以从组装加采样得到的结构中挑选一个

用-c或者-constraints加上约束

-seed可以省略掉，这样默认的种子是11

-constraints或者-c可以省略掉，表示不添加约束信息

CONSTRAINTS_FILE文件里面需要包含约束信息，例如可以加进DCA分析的信息:

    8 23 10
    9 22 10
    10 21 10

这里前面两列代表碱基的序号，例如8代表第8个碱基，23代表第23个碱基，最后一列是碱基之间的最小距离。

因此这里8 23 10就代笔第8个碱基和第23个碱基之间的最小的距离是10Å，9 22 10就代表第9个碱基和第22个碱基之间的最小距离是10Å。

<h4 id='聚类'>聚类</h4>
`nsp cluster -list <LIST_FILE> -k <NUMBER_OF_CLUSTERS>`

使用-list来设置需要对哪些结构进行聚类，用-k来设置聚类的数目。

经过一段时间的运行之后，会在屏幕上打印出聚类的结构。

LIST_FILE文件包含了要聚类的结构的名字:

    test.sample.1.pdb
    test.sample.2.pdb
    test.sample.3.pdb
    test.sample.4.pdb
    test.sample.5.pdb

<h4 id='打分'>RNA三级结构打分</h4>
1. 对单个结构打分

`3dRNAscore -s <PDB_FILE>`

2. 对多个结构打分

`3dRNAscore -s:l <LIST_FILE>`

LIST_FILE文件包含了要聚类的结构的名字:

    test.sample.1.pdb
    test.sample.2.pdb
    test.sample.3.pdb
    test.sample.4.pdb
    test.sample.5.pdb

<h4 id='RMSD'>计算RMSD</h4>
`nsp rmsd -pdb <PDB_FILE_1> <PDB_FILE_2>`

<h4 id='获取序列'>获取分子的序列</h4>
`nsp seq -pdb <PDB_FILE>`

####获取分子的碱基的个数
`nsp len -pdb <PDB_FILE>`

####截取分子中的指定的残基
`nsp sub -pdb <PDB_FILE> -num <FRAG1> <FRAG2> <FRAG3> <FRAG4> ...`

每个FRAG是指一个碱基段，格式为单个残基号`N`或者指定起点和终点的片段`BEGIN-END`。

例如`1`代表第一个残基，`4-11`代表第4到第11个碱基组成的片段

####去掉分子中多余的行，只留下ATOM行
`nsp rna -pdb <PDB_FILE>`