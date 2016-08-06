核酸分子二级结构及三级结构预测
Nucleic acid 2D and 3D Structure Prediction (NSP)

#安装
配置3dRNA预测的环境，在集群上里输入如下命令：
`source $HOME/../wangjian/wangjian.sh`

#用法
##二级结构预测
###用自由能最小方法预测
nsp ss_pred -seq \<SEQUENCE\>

###结合自由能最小方法以及DCA预测的DI值进行二级结构预测
nsp ss_dca -seq \<SEQUENCE\> -di \<DI_FILE\> [-k \<K\>]
k值是用来设置读取前k*L个DI值，如果k=1，就代表读取前L个，如果k=0.5，就代表读取前L/2个。

##三级结构预测
###组装
1. 组装
    nsp assemble -name \<JOB_NAME\> -seq \<SEQUENCE\> -ss "\<SECONDARY_STRUCTURE\>"
2. 组装+采样
    nsp assemble -name \<JOB_NAME\> -seq \<SEQUENCE\> -ss "\<SECONDARY_STRUCTURE\>" -sample -num \<NUMBER_OF_PREDICTIONS\>

###优化
    nsp mc3p -name \<JOB_NAME\> -seq \<SEQUENCE\> -ss "\<SECONDARY_STRUCTURE\>" -pdb \<INITIAL_PDB_FILE\> [-\<constraints|c\> \<CONSTRAINTS_FILE\>] [-seed \<SEED\>]
用-name设置名字，-seq设置序列，用-ss设置二级结构，用-seed设置种子，用-pdb设置起始结构，起始结构可以就用组装后得到的结构，也可以从组装加采样得到的结构中挑选一个，用-c或者-constraints加上约束。
-seed可以省略掉，这样默认的种子是11。
-constraints或者-c可以省略掉，表示不添加约束信息。
CONSTRAINTS_FILE文件里面需要包含约束信息，例如可以加进DCA分析的信息。
例子：
    8 23 10
    9 22 10
    10 21 10
这里前面两列代表碱基的序号，例如8代表第8个碱基，23代表第23个碱基，最后一列是碱基之间的最小距离。
因此这里8 23 10就代笔第8个碱基和第23个碱基之间的最小的距离是10Å，9 22 10就代表第9个碱基和第22个碱基之间的最小距离是10Å。

##结构聚类
    nsp cluster -list \<LIST_FILE\> -k \<NUMBER_OF_CLUSTERS\>
使用-list来设置需要对哪些结构进行聚类，用-k来设置聚类的数目。
经过一段时间的运行之后，会在屏幕上打印出聚类的结构
LIST_FILE文件包含了要聚类的结构的名字:
    test.sample.1.pdb
    test.sample.2.pdb
    test.sample.3.pdb
    test.sample.4.pdb
    test.sample.5.pdb
##RNA三级结构打分
1. 对单个结构打分
    3dRNAscore -s \<PDB_FILE\>
2. 对多个结构打分
    3dRNAscore -s:l \<LIST_FILE\>
LIST_FILE文件包含了要聚类的结构的名字:
    test.sample.1.pdb
    test.sample.2.pdb
    test.sample.3.pdb
    test.sample.4.pdb
    test.sample.5.pdb


