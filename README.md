核酸分子二级结构及三级结构预测
Nucleic acid 2D and 3D Structure Prediction (NSP)

#安装
配置3dRNA预测的环境，在集群上里输入如下命令：
`source /home/wangjian/wangjian.sh`

#用法
##二级结构预测
###用自由能最小方法预测
nsp ss_pred -seq <SEQUENCE>

###结合自由能最小方法以及DCA预测的DI值进行二级结构预测
nsp ss_dca -seq <SEQUENCE> -di <DI_FILE> [-k <K>]
k值是用来设置读取前k*L个DI值，如果k=1，就代表读取前L个，如果k=0.5，就代表读取前L/2个。

##三级结构预测
###组装
1. 组装
命令是nsp assemble，并设置-name和-seq和-ss这3个参数的值
例子：
`nsp assemble -name test–seq AAAAUUUU –ss ‘(((())))’`
会生成一个叫test.assemble.pdb的文件，这就是组装好的结构。
2. 组装+采样
如果要采样，就加上-sample，然后用-num来设置数目。
例子：
`nsp assemble -name test-seq AAAAUUUU –ss ‘(((())))’ –sample –num 100`
这样就会生成一个test.assemble.pdb的文件，以及test.sample.1.pdb, test.sample.2.pdb,
test.sample.3.pdb, ..., test.sample.1000.pdb这100个采样结构。

###优化
命令是nsp mc3p，用-name设置名字，-seq设置序列，用-ss设置二级结构，用-seed设置种子，用-pdb设置起始结构，起始结构可以就用组装后得到的结构，也可以从组装加采样得到的结构中挑选一个
例子：
`nsp mc3p –name test –seqAAAAUUUU –ss ‘(((())))’ –seed 123 –pdb test.assemble.pdb`
经过一段时间的模拟后会生成test.sample.123.pdb这个结构，这里的123代表设置的种子。还会生成test.mc.123.pdb这个文件，将这个文件拷贝到windows下面，用pymol打开即可看到模拟的过程。
-seed可以省略掉，这样默认的种子是11
例子：
nsp mc3p –name test –seqAAAAUUUU –ss ‘(((())))’ –pdb test.assemble.pdb
上面这句命令等价于
`nsp mc3p –name test –seqAAAAUUUU –ss ‘(((())))’ –seed 11 –pdb test.assemble.pdb`
经过一段时间的模拟后会生成test.sample.11.pdb以及test.mc.11.pdb两个文件。
也可以在命令中用-c或者-constraints加上约束
例子：
`nsp mc3p –name test –seqAAAACCCAAAUUUUCCCCCCUUU –ss ‘((((......)))).........’ –pdb test.assemble.pdb –constraintstest.constraints`
或者
`nsp mc3p –name test –seqAAAACCCAAAUUUUCCCCCCUUU –ss ‘((((......)))).........’ –pdb test.assemble.pdb –ctest.constraints`
经过一段时间的模拟后会生成test.sample.11.pdb以及test.mc.11.pdb两个文件。
这里的test.constraints是一个文件，文件里面需要包含约束信息，例如可以加进DCA分析的信息。
例子：
```
8 23 10
9 22 10
10 21 10
```
这里前面两列代表碱基的序号，例如8代表第8个碱基，23代表第23个碱基，最后一列是碱基之间的最小距离。
因此这里8 23 10就代笔第8个碱基和第23个碱基之间的最小的距离是10Å，9 22 10就代表第9个碱基和第22个碱基之间的最小距离是10Å。

##结构聚类
命令是nsp cluster，使用-list来设置需要对哪些结构进行聚类，用-k来设置聚类的数目。
例子：
`nsp cluster –list test.list–k 5`
经过一段时间的运行之后，会在屏幕上打印出聚类的结构
这里test.list是一个文件，文件包含了要聚类的结构的名字
例子：
```
test.sample.1.pdb
test.sample.2.pdb
test.sample.3.pdb
test.sample.4.pdb
test.sample.5.pdb
```
##RNA三级结构打分
命令是3dRNAscore
1. 对单个结构打分
例子
3dRNAscore –stest.sample.1.pdb
2. 对多个结构打分
3dRNAscore –s:l test.list
这里test.list是一个文件，文件包含了要聚类的结构的名字
例子：
```
test.sample.1.pdb
test.sample.2.pdb
test.sample.3.pdb
test.sample.4.pdb
test.sample.5.pdb
```


