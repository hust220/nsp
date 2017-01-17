# Nucleic acids 2D and 3D Structure Prediction (NSP)

Project:[https://github.com/hust220/nsp](https://github.com/hust220/nsp)

## Catalog
1. INSTALLATION
    * Installation requirements 
    * Compile and install
2. USAGE
    * Secondary structure prediction
    * Tertiary structure prediction
    * Clustering
    * RNA tertiary structure scoring
    * Calculating RMSD
    * Retrieving sequence

# INSTALLATION

---

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

# USAGE

Before use,you need to set the environment variable NSP to the folder where the template library is located. 

    export NSP=<PATH/OF/TEMPLATES/LIBRARY>

## DCA

    nsp dca -method <METHOD> -in <FASTA_FILE> -out <DI_FILE> -n <N>

## Secondary structure prediction

*   Free energy minimization method 

        nsp ssp_fe -seq <SEQUENCE>

*   Combining free energy minimization method and DI value of the DCA prediction for secondary structure prediction 
    
        nsp ssp_dca -seq <SEQUENCE> -di <DI_FILE> [-k <K>]

    k value is used to set to read the first K\*L DI values,if k=1,on behalf of the read before L;if k=0.5,on behalf of the read before L/2. 

*   Calculate MCC

        nsp mcc -nat "<SECONDARY_STRUCTURE_OF_NATIVE>" -pred "<SECONDARY_STRUCTURE_OF_PREDICTION>"

*   Calculate STY

        nsp sty -nat "<SECONDARY_STRUCTURE_OF_NATIVE>" -pred "<SECONDARY_STRUCTURE_OF_PREDICTION>"

*   Calculate PPV

        nsp ppv -nat "<SECONDARY_STRUCTURE_OF_NATIVE>" -pred "<SECONDARY_STRUCTURE_OF_PREDICTION>"

## Tertiary structure prediction

*   Assembly
    1.  assemble 

            nsp assemble -name <JOB_NAME> -seq <SEQUENCE> -ss "<SECONDARY_STRUCTURE>" -n <NUMBER_OF_PREDICTIONS>

    2.  assemble and sampling
    
            nsp sample -name <JOB_NAME> -seq <SEQUENCE> -ss "<SECONDARY_STRUCTURE>" -n <NUMBER_OF_PREDICTIONS>

*   Optimization

        nsp opt -name <JOB_NAME> -seq <SEQUENCE> -ss "<SECONDARY_STRUCTURE>" -init <INITIAL_PDB_FILE> [-<constraints|c> <CONSTRAINTS_FILE>] [-seed <SEED>]

    Set name with -name, set sequence with -seq,set secondary structure with -ss ,set seed with -seed.

    Set file for storing the optimized structure with -out.

    Set the start structure with -init,the starting structure may be the structure after assembly, and the structure can be selected from the structure which is assemble with the sample.

    Add constraints with -constranits or -c

    \-seed can be omitted and the default seed is 11.

    \-constraints or \-c also can be omitted on behalf of do not add constraint information.

    Constraints_file file inside the information of  constraints, such as DCA analysis can be added to the file.

        8 23 10
        9 22 10
        10 21 10

    The first two column represents the base sequence number, for example 8 represents the eighth base, 23 represents the twenty-third base, the last column is the minimum distance between the base.

    Therefore,8 23 10 represents the minimum distance between the eighth base and the twenty-third base is 10Å,9 22 10 represents the minimum distance between the ninth base and the twenty-two base is 10Å.

## Clustering

    nsp cluster -list <LIST_FILE> -k <NUMBER_OF_CLUSTERS>`

Set the structures that need to cluster with -list and set the amount with -k.

After a period of operation, the structure of the cluster will be printed out on the screen.

The list_file file contains the name of the structure to cluster:

## RNA tertiary structure scoring

1.  score a single structure

        nsp score -s <PDB_FILE>

2.  score multiple structures

        nsp score -s:l <LIST_FILE>

    The list_file file contains the name of the structure to cluster. 

## Calculating RMSD

    nsp rmsd -s <PDB_FILE_1> <PDB_FILE_2>

## Get the sequence of the molecule 

    nsp seq -s <PDB_FILE>

## Get the number of residues in the molecule

    nsp len -s <PDB_FILE>

## Get of specified residues in the molecul

    nsp sub -s <PDB_FILE> -num <FRAG1> <FRAG2> <FRAG3> <FRAG4> ...

Each frag refers to a base segment, the format is a single residue number n or the specified starting point and end point of the fragment  `begin-end.`
for example, 1 represents the first residue, and 4-11 represents a fragment of fourth to eleventh bases.

## Remove unwanted rows from the molecule, leaving only atom lines

    nsp rna -s <PDB_FILE>

    

    
