<h4 id="cluster">Clustering</h4>

    nsp cluster -list <LIST_FILE> -k <NUMBER_OF_CLUSTERS>

Set the structures that need to cluster with -list and set the amount with -k.

After a period of operation, the structure of the cluster will be printed out on the screen.

The list_file file contains the name of the structure to cluster:

<h4 id="rmsd">Calculating RMSD</h4>

    nsp rmsd <PDB_FILE_1> <PDB_FILE_2>

<h4 id="seq">Get the sequence of the molecule </h4>

    nsp seq <PDB_FILE>

<h4 id="ss">Get the secondary structure of the molecule </h4>

    nsp ss <PDB_FILE>

<h4 id="len">Get the number of residues in the molecule</h4>

    nsp len <PDB_FILE>

<h4 id="sub">Get of specified residues in the molecul</h4>

    nsp sub <PDB_FILE> -num <FRAG1> <FRAG2> <FRAG3> <FRAG4> ...

Each frag refers to a base segment, the format is a single residue number n or the specified starting point and end point of the fragment  `begin-end.`
for example, 1 represents the first residue, and 4-11 represents a fragment of fourth to eleventh bases.

<h4 id="rna">Only residues A, U, G and C are retained</h4>

    nsp rna <PDB_FILE>

<h4 id="rna-format">Only residues A, U, G and C containing the complete atoms are retained</h4>

    nsp rna <PDB_FILE> -format

    

    
