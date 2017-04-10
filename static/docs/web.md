## 3dRNA

[Click to visit the 3dRNA web server](http://biophy.hust.edu.cn/3dRNA)

There are four types of tasks:

1.  3dRNA (only assembly)

    This type of task is the most basic fragment assembly method.
    It constructs a tree for the secondary structure, and then decomposes the tree into many elements.
    The tree is called SST (secondary structure tree), and the element is called SSE (smallest secondary element).
    3dRNA find 3D templates for each SSE, and then assemble all the templates into an integrated 3D structure.

2.  3dRNA with sampling (assembly -> sampling)

    During the sampling process, we will re-select the 3D template for each SSE at each step and then
    assemble the re-selected templates into a new structure.

3.  3dRNA with optimization (assembly -> optimization)

    We designed a new optimization method based on the Monte Carlo simulation. After the assembly is complete,
    3dRNA will optimize the structure to make it closer to the native state.

4.  Optimization (optimization)

    The optimization procedure can optimize any structure provided by the user.

![](static/images/3drna-task-type.png)

#### 3dRNA

1.  Select the type of the task to be '3dRNA'.

2.  (Optional) Input the email address to recieve results.

3.  Input the sequence.

    Only the four characters "A U G C" are allowed in the sequence.

4.  Input the secondary structure. 

    Please provide the dot-bracket form of the secondary structure.

    Example:

    a)  (((...)))

    b)  ((((((.[[[))))))........]]]

5.  Select the number of predictions.

    The default number is 5.

    The user can choose the integer from 1 to 10.

6.  Submit

    After the task is submitted, the page will jump to the results page.

![](static/images/3drna-assembly.png)

#### 3dRNA with sampling

1.  Select the type of the task to be '3dRNA with sampling'.

2.  (Optional) Input the email address to recieve results.

3.  Input the sequence.

    Only the four characters "A U G C" are allowed in the sequence.

4.  Input the secondary structure. 

    Please provide the dot-bracket form of the secondary structure.

    Example:

    a)  (((...)))

    b)  ((((((.[[[))))))........]]]

5.  Select the number of predictions.

    The default number is 5.

    The user can choose the integer from 1 to 10.

6.  Select the number of samplings.

    The default number is 1000.

7.  Submit

    After the task is submitted, the page will jump to the results page.

![](static/images/3drna-assembly-sampling.png)

#### 3dRNA with optimization

1.  Select the type of the task to be '3dRNA with optimization'.

2.  (Optional) Input the email address to recieve results.

3.  Input the sequence.

    Only the four characters "A U G C" are allowed in the sequence.

4.  Input the secondary structure. 

    Please provide the dot-bracket form of the secondary structure.

    Example:

    a)  (((...)))

    b)  ((((((.[[[))))))........]]]

5.  Select the number of predictions.

    The default number is 5.

    The user can choose the integer from 1 to 10.

6.  Copy the contents of the constraints file into the text area or just upload the constraints file.

    The type of the constraints can be 'DCA' or 'Distances'.

    'DCA' means the direct information calculated by direct coupling analysis.

    'Distances' means the distances between residues.

7.  Submit

    After the task is submitted, the page will jump to the results page.

![](static/images/3drna-assembly-optimization.png)

#### Optimization

1.  Select the type of the task to be 'Optimization'.

2.  (Optional) Input the email address to recieve results.

3.  Input the secondary structure. 

    Please provide the dot-bracket form of the secondary structure.

    Example:

    a)  (((...)))

    b)  ((((((.[[[))))))........]]]

4.  Copy the contents of the structure file into the text area or just upload the structure file.

5.  Select the number of predictions.

    The default number is 5.

    The user can choose the integer from 1 to 10.

6.  Copy the contents of the constraints file into the text area or just upload the constraints file.

    The type of the constraints can be 'DCA' or 'Distances'.

    'DCA' means the direct information calculated by direct coupling analysis.

    'Distances' means the distances between residues.

7.  Submit

    After the task is submitted, the page will jump to the results page.

![](static/images/3drna-optimization.png)

## 3dRNAscore

[Click to visit the 3dRNAscore web server](http://biophy.hust.edu.cn/3dRNAscore)

1. Upload structure file

2. Submit

![](static/images/3drnascore-intro.png)

## DCA

[Click to visit the DCA web server](http://biophy.hust.edu.cn/DCA)

Direct Coupling Analysis (DCA) is a statistical inference framework used to infer direct co-evolutionary couplings among residue pairs in multiple sequence alignments, which aims at disentangling direct from indirect correlations. 

Users only need to select the input type, output type, and then provide the necessary input according to the selected input type, and finally submitted on it.

![](static/images/dca-intro.png)

There are three types of input, the first is to provide only the sequence, the second is to provide multiple sequence alignment, the third is to provide direct interactions.

## 3dRPC

3dRPC is an automatic program of building three-dimensional RNA-protein complex structures by docking RNA and protein structure.

The input of the 3dRPC server is the structure of protein and RNA in PDB format:

    ATOM 84706 N  MET V 1 -134.594 111.719 180.266
    ATOM 84707 CA MET V 1 -134.339 111.991 178.856
    ATOM 84708 C  MET V 1 -133.342 113.135 178.673
    ATOM 84709 O  MET V 1 -133.677 114.186 178.125
    ATOM 84710 CB MET V 1 -135.647 112.297 178.123
    ATOM 84711 CG MET V 1 -135.482 112.480 176.627
    ATOM 84712 SD MET V 1 -134.220 111.380 175.959
    ATOM 84713 CE MET V 1 -134.680 109.827 176.724

You can upload the structures by clicking the corresponding buttons. Then, you can
input the number of predictions. The 3dRPC server can generate corresponding
number of complexes. Finally, submit your job by clicking the "Submit" button. Once
a job is submitted to the 3dRPC server, a web page is displayed to show the rate of
progress and the predicted result if the job is finished. You can download the
predicted complex structures in PDB format.

![](static/images/3drpc-intro.png)

