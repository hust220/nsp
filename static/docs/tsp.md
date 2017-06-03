## Tertiary structure prediction

*   #### Assembly
        nsp assemble -name <JOB_NAME> -seq <SEQUENCE> -ss "<SECONDARY_STRUCTURE>" \
        -n <NUMBER_OF_PREDICTIONS>

*   #### Sampling
        nsp sample -name <JOB_NAME> -seq <SEQUENCE> -ss "<SECONDARY_STRUCTURE>" \
        -n <NUMBER_OF_PREDICTIONS>

*   #### Optimization

        nsp opt -name <JOB_NAME> \
        -seq <SEQUENCE> \
        -ss "<SECONDARY_STRUCTURE>" \
        -init <INITIAL_PDB_FILE> \
        [-dca <DCA_FILE>] \
        [-<constraints|c> <CONSTRAINTS_FILE>] \
        [-queue <QUEUE>] \
        [-seed <SEED>]

    <table>
    <tr><td>JOB_NAME</td><td>job name</td></tr>
    <tr><td>SEQUENCE</td><td>sequence<br>example: AAAAACCCCUUUUU</td></tr>
    <tr><td>SECONDARY_STRUCTURE</td><td>secondary structure with dot-bracket notation<br>example: (((((....)))))</td></tr>
    <tr><td>INITIAL_PDB_FILE</td><td>initial structure file with 'pdb' or 'cif' format</td></tr>
    <tr>
        <td>CONSTRAINTS_FILE</td>
        <td>
        constraints file<br>example:<br>
        8 23 10<br>
        9 22 10<br>
        10 21 10<br>
        The first two column represents the base sequence number, for example 8 represents the eighth base, 23 represents the twenty-third base, the last column is the minimum distance between the base.
        </td></tr>
    <tr><td>DCA_FILE</td><td>dca file</td></tr>
    <tr>
        <td>QUEUE</td>
        <td>
            Queue of optimization actions.<br>
            Optimization action example:
            <table>
                <tr><td>Simulated Annealing Monte Carlo simulation</td><td>samc:200000:20-0.01</td></tr>
                <tr><td>Monte Carlo simulation by always heating</td><td>heat:100000:20<br>heat:100000</td></tr>
                <tr><td>Monte Carlo simulation by always cooling</td><td>cool:100000:20<br>cool:100000</td></tr>
                <tr><td>Monte Carlo simulation by always warming</td><td>warm:100000:20<br>warm:100000</td></tr>
            </table>
            Queue example:<br>
            heat:30000+warm:100000+cool:1000000
        </td>
    </tr>
    <tr><td>SEED</td><td>default value of seed is 11</td></tr>
    </table>

*   <h4 id='traj-cluster'>Trajectory Clustering</h4>

        nsp traj cluster <TRAJECTORY_FILE> -[prefix|p] <PREFIX> -k <NUMBER_OF_CLUSTERS> -aa

    Example:

        nsp traj cluster example.traj.pdb -p example.cluster -k 10 -aa

    This command will generate 10 files:

        example.cluster.1.pdb
        example.cluster.2.pdb
        example.cluster.3.pdb
        example.cluster.4.pdb
        example.cluster.5.pdb
        example.cluster.6.pdb
        example.cluster.7.pdb
        example.cluster.8.pdb
        example.cluster.9.pdb
        example.cluster.10.pdb


