## Secondary structure prediction

*   <b id="fem">Free energy minimization method</b>

        nsp ssp_fe -seq <SEQUENCE>

*   <b id="fem_dca">Combining free energy minimization method and DI value of the DCA prediction for secondary structure prediction</b>
    
        nsp ssp_dca -seq <SEQUENCE> -di <DI_FILE> [-k <K>]

    k value is used to set to read the first K\*L DI values,if k=1,on behalf of the read before L;if k=0.5,on behalf of the read before L/2. 

*   <b id="ss_tree">Show secondary structure tree</b>

        nsp ss_tree -seq <SEQ> -ss "<SS>"

    <SEQ> is the sequence, <SS> is the secondary structure.

    Example:

        nsp ss_tree -seq AAAACCCCUUUU -ss "((((....))))"

*   <b id="mcc">Calculate MCC</b>

        nsp mcc -nat "<SECONDARY_STRUCTURE_OF_NATIVE>" -pred "<SECONDARY_STRUCTURE_OF_PREDICTION>"

*   <b id="sty">Calculate STY</b>

        nsp sty -nat "<SECONDARY_STRUCTURE_OF_NATIVE>" -pred "<SECONDARY_STRUCTURE_OF_PREDICTION>"

*   <b id="ppv">Calculate PPV</b>

        nsp ppv -nat "<SECONDARY_STRUCTURE_OF_NATIVE>" -pred "<SECONDARY_STRUCTURE_OF_PREDICTION>"
