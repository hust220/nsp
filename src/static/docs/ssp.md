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
