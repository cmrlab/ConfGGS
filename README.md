# ConfGGS
This program package (ConfGGS) provides an implementation of CHARMM Force Field parameters, AMBER Force Field parameters and OPLS-AA Force Field parameters in Gradient Gravitational Search Algorithm (J. Comput. Chem. 2015, 36 (14), 1060–1068) towards Conformational Search for the Building Blocks of Proteins and flexible protein fragments. You can find ConfGGS at https://github.com/cmrlab/ConfGGS.
These codes work on MATLAB using Windows Operating Systems.
To run the programs follow the steps:
Step1: Read initial structure (Input.csv)
Step2: Read Cartesian coordinates and Force Field parameters (CHARMM/AMBER/OPLS-AA) The CHARMM Force Field parameters have been included here for all bonds (cod_bond.csv), all angles (cod_angles.csv) and all dihedral angles (cod_dihedral.csv) including charge.csv
Step3: Random initialization of initial structure with all parameter setup (Population Size, N=50, lower bound=-10 and upper bound=10, dimension (dim)=Total number of atoms*3, sigma=0.001, max_it=200
Step4: Read reference coordinates (ref.csv)
Step5: Evaluation of fitness for the randomly initialized structure (finderror.m)
Step6: Computation of bonded terms (findEA.m, findEB.m, findED.m) and non-bonded terms (findES.m, findELJ.m)for each updated structure using CHARMM Force Field parameters and ultimately total energy (findETotal.m)
Step7: Check conformation with desired fitness value (0.009)
Step8: If yes, Go to Step: 11
Step9: Conformation Search with desired fitness using GGS Algorithm (main.m)
Step10: Update Conformational Data Base and go to Step: 7
Step11: Return near optimal solution with fitness value 0.009
Step12: Write the results (Best_Coordinates.csv)
Step13: Write Output.csv
Step14: Write Output_bond.csv, Output_angle.csv and Output_dihedral.csv
Step15: Plot fitness versus iteration steps
Step16: End
MAE and MSD Computation were carried out using ref_dihedral and Output_dihedral
Example: The protein fragment PQITL can also be considered as example mentioned in the MS for evaluation purpose. Both Input and Output folders with required files have also been uploaded.
Note: Any kind of assistance required to setup and run the program package “ConfGGS” may contact Computational Modelling Research Laboratory, School of Chemistry (Autonomous), Sambalpur University, INDIA Contact E-mail id: prabhat@suniv.ac.in 
