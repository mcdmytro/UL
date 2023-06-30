This project allows to:
- calculate UL on the number of signal events of the X states, namely X(4274), X(4685), 
X(4630), X(4500) and X(4700), for two decay channels: Ds2317 and Ds2460;
- estimate contribution into systematic uncertanties coming from parameters of the 
considered states. 

Usage: 
$ ./run_UL_calc.sh - Runs cross-section UL calculation for the parameters, which are 
                     currently written in config.h
$ ./process_all.sh - Runs cross-section UL calculation for all states, all mass and 
                     width variations. Exports results in readable format into
                     results.txt and calculates syst errors gain.

Included files:
- Data*.root    - input data files;
- config.h      - a config file that defines considered state, decay channel and 
                  variation of mass and width (for systematics study);
- datafiles.C   - imports datafile;
- cutdata.C     - applies cuts;
- X_states.C    - defines parameters of X states;
- UL.cxx        - a script that performs fit and exports its results into dat files;
- B_up_lik.dat  - fit output files;
- B_ulplot.C    - performs systematics smearing of a likelihood curve and calculates
                  N^UL and cross-section UL @ 90% CL;

- run_UL_calc.sh    - runs UL.cxx and B_ulplot.C in sequence;

- process_all.sh    - runs run_UL_calc.sh for all states, all mass and width variations
                      and exports results into results.txt and results_short.txt. Also
                      runs syst error gain calculation;
- results.txt       - readable results;
- results_short.txt - output of cross-section UL values only. Used for syst uncertanties
                      gain in errors_calc.
- errors_calc.C     - calculates syst errors gain based on the values in results_short.txt