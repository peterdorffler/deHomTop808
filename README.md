# deHomTop808
```
  __          __  __                      ______                   __       __      __     
 /\ \        /\ \/\ \                    /\__  _\                /'_ `\   /'__`\  /'_ `\   
 \_\ \     __\ \ \_\ \    ___     ___ ___\/_/\ \/   ___   _____ /\ \L\ \ /\ \/\ \/\ \L\ \  
 /'_` \  /'__`\ \  _  \  / __`\ /' __` __`\ \ \ \  / __`\/\ '__`\/_> _ <_\ \ \ \ \/_> _ <_ 
/\ \L\ \/\  __/\ \ \ \ \/\ \L\ \/\ \/\ \/\ \ \ \ \/\ \L\ \ \ \L\ \/\ \L\ \\ \ \_\ \/\ \L\ \
\ \___,_\ \____\\ \_\ \_\ \____/\ \_\ \_\ \_\ \ \_\ \____/\ \ ,__/\ \____/ \ \____/\ \____/
 \/__,_ /\/____/ \/_/\/_/\/___/  \/_/\/_/\/_/  \/_/\/___/  \ \ \/  \/___/   \/___/  \/___/ 
                                                            \ \_\                          
                                                             \/_/                          
```

An 808-line Matlab educational code for combined multi-scale topology optimisation and phasor-based dehomogenisation.

## Getting Started

The code is documented in the paper: ["Woldseth, R.V., Sigmund, O. & Jensen, P.D.L. An 808 line phasor-based dehomogenisation Matlab code for multi-scale topology optimisation. Struct Multidisc Optim 67, 205 (2024)."](https://doi.org/10.1007/s00158-024-03880-1).

The code was developed and tested using MATLAB, version R2023b, including MATLAB Image Processing Toolbox.

The code can also be executed without the MATLAB Image Processing Toolbox, but the behaviour may change, see paper for details.

The program is executed with the function ```deHomTop808()```.

Additional FE models are included in the repo:
- ```prepFEA_cant()```
- ```prepFEA_mbb()```
- ```prepFEA_db()```

Files for two-load bridge example includes:
- ```twoLoadBridge_80_48_Rank3_data.mat```
- ```getPas_2loadbridge.m```

Below is a Matlab code snippet of how to use and execute the code for both multi-scale topology optimisation, on-the-fly phasor-based dehomogenisation and post dehomogenisation.

### Matlab example
```
% Grid size
nelX = 60; nelY = 30;
% Volume fraction
volFrac = 0.3;
% Filter radius of thickness fields
rMin = 2;
% Relative thickness bounds
wMin = 0.1; wMax = 1.0;
% Dehomogenisation length-scale relative to element size
dMin = 0.2;
% Frequency of on-the-fly dehomogenisation
deHomFrq = 20;
% Post evaluation of dehomogenised result
eval = true;

%% Run multi-scale TO + dehomogenisation
[rhoPhys0,TO] = deHomTop808(nelX,nelY,volFrac,rMin,wMin,wMax,dMin,deHomFrq,eval); 

%% Re-run dehomogenisation with 0.5 dMin
rhoPhys1 = deHomTop808(nelX,nelY,volFrac,rMin,wMin,wMax,0.5*dMin,deHomFrq,eval,TO); 
```

## Help

Please send your comments or questions to: pdlj@dtu.dk

## Authors

This Matlab code was written by R. V. Woldseth, O. Sigmund and P. D. L. Jensen,
TopOpt Group, Department of Civil and Mechanical Engineering,
Technical University of Denmark,
DK-2800 Lyngby, Denmark.                                                

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

The authors acknowledge the financial support from the InnoTop VILLUM investigator project through the Villum Foundation and nTopology inc. Furthermore, the authors would like to express their gratitude to Dr. Federico Ferrari for valuable discussions during the preparation of this work.
