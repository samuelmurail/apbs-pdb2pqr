README for mg-auto sor APBS example
===================================

The example files included in this directory shows how to use the SOR option in the ELEC section of the input files. We will use 1FAS a potent acetylcholinesterase inhibitor from the green mamba snake venom at 1.9-A resolution.

The mg-auto method now includes an option to enable the use of the SOR method. To use the sor option, use mg-auto followed by the sor keyword as in the 1FAS.in file. This will allowed mg-auto whether to use full multigrid method or SOR with only one level. The deciding factor is primarily dependant in the linearity and grid size of the problem. 

Stopping criteria uses the residual with a standard 1.0E-9 error tolerenace or a predetermined number of iterations, whichever accurs first.

To run make sure that the apbs binary is in your path and enter in the command line:

apbs 1FAS-1.in

The output will be written out to the fas.dx file. 

input file| Description                                   | APBS Version| Results             |
----------|-----------------------------------------------|-------------|---------------------|
1FAS.in   |sequential, mg-auto sor, lpbe, srfm smol       | 1.4.1       | Global net ELEC energy = 										   2.6488841277737E+04


