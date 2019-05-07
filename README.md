laminarSMOKE
============
CFD solver for laminar reacting flows with detailed kinetic mechanisms based on OpenFOAM and [OpenSMOKE++ framework][1]

If you use laminarSMOKE for your publications, we kindly ask you to cite the following papers:

> Cuoci, A., Frassoldati, A., Faravelli, T., Ranzi, E., 
> Numerical modeling of laminar flames with detailed kinetics based on the operator-splitting method
> (2013) Energy and Fuels, 27 (12), pp. 7730-7753, DOI: 10.1021/ef4016334
 
> Cuoci, A., Frassoldati, A., Faravelli, T., Ranzi, E., 
> A computational tool for the detailed kinetic modeling of laminar flames: Application to C2H4/CH4 coflow flames
> (2013) Combustion and Flame, 160 (5), pp. 870-886, DOI: 10.1016/j.combustflame.2013.01.011
 
> Cuoci, A., Frassoldati, A., Faravelli, T., Ranzi, E., 
> OpenSMOKE++: An object-oriented framework for the numerical modeling of reactive systems with detailed kinetic mechanisms 
> (2015) Computer Physics Communications, 192, pp. 237-264, DOI: 10.1016/j.cpc.2015.02.014

Compulsory libraries
--------------------
- Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
- RapidXML (http://rapidxml.sourceforge.net/)
- Boost C++ (http://www.boost.org/)
- OpenSMOKE++ (provided with the current version of laminarSMOKE)

Optional libraries
------------------
- Intel MKL (https://software.intel.com/en-us/intel-mkl)
- ODEPACK (http://computation.llnl.gov/casc/odepack/odepack_home.html)
- DVODE (http://computation.llnl.gov/casc/odepack/odepack_home.html)
- DASPK (http://www.engineering.ucsb.edu/~cse/software.html)
- Sundials (http://computation.llnl.gov/casc/sundials/main.html)
- MEBDF (http://wwwf.imperial.ac.uk/~jcash/IVP_software/readme.html)
- RADAU (http://www.unige.ch/~hairer/software.html)

Optional libraries (under testing)
----------------------------------
- ISATLib (mauro.bracconi@polimi.it)

Compilation
-----------
Three different options are available to compile the code, according to the level of support for the solution of ODE systems. The ISATLib is needed only if you want to apply the In Situ Adaptive Tabulation (ISAT) technique.
1. Minimalist: no external, optional libraries are required. Only the native OpenSMOKE++ ODE solver can be used.
2. Minimalist + Intel MKL: only the native OpenSMOKE++ ODE solver can be used, but linear algebra operations are managed by the Intel MKL libraries
3. Complete: all the optional libraries are linked to the code, in order to have the possibility to work with different ODE solvers

1. Instructions to compile the Minimalist version
-------------------------------------------------
1. Open the `mybashrc.minimalist`, choose the version of OpenFOAM you are using (2.2, 2.3, 2.4, 3.0, 4.x, dev) and adjust the paths to the compulsory external libraries
2. Type: `source mybashrc.minimalist`
3. Go to Section 4

2. Instructions to compile the Minimalist+MKL version
-----------------------------------------------------
1. Open the `mybashrc.minimalist.mkl`, choose the version of OpenFOAM you are using (2.2, 2.3, 2.4, 3.0, 4.x, dev) and adjust the paths to the compulsory external libraries and the paths to the Intel MKL library
2. Type: `source mybashrc.minimalist.mkl`
3. Go to Section 4

3. Instructions to compile the Complete version
-----------------------------------------------------
1. Open the `mybashrc.complete`, choose the version of OpenFOAM you are using (2.2, 2.3, 2.4, 3.0, 4.x, dev) and adjust the paths to the compulsory external libraries and the Intel MKL library. You can choose the additional external libraries you want to add to laminarSMOKE, by modifying the `EXTERNAL_ODE_SOLVERS` variable: in particular `1` means that the support is requested, while `0` means that no support is requested. Obviously, for each requested library, you need to provide the correct path.
2. Type: `source mybashrc.complete`
3. Go to Section 4

4. Compile the libraries
-----------------------------------------------------
1. Compile the user-defined boundary condition library: from the `libs/boundaryConditionsOpenSMOKE++` folder type `wmake`
2. Compile the customized radiation library: from the `libs/radiationOpenSMOKE++` folder type `wmake`
3. Go to Section 5

5. Compile the solvers
-----------------------------------------------------
1. Compile the steady-state (accounting for buoyancy) solver: from the `solver/laminarBuoyantSimpleSMOKE` folder type `wmake`
2. Compile the steady-state solver: from the `solver/laminarSimpleSMOKE` folder type `wmake`
3. Compile the unsteady (accounting for buoyancy) solver: from the `solver/laminarBuoyantPimpleSMOKE` folder type `wmake`
4. Compile the unsteady solver: from the `solver/laminarPimpleSMOKE` folder type `wmake`
5. Compile the post-processor: from the `solver/laminarSMOKEpostProcessor` folder type `wmake`
6. Go to Section 6

6. Compile the CHEMKIN Pre-Processor
-----------------------------------------------------
1. Compile the CHEMKIN Pre-Processor utility: from the `solvers/openSMOKEppCHEMKINPreProcessor` folder type `wmake`

Preprocessing of CHEMKIN files
-----------------------------------------------------
In order to run a simulation with laminarSMOKE, a CHEMKIN mechanism (kinetics, thermodynamic and transport properties) has to be pre-processed using the `openSMOKEppCHEMKINPreProcessor` utility (see Section 6 above). 
Examples of mechanisms ready to be pre-processed are available in the `run/kinetic-mechanisms` folder. In particular, for in each mechanism folder you can find the three files corresponding to the CHEMKIN input (kinetics, thermodynamics, and transport properties) and an additional `input.dic` file, containing the instructions for the `openSMOKEppCHEMKINPreProcessor`.
In order to pre-process a kinetic mechanisms, the operations to carry out are very simple. As an example, for `POLIMI_H2CO_1412` mechanism:
1. Go to the `run/kinetic-mechanisms/POLIMI_H2CO_1412`
2. Type `openSMOKEppCHEMKINPreProcessor`
3. If everything works correctly, a `kinetics-POLIMI_H2CO_1412` folder will be created, including the preprocessed CHEMKIN files (in XML folder). This is the folder which has to be supplied to the `laminarSMOKE` solver.

Run your first case
-----------------------------------------------------
The folder `run/tutorials/ToroFlames/F3/` contains a simple test case (laminar coflow diffusion flame fed with hydrogen).

1. Unsteady simulation: Open the `laminarBuoyantPimpleSMOKE` folder, build the mesh using the `blockMesh` utility, and run the case using the `laminarBuoyantPimpleSMOKE` solver. Even if you are interested in steady state conditions, we strongly suggest to always start with unsteady calculations to create a reasonable first-guess solution for the application of the steady state solver. In this case, you can stop the unsteady simulation after 50 ms of physical time.

2. Steady state simuation: you can now move to the `laminarBuoyantSimpleSMOKE` folder. Copy the last time folder calculated by the unsteady solver (point 1 above), build the mesh using the `blockMesh` utility, and run the case using the `laminarBuoyantSimpleSMOKE` solver. In order to reach the steady state conditions, 5000-6000 iterations are enough.

[1]: https://www.opensmokepp.polimi.it/
