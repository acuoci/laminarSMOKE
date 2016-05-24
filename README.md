laminarSMOKE
============
CFD solver for laminar reacting flows with detailed kinetic mechanisms based on OpenFOAM and OpenSMOKE++ Library

Compulsory libraries
--------------------
- Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
- RapidXML (http://rapidxml.sourceforge.net/)
- Boost C++ (http://www.boost.org/)
- OpenSMOKE++ (alberto.cuoci@polimi.it)

Optional libraries
------------------
- Intel MKL (https://software.intel.com/en-us/intel-mkl)
- ISATLib (mauro.bracconi@polimi.it)
- ODEPACK (http://computation.llnl.gov/casc/odepack/odepack_home.html)
- DVODE (http://computation.llnl.gov/casc/odepack/odepack_home.html)
- DASPK (http://www.engineering.ucsb.edu/~cse/software.html)
- Sundials (http://computation.llnl.gov/casc/sundials/main.html)
- MEBDF (http://wwwf.imperial.ac.uk/~jcash/IVP_software/readme.html)
- RADAU (http://www.unige.ch/~hairer/software.html)

Compilation
-----------
Three different options are available to compile the code, according to the level of support for the solution of ODE systems. The ISATLib is needed only if you want to apply the In Situ Adaptive Tabulation (ISAT) technique.
1. Minimalist: no external, optional libraries are required. Only the native OpenSMOKE++ ODE solver can be used.
2. Minimalist + Intel MKL: only the native OpenSMOKE++ ODE solver can be used, but linear algebra operations are managed by the Intel MKL libraries
3. Complete: all the optional libraries are linked to the code, in order to have the possibility to work with different ODE solvers

1. Instructions to compile the Minimalist version
-------------------------------------------------
1. Open the `mybashrc.minimalist`, choose the version of OpenFOAM you are using (2.2, 2.3, 2.4, 3.0) and adjust the paths to the compulsory external libraries
2. Type: `source mybashrc.minimalist`
3. Go to Section 4

2. Instructions to compile the Minimalist+MKL version
-----------------------------------------------------
1. Open the `mybashrc.minimalist.mkl`, choose the version of OpenFOAM you are using (2.2, 2.3, 2.4, 3.0) and adjust the paths to the compulsory external libraries and the paths to the Intel MKL library
2. Type: `source mybashrc.minimalist.mkl`
3. Go to Section 4

3. Instructions to compile the Complete version
-----------------------------------------------------
1. Open the `mybashrc.complete`, choose the version of OpenFOAM you are using (2.2, 2.3, 2.4, 3.0) and adjust the paths to the compulsory external libraries and the Intel MKL library. You can choose the additional external libraries you want to add to laminarSMOKE, by modifying the `EXTERNAL_ODE_SOLVERS` variable: in particular `1` means that the support is requested, while `0` means that no support is requested. Obviously, for each requested library, you need to provide the correct path.
2. Type: `source mybashrc.complete`
3. Go to Section 4

4. Compile the solvers
-----------------------------------------------------
1. Compile the steady-state (accounting for buoyancy) solver: from the `solver/laminarBuoyantSimpleSMOKE` folder type `wmake`
2. Compile the steady-state solver: from the `solver/laminarSimpleSMOKE` folder type `wmake`
3. Compile the unsteady (accounting for buoyancy) solver: from the `solver/laminarBuoyantPimpleSMOKE` folder type `wmake`
4. Compile the unsteady solver: from the `solver/laminarPimpleSMOKE` folder type `wmake`
5. Compile the post-processor: from the `solver/laminarSMOKEpostProcessor` folder type `wmake`

Run your first case
-------------------
The folder `run/examples/Toro-Flames/F3` contains a simple test case (laminar coflow diffusion flame fed with hydrogen).

1. Unsteady simulation: Open the `laminarBuoyantPimpleSMOKE-Global` folder, build the mesh using the `blockMesh` utility, and run the case using the `laminarBuoyantPimpleSMOKE` solver. Even if you are interested in steady state conditions, we strongly suggest to always start with unsteady calculations to create a reasonable first-guess solution for the application of the steady state solver. In this case, you can stop the unsteady simulation after 50 ms of physical time.

2. Steady state simuation: you can now move to the `laminarBuoyantSimpleSMOKE-Global` folder. Copy the last time folder calculated by the unsteady solver (point 1 above), build the mesh using the `blockMesh` utility, and run the case using the `laminarBuoyantSimpleSMOKE` solver. In order to reach the steady state conditions, 5000-6000 iterations are enough.
