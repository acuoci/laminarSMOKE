# Validation case
The H2/N2 flames experimentally studied by Toro et al.  (Proceedings of the Combustion Institute, 30 (2005), p. 485-492), were chosen as a first test-case for the validation of the numerical code, because the fuel is H2, whose chemistry is well-known and can be described with a relatively small number of species. 

The investigated coflow flame is characterized by a fuel nozzle with an internal diameter of 9 mm, surrounded by an air-coflow annulus (i.d. 95 mm). The fuel stream is a mixture of 50% H2 and 50% N2 (by volume) at ambient temperature. Two different average fuel exit velocities are considered in the present work: 27 cm/s (Flame F2) and 50 cm/s (Flame F3). Measurements of temperature, H2, O2, H2O and OH mole fractions are available along the axis and in radial direction at 3, 10, 20 and 30 mm from the fuel nozzle.

The computational domain has a width of 95 mm and a length of 150 mm. Computations are performed on grids with two levels of refinement, containing respectively about 1,700 and 25,000 cells. The hydrogen–oxygen chemistry was described using the `POLIMI_H2CO12` kinetic mechanism.

# Procedure
## Step 1
We start our calculations using the coarse grid and the `laminarBuoyantPimpleSMOKE` unsteady solver (`01-initial-solution` folder). In order to ignite the flame, a spark is provided close to the fuel nozzle. The unsteady simulation can be carried out for a limited amount of time, because it will be used only as first-guess solution for the steady-state solver (see next step). For this particular case, 0.1 s are enough.
```sh
$ cd 01-initial-solution
$ blockMesh
$ laminarBuoyantPimpleSMOKE
```
## Step 2
At the end of Step 1 we have a first guess solution to apply the `laminarBuoyantSimpleSMOKE` steady-state solver. Thus, we can transfer the solution obtained in Step 1 into `02-steady-state` folder and run the steady-state solver. Only a small number of iterations is needed (4,000-5,000 are enough).
```sh
$ cd 02-steady-state
$ cp -r ../01-initial-solution/0.1 .
$ mv 0.1 0
$ rm -r 0.1/uniform
$ blockMesh
$ laminarBuoyantSimpleSMOKE
```
## Step 3
We are now ready to refine the mesh. From `03-fine-grid` folder we can remap the last solution obtained in Step 2 and run a steady state simulation on the new finer grid (about 25,000 cells). We can directly use the steady-state solver because we have a very good first-guess solution from the coarse grid. About 10,000 iterations are enough to reach the steady-state solution.
```sh
$ cd 03-fine-grid
$ blockMesh
$ mapFields ../02-steady-state/ -consistent -sourceTime 4000
$ laminarBuoyantSimpleSMOKE
```
## Step 4
We can now compare the numerical results with the experimental measurements, included in the `exp` folder. After running the `sample` utility, both for the coarse and the fine grids, the comparison can be directly carried out using gnuplot and the provided `comparison.plt` script.
