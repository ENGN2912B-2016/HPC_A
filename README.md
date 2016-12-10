#**HPC**

*Arjun Sharma*, *Aravind A.V.P.S*, *Sayan Samanta*

###**Course Project:** A Study of Phase Separation in Ternary Alloys

**__Summary__**

In the [paper](https://www.researchgate.net/publication/225799332_A_study_of_phase_separation_in_ternary_alloys), the authors have simulated a ternary compound (having 3 components namely A,B and C) which undergoes 'Spinodal Decomposition'. Spinodal Decomposition is a type of phase evolution where a thermodynamically unstable mixed system changes over time as the temperature is dropped. 

In such a case, typically the rich become richer and the poor becomes poorer in a sense that from the mixed alloy if there is a slight perturbation in composition, the elements of same type cluster together (note that such a behavior is not a very commonplace, where in general systems tend to get more random in stead of clustering among like atoms). The authors have modeled the free energy based on the Cahn- Hilliard equation and then have minimized it which resulted into the evolution equations. These time dependent evolution equations are then solved at each time to see how the 'phase separation/phase evolution' occurs.

By minimizing the free energy functional (cahn hilliard), we get PDE's in time. Since the intention is to model bulk properties, the boundary condition used is periodic boundary condition (which is that if we discretized the system into grid points starting from 0 to n in each direction, then the value at grid point -1 is the same as the value as grid point n). This allowed them to use Fourier Transforms (since the functions are now periodic).

On using semi-implicit finite difference method on the equation, we can convert the PDE's into a set of ODE's for each timestep that are to be solved. Now the advantage of using fourier transform is that we can easily calculate the differential in real space by taking the variables in reciprocal space (by fourier transforming it), solving the set of equation in the reciprocal space and then inverse fourier transforming them in real space. The values of the composition of each element in each time step can then be recorded and plotted. For this particular paper, they have also played with the interfacial energy parameter to observe it's effect on the morphology of the phases as they evolve with time.

The project involves use of  math libraries such as random number generator, complex numbers, fourier transforms and solving simultaneous linear equations. Also, the results can be visualized using a visualization tool such as gnuplot.

**__Structure of the Project__**

Both serial code and multi-thread code has been written  for 2D and 3D systems. In case fo 2D the system size is 512x512, however due to the large amount of data generated and associated time constraint, the 3D system size was kept at 128x128.

For each cases, the code was run for different parameters (interfacial energy), to duplicate the results of the paper mentioned. Each folder has a bash script for submitting a job at CCV, and also bash and gnuplot scripts to generate the animated gif and png images to visualize the data from the generated dataset. Each folder also has a CMakeLists file to build and make the executable. Typical execution steps (same steps apply to each folder) are given below.

**Compilation Running and Visualisations instruction**

```bash
$ cd <folder_name>/ #enter the particular folder
$ cmake CMakeLists.txt #build using cmake
$ make #generate the executable
```
This is a bit unreliable in CCV specially in the login node, we recommend that you log into an compute node by typing `$ interact` before running the cmake files.

On the executable is made, submit a job to CCV (also catch up some sleep, it is going to take some time)
```bash
$ sbatch ternary_<remaining_name>.sh #submit the job to CCV
```

Once the datafiles are generated, proceed to make the png files and the animated gif
```bash
$ cd output/
$ ./plotAnimation.sh #run the script file to generate gnuplot script files
$ gnuplot #open gnuplot
```

Inside GNUplot (note the filenames will differ a bit depending on 2D or 3D)
```bash
gnuplot> load "plotAnimation.gp" #generate the png files
gnuplot> load "plotAnimation_gif.gp #generate the aninmated gif
```

**Sample Visualisation**

![visual 2D](https://github.com/ENGN2912B/HPC_A/blob/master/timestep86000.png)
![visual 3D](https://github.com/ENGN2912B/HPC_A/blob/master/timestep79000.png)
