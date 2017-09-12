[![GPLv3](https://www.gnu.org/graphics/gplv3-88x31.png)](https://www.gnu.org/licenses/gpl.html)

# Foreword

* **pizza** is a high-performance numerical code that quasi-geostrophic and non-rotating convection in a 2-D annulus geometry. pizza solves for the Navier-Stokes equation including Coriolis force couple with a temperature equation under the Boussinesq approximation.  

* **pizza** uses Chebyshev polynomials in the radial direction and Fourier decomposition in the azimuthal direction. The time-stepping scheme relies on a semi-implicit [Crank-Nicolson]( https://en.wikipedia.org/wiki/Crank–Nicolson_method) for the linear terms of the equations and a [Adams-Bashforth](https://en.wikipedia.org/wiki/Linear_multistep_method) scheme for the non-linear terms and the Coriolis force.  

* **pizza** is written in Fortran and designed to be used on supercomputing clusters. It requires [CMake](https://cmake.org) [MPI](http://www.open-mpi.org/), [FFTW](http://www.fftw.org) and [LAPACK](http://www.netlib.org/lapack/) to be compiled and executed. Postprocessing functions written in python (requiring [matplotlib](http://matplotlib.org/) and [scipy](http://www.scipy.org/) are also provided to allow a useful data analysis.  

* **pizza** is a free software. It can be used, modified and redistributed under the terms of the [GNU GPL v3 licence](http://www.gnu.org/licenses/gpl-3.0.en.html).


# Quickly start using pizza

### 1) In order to check out the code, use the command

```sh
$ git clone https://github.com/magic-sph/pizza.git
```
or via SSH (it requires a public key):

```sh
$ git clone ssh://git@github.com/magic-sph/pizza.git
```

### 2) Set up your compiler and compile the code

Create a directory where the sources will be built

```sh
$ mkdir build
$ cd build
```
Set up your compilers

```sh
$ export FC=mpiifort
```
or 
```sh
$ export FC=mpif90
```

Compile and produce the executable 

```sh
$ cmake ..
$ make -j
```
The executable `pizza.exe` has been produced!

### 3) Go to the samples directory and check that everything is fine

```sh
$ cd pizza/samples
$ ./pizza_wizard.py --nranks 4 --mpicmd mpiexec
```

If everything is correctly set, all auto-tests should pass!

### 4) You're ready for a production run

```sh
$ mkdir run
$ cp path_to/pizza/build/pizza.exe .
$ cp path_to/pizza/samples/RotE1e3EkPump/input.nml .
```
    
Then change the input namelist to the setup you want and run the code:

```sh
$ mpiexec -n 4 ./pizza.exe input.nml
```

### 5) Data visualisation and postprocessing

a) Set-up your PYTHON environment ([ipython](http://ipython.org/), [scipy](http://www.scipy.org/) and [matplotlib](http://matplotlib.org/) are needed)

b) Add pizza to your PYTHONPATH

```sh
$ export PYTHONPATH=$PYTHONPATH:path_to/pizza/python
```

c) You can now import the python classes:

```python
python> from pizza import *
```

and use them to read time series, graphic files, movies, ...

```python
python> ts = PizzaTs(field='e_kin', all=True)
python> ...
