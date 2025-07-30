# NACA Euler solver
Database generator for flow over NACA 4 digit airfoils - compressible Euler equations (wikipedia link: https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics)#Euler_equations).

**_NOTE:_** The database is not in the repo since it is too big.
**_NOTE:_** I wrote the solver and the mesher in a CFD class, and the reports I wrote are included.

### Build and Run
To build the program, run this commend:
``` shell
IN_FILE=input.txt OUT_DIR=./results make main
```
* input.txt - input file
* ./results - output directory

### Usage
``` shell
main 'input file' 'output directory'
```

## Dependencies
The only dependencies are:
* C compiler
* make files
* SQLite3

### Installing SQLite3 on Ubuntu
``` shell
sudo apt install sqlite3
```

## Mesh explanation
For example, the mesh for
* NACA - 0012
* ni - 95
* nj - 48
* num_points_on_airfoil - 61
* delta_y - 0.003
* XSF - 1.1
* YSF - 1.1
* r - 0.001
* omega - 0.1

will be: ![0012 mesh](./0012%20mesh.png) 

### Airfoil position (indexes)
``` C
int i_LE  = (ni - 1) / 2;
int i_TEL = i_LE - num_points_on_airfoil / 2;
int i_TEU = i_LE + num_points_on_airfoil / 2;
int j_LE  = 0;
int j_TEL = 0;
int j_TEU = 0;
```

## Solver explanation
For the given mesh above and the following parameters:
* Mach_inf - 0.9
* angle_of_attack_deg - 0
* density - 1.225
* environment_pressure - 101325
* delta_t - 1e-5
* Gamma - 1.4
* epse - 0.06

The flow around the airfoil (Mach number) is: ![0012 mesh](./flow%200012%20Mach0.9.png) 

## Example input file
``` shell
NACA
0012

ni
95

nj
48

num_points_on_airfoil
61

delta_y
0.003

XSF
1.1

YSF
1.1

r
0.001

omega
0.1

Mach_inf
0.9

angle_of_attack_deg
0

density
1.225

environment_pressure
101325

delta_t
1e-5

Gamma
1.4

epse
0.06

max_iteration
4e4
```
Where

| Parameter | Description |
| :---: | :---: |
| NACA | 4 digit NACA airfoil |
| ni | number of points in the i (xi) direction |
| nj | number of points in the j (eta) direction |
| num_points_on_airfoil | number of points on the airfoil |
| delta_y | initial cell size in the j direction |
| XSF | cell size in the i direction multiplayer |
| YSF | cell size in the j direction multiplayer |
| r | parameter for mesh convergence |
| omega | parameter for mesh convergence |
| Mach_inf | mach number far away |
| angle_of_attack_deg | angle of attack of the airfoil in degrees |
| density | air density far away (height) |
| environment_pressure | air pressure far away (height) |
| delta_t | time step between two iteration |
| Gamma | heat capacity ratio |
| epse | smoothing coefficient for Beam & Warming method |
| max_iteration | maximal number of iteration for flow solver |

**_NOTE:_** the order does not meters.
**_NOTE:_** the time step can be bigger than the time step calculated using the CFL number since we seek a steady state solution.

## Automation - Ubuntu
In order to automate a lot of runs together you can use the 'automate.c' file

### Build and Run
Run the following command:
``` shell
make automat
```
This will create the folder 'auto' in which there will be a lot of 'input.txt' files and one 'command_to_run' in which the will be something like so:
``` shell
make build_and_link_main
./build/main ./auto/input0.txt ./auto/results 
./build/main ./auto/input1.txt ./auto/results 
./build/main ./auto/input2.txt ./auto/results 
./build/main ./auto/input3.txt ./auto/results 
.
.
.
./build/main ./auto/input206.txt ./auto/results 
make clean_main
```
You need to copy all the text in the file and past it in the terminal.

**_NOTE:_** This might work on windows but I am not sure.