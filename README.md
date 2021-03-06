# 1D Wave Equation Solver

*MATLAB code to numerically solve a 1D Wave Equation and compare to exact solution*

## Scenario
The applications of the Wave Equation are vast, from image processing to system vibrations. This code is directly applied to a perpendicularly-vibrating system linearly defined along the x-axis with some length L (so solving in the x-dimension between bounds 0 and L) for arbitrary time since t = 0s (in the solution, time bounds 0 and 0.5s are chosen). This is a driven system due to an external environmental applied force on the system. The code is automated, so you can change the parameters to observe the solution for a longer time interval or with finer step sizes.

The solution to the 1D Wave Equation presented here is very versatile: it can be applied to a plethora of perpendicularly-vibrating systems, from strings or cables to bridges! Also, the code used to solve the 1D Wave Equation is highly dynamic: the parameters can be set to control the x-dimension length L of the system, the time interval studied in the solution, and even the initial displacement of the system at t = 0s due to external stimuli!

## Boundary Conditions (Cauchy Conditions)
The system is:
- fixed at both ends ![equation](https://latex.codecogs.com/svg.image?%5CRightarrow%20u(x=0,%20t)%20=%20u(0,t)%20=%200,%20u(x=L,%20t)%20=%20u(L,%20t)%20=%200)
- initially at rest ![equation](https://latex.codecogs.com/svg.image?%5CRightarrow%20%5Cpartial%20u(x,%20t=0)/%5Cpartial%20t%20=%20%5Cpartial%20u/%5Cpartial%20t%20%5Crvert_%7Bt=0%7D%20=%200)
- initially displaced to a general position f(x) = sin(πx) ![equation](https://latex.codecogs.com/svg.image?%5CRightarrow%20u(x,%20t=0)%20=%20u(x,0)%20=%20f(x)%20=%20%5Csin(%5Cpi%20x))

_Note that the initial displacement f(x) of the system can be set to a different function to observe the system response to a different initial stimulus._

## Theoretical Underpinning of Numerical Method
The general n-dimension Wave Equation is\
![equation](https://latex.codecogs.com/svg.image?\frac{\partial&space;^2&space;u}{\partial&space;t^2}&space;=&space;a^2&space;\nabla^2&space;u)

In 1 dimension, the 1D Wave Equation is reduced to\
![equation](https://latex.codecogs.com/svg.image?%5Cfrac%7B%5Cpartial%20%5E2%20u%7D%7B%5Cpartial%20t%5E2%7D%20=%20a%5E2%20%5Cfrac%7B%5Cpartial%20%5E2%20u%7D%7B%5Cpartial%20x%5E2%7D)

Moreover, central difference approximations can be used to represent the second-order partial derivatives:
![equation](https://latex.codecogs.com/svg.image?%5Cfrac%7B%5Cpartial%20%5E2%20u%7D%7B%5Cpartial%20x%5E2%7D%20=%20%5Cfrac%7Bu(x-%5CDelta%20x,%20t)%20-%202%20u(x,t)%20&plus;%20u(x&plus;%20%5CDelta%20x,%20t)%7D%7B%5CDelta%20x%5E2%7D%20&plus;%20O(%5CDelta%20x%5E2)%20%5Capprox%20%5Cfrac%7Bu_%7Bi-1,%20j%7D-2u_%7Bi,j%7D%20&plus;%20u_%7Bi&plus;1,j%7D%7D%7B%5CDelta%20x%5E2%7D)
![equation](https://latex.codecogs.com/svg.image?%5Cfrac%7B%5Cpartial%20%5E2%20u%7D%7B%5Cpartial%20t%5E2%7D%20=%20%5Cfrac%7Bu(x,%20t-%5CDelta%20t)%20-%202%20u(x,t)%20&plus;%20u(x,%20t&plus;%20%5CDelta%20t)%7D%7B%5CDelta%20t%5E2%7D%20&plus;%20O(%5CDelta%20t%5E2)%20%5Capprox%20%5Cfrac%7Bu_%7Bi,%20j-1%7D-2u_%7Bi,j%7D%20&plus;%20u_%7Bi,j&plus;1%7D%7D%7B%5CDelta%20t%5E2%7D)

Dropping the error terms and substituting the difference equations into the 1D Wave Equation, the following computational molecule is obtained, where the Courant Number r = aΔt/Δx\
![equation](https://latex.codecogs.com/svg.image?u_%7Bi,j&plus;1%7D%20=%20r%5E2%20u_%7Bi-1,j%7D%20&plus;%20%202(1%20-%20r%5E2%20)u_%7Bi,j%7D%20&plus;%20r%5E2%20u_%7Bi&plus;1,j%7D%20-%20u_%7Bi,j-1%7D)

Numerically, given the desired x and t step-sizes (Δx = h and Δt = k respectively), the x-t plane can be partitioned into a mesh, with mesh points of coordinates (ih, jk) for i = {0, 1, ..., xmax/h} and j = {0, 1, ..., tmax/k}. The numerical solution to the 1D Wave Equation is a set of u-values associated with each mesh point which forms the solution function u(x,t). In the mesh array, j gives the row and i gives the column index.

The numerical solution in MATLAB uses the computational molecule applied to 2 subsequent rows (j and j-1 rows) to find the next row (j+1 row). This is an **iterative explicit mesh method**. Row 0 and the first & last columns of the mesh grid are given by boundary conditions. Also, using the boundary condition that the first derivative of u partial to time is zero, with a backward difference approximation at t=0, a row j=-1 is defined where each mesh point is equal to the mesh point of same i in the j=0 row above. These two rows (-1 and 0) given by the boundary conditions are the first iteration used to find row 1 with the computational molecule. Then rows 0 and 1 are used, etc ... No boundary condition is given for an upper j row, so this iterative method is effective as it can be used for any time length from 0s to arbitrary t.

## Analytical Solution
The analytical solution to the 1D Wave Equation that I calculated with Separation of Variables is:\
![equation](https://latex.codecogs.com/svg.image?u(x,t)%20=%20%5Csum_%7Bn=1%7D%5E%7B%5Cinfty%7D%20b_n%20%5Csin%20%5Cleft(%20%5Cfrac%7Bn%5Cpi%7D%7BL%7D%20x%5Cright)%20%5Ccos%20%5Cleft(%20%5Cfrac%7Ban%5Cpi%7D%7BL%7D%20t%20%5Cright))\
where\
![equation](https://latex.codecogs.com/svg.image?b_n%20=%20%5Cfrac%7B2%7D%7BL%7D%20%5Cint_0%5EL%20f(x)%20%5Csin%7B%5Cleft(%5Cfrac%7B2%5Cpi%7D%7B2L%7D%20nx%20%5Cright)%7D%20%5CLongrightarrow%20b_n%20=%20%5Cfrac%7B2%7D%7BL%7D%20%5Cint_0%5EL%20f(x)%20%5Csin%7B%5Cleft(%5Cfrac%7B%5Cpi%20nx%7D%7BL%7D%20%5Cright)%7D%20%5C,dx)

For the given boundary condition where f(x) = sin(πx), bn = 1 for n=1 but bn = 0 for integer n>1 as then bn is the integral of orthogonal functions. Hence, here the exact solution is u(x,t) = sin(πx)cos(πt). This is compared to the numerical solution to show the effectiveness of this numerical method.


<!--- equations displayed with https://editor.codecogs.com/ -- LATEX -> URL encoded --->
