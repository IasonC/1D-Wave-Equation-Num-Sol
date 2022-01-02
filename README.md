# 1D Wave Equation Solver

*MATLAB code to numerically solve a 1D Wave Equation and compare to exact solution*

## Scenario
The applications of the Wave Equation are vast, from image processing to system vibrations. The scenario for this particular solution is modelling the perpendicular displacement  of a power transmission cable of length L=1 over an arbitrary time interval from t=0s.

## Theoretical Underpinning of Numerical Method
The general n-dimension Wave Equation is\
<h1 align="center"> ![equation](https://latex.codecogs.com/svg.image?\frac{\partial&space;^2&space;u}{\partial&space;t^2}&space;=&space;a^2&space;\nabla^2&space;u) </h1>

In 1 dimension, the 1D Wave Equation is reduced to\
<h1 align="center"> ![equation](https://latex.codecogs.com/svg.image?%5Cfrac%7B%5Cpartial%20%5E2%20u%7D%7B%5Cpartial%20t%5E2%7D%20=%20a%5E2%20%5Cfrac%7B%5Cpartial%20%5E2%20u%7D%7B%5Cpartial%20x%5E2%7D) </h1>

Moreover, central difference approximations can be used to represent the second-order partial derivatives:
<p align="center">
    ![equation](https://latex.codecogs.com/svg.image?%5Cfrac%7B%5Cpartial%20%5E2%20u%7D%7B%5Cpartial%20x%5E2%7D%20=%20%5Cfrac%7Bu(x-%5CDelta%20x,%20t)%20-%202%20u(x,t)%20&plus;%20u(x&plus;%20%5CDelta%20x,%20t)%7D%7B%5CDelta%20x%5E2%7D%20&plus;%20O(%5CDelta%20x%5E2)%20%5Capprox%20%5Cfrac%7Bu_%7Bi-1,%20j%7D-2u_%7Bi,j%7D%20&plus;%20u_%7Bi&plus;1,j%7D%7D%7B%5CDelta%20x%5E2%7D)
    ![equation](https://latex.codecogs.com/svg.image?%5Cfrac%7B%5Cpartial%20%5E2%20u%7D%7B%5Cpartial%20t%5E2%7D%20=%20%5Cfrac%7Bu(x,%20t-%5CDelta%20t)%20-%202%20u(x,t)%20&plus;%20u(x,%20t&plus;%20%5CDelta%20t)%7D%7B%5CDelta%20t%5E2%7D%20&plus;%20O(%5CDelta%20t%5E2)%20%5Capprox%20%5Cfrac%7Bu_%7Bi,%20j-1%7D-2u_%7Bi,j%7D%20&plus;%20u_%7Bi,j&plus;1%7D%7D%7B%5CDelta%20t%5E2%7D)
</p>

Dropping the error terms and substituting the difference equations into the 1D Wave Equation, the following computational molecule is obtained:\
<p align="center">
    ![equation](https://latex.codecogs.com/svg.image?u_%7Bi,j&plus;1%7D%20=%20r%5E2%20u_%7Bi-1,j%7D%20&plus;%20%202(1%20-%20r%5E2%20)u_%7Bi,j%7D%20&plus;%20r%5E2%20u_%7Bi&plus;1,j%7D%20-%20u_%7Bi,j-1%7D)
</p>

Numerically, given the desired x and t step-sizes (h and k respectively), the x-t plane can be partitioned into a mesh, with mesh points of coordinates (ih, jk) for i = {0, 1, ..., xmax/h} and j = {0, 1, ..., tmax/k}. The numerical solution to the 1D Wave Equation is a set of u-values associated with each mesh point which forms the solution function u(x,t). In the mesh array, j gives the row and i gives the column index.

The numerical solution in MATLAB uses the computational molecule applied to 2 subsequent rows (j and j-1 rows) to find the next row (j+1 row). This is an **iterative explicit mesh method**. Row 0 and the first & last columns of the mesh grid are given by boundary conditions. Also, using the boundary condition that the first derivative of u partial to time is zero, with a backward difference approximation at t=0, a row j=-1 is defined where each mesh point is equal to the mesh point of same i in the j=0 row above. These two rows (-1 and 0) given by the boundary conditions are the first iteration used to find row 1 with the computational molecule. Then rows 0 and 1 are used, etc ... No boundary condition is given for an upper j row, so this iterative method is effective as it can be used for any time length from 0s to arbitrary t.
