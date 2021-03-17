# N-Body-Simulator
An object oriented approach to computing N body dynamics, specifically in the context of Astrophysics, and in star clusters.
A seperate mode to simulate Kozai Oscialltions is also present.

Requirements:
Numpy
Matplotlib

The following integration schemes are available:
1.Euler
2.Leapfrog KDK 
3.Leapfrog DKD

The use of Leapfrog is recommended in the context of this program, and Euler is to be treated as a fallback 
in case of any issues with the Leapfrog integrators.

The code is optimized by the use of Numpy. Nonetheless, the base language is Python
and thus the code is not designed with heavy computation in mind. 

Authors:
1.Sagar Ramchandani
2.Mykhailo Lobodin
