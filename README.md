##2D-Heat-Equation
================

As a final project for Computational Physics, I implemented the Crank Nicolson method for evolving partial differential equations and applied it to the two dimension heat equation.

In terms of stability and accuracy, Crank Nicolson is a very stable time evolution scheme as it is implicit. However, in two dimensions, implicit methods are considerably more complicated than their one dimensional counterparts.

In addition to the code I wrote (about which I am very open to feedback), I have uploaded the final writeup I submitted that has some additional figures of visualizations using heatmaps and analysis of errors using a variety of different methods.

Currently, the configuration allows for simulation on a 2d grid with fixed boundary conditions, so the natural extension would be to include time varying BC's. Also, at some point I would like to automate the graphing process, since I was tweaking bits of Python to create most of the figures.
