My Master's dissertation project, 'Monte Carlo Simulation of Leading Order Top Quark Production'. This uses C++ to create a Monte Carlo simulation of the collision of two protons to produce top quarks. Essentially, it is a non-solvable integral which I instead do using a computational method.

Files:
Monte Carlo Simulation of Leading Order Top Quark Production.pdf - The scientific paper I handed in for my dissertation
gluonmatrix.hh, histogram.hh, qqmatrix.hh - Header files called by the main 'mc.cc' file
mc.cc - The main C++ code to produce the results
param.cfg - An editable text file with some parameters that can be varied (particle energies, iterations, output type). This can be altered by non-technical users
parambac.cfg - param.cfg backup file
