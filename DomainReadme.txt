Included files:
domainDecomp.cpp - Main executable. Creates the matrices based on the given domain,
		   boundary conditions, PDE, and artificial boundary conditions from
		   splitting the domain. The two A matricies Ai and Aj are then decomposed
		   via LU and solution vectors are calculated. Afterwards, the residual of both
		   matrices is displayed on screen. A loop is then started which switches
		   the boundary conditions as discussed in class, which iterates 50 times
		   (although the residual stabilizes after the 2 iteration). Finally, the
		   two solution vectors xi and xj are written to files.

Matrix.h - Matrix class declarations
Matrix.cpp - Matrix class definitions

To run:
Add all files into a project, compile, and run domainDecomp.cpp