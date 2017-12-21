Included files:
jacobi.cpp -      Executable for the Jacobi Method. Reads in the file data.txt and displays the 
		  residual at each step (as well as the iteration count).
		  It then displays the solution matrix x as well as the residual.

gaussSeidel.cpp - Executable for the Gauss Seidel Method. Reads in the file data.txt and 
		  displays the residual at each step (as well as the iteration count).
		  It then displays the solution matrix x as well as the residual.

sor.cpp - 	  Executable for the SOR Method. Reads in the file data.txt and 
		  displays the residual at each step (as well as the iteration count).
		  It then displays the solution matrix x as well as the residual. It uses
		  1.2 as its relaxation factor.

Matrix.h - Matrix class declarations
Matrix.cpp - Matrix class definitions

To run:
Add Matrix.h and Matrix.cpp to a project and then add one of the method .cpp files at a time.
Because they all have a main() method, they can't be loaded together.