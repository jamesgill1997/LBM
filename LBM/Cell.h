#pragma once

class Cell
{
	friend class Grid;

public:

	// Constructor
	Cell();

	// Destructor
	~Cell();
	
	// Constructor providing position of cell and type
	Cell(int xpos, int ypos, int label, double gravity_lbm);

private:

	// Cell type 1 = fluid, 0 = solid
	int type;	

	// F values before streaming
	double * f = new double[9];		

	// F values after streaming
	double * fnew = new double[9];	

	// X and Y positions
	int x;				
	int y;		

	// X and Y velocities
	double ux;			
	double uy;	

	// Density
	double density;		

	// Force in x direction
	double forcex;		

	// Equilibrium function to compute f^eq
	double equilibrium(int v);

	// Update macros
	void updateMacro();

	// Collide
	void collide(double tau);

	// Compute force
	double force(double tau, int v);
};