#pragma once
#include <vector>
#include "Cell.h"

class Grid
{
public:

	// Constructor
	Grid();

	// Destructor
	~Grid();

	// Overloaded constructor with simulation setup parameters
	Grid(double width, double height, int resolution, double timestep,
		double re, double gravity);

	// Method for writing out velocity
	void writeOut();

	// Method for performing a single timestep
	void timestep();

private:

	// Lattice spacing (dimensionless units)
	double dx;

	// Timestep (dimensionless units
	double dt;

	// Relaxation time
	double tau;	
	
	// Viscosity
	double nu;

	// Reynolds number
	double reynolds;	

	// Array of cells representing the grid
	std::vector<Cell*> cells;

	// Number of cells in the X and Y directions
	int nx;
	int ny;

	// LBM methods
	void stream();
};