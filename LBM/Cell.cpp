#include "Cell.h"
#include "Constants.h"


// Default constructor
Cell::Cell()
{
}


// Destructor
Cell::~Cell()
{
}

// Custom constructor that takes arguments for initialisation
Cell::Cell(int xpos, int ypos, int label, double gravity_lbm)
{
	// Initialise position
	x = xpos;
	y = ypos;
	type = label;

	// Initialise macroscopic
	ux = 0.0;
	uy = 0.0;
	density = 1.0;

	// Forcing
	forcex = density * gravity_lbm;

	// Initialise the f values to their equilibrium values
	for (int v = 0; v < 9; v++)
	{
		f[v] = equilibrium(v);
		fnew[v] = equilibrium(v);
	}
}

// Method for equilibrium
double Cell::equilibrium(int i)
{
	double A, B;

	// Compute terms of the expansion
	A = (c[i][0] * ux) + (c[i][1] * uy);

	B = ((c[i][0] * c[i][0]) - (cs * cs)) * (ux * ux) +
		((c[i][1] * c[i][1]) - (cs * cs)) * (uy * uy) +
		2 * c[i][0] * c[i][1] * ux * uy;

	// Return freq
	return density * w[i] * (1 + (A / (cs * cs)) + (B / (2.0 * (cs * cs * cs * cs))));
}

// Method for updating macros
void Cell::updateMacro()
{
	// Reset current cell macroscopic to zero
	density = 0.0;
	ux = 0.0;
	uy = 0.0;
	
	for (int v = 0; v < 9; v++)
	{
		density += f[v];
		ux += f[v] * c[v][0];
		uy += f[v] * c[v][1];
	}

	// Add forcing to X momentum
	ux += 0.5 * forcex;

	// Divide by density to get velocity
	ux /= density;
	uy /= density;
}

// Method for BGK collision step
void Cell::collide(double tau)
{
	if (type != 1) return;

	// Conduct collision
	for (int v = 0; v < 9; v++)
	{
		fnew[v] += (1.0 / tau) * (equilibrium(v) - fnew[v]) + force(tau, v);
	}
}

// Method for computing Gup forcing
double Cell::force(double tau, int v)
{
	// Intermediate variable
	double beta_v = 0.0;

	// Compute the lattice forces based on Guo's forcing scheme
	double lambda_v = (1.0 - 0.5 * (1.0 / tau)) * (w[v] / (cs * cs));

	// Dot product (sum over d dimensions)
	beta_v = ((c[v][0] * ux) + (c[v][1] * uy) * (1.0 / (cs * cs)));

	// Compute force
	return (forcex * (c[v][0] * (1.0 + beta_v) - ux)) * lambda_v;
}