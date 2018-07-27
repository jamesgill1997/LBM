// D2Q9 2D Lattice Model

#include "Grid.h"
#include <iostream>


// Entry point
int main()
{
	// Create grid
	Grid grid(2.0, 1.0, 30, 0.001, 100, 0.08);

	// Run simulation for a certain number of timesteps
	for (int t = 0; t < 1000000; t++)
	{
		grid.timestep();
		std::cout << "t - " << t << 'r' << std::flush;
	}

	// Write out results
	grid.writeOut();

	return 0;
}
