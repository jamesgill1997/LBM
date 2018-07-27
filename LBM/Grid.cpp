#include "Grid.h"
#include "Constants.h"
#include <iostream>
#include <fstream>

// Default constructor
Grid::Grid()
{
}

// Destructor
Grid::~Grid()
{
}

// Custom constructor to take user arguements and initialise the member variables
// Resolution is the number of cells for l
Grid::Grid(double width, double height, int resolution, double timestep, double re, double gravity)
{
	// Populate member variables
	reynolds = re;
	dt = timestep;
	dx = 1.0 / static_cast<double>(resolution);

	// Assign number of cells
	nx = static_cast<int>(std::floor(static_cast<double>(resolution)* width));
	ny = static_cast<int>(std::floor(static_cast<double>(resolution)* height));

	// Viscosity in LBM units
	nu = ((1.0 / re) * (dt / dx) * (1.0 / dx));

	// Relaxation time 
	tau = 0.5 + (nu / (cs * cs));

	// Work out the gravity force in LBM units
	double grav_lbm = gravity * dt * (dt / dx);

	// Initialise the cells array
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			// If cell at top and bottom make it a solid, else wall
			if (j == 0 || j == ny - 1)
			{
				cells.push_back(new Cell(i, j, 0, 0));
			}
			else
			{

				cells.push_back(new Cell(i, j, 1, grav_lbm));
			}
		}
	}
}

// Method to write velocity to file
void Grid::writeOut()
{
	// Create a CSV file
	std::ofstream file;
	file.open("./velocity.csv", std::ios::out);
	if (file.is_open())
	{
		// Line 1 = nx
		file << nx << ',';

		// Line 2 = ny
		file << ny;

		// Line 3+ = velocity
		for (size_t idx = 0; idx < cells.size(); idx++)
			file << ',' << std::sqrt(cells[idx]->ux * cells[idx]->ux + cells[idx]->uy * cells[idx]->uy);

		// Close file
		file.close();
	}
}

// Method to perform timestep
void Grid::timestep()
{
	// Stream calls collide
	stream();

	// Last step is to swap f and fnew points for the next time step
	for (int idx = 0; idx < nx * ny; idx++)
	{
		double * tmp = cells[idx]->f;
		cells[idx]->f = cells[idx]->fnew;
		cells[idx]->fnew = tmp;
	}
}

// Method for stream step
void Grid::stream()
{
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			// Get current cell 1D index
			int currentCell = j + i * ny;

			// Loop over lattice directions
			for (int v = 0; v < 9; v++)
			{
				// Get coordinates of where the population will come from
				int i_adj = i - c[v][0];
				int j_adj = j - c[v][1];

				// Periodic Boundary Conditions
				if (i_adj < 0) i_adj = nx - 1;
				if (j_adj < 0) j_adj = ny - 1;
				if (i_adj > nx - 1) i_adj = 0;
				if (j_adj > ny - 1) j_adj = 0;

				// Compute 1D index of the neighbour
				int adjacentCell = j_adj + i_adj * ny;

				// Read value from f of adjacent cell and write to fnew of current cell
				if (cells[adjacentCell]->type == 0)
				{
					cells[currentCell]->fnew[v] = cells[currentCell]->f[vopp[v]];
				}
				else
				{
					// Regular streaming
					cells[currentCell]->fnew[v] = cells[adjacentCell]->f[v];
				}
			}

			// Update macroscopic
			cells[currentCell]->updateMacro();

			// Collide
			cells[currentCell]->collide(tau);
		}
	}
}