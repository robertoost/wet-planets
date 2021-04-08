using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;

public class CellMaster
{
    private Vector3 startLocation;       // Minimum location of grid
    private int x_size;                  // Width size of grid
    private int y_size;                  // Height size of grid
    private int z_size;                  // Depth size of grid
    private float cellSize;              // Size of cells within grid

    // Holds all cells in the grid
    CellGrid cells;

    // Constants
    private const float GRAV = -9.81f;
    public const float VISCOSITY = 0.89f;
    public const float FLUID_DENSITY = 1.00f;
    public const float AIR_DENSITY = 1.00f;
    public const float ATMOSPHERIC_PRESSURE = 101325f;

    // Start is called before the first frame update
    public CellMaster(Vector3 _startLocation, int _x_size, int _y_size, int _z_size, float _cellSize)
    {
        startLocation = _startLocation;
        x_size = _x_size;
        y_size = _y_size;
        z_size = _z_size;
        cellSize = _cellSize;

        cells = new CellGrid();
    }

    //2. Update the grid based on the marker particles (figure 4)
    // TODO: it seems like there might be too many air cells?
    public void updateGrid(Particle[] particles)
    {
        //set “layer” field of all cells to −1
        foreach ((int, int, int) key in cells.Keys)
        {
            Cell currentCell = cells[key];
            currentCell.layer = -1;
        }

        // update cells that currently have fluid in them by looping through particles
        for (int p = 0; p < particles.Length; p++)
        {
            // Get cell corresponding to particle
            (int i, int j, int k) = locationToCellIndex(particles[p].getPosition());
            Cell particleCell = cells[i, j, k];

            // If cell does not yet exist, create it
            if (!particleCell)
            {
                Vector3 cellLocation = cellIndexToLocation((i, j, k));
                if (withinBounds(cellLocation))
                {
                    particleCell = new Cell(Cell.CellType.FLUID, 0, cellLocation);
                    cells.addCell(particleCell, i, j, k);
                }
            } 
            // If cell exists, set cell type to fluid and layer to 0
            else if (particleCell.cellType != Cell.CellType.SOLID)
            {
                particleCell.cellType = Cell.CellType.FLUID;
                particleCell.layer = 0;
            }
        }

        // create a buffer zone around the fluid
        for(int layer = 1; layer <= 2; layer++)      // TODO: replace with for i = 1 to max(2, dkc f le)
        {
            // Temp storage for new cells, store both keys and cells themselves: TODO make more efficient?
            List<(int, int, int)> neighboursVisited = new List<(int, int, int)>();
            List<Cell> newCells = new List<Cell>();

            // Loop through cells existing either from previous timestep or cells that were added because they contain particles
            foreach ((int, int, int) key in cells.Keys)
            {
                // Select cell
                (int i, int j, int k) = key;
                Cell currentCell = cells[key];

                // If the cell is not a solid and in current layer visit all its neighbours, to expand grid
                if (currentCell.cellType != Cell.CellType.SOLID && currentCell.layer == layer - 1)
                {
                    // Create list of neighbouring cells.
                    (int, int, int)[] neighbourArray = {(i + 1, j, k), (i - 1, j, k), (i, j + 1, k), 
                        (i, j - 1, k), (i, j, k + 1), (i, j, k - 1) };

                    foreach((int, int, int) neighbourKey in neighbourArray)
                    {
                        // Check whether neighbour has already been visited
                        if (neighboursVisited.Contains(neighbourKey)) {
                            continue;
                        }

                        // Check whether neighbour cell already exists
                        Cell neighbour = cells[neighbourKey];
                        if (neighbour)
                        {
                            // If cell exists and is not yet marked as current fluid layer (layer == 0) and not solid
                            // make the cell air.
                            if(neighbour.layer == -1 && neighbour.cellType != Cell.CellType.SOLID)
                            {
                                neighbour.cellType = Cell.CellType.AIR;
                                neighbour.layer = layer;
                            }
                        }
                        else
                        {
                            // If cell does not yet exist, create it. If it is within simulation bounds, make it air
                            // otherwise make it solid.
                            Vector3 locationNeighbour = cellIndexToLocation(neighbourKey);
                            if (withinBounds(locationNeighbour))
                            {
                                neighbour = new Cell(Cell.CellType.AIR, layer, locationNeighbour);
                            }
                            else
                            {
                                neighbour = new Cell(Cell.CellType.SOLID, layer, locationNeighbour);
                            }

                            // Store cell and its key in the temporary arrays.
                            newCells.Add(neighbour);
                            neighboursVisited.Add(neighbourKey);
                        }
                    }

                }
            }

            // Commit temporay cells
            for(int c = 0; c < newCells.Count; c++)
            {
                (int i_neigh, int j_neigh, int k_neigh) = neighboursVisited[c];
                cells.addCell(newCells[c], i_neigh, j_neigh, k_neigh);
            }
        }

        // Find all unused cells (cells with layer -1)
        List<(int, int, int)> unusedCells = new List<(int, int, int)>();
        foreach ((int, int, int) key in cells.Keys)
        {
            Cell currentCell = cells[key];
            if (currentCell.layer == -1)
            {
                unusedCells.Add(key);
            }
        }

        // Remove unused cells.
        foreach ((int, int, int) key in unusedCells)
        {
            cells.removeCell(key);
        }
    }

    // 3. Advance the velocity field, u
    public void velocityUpdate(float timeStep)
    {
        //3a. Apply convection using a backwards particle trace.
        convection(timeStep);

        //3b. Apply external forces.
        externalForces(timeStep);

        //3c. Apply viscosity.
        viscosity(timeStep);

        //3d. Calculate the pressure to satisfy ∇·u=0.
        //3e. Apply the pressure.
        pressure(timeStep);

        //3f. Extrapolate fluid velocities into buffer zone.
        extrapolateVelocities(timeStep);

        //3g. Set solid cell velocities.
    }

    // Get velocity of particle at location location.
    public Vector3 getVelocity(Vector3 location)
    {
        // Return velocity of cell at grid coordinates.
        return cells[locationToCellIndex(location)].velocity;
    }

    // Get location in scene from cell index
    Vector3 cellIndexToLocation((int, int, int) cellIndex)
    {
        (int i, int j, int k) = cellIndex;
        float x_loc = startLocation.x + (float)(0.5 + i) * cellSize;
        float y_loc = startLocation.y + (float)(0.5 + j) * cellSize;
        float z_loc = startLocation.z + (float)(0.5 + k) * cellSize;
        Vector3 location = new Vector3(x_loc, y_loc, z_loc);
        return location;
    }

    // Get cell of particle at location location.
    (int, int, int) locationToCellIndex(Vector3 location)
    {
        // Find grid coordinates.
        int i = (int)((location.x - startLocation.x) / cellSize);
        int j = (int)((location.y - startLocation.z) / cellSize);
        int k = (int)((location.z - startLocation.y) / cellSize);

        // Return velocity of cell at grid coordinates.
        return (i, j, k);
    }

    // Check whether cell location is within simulation bounds
    bool withinBounds(Vector3 location)
    {
        float x_min = location[0] - 0.5f * cellSize;
        float y_min = location[1] - 0.5f * cellSize;
        float z_min = location[2] - 0.5f * cellSize;

        float x_max = location[0] + 0.5f * cellSize;
        float y_max = location[1] + 0.5f * cellSize;
        float z_max = location[2] + 0.5f * cellSize;

        return (x_min >= startLocation[0] && y_min >= startLocation[1] && z_min >= startLocation[2]
            && x_max <= (startLocation[0] + x_size) && y_max <= (startLocation[1] + y_size) && z_max <= (startLocation[2] + z_size));
    }

    // Replace all cells' tempvelocities with their current velocities.
    void commitVelocities() {
        foreach ((int, int, int) key in cells.Keys)
        {
            (int i, int j, int k) = key;
            Cell currentCell = cells[key];
            currentCell.velocity = currentCell.tempVelocity;
        }
    }

    // Find velocities invoked from convection.
    void convection(float timeStep)
    {
        foreach (Cell currentCell in cells)
        {
            // Extraction location of cell
            Vector3 locationCell = currentCell.location;
            float x = locationCell.x; float y = locationCell.y; float z = locationCell.z;
            
            // Find velocity
            Vector3 V = getVelocity(x, y, z);
            V = getVelocity(x + 0.5f * timeStep * V.x, y + 0.5f * timeStep * V.y, z + 0.5f * timeStep * V.z);

            // Update velocity of cell
            currentCell.tempVelocity = V;
        }

        commitVelocities();
    }

    // Get the interpolated velocity at a point in space.
    Vector3 getVelocity(float x, float y, float z)
    {
        Vector3 V;
        V.x = getInterpolatedValue(x / cellSize, y / cellSize - 0.5f, z / cellSize - 0.5f, 0);
        V.y = getInterpolatedValue(x / cellSize - 0.5f, y / cellSize, z / cellSize - 0.5f, 1);
        V.z = getInterpolatedValue(x / cellSize - 0.5f, y / cellSize - 0.5f, z / cellSize, 2);
        return V;
    }

    // Get an interpolated data value from the grid.
    // Some cells might not exist and will be skipped, see Figure 5 in the paper
    float getInterpolatedValue(float x, float y, float z, int index)
    {
        // Convert to integer to get the rounded down indices of the cell.
        (int i, int j, int k) = ((int)x, (int)y, (int)z);

        // Get the cells to interpolate between.
        Cell[] cellArray = {cells[i, j, k], cells[i + 1, j, k], cells[i, j + 1, k], cells[i + 1, j + 1, k],
                                  cells[i, j, k + 1], cells[i + 1, j, k + 1], cells[i, j + 1, k + 1], cells[i + 1, j + 1, k + 1] };
        
        // Get the summed velocity of every cell. If the cell doesn't exist, it won't be counted.
        float velocitySum = (i + 1 - x) * (j + 1 - y) * (k + 1 - z) * (cellArray[0] ? cellArray[0].velocity[index] : 0)
            + (x - i) * (j + 1 - y) * (k + 1 - z) * (cellArray[1] ? cellArray[1].velocity[index] : 0)
            + (i + 1 - x) * (y - j) * (k + 1 - z) * (cellArray[2] ? cellArray[2].velocity[index] : 0)
            + (x - i) * (y - j) * (k + 1 - z) * (cellArray[3] ? cellArray[3].velocity[index] : 0)
            + (i + 1 - x) * (j + 1 - y) * (z - k) * (cellArray[4] ? cellArray[4].velocity[index] : 0)
            + (x - i) * (j + 1 - y) * (z - k) * (cellArray[5] ? cellArray[5].velocity[index] : 0)
            + (i + 1 - x) * (y - j) * (z - k) * (cellArray[6] ? cellArray[6].velocity[index] : 0)
            + (x - i) * (y - j) * (z - k) * (cellArray[7] ? cellArray[7].velocity[index] : 0);

        // Get the summed weight for every cell. If a corresponding cell doesn't exist, don't count that weight.
        float weightSum = (i + 1 - x) * (j + 1 - y) * (k + 1 - z) * (cellArray[0] ? 1 : 0)
            + (x - i) * (j + 1 - y) * (k + 1 - z) * (cellArray[1] ? 1 : 0)
            + (i + 1 - x) * (y - j) * (k + 1 - z) * (cellArray[2] ? 1 : 0)
            + (x - i) * (y - j) * (k + 1 - z) * (cellArray[3] ? 1 : 0)
            + (i + 1 - x) * (j + 1 - y) * (z - k) * (cellArray[4] ? 1 : 0)
            + (x - i) * (j + 1 - y) * (z - k) * (cellArray[5] ? 1 : 0)
            + (i + 1 - x) * (y - j) * (z - k) * (cellArray[6] ? 1 : 0)
            + (x - i) * (y - j) * (z - k) * (cellArray[7] ? 1 : 0);

        // Find the weighted interpolated value
        return velocitySum / weightSum;
    }

    // Apply gravity to every cell's velocity.
    // TODO: Change to incorporate our planetary gravity.
    void externalForces(float timeStep)
    {
        foreach (Cell cell in cells)
        {
            cell.velocity += new Vector3(0, GRAV * timeStep, 0);
        }
    }

    // Calculates the viscosity of every liquid cell's fluid.
    void viscosity(float timeStep)
    {
        foreach ((int,int,int) key in cells.Keys)
        {
            (int i, int j, int k) = key;

            Cell currentCell = cells[key];

            // Create list of neighbouring cells.
            Cell[] cellArray = {cells[i + 1, j, k], cells[i - 1, j, k], 
                                cells[i, j + 1, k], cells[i, j - 1, k], 
                                cells[i, j, k + 1], cells[i, j, k - 1] };


            // TODO: Do we actually need to decrease the laplacian counter if a value is omitted?
            // calculate the laplacian by adding neighbouring velocities together.
            int laplacianCounter = 0;
            Vector3 laplacian = Vector3.zero;
            foreach (Cell cell in cellArray) {

                // If a cell exists, include it in the laplacian.
                if ( cell ) {
                    laplacian += cell.velocity;
                    laplacianCounter++;
                }
            }

            // Subtract current cell velocity multiplied by the amount of included neighbours.
            laplacian -= laplacianCounter * currentCell.velocity;
            
            // Set the viscous velicity as the temp velocity.
            currentCell.tempVelocity = currentCell.velocity + timeStep * VISCOSITY * laplacian;
        }

        commitVelocities();
    }

    void pressure(float timeStep)
    {
        int cellCount = cells.Count;

        // Copy the list of cell keys 
        (int, int, int)[] cellKeyList = { };
        for (int i = 0; i < 0; i++)
        {
            cells.Keys.CopyTo(cellKeyList, i);
        }

        // Create a dictionary of every cell with an associated index for their values in the A matrix.
        Dictionary<Cell, int> cellIndices = new Dictionary<Cell, int>();

        for (int i = 0; i < cellKeyList.Length; i++)
        {
            (int, int, int) cellKey = cellKeyList[i];
            Cell cell = cells[cellKey];
            cellIndices[cell] = i;
        }

        // Create a matrix A to store coefficients for every cell and their neighbours
        Matrix<float> A = CreateMatrix.Sparse<float>(cellCount, cellCount);

        // Create a vector B to solve the A matrix for P later on.
        Vector<float> B = CreateVector.Dense<float>(cellCount);

        foreach ((int, int, int) key in cellKeyList)
        {
            // Get the current key, cell, and cell index.
            (int i, int j, int k) = key;
            Cell currentCell = cells[key];
            int currentCellIndex = cellIndices[currentCell];

            Cell[] fluidNeighbours = {};    // Neighbours that are fluid cells
            int nonSolidCount = 0;          // Number of non-solid neighbours
            int airCount = 0;          // Number of non-solid neighbours

            // Retrieve all neighbours for the current cell.
            Cell xMax = cells[i + 1, j, k];
            Cell yMax = cells[i, j + 1, k];
            Cell zMax = cells[i, j, k + 1];
            Cell xMin = cells[i - 1, j, k];
            Cell yMin = cells[i, j - 1, k];
            Cell zMin = cells[i, j, k - 1];
            
            Cell[] neighbours = {xMax, xMin, yMax, yMin, zMax, zMin}; 

            // Loop through neighbouring cells to analyze their types.
            foreach (Cell neighbour in neighbours)
            {
                if (neighbour.cellType == Cell.CellType.SOLID ) {
                    continue;
                }
                
                // Count the non solid neighbours for A.
                nonSolidCount++;

                // If a cell is a fluid cell, add 1 at their respective collumn in the A matrix.
                if (neighbour.cellType == Cell.CellType.FLUID)
                {
                    int fluidCellIndex = cellIndices[neighbour];
                    A[currentCellIndex, fluidCellIndex] = 1;
                }
                else    // So in this case the cell type is AIR. Count for the B matrix.
                {
                    airCount++;
                }
            }

            // In matrix A, set the current cell to the negative amount of non-solid neighbours.
            A[currentCellIndex, currentCellIndex] = -nonSolidCount;


            // Calculate the divergence, setting to 0 for velocity components pointing into solid cells.
            float divergence = (xMax ? xMax.velocity[0] : 0) - (xMin && xMin.cellType != Cell.CellType.SOLID ? currentCell.velocity[0] : 0)
                + (yMax ? yMax.velocity[1] : 0) - (yMin && yMin.cellType != Cell.CellType.SOLID ? currentCell.velocity[1] : 0)
                + (zMax ? zMax.velocity[2] : 0) - (zMin && zMin.cellType != Cell.CellType.SOLID ? currentCell.velocity[2] : 0);


            // In matrix B, set a value based off the divergence in the velocity field for the current cell.
            B[currentCellIndex] = currentCell.density() * cellSize * divergence / timeStep - airCount * ATMOSPHERIC_PRESSURE;
        }

        //// Solve for the actual pressure.
        //Vector<float> pressure = A.QR().Solve(B);

        //int index = 0;
        //foreach ((int, int, int) key in cellKeyList)
        //{

        //    // Get the current key, cell, and cell index.
        //    (int i, int j, int k) = key;
        //    Cell currentCell = cells[key];
        //    int currentCellIndex = cellIndices[currentCell];
        //    index++;

        //    // Retrieve all neighbours for the current cell.
        //    Cell xMin = cells[i - 1, j, k];
        //    Cell yMin = cells[i, j - 1, k];
        //    Cell zMin = cells[i, j, k - 1];

        //    Cell[] neighbours = { xMin, yMin, zMin };
        //}
    }

    void extrapolateVelocities(float timeStep) {

    }
}
