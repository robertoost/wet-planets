using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Sparse;
using MathNet.Numerics.LinearAlgebra.Sparse.Linear;



[Serializable]
public class CellMaster
{
    // User chosen variables in SimulationMaster
    public Vector3 startLocation;       // Minimum location of grid
    public int x_size = 5;                  // Max Width size of grid
    public int y_size = 5;                  // Max Height size of grid
    public int z_size = 5;                  // Max Depth size of grid
    public float cellSize = 0.2f;              // Size of cells within grid

    public float cflTimestepConstant = 2;
    public float minDeltaTime = 0.02f;
    private int minBorderLayers = 2;            // Min number of border layers of type air/solid around fluid

    public float atmospheric_pressure = 101.325f;
    public float viscosity_const = 1.0016f;

    // Parameters for solving for the pressure vector, influence program speed
    public float absTolerance = 1e-30f;
    public float relTolerance = 1e-5f;

    // Holds all cells in the grid
    public CellGrid cells = new CellGrid();

    // Constants
    private const float GRAV = -9.81f;
    public const float FLUID_DENSITY = 1f;
    public const float AIR_DENSITY = 1f;

    //2. Update the grid based on the marker particles (figure 4)
    public void updateGrid(List<Particle> particles)
    {
        //set “layer” field of all cells to −1
        foreach ((int, int, int) key in cells.Keys)
        {
            Cell currentCell = cells[key];
            currentCell.layer = -1;
        }

        // update cells that currently have fluid in them by looping through particles
        for (int p = 0; p < particles.Count; p++)
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
        for(int layer = 1; layer <= Mathf.Max(minBorderLayers, Mathf.Ceil(cflTimestepConstant)); layer++)
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
                if (currentCell.layer == layer - 1)
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
                            if(neighbour.layer == -1)
                            {
                                if(neighbour.cellType != Cell.CellType.SOLID) { neighbour.cellType = Cell.CellType.AIR; }
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
        extrapolateVelocities();

        //3g. Set solid cell velocities.
        setSolidCells();
    }

    // Get velocity of particle at location location.
    public Vector3 getParticleVelocity(Vector3 location)
    {
        Cell cell = cells[locationToCellIndex(location)];
        Vector3 velocity = cell ? cell.velocity : Vector3.zero;
        // Return velocity of cell at grid coordinates.
        return velocity;
        //return getVelocity(location); // TODO: switch to this version
    }

    // Used for determining timestep
    public float getMaxVelocityMagnitude()
    {
        float maxVelocity = 0f;
        foreach ((int, int, int) key in cells.Keys)
        {
            float currentVelocity = cells[key].velocity.sqrMagnitude;
            maxVelocity = Mathf.Max(currentVelocity, maxVelocity);
        }
        return Mathf.Sqrt(maxVelocity);
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
        return (location.x > startLocation.x && location.y > startLocation.y && location.z > startLocation.z
            && location.x < (startLocation.x + x_size) && location.y < (startLocation.y + y_size) && location.z < (startLocation.z + z_size));
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
        foreach ((int, int, int) key in cells.Keys)
        {
            // Get cell
            (int i, int j, int k) = key;
            Cell currentCell = cells[key];

            // Extraction location of cell
            Vector3 locationCell = currentCell.location;

            // Find velocity
            Vector3 V = getVelocity(locationCell);
            V = getVelocity(locationCell + 0.5f * timeStep * V);

            // Get neighbours into which velocity components point
            Cell xMin = cells[i - 1, j, k];
            Cell yMin = cells[i, j - 1, k];
            Cell zMin = cells[i, j, k - 1];

            // Update velocity of cell
            currentCell.tempVelocity = new Vector3();
            if (xMin && (xMin.cellType == Cell.CellType.FLUID || currentCell.cellType == Cell.CellType.FLUID)) { currentCell.tempVelocity.x = V.x; }
            if (yMin && (yMin.cellType == Cell.CellType.FLUID || currentCell.cellType == Cell.CellType.FLUID)) { currentCell.tempVelocity.y = V.y; }
            if (zMin && (zMin.cellType == Cell.CellType.FLUID || currentCell.cellType == Cell.CellType.FLUID)) { currentCell.tempVelocity.z = V.z; }
            
        }

        commitVelocities();
    }

    // Get the interpolated velocity at a point in space.
    // TODO incorporate startLocation
    Vector3 getVelocity(Vector3 location)
    {
        Vector3 V;
        float x = location.x;   float y = location.y;   float z = location.z;
        V.x = getInterpolatedValue(x / cellSize, y / cellSize - 0.5f, z / cellSize - 0.5f, 0);
        V.y = getInterpolatedValue(x / cellSize - 0.5f, y / cellSize, z / cellSize - 0.5f, 1);
        V.z = getInterpolatedValue(x / cellSize - 0.5f, y / cellSize - 0.5f, z / cellSize, 2);
        return V;
    }

    // Get an interpolated data value from the grid.
    // Some cells might not exist and will be skipped, see Figure 5 in the paper
    float getInterpolatedValue(float x, float y, float z, int index)
    {
        int i = (int)x;
        int j = (int)y;
        int k = (int)z;

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
    // TODO: Change to incorporate our planetary gravity. In this case also check other neighbours.
    void externalForces(float timeStep)
    {
        foreach ((int, int, int) key in cells.Keys)
        {
            // Get cell
            (int i, int j, int k) = key;
            Cell currentCell = cells[key];

            // Only apply to velocity components that border fluid cells            
            Cell yMin = cells[i, j - 1, k];
            if (yMin && (yMin.cellType == Cell.CellType.FLUID || currentCell.cellType == Cell.CellType.FLUID))
            {
                currentCell.velocity += new Vector3(0, GRAV * timeStep, 0);
            }
        }
    }

    // Calculates the viscosity of every liquid cell's fluid.
    void viscosity(float timeStep)
    {
        foreach ((int,int,int) key in cells.Keys)
        {
            (int i, int j, int k) = key;

            Cell currentCell = cells[key];

            // Retrieve all neighbours for the current cell.
            Cell xMin = cells[i - 1, j, k];
            Cell yMin = cells[i, j - 1, k];
            Cell zMin = cells[i, j, k - 1];

            (int, int, int)[] neighbours = { (i + 1, j, k), (i, j + 1, k), (i, j, k + 1), (i - 1, j, k), (i, j - 1, k), (i, j, k - 1) };

            // calculate the laplacian by adding neighbouring velocities together.
            int laplacianCounter = 0;
            Vector3 laplacian = Vector3.zero;
            foreach ((int, int, int) neighbourKey in neighbours) {

                // Extract indices
                (int i_neigh, int j_neigh, int k_neigh) = neighbourKey;
                Cell neighbour = cells[neighbourKey];

                // If a cell exists include it in the laplacian, select only velocity components that point 
                // into fluid cells using borderFluidVelocityComponents.
                Vector3 relevantVelocityComponents = borderFluidVelocityComponents(i_neigh, j_neigh, k_neigh);
                if (neighbour && relevantVelocityComponents != new Vector3(0, 0, 0)) {
                    laplacian += relevantVelocityComponents;
                    laplacianCounter++;
                }
            }

            // Subtract current cell velocity multiplied by the amount of included neighbours.
            laplacian -= laplacianCounter * currentCell.velocity;   // TODO: FIX laplacianCounter

            // Get velocity from viscosity
            Vector3 newVelocity = currentCell.velocity + timeStep * viscosity_const * laplacian;

            // Set the viscous velicity as the temp velocity
            currentCell.tempVelocity.x = ((xMin && (xMin.cellType == Cell.CellType.FLUID || currentCell.cellType == Cell.CellType.FLUID)) ? newVelocity.x : currentCell.velocity.x);
            currentCell.tempVelocity.y = ((yMin && (yMin.cellType == Cell.CellType.FLUID || currentCell.cellType == Cell.CellType.FLUID)) ? newVelocity.y : currentCell.velocity.y);
            currentCell.tempVelocity.z = ((zMin && (zMin.cellType == Cell.CellType.FLUID || currentCell.cellType == Cell.CellType.FLUID)) ? newVelocity.z : currentCell.velocity.z);
        }

        commitVelocities();
    }

    // Get only velocity components of cell at given index which border fluid cells
    Vector3 borderFluidVelocityComponents(int i, int j, int k)
    {
        Cell cell = cells[i, j, k];
        Vector3 selectedVelocity = new Vector3(0f, 0f, 0f);

        if (xMinFluid(i, j, k)) { selectedVelocity.x = cell.velocity.x; }
        if (yMinFluid(i, j, k)) { selectedVelocity.y = cell.velocity.y; }
        if (zMinFluid(i, j, k)) { selectedVelocity.z = cell.velocity.z; }

        return selectedVelocity;
    }

    // Check whether cell in x velocity component is fluid
    bool xMinFluid(int i, int j, int k)
    {
        Cell xMin = cells[i - 1, j, k];
        return (xMin && xMin.cellType == Cell.CellType.FLUID);
    }

    // Check whether cell in y velocity component is fluid
    bool yMinFluid(int i, int j, int k)
    {
        Cell yMin = cells[i, j - 1, k];
        return (yMin && yMin.cellType == Cell.CellType.FLUID);
    }

    // Check whether cell in z velocity component is fluid
    bool zMinFluid(int i, int j, int k)
    {
        Cell zMin = cells[i, j, k - 1];
        return (zMin && zMin.cellType == Cell.CellType.FLUID);
    }


    // Calculate and apply pressure
    void pressure(float timeStep)
    {
        int cellCount = cells.Count;

        // Copy the list of cell keys 
        List<(int, int, int)> fluidCellKeyList = new List<(int, int, int)>();
        foreach ((int, int, int) key in cells.Keys)
        {
            Cell cell = cells[key];
            if (cell.cellType == Cell.CellType.FLUID)
            {
                fluidCellKeyList.Add(key);
            }
        }

        // Create a dictionary of every cell with an associated index for their values in the A matrix.
        Dictionary<Cell, int> cellIndices = new Dictionary<Cell, int>();
        int fluidCount = fluidCellKeyList.Count;
        for (int i = 0; i < fluidCount; i++)
        {
            (int, int, int) cellKey = fluidCellKeyList[i];
            Cell cell = cells[cellKey];
            cellIndices[cell] = i;
        }

        // Create a matrix A to store coefficients for every cell and their neighbours
        SparseRowMatrix A = new SparseRowMatrix(fluidCount, fluidCount);

        // Create a vector B to solve the A matrix for P later on.
        DenseVector B = new DenseVector(fluidCount);

        // TODO: “Ghost” pressure for solid walls?
        foreach ((int, int, int) key in fluidCellKeyList)
        {
            // Get the current key, cell, and cell index.
            (int i, int j, int k) = key;
            Cell currentCell = cells[key];

            // Only fluid cells are given a row.
            if (currentCell.cellType != Cell.CellType.FLUID) { continue; }

            int currentCellIndex = cellIndices[currentCell];        // Get corresponding row in matrix.
            Cell[] fluidNeighbours = {};                            // Neighbours that are fluid cells
            int nonSolidNeighCount = 0;                                  // Number of non-solid neighbours
            int fluidNeighCount = 0;

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
                nonSolidNeighCount++;

                // If a cell is a fluid cell, add 1 at their respective collumn in the A matrix.
                if (neighbour.cellType == Cell.CellType.FLUID)
                {
                    int fluidCellIndex = cellIndices[neighbour];
                    A.AddValue(currentCellIndex, fluidCellIndex, 1);
                    fluidNeighCount++;
                }
            }

            // In matrix A, set the current cell to the negative amount of non-solid neighbours.
            A.AddValue(currentCellIndex, currentCellIndex, -nonSolidNeighCount);

            // Calculate the divergence, setting to 0 for velocity components pointing into solid cells.
            float divergence = xMax.velocity[0] - (xMin.cellType != Cell.CellType.SOLID ? currentCell.velocity[0] : 0)
                                + yMax.velocity[1] - (yMin.cellType != Cell.CellType.SOLID ? currentCell.velocity[1] : 0)
                                + zMax.velocity[2] - (zMin.cellType != Cell.CellType.SOLID ? currentCell.velocity[2] : 0);

            // In matrix B, set a value based off the divergence in the velocity field for the current cell.
            int airNeighCount = nonSolidNeighCount - fluidNeighCount;
            float newB = currentCell.density() * cellSize * divergence / timeStep - airNeighCount * atmospheric_pressure;
            B.AddValue(currentCellIndex, newB);
        }

        // Solve for the actual pressure.
        // Create pressure vector that will be filled
        DateTime start = DateTime.UtcNow;
        DenseVector pressure = new DenseVector(fluidCount);

        // Initialize solver
        CGSolver solver = new CGSolver();
        DefaultLinearIteration iteration = (DefaultLinearIteration) solver.Iteration;
        iteration.SetParameters(relTolerance, absTolerance, 1e+5, 100000);
        solver.Iteration = iteration;

        // Initialize matrix that will be used in the preconditioner and initialize preconditioner
        SparseRowMatrix F = new SparseRowMatrix(fluidCount, fluidCount);
        ICCPreconditioner preconditioner = new ICCPreconditioner(F);
        preconditioner.Setup(A);

        // Couple preconditioner to solver
        solver.Preconditioner = preconditioner;

        // Solve for pressure
        pressure = (DenseVector) solver.Solve(A, B, pressure);
        DateTime end = DateTime.UtcNow;
        Debug.Log("time elapsed: " + (end - start));

        // Apply pressure to each cell (not to solid cells).
        foreach ((int, int, int) key in cells.Keys)
        {
            // Get the current key, cell, and cell index.
            (int i, int j, int k) = key;
            Cell currentCell = cells[key];

            // If celltype is solid, skip to next
            if(currentCell.cellType == Cell.CellType.SOLID) { continue; }

            // Get pressure of current cell.
            float cellPressure = (currentCell.cellType == Cell.CellType.FLUID ? (float) pressure.GetValue(cellIndices[currentCell]) : atmospheric_pressure);

            // Retrieve neighbours for the current cell.
            Cell xMin = cells[i - 1, j, k];
            Cell yMin = cells[i, j - 1, k];
            Cell zMin = cells[i, j, k - 1];

            // Calculate gradient of the pressure; only for velocity components bordering non-solid cells.
            Vector3 pressureGradient = new Vector3(0f,  0f, 0f);
            if (xMin && xMin.cellType != Cell.CellType.SOLID)
            {
                pressureGradient.x = cellPressure - (xMin.cellType == Cell.CellType.FLUID ? (float) pressure.GetValue(cellIndices[xMin]) : atmospheric_pressure);
            }
            if (yMin && yMin.cellType != Cell.CellType.SOLID)
            {
                pressureGradient.y = cellPressure - (yMin.cellType == Cell.CellType.FLUID ? (float) pressure.GetValue(cellIndices[yMin]) : atmospheric_pressure);
            }
            if (zMin && zMin.cellType != Cell.CellType.SOLID)
            {
                pressureGradient.z = cellPressure - (zMin.cellType == Cell.CellType.FLUID ? (float) pressure.GetValue(cellIndices[zMin]) : atmospheric_pressure);
            }

            // Apply to velocity
            currentCell.velocity -= timeStep * pressureGradient / (currentCell.density() * cellSize);
        }
    }

    // Extrapolate velocities to new cells
    void extrapolateVelocities() {
        // set layer field to 0 for fluid cells and −1 for non fluid cells
        foreach ((int, int, int) key in cells.Keys)
        {
            Cell currentCell = cells[key];
            currentCell.layer = ((currentCell.cellType == Cell.CellType.FLUID) ? 0 : -1);
        }

        for (int layer = 1; layer <= Mathf.Max(minBorderLayers, Mathf.Ceil(cflTimestepConstant)); layer++)
        {
            // Loop through non-fluid cells
            foreach ((int, int, int) key in cells.Keys)
            {
                (int i, int j, int k) = key;
                Cell currentCell = cells[key];
                if(currentCell.layer != -1) { continue; }   // Skip if the cell is a fluid cell

                // Retrieve all neighbours for the current cell.
                Cell xMax = cells[i + 1, j, k];
                Cell yMax = cells[i, j + 1, k];
                Cell zMax = cells[i, j, k + 1];
                Cell xMin = cells[i - 1, j, k];
                Cell yMin = cells[i, j - 1, k];
                Cell zMin = cells[i, j, k - 1];

                Cell[] neighbours = { xMax, xMin, yMax, yMin, zMax, zMin };

                // Check whether cell has a neighbour s.t. neighbour.layer == -1
                List<Cell> visitedNeighbours = new List<Cell>();
                foreach (Cell neighbour in neighbours)
                {
                    if(neighbour && neighbour.layer == layer - 1)
                    {
                        visitedNeighbours.Add(neighbour);
                    }
                }

                // If no such neighbour exists skip to next
                int numberOfVisited = visitedNeighbours.Count;
                if(numberOfVisited == 0) { continue; }

                // Average velocity of visitedNeighbours
                Vector3 totalVelocity = new Vector3(0, 0, 0);
                foreach(Cell neighbour in visitedNeighbours)
                {
                    totalVelocity += neighbour.velocity;
                }
                Vector3 averageVelocity = totalVelocity / numberOfVisited;

                // For velocity components of cell not bordering fluid cells set uj to the average 
                // of the neighbors of Cin which N.layer ==i−1
                if (!xMin || xMin.cellType != Cell.CellType.FLUID) {
                    currentCell.velocity.x = averageVelocity.x;
                }
                if (!yMin || yMin.cellType != Cell.CellType.FLUID)
                {
                    currentCell.velocity.y = averageVelocity.y;
                }
                if (!zMin || zMin.cellType != Cell.CellType.FLUID)
                {
                    currentCell.velocity.z = averageVelocity.z;
                }

                currentCell.layer = layer;
            }
        }
    }

    // Set velocities of solid cells
    void setSolidCells()
    {
        // Set velocity components that point into solid cells from liquid or air cells to zero.
        foreach ((int, int, int) key in cells.Keys)
        {
            // Get cell
            (int i, int j, int k) = key;
            Cell currentCell = cells[key];

            // Get neighbours into which velocity components point
            Cell xMin = cells[i - 1, j, k];
            Cell yMin = cells[i, j - 1, k];
            Cell zMin = cells[i, j, k - 1];

            // Set velocity component to zero if it points into solid neighbour.
            if(currentCell.cellType != Cell.CellType.SOLID)
            {
                if (xMin && (xMin.cellType == Cell.CellType.SOLID && currentCell.velocity.x < 0)) { currentCell.velocity.x = 0; }
                if (yMin && (yMin.cellType == Cell.CellType.SOLID && currentCell.velocity.y < 0)) { currentCell.velocity.y = 0; }
                if (zMin && (zMin.cellType == Cell.CellType.SOLID && currentCell.velocity.z < 0)) { currentCell.velocity.z = 0; }
            }
            else
            {
                if (xMin && (xMin.cellType != Cell.CellType.SOLID && currentCell.velocity.x > 0)) { currentCell.velocity.x = 0; }
                if (yMin && (yMin.cellType != Cell.CellType.SOLID && currentCell.velocity.y > 0)) { currentCell.velocity.y = 0; }
                if (zMin && (zMin.cellType != Cell.CellType.SOLID && currentCell.velocity.z > 0)) { currentCell.velocity.z = 0; }
            }
        }
    }
}
