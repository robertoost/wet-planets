using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;

public class CellMaster : MonoBehaviour
{
    public int x_size;
    public int y_size;
    public int z_size;
    public float cellSize;

    //private List<Cell> cells;
    CellGrid cells;

    private const float GRAV = -9.81f;
    public const float VISCOSITY = 0.89f;
    public const float FLUID_DENSITY = 1.00f;
    public const float ATMOSPHERIC_PRESSURE = 101325f;

    // Start is called before the first frame update
    void Start()
    {
        for(int x = 0; x < x_size; x++)
        {
            for(int y = 0; y < y_size; y++)
            {
                for(int z = 0; z < z_size; z++)
                {
                    // Initialize cell properties
                    Cell.CellType cellType = Cell.CellType.FLUID;
                    Vector3 velocity = new Vector3(0, 0, 0);
                    float pressure = -1;

                    // Initialize location
                    float x_loc = (float)(0.5 + x) * cellSize;
                    float y_loc = (float)(0.5 + y) * cellSize;
                    float z_loc = (float)(0.5 + z) * cellSize;
                    Vector3 location = new Vector3(x_loc, y_loc, z_loc);

                    // Create new cell
                    Cell newCell = new Cell(cellType, location, velocity, pressure, FLUID_DENSITY);
                    //cells.Add(newCell);
                    cells.addCell(newCell, x, y, z);
                }
            }
        }
    }

    // Update is called once per frame
    void FixedUpdate()
    {
        float timeStep = Time.fixedDeltaTime;

        //3a. Apply convection using a backwards particle trace.
        convection(timeStep);

        //3b. Apply external forces.
        externalForces(timeStep);

        //3c. Apply viscosity.
        viscosity(timeStep);

        //3d. Calculate the pressure to satisfy ∇·u=0.
        pressure(timeStep);

        //3e. Apply the pressure.

        //3f. Extrapolate fluid velocities into buffer zone.

        //3g. Set solid cell velocities.
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
        Matrix<float> A = SparseMatrixModule.zero<float>(cellCount, cellCount);

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
            B[currentCellIndex] = currentCell.density * cellSize * divergence / timeStep - airCount * ATMOSPHERIC_PRESSURE;
        }

        // Solve for the actual pressure.
        Vector<float> pressure = (Vector<float>) A.QR().Solve(B);

        
    }
}
