using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CellMaster : MonoBehaviour
{
    public int x_size;
    public int y_size;
    public int z_size;
    public float cellSize;

    //private List<Cell> cells;
    CellGrid cells;

    private const float GRAV = -9.81f;

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
                    Cell newCell = new Cell(cellType, location, velocity, pressure);
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
        findPressure(timeStep);

        //3e. Apply the pressure.

        //3f. Extrapolate fluid velocities into buffer zone.

        //3g. Set solid cell velocities.
    }

    void convection(float timeStep)
    {
        foreach ((int x, int y, int z) in cells.getKeys())
        {
            Cell currentCell = cells[x, y, z];
            // Extraction location of cell
            locationCell = currentCell.location;
            x = locationCell.x; y = locationCell.y; z = locationCell.z;
            
            // Find velocity
            Vector3 V = getVelocity(x, y, z);
            V = getVelocity(x + 0.5 * timeStep * V.x, y + 0.5 * timeStep * V.y, z + 0.5 * timeStep * V.z);

            // Update velocity of cell
            currentCell.tempVelocity = V;
        }

        // Commit velocities
        foreach ((int x, int y, int z) in cells.getKeys())
        {
            Cell currentCell = cells[x, y, z];
            currentCell.velocity = currentCell.tempVelocity;
        }
    }

    // Get the interpolated velocity at a point in space.
    Vector3 getVelocity(float x, float y, float z)
    {
        Vector3 V;
        V.x = getInterpolatedValue(x / cellSize, y / cellSize - 0.5, z / cellSize - 0.5, 0);
        V.y = getInterpolatedValue(x / cellSize - 0.5, y / cellSize, z / cellSize - 0.5, 1);
        V.z = getInterpolatedValue(x / cellSize - 0.5, y / cellSize - 0.5, z / cellSize, 2);
        return V;
    }

    // Get an interpolated data value from the grid.
    //TODO: Some cells might not exist and need to be skipped, see paper Figure 5
    float getInterpolatedValue(float x, float y, float z, int index)
    {
        int i = floor(x);
        int j = floor(y);
        int k = floor(z);

        Cell[] cellArray = {cells[i, j, k], cells[i + 1, j, k], cells[i, j + 1, k], cells[i + 1, j + 1, k],
                                  cells[i, j, k + 1], cells[i + 1, j, k + 1], cells[i, j + 1, k + 1], cells[i + 1, j + 1, k + 1] };

        float velocitySum = (i + 1 - x) * (j + 1 - y) * (k + 1 - z) * (cellArray[0]? cellArray[0].velocity[index] : 0)
            + (x - i) * (j + 1 - y) * (k + 1 - z) * (cellArray[1] ? cellArray[1].velocity[index] : 0)
            + (i + 1 - x) * (y - j) * (k + 1 - z) * (cellArray[2] ? cellArray[2].velocity[index] : 0)
            + (x - i) * (y - j) * (k + 1 - z) * (cellArray[3] ? cellArray[3].velocity[index] : 0)
            + (i + 1 - x) * (j + 1 - y) * (z - k) * (cellArray[4] ? cellArray[4].velocity[index] : 0)
            + (x - i) * (j + 1 - y) * (z - k) * (cellArray[5] ? cellArray[5].velocity[index] : 0)
            + (i + 1 - x) * (y - j) * (z - k) * (cellArray[6] ? cellArray[6].velocity[index] : 0)
            + (x - i) * (y - j) * (z - k) * (cellArray[7] ? cellArray[7].velocity[index] : 0);

        float weightSum = (i + 1 - x) * (j + 1 - y) * (k + 1 - z) * (cellArray[0] ? 1 : 0)
            + (x - i) * (j + 1 - y) * (k + 1 - z) * (cellArray[1] ? 1 : 0)
            + (i + 1 - x) * (y - j) * (k + 1 - z) * (cellArray[2] ? 1 : 0)
            + (x - i) * (y - j) * (k + 1 - z) * (cellArray[3] ? 1 : 0)
            + (i + 1 - x) * (j + 1 - y) * (z - k) * (cellArray[4] ? 1 : 0)
            + (x - i) * (j + 1 - y) * (z - k) * (cellArray[5] ? 1 : 0)
            + (i + 1 - x) * (y - j) * (z - k) * (cellArray[6] ? 1 : 0)
            + (x - i) * (y - j) * (z - k) * (cellArray[7] ? 1 : 0);

        return velocitySum / weightSum;
    }

    void externalForces(float timeStep)
    {
        foreach ((int x, int y, int z) in cells.getKeys())
        {
            cells[x, y, z].velocity += new Vector3(0, GRAV * timeStep, 0);
        }
    }

    void viscosity(float timeStep)
    {
    }

    void findPressure(float timeStep)
    {
    }
}
