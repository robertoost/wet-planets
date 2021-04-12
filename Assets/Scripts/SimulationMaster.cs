using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SimulationMaster : MonoBehaviour
{
    // User variables
    // public string particleTag;          // Tag used for particle GameObjects
    public Vector3 startLocation;       // Minimum location of grid
    public Vector3 particleFieldSize = new Vector3(2, 0.5f, 2);
    public float space_between;
    public GameObject particlePrefab;
    private List<Particle> particles = new List<Particle>();
    public CellMaster cellMaster;

    // Spawn particles at the very earliest moment
    void Awake()
    {
        SpawnParticles();
    }

    // Start is called before the first frame update
    void Start()
    {
        setFixedTimeStep(0f);       // Set variable time step for fixed update
        //setStartLocation();         // Shift start location to account for discretization
    }

    // public void DrawGizmos() {
    //     foreach(Cell cell in cells ) {
    //         bool isFluid = cell.cellType == Cell.CellType.FLUID;
    //         if ( !isFluid) {
    //             continue;
    //         }

    //         Gizmos.color = Color.cyan;
            
    //         Gizmos.DrawWireCube(cell.location, (Vector3.one * cellSize));
    //     }
    // }

    // private void OnDrawGizmosSelected() {
    //     foreach(Cell cell in cellMaster.cells ) {
    //         bool isFluid = cell.cellType == Cell.CellType.FLUID;
    //         if ( !isFluid) {
    //             continue;
    //         }

    //         Gizmos.color = Color.cyan;
            
    //         Gizmos.DrawWireCube(cell.location, (Vector3.one * cellMaster.cellSize));
    //     }
    // }

    void SpawnParticles() {

        float x_size = particleFieldSize.x;
        float y_size = particleFieldSize.y;
        float z_size = particleFieldSize.z;

        int numberParticlesX = (int) (x_size / space_between);
        int numberParticlesY = (int) (y_size / space_between);
        int numberParticlesZ = (int) (z_size / space_between);
        for (int x = 0; x < numberParticlesX; x++)
        {
            for (int y = 0; y < numberParticlesY; y++)
            {
                for(int z = 0; z< numberParticlesZ; z++)
                {
                    GameObject particleObject = Instantiate(particlePrefab, transform);
                    Particle particle = new Particle(particleObject);
                    particleObject.transform.position = new Vector3(startLocation.x + x * space_between,
                        startLocation.y + y * space_between,
                        startLocation.z + z * space_between);

                    particles.Add(particle);
                }
            }
        }
    }

    // Change time step based on maximum velocity of the particles
    void setFixedTimeStep(float maxVelocity)
    {
        float deltaTime = cellMaster.cflTimestepConstant * cellMaster.cellSize / maxVelocity;
        Time.fixedDeltaTime = Mathf.Min(deltaTime, cellMaster.minDeltaTime);
    }

    // Set grid at right position, has to be a bit shifted to account for discretization
    void setStartLocation()
    {
        cellMaster.startLocation = cellMaster.startLocation
            - new Vector3(0.5f * cellMaster.cellSize, 0.5f * cellMaster.cellSize, 0.5f * cellMaster.cellSize);
    }

    void FixedUpdate()
    {
        // Get timeStep for this update
        float timeStep = Time.fixedDeltaTime;

        //2. Update the grid based on the marker particles (figure 4)
        cellMaster.updateGrid(particles);

        // 3. Advance the velocity field, u
        cellMaster.velocityUpdate(timeStep);

        // 4. Move the particles through u for ∆t time
        // TODO:  Since∆tdoes not necessarily coincide with frameboundaries, the particles should not always be advanced by∆t,however.
        foreach (Particle particle in particles) {
            Vector3 velocityParticle = cellMaster.getParticleVelocity(particle.getPosition());
            particle.locationUpdate(timeStep, velocityParticle);
        }

        // 1. Change time step based on maximum velocity of particles
        float currentMaxVelocityMag = cellMaster.getMaxVelocityMagnitude();
        setFixedTimeStep(currentMaxVelocityMag);
        ;

    }
}
