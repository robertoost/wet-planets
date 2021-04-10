using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SimulationMaster : MonoBehaviour
{
    // User variables
    public string particleTag;          // Tag used for particle GameObjects
    public Vector3 startLocation;       // Minimum location of grid
    public int x_size;                  // Width size of grid
    public int y_size;                  // Height size of grid
    public int z_size;                  // Depth size of grid
    public float cellSize;              // Size of cells within grid

    // Constants
    public float viscosity = 1.0016f;
    public float atmospheric_pressure = 101.325f;

    private CellMaster cellMaster;
    private Particle[] particles;
    private int numberOfParticles;


    // Start is called before the first frame update
    void Start()
    {
        // Create cell master, responsible for updating grid and its velocities
        cellMaster = new CellMaster(startLocation, x_size, y_size, z_size, cellSize, viscosity, atmospheric_pressure);
        
        // Get all particle objects
        GameObject[] particleObjects = GameObject.FindGameObjectsWithTag(particleTag);
        numberOfParticles = particleObjects.Length;

        // Get Particle component from particle objects
        particles = new Particle[numberOfParticles];
        for (int i = 0; i < numberOfParticles; i++)
        {
            particles[i] = particleObjects[i].GetComponent<Particle>();
        }
    }

    void FixedUpdate()
    {
        // Get timeStep for this update
        float timeStep = Time.fixedDeltaTime;

        //2. Update the grid based on the marker particles (figure 4)
        cellMaster.updateGrid(particles);

        // 3. Advance the velocity field, u
        cellMaster.velocityUpdate(timeStep);

        // 4. Move the particles through u for ∆t time.
        for (int i = 0; i < particles.Length; i++)
        {
            Particle particle = particles[i];
            Vector3 velocityParticle = cellMaster.getParticleVelocity(particle.getPosition());
            particles[i].locationUpdate(timeStep, velocityParticle);
        }
    }
}
