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

    private CellMaster cellMaster;
    private Particle[] particles;
    private int numberOfParticles;


    // Start is called before the first frame update
    void Start()
    {
        cellMaster = new CellMaster(startLocation, x_size, y_size, z_size, cellSize);
        
        GameObject[] particleObjects = GameObject.FindGameObjectsWithTag(particleTag);
        numberOfParticles = particleObjects.Length;

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

        // 4. Move the particles throughufor ∆t time.
        for (int i = 0; i < particles.Length; i++)
        {
            Particle particle = particles[i];
            Debug.Log("particle loc " + particle.getPosition().ToString());
            Vector3 velocityParticle = cellMaster.getVelocity(particle.getPosition());
            particles[i].locationUpdate(timeStep, velocityParticle);
        }
    }
}
