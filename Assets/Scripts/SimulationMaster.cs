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

    // Start is called before the first frame update
    void Start()
    {
        setFixedTimeStep(0f);
        // Create cell master, responsible for updating grid and its velocities
        // cellMaster = new CellMaster(startLocation, x_size, y_size, z_size, cellSize, viscosity, atmospheric_pressure);
        
        // Get all particle objects
        // int numberOfParticles = particleObjects.Length;

        // Get Particle component from particle objects
        // particles = new Particle[numberOfParticles];
        // for (int i = 0; i < numberOfParticles; i++)
        // {
        //     particles[i] = particleObjects[i].GetComponent<Particle>();
        // }
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
    void Awake() {
        SpawnParticles();
    }

    void SpawnParticles() {
        // if ( particles.Count > 0 ) {
        //     foreach(Particle particle in particles) {
        //         GameObject.Destroy(particle.gameObject);
        //     }
        //     particles.Clear();
        // }

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
                    GameObject particleObject = Instantiate(particlePrefab);
                    Particle particle = new Particle(particleObject);
                    particleObject.transform.position = new Vector3(startLocation.x + x * space_between,
                        startLocation.y + y * space_between,
                        startLocation.z + z * space_between);

                    particles.Add(particle);
                }
            }
        }
        // throw new System.Exception("STOP");
    }

    void setFixedTimeStep(float maxVelocity)
    {
        float deltaTime = cellMaster.cflTimestepConstant * cellMaster.cellSize / maxVelocity;
        Time.fixedDeltaTime = Mathf.Min(deltaTime, cellMaster.minDeltaTime);
    }

    void FixedUpdate()
    {
        // Get timeStep for this update
        float timeStep = Time.fixedDeltaTime;


        //2. Update the grid based on the marker particles (figure 4)
        DateTime start = DateTime.UtcNow;
        cellMaster.updateGrid(particles);
        DateTime end = DateTime.UtcNow;
        TimeSpan timeDiff = end - start;
        Debug.Log("Update grid time elapsed " + timeDiff.ToString());

        // 3. Advance the velocity field, u
        start = DateTime.UtcNow;
        cellMaster.velocityUpdate(timeStep);
        end = DateTime.UtcNow;
        timeDiff = end - start;
        Debug.Log("Update velocity time elapsed " + timeDiff.ToString());

        // 4. Move the particles through u for ∆t time
        // TODO:  Since∆tdoes not necessarily coincide with frameboundaries, the particles should not always be advanced by∆t,however.
        start = DateTime.UtcNow;
        foreach (Particle particle in particles) {
            Vector3 velocityParticle = cellMaster.getParticleVelocity(particle.getPosition());
            particle.locationUpdate(timeStep, velocityParticle);
        }
        end = DateTime.UtcNow;
        timeDiff = end - start;
        Debug.Log("Update particles time elapsed " + timeDiff.ToString());

        float currentMaxVelocityMag = cellMaster.getMaxVelocityMagnitude();
        setFixedTimeStep(currentMaxVelocityMag);

    }
}
