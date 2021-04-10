using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SpawnParticles : MonoBehaviour
{
    public GameObject exampleParticle;
    public Vector3 startLocation;
    public float x_size;
    public float y_size;
    public float z_size;
    public float space_between;

    // Start is called before the first frame update
    void Awake()
    {
        int numberParticlesX = (int) (x_size / space_between);
        int numberParticlesY = (int) (y_size / space_between);
        int numberParticlesZ = (int) (z_size / space_between);
        for (int x = 0; x < numberParticlesX; x++)
        {
            for (int y = 0; y < numberParticlesY; y++)
            {
                for(int z = 0; z< numberParticlesZ; z++)
                {
                    if(x == 0 && y == 0 && z == 0)
                    {
                        exampleParticle.transform.position = startLocation;
                    }
                    else
                    {
                        GameObject newParticle = Instantiate(exampleParticle);
                        newParticle.transform.position = new Vector3(startLocation.x + x * space_between,
                            startLocation.y + y * space_between,
                            startLocation.z + z * space_between);
                    }
                }
            }
        }
    }
}
