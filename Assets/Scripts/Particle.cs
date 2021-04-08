using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Particle : MonoBehaviour
{
    public CellMaster cellMaster;

    public Particle() { }

    public void locationUpdate(float timeStep, Vector3 velocity)
    {
        transform.position += velocity * timeStep;
    }

    public Vector3 getPosition()
    {
        return transform.position;
    }
}
