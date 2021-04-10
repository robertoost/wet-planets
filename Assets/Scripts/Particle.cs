using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Particle
{
    public CellMaster cellMaster;
    public GameObject gameObject;
    public Particle(GameObject gameObject) { this.gameObject = gameObject; }

    public void locationUpdate(float timeStep, Vector3 velocity)
    {
        gameObject.transform.position += velocity * timeStep;
    }

    public Vector3 getPosition()
    {
        return gameObject.transform.position;
    }
}
