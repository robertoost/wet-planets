using System.Collections;
using System.Collections.Generic;
using UnityEngine;


public class CelestialObject : MonoBehaviour
{
    const float GRAV_CONST = 0.0001f;
    bool isFixed;
    public Rigidbody rb;
    float radius;
    float mass;
    Vector3 initialVelocity;
    public static List<CelestialObject> allCelestialObjects = new List<CelestialObject>();

    // On Awake, before any Start is executed, add this object to the list of all celestial objects.
    void Awake() {
        allCelestialObjects.Add(this);
    }

    // Start is called before the first frame update
    void Start()
    {
        rb = GetComponent<Rigidbody>();
        rb.velocity = initialVelocity;
    }

    void OnDestroy() {
        allCelestialObjects.Remove(this);
    }

    public void UpdateVelocity(float timeStep) {
        Vector3 totalAcceleration = CalculateTotalAcceleration();
        rb.velocity += totalAcceleration * timeStep;
    }

    public Vector3 CalculateAcceleration(CelestialObject other) {
        if (other == this) {
            return Vector3.zero;
        }

        // Calculate the acceleration using the other's difference as a force direction and square distance.
        Vector3 difference = other.rb.position - rb.position;
        Vector3 acceleration = difference.normalized * GRAV_CONST * other.mass / difference.sqrMagnitude;

        return acceleration;
    }

    public Vector3 CalculateTotalAcceleration() {
        Vector3 acceleration = Vector3.zero;

        // Calculate acceleration of this object using every other celestial object.
        foreach (CelestialObject other  in CelestialObject.allCelestialObjects) {
            
            // Ignore self.
            if (other == this) {
                continue;
            }

            // Add all acceleration values together to get the total acceleration for this object.
            acceleration += CalculateAcceleration(other);
        }
        
        return acceleration;
    }
}
