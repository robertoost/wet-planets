using UnityEngine;
public class CelestialSimulation : MonoBehaviour {
    void FixedUpdate() {
        foreach(CelestialObject celestialObject in CelestialObject.allCelestialObjects) {
            celestialObject.UpdateVelocity(Time.fixedDeltaTime);
        }
    }
}