using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

[CustomEditor(typeof(CelestialObject))]
public class CelestialObjectEditor : Editor
{
    CelestialObject celestialObject;
    bool showDebugInfo;
    public override void OnInspectorGUI () {
        DrawDefaultInspector ();

        EditorGUILayout.Space (10);
        EditorGUILayout.LabelField ("Debug", EditorStyles.boldLabel);
        showDebugInfo = EditorGUILayout.Foldout (showDebugInfo, "Debug info");
        if (showDebugInfo) {
            string[] gravityInfo = GetGravityInfo (celestialObject.transform.position, celestialObject);
            for (int i = 0; i < gravityInfo.Length; i++) {
                EditorGUILayout.LabelField (gravityInfo[i]);
            }
        }
    }
    void OnEnable () {
        celestialObject = (CelestialObject) target;
        showDebugInfo = EditorPrefs.GetBool (celestialObject.gameObject.name + nameof (showDebugInfo), false);
    }

    void OnDisable () {
        if (celestialObject) {
            EditorPrefs.SetBool (celestialObject.gameObject.name + nameof (showDebugInfo), showDebugInfo);
        }
    }

    static string[] GetGravityInfo (Vector3 point, CelestialObject ignore = null) {
        CelestialObject[] bodies = GameObject.FindObjectsOfType<CelestialObject> ();
        Vector3 totalAcc = Vector3.zero;

        // gravity
        var forceAndName = new List<FloatAndString> ();
        foreach (CelestialObject body in bodies) {
            if (body != ignore) {
                var offsetToBody = body.gameObject.transform.position - point;
                var sqrDst = offsetToBody.sqrMagnitude;
                float dst = Mathf.Sqrt (sqrDst);
                var dirToBody = offsetToBody / Mathf.Sqrt (sqrDst);
                var acceleration = CelestialObject.GRAV_CONST * body.gameObject.GetComponent<Rigidbody>().mass / sqrDst;
                totalAcc += dirToBody * acceleration;
                forceAndName.Add (new FloatAndString () { floatVal = acceleration, stringVal = body.gameObject.name });

            }
        }
        forceAndName.Sort ((a, b) => (b.floatVal.CompareTo (a.floatVal)));
        string[] info = new string[forceAndName.Count + 1];
        info[0] = $"acc: {totalAcc} (mag = {totalAcc.magnitude})";
        for (int i = 0; i < forceAndName.Count; i++) {
            info[i + 1] = $"acceleration due to {forceAndName[i].stringVal}: {forceAndName[i].floatVal}";
        }
        return info;
    }

    struct FloatAndString {
        public float floatVal;
        public string stringVal;
    }
}
