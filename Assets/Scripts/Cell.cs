using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Cell
{

    public enum CellType { FLUID, AIR, SOLID }

    public CellType cellType;
    public Vector3 location;
    public Vector3 velocity;
    public Vector3 tempVelocity;
    public float pressure;
    public float density;

    public Cell(CellType _cellType, Vector3 _location, Vector3 _velocity, float _pressure, float _density)
    {
        cellType = _cellType;
        location = _location;
        velocity = _velocity;
        pressure = _pressure;
        density = _density;
    }

    public static implicit operator bool(Cell cell) {
        return cell != null;
    }
}