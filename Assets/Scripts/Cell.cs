using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Cell
{

    public enum CellType { FLUID, AIR, SOLID }

    public CellType cellType;
    public int layer;
    public Vector3 location;
    public Vector3 velocity;
    public Vector3 tempVelocity;

    public Cell(CellType _cellType, int _layer, Vector3 _location)
    {
        cellType = _cellType;
        layer = _layer;
        location = _location;
        velocity = new Vector3(0f, 0f, 0f);
    }

    public static implicit operator bool(Cell cell) {
        return cell != null;
    }

    public float density()              // TODO: does density depend on number of particles?
    {
        if (cellType == CellType.AIR)
        {
            return CellMaster.AIR_DENSITY;
        } 
        else if (cellType == CellType.FLUID)
        {
            return CellMaster.FLUID_DENSITY;
        }
        else
        {
            return 99999f;
        }
    }
}