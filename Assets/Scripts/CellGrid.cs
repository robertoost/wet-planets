using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class  CellGrid
{
    Hashtable cells = new Hashtable();
    public Cell this[int x, int y, int z]
    {
        get { return (Cell) cells[(x, y, z)]; }
        set { cells[(x, y, z)] = value; }
    }

    //public Cell this[int i]
    //{
    //    get { return (Cell) cells[cells.Keys[i]]; }
    //    set { cells[cells.Keys[i]] = value; }
    //}

    //private void hashFunction(int x, int y, int z)
    //{
    //    return 541 * x + 79 * y + 31 * z;
    //}

    //private void invHashFunction(int hash)
    //{
    //    int rest1 = hash % 541;
    //    int x = floor(hash / 541);
    //    int z = rest1 % 79;
    //    int y = floor(rest1 / 79);
    //}

    //private Cell getCellVelocity(int x, int y, int z, int index)
    //{
    //    Cell cell = (Cell) cells[(x, y, z)];
    //    return (cell ? cell.velocity[index] : 0);
    //}

    public void addCell(Cell cell, int x, int y, int z)
    {
        cells.Add((x, y, z), cell);
    }

    public void removeCell(int x, int y, int z)
    {
        cells.Remove((x, y, z));
    }

    public ICollection getKeys()
    {
        return cells.Keys;
    }
}
