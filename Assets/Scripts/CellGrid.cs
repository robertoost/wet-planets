using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CellGrid : IEnumerable
{
    Dictionary<(int, int,int), Cell> cells = new Dictionary<(int, int,int), Cell>();

    // Indexer for a triple key.
    public Cell this[(int, int, int) key]
    {
        get {
            Cell cell;
            cells.TryGetValue(key, out cell);
            return cell;
        }
        set {
            cells[(key)] = value; 
        }
    }

    // Indexer for a set of indices.
    public Cell this[int x, int y, int z]
    {
        get { return this[(x, y, z)]; }
        set { this[(x, y, z)] = value; }
    }
    public int Count {
        get { return cells.Count; }
    }
    public ICollection Keys
    {
        get { return cells.Keys; }
    }

    //Loop over each cell.
    public IEnumerator<Cell> GetEnumerator()
    {
        foreach((int,int,int) key in cells.Keys)
        {
            yield return cells[key];
        }
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        throw new NotImplementedException();
    }

    public void addCell(Cell cell, int x, int y, int z)
    {
        cells.Add((x, y, z), cell);
    }

    public void removeCell((int, int, int) key)
    {
        cells.Remove(key);
    }
}
