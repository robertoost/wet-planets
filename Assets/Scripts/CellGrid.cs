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

    // // Iterate over all cells by their indices.
    // public IEnumerator<(int,int,int)> GetEnumerator()
    // {
    //     foreach ((int,int,int) key in cells.Keys) {
    //         yield return key;
    //     }
    // }

    IEnumerator IEnumerable.GetEnumerator()
    {
        throw new NotImplementedException();
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


}
