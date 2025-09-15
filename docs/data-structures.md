## Data structures
### Vectors and lists
Gmml2 makes frequent use of sequential data. The basic sequential data structures in C++ are the `std::vector` for variable length, and `std::array` for fixed length. This documentation makes use of notation like [1, 2, 3], which we will call a list (not to be confused with `std::list`). In the C++ context, it can represent either `std::vector<T>{1, 2, 3}` or `std::array<T, 3>{1, 2, 3}`. The common property is that these data structures can be accessed in constant-time given the index of an entry.
### Tabular data
Gmml2 was originally written in a style linking data together with pointers. It is in the process of being rewritten in a style where data is organized in tabular form, which enables better performance and is easier to copy. A common way of representing a table is the struct of vectors. It has a similar layout to a database table, with each vector representing a column.

For instance, we could define
```
struct AtomData {
  std::vector<std::string> names;
  std::vector<uint> numbers
};
```

Here, `names[n]` would be the name of the atom with index `n`, and `numbers[n]` would be the number of the same atom. When an entry is added to such a struct, it is important to push data into each of the vectors, in order to keep them the same size. For this reason it is recommended to let only a single function bear this responsibility, per data structure.

It is ill-adviced to delete entries from the middle of vectors, since this risks invalidating indices in use and will likely result in segmentation faults down the line. Prefer using a `std::vector<bool>` to flag entries as removed instead.

### Old style
The old style was called CentralDataStructure, or cds for short. Some parts of this code remain, and have not kept up to naming standards as they are planned to be rewritten. These are:

* Logic for finding rotatable bonds in residue linkages
* Carbohydrate creation
* Glycoprotein creation (mainly due to relying on the residue linkage logic)

There are also a few interface functions which convert from the pointer style data to the tabular style. These should disappear when the rewrite is finished.
