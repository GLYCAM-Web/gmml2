## Swig wrapping

Swig is a tool which can bridge C++ code to other languages such as Python.
Most linux distros have swig 4.0. Otherwise swig 4.0.2 must be installed from [their website](https://www.swig.org/download.html).

Swig wrapping adds the following dependencies

* `python3.9` (Version `3.9.12`)
* `python3.9-dev` (Version `3.9.12`)
* `swig` (Version `4.0.2`)

In order to wrap with Swig, use the `-w` flag of the make script, e.g `./make.sh -w`
