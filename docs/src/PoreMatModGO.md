# PoreMatModGO

To make `PoreMatMod.jl` more accessible, an express graphical interface, `PoreMatModGO`, is provided.  
The `PoreMatModGO` interface, built on `Pluto.jl`, enables interactive application of chemical substructure find/replace operations with minimal setup and no code provided by the user.

Follow the steps in [Getting Started](../manual/start).

Then, simply launch the GUI notebook via `Pluto`:

```julia
using PoreMatMod
PoreMatModGO()
```

The notebook may take several minutes to launch the first time, especially on Windows.

Prepare file inputs per the [manual](../manual/inputs) and load them into `PoreMatModGO` using the graphical interface.  
All [replacement](../manual/replace) modes are available with interactive visual previews and outputs are downloadable in `.cif`, `.xyz`, and `.vtk` file formats.
