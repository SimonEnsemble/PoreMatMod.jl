# PoreMatModGO

To make `PoreMatMod.jl` more accessible, an express graphical interface, `PoreMatModGO`, is provided.  
The `PoreMatModGO` interface, built on `Pluto.jl`, enables interactive application of chemical substructure find/replace operations with minimal setup and no code provided by the user.

Follow the steps in [Getting Started](../start).

Then, add the required dependencies and launch `Pluto.jl`:

```
pkg> add Pluto, PlutoUI, Bio3DView
julia> using Pluto
julia> Pluto.run()
```

`Pluto` should launch in a new tab in your default web browser.

Load `PoreMatModGO.jl` by copying and pasting the URL below into the "Open from file"
box on the `Pluto` main page.

[https://raw.githubusercontent.com/SimonEnsemble/PoreMatMod.jl/master/src/PoreMatModGO.jl](https://raw.githubusercontent.com/SimonEnsemble/PoreMatMod.jl/master/src/PoreMatModGO.jl?token=AD3TMGCN5XHK3VKRRIFAHV3BGBIE4)

Prepare file inputs per the [manual](../inputs) and load them into `PoreMatModGO` using the graphical interface.  
All [replacement](../replace) modes are available with interactive visual previews and outputs are downloadable in `.cif`, `.xyz`, and `.vtk` file formats.
