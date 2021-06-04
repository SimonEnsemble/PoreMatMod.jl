# MOFunGO

To make `MOFun.jl` more accessible to researchers without strong programming
backgrounds, an express graphical interface, `MOFunGO`, is provided.  The
`MOFunGO` interface, built on `Pluto.jl`, enables interactive application of
chemical substructure find/replace operations with minimal setup and no code
provided by the user.

Follow the steps in [Getting Started](../start).

Add the required libraries and launch `Pluto`:

```
pkg] add Pluto, PlutoUI, Bio3DView
julia> using Pluto
julia> Pluto.run()
```

`Pluto` should launch in a new tab in your default web browser.

Load `MOFunGO.jl` by copying and pasting the URL below into the "Open from file"
box on the `Pluto` main page.

[https://github.com/eahenle/MOFunGO/blob/main/MOFunGO.jl](https://github.com/eahenle/MOFunGO/blob/main/MOFunGO.jl)

Prepare file inputs per the [manual](../inputs) and load them into
`MOFunGO` using the graphical interface.  All [replacement](../replace) modes
are available with interactive visual previews and downloadable outputs in `.cif`,
`.xyz`, and `.vtk` file formats.
