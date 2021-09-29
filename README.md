![logo.JPG](logo.jpg)
---

`PoreMatMod.jl` is a [Julia](https://julialang.org/) package for (i) subgraph matching and (ii) modifying crystal structures such as metal-organic frameworks (MOFs).

Functioning as a "find-and-replace" tool on atomistic structure models of porous materials, `PoreMatMod.jl` is useful for:
:hammer: subgraph matching to filter databases of crystal structures

:hammer: create hypothetical crystal structure models of functionalized structures

:hammer: introduce defects into crystal structures

:hammer: repair artifacts of X-ray structure determination, such as missing hydrogen atoms, disorder, and guest molecules

See documentation [here](https://SimonEnsemble.github.io/PoreMatMod.jl/) and examples (in [Pluto](https://github.com/fonsp/Pluto.jl) notebooks) [here](https://github.com/SimonEnsemble/PoreMatMod.jl/tree/master/examples).

N.b. while `PoreMatMod.jl` was developed for MOFs and other porous crystalline materials, its find-and-replace operations can be applied to discrete, molecular structures as well by assigning an arbitrary unit cell.

| **Documentation** | **Build Status** | 
|:---:|:---:|
| [![Docs Badge](https://img.shields.io/badge/docs-latest-blue.svg)](https://SimonEnsemble.github.io/PoreMatMod.jl/) | [![Build Status](https://travis-ci.org/SimonEnsemble/PoreMatMod.jl.svg?branch=master)](https://app.travis-ci.com/github/SimonEnsemble/PoreMatMod.jl) |


## License
`PoreMatMod.jl` is licensed under the [MIT license](./LICENSE).

## Installation
`PoreMatMod.jl` is a registered pakcage and can be installed by entering the following the in package manager.

```
(v1.6) pkg> add PoreMatMod
```

## Documentation 
[![Docs Badge](https://img.shields.io/badge/docs-latest-blue.svg)](https://SimonEnsemble.github.io/PoreMatMod.jl/)
Please visit our [documentation pages](https://SimonEnsemble.github.io/PoreMatMod.jl/) for detailed instructions and examples.

## Citing
```latex
@misc{henle2021pmm,
    title={PoreMatMod.jl: Julia package for in silico post-synthetic modification of crystal structure models.},
    author={E. Adrian Henle and Nickolas Gantzler and Praveen K. Thallapally and Xiaoli Z. Fern and Cory M. Simon}
    year={2021}
    preprint={ChemRxiv}
}
```
## Contributing

`PoreMatMod.jl` is being actively developed and we encourage feature requests and community feedback [on GitHub](https://github.com/SimonEnsemble/PoreMatMod.jl/issues).
