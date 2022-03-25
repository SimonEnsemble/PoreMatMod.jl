![logo.JPG](logo.jpg)
---

`PoreMatMod.jl` is a [Julia](https://julialang.org/) package for (i) subgraph matching and (ii) modifying crystal structures such as metal-organic frameworks (MOFs).

Functioning as a "find-and-replace" tool on atomistic crystal structure models of porous materials, `PoreMatMod.jl` is useful for:

:hammer: subgraph matching to filter databases of crystal structures

:hammer: constructing hypothetical crystal structure models of functionalized structures

:hammer: introducing defects into crystal structures

:hammer: repairing artifacts of X-ray structure determination, such as missing hydrogen atoms, disorder, and guest molecules

N.b. while `PoreMatMod.jl` was developed for MOFs and other porous crystalline materials, its find-and-replace operations can be applied to discrete, molecular structures as well by assigning an arbitrary unit cell.

| **Documentation** | **Build Status** | **Test Coverage** |
|:---:|:---:|:---:|
| [![Docs Badge](https://img.shields.io/badge/docs-dev-blue.svg)](https://SimonEnsemble.github.io/PoreMatMod.jl/dev) | [![CI](https://github.com/SimonEnsemble/PoreMatMod.jl/actions/workflows/CI_build.yml/badge.svg)](https://github.com/SimonEnsemble/PoreMatMod.jl/actions/workflows/CI_build.yml) [![Docs](https://github.com/SimonEnsemble/PoreMatMod.jl/actions/workflows/doc_deployment.yml/badge.svg)](https://github.com/SimonEnsemble/PoreMatMod.jl/actions/workflows/doc_deployment.yml) [![weekly](https://github.com/SimonEnsemble/PoreMatMod.jl/actions/workflows/weekly.yml/badge.svg)](https://github.com/SimonEnsemble/PoreMatMod.jl/actions/workflows/weekly.yml) | [![codecov](https://codecov.io/gh/SimonEnsemble/PoreMatMod.jl/branch/master/graph/badge.svg?token=Z9VMLXS3U9)](https://codecov.io/gh/SimonEnsemble/PoreMatMod.jl) [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl) |


## Installation
`PoreMatMod.jl` is a registered Julia package and can be installed by entering the following line in the Julia REPL when in package mode (type `]` to enter package mode):

```
(v1.7) pkg> add PoreMatMod
```

## Documentation

Link to documentation [here](https://simonensemble.github.io/PoreMatMod.jl/dev).

## Gallery of examples

Link to examples [here](https://simonensemble.github.io/PoreMatMod.jl/dev/examples/) with raw [Pluto](https://github.com/fonsp/Pluto.jl) notebooks [here](https://github.com/SimonEnsemble/PoreMatMod.jl/tree/master/examples).

## Citing

If you found `PoreMatMod.jl` useful, please cite our paper in *J. Chem. Inf. Model.* (ACS Editors' Choice) [here](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01219) [preprint [here](https://chemrxiv.org/engage/chemrxiv/article-details/615cf5127d3da5dd7bee4a22)]. :point_down:

```latex
@article{Henle2022,
  doi = {10.1021/acs.jcim.1c01219},
  url = {https://doi.org/10.1021/acs.jcim.1c01219},
  year = {2022},
  month = jan,
  publisher = {American Chemical Society ({ACS})},
  volume = {62},
  number = {3},
  pages = {423--432},
  author = {E. Adrian Henle and Nickolas Gantzler and Praveen K. Thallapally and Xiaoli Z. Fern and Cory M. Simon},
  title = {{PoreMatMod}.jl: Julia Package for in Silico Postsynthetic Modification of Crystal Structure Models},
  journal = {Journal of Chemical Information and Modeling}
}
```
## Contributing

We encourage feature requests and feedback [on GitHub](https://github.com/SimonEnsemble/PoreMatMod.jl/issues).

## License
`PoreMatMod.jl` is licensed under the [MIT license](./LICENSE).

