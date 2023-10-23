# Mixed finite element method in Julia

*Mixed formulation of the finite element method for the Poisson equation*

This is an almost direct translation of the MATLAB implementation described in the ([Bahriawat, 2005](https://www.degruyter.com/document/doi/10.2478/cmam-2005-0016/html#APA)) and available [here](https://www.math.hu-berlin.de/~cc/cc_homepage/software/software.shtml). Currently only the **EBmfem** method is available. 

You can find the data structures of inputs and outputs for examples given in Sections 3.1 and 9.1. They are currently in the `data` folder, as well as images of the obtained fields (displacement and fluxes).

## License

This project is licensed under the MIT license, see the [LICENSE](https://github.com/jhortizo/julia-mfem-poisson/blob/main/LICENSE.md) file for license rights and limitations. 

The documents are licensed under [Creative Commons Attribution License](http://creativecommons.org/licenses/by/4.0/).

## Contributing guidelines

For the nature of this project at the moment, we do not expect a lot of contributions outside of our research team. Nevertheless, having some guidelines such as style (linting, name variables) is a good practice. Regarding commits, we are using emojis in the commit messages following [gitmoji](https://gitmoji.dev/)'s recommendations.

## Reference

Bahriawati, C. & Carstensen, C. (2005). Three Matlab Implementations of the Lowest-order Raviart-Thomas Mfem with a Posteriori Error Control. Computational Methods in Applied Mathematics, 5(4), 333-361. https://doi.org/10.2478/cmam-2005-0016