# SauterSchwab3D.jl

Sauter&amp;Schwab quadrature rules [1] for singular tetrahedron-tetrahedron, tetrahedron-triangle and triangle-triangle interaction integrals.

There are currently two version of the implemented. The frist uses standard 1D Gauss quadrature as the base quadrature to build the tensor-product higher dimensional rules. The second version uses 1D, 2D, 3D and 4D qudrature rules defined for simplices to build a simplex tensor-product quadrature rule of higher dimension. This can save up to an order of magnitude in the number of quadrature points needed to achieve the same accuracy in some cases.

The simplex quadrature rule were developed by Shunn &amp; Ham [2].

There is also a third version in the code but it is not recommended as it use an other kind of simplex quadrature rule (GrundmannMoeller), but unfortunately it has negative quadrature weights.

## Quadrature Rules

Details about the idea behind these quadrature rules can be found here [3].

If one is simply interested in the tabulated quadrature rules, then in ``SimplexTensorProductQuadrature.pdf`` are all quadrature rules for all the different case. The reference elements in the table are consistent with the reference elements used in this code, but they might differ from the reference elements in the literature.

## Code

There are three main file.

- ``pulled_back_integrals.jl`` Here, the actual quadrature rules are implemented for the induvidual subdomains.
- ``reorder_vertices.jl`` Mainly importend helper functions to detect the singularity type and to reorder the vertices to have the integration elements correctly aligned to effectively treat the singular kernel.
- ``sauterschwabintegral.jl`` Wrapper functions around the externel lower dimensional quadratre rules. Functions that will call to all relevant subdomain intergrals to compute a interaction integral. 'sautesrschwab_parametrized' takes two arguments, the integral kernel and a struct that encodes the type of integration and the lower dimensional quadrature rules.

## Examples

For each case, there is an example file that generetates a error convergence study and benchmark for the timing. The name convention is the following structure 'example_S_D.jl'.
| S  | Meaning            |
|----|--------------------|
| ct | Common Tetrahedron |
| cf | Common Face        |
| ce | Common Edge        |
| cv | Common Vertex      |
| pd | Positive Distance  |

| D    | Meaning                 |
|------|-------------------------|
| 3D   | Tetrahedron-Tetrahedron |
| 2.5D | Tetrahedron-Triangle    |
| 2D   | Triangle-Triangle       |

For example, the file 'example_cv_2.5D.jl' contains an example how to use the quadrature rule in the case, when a tetrahedron and a triangle share a vertex.

## References

- [1] S. Sauter and C. Schwab, ``Generating the Matrix Coefficients,`` Boundary Element Methods, pp 289-352, 2011

- [2] L. Shunn and F. Ham, ``Symmetric quadrature rules for tetrahedra based on a cubic close-packed lattice arrangement,`` in Journal of Computational and Applied Mathematics, vol. 236, Issue 17, pp. 4348-4364,2012.

- [3] C. Münger and K. Cools, ``Efficient and kernel-independent Evaluation of Singular Integrals in Volume Integral Equations,`` 2021 IEEE International Conference on Microwaves, Antennas, Communications and Electronic Systems (COMCAS), 2021, pp. 188-192.
