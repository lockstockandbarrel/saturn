# Code Modernization experiment

## Scientific-Subroutine-Package

This is an experiment looking at modernizing an early
IBM Fortran library, as discussed at

   https://fortran-lang.discourse.group/t/anecdotal-fortran/704/414 

It is based on early Fortran code that composed the IBM 360 Scientific
Subroutine Package, probably circa 1968.

A series of GNU/Linux commands (sed, vim, cut, expand, dos2unix, ...)
and the spag(1) refactoring software were applied to the files and
the code was restructured as an fpm(1) (Fortran Package Manager) package
in a github repository. Unit tests using the M_framework fpm(1) package
and using Grok to refine the code and create the unit tests is being
considered but for anyone curious what is here so far this took only a few
man-hours to accompish.

There are still several procedures where work storage has type mismatches,
and questions about how appropriate several quick changes were that were
made to make the code compile. This is temporarily being presented as 
part of the Discourse discussion and will be removed shortly, but if 
it is clear this can proceed and there is sufficient interest it will
hopefully re-appear as a proper github repository.

With no unit tests available so far so the only thing that can be said at 
this point is that everything compiles and all but five of the procedures
have been moved to a free-format module, but the experience has been
interesting so far, so I would not recommend anything be used in production
until it is clear about the pedigree and availability of the original 
version even though the code is extensively transformed; and that unit
and regression tests are included.

The topics touched on by the library are:

### Statistics

 + Probit analysis
 + Analysis of variance (factorial design)
 + Correlation analysis
 + Multiple linear regression
 + Stepwise regression
 + Polynomial regression
 + Canonical correlation
 + Factor analysis (principal components, varimax)
 + Discriminant analysis (many groups)
 + Time series analysis
 + Data screening and analysis
 + Nonparametric tests
 + Random number generation (uniform, normal)
 + Distribution functions

### Mathematics

 + Inversion
 + Eigenvalues and eigenvectors
 + Simultaneous linear algebraic equations
 + Transpositions
 + Matrix arithmetic (addition, product, etc.)
 + Matrix partitioning
 + Matrix tabulation and sorting of rows or columns
 + Elementary operations on rows or columns of matrices
 + Matrix factorization
 + Integration and differentiation of given or tabulated functions
 + Solution of systems of first-order differential equations
 + Fourier analysis of given or tabulated functions
 + Bessel and modified Bessel function evaluation
 + Gamma function evaluation
 + Jacobian elliptic functions
 + Elliptic, exponential, sine cosine, Fresnel integrals
 + Finding real roots of a given function
 + Finding real and complex roots of a real polynomial
 + Polynomial arithmetic (addition, division, etc.)
 + Polynomial evaluation, integration, differentiation
 + Chebyshev, Hermite, Laguerre, Legendre polynomials
 + Minimum of a function
 + Approximation, interpolation, and table construction

### References

 + http://www.ecalculations.com/
 + http://www.ebyte.it/library/codesnippets/IBM_System360_SSP.html
 + https://www.discourse.org/
 + https://fortran.uk/fortran-analysis-and-refactoring-with-plusfort/spag-fortran-code-restructuring/
