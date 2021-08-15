# DG_advection

Suite of C++ programs that solve adaptively, with Discontinuous Galerkin methods, advection dominated equations in 2D/3D, using a numerical flux approach. In particular, allowing different choices of `theta` parameter  in the weak form, we can recover the classical upwind scheme. The crucial point is that only normal components of the averages enter in the scheme. See paper by Brezzi-Marini-Süli-M3AS-2004 for an extended discussion. 

The correctness of the implementation is checked by observing the expected convergence rates on a sequence of uniformly refined meshes, for which we have rates `p+1` and `p` in L2 and H1, assuming a *smooth* solution, which we've built manually.

The parameters like theta, degree of the elements, refining strategy, and so on can be handled and tuned just by changing the `parameters.prm` file, without the need to recompile. This feature is provided by the `FunctionParser<dim>` class.