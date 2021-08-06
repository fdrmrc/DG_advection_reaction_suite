#DG Upwind [Work still in progress]


The convergence table obtained by building a smooth manufactured solution `u(x,y)=epx(x) sin(y)` and using:
- uniformly refined mesh at each cycle, rather than a broken gradient estimator is
- `fe_degree = 1`
- the r.h.s function, derived from `Function` base class has value `exp(x)*(cos(y)x - sin(y)y)`

`
cycle cells  dofs        L2             H1       
    0    64    256 9.776e-04    - 5.471e-02    - 
    1   256   1024 2.463e-04 1.99 2.744e-02 1.00 
    2  1024   4096 6.180e-05 1.99 1.373e-02 1.00 
    3  4096  16384 1.548e-05 2.00 6.871e-03 1.00 
    4 16384  65536 3.873e-06 2.00 3.436e-03 1.00 
    5 65536 262144 9.686e-07 2.00 1.718e-03 1.00 
`

Which confirms the theoretical fact that for uniformly refined meshes and smooth solutions one usually see EOC `p+1`, with `p` the polynomial degree of the finite element.