#include "../include/DG_upwind_tester.h"

#include <gtest/gtest.h>

#include <fstream>
#include <sstream>

using namespace dealii;

using DG_upwindTestTypes =
    ::testing::Types<std::integral_constant<int, 2>,
                     std::integral_constant<int, 3>>; // test with dim=1,2,3

using DG_upwind2DTester =
    DG_upwind_Tester<std::integral_constant<int, 2>>; // if u want to test in 2D only

// TYPED_TEST_CASE(DG_upwind_Tester, DG_upwindTestTypes);

TEST_F(DG_upwind2DTester, TestLinearElements)
{
    std::stringstream str;
    str
        << "  subsection AdvectionProblem<2>" << std::endl
        << "  set Advection coefficient       = 0" << std::endl
        << "  set Finite element degree       = 1" << std::endl
        << "  set Fun expression              = exp(x)*sin(y)" << std::endl
        << "  set Global refinement           = true " << std::endl
        << "  set Number of global refinement = 2" << std::endl
        << "  set Number of refinement cycles = 4" << std::endl
        << "  set Output filename             = DG_upwind" << std::endl
        << "  set Problem constants           =  " << std::endl
        << "  set Right hand side expression  = exp(x)*(x*cos(y) - y*sin(y)) +0.000000*exp(x)*sin(y)" << std::endl
        << "  set Theta                       = 0.5" << std::endl
        << "  set Use direct solver           = true" << std::endl
        << "end" << std::endl;
    parse_string(str.str());
    run();
    auto tmp = solution;
    Vector<double> L2_error_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(mapping, dof_handler, solution,
                                      *fun, L2_error_per_cell, QGauss<2>(fe->tensor_degree() + 1), VectorTools::L2_norm);
    const double L2_error = VectorTools::compute_global_error(triangulation, L2_error_per_cell,
                                                              VectorTools::L2_norm);
    ASSERT_LE(L2_error, 1e-4);
}
//
//
//
TEST_F(DG_upwind2DTester, TestQuadraticElements)
{
    std::stringstream str2;
    str2
        << "  subsection AdvectionProblem<2>" << std::endl
        << "  set Advection coefficient       = 0.0" << std::endl
        << "  set Finite element degree       = 2" << std::endl
        << "  set Fun expression              = exp(x)*sin(y)" << std::endl
        << "  set Global refinement           = true " << std::endl
        << "  set Number of global refinement = 2" << std::endl
        << "  set Number of refinement cycles = 3" << std::endl
        << "  set Output filename             = DG_upwind" << std::endl
        << "  set Problem constants           =  " << std::endl
        << "  set Right hand side expression  = exp(x)*(x*cos(y) - y*sin(y)) + 0.0*exp(x)*sin(y)" << std::endl
        << "  set Theta                       = 0.5" << std::endl
        << "  set Use direct solver           = true" << std::endl
        << "end" << std::endl;
    parse_string(str2.str());

    this->run();
    Vector<double> L2_error_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(mapping, dof_handler, solution,
                                      *fun, L2_error_per_cell, QGauss<2>(fe->tensor_degree() + 1), VectorTools::L2_norm);
    const double L2_error = VectorTools::compute_global_error(triangulation, L2_error_per_cell,
                                                              VectorTools::L2_norm);
    ASSERT_LE(L2_error, 1e-5);
}

TEST_F(DG_upwind2DTester, TestAdvectionCoefficient)
{
    std::stringstream str;
    str
        << "  subsection AdvectionProblem<2>" << std::endl
        << "  set Advection coefficient       = 2.0" << std::endl
        << "  set Finite element degree       = 2" << std::endl
        << "  set Fun expression              = exp(x)*sin(y)" << std::endl
        << "  set Global refinement           = true " << std::endl
        << "  set Number of global refinement = 2" << std::endl
        << "  set Number of refinement cycles = 3" << std::endl
        << "  set Output filename             = DG_upwind" << std::endl
        << "  set Problem constants           =  " << std::endl
        << "  set Right hand side expression  = exp(x)*(x*cos(y) - y*sin(y)) + 2.0*exp(x)*sin(y)" << std::endl
        << "  set Theta                       = 0.5" << std::endl
        << "  set Use direct solver           = true" << std::endl
        << "end" << std::endl;
    parse_string(str.str());

    this->run();
    Vector<double> L2_error_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(mapping, dof_handler, solution,
                                      *fun, L2_error_per_cell, QGauss<2>(fe->tensor_degree() + 1), VectorTools::L2_norm);
    const double L2_error = VectorTools::compute_global_error(triangulation, L2_error_per_cell,
                                                              VectorTools::L2_norm);

    ASSERT_LE(L2_error, 1e-4);
}
