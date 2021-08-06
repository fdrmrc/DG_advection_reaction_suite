/*
 * DG_upwind.h
 *
 *  Created on: 6 Aug 2021
 *      Author: marcofeder
 */

#ifndef INCLUDE_DG_UPWIND_H_
#define INCLUDE_DG_UPWIND_H_

// The first few files have already been covered in previous examples and will
// thus not be further commented on:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/fe/fe_dgq.h>
// This header is needed for FEInterfaceValues to compute integrals on
// interfaces:
#include <deal.II/fe/fe_interface_values.h>
//Solver
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/precondition_block.h>
// We are going to use gradients as refinement indicator.
#include <deal.II/numerics/derivative_approximation.h>
// Using using the mesh_loop from the MeshWorker framework
#include <deal.II/meshworker/mesh_loop.h>


#include <deal.II/base/convergence_table.h>


#include <iostream>
#include <fstream>

using namespace dealii;
template<int dim>
class AdvectionProblem {
public:
	AdvectionProblem();
	void run();

private: //define usual private members
	void setup_system();
	void assemble_system();
	void solve();
	void refine_grid();
	void output_results(const unsigned int cycle) const;

	Triangulation<dim> triangulation;
	const MappingQ1<dim> mapping;

	// Furthermore we want to use DG elements.
	const FE_DGQ<dim> fe;
	DoFHandler<dim> dof_handler;

	const QGauss<dim> quadrature;
	const QGauss<dim - 1> quadrature_face;

	SparsityPattern sparsity_pattern;
	SparseMatrix<double> system_matrix;

	Vector<double> solution;
	Vector<double> right_hand_side;
	mutable ConvergenceTable convergence_table;


};

#endif /* INCLUDE_DG_UPWIND_H_ */
