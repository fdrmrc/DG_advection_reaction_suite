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
#include <deal.II/lac/sparse_direct.h>
// We are going to use gradients as refinement indicator.
#include <deal.II/numerics/derivative_approximation.h>
// Using using the mesh_loop from the MeshWorker framework
#include <deal.II/meshworker/mesh_loop.h>


#include <deal.II/base/convergence_table.h>

//To enable parameter handling
#include <deal.II/base/function_parser.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parsed_convergence_table.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/symbolic_function.h>



#include <iostream>
#include <fstream>

// template <typename Integral>
// class DG_upwind_Tester;

//Simple struct used only for throwing an exception when theta parameter is not okay.
struct theta_exc{
	std::string message;
	theta_exc(std::string&& s):message{std::move(s)}{};
	const char* what() const { return message.c_str(); }

};


using namespace dealii;
template<int dim>
class AdvectionProblem : ParameterAcceptor{
public:
	AdvectionProblem();
	void run();
	void initialize_params(const std::string& filename);
	void parse_string(const std::string& parameters);

	//private: //define usual private members
protected:
	void setup_system();
	void assemble_system();
	void solve();
	void refine_grid();
	void output_results(const unsigned int cycle) const;
	void compute_error();

	Triangulation<dim> triangulation;
	const MappingQ1<dim> mapping;

	// Furthermore we want to use DG elements.
	std::unique_ptr<FE_DGQ<dim>> fe;
	DoFHandler<dim> dof_handler;


	SparsityPattern sparsity_pattern;
	SparseMatrix<double> system_matrix;

	Vector<double> solution;
	Vector<double> right_hand_side;
//	mutable ConvergenceTable convergence_table; //specified mutable as it is in the const-marked method output_results

	//Parameter handling
	FunctionParser<dim> exact_solution;
	FunctionParser<dim> rhs;
	std::unique_ptr<Functions::SymbolicFunction<dim>> fun;

	unsigned int fe_degree = 3; //high order only for convergence tests

	std::string fun_expression = "exp(x)*sin(y)";
	std::map<std::string, double> constants;
	std::string output_filename = "DG_upwind";
	ParsedConvergenceTable error_table;

	bool use_direct_solver = true;
	unsigned int n_refinement_cycles = 6;
	unsigned int n_global_refinements = 2;
	bool global_refinement = true;
	double theta = 0.5; //default is 0.5 so that I have classical upwind flux

	double advection_coefficient = 0.0; //no reaction term

	std::string rhs_expression = "exp(x)*(x*cos(y) - y*sin(y)) +" + std::to_string(advection_coefficient) + "*exp(x)*sin(y)";

	template <typename Integral>
	friend class DG_upwind_Tester;   
};

#endif /* INCLUDE_DG_UPWIND_H_ */
