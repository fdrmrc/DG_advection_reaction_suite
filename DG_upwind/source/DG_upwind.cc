/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2021 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Guido Kanschat, Texas A&M University, 2009
 *         Timo Heister, Clemson University, 2019
 *         Extension by:
 *         Marco Feder, SISSA, 2021
 */

#include "../include/DG_upwind.h"



//static ConvergenceTable convergence_table;


template<int dim>
class RightHandSide: public Function<dim> {
public:
	virtual double value(const Point<dim> &p,
			const unsigned int component = 0) const override;
};

template<int dim>
double RightHandSide<dim>::value(const Point<dim> &p,
		const unsigned int component) const {
	(void) component;
	return std::exp(p[0]) * (cos(p[1]) * p[0] - sin(p[1]) * p[1]);
}

template<int dim>
class Solution: public Function<dim> {
public:
	virtual double value(const Point<dim> &p,
			const unsigned int component = 0) const override;

//	virtual Tensor<1,dim> gradient(const Point<dim> & point,
//	           const unsigned int component = 0) const override; //need to provide this so that I can compute H1 norm, using VectorTools namespace
};

template<int dim>
double Solution<dim>::value(const Point<dim> &p, const unsigned int) const {

	return std::exp(p[0]) * sin(p[1]); //e^x sin(y)
}

//template<int dim>
//Tensor<1,dim> Solution<dim>::gradient(const Point<dim> & point,
//        const unsigned int component) const {
//
//	(void) component; //suppress warning about unused parameters
//	Tensor<1,dim> sol_grad;
//	const double etpx = std::exp(point[0]);
//	sol_grad[0] = etpx*sin(point[1]);
//	sol_grad[1] = etpx*cos(point[1]);
//	return sol_grad;
//
//}

template<int dim>
class BoundaryValues: public Function<dim> {
public:
	BoundaryValues() = default;
	virtual void value_list(const std::vector<Point<dim>> &points,
			std::vector<double> &values, const unsigned int component = 0) const
					override;
};

// Given the flow direction, the inflow boundary of the unit square $[0,1]^2$
// are the right and the lower boundaries.
template<int dim>
void BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points,
		std::vector<double> &values, const unsigned int component) const {
	(void) component;
	AssertIndexRange(component, 1);
	Assert(values.size() == points.size(),
			ExcDimensionMismatch(values.size(), points.size()));
	Solution<dim> exact_solution;
	for (unsigned int i = 0; i < values.size(); ++i) {
		values[i] = exact_solution.value(points[i]);
	}

}

// Finally, a function that computes and returns the wind field
// $\beta=\beta(\mathbf x)$. As explained in the introduction, we will use a
// rotational field around the origin in 2d. In 3d, we simply leave the
// $z$-component unset (i.e., at zero)
template<int dim>
Tensor<1, dim> beta(const Point<dim> &p) {
	Assert(dim >= 2, ExcNotImplemented());

	Tensor<1, dim> wind_field;
	wind_field[0] = -p[1];
	wind_field[1] = p[0];
//	if (wind_field.norm() > 1e-10)
//		wind_field /= wind_field.norm();

	return wind_field;
}

// The following objects are the scratch and copy objects we use in the call
// to MeshWorker::mesh_loop().
template<int dim>
struct ScratchData {
	ScratchData(const Mapping<dim> &mapping, const FiniteElement<dim> &fe,
			const Quadrature<dim> &quadrature,
			const Quadrature<dim - 1> &quadrature_face,
			const UpdateFlags update_flags = update_values | update_gradients
					| update_quadrature_points | update_JxW_values,
			const UpdateFlags interface_update_flags = update_values
					| update_gradients | update_quadrature_points
					| update_JxW_values | update_normal_vectors) :
			fe_values(mapping, fe, quadrature, update_flags), fe_interface_values(
					mapping, fe, quadrature_face, interface_update_flags) {
	}

	ScratchData(const ScratchData<dim> &scratch_data) :
			fe_values(scratch_data.fe_values.get_mapping(),
					scratch_data.fe_values.get_fe(),
					scratch_data.fe_values.get_quadrature(),
					scratch_data.fe_values.get_update_flags()), fe_interface_values(
					scratch_data.fe_interface_values.get_mapping(),
					scratch_data.fe_interface_values.get_fe(),
					scratch_data.fe_interface_values.get_quadrature(),
					scratch_data.fe_interface_values.get_update_flags()) {
	}

	FEValues<dim> fe_values;
	FEInterfaceValues<dim> fe_interface_values; //lavora sull'interfaccia di due celle.
};

struct CopyDataFace {
	FullMatrix<double> cell_matrix;
	std::vector<types::global_dof_index> joint_dof_indices;
};

struct CopyData {
	FullMatrix<double> cell_matrix;
	Vector<double> cell_rhs;
	std::vector<types::global_dof_index> local_dof_indices;
	std::vector<CopyDataFace> face_data;

	template<class Iterator>
	void reinit(const Iterator &cell, unsigned int dofs_per_cell) {
		cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
		cell_rhs.reinit(dofs_per_cell);

		local_dof_indices.resize(dofs_per_cell);
		cell->get_dof_indices(local_dof_indices);
	}
};

// We start with the constructor. The 1 in the constructor call of
// <code>fe</code> is the polynomial degree.
template<int dim>
AdvectionProblem<dim>::AdvectionProblem() :
		mapping(), dof_handler(triangulation){

	add_parameter("Finite element degree", fe_degree);
//	add_parameter("Exact solution expression", exact_solution_expression);
	add_parameter("Problem constants", constants);
	add_parameter("Output filename", output_filename);
	add_parameter("Use direct solver", use_direct_solver);
	add_parameter("Number of refinement cycles", n_refinement_cycles);
	add_parameter("Number of global refinement",n_global_refinements);
	add_parameter("Global refinement", global_refinement);
	add_parameter("Fun expression", fun_expression);

//
	this->prm.enter_subsection("Error table");
	error_table.add_parameters(this->prm);
	this->prm.leave_subsection();


}


template<int dim>
void AdvectionProblem<dim>::initialize(const std::string& filename){

	ParameterAcceptor::initialize(filename, "", ParameterHandler::Short);
}


template <int dim>
void
AdvectionProblem<dim>::parse_string(const std::string &parameters)
{
  ParameterAcceptor::prm.parse_input_from_string(parameters);
  ParameterAcceptor::parse_all_parameters();
}







template<int dim>
void AdvectionProblem<dim>::setup_system() {

	// first need to distribute the DoFs.
	if (!fe){
		fe              = std::make_unique<FE_DGQ<dim>>(fe_degree);
//		const auto vars = dim == 2 ? "x,y" : "x,y,z";
//		exact_solution.initialize(vars, exact_solution_expression, constants);
		fun = std::make_unique<Functions::SymbolicFunction<dim>>(fun_expression);
	}


	dof_handler.distribute_dofs(*fe);

	// To build the sparsity pattern for DG discretizations, we can call the
	// function analogue to DoFTools::make_sparsity_pattern, which is called
	// DoFTools::make_flux_sparsity_pattern:
	DynamicSparsityPattern dsp(dof_handler.n_dofs());
	DoFTools::make_flux_sparsity_pattern(dof_handler, dsp); //DG sparsity pattern generator
	sparsity_pattern.copy_from(dsp);

	// Finally, we set up the structure of all components of the linear system.
	system_matrix.reinit(sparsity_pattern);
	solution.reinit(dof_handler.n_dofs());
	right_hand_side.reinit(dof_handler.n_dofs());
}

//in the call to  MeshWorker::mesh_loop() we only need to specify what should happen on
// each cell, each boundary face, and each interior face. These three tasks
// are handled by the lambda functions inside the function below.

template<int dim>
void AdvectionProblem<dim>::assemble_system() {
	using Iterator = typename DoFHandler<dim>::active_cell_iterator;
	const BoundaryValues<dim> boundary_function;
	const RightHandSide<dim> rhs_function;

	const QGauss<dim> quadrature= fe->tensor_degree() +1;
	const QGauss<dim - 1> quadrature_face = fe->tensor_degree() +1;



	// This is the function that will be executed for each cell.
	const auto cell_worker = [&](const Iterator &cell,
			ScratchData<dim> &scratch_data, CopyData &copy_data) {
		const unsigned int n_dofs =
				scratch_data.fe_values.get_fe().n_dofs_per_cell();
		copy_data.reinit(cell, n_dofs);
		scratch_data.fe_values.reinit(cell);

		const auto &q_points = scratch_data.fe_values.get_quadrature_points();

		const FEValues<dim> &fe_v = scratch_data.fe_values;
		const std::vector<double> &JxW = fe_v.get_JxW_values();

		// We solve a homogeneous equation, thus no right hand side shows up in
		// the cell term.  What's left is integrating the matrix entries.
		for (unsigned int point = 0; point < fe_v.n_quadrature_points;
				++point) {
			auto beta_q = beta(q_points[point]);
			for (unsigned int i = 0; i < n_dofs; ++i) {
				for (unsigned int j = 0; j < n_dofs; ++j) {
					copy_data.cell_matrix(i, j) += -beta_q             // -\beta
					* fe_v.shape_grad(i, point)  // \nabla \phi_i
							* fe_v.shape_value(j, point) // \phi_j
							* JxW[point];                // dx
				}
				copy_data.cell_rhs(i) +=
						rhs_function.value(q_points[point]) // f(x_q)
						* fe_v.shape_value(i, point) //phi_i(x_q)
						* JxW[point]; //dx
			}
		}
	};

	// This is the function called for boundary faces and consists of a normal
	// integration using FEFaceValues. New is the logic to decide if the term
	// goes into the system matrix (outflow) or the right-hand side (inflow).
	const auto boundary_worker = [&](const Iterator &cell,
			const unsigned int &face_no, ScratchData<dim> &scratch_data,
			CopyData &copy_data) {
		scratch_data.fe_interface_values.reinit(cell, face_no);
		const FEFaceValuesBase<dim> &fe_face =
				scratch_data.fe_interface_values.get_fe_face_values(0);

		const auto &q_points = fe_face.get_quadrature_points();

		const unsigned int n_facet_dofs = fe_face.get_fe().n_dofs_per_cell();
		const std::vector<double> &JxW = fe_face.get_JxW_values();
		const std::vector<Tensor<1, dim>> &normals =
				fe_face.get_normal_vectors();

		std::vector<double> g(q_points.size());
		boundary_function.value_list(q_points, g);

		for (unsigned int point = 0; point < q_points.size(); ++point) {
			const double beta_dot_n = beta(q_points[point]) * normals[point];

			if (beta_dot_n > 0) {
				for (unsigned int i = 0; i < n_facet_dofs; ++i)
					for (unsigned int j = 0; j < n_facet_dofs; ++j)
						copy_data.cell_matrix(i, j) += fe_face.shape_value(i,
								point)   // \phi_i
						* fe_face.shape_value(j, point) // \phi_j
								* beta_dot_n                    // \beta . n
								* JxW[point];                   // dx
			} else
				for (unsigned int i = 0; i < n_facet_dofs; ++i)
					copy_data.cell_rhs(i) += -fe_face.shape_value(i, point) // \phi_i
					* /*exact_solution.value(q_points[point])*/ g[point]                     // g
							* beta_dot_n  // \beta . n
							* JxW[point]; // dx
		}
	};

	// This is the function called on interior faces. The arguments specify
	// cells, face and subface indices (for adaptive refinement). We just pass
	// them along to the reinit() function of FEInterfaceValues.
	const auto face_worker = [&](const Iterator &cell, const unsigned int &f,
			const unsigned int &sf, const Iterator &ncell,
			const unsigned int &nf, const unsigned int &nsf,
			ScratchData<dim> &scratch_data, CopyData &copy_data) {
		FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
		fe_iv.reinit(cell, f, sf, ncell, nf, nsf);
		const auto &q_points = fe_iv.get_quadrature_points();

		copy_data.face_data.emplace_back(); // @suppress("Method cannot be resolved")
		CopyDataFace &copy_data_face = copy_data.face_data.back();

		const unsigned int n_dofs = fe_iv.n_current_interface_dofs();
		copy_data_face.joint_dof_indices = fe_iv.get_interface_dof_indices();

		copy_data_face.cell_matrix.reinit(n_dofs, n_dofs);

		const std::vector<double> &JxW = fe_iv.get_JxW_values();
		const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

		for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint) {
			const double beta_dot_n = beta(q_points[qpoint]) * normals[qpoint];
			for (unsigned int i = 0; i < n_dofs; ++i)
				for (unsigned int j = 0; j < n_dofs; ++j)
					copy_data_face.cell_matrix(i, j) += fe_iv.jump(i, qpoint) // [\phi_i]
					* fe_iv.shape_value((beta_dot_n > 0), j, qpoint) // phi_j^{upwind}
							* beta_dot_n                          // (\beta . n)
							* JxW[qpoint];                                 // dx
		}
	};

	// The following lambda function will handle copying the data from the
	// cell and face assembly into the global matrix and right-hand side.
	//
	// While we would not need an AffineConstraints object, because there are
	// no hanging node constraints in DG discretizations, we use an empty
	// object here as this allows us to use its `copy_local_to_global`
	// functionality.
	const AffineConstraints<double> constraints;

	const auto copier = [&](const CopyData &c) {
		constraints.distribute_local_to_global(c.cell_matrix, c.cell_rhs,
				c.local_dof_indices, system_matrix, right_hand_side);

		for (auto &cdf : c.face_data) {
			constraints.distribute_local_to_global(cdf.cell_matrix,
					cdf.joint_dof_indices, system_matrix);
		}
	};

	ScratchData<dim> scratch_data(mapping, *fe, quadrature, quadrature_face);
	CopyData copy_data;

	// Here, we finally handle the assembly. We pass in ScratchData and
	// CopyData objects, the lambda functions from above, an specify that we
	// want to assemble interior faces once.
	MeshWorker::mesh_loop(dof_handler.begin_active(), dof_handler.end(),
			cell_worker, copier, scratch_data, copy_data,
			MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces
					| MeshWorker::assemble_own_interior_faces_once,
			boundary_worker, face_worker);
}

template<int dim>
void AdvectionProblem<dim>::solve() {

	if (use_direct_solver){

		SparseDirectUMFPACK system_matrix_inverse;
		system_matrix_inverse.initialize(system_matrix);
		system_matrix_inverse.vmult(solution, right_hand_side);
	}else{
		SolverControl solver_control(1000, 1e-15);
		SolverRichardson<Vector<double>> solver(solver_control);

		// Here we create the preconditioner,
		PreconditionBlockSSOR<SparseMatrix<double>> preconditioner;

		// then assign the matrix to it and set the right block size:
		preconditioner.initialize(system_matrix, fe->n_dofs_per_cell());

		// After these preparations we are ready to start the linear solver.
		solver.solve(system_matrix, solution, right_hand_side, preconditioner);

		std::cout << "  Solver converged in " << solver_control.last_step()
	    					<< " iterations." << std::endl;
	}
}




// We refine the grid according to a very simple refinement criterion, namely
// an approximation to the gradient of the solution.
template<int dim>
void AdvectionProblem<dim>::refine_grid() {



	// The <code>DerivativeApproximation</code> class computes the gradients to
	// float precision. This is sufficient as they are approximate and serve as
	// refinement indicators only.
	if (!global_refinement){
	Vector<float> gradient_indicator(triangulation.n_active_cells());

	// Now the approximate gradients are computed
	DerivativeApproximation::approximate_gradient(mapping, dof_handler,
			solution, gradient_indicator);

	// and they are cell-wise scaled by the factor $h^{1+d/2}$
	unsigned int cell_no = 0;
	for (const auto &cell : dof_handler.active_cell_iterators())
		gradient_indicator(cell_no++) *= std::pow(cell->diameter(),
				1 + 1.0 * dim / 2);

	// Finally they serve as refinement indicator.
	GridRefinement::refine_and_coarsen_fixed_number(triangulation,
			gradient_indicator, 0.3, 0.1);

	triangulation.execute_coarsening_and_refinement();
	}else{
	triangulation.refine_global(1); //just for testing on uniformly refined meshes
	}
}

// The output of this program consists of a vtk file of the adaptively
// refined grids and the numerical solutions. Finally, we also compute the
// L-infinity norm of the solution using VectorTools::integrate_difference().
template<int dim>
void AdvectionProblem<dim>::output_results(const unsigned int cycle) const {
	const std::string filename = "solution-" + std::to_string(cycle) + ".vtk";
	std::cout << "  Writing solution to <" << filename << ">" << std::endl;
	std::ofstream output(filename);

	DataOut<dim> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, "u", DataOut<dim>::type_dof_data);

	data_out.build_patches(mapping);

	data_out.write_vtk(output);

//	const unsigned int n_active_cells = triangulation.n_active_cells();
//	const unsigned int n_dofs         = dof_handler.n_dofs();

//
//	Vector<double> L2_error_per_cell(n_active_cells);
//	VectorTools::integrate_difference(mapping, dof_handler, solution,
//				Solution<dim>(), L2_error_per_cell, QGauss<dim>(fe->tensor_degree()+1), VectorTools::L2_norm);
//	const double L2_error = VectorTools::compute_global_error(triangulation, L2_error_per_cell,
//				VectorTools::L2_norm);
//
//	Vector<double> H1_error_per_cell(n_active_cells);
//	VectorTools::integrate_difference(mapping,dof_handler,solution,
//			  Solution<dim>(),H1_error_per_cell,QGauss<dim>(fe->tensor_degree()+1),VectorTools::H1_norm);
//	const double H1_error = VectorTools::compute_global_error(triangulation, H1_error_per_cell,
//			  VectorTools::H1_norm);
//
//
//
////
//	convergence_table.add_value("cycle", cycle);
//	convergence_table.add_value("cells", n_active_cells);
//	convergence_table.add_value("dofs", n_dofs);
//	convergence_table.add_value("L2", L2_error);
//	convergence_table.add_value("H1", H1_error);
	//	  convergence_table.add_value("Linfty", Linfty_error);







}



template<int dim>
void AdvectionProblem<dim>::compute_error(){
	error_table.error_from_exact(mapping, dof_handler, solution, *fun); //be careful: a FD approximation of the gradient is used to compute the H^1 norm if you're not relying on SymbolicFunction class
//	error_table.error_from_exact(mapping, dof_handler, solution, Solution<dim>()); //provided that Solution<dim> implements the Gradient function

}









//Usual run functions, running over several refinement cycles
template<int dim>
void AdvectionProblem<dim>::run() {
	for (unsigned int cycle = 0; cycle < n_refinement_cycles; ++cycle) {
		std::cout << "Cycle " << cycle << std::endl;

		if (cycle == 0) {
			GridGenerator::hyper_cube(triangulation);
			triangulation.refine_global(n_global_refinements);
		} else
			refine_grid();

		std::cout << "  Number of active cells:       "
				<< triangulation.n_active_cells() << std::endl;

		setup_system();

		std::cout << "  Number of degrees of freedom: " << dof_handler.n_dofs()
				<< std::endl;

		assemble_system();
		solve();
		compute_error();
		output_results(cycle);
	}


	std::cout <<"\n" <<"With ParsedConvergenceTable class: "<<"\n";
	error_table.output_table(std::cout);



	/*convergence_table.set_precision("L2", 3);
	convergence_table.set_precision("H1", 3);

	convergence_table.set_scientific("L2", true);
	convergence_table.set_scientific("H1", true);

	convergence_table.evaluate_convergence_rates(
			"L2", ConvergenceTable::reduction_rate_log2);

	convergence_table.evaluate_convergence_rates(
			"H1", ConvergenceTable::reduction_rate_log2);


	convergence_table.write_text(std::cout);*/


}

