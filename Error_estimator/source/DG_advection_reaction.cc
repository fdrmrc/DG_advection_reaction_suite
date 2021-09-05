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

#include "../include/DG_advection_reaction.h"

//template<int dim>
//class RightHandSide: public Function<dim> {
//public:
//	virtual double value(const Point<dim> &p,
//			const unsigned int component = 0) const override;
//};
//
//template<int dim>
//double RightHandSide<dim>::value(const Point<dim> &p,
//		const unsigned int component) const {
//	(void) component;
//	return std::exp(p[0]) * (cos(p[1]) * p[0] - sin(p[1]) * p[1]) + advection_coefficient*std::exp(p[0])*sin(p[1]);
//}

// template <int dim>
// class Solution : public Function<dim>
// {
// public:
// 	virtual double value(const Point<dim> &p,
// 						 const unsigned int component = 0) const override;

// 	//	virtual Tensor<1,dim> gradient(const Point<dim> & point,
// 	//	           const unsigned int component = 0) const override; //need to provide this so that I can compute H1 norm, using VectorTools namespace
// };

// template <int dim>
// double Solution<dim>::value(const Point<dim> &p, const unsigned int) const
// {

// 	return std::exp(p[0]) * sin(p[1]); //e^x sin(y)
// }

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

// template <int dim>
// class BoundaryValues : public Function<dim>
// {
// public:
// 	BoundaryValues() = default;
// 	virtual void value_list(const std::vector<Point<dim>> &points,
// 							std::vector<double> &values, const unsigned int component = 0) const
// 		override;
// };

// Given the flow direction, the inflow boundary of the unit square $[0,1]^2$
// are the right and the lower boundaries.
// template <int dim>
// void BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points,
// 									 std::vector<double> &values, const unsigned int component) const
// {
// 	(void)component;
// 	AssertIndexRange(component, 1);
// 	Assert(values.size() == points.size(),
// 		   ExcDimensionMismatch(values.size(), points.size()));
// 	Solution<dim> exact_solution;
// 	for (unsigned int i = 0; i < values.size(); ++i)
// 	{
// 		values[i] = exact_solution.value(points[i]);
// 	}
// }

// template <int dim>
// void BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points,
//                                      std::vector<double> &          values,
//                                      const unsigned int component) const
// {
//   (void)component;
//   AssertIndexRange(component, 1);
//   Assert(values.size() == points.size(),
//          ExcDimensionMismatch(values.size(), points.size()));
//   for (unsigned int i = 0; i < values.size(); ++i)
//     {
//       if (points[i](0) < 0.5)
//         values[i] = 1.;
//       else
//         values[i] = 0.;
//     }
// }

// Finally, a function that computes and returns the wind field
// $\beta=\beta(\mathbf x)$.
template <int dim>
Tensor<1, dim> beta(const Point<dim> &p)
{
	Assert(dim >= 2, ExcNotImplemented());
	(void)p;
	Tensor<1, dim> wind_field;
	wind_field[0] = 1.0;//;1.0;//-p[1];
	wind_field[1] = 1.0;//p[0];

	// const double r{std::sqrt(p[0]*p[0] + p[1]*p[1])};
	// wind_field[0] = -r*std::sin(std::atan2(p[1],p[0]));
	// wind_field[1] = +r*std::cos(std::atan2(p[1],p[0]));

	// if (wind_field.norm() > 1e-10)
	// 	wind_field /= wind_field.norm();

	return wind_field;
}

// The following objects are the scratch and copy objects we use in the call
// to MeshWorker::mesh_loop().
template <int dim>
struct ScratchData
{
	ScratchData(const Mapping<dim> &mapping, const FiniteElement<dim> &fe,
				const Quadrature<dim> &quadrature,
				const Quadrature<dim - 1> &quadrature_face,
				const UpdateFlags update_flags = update_values | update_gradients | update_quadrature_points | update_JxW_values,
				const UpdateFlags interface_update_flags = update_values | update_gradients | update_quadrature_points | update_JxW_values | update_normal_vectors) : fe_values(mapping, fe, quadrature, update_flags), fe_interface_values(mapping, fe, quadrature_face, interface_update_flags)
	{
	}

	ScratchData(const ScratchData<dim> &scratch_data) : fe_values(scratch_data.fe_values.get_mapping(),
																  scratch_data.fe_values.get_fe(),
																  scratch_data.fe_values.get_quadrature(),
																  scratch_data.fe_values.get_update_flags()),
														fe_interface_values(
															scratch_data.fe_interface_values.get_mapping(),
															scratch_data.fe_interface_values.get_fe(),
															scratch_data.fe_interface_values.get_quadrature(),
															scratch_data.fe_interface_values.get_update_flags())
	{
	}

	FEValues<dim> fe_values;
	FEInterfaceValues<dim> fe_interface_values; //lavora sull'interfaccia di due celle.
};

struct CopyDataFace
{
	FullMatrix<double> cell_matrix;
	std::vector<types::global_dof_index> joint_dof_indices;
	std::array<double, 2> values;
	std::array<unsigned int, 2> cell_indices;
};

struct CopyData
{
	FullMatrix<double> cell_matrix;
	Vector<double> cell_rhs;
	std::vector<types::global_dof_index> local_dof_indices;
	std::vector<CopyDataFace> face_data;

	double value;
	double value_estimator;
	unsigned int cell_index;

	FullMatrix<double> cell_mass_matrix;
	Vector<double> cell_mass_rhs;

	FullMatrix<double> cell_reconstruction_matrix;
	Vector<double> cell_reconstruction_rhs;

	template <class Iterator>
	void reinit(const Iterator &cell, unsigned int dofs_per_cell)
	{
		cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
		cell_mass_matrix.reinit(dofs_per_cell, dofs_per_cell);
		cell_reconstruction_matrix.reinit(dofs_per_cell, dofs_per_cell);

		cell_rhs.reinit(dofs_per_cell);
		cell_mass_rhs.reinit(dofs_per_cell);
		cell_reconstruction_rhs.reinit(dofs_per_cell); //TODO

		local_dof_indices.resize(dofs_per_cell);
		cell->get_dof_indices(local_dof_indices);
	}
};

template <int dim>
void get_function_jump(const FEInterfaceValues<dim> &fe_iv,
					   const Vector<double> &solution,
					   std::vector<double> &jump)
{
	const unsigned int n_q = fe_iv.n_quadrature_points;
	std::array<std::vector<double>, 2> face_values;
	jump.resize(n_q);
	for (unsigned int i = 0; i < 2; ++i)
	{
		face_values[i].resize(n_q);
		fe_iv.get_fe_face_values(i).get_function_values(solution,
														face_values[i]);
	}
	for (unsigned int q = 0; q < n_q; ++q)
		jump[q] = face_values[0][q] - face_values[1][q];
}

// We start with the constructor. The 1 in the constructor call of
// <code>fe</code> is the polynomial degree.
template <int dim>
AdvectionReaction<dim>::AdvectionReaction() : mapping(),
											  dof_handler(triangulation)/*,
											  dof_handler_continuous(triangulation),
											  fe_continuous(1)*/
{

	add_parameter("Finite element degree", fe_degree);
	add_parameter("Problem constants", constants);
	add_parameter("Output filename", output_filename);
	add_parameter("Use direct solver", use_direct_solver);
	add_parameter("Number of refinement cycles", n_refinement_cycles);
	add_parameter("Number of global refinement", n_global_refinements);
	add_parameter("Refinement", refinement);
	add_parameter("Exact solution expression", exact_solution_expression); //used only to test convegence rates
	add_parameter("Boundary conditions expression", boundary_conditions_expression); //used only to test convegence rates
	add_parameter("Theta", theta);
	add_parameter("Advection coefficient expression", advection_coefficient_expression);
	add_parameter("Right hand side expression", rhs_expression);

	//
	this->prm.enter_subsection("Error table");
	error_table.add_parameters(this->prm);
	this->prm.leave_subsection();
}

template <int dim>
void AdvectionReaction<dim>::initialize_params(const std::string &filename)
{

	ParameterAcceptor::initialize(filename, "", ParameterHandler::Short);
	if (theta < 0.0 || theta > 10.0 || std::abs(theta) < 1e-12)
	{
		throw(theta_exc("Theta parameter is not in a suitable range: see paper by Brezzi, Marini, Suli for an extended discussion"));
	}
}

template <int dim>
void AdvectionReaction<dim>::parse_string(const std::string &parameters)
{
	ParameterAcceptor::prm.parse_input_from_string(parameters);
	ParameterAcceptor::parse_all_parameters();
}

template <int dim>
void AdvectionReaction<dim>::setup_system()
{

	// first need to distribute the DoFs.
	if (!fe)
	{
		fe = std::make_unique<FE_DGQ<dim>>(fe_degree);
		const auto vars = dim == 2 ? "x,y" : "x,y,z";
		exact_solution.initialize(vars, exact_solution_expression, constants);
		rhs.initialize(vars, rhs_expression, constants);
		advection_coeff.initialize(vars,advection_coefficient_expression,constants);
		boundary_conditions.initialize(vars,boundary_conditions_expression,constants);
	}
	dof_handler.distribute_dofs(*fe);
	// dof_handler_continuous.distribute_dofs(fe_continuous);

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

template <int dim>
void AdvectionReaction<dim>::assemble_system()
{

	using Iterator = typename DoFHandler<dim>::active_cell_iterator;
	// const BoundaryValues<dim> boundary_function;
	//	const RightHandSide<dim> rhs_function;

	const QGauss<dim> quadrature = fe->tensor_degree() + 1;
	const QGauss<dim - 1> quadrature_face = fe->tensor_degree() + 1;

	// This is the function that will be executed for each cell.
	const auto cell_worker = [&](const Iterator &cell,
								 ScratchData<dim> &scratch_data, CopyData &copy_data)
	{
		FEValues<dim> fe_values_continuous(*fe,
										   quadrature,
										   update_values | update_gradients |
											   update_quadrature_points | update_JxW_values);

		const unsigned int n_dofs = scratch_data.fe_values.get_fe().n_dofs_per_cell();
		// std::cout << n_dofs << "sulla cella"
		// 		  << "\t" << fe_values_continuous.dofs_per_cell << "\n";
		copy_data.reinit(cell, n_dofs);
		scratch_data.fe_values.reinit(cell);

		const auto &q_points = scratch_data.fe_values.get_quadrature_points();

		const FEValues<dim> &fe_v = scratch_data.fe_values;
		const std::vector<double> &JxW = fe_v.get_JxW_values();

		// We solve a homogeneous equation, thus no right hand side shows up in
		// the cell term.  What's left is integrating the matrix entries.
		for (unsigned int point = 0; point < fe_v.n_quadrature_points;
			 ++point)
		{
			auto beta_q = beta(q_points[point]);
			for (unsigned int i = 0; i < n_dofs; ++i)
			{
				for (unsigned int j = 0; j < n_dofs; ++j)
				{
					copy_data.cell_matrix(i, j) += (-beta_q							 // -\beta
														* fe_v.shape_grad(i, point)	 // \nabla \phi_i
														* fe_v.shape_value(j, point) // \phi_j
													+
													advection_coeff.value(q_points[point]) *			 //gamma
														fe_v.shape_value(i, point)	 //phi_i
														* fe_v.shape_value(j, point) //phi_j
													) *
												   JxW[point]; // dx
				}
				copy_data.cell_rhs(i) +=
					/*rhs_function.value(q_points[point]) */
					rhs.value(q_points[point])	 // f(x_q)
					* fe_v.shape_value(i, point) //phi_i(x_q)
					* JxW[point];				 //dx
			}
		}
	};

	// This is the function called for boundary faces and consists of a normal
	// integration using FEFaceValues. New is the logic to decide if the term
	// goes into the system matrix (outflow) or the right-hand side (inflow).
	const auto boundary_worker = [&](const Iterator &cell,
									 const unsigned int &face_no, ScratchData<dim> &scratch_data,
									 CopyData &copy_data)
	{
		scratch_data.fe_interface_values.reinit(cell, face_no);
		const FEFaceValuesBase<dim> &fe_face =
			scratch_data.fe_interface_values.get_fe_face_values(0);

		const auto &q_points = fe_face.get_quadrature_points();

		const unsigned int n_facet_dofs = fe_face.get_fe().n_dofs_per_cell();
		const std::vector<double> &JxW = fe_face.get_JxW_values();
		const std::vector<Tensor<1, dim>> &normals =
			fe_face.get_normal_vectors();

		std::vector<double> g(q_points.size());
		exact_solution.value_list(q_points, g);

		for (unsigned int point = 0; point < q_points.size(); ++point)
		{
			const double beta_dot_n = beta(q_points[point]) * normals[point];

			if (beta_dot_n > 0)
			{
				for (unsigned int i = 0; i < n_facet_dofs; ++i)
					for (unsigned int j = 0; j < n_facet_dofs; ++j)
						copy_data.cell_matrix(i, j) += fe_face.shape_value(i,
																		   point)	   // \phi_i
													   * fe_face.shape_value(j, point) // \phi_j
													   * beta_dot_n					   // \beta . n
													   * JxW[point];				   // dx
			}
			else
				for (unsigned int i = 0; i < n_facet_dofs; ++i)
					copy_data.cell_rhs(i) += -fe_face.shape_value(i, point) // \phi_i
											 * g[point]						// g*/
											 * beta_dot_n					// \beta . n
											 * JxW[point];					// dx
		}
	};

	// This is the function called on interior faces. The arguments specify
	// cells, face and subface indices (for adaptive refinement). We just pass
	// them along to the reinit() function of FEInterfaceValues.
	const auto face_worker = [&](const Iterator &cell, const unsigned int &f,
								 const unsigned int &sf, const Iterator &ncell,
								 const unsigned int &nf, const unsigned int &nsf,
								 ScratchData<dim> &scratch_data, CopyData &copy_data)
	{
		FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
		fe_iv.reinit(cell, f, sf, ncell, nf, nsf);
		const auto &q_points = fe_iv.get_quadrature_points();

		copy_data.face_data.emplace_back(); // @suppress("Method cannot be resolved")
		CopyDataFace &copy_data_face = copy_data.face_data.back();

		const unsigned int n_dofs = fe_iv.n_current_interface_dofs();
		// std::cout << n_dofs << "sulla faccia" <<"\n";
		copy_data_face.joint_dof_indices = fe_iv.get_interface_dof_indices();

		copy_data_face.cell_matrix.reinit(n_dofs, n_dofs);

		const std::vector<double> &JxW = fe_iv.get_JxW_values();
		const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

		for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
		{
			const double beta_dot_n = beta(q_points[qpoint]) * normals[qpoint];
			for (unsigned int i = 0; i < n_dofs; ++i)
			{
				for (unsigned int j = 0; j < n_dofs; ++j)
				{
					copy_data_face.cell_matrix(i, j) += (beta(q_points[qpoint]) * normals[qpoint] * fe_iv.average(j, qpoint) * fe_iv.jump(i, qpoint) +
														 theta * std::abs(beta_dot_n) * fe_iv.jump(j, qpoint) * fe_iv.jump(i, qpoint)) *
														JxW[qpoint];
																				   // dx

					/*fe_iv.jump(i, qpoint)																		   // [\phi_i]
														*fe_iv.shape_value((beta_dot_n > 0), j, qpoint) // phi_j^{upwind}
														*JxW[qpoint];																				   // dx*/
				}
				// auto a  = fe_iv.shape_value(true, i,qpoint) - fe_iv.shape_value(false, i,qpoint);
				// std::cout << fe_iv.jump(i,qpoint) << "\t"<< a <<"\n"; //check that the jump do what it has to do!
			}
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

	const auto copier = [&](const CopyData &c)
	{
		constraints.distribute_local_to_global(c.cell_matrix, c.cell_rhs,
											   c.local_dof_indices, system_matrix, right_hand_side);

		for (auto &cdf : c.face_data)
		{
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
						  MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces | MeshWorker::assemble_own_interior_faces_once,
						  boundary_worker, face_worker);
}

template <int dim>
void AdvectionReaction<dim>::solve()
{

	if (use_direct_solver)
	{

		SparseDirectUMFPACK system_matrix_inverse;
		system_matrix_inverse.initialize(system_matrix);
		system_matrix_inverse.vmult(solution, right_hand_side);
	}
	else
	{
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
template <int dim>
void AdvectionReaction<dim>::refine_grid()
{

	if (refinement == "residual")
	{
		// for(auto& x: error_indicator_per_cell){
		// 	std::cout << x <<"\n";
		// }
		compute_local_projection_and_estimate();
		const double refinement_fraction = 0.3;
		GridRefinement::refine_and_coarsen_fixed_number(triangulation, error_indicator_per_cell, refinement_fraction, 0.1);
		triangulation.execute_coarsening_and_refinement();
		std::cout << error_indicator_per_cell.l2_norm() << "\n";
	}
	else if (refinement == "gradient")
	{

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
		GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
														  gradient_indicator, 0.25, 0.0);

		triangulation.execute_coarsening_and_refinement();
		std::cout << gradient_indicator.l2_norm() << "\n";
	}
	else if (refinement == "global")
	{
		triangulation.refine_global(1); //just for testing on uniformly refined meshes
	}
	else
	{
		Assert(false, ExcInternalError());
	}
}

// The output of this program consists of a vtk file of the adaptively
// refined grids and the numerical solutions. Finally, we also compute the
// L-infinity norm of the solution using VectorTools::integrate_difference().
template <int dim>
void AdvectionReaction<dim>::output_results(const unsigned int cycle) const
{
	const std::string filename = "solution-" + std::to_string(cycle) + ".vtk";
	std::cout << "  Writing solution to <" << filename << ">" << std::endl;
	std::ofstream output(filename);

	DataOut<dim> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, "u", DataOut<dim>::type_dof_data);

	data_out.build_patches(mapping);

	data_out.write_vtk(output);
}

template <int dim>
void AdvectionReaction<dim>::compute_error()
{
	error_table.error_from_exact(mapping, dof_handler, solution, exact_solution); //be careful: a FD approximation of the gradient is used to compute the H^1 norm if you're not relying on SymbolicFunction class
																				  //	error_table.error_from_exact(mapping, dof_handler, solution, Solution<dim>()); //provided that Solution<dim> implements the Gradient function
}

template <int dim>
double AdvectionReaction<dim>::compute_energy_norm()
{

	energy_norm_square_per_cell.reinit(triangulation.n_active_cells());

	using Iterator = typename DoFHandler<dim>::active_cell_iterator;
	//Working on a single cell K
	const auto cell_worker = [&](const Iterator &cell,
								 ScratchData<dim> &scratch_data, CopyData &copy_data)
	{
		const unsigned int n_dofs =
			scratch_data.fe_values.get_fe().n_dofs_per_cell();
		copy_data.reinit(cell, n_dofs);
		scratch_data.fe_values.reinit(cell);

		copy_data.cell_index = cell->active_cell_index();

		const auto &q_points = scratch_data.fe_values.get_quadrature_points();
		const FEValues<dim> &fe_v = scratch_data.fe_values;
		const std::vector<double> &JxW = fe_v.get_JxW_values();

		double error_square_norm{0.0};
		std::vector<double> sol_u(fe_v.n_quadrature_points);
		fe_v.get_function_values(solution, sol_u);

		for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
		{
			const double diff = (sol_u[point] - exact_solution.value(q_points[point]));
			error_square_norm += diff * diff * JxW[point];
		}
		copy_data.value = error_square_norm;
	};

	//Working on interior faces
	const auto face_worker = [&](const Iterator &cell,
								 const unsigned int &f,
								 const unsigned int &sf,
								 const Iterator &ncell,
								 const unsigned int &nf,
								 const unsigned int &nsf,
								 ScratchData<dim> &scratch_data,
								 CopyData &copy_data)
	{
		// const FEInterfaceValues<dim> &fe_iv = scratch_data.fe_iv.reinit(cell, f, sf, ncell, nf, nsf);

		FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
		fe_iv.reinit(cell, f, sf, ncell, nf, nsf);

		copy_data.face_data.emplace_back();
		CopyDataFace &copy_data_face = copy_data.face_data.back();
		copy_data_face.cell_indices[0] = cell->active_cell_index();
		copy_data_face.cell_indices[1] = ncell->active_cell_index();

		const auto &q_points = fe_iv.get_quadrature_points();
		const unsigned n_q_points = q_points.size();
		const std::vector<double> &JxW = fe_iv.get_JxW_values();
		std::vector<double> g(n_q_points);

		std::vector<double> jump(n_q_points);
		get_function_jump(fe_iv, solution, jump);

		const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

		double error_jump_square{0.0};
		for (unsigned int point = 0; point < n_q_points; ++point)
		{
			const double beta_dot_n = theta * std::abs(beta(q_points[point]) * normals[point]);
			error_jump_square += beta_dot_n * jump[point] * jump[point] * JxW[point];
		}

		copy_data.value = error_jump_square;
	};

	//Working on faces
	const auto boundary_worker = [&](const Iterator &cell,
									 const unsigned int &face_no,
									 ScratchData<dim> &scratch_data,
									 CopyData &copy_data)
	{
		scratch_data.fe_interface_values.reinit(cell, face_no);
		const FEFaceValuesBase<dim> &fe_fv = scratch_data.fe_interface_values.get_fe_face_values(0);
		const auto &q_points = fe_fv.get_quadrature_points();
		const unsigned n_q_points = q_points.size();
		const std::vector<double> &JxW = fe_fv.get_JxW_values();

		std::vector<double> g(n_q_points);
		// const BoundaryValues<dim> boundary_function;
		// boundary_function.value_list(q_points, g);

		std::vector<double> sol_u(n_q_points);
		fe_fv.get_function_values(solution, sol_u);

		const std::vector<Tensor<1, dim>> &normals = fe_fv.get_normal_vectors();

		double difference_norm_square = 0.;
		for (unsigned int point = 0; point < q_points.size(); ++point)
		{
			const double beta_dot_n = theta * std::abs(beta(q_points[point]) * normals[point]);
			const double diff = (boundary_conditions.value(q_points[point]) - sol_u[point]);
			difference_norm_square += beta_dot_n * diff * diff * JxW[point];
		}
		copy_data.value = difference_norm_square;
	};

	const auto copier = [&](const auto &copy_data)
	{
		if (copy_data.cell_index != numbers::invalid_unsigned_int)
		{
			// std::cout << copy_data.cell_index << "\n";
			energy_norm_square_per_cell[copy_data.cell_index] += copy_data.value;
		}
		for (auto &cdf : copy_data.face_data)
			for (unsigned int j = 0; j < 2; ++j)
				energy_norm_square_per_cell[cdf.cell_indices[j]] += cdf.values[j];
	};

	//   const UpdateFlags cell_flags =   update_gradients | update_quadrature_points | update_JxW_values;
	//   UpdateFlags face_flags = update_values | update_quadrature_points | update_JxW_values;
	// const QGauss<dim> quadrature_overintegration = fe->tensor_degree() + 1;
	// const QGauss<dim - 1> face_quadrature_overintegration = fe->tensor_degree() + 1;

	ScratchData<dim> scratch_data(mapping, *fe, QGauss<dim>{fe->tensor_degree() + 1},
								  QGauss<dim - 1>{fe->tensor_degree() + 1});

	CopyData copy_data;

	MeshWorker::mesh_loop(dof_handler.begin_active(),
						  dof_handler.end(),
						  cell_worker,
						  copier,
						  scratch_data,
						  copy_data,
						  MeshWorker::assemble_own_cells |
							  MeshWorker::assemble_own_interior_faces_once |
							  MeshWorker::assemble_boundary_faces,
						  boundary_worker,
						  face_worker);

	const double energy_error = std::sqrt(energy_norm_square_per_cell.l1_norm());
	return energy_error;
}

template <int dim>
void AdvectionReaction<dim>::compute_reconstruction_over_cell(const Iterator &cell)
{

	const QGauss<dim> quadrature = fe->tensor_degree() + 1;
	const QGauss<dim - 1> quadrature_face = fe->tensor_degree() + 1;
	std::cout << "sono entrato qui"
			  << "\n";

	// const unsigned int n_q_points      = quadrature.size();
	// const unsigned int n_face_q_points = quadrature_face.size();
	// const unsigned int dofs_per_cell = fe->n_dofs_per_cell();

	FEValues<dim> fe_values(*fe,
							quadrature,
							update_values | update_gradients |
								update_quadrature_points | update_JxW_values);

	std::cout << "sono entrato quo"
			  << "\n";
	FEInterfaceValues<dim> fe_iv(*fe,
								 quadrature_face,
								 update_values | update_quadrature_points |
									 update_normal_vectors |
									 update_JxW_values);

	// const FEFaceValuesBase<dim>& fe_fv = fe_iv.get_fe_face_values(0);

	fe_values.reinit(cell);
	const auto &q_points = fe_values.get_quadrature_points();
	std::cout << "sono entrato qua"
			  << "\n";
	const unsigned int n_q_points = q_points.size();
	const unsigned int dofs_per_cell = fe_values.dofs_per_cell;
	std::vector<Tensor<1, dim>> gradient_sol_u_at_quadrature_points(fe_values.n_quadrature_points);
	fe_values.get_function_gradients(solution, gradient_sol_u_at_quadrature_points);

	FullMatrix<double> cell_matrix_a(dofs_per_cell, dofs_per_cell);
	Vector<double> cell_rhs(dofs_per_cell);

	cell_matrix_a = 0.;
	cell_rhs = 0.;

	for (unsigned int point = 0; point < n_q_points; ++point)
	{
		for (unsigned int i = 0; i < dofs_per_cell; ++i)
		{
			for (unsigned int j = 0; j < dofs_per_cell; ++j)
			{
				cell_matrix_a(i, j) += beta(q_points[point]) *
									 fe_values.shape_grad(j, point) *
									 fe_values.shape_value(i, point) *
									 fe_values.JxW(point);
			}
			cell_rhs(i) += beta(q_points[point]) *
						   gradient_sol_u_at_quadrature_points[point] *
						   fe_values.shape_value(i, point) *
						   fe_values.JxW(point);
		}
	}

	for (const auto &face : cell->face_iterators())
	{
		if (!(face->at_boundary()))
		{ //not a boundary face
			std::cout << face->at_boundary() << "\n";
			fe_iv.reinit(cell, cell->face_iterator_to_index(face));
			const auto &q_points_face = fe_iv.get_quadrature_points();
			const unsigned n_face_q_points = q_points_face.size();
			const unsigned int n_dofs = fe_iv.n_current_interface_dofs();

			// std::vector<double> jump(n_face_q_points);
			// get_function_jump(fe_iv, solution, jump);

			const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();
			for (unsigned int point = 0; point < n_face_q_points; ++point)
			{
				const double beta_dot_n = beta(q_points[point]) * normals[point];
				for (unsigned int i = 0; i < n_dofs; ++i)
				{
					cell_rhs(i) += beta_dot_n *
								   // jump[point]*
								   fe_iv.shape_value(true, i, point) *
								   fe_iv.JxW(point);
				}
			}
		}
	}

	std::cout << "determinante è: " << cell_matrix_a.determinant() << "\n";
	// for(const auto& x : cell_matrix){
	// 	std::cout << x<<"\n";
	// }
	// std::cout << "\n";
	// FullMatrix<double> inverse(cell_matrix.m(), cell_matrix.m());
	// inverse.invert(cell_matrix);
	// Vector<double> u_reconstruct(cell_matrix.m());
	// inverse.vmult(u_reconstruct, cell_rhs); //A^{-1}*rhs = u_rec
}

// }

template <int dim>
void AdvectionReaction<dim>::compute_local_projection_and_estimate()
{

	//Working over one single cell
	using Iterator = typename DoFHandler<dim>::active_cell_iterator;
	error_indicator_per_cell.reinit(triangulation.n_active_cells());

	const auto cell_worker = [&](const Iterator &cell,
								 ScratchData<dim> &scratch_data, CopyData &copy_data)
	{
		const unsigned int n_dofs = scratch_data.fe_values.get_fe().n_dofs_per_cell();

		copy_data.reinit(cell, n_dofs);
		scratch_data.fe_values.reinit(cell);
		copy_data.cell_index = cell->active_cell_index();

		const auto &q_points = scratch_data.fe_values.get_quadrature_points();
		const unsigned n_q_points = q_points.size();

		const FEValues<dim> &fe_v = scratch_data.fe_values;
		const std::vector<double> &JxW = fe_v.get_JxW_values();

		std::vector<double> sol_u_at_quadrature_points(fe_v.n_quadrature_points);
		fe_v.get_function_values(solution, sol_u_at_quadrature_points);

		std::vector<Tensor<1, dim>> gradient_sol_u_at_quadrature_points(fe_v.n_quadrature_points);
		fe_v.get_function_gradients(solution, gradient_sol_u_at_quadrature_points);

		//Compute local L^2 projection of f-cU_h over the local finite element space
		for (unsigned int point = 0; point < n_q_points; ++point)
		{
			for (unsigned int i = 0; i < n_dofs; ++i)
			{
				for (unsigned int j = 0; j < n_dofs; ++j)
				{
					// std::cout << "Hey"<<"\n";
					copy_data.cell_mass_matrix(i, j) += fe_v.shape_value(i, point) * //phi_i(x_q)
														fe_v.shape_value(j, point) * //phi_j(x_q)
														JxW[point];					 // dx(x_q)

					// copy_data.cell_reconstruction_matrix(i, j) += beta(q_points[point]) *
					// 											  fe_v.shape_grad(j, point) *
					// 											  fe_v.shape_value(i, point) *
					// 											  JxW[point];
				}
				copy_data.cell_mass_rhs(i) +=
					(rhs.value(q_points[point])		  // f(x_q)
						 * fe_v.shape_value(i, point) //phi_i(x_q)
					 -
						
						advection_coeff.value(q_points[point])
					 	*fe_v.shape_value(i, point)				 //phi_i(x_q)
						 * sol_u_at_quadrature_points[point] //U_h(x_q)

					- beta(q_points[point])
					* gradient_sol_u_at_quadrature_points[point]
					*fe_v.shape_value(i,point)
					 ) *
					JxW[point]; //dx

				// copy_data.cell_reconstruction_rhs(i)+=beta(q_points[point])*
				// 										gradient_sol_u_at_quadrature_points[point]*
				// 										fe_v.shape_value(i,point)*
				// 										JxW[point];
			}
		}

		FullMatrix<double> inverse(fe_v.n_quadrature_points, fe_v.n_quadrature_points);
		inverse.invert(copy_data.cell_mass_matrix);
		Vector<double> proj(fe_v.n_quadrature_points); //projection of (f-c*U_h) on the local fe_space
		inverse.vmult(proj, copy_data.cell_mass_rhs);  //M^{-1}*rhs = proj
		// std::cout << "determinante è: " << copy_data.cell_reconstruction_matrix.determinant() << "\n";

		// std::cout << copy_data.cell_index <<"\n";
		// FullMatrix<double> test(fe_v.n_quadrature_points, fe_v.n_quadrature_points);
		// inverse.mmult(test,copy_data.cell_mass_matrix);
		// for(auto& x : test){
		// 	if(x>=1.2){std::cout << x <<"\n";};
		// }
		// std::cout << "determinante: "<<copy_data.cell_reconstruction_matrix.determinant() <<"\n";

		double square_norm_over_cell = 0.0;
		for (unsigned int point = 0; point < n_q_points; ++point)
		{
			const double diff = rhs.value(q_points[point]) - sol_u_at_quadrature_points[point] - beta(q_points[point])*gradient_sol_u_at_quadrature_points[point] - proj[point];
			square_norm_over_cell += diff * diff * JxW[point];
		}
		copy_data.value_estimator = square_norm_over_cell;
	};

	//Working over boundary face
	const auto boundary_worker = [&](const Iterator &cell,
									 const unsigned int &face_no,
									 ScratchData<dim> &scratch_data,
									 CopyData &copy_data)
	{
		scratch_data.fe_interface_values.reinit(cell, face_no);
		const FEFaceValuesBase<dim> &fe_fv = scratch_data.fe_interface_values.get_fe_face_values(0);
		const auto &q_points = fe_fv.get_quadrature_points();
		const unsigned n_q_points = q_points.size();
		const std::vector<double> &JxW = fe_fv.get_JxW_values();

		std::vector<double> g(n_q_points);
		exact_solution.value_list(q_points, g);

		std::vector<double> sol_u(n_q_points);
		fe_fv.get_function_values(solution, sol_u);

		const std::vector<Tensor<1, dim>> &normals = fe_fv.get_normal_vectors();

		double square_norm_over_bdary_face = 0.;
		for (unsigned int point = 0; point < q_points.size(); ++point)
		{
			const double beta_dot_n = beta(q_points[point]) * normals[point];

			if (beta_dot_n < 0) //\partial_{-T} \cap \partial_{- \Omega}
			{
				const double diff = std::abs(beta_dot_n) * (g[point] - sol_u[point]);
				square_norm_over_bdary_face += diff * diff * JxW[point];
			}
		}
		copy_data.value_estimator += square_norm_over_bdary_face;
	};

	//Working on interior faces
	const auto face_worker = [&](const Iterator &cell,
								 const unsigned int &f,
								 const unsigned int &sf,
								 const Iterator &ncell,
								 const unsigned int &nf,
								 const unsigned int &nsf,
								 ScratchData<dim> &scratch_data,
								 CopyData &copy_data)
	{
		FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
		fe_iv.reinit(cell, f, sf, ncell, nf, nsf);

		copy_data.face_data.emplace_back();
		CopyDataFace &copy_data_face = copy_data.face_data.back();
		copy_data_face.cell_indices[0] = cell->active_cell_index();
		copy_data_face.cell_indices[1] = ncell->active_cell_index();

		const auto &q_points = fe_iv.get_quadrature_points();
		const unsigned n_q_points = q_points.size();
		const unsigned int n_dofs = fe_iv.n_current_interface_dofs();

		copy_data.cell_reconstruction_rhs.reinit(n_dofs); //daghe n'ocio nel reinit di CopyData TODO

		const std::vector<double> &JxW = fe_iv.get_JxW_values();
		std::vector<double> g(n_q_points);

		std::vector<double> jump(n_q_points);
		get_function_jump(fe_iv, solution, jump);

		const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

		for (unsigned int point = 0; point < n_q_points; ++point)
		{
			const double beta_dot_n = beta(q_points[point]) * normals[point];
			for (unsigned int i = 0; i < n_dofs; ++i)
			{
				if (beta_dot_n < 0)
				{
					copy_data.cell_reconstruction_rhs(i) += beta_dot_n *
															jump[point] *
															fe_iv.shape_value(true, i, point) *
															JxW[point];
				}
			}
		}

		double error_jump_square{0.0};
		for (unsigned int point = 0; point < n_q_points; ++point)
		{
			const double beta_dot_n = beta(q_points[point]) * normals[point];
			if (beta_dot_n < 0)
			{
				error_jump_square += std::abs(beta_dot_n) * jump[point] * jump[point] * JxW[point];
			}
		}

		// copy_data.value = error_jump_square;
		copy_data_face.values[0] = error_jump_square;
		copy_data_face.values[1] = copy_data_face.values[0];
	};


	ScratchData<dim> scratch_data(mapping, *fe, QGauss<dim>{fe->tensor_degree() + 1},
								  QGauss<dim - 1>{fe->tensor_degree() + 1});

	const auto copier = [&](const auto &copy_data)
	{
		if (copy_data.cell_index != numbers::invalid_unsigned_int)
		{
			// std::cout << copy_data.cell_index << "\n";
			error_indicator_per_cell[copy_data.cell_index] += copy_data.value_estimator;
		}
		for (auto &cdf : copy_data.face_data)
		{
			for (unsigned int j = 0; j < 2; ++j)
			{
				error_indicator_per_cell[cdf.cell_indices[j]] += cdf.values[j];
			}
		}
	};

	// Here, we finally handle the assembly of the MASS matrix. We pass in ScratchData and
	// CopyData objects
	CopyData copy_data;
	MeshWorker::mesh_loop(dof_handler.begin_active(), dof_handler.end(),
						  cell_worker, copier, scratch_data, copy_data,
						  MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces | MeshWorker::assemble_own_interior_faces_once,
						  boundary_worker, face_worker);
}

//Usual run functions, running over several refinement cycles
template <int dim>
void AdvectionReaction<dim>::run()
{
	std::vector<double> energy_errors;
	std::vector<int> dofs_hist;
	for (unsigned int cycle = 0; cycle < n_refinement_cycles; ++cycle)
	{
		std::cout << "Cycle " << cycle << std::endl;

		if (cycle == 0)
		{
			GridGenerator::hyper_cube(triangulation);
			triangulation.refine_global(n_global_refinements);
		}
		else
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
		energy_errors.emplace_back(compute_energy_norm());
		dofs_hist.emplace_back(triangulation.n_active_cells());
	}

	error_table.output_table(std::cout);

	std::cout << "Error rate in energy norm"
			  << "\n";
	std::vector<double> rate_energy;
	for (unsigned int i = 0; i < n_refinement_cycles; ++i)
	{
		std::cout << "Cycle " << i << "\t" << energy_errors[i] << "\t";
		i >= 1 ? rate_energy.push_back(dim* std::log2(energy_errors[i - 1] / energy_errors[i]) / std::log2(dofs_hist[i] / dofs_hist[i - 1])) : rate_energy.push_back(.0);
		std::cout << rate_energy[i] << "\n";
	}
}

template class AdvectionReaction<2>;
template class AdvectionReaction<3>;
