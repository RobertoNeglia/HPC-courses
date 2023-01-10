#include "Poisson1D.hpp"

// Preprocessing step
void
Poisson1D::setup() {
  std::cout << "==============================================" << std::endl;

  // Build the mesh
  {
    std::cout << "Initializing the mesh" << std::endl;
    GridGenerator::subdivided_hyper_cube(mesh, N + 1, 0.0, 1.0, true);
    std::cout << "  Number of elements = " << mesh.n_active_cells() << std::endl;

    // Write the mesh to a file
    const std::string mesh_file_name = "mesh-" + std::to_string(N + 1) + ".vtk";
    GridOut           grid_out;
    std::ofstream     grid_out_file(mesh_file_name);
    grid_out.write_vtk(mesh, grid_out_file);
    std::cout << "  Mesh saved to " << mesh_file_name << std::endl;
  }

  std::cout << "-----------------------------------------" << std::endl;

  // Finite element space
  {
    std::cout << "Initializing the finite element space" << std::endl;

    // 1D FE are obtained with the FE_Q class (on higher dimensions, these represents
    // hexahedral finite elements: tetrahedral elements are obtained with FE_SimplexP)
    fe = std::make_unique<FE_Q<dim>>(r); // FE_Q: finite elements on squares

    std::cout << "  Degree                     = " << fe->degree << std::endl;
    std::cout << "  DoFs per cell              = " << fe->dofs_per_cell << std::endl;

    // Construct the quadrature formula
    quadrature = std::make_unique<QGauss<dim>>(r + 1); // Gauss-Legendre quadrature

    std::cout << "  Quadrature points per cell = " << quadrature->size() << std::endl;
  }

  std::cout << "-----------------------------------------" << std::endl;

  // Initialize the DoF handler
  // (Degree of Freedom, these are the control variables of the solution, u_j)
  {
    std::cout << "Initializing the DoF handler" << std::endl;
    dof_handler.reinit(mesh);

    // "Distribute" DoF
    // For a given FE space, initializes info on control variables (how many they are,
    // where they are, "global indices", etc)
    dof_handler.distribute_dofs(*fe);

    // Print number of unknowns
    std::cout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  std::cout << "-----------------------------------------" << std::endl;

  // Initialize linear algebra structures (system matrix, rhs, solution)
  {
    std::cout << "Initializing the linear system" << std::endl;

    std::cout << "  Initializing the sparsity pattern" << std::endl;
    // Sparsity pattern: data structure that indicates which entries of the matrix are
    // zero and which are not

    // DynamicSparsityPattern: this is stored in a memory- and access-inefficient way, but
    // fast to write
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    // Then the DSP is converted to a SparsityPattern (efficient, but can't be modified)
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    std::cout << "  Initializing the system matrix" << std::endl;
    // Use the sparsity pattern to initialize the matrix
    system_matrix.reinit(sparsity_pattern);

    std::cout << "  Initializing the system rhs" << std::endl;
    system_rhs.reinit(dof_handler.n_dofs());
    std::cout << "  Initializing the system solution" << std::endl;
    system_solution.reinit(dof_handler.n_dofs());
  }
}

void
Poisson1D::assemble() {
  std::cout << "==============================================" << std::endl;
  std::cout << "  Assembling the linear system" << std::endl;

  // This is n_loc (number of local DoFs for each element)
  const unsigned int dofs_per_cell = fe->dofs_per_cell;

  // Number of quadrature points for each elements
  const unsigned int n_q_points = quadrature->size();

  // FEValues object allows to compute basis functions, their derivatives, the
  // reference-to-current element mapping and its derivatives on all quadrature points of
  // all elements
  FEValues fe_values(*fe,
                     *quadrature,
                     // Specify what quantities we need FEValues to compute on quadrature points:
                     // - values of shape functions
                     // - derivatives of shape functions
                     // - position of quadrature points
                     // - product J_c(x_q)*w_q
                     update_values | update_gradients | update_quadrature_points |
                       update_JxW_values);

  // Local matrix and rhs vector that will be computer for each elements
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  // Vector to store global indices of the DoFs of the current element within the loop
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  // Reset global matrix and rhs
  system_matrix = 0.0;
  system_rhs    = 0.0;

    for (const auto &cell : dof_handler.active_cell_iterators()) { // loop over all elements
      // Reinitialize FEValues object on current element
      // Precomputes all quantities requested during construction (the flags) for all
      // quadrature nodes of the current element
      fe_values.reinit(cell);

      // Reset local matrix and rhs
      cell_matrix = 0.0;
      cell_rhs    = 0.0;

        for (unsigned int q = 0; q < n_q_points; q++) { // loop over quadrature points
            // Build the local matrix and local rhs (local contribution) of current cell
            // and current quadrature point
            for (unsigned int i = 0; i < dofs_per_cell; i++) {     // local rows
                for (unsigned int j = 0; j < dofs_per_cell; j++) { // local cols
                  // ::shape_grad(i,q) returns the gradient of the i-th basis function at
                  // the q-th quadrature node, already mapped on the physical element
                  // RMK: DONT HAVE TO DEAL WITH THE MAPPING, IS ALL HIDDEN INSIDE
                  // FEValues
                  cell_matrix(i, j) += diffusion_term.value(fe_values.quadrature_point(q)) // mu(x)
                                     * fe_values.shape_grad(i, q)                          // (I)
                                     * fe_values.shape_grad(j, q)                          // (II)
                                     * fe_values.JxW(q);                                   // (III)
                }

              cell_rhs(i) += forcing_term.value(fe_values.quadrature_point(q)) *
                             fe_values.shape_value(i, q) * fe_values.JxW(q);
            }
        }
      // Local matrix and vector are constructed

      // Retrieve the global indices of the DoFs of the current cell
      cell->get_dof_indices(dof_indices);

      // Sum them in the global matrix and vector
      system_matrix.add(dof_indices, cell_matrix);
      system_rhs.add(dof_indices, cell_rhs);
    }

  // A_hat is assembled (N_h + 2 x N_h + 2)

  // Boundary conditions
  {
    // Map that stores for each DoF corresponding to a Dirichlet bc the corresponding value
    // If u_i = b, does the mapping i -> b
    std::map<types::global_dof_index, double> boundary_values;

    // Boundary data function (always zero)
    Functions::ZeroFunction<dim> bc_function;

    // Boundary ids:
    //    - left boundary = 0
    //    - right boundary = 1
    // For each boundary id, store the corresponding boundary function
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;
    boundary_functions[0] = &bc_function;
    boundary_functions[1] = &bc_function;
    // RMK: 0 and 1 ARE NOT the coordinates of the extrema, but they're just ids (see
    // comment above)

    // this fills the boundary_values map
    VectorTools::interpolate_boundary_values(dof_handler, boundary_functions, boundary_values);

    // Modify the linear system by applying boundary conditions
    // Replace the equations for the boundary DoFs with the corresponding u_i = 0
    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, system_solution, system_rhs, true);
  }
}

void
Poisson1D::solve() {
  std::cout << "==============================================" << std::endl;

  // relative stopping criterion
  double        tol    = 1e-6 * system_rhs.l2_norm();
  unsigned int  max_it = 1000;
  SolverControl solver_control(max_it, tol);

  // System matrix is SPD, use Conjugate Gradient
  SolverCG<Vector<double>> solver(solver_control);

  std::cout << "  Solving the linear system" << std::endl;
  // Do not use any preconditioner, use the identity
  solver.solve(system_matrix, system_solution, system_rhs, PreconditionIdentity());
  std::cout << solver_control.last_step() << "  CG iterations" << std::endl;
}

void
Poisson1D::output() const {
  std::cout << "==============================================" << std::endl;

  // Manage writing the result to a file
  DataOut<dim> data_out;

  // Multiple variables can be written (if defined on the same mesh) to a single file
  // Each can be added eith add_data_vector method, passing the DoF handler and a name
  data_out.add_data_vector(dof_handler, system_solution, "solution");

  // Once all vectors have been inserted, this finalize the DataOut object
  data_out.build_patches();

  // Write output file in a vtk format (there are many more)
  const std::string out_file_name = "output-" + std::to_string(N + 1) + ".vtk";
  std::ofstream     output_file(out_file_name);
  data_out.write_vtk(output_file);

  std::cout << "Output written to " << out_file_name << std::endl;
  std::cout << "==============================================" << std::endl;
}