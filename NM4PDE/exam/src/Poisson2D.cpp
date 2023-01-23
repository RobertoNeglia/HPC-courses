#include "Poisson2D.hpp"

void
Poisson2D::setup() {
  pcout << "===============================================" << std::endl;

  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    Triangulation<dim> mesh_serial;

    {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(mesh_serial);

      std::ifstream mesh_in_file("../mesh/mesh-square-h" + std::to_string(1.0 / (N + 1)) +
                                 ".msh");
      grid_in.read_msh(mesh_in_file);
    }

    // Each parallel process now owns a copy of the whole mesh
    // I need now to split the mesh between them

    {
      // create a partition of the triangulation
      GridTools::partition_triangulation(mpi_size, mesh_serial);

      // describe what the triangulation contains
      const auto construction_data =
        TriangulationDescription::Utilities::create_description_from_triangulation(
          mesh_serial, MPI_COMM_WORLD);

      // create the triangulation for each process, with only owned elements
      mesh.create_triangulation(construction_data);
    }

    pcout << "  Number of elements = " << mesh.n_global_active_cells() << std::endl;

    // no writing to a file (doing so in parallel is complicated)
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    pcout << "Initializing the finite element space" << std::endl;

    // Construct the finite element object. Notice that we use the FE_SimplexP
    // class here, that is suitable for triangular (or tetrahedral) meshes.
    fe = std::make_unique<FE_SimplexP<dim>>(r);

    pcout << "  Degree                     = " << fe->degree << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell << std::endl;

    // Construct the quadrature formula of the appopriate degree of exactness.
    quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size() << std::endl;

    quadrature_boundary = std::make_unique<QGaussSimplex<dim - 1>>(r + 1);

    pcout << "  Quadrature points per boundary cell = " << quadrature_boundary->size()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    // Initialize the DoF handler with the mesh we constructed.
    dof_handler.reinit(mesh);

    // "Distribute" the degrees of freedom. For a given finite element space,
    // initializes info on the control variables (how many they are, where
    // they are collocated, their "global indices", ...).
    dof_handler.distribute_dofs(*fe);

    pcout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;

    locally_owned_dofs = dof_handler.locally_owned_dofs();
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    // We first initialize a "sparsity pattern", i.e. a data structure that
    // indicates which entries of the matrix are zero and which are different
    // from zero. To do so, we construct first a DynamicSparsityPattern (a
    // sparsity pattern stored in a memory- and access-inefficient way, but
    // fast to write) and then convert it to a SparsityPattern (which is more
    // efficient, but cannot be modified).
    pcout << "  Initializing the sparsity pattern" << std::endl;
    TrilinosWrappers::SparsityPattern sparsity_pattern(locally_owned_dofs,
                                                       MPI_COMM_WORLD);

    DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
    sparsity_pattern.compress();

    // Then, we use the sparsity pattern to initialize the system matrix
    pcout << "  Initializing the system matrix" << std::endl;
    system_matrix.reinit(sparsity_pattern);

    // Finally, we initialize the right-hand side and solution vectors.
    pcout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  }
}

void
Poisson2D::assemble() {
  pcout << "===============================================" << std::endl;

  pcout << "  Assembling the linear system" << std::endl;

  // Number of local DoFs for each element.
  const unsigned int dofs_per_cell = fe->dofs_per_cell;

  // Number of quadrature points for each element.
  const unsigned int n_q = quadrature->size();

  // FEValues instance. This object allows to compute basis functions, their
  // derivatives, the reference-to-current element mapping and its
  // derivatives on all quadrature points of all elements.
  FEValues<dim> fe_values(
    *fe,
    *quadrature,
    // Here we specify what quantities we need FEValues to compute on
    // quadrature points. For our test, we need:
    // - the values of shape functions (update_values);
    // - the derivative of shape functions (update_gradients);
    // - the position of quadrature points (update_quadrature_points);
    // - the product J_c(x_q)*w_q (update_JxW_values).
    update_values | update_gradients | update_quadrature_points | update_JxW_values);

  // Since we need to compute integrals on the boundary for Neumann conditions,
  // we also need a FEValues object to compute quantities on boundary edges
  // (faces).
  FEFaceValues<dim> fe_values_boundary(*fe,
                                       *quadrature_boundary,
                                       update_values | update_quadrature_points |
                                         update_JxW_values);

  // Local matrix and right-hand side vector. We will overwrite them for
  // each element within the loop.
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  // We will use this vector to store the global indices of the DoFs of the
  // current element within the loop.
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  // Reset the global matrix and vector, just in case.
  system_matrix = 0.0;
  system_rhs    = 0.0;

    for (const auto &cell : dof_handler.active_cell_iterators()) {
      if (!cell->is_locally_owned())
        continue;

      // Reinitialize the FEValues object on current element. This
      // precomputes all the quantities we requested when constructing
      // FEValues (see the update_* flags above) for all quadrature nodes of
      // the current cell.
      fe_values.reinit(cell);

      // We reset the cell matrix and vector (discarding any leftovers from
      // previous element).
      cell_matrix = 0.0;
      cell_rhs    = 0.0;

        for (unsigned int q = 0; q < n_q; ++q) {
            // Here we assemble the local contribution for current cell and
            // current quadrature point, filling the local matrix and vector.

            // Here we iterate over *local* DoF indices.
            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                  // FEValues::shape_grad(i, q) returns the gradient of the i-th
                  // basis function at the q-th quadrature node, already mapped
                  // on the physical element: we don't have to deal with the
                  // mapping, it's all hidden inside FEValues.
                  // Diffusion term
                  cell_matrix(i, j) +=
                    function_mu.value(fe_values.quadrature_point(q)) // mu(x)
                    * fe_values.shape_grad(i, q)                     // grad(phi_i)
                    * fe_values.shape_grad(j, q)                     // grad(phi_j)
                    * fe_values.JxW(q);                              // dx

                  // Reaction term
                  cell_matrix(i, j) +=
                    function_sigma.value(fe_values.quadrature_point(q)) // sigma(x)
                    * fe_values.shape_value(i, q)                       // phi_i
                    * fe_values.shape_value(j, q)                       // phi_j
                    * fe_values.JxW(q);                                 // dx
                }

              cell_rhs(i) += forcing_term.value(fe_values.quadrature_point(q)) // f
                           * fe_values.shape_value(i, q)                       // phi_i
                           * fe_values.JxW(q);                                 // dx
            }
        }

      // Neumann boundary condition
      // check if the element is on the boundary
        if (cell->at_boundary()) {
            // the element is on the boundary, but i dont know which one
            for (unsigned int face_number = 0; face_number < cell->n_faces();
                 face_number++) {
                if (cell->face(face_number)->at_boundary() &&
                    (cell->face(face_number)->boundary_id() == 1 ||
                     cell->face(face_number)->boundary_id() == 3)) {
                  // I'm on a boundary in which I must impose Neumann conditions
                  fe_values_boundary.reinit(cell, face_number);

                    for (unsigned int q = 0; q < quadrature_boundary->size();
                         q++) { // loop on quadrature nodes
                        for (unsigned int i = 0; i < dofs_per_cell; i++) {
                          cell_rhs(i) += function_q.value(
                                           fe_values_boundary.quadrature_point(q)) // q(x)
                                       * fe_values_boundary.shape_value(i, q) // phi_i(x)
                                       * fe_values_boundary.JxW(q); // Jc wq ( = dx)
                        }
                    }
              }
            }
      }

      // At this point the local matrix and vector are constructed: we
      // need to sum them into the global matrix and vector. To this end,
      // we need to retrieve the global indices of the DoFs of current
      // cell.
      cell->get_dof_indices(dof_indices);

      // Then, we add the local matrix and vector into the corresponding
      // positions of the global matrix and vector.
      system_matrix.add(dof_indices, cell_matrix);
      system_rhs.add(dof_indices, cell_rhs);
    }

  // communication between processes (for interfaces between partitions)
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

  // Dirichlet boundary conditions.
  {
    // We construct a map that stores, for each DoF corresponding to a
    // Dirichlet condition, the corresponding value. E.g., if the Dirichlet
    // condition is u_i = b, the map will contain the pair (i, b).
    std::map<types::global_dof_index, double> boundary_values;

    Functions::ZeroFunction<dim> zero_function;

    // Then, we build a map that, for each boundary tag, stores the
    // corresponding boundary function.
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;
    boundary_functions[0] = &zero_function;
    boundary_functions[2] = &zero_function;

    // interpolate_boundary_values fills the boundary_values map.
    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values);

    // Finally, we modify the linear system to apply the boundary
    // conditions. This replaces the equations for the boundary DoFs with
    // the corresponding u_i = 0 equations.
    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution, system_rhs, true);
  }
}

void
Poisson2D::solve() {
  pcout << "===============================================" << std::endl;

  // Here we specify the maximum number of iterations of the iterative solver,
  // and its tolerance.
  SolverControl solver_control(1000, 1e-6 * system_rhs.l2_norm());

  // Since the system matrix is symmetric and positive definite, we solve the
  // system using the conjugate gradient method.
  SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);

  TrilinosWrappers::PreconditionSSOR preconditioner;
  preconditioner.initialize(system_matrix,
                            TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));

  pcout << "  Solving the linear system" << std::endl;
  // We don't use any preconditioner for now, so we pass the identity matrix
  // as preconditioner.
  solver.solve(system_matrix, solution, system_rhs, preconditioner);
  pcout << "  " << solver_control.last_step() << " CG iterations" << std::endl;
}

void
Poisson2D::output() const {
  pcout << "===============================================" << std::endl;

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  TrilinosWrappers::MPI::Vector solution_ghosted(locally_owned_dofs,
                                                 locally_relevant_dofs,
                                                 MPI_COMM_WORLD);

  // communication of ghost cells between processes
  solution_ghosted = solution;

  // The DataOut class manages writing the results to a file.
  DataOut<dim> data_out;

  // It can write multiple variables (defined on the same mesh) to a single
  // file. Each of them can be added by calling add_data_vector, passing the
  // associated DoFHandler and a name.
  data_out.add_data_vector(dof_handler, solution_ghosted, "solution");

  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());

  data_out.add_data_vector(partitioning, "partitioning");

  // Once all vectors have been inserted, call build_patches to finalize the
  // DataOut object, preparing it for writing to file.
  data_out.build_patches();

  // no parallel output support on vtk files

  const std::string output_file_name = "output-h" + std::to_string(1.0 / (N + 1));

  DataOutBase::DataOutFilter data_filter(DataOutBase::DataOutFilterFlags(false, true));
  data_out.write_filtered_data(data_filter);
  data_out.write_hdf5_parallel(data_filter, output_file_name + ".h5", MPI_COMM_WORLD);

  std::vector<XDMFEntry> xdmf_entries({data_out.create_xdmf_entry(
    data_filter, output_file_name + ".h5", 0.0, MPI_COMM_WORLD)});
  data_out.write_xdmf_file(xdmf_entries, output_file_name + ".xdmf", MPI_COMM_WORLD);

  pcout << "Output written to " << output_file_name << std::endl;

  pcout << "===============================================" << std::endl;
}

double
Poisson2D::compute_error(const VectorTools::NormType &norm_type) const {
  FE_SimplexP<dim> fe_linear(1);
  MappingFE        mapping(fe_linear);

  // In general it is advisable to not use the same quadrature formula used to compute the
  // solution of the system, but to use a more precise formula for the error (in this case
  // with one more node wrt to the solution quadrature formula)
  const QGaussSimplex<dim> quadrature_error(r + 2);
  Vector<double>           error_per_cell(mesh.n_active_cells());

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    exact_solution,
                                    error_per_cell,
                                    quadrature_error,
                                    norm_type);

  const double error = VectorTools::compute_global_error(mesh, error_per_cell, norm_type);

  return error;
}