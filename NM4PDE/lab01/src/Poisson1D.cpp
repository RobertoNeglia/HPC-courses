#include "Poisson1D.hpp"

void
Poisson1D::setup() {
    // build the mesh
    GridGenerator::subdivided_hyper_cube(mesh, N + 1, 0.0, 1.0, true);
    std::cout << "Number of mesh elements: " << mesh.n_active_cells() << std::endl;
    {
        // write the mesh to a file
        const std::string mesh_file_name = "mesh-" + std::to_string(N + 1) + ".vtk";
        std::ofstream grid_out_file(mesh_file_name);
        GridOut grid_out;
        grid_out.write_vtk(mesh, grid_out_file);
    }

    // finite elements
    {
        fe = std::make_unique<FE_Q<dim>>(r);
        quadrature = std::make_unique<QGauss<dim>>(r + 1);
    }

    // initialize dof (degree of freedom -> these are the control variables, the u of the solution) handler
    {
        dof_handler.reinit(mesh);
        dof_handler.distribute_dofs(*fe);

        std::cout << "Number of DoFs: " << dof_handler.n_dofs() << std::endl;
    }

    // build linear algebra structure (just allocation)
    {
        DynamicSparsityPattern dsp(dof_handler.n_dofs());
        DoFTools::make_sparsity_pattern(dof_handler, dsp);
        sparsity_pattern.copy_from(dsp);

        system_matrix.reinit(sparsity_pattern);

        system_rhs.reinit(dof_handler.n_dofs());
        solution.reinit(dof_handler.n_dofs());
    }
}

void
Poisson1D::assemble() {
    std::cout << "Assembling the linear system" << std::endl;

    system_matrix = 0.0;
    system_rhs = 0.0;

    const unsigned int dofs_per_cell = fe->dofs_per_cell;
    FullMatrix<double> cell_matrix(dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    const unsigned int n_q_points = quadrature->size();

    FEValues<dim> fe_values(*fe, *quadrature, update_values | update_gradients | update_quadrature_points | update_JxW_values);



    // for (const auto& cell; dof_handler.active_cell_iterators()) {
    //     cell_matrix = 0.0;
    //     cell_rhs = 0.0;

    //     fe_values.reinit(cell);

    //     for (unsigned int q = 0; q < n_q_points; q++) {
    //         for (unsigned int i = 0; i < dofs_per_cell; i++) {
    //             for (unsigned int j = 0;j < dofs_per_cell; j++) {
    //                 cell_matrix(i, j) = ...;

    //             }

    //             cell_rhs(i) = ...;
    //         }
    //     }
    // }

}

void
Poisson1D::solve() {
}

void
Poisson1D::output() const {
}