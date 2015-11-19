// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/lac/generic_linear_algebra.h>

//#define USE_PETSC_LA
#ifndef USE_PETSC_LA
//#define USE_TRILINOS_SPARSITY_PATTERN
#endif
#include <fstream>
#include <ostream>

template<int dim>
void test()
{
  // Declare common variables
  MPI_Comm                                  mpi_communicator(MPI_COMM_WORLD);
  parallel::distributed::Triangulation<dim> triangulation(mpi_communicator);
  IndexSet                                  locally_owned_dofs;
  IndexSet                                  locally_relevant_dofs;
  const unsigned int   fe_degree    = 1;
  const unsigned int   n_components = 1; //dim+2;
  const FESystem<dim>  fe (FE_Q<dim>(fe_degree), n_components);
  DoFHandler<dim>      dof_handler(triangulation);

#ifdef USE_PETSC_LA
  typedef dealii::LinearAlgebraPETSc::MPI::SparseMatrix    SparseMatrixType;
  std::cout << "USE_PETSC_LA" << std::endl;
#else
  typedef dealii::LinearAlgebraTrilinos::MPI::SparseMatrix SparseMatrixType;
  std::cout << "USE_Trilinos_LA" << std::endl;
#endif
  SparseMatrixType             system_matrix;

  const unsigned int myid(Utilities::MPI::this_mpi_process (mpi_communicator));
  const bool I_am_host(myid == 0);
  unsigned int field_output_counter = 0;

  // Create mesh
  {
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(1);
    triangulation.begin_active()->set_refine_flag();
    triangulation.execute_coarsening_and_refinement();
  }

  // Write out mesh
  {
    DataOut<dim> data_out;
    data_out.attach_triangulation (triangulation);
    Vector<float> subdomain (triangulation.n_active_cells());
    std::fill (subdomain.begin(), subdomain.end(), myid);
    data_out.add_data_vector (subdomain, "subdomain");

    Vector<float> cell_index (triangulation.n_active_cells());
    {
      typename Triangulation<dim>::active_cell_iterator
      cell = triangulation.begin_active(),
      endc = triangulation.end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            cell_index[cell->active_cell_index()] = cell->active_cell_index();
          }
    }
    data_out.add_data_vector (cell_index, "cell_index");


    data_out.build_patches();

    const std::string output_tag = "solution-" +
                                   Utilities::int_to_string (field_output_counter, 4);
    const std::string slot_itag = ".slot-" + Utilities::int_to_string (myid, 4);

    std::ofstream output ((output_tag + slot_itag + ".vtu").c_str());
    data_out.write_vtu (output);

    if (I_am_host)
      {
        std::vector<std::string> filenames;
        for (unsigned int i=0;
             i<Utilities::MPI::n_mpi_processes (mpi_communicator);
             ++i)
          {
            filenames.push_back (output_tag +
                                 ".slot-" +
                                 Utilities::int_to_string (i, 4) +
                                 ".vtu");
          }
        std::ofstream master_output ((output_tag + ".pvtu").c_str());
        data_out.write_pvtu_record (master_output, filenames);
      }
  }

  // Setup system
  {
    dof_handler.clear();
    dof_handler.distribute_dofs (fe);
    locally_owned_dofs.clear();
    locally_owned_dofs = dof_handler.locally_owned_dofs();

    locally_relevant_dofs.clear();
    DoFTools::extract_locally_relevant_dofs (dof_handler,
                                             locally_relevant_dofs);

#ifdef USE_TRILINOS_SPARSITY_PATTERN
    TrilinosWrappers::SparsityPattern sparsity_pattern (locally_owned_dofs,
                                                        mpi_communicator);
    std::cout << "USE_TRILINOS_SPARSITY_PATTERN" << std::endl;
#else
    DynamicSparsityPattern sparsity_pattern (locally_relevant_dofs);
    std::cout << "USE_DynamicSparsityPattern" << std::endl;
#endif
    DoFTools::make_sparsity_pattern (dof_handler,
                                     sparsity_pattern,
                                     /*const ConstraintMatrix constraints = */ ConstraintMatrix(),
                                     /*const bool keep_constrained_dofs = */ true,
                                     myid);
#ifdef USE_TRILINOS_SPARSITY_PATTERN
    sparsity_pattern.compress();
    system_matrix.reinit (sparsity_pattern);
#else
    SparsityTools::distribute_sparsity_pattern (sparsity_pattern,
                                                dof_handler.n_locally_owned_dofs_per_processor(),
                                                mpi_communicator,
                                                locally_relevant_dofs);
    system_matrix.reinit (locally_owned_dofs,
                          locally_owned_dofs,
                          sparsity_pattern,
                          mpi_communicator);
#endif
  }

  // Pretend to assemble system, just write something into system matrix
  {
    const unsigned int dofs_per_cell = fe.dofs_per_cell;

    std::vector<types::global_dof_index> dof_indices (dofs_per_cell);
    std::vector<types::global_dof_index> dof_indices_neighbor (dofs_per_cell);

    std::vector<double> fake_etries_this(dofs_per_cell, 3.1416);
    std::vector<double> fake_etries_neighbor(dofs_per_cell, 2.71828);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices (dof_indices);
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              system_matrix.add (dof_indices[i], dofs_per_cell,
                                 & (dof_indices[0]), & (fake_etries_this[0]));
            }
          // Get dof_indices on neighboring cell.
          // We have several different cases to handle.
          for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
            {
              const bool face_is_at_boundary = cell->at_boundary (face_no);

              // Loop through all sub components of current face and extract neighboring DoFs.
              // The following statement means, if the face has children, loop through all of
              // its children; if not, treat itself as its only child.
              const unsigned int n_subfaces =
                cell->face (face_no)->has_children() ?
                cell->face (face_no)->n_children() : 1;
              for (unsigned int subface_no=0; subface_no < n_subfaces; ++subface_no)
                {
                  if (face_is_at_boundary)
                    {
                      // No neighboring cell exists
                      std::fill (dof_indices_neighbor.begin(),
                                 dof_indices_neighbor.end(),
                                 numbers::invalid_unsigned_int);
                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        {
                          system_matrix.add (dof_indices[i], dofs_per_cell,
                                             & (dof_indices[0]), & (fake_etries_this[0]));
                        }
                    }
                  else
                    {
                      const typename DoFHandler<dim>::cell_iterator
                      neighbor = cell->neighbor (face_no);

                      unsigned int neighbor_face_no = numbers::invalid_unsigned_int;
                      // Make sure not subtracting two unsigned integers.
                      const int this_cell_level = cell->level();
                      const int neighbor_active_cell_level = neighbor->level() + neighbor->has_children();
                      switch (this_cell_level-neighbor_active_cell_level)
                        {
                        case -1:
                        {
                          // This cell is one level coarser than neighbor. So we have several subfaces
                          // to deal with. The outer loop will take care this.
                          neighbor_face_no = cell->neighbor_of_neighbor(face_no);
                          cell->neighbor_child_on_subface (face_no, subface_no)->get_dof_indices (dof_indices_neighbor);
                          break;
                        }
                        case 0:
                        {
                          // This cell and neighbor cell are at same level
                          // Nothing to do.
                          break;
                        }
                        case 1:
                        {
                          // This cell is one level finer than neighbor. So we are now
                          // on a subface of neighbor cell.
                          neighbor_face_no = cell->neighbor_of_coarser_neighbor (face_no).first;
                          neighbor->get_dof_indices (dof_indices_neighbor);
                          break;
                        }
                        default:
                        {
                          Assert (false,
                                  ExcMessage ("Refinement level difference between face neighbor cells can't greater than 1."));
                          break;
                        }
                        } // End switch level difference

                      // Only take care of faces with hanging nodes.
                      if (this_cell_level != neighbor_active_cell_level)
                        {
                          deallog << cell->active_cell_index() << ", "
                                  << face_no <<  ", "
                                  << subface_no << std::endl;

                          for (unsigned int i=0; i<dofs_per_cell; ++i)
                            if (fe.has_support_on_face (i, face_no) == true)
                              {
                                system_matrix.add (dof_indices[i], dofs_per_cell,
                                                   & (dof_indices[0]), & (fake_etries_this[0]));
                                std::vector<types::global_dof_index> neighbor_effective_dof_indices;
                                std::vector<double> neighbor_effective_values;
                                for (unsigned int j=0; j<dofs_per_cell; ++j)
                                  if (fe.has_support_on_face (j, neighbor_face_no))
                                    {
                                      neighbor_effective_dof_indices.push_back(dof_indices_neighbor[j]);
                                      neighbor_effective_values.push_back(fake_etries_neighbor[j]);
                                    }
                                system_matrix.add (dof_indices[i], neighbor_effective_dof_indices.size(),
                                                   & (neighbor_effective_dof_indices[0]), & (neighbor_effective_values[0]));
                              }
                        }
                    } // End if (face_is_at_boundary) ... else.
                } // End for all subfaces
            } // End for all faces
        } // End for all locally owned active cells

    system_matrix.compress (VectorOperation::add);
  } // End assemble system,

  return;
}

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  // unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("2d");
  test<2>();
  deallog.pop();

  // deallog.push("3d");
  // test<3>();
  // deallog.pop();

  return (0);
}
