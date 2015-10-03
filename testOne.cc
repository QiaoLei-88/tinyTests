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



// recursively refine a 2d mesh

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/utilities.h>

#include <fstream>
#include <ostream>
#include <unistd.h>

template<int dim>

void test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  if (true)
    {
      if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
        deallog << "hyper_cube" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
      //dealii::Triangulation<dim>::none
      //limit_level_difference_at_vertices

      GridGenerator::hyper_cube(tr);
      tr.refine_global(1);

      unsigned int level = 0;
      while (tr.n_global_active_cells() < 20000/Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
        {
          if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
            {
              deallog << "refine_loop..." << std::endl;
              std::cerr <<"level: " << level
                        <<", n_global_active_cells: " << tr.n_global_active_cells()
                        << std::endl;
            }
          usleep(20);
          MPI_Barrier(MPI_COMM_WORLD);
          for (unsigned int my_id = 0; my_id<Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); ++my_id)
            {
              MPI_Barrier(MPI_COMM_WORLD);
              if (my_id == Utilities::MPI::this_mpi_process (MPI_COMM_WORLD))
                {
                  std::cerr <<"level: " << level
                            << ", myid: " << Utilities::MPI::this_mpi_process (MPI_COMM_WORLD)
                            << ", n_active_cells: " << tr.n_active_cells() << std::endl;
                  usleep(20);
                }
              MPI_Barrier(MPI_COMM_WORLD);
            }
          std::vector<bool> flags (tr.n_active_cells(), false);
          // refine one fifth of all cells each
          // time (but at least one).
          // note that only the own marked cells
          // will be refined.
          for (unsigned int i=0; i<tr.n_active_cells() / 5 + 1; ++i)
            {
              const unsigned int x = Testing::rand() % flags.size();
              flags[x] = true;
            }

          unsigned int index=0;
          for (typename Triangulation<dim>::active_cell_iterator
               cell = tr.begin_active();
               cell != tr.end(); ++cell, ++index)
            if (flags[index] && cell->is_locally_owned())
              {
                cell->set_refine_flag();
              }

          Assert (index == tr.n_active_cells(), ExcInternalError());
          tr.execute_coarsening_and_refinement ();

          if (myid == 0)
            {
              deallog << "#cells = " << tr.n_global_active_cells() << std::endl;
            }
          ++level;
        }
    }

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();

}
