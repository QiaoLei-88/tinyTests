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

#include <fstream>
#include <ostream>

template<int dim>


int main(int argc, char *argv[])
{
  // Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  // unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("2d");
  test<2>();
  deallog.pop();

  return (0);
}
