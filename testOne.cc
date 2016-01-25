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

#include <deal.II/base/utilities.h>

/*
 * Purpose:
 *   1. Combine separated file written by parallel processors;
 *   2. Sort data according to coordinate;
 *   3. Compute local Mach number and pressure coefficient.
 */

template<int dim>
class DataProcessor
{
public:
  DataProcessor(const std::string &file_name_prefix_in);
private:
  void run() const;

  const std::string    &file_name_prefix;
  const unsigned int    solution_component = dim + 2;
};


int main(int argc, char *argv[])
{
  // Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  // unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("2d");

  deallog.pop();

  return (0);
}

template<int dim>
DataProcessor<dim>::DataProcessor(const std::string &file_name_prefix_in)
  :
  file_name_prefix(file_name_prefix_in)
{}

template<int dim>
void
DataProcessor<dim>::run() const
{
  return;
}
