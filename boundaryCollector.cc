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
#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/base/point.h>
#include <boost/container/map.hpp>
/*
 * Purpose:
 *   1. Combine separated file written by parallel processors;
 *   2. Sort data according to coordinate;
 *   3. Compute local Mach number and pressure coefficient.
 */

template<int dim>
class DataProcessor
{
private:
  static const unsigned int n_components             = dim+2;
  static const unsigned int first_momentum_component = 0;
  static const unsigned int first_velocity_component = 0;
  static const unsigned int density_component        = dim;
  static const unsigned int energy_component         = dim+1;
  static const unsigned int pressure_component       = dim+1;

  struct DataStruct
  {
    Point<dim> position;
    std_cxx11::array<double, n_components> solution;
    double velocity_magnitude;
    double sound_speed;
    double local_Mach;
    double pressure_coefficient;
  };

  typedef boost::container::map<double, DataStruct> DataMap;

public:
  DataProcessor();

private:
  void run() const;

  const double gas_gamm;
};


int main(int argc, char *argv[])
{
  DataProcessor<2> data_processor;
  return (0);
}

template<int dim>
DataProcessor<dim>::DataProcessor()
  :
  gas_gamm(1.4)
{
  run();
}

template<int dim>
void
DataProcessor<dim>::run() const
{
  // Declare data storing object
  DataMap data_map;
  data_map.clear();
  double Mach = -1.0;

  // Gathering data written by all slots
  for ( unsigned int i_slot=0; i_slot<10000; ++i_slot)
    {
      // Open file
      const std::string file_name = "boundaryData.slot-"
                                    + Utilities::int_to_string (i_slot,4)
                                    + ".raw";
      std::ifstream file_in (file_name.c_str());
      if (!file_in)
        {
          std::cerr << "There is no file '" << file_name << "'.\n"
                    << "Exit!" << std::endl;
          break;
        }
      std::cerr << "Processing '" << file_name << "'." << std::endl;

      // Read first line for Mach number
      {
        double slot_Mach;
        file_in >> slot_Mach;
        AssertThrow (slot_Mach == Mach || Mach < 0.0,
                     ExcMessage("Referencing Mach number mismatch!"));
        Mach = slot_Mach;
      }

      // Go through all lines
      while (true)
        {
          double x_index;
          DataStruct data_struct;

          for (unsigned int id=0; id<dim; ++id)
            {
              file_in >> data_struct.position[id];
            }

          // Check file status after first reading try to avoid
          // trouble in ending white line.
          if (file_in.eof())
            {
              break;
            }

          for (unsigned int ic=0; ic<n_components; ++ic)
            {
              file_in >> data_struct.solution[ic];
            }

          data_struct.velocity_magnitude = 0.0;
          for (unsigned int ic=first_velocity_component, id=0; id<dim; ++ic, ++id)
            {
              data_struct.velocity_magnitude += data_struct.solution[ic] * data_struct.solution[ic];
            }
          data_struct.velocity_magnitude = std::sqrt(data_struct.velocity_magnitude);

          data_struct.sound_speed
            = std::sqrt (gas_gamm *
                         data_struct.solution[pressure_component] /
                         data_struct.solution[density_component]);
          data_struct.local_Mach = data_struct.velocity_magnitude /
                                   data_struct.sound_speed;

          data_struct.pressure_coefficient =
            2.0 * (data_struct.solution[pressure_component] - 1.0/gas_gamm) /
            (Mach * Mach);

          // Order data according to x-coordinate by default.
          x_index = data_struct.position[0];
          if (data_struct.position[dim-1] < 0.0)
            {
              // Order date circle form TE for foil case.
              // Z-up configuration in 3D is assumed.
              x_index = std::copysign(x_index, data_struct.position[dim-1]);
            }
          else if (data_struct.position[dim-1] > 0.9)
            {
              // Order date counter clock wise for GAMM Channel case
              x_index = -3.0 - x_index;
            }

          // Insert new entry
          data_map[x_index] = data_struct;
        }
      file_in.close();
    }

  // Write out collected data
  {
    std::ofstream file_out ("collectedBoundaryData.txt");

    file_out << "# Mach_{infty} = " << Mach << std::endl;

    {
      unsigned int column = 0;
      file_out << "#" << ++column << " x-index    ";
      file_out << "#" << ++column << " x           ";
      file_out << "#" << ++column << " y          ";
      if (dim == 3)
        {
          file_out << "#" << ++column << " z          ";
        }
      file_out << "#" << ++column << " u          ";
      file_out << "#" << ++column << " v          ";
      if (dim == 3)
        {
          file_out << "#" << ++column << " w          ";
        }
      file_out << "#" << ++column << " rho        ";
      file_out << "#" << ++column << " p          ";
      file_out << "#" << ++column << " |velocity| ";
      file_out << "#" << ++column << " v_sound    ";
      file_out << "#" << ++column << " Mach       ";
      file_out << "#" << ++column << " Cp         ";
      file_out << std::endl;
    }

    file_out.setf (std::ios::scientific);
    file_out.precision (6);

    for (typename DataMap::iterator it=data_map.begin(); it!=data_map.end(); ++it)
      {
        const unsigned int column_with = 13;
        file_out << std::setw (column_with) << it->first << ' ';
        const DataStruct &data_struct = it->second;
        for (unsigned int id=0; id<dim; ++id)
          {
            file_out << std::setw (column_with) << data_struct.position[id] << ' ';
          }
        for (unsigned int ic=0; ic<n_components; ++ic)
          {
            file_out << std::setw (column_with) << data_struct.solution[ic] << ' ';
          }
        file_out
            << std::setw (column_with) << data_struct.velocity_magnitude <<  ' '
            << std::setw (column_with) << data_struct.sound_speed <<  ' '
            << std::setw (column_with) << data_struct.local_Mach <<  ' '
            << std::setw (column_with) << data_struct.pressure_coefficient <<  ' '
            << std::endl;
      }

    file_out.close();
  }

  return;
}
