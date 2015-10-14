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

#include <deal.II/base/std_cxx11/bind.h>
#include <deal.II/grid/tria.h>
#include <iostream>


struct ClassA
{
  boost::signals2::signal<void ()>    SigA;
  boost::signals2::signal<void (int)> SigB;
};

struct ClassB
{
  virtual void PrintFoo()
  {
    std::cout << "ClassB::Foo" << std::endl;
  }
  virtual void PrintInt(int i)
  {
    std::cout << "ClassB::Bar: " << i << std::endl;
  }
};

struct ClassC : public ClassB
{
  virtual void PrintFoo()
  {
    std::cout << "ClassC::Foo" << std::endl;
  }
  virtual void PrintInt(int i)
  {
    std::cout << "ClassC::Bar: " << i << std::endl;
  }
};

int main()
{
  ClassA a;
  ClassB b, b2;

  ClassB *c1 = new ClassC;
  ClassC *c2 = new ClassC;

  a.SigA.connect(dealii::std_cxx11::bind(&ClassB::PrintFoo, &b));
  a.SigB.connect(std::bind(&ClassB::PrintInt, &b,  std::placeholders::_1));
  a.SigB.connect(std::bind(&ClassB::PrintInt, &b2, std::placeholders::_1));
  a.SigB.connect(std::bind(&ClassB::PrintInt, c1,  std::placeholders::_1));
  a.SigB.connect(std::bind(&ClassC::PrintInt, c2,  std::placeholders::_1));

  a.SigA();
  a.SigB(4);

  delete c1;
  delete c2;

  return (0);
}

