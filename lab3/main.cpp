#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

// Source term (right-hand side)
class Source : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 0;
  }
};

// Normal derivative (Neumann boundary condition)
class dUdN : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = exp(x[0]);
  }
};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    double dist = std::sqrt((x[0]+0.9) * (x[0]+0.9) + (x[1]) * (x[1]));
    return (dist < 0.2 + DOLFIN_EPS) && on_boundary;
  }
};


int main()
{
  // Create mesh and function space
  auto mesh = std::make_shared<Mesh>("moon.xml"); // Загружаем сетку
  // auto mesh = std::make_shared<Mesh>(
  //   UnitSquareMesh::create({{32, 32}}, CellType::Type::triangle));
  auto V = std::make_shared<Poisson::FunctionSpace>(mesh);

  // Define boundary condition
  auto u0 = std::make_shared<Constant>(1.0);
  auto boundary = std::make_shared<DirichletBoundary>();
  DirichletBC bc(V, u0, boundary);

  // Define variational forms
  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);
  auto f = std::make_shared<Source>();
  auto g = std::make_shared<dUdN>();
  L.f = f;
  L.g = g;

  // Compute solution
  Function u(V);
  solve(a == L, u, bc);

  // Save solution in VTK format
  File file("poisson.pvd");
  file << u;

  return 0;
}