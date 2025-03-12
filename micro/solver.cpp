#include <dolfin.h>
#include "NavierStokes.h"

using namespace dolfin;

int main()
{
    // Загрузка сетки и граничных меток
    auto mesh = std::make_shared<Mesh>("xmls/cockroach.xml");
    auto boundaries = std::make_shared<MeshFunction<std::size_t>>(mesh, "xmls/cockroach_facet_region.xml");

    // Создание функционального пространства
    auto W = std::make_shared<NavierStokes::FunctionSpace>(mesh);

    // Граничные условия
    auto u_in = std::make_shared<Constant>(1.0, 0.0, 0.0);
    auto zero = std::make_shared<Constant>(0.0, 0.0, 0.0);
    DirichletBC bc_wall2(W->sub(0), u_in, boundaries, 2);
    DirichletBC bc_wall3(W->sub(0), zero, boundaries, 3);
    DirichletBC bc_wall4(W->sub(0), zero, boundaries, 4);
    std::vector<const DirichletBC*> bcs = { &bc_wall2, &bc_wall3, &bc_wall4 };

    // Используем shared_ptr для u0
    auto u0 = std::make_shared<Function>(W->sub(0)->collapse());
    u0->vector()->zero();

    // Внешняя сила
    NavierStokes::LinearForm L(W);
    L.f = std::make_shared<Constant>(0.0, 0.0, 0.0);

    const double tol = 1e-6;
    const int maxiter = 20;
    double err = 1.0;
    int iter = 0;
    Function w(W);

    // Parameters params("solver_parameters");
    
    // params.add("linear_solver", "gmres");
    // params.add("preconditioner", "ilu");

    // Parameters krylov_params = params("krylov_solver");
    // krylov_params.add("relative_tolerance", 1e-6);
    // krylov_params.add("absolute_tolerance", 1e-10);
    // krylov_params.add("maximum_iterations", 1000);

    while (err > tol && iter < maxiter)
    {
        NavierStokes::BilinearForm a(W, W);
        a.u0 = u0;
    
        solve(a == L, w, bcs);
        // solve(a == L, w, bcs, params);
    
        Function u = w[0];
        auto diff_vector = std::make_shared<Vector>(*u.vector());
        *diff_vector -= *u0->vector();
        double diff = diff_vector->norm("l2");
    
        std::cout << "Iteration " << iter << ": error = " << diff << std::endl;
        err = diff;
        *u0 = u;
        iter++;
    }

    File ufile("velocity.pvd");
    ufile << w[0];
    File pfile("pressure.pvd");
    pfile << w[1];

    File mf("mf.pvd");
    mf << *boundaries;

    return 0;
}