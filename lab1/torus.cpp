#include <gmsh.h>


int main(int argc, char **argv) {
    gmsh::initialize();
    gmsh::model::add("torus");
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.3);

    // Радиусы тора
    double R = 10.0;
    double r = 4.;
    double delta = 1.;

    int outer_torus = gmsh::model::occ::addTorus(0, 0, 0, R, r);
    int inner_torus = gmsh::model::occ::addTorus(0, 0, 0, R, r - delta);

    gmsh::vectorpair outDimTags;
    std::vector<std::vector<std::pair<int, int>>> ov;
    gmsh::model::occ::cut({{3, outer_torus}}, {{3, inner_torus}}, outDimTags, ov, 3, true, false);

    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(3);

    gmsh::write("output/torus.msh");
    // gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}