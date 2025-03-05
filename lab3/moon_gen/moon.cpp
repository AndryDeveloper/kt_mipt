#include <gmsh.h>

int main(int argc, char **argv) {
    gmsh::initialize();
    gmsh::model::add("moon");
    // gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.05);
    gmsh::option::setNumber("Mesh.MshFileVersion", 2);

    double outer_r = 1.0;
    double inner_r = 0.8;
    double offset = 0.3;

    int outer_circle = gmsh::model::occ::addDisk(0, 0, 0, outer_r, outer_r);
    int inner_circle = gmsh::model::occ::addDisk(offset, 0, 0, inner_r, inner_r);

    gmsh::vectorpair outDimTags;
    std::vector<std::vector<std::pair<int, int>>> ov;
    gmsh::model::occ::cut({{2, outer_circle}}, {{2, inner_circle}}, outDimTags, ov);

    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(2);
    gmsh::write("../moon.msh");

    gmsh::finalize();
    return 0;
}