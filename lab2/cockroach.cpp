#include <iostream>
#include <cmath>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

using namespace std;


struct Box {
    double x_min, x_max, z_min, z_max, y_min, y_max;
    Box(double x_min, double x_max, double z_min, double z_max, double y_min, double y_max)
        : x_min(x_min), x_max(x_max), z_min(z_min), z_max(z_max), y_min(y_min), y_max(y_max) {}
};

struct Pair {
    double x; 
    double z;
};

Box get_box(int idx) {
    switch (idx) {
        case 0: return Box(12.8, 32.7, 31.5, 42, 0.323, 1.63);
        case 1: return Box(36.5, 40.2, 30.3, 44.5, 0.323, 1.63);
        case 2: return Box(46.4, 51.9, 28.9, 34.3, 0.323, 1.63);
        case 3: return Box(49.3, 68.2, 26, 42, 0.323, 1.63);
        case 4: return Box(12.8, 32.7, 6., 16.1, 0.323, 1.63);
        case 5: return Box(36.5, 40.2, 2.0, 17.2, 0.323, 1.63);
        case 6: return Box(46.4, 51.9, 13.2, 18.8, 0.323, 1.63);
        case 7: return Box(49.3, 68.2, 6, 21.9, 0.323, 1.63);
        default: return Box(0, 0, 0, 0, 0, 0); // Значение по умолчанию
    }
}

bool is_inside_box(int box_id, double pointX, double pointY, double pointZ) {
    Box box = get_box(box_id);
    return (pointX >= box.x_min && pointX <= box.x_max &&
            pointZ >= box.z_min && pointZ <= box.z_max &&
            pointY >= box.y_min && pointY <= box.y_max);
}

Pair get_center(Box box) {
    const double x_0 = 31.5;
    const double z_0 = 24.5;

    Pair corner1 = {box.x_min, box.z_min}; // (x_min, z_min)
    Pair corner2 = {box.x_max, box.z_max}; // (x_max, z_max)

    double distance1 = std::sqrt(std::pow(corner1.x - x_0, 2) + std::pow(corner1.z - z_0, 2));
    double distance2 = std::sqrt(std::pow(corner2.x - x_0, 2) + std::pow(corner2.z - z_0, 2));

    return (distance1 < distance2) ? corner1 : corner2;
}


// Класс расчётной точки
class CalcNode
{
// Класс сетки будет friend-ом точки
friend class CalcMesh;

protected:
    // Координаты
    double x;
    double y;
    double z;
    // Некая величина, в попугаях
    double smth;
    // Скорость
    double vx;
    double vy;
    double vz;

public:
    // Конструктор по умолчанию
    CalcNode() : x(0.0), y(0.0), z(0.0), smth(0.0), vx(0.0), vy(0.0), vz(0.0)
    {
    }

    // Конструктор с указанием всех параметров
    CalcNode(double x, double y, double z, double smth, double vx, double vy, double vz) 
            : x(x), y(y), z(z), smth(smth), vx(vx), vy(vy), vz(vz)
    {
    }

    // Метод отвечает за перемещение точки
    // Движемся время tau из текущего положения с текущей скоростью
    void move(double tau) {
        x += vx * tau;
        y += vy * tau;
        z += vz * tau;
    }
};

// Класс элемента сетки
class Element
{
// Класс сетки будет friend-ом и элемента тоже
// (и вообще будет нагло считать его просто структурой)
friend class CalcMesh;

protected:
    // Индексы узлов, образующих этот элемент сетки
    unsigned long nodesIds[4];
};

// Класс расчётной сетки
class CalcMesh
{
protected:
    // 3D-сетка из расчётных точек
    vector<CalcNode> nodes;
    vector<Element> elements;
    vector<int> laps[8];

public:
    // Конструктор сетки из заданного stl-файла
    CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints) {

        // Пройдём по узлам в модели gmsh
        nodes.resize(nodesCoords.size() / 3);
        for(unsigned int i = 0; i < nodesCoords.size() / 3; i++) {
            // Координаты заберём из gmsh
            double pointX = nodesCoords[i*3];
            double pointY = nodesCoords[i*3 + 1];
            double pointZ = nodesCoords[i*3 + 2];
            // Модельная скалярная величина распределена как-то вот так
            double smth = pow(pointX, 2)*sin(pointX) + pow(pointY, 2)*sin(pointY) + pow(pointZ, 2)*sin(pointZ);

            bool p = false;
            for (int j = 0; j < 8; j++){
                bool q = is_inside_box(j, pointX, pointY, pointZ);
                for (int k = 0; k < 8; k++){
                    if (k != j && is_inside_box(k, pointX, pointY, pointZ)){
                        q = false;
                    }
                }

                if (q){
                    laps[j].push_back(i);
                    nodes[i] = CalcNode(pointX, pointY, pointZ, smth, 50.0, 0.0, 0.0);
                    p = true;
                }
            }
            if (!p){
                nodes[i] = CalcNode(pointX, pointY, pointZ, smth, 50.0, 0.0, 0.0);
            }
        }
        std::cout << laps[0].size() << ' ' << 
                    laps[1].size() << ' ' <<
                    laps[2].size() << ' ' <<
                    laps[3].size() << ' ' <<
                    laps[4].size() << ' ' <<
                    laps[5].size() << ' ' <<
                    laps[6].size() << ' ' <<
                    laps[7].size() << ' ' << std::endl;

        // Пройдём по элементам в модели gmsh
        elements.resize(tetrsPoints.size() / 4);
        for(unsigned int i = 0; i < tetrsPoints.size() / 4; i++) {
            elements[i].nodesIds[0] = tetrsPoints[i*4] - 1;
            elements[i].nodesIds[1] = tetrsPoints[i*4 + 1] - 1;
            elements[i].nodesIds[2] = tetrsPoints[i*4 + 2] - 1;
            elements[i].nodesIds[3] = tetrsPoints[i*4 + 3] - 1;
        }
    }

    // Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    void doTimeStep(double tau, double step) {
        if (step < 60){
            for (int j = 0; j < 8; j++){
                Box box = get_box(j);
                Pair cnt = get_center(box);
    
                for (int idx : laps[j]){
                    nodes[idx].vx = -(nodes[idx].z - cnt.z) * 10 * cos(step / 3.) + 50;
                    nodes[idx].vz = (nodes[idx].x - cnt.x - tau*50*step) * 10 * cos(step / 3.);
                }
            }
            for(unsigned int i = 0; i < nodes.size(); i++) {
                nodes[i].smth = pow(nodes[i].x, 2)*sin(nodes[i].x) + pow(nodes[i].y, 2)*sin(nodes[i].y) + pow(nodes[i].z, 2)*sin(nodes[i].z);
            }
        }
        else {
            double x_0 = 31.5 + tau*50*step;
            double z_0 = 24.5;
            
            for(unsigned int i = 0; i < nodes.size(); i++) {

                nodes[i].vx += (nodes[i].x - x_0)*abs(nodes[i].x - x_0);
                nodes[i].vy += (nodes[i].y - 0.5)*abs(nodes[i].x - 0.5);
                nodes[i].vz += (nodes[i].z - z_0)*abs(nodes[i].x - z_0);
            }
            for(unsigned int i = 0; i < nodes.size(); i++) {
                nodes[i].smth *= ((nodes[i].x - x_0)*(nodes[i].x - x_0) + (nodes[i].y - 0.5)*(nodes[i].y - 0.5) + (nodes[i].z - z_0)*(nodes[i].z - z_0));
            }
        }
        // По сути метод просто двигает все точки
        for(unsigned int i = 0; i < nodes.size(); i++) {
            nodes[i].move(tau);
        }
    }

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number) {
        // Сетка в терминах VTK
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        // Точки сетки в терминах VTK
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        // Скалярное поле на точках сетки
        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("smth");

        // Векторное поле на точках сетки
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        // Обходим все точки нашей расчётной сетки
        for(unsigned int i = 0; i < nodes.size(); i++) {
            // Вставляем новую точку в сетку VTK-снапшота
            dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);

            // Добавляем значение векторного поля в этой точке
            double _vel[3] = {nodes[i].vx, nodes[i].vy, nodes[i].vz};
            vel->InsertNextTuple(_vel);

            // И значение скалярного поля тоже
            smth->InsertNextValue(nodes[i].smth);
        }

        // Грузим точки в сетку
        unstructuredGrid->SetPoints(dumpPoints);

        // Присоединяем векторное и скалярное поля к точкам
        unstructuredGrid->GetPointData()->AddArray(vel);
        unstructuredGrid->GetPointData()->AddArray(smth);

        // А теперь пишем, как наши точки объединены в тетраэдры
        for(unsigned int i = 0; i < elements.size(); i++) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId( 0, elements[i].nodesIds[0] );
            tetra->GetPointIds()->SetId( 1, elements[i].nodesIds[1] );
            tetra->GetPointIds()->SetId( 2, elements[i].nodesIds[2] );
            tetra->GetPointIds()->SetId( 3, elements[i].nodesIds[3] );
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        // Создаём снапшот в файле с заданным именем
        string fileName = "output/cockroach-step-" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
};

int main()
{
    // Шаг точек по пространству
    double h = 4.0;
    // Шаг по времени
    double tau = 0.003;

    const unsigned int GMSH_TETR_CODE = 4;

    // Теперь придётся немного упороться:
    // (а) построением сетки средствами gmsh,
    // (б) извлечением данных этой сетки в свой код.
    gmsh::initialize();
    gmsh::model::add("t13");

    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 1.);

    // Считаем STL
    try {
        gmsh::merge("cockroach.STL"); 
        // путь к файлу отсчитывается от точки запуска 
        // если вы собирали все в директории build как цивилизованные люди, 
        // то перейдите на уровень выше
        // и запускайте бинарник как ./build/tetr3d
    } catch(...) {
        gmsh::logger::write("Could not load STL mesh: bye!");
        gmsh::finalize();
        return -1;
    }

    // Восстановим геометрию
    double angle = 40;
    bool forceParametrizablePatches = false;
    bool includeBoundary = true;
    double curveAngle = 180;
    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary, forceParametrizablePatches, curveAngle * M_PI / 180.);
    gmsh::model::mesh::createGeometry();

    // Зададим объём по считанной поверхности
    std::vector<std::pair<int, int> > s;
    gmsh::model::getEntities(s, 2);
    std::vector<int> sl;
    for(auto surf : s) sl.push_back(surf.second);
    int l = gmsh::model::geo::addSurfaceLoop(sl);
    gmsh::model::geo::addVolume({l});

    gmsh::model::geo::synchronize();

    // Зададим мелкость желаемой сетки
    int f = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(f, "F", "4");
    gmsh::model::mesh::field::setAsBackgroundMesh(f);

    // Построим сетку
    gmsh::model::mesh::generate(3);

    // Теперь извлечём из gmsh данные об узлах сетки
    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    // И данные об элементах сетки тоже извлечём, нам среди них нужны только тетраэдры, которыми залит объём
    std::vector<std::size_t>* tetrsNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
    for(unsigned int i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] != GMSH_TETR_CODE)
            continue;
        tetrsNodesTags = &elementNodeTags[i];
    }

    if(tetrsNodesTags == nullptr) {
        cout << "Can not find tetra data. Exiting." << endl;
        gmsh::finalize();
        return -2;
    }

    cout << "The model has " <<  nodeTags.size() << " nodes and " << tetrsNodesTags->size() / 4 << " tetrs." << endl;

    // На всякий случай проверим, что номера узлов идут подряд и без пробелов
    for(int i = 0; i < nodeTags.size(); ++i) {
        // Индексация в gmsh начинается с 1, а не с нуля. Ну штош, значит так.
        assert(i == nodeTags[i] - 1);
    }
    // И ещё проверим, что в тетраэдрах что-то похожее на правду лежит.
    assert(tetrsNodesTags->size() % 4 == 0);

    // TODO: неплохо бы полноценно данные сетки проверять, да

    CalcMesh mesh(nodesCoord, *tetrsNodesTags);

    gmsh::finalize();

    mesh.snapshot(0);

    for(unsigned int step = 1; step < 70; step++) {
        mesh.doTimeStep(tau, step);
        mesh.snapshot(step);
    }

    return 0;
}
