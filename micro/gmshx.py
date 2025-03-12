import gmsh
import argparse
import os
import math

def main():
    parser = argparse.ArgumentParser(description="Генерация сетки для аэродинамических расчётов")
    parser.add_argument("filename", help="Имя STEP-файла (например, model.step)")
    args = parser.parse_args()

    gmsh.initialize()
    gmsh.model.add("AeroMesh")

    # Импорт и поворот модели ====================================================
    step_file = args.filename
    dimTags = gmsh.model.occ.importShapes(step_file, highestDimOnly=True)
    
    # Поворачиваем модель на 90 градусов вокруг оси X для выравнивания по OZ
    gmsh.model.occ.rotate(dimTags, 0, 0, 0, 1, 0, 0, math.pi/2)
    gmsh.model.occ.rotate(dimTags, 0, 0, 0, 0, 0, 1, math.pi/2)
    gmsh.model.occ.synchronize()

    # Получаем поверхности объекта после поворота
    object_surfaces = gmsh.model.getBoundary(dimTags, oriented=False)

    # Определяем новые габариты после поворота
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(*dimTags[0])
    L = max(xmax - xmin, ymax - ymin, zmax - zmin)

    # Параметры расчётной области
    padding = {
        'upstream': 3 * L,
        'downstream': 3 * L,
        'sides': 2 * L,
        'top': 2 * L,
        'bottom': 0.05 * L
    }

    # Создаём расчётную область
    box = gmsh.model.occ.addBox(
        xmin - padding['upstream'],
        ymin - padding['sides'],
        zmin - padding['bottom'],
        (xmax - xmin) + padding['upstream'] + padding['downstream'],
        (ymax - ymin) + 2 * padding['sides'],
        (zmax - zmin) + padding['top'] + padding['bottom']
    )

    # Вырезаем объект из расчётной области
    fluid_domain, _ = gmsh.model.occ.cut([(3, box)], dimTags, removeObject=True, removeTool=True)
    gmsh.model.occ.synchronize()

    # Настройка размера сетки ====================================================
    # base_size = L/20 # Базовый размер сетки
    
    # # Градиентное увеличение размера элементов от объекта
    # gmsh.option.setNumber("Mesh.CharacteristicLengthMin", base_size)
    # gmsh.option.setNumber("Mesh.CharacteristicLengthMax", base_size * 20)
    
    # # Поле для плавного изменения размера элементов
    # gmsh.model.mesh.field.add("Distance", 1)
    # gmsh.model.mesh.field.setNumbers(1, "FacesList", [s[1] for s in object_surfaces])
    # gmsh.model.mesh.field.add("Threshold", 2)
    # gmsh.model.mesh.field.setNumber(2, "InField", 1)
    # gmsh.model.mesh.field.setNumber(2, "SizeMin", base_size)
    # gmsh.model.mesh.field.setNumber(2, "SizeMax", base_size * 20)
    # gmsh.model.mesh.field.setNumber(2, "DistMin", 0.1*L)
    # gmsh.model.mesh.field.setNumber(2, "DistMax", 2*L)
    # gmsh.model.mesh.field.setAsBackgroundMesh(2)

    # Создаём физические группы
    
    # 1. Объём жидкости (3D)
    gmsh.model.addPhysicalGroup(3, [fluid_domain[0][1]], 1)
    gmsh.model.setPhysicalName(3, 1, "Fluid")

    # 2. Граничные поверхности (2D)
    all_surfaces = gmsh.model.getEntities(dim=2)
    
    # Идентификация границ
    external_tags = []
    ground_tags = []

    for surface in all_surfaces:
        com = gmsh.model.occ.getCenterOfMass(*surface)
        
        # Земля (Z = zmin - padding['bottom'])
        if abs(com[2] - (zmin - padding['bottom'])) < 1e-3:
            ground_tags.append(surface[1])
        
        # Внешние границы
        elif (
            abs(com[0] - (xmin - padding['upstream'])) < 1e-3 or      # Вход потока
            abs(com[0] - (xmax + padding['downstream'])) < 1e-3 or    # Выход
            abs(com[1] - (ymin - padding['sides'])) < 1e-3 or         # Левая стенка
            abs(com[1] - (ymax + padding['sides'])) < 1e-3 or         # Правая стенка
            abs(com[2] - (zmax + padding['top'])) < 1e-3              # Верхняя граница
        ):
            external_tags.append(surface[1])

    # Назначаем физические группы
    gmsh.model.addPhysicalGroup(2, external_tags, 2)
    gmsh.model.setPhysicalName(2, 2, "External")
    
    gmsh.model.addPhysicalGroup(2, ground_tags, 4)
    gmsh.model.setPhysicalName(2, 4, "Ground")
    
    obj_tags = [s[1] for s in object_surfaces]
    gmsh.model.addPhysicalGroup(2, obj_tags, 3)
    gmsh.model.setPhysicalName(2, 3, "Object")

    # Настройки генерации сетки
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    gmsh.option.setNumber("Mesh.Optimize", 2)
    gmsh.option.setNumber("Mesh.MeshSizeFactor", 2.)
    
    # Генерация и экспорт сетки
    gmsh.model.mesh.generate(3)
    output_dir = "meshs"
    os.makedirs(output_dir, exist_ok=True)
    msh_file = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(step_file))[0]}.msh")
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(msh_file)
    gmsh.finalize()

if __name__ == "__main__":
    main()