from __future__ import annotations

from itertools import cycle

from openalea.plantgl.all import (
    AxisRotated,
    Color3,
    Material,
    Scene,
    Shape,
    Translated,
    Viewer,
)


def build_scene(
    mtg,
    position=(0, 0, 0),
    orientation=0,
    leaf_material=None,
    stem_material=None,
    soil_material=None,
    colors=None,
):
    """
    Returns a plantgl scene from an mtg.
    """
    if not isinstance(mtg, list):
        mtg = [mtg]
    if not isinstance(position, list):
        position = [position]
    if not isinstance(orientation, list):
        orientation = [orientation]
    if colors is None:
        if leaf_material is None:
            leaf_material = Material(Color3(0, 180, 0))
        if stem_material is None:
            stem_material = Material(Color3(0, 130, 0))
        if soil_material is None:
            soil_material = Material(Color3(170, 85, 0))
            # colors = g.property('color')

    scene = Scene()

    def geom2shape(vid, mesh, scene, colors, position, orientation, shape_id=None):
        shape = None
        if shape_id is None:
            shape_id = vid
        if isinstance(mesh, list):
            for m in mesh:
                geom2shape(vid, m, scene, colors, position, orientation)
            return
        if mesh is None:
            return
        if isinstance(mesh, Shape):
            shape = mesh
            mesh = mesh.geometry
        label = labels.get(vid)
        is_green = greeness.get(vid)
        mesh = Translated(position, AxisRotated((0, 0, 1), orientation, mesh))
        if colors:
            shape = Shape(mesh, Material(Color3(*colors.get(vid, [0, 0, 0]))))
        elif not greeness:
            if not shape:
                shape = Shape(mesh)
        elif label.startswith("Stem") and is_green:
            shape = Shape(mesh, stem_material)
        elif label.startswith("Leaf") and is_green:
            shape = Shape(mesh, leaf_material)
        elif not is_green:
            shape = Shape(mesh, soil_material)
        shape.id = shape_id

        scene.add(shape)

    nump = []
    count = 0
    for i, (g, p, o) in enumerate(zip(cycle(mtg), position, cycle(orientation))):
        geometries = g.property("geometry")
        greeness = g.property("is_green")
        labels = g.property("label")

        for vid, mesh in geometries.items():
            geom2shape(vid, mesh, scene, colors, p, o, vid + count)
            nump.append(i)
        count += len(g)

    return scene, nump


def build_IC_scene(
    mtg1,
    mtg2,
    orientation=0,
    leaf_material=None,
    stem_material=None,
    soil_material=None,
    colors=None,
):
    """
    Returns a plantgl scene from 2 mtg, for intercropping.
    """
    if not isinstance(mtg1, list):
        mtg1 = [mtg1]
    if not isinstance(mtg2, list):
        mtg2 = [mtg2]
    if not isinstance(orientation, list):
        orientation = [orientation]
    if colors is None:
        if leaf_material is None:
            leaf_material = Material(Color3(0, 180, 0))
        if stem_material is None:
            stem_material = Material(Color3(0, 130, 0))
        if soil_material is None:
            soil_material = Material(Color3(170, 85, 0))
            # colors = g.property('color')

    scene = Scene()

    def geom2shape(vid, mesh, scene, colors, position, orientation, shape_id=None):
        shape = None
        if shape_id is None:
            shape_id = vid
        if isinstance(mesh, list):
            for m in mesh:
                geom2shape(vid, m, scene, colors, position, orientation)
            return
        if mesh is None:
            return
        if isinstance(mesh, Shape):
            shape = mesh
            mesh = mesh.geometry
        label = labels.get(vid)
        is_green = greeness.get(vid)
        mesh = Translated(position, mesh)  # AxisRotated((0,0,1),orientation, mesh))
        if colors:
            shape = Shape(mesh, Material(Color3(*colors.get(vid, [0, 0, 0]))))
        elif not greeness:
            if not shape:
                shape = Shape(mesh)
        elif label.startswith("Stem") and is_green:
            shape = Shape(mesh, stem_material)
        elif label.startswith("Leaf") and is_green:
            shape = Shape(mesh, leaf_material)
        elif not is_green:
            shape = Shape(mesh, soil_material)
        shape.id = shape_id

        scene.add(shape)

    nump = []
    count = 0

    for mtg in [mtg1, mtg2]:
        if mtg == mtg1:
            position = (0, 0, 0)

            for i, (g, p, o) in enumerate(
                zip(cycle(mtg), position, cycle(orientation))
            ):
                geometries = g.property("geometry")
                greeness = g.property("is_green")
                labels = g.property("label")

                for vid, mesh in geometries.items():
                    geom2shape(vid, mesh, scene, colors, p, o, vid + count)
                    nump.append(i)
                count += len(g)

        else:
            position = (1, 0, 0)
            ts = Translated(position, mtg2)
            sh = Shape(ts, Material(Color3(0, 100, 0)))
            scene.add(sh)

    return scene, nump


def display_scene(scene):
    """
    Display a plantgl scene
    """

    Viewer.display(scene)


def display_mtg(
    mtg, leaf_material=None, stem_material=None, soil_material=None, colors=None
):
    """
    Display a scene from a mtg and return the scene displayed.
    """

    scene = build_scene(
        mtg,
        leaf_material=leaf_material,
        stem_material=stem_material,
        soil_material=soil_material,
        colors=colors,
    )

    display_scene(scene)

    return scene
