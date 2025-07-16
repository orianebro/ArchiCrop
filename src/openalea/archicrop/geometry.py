"""Geometric primitives for cereals"""

from __future__ import annotations

import openalea.plantgl.all as pgl
from openalea.mtg.turtle import TurtleFrame

from .cereal_leaf import compute_leaf_area, leaf_mesh
from .growth import compute_organ_growth
from .internode import stem_mesh


def leaf_area_plant(g):
    S = 0
    for k,leaf in g.properties()["shape"].items():
        S += compute_leaf_area(leaf, g.properties()["visible_length"][k], g.properties()["mature_length"][k], g.properties()["shape_max_width"][k])
    return S
    # return sum(g.properties()["visible_leaf_area"].values())


def mesh_area(mesh):
    if mesh:
        sc = pgl.SurfComputer(pgl.Discretizer())
        mesh.apply(sc)
        return sc.surface
    return None


def _is_iterable(x):
    try:
        iter(x)
    except TypeError:
        return False
    return True


def as_tuples(pgl_3List, offset=0):
    """return pgl list of 3 numbers kind (indes3, vector3) as a list of python
    tuples
    """
    if not _is_iterable(offset):
        offset = [offset] * 3
    return [(i[0] + offset[0], i[1] + offset[1], i[2] + offset[2]) for i in pgl_3List]


def addSets(pglset1, pglset2, translate=(0, 0, 0)):
    """
    Create a new TriangleSet by addition of two existing ones
    if translate is not None, pglset2 is translated with vector translate
    """
    points = as_tuples(pglset1.pointList) + as_tuples(
        pglset2.pointList, offset=translate
    )
    index = as_tuples(pglset1.indexList) + as_tuples(
        pglset2.indexList, offset=len(pglset1.pointList)
    )
    return pgl.TriangleSet(points, index)


def compute_organ_geometry(element_node, classic=False):
    """compute geometry of Adel base elements (LeafElement and StemElement)
    element_node should be a mtg node proxy"""
    n = element_node
    geom = None

    if n.label.startswith("Leaf"):  # leaf element
        if n.visible_length > 0.0001:  # filter less than 0.001 mm leaves
            if n.shape is not None and n.srb is not None:
                geom = leaf_mesh(
                    n.shape,
                    n.mature_length,
                    n.shape_max_width,
                    n.visible_length,
                    n.srb,
                    n.srt,
                    # flipx allows x-> -x to place the shape along
                    #  with the tiller positioned with
                    # turtle.down()
                    flipx=True,
                    stem_diameter=n.stem_diameter,
                )
            if n.lrolled > 0:
                rolled = stem_mesh(
                    n.lrolled, n.lrolled, n.d_rolled, classic
                )
                if geom is None:
                    geom = rolled
                else:
                    geom = addSets(rolled, geom, translate=(0, 0, n.lrolled))
    elif n.label.startswith("Stem"):  # stem element
        geom = stem_mesh(n.length, n.visible_length, n.stem_diameter, n.stem_diameter, classic)

    return geom


class CerealsTurtle(pgl.PglTurtle):
    def __init__(self):
        super().__init__()
        self.context = {}

    def transform(self, mesh, face_up=False):
        x = self.getUp()
        z = pgl.Vector3(0, 0, 1) if face_up else self.getHeading()
        bo = pgl.BaseOrientation(x, z ^ x)
        matrix = pgl.Transform4(bo.getMatrix())
        matrix.translate(self.getPosition())
        # print 'Position ', turtle.getPosition()
        return mesh.transform(matrix)

    def getFrame(self):
        pos = self.getPosition()
        Head = self.getHeading()
        Up = self.getUp()
        return {"Position": pos, "Head": Head, "Up": Up}

    def setFrame(self, frame):
        self.move(frame["Position"])
        self.setHead(frame["Head"], frame["Up"])


class CerealsVisitor:
    """Performs geometric interpretation of mtg nodes"""

    def __init__(self, classic):
        self.classic = classic

    def __call__(self, g, v, turtle):
        # 1. retrieve the node
        n = g.node(v)

        # Go to plant position if first plant element
        if n.parent() is None:
            turtle.move(0, 0, 0)
            # initial position to be compatible with canMTG positioning
            turtle.setHead(0, 0, 1, -1, 0, 0)

        # incline turtle at the base of stems,
        if n.label.startswith("Stem"):
            azim = float(n.azimuth) if n.azimuth else 0.0
            if azim:
                # print 'node', n._vid, 'azim ', azim
                turtle.rollL(azim)

        if n.label.startswith("Leaf") or n.label.startswith("Stem"):  # noqa: SIM102
            # update geometry of elements
            if n.length > 0:
                mesh = compute_organ_geometry(n, self.classic)
                if mesh:  # To DO : reset to None if calculated so ?
                    n.geometry = turtle.transform(mesh)
                    n.anchor_point = turtle.getPosition()
                    n.heading = turtle.getHeading()

        # 3. Update the turtle and context
        turtle.setId(v)
        if n.label.startswith("Stem"):
            if n.visible_length > 0:  
                turtle.f(n.visible_length)
            turtle.context.update({"top": turtle.getFrame()})
        if n.label.startswith("Leaf"):  # noqa: SIM102
            if n.lrolled > 0:
                turtle.f(n.lrolled)
                turtle.context.update({"top": turtle.getFrame()})


def mtg_interpreter(g, classic=False):
    """
    Compute/update the geometry on each node of the MTG using Turtle geometry.
    """
    # BUG : sub_mtg mange le vertex plant => on perd la plante !
    # plants = g.component_roots_at_scale(g.root, scale=1)
    # nplants = g.nb_vertices(scale=1)
    # gt = MTG()

    # for plant in plants:
    #   gplant = g.sub_mtg(plant)
    turtle = CerealsTurtle()
    visitor = CerealsVisitor(classic)

    scene = TurtleFrame(g, visitor=visitor, turtle=turtle, gc=False, all_roots=True)  # noqa: F841

    return g

### for growing organs

def compute_growing_organ_geometry(element_node, classic=False):
    """compute geometry of Adel base elements (LeafElement and StemElement)
    element_node should be a mtg node proxy"""

    n = element_node
    geom = None

    if n.label.startswith("Leaf"):  # leaf element
        if n.visible_length > 0.0001:  # filter less than 0.001 mm leaves
            if n.shape is not None and n.srb is not None:
                # stem_diameter=0
                geom = leaf_mesh(
                    leaf=n.shape,
                    L_shape=n.mature_length,
                    Lw_shape=n.shape_max_width,
                    length=n.visible_length,
                    s_base=n.srb,
                    s_top=n.srt,
                    # flipx allows x-> -x to place the shape along
                    # with the tiller positioned with
                    # turtle.down()
                    flipx=True,
                    inclination=n.inclination, 
                    stem_diameter=n.stem_diameter,
                )

            ### Add a mesh that is senescent
            if n.senescent_length > 0.0001: # filter less than 0.001 mm leaves  # noqa: SIM102
                if n.srt is not None: # and n.shape is not None
                    geom_senescent = leaf_mesh(
                        leaf=n.shape,
                        L_shape=n.mature_length,
                        Lw_shape=n.shape_max_width,
                        length=n.visible_length,
                        s_base=n.srt,
                        s_top=1.,
                        # flipx allows x-> -x to place the shape along
                        # with the tiller positioned with
                        # turtle.down()
                        flipx=True,
                        inclination=n.inclination, 
                        stem_diameter=n.stem_diameter,
                    )
                    n.geometry_senescent = geom_senescent


            if n.lrolled > 0:
                rolled = stem_mesh(
                    n.lrolled, n.lrolled, n.d_rolled, classic
                )
                if geom is None:
                    geom = rolled
                else:
                    geom = addSets(rolled, geom, translate=(0, 0, n.lrolled))
    elif n.label.startswith("Stem"):  # stem element
        geom = stem_mesh(length=n.length, visible_length=n.visible_length, stem_diameter=n.stem_diameter, classic=classic)

    return geom


class CerealsVisitorConstrained(CerealsVisitor):
    def __init__(self, classic):
        super().__init__(classic)

    def __call__(self, g, v, turtle, time, growth):
        
        # 1. retrieve the node
        geoms_senesc = g.property("geometry_senescent")
        n = g.node(v)

        # Go to plant position if first plant element
        if n.parent() is None:
            turtle.move(0, 0, 0)
            # initial position to be compatible with canMTG positioning
            turtle.setHead(0, 0, 1, -1, 0, 0)

        # Manage inclination of tiller
        if g.edge_type(v) == "+" and not n.label.startswith("Leaf"):
            # axis_id = g.complex(vid); g.property('insertion_angle')
            # TODO : vary as function of age and species(e.g. rice)
            # print(n.label, n.visible_length, n.tiller_angle)
            angle = 2*n.tiller_angle if g.order(v) == 1 else n.tiller_angle 
            turtle.down(angle)
            turtle.elasticity = n.gravitropism_coefficient 
            turtle.tropism = (0, 0, 1)


        # incline turtle at the base of stems,
        if n.label.startswith("Stem"):
            azim = float(n.azimuth) if n.azimuth else 0.0
            if azim:
                # print 'node', n._vid, 'azim ', azim
                turtle.rollL(azim)

        # Try to build the mesh from GC rather than Cylinder
        if n.label.startswith("Leaf") or n.label.startswith("Stem"):
            # update geometry of elements
            # if n.length > 0:
            # print(v)
            # nb_of_growing_internodes, nb_of_growing_leaves = growing_organs(g, time)
            growth = compute_organ_growth(v, n, time, growth, self.classic)
            mesh = compute_growing_organ_geometry(n, self.classic)
            if mesh:  # To DO : reset to None if calculated so ?
                n.geometry = turtle.transform(mesh)
                n.anchor_point = turtle.getPosition()
                n.heading = turtle.getHeading()
            if v in geoms_senesc and geoms_senesc[v] is not None:
                geoms_senesc[v] = turtle.transform(geoms_senesc[v])

        # 3. Update the turtle and context
        turtle.setId(v)
        if n.label.startswith("Stem"):
            if n.visible_length > 0:
                turtle.f(n.visible_length)
            turtle.context.update({"top": turtle.getFrame()})
        if n.label.startswith("Leaf"):  # noqa: SIM102
            if n.lrolled > 0:
                turtle.f(n.lrolled)
                turtle.context.update({"top": turtle.getFrame()})

        return growth
    
