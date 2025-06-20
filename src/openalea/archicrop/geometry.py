"""Geometric primitives for cereals"""

from __future__ import annotations

from math import cos, pi, radians, sin

import numpy as np
from scipy.integrate import trapezoid

import openalea.plantgl.all as pgl
from openalea.mtg.turtle import TurtleFrame

from . import fitting

# Deactivate arrange
#LEAF_MOVE = False

def blade_elt_area(s, r, Lshape=1, Lwshape=1, sr_base=0, sr_top=1):
    """surface of a blade element, positioned with two relative curvilinear absisca"""

    S = 0
    sr_base = min([1, max([0, sr_base])])
    sr_top = min([1, max([sr_base, sr_top])])
    sre = [sr for sr in zip(s, r) if (sr_base < sr[0] < sr_top)]
    if len(sre) > 0:
        se, re = list(zip(*sre))
        snew = [sr_base, *list(se), sr_top]
        rnew = [np.interp(sr_base, s, r), *list(re), np.interp(sr_top, s, r)]
    else:
        snew = [sr_base, sr_top]
        rnew = [np.interp(sr_base, s, r), np.interp(sr_top, s, r)]

    S = trapezoid(rnew, snew) * Lshape * Lwshape

    return S  # noqa: RET504


def form_factor(leaf):
    _, _, s, r = leaf
    return blade_elt_area(s, r, 1, 1, 0, 1)


def compute_leaf_area(leaf, length=1, mature_length=1, width_max=1, form_factor=None):
    """
    Leaf area as a function of length
    -------

    Parameters
    ----------
    leaf: x,y,s,r coordinate describing mature leaf shape
    length: current length of the leaf
    mature_length: the length of the leaf once mature
    width_max: maximal width of the mature leaf
    form_factor: (optional) the form_factor of the mature leaf (if known), used to avoid integration

    Returns
    -------
    the area of the leaf corresponding to the distal part up to length (computed with trapeze aera)
    """
    if length >= mature_length and form_factor is not None:
        return length * width_max * form_factor
    _, _, s, r = leaf
    sr_b = 1 - float(length) / mature_length
    return blade_elt_area(s, r, mature_length, width_max, sr_base=sr_b, sr_top=1)


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


def arrange_leaf(leaf, stem_diameter=0, inclination=1, relative=True):
    """Arrange a leaf to be placed along a stem with a given inclination.

    Args:
        leaf: a x, y, s, r tuple describing leaf shape
        stem_diameter: the diameter of the sem at the leaf insertion point
        inclination: if relative=False, the leaf basal inclination (deg). A
        multiplier to leaf basal inclination angle otherwise
        relative: (bool) controls the meaning of inclination parameter

    Returns:
        a modified x, y, s, r tuple

    """
    #if not LEAF_MOVE:
    #    return leaf

    x, y, s, r = list(map(np.array, leaf))
    if relative and inclination == 1:
        x1, y1 = x, y
    else:
        basal_inclination = pgl.angle((x[1] - x[0], y[1] - y[0]), (0, 1))

        if relative:
            angle = inclination * basal_inclination
            angle = min(pi, angle)
        else:
            angle = radians(inclination)

        rotation_angle = basal_inclination - angle

        # rotation of the midrib
        cos_a = cos(rotation_angle)
        sin_a = sin(rotation_angle)

        x1 = x[0] + cos_a * x - sin_a * y
        y1 = y[0] + sin_a * x + cos_a * y
    leaf = x1 + stem_diameter / 2.0, y1, s, r

    return leaf  # noqa: RET504


def leaf_mesh(
    leaf,
    L_shape=1,
    Lw_shape=1,
    length=1,
    s_base=0,
    s_top=1,
    flipx=False,
    twist=0,
    volume=0,
    stem_diameter=0,
    inclination=1,
    relative=True,
):
    """Compute mesh for a leaf element along a scaled leaf shape

    Args:
        leaf: a x,y,s,r tuple describing leaf shape
        L_shape: length of the shape
        Lw_shape: width of the shape
        length: the total visible length to be meshed
        s_base: normalised position (on the shape) of the start of the element
        s_top: normalised position (on the shape) of the end of the element
        flipx: should leaf shape be flipped ?
        twist:
        volume: float value of the thickness of the leaf.
              Default is 0. Else it indicates the relative depth of the leaf
              along the midrib.
        stem_diameter: the diameter of the sem at the leaf insertion point
        inclination: if relative=False, the leaf basal inclination (deg). A
        multiplier to leaf basal inclination angle otherwise
        relative: (bool) controls the meaning of inclination parameter

    Returns:
        a PlanGl mesh representing the element
    """
    shape = arrange_leaf(
        leaf,
        stem_diameter=stem_diameter / L_shape,
        inclination=inclination,
        relative=relative,
    )

    # flip to position leaves along tiller emitted
    if flipx:
        # to position leaves along tiller emitted
        shape = (-shape[0],) + shape[1:]

    mesh = fitting.mesh4(
        shape, L_shape, length, s_base, s_top, Lw_shape, twist=twist, volume=volume
    )

    if mesh:
        pts, ind = mesh
        mesh = None if len(ind) < 1 else fitting.plantgl_shape(pts, ind)
    else:
        if length > 0:
            msg = f"Error : no mesh. s_base = {s_base}, s_top = {s_top}, length = {length}"
            raise ValueError(msg)
        mesh = None

    return mesh


# Meshing function for StemElements


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


def slim_cylinder(length, radius_base, radius_top):
    """
    Try to construct a cylinder with a low number of triangles.
    """
    rb, rt = radius_base, radius_top
    a1, a2, a3 = 0, 2 * pi / 3.0, 4 * pi / 3.0
    r = rb
    p1 = (r * cos(a1), r * sin(a1), 0)
    p2 = (r * cos(a2), r * sin(a2), 0)
    p3 = (r * cos(a3), r * sin(a3), 0)
    r = rt
    q1 = (r * cos(a1 + pi), r * sin(a1 + pi), length)
    q2 = (r * cos(a2 + pi), r * sin(a2 + pi), length)
    q3 = (r * cos(a3 + pi), r * sin(a3 + pi), length)

    return pgl.TriangleSet(
        [p1, p2, p3, q1, q2, q3],
        [
            (2, 1, 0),
            (3, 4, 5),
            (0, 5, 4),
            (0, 4, 2),
            (2, 4, 3),
            (3, 1, 2),
            (1, 3, 5),
            (5, 0, 1),
        ],
    )


def stem_mesh(length, visible_length, stem_diameter, classic=False, slices=24):  # noqa: ARG001
    """Compute mesh for a stem element
    - classic indicates
    """

    if classic:
        solid = True
        # 6 is the minimal number of slices for a correct computation of star
        #  (percentage error lower than 5)
        slices = 10
        stem = pgl.Tapered(
            stem_diameter / 2.0,
            stem_diameter / 2.0,
            pgl.Cylinder(1.0, visible_length, solid, slices),
        )
        tessel = pgl.Tesselator()
        stem.apply(tessel)
        mesh = tessel.triangulation
    else:
        mesh = slim_cylinder(visible_length, stem_diameter / 2.0, stem_diameter / 2.0)

    return mesh


def compute_element(element_node, classic=False):
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
                mesh = compute_element(n, self.classic)
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
