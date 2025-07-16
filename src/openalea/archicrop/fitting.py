from __future__ import annotations

from heapq import *  # noqa: F403

import numpy as np
from scipy.integrate import simpson, trapezoid
from scipy.interpolate import splev, splprep

import openalea.plantgl.all as pgl
from openalea.plantgl.all import Vector3


def curvilinear_abscisse(x, y, z=None):
    """Curvilinear abcissa along a polyline"""
    s = np.zeros(len(x))
    if z is None:
        s[1:] = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)
    else:
        s[1:] = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2 + np.diff(z) ** 2)
    return s.cumsum()


def fit2(x, y, s, r):
    """
    x, y: 2d points
    s: curvilinear abscisse
    r: leaf radius

    Algo:
        1.1 fit x, y
        1.2 discretize the splines (high: 100)
        1.3 compute curvilinear abscisse
        1.4 normalize :
          + compute length and divide x, y by length

        2.1 fit r, s
        2.2 discretize r, s (high 500)
        2.3 normalize :
          + compute max r, s and divide r, s
        2.4 compute leaf surface
        2.5 r_t1-> found r(s): interp(s_t1,s_t2, r_t2)

        3. fit x(t1), y(t1), r(t1)
        4. return (x, y, r) and surface value
    """
    # spline parameters
    smooth = len(x) - np.sqrt(2 * len(x))  # smoothness parameter
    smooth = 0.01

    k = 3  # spline order
    nest = -1  # estimate of number of knots needed (-1 = maximal)

    try:
        tckp, u = splprep([x, y], s=smooth, k=k, nest=nest)
    except:  # noqa: E722
        try:
            tckp, u = splprep([x, y])
        except:  # noqa: E722
            tckp, u = splprep([x, y], k=1)

    xnew, ynew = splev(np.linspace(0, 1, 100), tckp)
    snew = curvilinear_abscisse(xnew, ynew)
    # 1.4
    length = snew.max()
    snew /= length
    xnew /= length
    ynew /= length

    # 2.1
    try:
        tckp2, v = splprep([s, r], s=smooth / 100.0, k=k, nest=nest)
    except:  # noqa: E722
        tckp2, v = splprep([s, r])

    snew2, rnew2 = splev(np.linspace(0, 1, 200), tckp2)

    snew2 = snew2 / snew2.max()
    rnew2 = np.maximum(rnew2, 0) / rnew2.max()

    # 2.4 leaf surface (integral r ds)
    leaf_surface = simpson(rnew2, snew2)

    # 2.5 Compute rnew
    rnew = np.interp(snew, snew2, rnew2)
    return (xnew, ynew, snew, rnew), leaf_surface


def discretize(leaf_spline, nb_polygones, length_max, radius_max):
    """
    Compute a mesh from a leaf represented by a 3d spline containing
    x, y and radius.

    Final radius have to be null.
    """
    return partial_leaf(leaf_spline, nb_polygones, length_max, length_max, radius_max)


def partial_leaf(leaf_spline, nb_polygones, length_max, length, radius_max):
    """
    Compute a mesh from a leaf represented by a 3d spline containing
    x, y and radius.
    Final radius have to be null.
    The leaf is obtain by evaluating the central curve on a domain
    [0, l/lmax] and the radius on a domain [1-l/lmax, 1].

    """
    if length <= 0.0:
        return None
    
    length = min(length, length_max)

    param = float(length) / length_max

    xf, yf, rnull = splev(np.linspace(0, param, nb_polygones), leaf_spline)
    xnull, ynull, rf = splev(np.linspace(1 - param, 1, nb_polygones), leaf_spline)
    rf = np.where(rf > 0.0, rf, np.zeros(len(rf)))
    rf[-1] = 0.0

    xf *= length_max
    yf *= length_max
    rf *= radius_max

    # build a mesh
    n = len(xf)
    points = list(zip(xf, yf, -rf / 2.0))
    points.extend(list(zip(xf, yf, rf / 2.0)))

    ind = np.array(range(n - 2))
    indices = list(zip(ind, ind + n, ind + (n + 1)))
    indices.extend(list(zip(ind, ind + (n + 1), ind + 1)))
    # add only one triangle at the end !!
    indices.append((n - 2, 2 * n - 2, n - 1))

    return plantgl_shape(points, indices)


def leaf_to_mesh_2d(x, y, r, twist_start=0, twist_end=0):
    n = len(x)
    theta = np.linspace(np.radians(twist_start), np.radians(twist_end), n)
    points = list(
        zip(x, -r / 2.0 * abs(np.cos(theta)), y + abs(np.sin(theta)) * r / 2.0)
    )
    points.extend(
        list(zip(x, r / 2.0 * abs(np.cos(theta)), y - abs(np.sin(theta)) * r / 2.0))
    )

    ind = np.array(range(n - 2)) if n > 2 else np.array([0])
    indices = list(zip(ind, ind + n, ind + (n + 1)))
    indices.extend(list(zip(ind, ind + (n + 1), ind + 1)))

    # add only one triangle at the end !!
    if n < 2:
        assert len(points) == 4
        if r[-1] < 0.001:
            indices = indices[0:1]
    elif r[-1] < 0.001:
        indices.append((n - 2, 2 * n - 2, 2 * n - 1))
    else:
        indices.append((n - 2, 2 * n - 2, 2 * n - 1))
        indices.append((n - 2, 2 * n - 1, n - 1))

    return points, indices


def leaf_to_mesh(x, z, w, twist_start=0, twist_end=0, volume=0.1):
    """
    Convert a leaf polyline (x, z position) with this width and the angle
    twist and volume to a mesh.

    Args:
        x: x positions of the leaf polyline
        z: z positions of the leaf polyline
        w: width values at each point of the leaf polyline
        twist_start: angle in degree
        twist_end: angle in degree
        volume: floaf value of the thickness of the leaf.

    Returns:
        vertices, faces

        vertices is a list of the point position (x, y, z) of the mesh
        faces is a list of tuple of 3 indice value (vertices)
    """

    if volume == 0:
        return leaf_to_mesh_2d(x, z, w, twist_start, twist_end)
    n = len(x)

    theta = np.linspace(np.radians(twist_start), np.radians(twist_end), n)

    def rotate(arr, theta):
        arr[:, 1] = arr[:, 1] * np.cos(theta)
        arr[:, 2] = arr[:, 1] * np.sin(theta) + arr[:, 2]
        return arr

    def normalize(v):
        norm = np.linalg.norm(v)
        if norm == 0:
            return v
        return v / norm

    y = [0] * n
    pts = np.array(list(zip(x, y, z)))

    # Compute volume vector
    v = []
    for x0, z0, x1, z1 in zip(x[:-1], z[:-1], x[1:], z[1:]):
        vv = np.array([-(z1 - z0), 0, x1 - x0])
        v.append(normalize(vv))
    v.append((0, 0, 0))
    v = np.array(v) * volume

    # Compute width leaf
    width = np.zeros_like(pts)
    width[:, 1] = np.array(w) / 2.0
    width[-1, 1] = 0

    # Computes vertices of the mesh
    vertices = (
        list(rotate(pts - width + v, theta))
        + list(rotate(pts + width + v, theta))
        + list(rotate(pts + width - v, theta))
        + list(rotate(pts - width - v, theta))
    )

    # Computes faces of the mesh
    ind = np.array(range(n - 2)) if n > 2 else np.array([0])

    faces = list(zip(ind, ind + n, ind + n + 1))
    faces.extend(list(zip(ind, ind + n + 1, ind + 1)))

    faces.extend(list(zip(ind + n, ind + 2 * n, ind + 2 * n + 1)))
    faces.extend(list(zip(ind + n, ind + 2 * n + 1, ind + n + 1)))

    faces.extend(list(zip(ind + 2 * n, ind + 3 * n, ind + 3 * n + 1)))
    faces.extend(list(zip(ind + 2 * n, ind + 3 * n + 1, ind + 2 * n + 1)))

    faces.extend(list(zip(ind + 3 * n, ind, ind + 1)))
    faces.extend(list(zip(ind + 3 * n, ind + 1, ind + 3 * n + 1)))

    # close extremity of the leaf
    faces.append((n - 2, 2 * n - 2, n - 1))
    faces.append((2 * n - 2, 3 * n - 2, n - 1))
    faces.append((3 * n - 2, 4 * n - 2, n - 1))
    faces.append((4 * n - 2, n - 2, n - 1))

    # close base of the leaf
    faces.append((0, n, 2 * n))
    faces.append((2 * n, 3 * n, 0))

    return vertices, faces


def leaf_element(leaf, length_max=1, length=1, s_base=0, s_top=1, radius_max=1):
    def insert_values(a, values):
        l = a.tolist()  # noqa: E741
        l.extend(values)
        return np.unique(l)

    s_base = min(s_base, s_top, 1.0)
    s_top = max(s_base, s_top, 0.0)

    if length <= 0.0:
        return None
    
    length = min(length, length_max)

    x, y, s, r = leaf
    # just scale leaf case
    if length == length_max and s_base == 0 and s_top == 1:
        x *= length_max
        y *= length_max
        s *= length_max
        r *= radius_max
        return x, y, s, r

    # 1. compute s_xy and s_r for length vs length_max

    # force the leaf width to zero at the top
    r[-1] = 0

    param = float(length) / length_max
    # n = len(s)
    # add param in s_xy
    sp = insert_values(s, [1 - param, param])
    s_xy = np.compress(sp <= param, sp)
    s_r = np.compress(sp >= (1 - param), sp)
    s_xy = np.union1d(s_xy, s_r - s_r[0])
    s_r = s_xy + s_r[0]
    # add s_base and s_top in s

    s_base *= param
    s_top *= param
    s_valid = insert_values(s_xy, [s_base, s_top])
    s_valid = np.compress(s_valid <= s_top, s_valid)
    s_valid = np.compress(s_valid >= s_base, s_valid)

    # delete small intervals COMMIT from  Here !
    eps = (s_top - s_base) / (len(s) * 2)
    ds = s_valid[1:] - s_valid[:-1]
    error = ds >= eps
    s_val = np.compress(error, s_valid[1:])
    s_val = insert_values(s_val, [s_valid[0], s_valid[-1]])

    # find param and 1-param in s...
    # build xf, yf, rf with the same nb of points

    s_xyf = s_val
    s_rf = s_val + s_r[0]

    xf = np.interp(s_xyf, s, x)
    yf = np.interp(s_xyf, s, y)
    rf = np.interp(s_rf, s, r)

    cond = rf <= 0
    if cond.any():
        # delete points with negative radius
        index = cond.searchsorted(True)
        xf = xf[: index + 1]
        yf = yf[: index + 1]
        rf = rf[: index + 1]
        rf[-1] = 0.0

    xf *= length_max
    yf *= length_max
    rf *= radius_max

    return xf, yf, s_val, rf


def mesh4(leaf, length_max, length, s_base, s_top, radius_max, twist=0, volume=0.1):
    xf, yf, s_val, rf = leaf_element(
        leaf, length_max, length, s_base, s_top, radius_max
    )

    if len(xf) < 2:
        # All the radius are negative or null.
        # Degenerated element.
        return [], []

    pts, ind = leaf_to_mesh(
        xf,
        yf,
        rf,
        twist_start=twist * min(s_val),
        twist_end=twist * max(s_val),
        volume=volume,
    )

    # def _surf(ind, pts):
    #     from openalea.plantgl.all import norm, cross, Vector3
    #     A, B, C = [Vector3(pts[i]) for i in ind]
    #     return norm(cross(B - A, C - A)) / 2.0
    #
    # ind = [id for id in ind if _surf(id, pts) > 1e-6]

    return pts, ind


def write_smf(filename, points, indices):
    """
    Write a smf file for QSlim software.
    See http://people.scs.fsu.edu/~burkardt/data/smf/smf.html
    """

    header = """#$SMF 1.0
#$vertices %d
#faces %d
#

""" % (len(points), len(indices))  # noqa: UP031

    pts_str = "\n".join(["v %f %f %f" % (pt[0], pt[1], pt[2]) for pt in points])  # noqa: UP031
    ind_str = "\n".join(["f %d %d %d" % (ind[0] + 1, ind[1] + 1, ind[2] + 1) for ind in indices])  # noqa: UP031

    f = open(filename, "w")  # noqa: PTH123, SIM115
    f.write(header)
    f.write(pts_str)
    f.write("\n")
    f.write(ind_str)

    f.close()


def read_smf(filename):
    """
    Read a smf file containing a mesh.
    Returns the point set and the face set.
    """
    f = open(filename)  # noqa: PTH123, SIM115
    points = []
    indices = []

    for line in f:
        line.strip()
        if not line:
            continue
        if line[0] == "v":
            pts = line.split()[1:]
            points.append([float(pt) for pt in pts])
        elif line[0] in ["f", "t"]:
            ind = line.split()[1:]
            indices.append([int(i) - 1 for i in ind])
        else:
            continue

    f.close()

    return points, indices


def plantgl_shape(points, indices):
    indices = [list(map(int, index)) for index in indices]
    return pgl.TriangleSet(points, indices, normalPerVertex=False)


def fit3(x, y, s, r, nb_points):
    leaf, leaf_surface = fit2(x, y, s, r)  # use here simpson as leaf is a smooth shape
    xn, yn, sn, rn = simplify(leaf, nb_points)
    return xn, yn, sn, rn



def max_distance(pts, line):
    d_line = line.__normSquared__()
    max_dist = 0.0
    index = 0
    for i, pt in enumerate(pts):
        d = (Vector3(pt) ^ line).__normSquared__()
        if d > max_dist:
            max_dist = d
            index = i

    return index, max_dist / d_line


def distance(pt, p0, p1):
    line = p1 - p0
    length = line.__normSquared__()
    d = ((pt - p0) ^ line).__normSquared__()
    d /= length
    return d


def cost(polyline, nb_points):
    nb_points += 2
    n = len(polyline)
    pts = list(polyline)
    _cost = []

    sibling = [[i - 1, i + 1] for i in range(n)]

    # compute the cost for each points
    heap_cost = []
    _cost = {}

    for i in range(1, n - 1):
        d = distance(polyline[i], polyline[i - 1], polyline[i + 1])
        _cost[i] = d
        heappush(heap_cost, [d, i])  # noqa: F405

    while len(heap_cost) > nb_points - 2:
        d, i = heappop(heap_cost)  # noqa: F405

        assert pts[i] is not None
        pts[i] = None

        # update i-1 and i+1 distance
        il, ir = sibling[i]

        if il != 0:
            sibling[il][1] = ir
            ill = sibling[il][0]
            dl = distance(pts[il], pts[ill], pts[ir])
            old_dl = _cost[il]
            heap_index = heap_cost.index([old_dl, il])
            _cost[il] = dl
            del heap_cost[heap_index]
            heappush(heap_cost, [dl, il])  # noqa: F405
        if ir != n - 1:
            sibling[ir][0] = il
            irr = sibling[ir][1]
            dr = distance(pts[ir], pts[il], pts[irr])
            old_dr = _cost[ir]
            heap_index = heap_cost.index([old_dr, ir])
            _cost[ir] = dr
            del heap_cost[heap_index]
            heappush(heap_cost, [dr, ir])  # noqa: F405
    # bug in the algorithm...
    pts[1] = None
    pts[-2] = None
    return pts



def simplify(leaf, nb_points, scale_radius=True):
    """ " Simplify a 2d polyline up to nb_points
    Parameters
    ----------
    leaf : a x, y, s, r tuple of iterables representing the 2d polyline to be simplified
    nb_points : the number of points along the polyline after simplification
    scale_radius : (bool, default=True)  should radius be scaled after simplification to ensure
        polygonial_area(simplified_leaf) = smooth_area(input_leaf) ?

    Returns
    -------
    a x, y, s, r tuple of iterables for the simplified 2d polyline
    """
    xn, yn, sn, rn = leaf

    points = [pgl.Vector3(*pt) for pt in zip(xn, rn, yn)]
    pts = [_f for _f in cost(points, nb_points) if _f]
    coords = ((pt.x, pt.y, pt.z) for pt in pts)
    x, r, y = list(map(np.array, zip(*coords)))
    s = curvilinear_abscisse(x, y)
    # keep smax similar to sn
    adj = max(sn) / max(s)
    s *= adj
    x *= adj
    y *= adj

    if scale_radius:
        leaf_surface = simpson(rn, sn)  # use here simpson as leaf is a smooth shape
        new_surface = trapezoid(
            r, s
        )  # use here trapezoid as the new surface is a polygonial shape
        scale_r = leaf_surface / new_surface
        r *= scale_r

    return x, y, s, r

