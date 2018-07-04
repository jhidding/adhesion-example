import numpy as np
import cft
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import collections
from scipy.spatial import ConvexHull


def lagrangian_vertices(box, pot, t):
    """Get a grid of vertices as described in `box` and lift them to the
    given potential `pot` at time `t`."""
    q = np.indices(box.shape) * box.L/box.N - box.L/2
    z = np.sum(q**2, axis=0) - 2 * t * pot
    return np.concatenate([q, np.expand_dims(z, 0)], 0).reshape(box.dim+1, -1).T


def plot_regular_triangulation(ch, selection, xlim=None, ylim=None, ax=None):
    """Plot the regular triangulation.
    
    :param ch: Convex Hull (as returned by `scipy.spatial.ConvexHull`.)
    :param selection: indices of selected simplices.
    :param xlim: xlim of the axes
    :param ylim: ylim of the axes
    :param ax: the Axes instance to plot to, if left emtpy, a new figure
    is created.
    :return: None
    """
    x = ch.points[:,0]
    y = ch.points[:,1]
    t = ch.simplices[selection]

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
    ax.set_aspect('equal')
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    ax.triplot(x, y, t, linewidth=1.0)
    

def delaunay_areas(box, ch, selection):
    """Compute areas of selected simplices.
    
    :param ch: Convex Hull
    :param selection: indices of selected simplices.
    :return: ndarray of areas, having same shape as `selection`.
    """
    a, b, c = ch.points[ch.simplices[selection]][:,:,:box.dim].transpose([1,0,2])
    return np.abs(np.cross(a - b, c - b) / 2)


def delaunay_class(box, ch, selection, threshold):
    """Compute the classification of each simplex using the given threshold.
    
    :param ch: Convex hull
    :param selection: indices of selected simplices.
    :return: Number of edges with square length longer than the threshold for
    each simplex. Array of integer of same shape as selection.
    
    ======  =========
    number  class
    ======  =========
    0       void
    1       curto-parabolic point
    2       filament
    3       node
    ======  =========
    """
    pts = ch.points[ch.simplices[selection]][:,:,:box.dim]
    dists2 = ((pts - np.roll(pts, 1, axis=1))**2).sum(axis=2)
    return np.where(dists2 > threshold, 1, 0).sum(axis=1)


def voronoi_points(ch, selection):
    """Compute the dual vertices of the selected simplices."""
    return - ch.equations[selection][:,:2] / ch.equations[selection][:,2][:,None] / 2


def edges(box, ch, valid):
    """List all edges in the convex hull."""
    nb = np.zeros(shape=(len(valid), 2*(box.dim+1)), dtype=int)
    nb[:,1::2] = ch.neighbors[valid]  # neighbours index into simplices and equations
    nb[:,0::2] = valid[:,None]        # so does `valid`, we intersperse them to create pairs
    return np.unique(np.sort(nb.reshape([-1, 2]), axis=1), axis=0)


def edge_points(ch, edges):
    """Compute the dual vertices for all edges."""
    save = np.seterr(invalid = 'ignore', divide = 'ignore')
    pts = - ch.equations[:,:2] / ch.equations[:,2][:,None] / 2
    np.seterr(**save)
    return pts[edges]


def edge_length(ch, edges):
    """Get the length of each edge (in the Delaunay triangulation)."""
    # find the points common to both simplices, should always be two points
    # this operation performs a loop in Python TODO: find numpy expression
    edge_verts = np.array([np.intersect1d(x[0], x[1]) for x in ch.simplices[edges]])
    return np.sqrt(np.sum((ch.points[edge_verts][:,1,:2] - ch.points[edge_verts][:,0,:2])**2, axis=1))


def plot_power_diagram(box, ch, valid, xlim, ylim, ax=None, point_scale=10):
    """Plot the power diagram."""
    m_edges = edges(box, ch, valid)
    m_edge_lengths = edge_length(ch, m_edges)
    m_edge_points = edge_points(ch, m_edges)

    edge_sel = np.where(m_edge_lengths > np.sqrt(2)*box.res)[0]
    
    lc_grid = collections.LineCollection(m_edge_points, linewidths=1.0, color='#888888')
    lc = collections.LineCollection(m_edge_points[edge_sel], linewidths=m_edge_lengths[edge_sel], color='maroon')
    lc.set_capstyle('round')

    X = voronoi_points(ch, valid)
    mass = delaunay_areas(box, ch, valid)
    big_points = np.where(delaunay_class(box, ch, valid, threshold=1.0) > 2)

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
    ax.add_collection(lc_grid)
    ax.add_collection(lc)
    ax.set_aspect('equal')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.scatter(X[big_points,0], X[big_points,1], s=mass[big_points]**1.5 * point_scale * 1.25, zorder=2, c='black', alpha=0.5)
    ax.scatter(X[big_points,0], X[big_points,1], s=mass[big_points]**1.5 * point_scale, zorder=4, alpha=0.5)

    
def get_convex_hull(box, pot_0, t):
    """Get the Convex hull, selection and valid simplices for a given box,
    potential and time.
    
    The selection checks for each simplex that it is pointed upward. We do not whish
    to get non-physical simplices that face downward, as they are a remnant of using
    the Convex Hull directly to do our computation.
    
    The 'valid' simplices are those that are in the selection and also have no neighbours
    outside the selection."""
    pts = lagrangian_vertices(box, pot_0, t)
    ch = ConvexHull(pts)
    selection = np.where(np.dot(ch.equations[:,0:3], [0, 0, -1]) > 0.00001)[0]
    valid = selection[np.where(np.all(np.isin(ch.neighbors[selection], selection), axis=1))[0]]
    return ch, selection, valid


def plot_image(box, data, title=None, ax=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    ax.set_title(title)
    ax.imshow(data, extent=[-box.L/2, box.L/2, -box.L/2, box.L/2], interpolation='nearest', cmap='RdYlBu')