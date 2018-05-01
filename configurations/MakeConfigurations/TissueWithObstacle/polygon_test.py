import numpy as np

def inside_polygon(x, y, points):
    """
    Return True if a coordinate (x, y) is inside a polygon defined by
    a list of verticies [(x1, y1), (x2, x2), ... , (xN, yN)].

    Reference: http://www.ariel.com.au/a/python-point-int-poly.html
    """
    n = points.shape[0]
    inside = False
    p1x, p1y = points[0,:]
    for i in range(1, n + 1):
        p2x, p2y = points[i % n,:]
        if y > np.minimum(p1y, p2y):
            if y <= np.maximum(p1y, p2y):
                if x <= np.maximum(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y
    return inside

def area(points):
    i = np.arange(points.shape[0])
    j = np.roll(i,1)
    return 0.5*np.abs(np.sum(points[i,0]*points[j,1]-points[i,1]*points[j,0]))

def perim(points):
    i = np.arange(points.shape[0])
    j = np.roll(i,1)
    x, y = points[:,0], points[:,1]
    return np.sum(np.sqrt((x[i]-x[j])**2 + (y[i]-y[j])**2))

def normals(points):
    i = np.arange(points.shape[0])
    j = np.roll(i,1)
    k = np.roll(i,-1)
    x, y = points[:,0], points[:,1]
    nx = -((x[j]-x[i]) + (x[k]-x[i]))
    ny = -((y[j]-y[i]) + (y[k]-y[i]))
    nz = np.zeros(points.shape[0])
    len_n = np.sqrt(nx*nx + ny*ny)
    return np.vstack((nx/len_n,ny/len_n,nz)).T
