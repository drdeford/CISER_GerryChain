import math, random


# Data conventions: A point is a pair of floats (x, y). A circle is a triple of floats (center x, center y, radius).

# Returns the smallest circle that encloses all the given points. Runs in expected O(n) time, randomized.
# Input: A sequence of pairs of floats or ints, e.g. [(0,5), (3.1,-2.7)].
# Output: A triple of floats representing a circle.
# Note: If 0 points are given, None is returned. If 1 point is given, a circle of radius 0 is returned.
# 
# Initially: No boundary points known
def make_circle(points):
	# Convert to float and randomize order
	shuffled = [(float(x), float(y)) for (x, y) in points]
	random.shuffle(shuffled)
	
	# Progressively add points to circle or recompute circle
	c = None
	for (i, p) in enumerate(shuffled):
		if c is None or not is_in_circle(c, p):
			c = _make_circle_one_point(shuffled[ : i + 1], p)
	return c


# One boundary point known
def _make_circle_one_point(points, p):
	c = (p[0], p[1], 0.0)
	for (i, q) in enumerate(points):
		if not is_in_circle(c, q):
			if c[2] == 0.0:
				c = make_diameter(p, q)
			else:
				c = _make_circle_two_points(points[ : i + 1], p, q)
	return c


# Two boundary points known
def _make_circle_two_points(points, p, q):
	circ = make_diameter(p, q)
	left = None
	right = None
	px, py = p
	qx, qy = q
	
	# For each point not in the two-point circle
	for r in points:
		if is_in_circle(circ, r):
			continue
		
		# Form a circumcircle and classify it on left or right side
		cross = _cross_product(px, py, qx, qy, r[0], r[1])
		c = make_circumcircle(p, q, r)
		if c is None:
			continue
		elif cross > 0.0 and (left is None or _cross_product(px, py, qx, qy, c[0], c[1]) > _cross_product(px, py, qx, qy, left[0], left[1])):
			left = c
		elif cross < 0.0 and (right is None or _cross_product(px, py, qx, qy, c[0], c[1]) < _cross_product(px, py, qx, qy, right[0], right[1])):
			right = c
	
	# Select which circle to return
	if left is None and right is None:
		return circ
	elif left is None:
		return right
	elif right is None:
		return left
	else:
		return left if (left[2] <= right[2]) else right


def make_circumcircle(p0, p1, p2):
	# Mathematical algorithm from Wikipedia: Circumscribed circle
	ax, ay = p0
	bx, by = p1
	cx, cy = p2
	ox = (min(ax, bx, cx) + max(ax, bx, cx)) / 2.0
	oy = (min(ay, by, cy) + max(ay, by, cy)) / 2.0
	ax -= ox; ay -= oy
	bx -= ox; by -= oy
	cx -= ox; cy -= oy
	d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2.0
	if d == 0.0:
		return None
	x = ox + ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
	y = oy + ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
	ra = math.hypot(x - p0[0], y - p0[1])
	rb = math.hypot(x - p1[0], y - p1[1])
	rc = math.hypot(x - p2[0], y - p2[1])
	return (x, y, max(ra, rb, rc))


def make_diameter(p0, p1):
	cx = (p0[0] + p1[0]) / 2.0
	cy = (p0[1] + p1[1]) / 2.0
	r0 = math.hypot(cx - p0[0], cy - p0[1])
	r1 = math.hypot(cx - p1[0], cy - p1[1])
	return (cx, cy, max(r0, r1))


_MULTIPLICATIVE_EPSILON = 1 + 1e-14

def is_in_circle(c, p):
	return c is not None and math.hypot(p[0] - c[0], p[1] - c[1]) <= c[2] * _MULTIPLICATIVE_EPSILON


# Returns twice the signed area of the triangle defined by (x0, y0), (x1, y1), (x2, y2).
def _cross_product(x0, y0, x1, y1, x2, y2):
	return (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)


def _continuous_area(geo):
    """returns geo.area"""
    
    return geo.area

def _continuous_perimeter(geo):
    """returns geo.length"""
    
    return geo.length
from math import pi, sqrt
from shapely.geometry.polygon import orient

def perimeter(geo, geo_cell = None, gross = False):
    """
    Return perimeters of geometries in GeoSeries as Series of floats.
    
    Keyword arguments:
        geo -- GeoSeries or GeoDataFrame
        geo_cell -- GeoSeries or GeoDataFrame representing units used to build
            geo (the "container"); does not have to nest cleanly
        gross -- Calculate perimeter of geo by measuringstraight line distance 
            between nodes of geo_cell on boundary of geo (see Schwartzberg 1966)
            
    This function calculates continuous or discrete perimeter. 
    
    Continuous (Euclidean) perimeter is calculated if only geo argument is 
    provided. Currently this function just returns GeoSeries.length. 
    Future improvements could include:
        
        * Checking for lat-long coordinate system and performing geodetic
        measurement
        * Determining appropriate local CRS (most likely a State Plane or UTM
        zone) and performing calculation in that CRS.
        
    NOT YET OPERATIONALIZED: Discrete perimeter is calculated if a second
    geographic argument is provided that represents the "cells" or "building 
    blocks" of the first, larger geography.
    """

    if geo_cell is None:
        # Continuous perimeter
        return _continuous_perimeter(geo)
    elif gross:
        return _gross_perimeter(geo, geo_cell)
    else:
        return _discrete_perimeter(geo, geo_cell)

def area(geo, geo_cell = None, convex_hull = False):
    """
    Return areas of geometries in GeoSeries as Series of floats.
    
    Keyword arguments:
        geo -- GeoSeries or GeoDataFrame
        geo_cell -- GeoSeries or GeoDataFrame representing units used to build
            geo (the "container"); does not have to nest cleanly
        convex_hull -- Calculate area of convex hull of geo
        
    This function calculates continuous or area. 
    
    Continuous (Euclidean) area is calculated if only geo argument is 
    provided. Currently this function just returns GeoSeries.area. 
    Future improvements could include:
        
        * Checking for lat-long coordinate system and performing geodetic
        measurement
        * Determining appropriate local CRS (most likely a State Plane or UTM
        zone) and performing calculation in that CRS.
        
    NOT YET OPERATIONALIZED: Discrete area is calculated if a second
    geographic argument is provided that represents the "cells" or "building 
    blocks" of the first, larger geography.
    """

    if geo_cell is None:
        # Continuous area
        if convex_hull:
            return _continuous_area(geo.convex_hull)
        else:        
            return _continuous_area(geo)
    else:
        return _discrete_area(geo)
    
def polsby_popper_shape(geo, geo_cell = None):
    """
    Returns Polsby-Popper (1991) compactness of geo as float
    
    Keyword arguments:
        geo -- GeoSeries or GeoDataFrame
        geo_cell -- GeoSeries or GeoDataFrame representing units used to build
            geo (the "container"); does not have to nest cleanly
    """
    
    return 4 * pi * area(geo, geo_cell) / (perimeter(geo, geo_cell) ** 2)

def schwartzberg(geo, inverse = True, geo_cell = None):
    """
    Returns Schwartzberg (1966) compactness of geo as float
    
    Keyword arguments:
        geo -- GeoSeries or GeoDataFrame
        inverse -- Boolean, return the original (1 to infinity) Schwartzberg,
            or the inverse (0 to 1); most analysts use the inverse, so this
            function defaults to inverse = True
        geo_cell -- GeoSeries or GeoDataFrame representing units used to build
            geo (the "container"); does not have to nest cleanly
    """
    
    schw = polsby_popper(geo, geo_cell) ** -0.5
    
    if inverse:
        return schw ** -1
    else:
        return schw

def c_hull_ratio(geo):
    
    return area(geo) / area(geo, convex_hull = True)

def reock(geo):
    """
    Returns Reock (1961) compactness of geo as float
    
    Keyword arguments:
        geo -- GeoSeries or GeoDataFrame
    """
    
    # make_circle returns x, y, r of circle. Use r to calculate area.
    mbc_area = geo.convex_hull.apply(lambda x: pi * make_circle(list(x.exterior.coords))[2] ** 2)
    return geo.area / mbc_area

def _polar_moment_of_area(pts, centroid = None):
    """Returns the polar moment of area of a shape
    
    Keyword arguments:
        pts -- Iterable of tuples containing xy coordinate pairs of points
            that define a closed linear ring (i.e. a polygon boundary)
    
    The polar moment of area is the sum of the second moment of area with 
    respect to the x-axis (I_x) and the y-axis (I_y). The moments are 
    calculated with respect to the centroid of the shape.
    
    Adapted from https://leancrew.com/all-this/2018/01/python-module-for-section-properties/
    """
    
    x = [ c[0] for c in pts ]
    y = [ c[1] for c in pts ]

    s = 0
    sx = sy = 0
    sxx = syy = 0 #sxy = 0

    for i in range(len(pts) - 1):
        s += x[i]*y[i+1] - x[i+1]*y[i]
        sx += (x[i] + x[i+1])*(x[i]*y[i+1] - x[i+1]*y[i])
        sy += (y[i] + y[i+1])*(x[i]*y[i+1] - x[i+1]*y[i])
        sxx += (y[i]**2 + y[i]*y[i+1] + y[i+1]**2)*(x[i]*y[i+1] - x[i+1]*y[i])
        syy += (x[i]**2 + x[i]*x[i+1] + x[i+1]**2)*(x[i]*y[i+1] - x[i+1]*y[i])
        # sxy += (x[i]*y[i+1] + 2*x[i]*y[i] + 2*x[i+1]*y[i+1] + x[i+1]*y[i])*(x[i]*y[i+1] - x[i+1]*y[i])

    a = s/2
    if centroid:
        cx, cy = centroid[0], centroid[1]
    else:
        cx, cy = sx/(6*a), sy/(6*a)

    mi = (sxx/12 - a*cy**2) + (syy/12 - a*cx**2)  

    return mi

def _mass_moment_of_inertia(geo, geo_cell, wt):
    """
    Returns the mass moment of inertia of a shape
    
    Keyword arguments:
        geo -- GeoSeries or GeoDataFrame
        geo_cell -- GeoSeries or GeoDataFrame representing units used to build
            geo (the "container"); does not have to nest cleanly
        wt -- Str: Name of column in geo_cell to use as the weight of the unit 
            (e.g. population count) in the moment of inertia calculation
    
    Each geo_cell is assigned to a container geo based on a representative point,
    a point guaranteed to be in geo_cell. The mass of each geo_cell is assumed
    to fall at its centroid. If geo_cell does not nest cleanly (that is, it 
    overlaps neighboring geos), the mass is is nonetheless assigned entirely to
    one container geo.
    """

    # Copy geo_cell so that it is not modified during function    
    tmp = geo_cell[[wt, geo_cell.geometry.name]].copy()
    
    # Calculate area of each cell, used later
    tmp["dA"] = geo_cell.area
    
    # Assign container ID using representative point, since centroid of small
    # geometry may fall outside the geometry, and therefore outside the container
    # for small geometries near edge of container
    tmp.geometry = geo_cell.representative_point()
    tmp.crs = geo_cell.crs # For some reason representative_point() destroys CRS, reassign it
    tmp = gpd.sjoin(tmp, geo[[geo.geometry.name]]) # After sjoin, order of GeoDataFrame changes. Index is preserved.

    # Replace representative point with centroid, which is needed for MI
    tmp.geometry = geo_cell.centroid

    # Calculate population-weighted centroid of each container geo
    # index_right refers to index of each feature in container
    cx = tmp.groupby("index_right").apply(lambda row: (row.geometry.x * row[wt]).sum() / row[wt].sum()).rename("cx")
    cy = tmp.groupby("index_right").apply(lambda row: (row.geometry.y * row[wt]).sum() / row[wt].sum()).rename("cy")
    tmp = tmp.merge(cx, on = "index_right").merge(cy, on = "index_right")
    
    # Calculate moment of inertia of each container geo
    tmp["centroid_distance"] = ((tmp.geometry.x - tmp["cx"])**2 + (tmp.geometry.y - tmp["cy"])**2).apply(sqrt)
    mi = tmp.groupby("index_right").apply(lambda row: (row.centroid_distance**2 * row[wt]).sum()).rename("mi")
    
    # Calculate radius of circle with same area as geo
    circle_mi = tmp.groupby("index_right").apply(lambda row: row[wt].sum() * row["dA"].sum() / (2 * pi))
    
    # Calculate compactness as ratio of MI of circle with same area and equal mass distribution to MI of container
    c_mi = circle_mi / mi
    c_mi.rename("c_mi", inplace = True)
      
    return c_mi
    
def moment_of_inertia(geo, geo_cell = None, wt = 1):
    """
    Returns moment of inertia shape index (MacEachren 1985; Li, et al. 2013) of geo as float
    
    Keyword arguments:
        geo -- GeoSeries or GeoDataFrame
        geo_cell -- GeoSeries or GeoDataFrame representing units used to build
            geo (the "container"); does not have to nest cleanly
        wt -- Int, float, or str: The weight of the unit (e.g. population
            count) in the moment of inertia calculation; this is either a
            number applied as the weight of every unit, or the name of a 
            column in geo_cell with weights; if you want the units weighted
            by area, do not supply values for geo_cell or wt
            
    The moment of inertia of a shape changes as the shape is translated in
    coordinate space. The moment of inertia about the centroid is the 
    minimum moment of inertia for a given shape. The MI of a circle of the
    same area as the shape is divided by the MI of the shape to convert to a 
    "shape index".
    
    When geo_cell is missing, the value calculated is knowsn as the "area 
    moment of inertia" or "second moment of area". It is calculated from the
    polygon coordinates. As with other compactness measures, the final value 
    varies from 0 (least compact) to 1 (most compact, a circle).
    
    When geo_cell is supplied, the value calculated is the "mass moment of
    inertia" or, simply, the "moment of inertia". In this case each cell
    is weighted by some value, specified by the parameter wt. For demographic
    data, the "mass" will usually be the cell population. if the mass is the 
    cell area, this is equivalent to the second moment of inertia, and geo_cell
    should just be omitted. Importanly, in calculating the shape index, the 
    reference circle has the same area of the shape *and uniform density*.
    Therefore, unlike most other compactness measures, this shape index is *not* 
    constrained to be between 0 and 1. If the shape has mass (population)
    concentrated near the center of the shape, the shape can be more compact
    than the reference circle, and the return value will be greater than 1.
    """
    
    if geo_cell is None:

        # Initialize list to hold moments of inertia
        moments = []
        
        # Iterate geometry, which may be a multipolygon
        for geom in geo.geometry:
            
            # Initialize moment of inertia
            mi = 0
            
            # Find centroid of geom, possibly multipolygon
            centroid = (geom.centroid.x, geom.centroid.y)
            
            # geom will contain one or more polygons
            for poly in geom:

                # Give polygon positive orientation (clockwise(?)) so moment
                # of inertia is positive and holes are negative
                poly = orient(poly, sign = 1)
                                
                # Add moment of inertia for exterior points of polygon
                mi += _polar_moment_of_area(poly.exterior.coords, centroid)           
                for interior in poly.interiors:
                    
                    # Subtract moment of inertia for each hole in polygon
                    mi += _polar_moment_of_area(interior.coords, centroid)
                        
            # Convert to shape index by comparing to circle of equal area
            compactness_mi = geom.area**2 / (2 * pi * mi)
            
            moments.append(compactness_mi)
            
        return pd.Series(moments, index = geo.index)
    
    else:
        
        # wt is int or float, repeat value in "wt" column
        if type(wt) in [int, float, tuple, list, pd.Series]:
            geo_cell["wt"] = wt
            wt = "wt"

        return _mass_moment_of_inertia(geo, geo_cell, wt)
    
        
