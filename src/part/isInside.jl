# functions for determining whether a point lies within a quadrilateral

include("ParticleTypes.jl")  # loads Point2D type

"""
  isInside(polygon, n, p)

Returns true if the point p lies inside the polygon with n vertices
"""
function isInside(polygon::Vector{Point2D{T}}, n::Int, p::Point2D{T}; extreme = Point2D(100000.0, p[2])) where T<:Real
  # There must be at least 3 vertices in polygon
  if (n < 3) return false end

  # Create a point for line segment from p to infinite
  

  # Count intersections of the above line with sides of polygon
  count = 0
  for i=1:n
    next = mod(i,n)+1
    # Check if the line segment from 'p' to 'extreme' intersects
    # with the line segment from 'polygon[i]' to 'polygon[next]'
    if (doIntersect(polygon[i], polygon[next], p, extreme))
      # If the point 'p' is colinear with line segment 'i-next',
      # then check if it lies on segment. If it lies, return true,
      # otherwise false
      if (orientation(polygon[i], p, polygon[next]) == 0)
         return onSegment(polygon[i], p, polygon[next])
      end
      count += 1
    end
  end

  # Return true if count is odd, false otherwise
  if isodd(count)
    return true
  else
    return false
  end
end

"""
  doIntersect(p1,q1,p2,q2) -> bool

Checks whether two lines, defined by (p1,q1) and (p2,q2) where pi,qi are 2-arrays.

Algorithm explained in https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/.

Note we have NOT implemented the special case for overlapping line segments, as that will never be the case in finite element applications
"""
function doIntersect(p1::Point2D{T},q1::Point2D{T},p2::Point2D{T},q2::Point2D{T}) where T<:Real
  # Find the four orientations needed for general case
  o1 = orientation(p1, q1, p2)
  o2 = orientation(p1, q1, q2)
  o3 = orientation(p2, q2, p1)
  o4 = orientation(p2, q2, q1)

  # General case
  if (o1 != o2 && o3 != o4)
    return true
  else
    return false
  end
end
function doIntersect(p1,q1,p2,q2)
  # Find the four orientations needed for general case
  o1 = orientation(p1, q1, p2)
  o2 = orientation(p1, q1, q2)
  o3 = orientation(p2, q2, p1)
  o4 = orientation(p2, q2, q1)

  # General case
  if (o1 != o2 && o3 != o4)
    return true
  else
    return false
  end
end

"""
  orientation(p,q,r) -> Int

Computes orientation of 3 (ordered) points p,q,r. where each is a 2-array
2 -> counterclockwise (CCW)
1 -> clockwise        (CW)
0 -> collinear

Algorithm is explained in https://www.geeksforgeeks.org/orientation-3-ordered-points/
"""
function orientation(p::Point2D{T},q::Point2D{T},r::Point2D{T}) where T<:Real
  val = (q.y-p.y)*(r.x-q.x) - (q.x-p.x)*(r.y-q.y);

  if (val == 0) return 0 end  # colinear

  return (val > 0) ? 1 : 2 # clock- or counterclock-wise
end

"""
  onSegment(p,q,r) -> bool

Given three colinear points p, q, r, the function checks if
point q lies on line segment 'pr'

Note that it ASSUMES that p,r,q are co-linear!
"""
function onSegment(p::Point2D{T},q::Point2D{T},r::Point2D{T}) where T<:Real
  maxX = (p.x >= r.x ? p.x : r.x)
  minX = (p.x >= r.x ? r.x : p.x)
  maxY = (p.y >= r.y ? p.y : r.y)
  minY = (p.y >= r.y ? p.y : r.y)

  if (q.x <= maxX && q.x >= minX && q.y <= maxY && q.y >= minY)
      return true
  else
    return false
  end
end