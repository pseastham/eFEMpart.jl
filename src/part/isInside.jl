using BenchmarkTools

include("ParticleTypes.jl")

"""
  isInside(polygon, n, p)

Returns true if the point p lies inside the polygon with n vertices
"""
function isInside(polygon::Vector{Point2D}, n::Int, p::Point2D; extreme = Point2D(100000.0, p[2]))
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
function doIntersect(p1::Point2D,q1::Point2D,p2::Point2D,q2::Point2D)
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
function orientation(p::Point2D,q::Point2D,r::Point2D)
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
function onSegment(p::Point2D,q::Point2D,r::Point2D)
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

# ---
# NEED TESTS
# ---

function onSegment_TEST()
  p = Point2D(0.0,0.0)
  r = Point2D(1.0,1.0)
  q = Point2D(0.5,0.5)

  #@code_warntype onSegment(p,q,r)
  @btime onSegment($p,$q,$r)

  return  onSegment(p,q,r)
end

function orientation_TEST()
  p   = Point2D(0.0,0.0)
  q   = Point2D(1.0,0.0)
  r   = Point2D(1.0,1.0)
  mid = Point2D(0.5,0.0)

  #@code_warntype orientation(p,q,r)
  @btime orientation($p,$q,$r)
  
  # CCW
  println("should be 2: ",orientation(p,q,r))

  # CW
  println("should be 1: ",orientation(r,q,p))

  # colinear
  println("should be 0: ",orientation(p,q,mid))

  nothing
end

function  doIntersect_TEST()
  p1 = Point2D(0.0,0.0)
  q1 = Point2D(1.0,0.0)
  p2 = Point2D(0.5,0.5)
  q2 = Point2D(0.5,-0.5)

  #@code_warntype doIntersect(p1,q1,p2,q2)
  @btime doIntersect($p1,$q1,$p2,$q2)

  nothing
end

function isInside_Point2D_TEST()
  n = 4
  p1 = Point2D(0.0,0.0)
  p2 = Point2D(1.0,0.0)
  p3 = Point2D(1.0,1.0)
  p4 = Point2D(0.0,1.0)
  p  = Point2D(rand(),rand())
  myextreme = [100000.0, p[2]]
  
  polygon = [p1,p2,p3,p4]

  #@code_warntype isInside(polygon, n, p)
  @btime isInside($polygon, $n, $p; extreme=$myextreme)

  return isInside(polygon,n,p)
end
function isInside_TEST()
  n = 4
  p1 = [0.0,0.0]
  p2 = [1.0,0.0]
  p3 = [1.0,1.0]
  p4 = [0.0,1.0]
  p = [50+rand(),rand()]
  myextreme = [100000.0, p[2]]
  
  polygon = [p1,p2,p3,p4]

  #@code_warntype isInside(polygon, n, p)
  @btime isInside($polygon, $n, $p; extreme=$myextreme)

  return isInside(polygon,n,p)
end