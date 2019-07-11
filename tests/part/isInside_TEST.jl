using BenchmarkTools

include("../../src/part/isInside.jl")  # loads functions to be tested

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