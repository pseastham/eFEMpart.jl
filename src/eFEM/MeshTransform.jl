# File loaded into eFEMpart

"""
  linearInterp()

linear interpolates the point x0 which lies in (x1,x2) onto the domain (s1,s2)
"""
function linearInterp(x0,x1,x2,s1,s2)
  A = [x1 1; x2 1]; rhs=[s1,s2]
  mbtemp = A\rhs
  m = mbtemp[1]
  b = mbtemp[2]

  s0 = m*x0 + b

  return s0
end

"""
  onSegment(p,q,r) -> bool

Given three colinear points p, q, r, the function checks if
point q lies on line segment 'pr'
"""
function onSegment(p,q,r)
  if (q[1] <= maximum([p[1], r[1]]) && q[1] >= minimum([p[1],r[1]]) &&
          q[2] <= maximum([p[2], r[2]]) && q[2] >= minimum([p[2], r[2]]))
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
function orientation(p,q,r)
    val = (q[2] - p[2]) * (r[1] - q[1]) -
              (q[1] - p[1]) * (r[2] - q[2]);

    if (val == 0) return 0 end  # colinear

    return (val > 0) ? 1 : 2 # clock or counterclock wise
end


"""
  doIntersect(p1,q1,p2,q2) -> bool

Checks whether two lines, defined by (p1,q1) and (p2,q2) where pi,qi are 2-arrays.

Algorithm explained in https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/.

Note we have not implemented the special case for overlapping line segments, as that will never be the case in finite element applications
"""
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
  isInside(polygon, n, p)

Returns true if the point p lies inside the polygon with n vertices
"""
function isInside(polygon, n, p)
  # There must be at least 3 vertices in polygon
  if (n < 3)  return false end

  # Create a point for line segment from p to infinite
  extreme = [100000.0, p[2]]

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
  reverseQuadMap(quad,point) -> (s,t)

Computes the (s,t) computational domain point for an arbitrary point (x,y) in the physical element quad

For use with meshTransform
"""
function reverseQuadMap(quad,point)
  PRINT = false   # for testing

  x0 = point[1];   y0 = point[2]
  x1 = quad[1][1]; y1 = quad[1][2]
  x2 = quad[2][1]; y2 = quad[2][2]
  x3 = quad[3][1]; y3 = quad[3][2]
  x4 = quad[4][1]; y4 = quad[4][2]

  if PRINT
    println("(x0,y0) = (",x0,",",y0,")")
    println()
    println("(x1,y1) = (",x1,",",y1,")")
    println("(x2,y2) = (",x2,",",y2,")")
    println("(x3,y3) = (",x3,",",y3,")")
    println("(x4,y4) = (",x4,",",y4,")")
    println()
  end

  # define ci
  c1 =  x1 + x2 + x3 + x4 - 4*x0
  c5 =  y1 + y2 + y3 + y4 - 4*y0

  c2 = -x1 + x2 + x3 - x4
  c6 = -y1 + y2 + y3 - y4

  c3 = -x1 - x2 + x3 + x4
  c7 = -y1 - y2 + y3 + y4

  c4 =  x1 - x2 + x3 - x4
  c8 =  y1 - y2 + y3 - y4

  if PRINT
    println("c1 = ",c1)
    println("c2 = ",c2)
    println("c3 = ",c3)
    println("c4 = ",c4)
    println("c5 = ",c5)
    println("c6 = ",c6)
    println("c7 = ",c7)
    println("c8 = ",c8)
    println()
  end

  s = 0.0
  t = 0.0

  # case 1
  if abs(c4) < 1e-12 && abs(c8) < 1e-12
    #println("case 1")
    A = [c2 c3; c6 c7]; b = [-c1;-c5]

    x = A\b
    s = x[1]; t = x[2]

  # case 2A -- NOT CHECKED YET
  elseif abs(c4) < 1e-12
    #println("case 2A")
    if PRINT println("Case 2A!") end

    # c2 != 0
    if abs(c2) > 1e-12
      A = -c3*c8
      B = c2*c7 - c3*c6
      C = c2*c5 - c1*c6
      det = B^2 - 4*A*C

      println("A = ",A)
      println("B = ",B)
      println("C = ",C)

      if abs(A) < 1e-12
        t = -C/B

      else
        if PRINT println("det = ",det) end
        t1 = 0.5*(-B + sqrt(det))/A
        t2 = 0.5*(-B - sqrt(det))/A
        println("t1 = ",t1)
        println("t2 = ",t2)
        println()

        t = t1
      end

      s = -(c1+c3*t)/c2

    # c2==0
    else
      A = c3*c8
      B = c3*c6 - c2*c7 + c1*c8
      C = c1*c6 - c2*c5
      det = B^2 - 4*A*C

      println("A = ",A)
      println("B = ",B)
      println("C = ",C)

      if PRINT println("det = ",det) end

      t1 = 0.5*(-B + sqrt(det))/A
      t2 = 0.5*(-B - sqrt(det))/A
      println("t1 = ",t1)
      println("t2 = ",t2)
      println()

      s1 = -(c5+c7*t1)/(c6+c8*t1)
      s2 = -(c5+c7*t2)/(c6+c8*t2)

      println((c5+c7*t2))

      println("s1 = ",s1)
      println("s2 = ",s2)

      if (s1>1.0 || s1<-1.0)
        t = t2
        s = s2
      else
        t = t1
        s = s1
      end
    end
  # case 2B ---- NOT CHECKED YET
  elseif abs(c8) < 1e-12
    #println("case 2B")
    if PRINT println("Case 2B!") end
    # c6 != 0
    if abs(c6) > 1e-12
      A = -c7*c4
      B = c6*c3 - c7*c2
      C = c6*c1 - c5*c2
      det = B^2 - 4*A*C

      println("A = ",A)
      println("B = ",B)
      println("C = ",C)

      if abs(A) < 1e-12
        t = -C/B

      else
        if PRINT println("det = ",det) end
        t1 = 0.5*(-B + sqrt(det))/A
        t2 = 0.5*(-B - sqrt(det))/A
        println("t1 = ",t1)
        println("t2 = ",t2)
        println()

        t = t1
      end

      s = -(c1+c3*t)/c2

    # c6==0
    else
      A = c7*c4
      B = c7*c2 - c6*c3 + c5*c4
      C = c5*c2 - c6*c1
      det = B^2 - 4*A*C

      println("A = ",A)
      println("B = ",B)
      println("C = ",C)

      if PRINT println("det = ",det) end

      t1 = 0.5*(-B + sqrt(det))/A
      t2 = 0.5*(-B - sqrt(det))/A
      println("t1 = ",t1)
      println("t2 = ",t2)
      println()

      s1 = -(c1+c3*t1)/(c2+c4*t1)
      s2 = -(c1+c3*t2)/(c2+c4*t2)

      println((c1+c3*t2))

      println("s1 = ",s1)
      println("s2 = ",s2)

      if (s1>1.0 || s1<-1.0)
        t = t2
        s = s2
      else
        t = t1
        s = s1
      end
    end
  # case 3 ----- NOT CHECKED YET
  else
    #println("case 3")
    if PRINT println("Case 3!") end
    # define A,B,C
    A = c4*c7 - c3*c8                   # t^2
    B = c4*c5 - c3*c6 + c2*c7 - c1*c8   # t
    C = c2*c5 - c1*c6                   # 1

    # solve for solution t
    if abs(A)<1e-12
      t = -C/B
    else
      det = sqrt(B^2 - 4*A*C)
      t1 = 0.5*(-B + det)/A
      t2 = 0.5*(-B - det)/A

      if t1 > 1.0 || t1 < -1.0
        t = t2
      else
        t = t1
      end
    end

    # solve for solution s
    s = -(c1+c3*t)/(c2+c4*t)
  end

  # return s,t
  if PRINT println(s," ",t) end
  return s,t
end

"""
  isInsideDomain(xyA,cmA,point) -> bool

Transforms finite element solutions from domain A to point xyB, returns solution Uinterp at point xyB

Can handle 4 point (linear) and 9 point (quadratic) quadrilateral elements
"""
function isInsideDomain(xyA,cmA,point)
  nptsA   = size(xyA,1)
  nelA    = size(cmA,1)
  npe     = size(cmA,2); npe == 4 ? 4 : 8 # adjust for only 8 bdry pts
  order   = (npe==4 ? 1 : 2)
  Uinterp = 0.0
  foundInside = false
  foundVal = false

  # check over all nodes in Q1 for Q2 match
  #nodeTOL = 1e0
  #for ti=1:nptsA
  #    if sqrt((xyA[ti,1]-point[1])^2 + (xyA[ti,2]-point[2])^2) < nodeTOL
  #      println("found!")
  #      Uinterp = UA[ti]
  #      foundVal = true
  #      break
  #    end
  #end

  # loop over elements in A to find which one contains above point
  if foundVal==false
    for el=1:nelA
      if foundInside == false
        c = cmA[el,1:4]
        p1 = [xyA[c[1],1],xyA[c[1],2]]
        p2 = [xyA[c[2],1],xyA[c[2],2]]
        p3 = [xyA[c[3],1],xyA[c[3],2]]
        p4 = [xyA[c[4],1],xyA[c[4],2]]
        quad  = [p1,p2,p3,p4]

        # check whether point is in element
        if isInside(quad,4,point)
          return true
        # check whether point is on side of element
        else
          if onSegment(p1,point,p2)
            return true
          elseif onSegment(p2,point,p3)
            return true
          elseif onSegment(p3,point,p4)
            return true
          elseif onSegment(p4,point,p1)
            return true
          end     
        end
      end
    end
  end

  return false
end


"""
  pointTransform(xyA,cmA,UA,point) -> Uinterp

Transforms finite element solutions from domain A to point xyB, returns solution Uinterp at point xyB

Can handle 4 point (linear) and 9 point (quadratic) quadrilateral elements
"""
function pointTransform(xyA,cmA,UA,point)
  nptsA   = size(xyA,1)
  nelA    = size(cmA,1)
  npe     = size(cmA,2); npe == 4 ? 4 : 8 # adjust for only 8 bdry pts
  order   = (npe==4 ? 1 : 2)
  Uinterp = 0.0
  foundInside = false
  foundVal = false

  # check over all nodes in Q1 for Q2 match
  #nodeTOL = 1e0
  #for ti=1:nptsA
  #    if sqrt((xyA[ti,1]-point[1])^2 + (xyA[ti,2]-point[2])^2) < nodeTOL
  #      println("found!")
  #      Uinterp = UA[ti]
  #      foundVal = true
  #      break
  #    end
  #end

  # loop over elements in A to find which one contains above point
  if foundVal==false
    for el=1:nelA
      if foundInside == false
        c = cmA[el,1:4]
        p1 = [xyA[c[1],1],xyA[c[1],2]]
        p2 = [xyA[c[2],1],xyA[c[2],2]]
        p3 = [xyA[c[3],1],xyA[c[3],2]]
        p4 = [xyA[c[4],1],xyA[c[4],2]]
        quad  = [p1,p2,p3,p4]

        # check whether point is in element
        if isInside(quad,4,point)
          s,t = reverseQuadMap(quad,point)

          #println("s=",round(s,2))
          #println("t=",round(t,2))

          phi,_,_ = shape2D(s,t,1)

          x = shapeEval([p1[1],p2[1],p3[1],p4[1]],phi)
          y = shapeEval([p1[2],p2[2],p3[2],p4[2]],phi)

          maxerr = maximum([abs(x-point[1]),abs(y-point[2])])
          if maxerr > 1e-2
            error("(s,t) error is too large, = ",maxerr)
          end
          #println("good transform")
          Uinterp = shapeEval(UA[c],phi)
          foundInside = true
        # check whether point is on side of element
        else
          if onSegment(p1,point,p2)
            #println("on Segment!")
          elseif onSegment(p2,point,p3)
            #println("on Segment!")
          elseif onSegment(p3,point,p4)
            #println("on Segment!")
          elseif onSegment(p4,point,p1)
            #println("on Segment!")
          end     
        end
      end
    end
  end

  if !(foundInside || foundVal)
    errorstr = string("point (",point[1],",",point[2],") not inside mesh")
    error(errorstr)
  end

  return Uinterp
end
# point transform for mesh
function pointTransform(mesh::T,UA,point) where T<:AbstractMesh
  # transform meshA to xy and cm
  xyA = NodesToArray(mesh.xy)
  cmA = ElementsToArray(mesh.cm)

  Uinterp = pointTransform(xyA,cmA,UA,point)
  return Uinterp
end

"""
  meshTransform(xyA,cmA,UA,xyB) -> UB

Transforms finite element solutions from domain A to domain B, returns solutions on UB

Can handle 4 point (linear) and 9 point (quadratic) quadrilateral elements
"""
function meshTransform(xyA,cmA,UA,xyB)
  nptsB = size(xyB,1)
  UB    = zeros(nptsB)

  for i=1:nptsB
    point = [xyB[i,1],xyB[i,2]]
    UB[i] = pointTransform(xyA,cmA,UA,point)
  end

  return UB
end

# transform mesh using mesh objects and points
# xyB is in Nx2 array where xyB[:,1] = x-vals, xyB[:,2] = y-vals
function meshTransform(meshA::T,UA,xyB) where T<:AbstractMesh
  # transform meshA to xy and cm
  xyA = NodesToArray(meshA.xy)
  cmA = ElementsToArray(meshA.cm)

  # call meshTransform
  UB = meshTransform(xyA,cmA,UA,xyB)

  # return interpolated solutions
  return UB
end

# transform mesh using mesh objects
function meshTransform(meshA::T1,UA,meshB::T2) where 
                  {T1<:AbstractMesh, T2<:AbstractMesh}
  # transform meshA to xy and cm
  xyA = NodesToArray(meshA.xy)
  cmA = ElementsToArray(meshA.cm)

  # transform meshB to xy
  xyB = NodesToArray(meshB.xy)

  # call meshTransform
  UB = meshTransform(xyA,cmA,UA,xyB)

  # return interpolated solutions
  return UB
end

"""
  car2pol(x,y) -> r,theta

Transforms cartesian coordinates into polar
"""
function car2pol(x,y)
  r = sqrt(x^2 + y^2)
  theta = atan2(y,x)

  return r,theta
end


"""
  pol2car(r,theta) -> x,y

Transforms polar coordinates into cartesian
"""
function pol2car(r,theta)
    x = r*cos(theta)
    y = r*sin(theta)

  return x,y
end
