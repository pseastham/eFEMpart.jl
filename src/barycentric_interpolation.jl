# Computes barycentric weights, can be used for value interpolation
# from finite element mesh

#   defining quadratic function
quadratic(a, sqr_term, b) = (-b + sqr_term) / 2a

#   finding positive and negative components of quadratic function
function quadratic2(a::T, b::T, c::T) where T<:Real
    sqr_term = sqrt(b^2-4a*c)
    r1 = quadratic(a, sqr_term, b)
    r2 = quadratic(a, -sqr_term, b)

    r1, r2  # the two roots
end

function computeBaryWeights!(w::Vector{T},x::Vector{T}, y::Vector{T}, p::Vector{T}) where T<:Real
    #   Computing the Barycentric weights
    #

    #   INPUTS:
    #   x  = (x1, x2, x3, x4) [vector of x values]
    #   y  = (y1, y2, y3, y4) [vector of y values]
    #   z  = (f(x1, y1), f(x2, y2), f(x3, y3), f(x4, y4)) [vector of z values]
    #   xi = x-value of point we want to interpolate
    #   yi = y-value of point we want to interpolate

    #   OUTPUT (IN PLACE):
    #   w  = (w1, w2, w3, w4) [vector of barycentric weights for sanity check]

    # dummy variable
    bad = -1000.0

    a = [x[1] - p[1], y[1]-p[2], 0.0]
    b = [x[2] - x[1], y[2]-y[1], 0.0]
    c = [x[4] - x[1], y[4]-y[1], 0.0]
    d = [x[1] - x[2] - x[4] + x[3], y[1] - y[2] - y[4] + y[3], 0.0]

    # FINDING μ terms
    A = cross(c,d)
    B = cross(c,b)+cross(a,d)
    C = cross(a,b)

    if (abs(A[3]) < 1.0e-10) # if A=0, we have a line
        u1 = -C[3]/B[3]
        u2 = u1
    else
        sqr_term = sqrt(B[3]^2-4*A[3]*C[3])
        if (sqr_term > zero(T)) # only looking at real solutions
            u1, u2 = quadratic2(A[3], B[3], C[3])
        else #complex solutions don't matter
            u1 = bad
            u2 = u1
        end
    end

    μ = bad
    if (u1 >= zero(T) && u1 <= one(T))
        μ = u1
    end
    if (u2 >= zero(T) && u2 <= one(T))
        μ = u2
    end

    # FINDING λ terms
    A = cross(b,d)
    B = cross(b,c)+cross(a,d)
    C = cross(a,c)

    if (abs(A[3]) < 1.0e-10) # if A=0, we have a line
        u3 = -C[3]/B[3]
        u4 = u3
    else
        sqr_term = sqrt(B[3]^2-4*A[3]*C[3])
        if (sqr_term > zero(T)) # only looking at real solutions
            u3, u4 = quadratic2(A[3], B[3], C[3])
        else #complex solutions don't matter
            u3 = bad
            u4 = u3
        end
    end

    λ = bad
    if (u3 >= zero(T) && u3 <= one(T))
        λ = u3
    end
    if (u4 >= zero(T) && u4 <= one(T))
        λ = u4
    end

    # barycentric weights
    w[1] = (1-μ)*(1-λ)
    w[2] = λ*(1-μ)
    w[3] = λ*μ
    w[4] = (1-λ)*μ

    return w #returning vector of weights
end
# w, a, b, c input just to avoid re-allocation
function computeBaryWeights!(w::Vector{T},polygon::Vector{Point2D{T}}, p::Point2D{T},
                                a::Vector{T},b::Vector{T},c::Vector{T},d::Vector{T}) where T<:Real
    #   Computing the Barycentric weights
    #

    #   INPUTS:
    #   x  = (x1, x2, x3, x4) [vector of x values]
    #   y  = (y1, y2, y3, y4) [vector of y values]
    #   z  = (f(x1, y1), f(x2, y2), f(x3, y3), f(x4, y4)) [vector of z values]
    #   xi = x-value of point we want to interpolate
    #   yi = y-value of point we want to interpolate

    #   OUTPUT (IN PLACE):
    #   w  = (w1, w2, w3, w4) [vector of barycentric weights for sanity check]

    # dummy variable
    bad = -1000.0

    a[1] = polygon[1].x - p.x
    a[2] = polygon[1].y-p.y
    a[3] = 0.0
    b[1] = polygon[2].x - polygon[1].x
    b[2] = polygon[2].y-polygon[1].y
    b[3] = 0.0
    c[1] = polygon[4].x - polygon[1].x
    c[2] = polygon[4].y-polygon[1].y
    c[3] = 0.0
    d[1] = polygon[1].x - polygon[2].x - polygon[4].x + polygon[3].x
    d[2] = polygon[1].y - polygon[2].y - polygon[4].y + polygon[3].y
    d[3] = 0.0

    # FINDING μ terms
    A1 = c[2]*d[3]-c[3]*d[2]; B1 = c[2]*b[3]-c[3]*b[2] + a[2]*d[3]-a[3]*d[2]; C1 = a[2]*b[3]-a[3]*b[2]
    A2 = c[3]*d[1]-c[1]*d[3]; B2 = c[3]*b[1]-c[1]*b[3] + a[3]*d[1]-a[1]*d[3]; C2 = a[3]*b[1]-a[1]*b[3]
    A3 = c[1]*d[2]-c[2]*d[1]; B3 = c[1]*b[2]-c[2]*b[1] + a[1]*d[2]-a[2]*d[1]; C3 = a[1]*b[2]-a[2]*b[1]

    #A = cross(c,d)
    #B = cross(c,b) + cross(a,d)
    #C = cross(a,b)

    #if (abs(A[3]) < 1.0e-10) # if A=0, we have a line
    if (abs(A3) < 1.0e-10) # if A=0, we have a line        
        #u1 = -C[3]/B[3]
        u1 = -C3/B3
        u2 = u1
    else
        #sqr_term = sqrt(B[3]^2-4*A[3]*C[3])
        sqr_term = sqrt(B3^2-4*A3*C3)
        if (sqr_term > zero(T)) # only looking at real solutions
            #u1, u2 = quadratic2(A[3], B[3], C[3])
            u1, u2 = quadratic2(A3, B3, C3)
        else #complex solutions don't matter
            u1 = bad
            u2 = u1
        end
    end

    μ = bad
    if (u1 >= zero(T) && u1 <= one(T))
        μ = u1
    end
    if (u2 >= zero(T) && u2 <= one(T))
        μ = u2
    end

    # FINDING λ terms
    A1 = b[2]*d[3]-b[3]*d[2]; B1 = b[2]*c[3]-b[3]*c[2] + a[2]*d[3]-a[3]*d[2]; C1 = a[2]*c[3]-a[3]*c[2]
    A2 = b[3]*d[1]-b[1]*d[3]; B2 = b[3]*c[1]-b[1]*c[3] + a[3]*d[1]-a[1]*d[3]; C2 = a[3]*c[1]-a[1]*c[3]
    A3 = b[1]*d[2]-b[2]*d[1]; B3 = b[1]*c[2]-b[2]*c[1] + a[1]*d[2]-a[2]*d[1]; C3 = a[1]*c[2]-a[2]*c[1]

    #A = cross(b,d)
    #B = cross(b,c) + cross(a,d)
    #C = cross(a,c)

    #if (abs(A[3]) < 1.0e-10) # if A=0, we have a line
    if (abs(A3) < 1.0e-10) # if A=0, we have a line
        u3 = -C3/B3
        u4 = u3
    else
        #sqr_term = sqrt(B[3]^2-4*A[3]*C[3])
        sqr_term = sqrt(B3^2-4*A3*C3)
        if (sqr_term > zero(T)) # only looking at real solutions
            #u3, u4 = quadratic2(A[3], B[3], C[3])
            u3, u4 = quadratic2(A3, B3, C3)
        else #complex solutions don't matter
            u3 = bad
            u4 = u3
        end
    end

    λ = bad
    if (u3 >= zero(T) && u3 <= one(T))
        λ = u3
    end
    if (u4 >= zero(T) && u4 <= one(T))
        λ = u4
    end

    # barycentric weights
    w[1] = (1-μ)*(1-λ)
    w[2] = λ*(1-μ)
    w[3] = λ*μ
    w[4] = (1-λ)*μ

    return w #returning vector of weights
end
function computeBaryWeights(x::Vector{T}, y::Vector{T}, p::Vector{T}) where T<:Real
    w = zeros(T,4)
    computeBaryWeights!(w::Vector{T},x::Vector{T}, y::Vector{T}, p::Vector{T})

    return w
end