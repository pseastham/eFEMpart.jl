# interface to C++ library "figtree" from
#   Vlad I. Morariu, Balaji Vasan Srinivasan, Vikas C. Raykar, 
#   Ramani Duraiswami, and Larry S. Davis. Automatic online tuning for 
#   fast Gaussian summation. Advances in Neural Information Processing 
#   Systems (NIPS), 2008. 
# and website that hosts code is 
#   http://www.umiacs.umd.edu/~morariu/figtree/ 

using LinearAlgebra, Plots
gr(size=(900,300))

#include("../../src/part/fgt.jl")  # loads in functions to be tests

function fgtTest()
    d = 2                           # dimension of data
    M = 1_000                      # number of targets 
    N = 2
    #N = 1_000                      # number of sources
    h = 0.8                         # bandwidth
    ε = 1e-2                        # max tolerance
    x = [-0.5,0.5]
    #x = rand(N*d)                   # source vectors, reshape'd to 1D array
    y = rand(M*d)                   # target vectors, reshape'd to 1D array
    W = 1                           # number of sources being evaluated 
    q=[0.5,1.0]
    #q = rand(N*W)                   # weights for sources
    v1 = rand(M)                     # output values array
    v2 = rand(M)                    # output 2
    #v3 = rand(M)                    # output 3

    temp = norm(v1-v2)

    fgt!(v1,d,M,N,h,ε,x,y,q,W)
    dgt!(v1,d,M,N,h,ε,x,y,q,W)
    #mydgt!(v3,d,M,N,h,ε,x,y,q,W)

    println("fast (C++) transform: ")
    @time fgt!(v1,d,M,N,h,ε,x,y,q,W)
    println("slow (C++) transform: ")
    @time dgt!(v2,d,M,N,h,ε,x,y,q,W)
    #println("slow (Julia) transform: ")
    #@time mydgt!(v3,d,M,N,h,ε,x,y,q,W)
    println()

    println("(check that two algorithms give same result)")
    println("initial norm: ",temp)
    println("later norm:   ",norm(v1-v2))
end

function fgt1Dtest()
    d = 1                           # dimension of data
    M = 10_000                      # number of targets 
    #N = 2
    N = 10_000                      # number of sources
    h = 0.2                         # bandwidth
    ε = 1e-2                        # max tolerance
    #x = [-0.5,0.5]
    x = rand(N*d)                   # source vectors, reshape'd to 1D array
    y = collect(range(-2,stop=2,length=M))                   # target vectors, reshape'd to 1D array
    W = 1                           # number of sources being evaluated 
    #q=[0.5,1.0]
    q = rand(N*W)                   # weights for sources
    v1 = rand(M)                     # output values array
    v2 = rand(M)                    # output 2
    v3 = rand(M)                    # output 3

    @time fgt!(v1,d,M,N,h,ε,x,y,q,W)
    @time dgt!(v2,d,M,N,h,ε,x,y,q,W)
    @time mydgt!(v3,d,M,N,h,ε,x,y,q,W)

    p1 = plot(y,v1,label="fgt (c++)")
    p2 = plot(y,v2,label="dgt (c++)")
    p3 = plot(y,v3,label="mydgt (julia)")

    l = @layout [a a a]
    p = plot(p1,p2,p3,layout = l)
    gui(p)
end