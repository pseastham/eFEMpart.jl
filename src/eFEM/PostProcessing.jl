# Code to be loaded into eFEMpart

function computeStress(mesh::FluidMesh,fluidsol::FluidSolution,param::P) where P<:AbstractConstantParameter
    Nnodes = length(mesh.xy)
    Nel    = length(mesh.cm)
    σxx = zeros(Float64,Nnodes)
    σxy = zeros(Float64,Nnodes)
    σyy = zeros(Float64,Nnodes)

    xN = zeros(Float64,9)
    yN = zeros(Float64,9)

    nodeCounter = zeros(Int,Nnodes)

    # compute pressure at Q2 mesh
    xyArr = NodesToArray(mesh.xy)
    xypArr = NodesToArray(mesh.xyp)
    cmpArr = ElementsToArray(mesh.cmp)
    p = meshTransform(xypArr,cmpArr,fluidsol.p,xyArr)

    s = [-1.0,1.0,1.0,-1.0,0.0,1.0,0.0,-1.0,0.0]
    t = [-1.0,-1.0,1.0,1.0,-1.0,0.0,1.0,0.0,0.0]

    for el=1:Nel
        c = mesh.cm[el].NodeList
        xN = [mesh.xy[i].x for i in c]
        yN = [mesh.xy[i].y for i in c]
        
        for ti=1:9
            phi,dphidx,dphidy,jac = derivShape2D(s[ti],t[ti],xN,yN,2)

            # interpolate values at gauss points
            dudxg = shapeEval(fluidsol.u[c],dphidx)/jac
            dudyg = shapeEval(fluidsol.u[c],dphidy)/jac
            dvdxg = shapeEval(fluidsol.v[c],dphidx)/jac
            dvdyg = shapeEval(fluidsol.v[c],dphidy)/jac

            σxx[c[ti]] += -p[c[ti]] + 2*param.μ*dudxg
            σxy[c[ti]] += param.μ*(dudyg + dvdxg)
            σyy[c[ti]] += -p[c[ti]] + 2*param.μ*dvdyg

            nodeCounter[c[ti]] += 1
        end
    end
    
    # divide through by number of times contributed (finds average at each node from surrounding elements)
    for ti=1:Nnodes
        σxx[ti] = σxx[ti]/nodeCounter[ti]
        σxy[ti] = σxy[ti]/nodeCounter[ti]
        σyy[ti] = σyy[ti]/nodeCounter[ti]
    end

    return σxx, σxy, σyy
end

function computeStress(mesh::FluidMesh,fluidsol::FluidSolution,param::P) where P<:AbstractVariableParameter
    Nnodes = length(mesh.xy)
    Nel    = length(mesh.cm)
    σxx = zeros(Float64,Nnodes)
    σxy = zeros(Float64,Nnodes)
    σyy = zeros(Float64,Nnodes)

    xN = zeros(Float64,9)
    yN = zeros(Float64,9)

    nodeCounter = zeros(Int,Nnodes)

    # compute pressure at Q2 mesh
    xyArr = NodesToArray(mesh.xy)
    xypArr = NodesToArray(mesh.xyp)
    cmpArr = ElementsToArray(mesh.cmp)
    p = meshTransform(xypArr,cmpArr,fluidsol.p,xyArr)

    s = [-1.0,1.0,1.0,-1.0,0.0,1.0,0.0,-1.0,0.0]
    t = [-1.0,-1.0,1.0,1.0,-1.0,0.0,1.0,0.0,0.0]

    for el=1:Nel
        c = mesh.cm[el].NodeList
        xN = [mesh.xy[i].x for i in c]
        yN = [mesh.xy[i].y for i in c]
        
        for ti=1:9
            phi,dphidx,dphidy,jac = derivShape2D(s[ti],t[ti],xN,yN,2)

            # interpolate values at gauss points
            dudxg = shapeEval(fluidsol.u[c],dphidx)/jac
            dudyg = shapeEval(fluidsol.u[c],dphidy)/jac
            dvdxg = shapeEval(fluidsol.v[c],dphidx)/jac
            dvdyg = shapeEval(fluidsol.v[c],dphidy)/jac
            μg    = shapeEval(param.μ[c],phi)

            σxx[c[ti]] += -p[c[ti]] + 2*μg*dudxg
            σxy[c[ti]] += μg*(dudyg + dvdxg)
            σyy[c[ti]] += -p[c[ti]] + 2*μg*dvdyg

            nodeCounter[c[ti]] += 1
        end
    end
    
    # divide through by number of times contributed (finds average at each node from surrounding elements)
    for ti=1:Nnodes
        σxx[ti] = σxx[ti]/nodeCounter[ti]
        σxy[ti] = σxy[ti]/nodeCounter[ti]
        σyy[ti] = σyy[ti]/nodeCounter[ti]
    end

    return σxx, σxy, σyy
end

function computeDerivative(mesh::M,u::Vector{Float64},var::Symbol) where M<:AbstractMesh
    Nnodes = length(mesh.xy)
    Nel    = length(mesh.cm)
    dudi = zeros(Float64,Nnodes)

    xN = zeros(Float64,9)
    yN = zeros(Float64,9)

    nodeCounter = zeros(Int,Nnodes)

    s = [-1.0,1.0,1.0,-1.0,0.0,1.0,0.0,-1.0,0.0]
    t = [-1.0,-1.0,1.0,1.0,-1.0,0.0,1.0,0.0,0.0]

    for el=1:Nel
        c = mesh.cm[el].NodeList
        xN = [mesh.xy[i].x for i in c]
        yN = [mesh.xy[i].y for i in c]
        
        for ti=1:9
            phi,dphidx,dphidy,jac = derivShape2D(s[ti],t[ti],xN,yN,2)

            # interpolate values at gauss points
            if var == :dx
                dudig = shapeEval(u[c],dphidx)/jac
            elseif var == :dy
                dudig = shapeEval(u[c],dphidy)/jac
            else
                error("variable ain't right")
            end

            dudi[c[ti]] += dudig
            nodeCounter[c[ti]] += 1
        end
    end
    
    # divide through by number of times contributed (finds average at each node from surrounding elements)
    for ti=1:Nnodes
        dudi[ti] = dudi[ti]/nodeCounter[ti]
    end

    return dudi
end

function computeExtension(mesh::FluidMesh,fluidsol::FluidSolution,param)
    nNodes = length(mesh.xy)

    σxx,σxy,σyy = computeStress(mesh,fluidsol,param)

    λField = zeros(Float64,nNodes)
    evecXField = zeros(Float64,nNodes)
    evecYField = zeros(Float64,nNodes)

    for ti=1:nNodes
        # construct stress tensor
        T = [σxx[ti] σxy[ti]; σxy[ti] σyy[ti]]

        # find eigenvalues and eigenvectors of stress tensor T
        F = eigfact(T)

        # find maximum eigenvalue
        λField[ti]     = abs(maximum(abs.(F[:values])))
    end

    return λField
end

function computeCompressibility(mesh::M,fSol::FluidSolution) where M<:AbstractMesh
    u = fSol.u; v=fSol.v

    dudx = computeDerivative(mesh,fSol.u,:dx)
    dvdy = computeDerivative(mesh,fSol.v,:dy)

    compress = dudx .+ dvdy

    return compress
end