function update_particles_noFEMcl!(mesh,u,v,pList,wList,param,data,cl)
    pUarr,pVarr = compute_particle_velocity_noFEMcl(mesh,u,v,pList,wList,param,data,cl)
    for ti=1:data.n_particles
        pList[ti].pos.x += pUarr[ti]*param.Δt
        pList[ti].pos.y += pVarr[ti]*param.Δt
    end

    nothing
end

# does not use cell lists
function compute_particle_velocity_noFEMcl(mesh,u,v,pList,wList,param,data,cl)
    StokesParticles.compute_gravity_force!(data.gfX,data.gfY,data.n_particles,param.G)
    compute_seepage_force_nocl!(data.sfX,data.sfY,mesh,u,v,pList,data)
    StokesParticles.compute_cohesion_force_cl!(data.cfX,data.cfY,cl,pList,param.sC,param.ϵ)
    StokesParticles.compute_adhesion_force!(data.afX,data.afY,pList,wList,data.pointOnWall,param.sA,param.ϵ)

    pUarr = data.gfX + data.cfX + data.afX + data.sfX
    pVarr = data.gfY + data.cfY + data.afY + data.sfY

    return pUarr,pVarr
end