function plot_state( prob, state, t;
    var=:omega, xlims=:auto, ylims=:auto, clims=:auto, clevs=30)

    if (xlims==:auto) & (ylims==:auto)
        #Would be nice to add vel mag
        if var == :omega
            plot_ω(state, prob.model.grid, lev=1, clims=clims, xlims=xlims,
                ylims=ylims, clevs=clevs )
            plot_bodies(state.xb)
        elseif var == :vel
            plt = plot(layout=(2,1))
            plot_u(state, prob.model.grid, lev=1, clims=clims, xlims=xlims,
                ylims=ylims, plt=plt, clevs=clevs )
            plot_bodies(prob.model.bodies, plt=plt )
        elseif var == :psi
            plot_ψ(state, prob.model.grid, lev=glev, clims=clims, xlims=xlims,
                ylims=ylims, clevs=clevs )
            plot_bodies(prob.model.bodies)
        end
    else
        if var == :vel
            plt = plot(layout=(2,1))
        end

        for glev = prob.model.grid.mg : -1 : 1

            #Would be nice to add vel mag
            if var == :omega
                plot_ω(state, prob.model.grid, lev=glev, clims=clims, xlims=xlims,
                    ylims=ylims, clevs=clevs )
                if glev==1
                    plot_bodies(state.xb)
                end
            elseif var == :vel
                plot_u(state, prob.model.grid, lev=glev, clims=clims, xlims=xlims,
                    ylims=ylims, plt=plt, clevs=clevs )
                if glev==1
                    plot_bodies(prob.model.bodies, plt=plt )
                end
            elseif var == :psi
                plot_ψ(state, prob.model.grid, lev=glev, clims=clims, xlims=xlims,
                    ylims=ylims, clevs=clevs )
                if glev==1
                    plot_bodies(prob.model.bodies)
                end
            end

        end
    end


end


function plot_ω(state::IBState, grid::T; lev=1,
    xlims=:auto, ylims=:auto, clims=:auto,
    colorbar=:false, framestyle=:box, clevs=30) where T <:Grid

    h = grid.h
    len = grid.len

    fac = 2.0^(Float64(lev-1))
    δ = h * fac
    xlen = len*fac

    ylen = xlen*(grid.ny/grid.nx)
    offx = fac * len/2.0 - len/2.0 + grid.offx
    offy = fac * (grid.ny*h)/2.0 - (grid.ny*h)/2.0 + grid.offy

    x = range(-offx+δ, xlen-offx-δ, length=grid.nx-1)
    y = range(-offy+δ, ylen-offy-δ, length=grid.ny-1)

    ω = reshape( state.Γ[:, lev], grid.nx-1, grid.ny-1 )' / δ^2.0

    if xlims==:auto
        xlims=(x[1], x[end])
    end

    if ylims==:auto
        ylims=(y[1], y[end])
    end

    if clims ≠ :auto
        ω[ ω .>= clims[2] ] .= clims[2]
        ω[ ω .<= clims[1] ] .= clims[1]
    end

    display(contourf!(x,y, ω, c=cgrad(:seaborn_icefire_gradient),#c=cgrad(:RdBu_11),#c=cgrad(:diverging_gkr_60_10_c40_n256), #
        colorbar=colorbar, lw=0, levels=clevs, aspect_ratio=:equal,
        clims=clims, xlims=xlims, ylims=ylims, framestyle=framestyle))
end


function plot_bodies(xbv; plt=missing)

    if ismissing(plt)
        for j = 1:length(xbv)
            xb = xbv[j]

            #anoying hack: find interior point of body for fill
            interior = sum(xb[:,2])/length(xb[:,2])

            display(plot!(xb[:, 1], xb[:, 2], lw=3, fillrange=interior,
            fillcolor=:gray, legend=false, color=:gray))
        end
    else
        for j = 1:length(bodies)
            xb = xbv[j]
            for jj = 1 : length(plt)
                #anoying hack: find interior point of body for fill
                interior = sum(xb[:,2])/length(xb[:,2])

                display(plot!(plt[jj],xb[:, 1], xb[:, 2], lw=3,
                fillrange=interior, fillcolor=:gray, legend=false, color=:gray))
            end
        end
    end
end


function plot_u(state::IBState, grid::T; lev=1,
    xlims=:auto, ylims=:auto, clims=:auto,
    colorbar=:false, framestyle=:box, plt=missing, clevs=30) where T <:Grid
    nx, ny = grid.nx, grid.ny
    nu = ny * (nx+1); nv = nx * (ny+1);
    # hc = grid.h*2^(lev-1)
    # x = hc*(1:nx-1) .- grid.offx
    # y = hc*(1:ny) .- grid.offy

    h = grid.h
    len = grid.len

    fac = 2.0^(Float64(lev-1))
    δ = h * fac
    xlen = len*fac

    ylen = xlen*(grid.ny/grid.nx)
    offx = fac * len/2.0 - len/2.0 + grid.offx
    offy = fac * (grid.ny*h)/2.0 - (grid.ny*h)/2.0 + grid.offy

    #-- x vels
    x = range(-offx, xlen-offx, length=grid.nx+1)
    y = range(-offy+δ/2.0, ylen-offy-δ/2.0, length=grid.ny)
    qx = state.q[1:nu, lev] .+ state.q0[1:nu, lev]

    #for some reason the top left velocity on grid 1 is erroneous.
    #extract interior velocities....
    #TODO: look into this...
    u = (reshape( qx, length(x), length(y) ) / δ)'
    u = u[2:end-1, 2:end-1]
    x = x[2:end-1]
    y = y[2:end-1]

    if xlims==:auto
        xlims=(x[1], x[end])
    end

    if ylims==:auto
        ylims=(y[1], y[end])
    end

    if ismissing(plt)
        display(
            contourf!(x, y, u, c=cgrad(:seaborn_icefire_gradient),#c=cgrad(:temperaturemap), #c=cgrad(:diverging_gkr_60_10_c40_n256),#c=cgrad(:blackbody),#c=cgrad(:temperaturemap),#,#
            colorbar=:true, lw=0, levels=clevs, aspect_ratio=:equal,
            clims=clims, xlims=xlims, ylims=ylims, framestyle=framestyle,
            legend=:false, widen=false)
            )
    else

        display(
            contourf!(plt[1], x, y, u, c=cgrad(:seaborn_icefire_gradient),#c=cgrad(:temperaturemap), #c=cgrad(:diverging_gkr_60_10_c40_n256),#c=cgrad(:blackbody),#c=cgrad(:temperaturemap),#,#
            colorbar=:true, lw=0, levels=clevs, aspect_ratio=:equal,
            clims=clims, xlims=xlims, ylims=ylims, framestyle=framestyle,
            legend=:false, widen=false)
            )
    end

    #-- y vels
    x = range(-offx+δ/2.0, xlen-offx-δ/2.0, length=grid.nx)
    y = range(-offy, ylen-offy, length=grid.ny+1)
    qy = state.q[nu.+(1:nv), lev] .+ state.q0[nu.+(1:nv), lev]
    v = (reshape( qy, length(x), length(y) ) / δ)'

    # vplot = contourf!(x, y, v, c=cgrad(:seaborn_icefire_gradient),#c=cgrad(:temperaturemap), #c=cgrad(:diverging_gkr_60_10_c40_n256),#c=cgrad(:blackbody),#c=cgrad(:temperaturemap),#,#
    #     colorbar=:true, lw=0, levels=30, aspect_ratio=:equal,
    #     clims=clims, xlims=xlims, ylims=ylims, framestyle=framestyle,
    #     legend=:false)

    if xlims==:auto
        xlims=(x[1], x[end])
    end

    if ylims==:auto
        ylims=(y[1], y[end])
    end

    if ismissing(plt)
        display(
            contourf!(x, y, v, c=cgrad(:seaborn_icefire_gradient),#c=cgrad(:temperaturemap), #c=cgrad(:diverging_gkr_60_10_c40_n256),#c=cgrad(:blackbody),#c=cgrad(:temperaturemap),#,#
            colorbar=:true, lw=0, levels=clevs, aspect_ratio=:equal,
            clims=clims, xlims=xlims, ylims=ylims, framestyle=framestyle,
            legend=:false, widen=false)
            )
    else
        display(
            contourf!(plt[2], x, y, v, c=cgrad(:seaborn_icefire_gradient),#c=cgrad(:temperaturemap), #c=cgrad(:diverging_gkr_60_10_c40_n256),#c=cgrad(:blackbody),#c=cgrad(:temperaturemap),#,#
            colorbar=:true, lw=0, levels=clevs, aspect_ratio=:equal,
            clims=clims, xlims=xlims, ylims=ylims, framestyle=framestyle,
            legend=:false, widen=false)
            )
    end

end


function plot_ψ(state::IBState, grid::T; lev=1,
    xlims=:auto, ylims=:auto, clims=:auto,
    colorbar=:false, framestyle=:box, clevs=30) where T <:Grid

    h = grid.h
    len = grid.len

    fac = 2.0^(Float64(lev-1))
    δ = h * fac
    xlen = len*fac

    ylen = xlen*(grid.ny/grid.nx)
    offx = fac * len/2.0 - len/2.0 + grid.offx
    offy = fac * (grid.ny*h)/2.0 - (grid.ny*h)/2.0 + grid.offy

    x = range(-offx+δ, xlen-offx-δ, length=grid.nx-1)
    y = range(-offy+δ, ylen-offy-δ, length=grid.ny-1)

    ψ = reshape( state.ψ[:, lev], length(x), length(y) )'
    display(contourf!(x, y, ψ, c=cgrad(:seaborn_icefire_gradient),
        colorbar=:true, lw=1, levels=clevs, aspect_ratio=:equal,
        framestyle=framestyle, clims=clims, xlims=xlims, ylims=ylims,
        background_color=:transparent, foreground_color=:transparent))
end

using .IBPM.Quantities

mutable struct ParticlePlot
    framerate::Int64
    particlePos::Any
    particlePosEdgesX::Any
    particlePosEdgesY::Any
    discretization::NTuple{2, Int64}
    distance0X::Float64
    distance0Y::Float64
    xstate::Any
    ystate::Any
    function ParticlePlot(framerate, discretization, grid)
        pp = new()

        nx, ny = grid.nx, grid.ny
        nu = ny * (nx+1); nv = nx * (ny+1);
        # hc = grid.h*2^(lev-1)
        # x = hc*(1:nx-1) .- grid.offx
        # y = hc*(1:ny) .- grid.offy
        lev = 1
        h = grid.h
        len = grid.len

        fac = 2.0^(Float64(lev-1))
        δ = h * fac
        xlen = len*fac

        ylen = xlen*(grid.ny/grid.nx)
        offx = fac * len/2.0 - len/2.0 + grid.offx
        offy = fac * (grid.ny*h)/2.0 - (grid.ny*h)/2.0 + grid.offy

        #-- x vels
        x = range(-offx, xlen-offx, length=grid.nx+1)
        y = range(-offy+δ/2.0, ylen-offy-δ/2.0, length=grid.ny)
        # qx = state.q[1:nu, lev] .+ state.q0[1:nu, lev]

        # #for some reason the top left velocity on grid 1 is erroneous.
        # #extract interior velocities....
        # #TODO: look into this...
        # u = (reshape( qx, length(x), length(y) ) / δ)'
        # u = u[2:end-1, 2:end-1]
        # println(size(u))

        pp.xstate = x[2:end-1]
        pp.ystate = y[2:end-1]

        nParticles = discretization[1] * discretization[2]
        pp.framerate = framerate
        pp.particlePos = Array{Float64}(undef, 2, nParticles)
        xlen = x[end] - x[1]
        ylen = y[end] - y[1]
        iter = 1
        xdistance = xlen / (discretization[1] - 1)
        ydistance = ylen / (discretization[2] - 1)
        pp.distance0X = xdistance
        pp.distance0Y = ydistance
        for i in 0:discretization[1]-1
            xpos = i * xdistance + pp.xstate[1]
            for j in 0:discretization[2]-1
                ypos = j * ydistance + pp.ystate[1]
                pp.particlePos[1,iter] = xpos
                pp.particlePos[2,iter] = ypos
                iter += 1
            end
        end
        pp.discretization = discretization

        pp.particlePosEdgesY = Array{Float64}(undef, 2, 2*discretization[2])
        pp.particlePosEdgesX = Array{Float64}(undef, 2, 2*discretization[1]- 4)
        iter = 1
        for j in 0:discretization[2]-1
            ypos = j * ydistance + pp.ystate[1]
            pp.particlePosEdgesY[1, iter] = pp.xstate[1]
            pp.particlePosEdgesY[2, iter] = ypos
            iter += 1
            pp.particlePosEdgesY[1, iter] = pp.xstate[end]
            pp.particlePosEdgesY[2, iter] = ypos
            iter += 1
        end
        iter = 1
        for i in 1:discretization[1]-2
            xpos = i * xdistance + pp.xstate[1]
            pp.particlePosEdgesX[1, iter] = xpos
            pp.particlePosEdgesX[2, iter] = pp.ystate[1]
            iter += 1
            pp.particlePosEdgesX[1, iter] = xpos
            pp.particlePosEdgesX[2, iter] = pp.ystate[end]
            iter += 1
        end
        
        return pp
    end
end

function particleUpdate!(dt, particlePlot::ParticlePlot, state)
    x_vel = x_velocity(state)[:, :, 1][2:end, :]
    y_vel = y_velocity(state)[:, :, 1][:, 2:end]
    # println(size(x_vel))
    # println(size(y_vel))
    # println(x_vel[1])
    # println(x_vel[200])
    # println(x_vel[300])
    # println()
    
    # println(size(particlePlot.particlePos))
    for (i, particle) in enumerate(eachcol(particlePlot.particlePos))
        # println(particle)
        xpos = particle[1]
        ypos = particle[2]
        xind = findfirst(x -> x >= xpos, particlePlot.xstate)
        yind = findfirst(y -> y >= ypos, particlePlot.ystate)
        # if isnothing(xind)
        #     xind = 1
        # end
        # if isnothing(yind)
        #     yind = 1
        # end
        if isnothing(xind) || isnothing(yind)
            # randnum = rand(1:size(particlePlot.particlePosEdges)[2])
            # particlePlot.particlePos[1, i] = particlePlot.particlePosEdges[1, randnum]
            # particlePlot.particlePos[2, i] = particlePlot.particlePosEdges[2, randnum]
            # first iterate through all edge particle position
            # then iterate through all current particles
            # calculate min distance between edge particle and any current particle
            # if min distance is larger than threshold, put new particle at that edge particle position
            # else find a particle position that does
            # TODO: what to do if no potential particle position?
            found = false
            for (j, particleEdge) in enumerate(eachcol(particlePlot.particlePosEdgesY))
                xposEdge = particleEdge[1]
                yposEdge = particleEdge[2]
                mindist = 10 #Arbitray big number cause lazy
                mindistx = 0
                mindisty = 0
                for particle2 in eachcol(particlePlot.particlePos)
                    xpos2 = particle2[1]
                    ypos2 = particle2[2]
                    xdist = abs(xposEdge - xpos2)
                    ydist = abs(yposEdge - ypos2)
                    dist = sqrt(xdist^2 + ydist^2)
                    if dist < mindist
                        mindist = dist
                        mindistx = xdist
                        mindisty = ydist
                    end
                end
                if mindistx >= particlePlot.distance0X
                    particlePlot.particlePos[1, i] = particlePlot.particlePosEdgesY[1, j]
                    particlePlot.particlePos[2, i] = particlePlot.particlePosEdgesY[2, j]
                    break
                end
                if mindisty >= particlePlot.distance0Y
                    particlePlot.particlePos[1, i] = particlePlot.particlePosEdgesY[1, j]
                    particlePlot.particlePos[2, i] = particlePlot.particlePosEdgesY[2, j]
                    break
                end
            end
            if !found
                for (j, particleEdge) in enumerate(eachcol(particlePlot.particlePosEdgesX))
                    xposEdge = particleEdge[1]
                    yposEdge = particleEdge[2]
                    mindist = 10 #Arbitray big number cause lazy
                    mindistx = 0
                    mindisty = 0
                    for particle2 in eachcol(particlePlot.particlePos)
                        xpos2 = particle2[1]
                        ypos2 = particle2[2]
                        xdist = abs(xposEdge - xpos2)
                        ydist = abs(yposEdge - ypos2)
                        dist = sqrt(xdist^2 + ydist^2)
                        if dist < mindist
                            mindist = dist
                            mindistx = xdist
                            mindisty = ydist
                        end
                    end
                    if mindistx >= particlePlot.distance0X
                        particlePlot.particlePos[1, i] = particlePlot.particlePosEdgesX[1, j]
                        particlePlot.particlePos[2, i] = particlePlot.particlePosEdgesX[2, j]
                        break
                    end
                    if mindisty >= particlePlot.distance0Y
                        particlePlot.particlePos[1, i] = particlePlot.particlePosEdgesX[1, j]
                        particlePlot.particlePos[2, i] = particlePlot.particlePosEdgesX[2, j]
                        break
                    end
                end
            end
            continue
        end
        xvel = x_vel[xind, yind]
        yvel = y_vel[xind, yind]
        newx = xpos + dt*xvel
        newy = ypos + dt*yvel
        particlePlot.particlePos[1, i] = newx
        particlePlot.particlePos[2, i] = newy
    end
end

# function newParticle(particlePlot::ParticlePlot)
#     particlePlot.particlePos
# end

function plotParticles(dt, particlePlot::ParticlePlot, state, grid)
    scatter!(particlePlot.particlePos[1,:], particlePlot.particlePos[2,:])
    dt = 1/particlePlot.framerate
    particleUpdate!(dt, particlePlot, state)
end