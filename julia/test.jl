module NBody

using WriteVTK

# Profile.init(delay=0.01)

# Constants
const dt = 0.01
const tstop = 10*dt
const particles = 100000
const cutoff = 0.1/∛particles

type Body
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    mass::Float64
    name::Int64
end

immutable Pair
    m::Int32
    n::Int32
end

immutable Bin
    nx::Int
    ny::Int
    nz::Int
end

function get_pairs(bodies)
    ngrid = floor(Integer, 0.2/cutoff)
    nx = ngrid
    ny = ngrid
    nz = ngrid
    particles_in_bin = [ Array(Int,0) for x=1:nx, y=1:ny, z=1:nz]
    # @printf("particles_in_bin = %s\n", particles_in_bin)
    bin_of_particle = Array(Any, length(bodies))
    @printf("BINNING")
    @time for n in 1:length(bodies)
        try
            body_nx = max(1,min(nx, 1+floor(Int64, bodies[n].x * nx)))
            body_ny = max(1,min(ny, 1+floor(Int64, bodies[n].y * ny)))
            body_nz = max(1,min(nz, 1+floor(Int64, bodies[n].z * nz)))
            bin_of_particle[n] = Bin(body_nx, body_ny, body_nz)
            push!(particles_in_bin[body_nx, body_ny, body_nz],n)
            if 1 < body_nx
                push!(particles_in_bin[body_nx-1, body_ny, body_nz],n)
            end
            if body_nx < nx
                push!(particles_in_bin[body_nx+1, body_ny, body_nz],n)
            end
            if 1 < body_ny
                push!(particles_in_bin[body_nx, body_ny-1, body_nz],n)
            end
            if body_ny < ny
                push!(particles_in_bin[body_nx, body_ny+1, body_nz],n)
            end
            if 1 < body_nz
                push!(particles_in_bin[body_nx, body_ny, body_nz-1],n)
            end
            if body_nz < nz
                push!(particles_in_bin[body_nx, body_ny, body_nz+1],n)
            end
        catch
            @printf("bodies[%d] - %s\n", n, bodies[n])
            @printf("%d %d %d\n", nx, ny, nz)
        end
    end

    @printf("PAIRING")
    pairs = @time @parallel (append!) for n in 1:length(bodies)
        pairs_for_body(n, bin_of_particle, particles_in_bin, bodies)
    end
    pairs
end

function pairs_for_body(n, bin_of_particle, particles_in_bin, bodies)
        pares = Array(Pair, 0)
        bin = bin_of_particle[n]
        for m in particles_in_bin[bin.nx, bin.ny, bin.nz]
            if m < n
                dx = bodies[m].x - bodies[n].x
                dy = bodies[m].y - bodies[n].y
                dz = bodies[m].z - bodies[n].z
                dsq = dx^2 + dy^2 + dz^2
                if dsq < cutoff
                    push!(pares,Pair(m,n))
                end
            end
        end
        pares
    end

function advance!(pairs, bodies, dt)

    for pair in pairs
        i, j = pair.m, pair.n

        dx = bodies[j].x - bodies[i].x
        dy = bodies[j].y - bodies[i].y
        dz = bodies[j].z - bodies[i].z
        dsq = dx^2 + dy^2 + dz^2
        distance = sqrt(dsq)
        mag = dt / (dsq * distance)

        bodies[i].vx -= dx * bodies[j].mass * mag
        bodies[i].vy -= dy * bodies[j].mass * mag
        bodies[i].vz -= dz * bodies[j].mass * mag

        bodies[j].vx += dx * bodies[i].mass * mag
        bodies[j].vy += dy * bodies[i].mass * mag
        bodies[j].vz += dz * bodies[i].mass * mag
    end

    for b in bodies
        b.x += dt * b.vx
        b.y += dt * b.vy
        b.z += dt * b.vz
        if b.x < 0 || 1 < b.x
            b.vx = -b.vx
        end
        if b.y < 0 || 1 < b.y
            b.vy = -b.vy
        end
        if b.z < 0 || 1 < b.z
            b.vz = -b.vz
        end
    end

    for b in bodies
        b.x += dt * b.vx
        b.y += dt * b.vy
        b.z += dt * b.vz
    end
    bodies

end

function RandomBody(name)
    Body(rand(),
         rand(),
         rand(),
         rand(),
         rand(),
         rand(),
         0.00001,
         name,
         )
end

function nbody(tstop::Float64)

    @printf("INITIALIZING")
    @time bodies = [ RandomBody(name) for name in 1:particles]

    x = collect(1:10)
    y = collect(1:10)
    z = collect(1:10)

    pvd = paraview_collection("my_pvd_file")

    no_cells = convert(Array{WriteVTK.MeshCell,1}, [])

    # pairs = Tuple{Int64,Int64}[]
    # for i=1:length(bodies)
    #     for j=i+1:length(bodies)
    #         push!(pairs, (i,j))
    #     end
    # end

    TIME = 0
    i = 0
    while ( TIME < tstop)
        walltime = time()
        @printf(" at time %s\n", TIME)
        pairs = get_pairs(bodies)
        @printf("ADVANCING")
        @time bodies = advance!(pairs, bodies, dt)

        @printf("outputting")
        @time begin
            vtkfile = vtk_grid("my_vtk_file_$(i)", [ [b.x, b.y, b.z][ii] for b in bodies, ii in 1:3], no_cells)
            vtk_save(vtkfile)
            collection_add_timestep(pvd, vtkfile, i*1.0)
        end
        TIME += dt
        i += 1
        @printf(" ELAPSED:::  %s", (time()-walltime))
    end
    vtk_save(pvd)
end

nbody(tstop)
# Profile.print()

end # module
