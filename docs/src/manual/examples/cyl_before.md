# Example: 2D Uniform flow around a stationary cylinder (Old UI)
 
In this example we will solve the problem of a uniform freestream of 1 m/s flowing around a 2 dimentional cylinder of radius 0.5 m over the time interval from 0 to 10 seconds.
 
The general workflow is to define the problem, choose desired state information to save, solve the problem, then analyze the solution. The full code for solving this example problem is:

```julia
include("../src/ibpm.jl")
using .ibpm
using FileIO

boundary = (-1.0, 3.0, -2.0, 2.0)
Re = 100.0
type = :cylinder
lengthscale = 0.5
motion = :static
center = [0.0; 0.0]
body = [(type=type, lengthscale=lengthscale, motion=motion, center=center)]

freestream = (Ux=t->t^0.0,)
T = 1.0
Δx = 0.02
Δt = 0.004
mg = 5

svfc = [ (t,state)->state.CL; (t,state)->state  ]
svti = [ 0.0, T/30.0 ]
svty = [ Vector{Float64}; Any ]

save_info = ( save_fcns = svfc, save_types = svty, save_times = svti )

prob, data, runtime = IBPM_advance( Re, boundary, body, freestream,
    Δx=Δx, Δt=Δt, T=T, save_info=save_info )

using Plots
    
ibpm.plot_state( prob, data[2].val[end], data[2].t[end], var=:omega,
    xlims=(-4.0, 10.0), ylims=(-3.0, 3.0), clims=(-5.0, 5.0), clevs=40)

anim = @animate for j = 1 : length(data[2].t)
            ibpm.plot_state( prob, data[2].val[j], data[2].t[j], var=:omega,
            xlims=(-4.0, 10.0), ylims=(-3.0, 3.0), clims=(-5.0, 5.0), clevs=40 )
        end
gif(anim, "cyl_100.gif", fps=10)

```

## Step 1: Defining the Problem Parameters

The problem is defined by:
- Domain
   - Left boundary
   - Right boundary
   - Upper boundary
   - Lower boundary
- Reynolds Number
- Rigid Bodies
   - Position
   - Size
   - Static or motion (can be function of time)
- Freestream
   - X axis velocity (can be function of time)
   - Y axis velocity (can be funciton of time)
   - Inclination
- Total simulation time
- Grid size
- Time step size
- Number of sub domains

We can initialize these parameters for later use in the solver:

```julia
boundary = (-1.0, 3.0, -2.0, 2.0)
Re = 100.0
type = :cylinder
lengthscale = 0.5
motion = :static
center = [0.0; 0.0]
body = [(type=type, lengthscale=lengthscale, motion=motion, center=center)]
```

There are a number of optional arguments that can be specified when inputted into the solver:

```julia
freestream = (Ux=t->t^0.0,)
T = 1.0 
Δx = 0.02   
Δt = 0.004
mg = 5
```

Finally, the desired saved calculated quantities can be optionally defined:

```julia
svfc = [ (t,state)->state.CL; (t,state)->state  ]
svti = [ 0.0, T/30.0 ]
svty = [ Vector{Float64}; Any ]

save_info = ( save_fcns = svfc, save_types = svty, save_times = svti )
```

## Step 2: Solve the problem

Once all required arguments are defined, they can be used to call the solver method IBPM\_advance() See [IBPM\_advance](../../files/ibpm.md).:

```julia
prob, data, runtime = IBPM_advance( Re, boundary, body, freestream, 
    Δx=Δx, Δt=Δt, T=T, save_info=save_info )
```

This solver returns three things:
- prob::IBProblem : for use in plotting
- data::Array : saved data for use in analysis / plotting
- runtime::Float64 : time elapsed while running simulation

## Step 3: Plot / Analyze the Solution

The data returned from the solver can be used with inbuilt plotting functionality at one time step (in this case the final time step):
```julia
ibpm.plot_state( prob, data[2].val[end], data[2].t[end], var=:omega,
        xlims=(-4.0, 10.0), ylims=(-3.0, 3.0), clims=(-5.0, 5.0), clevs=40)
```
or over a range of times and displayed as a .gif:
```julia
anim = @animate for j = 1 : length(data[2].t)
            ibpm.plot_state( prob, data[2].val[j], data[2].t[j], var=:omega,
            xlims=(-4.0, 10.0), ylims=(-3.0, 3.0), clims=(-5.0, 5.0), clevs=40 )
        end
gif(anim, "cyl_100.gif", fps=10)
```