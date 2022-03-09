# Example: 2D Uniform flow around a stationary cylinder (New UI)

In this example we will solve the problem of a uniform freestream of 1 m/s flowing around a 2 dimentional cylinder of radius 0.5 m over the time interval from 0 to 10 seconds.
 
The general workflow is to define the problem, choose desired state information to save, solve the problem, then analyze the solution. The full code for solving this example problem is:

```julia
using .IBPM
using .IBPM.Quantities: lift_coef

Δx = 0.02
boundary = ((-1.0, 3.0), (-2.0, 2.0))
gridcount = 5
grid = MultiGrid(Δx, boundary; mg=gridcount)

cylinder = make_cylinder(0.5, Δx, (0.0, 0.0), motion=Static())
bodies = [cylinder]

T = 10.0
Δt = 0.004
times = range(0, T, step=Δt)

Re = 100.0
freestream = t -> (1.0, 0.0)

problem = IBProblem(grid, bodies, times; Re, freestream)

lift = StateData(state -> lift_coef(state)[1])
states = StateData(saveat=LinRange(0, T, 31))
save_data = [lift, states]

solve!(save_data, problem)

anim = @animate for (t, state) in zip(states.t, states)
    print("t = ", t, '\r')
    IBPM.plot_state(
        problem, state, t, var=:omega,
        xlims=(-4.0, 10.0), ylims=(-3.0, 3.0), clims=(-5.0, 5.0), clevs=40
    )
end

gif(anim, "cyl_100.gif", fps=10)
```

## Step 1: Defining the Problem Parameters

The problem is defined by:
- Grid
   - Domain boundaries
   - Grid spacing
   - (For a MultiGrid, number of grids)
- Rigid Bodies
   - Position
   - Size
   - Motion type (optional)
- Simulation times
- Reynolds Number
- Freestream
   - X axis velocity (can be function of time)
   - Y axis velocity (can be funciton of time)
   - Inclination (optional)

### Grid

For this example, we are using the MultiGrid method for discretizing the domain. For other grid choices, see #todo#
The multigrid is optionally defined using the number of grids. By default the multigrid only uses one grid.
```julia
Δx = 0.02
boundary = ((-1.0, 3.0), (-2.0, 2.0))
gridcount = 5
grid = MultiGrid(Δx, boundary; mg=gridcount)
```

### Bodies

Bodies are defined individually. The definition of a 2 dimensional cylinder necessitates the radius and position. 
The motion of the cylinder over time can be optionally defined, being either static or in motion.
See #motion link# for options on body motion. The motion is set to static by default.

```julia
cylinder = make_cylinder(0.5, Δx, (0.0, 0.0), motion=Static())
```

Once all bodies are defined, in this case the one cylinder, they can be arranged in an array of bodies.

```julia
bodies = [cylinder]
```

### Simulation times

The simulation times needs to be defined for the problem. To do this, simply use the user defined total 
simulation time and time step.

```julia
T = 10.0
Δt = 0.004
times = range(0, T, step=Δt)
```

### Reynolds Number, Freestream

Optional definitions of Reynolds number and freestream velocity. Freestream must be defined as a function of time. For options on freestream definition, see #freestream#.

```julia
Re = 100.0
freestream = t -> (1.0, 0.0)
```

### Creating the Problem

The problem can be defined once all of the above have been defined. 
```julia
problem = IBProblem(grid, bodies, times; Re, freestream)
```

## Step 2: Define Save Data
Save data must be predefined before being solved. This is done by creating instances of SaveData. There are preexisting supported calculateable
quantities built into module functionality, but user defined functions can be defined, calculated and saved. For quantities that can have solutions for
individual bodies (for example lift and drag), access the solution for the desired body by accessing the respective index when creating the array of bodies above.
```julia
lift = StateData(state -> lift_coef(state)[1])
drag = StateData(state -> drag_coef(state)[1])
states = StateData(saveat=LinRange(0, T, 31))
save_data = [lift, drag, states]
```

## Step 3: Solve the problem

Once all required arguments are defined, they can be used to call the solve! method:

```julia

solve!(save_data, problem)
```

This solver does not return anything, but rather modifies the array of inputted saved quantities save_data depending on how each individual StateData is defined.

## Step 4: Plot / Analyze the Solution

The data can be accessed by referencing each specific desired calcualted quantitiy. For example, the lift coefficient at the third time step is:
```julia
lift[3]
```
or the 2 dimentional state defined by a 2 dimentional matrix at the last time step is:
```julia
state[end]
```
To plot:
```julia
anim = @animate for (t, state) in zip(states.t, states)
    print("t = ", t, '\r')
    IBPM.plot_state(
        problem, state, t, var=:omega,
        xlims=(-4.0, 10.0), ylims=(-3.0, 3.0), clims=(-5.0, 5.0), clevs=40
    )
end
```