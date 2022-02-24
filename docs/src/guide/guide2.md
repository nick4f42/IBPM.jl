# Example 1: 2D Uniform flow around a stationary cylinder (New UI)

In this example we will solve the problem of a uniform freestream of 1 m/s flowing around a 2 dimentional cylinder of radius 0.5 m over the time interval from 0 to 10 seconds.
 
The general workflow is to define the problem, choose desired state information to save, solve the problem, then analyze the solution. The full code for solving this example problem is:

```@example
using IBPM.Quantities: lift_coef

T = 1.0
Δx = 0.02
Δt = 0.004
Re = 100

boundary = ((-1.0, 3.0), (-2.0, 2.0))

grid = MultiGrid(Δx, boundary; grids=5)

body1 = Bodies.cylinder(0.5, Δx, (0.0, 0.0), motion=Static())
bodies = [body1]

freestream = Freestream((1.0, 0.0))

problem = IBProblem(grid, bodies, Δt; freestream, Re)

lift = StateData(lift_coeff)
states = StateData()
save_data = [lift, states]

time_range = 0:Δt:T

solve!(save_data, problem, time_range)
```

# Step 1: Defining the Problem Parameters

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

We can initialize these parameters and use them to define the problem:

```@example
T = 1.0
Δx = 0.02
Δt = 0.004
Re = 100

boundary = ((-1.0, 3.0), (-2.0, 2.0))

grid = MultiGrid(Δx, boundary; grids=5)

body1 = Bodies.cylinder(0.5, Δx, (0.0, 0.0), motion=Static())
bodies = [body1]

freestream = Freestream((1.0, 0.0))

problem = IBProblem(grid, bodies, Δt; freestream, Re)
```

The desired saved calculated quantities can be optionally defined:

```@example
lift = StateData(lift_coeff)
states = StateData()
save_data = [lift, states]
```

# Step 2: Solve the problem

Once all required arguments are defined, they can be used to call the solve! method for a user defined time range or by default based on Δt.:

```@example
time_range = 0:Δt:T #Interpolation maybe?

solve!(save_data, problem, time_range)
```

This solver does not return anything, but rather modifies the array of inputted saved quantities save_data and stores all data for each time step.

# Step 3: Plot / Analyze the Solution

The data can be accessed by referencing each specific desired calcualted quantitiy. For example, the lift coefficient at the third time step is:
```@example
lift[3]
```
or the 2 dimentional state defined by a 2 dimentional matrix at the last time step is:
```@example
state[end]
```
To plot:
```@example
#TODO
```