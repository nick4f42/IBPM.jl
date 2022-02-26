include("../src/IBPM.jl")
using .IBPM
using FileIO #For saving data as a jld2 file

#Below, the Reynolds #, boundary of the finest domain, and body must be specified.
#all other variables have default values and needn't be provided the user

#--necessary variables
    boundary = (-4.0, 4.023, -2.0148, 1.99) #left, right, bottom, and top of domain
    Re = 200.0 #Reynolds #
    # specify body as a vector of named tuples with keys type, lengthscale,
    #center, motion.
    #Type and lengthscale must be specified. The others have defaults associated
    #with a stationary body centered at (x,y)=(0.0,0.0)
    type = :plate #:cylinder, :plate are supported
    lengthscale = 1.0 #key lengthscale. e.g., for cylinder is radius. Supports Float64
    motion = :movinggrid #type of body motion. supports
                     #:static (default)
                     #:movinggrid
                     #:movingbody
        #NOTE: if not :static, prescribed motion is given as a function of time
        #For :movinggrid, precribe as named tuple with keys U (x velocity),
        #V (y velocity), \theta (rotational motion),
        #\dot{\theta} (rotational velocity), xc & yc (coords for center of rot)
        #DEFAULT: (U=(t -> 1.0), V=(t->t->0.0), \theta=(t->0.0),
        #\dot{\theta}=(t->0.0), xc=0.0, yc=0.0  )
    mfcn = (U=t -> sin(t), V=t -> 0.0, θ=t-> 0.0, θ̇=t->0.0)
    center = [0.0; 0.5] #body CoM is centered here. default: [0.0; 0.0]
                        #displace in y to trigger asymmetry
    spec = (α = 90.0*π/180.0,) #:plate type supports AoA within spec variable
    body = [(type=type, lengthscale=lengthscale, motion=motion, spec=spec,
        motionfcn=mfcn, center=center)]
#--

#--optional variables
    freestream = (Ux=t->0.0,) #freestream conditions
                                #can be provided as constants or
                                #functions of time
                                #(default: (Ux=t->0.0,
                                #           Uy=t->0.0,
                                #           inclination=t->0.0))

    T = 2.0π #final time to run to (default = 20.0*dt)

    #simulation parameters (these are all optional)
    Δx = 0.02 #default (==missing) gives a grid Re of 2
    Δt = 0.002 #default (==missing) aims for a CFL of 0.1 with a
              #fairly conservative safety factor on max vel
    mg=3      #Number of sub-domains. Default is 5
#--

#--save_info (optional)
    #gives the code information for a data structure to return
    #user provides as a Named Tuple with three keys: save_fcns, save_times,
    #   save_types
    #default (if save_info is unspecified): the code returns necessary
    #   information for a restart.

    #Even if save_info is provided, only the save_fcns key is necessary.

    #vector of functions that specify which data to save.
    #default: save full state at the final time
    #NOTE: if you want to plot the flowfield later, you have to save the full
    #state as (t,state)-> state
    svfc = [ (t,state)->state.CL; (t,state)->state  ]

    #corresponding vector of save instances (1 per function)
    #each can be prescribed as scalar Float or array.
    #default: save at the final time instance T for each save variable
    #if any of the entries is either 0.0 or missing, that variable will be saved
    #every timestep
    svti = [ 0.0, T/30.0 ]
    #types of each variable
    #default: Any for each save variable
    svty = [ Vector{Float64}; Any ]

    #Store in Named Tuple
    save_info = ( save_fcns = svfc, save_types = svty, save_times = svti )

#--

#--run the simulation based on user data
    prob, data, runtime = IBPM_advance( Re, boundary, body, freestream,
        Δx=Δx, Δt=Δt, T=T, save_info=save_info )

        #Uncomment to save data as a JLD2 file (caution, could be large!)
        # FileIO.save("examples/output.jld2",  "data", data, "prob", prob)
#--


#--after running, you can use the data struct to plot/analyze quantities.
    using Plots

    #Looking at timetraces of lift or drag is a good starting point. Lift was
    #the first variable saved in our data struct so we can plot that as...
    display(
        plot( data[1].t, [data[1].val[j][1] for j in 1:length(data[1].t)],
        legend=:false )
        )

    #Plotting the flowfield is trickier, so there's more built-in support for
    #that. First off, we can plot a snapshot of the flowfield in terms of either
    #vorticity (var=:omega, default), x and y velocity (var=:vel), or the
    #streamfunction (var=:psi). Let's do that at the final time:
    IBPM.plot_state( prob, data[2].val[end], data[2].t[end], var=:omega,
        xlims=(-4.0, 4.0), ylims=(-3.0, 3.0), clims=(-5.0, 5.0), clevs=40)

    #We can also use Julia's handy @animate macro to make a gif of the vorticity
    # field from our save data:
    anim = @animate for j = 1 : length(data[2].t)
                IBPM.plot_state( prob, data[2].val[j], data[2].t[j], var=:omega,
                xlims=(-6.0, 6.0), ylims=(-3.0, 3.0), clims=(-5.0, 5.0), clevs=40 )
            end
    #To save the animation, use gif():
    gif(anim, "examples/osc_plate_movinggrid.gif", fps=10)
#--
