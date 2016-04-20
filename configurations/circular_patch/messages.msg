INFO: 2016-Jan-03 17:11:20: Messages will be sent to messages.msg.
INFO: 2016-Jan-03 17:11:20: Simulation box is fixed with size (lx,ly,lz) = (1000,1000,10).
INFO: 2016-Jan-03 17:11:20: Reading particle coordinates from file: patch.dat
INFO: 2016-Jan-03 17:11:20: Column 1 of input file is : id.
INFO: 2016-Jan-03 17:11:20: Column 2 of input file is : x.
INFO: 2016-Jan-03 17:11:20: Column 3 of input file is : y.
INFO: 2016-Jan-03 17:11:20: Column 4 of input file is : nx.
INFO: 2016-Jan-03 17:11:20: Column 5 of input file is : ny.
INFO: 2016-Jan-03 17:11:20: Column 6 of input file is : nvx.
INFO: 2016-Jan-03 17:11:20: Column 7 of input file is : nvy.
INFO: 2016-Jan-03 17:11:20: Column 8 of input file is : nvz.
INFO: 2016-Jan-03 17:11:20: Column 9 of input file is : area.
INFO: 2016-Jan-03 17:11:20: Read data for 1000 particles.
INFO: 2016-Jan-03 17:11:20: Generated group 'all' containing all particles.
INFO: 2016-Jan-03 17:11:20: There are 1 distinct particle types in the system.
INFO: 2016-Jan-03 17:11:20: Periodic flag for the system is set to false.
INFO: 2016-Jan-03 17:11:20: Finished reading system coordinates from patch.dat.
INFO: 2016-Jan-03 17:11:20: Created cell list with 162450 cells.
INFO: 2016-Jan-03 17:11:20: Each cell has dimensions (3.5087719298245612,3.5087719298245612,5).
INFO: 2016-Jan-03 17:11:20: Using cell lists for neighbour list builds.
INFO: 2016-Jan-03 17:11:20: Neighbour list will also build faces.
INFO: 2016-Jan-03 17:11:20: Neighbour list.  Setting contact distance to 2.8.
INFO: 2016-Jan-03 17:11:20: Neighbour list.  Setting maximum face perimeter to 30.0.
INFO: 2016-Jan-03 17:11:20: Created neighbour list with cutoff 3 and padding distance 0.5.
WARNING: 2016-Jan-03 17:11:20: Constraint. Maximum number of iterations has not been set. Assuming 100
WARNING: 2016-Jan-03 17:11:20: Constraint. Tolerance has not been set. Assuming 1e-6.
WARNING: 2016-Jan-03 17:11:20: Constraint. No particle groups set. Assuming "all".
INFO: 2016-Jan-03 17:11:20: Planar constraint. Setting plane dimensions to lx = 1000 and ly = 1000.
INFO: 2016-Jan-03 17:11:20: Plane constraint. Particle motion will not be affected by the box boundaries.
INFO: 2016-Jan-03 17:11:20: Added  : plane to the list of constraints.
INFO: 2016-Jan-03 17:11:20: Adding constraint of type plane.
INFO: 2016-Jan-03 17:11:20: Parameter unlimited for constraint plane is set to .
WARNING: 2016-Jan-03 17:11:20: Value object. No minimum value set. Assuming 0.
WARNING: 2016-Jan-03 17:11:20: Value object. No maximum value set. Assuming 1.0.
WARNING: 2016-Jan-03 17:11:20: Value object. No number of steps set. Assuming 1000.
INFO: 2016-Jan-03 17:11:20: Using constant value object.
INFO: 2016-Jan-03 17:11:20: Global area stiffness (K) for vertex-particle pair potential is set to 5.0.
INFO: 2016-Jan-03 17:11:20: Global perimeter stiffness (gamma) for vertex-particle pair potential is set to 0.0.
INFO: 2016-Jan-03 17:11:20: Global contact stiffness (lambda) for vertex-particle pair potential is set to 0.0.
INFO: 2016-Jan-03 17:11:20: Added pair potential : vp to the list of pair interactions.
INFO: 2016-Jan-03 17:11:20: Pair potential vp has neighbour list. Neighbour list updates will be performed during the simulation.
INFO: 2016-Jan-03 17:11:20: Added vp to the list of pair potentials.
INFO: 2016-Jan-03 17:11:20: Parameter K for pair potential vp is set to 5.0.
INFO: 2016-Jan-03 17:11:20: Parameter gamma for pair potential vp is set to 0.0.
INFO: 2016-Jan-03 17:11:20: Parameter lambda for pair potential vp is set to 0.0.
INFO: 2016-Jan-03 17:11:20: Added logger. Logged quantities will be sent to cell.log.
INFO: 2016-Jan-03 17:11:20: Logs will be produced every 10 time steps.
INFO: 2016-Jan-03 17:11:20: Adding log quantity : area.
INFO: 2016-Jan-03 17:11:20: Adding log quantity : avg_perim.
INFO: 2016-Jan-03 17:11:20: Adding log quantity : size.
INFO: 2016-Jan-03 17:11:20: Adding log quantity : vp_energy.
INFO: 2016-Jan-03 17:11:20: Adding logger to file cell.log.
WARNING: 2016-Jan-03 17:11:20: Value object. No minimum value set. Assuming 0.
WARNING: 2016-Jan-03 17:11:20: Value object. No maximum value set. Assuming 1.0.
WARNING: 2016-Jan-03 17:11:20: Value object. No number of steps set. Assuming 1000.
INFO: 2016-Jan-03 17:11:20: Using constant value object.
WARNING: 2016-Jan-03 17:11:20: No group has been set. Assuming group 'all'.
WARNING: 2016-Jan-03 17:11:20: Brownian dynamics is a legacy integrator keep for backwards compatibility. It will be removed in later releaser. Please consider using brawnian_pos and brownian_align instead.
INFO: 2016-Jan-03 17:11:20: Brownian dynamics integrator. Setting magnitude of active velocity to 0.0.
INFO: 2016-Jan-03 17:11:20: Brownian dynamics integrator. Setting rotational diffusion rate to 0.0.
INFO: 2016-Jan-03 17:11:20: Brownian dynamics integrator. Setting mobility to 1.0.
WARNING: 2016-Jan-03 17:11:20: Brownian dynamics integrator. Rotational mobility not set. Using default value 1.
INFO: 2016-Jan-03 17:11:20: Brownian dynamics integrator. Setting random number generator seed to 1.
WARNING: 2016-Jan-03 17:11:20: Brownian dynamics integrator. Assuming polar order parameter.
INFO: 2016-Jan-03 17:11:20: Adding integrator of type brownian.
INFO: 2016-Jan-03 17:11:20: Parameter dt for integrator brownian is set to 0.001.
INFO: 2016-Jan-03 17:11:20: Parameter mu for integrator brownian is set to 1.0.
INFO: 2016-Jan-03 17:11:20: Parameter nu for integrator brownian is set to 0.0.
INFO: 2016-Jan-03 17:11:20: Parameter seed for integrator brownian is set to 1.
INFO: 2016-Jan-03 17:11:20: Parameter v0 for integrator brownian is set to 0.0.
INFO: 2016-Jan-03 17:11:20: Starting simulation run for 10000 steps.
INFO: 2016-Jan-03 17:11:27: Built neighbour list 291 time. Average number of steps between two builds : 34.364261168384878.
INFO: 2016-Jan-03 17:11:27: Dump type set to vtp
INFO: 2016-Jan-03 17:11:27: Dump will be sent to the file with base name cell with extension vtp.
INFO: 2016-Jan-03 17:11:27: Dumping will start after 0 time steps.
INFO: 2016-Jan-03 17:11:27: Dump will be produced every 500 time steps.
INFO: 2016-Jan-03 17:11:27: Each dumped time step will be stored in a separate file labelled by the time step.
INFO: 2016-Jan-03 17:11:27: Adding dump to file cell.
INFO: 2016-Jan-03 17:11:27: Parameter freq for dump cell is set to 500.
INFO: 2016-Jan-03 17:11:27: Parameter multi for dump cell is set to .
INFO: 2016-Jan-03 17:11:27: Parameter start for dump cell is set to 0.
INFO: 2016-Jan-03 17:11:27: Parameter type for dump cell is set to vtp.
INFO: 2016-Jan-03 17:11:27: Dump type set to vtp
INFO: 2016-Jan-03 17:11:27: Dump will be sent to the file with base name cell_dual with extension vtp.
INFO: 2016-Jan-03 17:11:27: Dumping will start after 0 time steps.
INFO: 2016-Jan-03 17:11:27: Dump will be produced every 500 time steps.
INFO: 2016-Jan-03 17:11:27: Each dumped time step will be stored in a separate file labelled by the time step.
INFO: 2016-Jan-03 17:11:27: VTP file will contain dual mesh.
INFO: 2016-Jan-03 17:11:27: Adding dump to file cell_dual.
INFO: 2016-Jan-03 17:11:27: Parameter dual for dump cell_dual is set to .
INFO: 2016-Jan-03 17:11:27: Parameter freq for dump cell_dual is set to 500.
INFO: 2016-Jan-03 17:11:27: Parameter multi for dump cell_dual is set to .
INFO: 2016-Jan-03 17:11:27: Parameter start for dump cell_dual is set to 0.
INFO: 2016-Jan-03 17:11:27: Parameter type for dump cell_dual is set to vtp.
INFO: 2016-Jan-03 17:11:27: Dump type set to full
INFO: 2016-Jan-03 17:11:27: Dump will be sent to the file with base name cell with extension dat.
INFO: 2016-Jan-03 17:11:27: Dumping will start after 0 time steps.
INFO: 2016-Jan-03 17:11:27: Dump will be produced every 2500 time steps.
INFO: 2016-Jan-03 17:11:27: Each dumped time step will be stored in a separate file labelled by the time step.
INFO: 2016-Jan-03 17:11:27: Include info header into the dump file.
INFO: 2016-Jan-03 17:11:27: Adding dump to file cell.
INFO: 2016-Jan-03 17:11:27: Parameter a0 for dump cell is set to .
INFO: 2016-Jan-03 17:11:27: Parameter boundary for dump cell is set to .
INFO: 2016-Jan-03 17:11:27: Parameter cell_area for dump cell is set to .
INFO: 2016-Jan-03 17:11:27: Parameter cell_perim for dump cell is set to .
INFO: 2016-Jan-03 17:11:27: Parameter cont_num for dump cell is set to .
INFO: 2016-Jan-03 17:11:27: Parameter coordinate for dump cell is set to .
INFO: 2016-Jan-03 17:11:27: Parameter freq for dump cell is set to 2500.
INFO: 2016-Jan-03 17:11:27: Parameter header for dump cell is set to .
INFO: 2016-Jan-03 17:11:27: Parameter multi for dump cell is set to .
INFO: 2016-Jan-03 17:11:27: Parameter start for dump cell is set to 0.
INFO: 2016-Jan-03 17:11:27: Parameter type for dump cell is set to full.
INFO: 2016-Jan-03 17:11:27: Parameter velocity for dump cell is set to .
INFO: 2016-Jan-03 17:11:27: Population. Applying population to group all.
INFO: 2016-Jan-03 17:11:27: Population. Population control frequency set to 1500.
WARNING: 2016-Jan-03 17:11:27: Cell population control. No random number generator seed specified. Using default 0.
INFO: 2016-Jan-03 17:11:27: Cell population control. Setting division rate 0.01.
INFO: 2016-Jan-03 17:11:27: Cell population control. Setting death rate 0.0.
WARNING: 2016-Jan-03 17:11:27: Cell population control. No split distance set. Using default 0.5.
INFO: 2016-Jan-03 17:11:27: Cell population control. Maximum native cell area set to 7.0.
WARNING: 2016-Jan-03 17:11:27: Cell population control. No maximum cell age set. Using default 10,000.
WARNING: 2016-Jan-03 17:11:27: Cell population control. No split of the offspring separation set. Assuming 0.5.
INFO: 2016-Jan-03 17:11:27: Cell population control. Growth rate set to 0.005.
WARNING: 2016-Jan-03 17:11:27: Cell population control. No growth proability set. Assuming 0.5.
WARNING: 2016-Jan-03 17:11:27: Cell population control. Contact distance and neigbour list distance will not be rescaled.
INFO: 2016-Jan-03 17:11:27: Adding population of type cell.
INFO: 2016-Jan-03 17:11:27: Parameter death_rate for population cell is set to 0.0.
INFO: 2016-Jan-03 17:11:27: Parameter division_rate for population cell is set to 0.01.
INFO: 2016-Jan-03 17:11:27: Parameter freq for population cell is set to 1500.
INFO: 2016-Jan-03 17:11:27: Parameter group for population cell is set to all.
INFO: 2016-Jan-03 17:11:27: Parameter growth_rate for population cell is set to 0.005.
INFO: 2016-Jan-03 17:11:27: Parameter max_area for population cell is set to 7.0.
INFO: 2016-Jan-03 17:11:27: Starting simulation run for 1000000 steps.
