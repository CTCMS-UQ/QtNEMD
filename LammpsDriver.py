#!/usr/bin/python

""" This is a template class to interface with whatever MD backend code we decide to 
    use. It should only include the broadest outline of functionality, which can be 
    overridden to support, e.g. LAMMPS or Debra's MD code."""

import lammps

class MDInterface:
    def __init__(self):
        # LAMMPS instance. Don't print output to the terminal
        lmp_args = ["-log", "none", "-screen", "none"]
        self._lmp = lammps.lammps(cmdargs=lmp_args)

        # Initialise user-variables with default values
        # These parameters require a restart when they change
        self._xmax = 4.0
        self._ymax = 4.0
        self._zmax = 4.0
        self._rho = 0.8442

        # These parameters do not require a restart when they change
        self._sigma = 1.0
        self._eps= 1.0

        self._flowrate = 0.001
        self._temp = 1.0
        self._do_nemd = False

        # LJ cutoff parameter. This is just to aid in plotting, probably don't want to change this
        self._rcut = 2.5
        
        # Initialise the thermodynamic output parameters to None until they are calculated
        self._npart = None
        self._vol = None

        self._box_bounds = (None, None, None)

        self._x = None
        self._y = None
        self._z = None

        self._g2_arrays = None

        # Finally, pass these values into the simulation
        self.reset_and_update_parameters()

###################### Getters and setters (properties) #######################

    # Box coords/initial conditions. Need to reset the simulation if these change
    @property
    def xmax(self):
        return(self._xmax)
    @xmax.setter
    def xmax(self, value):
        self._xmax = value
        self.reset_and_update_parameters()
    @property
    def ymax(self):
        return(self._ymax)
    @ymax.setter
    def ymax(self, value):
        self._ymax = value
        self.reset_and_update_parameters()
    @property
    def zmax(self):
        return(self._zmax)
    @zmax.setter
    def zmax(self, value):
        self._ymax = value
        self.reset_and_update_parameters()
    @property
    def reduced_density(self):
        return(self._rho)
    @reduced_density.setter
    def reduced_density(self, value):
        self._rho=value
        self.reset_and_update_parameters()

    # Lennard-Jones parameters. These do not require a restart when they change
    @property
    def sigma(self):
        return(self._sigma)
    @sigma.setter
    def sigma(self, value):
        self._sigma=value
        self.update_parameters()
    @property
    def eps(self):
        return(self._eps)
    @eps.setter
    def eps(self, value):
        self._eps=value
        self.update_parameters()

    # Misc simulation parameters. Also don't require a restart
    @property
    def flowrate(self):
        return(self._flowrate)
    @flowrate.setter
    def flowrate(self, value):
        self._flowrate = value
        self.update_parameters()
    @property
    def temp(self):
        return(self._lmp.get_thermo("temp"))
    @temp.setter
    def temp(self, value):
        self._temp = value
        self.update_parameters()
    
    # Get the status of whether we're doing nonequilibrium calculations. No setter, since this is handled by toggle_nemd()
    @property
    def do_nemd(self):
        return(self._do_nemd)

    # Atomic coordinates. These should be read-only to external classes
    @property
    def x(self):
        coords = self._lmp.numpy.extract_atom("x")
        return(coords[:,0])
    @property
    def y(self):
        coords = self._lmp.numpy.extract_atom("x")
        return(coords[:,1])
    @property
    def z(self):
        coords = self._lmp.numpy.extract_atom("x")
        return(coords[:,2])
    @property
    def box_bounds(self):
        return(self._box_bounds)
    
    @property
    def npart(self):
        return(self._lmp.get_natoms())
    
    @property
    def vol(self):
        return(self._lmp.get_thermo("vol"))

###################### Callable methods ###########################
    def setup(self):
        # This function gets called at the start of the program's run, as well as whenever we restart the simulation.
        self.reset_and_update_parameters()
        self._lmp.command("run 0")

        # Now update all of our thermodynamic parameters with values from LAMMPS
        self.update_from_lammps()

    def update_from_lammps(self):
        # Update this class's parameters with the values from LAMMPS. This doesn't change the underlying simulation
        self._npart = self._lmp.get_natoms()
        self._temp = self._lmp.get_thermo("temp")
        self._vol = self._lmp.get_thermo("vol")

        # And particle positions
        coords = self._lmp.numpy.extract_atom("x")
        self._x = coords[:,0]
        self._y = coords[:,1]
        self._z = coords[:,2]
        self._box_bounds = (self._lmp.get_thermo("lx"), self._lmp.get_thermo("ly"), self._lmp.get_thermo("lz"))

        # Finally, grab the g(2) radial distribution function
        self._g2_arrays = self._lmp.numpy.extract_compute("g2", lammps.LMP_STYLE_GLOBAL, lammps.LMP_TYPE_ARRAY)

    # g(2) radial distribution function. This returns into two arrays for r and g(2)(r), indexed in a dict
    def g2_compute(self):
        return({'r': self._g2_arrays[:,0], 'g2': self._g2_arrays[:,1]})

    def reset_and_update_parameters(self):
        # Reset the simulation and initialise the parameters
        self._lmp.command(f"clear")

        self._lmp.command(f"suffix            opt")
        self._lmp.command(f"units		        lj")
        self._lmp.command(f"atom_style	    atomic")
        #self._lmp.command(f"dimension	        2")
        self._lmp.command(f"lattice		    fcc {self._rho}")
        self._lmp.command(f"region		    box prism 0 {self._xmax} 0 {self._ymax} -{self._zmax}  {self._zmax} 0 0 0")
        self._lmp.command(f"create_box	    2 box")
        self._lmp.command(f"create_atoms	    1 box")
        self._lmp.command(f"mass		        * 1.0")
        self._lmp.command(f"velocity	        all create 1.44 87287 loop geom")
        self._lmp.command(f"region		    slice block 4 6 INF INF INF INF")
        self._lmp.command(f"set		        region slice type 2")
        self._lmp.command(f"pair_style	    lj/cut {self._rcut}")
        self._lmp.command(f"pair_coeff	    * * {self._eps} {self._sigma}")
        self._lmp.command(f"neighbor	        0.3 bin")
        self._lmp.command(f"neigh_modify	    delay 0 every 1")
        self._lmp.command("compute            g2 all rdf 100")
        if self.do_nemd:
            self._lmp.command("compute        sllodtemp all temp/deform")
            self._lmp.command("thermo_modify  temp sllodtemp")
            self._lmp.command(f"fix		    1 all nvt/sllod temp {self._temp} {self._temp} 1.0 tchain 1")
            self._lmp.command(f"fix		    2 all deform 1 xy erate {self._flowrate} remap v")
        else:
            self._lmp.command(f"fix		    1 all nvt temp {self._temp} {self._temp} 1.0 tchain 1")


    def update_parameters(self):
        # Only update the parameters which don't require a reset
        if self.do_nemd:
            self._lmp.command("unfix 1")
            self._lmp.command("unfix 2")
            self._lmp.command("uncompute sllodtemp")
            
            self._lmp.command(f"pair_coeff * * {self._eps} {self._sigma}")

            self._lmp.command("compute sllodtemp all temp/deform")
            self._lmp.command("thermo_modify temp sllodtemp")

            self._lmp.command(f"fix 1 all nvt/sllod temp {self._temp} {self._temp} 1.0 tchain 1")
            self._lmp.command(f"fix 2 all deform 1 xy erate {self._flowrate} remap v")
        else:
            self._lmp.command("unfix 1")
            self._lmp.command(f"pair_coeff * * {self._eps} {self._sigma}")
            self._lmp.command(f"fix 1 all nvt temp {self._temp} {self._temp} 1.0 tchain 1")

    def toggle_nemd(self):
        # Toggle NEMD field on or off. This needs to be a separate function to update_parameters() since
        # different fixes are defined (and therefore different fixes need to be deleted) depending on
        # whether or not we're doing NEMD. E.g. if we just do update_parameters() after switching
        # do_nemd to True then it will try to do "unfix 2" which will not be defined (since we we're
        # previously doing equilibrium simulations) and throw an exception.
        if self._do_nemd:
            self._lmp.command("unfix 1")
            self._lmp.command("unfix 2")
            self._lmp.command("uncompute sllodtemp")
            self._do_nemd = False

            self._lmp.command(f"fix 1 all nvt temp {self._temp} {self._temp} 1.0 tchain 1")
        else:
            self._lmp.command("unfix 1")
            self._do_nemd = True

            self._lmp.command("compute sllodtemp all temp/deform")
            self._lmp.command("thermo_modify temp sllodtemp")

            self._lmp.command(f"fix 1 all nvt/sllod temp {self._temp} {self._temp} 1.0 tchain 1")
            self._lmp.command(f"fix 2 all deform 1 xy erate {self._flowrate} remap v")

    def format_params(self):
        """ Make a format string for the input which can be either written to file or displayed in Qt."""
        
        # Write the input parameters, making sure to preserve whitespace and newlines
        param_str  =        f"""units lj"""
        param_str +=        f"""\natom_style atomic"""
        param_str +=        f"""\ndimension 3"""
        param_str +=        f"""\nlattice fcc {self._rho}"""
        param_str +=        f"""\nregion box prism 0 {self._xmax} 0 {self._ymax} -{self._zmax}  {self._zmax} 0 0 0"""
        param_str +=        f"""\ncreate_box 2 box"""
        param_str +=        f"""\ncreate_atoms 1 box"""
        param_str +=        f"""\nmass * 1.0"""
        param_str +=        f"""\nvelocity all create 1.44 87287 loop geom"""
        param_str +=        f"""\nregion slice block 4 6 INF INF INF INF"""
        param_str +=        f"""\nset region slice type 2"""
        param_str +=        f"""\npair_style lj/cut 2.5"""
        param_str +=        f"""\npair_coeff * * {self._eps} {self._sigma} 1.0"""
        param_str +=        f"""\nneighbor 0.3 bin"""
        param_str +=        f"""\nneigh_modify 0 every 1"""
        if self.do_nemd:
            param_str +=    f"""\ncompute all temp/deform"""
            param_str +=    f"""\nthermo_modify temp sllodtemp"""
            param_str +=    f"""\nfix 1 all nvt/sllod temp {self._temp} {self._temp} 1.0 tchain 1"""
            param_str +=    f"""\nfix 2 all deform 1 xy erate {self._flowrate} remap v"""
        else:
            param_str +=    f"""\nfix 1 all nvt temp {self._temp} {self._temp} 1.0 tchain 1"""



        return(param_str)

    def run(self, nsteps):
        # First, advance the simulation
        if(nsteps >= 0):
            self._lmp.command(f"run {nsteps}")
        else:
            raise(ValueError)
        
        # Now update all of our thermodynamic parameters with values from LAMMPS
        self.update_from_lammps()
