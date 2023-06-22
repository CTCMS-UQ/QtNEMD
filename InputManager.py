#!/usr/bin/python3
class InputManager:
    """ Overgrown struct to hold user input for the MD routines."""
    def __init__(self):
        # Initialise with default values
        # Particle data
        self.xmax = 4.0
        self.ymax = 4.0
        self.zmax = 4.0
        self.reduced_density = 0.8442
        self.sigma = 1.0
        self.eps = 1.0

        # Simulation parameters
        self.flowrate = 0.1
        self.temp = 1.0
        self.do_nemd = False

        # LJ cutoff parameter. This is just to aid in plotting, probably don't want to change this
        self.rcut = 2.5

    def reset_and_update(self, lmp):
        # Reset the simulation and initialise the parameters
        lmp.command(f"clear")

        lmp.command(f"suffix            opt")
        lmp.command(f"units		        lj")
        lmp.command(f"atom_style	    atomic")
        lmp.command(f"dimension	        3")
        lmp.command(f"lattice		    fcc {self.reduced_density}")
        lmp.command(f"region		    box prism 0 {self.xmax} 0 {self.ymax} -{self.zmax}  {self.zmax} 0 0 0")
        lmp.command(f"create_box	    2 box")
        lmp.command(f"create_atoms	    1 box")
        lmp.command(f"mass		        * 1.0")
        lmp.command(f"velocity	        all create 1.44 87287 loop geom")
        lmp.command(f"region		    slice block 4 6 INF INF INF INF")
        lmp.command(f"set		        region slice type 2")
        lmp.command(f"pair_style	    lj/cut {self.rcut}")
        lmp.command(f"pair_coeff	    * * {self.eps} {self.sigma}")
        lmp.command(f"neighbor	        0.3 bin")
        lmp.command(f"neigh_modify	    delay 0 every 1")
        lmp.command("compute            g2 all rdf 100")
        if self.do_nemd:
            lmp.command("compute        sllodtemp all temp/deform")
            lmp.command("thermo_modify  temp sllodtemp")
            lmp.command(f"fix		    1 all nvt/sllod temp {self.temp} {self.temp} 1.0 tchain 1")
            lmp.command(f"fix		    2 all deform 1 xy erate {self.flowrate} remap v")
        else:
            lmp.command(f"fix		    1 all nvt temp {self.temp} {self.temp} 1.0 tchain 1")

    def update_parameters(self, lmp):
        # Only update the parameters which don't require a reset
        if self.do_nemd:
            lmp.command("unfix 1")
            lmp.command("unfix 2")
            lmp.command("uncompute sllodtemp")
            
            lmp.command(f"pair_coeff * * {self.eps} {self.sigma}")

            lmp.command("compute sllodtemp all temp/deform")
            lmp.command("thermo_modify temp sllodtemp")

            lmp.command(f"fix 1 all nvt/sllod temp {self.temp} {self.temp} 1.0 tchain 1")
            lmp.command(f"fix 2 all deform 1 xy erate {self.flowrate} remap v")
            lmp.command(f"fix 3 all enforce2d")
        else:
            lmp.command("unfix 1")
            lmp.command(f"pair_coeff * * {self.eps} {self.sigma}")
            lmp.command(f"fix 1 all nvt temp {self.temp} {self.temp} 1.0 tchain 1")
            lmp.command(f"fix 3 all enforce2d")

    def toggle_nemd(self, lmp):
        # Toggle NEMD field on or off. This needs to be a separate function to update_parameters() since
        # different fixes are defined (and therefore different fixes need to be deleted) depending on
        # whether or not we're doing NEMD. E.g. if we just do update_parameters() after switching
        # do_nemd to True then it will try to do "unfix 2" which will not be defined (since we we're
        # previously doing equilibrium simulations) and throw an exception.
        if self.do_nemd:
            lmp.command("unfix 1")
            lmp.command("unfix 2")
            lmp.command("uncompute sllodtemp")
            self.do_nemd = False

            lmp.command(f"fix 1 all nvt temp {self.temp} {self.temp} 1.0 tchain 1")
            lmp.command(f"fix 3 all enforce2d")            
        else:
            lmp.command("unfix 1")
            self.do_nemd = True

            lmp.command("compute sllodtemp all temp/deform")
            lmp.command("thermo_modify temp sllodtemp")

            lmp.command(f"fix 1 all nvt/sllod temp {self.temp} {self.temp} 1.0 tchain 1")
            lmp.command(f"fix 2 all deform 1 xy erate {self.flowrate} remap v")
            lmp.command(f"fix 3 all enforce2d")

    def format_params(self):
        """ Make a format string for the input which can be either written to file or displayed in Qt."""
        
        # Write the input parameters, making sure to preserve whitespace and newlines
        param_str  =        f"""units lj"""
        param_str +=        f"""\natom_style atomic"""
        param_str +=        f"""\ndimension 2"""
        param_str +=        f"""\nlattice sq2 {self.reduced_density}"""
        param_str +=        f"""\nregion box prism 0 {self.xmax} 0 {self.ymax} -{self.zmax}  {self.zmax} 0 0 0"""
        param_str +=        f"""\ncreate_box 2 box"""
        param_str +=        f"""\ncreate_atoms 1 box"""
        param_str +=        f"""\nmass * 1.0"""
        param_str +=        f"""\nvelocity all create 1.44 87287 loop geom"""
        param_str +=        f"""\nregion slice block 4 6 INF INF INF INF"""
        param_str +=        f"""\nset region slice type 2"""
        param_str +=        f"""\npair_style lj/cut 2.5"""
        param_str +=        f"""\npair_coeff * * {self.eps} {self.sigma} 1.0"""
        param_str +=        f"""\nneighbor 0.3 bin"""
        param_str +=        f"""\nneigh_modify 0 every 1"""
        if self.do_nemd:
            param_str +=    f"""\ncompute all temp/deform"""
            param_str +=    f"""\nthermo_modify temp sllodtemp"""
            param_str +=    f"""\nfix 1 all nvt/sllod temp {self.temp} {self.temp} 1.0 tchain 1"""
            param_str +=    f"""\nfix 2 all deform 1 xy erate {self.flowrate} remap v"""
        else:
            param_str +=    f"""\nfix 1 all nvt temp {self.temp} {self.temp} 1.0 tchain 1"""

        param_str +=        f"""\nfix 3 all enforce2d\n"""


        return(param_str)
