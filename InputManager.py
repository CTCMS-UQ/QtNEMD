#!/usr/bin/python3
class InputManager:
    """ Overgrown struct to hold user input for the MD routines."""
    def __init__(self):
        # Initialise with default values
        # Particle data
        self.xmax = 10
        self.ymax = 10
        self.zmax = 0.5
        self.reduced_density = 0.8442
        self.sigma = 1.0
        self.eps = 1.0

        # Simulation parameters
        self.flowrate = 0.01
        self.temp = 1.0
        self.do_nemd = False

    def reset_and_update(self, lmp):
        # Reset the simulation and initialise the parameters
        lmp.command(f"clear")

        lmp.command(f"suffix            opt")
        lmp.command(f"units		        lj")
        lmp.command(f"atom_style	    atomic")
        lmp.command(f"dimension	        2")
        lmp.command(f"lattice		    sq2 {self.reduced_density}")
        lmp.command(f"region		    box prism 0 {self.xmax} 0 {self.ymax} -{self.zmax}  {self.zmax} 0 0 0")
        lmp.command(f"create_box	    2 box")
        lmp.command(f"create_atoms	    1 box")
        lmp.command(f"mass		        * 1.0")
        lmp.command(f"velocity	        all create 1.44 87287 loop geom")
        lmp.command(f"region		    slice block 4 6 INF INF INF INF")
        lmp.command(f"set		        region slice type 2")
        lmp.command(f"pair_style	    lj/cut 2.5")
        lmp.command(f"pair_coeff	    * * {self.eps} {self.sigma} 1.0")
        lmp.command(f"neighbor	        0.3 bin")
        lmp.command(f"neigh_modify	    delay 0 every 1")
        if self.do_nemd:
            lmp.command(f"fix		        1 all nvt/sllod temp {self.temp} {self.temp} 1.0 tchain 1")
            lmp.command(f"fix		        2 all deform 1 xy erate {self.flowrate} remap v")
        else:
            lmp.command(f"fix		        1 all nvt temp {self.temp} {self.temp} 1.0 tchain 1")

        lmp.command(f"fix               3 all enforce2d")

    def update_parameters(self, lmp):
        # Only update the parameters which don't require a reset
        if self.do_nemd:
            lmp.command("unfix 1")
            lmp.command("unfix 2")
            lmp.command("unfix 3")
            
            lmp.command(f"pair_coeff * * {self.eps} {self.sigma}")
            lmp.command(f"fix 1 all nvt/sllod temp {self.temp} {self.temp} 1.0 tchain 1")
            lmp.command(f"fix 2 all deform 1 xy erate {self.flowrate} remap v")
            lmp.command(f"fix 3 all enforce2d")
        else:
            lmp.command("unfix 1")
            lmp.command("unfix 3")
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
            lmp.command("unfix 3")
            self.do_nemd = False

            lmp.command(f"fix 1 all nvt temp {self.temp} {self.temp} 1.0 tchain 1")
            lmp.command(f"fix 3 all enforce2d")            
        else:
            lmp.command("unfix 1")
            lmp.command("unfix 3")
            self.do_nemd = True

            lmp.command(f"fix 1 all nvt/sllod temp {self.temp} {self.temp} 1.0 tchain 1")
            lmp.command(f"fix 2 all deform 1 xy erate {self.flowrate} remap v")
            lmp.command(f"fix 3 all enforce2d")

#    def format_params(self):
#        """ Make a format string which can be either written to file or displayed in Qt."""
#        
#        # Write the input parameters, making sure to preserve whitespace and newlines
#        param_str  = f"""{self.tr:.3f} {self.drw:.3f} {self.drf:.3f} {self.delta:.3f} {self.latt} {self.npart} 25"""
#        param_str += f"""\n{self.fe0:.3f} {self.rcut:.3f} 1 {self.kh:.3f} {self.nprint} """
#        param_str += f"""\n{self.mix:.3f} {self.eps1:.3f} {self.eps2:.3f} {self.qvol:.3f}"""
#        param_str += f"""\n{self.kf:.3f} {self.r0:.3f} {self.limol:.3f} {self.yzdivx:.3f}"""
#        param_str += f"""\n{self.dxxdiv:.3f}"""
#        param_str += f"""\n{self.ntype} {self.non} {self.ngaus} {self.e0:.3f}"""
#        param_str += f"""\n{self.nplot} {self.maxtau} {self.eqtim} {self.ncyc}"""
#
#        # Now write the header
#        param_str += """\nTR,DRW,DRF,DELTA,LATT,NPART,NLP
#FE0,RCUT,NLAYER,KH,NPRINT
#MIX, EPS1, EPS2, QVOL
#KF,R0,LIMOL,YZDIVX
#DXXDIV
#NTYPE,NON,NGAUS,E0
#NPLOT,MAXTAU,EQTIM,NCYC"""
#        
#        return(param_str)
