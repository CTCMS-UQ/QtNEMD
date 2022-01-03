#!/usr/bin/python

""" This is a template class to interface with whatever MD backend code we decide to 
    use. It should only include the broadest outline of functionality, which can be 
    overridden to support, e.g. LAMMPS or Debra's MD code."""

import TTCF

class MDInterface:
    def __init__(self):
        self._tr     = 1.0
        self._drw    = 1.0
        self._drf    = 0.5
        self._delta  = 5.0E-003
        self._latt   = 1
        self._npart  = 500

        self._fe0    = 1.0
        self._rcut   = 2.5
        self._kh     = 1
        self._nprint = 100

        self._mix    = 0.2
        self._eps1   = 1
        self._eps2   = 1
        self._qvol   = 100

        self._kf     = 1.0
        self._r0     = 1.0
        self._limol  = 1
        self._yzdivx = 1.0

        self._dxxdiv = 1.0

        self._ntype  = 1
        self._non    = 0
        self._ngaus  = 4
        self._e0     = 1.0

        self._nplot  = 1
        self._maxtau = 1000
        self._eqtim  = 1000
        self._ncyc   = 1

        # Lastly, set the internal flag determining whether to turn on the NEMD field 
        self._do_nemd = False
        self._iflag = 0
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
    def yzdivx(self):
        return(self._yzdivx)
    @yzdivx.setter
    def yzdivx(self, value):
        self._yzdivx = value
        self.reset_and_update_parameters()
    @property
    def reduced_density(self):
        return(self._drf)
    @reduced_density.setter
    def reduced_density(self, value):
        self._drf=value
        self.reset_and_update_parameters()

    @property
    def npart(self):
        return(self._npart)
    @npart.setter
    def npart(self, value):
        self._npart=value
        self.reset_and_update_parameters()

    # WCA parameters. These do not require a restart when they change
    @property
    def eps(self):
        return(self._kf)
    @eps.setter
    def eps(self, value):
        self._kf=value
        self.update_parameters()

    # Misc simulation parameters. Also don't require a restart
    @property
    def fieldstrength(self):
        return(self._fe0)
    @fieldstrength.setter
    def fieldstrength(self, value):
        self._fe0 = value
        self.update_parameters()
    @property
    def temp(self):
        return(TTCF.averg.temp)
    @temp.setter
    def temp(self, value):
        self._tr = value
        self.update_parameters()
    
    # Get the status of whether we're doing nonequilubrium calculations. No setter, since this is handled by toggle_nemd()
    @property
    def do_nemd(self):
        return(self._do_nemd)

    # Atomic coordinates. These should be read-only to external classes
    @property
    def x(self):
        return(TTCF.coord.x)
    @property
    def y(self):
        return(TTCF.coord.y)
    @property
    def z(self):
        return(TTCF.coord.z)
    @property
    def box_bounds(self):
        return(self._box_bounds)
    
    @property
    def vol(self):
        return(TTCF.nopart.npart / TTCF.parm.drf)

###################### Callable methods ###########################
    def setup(self):
        # This function gets called at the start of the program's run, as well as whenever we restart the simulation.
        self.reset_and_update_parameters()
        TTCF.setup()
        TTCF.md(1,0)

        # Now update all of our thermodynamic parameters with values from LAMMPS
        self.get_params_from_TTCF()

    def get_params_from_TTCF(self):
        # Update this class's parameters with the values from LAMMPS. This doesn't change the underlying simulation
        self._temp = TTCF.averg.temp
        self._vol = TTCF.nopart.npart / TTCF.parm.drf
        # And particle positions
        self._x = TTCF.coord.x
        self._y = TTCF.coord.y
        self._z = TTCF.coord.z
        self._box_bounds = (TTCF.parm.cubex, TTCF.parm.cubey, TTCF.parm.cubez)

        # Finally, grab the g(2) radial distribution function
        #self._g2_arrays = self._lmp.numpy.extract_compute("g2", lammps.LMP_STYLE_GLOBAL, lammps.LMP_TYPE_ARRAY)

    # g(2) radial distribution function. This returns into two arrays for r and g(2)(r), indexed in a dict
    def g2_compute(self):
      # TODO
      pass

    def reset_and_update_parameters(self):
        # Reset the simulation and initialise the parameters
        TTCF.inener.tr      = self._tr      
        TTCF.parm.drw       = self._drw    
        TTCF.parm.drf       = self._drf    
        TTCF.simul.delta    = self._delta  
        TTCF.iparm.latt     = self._latt   
        TTCF.nopart.npart   = self._npart   
                                        
        TTCF.parm.fe0       = self._fe0    
        TTCF.parm.rcut      = self._rcut   
        TTCF.parm.kh        = self._kh     
        TTCF.iparm.nprint   = self._nprint 
                                          
        TTCF.parm.mix       = self._mix    
        TTCF.parm.eps1      = self._eps1   
        TTCF.parm.eps2      = self._eps2   
        TTCF.coord.qvol     = self._qvol   
                                          
        TTCF.parm.kf        = self._kf     
        TTCF.parm.r0        = self._r0     
        TTCF.nopart.limol   = self._limol   
        TTCF.parm.yzdivx    = self._yzdivx 
                                         
        TTCF.parm.dxxdiv    = self._dxxdiv 
        TTCF.iparm.ntype    = self._ntype  
                                          
        TTCF.iparm.non      = self._non    
        TTCF.iparm.ngaus    = self._ngaus  
        TTCF.inener.e0      = self._e0     
                                       
        TTCF.iparm.nplot    = self._nplot  
        TTCF.iparm.maxtau   = self._maxtau 
        TTCF.iparm.eqtim    = self._eqtim  
        TTCF.iparm.ncyc     = self._ncyc
        TTCF.setup()

    def update_parameters(self):
        # Only update the parameters which don't require a reset
        TTCF.inener.tr      = self._tr      
        TTCF.parm.drw       = self._drw    
        TTCF.parm.drf       = self._drf    
        TTCF.simul.delta    = self._delta  
        TTCF.iparm.latt     = self._latt   
        TTCF.nopart.npart   = self._npart   
                                        
        TTCF.parm.fe0       = self._fe0    
        TTCF.parm.rcut      = self._rcut   
        TTCF.parm.kh        = self._kh     
        TTCF.iparm.nprint   = self._nprint 
                                          
        TTCF.parm.mix       = self._mix    
        TTCF.parm.eps1      = self._eps1   
        TTCF.parm.eps2      = self._eps2   
        TTCF.coord.qvol     = self._qvol   
                                          
        TTCF.parm.kf        = self._kf     
        TTCF.parm.r0        = self._r0     
        TTCF.nopart.limol   = self._limol   
        TTCF.parm.yzdivx    = self._yzdivx 
                                         
        TTCF.parm.dxxdiv    = self._dxxdiv 
        TTCF.iparm.ntype    = self._ntype  
                                          
        TTCF.iparm.non      = self._non    
        TTCF.iparm.ngaus    = self._ngaus  
        TTCF.inener.e0      = self._e0     
                                       
        TTCF.iparm.nplot    = self._nplot  
        TTCF.iparm.maxtau   = self._maxtau 
        TTCF.iparm.eqtim    = self._eqtim  
        TTCF.iparm.ncyc     = self._ncyc      

    def toggle_nemd(self):
        if(self._do_nemd):
          self._do_nemd = False
          self._iflag = 0
          self.update_parameters()
        else:
          self._do_nemd = True
          self._iflag = 1
          self.update_parameters()

    def format_params(self):
        """ Make a format string for the input which can be either written to file or displayed in Qt."""
        
        # Write the input parameters, making sure to preserve whitespace and newlines
        # Write the input parameters, making sure to preserve whitespace and newlines
        param_str  =   f"""{self._tr:.3f}  {self._drw:.3f}  {self._drf:.3f}   {self._delta:.3f} {self._latt} {self._npart} 25"""
        param_str += f"""\n{self._fe0:.3f} {self._rcut:.3f} 1 {self._kh:.3f}  {self._nprint} """
        param_str += f"""\n{self._mix:.3f} {self._eps1:.3f} {self._eps2:.3f}  {self._qvol:.3f}"""
        param_str += f"""\n{self._kf:.3f}  {self._r0:.3f}   {self._limol:.3f} {self._yzdivx:.3f}"""
        param_str += f"""\n{self._dxxdiv:.3f}"""
        param_str += f"""\n{self._ntype} {self._non}    {self._ngaus} {self._e0:.3f}"""
        param_str += f"""\n{self._nplot} {self._maxtau} {self._eqtim} {self._ncyc}"""

        # Now write the header
        param_str += """\nTR,DRW,DRF,DELTA,LATT,NPART,NLP
FE0,RCUT,NLAYER,KH,NPRINT
MIX, EPS1, EPS2, QVOL
KF,R0,LIMOL,YZDIVX
DXXDIV
NTYPE,NON,NGAUS,E0
NPLOT,MAXTAU,EQTIM,NCYC"""
        return(param_str)

    def run(self, nsteps):
        # First, advance the simulation
        if(nsteps >= 0):
            TTCF.md(nsteps, self._iflag)
        else:
            raise(ValueError)
        
        # Now update all of our thermodynamic parameters with values from LAMMPS
        self.get_params_from_TTCF()