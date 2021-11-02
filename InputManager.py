#!/usr/bin/python3
import TTCF

import numpy as np

class InputManager:
    """ Overgrown struct to hold user input for the MD routines."""
    def __init__(self):
        # Initialise with default values
        self.tr     = 1.0
        self.drw    = 1.0
        self.drf    = 0.5
        self.delta  = 5.0E-003
        self.latt   = 1
        self.npart  = 500

        self.fe0    = 1.0
        self.rcut   = 2.5
        self.kh     = 1
        self.nprint = 100

        self.mix    = 0.2
        self.eps1   = 1
        self.eps2   = 1
        self.qvol   = 100

        self.kf     = 1.0
        self.r0     = 1.0
        self.limol  = 1
        self.yzdivx = 1.0

        self.dxxdiv = 1.0

        self.ntype  = 1
        self.non    = 0
        self.ngaus  = 4
        self.e0     = 1.0

        self.nplot  = 1
        self.maxtau = 1000
        self.eqtim  = 1000
        self.ncyc   = 1

        # Lastly, set the internal flag determining whether to turn on the NEMD field and zero our
        # timestep value
        self.do_nemd = False
        self.tau = 0
        
        # Now store the variables in their corresponding Fortran blocks
        self.pass_to_fortran()

    def pass_to_fortran(self):
        """ Passes the current value of input parameters to their corresponding Fortran blocks. """

        TTCF.inener.tr      = self.tr      
        TTCF.parm.drw       = self.drw    
        TTCF.parm.drf       = self.drf    
        TTCF.simul.delta    = self.delta  
        TTCF.iparm.latt     = self.latt   
        TTCF.nopart.npart   = self.npart   
                                        
        TTCF.parm.fe0       = self.fe0    
        TTCF.parm.rcut      = self.rcut   
        TTCF.parm.kh        = self.kh     
        TTCF.iparm.nprint   = self.nprint 
                                          
        TTCF.parm.mix       = self.mix    
        TTCF.parm.eps1      = self.eps1   
        TTCF.parm.eps2      = self.eps2   
        TTCF.coord.qvol     = self.qvol   
                                          
        TTCF.parm.kf        = self.kf     
        TTCF.parm.r0        = self.r0     
        TTCF.nopart.limol   = self.limol   
        TTCF.parm.yzdivx    = self.yzdivx 
                                         
        TTCF.parm.dxxdiv    = self.dxxdiv 
        TTCF.iparm.ntype    = self.ntype  
                                          
        TTCF.iparm.non      = self.non    
        TTCF.iparm.ngaus    = self.ngaus  
        TTCF.inener.e0      = self.e0     
                                       
        TTCF.iparm.nplot    = self.nplot  
        TTCF.iparm.maxtau   = self.maxtau 
        TTCF.iparm.eqtim    = self.eqtim  
        TTCF.iparm.ncyc     = self.ncyc   
    
    def format_params(self):
        """ Make a format string which can be either written to file or displayed in Qt."""
        
        # Write the input parameters, making sure to preserve whitespace and newlines
        param_str  = f"""{self.tr:.3f} {self.drw:.3f} {self.drf:.3f} {self.delta:.3f} {self.latt} {self.npart} 25"""
        param_str += f"""\n{self.fe0:.3f} {self.rcut:.3f} 1 {self.kh:.3f} {self.nprint} """
        param_str += f"""\n{self.mix:.3f} {self.eps1:.3f} {self.eps2:.3f} {self.qvol:.3f}"""
        param_str += f"""\n{self.kf:.3f} {self.r0:.3f} {self.limol:.3f} {self.yzdivx:.3f}"""
        param_str += f"""\n{self.dxxdiv:.3f}"""
        param_str += f"""\n{self.ntype} {self.non} {self.ngaus} {self.e0:.3f}"""
        param_str += f"""\n{self.nplot} {self.maxtau} {self.eqtim} {self.ncyc}"""

        # Now write the header
        param_str += """\nTR,DRW,DRF,DELTA,LATT,NPART,NLP
FE0,RCUT,NLAYER,KH,NPRINT
MIX, EPS1, EPS2, QVOL
KF,R0,LIMOL,YZDIVX
DXXDIV
NTYPE,NON,NGAUS,E0
NPLOT,MAXTAU,EQTIM,NCYC"""
        
        return(param_str)

    def read_from_file(self, input_file):
        # Read the input file line by line and assign to internal variables
        pass
        with open(input_file, "r") as ifp:
            # Manually iterate through the data, since each line needs to go into a different set of 
            # variables
            data = next(ifp).split()
            tr,drw,drf,delta,latt,npart,nlp = data

            data = next(ifp).split()
            fe0,rcut,nlayer,kh,nprint = data

            data = next(ifp).split()
            mix, eps1, eps2, qvol = data

            data = next(ifp).split()
            kf,r0,limol,xydivz = data

            dxxdiv = float(next(ifp))

            data = next(ifp).split()
            ntype,non,ngaus,e0 = data

            data = next(ifp).split()
            nplot,maxtau,eqtim,ncyc = data

        
        # now cast the data we read from the file into the proper numpy types
        self.tr, self.drw, self.drf, self.delta = np.array([tr, drw, drf, delta]).astype(np.double)
        self.latt, self.npart = np.array([latt, npart]).astype(np.int)
        self.fe0, self.rcut, self.kh = np.array([fe0, rcut, kh]).astype(np.double)
        self.nprint = np.int(nprint)
        self.mix, self.eps1, self.eps2, self.qvol = np.array([mix, eps1, eps2, qvol]).astype(np.double)

        self.kf, self.r0, self.limol, self.xydivz = np.array([kf,r0,limol,xydivz]).astype(np.double)
        self.dxxdiv = np.double(dxxdiv)
        self.ntype, self.non, self.ngaus = np.array([ntype,non,ngaus]).astype(np.int)
        self.e0 = np.double(e0)
        self.nplot, self.maxtau, self.eqtim, self.ncyc = np.array([nplot,maxtau,eqtim,ncyc]).astype(np.int)
        self.pass_to_fortran()
    
if __name__ == "__main__":
    params = InputManager()
    params.format_params()
