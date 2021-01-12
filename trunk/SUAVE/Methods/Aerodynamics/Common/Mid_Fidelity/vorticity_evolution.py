## @ingroup Methods-Aerodynamics-Common-Mid_Fidelity
# vorticity_evolution.py
# 
# Created:  Dec 2020, R. Erhard
#           

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# package imports 
import numpy as np 
from SUAVE.Core import Data



# This script is used for the viscous Vortex Particle Method (vVPM).
# It consists of the evolution of vorticity over one time step.


def vorticity_evolution():
    
    # sigma: smoothing radius (core size), for simplicity, all core sizes are constant for all filaments at all times
    
    alpha_dot_grad_u = (-1/(4*np.pi)) * sum( (((xp-xq)**2 + 2.5*sigma**2)/((xp-xq)*sigma**2)**2.5)) * np.cross(ap, aq) + 3*(()/())*(np.cross((xp-xq),aq)) * (xp-xq)  )
    nu_del2_omega = (-1/(4*np.pi)) * sum( 105*nu*(sigma**4/(((xp-xq)**2+sigma**2)**4.5) ) * (vol_p*aq - vol_q*ap) )
    
    dalpha_dt = alpha_dot_grad_u + nu_del2_omega
    
    
    return