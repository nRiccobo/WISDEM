import copy
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import PchipInterpolator
from openmdao.api import ExplicitComponent, Group
from wisdem.ccblade.ccblade_component import CCBladeLoads, AeroHubLoads
import wisdem.ccblade._bem as _bem
import wisdem.commonse.utilities as util
from wisdem.commonse.akima import Akima
from wisdem.commonse import gravity
from wisdem.commonse.csystem import DirectionVector
from wisdem.rotorse import RPM2RS, RS2RPM
import wisdem.pyframe3dd.pyframe3dd as pyframe3dd

class BladeCurvature(ExplicitComponent):
    # OpenMDAO component that computes the 3D curvature of the blade
    def initialize(self):
        self.options.declare('analysis_options')

    def setup(self):
        n_span = self.options['analysis_options']['blade']['n_span']

        # Inputs
        self.add_input('r',         val=np.zeros(n_span), units='m',      desc='location in blade z-coordinate')
        self.add_input('precurve',  val=np.zeros(n_span), units='m',      desc='location in blade x-coordinate')
        self.add_input('presweep',  val=np.zeros(n_span), units='m',      desc='location in blade y-coordinate')
        self.add_input('precone',   val=0.0,              units='deg',    desc='precone angle')

        # Outputs
        self.add_output('3d_curv',  val=np.zeros(n_span),units='deg',    desc='total cone angle from precone and curvature')
        self.add_output('x_az',     val=np.zeros(n_span), units='m',      desc='location of blade in azimuth x-coordinate system')
        self.add_output('y_az',     val=np.zeros(n_span), units='m',      desc='location of blade in azimuth y-coordinate system')
        self.add_output('z_az',     val=np.zeros(n_span), units='m',      desc='location of blade in azimuth z-coordinate system')
        self.add_output('s',        val=np.zeros(n_span), units='m',      desc='cumulative path length along blade')

    def compute(self, inputs, outputs):

        r = inputs['r']
        precurve = inputs['precurve']
        presweep = inputs['presweep']
        precone = inputs['precone']

        n = len(r)
        dx_dx = np.eye(3*n)

        x_az, x_azd, y_az, y_azd, z_az, z_azd, cone, coned, s, sd = _bem.definecurvature_dv2(r, dx_dx[:, :n],
                                                                                             precurve, dx_dx[:, n:2*n],
                                                                                             presweep, dx_dx[:, 2*n:],
                                                                                             0.0, np.zeros(3*n))

        totalCone = precone + np.degrees(cone)
        s = r[0] + s

        outputs['3d_curv'] = totalCone
        outputs['x_az'] = x_az
        outputs['y_az'] = y_az
        outputs['z_az'] = z_az
        outputs['s'] = s

class TotalLoads(ExplicitComponent):
    # OpenMDAO component that takes as input the rotor configuration (tilt, cone), the blade twist and mass distributions, and the blade aerodynamic loading, and computes the total loading including gravity and centrifugal forces
    def initialize(self):
        self.options.declare('analysis_options')

    def setup(self):
        n_span = self.options['analysis_options']['blade']['n_span']

        # Inputs
        self.add_input('r',                 val=np.zeros(n_span),   units='m',      desc='radial positions along blade going toward tip')
        self.add_input('aeroloads_Px',      val=np.zeros(n_span),   units='N/m',    desc='distributed loads in blade-aligned x-direction')
        self.add_input('aeroloads_Py',      val=np.zeros(n_span),   units='N/m',    desc='distributed loads in blade-aligned y-direction')
        self.add_input('aeroloads_Pz',      val=np.zeros(n_span),   units='N/m',    desc='distributed loads in blade-aligned z-direction')
        self.add_input('aeroloads_Omega',   val=0.0,                units='rpm',    desc='rotor rotation speed')
        self.add_input('aeroloads_pitch',   val=0.0,                units='deg',    desc='pitch angle')
        self.add_input('aeroloads_azimuth', val=0.0,                units='deg',    desc='azimuthal angle')
        self.add_input('theta',             val=np.zeros(n_span),   units='deg',    desc='structural twist')
        self.add_input('tilt',              val=0.0,                units='deg',    desc='tilt angle')
        self.add_input('3d_curv',           val=np.zeros(n_span),   units='deg',    desc='total cone angle from precone and curvature')
        self.add_input('z_az',              val=np.zeros(n_span),   units='m',      desc='location of blade in azimuth z-coordinate system')
        self.add_input('rhoA',              val=np.zeros(n_span),   units='kg/m',   desc='mass per unit length')
        self.add_input('dynamicFactor',     val=1.0,                                desc='a dynamic amplification factor to adjust the static deflection calculation')

        # Outputs
        self.add_output('Px_af', val=np.zeros(n_span), desc='total distributed loads in airfoil x-direction')
        self.add_output('Py_af', val=np.zeros(n_span), desc='total distributed loads in airfoil y-direction')
        self.add_output('Pz_af', val=np.zeros(n_span), desc='total distributed loads in airfoil z-direction')


    def compute(self, inputs, outputs):

        dynamicFactor = inputs['dynamicFactor']
        r             = inputs['r']
        theta         = inputs['theta']
        tilt          = inputs['tilt']
        totalCone     = inputs['3d_curv']
        z_az          = inputs['z_az']
        rhoA          = inputs['rhoA']

        # keep all in blade c.s. then rotate all at end

        # --- aero loads ---
        P_a = DirectionVector(0, 0, 0)
        P_a.x, P_a.y, P_a.z = inputs['aeroloads_Px'], inputs['aeroloads_Py'], inputs['aeroloads_Pz']

        # --- weight loads ---
        # yaw c.s.
        weight = DirectionVector(0.0, 0.0, -rhoA*gravity)
        P_w = weight.yawToHub(tilt).hubToAzimuth(inputs['aeroloads_azimuth']).azimuthToBlade(totalCone)

        # --- centrifugal loads ---
        # azimuthal c.s.
        Omega = inputs['aeroloads_Omega']*RPM2RS
        load = DirectionVector(0.0, 0.0, rhoA*Omega**2*z_az)
        P_c = load.azimuthToBlade(totalCone)

        # --- total loads ---
        P = P_a + P_w + P_c

        # rotate to airfoil c.s.
        P = P.bladeToAirfoil(theta + inputs['aeroloads_pitch'])

        outputs['Px_af'] = dynamicFactor * P.x
        outputs['Py_af'] = dynamicFactor * P.y
        outputs['Pz_af'] = dynamicFactor * P.z




class RunFrame3DD(ExplicitComponent):
    def initialize(self):
        self.options.declare('analysis_options')
        self.options.declare('pbeam',default=False) # Recover old pbeam c.s. and accuracy

    def setup(self):
        blade_init_options = self.options['analysis_options']['blade']
        self.n_span = n_span = blade_init_options['n_span']
        self.n_freq = n_freq = blade_init_options['n_freq']

        # Locations of airfoils in global c.s.
        self.add_input('x_az',     val=np.zeros(n_span), units='m',      desc='location of blade in azimuth x-coordinate system (from LE to TE)')
        self.add_input('y_az',     val=np.zeros(n_span), units='m',      desc='location of blade in azimuth y-coordinate system (from PS to SS)')
        self.add_input('z_az',     val=np.zeros(n_span), units='m',      desc='location of blade in azimuth z-coordinate system (from root to tip)')
        self.add_input('theta',    val=np.zeros(n_span),   units='deg',    desc='structural twist')
        
        # all inputs/outputs in airfoil coordinate system
        self.add_input('Px_af', val=np.zeros(n_span), desc='distributed load (force per unit length) in airfoil x-direction')
        self.add_input('Py_af', val=np.zeros(n_span), desc='distributed load (force per unit length) in airfoil y-direction')
        self.add_input('Pz_af', val=np.zeros(n_span), desc='distributed load (force per unit length) in airfoil z-direction')

        self.add_input('xu_strain_spar',    val=np.zeros(n_span), desc='x-position of midpoint of spar cap on upper surface for strain calculation')
        self.add_input('xl_strain_spar',    val=np.zeros(n_span), desc='x-position of midpoint of spar cap on lower surface for strain calculation')
        self.add_input('yu_strain_spar',    val=np.zeros(n_span), desc='y-position of midpoint of spar cap on upper surface for strain calculation')
        self.add_input('yl_strain_spar',    val=np.zeros(n_span), desc='y-position of midpoint of spar cap on lower surface for strain calculation')
        self.add_input('xu_strain_te',      val=np.zeros(n_span), desc='x-position of midpoint of trailing-edge panel on upper surface for strain calculation')
        self.add_input('xl_strain_te',      val=np.zeros(n_span), desc='x-position of midpoint of trailing-edge panel on lower surface for strain calculation')
        self.add_input('yu_strain_te',      val=np.zeros(n_span), desc='y-position of midpoint of trailing-edge panel on upper surface for strain calculation')
        self.add_input('yl_strain_te',      val=np.zeros(n_span), desc='y-position of midpoint of trailing-edge panel on lower surface for strain calculation')

        self.add_input('r',     val=np.zeros(n_span), units='m',        desc='locations of properties along beam')
        self.add_input('A',     val=np.zeros(n_span), units='m**2',     desc='airfoil cross section material area')
        self.add_input('EA',    val=np.zeros(n_span), units='N',        desc='axial stiffness')
        self.add_input('EIxx',  val=np.zeros(n_span), units='N*m**2',   desc='edgewise stiffness (bending about :ref:`x-axis of airfoil aligned coordinate system <blade_airfoil_coord>`)')
        self.add_input('EIyy',  val=np.zeros(n_span), units='N*m**2',   desc='flapwise stiffness (bending about y-axis of airfoil aligned coordinate system)')
        self.add_input('EIxy',  val=np.zeros(n_span), units='N*m**2',   desc='coupled flap-edge stiffness')
        self.add_input('GJ',    val=np.zeros(n_span), units='N*m**2',   desc='torsional stiffness (about axial z-direction of airfoil aligned coordinate system)')
        self.add_input('rhoA',  val=np.zeros(n_span), units='kg/m',     desc='mass per unit length')
        self.add_input('rhoJ',  val=np.zeros(n_span), units='kg*m',     desc='polar mass moment of inertia per unit length')
        self.add_input('x_ec',  val=np.zeros(n_span), units='m',        desc='x-distance to elastic center from point about which above structural properties are computed (airfoil aligned coordinate system)')
        self.add_input('y_ec',  val=np.zeros(n_span), units='m', desc='y-distance to elastic center from point about which above structural properties are computed')

        # outputs
        n_freq2 = int(n_freq/2)
        self.add_output('root_F', np.zeros(3), units='N',   desc='Blade root forces in blade c.s.')
        self.add_output('root_M', np.zeros(3), units='N*m', desc='Blade root moment in blade c.s.')
        self.add_output('flap_mode_shapes', np.zeros((n_freq2,5)), desc='6-degree polynomial coefficients of mode shapes in the flap direction (x^2..x^6, no linear or constant term)')
        self.add_output('edge_mode_shapes', np.zeros((n_freq2,5)), desc='6-degree polynomial coefficients of mode shapes in the edge direction (x^2..x^6, no linear or constant term)')
        self.add_output('all_mode_shapes', np.zeros((n_freq,5)), desc='6-degree polynomial coefficients of mode shapes in the edge direction (x^2..x^6, no linear or constant term)')
        self.add_output('flap_mode_freqs', np.zeros(n_freq2), units='Hz', desc='Frequencies associated with mode shapes in the flap direction')
        self.add_output('edge_mode_freqs', np.zeros(n_freq2), units='Hz', desc='Frequencies associated with mode shapes in the edge direction')
        self.add_output('freqs',            val=np.zeros(n_freq),  units='Hz', desc='ration of 2nd and 1st natural frequencies, should be ratio of edgewise to flapwise')
        self.add_output('freq_distance',    val=0.0,              desc='ration of 2nd and 1st natural frequencies, should be ratio of edgewise to flapwise')
        self.add_output('dx',               val=np.zeros(n_span), units='m', desc='deflection of blade section in airfoil x-direction')
        self.add_output('dy',               val=np.zeros(n_span), units='m', desc='deflection of blade section in airfoil y-direction')
        self.add_output('dz',               val=np.zeros(n_span), units='m', desc='deflection of blade section in airfoil z-direction')
        self.add_output('strainU_spar',     val=np.zeros(n_span), desc='strain in spar cap on upper surface at location xu,yu_strain with loads P_strain')
        self.add_output('strainL_spar',     val=np.zeros(n_span), desc='strain in spar cap on lower surface at location xl,yl_strain with loads P_strain')
        self.add_output('strainU_te',       val=np.zeros(n_span), desc='strain in trailing-edge panels on upper surface at location xu,yu_te with loads P_te')
        self.add_output('strainL_te',       val=np.zeros(n_span), desc='strain in trailing-edge panels on lower surface at location xl,yl_te with loads P_te')

        
    def compute(self, inputs, outputs):

        # Unpack inputs
        r     = inputs['r']
        x_az  = inputs['x_az']
        y_az  = inputs['y_az']
        z_az  = inputs['z_az']
        theta = inputs['theta']
        x_ec  = inputs['x_ec']
        y_ec  = inputs['y_ec']
        A     = inputs['A']
        rhoA  = inputs['rhoA']
        rhoJ  = inputs['rhoJ']
        GJ    = inputs['GJ']
        EA    = inputs['EA']
        EIxx  = inputs['EIxx']
        EIyy  = inputs['EIyy']
        EIxy  = inputs['EIxy']
        Px_af = inputs['Px_af']
        Py_af = inputs['Py_af']
        Pz_af = inputs['Pz_af']
        xu_strain_spar = inputs['xu_strain_spar']
        xl_strain_spar = inputs['xl_strain_spar']
        yu_strain_spar = inputs['yu_strain_spar']
        yl_strain_spar = inputs['yl_strain_spar']
        xu_strain_te   = inputs['xu_strain_te']
        xl_strain_te   = inputs['xl_strain_te']
        yu_strain_te   = inputs['yu_strain_te']
        yl_strain_te   = inputs['yl_strain_te']
        #np.savez('nrel5mw_test.npz',r=r,x_az=x_az,y_az=y_az,z_az=z_az,theta=theta,x_ec=x_ec,y_ec=y_ec,A=A,rhoA=rhoA,rhoJ=rhoJ,GJ=GJ,EA=EA,EIxx=EIxx,EIyy=EIyy,EIxy=EIxy,Px_af=Px_af,Py_af=Py_af,Pz_af=Pz_af,xu_strain_spar=xu_strain_spar,xl_strain_spar=xl_strain_spar,yu_strain_spar=yu_strain_spar,yl_strain_spar=yl_strain_spar,xu_strain_te=xu_strain_te,xl_strain_te=xl_strain_te,yu_strain_te=yu_strain_te,yl_strain_te=yl_strain_te)
        
        # Determine principal C.S. (with swap of x, y for profile c.s.)
        # Can get to Hansen's c.s. from Precomp's c.s. by rotating around z -90 deg, then y by 180 (swap x-y)
        EIxx_cs , EIyy_cs = EIyy.copy() , EIxx.copy()
        x_ec_cs , y_ec_cs = y_ec.copy() , x_ec.copy()
        EIxy_cs = EIxy.copy()

        # translate to elastic center
        EIxx_cs -= y_ec_cs**2 * EA
        EIyy_cs -= x_ec_cs**2 * EA
        EIxy_cs -= x_ec_cs * y_ec_cs * EA

        # get rotation angle
        alpha = 0.5*np.arctan2(2*EIxy_cs, EIyy_cs-EIxx_cs)

        # get moments and positions in principal axes
        EI11 = EIxx_cs - EIxy_cs*np.tan(alpha)
        EI22 = EIyy_cs + EIxy_cs*np.tan(alpha)
        ca   = np.cos(alpha)
        sa   = np.sin(alpha)
        def rotate(x,y):
            x2 =  x*ca + y*sa
            y2 = -x*sa + y*ca
            return x2, y2
        
        # Now store alpha for later use in degrees
        alpha = np.rad2deg(alpha)
        
        # Frame3dd call
        # ------- node data ----------------
        n     = len(z_az)
        rad   = np.zeros(n) # 'radius' of rigidity at node- set to zero
        inode = 1 + np.arange(n) # Node numbers (1-based indexing)
        if self.options['pbeam']:
            nodes = pyframe3dd.NodeData(inode, np.zeros(n), np.zeros(n), r, rad)
            L     = np.diff(r)
        else:
            nodes = pyframe3dd.NodeData(inode, x_az, y_az, z_az, rad)
            L     = np.sqrt(np.diff(x_az)**2 + np.diff(y_az)**2 + np.diff(z_az)**2)
        # -----------------------------------

        # ------ reaction data ------------
        # Pinned at root
        rnode = np.array([1])
        rigid = np.array([1e16])
        reactions = pyframe3dd.ReactionData(rnode, rigid, rigid, rigid, rigid, rigid, rigid, float(rigid))
        # -----------------------------------

        # ------ frame element data ------------
        elem = np.arange(1, n) # Element Numbers
        N1   = np.arange(1, n) # Node number start
        N2   = np.arange(2, n+1) # Node number finish

        E   = EA   / A
        rho = rhoA / A
        J   = rhoJ / rho
        G   = GJ   / J
        if self.options['pbeam']:
            # Use airfoil c.s.
            Ix  = EIyy / E
            Iy  = EIxx / E 
        else:
            # Will further rotate to principle axes
            Ix  = EI11 / E
            Iy  = EI22 / E

        # Have to convert nodal values to find average at center of element
        Abar,_   = util.nodal2sectional(A)
        Ebar,_   = util.nodal2sectional(E)
        rhobar,_ = util.nodal2sectional(rho)
        Jbar,_   = util.nodal2sectional(J)
        Gbar,_   = util.nodal2sectional(G)
        Ixbar,_  = util.nodal2sectional(Ix)
        Iybar,_  = util.nodal2sectional(Iy)

        # Angle of element principal axes relative to global coordinate system
        if self.options['pbeam']:
            # Work in airfoil c.s. for both global and local c.s.
            roll = np.zeros(n-1)
        else:
            # Global c.s. is blade, local element c.s. is airfoil (twist + principle rotation)
            roll,_ = util.nodal2sectional(theta + alpha)
            
        Asx = Asy = 1e-6*np.ones(elem.shape) # Unused when shear=False
        elements   = pyframe3dd.ElementData(elem, N1, N2, Abar, Asx, Asy, Jbar, Ixbar, Iybar, Ebar, Gbar, roll, rhobar)
        # -----------------------------------

        # ------ options ------------
        shear   = False # If not false, have to compute Asx or Asy
        geom    = (not self.options['pbeam']) # Must be true for spin-stiffening
        dx      = -1 # Don't need stress changes within element for now
        options = pyframe3dd.Options(shear, geom, dx)
        # -----------------------------------

        # initialize frame3dd object
        blade = pyframe3dd.Frame(nodes, reactions, elements, options)

        # ------- enable dynamic analysis ----------
        Mmethod = 1 # 1= Subspace-Jacobi iteration, 2= Stodola (matrix iteration) method
        lump    = 0 # 0= consistent mass matrix, 1= lumped mass matrix
        tol     = 1e-9 # frequency convergence tolerance
        shift   = 0.0 # frequency shift-factor for rigid body modes, make 0 for pos.def. [K]
        blade.enableDynamics(self.n_freq, Mmethod, lump, tol, shift)
        # ----------------------------

        # ------ load case 1, blade 1 ------------
        # trapezoidally distributed loads- already has gravity, centrifugal, aero, etc.
        gx = gy = gz = 0.0
        load = pyframe3dd.StaticLoadCase(gx, gy, gz)

        if not self.options['pbeam']:
            # Have to further move the loads into principle directions
            P = DirectionVector(Px_af, Py_af, Pz_af).bladeToAirfoil(alpha)
            Px_af = P.x
            Py_af = P.y
            Pz_af = P.z
            
        Px, Py, Pz = Pz_af, Py_af, -Px_af # switch to local c.s.
        xx1 = xy1 = xz1 = np.zeros(n-1)
        xx2 = xy2 = xz2 = L - 1e-6  # subtract small number b.c. of precision
        wx1 = Px[:-1]
        wx2 = Px[1:]
        wy1 = Py[:-1]
        wy2 = Py[1:]
        wz1 = Pz[:-1]
        wz2 = Pz[1:]
        load.changeTrapezoidalLoads(elem, xx1, xx2, wx1, wx2, xy1, xy2, wy1, wy2, xz1, xz2, wz1, wz2)
        blade.addLoadCase(load)
        
        # Debugging
        #blade.write('blade.3dd')

        # run the analysis
        displacements, forces, reactions, internalForces, mass, modal = blade.run()

        # For now, just 1 load case and blade
        iCase = 0

        # Displacements in global (blade) c.s.
        dx = displacements.dx[iCase,:]
        dy = displacements.dy[iCase,:]
        dz = displacements.dz[iCase,:]
        
        # Mode shapes and frequencies
        freq_x, freq_y, mshapes_x, mshapes_y = util.get_mode_shapes(r, modal.freq, modal.xdsp, modal.ydsp, modal.zdsp, modal.xmpf, modal.ympf, modal.zmpf)

        # shear and bending, one per element (convert from local to global c.s.)
        Fz = np.r_[-forces.Nx[iCase,0],  forces.Nx[iCase, 1::2]]
        Vy = np.r_[-forces.Vy[iCase,0],  forces.Vy[iCase, 1::2]]
        Vx = np.r_[ forces.Vz[iCase,0], -forces.Vz[iCase, 1::2]]

        Tz = np.r_[-forces.Txx[iCase,0],  forces.Txx[iCase, 1::2]]
        My = np.r_[-forces.Myy[iCase,0],  forces.Myy[iCase, 1::2]]
        Mx = np.r_[ forces.Mzz[iCase,0], -forces.Mzz[iCase, 1::2]]

        def strain(xu, yu, xl, yl):
            # use profile c.s. to use Hansen's notation
            xuu, yuu = yu, xu
            xll, yll = yl, xl

            # convert to principal axes, unless already there
            if self.options['pbeam']:
                M1,M2 = rotate(My, Mx)
            else:
                M1,M2 = My,Mx

            # compute strain
            x,y = rotate(xuu, yuu)
            strainU = -(M1/EI11*y - M2/EI22*x + Fz/EA)  # negative sign because Hansen c3 is opposite of Precomp z

            x,y = rotate(xll, yll)
            strainL = -(M1/EI11*y - M2/EI22*x + Fz/EA)

            return strainU, strainL
        
        # ----- strain -----
        strainU_spar, strainL_spar = strain(xu_strain_spar, yu_strain_spar, xl_strain_spar, yl_strain_spar)
        strainU_te,   strainL_te   = strain(xu_strain_te,  yu_strain_te,    xl_strain_te,   yl_strain_te)

        # Store outputs
        outputs['root_F'] = -1.0 * np.array([reactions.Fx.sum(), reactions.Fy.sum(), reactions.Fz.sum()])
        outputs['root_M'] = -1.0 * np.array([reactions.Mxx.sum(), reactions.Myy.sum(), reactions.Mzz.sum()])
        outputs['freqs'] = modal.freq
        outputs['edge_mode_shapes'] = mshapes_y
        outputs['flap_mode_shapes'] = mshapes_x
        # Dense numpy command that interleaves and alternates flap and edge modes
        outputs['all_mode_shapes'] = np.c_[mshapes_x, mshapes_y].flatten().reshape((self.n_freq,5))
        outputs['edge_mode_freqs']  = freq_y
        outputs['flap_mode_freqs']  = freq_x
        outputs['freq_distance']    = freq_y[0] / freq_x[0]
        outputs['dx'] = dx
        outputs['dy'] = dy
        outputs['dz'] = dz
        outputs['strainU_spar'] = strainU_spar
        outputs['strainL_spar'] = strainL_spar
        outputs['strainU_te'] = strainU_te
        outputs['strainL_te'] = strainL_te

        
class TipDeflection(ExplicitComponent):
    # OpenMDAO component that computes the blade deflection at tip in yaw x-direction
    def setup(self):
        # Inputs
        self.add_input('dx_tip',        val=0.0,    units='m',      desc='deflection at tip in blade x-direction')
        self.add_input('dy_tip',        val=0.0,    units='m',      desc='deflection at tip in blade y-direction')
        self.add_input('dz_tip',        val=0.0,    units='m',      desc='deflection at tip in blade z-direction')
        #self.add_input('theta_tip',     val=0.0,    units='deg',    desc='twist at tip section')
        self.add_input('pitch_load',    val=0.0,    units='deg',    desc='blade pitch angle')
        self.add_input('tilt',          val=0.0,    units='deg',    desc='tilt angle')
        self.add_input('3d_curv_tip',   val=0.0,    units='deg',    desc='total coning angle including precone and curvature')
        self.add_input('dynamicFactor', val=1.0,                    desc='a dynamic amplification factor to adjust the static deflection calculation') #)
        # Outputs
        self.add_output('tip_deflection', val=0.0,  units='m',      desc='deflection at tip in yaw x-direction')

    def compute(self, inputs, outputs):

        dx            = inputs['dx_tip']
        dy            = inputs['dy_tip']
        dz            = inputs['dz_tip']
        pitch         = inputs['pitch_load'] #+ inputs['theta_tip']
        azimuth       = 180.0 # The blade is assumed in front of the tower, although the loading may correspond to another azimuthal position
        tilt          = inputs['tilt']
        totalConeTip  = inputs['3d_curv_tip']
        dynamicFactor = inputs['dynamicFactor']

        dr = DirectionVector(dx, dy, dz)

        delta = dr.airfoilToBlade(pitch).bladeToAzimuth(totalConeTip).azimuthToHub(azimuth).hubToYaw(tilt)

        outputs['tip_deflection'] = dynamicFactor * delta.x

        
class DesignConstraints(ExplicitComponent):
    # OpenMDAO component that formulates constraints on user-defined maximum strains, frequencies   
    def initialize(self):
        self.options.declare('analysis_options')
        self.options.declare('opt_options')

    def setup(self):
        blade_init_options = self.options['analysis_options']['blade']
        self.n_span = n_span = blade_init_options['n_span']
        self.n_freq = n_freq = blade_init_options['n_freq']
        n_freq2 = int(n_freq/2)
        self.opt_options   = opt_options   = self.options['opt_options']
        self.n_opt_spar_cap_ss = n_opt_spar_cap_ss = opt_options['optimization_variables']['blade']['structure']['spar_cap_ss']['n_opt']
        self.n_opt_spar_cap_ps = n_opt_spar_cap_ps = opt_options['optimization_variables']['blade']['structure']['spar_cap_ps']['n_opt']
        # Inputs strains
        self.add_input('strainU_spar',     val=np.zeros(n_span), desc='strain in spar cap on upper surface at location xu,yu_strain with loads P_strain')
        self.add_input('strainL_spar',     val=np.zeros(n_span), desc='strain in spar cap on lower surface at location xl,yl_strain with loads P_strain')

        self.add_input('min_strainU_spar', val=0.0, desc='minimum strain in spar cap suction side')
        self.add_input('max_strainU_spar', val=0.0, desc='minimum strain in spar cap pressure side')
        self.add_input('min_strainL_spar', val=0.0, desc='maximum strain in spar cap suction side')
        self.add_input('max_strainL_spar', val=0.0, desc='maximum strain in spar cap pressure side')
        
        self.add_input('s',                     val=np.zeros(n_span),       desc='1D array of the non-dimensional spanwise grid defined along blade axis (0-blade root, 1-blade tip)')
        self.add_input('s_opt_spar_cap_ss',         val=np.zeros(n_opt_spar_cap_ss),desc='1D array of the non-dimensional spanwise grid defined along blade axis to optimize the blade spar cap suction side')
        self.add_input('s_opt_spar_cap_ps',         val=np.zeros(n_opt_spar_cap_ss),desc='1D array of the non-dimensional spanwise grid defined along blade axis to optimize the blade spar cap suction side')

        # Input frequencies
        self.add_input('rated_Omega', val=0.0,                units='rpm', desc='rotor rotation speed at rated')
        self.add_input('delta_f',     val=1.1,                             desc='minimum frequency margin')
        self.add_input('flap_mode_freqs', np.zeros(n_freq2), units='Hz', desc='Frequencies associated with mode shapes in the flap direction')
        self.add_input('edge_mode_freqs', np.zeros(n_freq2), units='Hz', desc='Frequencies associated with mode shapes in the edge direction')

        # Outputs
        # self.add_output('constr_min_strainU_spar',     val=np.zeros(n_opt_spar_cap_ss), desc='constraint for minimum strain in spar cap suction side')
        self.add_output('constr_max_strainU_spar',     val=np.zeros(n_opt_spar_cap_ss), desc='constraint for maximum strain in spar cap suction side')
        # self.add_output('constr_min_strainL_spar',     val=np.zeros(n_opt_spar_cap_ps), desc='constraint for minimum strain in spar cap pressure side')
        self.add_output('constr_max_strainL_spar',     val=np.zeros(n_opt_spar_cap_ps), desc='constraint for maximum strain in spar cap pressure side')
        self.add_output('constr_flap_f_margin',      val=np.zeros(n_freq2),             desc='constraint on flap blade frequency such that ratio of 3P/f is above or below delta with constraint <= 0')
        self.add_output('constr_edge_f_margin',      val=np.zeros(n_freq2),             desc='constraint on edge blade frequency such that ratio of 3P/f is above or below delta with constraint <= 0')

    def compute(self, inputs, outputs):
        
        # Constraints on blade strains
        s               = inputs['s']
        s_opt_spar_cap_ss   = inputs['s_opt_spar_cap_ss']
        s_opt_spar_cap_ps   = inputs['s_opt_spar_cap_ps']
        
        strainU_spar     = inputs['strainU_spar']
        strainL_spar     = inputs['strainL_spar']
        # min_strainU_spar = inputs['min_strainU_spar']
        if inputs['max_strainU_spar'] == np.zeros_like(inputs['max_strainU_spar']):
            max_strainU_spar =  np.ones_like(inputs['max_strainU_spar'])
        else:
            max_strainU_spar = inputs['max_strainU_spar']
        # min_strainL_spar = inputs['min_strainL_spar']
        if inputs['max_strainL_spar'] == np.zeros_like(inputs['max_strainL_spar']):
            max_strainL_spar =  np.ones_like(inputs['max_strainL_spar'])
        else:
            max_strainL_spar = inputs['max_strainL_spar']

        # outputs['constr_min_strainU_spar'] = abs(np.interp(s_opt_spar_cap_ss, s, strainU_spar)) / abs(min_strainU_spar)
        outputs['constr_max_strainU_spar'] = abs(np.interp(s_opt_spar_cap_ss, s, strainU_spar)) / max_strainU_spar
        # outputs['constr_min_strainL_spar'] = abs(np.interp(s_opt_spar_cap_ps, s, strainL_spar)) / abs(min_strainL_spar)
        outputs['constr_max_strainL_spar'] = abs(np.interp(s_opt_spar_cap_ps, s, strainL_spar)) / max_strainL_spar

        # Constraints on blade frequencies
        threeP = 3. * inputs['rated_Omega'] / 60. # TODO: CHange this to nBlades
        flap_f = inputs['flap_mode_freqs']
        edge_f = inputs['edge_mode_freqs']
        delta  = inputs['delta_f']
        outputs['constr_flap_f_margin'] = np.array( [min([threeP-(2-delta)*f, delta*f-threeP]) for f in flap_f] ).flatten()
        outputs['constr_edge_f_margin'] = np.array( [min([threeP-(2-delta)*f, delta*f-threeP]) for f in edge_f] ).flatten()
        
        
class RotorLoadsDeflStrains(Group):
    # OpenMDAO group to compute the blade elastic properties, deflections, and loading
    def initialize(self):
        self.options.declare('analysis_options')
        self.options.declare('opt_options')
    def setup(self):
        analysis_options = self.options['analysis_options']
        opt_options     = self.options['opt_options']

        # Load blade with rated conditions and compute aerodynamic forces
        promoteListAeroLoads =  ['r', 'theta', 'chord', 'Rtip', 'Rhub', 'hub_height', 'precone', 'tilt', 'airfoils_aoa', 'airfoils_Re', 'airfoils_cl', 'airfoils_cd', 'airfoils_cm', 'nBlades', 'rho', 'mu', 'Omega_load','pitch_load']
        # self.add_subsystem('aero_rated',        CCBladeLoads(analysis_options = analysis_options), promotes=promoteListAeroLoads)
        self.add_subsystem('aero_gust',         CCBladeLoads(analysis_options = analysis_options), promotes=promoteListAeroLoads)
        # self.add_subsystem('aero_storm_1yr',    CCBladeLoads(analysis_options = analysis_options), promotes=promoteListAeroLoads)
        # self.add_subsystem('aero_storm_50yr',   CCBladeLoads(analysis_options = analysis_options), promotes=promoteListAeroLoads)
        # Add centrifugal and gravity loading to aero loading
        promotes=['tilt','theta','rhoA','z','totalCone','z_az']
        self.add_subsystem('curvature',         BladeCurvature(analysis_options = analysis_options),  promotes=['r','precone','precurve','presweep','3d_curv','x_az','y_az','z_az'])
        promoteListTotalLoads = ['r', 'theta', 'tilt', 'rhoA', '3d_curv', 'z_az', 'aeroloads_Omega', 'aeroloads_pitch']
        # self.add_subsystem('tot_loads_rated',       TotalLoads(analysis_options = analysis_options),      promotes=promoteListTotalLoads)
        self.add_subsystem('tot_loads_gust',        TotalLoads(analysis_options = analysis_options),      promotes=promoteListTotalLoads)
        # self.add_subsystem('tot_loads_storm_1yr',   TotalLoads(analysis_options = analysis_options),      promotes=promoteListTotalLoads)
        # self.add_subsystem('tot_loads_storm_50yr',  TotalLoads(analysis_options = analysis_options),      promotes=promoteListTotalLoads)
        promoteListFrame3DD = ['x_az','y_az','z_az','theta','r','A','EA','EIxx','EIyy','EIxy','GJ','rhoA','rhoJ','x_ec','y_ec','xu_strain_spar','xl_strain_spar','yu_strain_spar','yl_strain_spar','xu_strain_te','xl_strain_te','yu_strain_te','yl_strain_te']
        self.add_subsystem('frame',     RunFrame3DD(analysis_options = analysis_options),      promotes=promoteListFrame3DD)
        self.add_subsystem('tip_pos',   TipDeflection(),                                  promotes=['tilt','pitch_load'])
        self.add_subsystem('aero_hub_loads', AeroHubLoads(analysis_options = analysis_options), promotes = promoteListAeroLoads)
        self.add_subsystem('constr',    DesignConstraints(analysis_options = analysis_options, opt_options = opt_options))

        # Aero loads to total loads
        # self.connect('aero_rated.loads_Px',     'tot_loads_rated.aeroloads_Px')
        # self.connect('aero_rated.loads_Py',     'tot_loads_rated.aeroloads_Py')
        # self.connect('aero_rated.loads_Pz',     'tot_loads_rated.aeroloads_Pz')
        self.connect('aero_gust.loads_Px',      'tot_loads_gust.aeroloads_Px')
        self.connect('aero_gust.loads_Py',      'tot_loads_gust.aeroloads_Py')
        self.connect('aero_gust.loads_Pz',      'tot_loads_gust.aeroloads_Pz')
        # self.connect('aero_storm_1yr.loads_Px', 'tot_loads_storm_1yr.aeroloads_Px')
        # self.connect('aero_storm_1yr.loads_Py', 'tot_loads_storm_1yr.aeroloads_Py')
        # self.connect('aero_storm_1yr.loads_Pz', 'tot_loads_storm_1yr.aeroloads_Pz')
        # self.connect('aero_storm_50yr.loads_Px', 'tot_loads_storm_50yr.aeroloads_Px')
        # self.connect('aero_storm_50yr.loads_Py', 'tot_loads_storm_50yr.aeroloads_Py')
        # self.connect('aero_storm_50yr.loads_Pz', 'tot_loads_storm_50yr.aeroloads_Pz')

        # Total loads to strains
        self.connect('tot_loads_gust.Px_af', 'frame.Px_af')
        self.connect('tot_loads_gust.Py_af', 'frame.Py_af')
        self.connect('tot_loads_gust.Pz_af', 'frame.Pz_af')

        # Blade distributed deflections to tip deflection
        self.connect('frame.dx', 'tip_pos.dx_tip', src_indices=[-1])
        self.connect('frame.dy', 'tip_pos.dy_tip', src_indices=[-1])
        self.connect('frame.dz', 'tip_pos.dz_tip', src_indices=[-1])
        self.connect('3d_curv',  'tip_pos.3d_curv_tip', src_indices=[-1])

        # Strains from frame3dd to constraint
        self.connect('frame.strainU_spar', 'constr.strainU_spar')
        self.connect('frame.strainL_spar', 'constr.strainL_spar')
        self.connect('frame.flap_mode_freqs', 'constr.flap_mode_freqs')
        self.connect('frame.edge_mode_freqs', 'constr.edge_mode_freqs')

         
