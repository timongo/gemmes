#!/usr/bin/python

# Built-in
import sys

# Common
import numpy as np
from scipy.integrate import ode

# Library-specific
if sys.version[0]=='2':
    import _utils
    import _plot
else:
    try:
        import gemmes._utils as _utils
        import gemmes._plot as _plot
    except ModuleNotFoundError:
        import _utils
        import _plot
    except ImportError:
        import _utils
        import _plot

__all__ = ['GemmesIntegrator']

class GemmesIntegrator(object):
    """
    Integrator for the dynamical system of Gemmes
    Financial impacts of climate change mitigation policies and their macroeconomic implications: a stock-flow consistent approach
    https://www.tandfonline.com/doi/full/10.1080/14693062.2019.1698406?scroll=top&needAccess=true
    
    In order to ease comparisons, we choose to name variables according to 
    their name in the paper, rather than their corresponding name in the
    coping (first) paper. For instance, the growth rate of the population 
    is names deltaN instead of q, and the limit for the sigmoid is 
    Nbar instead of PN    
    """

    def __init__(self, name='',
                 alpha=0.02, # growth rate of productivity
                 deltaN=0.0305, # leading growth rate of workforce
                 Nbar=7.055925493, # maximum population in the logistic evolution
                 nu=2.7, # Constant capital-to-output ratio
                 pi1=0., # damage function parameter
                 pi2=0.00236, # damage function parameter
                 pi3=0.0000819, # damage function parameter
                 zeta3=6.754, # damage function parameter
                 fK=0., # (paper = 1./3) fraction of environmental damage allocated to the stock of capital
                 delta=0.04, # capital depreciation rate
                 sA=0, # Fraction of abatement costs that are subsidized
                 apC=0.02, # carbon price parameter
                 bpC=0., # carbon price parameter
                 conv10to15=1.160723971/1000., # conversion factor
                 deltagsigma=-0.001, # dynamics of emissivity
                 eta=0.5, # relaxation parameter of inflation 
                 etar=10., # relaxation parameter of the interest rate
                 mu=1.3, # markup of prices over the average cost
                 omitted=0.3, # offset for the production cost in the inflation
                 rstar=0.01, # Long-term interest rate target of the economy
                 phitaylor=0.5, # parameter characterizing the reactivity of the monetary policy
                 istar=0.02, #interest rate targeted by the monetary policy
                 srep=0.1, # Fraction of the outstanding debt repaid yearly
                 S=3.1, # Climate sensitivity
                 gammastar=0.0176, # Heat exchange coefficient between temperature layers
                 F2CO2=3.6813, # Change in radiative forcing resulting from doubling of CO2
                 CATpind=588., # CO2 preindustrial concentration in atmosphere
                 CUPpind=360., # CO2 preindustrial concentration in upper layer of ocean and biosphere
                 CLOpind=1720., # CO2 preindustrial concentration in bottom layer of the ocean
                 Fexo0=0.5, # Initial value of the exogenous radiative forcing
                 Fexo1=1., # value of the exogenous radiative forcing in 2100
                 phi12=0.0239069, # Transfer coefficient for carbon from the atmosphere to the upper ocean
                 phi23=0.0013409, # Transfer coefficient for carbon from the upper ocean/biosphere to the lower ocean
                 deltapbs=-.00505076337, # Exogenous growth rate of the back-stop technology price
                 theta=2.6, # parameter of the abatement cost function
                 phi0=-.291535421, # Constant of the short-term Philips curve
                 phi1=.468777035, # Slope of the short-term Philips curve
                 kappa0=.031774, # Constant of the investment function
                 kappa1=.575318414, # Slope of the investment function
                 kappamin=0., # Minimum of the investment function
                 kappamax=0.3, # Maximum of the investment function
                 div0=0.0512863, # Constant of the dividend function
                 div1=0.4729, # Slope of the dividend function
                 divmin=0., # Minimum of the dividend function
                 divmax=1.,# Maximum of the dividend function
                 C=49.75817526819656, # Heat capacity of the atmosphere, biosphere and upper ocean
                 C0=3.52, # Heat capacity of the deeper ocean
                 deltaEland=-.0220096, # Growth rate of land-use change CO2 of emissions
                 CR0=0.17, # Constant of the leverage function
                 CRlev=.1,# slope of the leverage function
                 #Initial conditions
                 CO2AT_ini=851.,
                 CO2UP_ini=460.,
                 CO2LO_ini=1740.,
                 d_ini=1.53282171,
                 Eind_ini=35.85,
                 Eland_ini=2.6,
                 gsigma_ini=-0.0105,
                 pbs_ini=547.2220801465,
                 n_ini=0.03,
                 N_ini=4.825484061,
                 T_ini=0.85,
                 T0_ini=0.0068,
                 Y_ini=59.7387,
                 lambda_ini=0.674828,
                 omega_ini=0.5782323,
                 r_ini=0.10627649250000004,
                 # Options
                 Damage = 'No', # 'No', 'Extreme'
                 # Integration parameters
                 dt=0.5, # half a year time step
                 tmax=84., # up to 2100
                 method='dop853'):
        """
        Initialization
        """

        self.name = name # Name of the run, useful for plotting multiple runs

        self.alpha = alpha # growth rate of productivity
        self.deltaN = deltaN # leading growth rate of workforce
        self.Nbar = Nbar # maximum population in the logistic evolution
        self.nu = nu # Constant capital-to-output ratio
        self.pi1 = pi1 # damage function parameter
        self.pi2 = pi2 # damage function parameter
        self.pi3 = pi3 # damage function parameter
        self.zeta3 = zeta3 # damage function parameter
        self.fK = fK # fraction of environmental damage allocated to the stock of capital
        self.delta = delta # capital depreciation rate
        self.sA = sA # Fraction of abatement costs that are subsidized
        self.apC = apC # carbon price parameter
        self.bpC = bpC # carbon price parameter
        self.conv10to15 = conv10to15 # conversion factor
        self.convCO2toC = 1./3.666 # conversion from tCO2 to tC
        self.deltagsigma = deltagsigma # dynamics of emissivity
        self.eta = eta # relaxation parameter of inflation 
        self.etar = etar # relaxation parameter of the interest rate
        self.mu = mu # markup of prices over the average cost
        self.omitted = omitted # offset for the production cost in the inflation
        self.rstar = rstar # Long-term interest rate target of the economy
        self.phitaylor = phitaylor # parameter characterizing the reactivity of the monetary policy
        self.istar = istar #interest rate targeted by the monetary policy
        self.srep = srep # Fraction of the outstanding debt repaid yearly
        self.gammastar = gammastar # Heat exchange coefficient between temperature layers
        self.F2CO2 = F2CO2 # Change in radiative forcing resulting from doubling of CO2
        self.S = S # Climate sensitivity to doubling of CO2
        self.rho = self.F2CO2/self.S # Radiative feedback parameter
        self.CATpind = CATpind # CO2 preindustrial concentration in atmosphere
        self.CUPpind = CUPpind # CO2 preindustrial concentration in upper layer of ocean and biosphere
        self.CLOpind = CLOpind # CO2 preindustrial concentration in bottom layer of the ocean
        self.Fexo0 = Fexo0 # Initial value of the exogenous radiative forcing
        self.Fexo1 = Fexo1 # value of the exogenous radiative forcing in 2100
        self.phi12 = phi12 # Transfer coefficient for carbon from the atmosphere to the upper ocean
        self.phi23 = phi23 # Transfer coefficient for carbon from the upper ocean/biosphere to the lower ocean

        # Phi transfer matrix
        CATUP = self.CATpind/self.CUPpind
        CUPLO = self.CUPpind/self.CLOpind
        self.phimat = np.zeros((3,3))
        self.phimat[0,0] = -self.phi12
        self.phimat[0,1] = self.phi12*CATUP
        self.phimat[1,0] = self.phi12
        self.phimat[1,1] = -self.phi12*CATUP - self.phi23
        self.phimat[1,2] = self.phi23*CUPLO
        self.phimat[2,1] = self.phi23
        self.phimat[2,2] = -self.phi23*CUPLO

        self.deltapbs = deltapbs # Exogenous growth rate of the back-stop technology price
        self.theta = theta # parameter of the abatement cost function
        self.phi0 = phi0 # Constant of the short-term Philips curve
        self.phi1 = phi1 # Slope of the short-term Philips curve
        self.kappa0 = kappa0 # Constant of the investment function
        self.kappa1 = kappa1 # Slope of the investment function
        self.kappamin = kappamin # Minimum of the investment function
        self.kappamax = kappamax # Maximum of the investment function
        self.div0 = div0 # Constant of the dividend function
        self.div1 = div1 # Slope of the dividend function
        self.divmin = divmin # Minimum of the dividend function
        self.divmax = divmax# Maximum of the dividend function
        self.C = C # Heat capacity of the atmosphere, biosphere and upper ocean
        self.C0 = C0 # Heat capacity of the deeper ocean
        self.deltaEland = deltaEland # Growth rate of land-use change CO2 of emissions
        self.CR0 = CR0 # Constant of the leverage function
        self.CRlev = CRlev# slope of the leverage function
        # Initial conditions
        self.CO2AT_ini = CO2AT_ini
        self.CO2UP_ini = CO2UP_ini
        self.CO2LO_ini = CO2LO_ini
        self.d_ini = d_ini
        self.Eind_ini = Eind_ini
        self.Eland_ini = Eland_ini
        self.gsigma_ini = gsigma_ini
        self.pbs_ini = pbs_ini
        self.n_ini = n_ini
        self.N_ini = N_ini
        self.T_ini = T_ini
        self.T0_ini = T0_ini
        self.Y_ini = Y_ini
        self.lambda_ini = lambda_ini
        self.omega_ini = omega_ini
        self.r_ini = r_ini
        # Options
        self.Damage = Damage
        # Integration parameters
        self.dt = dt # half a year time step
        self.tmax = tmax # up to 2100
        self.method=method # integration algorithm
        # Store Param class to facilitate handling / understanding
        self.params = _utils.Params()

        # Initilaze the solution
        self.sol = None
        
        # Variables for plotting
        self.eco_vars = ['omega','lambda','d','N','rCB','pi','piK']
        self.gdp_vars = ['Y','Y0','g0','i']
        self.emi_vars = ['Eind','Eland','sigma','gsigma','n']
        self.clim_vars = ['T','T0','CO2at','CO2up','CO2lo']
        self.pcar_vars = ['pC','pbs']
        self.dam_vars = ['A','DamY','DamK','Dam']

    def Phi(self,lam):
        """
        Philips curve
        """
        return self.phi0 + self.phi1*lam

    def Taylor(self,i):
        """
        Interest rate as a function of inflation
        """
        
        r = self.rstar - self.phitaylor*self.istar + (1.+self.phitaylor)*i
        try:
            return min(max(r,0),1)
        except ValueError:
            r[r<0.] = 0.
            r[r>1.] = 1.
            return r
    
    def Kappa(self,pi):
        """
        Investment function
        """
        kappa = self.kappa0 + self.kappa1*pi
        try:
            if kappa < self.kappamin:
                return self.kappamin
            elif kappa > self.kappamax:
                return self.kappamax
            else:
                return kappa
        except ValueError:
            kappa[kappa<self.kappamin] = self.kappamin
            kappa[kappa>self.kappamax] = self.kappamax
            return kappa

    def Delta(self,pi):
        """
        Dividend function
        """
        div = self.div0 + self.div1*pi
        try:
            if div < self.divmin:
                return self.divmin
            elif div > self.divmax:
                return self.divmax
            else:
                return div
        except ValueError:
            div[div<self.divmin] = self.divmin
            div[div>self.divmax] = self.divmax
            return div

    def Tau(self,lev):
        """
        Leverage function
        """
        CR = self.CR0 + self.CRlev*lev
        
        try:
            return min(CR,1)
        except ValueError:
            CR[CR<0.] = 0.
            CR[CR>1.] = 1.
            return CR

    def InitialConditions(self):
        """
        Returns
        sigma = Eind/((1-n)*Y)
        pC = root of n = min((pC/pbs)**(1/(theta-1)),1)
        etc
        """

        t2016 = 1.
        t2100 = 2100. - 2016.

        # Carbon price
        pC = self.pbs_ini*self.n_ini**(self.theta - 1.)

        # emissivity
        sigma = self.Eind_ini/((1-self.n_ini)*self.Y_ini)
        
        # Abatement cost
        A = sigma/1000.*self.pbs_ini*self.n_ini**self.theta/self.theta

        # Temperature damage
        if self.Damage == 'No':
            Dam = 0.
        elif self.Damage == 'Q':
            Dam = 1. - 1./(1 + self.pi1*self.T_ini
                         + self.pi2*self.T_ini**2)
        elif self.Damage == 'Extreme':
            Dam = 1. - 1./(1 + self.pi1*self.T_ini
                         + self.pi2*self.T_ini**2
                         + self.pi3*self.T_ini**self.zeta3)
        DamK = self.fK*Dam
        DamY = 1. - (1.-Dam)/(1.-DamK)
        deltaD = (self.delta + DamK)

        # Total cost of climate change
        TC = (1-DamY)*(1-A)
        # Capital
        Y0 = self.Y_ini/TC
        K = Y0*self.nu
        L = self.lambda_ini*self.N_ini
        a = Y0/L

        # Initial price index
        p = 1.
        # production cost
        c = self.omega_ini
        # inflation
        i = self.eta*(self.mu*(c+self.omitted)-1.)
        # Central Bank interest rate
        rCB = self.Taylor(i)

        # Net transfers between public and private sectors
        Transfers = self.sA*A*Y0 - pC*self.conv10to15*self.Eind_ini

        # wages
        w = self.omega_ini*p*self.Y_ini/L
        # Nominal debt
        D = self.d_ini*p*self.Y_ini

        Pi = p*self.Y_ini - w*L - rCB*D - deltaD*p*K + p*Transfers
        piK = Pi/(p*K)
        pi = Pi/(p*self.Y_ini)

        Find = self.F2CO2*np.log(self.CO2AT_ini/self.CATpind)/np.log(2.)
        Fexo = self.Fexo0

        # leverage
        leverage = self.d_ini*TC/self.nu
        CR = self.Tau(leverage)

        # Economic growth rate
        g0 = ((1.-CR)*(1.25*self.Kappa(pi)*TC/self.nu - deltaD)
             +CR*(piK - self.Delta(piK) - TC*self.srep*self.d_ini/self.nu))

        return K,D,w,a,p,sigma,pC,Y0,g0,i,A,pi,piK,Dam,DamK,DamY,Fexo,Find,rCB
        
    def SystemFun(self,t,u):
        """
        u[0] = K
        u[1] = N
        u[2] = D
        u[3] = w
        u[4] = a
        u[5] = p
        u[6] = Eland
        u[7] = sigma
        u[8] = gsigma
        u[9:12] = CO2-e concentration (AT, UP, LO)
        u[12] = T
        u[13] = T0
        u[14] = pbs
        u[15] = pC
        """

        t2016 = 1.
        t2100 = 2100. - 2016.

        K = u[0]
        N = u[1]
        D = u[2]
        w = u[3]
        a = u[4]
        p = u[5]
        Eland = u[6]
        sigma = u[7]
        gsigma = u[8]
        CO2 = np.array(u[9:12])
        T = u[12]
        T0 = u[13]
        pbs = u[14]
        pC = u[15]

        # saturation of pC
        #expo = (self.theta-1.)/self.theta
        #pC = pbs*min(pC/pbs,0.95*(self.theta/(sigma*pbs))**expo)

        # Abatement cost
        n = min((pC/((1.-self.sA)*pbs))**(1./(self.theta-1.)),1)
        A = sigma/1000.*pbs*n**self.theta/self.theta

        # Temperature damage
        if self.Damage == 'No':
            Dam = 0.
        elif self.Damage == 'Q':
            Dam = 1. - 1./(1 + self.pi1*T
                         + self.pi2*T**2)
        elif self.Damage == 'Extreme':
            Dam = 1. - 1./(1 + self.pi1*T
                         + self.pi2*T**2
                         + self.pi3*T**self.zeta3)
        DamK = self.fK*Dam
        DamY = 1. - (1.-Dam)/(1.-DamK)
        deltaD = (self.delta + DamK)

        # Total cost of climate change
        TC = (1-DamY)*(1-A)
        # GDP
        Y0 = K/self.nu
        Y = TC*Y0
        # Workforce
        L = Y0/a
        # Wage share
        omega = w*L/(p*Y)
        # employment parameter
        lam = L/N
        # Debt ratio
        d = D/(p*Y)
        # Industrial emission
        Eind = Y0*sigma*(1.-n)

        # inflation
        i = self.eta*(self.mu*(omega+self.omitted) - 1.)
        # central bank interest rate
        rCB = self.Taylor(i)
        # Net transfers between public and private sectors
        Transfers = self.sA*A*Y0 - pC*self.conv10to15*Eind
        # profits
        Pi = p*Y - w*L - rCB*D - deltaD*p*K + p*Transfers
        # piK = (TC*(1.-omega-rCB*d) - pC*self.conv10to15*sigma*(1-n) + self.sA*A)/self.nu - deltaD
        # pi = self.nu/TC*piK
        piK = Pi/(p*K)
        pi = Pi/(p*Y)
        # leverage
        leverage = d*TC/self.nu
        CR = self.Tau(leverage)
        # Population growth
        beta = self.deltaN*(1.-N/self.Nbar)
        # Temperature change
        CO2AT = CO2[0]
        Find = self.F2CO2*np.log(CO2AT/self.CATpind)/np.log(2.)
        Fexo = min(self.Fexo0 + (self.Fexo1-self.Fexo0)*t/84.,self.Fexo1)
        F = Find + Fexo
        E = Eind + Eland

        # Investment
        Id = 1.25*self.Kappa(pi)*Y
        Pir = Pi - self.Delta(piK)*p*K
        I = CR*(Pir/p + deltaD*K - self.srep*D/p)+(1.-CR)*Id
        # Debt demand
        Dd = p*Id - Pir + self.srep*D - deltaD*p*K

        # Computation of results
        # Kdot
        Kdot = I - deltaD*K
        # Ndot
        Ndot = beta*N
        # Ddot
        Ddot = (1.-CR)*Dd - self.srep*D
        # wdot
        wdot = self.Phi(lam)*w
        # adot
        adot = a*self.alpha
        # pdot
        pdot = i*p
        # Elanddot
        Elanddot = self.deltaEland*Eland
        # sigmadot
        sigmadot = gsigma*sigma
        # gsigmadot
        gsigmadot = self.deltagsigma*gsigma
        # CO2dot
        CO2dot = np.array([E*self.convCO2toC,0,0]) + np.dot(self.phimat,CO2)            
        # Tdot
        Tdot = (F - self.rho*T - self.gammastar*(T-T0))/self.C
        # T0dot
        T0dot = self.gammastar*(T-T0)/self.C0            
        # pbsdot
        pbsdot = pbs*self.deltapbs
        # pCdot
        pCdot = pC*(self.apC + self.bpC/(t+t2016))

        # Auxiliary variables computation
        # Economic growth rate
        g0 = ((1.-CR)*(1.25*self.Kappa(pi)*TC/self.nu - deltaD)
             +CR*(piK - self.Delta(piK) - TC*self.srep*d/self.nu))
        
        self.omega = omega
        self.lam = lam
        self.d = d
        self.Y0 = Y0
        self.Y = Y
        self.g0 = g0
        self.Eind = Eind
        self.i = i
        self.A = A
        self.n = n
        self.pi = pi
        self.piK = piK
        self.Dam = Dam
        self.DamK = DamK
        self.DamY = DamY
        self.Fexo = Fexo
        self.Find = Find
        self.rCB = rCB

        return [Kdot,
                Ndot,
                Ddot,
                wdot,
                adot,
                pdot,
                Elanddot,
                sigmadot,
                gsigmadot,
                CO2dot[0], # dCO2/dt (AT)
                CO2dot[1], # dCO2/dt (UP)
                CO2dot[2], # dCO2/dt (LO)
                Tdot,
                T0dot,
                pbsdot,
                pCdot] 


    def Solve(self, plot=True, fignumber=1, verb=-1):
        """
        Solve the system
        """

        (K_ini,
         D_ini,
         w_ini,
         a_ini,
         p_ini,
         sigma_ini,
         pC_ini,
         Y0_ini,
         g0_ini,
         i_ini,
         A_ini,
         pi_ini,
         piK_ini,
         Dam_ini,
         DamK_ini,
         DamY_ini,
         Fexo_ini,
         Find_ini,
         rCB_ini) = self.InitialConditions()
        
        U_ini = [K_ini,
                 self.N_ini,
                 D_ini,
                 w_ini,
                 a_ini,
                 p_ini,
                 self.Eland_ini,
                 sigma_ini,
                 self.gsigma_ini,
                 self.CO2AT_ini,
                 self.CO2UP_ini,
                 self.CO2LO_ini,
                 self.T_ini,
                 self.T0_ini,
                 self.pbs_ini,
                 pC_ini,
                 self.omega_ini,
                 self.lambda_ini,
                 self.d_ini,
                 Y0_ini,
                 self.Y_ini,
                 g0_ini,
                 self.Eind_ini,
                 i_ini,
                 A_ini,
                 self.n_ini,
                 pi_ini,
                 piK_ini,
                 Dam_ini,
                 DamK_ini,
                 DamY_ini,
                 Fexo_ini,
                 Find_ini,
                 rCB_ini]

        U_ini_sys = [K_ini,
                     self.N_ini,
                     D_ini,
                     w_ini,
                     a_ini,
                     p_ini,
                     self.Eland_ini,
                     sigma_ini,
                     self.gsigma_ini,
                     self.CO2AT_ini,
                     self.CO2UP_ini,
                     self.CO2LO_ini,
                     self.T_ini,
                     self.T0_ini,
                     self.pbs_ini,
                     pC_ini]

        nts = int(self.tmax/self.dt)+1
        t = np.zeros(nts)
        U = np.zeros((nts,34))
        U[0,:] = U_ini

        if self.method!='rk4':
            system = ode(self.SystemFun).set_integrator(self.method)
            system.set_initial_value(U_ini_sys,0.)        
        
            k=0
            while system.successful() and k < nts-1:
                k=k+1
                system.integrate(system.t + self.dt)
                t[k] = system.t
                U[k,:16] = system.y[:]
                U[k,16] = self.omega
                U[k,17] = self.lam
                U[k,18] = self.d
                U[k,19] = self.Y0
                U[k,20] = self.Y
                U[k,21] = self.g0
                U[k,22] = self.Eind
                U[k,23] = self.i 
                U[k,24] = self.A
                U[k,25] = self.n
                U[k,26] = self.pi
                U[k,27] = self.piK
                U[k,28] = self.Dam
                U[k,29] = self.DamK
                U[k,30] = self.DamY
                U[k,31] = self.Fexo
                U[k,32] = self.Find
                U[k,33] = self.rCB
        else:
            u = U_ini_sys
            k=0
            while k < nts-1:
                k=k+1
                u = self.rk4(self.SystemFun,t[k-1],u)
                t[k] = t[k-1]+self.dt
                U[k,:16] = u[:]
                U[k,16] = self.omega
                U[k,17] = self.lam
                U[k,18] = self.d
                U[k,19] = self.Y0
                U[k,20] = self.Y
                U[k,21] = self.g0
                U[k,22] = self.Eind
                U[k,23] = self.i 
                U[k,24] = self.A
                U[k,25] = self.n
                U[k,26] = self.pi
                U[k,27] = self.piK
                U[k,28] = self.Dam
                U[k,29] = self.DamK
                U[k,30] = self.DamY
                U[k,31] = self.Fexo
                U[k,32] = self.Find
                U[k,33] = self.rCB
            
        # Keep track of the solution for later use (plot, treatment...)
        self.sol = {'t':t,'U':U}

        # Optional plotting right away after computation
        if plot:
            self.plot(fignumber=fignumber)

        return t,U

    def rk4(self,f,t,u):

        dt = self.dt

        k1 = np.array(f(t,u))
        k2 = np.array(f(t+0.5*dt,u+k1*0.5*dt))
        k3 = np.array(f(t+0.5*dt,u+k2*0.5*dt))
        k4 = np.array(f(t+dt,u+k3*dt))
        
        return u+(k1+2.*k2+2.*k3+k4)*dt/6.
        
    
    def plot(self, vars = None, fs=None, dmargin=None, draw=True, fignumber=1):
        """
        Plot the solution
        """
        
        if vars is None:
            vars = self.eco_vars+self.gdp_vars

        assert self.sol is not None, "System not solved yet ! "
        lax = _plot.plot_basic(self, vars=vars, fs=fs, dmargin=dmargin, draw=draw,fignumber=fignumber)
        return lax


    def compare(self, lGemInt, Type='basic',
               fs=None, dmargin=None, draw=True):
        assert self.sol is not None, "System not solved yet !"
        # Prepare inputs
        if type(lGemInt) is not list:
            lGemInt = [lGemInt]
        msg = "Some systems not solved yet ! "
        assert all([ss.sol is not None for ss in lGemInt]), msg

        if Type=='basic':
            lax = _plot.plot_basic([self]+lGemInt, fs=fs, dmargin=dmargin, draw=draw)
        return lax

