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


__all__ = ['GemmesIntegrator']

class GemmesIntegrator(object):
    """
    Integrator for the dynamical system of Gemmes
    Coping with the Collapse: A Stock-Flow Consistent Monetary Macrodynamics of Global Warming
    https://www.sciencedirect.com/science/article/pii/S0921800916309569
    """

    def __init__(self, name='',
                 alpha=0.02,
                 q=0.0305,
                 PN=7.056,
                 nu=2.7,
                 pi1=0.,
                 pi2=0.00236,
                 pi3=0.00000507,
                 zeta3=6.754,
                 fK=1./3,
                 delta=0.04,
                 r=0.03,
                 apC=0.02,
                 bpC=0.,
                 conv10to15=1.1607/1000,
                 deltagsigma=-0.001,
                 eta=0.5,
                 mu=1.3,
                 rstar=0.01,
                 rtaylor=0.5,
                 istar=0.02,
                 S=3.1,
                 gammastar=0.0176,
                 F2CO2=3.681,
                 CATpind=588.,
                 CUPpind=360.,
                 CLOpind=1720.,
                 Fexo0=0.5,
                 Fexo1=1.,
                 phi12=0.024,
                 phi23=0.001,
                 deltapbs=-0.005,
                 theta=2.6,
                 phi0=-0.292,
                 phi1=0.469,
                 kappa0=0.0318,
                 kappa1=0.575,
                 kappamin=0.,
                 kappamax=0.3,
                 div0=0.138,
                 div1=0.473,
                 divmin=0.,
                 divmax=0.3,
                 C=49.76,
                 C0=3.52,
                 #Initial conditions
                 CO2AT_ini=851.,
                 CO2UP_ini=460.,
                 CO2LO_ini=1740.,
                 d_ini=1.53,
                 deltaEland=-0.022,
                 Eind_ini=35.85,
                 Eland_ini=2.6,
                 gsigma_ini=-0.0152,
                 pbs_ini=547.22,
                 n_ini=0.03,
                 N_ini=4.83,
                 T_ini=0.85,
                 T0_ini=0.0068,
                 Y_ini=59.74,
                 lambda_ini=0.675,
                 omega_ini=0.578,
                 # Integration parameters
                 dt=0.5, # half a year time step
                 tmax=84., # up to 2100
                 nofeedback=False):

        """
        Initialization
        """

        self.name = name # Name of the run, useful for plotting multiple runs

        self.alpha = alpha # Constant growth rate of labor productivity
        self.q = q # Speed of growth of the workforce dynamics
        self.PN = PN # Upper limit of the workforce dynamics in billions
        self.nu = nu # Constant capital-to-output ratio
        if nofeedback:
            self.pi1 = 0. # Damage function parameter
            self.pi2 = 0. # Damage function parameter
            self.pi3 = 0. # Damage function parameter
            self.Afac = 0.
            self.nfac = 1.
        else:
            self.pi1 = pi1 # Damage function parameter
            self.pi2 = pi2 # Damage function parameter
            self.pi3 = pi3 # Damage function parameter
            self.Afac = conv10to15 # conversion factor due to weird units
            self.nfac = 1.
        self.zeta3 = zeta3 # Damage function parameter
        self.fK = fK # Fraction of environmental damage allocated to the stock of capital
        self.delta = delta # Depreciation rate of capital
        self.r = r # Short-term interest rate of the economy
        self.apC = apC # Carbon price evolution parameter
        self.bpC = bpC # Carbon price evolution parameter
        self.conv10to15 = conv10to15 # Scaling factor with trillion dollars for the carbon price
        self.convCO2toC = 1./3.66667 # conversion from tCO2 to tC
        self.deltagsigma = deltagsigma # Variation rate of the growth of emission intensity
        self.eta = eta # Relaxation parameter of the inflation
        self.mu = mu # "markup"
        self.rstar = rstar # Parameter for the Taylor function defining the interest rate as a function of inflation
        self.rtaylor = rtaylor # Parameter for the Taylor function defining the interest rate as a function of inflation
        self.istar = istar # Parameter for the Taylor function defining the interest rate as a function of inflation
        self.gammastar = gammastar # Heat exchange coefficient between temperature layers
        self.F2CO2 = F2CO2 # Change in the radiative forcing resulting from
        # a doubling of CO2-e concentration wrt the preindustrial period
        self.S = S # Climate sensitivity to doubling of CO2
        self.rho = self.F2CO2/self.S # Radiative feedback parameter
        self.CATpind = CATpind # CO2-e preindustrial concentration (atmosphere)
        self.CUPpind = CUPpind # CO2-e preindustrial concentration (biosphere + upper ocean)
        self.CLOpind = CLOpind # CO2-e preindustrial concentration (deep ocean)
        self.Fexo0 = Fexo0 # Parameter of exogenous radiative forcing
        self.Fexo1 = Fexo1 # Parameter of exogenous radiative forcing
        self.deltaEland = deltaEland # Growth rate of land use change CO2-e emissions
        self.phi12 = phi12 # Transfer coefficient for carbon from the atmosphere to
        # the upper ocean/biosphere
        self.phi23 = phi23 # Transfer coefficient for carbon from the upper ocean/
        # biosphere to the lower ocean

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

        self.deltapbs = deltapbs # Variation rate of land use change CO2-e emissions
        self.theta = theta # Parameter of the abatement cost function
        self.phi0 = phi0 # Constant term of the Philips curve
        self.phi1 = phi1 # Slope of the Philips curve
        self.kappa0 = kappa0 # Constant term of the investment function
        self.kappa1 = kappa1 # Slope of the investment function
        self.kappamin = kappamin # Lower boundary of the investment function
        self.kappamax = kappamax # Upper boundary of the investment function
        self.div0 = div0 # Constant term of the dividend function
        self.div1 = div1 # Slope of the dividend function
        self.divmin = divmin # Lower boundary of the dividend function
        self.divmax = divmax # Upper boundary of the dividend function
        self.C = C # Heat capacity of the atmosphere, biosphere and upper ocean
        self.C0 = C0 # Heat capacity of the deeper ocean

        # Initial conditions
        self.CO2AT_ini = CO2AT_ini 
        self.CO2UP_ini = CO2UP_ini
        self.CO2LO_ini = CO2LO_ini
        self.d_ini = d_ini
        self.Eind_ini = Eind_ini
        self.Eland_ini = Eland_ini
        self.gsigma_ini = gsigma_ini
        self.pbs_ini = pbs_ini
        if nofeedback:
            self.n_ini = 0.
        else:
            self.n_ini = n_ini
        self.N_ini = N_ini
        self.T_ini = T_ini
        self.T0_ini = T0_ini
        self.Y_ini = Y_ini
        self.lambda_ini = lambda_ini
        self.omega_ini = omega_ini

        self.dt = dt
        self.tmax = tmax

        # Store Param class to facilitate handling / understanding
        self.params = _utils.Params()

        # Initilaze the solution
        self.sol = None
        
        # Variables for plotting
        self.eco_vars = ['omega','lambda','d','N','pi']
        self.gdp_vars = ['Y','Y0','g','i']
        self.emi_vars = ['Eind','Eland','sigma','gsigma','n']
        self.clim_vars = ['T','T0','CO2at','CO2up','CO2lo']
        self.pcar_vars = ['pC','pbs']
        self.dam_vars = ['A','DY','DK','D']

    def Phi(self,x):
        """
        Philips curve
        """
        # return -0.04006410256410257 + 6.410256410256412e-05/(1.-x)**2 # Grasselli
        return self.phi0 + self.phi1*x

    def Taylor(self,i):
        """
        Interest rate as a function of inflation
        """
        
        r = self.rstar - self.rtaylor*self.istar + (1.+self.rtaylor)*i
        try:
            return min(max(r,0),1)
        except ValueError:
            r[r<0.] = 0.
            r[r>1.] = 1.
            return r
    
    def Kappa(self,x):
        """
        Investment function
        """
        kappa = self.kappa0 + self.kappa1*x
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

    def Kappainv(self,kappa):
        """
        inverse of kappa
        """
        
        if kappa>self.kappamax:
            return None
        elif kappa<self.kappamin:
            return None
        else:
            return (kappa - self.kappa0)/self.kappa1

    def Delta(self,x):
        """
        Dividend function
        """
        div = self.div0 + self.div1*x
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

    def Deltainv(self,delta):
        """
        inverse of delta
        """
        
        if delta>self.deltamax:
            return None
        elif delta<self.deltamin:
            return None
        else:
            return (delta - self.delta0)/self.delta1

    def Solve(self, plot=True, fignumber=1, verb=-1):

        def InitialConditions(omega,d,Y,pbs,T,n,Eind):
            """
            Returns
            sigma = ((1- DY)/Y)/((1-n)/Eind + pbs*(1-DY)/Y)
            pC = root of n = min((pC/pbs)**(1/(theta-1)),1)
            """

            D = 1. - 1./(1 + self.pi1*T
                         + self.pi2*T**2
                         + self.pi3*T**self.zeta3)
            DK = self.fK*D
            DY = 1. - (1.-D)/(1.-DK)
            deltaD = (self.delta + DK)
            
            sigma = ((1. - DY)/Y)/((1.-n)/Eind + n**self.theta/self.theta*pbs*(1.-DY)/Y)
            A = self.Afac*sigma*pbs*n**self.theta/self.theta
            pC = pbs*n**(self.theta - 1.)

            TotalCost = (1-DY)*(1-A)
            # i = self.eta*(self.mu*(omega+0.3) - 1.)
            # r = self.Taylor(i)
            r = self.r
            pi = (1. - omega - r*d
                  -(pC*self.conv10to15*sigma*(1-n) + deltaD*self.nu)/TotalCost)
            # Unitary cost of production
            c = omega + r*d + self.Delta(pi) + (pC*self.conv10to15*sigma*(1-n) + deltaD*self.nu)/TotalCost
            #c = omega + self.Delta(pi) + (pC*self.conv10to15*sigma*(1-n) + deltaD*self.nu)/TotalCost
            # Inflation
            i = self.eta*(self.mu*c - 1.)

            g = self.Kappa(pi)*TotalCost/self.nu - deltaD 

            Y0 = Y/((1.-DY)*(1.-A))

            return sigma,pC,i,A,pi,D,DK,DY,g,Y0
            
        def f(t,u):
            """
            u[0] = omega, the share of work in the GDP
            u[1] = lambda, the work rate
            u[2] = d, the debt ratio
            u[3] = N, the population
            u[4] = T, the temperature of the atmosphere
            u[5] = T0, the temperature of the deep ocean
            u[6] = Y, the GDP
            u[7] = sigma, the emission intensity
            u[8] = gsigma, the rate of increase of emission intensity
            u[9:12] = CO2-e concentration (AT, UP, LO)
            u[12] = Eland, the land-use emissions
            u[13] = pC, the carbon price
            u[14] = pbs, the backstop technology price
            """

            omega = u[0]
            lam = u[1]
            d = u[2]
            N = u[3]
            T = u[4]
            T0 = u[5]
            Y = u[6]
            sigma = u[7]
            gsigma = u[8]
            CO2 = np.array(u[9:12])
            Eland = u[12]
            pC = u[13]
            pbs = u[14]
            
            # saturation of pC
            expo = (self.theta-1.)/self.theta
            pC = pbs*min(pC/pbs,0.95*(self.theta/(sigma*pbs))**expo)

            # Abatement cost
            n = self.nfac*min((pC/pbs)**(1./(self.theta-1.)),1)
            A = self.Afac*sigma*pbs*n**self.theta/self.theta
            # Temperature damage
            D = 1. - 1./(1 + self.pi1*T
                         + self.pi2*T**2
                         + self.pi3*T**self.zeta3)
            DK = self.fK*D
            DY = 1. - (1.-D)/(1.-DK)
            deltaD = (self.delta + DK)
            # Total cost of climate change
            TotalCost = (1-DY)*(1-A)
            # i = self.eta*(self.mu*(omega+0.3) - 1.)
            # Profit rate
            # r = self.Taylor(i)
            r = self.r
            pi = (1. - omega - r*d
                  -(pC*self.conv10to15*sigma*(1-n) + deltaD*self.nu)/TotalCost)
            # Unitary cost of production
            c = omega + r*d + self.Delta(pi) + (pC*self.conv10to15*sigma*(1-n) + deltaD*self.nu)/TotalCost
            # c = omega + self.Delta(pi) + (pC*self.conv10to15*sigma*(1-n) + deltaD*self.nu)/TotalCost
            # Inflation
            i = self.eta*(self.mu*c - 1.)
            # Economic growth rate
            g = self.Kappa(pi)*TotalCost/self.nu - deltaD
            # Population growth
            beta = self.q*(1.-N/self.PN)
            # Temperature change
            CO2AT = CO2[0]
            Find = self.F2CO2*np.log(CO2AT/self.CATpind)/np.log(2.)
            t2016 = 1.
            t2100 = 2100. - 2016.
            F = Find + self.Fexo0 + (self.Fexo1-self.Fexo0)*(t - t2016)/(t2100 - t2016)
            Y0 = Y/((1.-DY)*(1.-A))
            Eind = Y0*sigma*(1-n)
            E = Eind + Eland
            CO2dot = np.array([E*self.convCO2toC,0,0]) + np.dot(self.phimat,CO2)

            self.Eind = Eind
            self.i = i
            self.A = A
            self.pi = pi
            self.D = D
            self.DK = DK
            self.DY = DY
            self.g = g
            self.n = n
            self.Y0 = Y0
            
            return [omega*(self.Phi(lam) - i - self.alpha), #domega/dt
                    lam*(g - self.alpha - beta), # dlam/dt
                    #(-d*(g + i) # dd/dt
                    #  +self.Kappa(pi) + self.Delta(pi) - pi
                    #  -self.nu*deltaD/TotalCost),
                    (-d*(g + i - self.r) # dd/dt
                      +omega + self.Kappa(pi) + self.Delta(pi) - 1
                      +pC*self.conv10to15*sigma*(1-n)),
                    N*beta, #dN.dt
                    (F - self.rho*T - self.gammastar*(T-T0))/self.C, #dT/dt
                    self.gammastar*(T-T0)/self.C0, # dT0/dt
                    Y*g, #dY/dt
                    sigma*gsigma, #dsigma/dt
                    gsigma*self.deltagsigma, #dgsigma/dt
                    CO2dot[0], # dCO2/dt (AT)
                    CO2dot[1], # dCO2/dt (UP)
                    CO2dot[2], # dCO2/dt (LO)
                    Eland*self.deltaEland, #dEland/dt
                    pC*(self.apC + self.bpC/(t+t2016)), #dpC/dt
                    pbs*self.deltapbs] # dpbs/dt


        var_ini = InitialConditions(self.omega_ini,
                                    self.d_ini,
                                    self.Y_ini,
                                    self.pbs_ini,
                                    self.T_ini,
                                    self.n_ini,
                                    self.Eind_ini)

        sigma_ini,pC_ini,i_ini,A_ini,pi_ini,D_ini,DK_ini,DY_ini,g_ini,Y0_ini = var_ini
        
        U_ini = [self.omega_ini,
                 self.lambda_ini,
                 self.d_ini,
                 self.N_ini,
                 self.T_ini,
                 self.T0_ini,
                 self.Y_ini,
                 sigma_ini,
                 self.gsigma_ini,
                 self.CO2AT_ini,
                 self.CO2UP_ini,
                 self.CO2LO_ini,
                 self.Eind_ini,
                 self.Eland_ini,
                 pC_ini,
                 self.pbs_ini,
                 i_ini,
                 A_ini,
                 pi_ini,
                 D_ini,
                 DK_ini,
                 DY_ini,
                 g_ini,
                 self.n_ini,
                 Y0_ini]

        U_ini_sys = [self.omega_ini,
                     self.lambda_ini,
                     self.d_ini,
                     self.N_ini,
                     self.T_ini,
                     self.T0_ini,
                     self.Y_ini,
                     sigma_ini,
                     self.gsigma_ini,
                     self.CO2AT_ini,
                     self.CO2UP_ini,
                     self.CO2LO_ini,
                     self.Eland_ini,
                     pC_ini,
                     self.pbs_ini]
        
        system = ode(f).set_integrator('dop853', verbosity=verb)
        system.set_initial_value(U_ini_sys,self.dt)
        
        nts = int(self.tmax/self.dt)
        t = np.zeros(nts)
        U = np.zeros((nts,25))
        U[0,:] = U_ini
        #
        k=0
        while system.successful() and k < nts-1:
            k=k+1
            system.integrate(system.t + self.dt)
            t[k] = system.t
            U[k,:12] = system.y[:12]
            U[k,12] = self.Eind
            U[k,13:16] = system.y[12:15]
            U[k,16] = self.i
            U[k,17] = self.A
            U[k,18] = self.pi
            U[k,19] = self.D
            U[k,20] = self.DK
            U[k,21] = self.DY
            U[k,22] = self.g
            U[k,23] = self.n
            U[k,24] = self.Y0
            
        # Keep track of the solution for later use (plot, treatment...)
        self.sol = {'t':t,'U':U}

        # Optional plotting right away after computation
        if plot:
            self.plot(fignumber=fignumber)

        return t,U

    def plot(self, vars = eco_vars+gdp_vars, fs=None, dmargin=None, draw=True, fignumber=1):
        """
        Plot the solution
        """
        
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


class GemmesIntegrator_V2(object):
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
                 Nbar=7.056, # maximum population in the logistic evolution
                 nu=2.7, # Constant capital-to-output ratio
                 pi1=0., # damage function parameter
                 pi2=0.00236, # damage function parameter
                 pi3=0.0000819, # damage function parameter
                 zeta3=6.754, # damage function parameter
                 fK=1./3, # fraction of environmental damage allocated to the stock of capital
                 delta=0.04, # capital depreciation rate
                 sA=0, # Fraction of abatement costs that are subsidized
                 apC=0.02, # carbon price parameter
                 bpC=0., # carbon price parameter
                 conv10to15=1.1607/1000, # conversion factor
                 deltagsigma=-0.001, # dynamics of emissivity
                 eta=0.192, # relaxation parameter of inflation 
                 etar=10., # relaxation parameter of the interest rate
                 mu=1.875, # markup of prices over the average cost
                 rstar=0.01, # Long-term interest rate target of the economy
                 phitaylor=0.5, # parameter characterizing the reactivity of the monetary policy
                 istar=0.02, #interest rate targeted by the monetary policy
                 srep=0.1, # Fraction of the outstanding debt repaid yearly
                 S=3.1,
                 gammastar=0.0176, # Heat exchange coefficient between temperature layers
                 F2CO2=3.681, # Change in radiative forcing resulting from doubling of CO2
                 CATpind=588., # CO2 preindustrial concentration in atmosphere
                 CUPpind=360., # CO2 preindustrial concentration in upper layer of ocean and biosphere
                 CLOpind=1720., # CO2 preindustrial concentration in bottom layer of the ocean
                 Fexo0=0.5, # Initial value of the exogenous radiative forcing
                 Fexo1=1., # value of the exogenous radiative forcing in 2100
                 phi12=0.024, # Transfer coefficient for carbon from the atmosphere to the upper ocean
                 phi23=0.001, # Transfer coefficient for carbon from the upper ocean/biosphere to the lower ocean
                 deltapbs=-0.005, # Exogenous growth rate of the back-stop technology price
                 theta=2.6, # parameter of the abatement cost function
                 phi0=-0.292, # Constant of the short-term Philips curve
                 phi1=0.469, # Slope of the short-term Philips curve
                 kappa0=0.0397, # Constant of the investment function
                 kappa1=0.719, # Slope of the investment function
                 kappamin=0., # Minimum of the investment function
                 kappamax=0.3, # Maximum of the investment function
                 div0=0.051, # Constant of the dividend function
                 div1=0.473, # Slope of the dividend function
                 divmin=0., # Minimum of the dividend function
                 divmax=1.,# Maximum of the dividend function
                 C=49.76, # Heat capacity of the atmosphere, biosphere and upper ocean
                 C0=3.52, # Heat capacity of the deeper ocean
                 deltaEland=-0.022, # Growth rate of land-use change CO2 of emissions
                 CR0=0.17, # Constant of the leverage function
                 CRlev=.1,# slope of the leverage function
                 #Initial conditions
                 CO2AT_ini=851.,
                 CO2UP_ini=460.,
                 CO2LO_ini=1740.,
                 d_ini=1.53,
                 Eind_ini=35.85,
                 Eland_ini=2.6,
                 gsigma_ini=-0.0152,
                 pbs_ini=547.22,
                 n_ini=0.03,
                 N_ini=4.83,
                 T_ini=0.85,
                 T0_ini=0.0068,
                 Y_ini=59.74,
                 lambda_ini=0.675,
                 omega_ini=0.578,
                 r_ini=0.024,
                 # Integration parameters
                 dt=0.5, # half a year time step
                 tmax=84., # up to 2100
                 nofeedback=False):
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
        self.convCO2toC = 1./3.66667 # conversion from tCO2 to tC
        self.deltagsigma = deltagsigma # dynamics of emissivity
        self.eta = eta # relaxation parameter of inflation 
        self.etar = etar # relaxation parameter of the interest rate
        self.mu = mu # markup of prices over the average cost
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
        # Integration parameters
        self.dt = dt # half a year time step
        self.tmax = tmax # up to 2100
        self.nofeedback = nofeedback
        # Store Param class to facilitate handling / understanding
        self.params = _utils.Params_V2()

        # Initilaze the solution
        self.sol = None
        
        # Variables for plotting
        self.eco_vars = ['omega','lambda','d','N','r','piK']
        self.gdp_vars = ['Y','Y0','g','i']
        self.emi_vars = ['Eind','Eland','sigma','gsigma','n']
        self.clim_vars = ['T','T0','CO2at','CO2up','CO2lo']
        self.pcar_vars = ['pC','pbs']
        self.dam_vars = ['A','DY','DK','D']

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


    def Solve(self, plot=True, fignumber=1, verb=-1):

        def InitialConditions(omega,d,r,Y,pbs,T,n,Eind):
            """
            Returns
            sigma = ((1- DY)/Y)/((1-n)/Eind + pbs*(1-DY)/Y)
            pC = root of n = min((pC/pbs)**(1/(theta-1)),1)
            """

            D = 1. - 1./(1 + self.pi1*T
                         + self.pi2*T**2
                         + self.pi3*T**self.zeta3)
            DK = self.fK*D
            DY = 1. - (1.-D)/(1.-DK)
            deltaD = (self.delta + DK)
            
            sigma = ((1. - DY)/Y)/((1.-n)/Eind + n**self.theta/self.theta*pbs*(1.-DY)/Y)
            A = sigma*pbs*n**self.theta/self.theta
            pC = pbs*n**(self.theta - 1.)

            TC = (1-DY)*(1-A)
            i = self.eta*(self.mu*omega - 1.)
            rCB = self.Taylor(i)
            self.TCi = TC
            self.omegai = omega
            self.ri = r
            self.di = d
            self.pCi = pC
            self.sigmai = sigma
            self.ni = n
            self.sAi = self.sA
            self.Ai = A
            self.nui = self.nu
            self.deltaDi = deltaD
            piK = (TC*(1.-omega-r*d) - 0.*pC*sigma*(1-n) + self.sA*A)/self.nu - deltaD
            leverage = d*TC/self.nu
            CR = self.Tau(leverage)
            # Economic growth rate
            g = ((1.-CR)*(self.Kappa(piK)*TC/self.nu - deltaD)
                 +CR*(piK - self.Delta(piK) - TC*self.srep*d/self.nu))

            Y0 = Y/TC

            return sigma,pC,i,A,piK,D,DK,DY,g,Y0
        
        def f(t,u):
            """
            u[0] = omega, the share of work in the GDP
            u[1] = lambda, the work rate
            u[2] = d, the debt ratio
            u[3] = N, the population
            u[4] = r, the interest rate
            u[5] = T, the temperature of the atmosphere
            u[6] = T0, the temperature of the deep ocean
            u[7] = Y, the GDP
            u[8] = sigma, the emission intensity
            u[9] = gsigma, the rate of increase of emission intensity
            u[10:13] = CO2-e concentration (AT, UP, LO)
            u[13] = Eland, the land-use emissions
            u[14] = pC, the carbon price
            u[15] = pbs, the backstop technology price
            """

            omega = u[0]
            lam = u[1]
            d = u[2]
            N = u[3]
            r = u[4]
            T = u[5]
            T0 = u[6]
            Y = u[7]
            sigma = u[8]
            gsigma = u[9]
            CO2 = np.array(u[10:13])
            Eland = u[13]
            pC = u[14]
            pbs = u[15]
            
            # saturation of pC
            expo = (self.theta-1.)/self.theta
            pC = pbs*min(pC/pbs,0.95*(self.theta/(sigma*pbs))**expo)

            # Abatement cost
            n = min((pC/pbs)**(1./(self.theta-1.)),1)
            A = sigma*pbs*n**self.theta/self.theta

            # Temperature damage
            D = 1. - 1./(1 + self.pi1*T
                         + self.pi2*T**2
                         + self.pi3*T**self.zeta3)
            DK = self.fK*D
            DY = 1. - (1.-D)/(1.-DK)
            deltaD = (self.delta + DK)
            # Total cost of climate change
            TC = (1-DY)*(1-A)
            # inflation
            i = self.eta*(self.mu*omega - 1.)
            # central bank interest rate
            rCB = self.Taylor(i)
            # return on assets piK
            piK = (TC*(1.-omega-r*d) - 0.*pC*sigma*(1-n) + self.sA*A)/self.nu - deltaD
            # leverage
            leverage = d*TC/self.nu
            CR = self.Tau(leverage)
            # Economic growth rate
            g = ((1.-CR)*(self.Kappa(piK)*TC/self.nu - deltaD)
                 +CR*(piK - self.Delta(piK) - TC*self.srep*d/self.nu))
            # Population growth
            beta = self.deltaN*(1.-N/self.Nbar)
            # Temperature change
            CO2AT = CO2[0]
            Find = self.F2CO2*np.log(CO2AT/self.CATpind)/np.log(2.)
            t2016 = 1.
            t2100 = 2100. - 2016.
            F = Find + self.Fexo0 + (self.Fexo1-self.Fexo0)*(t - t2016)/(t2100 - t2016)
            Y0 = Y/TC
            Eind = Y0*sigma*(1-n)
            E = Eind + Eland
            CO2dot = np.array([E*self.convCO2toC,0,0]) + np.dot(self.phimat,CO2)

            self.Eind = Eind
            self.i = i
            self.A = A
            self.piK = piK
            self.D = D
            self.DK = DK
            self.DY = DY
            self.g = g
            self.n = n
            self.Y0 = Y0
            
            return [omega*(self.Phi(lam) - i - self.alpha), #domega/dt
                    lam*(g - self.alpha - beta), # dlam/dt
                    (-d*(g + i + CR*self.srep) # dd/dt
                     +(1.-CR)*self.Kappa(piK)
                     -self.nu*((1.-CR)*deltaD + piK - self.Delta(piK))/TC),
                    N*beta, #dN.dt
                    (rCB-r)/self.etar, # dr/dt
                    (F - self.rho*T - self.gammastar*(T-T0))/self.C, #dT/dt
                    self.gammastar*(T-T0)/self.C0, # dT0/dt
                    Y*g, #dY/dt
                    sigma*gsigma, #dsigma/dt
                    gsigma*self.deltagsigma, #dgsigma/dt
                    CO2dot[0], # dCO2/dt (AT)
                    CO2dot[1], # dCO2/dt (UP)
                    CO2dot[2], # dCO2/dt (LO)
                    Eland*self.deltaEland, #dEland/dt
                    pC*(self.apC + self.bpC/(t+t2016)), #dpC/dt
                    pbs*self.deltapbs] # dpbs/dt

        var_ini = InitialConditions(self.omega_ini,
                                    self.d_ini,
                                    self.r_ini,
                                    self.Y_ini,
                                    self.pbs_ini,
                                    self.T_ini,
                                    self.n_ini,
                                    self.Eind_ini)

        sigma_ini,pC_ini,i_ini,A_ini,piK_ini,D_ini,DK_ini,DY_ini,g_ini,Y0_ini = var_ini
        
        U_ini = [self.omega_ini,
                 self.lambda_ini,
                 self.d_ini,
                 self.N_ini,
                 self.r_ini,
                 self.T_ini,
                 self.T0_ini,
                 self.Y_ini,
                 sigma_ini,
                 self.gsigma_ini,
                 self.CO2AT_ini,
                 self.CO2UP_ini,
                 self.CO2LO_ini,
                 self.Eind_ini,
                 self.Eland_ini,
                 pC_ini,
                 self.pbs_ini,
                 i_ini,
                 A_ini,
                 piK_ini,
                 D_ini,
                 DK_ini,
                 DY_ini,
                 g_ini,
                 self.n_ini,
                 Y0_ini]

        U_ini_sys = [self.omega_ini,
                     self.lambda_ini,
                     self.d_ini,
                     self.N_ini,
                     self.r_ini,
                     self.T_ini,
                     self.T0_ini,
                     self.Y_ini,
                     sigma_ini,
                     self.gsigma_ini,
                     self.CO2AT_ini,
                     self.CO2UP_ini,
                     self.CO2LO_ini,
                     self.Eland_ini,
                     pC_ini,
                     self.pbs_ini]
        
        system = ode(f).set_integrator('dop853', verbosity=verb)
        system.set_initial_value(U_ini_sys,self.dt)
        
        nts = int(self.tmax/self.dt)
        t = np.zeros(nts)
        U = np.zeros((nts,26))
        U[0,:] = U_ini
        #
        k=0
        while system.successful() and k < nts-1:
            k=k+1
            system.integrate(system.t + self.dt)
            t[k] = system.t
            U[k,:13] = system.y[:13]
            U[k,13] = self.Eind
            U[k,14:17] = system.y[13:16]
            U[k,17] = self.i
            U[k,18] = self.A
            U[k,19] = self.piK
            U[k,20] = self.D
            U[k,21] = self.DK
            U[k,22] = self.DY
            U[k,23] = self.g
            U[k,24] = self.n
            U[k,25] = self.Y0
            
        # Keep track of the solution for later use (plot, treatment...)
        self.sol = {'t':t,'U':U}

        # Optional plotting right away after computation
        if plot:
            self.plot(fignumber=fignumber)

        return t,U

    def plot(self, vars = eco_vars+gdp_vars, fs=None, dmargin=None, draw=True, fignumber=1):
        """
        Plot the solution
        """
        
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

