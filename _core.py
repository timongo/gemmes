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
    import gemmes._utils as _utils
    import gemmes._plot as _plot


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
                 deltagsigma=-0.001,
                 eta=0.5,
                 m=1.875,
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
                 C=55,
                 C0=3.52,
                 #Initial conditions
                 CO2AT_ini=851.,
                 CO2UP_ini=460.,
                 CO2LO_ini=1740.,
                 d_ini=1.53,
                 deltaEland=-0.022,
                 Eind_ini=38.85,
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
                 tmax=84.): # up to 2100

        """
        Initialization
        """
        self.name = name # Name of the run, useful for plotting multiple runs

        self.alpha = alpha # Constant growth rate of labor productivity
        self.q = q # Speed of growth of the workforce dynamics
        self.PN = PN # Upper limit of the workforce dynamics in billions
        self.nu = nu # Constant capital-to-output ratio
        self.pi1 = pi1 # Damage function parameter
        self.pi2 = pi2 # Damage function parameter
        self.pi3 = pi3 # Damage function parameter
        self.zeta3 = zeta3 # Damage function parameter
        self.fK = fK # Fraction of environmental damage allocated to the stock of capital
        self.delta = delta # Depreciation rate of capital
        self.r = r # Short-term interest rate of the economy
        self.apC = apC # Carbon price evolution parameter
        self.bpC = bpC # Carbon price evolution parameter
        self.deltagsigma = deltagsigma # Variation rate of the growth of emission intensity
        self.eta = eta # Relaxation parameter of the inflation
        self.m = m # "markup"
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
        self.n_ini = n_ini
        self.N_ini = N_ini
        self.T_ini = T_ini
        self.T0_ini = T0_ini
        self.Y_ini = Y_ini
        self.lambda_ini = lambda_ini
        self.omega_ini = omega_ini

        self.dt = dt
        self.tmax = tmax

        # Store Param class to facailitate handling / understanding
        self.params = _utils.Params()

        # Initilaze the solution
        self.sol = None
        
    def Phi(self,x):
        """
        Philips curve
        """
        return self.phi0 + self.phi1*x

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
            div[div>self.divmax] = self.kappamax
            return div

    def Solve(self, plot=True, verb=-1):

        def InitialConditions(Y,pbs,T,n,Eind):
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
            
            sigma = ((1. - DY)/Y)/((1.-n)/Eind + n**self.theta/self.theta*pbs*(1.-DY)/Y)
            pC = pbs*n**(self.theta - 1.)

            return sigma,pC
            
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
            
            # Abaatement cost
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
            TotalCost = (1-DY)*(1-A)
            # Profit rate
            pi = (1. - omega - self.r*d
                  -(pC*sigma*(1-n) + deltaD*self.nu)/TotalCost)
            # Inflation
            c = omega + self.r*d + self.Delta(pi) + self.nu*deltaD/TotalCost
            i = self.eta*(self.m*c - 1.)
            # Economic growth rate
            g = self.Kappa(pi)*TotalCost/self.nu - deltaD
            # Population growth
            beta = self.q*(1-N/self.PN)
            # Temperature change
            CO2AT = CO2[0]
            Find = self.F2CO2*np.log(CO2AT/self.CATpind)/np.log(2.)
            t2016 = 1.
            t2100 = 2100 - 2016
            F = Find + self.Fexo0 + (self.Fexo1-self.Fexo0)*(t - t2016)/(t2100 - t2016)
            Y0 = Y/((1.-DY)*(1.-A))
            Eind = Y0*sigma*(1-n)
            E = Eind + Eland
            CO2dot = np.array([E,0,0]) + np.dot(self.phimat,CO2)

            self.other = (i,g,pi)
            
            return [omega*(self.Phi(lam) - i - self.alpha),
                    lam*(g - self.alpha - beta),
                    (-d*(g + i)
                     +self.Kappa(pi) + self.Delta(pi) - pi
                     -self.nu*deltaD/TotalCost),
                    N*beta,
                    (F - self.rho*T - self.gammastar*(T-T0))/self.C,
                    self.gammastar*(T-T0)/self.C0,
                    Y*g,
                    sigma*gsigma,
                    gsigma*self.deltagsigma,
                    CO2dot[0],
                    CO2dot[1],
                    CO2dot[2],
                    Eland*self.deltaEland,
                    pC*(self.apC + self.bpC/t),
                    pbs*self.deltapbs]


        sigma_ini,pC_ini = InitialConditions(self.Y_ini,
                                             self.pbs_ini,
                                             self.T_ini,
                                             self.n_ini,
                                             self.Eind_ini)
        
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
                 self.Eland_ini,
                 pC_ini,
                 self.pbs_ini]
        
        system = ode(f).set_integrator('dop853', verbosity=verb)
        system.set_initial_value(U_ini,self.dt)

        t = [self.dt]
        U = [U_ini]
        dudt = f(self.dt,U_ini)
        other = [self.other]
        while system.successful() and system.t < self.tmax:
            system.integrate(system.t + self.dt)
            t.append(system.t)
            U.append(system.y)
            other.append(self.other)
            
        U = np.array(U)
        other = np.array(other)
        t = np.array(t)
        
        # Keep track of the solution for later use (plot, treatment...)
        self.sol = {'U':U, 't':t, 'other':other}

        # Optional plotting right away after computation
        if plot:
            self.plot()

        return t,U,other

    def plot(self, Type='basic', fs=None, dmargin=None, draw=True):
        """ Plot the solution """
        assert self.sol is not None, "System not solved yet ! "
        if Type=='basic':
            lax = _plot.plot_basic(self, fs=fs, dmargin=dmargin, draw=draw)
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
