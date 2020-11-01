''' Written by Eric Deal, 10.2020. Script should be run using python 3.'''

import matplotlib.pyplot as plt, numpy as np, scipy.special as ss

class glacial_profile_model(object):
    '''Script to accompany the paper 'The Sliding Ice Incision Model: A new approach to understanding glacial landscape evolution' be Eric Deal
    and Guenther Prasicek. The parameter names here should match those given in the paper.'''

    def __init__(self,
        # system size and boundary conditions
        U=1e-3, # uplift rate in m/yr
        L=50e3, # profile length in m
        # climatic forcing, P and beta must be specified for the erodability coefficients, either ELA is specified, or To is specified. 
        # If To is specified then ELA is determined by the combination of P, beta and To, if ELA is specified, then To is adjusted to give the
        # chosen ELA with the given P and beta.
        delta=2e-3, # solid preciipitation lapse rate in 1/yr
        P=1, # mean rainfall rate in m/yr
        To=None, # Elevation where 100% of precipitation is converted to glacial ice in m
        ELA=None, # Equilibrium Line Altitude (elevation, h, where P - delta*(To-h) = 0) in m
        # fluvial parameters
        K=1e-6, # fluvial erodibility coefficient (units depend on n and m)
        fluvial_slope_exp=1, # fluvial slope exponent, n
        # glacial parameters
        ce=1e-4, # glacial erodibility coefficient (units depend on ell)
        ice_sliding_exp=7/9, # ice sliding exponent, ell
        # plot profile or not
        plot=False):

        if plot is True: plt.figure(1, figsize=(10,6))

        ### set variables and parameters (all units in meters and years unless otherwise specified)
        # number of points to calculate solution on between 0 and L
        self.N = 300

        # physical constraints
        self.U, self.L, self.P, self.delta = (U, L, P, delta)
        # if neither ELA or To is user defined
        if ELA is None and To is None: self.ELA, self.To = 1000, 1000 + self.P/self.delta
        # if To is user defined, but not ELA
        elif ELA is None and To is not None: self.ELA, self.To = To - self.P/self.delta, To
        # if ELA is user defined, then To is ignored, whether it was user defined or not.
        else: self.ELA, self.To = ELA, ELA + self.P/self.delta
        
        # fluvial exponents (the model requires that h*m/n=1)
        self.h_exp = 2 # fluvial hack's exponent
        self.omega_f = 1/2 # fluvial channel width exponent (just for plotting, value does not influence model)
        self.n = fluvial_slope_exp # slope exponent in SPIM
        self.m = self.n/2 # flux exponent in SPIM

        # fluvial constants
        self.xo = 1000 # critical catchment area where hillslopes take over, max elevation is assumed to be proportional to fluvial elevation at xo
        self.kc = 2 # fluvial Hack's coefficient
        self.kw = 1 # fluvial channel width coefficient (just for plotting, value does not influence model)
        self.K = K # fluvial erosion coefficient
        self.ks = (self.U / ((self.kc * self.P)**self.m * self.K))**(1/self.n) # fluvial steepness index

        # fluvial profile solution
        self.zf = lambda x: self.ks*np.log(self.L/x)
        self.Wf = lambda x: self.kw * x**self.omega_f # fluvial channel width (m)
        self.qf = lambda x: (self.P * self.kc * x**self.h_exp)**(1-self.omega_f) / self.kw # fluvial water flux density (m^2/yr)

        # glacial exponents
        self.eta = 2 # glacial Hack’s exponent
        self.phi = .5 # network branching parameter
        self.omega_g = 1/4 # glacial channel width exponent
        self.ell = ice_sliding_exp # glacial sliding erosion exponent
        self.psi = 3. # Glen’s flow law exponent
        self.gamma = 2/3 # ice flux approximation exponent
        self.kappa = self.psi / (1 + self.gamma*(self.psi - 1))
        self.lamb = self.gamma*(self.psi - 1) / (1 + self.gamma*(self.psi - 1))
        self.nu = self.ell * self.kappa # slope exponent in SIIM
        self.mu = self.ell * self.lamb * (1 - self.omega_g) # flux exponent in SIIM
        self.c1 = 1 - self.eta * self.mu / self.nu # ice accumulation exponent
        self.c2 = 1 - self.mu / self.nu # ice loss exponent
        self.r = (self.nu - self.eta * self.mu) / (self.nu + self.mu)
        # ice flow constants
        self.ro, self.g, self.cd, self.cb = (910.0, 9.81, 2.39e-24, 1.7e-19)
        self.fs, self.fd = (self.ro*self.g)**self.psi * self.cb * 31556926.0, (self.ro*self.g)**self.psi * self.cd * 31556926.0
        # glacial constants
        self.D = (1 - self.phi) * self.eta # divergence operator
        self.sigma = self.D / (1 + self.D) # ratio of relief above ELA to glacial relief
        self.cc = 2 # glacial hack's coefficient
        self.ca = 1/2 # ice flux approximation coefficient
        self.cw = 50 # glacial channel width coefficient
        self.ce = ce # ice sliding erosion coefficient
        self.cp = self.ca * (self.cb/self.cd)**((1-self.gamma)/(self.psi-1))
        self.ct = (self.cp * self.fs**(1/(self.psi-1)))**(self.lamb/self.gamma)
        self.C = self.ce*(self.ct/self.cw**self.lamb)**self.ell
        self.cs = (self.U * ss.beta(self.c1, self.c2)**self.nu / (self.C * (self.sigma*self.delta*self.cc)**self.mu)) ** (1 / (self.mu + self.nu))
        self.ck = (self.ks / self.cs) / (self.r * (self.sigma - 1))

        # Pure fluvial profile (this is the case if the max topography does not exceed the ELA)
        if self.ks < self.ELA / np.log(self.L/self.xo):
            self.xf = np.linspace(self.xo, self.L, self.N) # x_axis for plotting fluvial profile
            self.Hmax = 0 # maximum ice thickness
            self.xt = 0 # location of glacial terminus
            self.ho = self.zf(self.xo) # maximum profile elevation

            if plot == True:
                plt.plot(self.xf/1000, self.zf(self.xf), 'k--', lw=1, label='Fluvial profile')

        # Purely glacial
        elif self.cs >= ((1+self.D)*self.ELA / self.L**self.r):
            self.t = self.nu / (self.nu + self.mu)
            self.theta = self.mu / (self.nu + self.mu)

            # function to help calculate max elevation iteratively
            def calc_ho(hinput):
                xt = (1 - self.ELA / hinput) / self.sigma * self.L
                return self.cs * self.L**self.theta * xt**(self.r - self.theta) * (ss.betainc(self.c1, self.c2, self.L/xt) * ss.beta(self.c1, self.c2))**self.t

            # calculate max elevation iteratively
            hin, count, diff_h = (10 * self.ELA * (1 + self.D), 0, 1e4)
            while (np.abs(diff_h) > .1) and (count < 100):
                diff_h = hin - calc_ho(hin)
                hin -= .3*diff_h
                count += 1

            if count >= 99: print('warning, did not converge')

            self.xt = (1 - self.ELA / hin) / self.sigma * self.L # location of glacial terminus along profile (m)
            self.ho = self.cs * self.L**self.theta * self.xt**(self.r - self.theta) * ss.betainc(self.c1, self.c2, self.L/self.xt) # maximum profile elevation (m)

            # glacial profile solution
            self.h = lambda x: self.cs * self.L**self.theta * self.xt**(self.r-self.theta) * ss.betainc(self.c1, self.c2, self.L/self.xt) * (1 - ss.betainc(self.c1, self.c2, x/self.xt)/ss.betainc(self.c1, self.c2, self.L/self.xt))
            # mean elevation of glacier over glacial profile
            self.Sgmean = self.ho / self.xt
            self.I = lambda x: self.sigma * self.delta * self.Sgmean * (self.xt - x) # Effective ice accumulation rate (m/yr)
            self.qg = lambda x: (self.I(x) * self.cc * x**self.eta)**(1-self.omega_g) / self.cw # Ice flux density (m^2/yr)
            self.H = lambda x: self.cp * (self.ce / self.U)**(self.gamma/self.ell) * self.qg(x) ** self.gamma # Ice thickness at steady state (m)
            self.Wg = lambda x: self.cw * x**self.omega_g # glacial channel width (m)
            self.zg = lambda x: self.h(x) - self.H(x) # subglacial bedrock elevation (m)
            self.xg = np.linspace(1, self.L, self.N) # glacial profile x-axis (m)
            self.Hmax = np.max(self.H(self.xg)) # maximum ice thickness

            if plot == True:
                plt.plot(self.xg/1000, self.zf(self.xg), ':', c='gray', lw=1, label='Fluvial profile cont.')
                plt.plot(self.xg/1000, self.zg(self.xg), c='k', lw=1, label='Glacial profile')
                plt.plot(self.xg/1000, self.h(self.xg), c='c', lw=1, label='Ice surface')
                plt.fill_between(self.xg/1000, self.zg(self.xg), self.h(self.xg), alpha=.3, color='c')
                
        # mixed glacial-fluvial
        else:
            self.xt = (self.ck * np.real(ss.lambertw((self.L * np.exp(-self.ELA / self.ks))**self.r / self.ck)))**(1/self.r) # glacial terminus location along the profile

            # glacial profile solution
            self.ho = self.ELA + self.cs * self.sigma * self.xt**self.r # maximum elevation
            self.Sgmean = (self.ho - self.zf(self.xt)) / self.xt # mean slope over glacier
            self.h = lambda x: self.ELA + self.cs * self.xt**self.r * (self.sigma - ss.betainc(self.c1, self.c2, x/self.xt)) # Ice surface elevation as a function of x
            self.I = lambda x: self.sigma * self.delta * self.Sgmean * (self.xt - x) # Effective ice accumulation rate (m/yr)
            self.qg = lambda x: (self.I(x) * self.cc * x**self.eta)**(1-self.omega_g) / self.cw # Ice flux density (m^2/yr)
            self.H = lambda x: self.cp * (self.ce / self.U)**(self.gamma/self.ell) * self.qg(x) ** self.gamma # Ice thickness at steady state (m)
            self.Wg = lambda x: self.cw * x**self.omega_g # glacial channel width (m)
            self.zg = lambda x: self.h(x) - self.H(x) # subglacial bedrock elevation (m)
            self.xf = np.linspace(self.xt, self.L, self.N) # fluvial profile x-axis (m)
            self.xg = np.linspace(1, self.xt, self.N) # glacial profile x-axis (m)
            self.Hmax = np.max(self.H(self.xg)) # maximum ice thickness

            if plot == True:
                plt.plot(self.xf/1000, self.zf(self.xf), '-', c='gray', lw=1, label='Fluvial profile')
                plt.plot(self.xg/1000, self.zf(self.xg), ':', c='gray', lw=1, label='Fluvial profile cont.')
                plt.plot(self.xg/1000, self.zg(self.xg), c='k', lw=1, label='Glacial profile')
                plt.plot(self.xg/1000, self.h(self.xg), c='c', lw=1, label='Ice surface')
                plt.fill_between(self.xg/1000, self.zg(self.xg), self.h(self.xg), alpha=.3, color='c')

        if plot is True:
            plt.plot((0, self.L/1e3), (self.ELA, self.ELA), 'r:', lw=0.5, label='ELA')
            plt.xlabel('Along profle distance (km)')
            plt.ylabel('Elevation above baselevel')
            plt.xlim(0,self.L/1e3)
            plt.ylim(0,self.ho+100)
            plt.legend(ncol=2)


gp = glacial_profile_model(U=1e-3, ELA=1000, P=3, delta=3e-3, plot=True)
print('Maximum profile elevation: %i m\nMaximum ice thickness: %i m\nGlacial terminus: %0.1f km' % (gp.ho, gp.Hmax, gp.xt/1e3))

plt.figure(2)
plt.plot(gp.xg/1e3, gp.qg(gp.xg), 'b', label='Ice flux density')
plt.loglog(gp.xf/1e3, gp.qf(gp.xf), 'b--', label='Water flux density')
plt.xlim(0, gp.L/1e3)
plt.ylabel('Flux density ($km^2/yr$)')
plt.xlabel('Along profile distance (km)')
plt.legend()

plt.figure(3)
plt.plot(gp.xg/1e3, gp.Wg(gp.xg), 'g', label='Glacial channel width (m)')
plt.loglog(gp.xf/1e3, gp.Wf(gp.xf), 'g--', label='Fluvial channel width (m)')
plt.xlim(0, gp.L/1e3)
plt.xlabel('Along profile distance (km)')
plt.ylabel('Channel thickness')
plt.legend()

plt.figure(4)
plt.plot(gp.xg/1e3, gp.H(gp.xg), 'c', label='Ice thickness')
plt.xlim(0, gp.xt/1e3 + 1)
plt.xlabel('Along profile distance (km)')
plt.ylabel('Ice thicknes (m)')
