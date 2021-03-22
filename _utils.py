#!/usr/bin/python

class Params(object):
    """
    Class to facilitate parameters handling (symbols, names, meaning,
    index, units..), in the form of an enhanced dictionnary  
    Necessary because of the new variable, the interest rate, in the second version
    """

    def __init__(self):
        u = [None for ii in range(0,34)]

        u[0] = {'var':'K',
                'meaning':'Capital stock',
                'symbol':r'$K$',
                'units':'prod. unit',
                'plot':'lin',
                'ylim':'',
                'index':0}

        u[1] = {'var':'N',
                'meaning':'Workforce population',
                'symbol':r'$N$',
                'units':'persons',
                'plot':'lin',
                'ylim':'',
                'index':1}

        u[2] = {'var':'D',
                'meaning':'Debt',
                'symbol':r'$D$',
                'units':'money',
                'plot':'lin',
                'ylim':'',
                'index':2}

        u[3] = {'var':'w',
                'meaning':'wages',
                'symbol':r'$w$',
                'units':'money/person',
                'plot':'lin',
                'ylim':'',
                'index':3}

        u[4] = {'var':'a',
                'meaning':'productivity',
                'symbol':r'$a$',
                'units':'(prod. unit)/person',
                'plot':'lin',
                'ylim':'',
                'index':4}

        u[5] = {'var':'p',
                'meaning':'price index',
                'symbol':r'$p$',
                'units':'money/(prod. unit)',
                'plot':'lin',
                'ylim':'',
                'index':5}

        u[6] = {'var':'Eland',
                'meaning':'Land emissions',
                'symbol':r'$E_\mathrm{land}$',
                'units':'GtCO2/year',
                'plot':'lin',
                'ylim':'',
                'index':6}

        u[7] = {'var':'sigma',
                'meaning':'emission intensity',
                'symbol':r'$\sigma$',
                'units':'tCO2/(1000(prod. unit))',
                'plot':'lin',
                'ylim':'',
                'index':7}

        u[8] = {'var':'gsigma',
                'meaning':'rate of increase of emission intensity',
                'symbol':r'$g_{\sigma}$',
                'units':'1/year',
                'plot':'lin',
                'ylim':'',
                'index':8}

        u[9] = {'var':'CO2at',
                'meaning':'CO2-eq concentration in atmosphere',
                'symbol':r'$CO_2^{AT}$',
                'units':'GtC',
                'plot':'lin',
                'ylim':'',
                 'index':9}

        u[10] = {'var':'CO2up',
                 'meaning':'CO2-eq concentration in upper ocean',
                 'symbol':r'$CO_2^{UP}$',
                 'units':'GtC',
                 'plot':'lin',
                 'ylim':'',
                 'index':10}

        u[11] = {'var':'CO2lo',
                 'meaning':'CO2-eq concentration in lower ocean',
                 'symbol':r'$CO_2^{LO}$',
                 'units':'GtC',
                 'plot':'lin',
                 'ylim':'',
                 'index':11}


        u[12] = {'var':'T',
                'meaning':'temperature of atmosphere',
                'symbol':r'$T$',
                'units':r'Celsius',
                'plot':'lin',
                'ylim':'',
                'ylim':'',
                'index':12}

        u[13] = {'var':'T0',
                'meaning':'temperature of deep ocean',
                'symbol':r'$T_0$',
                'units':'Celsius',
                'plot':'lin',
                'ylim':'',
                'index':13}

        u[14] = {'var':'pbs',
                 'meaning':'price of backstop technology',
                 'symbol':r'$p_\mathrm{bs}$',
                 'units':'(prod. unit)/tC',
                 'plot':'lin',
                 'ylim':'',
                 'index':14}

        u[15] = {'var':'pC',
                 'meaning':'carbon price',
                 'symbol':r'$p_C$',
                 'units':'(prod. Unit)/tC',
                 'plot':'lin',
                 'ylim':'',
                 'index':15}

        u[16] = {'var':'omega',
                'meaning':'share of labour in GDP',
                'symbol':r'$\omega$',
                'units':'adim.',
                'plot':'lin',
                'ylim':'0,1',
                'index':16}

        u[17] = {'var':'lambda',
                'meaning':'employment rate',
                'symbol':r'$\lambda$',
                'units':'adim.',
                'plot':'lin',
                'ylim':'0,1',
                'index':17}

        u[18] = {'var':'d',
                'meaning':'private debt ratio',
                'symbol':r'$d$',
                'units':'adim.',
                'plot':'lin',
                'ylim':'',
                'index':18}

        u[19] = {'var':'Y0',
                 'meaning':'output before abatement and damage',
                 'symbol':r'$Y_0$',
                 'units':'prod. unit',
                 'plot':'lin',
                 'ylim':'',
                 'index':19}

        u[20] = {'var':'Y',
                'meaning':'GDP',
                'symbol':r'$Y$',
                'units':'(prod. unit)/year',
                'plot':'lin',
                'ylim':'',
                'index':20}

        u[21] = {'var':'g0',
                 'meaning':'real GDP growth rate of Y0',
                 'symbol':r'$g_0$',
                 'units':'1/year',
                 'plot':'lin',
                 'ylim':'',
                 'index':21}

        u[22] = {'var':'Eind',
                 'meaning':'industrial CO2-eq emissions',
                 'symbol':r'$E_\mathrm{ind}$',
                 'units':'GtCO2/year',
                 'plot':'lin',
                 'ylim':'',
                 'index':22}

        u[23] = {'var':'i',
                 'meaning':'inflation',
                 'symbol':r'$i$',
                 'units':'1./year',
                 'plot':'lin',
                 'ylim':'',
                 'index':23}

        u[24] = {'var':'A',
                 'meaning':'abatement',
                 'symbol':r'$A$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':24}

        u[25] = {'var':'n',
                 'meaning':'reduction emission fraction',
                 'symbol':r'$n$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':25}
        
        u[26] = {'var':'pi',
                 'meaning':'profit share',
                 'symbol':r'$\pi$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':26}

        u[27] = {'var':'piK',
                 'meaning':'return on assets',
                 'symbol':r'$\pi_K$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':27}

        u[28] = {'var':'Dam',
                 'meaning':'damage',
                 'symbol':r'$\mathrm{Dam}$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':28}

        u[29] = {'var':'DamK',
                 'meaning':'capital damage',
                 'symbol':r'$\mathrm{Dam}_K$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':29}

        u[30] = {'var':'DamY',
                 'meaning':'output damage',
                 'symbol':r'$\mathrm{Dam}_Y$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':30}

        u[31] = {'var':'Fexo',
                'meaning':'exogenous radiative forcing',
                'symbol':r'$F_\mathrm{exo}$',
                'units':'W/m^2',
                'plot':'lin',
                'ylim':'',
                'index':31}

        u[32] = {'var':'Find',
                'meaning':'Industrial radiative forcing',
                'symbol':r'$F_\mathrm{ind}$',
                'units':'W/m^2',
                'plot':'lin',
                'ylim':'',
                'index':32}

        u[33] = {'var':'rCB',
                 'meaning':'central bank interest rate',
                 'symbol':r'$r_\mathrm{CB}$',
                 'units':'1./year',
                 'plot':'lin',
                 'ylim':'',
                 'index':33}

        self._lparams = u
        self._version = 'v2'

    @property
    def lparams(self):
        return self._lparams

    @property
    def version(self):
        return self._version

    def __getitem__(self, ind):
        return self.lparams[ind]

    def get(self, key=None, val=None, unique=True):
        assert not (key is None and val is not None)
        if key is None and val is None:
            out = self.lparams
        
        elif val is None:
            out = [p[key] for p in self.lparams]

        else:
            out = [p for p in self.lparams if p[key]==val]
        return out
