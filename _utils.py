#!/usr/bin/python

class Params(object):
    """ Class to facilitate parameters handling (symbols, names, meaning,
    index, units..), in the form of an enhanced dictionnary  """

    def __init__(self):
        u = [None for ii in range(0,25)]
        u[0] = {'var':'omega',
                'meaning':'share of labour in GDP',
                'symbol':r'$\omega$',
                'units':'adim.',
                'plot':'lin',
                'ylim':'0,1',
                'index':0}

        u[1] = {'var':'lambda',
                'meaning':'employment rate',
                'symbol':r'$\lambda$',
                'units':'adim.',
                'plot':'lin',
                'ylim':'0,1',
                'index':1}

        u[2] = {'var':'d',
                'meaning':'private debt ratio',
                'symbol':r'$d$',
                'units':'adim.',
                'plot':'lin',
                'ylim':'',
                'index':2}

        u[3] = {'var':'N',
                'meaning':'total working force',
                'symbol':r'$N$',
                'units':'person',
                'plot':'lin',
                'ylim':'',
                'index':3}

        u[4] = {'var':'T',
                'meaning':'temperature of atmosphere',
                'symbol':r'$T$',
                'units':r'Celsius',
                'plot':'lin',
                'ylim':'',
                'ylim':'',
                'index':4}

        u[5] = {'var':'T0',
                'meaning':'temperature of deep ocean',
                'symbol':r'$T_0$',
                'units':'Celsius',
                'plot':'lin',
                'ylim':'',
                'index':5}

        u[6] = {'var':'Y',
                'meaning':'GDP',
                'symbol':r'$Y$',
                'units':'(prod. unit)/year',
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

        u[12] = {'var':'Eind',
                 'meaning':'industrial CO2-eq emissions',
                 'symbol':r'$E_\mathrm{ind}$',
                 'units':'GtCO2/year',
                 'plot':'lin',
                 'ylim':'',
                 'index':12}

        u[13] = {'var':'Eland',
                 'meaning':'land use CO2-eq emissions',
                 'symbol':r'$E_\mathrm{land}$',
                 'units':'GtCO2/year',
                 'plot':'lin',
                 'ylim':'',
                 'index':13}

        u[14] = {'var':'pC',
                 'meaning':'carbon price',
                 'symbol':r'$p_C$',
                 'units':'(prod. Unit)/tC',
                 'plot':'lin',
                 'ylim':'',
                 'index':14}

        u[15] = {'var':'pbs',
                 'meaning':'price of backstop technology',
                 'symbol':r'$p_\mathrm{bs}$',
                 'units':'(prod. unit)/tC',
                 'plot':'lin',
                 'ylim':'',
                 'index':15}

        u[16] = {'var':'i',
                 'meaning':'inflation',
                 'symbol':r'$i$',
                 'units':'1./year',
                 'plot':'lin',
                 'ylim':'',
                 'index':16}

        u[17] = {'var':'A',
                 'meaning':'abatement',
                 'symbol':r'$A$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':17}

        u[18] = {'var':'pi',
                 'meaning':'profit share',
                 'symbol':r'$\pi$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':18}

        u[19] = {'var':'D',
                 'meaning':'damage',
                 'symbol':r'$D$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':19}

        u[20] = {'var':'DK',
                 'meaning':'capital damage',
                 'symbol':r'$D^K$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':20}

        u[21] = {'var':'DY',
                 'meaning':'output damage',
                 'symbol':r'$D^Y$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':21}
        
        u[22] = {'var':'g',
                 'meaning':'real GDP growth rate',
                 'symbol':r'$g$',
                 'units':'1/year',
                 'plot':'lin',
                 'ylim':'',
                 'index':22}

        u[23] = {'var':'n',
                 'meaning':'reduction emission fraction',
                 'symbol':r'$n$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':23}
        
        u[24] = {'var':'Y0',
                 'meaning':'output before abatement and damage',
                 'symbol':r'$Y_0$',
                 'units':'prod. unit',
                 'plot':'lin',
                 'ylim':'',
                 'index':24}


        self._lparams = u
        self._version = 'v1'

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


class Params_V2(object):
    """ Class to facilitate parameters handling (symbols, names, meaning,
    index, units..), in the form of an enhanced dictionnary  
    Necessary because of the new variable, the interest rate, in the second version
    """

    def __init__(self):
        u = [None for ii in range(0,27)]
        u[0] = {'var':'omega',
                'meaning':'share of labour in GDP',
                'symbol':r'$\omega$',
                'units':'adim.',
                'plot':'lin',
                'ylim':'0,1',
                'index':0}

        u[1] = {'var':'lambda',
                'meaning':'employment rate',
                'symbol':r'$\lambda$',
                'units':'adim.',
                'plot':'lin',
                'ylim':'0,1',
                'index':1}

        u[2] = {'var':'d',
                'meaning':'private debt ratio',
                'symbol':r'$d$',
                'units':'adim.',
                'plot':'lin',
                'ylim':'',
                'index':2}

        u[3] = {'var':'N',
                'meaning':'total working force',
                'symbol':r'$N$',
                'units':'person',
                'plot':'lin',
                'ylim':'',
                'index':3}

        u[4] = {'var':'r',
                'meaning':'interest rate',
                'symbol':r'$r$',
                'units':'1./year',
                'plot':'lin',
                'ylim':'',
                'index':4}

        u[5] = {'var':'T',
                'meaning':'temperature of atmosphere',
                'symbol':r'$T$',
                'units':r'Celsius',
                'plot':'lin',
                'ylim':'',
                'ylim':'',
                'index':5}

        u[6] = {'var':'T0',
                'meaning':'temperature of deep ocean',
                'symbol':r'$T_0$',
                'units':'Celsius',
                'plot':'lin',
                'ylim':'',
                'index':6}

        u[7] = {'var':'Y',
                'meaning':'GDP',
                'symbol':r'$Y$',
                'units':'(prod. unit)/year',
                'plot':'lin',
                'ylim':'',
                'index':7}

        u[8] = {'var':'sigma',
                'meaning':'emission intensity',
                'symbol':r'$\sigma$',
                'units':'tCO2/(1000(prod. unit))',
                'plot':'lin',
                'ylim':'',
                'index':8}

        u[9] = {'var':'gsigma',
                'meaning':'rate of increase of emission intensity',
                'symbol':r'$g_{\sigma}$',
                'units':'1/year',
                'plot':'lin',
                'ylim':'',
                'index':9}

        u[10] = {'var':'CO2at',
                'meaning':'CO2-eq concentration in atmosphere',
                'symbol':r'$CO_2^{AT}$',
                'units':'GtC',
                'plot':'lin',
                'ylim':'',
                 'index':10}

        u[11] = {'var':'CO2up',
                 'meaning':'CO2-eq concentration in upper ocean',
                 'symbol':r'$CO_2^{UP}$',
                 'units':'GtC',
                 'plot':'lin',
                 'ylim':'',
                 'index':11}

        u[12] = {'var':'CO2lo',
                 'meaning':'CO2-eq concentration in lower ocean',
                 'symbol':r'$CO_2^{LO}$',
                 'units':'GtC',
                 'plot':'lin',
                 'ylim':'',
                 'index':12}

        u[13] = {'var':'Eind',
                 'meaning':'industrial CO2-eq emissions',
                 'symbol':r'$E_\mathrm{ind}$',
                 'units':'GtCO2/year',
                 'plot':'lin',
                 'ylim':'',
                 'index':13}

        u[14] = {'var':'Eland',
                 'meaning':'land use CO2-eq emissions',
                 'symbol':r'$E_\mathrm{land}$',
                 'units':'GtCO2/year',
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

        u[16] = {'var':'pbs',
                 'meaning':'price of backstop technology',
                 'symbol':r'$p_\mathrm{bs}$',
                 'units':'(prod. unit)/tC',
                 'plot':'lin',
                 'ylim':'',
                 'index':16}

        u[17] = {'var':'i',
                 'meaning':'inflation',
                 'symbol':r'$i$',
                 'units':'1./year',
                 'plot':'lin',
                 'ylim':'',
                 'index':17}

        u[18] = {'var':'A',
                 'meaning':'abatement',
                 'symbol':r'$A$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':18}

        u[19] = {'var':'piK',
                 'meaning':'return on assets',
                 'symbol':r'$\pi_K$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':19}

        u[20] = {'var':'D',
                 'meaning':'damage',
                 'symbol':r'$D$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':20}

        u[21] = {'var':'DK',
                 'meaning':'capital damage',
                 'symbol':r'$D_K$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':21}

        u[22] = {'var':'DY',
                 'meaning':'output damage',
                 'symbol':r'$D_Y$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':22}
        
        u[23] = {'var':'g',
                 'meaning':'real GDP growth rate',
                 'symbol':r'$g$',
                 'units':'1/year',
                 'plot':'lin',
                 'ylim':'',
                 'index':23}

        u[24] = {'var':'n',
                 'meaning':'reduction emission fraction',
                 'symbol':r'$n$',
                 'units':'adim.',
                 'plot':'lin',
                 'ylim':'',
                 'index':24}
        
        u[25] = {'var':'Y0',
                 'meaning':'output before abatement and damage',
                 'symbol':r'$Y_0$',
                 'units':'prod. unit',
                 'plot':'lin',
                 'ylim':'',
                 'index':25}

        u[26] = {'var':'rCB',
                'meaning':'central bank interest rate',
                'symbol':r'$r_\mathrm{CB}$',
                'units':'1./year',
                'plot':'lin',
                'ylim':'',
                'index':26}

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
