#!/usr/bin/python


class Params(object):
    """ Class to facilitate parameters handling (symbols, names, meaning,
    index, units..), in the form of an enhanced dictionnary  """

    def __init__(self):
        u = [None for ii in range(0,15)]
        u[0] = {'var':'omega',
                'meaning':'share of labour in GDP',
                'symbol':r'$\omega$',
                'units':'adim.'}

        u[1] = {'var':'llambda',
                'meaning':'employment rate',
                'symbol':r'$\lambda$',
                'units':'adim.'}

        u[2] = {'var':'d',
                'meaning':'private debt ratio',
                'symbol':r'$d$',
                'units':'adim.'}

        u[3] = {'var':'N',
                'meaning':'total population',
                'symbol':r'$N$',
                'units':'person'}

        u[4] = {'var':'T',
                'meaning':'temperature of atmosphere',
                'symbol':r'$T$',
                'units':r'Celsius'}

        u[5] = {'var':'T0',
                'meaning':'temperature of deep ocean',
                'symbol':r'$T_0$',
                'units':'Celsius'}

        u[6] = {'var':'Y',
                'meaning':'GDP',
                'symbol':r'$Y$',
                'units':'money'}

        u[7] = {'var':'sigma',
                'meaning':'emission intensity',
                'symbol':r'$\sigma$',
                'units':'GtC/year/money'}

        u[8] = {'var':'gsigma',
                'meaning':'rate of increase of emission intensity',
                'symbol':r'$g_{\sigma}$',
                'units':'1/year'}

        u[9] = {'var':'',
                 'meaning':'CO2-eq concentration in atmosphere',
                 'symbol':r'$CO_2^{AT}$',
                 'units':'GtC'}

        u[10] = {'var':'',
                 'meaning':'CO2-eq concentration in upper ocean',
                 'symbol':r'$CO_2^{UP}$',
                 'units':'GtC'}

        u[11] = {'var':'',
                 'meaning':'CO2-eq concentration in lower ocean',
                 'symbol':r'$CO_2^{LO}$',
                 'units':'GtC'}

        u[12] = {'var':'Eland',
                 'meaning':'land_use CO2-eq emissions',
                 'symbol':r'$E_{land}$',
                 'units':'GtC/year'}

        u[13] = {'var':'pC',
                'meaning':'carbon price',
                'symbol':r'$p_C$',
                'units':'money'}

        u[14] = {'var':'pbs',
                 'meaning':'price of backstop technology',
                 'symbol':r'$p_{bs}$',
                 'units':'money'}

        self._lparams = u

    @property
    def lparams(self):
        return self._lparams

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
