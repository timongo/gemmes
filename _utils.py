#!/usr/bin/python


class Params(object):
    """ Class to facilitate parameters handling (symbols, names, meaning,
    index, units..), in the form of an enhanced dictionnary  """

    def __init__(self):
        u = [None for ii in range(0,15)]
        u[0] = {'var':'omega',
                'meaning':'share of labour in GDP',
                'symbol':r'$\omega$',
                'units':'adim.',
                'plotmethod':'lin'}

        u[1] = {'var':'llambda',
                'meaning':'employment rate',
                'symbol':r'$\lambda$',
                'units':'adim.',
                'plotmethod':'lin'}

        u[2] = {'var':'d',
                'meaning':'private debt ratio',
                'symbol':r'$d$',
                'units':'adim.',
                'plotmethod':'lin'}

        u[3] = {'var':'N',
                'meaning':'total population',
                'symbol':r'$N$',
                'units':'person',
                'plotmethod':'lin'}

        u[4] = {'var':'T',
                'meaning':'temperature of atmosphere',
                'symbol':r'$T$',
                'units':r'Celsius',
                'plotmethod':'lin'}

        u[5] = {'var':'T0',
                'meaning':'temperature of deep ocean',
                'symbol':r'$T_0$',
                'units':'Celsius',
                'plotmethod':'lin'}

        u[6] = {'var':'Y',
                'meaning':'GDP',
                'symbol':r'$Y$',
                'units':'money',
                'plotmethod':'log'}

        u[7] = {'var':'sigma',
                'meaning':'emission intensity',
                'symbol':r'$\sigma$',
                'units':'GtC/year/money',
                'plotmethod':'lin'}

        u[8] = {'var':'gsigma',
                'meaning':'rate of increase of emission intensity',
                'symbol':r'$g_{\sigma}$',
                'units':'1/year',
                'plotmethod':'lin'}

        u[9] = {'var':'',
                 'meaning':'CO2-eq concentration in atmosphere',
                 'symbol':r'$CO_2^{AT}$',
                 'units':'GtC',
                'plotmethod':'log'}

        u[10] = {'var':'',
                 'meaning':'CO2-eq concentration in upper ocean',
                 'symbol':r'$CO_2^{UP}$',
                 'units':'GtC',
                 'plotmethod':'log'}

        u[11] = {'var':'',
                 'meaning':'CO2-eq concentration in lower ocean',
                 'symbol':r'$CO_2^{LO}$',
                 'units':'GtC',
                 'plotmethod':'log'}

        u[12] = {'var':'Eland',
                 'meaning':'land use CO2-eq emissions',
                 'symbol':r'$E_{land}$',
                 'units':'GtC/year',
                 'plotmethod':'lin'}

        u[13] = {'var':'Eind',
                 'meaning':'industrial CO2-eq emissions',
                 'symbol':r'$E_{ind}$',
                 'units':'GtC/year',
                 'plotmethod':'lin'}

        u[14] = {'var':'pC',
                'meaning':'carbon price',
                'symbol':r'$p_C$',
                'units':'money',
                'plotmethod':'lin'}

        u[15] = {'var':'pbs',
                 'meaning':'price of backstop technology',
                 'symbol':r'$p_{bs}$',
                 'units':'money',
                 'plotmethod':'lin'}

        u[16] = {'var':'pbs',
                 'meaning':'price of backstop technology',
                 'symbol':r'$p_{bs}$',
                 'units':'money',
                 'plotmethod':'lin'}

        u[17] = {'var':'i',
                 'meaning':'inflation',
                 'symbol':r'$i$',
                 'units':'adim.',
                 'plotmethod':'lin'}



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
