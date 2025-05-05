# Parameters structure
# Edit values here
class paramsStruct:
    def __init__(self, **kwargs):
        ## GEOMETRY AND BOUNDARY CONDITIONS
        self.Minf = kwargs.get('Minf', 0.8)
        self.isWall = kwargs.get('isWall', 0)
        self.foil = kwargs.get('foil', 'NACA0010')

        ## FLUID PARAMETERS
        self.gamma = kwargs.get('gamma', 1.4)
        
        ## MESH PARAMETERS
        self.kConst = kwargs.get('kConst', 3)
        self.ds = kwargs.get('ds', 0.004)
        self.jMax = kwargs.get('jMax', 217)
        self.kMax = kwargs.get('kMax', 71)
        self.jLE = kwargs.get('jLE', 109)
        self.jTE = kwargs.get('jTE', 73)
        self.xSF = kwargs.get('xSF', 1.12)
        self.ySF = kwargs.get('ySF', 1.12 if self.isWall == 0 else 1.0)
        self.dxdy = kwargs.get('dxdy', 1.0)
        self.gridType = kwargs.get('gridType', 'O')
        self.gridGenType = kwargs.get('gridGenType', 'Steger-Sorenson')

        ## RELAXATION PARAMETERS
        self.method = kwargs.get('method', 'PJ')
        self.wSLOR = kwargs.get('wSLOR', 1.89)
        self.iterMax = kwargs.get('iterMax', 1000)
        self.convCriteria = kwargs.get('convCriteria', 0.001)

    @classmethod
    def create(cls, **kwargs):
        return cls(**kwargs)