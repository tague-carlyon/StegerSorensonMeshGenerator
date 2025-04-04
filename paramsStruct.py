# Parameters structure
# Edit values here
class Params:
    def __init__(self):
        ## GEOMETRY AND BOUNDARY CONDITIONS
        # freestream Mach number
        self.Minf = 0.8
        # boolean whether upper boundary is wall or freestream
        # 1:  wall boundary
        # 0:  free stream
        self.isWall = 0
        # airfoil thickness
        self.th = 0.10
        # foil profile
        # options are 'NACA0010', 'BICONVEX'
        self.foil = 'NACA0010'

        ## FLUID PARAMETERS
        # gas constant
        # self.gamma = 1.4 for regular air
        # self.gamma = 1 for lsd
        self.gamma = 1.4
        
        ## MESH PARAMETERS
        # number of equally spaced grid points on all sides of airfoil
        self.kConst = 3
        # total grid points in x direction
        self.jMax = 217
        # total grid points in y direction
        self.kMax = 71
        # index at start of airfoil
        self.jLE = 109
        # index at end of airfoil
        self.jTE = 73
        # stretching factor in x direction
        self.xSF = 1.12
        # stretching factor in y direction (depends on type of upper BC)
        self.ySF = 1.12 if self.isWall==0 else 1.0
        # grid ratio at airfoil surface
        self.dxdy = 1.0
        self.gridType = 'O'
        # grid generation type 
        # options are 'TTM' (Thompson, Thames, Mastin), 'LSD' (linear stretching), 'UDM' (uniform)
        self.gridGenType = 'Elliptic'

        ## RELAXATION PARAMETERS
        # relaxation method
        # oprions are 'PJ', 'GS', 'SLOR', 'ADI'
        self.method = 'SLOR'
        # SLOR relaxation parameter
        self.wSLOR = 1.89
        # maximum number of iterations
        self.iterMax = 1000
        # ratio of current to initial residual L2 norm at convergence
        self.convCriteria = 0.001