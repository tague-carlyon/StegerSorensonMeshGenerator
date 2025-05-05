import numpy as np
import matplotlib.pyplot as plt

class meshStruct:
    """
    A class representing a mesh structure.

    Attributes:
        params (dict): The dictionary containing configuration or initialization data.
        meshXandYs (numpy.ndarray): A NumPy array initialized to zeros with the same size as the mesh.
    """

    def __init__(self, params):
        """
        Initializes an instance of the meshStruct class.

        Args:
            params (dict): A dictionary containing configuration or initialization 
                           data for the mesh structure. This may include parameters 
                           such as mesh size, dimensions, or other relevant settings.
        """
        self.params = params
        self.jMax = params.jMax
        self.kMax = params.kMax
        self.jLE = params.jLE
        self.jTE = params.jTE
        #                       xi       , eta
        self.meshXs = np.zeros((self.jMax, self.kMax))
        self.meshYs = np.zeros((self.jMax, self.kMax))
        self.xSF = params.xSF
        self.ySF = params.ySF
        self.dxdy = params.dxdy
        self.ds = 0.004
        self.gridType = params.gridType
        self.etas = np.zeros((self.jMax, self.kMax))
        self.xis = np.zeros((self.jMax, self.kMax))
        for i in range(self.jMax):
            self.etas[i, :] = np.linspace(0, self.kMax-1, self.kMax)
        for i in range(self.kMax):
            self.xis[:, i] = np.linspace(0, self.jMax-1, self.jMax)
        
    def assignInternalConditions(self):
        """
        Assigns boundary conditions to the mesh structure.
        """
        if self.gridType == 'O':
            js = np.linspace(0, self.jMax, self.jLE)
            # create x array for airfoil points with cosine spacing
            xs = 0.5-0.5*np.cos((js-1) / (self.jMax-1)*np.pi)
            ys = np.zeros((2, self.jLE))
            # create airfoil
            if 'NACA' in self.params.foil and self.params.foil[4:6] == '00':
                # constant needed for NACA 4-digit airfoil generation
                xint = 1.0
                # get the NACA thickness
                th = float(self.params.foil[-2:])/100
                # y locations of for the top of the airfoil
                ys[0, :] = 5 * th * (0.2969 * np.sqrt(xs * xint) - \
                            0.1260 * xs * xint - 0.3516 * (xs * xint)**2 + \
                            0.2843 * (xs * xint)**3 - 0.1015 * (xs * xint)**4)
                # y locations of for the bottom of the airfoil
                ys[1, :] = -ys[0, :]
                
                # assign bottom of the airfoil to the mesh
                self.meshXs[:self.jLE-1, 0] = xs[1:][::-1]
                self.meshYs[:self.jLE-1, 0] = ys[1, 1:][::-1]
                # assign top of the airfoil to the mesh
                self.meshXs[self.jLE:, 0] = xs[1:]
                self.meshYs[self.jLE:, 0] = ys[0, 1:]
                
                # fix trailing edge of airfoil
                self.meshYs[0, 0] = 0.0
                self.meshYs[-1, 0] = 0.0
                self.meshYs[1, 0] = 0.5*(self.meshYs[1, 0] + 0.25 * (self.meshYs[2, 0] + self.meshYs[0, 0]))
                self.meshYs[-2, 0] = 0.5*(self.meshYs[self.jMax-2, 0] + 0.25 * (self.meshYs[self.jMax-1, 0] + self.meshYs[self.jMax-3, 0]))

            else:
                raise ValueError('Invalid airfoil type.')

    def linearStretchMeshToWall(self):
        """
        Linearly stretches the mesh structure to the wall.
        """
        if self.gridType == 'O':
            self.meshXs[0, 1] = self.meshXs[0, 0] + self.ds

            # set trailing edge peice
            for k in range(2, self.kMax):
                self.meshXs[0, k] = self.meshXs[0, k-1] + (self.meshXs[0, k-1] - self.meshXs[0, k-2]) * self.xSF
                self.meshYs[0, k] = 0.0
                self.meshXs[-1, k] = self.meshXs[0, k]
                self.meshYs[-1, k] = 0.0

            # assign wall boundary conditions
            XMAX = self.meshXs[0, -1]
            js = np.linspace(1, self.jMax, self.jMax, dtype=int)
            self.meshXs[:, -1] = -XMAX * np.cos((self.jLE-js) / (self.jLE-1)*np.pi)
            self.meshYs[:, -1] = -XMAX * np.sin((self.jLE-js) / (self.jLE-1)*np.pi)

            # linear stretch to wall using matrix operations
            percentages = (self.meshXs[0, 1:-1] - self.meshXs[0, 0]) / (self.meshXs[0, -1] - self.meshXs[0, 0])
            self.meshXs[:, 1:-1] = self.meshXs[:, [0]] + percentages * (self.meshXs[:, [-1]] - self.meshXs[:, [0]])
            self.meshYs[:, 1:-1] = self.meshYs[:, [0]] + percentages * (self.meshYs[:, [-1]] - self.meshYs[:, [0]])

    def relaxMesh(self):
        """
        Relaxes the mesh structure.
        """
        Res = self.params.convCriteria + 1.0

        dx = np.zeros((self.jMax-2, self.kMax-2))
        dy = np.zeros((self.jMax-2, self.kMax-2))

        currIter = 0

        while Res > self.params.convCriteria:

            currIter += 1
            if currIter % 100 == 0:
                self.plotMesh()

            resx, resy, alpha, beta, gamma = self.computeResidual()

            Res = np.linalg.norm(resx + resy, ord=2)
            print("Residual for meshGen: ", Res)

            match self.params.method:
                case 'PJ':
                    # calculate the relaxation factors
                    dx = -resx/(-2 * alpha - 2 * gamma)
                    dy = -resy/(-2 * alpha - 2 * gamma)
            
            # update the mesh
            self.meshXs[1:-1, 1:-1] += dx
            self.meshYs[1:-1, 1:-1] += dy


    def computeResidual(self):
        
        phi = np.zeros((self.jMax, self.kMax-2))
        psi = np.zeros((self.jMax, self.kMax-2))

        resx = np.zeros((self.jMax, self.kMax-2))
        resy = np.zeros((self.jMax, self.kMax-2))

        x_xi = np.zeros((self.jMax, self.kMax))
        y_xi = np.zeros((self.jMax, self.kMax))
        
        x_eta = np.zeros((self.jMax, self.kMax-2))
        y_eta = np.zeros((self.jMax, self.kMax-2))

        x_xixi = np.zeros((self.jMax, self.kMax))
        y_xixi = np.zeros((self.jMax, self.kMax))
        
        x_xieta = np.zeros((self.jMax, self.kMax-2))
        y_xieta = np.zeros((self.jMax, self.kMax-2))

        x_xi[1:-1, :] = (self.meshXs[2:, 1:-1] - self.meshXs[:-2, 1:-1]) / 2
        x_xi[0, :] = (self.meshXs[1, 1:-1] - self.meshXs[-1, 1:-1]) / 2
        x_xi[-1, :] = x_xi[0, :]

        
        x_xixi[1:-1, :] = (self.meshXs[:-2, :] - 2 * self.meshXs[1:-1, :] + self.meshXs[2:, :])
        x_xixi[0, :] = (self.meshXs[-1, :] - 2 * self.meshXs[0, :] + self.meshXs[1, :])
        x_xixi[-1, :] = x_xixi[0, :]

        x_xieta[:, :] = (x_xi[:, 2:] - x_xi[:, :-2]) / 2

        x_eta = (self.meshXs[:, 2:] - self.meshXs[:, :-2]) / 2
        x_ee = (self.meshXs[:, :-2] - 2 * self.meshXs[:, 1:-1] + self.meshXs[:, 2:])    

        y_xi[1:-1, :] = (self.meshYs[2:, 1:-1] - self.meshYs[:-2, 1:-1]) / 2
        y_xi[0, :] = (self.meshYs[1, 1:-1] - self.meshYs[-1, 1:-1]) / 2
        y_xi[-1, :] = y_xi[0, :]
        
        y_xixi[1:-1, :] = (self.meshYs[:-2, :] - 2 * self.meshYs[1:-1, :] + self.meshYs[2:, :])
        y_xixi[0, :] = (self.meshYs[-1, :] - 2 * self.meshYs[0, :] + self.meshYs[1, :])
        y_xixi[-1, :] = y_xixi[0, :]

        y_xieta[:, :] = (y_xi[:, 2:] - y_xi[:, :-2]) / 2

        y_eta = (self.meshYs[:, 2:] - self.meshYs[:, :-2]) / 2
        y_ee = (self.meshYs[:, :-2] - 2 * self.meshYs[:, 1:-1] + self.meshYs[:, 2:])    

        match self.params.gridGenType:
            case 'Elliptic':
                x_eta = (self.meshXs[:, 2:] - self.meshXs[:, :-2]) / 2
                y_eta = (self.meshYs[:, 2:] - self.meshYs[:, :-2]) / 2
            case 'Steger-Sorenson':
                #y_eta[:, 0] = np.sign(y_eta[:, 0]) * np.abs(self.ds * x_xi[:, 0] / np.sqrt(x_xi[:, 0]**2 + y_xi[:, 0]**2))
                #x_eta[:, 0] = np.sign(x_eta[:, 0]) * np.abs(self.ds * y_xi[:, 0] / np.sqrt(x_xi[:, 0]**2 + y_xi[:, 0]**2))

                #x_ee[:, 0] = 0.5 * (7 * self.meshXs[1:-1, 0] + 8 * self.meshXs[1:-1, 1] - self.meshXs[1:-1, 2]) - 3 * x_eta[:, 0]
                #y_ee[:, 0] = 0.5 * (7 * self.meshYs[1:-1, 0] + 8 * self.meshYs[1:-1, 1] - self.meshYs[1:-1, 2]) - 3 * y_eta[:, 0]
                True

        alpha = x_eta**2 + y_eta**2 

        beta = x_xi * x_eta + y_xi * y_eta
        
        gamma = x_xi**2 + y_xi**2

        match self.params.gridGenType:
            case 'Elliptic':
                resx = (alpha * x_xixi + beta * x_xieta + gamma * x_ee)
                resy = (alpha * y_xixi + beta * y_xieta + gamma * y_ee)

            case 'Steger-Sorenson':
                
                expa = np.exp(-self.etas[1:-1, 1:-1])
                expb = np.exp(-self.etas[1:-1, 1:-1])
                J = x_xi * y_eta - x_eta * y_xi

                Rx = -J[:, 0] ** 2 * (alpha[:, 0] * x_xixi[:, 0] - 2 * beta[:, 0] * x_xieta[:, 0] + gamma[:, 0] * x_ee[:, 0])
                Ry = -J[:, 0] ** 2 * (alpha[:, 0] * y_xixi[:, 0] - 2 * beta[:, 0] * y_xieta[:, 0] + gamma[:, 0] * y_ee[:, 0])
                # P0 and Q0 are defined at the wall of the airfoil
                P0 = J[:, 0] * (y_eta[:, 0] * Rx - x_eta[:, 0] * Ry)
                Q0 = J[:, 0] * (y_xi[:, 0] * Rx + x_xi[:, 0] * Ry)

                for i in range(0, self.jMax-2):
                    phi[i, :] = P0[i] * expa[i, :]
                    psi[i, :] = Q0[i] * expb[i, :]

                # calculate the residuals
                resx = (alpha * x_xixi - 2 * beta * x_xieta + gamma * x_ee) + J ** 2 * (phi * x_xi + psi * y_xi)
                resy = (alpha * y_xixi - 2 * beta * y_xieta + gamma * y_ee) + J ** 2 * (phi * y_xi + psi * x_xi)
                print("Shape of resx:", resx.shape)

        return resx, resy, alpha, beta, gamma

    def plotMesh(self):
        """
        Plots the mesh structure.
        """
        # plot the mesh
        fig, ax = plt.subplots()
        for x in range(self.jMax):
            ax.plot(self.meshXs[x, :], self.meshYs[x, :])
        for y in range(self.kMax):
            ax.plot(self.meshXs[:, y], self.meshYs[:, y])
        ax.set_aspect('equal', 'box')
        plt.show()