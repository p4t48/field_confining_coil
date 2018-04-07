"""
This file contains the Geometry class. 


"""


import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import chain
from mpl_toolkits.mplot3d import Axes3D


def find_nearest(array,value):
    """ 
    Returns the location of the element in an array which is closest to a given value. 
    """
    
    idx = (np.abs(array-value)).argmin()
    return idx

def find_furthest(array,value):
    """ 
    Returns the location of the element in an array which is furthest to a given value. 
    """
    
    idx = (np.abs(array-value)).argmax()
    return idx


class Geometry:
    
    def __init__(self, dataFile):
        self.dataFile = dataFile
        self.fileName = re.search('(.+).csv', self.dataFile).group(1)
        self.workingFolder = re.search('(.+)/(.+)/', self.dataFile).group(0)
        self.points2D = []
        self.points = []


    def BoundaryPoints(self):
        """ 
        The method extracts points from a COMSOL generated file. 
        The points are the boundary points of the isopotentials drawn with COMSOL. 
        """

        data = pd.read_csv(self.dataFile, skiprows=8, header=None, names=["x","y","IsoLevel","Color","Radius"])
        
        # Sections of the coil such that I can sort them and get a continous array of consecutive loops.
        leftSol = []
        rightSol = []
        topRightFlux = []
        topLeftFlux = []
        bottomLeftFlux = []
        bottomRightFlux = []

        # Each isopotential value has four individual loops which need to be extracted
        for value in np.unique(data['IsoLevel'].values):

            # Start with one given isopotential value
            extracted = data[data['IsoLevel'] == value]

            # The color map set in COMSOL helps to sort the points
            posCol = extracted[extracted['Color'] > 0]
            negCol = extracted[extracted['Color'] < 0]

            # Data for the loops in the solenoid parts of the coil
            xPosPosCol = posCol[posCol['x'] > 0]
            xNegNegCol = negCol[negCol['x'] < 0]
            
            # Data for the loops in the flux return parts of the coil
            xNegPosCol = posCol[posCol['x'] < 0]            
            xPosNegCol = negCol[negCol['x'] > 0]

            # X closest, furthest to zero  criterion on left right parts is more stable
            for i in [xNegNegCol, xPosPosCol]:
                temp = find_nearest(i['x'],0)
                temp = i.loc[temp]
                minX = [temp[0], temp[1]]

                temp = find_furthest(i['x'],0)
                temp = i.loc[temp]
                maxX = [temp[0], temp[1]]


                if i[i['x'] < 0]['x'].any() :
                    leftSol.append([minX, maxX])
                else:
                    rightSol.append([minX, maxX])

            # Y closest, furthest to zero criterion on top bottom parts is more stable
            for i in [xNegPosCol, xPosNegCol]:
                temp = find_nearest(i['y'],0)
                temp = i.loc[temp]
                minY = [temp[0], temp[1]]

                temp = find_furthest(i['y'],0)
                temp = i.loc[temp]
                maxY = [temp[0], temp[1]]

                if i[i['x'] > 0]['x'].any() and i[i['y'] > 0]['y'].any():
                    topRightFlux.append([minY, maxY])
                elif i[i['x'] < 0]['x'].any() and i[i['y'] > 0]['y'].any():
                    topLeftFlux.append([minY, maxY])
                elif i[i['x'] < 0]['x'].any() and i[i['y'] < 0]['y'].any():
                    bottomLeftFlux.append([minY, maxY])
                else:
                    bottomRightFlux.append([minY, maxY])
       
        # Invert such that the coil gets built up in counter clockwise orientation
        bottomLeftFlux = bottomLeftFlux[::-1]
        bottomRightFlux = bottomRightFlux[::-1]
        topRightFlux = topRightFlux[::-1]
        topLeftFlux = topLeftFlux[::-1]
        self.points2D = list(chain.from_iterable([topLeftFlux, leftSol, bottomLeftFlux, bottomRightFlux, rightSol, topRightFlux]))

    def Build3DGeom(self, zHeight):
        """ 
        This method uses the previously extracted points (in 2D) and extends the geometry to 3D. 
        The points are grouped into groups of four which build up one single loop. 
        """

        points2D = list(chain.from_iterable(self.points2D))
        points2D = np.array(points2D)
        bottom = np.hstack((points2D, 0*np.ones((points2D.shape[0],1))))
        top = np.hstack((points2D, zHeight*np.ones((points2D.shape[0],1))))

        for i in range(int((len(bottom)-1)/2)):
            # The five points build a disconnected loop coil (for winding sense use Arrows3D)
#            loop = [bottom[2*i], bottom[2*i+1], top[2*i+1], top[2*i], bottom[2*i]] # Disconnected loop
#            loop = [top[2*i], bottom[2*i], bottom[2*i+1], top[2*i+1], top[2*i], top[2*i+2]] # connected_b
#            loop = [top[2*i+1], top[2*i], bottom[2*i], bottom[2*i+1], top[2*i+1], top[2*i+3]] # connected_a
            loop = [bottom[2*i], bottom[2*i+1], top[2*i+1], top[2*i], bottom[2*i+2]] # connected_bc
            self.points.append(loop)

        self.points = np.array(self.points)

    def Plot2D(self):
        """ 
        Visualize the extracted points to check if everything is ok. 
        """

        test = list(chain.from_iterable(self.points2D))
        test = np.array(test)
        test = np.column_stack(test)

        x = test[0]
        y = test[1]

        plt.plot(x,y,'bo')
        plt.show()
 
    def Plot3D(self):
        """ 
        Visualize the 3D geometry to check if everything is ok. 
        """

        test = list(chain.from_iterable(self.points))
        test = np.array(test)
        test = np.column_stack(test)

        x = test[0]
        y = test[1]
        z = test[2]

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.scatter(x,y,z)
        plt.show()


    def Arrows3D(self):
        """ 
        In order to check the winding direction of the current loop, plot the wires as arrows. 
        """

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        x, y, z = [], [], []
        u, v, w = [], [], []

        for loop in self.points:
            location = loop[1:5]
            direction = loop[1:5] - loop[0:4]
            location = np.column_stack(location)
            direction = np.column_stack(direction)

            x = location[0]
            y = location[1]
            z = location[2]

            u = direction[0]
            v = direction[1]
            w = direction[2]

            ax.quiver(x,y,z,u,v,w, length=0.05, arrow_length_ratio=0.1)
        plt.show()

    def Draw3D(self):
        """ 
        Draws a 3D squeleton of the coil as illustration of the coil geometry. 
        """
        

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        point  = np.array([0.02, 0.02, 0.175])
        normal = np.array([0, 0, 0.1])
        d = -point.dot(normal)
        # create x,y
        xx, yy = np.meshgrid(np.arange(-0.2,0.2,0.01), np.arange(-0.2,0.2,0.01))

        # calculate corresponding z
        z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]

        # plot the surface
        ax.plot_surface(xx, yy, z, alpha=0.5)
        ax.hold(True)
        
        x, y, z = [], [], []

        for loop in self.points:
            location = loop[0:5]
            location = np.column_stack(location)

            x = location[0]
            y = location[1]
            z = location[2]

            ax.plot(x, y, z, 'k')
        plt.axis('off')
        plt.show()

        
    def SquareBoundaryPoints(self):
        """ 
        The method extracts points from a COMSOL generated file for a square coil geometry. 
        The points are the boundary points of the isopotentials drawn with COMSOL. 
        """
        
        data = pd.read_csv(self.dataFile, skiprows=8, header=None, names=["x","y","IsoLevel","Color","Radius"])
     
        rightTrapezoid = []
        leftTrapezoid = []
        topTrapezoid = []
        bottomTrapezoid = []

        # Each isopotential value has four individual loops which need to be extracted
        for value in np.unique(data['IsoLevel'].values):

            # Start with one given isopotential value
            extracted = data[data['IsoLevel'] == value]

            # The color map set in COMSOL helps to sort the points
            posCol = extracted[extracted['Color'] > 0.5]
            negCol = extracted[extracted['Color'] < -0.5]

            zeroNegCol = extracted[extracted['Color'] < 0.5]
            zeroCol = zeroNegCol[zeroNegCol['Color'] > -0.5]

            yPosZeroCol = zeroCol[zeroCol['y'] > 0]
            yNegZeroCol = zeroCol[zeroCol['y'] < 0]

            # X closest, furthest to zero  criterion on left right parts is more stable
            for i in [posCol, negCol]:
                temp = find_nearest(i['x'],0)
                temp = i.loc[temp]
                minX = [temp[0], temp[1]]

                temp = find_furthest(i['x'],0)
                temp = i.loc[temp]
                maxX = [temp[0], temp[1]]


                if i[i['x'] < 0]['x'].any() :
                    leftTrapezoid.append([minX, maxX])
                else:
                    rightTrapezoid.append([minX, maxX])


            for i in [yPosZeroCol, yNegZeroCol]:
                temp = find_nearest(i['y'],0)
                temp = i.loc[temp]
                minY = [temp[0], temp[1]]

                temp = find_furthest(i['y'],0)
                temp = i.loc[temp]
                maxY = [temp[0], temp[1]]

 
                if i[i['y'] < 0]['y'].any() :
                    bottomTrapezoid.append([minY, maxY])
                else:
                    topTrapezoid.append([minY, maxY])

        leftTrapezoid = leftTrapezoid[::-1]
        rightTrapezoid = rightTrapezoid[::-1]
        self.points2D = list(chain.from_iterable([leftTrapezoid, bottomTrapezoid, rightTrapezoid, topTrapezoid]))


    def HexagonalBoundaryPoints(self):
        """ The method extracts points from a COMSOL generated file for a hexagonal coil geometry. The points are the boundary points of the isopotentials drawn with COMSOL. """
        
        data = pd.read_csv(self.dataFile, skiprows=8, header=None, names=["x","y","IsoLevel","Color","Radius"])
     
        topRightTrapezoid = []
        topLeftTrapezoid = []
        bottomRightTrapezoid = []
        bottomLeftTrapezoid = []
        topTrapezoid = []
        bottomTrapezoid = []

        # Each isopotential value has four individual loops which need to be extracted
        for value in np.unique(data['IsoLevel'].values):

            # Start with one given isopotential value
            extracted = data[data['IsoLevel'] == value]

            # The color map set in COMSOL helps to sort the points
            posCol = extracted[extracted['Color'] > 0.3]
            negCol = extracted[extracted['Color'] < -0.3]

            xNegPosCol = posCol[posCol['x'] < 0]
            xPosPosCol = posCol[posCol['x'] > 0]
            xPosNegCol = negCol[negCol['x'] > 0]
            xNegNegCol = negCol[negCol['x'] < 0]

            zeroNegCol = extracted[extracted['Color'] < 0.3]
            zeroCol = zeroNegCol[zeroNegCol['Color'] > -0.3]

            yPosZeroCol = zeroCol[zeroCol['y'] > 0]
            yNegZeroCol = zeroCol[zeroCol['y'] < 0]

            # X closest, furthest to zero  criterion on left right parts is more stable
            for i in [xNegPosCol, xPosPosCol, xPosNegCol, xNegNegCol]:
                temp = find_nearest(i['y'],0)
                temp = i.loc[temp]
                minX = [temp[0], temp[1]]

                temp = find_furthest(i['y'],0)
                temp = i.loc[temp]
                maxX = [temp[0], temp[1]]


                if i[i['x'] < 0]['x'].any() and i[i['Color'] < 0]['Color'].any():
                    bottomLeftTrapezoid.append([minX, maxX])
                if i[i['x'] < 0]['x'].any() and i[i['Color'] > 0]['Color'].any():
                    topLeftTrapezoid.append([minX, maxX])
                if i[i['x'] > 0]['x'].any() and i[i['Color'] < 0]['Color'].any():
                    topRightTrapezoid.append([minX, maxX])
                else:
                    bottomRightTrapezoid.append([minX, maxX])


            for i in [yPosZeroCol, yNegZeroCol]:
                temp = find_nearest(i['y'],0)
                temp = i.loc[temp]
                minY = [temp[0], temp[1]]

                temp = find_furthest(i['y'],0)
                temp = i.loc[temp]
                maxY = [temp[0], temp[1]]

 
                if i[i['y'] < 0]['y'].any() :
                    bottomTrapezoid.append([minY, maxY])
                else:
                    topTrapezoid.append([minY, maxY])

        topLeftTrapezoid = topLeftTrapezoid[::-1]
        bottomTrapezoid = bottomTrapezoid[::-1]
        topRightTrapezoid = topRightTrapezoid[::-1]
        self.points2D = list(chain.from_iterable([topTrapezoid, topLeftTrapezoid, bottomLeftTrapezoid, bottomTrapezoid, bottomRightTrapezoid, topRightTrapezoid]))

    def SquareSimpleMisalignment(self, shift, coilWidth):
        """ Misaligns the outer wire locations of the 2D points by "shift". """

        topTriangle, bottomTriangle, rightTriangle, leftTriangle = [], [], [], []
        crit = coilWidth/2 - 0.001

        for k in [topTriangle, bottomTriangle, rightTriangle, leftTriangle]:

            if k == topTriangle:
                truth =  [x[1][1] > crit for x in self.points2D]
            elif k == bottomTriangle:
                truth =  [x[1][1] < -crit for x in self.points2D]
            elif k == rightTriangle:
                truth =  [x[1][0] > crit for x in self.points2D]
            elif k == leftTriangle:
                truth =  [x[1][0] < -crit for x in self.points2D]

            for i,j in enumerate(truth):
                if j == True:
                    k.append(self.points2D[i])

            if k == topTriangle:
                for i,j in enumerate(k):
                    k[i][1][0] = j[1][0] + shift

            if k == bottomTriangle:
                for i,j in enumerate(k):
                    k[i][1][0] = j[1][0] - shift

            if k == rightTriangle:
                for i,j in enumerate(k):
                    k[i][1][1] = j[1][1] - shift

            if k == leftTriangle:
                for i,j in enumerate(k):
                    k[i][1][1] = j[1][1] + shift

        self.points2D = list(chain.from_iterable([leftTriangle, bottomTriangle, rightTriangle, topTriangle]))


    def SquareFrontPanelMisalignment(self, shift, coilWidth):
        """ Shifts the position of the front panel of the 2D points by "shift". """

        topTriangle, bottomTriangle, rightTriangle, leftTriangle = [], [], [], []
        crit = coilWidth/2 - 0.001

        for k in [topTriangle, bottomTriangle, rightTriangle, leftTriangle]:

            if k == topTriangle:
                truth =  [x[1][1] > crit for x in self.points2D]
            elif k == bottomTriangle:
                truth =  [x[1][1] < -crit for x in self.points2D]
            elif k == rightTriangle:
                truth =  [x[1][0] > crit for x in self.points2D]
            elif k == leftTriangle:
                truth =  [x[1][0] < -crit for x in self.points2D]

            for i,j in enumerate(truth):
                if j == True:
                    k.append(self.points2D[i])

            if k == topTriangle:
                for i,j in enumerate(k):
                    k[i][1][1] = j[1][1] + shift

            if k == bottomTriangle:
                for i,j in enumerate(k):
                    k[i][1][1] = j[1][1] - shift

            if k == rightTriangle:
                for i,j in enumerate(k):
                    k[i][1][0] = j[1][0] + shift

            if k == leftTriangle:
                for i,j in enumerate(k):
                    k[i][1][0] = j[1][0] - shift

        self.points2D = list(chain.from_iterable([leftTriangle, bottomTriangle, rightTriangle, topTriangle]))

    def PCBThickness(self, thickness, coilWidth):
        """ Misaligns the 2D points by "thickness/sqrt(2)" to simulate the PCB thickness. """

        topTriangle, bottomTriangle, rightTriangle, leftTriangle = [], [], [], []
        shift = thickness/np.sqrt(2)
        crit = coilWidth/2 - 0.001

        for k in [topTriangle, bottomTriangle, rightTriangle, leftTriangle]:

            if k == topTriangle:
                truth =  [x[1][1] > crit for x in self.points2D]
            elif k == bottomTriangle:
                truth =  [x[1][1] < -crit for x in self.points2D]
            elif k == rightTriangle:
                truth =  [x[1][0] > crit for x in self.points2D]
            elif k == leftTriangle:
                truth =  [x[1][0] < -crit for x in self.points2D]

            for i,j in enumerate(truth):
                if j == True:
                    k.append(self.points2D[i])

            if k == topTriangle:
                for i,j in enumerate(k):
                    k[i][0][1] = j[0][1] + shift

            if k == bottomTriangle:
                for i,j in enumerate(k):
                    k[i][0][1] = j[0][1] - shift

            if k == rightTriangle:
                for i,j in enumerate(k):
                    k[i][0][0] = j[0][0] + shift

            if k == leftTriangle:
                for i,j in enumerate(k):
                    k[i][0][0] = j[0][0] - shift

        self.points2D = list(chain.from_iterable([leftTriangle, bottomTriangle, rightTriangle, topTriangle]))

    def BuildSqSolenoid(self, nWindings, loopSpacing, sideLength):
        """ Defines points for a square cross section solenoid composed of n loops and 'spacing' between the separate loops """
    
        start = -nWindings*loopSpacing/2
        self.points = np.array([[sideLength, start, sideLength],[-sideLength, start, sideLength],[-sideLength, start, -sideLength],[sideLength, start, -sideLength],[sideLength, start, sideLength]])
        x1 = np.array([sideLength, start + loopSpacing, sideLength])
        x2 = np.array([-sideLength, start + loopSpacing, sideLength])
        x3 = np.array([-sideLength, start + loopSpacing, -sideLength])
        x4 = np.array([sideLength, start + loopSpacing, -sideLength])
        x5 = np.array([sideLength, start + loopSpacing, sideLength])

        for i in range(nWindings):
            self.points = np.vstack((self.points, x1, x2, x3, x4, x5))
            x1[1] += loopSpacing
            x2[1] += loopSpacing
            x3[1] += loopSpacing
            x4[1] += loopSpacing
            x5[1] += loopSpacing

        size = self.points.size
        self.points = self.points.reshape((size/15, 5, 3))
