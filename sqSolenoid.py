import re
from os.path import isfile
import glob
import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from itertools import chain
from mpl_toolkits.mplot3d import Axes3D
from point_extraction import *


class Solenoid(Geometry):

    def __init__(self, dataFile):
        Geometry.__init__(self, dataFile)
        self.fieldVals = []
        self.gradVals = []
        self.fieldGrid = np.array([])
        self.gradientGrid = np.array([])

    def Bfield(self, current, loc, pk, pk1):
        """ Calculates the B-field of a straight current carrying wire with the formula (9) of J. D. Hanson and S. P. Hirshman """

        m04Pi = 10**(-7)
        Ri = loc - pk
        RiNorm = np.linalg.norm(Ri)
        Rf = loc - pk1
        RfNorm = np.linalg.norm(Rf)
        
        geom = np.cross(Ri,Rf) * ( (RiNorm + RfNorm) / (RiNorm*RfNorm*(RiNorm*RfNorm + np.dot(Ri, Rf))) )
    
        return m04Pi * current * geom

    def loopField(self, coilCurrent, loc, loop):
        """ Calculates the B field components at a given point loc for one loop of the coil."""

        singleLoopField = np.array([0.0,0.0,0.0])
        wireFields = np.array([self.Bfield(coilCurrent, loc, loop[i], loop[i+1]) for i in range(len(loop)-1)])

        # Sum up all the contributions from the wires composing the loop
        for i in wireFields:
            singleLoopField += i

        return singleLoopField

    def CalcB(self, coilCurrent, xPos, yIni, yFin, zPos, nBy, grid):
        """ Calculates the B field components along a given axis for a coil geometry defined by the points  """

        # Fields produced by all the loops at the evaluation points
        evalFields = np.zeros(shape=(nBy,6))

        # Spacing between the evaluation points
        evalSpacing = (yFin - yIni) / float(nBy)

        # Calculates fields at all the specified evaluation points
        for k in range(nBy):
            # Initialize arrays to store computed data
            loc = np.array([xPos, yIni + k * evalSpacing, zPos])
            singleEvalField = np.array([0.0,0.0,0.0])

            wireField = np.array([self.loopField(coilCurrent, loc, loop) for loop in self.points])

            for i in wireField:
                singleEvalField += i

            # Save (x,y,z,Bx,By,Bz) for later field maps
            evalFields[k] = np.append((xPos, yIni + k * evalSpacing, zPos), singleEvalField)

        # Data for simple line plot or contour plot
        if grid == False:
            self.fieldVals = evalFields
        elif grid == True:
            if self.fieldGrid.size:
                self.fieldGrid = np.vstack((self.fieldGrid, evalFields))
            else:
                self.fieldGrid = evalFields

    def CalcdB(self, coilCurrent, xPos, yIni, yFin, zPos, nBy, grid):
        """ Calculates the gradients dBx/dz, dBy/dz, dBz/dz along the solenoid axis """

        # Fields produced by all the loops at the evaluation points
        evalGradients = np.zeros(shape=(nBy,6))

        # Spacing between the evaluation points
        evalSpacing = (yFin - yIni) / float(nBy)

        h = 0.00001
        for k in range(nBy):
            # Location at which the field should be evaluated min, max <-> -+ h
            locMin = np.array([xPos, yIni + k * evalSpacing - h, zPos])
            locMax = np.array([xPos, yIni + k * evalSpacing + h, zPos])

            # Stores the field at one give evaluation point (-+ h)
            evalPointFieldMin = np.array([0.0,0.0,0.0])
            evalPointFieldMax = np.array([0.0,0.0,0.0])

            wireFieldMin = np.array([self.loopField(coilCurrent, locMin, loop) for loop in self.points])
            wireFieldMax = np.array([self.loopField(coilCurrent, locMax, loop) for loop in self.points])

            for i in wireFieldMin:
                evalPointFieldMin += i
            for i in wireFieldMax:
                evalPointFieldMax += i

            # Calculate the gradient with symmetric derivative (f(x + h) - f(x - h))/2h
            gradient = (evalPointFieldMax - evalPointFieldMin)/(2*h)

            # Save (x,y,z,dBx/dz,dBy/dz,dBz/dz) for later field maps
            evalGradients[k] = np.append((xPos, yIni + k*evalSpacing, zPos), gradient)
 
          
        # Data for simple line or contour plot
        if grid == False:
            self.gradVals = evalGradients
        elif grid == True:
            if self.gradientGrid.size:
                self.gradientGrid = np.vstack((self.gradientGrid, evalGradients))
            else:
                self.gradientGrid = evalGradients


    def PlotB(self, coilCurrent, xPos, yIni, yFin, zPos, nBy):
        """ Plots B components as a function of z through (x,y) in the solenoid """

        self.CalcB(coilCurrent, xPos, yIni, yFin, zPos, nBy, False)

        # Get z field components in uT without the coordinates (x,y,z,Bx,By,Bz)
        fieldComponents = np.column_stack(self.fieldVals[:,-3:]*10**6) # In uT

        yVals = np.linspace(yIni, yFin, nBy)

        # Plot Bx, By, Bz along the centered z axis
        plt.subplot(3,1,1)
        plt.xlim([yIni,yFin])
        plt.title('Magnetic field components', size=28)
        plt.plot(yVals, fieldComponents[0]*1000) # In nT
        plt.ylabel('Bx (nT)', size=20)

        plt.subplot(3,1,2)
        plt.xlim([yIni,yFin])
        plt.plot(yVals, fieldComponents[1])
        plt.ylabel('By (nT)', size=20)

        plt.subplot(3,1,3)
        plt.xlim([yIni,yFin])
        plt.plot(yVals, fieldComponents[2]*1000) # In nT
        plt.ylabel('Bz (uT)', size=20)
        plt.xlabel('Position, y (m)', size=20)

        plt.subplots_adjust(hspace=0.4)
        plt.show()
        plt.clf()

    def PlotdB(self, coilCurrent, xPos, yIni, yFin, zPos, nBy):
        """ Plots dB components as a function of z through (x,y) in the solenoid """

        self.CalcdB(coilCurrent, xPos, yIni, yFin, zPos, nBy, False)

        # Get z field components in uT without the coordinates (x,y,z,Bx,By,Bz)
        gradientComponents = np.column_stack(self.gradVals[:,-3:]*10**4) # In uT/cm

        yVals = np.linspace(yIni, yFin, nBy)

        # Plot Bx, By, Bz along the centered z axis
        plt.subplot(3,1,1)
        plt.xlim([yIni,yFin])
        plt.title('Magnetic gradient components', size=28)
        plt.plot(yVals, gradientComponents[0])
        plt.ylabel('dBx/dz (uT/cm)', size=20)

        plt.subplot(3,1,2)
        plt.xlim([yIni,yFin])
        plt.plot(yVals, gradientComponents[1])
        plt.ylabel('dBy/dz (uT/cm)', size=20)

        plt.subplot(3,1,3)
        plt.xlim([yIni,yFin])
        plt.plot(yVals, gradientComponents[2])
        plt.ylabel('dBz/dz (uT/cm)', size=20)
        plt.xlabel('Position, y (m)', size=20)

        plt.subplots_adjust(hspace=0.4)
        plt.show()
        plt.clf()

    def CalcBGrid(self, coilCurrent, xVals, yIni, yFin, zVals, nBy):
        """ Calculates a grid of field values for a countour plot """

        # Make sure the grid is empty before recalculating
        self.fieldGrid = np.array([])

        for i in xVals:
            for k in zVals:
                self.CalcB(coilCurrent, i, yIni, yFin, k, nBy, True)

    def CalcdBGrid(self, coilCurrent, xVals, yIni, yFin, zVals, nBy):
        """ Calculates a grid of gradient values for a countour plot """

        # Make sure the grid is empty before recalculating
        self.gradientGrid = np.array([])

        for i in xVals:
            for k in zVals:
                self.CalcdB(coilCurrent, i, yIni, yFin, k, nBy, True)


    def ContourB(self, coilCurrent, xVals, yIni, yFin, zVals, nBy):
        """ Makes the contour plot of the calculated field values over the given region. """

        self.CalcBGrid(coilCurrent, xVals, yIni, yFin, zVals, nBy)

        minX = min(xVals)
        maxX = max(xVals)
        x = self.fieldGrid[:,0]
        y = self.fieldGrid[:,1]
        z = self.fieldGrid[:,2]
        Bx = self.fieldGrid[:,3]*10**6 # In uT
        By = self.fieldGrid[:,4]*10**6 # In uT
        Bz = self.fieldGrid[:,5]*10**6 # In uT

        # Define the grid
        yi = np.linspace(min(y),max(y),400)
        xi = np.linspace(min(x),max(x),400)
        Bxi = griddata(x, y, Bx, xi, yi, interp='linear')
        Byi = griddata(x, y, By, xi, yi, interp='linear')
        Bzi = griddata(x, y, Bz, xi, yi, interp='linear')

        plt.subplot(3,1,1)
        plt.contourf(xi, yi, Bxi, 50, cmap='cubehelix')
        plt.title('Flux return field maps', size=20)
        plt.ylabel('y (m)', size=20)
        cbarx = plt.colorbar()
        cbarx.ax.set_ylabel('Bx ($\mu$T)')
        plt.scatter(x, y, marker='o', s=0.1)
        plt.ylim(yIni,yFin)
        plt.xlim(minX,maxX)

        plt.subplot(3,1,2)
        plt.contourf(xi, yi, Byi, 50, cmap='cubehelix')
        plt.ylabel('y (m)', size=20)
        cbary = plt.colorbar()
        cbary.ax.set_ylabel('By ($\mu$T)')
        plt.scatter(x, y, marker='o', s=0.1)
        plt.ylim(yIni,yFin)
        plt.xlim(minX,maxX)

        plt.subplot(3,1,3)
        plt.contourf(xi, yi, Bzi, 50, cmap='cubehelix')
        plt.ylabel('y (m)', size=20)
        plt.xlabel('Transverse direction, x (m)', size=20)
        cbarz = plt.colorbar()
        cbarz.ax.set_ylabel('Bz ($\mu$T)')
        plt.scatter(x, y, marker='o', s=0.1)
        plt.ylim(yIni,yFin)
        plt.xlim(minX,maxX)

        plt.savefig("%s_B.pdf" % self.fileName)
        plt.clf()

    def ContourBi(self, coilCurrent, xVals, yIni, yFin, zVals, nBy):
        """ Makes the contour plot of the calculated Bi (i.e. x, y) field values over the given region. """

        self.CalcBGrid(coilCurrent, xVals, yIni, yFin, zVals, nBy)

        minX = min(xVals)
        maxX = max(xVals)
        x = self.fieldGrid[:,0]
        y = self.fieldGrid[:,1]
        z = self.fieldGrid[:,2]
        Bx = self.fieldGrid[:,3]*10**6 # In uT
        By = self.fieldGrid[:,4]*10**6 # In uT

        # Define the grid
        yi = np.linspace(min(y),max(y),400)
        xi = np.linspace(min(x),max(x),400)
        Bxi = griddata(x, y, Bx, xi, yi, interp='linear')
        Byi = griddata(x, y, By, xi, yi, interp='linear')

        #levels = np.arange(-35, 35, 0.5)

        plt.contourf(xi, yi, Bxi, 50, cmap='cubehelix', extend='both')
        plt.title('Prototype coil Bx map', size=20)
        plt.ylabel('y (m)', size=20)
        plt.xlabel('x (m)', size=20)
        cbarx = plt.colorbar()
        cbarx.ax.set_ylabel('Bx ($\mu$T)')
        plt.scatter(x, y, marker='o', s=0.1)
        plt.ylim(yIni,yFin)
        plt.xlim(minX,maxX)

        plt.savefig("%s_Bx.pdf" % self.fileName)
        plt.clf()

        plt.contourf(xi, yi, Byi, 50, cmap='cubehelix', extend='both')
        plt.title('Prototype coil By map', size=20)
        plt.ylabel('y (m)', size=20)
        plt.xlabel('x (m)', size=20)
        cbarx = plt.colorbar()
        cbarx.ax.set_ylabel('By ($\mu$T)')
        plt.scatter(x, y, marker='o', s=0.1)
        plt.ylim(yIni,yFin)
        plt.xlim(minX,maxX)

        plt.savefig("%s_By.pdf" % self.fileName)
        plt.clf()


    def ContourdB(self, coilCurrent, xVals, yIni, yFin, zVals, nBy):
        """ Makes the contour plot of the calculated field values over the given region. """

        self.CalcdBGrid(coilCurrent, xVals, yIni, yFin, zVals, nBy)

        minX = min(xVals)
        maxX = max(xVals)
        x = self.gradientGrid[:,0]
        y = self.gradientGrid[:,1]
        z = self.gradientGrid[:,2]
        dBx = self.gradientGrid[:,3]*10**4 # In uT/cm
        dBy = self.gradientGrid[:,4]*10**4 # In uT/cm
        dBz = self.gradientGrid[:,5]*10**4 # In uT/cm

        # Define the grid
        yi = np.linspace(min(y),max(y),400)
        xi = np.linspace(min(x),max(x),400)
        dBxi = griddata(x, y, dBx, xi, yi, interp='linear')
        dByi = griddata(x, y, dBy, xi, yi, interp='linear')
        dBzi = griddata(x, y, dBz, xi, yi, interp='linear')

        plt.subplot(3,1,1)
        plt.contourf(xi, yi, dBxi, 50, cmap='cubehelix')
        plt.title('Flux return gradient maps', size=20)
        plt.ylabel('y (m)', size=20)
        cbarx = plt.colorbar()
        cbarx.ax.set_ylabel('dBx/dy ($\mu$T/cm)')
        plt.scatter(x, y, marker='o', s=0.1)
        plt.ylim(yIni,yFin)
        plt.xlim(minX,maxX)

        plt.subplot(3,1,2)
        plt.contourf(xi, yi, dByi, 50, cmap='cubehelix')
        plt.ylabel('y (m)', size=20)
        cbary = plt.colorbar()
        cbary.ax.set_ylabel('dBy/dy ($\mu$T/cm)')
        plt.scatter(x, y, marker='o', s=0.1)
        plt.ylim(yIni,yFin)
        plt.xlim(minX,maxX)

        plt.subplot(3,1,3)
        plt.contourf(xi, yi, dBzi, 50, cmap='cubehelix')
        plt.ylabel('y (m)', size=20)
        plt.xlabel('Transverse direction, x (m)', size=20)
        cbarz = plt.colorbar()
        cbarz.ax.set_ylabel('dBz/dy ($\mu$T/cm)')
        plt.scatter(x, y, marker='o', s=0.1)
        plt.ylim(yIni,yFin)
        plt.xlim(minX,maxX)

        plt.savefig("%s_dB.pdf" % self.fileName)
        plt.clf()

    def MeandB(self, radius, solCenter, dB, current, shift):
        """ Calculates mean, median, std of a given region of interest (position of magnetometer) in the coil.  """
     
        outputFile = self.workingFolder + "%1.4f_A_connected_b_shifted.txt" % current

        # Create header of output file
        if not isfile(outputFile):
            f = open(outputFile,'w')
            f.write("#Inner(m)\tMisalignment(m)\tMean_dBy/dy(nT/cm)\tMedian_dBy/dy(nT/cm)\tStd(nT/cm)\tNumber_of_values\n")
            f.close()

        # Find the gradient values in the square region of interest at the center of the solenoid
        # Have to use logical_and to compare arrays of booleans
        xConditionLow = solCenter - radius < self.gradientGrid[:,0]
        xConditionHigh = self.gradientGrid[:,0] < solCenter + radius
        xCondition = np.logical_and(xConditionLow, xConditionHigh)

        yConditionLow = -radius < self.gradientGrid[:,1]
        yConditionHigh = self.gradientGrid[:,1] < radius
        yCondition = np.logical_and(yConditionLow,yConditionHigh)

        condition = np.logical_and(xCondition,yCondition)

        # Extract values of interest from the gradient values in the defined region
        meanGrad = np.mean(self.gradientGrid[condition,dB]*10**7) # In nT/cm
        medianGrad = np.median(self.gradientGrid[condition,dB]*10**7) # In nT/cm
        stdGrad = np.std(self.gradientGrid[condition,dB]*10**7) # In nT/cm
        lenGrad = len(self.gradientGrid[condition,dB])
        innerSq = float(re.search(r'_outer_D_inner_([-+]?\d*\.\d+|\d+)',self.dataFile).group(1))

        # Output to file in the same folder as data
        f = open(outputFile, 'a')
        f.write('%f\t%f\t%f\t%f\t%f\t%d\n'%(innerSq,shift,meanGrad,medianGrad,stdGrad,lenGrad))
        f.close()

    def ROI(self, zValue, connection, current):

        # Format data from region of interest to DataFrame (simpler export)
        gradientData = self.gradientGrid
        fieldData = self.fieldGrid

        gradientData[:,[3,4,5]] *= 10**7 # In nT/cm
        fieldData[:,[3,4,5]] *= 10**6 # In uT

        data = np.concatenate((gradientData, fieldData[:,[3,4,5]]), axis=1)
        df = pd.DataFrame.from_records(data, columns=["x","y","z","dBx/dy(nT/cm)","dBy/dy(nT/cm)","dBz/dy(nT/cm)","Bx(uT)","By(uT)","Bz(uT)"])

        # Pandas export to csv
        outputFile = self.workingFolder + "SqSolenoid_ROI_connected_%s_z_%1.3fm_current_%1.3fA.csv" % (connection, zValue, current)
        df.to_csv(outputFile, index=False, sep='\t')
        

# Set parameters for which to calculate the maps.

current = 0.1
coilWidth = 0.2

""" Calulate for simple square solenoid. 

dataFile = "01082016/square_coil_50_iso_25_outer/square_coil_50_iso_25_outer_D_inner_0.001.csv"
s = Solenoid(dataFile)
#s.SquareBoundaryPoints()
#s.Build3DGeom(coilWidth)

nWindings = 500
loopSpacing = 0.001
sideLength = 0.05
nPoints = 50

s.BuildSqSolenoid(nWindings, loopSpacing, sideLength)
s.Plot3D()
s.Arrows3D()
s.Draw3D()

ROILength = sideLength - 0.005
xVals = np.linspace(-ROILength,ROILength, num=nPoints)


for zValue in np.arange(-ROILength, ROILength+0.005, 0.005):
    print(zValue)
    s.CalcdBGrid(current, xVals, -ROILength, ROILength, [zValue], nPoints)
    s.CalcBGrid(current, xVals, -ROILength, ROILength, [zValue], nPoints)
    s.ROI(zValue, 'disc', current)    
"""

""" Calculate for one geometry only.  """

dataFile = "coil_geometries/square_coil_iso_100_outer_100/square_coil_iso_100_outer_100_inner_0.001.csv"


s = Solenoid(dataFile)

s.SquareBoundaryPoints()
#s.SquareFrontPanelMisalignment(0.007, coilWidth)
#s.PCBThickness(0.0009, coilWidth)

s.Plot2D()

s.Build3DGeom(coilWidth)
s.Plot3D()
s.Arrows3D()
s.Draw3D()


xMin = 0.05
xMax = 0.12
nPoints = 100
xVals = np.linspace(xMin,xMax, num=nPoints)

"""
for zValue in np.arange(0.01, coilWidth-0.005, 0.01):
    print(zValue)
    s.CalcdBGrid(current, xVals, -0.04, 0.04, [zValue], nPoints)
    s.CalcBGrid(current, xVals, -0.04, 0.04, [zValue], nPoints)
    s.ROI(zValue, 'bc', current)
"""

 
"""

    xMin = 0.13 + k
    xMax = 0.2 + k 
    nPoints = 100
    xVals = np.linspace(xMin,xMax, num=nPoints)

#    for zValue in (0.25/2, 0.25/2-0.05, 0.25/2+0.05, 0.01, 0.24):
    s.CalcdBGrid(current, xVals, -0.13, 0.13, [coilWidth/2], nPoints)
    s.CalcBGrid(current, xVals, -0.13, 0.13, [coilWidth/2], nPoints)
    s.ROI(coilWidth/2, 'bc', current, k)
"""
