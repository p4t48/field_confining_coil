from sqSolenoid import *
import numpy as np

# Set parameters for which to calculate the maps.

current = 0.1
coilWidth = 0.2


""" Calculate for one geometry only.  """

dataFile = "coil_geometries/square_coil_iso_100_outer_100/square_coil_iso_100_outer_100_inner_0.001.csv"


s = Solenoid(dataFile)

s.SquareBoundaryPoints()
s.SquareFrontPanelMisalignment(0.007, coilWidth)
s.PCBThickness(0.0009, coilWidth)

s.Plot2D()
#s.Build3DGeom(coilWidth)
#s.Plot3D()
#s.Arrows3D()
#s.Draw3D()


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
