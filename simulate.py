import sqSolenoid as sq

# Set parameters for which to calculate the maps.

current = 0.1
coilWidth = 0.2


""" Calculate for one geometry only.  """

dataFile = "coil_geometries/hexagonal_coil_iso_50_outer_100/hexagonal_coil_iso_50_outer_100_inner_0.01.csv"


s = sq.Solenoid(dataFile)

s.HexagonalBoundaryPoints()
s.Plot2D()
