import sqSolenoid as sq

# Set parameters for which to calculate the maps.

current = 0.1
coilWidth = 0.2


""" Calculate for one geometry only.  """

dataFile = "coil_geometries/square_coil_iso_100_outer_100/square_coil_iso_100_outer_100_inner_0.001.csv"


s = sq.Solenoid(dataFile)

s.SquareBoundaryPoints()
s.Plot2D()
