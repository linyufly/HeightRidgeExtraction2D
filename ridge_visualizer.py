# import shapefile
import csv
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

kDensityMapFile = '/home/linyufly/GitHub/TrajectoryDensityClustering/density_map.csv'
kRidgeFile = 'ridges.txt'
# kRidgeShapeFilePrefix = 'ridges'
kLengthThreshold = 50

def main():
  with open(kDensityMapFile, 'r') as csv_file:
    table = [row for row in csv.reader(csv_file, delimiter = ',')]

  np_table = np.array([[float(element) for element in row] for row in table])

  fig, ax = plt.subplots()

  im = plt.imshow(np_table, cmap = matplotlib.cm.RdBu_r, vmin = np_table.min(), vmax = np_table.max(), extent = [-105.000000, 0.000000, 0.000000, 65.000000])
  im.set_interpolation('bilinear')
  # cb = fig.colorbar(im)

  # plt.show()

  reader = open(kRidgeFile, 'r')

  num_traj = int(reader.readline())

  # writer = shapefile.Writer(shapefile.POLYLINE)
  # writer.field('ID', 'C', '20')

  for count in range(num_traj):
    print 'count: {0}', count

    num_points = int(reader.readline())
    curr_traj = []

    for p in range(num_points):
      curr_traj.append([float(num) for num in reader.readline().split()])

    print curr_traj
    print len(curr_traj)

    if len(curr_traj) >= kLengthThreshold:
      # writer.line(parts = [curr_traj])
      # writer.record('line {0}', count)

      x = np.array([point[0] for point in curr_traj])
      y = np.array([point[1] for point in curr_traj])

      plt.plot(x, y)

  # writer.save(kRidgeShapeFilePrefix)

  reader.close()

  plt.show()

if __name__ == '__main__':
  main()

