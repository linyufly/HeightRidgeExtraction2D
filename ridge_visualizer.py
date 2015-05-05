import shapefile

kRidgeFile = 'ridges.txt'
kRidgeShapeFilePrefix = 'ridges'
kLengthThreshold = 0

def main():
  reader = open(kRidgeFile, 'r')

  num_traj = int(reader.readline())

  writer = shapefile.Writer(shapefile.POLYLINE)
  writer.field('ID', 'C', '20')

  for count in range(num_traj):
    print 'count: {0}', count

    num_points = int(reader.readline())
    curr_traj = []

    for p in range(num_points):
      curr_traj.append([float(num) for num in reader.readline().split()])

    print curr_traj
    print len(curr_traj)

    if len(curr_traj) >= kLengthThreshold:
      writer.line(parts = [curr_traj])
      writer.record('line {0}', count)

  writer.save(kRidgeShapeFilePrefix)

  reader.close()

if __name__ == '__main__':
  main()

