// Author: Mingcheng Chen (linyufly@gmail.com)

#include "height_ridge_extractor_2d.h"

#include "util.h"

#include <cstdio>

#include <algorithm>
#include <vector>

namespace {

const char *kDensityMapFile = "/home/linyufly/GitHub/TrajectoryDensityClustering/density_map.txt";
const char *kRidgeFile = "ridges.txt";

void extract_ridges_test() {
  printf("extract_ridges_test {\n");

  FILE *fin = fopen(kDensityMapFile, "r");

  double x_min, x_max, y_min, y_max;
  int nx, ny;

  fscanf(fin, "%lf %lf %lf %lf", &x_min, &x_max, &y_min, &y_max);
  fscanf(fin, "%d %d", &nx, &ny);

  double **density_map = create_matrix<double>(nx, ny);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      fscanf(fin, "%lf", &density_map[x][y]);
    }
  }

  fclose(fin);

  std::vector<std::vector<std::pair<double, double> > > ridges;

  int dimensions[] = {nx, ny};
  double spacing[] = {(x_max - x_min) / (nx - 1), (y_max - y_min) / (ny - 1)};
  double origin[] = {x_min, y_min};
  HeightRidgeExtractor2D::extract_ridges(dimensions, spacing, origin, density_map, &ridges);

  FILE *fout = fopen(kRidgeFile, "w");

  fprintf(fout, "%d\n", static_cast<int>(ridges.size()));
  for (int r = 0; r < ridges.size(); r++) {
    fprintf(fout, "%d\n", static_cast<int>(ridges[r].size()));
    for (int p = 0; p < ridges[r].size(); p++) {
      fprintf(fout, "%lf %lf\n", ridges[r][p].first, ridges[r][p].second);
    }
  }

  fclose(fout);

  printf("} extract_ridges_test\n");
}

}

int main() {
  extract_ridges_test();

  return 0;
}

