// Author: Mingcheng Chen (linyufly@gmail.com)

#include "height_ridge_extractor_2d.h"

#include "math.h"
#include "util.h"

#include <vtkMath.h>

namespace {

double partial_x(double **scalar_field, int nx, int ny, int x, int y, double dx, double dy) {
  int x_1 = (x > 0) ? x - 1 : x;
  int x_2 = (x < nx - 1) ? x + 1 : x;
  return (scalar_field[x_2][y] - scalar_field[x_1][y]) / (dx * (x_2 - x_1));
}

double partial_y(double **scalar_field, int nx, int ny, int x, int ny, double dx, double dy) {
  int y_1 = (y > 0) ? y - 1 : y;
  int y_2 = (y < ny - 1) ? y + 1 : y;
  return (scalar_field[x][y_2] - scalar_field[x][y_1]) / (dy * (y_2 - y_1));
}

double partial_xx(double **scalar_field, int nx, int ny, int x, int y, double dx, double dy) {
  int x_1 = (x > 0) ? x - 1 : x;
  int x_2 = (x < nx - 1) ? x + 1 : x;
  return (partial_x(scalar_field, nx, ny, x_2, y, dx, dy)
         - partial_x(scalar_field, nx, ny, x_1, y, dx, dy)) / (dx * (x_2 - x_1));
}

double partial_yy(double **scalar_field, int nx, int ny, int x, int y, double dx, double dy) {
  int y_1 = (y > 0) ? y - 1 : y;
  int y_2 = (y < ny - 1) ? y + 1 : y;
  return (partial_y(scalar_field, nx, ny, x, y_2, dx, dy)
         - partial_y(scalar_field, nx, ny, x, y_1, dx, dy)) / (dy * (y_2 - y_1));
}

double partial_xy(double **scalar_field, int nx, int ny, int x, int y, double dx, double dy) {
  int x_1 = (x > 0) ? x - 1 : x;
  int x_2 = (x < nx - 1) ? x + 1 : x;
  return (partial_y(scalar_field, nx, ny, x_2, y, dx, dy)
         - partial_y(scalar_field, nx, ny, x_1, y, dx, dy)) / (dx * (x_2 - x_1));
}

}

void HeightRidgeExtractor2D::extract_ridges(
    int *dimensions, double *spacing, double *origin,
    double **scalar_field,
    std::vector<std::pair<double, double> > *ridges) {
  int nx = dimensions[0];
  int ny = dimensions[1];

  double dx = spacing[0];
  double dy = spacing[1];

  double ***grad = create_3d_array<double>(nx, ny, 2);
  double **eig_val = create_matrix<double>(nx, ny);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      grad[x][y][0] = partial_x(scalar_field, nx, ny, x, y, dx, dy);
      grad[x][y][1] = partial_y(scalar_field, nx, ny, x, y, dx, dy);
    }
  }

  double ***e2 = create_3d_array<double>(nx, ny, 2);
  double **hessian = create_matrix<double>(2, 2);
  double **v = create_matrix<double>(2, 2);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      hessian[0][0] = partial_xx(scalar_field, nx, ny, x, y, dx, dy);
      hessian[0][1] = hessian[1][0]
                    = partial_xy(scalar_field, nx, ny, x, y, dx, dy);
      hessian[1][1] = partial_yy(scalar_field, nx, ny, x, y, dx, dy);

      double w[2];
      vtkMath::JacobiN(hessian, 2, w, v);

      eig_val[x][y] = w[1];
      e2[x][y][0] = v[0][1];
      e2[x][y][1] = v[1][1];
    }
  }

  for (int x = 0; x + 1 < nx; x++) {
    for (int y = 0; y + 1 < ny; y++) {
      double directions[4][2];
      int count = 0;
      for (int dx = 0; dx < 2; dx++) {
        for (int dy = 0; dy < 2; dy++) {
          directions[count][0] = e2[x + dx][y + dy][0];
          directions[count][1] = e2[x + dx][y + dy][1];
        }
      }


    }
  }

  delete_matrix(eig_val);
  delete_matrix(v);
  delete_matrix(hessian);
  delete_3d_array(e2);
  delete_3d_array(grad);
}
