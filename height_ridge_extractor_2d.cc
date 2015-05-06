// Author: Mingcheng Chen (linyufly@gmail.com)

#include "height_ridge_extractor_2d.h"

#include "math.h"
#include "util.h"

#include <algorithm>
#include <set>

#include <vtkMath.h>

namespace {

const double kEigenValueThreshold = -25000.0;

double linear_interpolation(double s1, double s2, double x1, double x2) {
  return (x2 - x1) * -s1 / (s2 - s1) + x1;
}

double dot_product(double x1, double y1, double x2, double y2) {
  return x1 * x2 + y1 * y2;
}

double dot_product(double *v1, double *v2) {
  return v1[0] * v2[0] + v1[1] * v2[1];
}

double partial_x(double **scalar_field, int nx, int ny, int x, int y, double dx, double dy) {
  int x_1 = (x > 0) ? x - 1 : x;
  int x_2 = (x < nx - 1) ? x + 1 : x;
  return (scalar_field[x_2][y] - scalar_field[x_1][y]) / (dx * (x_2 - x_1));
}

double partial_y(double **scalar_field, int nx, int ny, int x, int y, double dx, double dy) {
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
    std::vector<std::vector<std::pair<double, double> > > *ridges) {
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

      /// DEBUG ///
      // printf("%lf %lf\n%lf %lf\n--\n", hessian[0][0], hessian[0][1], hessian[1][0], hessian[1][1]);

      double w[2];
      vtkMath::JacobiN(hessian, 2, w, v);

      /// DEBUG ///
      // printf("w: %lf, %lf\n", w[0], w[1]);
      // printf("v[1]: %lf, %lf\n", v[0][1], v[1][1]);

      eig_val[x][y] = w[1];
      e2[x][y][0] = v[0][1];
      e2[x][y][1] = v[1][1];
    }
  }

  int ***edge_point = create_3d_array<int>(nx, ny, 2);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int d = 0; d < 2; d++) {
        edge_point[x][y][d] = -1;
      }
    }
  }

  std::vector<std::pair<double, double> > edge_points;
  std::set<std::pair<int, int> > links;

  for (int x = 0; x + 1 < nx; x++) {
    for (int y = 0; y + 1 < ny; y++) {
      double **directions = create_matrix<double>(4, 2);
      int count = 0;

      bool trimmed = false;
      for (int dx = 0; dx < 2; dx++) {
        for (int dy = 0; dy < 2; dy++) {
          directions[count][0] = e2[x + dx][y + dy][0];
          directions[count][1] = e2[x + dx][y + dy][1];
          count++;

          if (eig_val[x + dx][y + dy] > kEigenValueThreshold) {
            trimmed = true;
            break;
          }
        }

        if (trimmed) {
          break;
        }
      }

      if (trimmed) {
        continue;
      }

      double *pca = principal_component(directions, count, 2);

      double local_e2[2][2][2];
      for (int dx = 0; dx < 2; dx++) {
        for (int dy = 0; dy < 2; dy++) {
          local_e2[dx][dy][0] = e2[x + dx][y + dy][0];
          local_e2[dx][dy][1] = e2[x + dx][y + dy][1];

          if (dot_product(local_e2[dx][dy], pca) < 0.0) {
            local_e2[dx][dy][0] *= -1.0;
            local_e2[dx][dy][1] *= -1.0;
          }
        }
      }

      delete [] pca;
      delete_matrix(directions);

      int local_indices[4][2] = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};
      double scalars[4];

      for (int c = 0; c < 4; c++) {
        int dx = local_indices[c][0];
        int dy = local_indices[c][1];
        scalars[c] = dot_product(local_e2[dx][dy], grad[x + dx][y + dy]);
      }

      std::vector<int> intersections;

      for (int curr = 0; curr < 4; curr++) {
        int next = (curr + 1) % 4;
        if (scalars[curr] < 0.0 && scalars[next] >= 0.0
            || scalars[curr] >= 0.0 && scalars[next] < 0.0) {
          double x_1 = origin[0] + spacing[0] * (x + local_indices[curr][0]);
          double y_1 = origin[1] + spacing[1] * (y + local_indices[curr][1]);

          double x_2 = origin[0] + spacing[0] * (x + local_indices[next][0]);
          double y_2 = origin[1] + spacing[1] * (y + local_indices[next][1]);

          double x_itx = linear_interpolation(scalars[curr], scalars[next], x_1, x_2);
          double y_itx = linear_interpolation(scalars[curr], scalars[next], y_1, y_2);

          /// DEBUG ///
          // printf("%lf %lf, (x1, y1): (%lf, %lf), (x2, y2): (%lf, %lf), (x_itx, y_itx): (%lf, %lf)\n",
          //        scalars[curr], scalars[next], x_1, y_1, x_2, y_2, x_itx, y_itx);

          int dx_1 = local_indices[curr][0];
          int dy_1 = local_indices[curr][1];
          int dx_2 = local_indices[next][0];
          int dy_2 = local_indices[next][1];

          int located_x, located_y, located_d;

          if (dx_1 == dx_2) {
            located_x = x + dx_1;
            located_y = y;
            located_d = 1;
          } else {
            if (dy_1 != dy_2) {
              printf("Error: Invalid dy_1 and dy_2.\n");
              exit(0);
            }

            located_x = x;
            located_y = y + dy_1;
            located_d = 0;
          }

          if (edge_point[located_x][located_y][located_d] == -1) {
            edge_point[located_x][located_y][located_d] = edge_points.size();
            edge_points.push_back(std::pair<double, double>(x_itx, y_itx));
          }

          intersections.push_back(edge_point[located_x][located_y][located_d]);
        }
      }

      if (intersections.size() % 2 != 0) {
        printf("Error: Invalid number of intersections: %d.\n", static_cast<int>(intersections.size()));

        printf("%.20lf %.20lf %.20lf %.20lf\n", scalars[0], scalars[1], scalars[2], scalars[3]);
        exit(0);
      }

      for (int i = 0; i < intersections.size(); i += 2) {
        links.insert(std::pair<int, int>(intersections[i], intersections[i + 1]));
      }
    }
  }

  /// DEBUG ///
  printf("after construction.\n");

  std::vector<std::vector<int> > point_link(edge_points.size());
  for (std::set<std::pair<int, int> >::iterator itr = links.begin();
       itr != links.end(); ++itr) {
    point_link[itr->first].push_back(itr->second);
    point_link[itr->second].push_back(itr->first);
  }

  std::vector<bool> used(edge_points.size());
  std::fill(used.begin(), used.end(), false);

  /// DEBUG ///
  printf("before connecting.\n");

  for (int p = 0; p < edge_points.size(); p++) {
    if (used[p]) {
      continue;
    }

    if (point_link[p].size() == 1) {
      std::vector<std::pair<double, double> > path;
      for (int curr = p; curr != -1; ) {
        path.push_back(std::pair<double, double>(edge_points[curr].first, edge_points[curr].second));
        used[curr] = true;

        int new_curr = -1;
        for (int idx = 0; idx < point_link[curr].size(); idx++) {
          int next = point_link[curr][idx];
          if (!used[next]) {
            new_curr = next;
            break;
          }
        }

        curr = new_curr;
      }

      ridges->push_back(path);
    }
  }

  /// DEBUG ///
  printf("curr size = %d\n", static_cast<int>(ridges->size()));
  printf("check point 1\n");

  for (int p = 0; p < edge_points.size(); p++) {
    if (used[p]) {
      continue;
    }

    std::vector<std::pair<double, double> > path;
    for (int curr = p; curr != -1; ) {
      path.push_back(std::pair<double, double>(edge_points[curr].first, edge_points[curr].second));
      used[curr] = true;

      int new_curr = -1;
      for (int idx = 0; idx < point_link[curr].size(); idx++) {
        int next = point_link[curr][idx];
        if (!used[next]) {
          new_curr = next;
          break;
        }
      }

      curr = new_curr;
    }

    ridges->push_back(path);
  }

  /// DEBUG ///
  printf("curr size = %d\n", static_cast<int>(ridges->size()));
  printf("check point 2\n");

  int longest = 0;
  for (int r = 0; r < ridges->size(); r++) {
    longest = std::max(longest, static_cast<int>((*ridges)[r].size()));
  }

  printf("longest = %d\n", longest);

  delete_matrix(eig_val);
  delete_matrix(v);
  delete_matrix(hessian);
  delete_3d_array(edge_point);
  delete_3d_array(e2);
  delete_3d_array(grad);
}

