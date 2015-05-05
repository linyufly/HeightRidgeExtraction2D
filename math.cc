// Author: Mingcheng Chen (linyufly@gmail.com)

#include "math.h"

#include "util.h"

#include <vtkMath.h>

double **matrix_matrix_multiplication(double **a, double **b, int m, int n, int p) {
  double **c = create_matrix<double>(m, p);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < p; j++) {
      c[i][j] = 0.0;
      for (int k = 0; k < n; k++) {
        c[i][j] += a[i][k] * b[k][j];
      }
    }
  }
  return c;
}

double **transpose(double **a, int n, int m) {
  double **t = create_matrix<double>(m, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      t[j][i] = a[i][j];
    }
  }
  return t;
}

double *principal_component(double **vectors, int n, int k) {
  double **vec_trans = transpose(vectors, n, k);
  double **symm = matrix_matrix_multiplication(vec_trans, vectors, k, n, k);

  double *eigen_values = new double[k];
  double **eigen_vectors = create_matrix<double>(k, k);

  vtkMath::JacobiN(symm, k, eigen_values, eigen_vectors);

  double *result = new double[k];
  for (int i = 0; i < k; i++) {
    result[i] = eigen_vectors[i][0];
  }

  delete [] eigen_values;

  delete_matrix(eigen_vectors);
  delete_matrix(vec_trans);
  delete_matrix(symm);

  return result;
}
