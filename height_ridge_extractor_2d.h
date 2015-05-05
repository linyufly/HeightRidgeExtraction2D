// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef HEIGHT_RIDGE_EXTRACTOR_2D_H_
#define HEIGHT_RIDGE_EXTRACTOR_2D_H_

#include <algorithm>
#include <vector>

class HeightRidgeExtractor2D {
 public:
  static void extract_ridges(
      int *dimensions, double *spacing, double *origin,
      double **scalar_field,
      std::vector<std::pair<double, double> > *ridges);
};

#endif  // HEIGHT_RIDGE_EXTRACTOR_2D_H_

