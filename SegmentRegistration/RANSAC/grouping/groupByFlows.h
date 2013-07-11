

#include <iostream>
#include <map>
#include <vector>
#include <time.h>

#include "MeanShift/MeanShiftCluster.h"
#include "armadillo/include/armadillo"
using namespace arma;
// |
// | :: mat
// | :: vec

namespace groupsac  {
namespace grouping  {

mat weightedL2distance(const mat & x1s, const mat & x2)
{
  mat weights = "1 1 10 10";
  int x1_num = x1s.n_cols;
  return sum( weights*pow( (repmat(x2,1,x1_num) - x1s),2 ) );    //weighted dist squared
}

/// use Optical Flow based segmentation to group the data points for GroupSAC
bool groupByFlows(const mat & xs1, const mat & xs2, double bandwidth, int & seg_num, mat & vis_map, mat & clustCent, bool bVerbose = false)
{
  using namespace std;
  //% output:
  //% seg_num -   the total number of segments
  //% vis_map -   the row is the point index, the column is the segment index
  //% clustCent - the cluster centers

  clock_t start = clock();

  //%% generate the flow groups
  mat x_dists = xs1.row(0) - xs2.row(0);
  mat y_dists = xs1.row(1) - xs2.row(1);
  mat xs = xs1;
  xs = join_cols(xs, x_dists);
  xs = join_cols(xs, y_dists);

  map< int,vector<int> > cluster2dataCell;
  vector<int> data2cluster;
  MeanShiftCluster(xs, bandwidth, weightedL2distance, clustCent, data2cluster, cluster2dataCell);
  
  clock_t end = clock();
  if (bVerbose)
  {
      cout << "clustering takes "<<  ((double)(end-start))/CLOCKS_PER_SEC << " seconds" << endl;
  }
  
  seg_num = clustCent.n_cols;
  
  vis_map = zeros(cluster2dataCell.size(), x_dists.n_cols);
  // Build vis_map :
  for (int j=0; j < cluster2dataCell.size(); ++j)
  {
    for (int i = 0; i < cluster2dataCell[j].size(); ++i)
    {
      vis_map(j,cluster2dataCell[j][i]) = 1;
    }
  }
  
  return true; //TODO make a better return
}
}; // namespace grouping
}; // namespace groupsac
