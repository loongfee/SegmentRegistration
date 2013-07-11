#include <map>
#include <Eigen/Eigen>
#include <vector>
#include "lineConstant.h"
#include "func.h"
using namespace std;

class SpectralMatching
{
public:
	//each node in the graph is a possible line matching pair in the source and reference image
	struct Node{
		unsigned int srcID;//the index of line in the source image
		unsigned int refID;//the index of line in the reference image
	};

	// Specifies a vector of nodes.
	typedef std::vector<Node> Nodes_list;

	struct CompareL {
		bool operator() (const double& lhs, const double& rhs) const
		{return lhs>rhs;}
	};
	typedef  std::multimap<double,unsigned int,CompareL> EigenMAP;
	struct CompareS {
		bool operator() (const double& lhs, const double& rhs) const
		{return lhs<rhs;}
	};
	typedef  std::multimap<double,unsigned int,CompareS> DISMAP;
    SpectralMatching(){};
    void Matching(const std::vector<cvline_polar> &dataSource, const std::vector<cvline_polar> &dataReference,
    		std::vector<unsigned int> &matchResult);
	void setParameters(Eigen::VectorXd parameters);
	cvline_polar transform(const cvline_polar& l);
private:
    /* Build the symmetric non-negative adjacency matrix M, whose nodes are the potential assignments a = (i_l, j_r)
  * and whose weights on edges measure the agreements between pairs of potential assignments. That is where the pairwise
  * constraints are applied(c.f. A spectral technique for correspondence problems using pairwise constraints, M.Leordeanu).
  */
    void BuildAdjacencyMatrix_(const std::vector<cvline_polar> &dataSource, const std::vector<cvline_polar> &dataReference);
    /* Get the final matching from the principal eigenvector.
    */
    void MatchingResultFromPrincipalEigenvector_(const std::vector<cvline_polar> &linesInLeft,const std::vector<cvline_polar> &linesInRight,
    		std::vector<unsigned int > &matchResult);
	double dist(const cvline_polar& l1, const cvline_polar& l2);
	fPoint dist2(const cvline_polar& l1, const cvline_polar& l2);
    double globalRotationAngle_;//the approximate global rotation angle between image pairs
	Eigen::VectorXd m_parameters;
    /*construct a map to store the principal eigenvector and its index.
     *each pair in the map is in this form (eigenvalue, index);
     *Note that, we use eigenvalue as key in the map and index as their value.
     *This is because the map need be sorted by the eigenvalue rather than index
     *for our purpose.
      */
    EigenMAP eigenMap_;
    Nodes_list nodesList_;//save all the possible matched line pairs
    double minOfEigenVec_;//the acceptable minimal value in the principal eigen vector;
};