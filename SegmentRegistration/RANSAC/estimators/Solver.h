#ifndef GROUPSAC_ESTIMATORS_SOLVER_H
#define GROUPSAC_ESTIMATORS_SOLVER_H

#include <vector>
using namespace std;

namespace groupsac  {
namespace estimators  {

// General Solver class.
// -- Compute a list of model that fit the candidates data.
template<typename T, typename Model>
class Solver
{
  /// Define the minimum number of candidates to compute a model
  enum { MINIMUM_SAMPLES = 0 };
public :
  typedef Model ModelType;
  /**
    * Find the solutions that fit the candidates data.
    *
    * \param[in] candidates (The input data).
    * \param[in] model (The model(s) that fit the data).
    *
    * \return True if success, else false.
    */
  virtual bool solve(const T & candidates, vector<Model> & model)const =0;

  /**
  * How many points are required at least to compute a model
  * \return the internal MINIMUM_SAMPLES value.
  */
  virtual int get_MINIMUM_SAMPLES() const =0;
};

}; // namespace estimators
}; // namespace groupsac

#endif //GROUPSAC_ESTIMATORS_SOLVER_H
