
#pragma once

#include <limits>
#include <iterator>
#include "non_randomness.h"

namespace groupsac  {
namespace ransac  {
namespace prosac {

using namespace std;

class PROSAC_handler
{
public:
  PROSAC_handler(int min_sample_num, int rounds_to_equal, const std::vector<int> & ordering)
  {
    // initialization
  	PR.exceeded = false;         
	  PR.n = min_sample_num-1; // the index of the last element in the current subset
	  PR.n_rounds = 0;         // the rounds used in the current subset
    PR.max_round = std::numeric_limits<int>::max();
		PR.T_n = 0;

    int T_N = rounds_to_equal;
    PR.T_n_ratio = T_N / (double)nchoosek(ordering.size(),min_sample_num);

    m_vec_ordering = ordering;
  }

  vector<int> fun_candidates(int nbCandidates, int round)
  {
    // if the budget in the current subset is used up
		if (PR.n_rounds >= PR.max_round)
    {
			PR.n_rounds = 0;
    }
		
		// start of a new subset
		if (PR.n_rounds == 0)
    {
			if (PR.n != m_vec_ordering.size())
      {
				PR.n = PR.n + 1;
        if ( PR.n > m_vec_ordering.size())
          PR.n = PR.n -1;
				int T_n_1 = PR.T_n;
				PR.T_n = PR.T_n_ratio * nchoosek(PR.n, nbCandidates);
				PR.max_round = ceil((double)(PR.T_n - T_n_1));
      }
			else
      {
				PR.max_round = std::numeric_limits<int>::max(); // always use the entire data set
				PR.exceeded = true;    // when exceeeded, does not have sample from the last element
      }
    }
		
		++PR.n_rounds; // One more iteration done.
		
    std::vector<int> vec_candidates;
    for (int i=0; i < PR.n; ++i) //		candidates = ordering(1:PR.n);
    {
      vec_candidates.push_back( m_vec_ordering[i] );
    }
    return vec_candidates;
  }

private:
  //-- Temporary function (Must be defined more generally).
  //---- (not yet done) Deterministic sampler
  //---- x Random sampler
  static vector<int> rand_sampler(const int nbPutatives, int MININUM_SAMPLES)
  {
    vector<int> array;
    int i = 0;
    std::set<int> sample_set;
    while (i < MININUM_SAMPLES)
    {
      int random_value_in_range = rand() % nbPutatives;
      if (sample_set.insert(random_value_in_range).second)
      {
        array.push_back(random_value_in_range);
        ++i;
      }
    }
    return array;
  }

  // pick up {sample_num} samples from {candidates}
  static vector<int> draw_samples(int candidatesSize, int sample_num)
  {
    vector<int> idx_rand = rand_sampler(candidatesSize,sample_num);
    return idx_rand;
  }

public:
  //Prosac sampler
  template <typename T>
  T sampler(const T & candidates, int MININUM_SAMPLES) const
  {
    vector<int> sampled;
    if (!PR.exceeded)
    {
			// draw samples from 1 to n-1
			sampled = draw_samples(candidates.size()-1, MININUM_SAMPLES-1);
			
			// the n-th point is mandatory
			assert(m_vec_ordering[PR.n-1] == candidates[candidates.size()-1]);
      sampled .push_back(candidates.size()-1); // add last point
    }
		else
    {
			sampled = draw_samples(candidates.size(), MININUM_SAMPLES);
    }
    return sampled;
  }

   bool fun_termination
    (
      vector<int> & best_inliers,
      const vector<int> & inliers,
      int round,
      int max_rounds,
      double logConfidence,
      int iSolverMinimalSamples,
      int iPutativesNumber
    )
  {
    const int min_sample_num = iSolverMinimalSamples;
    const double l1mp = logConfidence;
    // first update the minimal inlier we need for non-randomness if we start a new group configration
    if (PR.n_rounds == 0)
    {
      min_inlier = non_randomness(min_sample_num, PR.n);
      rounds_needed = ransac_rounds_needed(std::numeric_limits<int>::max(), min_sample_num, l1mp, PR.n, best_inliers.size());
      cout << "n: "               << PR.n << endl
           << "round: "           << PR.n_rounds << endl
           << "inliers :"         << best_inliers.size() << endl
           << "non_randomness: "  << min_inlier << endl
           << "rounds needed: "   << rounds_needed << endl
           << "allowed: "         << PR.max_round << endl;
    }

    if (inliers.size() > best_inliers.size())
    {
      best_inliers = inliers;
      rounds_needed = ransac_rounds_needed(std::numeric_limits<int>::max(), min_sample_num, l1mp, PR.n, inliers.size());
      cout << "n: "               << PR.n << endl
           << "round: "           << PR.n_rounds << endl
           << "inliers :"         << best_inliers.size() << endl
           << "non_randomness: "  << min_inlier << endl
           << "rounds needed: "   << rounds_needed << endl
           << "allowed: "         << PR.max_round << endl;
    }

    bool bTerminate = (inliers.size() > min_inlier && round >= rounds_needed) || round > max_rounds;

    if (bTerminate)
    {
      cout << "ready to terminate when non_randomness: " << endl;
      copy(inliers.begin(), inliers.end(), ostream_iterator<int>(cout, " "));
      cout << endl
           << "rounds needed: " << rounds_needed << endl
           << "min_inlier: " << min_inlier << endl;
    }
    return bTerminate;
  }

private:
 struct
  {
    bool exceeded;
    int n;
    int n_rounds;
    int max_round;
    int T_n;
    double T_n_ratio;
  } PR;

 std::vector<int> m_vec_ordering;

 int min_inlier;
 int rounds_needed;
};

}; // namespace prosac
}; // namespace ransac
}; // namespace groupsac
