#include "../struct/support.hpp"

#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <vector>

/**
  Auxiliary function for removing (erasing) elements from an (unordered)
  map that satisfy a certain predicate. This functionality is *obsolete*
  in case the 'uniform container erasure' extensions are available. Note
  that the predicate will only be evaluated for the keys of the map.

  \param map Map
  \param pred Predicate
*/


int main(int, char**)
{
  Support S1;

  Support S2(0.001, 0.0004, 0.025);

  std::cout << S2 << std::endl;

  if (S2.min_pvalue() < 0.05)
  {
    std::cout << "min-pvalue below 0.05" << std::endl;
  }

}
