#include "../struct/interval.hpp"

#include <algorithm>
#include <iostream>
#include <unordered_map>

/**
  Auxiliary function for removing (erasing) elements from an (unordered)
  map that satisfy a certain predicate. This functionality is *obsolete*
  in case the 'uniform container erasure' extensions are available. Note
  that the predicate will only be evaluated for the keys of the map.

  \param map Map
  \param pred Predicate
*/

template<typename Map, typename Pred> void erase_if(Map& map, Pred pred)
{
  for(auto it = map.begin(); it!= map.end(); )
  {
    if(pred(it->first))
      it = map.erase(it);
    else
      ++it;
  }
}

int main(int, char**)
{
  Interval I1;
  Interval I2(2, 3);
  Interval I3(2, 7);
  Interval I4(3, 5);

  std::unordered_map<Interval, bool> interval_map;
  interval_map[I1] = true;
  interval_map[I2] = false;
  interval_map[I3] = true;
  interval_map[I4] = true;

  for(auto&& pair : interval_map)
    std::cout << pair.first << ": " << pair.second << "\n";

  std::cout << "Before `erase_if`, there are "
            << interval_map.size()
            << " intervals"
            << std::endl;

  erase_if(interval_map,
           [] (const Interval& I)
           {
             return I.start() == 2;
           }
  );

  for(auto&& pair : interval_map)
    std::cout << pair.first << ": " << pair.second << "\n";

  std::cout << "After `erase_if`, there are "
            << interval_map.size()
            << " intervals"
            << std::endl;
}
