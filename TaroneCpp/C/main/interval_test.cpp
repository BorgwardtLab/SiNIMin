#include "../struct/interval.hpp"

#include <iostream>
#include <unordered_map>

int main(int, char**)
{
  Interval I1;
  Interval I2(2, 3);

  std::unordered_map<Interval, bool> interval_map;
  interval_map[I2] = true;
  interval_map[I1] = false;

  for(auto&& pair : interval_map)
    std::cout << pair.first << ": " << pair.second << "\n";
}
