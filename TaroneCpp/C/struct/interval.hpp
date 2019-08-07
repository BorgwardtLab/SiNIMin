#ifndef INTERVAL_HPP__
#define INTERVAL_HPP__

#include <functional>

class Interval
{
public:
  using value_type = int;

  Interval() = default;

  Interval(value_type start, value_type end)
    : _start(start)
    , _end(end)
  {
  }

  value_type start() const { return _start; }
  value_type end() const { return _end; }

private:
  value_type _start = 0;
  value_type _end = 0;
};

namespace std
{
  template<> struct hash<Interval>
  {
    using value_type = Interval::value_type;

    std::size_t operator()(const Interval& interval) const
    {
      return hash<value_type>()(interval.start()) ^ hash<value_type>()(interval.end());
    }
  };
}

#endif
