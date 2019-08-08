#ifndef INTERVAL_HPP__
#define INTERVAL_HPP__

#include <functional>
#include <iostream>

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

  value_type start() const noexcept { return _start; }
  value_type end() const noexcept { return _end; }

  bool contains(value_type x) const noexcept
  {
    return _start <= x && _x <= _end;
  }

  bool operator==(const Interval& other) const noexcept
  {
    return _start == other._start && _end == other._end;
  }

  bool operator!=(const Interval& other) const noexcept
  {
    return !this->operator==(other);
  }

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

  ostream& operator<<(ostream& o, const Interval& interval)
  {
    o << "[" << interval.start() << ":" << interval.end() << "]";
    return o;
  }
}

#endif
