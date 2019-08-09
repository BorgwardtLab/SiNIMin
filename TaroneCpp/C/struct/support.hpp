#ifndef SUPPORT_HPP__
#define SUPPORT_HPP__

#include <functional>
#include <iostream>


class Support
{
  public:

    // Default constructor.
    Support() = default;

    Support(double min_pv, double envelope, double pval)
    : _min_pv(min_pv)
    , _envelope(envelope)
    , _pval(pval)
    {
    }

    // Functions to print different descriptors.
    double min_pvalue() const noexcept { return _min_pv; }
    double envelope() const noexcept { return _envelope; }
    double pvalue() const noexcept { return _pval; }


  private:

    double _min_pv;
    double _envelope;
    double _pval;

};


namespace std
{
  ostream& operator<<(ostream& o, const Support& support)
    {
      o << support.min_pvalue() << ", " << support.envelope() << ", ";
      o << support.pvalue();
      return o;
    }
}

#endif
