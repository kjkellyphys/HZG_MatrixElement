#ifndef _LORENTZ_VECTOR_HH_
#define _LORENTZ_VECTOR_HH_

#include <cmath>

class LorentzVector
{
 public:

  double px,py,pz,E;

  LorentzVector(): px(0), py(0), pz(0), E(0) {}
  LorentzVector(double p1, double p2, double p3, double p0): px(p1), py(p2), pz(p3), E(p0) {}
  LorentzVector(const LorentzVector& p): px(p.px), py(p.py), pz(p.pz), E(p.E) {}
  double p()   const {return sqrt(px*px + py*py + pz*pz);}
  double pt()  const {return sqrt(px*px + py*py);}
  double y()   const {return 0.5*log((E + pz)/(E - pz));}
  double Et()  const {return E/p()*pt();}
  double eta() const {return 0.5*log((p() + pz)/(p() - pz));}
  double phi() const
  {
    double r = atan2(py,px);
    if(r < 0)
      r += 2*M_PI;
    return r;
  }
  void add(LorentzVector v)
  {
    px += v.px;
    py += v.py;
    pz += v.pz;
    E  += v.E;
  }
  void subtract(LorentzVector v)
  {
    px -= v.px;
    py -= v.py;
    pz -= v.pz;
    E  -= v.E;
  }
  bool isEqual(LorentzVector v)
  {
    return px == v.px && py == v.py && pz == v.pz && E == v.E;
  }
};

#endif
