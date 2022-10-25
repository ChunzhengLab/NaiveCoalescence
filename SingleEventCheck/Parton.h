#include <vector>
#include "TVector2.h"
using namespace std;

class Parton {
 private:
  unsigned int serial_number;
  float baryon_number;
  TVector2 position_vector;

 public:
  Parton() {
    serial_number = 9999;
    baryon_number = -999;
    position_vector.SetX(-999);
    position_vector.SetY(-999);
  };
  void SetSerialNumber(unsigned int n) { serial_number = n; }
  void SetBaryonNumber(float bn) { baryon_number = bn; }
  void SetPositionVector(float x, float y) { position_vector.Set(x, y); }
  void SetPositionVector(TVector2 pv) { position_vector = pv; }

  unsigned int GetSerialNumber() { return serial_number; }
  float GetBaryonNumber() { return baryon_number; }
  TVector2 GetPositionVector() { return position_vector; }
};