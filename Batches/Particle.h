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

class Hardon {
 protected:
  unsigned int serial_number;
  int baryon_number;
  TVector2 position_vector;
  float mean_distance;
  vector<unsigned int> vec_partons_serial_number;
  vector<float> vec_partons_baryon_number;
  vector<TVector2> vec_partons_position;
 public:
  Hardon() {
    serial_number = -999;
    baryon_number = 9999;
    position_vector.SetX(-999);
    position_vector.SetY(-999);
    mean_distance = -999;
    vector<unsigned int>().swap(vec_partons_serial_number);
    vector<float>().swap(vec_partons_baryon_number);
    vector<TVector2>().swap(vec_partons_position);
  }
  void SetSerialNumber(unsigned int n) { serial_number = n; }
  void SetBaryonNumber(int bn) { baryon_number = bn; }
  void SetPositionVector(float x, float y) { position_vector.Set(x, y); }
  void SetPositionVector(TVector2 pv) { position_vector = pv; }
  void SetMeanDistance(float dis) { mean_distance = dis; }
  void SetVecPartonsSerialNumber(vector<unsigned int> sm) { vec_partons_serial_number.assign(sm.begin(), sm.end()); }
  void SetVecPartonsBaryonNumber(vector<float> bn) { vec_partons_baryon_number.assign(bn.begin(), bn.end()); }
  void SetVecPartonsPosition(vector<TVector2> vecpv) { vec_partons_position.assign(vecpv.begin(), vecpv.end()); }

  unsigned int GetSerialNumber() { return serial_number; }
  int GetBaryonNumber() { return baryon_number; }
  TVector2 GetPositionVector() { return position_vector; }
  float GetDistance() { return mean_distance; }
  vector<unsigned int> GetVecPartonsSerialNumber() { return vec_partons_serial_number; }
  vector<float> GetVecPartonsBaryonNumber() { return vec_partons_baryon_number; }
  vector<TVector2> GetVecPartonsPosition() { return vec_partons_position; }
};

class HardonCandidate : public Hardon {
 public:
  friend bool operator<(const HardonCandidate& a, const HardonCandidate& b) { return a.mean_distance < b.mean_distance; }
};

