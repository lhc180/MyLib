/*******************************/
// Segmental SNRを求めるクラス
/******************************/

#ifndef SNR_H
#define SNR_H

#include <itpp/itbase.h>
#include <itpp/itcomm.h>

class cSNR
{
protected:
  double sumSnr;
  int nFrame;

public:
  cSNR();
  // ~cSNR --- default destructor
  
  void reset();
  // SNRを返す
  double calc(itpp::vec &x, itpp::vec &y);
  double calc(itpp::ivec &x, itpp::ivec &y);
  double calc(double x[], double y[], int sSamples);
  double calc(std::vector<double> &x, std::vector<double> &y);

  // これまで計算してきたSNRをフレーム数で割ってsegmental SNRを返す
  double getSegmentalSNR();
};

#endif
