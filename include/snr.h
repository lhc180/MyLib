/*******************************/
// Segmental SNR�����߂�N���X
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
  // SNR��Ԃ�
  double calc(itpp::vec &x, itpp::vec &y);
  double calc(itpp::ivec &x, itpp::ivec &y);
  double calc(double x[], double y[], int sSamples);
  double calc(std::vector<double> &x, std::vector<double> &y);

  // ����܂Ōv�Z���Ă���SNR���t���[�����Ŋ�����segmental SNR��Ԃ�
  double getSegmentalSNR();
};

#endif
