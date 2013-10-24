#ifndef DISTORTION_H
#define DISTORTION_H

#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include "mymatrix.h"

/************************************************************************************
 * distortion.h
 *   
 * �c�ݎړx���܂Ƃ߂��w�b�_
 *
 * Contents:
 *   cSNR, calcMSE
 *
 * Last Updated: <2013/10/18 18:17:56 from Yoshitos-iMac.local by yoshito>
 ************************************************************************************/

namespace mylib{
  class cSNR
  {
  protected:
    double sumSnr;
    int nFrame;

  public:
    cSNR(){
      reset();
    }
    // ~cSNR --- default destructor
  
    void reset()
    {
      sumSnr = 0.0;
      nFrame = 0;
    }
    // SNR��Ԃ�
    double calc(itpp::vec &x, itpp::vec &y);
    double calc(itpp::ivec &x, itpp::ivec &y);
    double calc(double x[], double y[], int sSamples);
    double calc(std::vector<double> &x, std::vector<double> &y);

    // ����܂Ōv�Z���Ă���SNR���t���[�����Ŋ�����segmental SNR��Ԃ�
    double getSegmentalSNR();
  };

  // system()�֐�����PQevalAudio���Ăяo���āAODG��Ԃ�
  double CalcPEAQ(const char* referenceName, const char *objectName);

  // wav�t�@�C������snr�����߂�
  // �Ȓ��͒Z�����ɍ��킹��
  double CalcSNR(const char *xFileName, const char *yFileName, const int frameLength);
  // The other name
  double CalcSDR(const char *xFileName, const char *yFileName, const int frameLength);

  // mean square error
  double CalcMSE(const std::vector<double> &original, const std::vector<double> &object);

  // ----- 2D version ------
  
  double CalcMSE(const vec_2D &original, const vec_2D &object);
 
  // u_char version
  double CalcMSE(const Vector_2D< u_char > &original, const Vector_2D< u_char > &object);

  double CalcMSE(const itpp::Mat< u_char > &original, const itpp::Mat< u_char > &object);
  
  // peak is the maximum possible value of pixel.
  double CalcPSNR(const vec_2D &original, const vec_2D &object, double peak);
  
  // u_char version
  double CalcPSNR(const Vector_2D< u_char > &original, const Vector_2D< u_char > &object,
                         double peak);

  // u_char version
  double CalcPSNR(const itpp::Mat< u_char > &original, const itpp::Mat< u_char > &object,
                  double peak);

  
  
} // namespace mylib

#endif
