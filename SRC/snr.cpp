#include <cassert>
#include <cmath>
#include "../include/snr.h"

const double DELTA = pow(1.0,-10000); // escape from 0 division.

cSNR::cSNR()
{
  reset();
}

void cSNR::reset()
{
  sumSnr = 0.0;
  nFrame = 0;
}

double cSNR::calc(itpp::vec &x, itpp::vec &y)
{
  assert(x.size() == y.size());

  double signal = 0.0;
  double noise = 0.0;
  
  for(int n = 0; n < x.size(); n++){
    signal += pow(x[n],2);
    noise += pow(x[n]-y[n],2);
  }

  if(noise == 0.0){
    noise = DELTA;
  }
  
  double snr = 10.0*log10(1+signal/noise); // 1ÇÕâΩÇ©ÇÃï‚ê≥ÇÃÇΩÇﬂ

  // ÉÅÉìÉoïœêîÇ…â¡Ç¶ÇÈ
  sumSnr += snr;
  nFrame++;

  return snr;
}

double cSNR::calc(itpp::ivec &x, itpp::ivec &y)
{
  itpp::vec xDouble = itpp::to_vec(x);
  itpp::vec yDouble = itpp::to_vec(y);

  double snr = calc(xDouble,yDouble);

  return snr;
}

double cSNR::calc(double x[], double y[], int nSamples)
{
  double signal = 0.0;
  double noise = 0.0;
  
  for(int n = 0; n < nSamples; n++){
    signal += pow(x[n],2);
    noise += pow(x[n]-y[n],2);
  }

  if(noise == 0.0){
    noise = DELTA;
  }
  
  double snr = 10.0*log10(1+signal/noise); // 1ÇÕâΩÇ©ÇÃï‚ê≥ÇÃÇΩÇﬂ

  // ÉÅÉìÉoïœêîÇ…â¡Ç¶ÇÈ
  sumSnr += snr;
  nFrame++;

  return snr;
}

double cSNR::calc(std::vector<double> &x, std::vector<double> &y)
{
  assert(x.size() == y.size());

  double signal = 0.0;
  double noise = 0.0;
  
  for(int n = 0; n < static_cast<int>(x.size()); n++){
    signal += pow(x[n],2);
    noise += pow(x[n]-y[n],2);
  }

  if(noise == 0.0){
    noise = DELTA;
  }
  
  double snr = 10.0*log10(1+signal/noise); // 1ÇÕâΩÇ©ÇÃï‚ê≥ÇÃÇΩÇﬂ

  // ÉÅÉìÉoïœêîÇ…â¡Ç¶ÇÈ
  sumSnr += snr;
  nFrame++;

  return snr;
}

double cSNR::getSegmentalSNR()
{
  double segmentalSNR = sumSnr / static_cast<double>(nFrame);

  return segmentalSNR;
}

