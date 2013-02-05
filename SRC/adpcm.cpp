/***************************************/
// delay free ADPCM の定義
// このソースコードは以下の文献で提案されたADPCMを実装してある
// DELAY-FREE AUDIO CODING BASED ON ADPCM AND ERROR FEEDBACK
// 変数や配列の名前は論文を参照
/****************************************/
/*************to do list****************/
// noise shaperの入力のxHatをeHatかpにしてみる
/***************************************/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include "../include/adpcm.h"
#include <itpp/itbase.h>

const double MU_TILDE = 0.002;	// base step size
const double L_MIN = 0.000000001;	// 式(7)で分母が0にならないようにする値
const double EPSILON = pow(10,-8);


template<typename type>
inline type maxOf(type x1, type x2)
{
  if(x1 > x2){
    return x1;
  }
  else{
    return x2;
  }
}

/**********adpcmDelayFree class*****************/

void adpcmDelayFree::encoding(const itpp::vec &inputData, itpp::vec &encodedData)
{
  encodedData.set_size(inputData.size());
  
  double v = 1.0;
  double e,eHat;
  double q = 0.0;
  double xHat = 0.0;
  double p = 0.0;
  
  start = true;
  
  for(int n = 0; n < inputData.size(); n++){
    e = inputData[n] - p;
    e /= v;

    /**********ここで量子化してqを得る**************/
   
    q = Quantizer.quantize(e);
    /*******************************************/
  
    eHat = q*v;

    xHat = eHat + p;

    /***************vの更新*************************/
      
    v = envelopeEstimation(v, eHat);
    
    /***********************************************/

    /**************pの更新**************************/
    
    p = predictor(xHat);

    /***************************************************/

    encodedData[n] = q;
    if(start){
      start = false;
    }
  }
}
  
// indexのベクトルを返す
// encodedDataがivecなのに注意
void adpcmDelayFree::encoding(const itpp::vec &inputData, itpp::ivec &encodedData)
{

  encodedData.set_size(inputData.size());


  double v = 1.0;
  double e,eHat;
  double q = 0.0;
  double xHat = 0.0;
  double p = 0.0;
  
  start = true;
  
  for(int n = 0; n < inputData.size(); n++){
    e = inputData[n] - p;
    e /= v;

    /**********ここで量子化してqを得る**************/
    int index = Quantizer.quantizeReturnIndex(e);

    q = Quantizer.quantize(e);
    /*******************************************/
  
    eHat = q*v;

    xHat = eHat + p;

    /***************vの更新*************************/
      
    v = envelopeEstimation(v, eHat);
    
    /***********************************************/

    /**************pの更新**************************/
    
    p = predictor(xHat);

    /***************************************************/

    encodedData[n] = index;
    if(start){
      start = false;
    }
  }
  

}


void adpcmDelayFree::decoding(const itpp::vec &encodedData, itpp::vec &decodedData)
{
  double eHat;
  double v = 1.0;
  double p = 0.0;
  double xHat;

  decodedData.set_size(encodedData.size());

  start = true;
  
  // std::cout << "adpcmDelayFree::decoding is called.\n";
  

  for(int n = 0; n < encodedData.size(); n++){
    eHat = encodedData[n] * v;
    
    xHat = eHat + p;

    // std::cout << "## eHat = " << eHat << std::endl;
    // vの更新
    v = envelopeEstimation(v,eHat);

    // std::cout << "## v = " << v << '\n';
    
    // pの更新
    p = predictor(xHat);

    //std::cout << "## p = " << p << '\n';
    

    decodedData[n] = xHat;
    
    if(start){
      start = false;
    }
  }
}
    
// indexのベクトルから復号データを作る
void adpcmDelayFree::decoding(const itpp::ivec &encodedData, itpp::vec &decodedData)
{
  decodedData.set_size(encodedData.size());

  double eHat;
  double v = 1.0;
  double p = 0.0;
  double xHat;

  start = true;

  //  std::cout << "adpcmDelayFree::decoding is called.\n";
  
  for(int n = 0; n < encodedData.size(); n++){
    double encoded = Quantizer.dequantize(encodedData[n]);
    eHat = encoded * v;
    
    xHat = eHat + p;

    // std::cout << "## eHat = " << eHat << std::endl;
    // vの更新
    v = envelopeEstimation(v,eHat);

    // std::cout << "## v = " << v << '\n';
    
    // pの更新
    p = predictor(xHat);

    decodedData[n] = xHat;
    
    if(start){
      start = false;
    }
  }
}
    

double adpcmDelayFree::envelopeEstimation(const double v, const double eHat)
{
  double lambda;
  static double vMin = L_MIN;
  double vOut;

  // 初期化
  if(start == true){
    vMin = L_MIN;
    // start = false;
  }

  if(pow(eHat,2) > pow(v,2)){
    lambda = lambdaAT;
  }
  else{
    lambda = lambdaRT;
  }
      
  vOut = sqrt(maxOf(pow(vMin,2), lambda*pow(eHat,2) + (1-lambda)*pow(v,2)));
  // if(v < vMin){
  //     vMin = v;
  //   }

  

  return vOut;
}
   

  // 
double adpcmDelayFree::predictor(const double xHat) 
{
  double p;
  itpp::vec f(delay+1),b(delay+1);
  itpp::vec mu(delay+1),l(delay+1); // テキストと係数を合わせるため+1してある
  // 0番目は捨て
  static itpp::vec b1Delayed(delay+1);
  
  // 初期化
  if(start == true){
    b1Delayed.zeros();
    kappa.zeros();
    l.zeros();
    // start = false;
  }

  f[0] = xHat;
  b[0] = xHat;
  
  
  for(int m = 1; m <= delay; m++){
    f[m] = f[m-1] + kappa[m]*b1Delayed[m-1]; // ## b1Delayed
    b[m] = b1Delayed[m-1] + kappa[m]*f[m-1]; // ## b1Delayed
  }

  // std::cout << "## f = " << f << std::endl;
  // std::cout << "## b = " << b << std::endl;

  // pを求める
  /* p = 0.0;
     for(int m = 1; m <= delay; m++){
     p += kappa[m]*b1Delayed[m-1];	// ## b1Delayed
     }*/
  
  p = xHat - f[delay];		// ##


  // kappaを更新する
  // 論文の式(6)におけるbm(n)は間違えの可能性があるので、bm-1(n)としてとりあえずやる
  for(int m = 1; m <= delay; m++){
    l[m] = (MU_TILDE)*l[m] + (1-MU_TILDE)*(pow(f[m-1],2)+pow(b1Delayed[m-1],2)); // ## b1Delayed
    // std::cout << l[m] << " ";
    if(l[m] < L_MIN){
      l[m] = L_MIN;
      // std::cout << "sgn(l) = " << itpp::sgn(l[m]) << std::endl;
    }
    mu[m] = MU_TILDE / (l[m]);	// この辺はテキストと変えた
    kappa[m] = kappa[m] - mu[m]*(f[m]*b1Delayed[m-1] + f[m-1]*b[m]); // ## 1項目b1Delayed
    if(fabs(kappa[m]) > 1 - EPSILON){
      kappa[m] = itpp::sgn(kappa[m])*(1-EPSILON);
    }

  }
  
  // std::cout << "## kappa = " << kappa << std::endl;
  // std::cout << "## l = " << l << std::endl;

 

  for(int m = 0; m <= delay; m++){
    b1Delayed[m] = b[m];
  }

  return p;
    
}

// double adpcmDelayFree::quantize(const double input)
// {
//   int index = 0;
//   double minimumError = 100.0;
//   double error;
  
//   for(int i = 0; i < quantizeVec.size(); i++){
//     error = pow(input - quantizeVec[i],2);
//     if(error < minimumError){
//       index = i;
//       minimumError = error;
//     }
//   }

//   return quantizeVec[index];
// }

// int adpcmDelayFree::quantizeReturnIndex(const double input)
// {
//   int index = 0;
//   double minimumError = 100.0;
//   double error;
  
//   for(int i = 0; i < quantizeVec.size(); i++){
//     error = pow(input - quantizeVec[i],2);
//     if(error < minimumError){
//       index = i;
//       minimumError = error;
//     }
//   }

//   return index;
// }
/************end of adpcmDelayFree class*****************/

/*************adpcmWithNoiseShaper class*****************/
// このソースコードは以下の文献で提案されたADPCMを実装してある
// DELAY-FREE AUDIO CODING BASED ON ADPCM AND ERROR FEEDBACK
// noise shaperを継承を使って実装する
// 変数や配列の名前は論文を参照
//
// GAL algorithmによる更新式は以下の論文を参照した
// Simulation and Performance Analysis of Adaptive Filtering Algorithms in Noise Cancellation



const double BETA = 0.5;	// noise shaperで用いる定数

void adpcmWithNoiseShaper::encoding(const itpp::vec &inputData, itpp::vec &encodedData)
{
  encodedData.set_size(inputData.size());
  
  double v = 1.0;
  double e,eHat;
  double q = 0.0;
  double xHat = 0.0;
  double p = 0.0;
  double s = 0.0;
  double u = 0.0;
  double x;
  
  start = true;
  
  for(int n = 0; n < inputData.size(); n++){
    u = inputData[n];
    
    x = u + s;

    e = x - p;
    e /= v;

    /**********ここで量子化してqを得る**************/
    q = Quantizer.quantize(e);
    /*******************************************/
  
    eHat = q*v;

    xHat = eHat + p;

    /***************vの更新*************************/
      
    v = envelopeEstimation(v, eHat);
    // std::cout << "## vOut = " << v << std::endl;
    /***********************************************/

    /**************pの更新**************************/
    
    p = predictor(xHat);

    /***************************************************/

    /***************noise shaper*********************/
    s = noiseShaper(u, x, xHat);

    encodedData[n] = q;
    
    if(start){
      start = false;
    }
  }
}
 
// indexのベクトルを返す
void adpcmWithNoiseShaper::encoding(const itpp::vec &inputData, itpp::ivec &encodedData)
{
  encodedData.set_size(inputData.size());
  
  double v = 1.0;
  double e,eHat;
  double q = 0.0;
  double xHat = 0.0;
  double p = 0.0;
  double s = 0.0;
  double u = 0.0;
  double x;
  
  start = true;
  
  for(int n = 0; n < inputData.size(); n++){
    u = inputData[n];
    
    x = u + s;

    e = x - p;
    e /= v;

    /**********ここで量子化してqを得る**************/
    int index = Quantizer.quantizeReturnIndex(e);

    q = Quantizer.quantize(e);
    /*******************************************/
  
    eHat = q*v;

    xHat = eHat + p;

    /***************vの更新*************************/
      
    v = envelopeEstimation(v, eHat);
    

    /***********************************************/

    /**************pの更新**************************/
    
    p = predictor(xHat);

    /***************************************************/

    /***************noise shaper*********************/
    s = noiseShaper(u, x, xHat);

    encodedData[n] = index;
    
    if(start){
      start = false;
    }
  }
}
 
  

double adpcmWithNoiseShaper::noiseShaper(const double u, const double x, const double xHat)
{
  static itpp::vec xVec(hs.size()),xHatVec(hs.size()); // 0は現在の信号
  itpp::vec f(delay), b(delay);
  static itpp::vec b1Delayed(delay);
  
  if(start){
    xVec.zeros();
    xHatVec.zeros();
    b1Delayed.zeros();
  }

  xVec[0] = x;
  xHatVec[0] = xHat;

  f[0] = xHat - u;
  b[0] = BETA * (xHat - u);

  //  std::cout << "## kappa = " << kappa << std::endl;
  // ここは+-逆かも
  for(int m = 0; m < delay; m++){
    f[m] = f[m-1] - kappa[m]*b1Delayed[m-1];
    b[m] = BETA * (b1Delayed[m-1] - kappa[m]*f[m-1]);
  }

  //  std::cout << "## f = " << f << std::endl;
  // std::cout << "## b = " << b << std::endl;

  double sa = 0.0;
  for(int m = 1; m <= delay; m++){
    sa += kappa[m] * b1Delayed[m-1];
  }
  
  double ss = 0.0;
  for(int l = 1; l < hs.size(); l++){
    ss += hs[l] * (xHatVec[l] - xVec[l]);
  }

  double s = sa - ss;

  // 係数の更新
  for(int m = 0; m < delay; m++){
    b1Delayed[m] = b[m];
  }
  
  for(int l = hs.size() - 1; l >= 1; l--){
    xVec[l] = xVec[l-1];
    xHatVec[l] = xHatVec[l-1];
  }

  return s;
}

  /*****************end of adpcmWithNoiseShaper class***********/

  /****************adpcmLMS class************************/

  const double MU = 0.1;	// base step size

void adpcmLMS::encoding(const itpp::vec &inputData, itpp::vec &encodedData)
{
  encodedData.set_size(inputData.size());
  
  double v = 1.0;
  double e,eHat;
  double q = 0.0;
  double xHat = 0.0;
  double p = 0.0;
  
  start = true;
  
  for(int n = 0; n < inputData.size(); n++){
    e = inputData[n] - p;
    e /= v;

    /**********ここで量子化してqを得る**************/
    q = Quantizer.quantize(e);
    /*******************************************/
  
    eHat = q*v;

    xHat = eHat + p;

    /***************vの更新*************************/
      
    v = envelopeEstimation(v, eHat);
    
    /***********************************************/

    /**************pの更新**************************/
    
    p = predictor(xHat,eHat);

    /***************************************************/

    encodedData[n] = q;

    if(start){
      start = false;
    }
  }
}

// indexのベクトルを返す
void adpcmLMS::encoding(const itpp::vec &inputData, itpp::ivec &encodedData)
{
  encodedData.set_size(inputData.size());
  
  double v = 1.0;
  double e,eHat;
  double q = 0.0;
  double xHat = 0.0;
  double p = 0.0;
  
  start = true;
  
  for(int n = 0; n < inputData.size(); n++){
    e = inputData[n] - p;
    e /= v;

    /**********ここで量子化してqを得る**************/
    int index = Quantizer.quantizeReturnIndex(e);
    q = Quantizer.quantize(e);
    /*******************************************/
  
    eHat = q*v;

    xHat = eHat + p;

    /***************vの更新*************************/
      
    v = envelopeEstimation(v, eHat);
    
    /***********************************************/

    /**************pの更新**************************/
    
    p = predictor(xHat,eHat);

    /***************************************************/

    encodedData[n] = index;

    if(start){
      start = false;
    }
  }
}
  
void adpcmLMS::decoding(const itpp::vec &encodedData, itpp::vec &decodedData)
{
  double eHat;
  double v = 1.0;
  double p = 0.0;
  double xHat;

  start = true;
  
  for(int n = 0; n < encodedData.size(); n++){
    eHat = encodedData[n] * v;
    
    xHat = eHat + p;

    // std::cout << "## eHat = " << eHat << std::endl;
    // vの更新
    v = envelopeEstimation(v,eHat);

    // std::cout << "## v = " << v << '\n';
    
    // pの更新
    p = predictor(xHat,eHat);

    decodedData[n] = xHat;
    
    if(start){
      start = false;
    }
  }
}

void adpcmLMS::decoding(const itpp::ivec &encodedData, itpp::vec &decodedData)
{
  double eHat;
  double v = 1.0;
  double p = 0.0;
  double xHat;

  start = true;
  
  for(int n = 0; n < encodedData.size(); n++){
    eHat = Quantizer.dequantize(encodedData[n]) * v;
    
    xHat = eHat + p;

    // std::cout << "## eHat = " << eHat << std::endl;
    // vの更新
    v = envelopeEstimation(v,eHat);

    // std::cout << "## v = " << v << '\n';
    
    // pの更新
    p = predictor(xHat,eHat);

    decodedData[n] = xHat;
    
    if(start){
      start = false;
    }
  }
}
    

   

// 
double adpcmLMS::predictor(const double xHat, const double eHat) 
{
  static itpp::vec xHatVec(delay+1), eHatVec(delay+1), h(delay+1);
  double p;

  if(start == true){
    xHatVec.zeros();
    eHatVec.zeros();
    h.zeros();
    // start = false;
  }
  
  xHatVec[0] = xHat;
  eHatVec[0] = eHat;


  for(int m = 1; m <= delay; m++){
    h[m] = h[m] + MU * eHatVec[m] * xHatVec[m];
  }

  p = 0.0;
  for(int m = 1; m <= delay; m++){
    p += h[m] * xHatVec[m];
  }



  for(int m = delay; m >= 1; m--){
    xHatVec[m] = xHatVec[m-1];
    eHatVec[m] = eHatVec[m-1];
  }

  return p;
}

  /**************end of adpcmLMS class**************/

