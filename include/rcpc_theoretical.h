#ifndef RCPC_THEORETICAL_H
#define RCPC_THEORETICAL_H

#include<cstdlib>
#include<assert.h>
#include<cmath>
/**************************************
 * RCPCの理論値を格納したクラス
 *
 * パンクチャリング周期（符号化率の分子）と符号化率の分母をコンストラクタで与える。
 * getBERの引数で所望のchannel SNRを与える
 *------------------------------------------
 *thrRCPC
 *
 *       コンストラクタ -- 引数で
 * 
 *
 *************************************/

class thrRCPC{
public:
  thrRCPC():puncturing_period(1), codingRate_denom(2){}			// デフォルトコンストラクタ
  thrRCPC(const unsigned int p, const unsigned int q)
  {
    set(p,q);
  } // コンストラクタ
  ~thrRCPC()
  {
    puncturing_period = 0;
    codingRate_denom = 0;
    // numEbOverN0 = 0;
  }
					       // デストラクタ
  // デフォルトコピーコンストラクタ

  void set(const unsigned int p, const unsigned int q)
  {
    puncturing_period = p;
    codingRate_denom = q;
  }

  unsigned int getPuncturingPeriod() const
  {
    return puncturing_period;
  }

  unsigned int getCodingRate_denom() const
  {  
    return codingRate_denom;
  }

  // unsigned int getNumEbOverN0()
//   {
//     return numEbOverN0;
//   }

  // EsN0は真値
  double getBER(double EsN0) const
  {
    double ber = calcBER(EsN0);

    return ber;
  }

private:
  unsigned int puncturing_period; // パンクチャリング周期、符号化率の分子
  unsigned int codingRate_denom; // 符号化率の分母
  // unsigned int numEbOverN0;	 // berの要素数
  
  double calcBER(double EsN0) const;

};


  /******************************************
 setBER -- RCPC符号の理論値をセットする関数
 引数:
      ber -- 理論値を格納する配列
      puncturing_period -- パンクチャリング周期
      codingRate_denom -- 符号化率の分母
 注意:
      符号化率はpuncturing_period/codingRate_denomとなる
 ******************************************/
inline double thrRCPC::calcBER(double snr) const
{
  unsigned int d_free;		// 最小自由距離
  unsigned int *beta;		// β
  unsigned int num_d;		// Σで計算する回数
  
  // 各puncturing_periodによってパラメータを設定する。
  if(puncturing_period == 1 && codingRate_denom == 2){
    d_free = 5;
    num_d = 10;
    beta = new unsigned int[num_d];
    beta[0] = 1;
    beta[1] = 4;
    beta[2] = 12;
    beta[3] = 32;
    beta[4] = 80;
    beta[5] = 192;
    beta[6] = 448;
    beta[7] = 1024;
    beta[8] = 2304;
    beta[9] = 5120;
  }
  else if(puncturing_period == 2){
    switch(codingRate_denom){
    case 3:
      d_free = 3;
      num_d = 8;
      beta = new unsigned int[num_d];
      beta[0] = 1;
      beta[1] = 10;
      beta[2] = 54;
      beta[3] = 226;
      beta[4] = 856;
      beta[5] = 3072;
      beta[6] = 10647;
      beta[7] = 35998;
      break;
    
    default:
      exit(0);
      break;
    }
  }
  else if(puncturing_period == 3){
    switch(codingRate_denom){
    case 5:
      d_free = 4;
      num_d = 7;
      beta = new unsigned int[num_d];
      beta[0] = 4;
      beta[1] = 32;
      beta[2] = 104;
      beta[3] = 312;
      beta[4] = 961;
      beta[5] = 2861;
      beta[6] = 8278;
      break;

    case 4:
      d_free = 3;
      num_d = 6;
      beta = new unsigned int[num_d];
      beta[0] = 15;
      beta[1] = 104;
      beta[2] = 540;
      beta[3] = 2536;
      beta[4] = 11302;
      beta[5] = 48638;
      break;

    default:
      exit(0);
      break;
    }
  }
  else if(puncturing_period == 4){
    switch(codingRate_denom){
    case 5:
      d_free = 2;
      num_d = 7;
      beta = new unsigned int[num_d];
      beta[0] = 3;
      beta[1] = 54;
      beta[2] = 387;
      beta[3] = 2299;
      beta[4] = 12679;
      beta[5] = 66708;
      beta[6] = 339837;
      break;

    case 6:
      d_free = 4;
      num_d = 7;
      beta = new unsigned int[num_d];
      beta[0] = 66;
      beta[1] = 0;
      beta[2] = 1108;
      beta[3] = 0;
      beta[4] = 13836;
      beta[5] = 0;
      beta[6] = 154966;
      break;

    case 7:
      d_free = 4;
      num_d = 7;
      beta = new unsigned int[num_d];
      beta[0] = 4;
      beta[1] = 19;
      beta[2] = 72;
      beta[3] = 224;
      beta[4] = 651;
      beta[5] = 1850;
      beta[6] = 5090;
      break;

    default:
      exit(1);
      break;
    }
  }
  else if(puncturing_period == 5){
    switch(codingRate_denom){
    case 6:
      d_free = 2;
      num_d = 6;
      beta = new unsigned int[num_d];
      beta[0] = 6;
      beta[1] = 116;
      beta[2] = 1004;
      beta[3] = 7183;
      beta[4] = 47588;
      beta[5] = 300505;
      break;
      
    case 7:
      d_free = 3;
      num_d = 7;
      beta = new unsigned int[num_d];
      beta[0] = 6;
      beta[1] = 92;
      beta[2] = 375;
      beta[3] = 2148;
      beta[4] = 8448;
      beta[5] = 38025;
      beta[6] = 150243;
      break;
      
    case 8:
      d_free = 4;
      num_d = 7;
      beta = new unsigned int[num_d];
      beta[0] = 22;
      beta[1] = 85;
      beta[2] = 292;
      beta[3] = 1043;
      beta[4] = 3217;
      beta[5] = 10324;
      beta[6] = 31724;
      break;

    case 9:
      d_free = 4;
      num_d = 7;
      beta = new unsigned int[num_d];
      beta[0] = 2;
      beta[1] = 21;
      beta[2] = 66;
      beta[3] = 205;
      beta[4] = 578;
      beta[5] = 1551;
      beta[6] = 4290;
      break;

    default:
      exit(0);
      break;
    }
  }
  else if(puncturing_period == 6){
    switch(codingRate_denom){
    case 7:
      d_free = 2;
      num_d = 6;
      beta = new unsigned int[num_d];
      beta[0] = 13;
      beta[1] = 252;
      beta[2] = 2308;
      beta[3] = 18083;
      beta[4] = 133171;
      beta[5] = 938718;
      break;
      
    case 8:
      d_free = 3;
      num_d = 6;
      beta = new unsigned int[num_d];
      beta[0] = 30;
      beta[1] = 208;
      beta[2] = 1080;
      beta[3] = 5252;
      beta[4] = 24542;
      beta[5] = 110940;
      break;

    case 9:
      d_free = 3;
      num_d = 7;
      beta = new unsigned int[num_d];
      beta[0] = 2;
      beta[1] = 61;
      beta[2] = 170;
      beta[3] = 892;
      beta[4] = 2698;
      beta[5] = 11072;
      beta[6] = 35954;
      break;
      
    case 10:
      d_free = 4;
      num_d = 7;
      beta = new unsigned int[num_d];
      beta[0] = 8;
      beta[1] = 64;
      beta[2] = 208;
      beta[3] = 624;
      beta[4] = 1946;
      beta[5] = 6010;
      beta[6] = 18140;
      break;

    case 11:
      d_free = 4;
      num_d = 7;
      beta = new unsigned int[num_d];
      beta[0] = 2;
      beta[1] = 19;
      beta[2] = 68;
      beta[3] = 191;
      beta[4] = 530;
      beta[5] = 1516;
      beta[6] = 3965;
      break;

    default:
      exit(0);
      break;
    }
  }
  else if(puncturing_period == 7){
    switch(codingRate_denom){
    case 8:
      d_free = 2;
      num_d = 5;
      beta = new unsigned int[num_d];
      beta[0] = 28;
      beta[1] = 467;
      beta[2] = 4606;
      beta[3] = 39748;
      beta[4] = 323416;
      break;
      
    case 9:
      d_free = 2;
      num_d = 6;
      beta = new unsigned int[num_d];
      beta[0] = 4;
      beta[1] = 46;
      beta[2] = 540;
      beta[3] = 1838;
      beta[4] = 14334;
      beta[5] = 53728;
      break;

    case 10:
      d_free = 3;
      num_d = 7;
      beta = new unsigned int[num_d];
      beta[0] = 8;
      beta[1] = 130;
      beta[2] = 540;
      beta[3] = 1786;
      beta[4] = 9668;
      beta[5] = 37046;
      beta[6] = 139168;
      break;

    case 11:
      d_free = 4;
      num_d = 7;
      beta = new unsigned int[num_d];
      beta[0] = 36;
      beta[1] = 189;
      beta[2] = 340;
      beta[3] = 1866;
      beta[4] = 6416;
      beta[5] = 18464;
      beta[6] = 65900;
      break;

    case 12:
      d_free = 4;
      num_d = 7;
      beta = new unsigned int[num_d];
      beta[0] = 8;
      beta[1] = 48;
      beta[2] = 160;
      beta[3] = 522;
      beta[4] = 1504;
      beta[5] = 4650;
      beta[6] = 13647;
      break;

    case 13:
      d_free = 4;
      num_d = 8;
      beta = new unsigned int[num_d];
      beta[0] = 2;
      beta[1] = 17;
      beta[2] = 70;
      beta[3] = 206;
      beta[4] = 528;
      beta[5] = 1419;
      beta[6] = 3934;
      beta[7] = 10346;
      break;

    default:
      exit(1);
      break;
    }
  }
  else{
    std::cout << "Uncorrect value.\n";
    exit(1);
  }

    
  double probability = 0.0;
    
  for(unsigned int d=0;d<num_d;d++){
    double temp = (1.0/2.0)* erfc( sqrt( static_cast<double>(d_free+d) *  snr) );
    probability += static_cast<double>(beta[d])*temp;
  }
  probability /= static_cast<double>(puncturing_period);

  delete[] beta;
  beta = NULL;

  return probability;
}


#endif
