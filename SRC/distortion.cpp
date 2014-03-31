#include <cassert>
#include "../include/wav.h"
#include "../include/distortion.h"

namespace mylib{
  const std::string PQ_COMMAND("PQevalAudio");
  const std::string TEMP_FILENAME("_temp.dat");

  const double DELTA = pow(1.0,-10000); // escape from 0 division.

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
  
    double snr = 10.0*log10(1+signal/noise); // 1は何かの補正のため

    // メンバ変数に加える
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
  
    double snr = 10.0*log10(1+signal/noise); // 1は何かの補正のため

    // メンバ変数に加える
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
  
    double snr = 10.0*log10(1+signal/noise); // 1は何かの補正のため

    // メンバ変数に加える
    sumSnr += snr;
    nFrame++;

    return snr;
  }

  double cSNR::getSegmentalSNR()
  {
    double segmentalSNR = sumSnr / static_cast<double>(nFrame);

    return segmentalSNR;
  }


  // PQevalAudioがインストールしてあることがマスト！
  double CalcPEAQ(const char *reference, const char *object)
  {
    std::string referenceName(reference);
    std::string objectName(object);

    std::string pqCommand(PQ_COMMAND.c_str());

    pqCommand += ' ';
    pqCommand += referenceName;
    pqCommand += ' ';
    pqCommand += objectName;
    pqCommand += " > ";
    pqCommand += TEMP_FILENAME;
  
    // コマンド実行
    if(system(pqCommand.c_str()) != 0){
      std::cerr << "Error at calcPEAQ.\n"
                << "Can not execute \"" << pqCommand << "\"\n";
      exit(1);
    }

    // これでtemp.datが出来て、中にPEAQ ODGの値がある
    std::ifstream tempFile(TEMP_FILENAME.c_str());
    assert(!tempFile.fail());

    tempFile.seekg(-7,std::ios::end);
  
    double odg;
    tempFile >> odg;

    tempFile.close();

    // 作成したtemp.datを削除
    std::string rmCommand("rm ");
  
    rmCommand += TEMP_FILENAME;

    if(system(rmCommand.c_str()) != 0){
      std::cerr << "Error at calcPEAQ.\n"
                << "Can not execute \"" << rmCommand << "\"\n";
      exit(1);
    }

    return odg;
  }
  
  // WAVファイル同士のSNRを求める
  double CalcSNR(const char *xFileName, const char *yFileName, const int frameLength)
  {
    cWavRead xFile(xFileName);
    cWavRead yFile(yFileName);

    assert(xFile.channels() == yFile.channels());
    assert(xFile.bitsPerSample() == yFile.bitsPerSample());
    assert(xFile.samplingFreq() == yFile.samplingFreq());

    cSNR SNR;
  
    bool bFinalFrame = false;
    while(1){
      int xFrameLength, yFrameLength;
      itpp::ivec x, y;
      xFrameLength = xFile.getWave(x, frameLength*xFile.channels());
      yFrameLength = yFile.getWave(y, frameLength*yFile.channels());
      // frameLength分取れなかったら最後のフレームということ
      if(frameLength*static_cast<int>(xFile.channels()) != xFrameLength){
        bFinalFrame = true;
      }
      if(frameLength*static_cast<int>(yFile.channels()) != yFrameLength){
        bFinalFrame = true;
      }
      // 長さが短い方に合わせる
      if(xFrameLength > yFrameLength){
        x = x.left(yFrameLength);
      }
      else if(xFrameLength < yFrameLength){
        y = y.left(xFrameLength);
      }

      SNR.calc(x,y);

      if(bFinalFrame){
        break;
      }
    } // while 1

    return SNR.getSegmentalSNR();
  }

  // Just the other name
  double CalcSDR(const char *xFileName, const char *yFileName, const int frameLength)
  {
    return CalcSNR(xFileName, yFileName, frameLength);
  }

  // mean square error
  double CalcMSE(const std::vector<double> &original, const std::vector<double> &object)
  {
    assert(original.size() == object.size());

    double mse = 0.0;
    for(int i = 0; i < static_cast<int>(original.size()); i++){
      mse += pow(object[i] - original[i], 2);
    }

    mse /= original.size();

    return mse;
  }

  double CalcMSE(const vec_2D &original, const vec_2D &object)
  {
    assert(original.is_rectangular() && object.is_rectangular());
    assert(original.size_rows() == object.size_rows() && original.size_cols() == object.size_cols());

    int vSize = original.size_rows();
    int hSize = original.size_cols();
    double mse = 0.0;
    for(int v = 0; v < vSize; v++)
      {
        for(int h = 0; h < hSize; h++)
          {
            mse += (original(v,h) - object(v,h)) * (original(v,h) - object(v,h));
          }
      }
    mse /= static_cast<double>(vSize * hSize);
    return mse;
  } 

  double CalcMSE(const Vector_2D< u_char > &original, const Vector_2D< u_char > &object)
  {
    assert(original.is_rectangular() && object.is_rectangular());
    assert(original.size_rows() == object.size_rows() && original.size_cols() == object.size_cols());

    int vSize = original.size_rows();
    int hSize = original.size_cols();
    double mse = 0.0;
    for(int v = 0; v < vSize; v++)
      {
        for(int h = 0; h < hSize; h++)
          {
            u_char diff = (original(v,h) >= object(v,h)) ? (original(v,h) - object(v,h)) :
              (object(v,h) - original(v,h));
            mse += static_cast<double>(diff*diff);
          }
      }
    mse /= static_cast<double>(vSize * hSize);

    return mse;
  }

  double CalcMSE(const itpp::Mat< u_char > &original, const itpp::Mat< u_char > &object)
  {
    assert(original.rows() == object.rows() && original.cols() == object.cols());

    int vSize = original.rows();
    int hSize = original.cols();
    double mse = 0.0;
    for(int v = 0; v < vSize; v++)
      {
        for(int h = 0; h < hSize; h++)
          {
            u_char diff = (original(v,h) >= object(v,h)) ? (original(v,h) - object(v,h)) :
              (object(v,h) - original(v,h));
            mse += static_cast<double>(diff*diff);
          }
      }
    mse /= static_cast<double>(vSize * hSize);

    return mse;
  }
  
  double CalcPSNR(const vec_2D &original, const vec_2D &object, double peak)
  {
    assert(original.is_rectangular() && object.is_rectangular());
    assert(original.size_rows() == object.size_rows() && original.size_cols() == object.size_cols());

    int vSize = original.size_rows();
    int hSize = original.size_cols();
    double mse = 0.0;
    for(int v = 0; v < vSize; v++)
      {
        for(int h = 0; h < hSize; h++)
          {
            assert(original(v,h) <= peak && object(v,h) <= peak);
            mse += (original(v,h) - object(v,h)) * (original(v,h) - object(v,h));
          }
      }
    mse /= static_cast<double>(vSize * hSize);
    if(mse == 0)
      {
        mse = 1.0/INFTY;
      }
    double psnr = 10 * log10( peak * peak / mse);
    return psnr;
  }

  double CalcPSNR(const Vector_2D< u_char > &original, const Vector_2D< u_char > &object,
                  double peak)
  {
    assert(original.is_rectangular() && object.is_rectangular());
    assert(original.size_rows() == object.size_rows() && original.size_cols() == object.size_cols());

    int vSize = original.size_rows();
    int hSize = original.size_cols();
    double mse = 0.0;
    for(int v = 0; v < vSize; v++)
      {
        for(int h = 0; h < hSize; h++)
          {
            assert(original(v,h) <= peak && object(v,h) <= peak);
            u_char diff = (original(v,h) >= object(v,h)) ? (original(v,h) - object(v,h)) :
              (object(v,h) - original(v,h));
            mse += static_cast<double>(diff*diff);
          }
      }
    mse /= static_cast<double>(vSize * hSize);
    if(mse == 0)
      {
        mse = 1.0 / INFTY;
      }
    double psnr = 10 * log10( peak * peak / mse);
    return psnr;
  }

  // u_char version
  double CalcPSNR(const itpp::Mat< u_char > &original, const itpp::Mat< u_char > &object,
                  double peak)
  {
    assert(original.rows() == object.rows() && original.cols() == object.cols());

    int vSize = original.rows();
    int hSize = original.cols();
    double mse = 0.0;
    for(int v = 0; v < vSize; v++)
      {
        for(int h = 0; h < hSize; h++)
          {
            assert(original(v,h) <= peak && object(v,h) <= peak);
            u_char diff = (original(v,h) >= object(v,h)) ? (original(v,h) - object(v,h)) :
              (object(v,h) - original(v,h));
            mse += static_cast<double>(diff*diff);
          }
      }
    mse /= static_cast<double>(vSize * hSize);
    if(mse == 0)
      {
        mse = 1.0 / INFTY;
      }
    double psnr = 10 * log10( peak * peak / mse);
    return psnr;
  }
  
}

