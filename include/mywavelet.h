#ifndef MYWAVELET_H
#define MYWAVELET_H

#include <iostream>
#include <cstdlib>
#include <gsl/gsl_wavelet2d.h>
#include <itpp/itbase.h>
#include "myutl.h"
#include "mymatrix.h"

/************************************************************************************
 * Wavelet
 *
 * ウェーブレット変換と変換領域を1次元配列にするための関数を含む
 *
 * Contents:
 *   Wavelet2D
 *   To1D_forWVT
 *   To2D_forWVT
 *
 * Last Updated: <2014/04/15 22:36:16 from Okauchi.local by yoshito>
 ************************************************************************************/

// =============== To Do List ===============
// -- ヘッダと実装ファイルの分割
// ==========================================


namespace mylib{
  enum WaveletType{
    WAVELET_DAUBECHIES,
    WAVELET_DAUBECHIES_CENTERED,
    WAVELET_HAAR,
    WAVELET_HAAR_CENTERED,
    WAVELET_BSPLINE,
    WAVELET_BSPLINE_CENTERED
  };

  class Wavelet2D
  {
  private:
    gsl_wavelet *wavelet_;

    void Alloc(int k, WaveletType waveletType)
    {
      switch (waveletType){
      case WAVELET_DAUBECHIES:
        wavelet_ = gsl_wavelet_alloc(gsl_wavelet_daubechies, k);
        break;

      case WAVELET_DAUBECHIES_CENTERED:
        wavelet_ = gsl_wavelet_alloc(gsl_wavelet_daubechies_centered, k);
        break;

      case WAVELET_HAAR:
        wavelet_ = gsl_wavelet_alloc(gsl_wavelet_haar, k);
        break;

      case WAVELET_HAAR_CENTERED:
        wavelet_ = gsl_wavelet_alloc(gsl_wavelet_haar_centered, k);
        break;

      case WAVELET_BSPLINE:
        wavelet_ = gsl_wavelet_alloc(gsl_wavelet_bspline, k);
        break;

      case WAVELET_BSPLINE_CENTERED:
        wavelet_ = gsl_wavelet_alloc(gsl_wavelet_bspline_centered, k);
        break;

      default:
        std::cout << "Error: Undefined wavelet is used in \"Wavelet2D\"." << std::endl;
        exit(1);
      }
    }
    void Free()
    {
      gsl_wavelet_free(wavelet_);
    }

  public:
    explicit Wavelet2D(int k = 2, WaveletType waveletType = WAVELET_HAAR_CENTERED)
    {
      Alloc(k, waveletType);
    }

    virtual ~Wavelet2D()
    {
      Free();
    }

    // virtual void Clear();
    itpp::mat Forward(const itpp::mat &input);
    itpp::mat Inverse(const itpp::mat &input);
    gsl_wavelet GslWavelet() const
    {
      return *wavelet_;
    }
  
  };
  
  /************************************************************************************
   * To1D_forWVT -- ウェーブレット変換用に左上から2の累乗サイズのブロックずつとっていく
   * 
   * Arguments:
   *   input -- 二次元配列
   ************************************************************************************/
  template <typename kind>
  inline std::vector< kind > To1D_forWVT(const Vector_2D< kind > &input)
  {
    int height = input.Height();
    int width = input.Width();
    assert(Radix2(height) && Radix2(width) && (height == width) && input.is_rectangular());
    
    if (height == 1){
      return input.To1D();
    } // if height
    else {
      DividePixels< kind > quarterPixels(input, height/2);
      std::vector< kind > output = To1D_forWVT(quarterPixels.GetBlock(0,0)); // 再帰関数

      Vector_2D< kind > temp = quarterPixels.GetBlock(0,1);
      output = Concat(output, temp.To1D());
    
      temp = quarterPixels.GetBlock(1,0);
      output = Concat(output, temp.To1D());

      temp = quarterPixels.GetBlock(1,1);
      output = Concat(output, temp.To1D());

      return output;
    } // else height
  }

  /************************************************************************************
   * To2D_forWVT -- To1D_forWVTによって1次元にされたものを2次元に直す
   * 
   * Arguments:
   *   input -- 1次元配列
   ************************************************************************************/
  template <typename kind>
  inline Vector_2D< kind > To2D_forWVT(const std::vector< kind > &input)
  {
    int sizeVec = input.size();
    if (sizeVec == 1){
      return To2D(input, 1);
    } // if size                               
    else{
      int side = static_cast< int >(sqrt(sizeVec));
      assert(Radix2(side));
      CompPixels< kind > output(side, side);
      Vector_2D< kind > temp2D = To2D_forWVT(Mid(input, 0, sizeVec/4)); // 再帰関数
      assert(temp2D.Width() == side/2);
      output(0, 0, temp2D);
      
      std::vector< kind > temp = Mid(input, sizeVec/4, sizeVec/4);
      output(0, side/2, To2D(temp, side/2));

      temp = Mid(input, 2*sizeVec/4, sizeVec/4);
      output(side/2, 0, To2D(temp, side/2));

      temp = Mid(input, 3*sizeVec/4, sizeVec/4);
      output(side/2, side/2, To2D(temp, side/2));

      return output.Get();
    } // else size
  }

  /************************************************************************************
   * NumAreas -- ウェーブレット変換された領域が何個のエリア分割数
   * ただし右下右上左下は同じエリア扱い
   * 
   * Arguments:
   *   waveletSize -- wavelet変換領域の1辺の長さ
   *
   * Return Value:
   *   int -- 分割数
   ************************************************************************************/
  inline int NumAreas(int waveletSize)
  {
    return log2(static_cast< double >(waveletSize)) + 1;
  }

  // 上で分けたエリアの番号の総ピクセル数
  inline int NumPix(int area)
  {
    int startPos = static_cast< int >(pow(pow(2, area-1),2));
    int numPix = static_cast< int >(pow(pow(2, area),2)) - startPos; // 正方形から左上1/4を抜いた部分

    return numPix;
  }

  
} // end of mylib

#endif
