#ifndef MYWAVELET_H
#define MYWAVELET_H

#include <iostream>
#include <cstdlib>
#include <gsl/gsl_wavelet2d.h>
#include "myutl.h"
#include "mymatrix.h"

// 2次元ウェーブレット変換のクラス
// 1次元は作らない

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
  WaveletType waveletType_;
  gsl_wavelet *wavelet_;
  

  void Alloc();
  void Free();

public:
  explicit Wavelet2D(WaveletType waveletType = WAVELET_HAAR_CENTERED): waveletType_(waveletType)
  {
    Alloc();
  }

  virtual ~Wavelet2D()
  {
    Free();
  }

  // virtual void Clear();

  virtual mylib::vec_2D Forward(const mylib::vec_2D &input);
  virtual mylib::vec_2D Inverse(const mylib::vec_2D &input);
};

inline void Wavelet2D::Alloc()
{


  switch (waveletType_){
  case WAVELET_DAUBECHIES:
    wavelet_ = gsl_wavelet_alloc(gsl_wavelet_daubechies, 4);
    break;

  case WAVELET_DAUBECHIES_CENTERED:
    wavelet_ = gsl_wavelet_alloc(gsl_wavelet_daubechies_centered, 4);
    break;

  case WAVELET_HAAR:
    wavelet_ = gsl_wavelet_alloc(gsl_wavelet_haar, 2);
    break;

  case WAVELET_HAAR_CENTERED:
    wavelet_ = gsl_wavelet_alloc(gsl_wavelet_haar_centered, 2);
    break;

  case WAVELET_BSPLINE:
    wavelet_ = gsl_wavelet_alloc(gsl_wavelet_bspline, 103);
    break;

  case WAVELET_BSPLINE_CENTERED:
    wavelet_ = gsl_wavelet_alloc(gsl_wavelet_bspline_centered, 103);
    break;

  default:
    std::cout << "Error: Undefined wavelet is used in \"Wavelet2D\"." << std::endl;
    exit(1);
  }


}

inline void Wavelet2D::Free()
{
  gsl_wavelet_free(wavelet_);
}

inline mylib::vec_2D Wavelet2D::Forward(const mylib::vec_2D &input)
{
  int height = input.size_rows();
  int width = input.size_cols();
  assert(mylib::Radix2(height) && mylib::Radix2(width));
  
  int workspaceSize = mylib::maxOf(width, height);
  gsl_wavelet_workspace *workspace = gsl_wavelet_workspace_alloc(workspaceSize);

  // gsl_matrixにinputを格納
  gsl_matrix *inputGslMatrix = gsl_matrix_alloc(height, width);
  for (int row = 0; row < height; ++row){
    for (int col = 0; col < width; ++col){
      gsl_matrix_set(inputGslMatrix, row, col, input(row,col));
    } // for col
  } // for row

  int status = gsl_wavelet2d_nstransform_matrix_forward(wavelet_, inputGslMatrix, workspace);
  if (status != GSL_SUCCESS){
    std::cout << "Error in Wavelete2D::Forward." << std::endl;
    std::cout << "Height and width of matrix are not equal or insufficient workspace is provided." << std::endl;
    exit(1);
  }

  mylib::vec_2D output(height, width);
  for (int row = 0; row < height; ++row){
    for (int col = 0; col < width; ++col){
      output(row, col) = gsl_matrix_get(inputGslMatrix, row, col);
    } // for col
  } // for row

  gsl_matrix_free(inputGslMatrix);
  gsl_wavelet_workspace_free(workspace);

  return output;
}

inline mylib::vec_2D Wavelet2D::Inverse(const mylib::vec_2D &input)
{
  int height = input.size_rows();
  int width = input.size_cols();
  assert(mylib::Radix2(height) && mylib::Radix2(width));
  
  int workspaceSize = mylib::maxOf(width, height);
  gsl_wavelet_workspace *workspace = gsl_wavelet_workspace_alloc(workspaceSize);

  // gsl_matrixにinputを格納
  gsl_matrix *inputGslMatrix = gsl_matrix_alloc(height, width);
  for (int row = 0; row < height; ++row){
    for (int col = 0; col < width; ++col){
      gsl_matrix_set(inputGslMatrix, row, col, input(row,col));
    } // for col
  } // for row

  int status = gsl_wavelet2d_nstransform_matrix_inverse(wavelet_, inputGslMatrix, workspace);
  if (status != GSL_SUCCESS){
    std::cout << "Error in Wavelete2D::Forward." << std::endl;
    std::cout << "Height and width of matrix are not equal or insufficient workspace is provided." << std::endl;
    exit(1);
  }

  mylib::vec_2D output(height, width);
  for (int row = 0; row < height; ++row){
    for (int col = 0; col < width; ++col){
      output(row, col) = gsl_matrix_get(inputGslMatrix, row, col);
    } // for col
  } // for row

  gsl_matrix_free(inputGslMatrix);
  gsl_wavelet_workspace_free(workspace);

  return output;
}


#endif
