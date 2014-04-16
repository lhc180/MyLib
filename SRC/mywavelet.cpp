#include "../include/mywavelet.h"

namespace mylib {

  itpp::mat Wavelet2D::Forward(const itpp::mat &input)
  {
    int height = input.rows();
    int width = input.cols();
    assert(mylib::Radix2(height) && mylib::Radix2(width));
  
    int workspaceSize = std::max(width, height);
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

    itpp::mat output(height, width);
    for (int row = 0; row < height; ++row){
      for (int col = 0; col < width; ++col){
        output(row, col) = gsl_matrix_get(inputGslMatrix, row, col);
      } // for col
    } // for row

    gsl_matrix_free(inputGslMatrix);
    gsl_wavelet_workspace_free(workspace);

    return output;
  }

  itpp::mat Wavelet2D::Inverse(const itpp::mat &input)
  {
    int height = input.rows();
    int width = input.cols();
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

    itpp::mat output(height, width);
    for (int row = 0; row < height; ++row){
      for (int col = 0; col < width; ++col){
        output(row, col) = gsl_matrix_get(inputGslMatrix, row, col);
      } // for col
    } // for row

    gsl_matrix_free(inputGslMatrix);
    gsl_wavelet_workspace_free(workspace);

    return output;
  }
  
  
} // mylib
