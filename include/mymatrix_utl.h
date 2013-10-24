#ifndef MYMATRIX_UTL_H
#define MYMATRIX_UTL_H

/*************************************/
// Vector_2Dとitppの互換性を持たせるための関数群
/************************************/

#include "./mymatrix.h"
#include "./mybinmat.h"
#include <itpp/itbase.h>
#include <cassert>

namespace mylib
{
  // 1である要素番号のみを保持
  inline Vector_2D<int> compressGF2matToVector_2D(const itpp::GF2mat &inputMat)
  {
    Vector_2D<int> mat(inputMat.rows());

    for(int row = 0; row < inputMat.rows(); row++){
      for(int col = 0; col < inputMat.cols(); col++){
	if(inputMat(row,col) == 1){
	  mat(row).push_back(col);
	}
      }
    }

    return mat;
  }

  inline ivec_2D CompressBinMat(const BinMat& input)
  {
    ivec_2D output(input.Rows());

    for (int r = 0; r < input.Rows(); ++r){
      for (int c = 0; c < input.Cols(); ++c){
        if (input.Get(r, c) == 1){
          output(r).push_back(c);
        }
      } // for c
    } // for r

    return output;
  }

  // itpp::bvecと1である要素番号が格納されたmylib::Vector_2D<int>の掛け算
  // ただし、Vector_2Dの方は転置した形で格納されていなければならない
  // つまり、bvec * Mat の場合は、bvec * Mat_transposeを入れる
  // 同様に、bvec * Mat_transposeは bvec * Matとする
  inline itpp::bvec times(const itpp::bvec &in_bvec, const Vector_2D<int> &in_mat)
  {
    itpp::bvec out(in_mat.size_rows());
    
    for(int row = 0; row < out.size(); row++){
      itpp::bin sum = 0;
      for(int col = 0; col < in_mat.size_cols(row); col++){
	int nElem = in_mat(row,col); // 1である要素番号
	// assert(nElem < in_bvec.size()); // ##
	sum += in_bvec[nElem] * itpp::bin(1);
      }	// for col
      out[row] = sum;
    } // for row

    return out;
  }

  inline BinVec Times(const BinVec& v, const ivec_2D& m)
  {
    BinVec output(m.size_rows());

    for (int r = 0; r < m.size_rows(); ++r){
      itpp::bin sum = 0;
      for (int c = 0; c < m.size_cols(r); ++c){
        int index = m(r, c);    // 1である要素番号
        // assert(index < v.Size());
        sum += v.Get(index) * itpp::bin(1);
      } // for c
      output.Set(r, sum);
    } // for r

    return output;
  }

} // end of namespace mylib
#endif
