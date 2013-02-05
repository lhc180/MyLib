#ifndef MYMATRIX_H
#define MYMATRIX_H

/***********************************/
// std::vectorによる二次元vectorを使いやすいようにクラス化する
// 各行の列数は可変
/**********************************/

#include <vector>
#include <cassert>
#include <complex>
#include "myutl.h"
// #include <itpp/itbase.h>


namespace mylib
{

  template<typename kind>
  class Vector_2D{
  protected:
    std::vector< std::vector<kind> > vMat_; // 行番目のnElem列目の要素値

    void check_rowRange(int row) const{
      assert((row >= 0) && (row < static_cast<int>(vMat_.size())));
    }

    void check_colRange(int row, int col) const{
      assert((col >= 0) && (col < static_cast<int>(vMat_[row].size())));
    }
    
  public:
    // constructor
    Vector_2D(int rows = 0, int cols = 0): vMat_(rows, std::vector< kind>(cols,0))
    { }
    
    // destructor
    ~Vector_2D(){
      clear();
    }

    // copy constructor is default

    void set_size(int rows = 0, int cols = 0){
      assert(rows >= 0 && cols >= 0);
      vMat_.resize(rows);
    
      for(int i = 0; i < rows; i++){
	vMat_[i].resize(cols);
      }
    }
  
    void resize(int rows = 0, int cols = 0){
      set_size(rows, cols);
    }
   
    // 単独の行の列数を設定
    void set_colSize(int row, int cols = 0){
      check_rowRange(row);
      assert(cols >= 0);
      vMat_[row].resize(cols);
    }

    // 全ての行の列数が等しいかどうか
    bool is_rectangular() const{
      
      if(vMat_.size() == 0){
	return true;
      }
      else{
	bool result = true;
	int defaultCols = vMat_[0].size();
	for(int row = 1; row < static_cast<int>(vMat_.size()); row++){
	  if(defaultCols != static_cast<int>(vMat_[row].size())){
	    result = false;
	    break;
	  }
	} // for row
	return result;
      }
    }

    void clear(){
      vMat_.clear();
    }

    void zeros(){
      for(int row = 0; row < vMat_.size(); row++){
	for(int col = 0; col < vMat_[row].size(); col++){
	  vMat_[row][col] = 0;
	} // for col
      }	  // for row
    }

    void ones(){
      for(int row = 0; row < vMat_.size(); row++){
	for(int col = 0; col < vMat_[row].size(); col++){
	  vMat_[row][col] = 1;
	} // for col
      }	  // for row
    } 
    
    int size_rows() const{
      return vMat_.size();
    }

    int Height() const{
      return vMat_.size();
    }
    
    // 指定した行番号の列数を返す
    int size_cols(int row = 0) const{
      check_rowRange(row);
      return vMat_[row].size();
    }

    int Width(int row = 0) const{
      check_rowRange(row);
      return vMat_[row].size();
    }
    
    // row行col列目の要素へのアクセス
    const kind &operator()(int row, int col) const{
      check_rowRange(row);
      check_colRange(row, col);
      return vMat_[row][col];
    }
    kind &operator()(int row, int col){
      check_rowRange(row);
      check_colRange(row, col);
      return vMat_[row][col];
    }

    // row行のvectorへのアクセス
    const std::vector<kind> &operator()(int row) const{
      check_rowRange(row);
      return vMat_[row];
    }
    std::vector<kind> &operator()(int row){
      check_rowRange(row);
      return vMat_[row];
    }

    // num*colsの行列を追加する
    void add_rows(int cols = 0, int num = 1){
      assert(num >= 1);
      for(int i = 0; i < num; i++){
	vMat_.push_back(std::vector<kind>(cols));
      }
    }
  
    void add_row(const std::vector<kind> &input)
    {
      vMat_.push_back(input);
    }
        
    // ある行に列を追加する
    void add_cols(int row, kind elem = 0, int num = 1){
      check_rowRange(row);
      assert(num >= 1);
      for(int i = 0; i < num; i++){
	vMat_[row].push_back(elem);
      }
    }
    
    // ある行にvectorをセットする
    // void set_row(int row, const std::vector< kind > &input){
    //   check_rowRange(row);
    //   vMat_[row] = input;
    // }

    
    // 最後の行に要素を追加する
    // void push_back(kind elem){
    //   int lastRow = vMat_.size() - 1;
    //   assert(lastRow >= 0);

    //   vMat_[lastRow].push_back(elem);
    // }

    // push_backはmat(i).push_back(a)でいける

    // Add element to row.
    void push_back(int row, kind elem)
    {
      check_rowRange(row);
      vMat_[row].push_back(elem);
    }

    // Get sub-Vector_2D form rows r1 to r2 and columns c1 to c2.
    // -1 as r2 and c2 indecates the last row and column, respectively.
    Vector_2D< kind > get(int r1, int r2, int c1, int c2) const
    {
      is_rectangular();
      check_rowRange(r1);
      check_colRange(r1, c1);
          
      int startRow = r1;
      int lastRow  = (r2 == -1) ? (size_rows()-1) : r2;
      int numRows  = lastRow - startRow + 1;
      int startCol = c1;
      int lastCol  = (c2 == -1) ? (size_cols(lastRow)-1) : c2;
      int numCols  = lastCol - startCol + 1;
      
      Vector_2D< kind > output(numRows);
      for(int r = 0; r < numRows; r++)
        {
          output(r) = getMid(vMat_[startRow + r], startCol, numCols);
        } // for r

      return output;
      
    }

    // transform to std::vector.
    // The order is (0,0) -> (0,size_cols(0)-1) -> (1, 0) ~~~ (size_rows()-1, size_cols()-1)
    std::vector< kind > to_1D() const
    {
      std::vector< kind > output(0);
      
      for(int row = 0; row < vMat_.size(); row++)
        {
          for(int col = 0; col < vMat_[row].size(); col++)
            {
              output.push_back(vMat_[row][col]);
            }
        }
      
      return output;
    }

  };

  
  typedef Vector_2D<double> vec_2D;
  typedef Vector_2D<int> ivec_2D;
  typedef Vector_2D< std::complex<double> > cvec_2D;

  
  template< typename kind >
  inline Vector_2D< kind > to_2D(const std::vector< kind > &input, int numCols)
  {
    int numRows = input.size() / numCols + 
      ((input.size() % numCols == 0) ? 0 : 1);
    
    Vector_2D< kind > output(numRows);
    for(int i = 0; i < input.size(); i++)
      {
        int row = i / numCols;
        output(row).push_back(input[i]);
      }
    
    return output;
  }

  template< typename kind >
  inline Vector_2D< double > toDouble(const Vector_2D< kind > &input)
  {
    Vector_2D< double > output(input.size_rows(), 0);

    for(int row = 0; row < input.size_rows(); row++)
      {
        for(int col = 0; col < input.size_cols(row); col++)
          {
            output(row).push_back(static_cast< double >(input(row, col)));
          } // for col
      }     // for row
    
    return output;
  }

  template< typename kind >
  inline Vector_2D< u_char > toU_char(const Vector_2D< kind > &input)
  {
    Vector_2D< u_char > output(input.size_rows(), 0);

    for(int row = 0; row < input.size_rows(); row++)
      {
        for(int col = 0; col < input.size_cols(row); col++)
          {
            output(row).push_back(static_cast< u_char >(input(row, col)));
          } // for col
      }     // for row
    
    return output;
  }

  inline vec_2D Abs(const vec_2D &input)
  {
    vec_2D output(input.size_rows(), 0);

    for (int row = 0; row < input.size_rows(); ++row){
      output(row) = Abs(input(row));
    } // for row

    return output;
  }

  template <typename kind>
  inline kind Max(const Vector_2D< kind > &input)
  {
    kind max = input(0,0);

    for (int i = 0; i < input.size_rows(); ++i){
      int temp = Max(input(i));
      if(max < temp){
        max = temp;
      }
    } // for i
    return max;
  }

  template <typename kind>
  inline kind Min(const Vector_2D< kind > &input)
  {
    kind min = input(0,0);

    for (int i = 0; i < input.size_rows(); ++i){
    int temp = Min(input(i));
      if(min > temp){
        min = temp;
      }
    } // for curInput
    return min;
  }
  
} // end of namespace mylib



template <typename kind>
inline std::ostream& operator << (std::ostream& outFile, const mylib::Vector_2D< kind >& vec2D)
{
  for (int row = 0; row < vec2D.size_rows(); ++row)
    {
      outFile << "[ ";
      for (int col = 0; col < vec2D.size_cols(row); ++col)
        {
          outFile << vec2D(row, col) << ' ';
        } // for col
      outFile << "]" << std::endl;
    } // for row

  return outFile;
}

#endif
