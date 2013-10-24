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
    Vector_2D(int rows = 0, int cols = 0, kind val = 0): vMat_(rows, std::vector< kind>(cols,val))
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
      // check_rowRange(row);
      // check_colRange(row, col);
      return vMat_[row][col];
    }
    kind &operator()(int row, int col){
      // check_rowRange(row);
      // check_colRange(row, col);
      return vMat_[row][col];
    }

    // row行のvectorへのアクセス
    const std::vector<kind> &operator()(int row) const{
      // check_rowRange(row);
      return vMat_[row];
    }
    std::vector<kind> &operator()(int row){
      // check_rowRange(row);
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
      // check_rowRange(row);
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
    std::vector< kind > To1D() const
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
  inline Vector_2D< kind > To2D(const std::vector< kind > &input, int numCols)
  {
    int numRows = input.size() / numCols + 
      ((input.size() % numCols == 0) ? 0 : 1);
    
    Vector_2D< kind > output(numRows, numCols);

    for(int i = 0, size = input.size(); i < size; i++){
      int row = i / numCols;
      int col = i % numCols;
      output(row, col) = input[i];
    } // for i
    
    return output;
  }

  template< typename kind >
  inline Vector_2D< int > ToInt(const Vector_2D< kind > &input)
  {
    Vector_2D< int > output(input.size_rows(), 0);

    for(int row = 0; row < input.size_rows(); row++)
      {
        output(row).resize(input.size_cols(row));
        for(int col = 0; col < input.size_cols(row); col++)
          {
            output(row, col) = static_cast< int >(input(row, col));
          } // for col
      }     // for row
    
    return output;
  }

  
  template< typename kind >
  inline Vector_2D< double > ToDouble(const Vector_2D< kind > &input)
  {
    Vector_2D< double > output(input.size_rows(), 0);

    for(int row = 0; row < input.size_rows(); row++)
      {
        output(row).resize(input.size_cols(row));
        for(int col = 0; col < input.size_cols(row); col++)
          {
            output(row, col) = static_cast< double >(input(row, col));
          } // for col
      }     // for row
    
    return output;
  }

  template< typename kind >
  inline Vector_2D< u_char > ToU_char(const Vector_2D< kind > &input)
  {
    Vector_2D< u_char > output(input.size_rows(), 0);

    // std::cout << "## called" << std::endl;
    
    for(int row = 0; row < input.size_rows(); row++)
      {
        output(row).resize(input.size_cols(row));
        for(int col = 0; col < input.size_cols(row); col++)
          {
            kind temp = input(row, col);
            if (temp > 255){
              temp = 255;
            } // if
            if(temp < 0){
              temp = 0;
            }
            output(row, col) = static_cast< u_char >(temp);
          } // for col
      }     // for row
    
    return output;
  }

  template < typename kind >
  inline itpp::Mat< kind > ToItppMat(const mylib::Vector_2D< kind >& input)
  {
    assert(input.is_rectangular());
    itpp::Mat< kind > output(input.Height(), input.Width());
    
    for (int r = 0, rows = input.Height(); r < rows; ++r){
      std::vector< kind > temp = input(r);
      output.set_row(r, ToItppVec(temp));
    } // for r

    return output;
  }

  template < typename kind >
  inline mylib::Vector_2D< kind > ToVector_2D(const itpp::Mat< kind >& input)
  {
    mylib::Vector_2D< kind > output(input.rows(), input.cols());

    for (int r = 0, rows = input.rows(); r < rows ; ++r){
      itpp::Vec< kind > temp = input.get_row(r);
      output(r) = ToStdVector(temp);
    } // for r
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

  template < typename kind >
  inline Vector_2D< kind > LevelShift(const Vector_2D< kind > &input, kind level)
  {
    Vector_2D< kind > output(input.size_rows(), 0);
  
    for(int row = 0; row < input.size_rows(); row++)
      {
        output(row) = LevelShift(input(row), level);
      }
  
    return output;
  }

  template<typename kind>
  inline Vector_2D< kind > Divide(const Vector_2D< kind > &input,
                                  const Vector_2D< kind > &qTable)
  {
    assert(input.size_rows() == qTable.size_rows());

    Vector_2D< kind > output(input.size_rows());
    for(int row = 0; row < input.size_rows(); row++)
      {
        output(row) = Divide(input(row), qTable(row));
      }

    return output;
  }

  template<typename kind>
  inline Vector_2D< kind > Multiply(const Vector_2D< kind > &input,
                                    const Vector_2D< kind > &qTable)
  {
    assert(input.size_rows() == qTable.size_rows());

    Vector_2D< kind > output(input.size_rows());
    for(int row = 0; row < input.size_rows(); row++)
      {
        output(row) = Multiply(input(row), qTable(row));
      }

    return output;
  }

  // Divides BMP pixel data into square areas.
  template< typename kind >
  class DividePixels
  {
  protected:
    Vector_2D< kind > pixels_;
    int numPix_;
    mutable int totalHsize_;
    mutable int totalVsize_;
    mutable int numHBlocks_;
    mutable int numVBlocks_;
  
  public:
    DividePixels(const Vector_2D< kind >  &pixels, 
                 int numPix): pixels_(pixels), numPix_(numPix)
    {
      assert(pixels_.is_rectangular());

      totalHsize_ = pixels_.size_cols(0);
      totalVsize_ = pixels_.size_rows();
        
      numHBlocks_ = totalHsize_ / numPix_;
      if(totalHsize_ % numPix_ != 0)
        {
          numHBlocks_++;
        }
    
      numVBlocks_ = totalVsize_ / numPix_;
      if(totalVsize_ % numPix_ != 0)
        {
          numVBlocks_++;
        }
    }
  
    Vector_2D< kind > GetBlock(int Vblock, int Hblock) const
    {
      assert(Vblock >= 0 && Vblock < numVBlocks_);
      assert(Hblock >= 0 && Hblock < numHBlocks_);
    
      int startRow = Vblock * numPix_;
      int startCol = Hblock * numPix_;
    
      Vector_2D< kind > output(numPix_, numPix_);
      for(int v = 0; v < numPix_; v++)
        {
          int row = startRow + v;
          if(startRow + v >= totalHsize_)
            {
              row = totalHsize_ - 1;
            }
          for(int h = 0; h < numPix_; h++)
            {
              int col = startCol + h;
              if(startCol + h >= totalHsize_)
                {
                  col = totalVsize_ - 1;
                }
              output(v, h) = pixels_(row, col);
            }
        }
        
      return output;
    }

    Vector_2D< kind > operator()(int vBlock, int hBlock) const
    {
      return GetBlock(vBlock, hBlock);
    }
    
    int Hblocks() const
    {
      return numHBlocks_;
    }
    int Vblocks() const
    {
      return numVBlocks_;
    }
  };

  // Compose Vector_2D divided into blocks.
  template< typename kind >
  class CompPixels
  {
  protected:
    Vector_2D< kind > pixels_;
    int totalHsize_;
    int totalVsize_;
  
  public:
    CompPixels(int totalVsize, int totalHsize):
      pixels_(totalVsize, totalHsize), totalHsize_(totalHsize), totalVsize_(totalVsize)
    {
    }
  
    void operator()(int startRow, int startCol, const Vector_2D< kind > &input)
    {
      assert(startRow >= 0 && startRow < totalVsize_);
      assert(startCol >= 0 && startCol < totalHsize_);
      assert(input.is_rectangular());
    
      for(int v = 0; v < input.size_rows(); v++)
        {
          int row = startRow + v;
          if(row >= totalVsize_)
            {
              break;
            }
          for(int h = 0; h < input.size_cols(0); h++)
            {
              int col = startCol + h;
              if(col >= totalHsize_)
                {
                  break;
                }
              pixels_(row, col) = input(v, h);
            }
        }
    }
  
    Vector_2D< kind > Get()
    {
      return pixels_;
    }
  };

  
  
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
