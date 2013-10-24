#ifndef MYBINMAT_H
#define MYBINMAT_H

#include "mybinvec.h"

namespace mylib
{
  class BinMat
  {
  private:
    BinVec data_;
    int rows_;
    int cols_;
    int totalSize_;

    void CheckRange(int row, int col) const;
    void CheckRowRange(int row) const;
    void CheckColRange(int col) const;
    void Init(int rows, int cols); // 値の代入のみ    
    void Alloc(int rows, int cols); // 値の代入とdata_の確保
    void Free();

  public:
    // Constructor
    BinMat(): data_(0), rows_(0), cols_(0), totalSize_(0)
    { }
    explicit BinMat(int rows, int cols = 0, const itpp::bin &initValue = 0)
    {
      Alloc(rows, cols);

      if (initValue == itpp::bin(0)){
        Zeros();
      }
      else{
        Ones();
      }
    }
    
    virtual ~BinMat()
    {
      Free();
    }

    BinMat(const BinMat& oldBinMat);
    BinMat(const itpp::bmat& oldBmat);
        
    int Rows() const { return rows_; }
    int Height() const { return rows_; }

    int Cols() const { return cols_; }
    int Width() const { return cols_; }

    void SetSize(int rows, int cols){
      Free();
      Alloc(rows, cols);
    }
    void Resize(int rows, int cols){
      Free();
      Alloc(rows, cols);
    }

    void Zeros() { data_.Zeros(); }
    void Ones() { data_.Ones(); }

    void Delete(){ Free(); }

    // Accessor
    // -- Getter
    itpp::bin Get(int row, int col) const;
    
    // -- Setter
    void Set(int row, int col, const itpp::bin& input);

    itpp::bmat Bmatify() const;

    BinVec GetRow(int r) const;
    BinMat GetRows(int r1, int r2) const;
    BinVec GetCol(int c) const;
    BinMat GetCols(int c1, int c2) const;

    void SetRow(int r, const BinVec& input);
    void SetRow(int r, const itpp::bvec& input);
    void SetCol(int c, const BinVec& input);
    void SetCol(int c, const itpp::bvec& input);
    
    void SwapRows(int r1, int r2);
    void SwapCols(int c1, int c2);
    
    // 下側にくっつける系
    void AppendRow(const BinVec& input);
    void AppendRow(const itpp::bvec& input);
    void AppendRows(const BinMat& input);
    BinMat ConcatVertic(const BinMat& input) const;
        
    // 右側にくっつける系
    void AppendCol(const BinVec& input);
    void AppendCol(const itpp::bvec& input);
    void AppendCols(const BinMat& input);
    BinMat ConcatHorizon(const BinMat& input) const;
        
    void InsRow(int r, const BinVec& input);
    void InsRow(int r, const itpp::bvec& input);
    void InsCol(int c, const BinVec& input);
    void InsCol(int c, const itpp::bvec& input);
    void DelRow(int r);
    void DelRows(int r1, int r2);
    void DelCol(int c);
    void DelCols(int c1, int c2);

    // r1 += r2
    void AddRows(int r1, int r2);

    void PermuteRows(const itpp::ivec& perm, bool inverse = false);
    void PermuteRows(const std::vector< int >& perm, bool inverse = false);
    void PermuteCols(const itpp::ivec& perm, bool inverse = false);
    void PermuteCols(const std::vector< int >& perm, bool inverse = false);
        
    BinMat Get(int r1, int r2, int c1, int c2) const;

    BinMat Transpose() const;

    bool IsZero() const;

    // 代入演算子
    BinMat operator = (const BinMat& input);
    BinMat operator = (const itpp::bmat& input);

    friend BinMat operator +(const BinMat& m1, const BinMat& m2);
    friend BinMat operator +(const BinMat& m, const itpp::bin& b);
    friend BinMat operator +(const itpp::bin& b, const BinMat& m);

    friend BinMat operator *(const BinMat& m1, const BinMat& m2);
    friend BinMat operator *(const BinMat& m, const itpp::bin& b);
    friend BinMat operator *(const itpp::bin& b, const BinMat& m);
    friend BinVec operator *(const BinVec& v, const BinMat& m); // v * m (output row vector)
    friend BinVec operator *(const BinMat& m, const BinVec& v); // m * v (output column vector)
    
  };

  inline BinMat::BinMat(const BinMat& oldBinMat)
  {
    Init(oldBinMat.rows_, oldBinMat.cols_);

    data_ = oldBinMat.data_;
  }

  inline BinMat BinMat::operator = (const BinMat& oldBinMat)
  {
    Init(oldBinMat.rows_, oldBinMat.cols_);
    
    data_ = oldBinMat.data_;

    return (*this);
  }
  
  
  inline BinMat::BinMat(const itpp::bmat& oldBmat)
  {
    Init(oldBmat.rows(), oldBmat.cols());

    itpp::bvec temp(0);
    for (int i = 0; i < oldBmat.rows(); ++i){
      temp = itpp::concat(temp, oldBmat.get_row(i));
    }

    data_ = temp;
  }

  inline BinMat BinMat::operator = (const itpp::bmat& oldBmat)
  {
    Init(oldBmat.rows(), oldBmat.cols());

    itpp::bvec temp(0);
    for (int i = 0; i < oldBmat.rows(); ++i){
      temp = itpp::concat(temp, oldBmat.get_row(i));
    }

    data_ = temp;

    return (*this);
  }
  
  inline void BinMat::CheckRange(int row, int col) const
  {
    assert(row >= 0 && row < rows_);
    assert(col >= 0 && col < cols_);
  }

  inline void BinMat::CheckRowRange(int row) const
  {
    assert(row >= 0 && row < rows_);
  }

  inline void BinMat::CheckColRange(int col) const
  {
    assert(col >= 0 && col < cols_);
  }

  inline void BinMat::Init(int rows, int cols)
  {
    assert(rows >= 0 && cols >= 0);

    rows_ = rows;
    cols_ = cols;
    totalSize_ = rows * cols;
  }
  
  inline void BinMat::Alloc(int rows, int cols)
  {
    Init(rows, cols);
    
    data_.SetSize(totalSize_);
  }

  inline void BinMat::Free()
  {
    data_.Delete();
    rows_ = 0;
    cols_ = 0;
    totalSize_ = 0;
  }

  inline itpp::bin BinMat::Get(int row, int col) const
  {
    CheckRange(row, col);

    return data_.Get(row*cols_ + col);
  }

  inline void BinMat::Set(int row, int col, const itpp::bin& input)
  {
    CheckRange(row, col);

    data_.Set(row*cols_ + col, input);
  }

  inline itpp::bmat BinMat::Bmatify() const
  {
    itpp::bmat output(rows_, cols_);

    data_.Begin();
    for (int r = 0; r < rows_; ++r){
      for (int c = 0; c < cols_; ++c){
        output(r, c) = data_.Current();
        data_.Next();
      } // for c
    } // for r

    return output;
  }
  
  inline BinVec BinMat::GetRow(int r) const
  {
    assert(r >= 0 && r < rows_);

    BinVec output = data_.Mid(r*cols_, cols_);
    
    return output;
  }

  inline BinMat BinMat::GetRows(int r1, int r2) const
  {
    if(r2 == -1){
      r2 = rows_ - 1;
    }
    CheckRowRange(r1);
    CheckRowRange(r2);
    assert(r2 > r1);

    int numRows = r2 - r1 + 1;

    BinMat output(numRows, cols_);
    output.data_ = data_.Mid(r1*cols_, numRows*cols_);
    
    return output;
  }

  inline BinVec BinMat::GetCol(int c) const
  {
    assert(c >= 0 && c < cols_);

    BinVec output(rows_);

    for (int r = 0; r < rows_; ++r){
      output.Set(r, data_.Get(r*cols_ + c));
    } // for r

    return output;
  }

  inline BinMat BinMat::GetCols(int c1, int c2) const
  {
    if(c2 == -1){
      c2 = cols_ - 1;
    }
    CheckColRange(c1);
    CheckColRange(c2);
    assert(c2 > c1);

    int numCols = c2 - c1 + 1;

    BinMat output(rows_, numCols);
    for (int r = 0; r < rows_; ++r){
      for (int c = 0; c < numCols; ++c){
        output.data_.Set(r*numCols +c, data_.Get(r*cols_ + c1 + c));
      } // for c
    } // for r

    return output;
  }

  inline void BinMat::SetRow(int r, const mylib::BinVec &input)
  {
    CheckRowRange(r);
    assert(cols_ == input.Size());
    
    input.Begin();
    for (int c = 0; c < cols_; ++c, input.Next()){
      data_.Set(r*cols_ + c, input.Current());
    } // for c
  }

  inline void BinMat::SetRow(int r, const itpp::bvec &input)
  {
    CheckRowRange(r);
    assert(cols_ == input.size());

    for (int c = 0; c < cols_; ++c){
      data_.Set(r*cols_ + c, input[c]);
    } // for c
  }

  inline void BinMat::SetCol(int c, const mylib::BinVec &input)
  {
    CheckColRange(c);
    assert(rows_ == input.Size());

    input.Begin();
    for (int r = 0; r < rows_; ++r, input.Next()){
      data_.Set(r*cols_ + c, input.Current());
    } // for r
  }

  inline void BinMat::SetCol(int c, const itpp::bvec &input)
  {
    CheckColRange(c);
    assert(rows_ == input.size());

    for (int r = 0; r < rows_; ++r){
      data_.Set(r*cols_ + c, input[r]);
    } // for r
  }
    
  inline void BinMat::SwapRows(int r1, int r2)
  {
    CheckRowRange(r1);
    CheckRowRange(r2);

    BinVec temp = GetRow(r1);
    SetRow(r1,GetRow(r2));
    SetRow(r2, temp);
  }

  inline void BinMat::SwapCols(int c1, int c2)
  {
    CheckColRange(c1);
    CheckColRange(c2);

    BinVec temp = GetCol(c1);
    SetCol(c1, GetCol(c2));
    SetCol(c2, temp);
  }

  inline void BinMat::AppendRow(const mylib::BinVec &input)
  {
    assert(cols_ == input.Size());

    BinMat temp = *this;
    Resize(rows_+1, cols_);

    for (int r = 0; r < temp.rows_; ++r){
      SetRow(r, temp.GetRow(r));
    } // for r

    SetRow(rows_-1, input);
  }

  inline void BinMat::AppendRow(const itpp::bvec &input)
  {
    assert(cols_ == input.size());

    BinMat temp = *this;
    Resize(rows_+1, cols_);

    for (int r = 0; r < temp.rows_; ++r){
      SetRow(r, temp.GetRow(r));
    } // for r

    SetRow(rows_-1, input);
  }
  
  inline void BinMat::AppendCol(const mylib::BinVec &input)
  {
    assert(rows_ == input.Size());
    
    BinMat temp = *this;
    Resize(rows_, cols_ + 1);

    for (int c = 0; c < temp.cols_; ++c){
      SetCol(c, temp.GetCol(c));
    } // for c

    SetCol(cols_-1, input);
  }

  inline void BinMat::AppendCol(const itpp::bvec &input)
  {
    assert(rows_ == input.size());

    BinMat temp = *this;
    Resize(rows_, cols_ + 1);

    for (int c = 0; c < temp.cols_; ++c){
      SetCol(c, temp.GetCol(c));
    } // for c

    SetCol(cols_-1, input);
  }

  inline void BinMat::InsRow(int r, const mylib::BinVec &input)
  {
    assert(r >= 0 && r <= rows_);
    
    BinMat temp(*this);
    Resize(rows_+1, cols_);
    for (int i = 0; i < r; ++i){
      SetRow(i, temp.GetRow(i));
    } // for i

    SetRow(r, input);           // Insert here

    for (int i = r + 1; i < rows_; ++i){
      SetRow(i , temp.GetRow(i-1));
    } // for i
    
  }

  inline void BinMat::InsRow(int r, const itpp::bvec& input)
  {
    assert(r >= 0 && r <= rows_);
    
    BinMat temp(*this);
    Resize(rows_+1, cols_);
    for (int i = 0; i < r; ++i){
      SetRow(i, temp.GetRow(i));
    } // for i

    SetRow(r, input);           // Insert here

    for (int i = r + 1; i < rows_; ++i){
      SetRow(i , temp.GetRow(i-1));
    } // for i
    
  }  

  inline void BinMat::InsCol(int c, const mylib::BinVec &input)
  {
    assert(c >= 0 && c <= cols_);

    BinMat temp(*this);
    Resize(rows_, cols_+1);
    for (int i = 0; i < c; ++i){
      SetCol(i, temp.GetCol(i));
    } // for i

    SetCol(c, input);

    for (int i = c + 1; c < cols_; ++i){
      SetCol(i, temp.GetCol(i-1));
    } // for i
  }

  inline void BinMat::InsCol(int c, const itpp::bvec &input)
  {
    assert(c >= 0 && c <= cols_);

    BinMat temp(*this);
    Resize(rows_, cols_+1);
    for (int i = 0; i < c; ++i){
      SetCol(i, temp.GetCol(i));
    } // for i

    SetCol(c, input);           // Insert here

    for (int i = c + 1; c < cols_; ++i){
      SetCol(i, temp.GetCol(i-1));
    } // for i
  }

  inline void BinMat::DelRow(int r)
  {
    CheckRowRange(r);

    BinMat temp(*this);
    Resize(rows_-1, cols_);

    for (int i = 0; i < r; ++i){
      SetRow(i, temp.GetRow(i));
    } // for i

    for (int i = r; i < rows_; ++i){
      SetRow(i, temp.GetRow(i+1));
    } // for i
    
  }

  inline void BinMat::DelRows(int r1, int r2)
  {
    if (r2 == -1){
      r2 = rows_-1;
    }
    assert(r1 < r2);
    CheckRowRange(r1);
    CheckRowRange(r2);

    int totalRows = rows_ - (r2 - r1 + 1);
    BinMat temp(*this);
    Resize(totalRows, cols_);

    for (int i = 0; i < r1; ++i){
      SetRow(i, temp.GetRow(i));
    } // for i

    for (int i = r1; i < rows_; ++i){
      SetRow(i, temp.GetRow(i+r2+1));
    } // for i
  }

  inline void BinMat::DelCol(int c)
  {
    CheckColRange(c);
    
    BinMat temp(*this);
    Resize(rows_, cols_-1);

    for (int i = 0; i < c; ++i){
      SetCol(i, temp.GetCol(i));
    } // for i

    for (int i = c; i < cols_; ++i){
      SetCol(i, temp.GetCol(i+1));
    } // for i
  }

  inline void BinMat::DelCols(int c1, int c2)
  {
    if (c2 == -1){
      c2 = cols_-1;
    }
    assert(c1 < c2);
    CheckColRange(c1);
    CheckColRange(c2);

    int totalCols = cols_ - (c2 - c1 +1);
    BinMat temp(*this);
    Resize(rows_, totalCols);

    for (int i = 0; i < c1; ++i){
      SetCol(i, temp.GetCol(i));
    } // for i

    for (int i = c1; i < cols_; ++i){
      SetCol(i, temp.GetCol(i+c2+1));
    } // for i
  }

  inline void BinMat::AddRows(int r1, int r2)
  {
    CheckRowRange(r1);
    CheckRowRange(r2);

    BinVec row1 = GetRow(r1);
    BinVec row2 = GetRow(r2);

    SetRow(r1, row1 + row2);
  }

  inline void BinMat::PermuteRows(const itpp::ivec& perm, bool inverse)
  {
    assert(rows_ == perm.size());

    BinMat temp(*this);

    if (inverse){
      for (int i = 0; i < rows_; ++i){
        SetRow(perm[i], temp.GetRow(i));
      } // for i
    }
    else{
      for (int i = 0; i < rows_; ++i){
        SetRow(i, temp.GetRow(perm[i]));
      } // for i
    }
  }

  inline void BinMat::PermuteRows(const std::vector< int >& perm, bool inverse)
  {
    assert(rows_ == static_cast< int >(perm.size()));

    BinMat temp(*this);

    if (inverse){
      for (int i = 0; i < rows_; ++i){
        SetRow(perm[i], temp.GetRow(i));
      } // for i
    }
    else{
      for (int i = 0; i < rows_; ++i){
        SetRow(i, temp.GetRow(perm[i]));
      } // for i
    }
  }

  
  inline void BinMat::PermuteCols(const itpp::ivec& perm, bool inverse)
  {
    assert(cols_ == perm.size());

    BinMat temp(*this);

    if (inverse){
      for (int i = 0; i < cols_; ++i){
        SetCol(perm[i], temp.GetCol(i));
      } // for i
    }
    else{
      for (int i = 0; i < cols_; ++i){
        SetCol(i, temp.GetCol(perm[i]));
      } // for i
    }
  }

  inline void BinMat::PermuteCols(const std::vector< int >& perm, bool inverse)
  {
    assert(cols_ == static_cast< int >(perm.size()));

    BinMat temp(*this);

    if (inverse){
      for (int i = 0; i < cols_; ++i){
        SetCol(perm[i], temp.GetCol(i));
      } // for i
    }
    else{
      for (int i = 0; i < cols_; ++i){
        SetCol(i, temp.GetCol(perm[i]));
      } // for i
    }
  }

  
  inline BinMat BinMat::Get(int r1, int r2, int c1, int c2) const
  {
    if (r2 == -1){
      r2 = rows_ - 1;
    }
    if (c2 == -1){
      c2 = cols_ - 1;
    }
    assert(r1 < r2);
    assert(c1 < c2);
    CheckRange(r1, c1);
    CheckRange(r2, c2);

    int numRows = r2 -r1 + 1;
    int numCols = c2 -c1 + 1;

    BinMat output(numRows, numCols);
    for (int r = 0; r < numRows; ++r){
      for (int c = 0; c < numCols; ++c){ 
        output.data_.Set(r*numCols + c, data_.Get((r+r1)*cols_ + c1 + c));
      }                         // for c
    }                           // for r

    return output;
  }

  inline BinMat BinMat::Transpose() const
  {
    BinMat output(cols_, rows_);

    for (int i = 0; i < rows_; ++i){
      BinVec temp = GetRow(i);
      output.SetCol(i, temp);
    } // for i
    return output;
  }

  inline bool BinMat::IsZero() const
  {
    return data_.IsZero();
  }
  
  inline void BinMat::AppendCols(const BinMat& input)
  {
    assert(input.rows_ == rows_);

    for (int i = 0; i < input.cols_; ++i){
      AppendCol(input.GetCol(i));
    } // for i
  }

  inline void BinMat::AppendRows(const BinMat& input)
  {
    assert(input.cols_ == cols_);

    for (int i = 0; i < input.rows_; ++i){
      AppendRow(input.GetRow(i));
    } // for i
  }

  inline BinMat BinMat::ConcatVertic(const BinMat& input) const
  {
    assert(input.cols_ == cols_);

    BinMat output(rows_ + input.rows_, cols_);

    for (int i = 0; i < rows_; ++i){
      output.SetRow(i, GetRow(i));
    } // for i
    for (int i = 0; i < input.rows_; ++i){
      output.SetRow(i+rows_, input.GetRow(i));
    } // for i

    return output;
  }

  inline BinMat BinMat::ConcatHorizon(const BinMat& input) const
  {
    assert(input.rows_ == rows_);

    BinMat output(rows_, cols_ + input.cols_);

    for (int i = 0; i < cols_; ++i){
      output.SetCol(i, GetRow(i));
    } // for i
    for (int i = 0; i < input.cols_; ++i){
      output.SetCol(i+cols_, input.GetCol(i));
    } // for i

    return output;
  }
  
  // 単位行列を作る
  inline BinMat DenseEye(int size)
  {
    assert(size > 0);

    BinMat output(size, size);
    output.Zeros();
    for (int i = 0; i < size; ++i){
      output.Set(i, i, 1);
    } // for i

    return output;
  }

  inline BinMat operator +(const BinMat& m1, const BinMat& m2)
  {
    assert(m1.rows_ == m2.rows_ && m1.cols_ == m2.cols_);

    BinMat output;
    output.Init(m1.rows_, m1.cols_);
    output.data_ = m1.data_ + m2.data_;
    
    return output;
  }

  inline BinMat operator +(const BinMat& m, const itpp::bin& b)
  {
    BinMat output;
    output.Init(m.rows_, m.cols_);
    output.data_ = m.data_ + b;

    return output;
  }

  inline BinMat operator +(const itpp::bin& b, const BinMat& m)
  {
    BinMat output;
    output.Init(m.rows_, m.cols_);
    output.data_ = b + m.data_;

    return output;
  }

  inline BinMat operator *(const BinMat& m1, const BinMat& m2)
  {
    assert(m1.cols_ == m2.rows_);
    BinMat output(m1.rows_, m2.cols_);

    for (int r = 0; r < m1.rows_; ++r){
      for (int c = 0; c < m2.cols_; ++c){
        output.Set(r, c, m1.GetRow(r) * m2.GetCol(c));
      } // for c
    } // for r

    return output;
  }

  inline BinMat operator *(const BinMat& m, const itpp::bin& b)
  {
    BinMat output;
    output.Init(m.rows_, m.cols_);
    output.data_ = m.data_ * b;
    
    return output;
  }

  inline BinMat operator *(const itpp::bin& b, const BinMat& m)
  {
    BinMat output;
    output.Init(m.rows_, m.cols_);
    output.data_ = b * m.data_;

    return output;
  }

  inline BinVec operator *(const BinVec& v, const BinMat& m)
  {
    assert(v.Size() == m.rows_);
    BinVec output(m.cols_);

    for (int i = 0; i < m.cols_; ++i){
      output.Set(i, v * m.GetCol(i));
    } // for i
    return output;
  }

  inline BinVec operator *(const BinMat& m, const BinVec& v)
  {
    assert(m.cols_ == v.Size());
    BinVec output(m.rows_);

    for (int i = 0; i < m.rows_; ++i){
      output.Set(i, m.GetRow(i) * v);
    } // for i
    return output;
  }
  
} // end of namespace mylib

inline std::ostream& operator << (std::ostream& outFile, const mylib::BinMat& binmat)
{
  for (int r = 0; r < binmat.Rows(); ++r){
    outFile << "[ ";
    for (int c = 0; c < binmat.Cols(); ++c){
      outFile << binmat.Get(r, c) << ' ';
    } // for c
    outFile << ']' << std::endl;
  } // for r

  return outFile;
}

#endif


