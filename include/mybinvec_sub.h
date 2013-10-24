// オリジナルのバイナリベクトルのクラス
// u_charの各ビットを1要素として持つ
// itpp::binやbvecが使える環境を前提
// 
// もしオリジナルバイナリクラスを作ればtemplateを使ってitppバージョンか
// オリジナルバージョンか切り替えできるようにする。
//
// イテレータも実装
// 走査をするならイテレータ使うと早くなる
//
// ゲッタとセッタは明示的に行う
// つまり、operator[]は実装しない
// 
// Last Updated: <2013/02/12 21:20:47 from Yoshitos-iMac.local by yoshito>
//
// ---------- ToDo List -----------
// -- イテレータの実装
// -- 各種メンバ関数
// -- オペレータ


#ifndef MYBINVEC_H
#define MYBINVEC_H

#include <itpp/itbase.h>
#include <cstdlib>
#include <cmath>
#include <myutl.h>

static const int SIZE_UCHAR = 8*sizeof(u_char);

namespace mylib
{
  class BinVec
  {
  private:
    int     size_;
    u_char *data_;

    void CheckRange(int size) const;
    void Alloc(int size);
    void Free();
        
  public:
    // Constructor
    BinVec():size_(0), data_(0)
    {
      assert(SIZE_UCHAR == 8);
    }
    explicit BinVec(int size, const itpp::bin &initValue = 0)
    {
      assert(SIZE_UCHAR == 8);

      Alloc(size);

      if(initValue == itpp::bin(0)){
        Zeros();
      }
      else{
        Ones();
      }
            
    }
    
    BinVec(const BinVec &oldBinVec); // Copy Constructor
    BinVec(const itpp::bvec &oldBvec);

    ~BinVec()                   // Destructor
    {
      Free();
    }
    
    int Length() const{
      return size_;
    }      
    int Size() const{
      return size_;
    }

    void SetSize(int size){
      Free();
      Alloc(size);
    }
    void Resize(int size){
      Free();
      Alloc(size);
    }
    
    void Zeros();               // 要素を全て0にする
                                // deleteするわけではない

    void Ones();                // 要素を全て1にする

    void Delete()
    {
      Free();
    }
    
    // Accessor
    // itpp::bin operator[](int i) const;
    // itpp::bin &operator[](int i); --- セッタは明示的に行わなければいけない
    // const itpp::bin &at(int i) const;
    // itpp::bin &at(int i );
    itpp::bin Get(int i) const;
    void Set(int i, const itpp::bin& input);

    itpp::bvec Bvecify() const;
    BinVec Right(int n) const;
    BinVec Left(int n) const;
    BinVec Mid(int start, int n) const;
    BinVec Get(int i1, int i2) const;

    void Swap(int i1, int i2);
    void Append(const BinVec& input);
    void Append(const itpp::bvec& input);
    void Append(const itpp::bin& input);
    void Ins(int i, itpp::bin& input);
    void Ins(int i, BinVec& input);
    void Ins(int i, itpp::bvec& input);
    void Del(int i);
    void Del(int i1, int i2);
    
    // 代入演算子
    BinVec operator = (const BinVec& input);
    BinVec operator = (const itpp::bvec& input);

    friend class BinVecIterator;
  };

  // 反復子
  class BinVecIterator
  {
  private:
    u_char* curData_;
    int     curBits_;
    int     current_;
    int     size_;

  public:
    explicit BinVecIterator(const BinVec& binVec)
    {
      Begin(binVec);
    }
    virtual ~BinVecIterator()
    {
      curData_ = 0;
    }
    void      Begin(const BinVec& binVec);
    bool      End();
    itpp::bin Current();

    void Next();
    void Next(int n);
            
  };
  
  // inline BinVec& operator ++(BinVec& obj);
  
  
  inline void BinVec::CheckRange(int size) const
  {
    assert(size >= 0);
    assert(size < size_);
  }
  
  inline void BinVec::Alloc(int size)
  {
    assert(size >= 0);

    size_ = size;
    data_ = new u_char[(size-1)/SIZE_UCHAR + 1];
  }

  inline void BinVec::Free()
  {
    delete[] data_;
    data_ = 0;
    size_ = 0;
    
  }

  // Copy Constructor
  inline BinVec::BinVec(const BinVec& oldBinVec)
  {
    Alloc(oldBinVec.size_);

    memcpy(data_, oldBinVec.data_, ((size_-1)/SIZE_UCHAR + 1)*sizeof(u_char));
  }

  inline BinVec::BinVec(const itpp::bvec &oldBvec)
  {
    Alloc(oldBvec.size());

    BinVecIterator ite(*this);
    int            i = 0;
    while (!ite.End()){
      Set(i, oldBvec[i]);
      i++;
      ite.Next();
    }                           // while ite
    
  }

  inline BinVec BinVec::operator = (const BinVec& oldBinVec)
  {
    Free();
    Alloc(oldBinVec.size_);

    memcpy(data_, oldBinVec.data_, ((size_-1)/SIZE_UCHAR + 1)*sizeof(u_char));
        
    return (*this);
  }

  inline BinVec BinVec::operator = (const itpp::bvec& oldBvec)
  {
    Free();
    Alloc(oldBvec.size());

    BinVecIterator ite(*this);
    int            i = 0;
    while (!ite.End()){
      Set(i, oldBvec[i]);
      i++;
      ite.Next();
    }                           // while ite

    return (*this);
  }
  
  inline void BinVec::Zeros()
  {
    for (int i = 0; i < (size_-1)/SIZE_UCHAR + 1; ++i)
      {
        data_[i] = 0x00;
      }                         // for i
  }

  inline void BinVec::Ones()
  {
    for (int i = 0; i < (size_-1)/SIZE_UCHAR + 1; ++i)
      {
        data_[i] = 0xff;
      }                         // for i
  }
  
  inline itpp::bin BinVec::Get(int i) const
  {
    CheckRange(i);
    return (data_[i/SIZE_UCHAR] >> (i%SIZE_UCHAR)) & 0x1;
  }

  inline void BinVec::Set(int i, const itpp::bin &input)
  {
    CheckRange(i);
    if(input == itpp::bin(1)){
      data_[i/SIZE_UCHAR] |= (0x1 << (i%SIZE_UCHAR));
    }
    else{
      data_[i/SIZE_UCHAR] &= (0xff ^ (1 << (i%SIZE_UCHAR)));
    }
  }

  inline itpp::bvec BinVec::Bvecify() const
  {
    itpp::bvec output(size_);

    int i = 0;
    for (BinVecIterator ite(*this); !ite.End(); ite.Next(), ++i)
      {
        output[i] = ite.Current();
      } // for ite
    
    return output;
  }

  inline BinVec BinVec::Left(int n) const
  {
    CheckRange(n);
    BinVec output(n);

    BinVecIterator ite(*this);
    for (int i = 0; i < n; ++i, ite.Next())
      {
        output.Set(i, ite.Current());
      } // for i

    return output;
  }

  inline BinVec BinVec::Right(int n) const
  {
    CheckRange(n);
    BinVec output(n);
    BinVecIterator ite(*this);
    int startPos = size_ - n;
    for (int i = 0; i < startPos; ++i){
      ite.Next();
    } // for i
    
    for (int i = 0; i < n; ++i, ite.Next())
      {
        assert(!ite.End());
        output.Set(i, ite.Current());
      } // for i

    return output;
  }

  inline BinVec BinVec::Mid(int start, int n) const
  {
    CheckRange(start);
    CheckRange(start + n - 1);
    
    BinVec output(n);
    BinVecIterator ite(*this);
    for (int i = 0; i < start; ++i){
      ite.Next();
    } // for i

    
    for (int i = 0; i < n; ++i, ite.Next()){
      assert(!ite.End());
      output.Set(i, ite.Current());
    } // for i
    
    return output;
  }

  inline BinVec BinVec::Get(int i1, int i2) const
  {
    if (i2 == -1){
      i2 = size_ - 1;
    }
    CheckRange(i1);
    CheckRange(i2);
    assert(i1 < i2);

    BinVec output(i2 - i1 + 1);
        
    for (int i = i1, k = 0; i <= i2; ++i, ++k){
      output.Set(k , Get(i));
    } // for i

    return output;
  }

  inline void BinVec::Swap(int i1, int i2)
  {
    CheckRange(i1);
    CheckRange(i2);

    itpp::bin temp = (*this).Get(i1);
    Set(i1, (*this).Get(i2));
    Set(i2, temp);
  }

  inline void BinVec::Append(const mylib::BinVec& input)
  {
    int totalSize = size_ + input.size_;
    BinVec temp(*this);

    Resize(totalSize);

    int i = 0;
    for (BinVecIterator ite(temp); !ite.End(); ite.Next(), ++i){
      Set(i, ite.Current());
    } // for ite
      
    for (BinVecIterator ite(input); !ite.End(); ite.Next(),++i){
      Set(i, ite.Current());
    } // for ite
  }

  inline void BinVec::Append(const itpp::bvec& input)
  {
    int totalSize = size_ + input.size();
    BinVec temp(*this);

    Resize(totalSize);

    int i = 0;
    for (BinVecIterator ite(temp); !ite.End(); ite.Next(), ++i){
      Set(i, ite.Current());
    } // for ite
      
    for (int k = 0; k < input.size(); ++k, ++i){
      Set(i, input[k]);
    } // for ite
  }

  inline void BinVec::Append(const itpp::bin& input)
  {
    int totalSize = size_ + 1;
    BinVec temp(*this);

    Resize(totalSize);

    int i = 0;
    for (BinVecIterator ite(temp); !ite.End(); ite.Next(), ++i){
      Set(i, ite.Current());
    } // for ite
    Set(totalSize - 1, input);
        
  }

  inline void BinVec::Ins(int index, itpp::bin& input)
  {
    assert(index >= 0 && index <= size_);
    
    int totalSize = size_ + 1;
    BinVec temp(*this);

    Resize(totalSize);
    BinVecIterator ite(temp);
    for (int i = 0; i < index; ++i, ite.Next()){
      Set(i, ite.Current());
    } // for i

    Set(index, input);          // Insert here

    for (int i = index + 1; !ite.End(); ++i, ite.Next()){
      Set(i, ite.Current());
    } // for i
    
  }

  inline void BinVec::Ins(int index, BinVec& input)
  {
    assert(index >= 0 && index <= size_);

    int totalSize = size_ + input.size_;
    BinVec temp(*this);

    Resize(totalSize);
    BinVecIterator ite(temp);
    for (int i = 0; i < index; ++i, ite.Next()){
      Set(i, ite.Current());
    } // for i


    BinVecIterator iteObj(input);
    for (int i = index; !iteObj.End(); ++i, iteObj.Next()){
      Set(i, iteObj.Current());
    } // for i

    for (int i = index + input.size_; !ite.End(); ++i, ite.Next()){

      Set(i, ite.Current());
    } // for i
  }

  inline void BinVec::Ins(int index, itpp::bvec& input)
  {
    assert(index >= 0 && index <= size_);

    int totalSize = size_ + input.size();
    BinVec temp(*this);

    Resize(totalSize);
    BinVecIterator ite(temp);
    for (int i = 0; i < index; ++i, ite.Next()){
      Set(i, ite.Current());
    } // for i

    // Insert here
    for (int i = 0; i < input.size(); ++i){
      Set(i + index, input[i]);
    } // for i

    for (int i = index + input.size(); !ite.End(); ++i, ite.Next()){
      Set(i, ite.Current());
    } // for i
  }

  inline void BinVec::Del(int index)
  {
    CheckRange(index);

    int totalSize = size_ - 1;
    BinVec temp(*this);
    
    Resize(totalSize);
    BinVecIterator ite(temp);
    for (int i = 0; i < index; ++i, ite.Next()){
      Set(i, ite.Current());
    } // for i

    ite.Next();                 // Delete here
    
    for (int i = index; !ite.End(); ++i, ite.Next()){
      Set(i, ite.Current());
    } // for i
  }

  inline void BinVec::Del(int i1, int i2)
  {
    if(i2 == -1){
      i2 = size_ - 1;
    }
    CheckRange(i1);
    CheckRange(i2);
    assert(i1 < i2);
    
    int totalSize = size_ - (i2 - i1 + 1);
    BinVec temp(*this);

    Resize(totalSize);
    BinVecIterator ite(temp);
    for (int i = 0; i < i1; ++i, ite.Next()){
      Set(i, ite.Current());
    } // for i

    ite.Next(i2 - i1 + 1);      // Delete here

    for (int i = i1; !ite.End(); ++i, ite.Next()){
      Set(i, ite.Current());
    } // for i
  }
  
  inline void BinVecIterator::Begin(const mylib::BinVec &binVec)
  {
    curData_ = &(binVec.data_[0]);
    curBits_ = 0;
    current_ = 0;
    size_    = binVec.size_;
  }

  inline bool BinVecIterator::End()
  {
    if(current_ >= size_){
      return true;
    }
    else{
      return false;
    }
  }

  inline itpp::bin BinVecIterator::Current()
  {
    return ((*curData_) >> curBits_) & 1;
  }

  inline void BinVecIterator::Next()
  {
    ++current_;
    ++curBits_;
    if(curBits_ >= SIZE_UCHAR){
      curBits_ = 0;
      ++curData_;
    }
    
  }

  inline void BinVecIterator::Next(int n)
  {
    assert(n > 0);
    for (int i = 0; i < n; ++i){
      Next();
    }
  }
  
} // end of mylib

inline std::ostream& operator << (std::ostream& outFile, const mylib::BinVec& binvec)
{
  mylib::BinVecIterator ite(binvec);
  outFile << "[ ";
  while (!ite.End()){
    outFile << ite.Current() << ' ';
    ite.Next();
  } // while ite
  outFile << ']';
  
  return outFile;
}

#endif
