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
// Last Updated: <2013/02/14 17:12:56 from Yoshitos-iMac.local by yoshito>
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

    struct {
      mutable u_char* curData_;
      mutable int     curBits_;
      mutable int     current_;
    } Ite_;
          
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

    bool IsZero() const;
        
    // 代入演算子
    BinVec operator = (const BinVec& input);
    BinVec operator = (const itpp::bvec& input);

    // 二項演算子
    friend BinVec operator +(const BinVec& v1, const BinVec& v2);
    friend BinVec operator +(const BinVec& v, const itpp::bin& b);
    friend BinVec operator +(const itpp::bin& b, const BinVec& v);
    
    friend itpp::bin operator *(const BinVec& v1, const BinVec& v2); // inner product
    friend BinVec operator *(const BinVec& v, const itpp::bin& b); // scalar product
    friend BinVec operator *(const itpp::bin& b, const BinVec& v);

    // イテレータ用関数++++++++++++++++++++++
    void      Begin() const;
    bool      End() const;
    itpp::bin Current() const;

    void Next() const;
    void Next(int n) const;
    // ------------------------------------
  };

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

    int            i = 0;
    for (Begin(); !End(); ++i, Next()){
      Set(i, oldBvec[i]);
    } // for ite
    
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

    int            i = 0;
    for (Begin(); !End(); ++i, Next()){
      Set(i, oldBvec[i]);
    } // for ite

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
    for (Begin(); !End(); Next(), ++i)
      {
        output[i] = Current();
      } // for ite
    
    return output;
  }

  inline BinVec BinVec::Left(int n) const
  {
    CheckRange(n);
    BinVec output(n);

    Begin();
    for (int i = 0; i < n; ++i, Next())
      {
        output.Set(i, Current());
      } // for i

    return output;
  }

  inline BinVec BinVec::Right(int n) const
  {
    CheckRange(n);
    BinVec output(n);
    Begin();
    int startPos = size_ - n;
    for (int i = 0; i < startPos; ++i){
      Next();
    } // for i
    
    for (int i = 0; i < n; ++i, Next())
      {
        assert(!End());
        output.Set(i, Current());
      } // for i

    return output;
  }

  inline BinVec BinVec::Mid(int start, int n) const
  {
    CheckRange(start);
    CheckRange(start + n - 1);
    
    BinVec output(n);
    Begin();
    for (int i = 0; i < start; ++i){
      Next();
    } // for i

    
    for (int i = 0; i < n; ++i, Next()){
      assert(!End());
      output.Set(i, Current());
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
    for (temp.Begin(); !temp.End(); temp.Next(), ++i){
      Set(i, temp.Current());
    } // for ite
      
    for (input.Begin(); !input.End(); input.Next(),++i){
      Set(i, input.Current());
    } // for ite
  }

  inline void BinVec::Append(const itpp::bvec& input)
  {
    int totalSize = size_ + input.size();
    BinVec temp(*this);

    Resize(totalSize);

    int i = 0;
    for (temp.Begin(); !temp.End(); temp.Next(), ++i){
      Set(i, temp.Current());
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
    for (temp.Begin(); !temp.End(); temp.Next(), ++i){
      Set(i, temp.Current());
    } // for ite
    Set(totalSize - 1, input);
        
  }

  inline void BinVec::Ins(int index, itpp::bin& input)
  {
    assert(index >= 0 && index <= size_);
    
    int totalSize = size_ + 1;
    BinVec temp(*this);

    Resize(totalSize);
    temp.Begin();
    for (int i = 0; i < index; ++i, temp.Next()){
      Set(i, temp.Current());
    } // for i

    Set(index, input);          // Insert here

    for (int i = index + 1; !temp.End(); ++i, temp.Next()){
      Set(i, temp.Current());
    } // for i
    
  }

  inline void BinVec::Ins(int index, BinVec& input)
  {
    assert(index >= 0 && index <= size_);

    int totalSize = size_ + input.size_;
    BinVec temp(*this);

    Resize(totalSize);
    temp.Begin();
    for (int i = 0; i < index; ++i, temp.Next()){
      Set(i, temp.Current());
    } // for i


    input.Begin();
    for (int i = index; !input.End(); ++i, input.Next()){
      Set(i, input.Current());
    } // for i

    for (int i = index + input.size_; !temp.End(); ++i, temp.Next()){
      Set(i, temp.Current());
    } // for i
  }

  inline void BinVec::Ins(int index, itpp::bvec& input)
  {
    assert(index >= 0 && index <= size_);

    int totalSize = size_ + input.size();
    BinVec temp(*this);

    Resize(totalSize);
    temp.Begin();
    for (int i = 0; i < index; ++i, temp.Next()){
      Set(i, temp.Current());
    } // for i

    // Insert here
    for (int i = 0; i < input.size(); ++i){
      Set(i + index, input[i]);
    } // for i

    for (int i = index + input.size(); !temp.End(); ++i, temp.Next()){
      Set(i, temp.Current());
    } // for i
  }

  inline void BinVec::Del(int index)
  {
    CheckRange(index);

    int totalSize = size_ - 1;
    BinVec temp(*this);
    
    Resize(totalSize);
    temp.Begin();
    for (int i = 0; i < index; ++i, temp.Next()){
      Set(i, temp.Current());
    } // for i

    temp.Next();                 // Delete here
    
    for (int i = index; !temp.End(); ++i, temp.Next()){
      Set(i, temp.Current());
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
    temp.Begin();
    for (int i = 0; i < i1; ++i, temp.Next()){
      Set(i, temp.Current());
    } // for i

    temp.Next(i2 - i1 + 1);      // Delete here

    for (int i = i1; !temp.End(); ++i, temp.Next()){
      Set(i, temp.Current());
    } // for i
  }

  inline bool BinVec::IsZero() const
  {
    for (int i = 0; i < (size_-1)/SIZE_UCHAR; ++i){
      if (data_[i] != 0x00){
        return false;
      }
    } // for i

    // 最後の配列の上位ビットは1の可能性があるが、size_を超えていたら関係ないので
    // 途中まで見る
    for (int i = 0; i < size_ % SIZE_UCHAR; ++i){
      if (((data_[(size_-1)/SIZE_UCHAR] >> i) & 1) != 0x00){
        return false;
      }
    } // for i
    
    return true;
  }
  
  inline void BinVec::Begin() const
  {
    Ite_.curData_ = &(data_[0]);
    Ite_.curBits_ = 0;
    Ite_.current_ = 0;
  }

  inline bool BinVec::End() const
  {
    if(Ite_.current_ >= size_){
      return true;
    }
    else{
      return false;
    }
  }

  inline itpp::bin BinVec::Current() const
  {
    return (*(Ite_.curData_) >> Ite_.curBits_) & 1;
  }

  inline void BinVec::Next() const
  {
    ++Ite_.current_;
    ++Ite_.curBits_;
    if(Ite_.curBits_ >= SIZE_UCHAR){
      Ite_.curBits_ = 0;
      ++Ite_.curData_;
    }
    
  }

  inline void BinVec::Next(int n) const
  {
    assert(n > 0);
    for (int i = 0; i < n; ++i){
      Next();
    }
  }

  inline BinVec operator +(const BinVec& v1, const BinVec& v2)
  {
    assert(v1.size_ == v2.size_);
    BinVec output(v1.size_);

    for (int i = 0; i < (v1.size_-1)/SIZE_UCHAR + 1; ++i){
      output.data_[i] = v1.data_[i] ^ v2.data_[i];
    } // for i
    
    return output;
  }

  inline BinVec operator +(const BinVec& v, const itpp::bin& b)
  {
    BinVec output(v.size_);

    if (b == 1){
      for (int i = 0; i < (v.size_-1)/SIZE_UCHAR + 1; ++i){
        output.data_[i] = v.data_[i] ^ 0xff;
      } // for i
    }
    else{
      output = v;
    }
    
    return output;
  }

  inline BinVec operator +(const itpp::bin& b, const BinVec& v)
  {
    BinVec output(v.size_);
    
    if (b == 1){
      for (int i = 0; i < (v.size_-1)/SIZE_UCHAR + 1; ++i){
        output.data_[i] = 0xff ^ v.data_[i];
      } // for i
    }
    else{
      output = v;
    }
   
    return output;
  }

  inline itpp::bin operator *(const BinVec& v1, const BinVec& v2)
  {
    assert(v1.size_ == v2.size_);
    itpp::bin output(0);

    v1.Begin();
    v2.Begin();
    for (int i = 0; i < v1.size_; ++i, v1.Next(), v2.Next()){
      output += v1.Current() & v2.Current();
    } // for i

    return output;
  }

  inline BinVec operator *(const BinVec& v, const itpp::bin& b)
  {
    BinVec output(v.size_);

    if (b == 1){
      for (int i = 0; i < (v.size_-1)/SIZE_UCHAR + 1; ++i){
        output.data_[i] = v.data_[i] & 0xff;
      } // for i
    }
    else{
      for (int i = 0; i < (v.size_-1)/SIZE_UCHAR + 1; ++i){
        output.data_[i] = v.data_[i] & 0x00;
      } // for i
    }

    return output;
  }

  inline BinVec operator *(const itpp::bin& b, const BinVec& v)
  {
    BinVec output(v.size_);

    if (b == 1){
      for (int i = 0; i < (v.size_-1)/SIZE_UCHAR + 1; ++i){
        output.data_[i] =  0xff & v.data_[i];
      } // for i
    }
    else{
      for (int i = 0; i < (v.size_-1)/SIZE_UCHAR + 1; ++i){
        output.data_[i] =  0x00 & v.data_[i];
      } // for i
    }

    return output;
  }

  inline itpp::bvec Bvecify(const BinVec& input)
  {
    itpp::bvec output(input.Size());
    
    int i = 0;
    for (input.Begin(); !input.End(); input.Next(), ++i)
      {
        output[i] = input.Current();
      } // for ite
    
    return output;

  }
  
} // end of mylib

inline std::ostream& operator << (std::ostream& outFile, const mylib::BinVec& binvec)
{

  outFile << "[ ";
  for (binvec.Begin();!binvec.End(); binvec.Next()){
    outFile << binvec.Current() << ' ';
  } // while ite
  outFile << ']';
  
  return outFile;
}



#endif
