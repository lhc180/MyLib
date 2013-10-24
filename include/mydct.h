#ifndef MYDCT_H
#define MYDCT_H

#include <vector>
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include "myutl.h"
#include "mymatrix.h"
// #include "myjpeg.h"

namespace mylib{

  const double SQRT2 = 1.41421356;                    // 2の平方根
  const double DIS_SQRT2 = 1.0 / SQRT2;               // 2の平方根の逆数
  const double PI_DIV16 = 3.14159265 / 16;    // 円周率/16

  
  const double DIV_SQRT2 = 1.0/sqrt(2.0);

  template < typename kind >
  class Dct
  {
  protected:
    u_int size_;
    itpp::Mat< double > cosT_;
    bool setDone_;
    
  public:  
    explicit Dct(u_int size)
    {
      SetSize(size);
    }
  
    ~Dct()
    {
      Clear();
    }
  
    void SetSize(u_int size)
    {
      assert((size > 0) && ((size & (size-1)) == 0));
  
      size_ = size;
      cosT_.set_size(size, size);
      double pi_2N = M_PI / (2.0*static_cast< double >(size));
  
      for(u_int k = 0; k < size; k++)
        {
          for(u_int n = 0; n < size; n++)
            {
              cosT_(k, n) = cos((2*n + 1) * k * pi_2N);
            } // for n
        }   // for k
  
      setDone_ = true;
    }
  
    std::vector< kind > Forward(const std::vector< kind > &input, kind shift = 0)
    {
      assert(setDone_ && (input.size() == size_));
  
      std::vector< kind > output(size_);
  
      double coeff = sqrt(2.0/static_cast< double >(size_));
  
      for(u_int k = 0; k < size_; k++)
        {
          double ck = (k == 0) ? DIS_SQRT2 : 1;
          double sum = 0;
      
          for(u_int n = 0; n < size_; n++)
            {
              sum += static_cast< double >(input[n]) * cosT_(k, n);
            } // for n
          output[k] = static_cast< kind >(coeff * ck * sum + shift);
        } // for k

      return output;
    }
    
    mylib::Vector_2D< kind > Forward(const mylib::Vector_2D< kind > &input, kind shift = 0) // 2D
    {
      assert(input.is_rectangular());
      assert(input.size_rows() == input.size_cols(0));
      assert(setDone_ && (static_cast<u_int>(input.size_rows()) == size_));

      mylib::Vector_2D< kind > output(size_, size_);
    
      double coeff = 2.0/static_cast< double >(size_);

      for(u_int k = 0; k < size_; k++)
        {
          double ck = (k == 0) ? DIS_SQRT2 : 1;
          for(u_int l = 0; l < size_; l++)
            {
              double cl = (l == 0) ? DIS_SQRT2 : 1;
              double sum = 0;          
              for(u_int m = 0; m < size_; m++)
                {
                  for(u_int n = 0; n < size_; n++)
                    {
                      sum += static_cast< double >(input(m, n)) * cosT_(l, n) * cosT_(k, m);
                    } // for n
                }     // for m
              output(k, l) = static_cast< kind >(coeff * cl * ck * sum + shift);
            } // for l
        }     // for k
      return output;
    }

    itpp::Mat< kind > Forward(const itpp::Mat< kind > &input, kind shift = 0) // 2D
    {

      assert(input.rows() == input.cols());
      assert(setDone_ && (static_cast<u_int>(input.rows()) == size_));

      itpp::Mat< kind > output(size_, size_);
    
      double coeff = 2.0/static_cast< double >(size_);

      for(u_int k = 0; k < size_; k++)
        {
          double ck = (k == 0) ? DIS_SQRT2 : 1;
          for(u_int l = 0; l < size_; l++)
            {
              double cl = (l == 0) ? DIS_SQRT2 : 1;
              double sum = 0;          
              for(u_int m = 0; m < size_; m++)
                {
                  for(u_int n = 0; n < size_; n++)
                    {
                      sum += static_cast< double >(input(m, n)) * cosT_(l, n) * cosT_(k, m);
                    } // for n
                }     // for m
              output(k, l) = static_cast< kind >(coeff * cl * ck * sum + shift);
            } // for l
        }     // for k
      return output;
    }

    
    // shiftに値を入れると一緒にシフトも出来るアイデア
    std::vector< kind > Inverse(const std::vector< kind > &input, kind shift = 0)
    {
      assert(setDone_ && (input.size() == size_));
  
      std::vector< kind > output(size_);
  
      double coeff = sqrt(2.0/size_);
  
      for(u_int n = 0; n < size_; n++)
        {
          double sum = 0;
      
          for(u_int k = 0; k < size_; k++)
            {
              double ck = (k == 0) ? DIS_SQRT2 : 1;
              sum += ck * static_cast< double >(input[k]) * cosT_(k, n);
            } // for k
          output[n] = static_cast< kind >(coeff * sum);
        } // for n

      return output;
    }
    
    mylib::Vector_2D< kind > Inverse(const mylib::Vector_2D< kind > &input, kind shift = 0) // 2D
    {
      assert(input.is_rectangular());
      assert(input.size_rows() == input.size_cols(0));
      assert(setDone_ && (static_cast<u_int>(input.size_rows()) == size_));

      mylib::Vector_2D< kind > output(size_, size_);
  
      double coeff = 2.0/static_cast< double >(size_);
      
      for(u_int m = 0; m < size_; m++){
        for(u_int n = 0; n < size_; n++){
          double sum = 0;
          
          for(u_int k = 0; k < size_; k++){
            double ck = (k == 0) ? DIS_SQRT2 : 1.0; // ##
            for(u_int l = 0; l < size_; l++){
              double cl = (l == 0) ? DIS_SQRT2 : 1.0; // ##
              sum += ck * cl * static_cast< double >(input(k, l)) * cosT_(l, n) * cosT_(k, m);
            } // for l
          }     // for k
          output(m, n) = static_cast< kind >(coeff * sum + shift); 
        } // for n
      }     // for m

      return output;
    }

    
    itpp::Mat< kind > Inverse(const itpp::Mat< kind > &input, kind shift = 0) // 2D
    {
      assert(input.rows() == input.cols());
      assert(setDone_ && (static_cast<u_int>(input.rows()) == size_));
      
      itpp::Mat< kind > output(size_, size_);
  
      double coeff = 2.0/static_cast< double >(size_);
      
      for(u_int m = 0; m < size_; m++){
        for(u_int n = 0; n < size_; n++){
          double sum = 0;
          
          for(u_int k = 0; k < size_; k++){
            double ck = (k == 0) ? DIS_SQRT2 : 1.0; // ##
            for(u_int l = 0; l < size_; l++){
              double cl = (l == 0) ? DIS_SQRT2 : 1.0; // ##
              sum += ck * cl * static_cast< double >(input(k, l)) * cosT_(l, n) * cosT_(k, m);
            } // for l
          }     // for k
          output(m, n) = static_cast< kind >(coeff * sum + shift); 
        } // for n
      }     // for m

      return output;
    }
    
    void Clear()
    {
      cosT_.clear();
    }
  };

  // itpp::dctを使った2次元DCT
  inline itpp::mat Dct2D(const itpp::mat& input)
  {
    itpp::mat output1(input.rows(), input.cols());
    
    for (int i = 0; i < input.rows(); ++i){
      itpp::vec temp = input.get_row(i);
      output1.set_row(i,itpp::dct(temp));
    } // for i

    itpp::mat output2(input.rows(), input.cols());

    for (int i = 0; i < input.cols(); ++i){
      itpp::vec temp = output1.get_col(i);
      output2.set_col(i, itpp::dct(temp));
    } // for i

    return output2;
  }

  inline itpp::mat Idct2D(const itpp::mat& input)
  {
    itpp::mat output1(input.rows(), input.cols());
    
    for (int i = 0; i < input.rows(); ++i){
      itpp::vec temp = input.get_row(i);
      output1.set_row(i,itpp::idct(temp));
    } // for i

    itpp::mat output2(input.rows(), input.cols());

    for (int i = 0; i < input.cols(); ++i){
      itpp::vec temp = output1.get_col(i);
      output2.set_col(i, itpp::idct(temp));
    } // for i

    return output2;
  }

  inline mylib::vec_2D Dct2D(const mylib::vec_2D& input)
  {
    itpp::mat input_itpp = ToItppMat(input);
    itpp::mat output_itpp = Dct2D(input_itpp);
    mylib::vec_2D output = ToVector_2D(output_itpp);
    return output;
  }

  inline mylib::vec_2D Idct2D(const mylib::vec_2D& input)
  {
    itpp::mat input_itpp = ToItppMat(input);
    itpp::mat output_itpp = Idct2D(input_itpp);
    mylib::vec_2D output = ToVector_2D(output_itpp);
    return output;
  }
  
} // end of namespace mylib

#endif
