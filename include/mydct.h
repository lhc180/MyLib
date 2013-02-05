#ifndef MYDCT_H
#define MYDCT_H

#include <vector>
#include <itpp/itbase.h>
#include "myutl.h"
#include "mymatrix.h"

const double DIV_SQRT2 = 1.0/sqrt(2.0);

template<typename kind>
class cDct
{
protected:
  u_int size_;
  mylib::Vector_2D< double > cosT_;
  bool setDone_;
    
public:
  cDct(): setDone_(false)
  { 
  }
      
  cDct(u_int size)
  {
    setSize(size);
  }
  
  ~cDct()
  {
    clear();
  }
  
  void setSize(u_int size)
  {
    assert((size > 0) && ((size & (size-1)) == 0));
  
    size_ = size;
    cosT_.set_size(size, size);
    double pi_2N = M_PI / (2.0*size);
  
    for(u_int k = 0; k < size; k++)
      {
        for(u_int n = 0; n < size; n++)
          {
            cosT_(k, n) = cos((2*n + 1) * k * pi_2N);
          } // for n
      }   // for k
  
    setDone_ = true;
  }
  
  std::vector< kind > forward(const std::vector< kind > &input)
  {
    assert(setDone_ && (input.size() == size_));
  
    std::vector< kind > output(size_);
  
    double coeff = sqrt(2.0/size_);
  
    for(u_int k = 0; k < size_; k++)
      {
        double ck = (k == 0) ? DIV_SQRT2 : 1;
        double sum = 0;
      
        for(u_int n = 0; n < size_; n++)
          {
            sum += static_cast< double >(input[n]) * cosT_(k, n);
          } // for n
        output[k] = static_cast< kind >(coeff * ck * sum);
      } // for k

    return output;
  }

  itpp::Vec< kind > forward(const itpp::Vec< kind > &input)
  {
    assert(setDone_ && (input.size() == size_));
  
    itpp::Vec< kind > output(size_);
  
    double coeff = sqrt(2.0/size_);
  
    for(u_int k = 0; k < size_; k++)
      {
        double ck = (k == 0) ? DIV_SQRT2 : 1;
        double sum = 0;
      
        for(u_int n = 0; n < size_; n++)
          {
            sum += static_cast< double >(input[n]) * cosT_(k, n);
          } // for n
        output[k] = static_cast< kind >(coeff * ck * sum);
      } // for k
  
    return output;
  }

  mylib::Vector_2D< kind > forward(const mylib::Vector_2D< kind > &input) // 2D
  {
    assert(input.is_rectangular());
    assert(input.size_rows() == input.size_cols(0));
    assert(setDone_ && (static_cast<u_int>(input.size_rows()) == size_));

    mylib::Vector_2D< kind > output(size_, size_);
    
    double coeff = 2.0/size_;

    for(u_int k = 0; k < size_; k++)
      {
        double ck = (k == 0) ? DIV_SQRT2 : 1;
        for(u_int l = 0; l < size_; l++)
          {
            double cl = (l == 0) ? DIV_SQRT2 : 1;
            double sum = 0;          
            for(u_int m = 0; m < size_; m++)
              {
                for(u_int n = 0; n < size_; n++)
                  {
                    sum += static_cast< double >(input(m, n)) * cosT_(l, n) * cosT_(k, m);
                  } // for n
              }     // for m
            output(k, l) = static_cast< kind >(coeff * cl * ck * sum);
          } // for l
      }     // for k
    return output;
  }

  std::vector< kind > inverse(const std::vector< kind > &input)
  {
    assert(setDone_ && (input.size() == size_));
  
    std::vector< kind > output(size_);
  
    double coeff = sqrt(2.0/size_);
  
    for(u_int n = 0; n < size_; n++)
      {
        double sum = 0;
      
        for(u_int k = 0; k < size_; k++)
          {
            double ck = (k == 0) ? DIV_SQRT2 : 1;
            sum += ck * static_cast< double >(input[k]) * cosT_(k, n);
          } // for k
        output[n] = static_cast< kind >(coeff * sum);
      } // for n

    return output;
  }

  itpp::Vec< kind > inverse(const itpp::Vec< kind > &input)
  {
    assert(setDone_ && (input.size() == size_));
  
    itpp::Vec< kind > output(size_);
  
    double coeff = sqrt(2.0/size_);
  
    for(u_int n = 0; n < size_; n++)
      {
        double sum = 0;
      
        for(u_int k = 0; k < size_; k++)
          {
            double ck = (k == 0) ? DIV_SQRT2 : 1;
            sum += ck * static_cast< double >(input[k]) * cosT_(k, n);
          } // for k
        output[n] = static_cast< kind >(coeff * sum);
      } // for n

    return output;
  }

  mylib::Vector_2D< kind > inverse(const mylib::Vector_2D< kind > &input) // 2D
  {
    assert(input.is_rectangular());
    assert(input.size_rows() == input.size_cols(0));
    assert(setDone_ && (static_cast<u_int>(input.size_rows()) == size_));

    mylib::Vector_2D< kind > output(size_, size_);
  
    double coeff = 2.0/size_;

    for(u_int m = 0; m < size_; m++)
      {
        for(u_int n = 0; n < size_; n++)
          {
            double sum = 0;          
            for(u_int k = 0; k < size_; k++)
              {
                double ck = (k == 0) ? DIV_SQRT2 : 1;
                for(u_int l = 0; l < size_; l++)
                  {
                    double cl = (l == 0) ? DIV_SQRT2 : 1;
                    sum += ck * cl * static_cast< double >(input(k, l)) * cosT_(l, n) * cosT_(k, m);
                  } // for l
              }     // for k
            output(m, n) = static_cast< kind >(coeff * sum);
          } // for n
      }     // for m

    return output;
  }
  
  void clear()
  {
    cosT_.clear();
  }
};

#endif
