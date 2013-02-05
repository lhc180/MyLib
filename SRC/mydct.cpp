#include <cmath>
#include "../include/mydct.h"

const double DIV_SQRT2 = 1.0/sqrt(2.0);

template< typename kind >
void cDct< kind >::setSize(u_int size)
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

template< typename kind >
std::vector< kind > cDct< kind >::forward(const std::vector< kind > &input)
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

template< typename kind >
itpp::Vec< kind > cDct< kind >::forward(const itpp::Vec< kind > &input)
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

// 2D
template< typename kind >
mylib::Vector_2D< kind > cDct< kind >::forward(const mylib::Vector_2D< kind > &input)
{
  assert(input.is_rectangular());
  assert(input.size_rows() == input.size_cols(0));
  assert(setDone_ && (input.size_rows() == size_));

  mylib::Vector_2D< kind > output(size_, size_);
    
  double coeff = 2.0/size_;

  for(u_int j = 0; j < size_; j++)
    {
      double cj = (j == 0) ? DIV_SQRT2 : 1;
      for(u_int k = 0; k < size_; k++)
        {
          double ck = (k == 0) ? DIV_SQRT2 : 1;
          double sum = 0;          
          for(u_int m = 0; m < size_; m++)
            {
              for(u_int n = 0; n < size_; n++)
                {
                  sum += input(m, n) * cosT_(k, n) * cosT_(j, m);
                } // for n
            }     // for m
          output(j, k) = static_cast< kind >(coeff * ck * cj * sum);
        } // for k
    }     // for j
  return output;
}


template< typename kind >
std::vector < kind > cDct< kind >::inverse(const std::vector< kind > &input)
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
  
template< typename kind >
itpp::Vec < kind > cDct< kind >::inverse(const itpp::Vec< kind > &input)
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
      
// 2D
template< typename kind >
mylib::Vector_2D< kind > cDct< kind >::inverse(const mylib::Vector_2D< kind > &input)
{
  assert(input.is_rectangular());
  assert(input.size_rows() == input.size_cols(0));
  assert(setDone_ && (input.size_rows() == size_));

  mylib::Vector_2D< kind > output(size_, size_);
  
  double coeff = 2.0/size_;

  for(u_int m = 0; m < size_; m++)
    {
      for(u_int n = 0; n < size_; n++)
        {
          double sum = 0;          
          for(u_int j = 0; j < size_; j++)
            {
              double cj = (j == 0) ? DIV_SQRT2 : 1;
              for(u_int k = 0; k < size_; k++)
                {
                  double ck = (k == 0) ? DIV_SQRT2 : 1;
                  sum += cj * ck * input(m, n) * cosT_(k, n) * cosT_(j, m);
                } // for k
            }     // for j
          output(m, n) = static_cast< kind >(coeff * sum);
        } // for n
    }     // for m

  return output;
}
