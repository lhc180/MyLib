#ifndef MYIMG_H
#define MYIMG_H

#include <vector>
#include <cmath>
#include "mymatrix.h"
#include "myutl.h"
#include "mybmp.h"
#include "mytiff.h"
#include "mydct.h"
#include "myjpeg.h"
#include "mywavelet.h"

static const int IMG_DCT_SIZE      = 8;
static const int IMG_DCT_SIZE2     = 64;
static const int ZIGZAG_TABLE_SIZE = 64;

// for zigzag sequence.
static const int ZIGZAG_TABLE[64] = {
  0,  1,  8, 16,  9,  2,  3, 10,
  17, 24, 32, 25, 18, 11,  4,  5,
  12, 19, 26, 33, 40, 48, 41, 34,
  27, 20, 13,  6,  7, 14, 21, 28,
  35, 42, 49, 56, 57, 50, 43, 36,
  29, 22, 15, 23, 30, 37, 44, 51,
  58, 59, 52, 45, 38, 31, 39, 46,
  53, 60, 61, 54, 47, 55, 62, 63
};



template<typename kind>
inline std::vector< kind > forwardZigzag(const std::vector< kind > &input)
{
  assert(input.size() == 64);
  
  std::vector< kind > output(64);
  
  for(int i = 0; i < 64; i++)
    {
      output[i] = input[ ZIGZAG_TABLE[i] ];
    }
  
  return output;
}


template<typename kind>
inline std::vector< kind > inverseZigzag(const std::vector< kind > &input)
{
  assert(input.size() == 64);
  
  std::vector< kind > output(64);

  for(int i = 0; i < 64; i++)
    {
      output[ ZIGZAG_TABLE[i] ] = input[i];
    }
  
  return output;
}

inline std::vector< double > levelShift(const std::vector < double > &input, int level)
{
  std::vector< double > output(input.size());
  for(int i = 0; i < input.size(); i++)
    {
      output[i] = input[i] + level;
    }

  return output;
  
}

inline mylib::Vector_2D< double > levelShift(const mylib::Vector_2D< double > &input, int level)
{
  mylib::Vector_2D< double > output(input.size_rows(), 0);
  
  for(int row = 0; row < input.size_rows(); row++)
    {
      output(row) = levelShift(input(row), level);
    }
  
  return output;
}

template<typename kind>
inline std::vector< kind > divideByQTable(const std::vector< kind > &input, const std::vector< kind > &qTable)
{
  assert(input.size() == qTable.size());
  
  std::vector< kind > output(input.size());
  for(int i = 0; i < input.size(); i++)
    {
      output[i] = input[i] / qTable[i];
    }

  return output;
}

template<typename kind>
inline mylib::Vector_2D< kind > divideByQTable(const mylib::Vector_2D< kind > &input,
                                               const mylib::Vector_2D< kind > &qTable)
{
  assert(input.size_rows() == qTable.size_rows());

  mylib::Vector_2D< kind > output(input.size_rows());
  for(int row = 0; row < input.size_rows(); row++)
    {
      output(row) = divideByQTable(input(row), qTable(row));
    }

  return output;
}

template<typename kind>
inline std::vector< kind > multiplyByQTable(const std::vector< kind > &input, const std::vector< u_int > &qTable)
{
  assert(input.size() == qTable.size());
  
  std::vector< kind > output(input.size());
  for(int i = 0; i < input.size(); i++)
    {
      output[i] = input[i] * qTable[i];
    }

  return output;
}

template<typename kind>
inline mylib::Vector_2D< kind > multiplyByQTable(const mylib::Vector_2D< kind > &input,
                                               const mylib::Vector_2D< u_int > &qTable)
{
  assert(input.size_rows() == qTable.size_rows());

  mylib::Vector_2D< kind > output(input.size_rows());
  for(int row = 0; row < input.size_rows(); row++)
    {
      output(row) = multiplyByQTable(input(row), qTable(row));
    }

  return output;
}

inline double calcMSE(const mylib::vec_2D &original, const mylib::vec_2D &object)
{
  assert(original.is_rectangular() && object.is_rectangular());
  assert(original.size_rows() == object.size_rows() && original.size_cols() == object.size_cols());

 int vSize = original.size_rows();
  int hSize = original.size_cols();
  double mse = 0.0;
  for(int v = 0; v < vSize; v++)
    {
      for(int h = 0; h < hSize; h++)
        {
          mse += (original(v,h) - object(v,h)) * (original(v,h) - object(v,h));
        }
    }
  mse /= static_cast<double>(vSize * hSize);
  return mse;
} 

// u_char version
inline double calcMSE(const mylib::Vector_2D< u_char > &original, const mylib::Vector_2D< u_char > &object)
{
  assert(original.is_rectangular() && object.is_rectangular());
  assert(original.size_rows() == object.size_rows() && original.size_cols() == object.size_cols());

  int vSize = original.size_rows();
  int hSize = original.size_cols();
  double mse = 0.0;
  for(int v = 0; v < vSize; v++)
    {
      for(int h = 0; h < hSize; h++)
        {
          u_char diff = (original(v,h) >= object(v,h)) ? (original(v,h) - object(v,h)) :
            (object(v,h) - original(v,h));
          mse += static_cast<double>(diff*diff);
        }
    }
  mse /= static_cast<double>(vSize * hSize);

  return mse;
}

// peak is the maximum possible value of pixel.
inline double calcPSNR(const mylib::vec_2D &original, const mylib::vec_2D &object, double peak)
{
  assert(original.is_rectangular() && object.is_rectangular());
  assert(original.size_rows() == object.size_rows() && original.size_cols() == object.size_cols());

  int vSize = original.size_rows();
  int hSize = original.size_cols();
  double mse = 0.0;
  for(int v = 0; v < vSize; v++)
    {
      for(int h = 0; h < hSize; h++)
        {
          assert(original(v,h) <= peak && object(v,h) <= peak);
          mse += (original(v,h) - object(v,h)) * (original(v,h) - object(v,h));
        }
    }
  mse /= static_cast<double>(vSize * hSize);
  if(mse == 0)
    {
      mse = 1.0/mylib::INFTY;
    }
  double psnr = 10 * log10( peak * peak / mse);
  return psnr;
}

// u_char version
inline double calcPSNR(const mylib::Vector_2D< u_char > &original, const mylib::Vector_2D< u_char > &object,
                       double peak)
{
  assert(original.is_rectangular() && object.is_rectangular());
  assert(original.size_rows() == object.size_rows() && original.size_cols() == object.size_cols());

  int vSize = original.size_rows();
  int hSize = original.size_cols();
  double mse = 0.0;
  for(int v = 0; v < vSize; v++)
    {
      for(int h = 0; h < hSize; h++)
        {
          assert(original(v,h) <= peak && object(v,h) <= peak);
          u_char diff = (original(v,h) >= object(v,h)) ? (original(v,h) - object(v,h)) :
            (object(v,h) - original(v,h));
          mse += static_cast<double>(diff*diff);
        }
    }
  mse /= static_cast<double>(vSize * hSize);
  if(mse == 0)
    {
      mse = 1.0 / mylib::INFTY;
    }
  double psnr = 10 * log10( peak * peak / mse);
  return psnr;
}


// Divides BMP pixel data into square areas.
template< typename kind >
class cDividePixels
{
protected:
  mylib::Vector_2D< kind > pixels_;
  int numPix_;
  int totalHsize_;
  int totalVsize_;
  int numHBlocks_;
  int numVBlocks_;
  
public:
  cDividePixels(const mylib::Vector_2D< kind >  &pixels, 
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
  
  mylib::Vector_2D< kind > getBlock(int Vblock, int Hblock)
  {
    assert(Vblock >= 0 && Vblock < numVBlocks_);
    assert(Hblock >= 0 && Hblock < numHBlocks_);
    
    int startRow = Vblock * numPix_;
    int startCol = Hblock * numPix_;
    
    mylib::Vector_2D< kind > output(numPix_, numPix_);
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
  
  int Hblocks()
  {
    return numHBlocks_;
  }
  int Vblocks()
  {
    return numVBlocks_;
  }
};

// Compose Vector_2D divided into blocks.
template< typename kind >
class cCompPixels
{
protected:
  mylib::Vector_2D< kind > pixels_;
  int totalHsize_;
  int totalVsize_;
  
public:
  cCompPixels(int totalVsize, int totalHsize):
    pixels_(totalVsize, totalHsize), totalHsize_(totalHsize), totalVsize_(totalVsize)
  {
  }
  
  void operator()(int startRow, int startCol, const mylib::Vector_2D< kind > &input)
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
  
  mylib::Vector_2D< kind > get()
  {
    return pixels_;
  }
};

#endif
