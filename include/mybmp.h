/*****************************************
 * This code provides interface of BMP file handler.
 * Suited for C++ vector.
 *****************************************/

#ifndef MYBMP_H
#define MYBMP_H

#include <string>
#include <vector>
#include "mymatrix.h"
#include "myutl.h"
#include "stbi_load.h"
#include "stbi_write.h"

class BmpReader
{
protected:
  u_char* pixels_;
  int     width_;
  int     height_;
  int     reqComp_;
  bool    setDone_;

public:
  BmpReader(): setDone_(false)
  { }

  BmpReader(const char* filename, int reqComp = 0)
  {
    Set(filename, reqComp);
  }

  ~BmpReader()
  {
    Close();
  }

  void Set(const char* filename, int reqComp = 0) // comp = 1: grey,
                            //        2: grey, alpha
                            //        3: red, green, blue
                            //        4: red, green, blue, alpha
  {
    assert((reqComp >= 0) && (reqComp <= 4));
    reqComp_ = reqComp;
    int numComp;
    pixels_  = stbi_load(filename, &width_, &height_, &numComp, reqComp_);
    if(reqComp == 0){
      reqComp_ = numComp;
    }
    setDone_ = true;
  }
  
  mylib::Vector_2D< u_char > GetGray();  
  std::vector < mylib::Vector_2D< u_char > > Get();
  
  void Close()
  {
    stbi_image_free(pixels_);
  }

  int Width()
  { return width_; }
  
  int Height()
  { return height_; }

  int Components()
  { return reqComp_; }
};

// For making BMP files, we do not define a class, but functions.

// for grey scale
bool MakeBmp(const char* filename,
                 const mylib::Vector_2D< u_char > &pixelVec);

// for color scale
bool MakeBmp(const char* filename,
             const std::vector< mylib::Vector_2D< u_char > > &pixelVec);

// 2次元配列にRGBRGB...の順に入れてある必要がある
bool MakeBmpRGB(const char* filename,
                  const mylib::Vector_2D< u_char > &pixelVec);

#endif










