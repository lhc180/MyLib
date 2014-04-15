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

namespace mylib{
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

    // intかu_char
    template < typename kind >
    itpp::Mat< kind > GetGray();
    
    template < typename kind >
    std::vector < itpp::Mat< kind > > Get();
  
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
  template < typename kind >
  bool MakeBmp(const char* filename,
               const itpp::Mat< kind > &pixelVec);

  // for color scale
  template < typename kind >
  bool MakeBmp(const char* filename,
               const std::vector< itpp::Mat< kind > > &pixelVec);

  // 2次元配列にRGBRGB...の順に入れてある必要がある
  template < typename kind >
  bool MakeBmpRGB(const char* filename,
                  const itpp::Mat< kind > &pixelVec);

  template < typename kind >
  std::vector< kind > Rgb2Ycc(const std::vector< kind >& rgb);
  template < typename kind >
  std::vector< kind > Ycc2Rgb(const std::vector< kind >& ycc);
  template < typename kind >
  itpp::Vec< kind > Rgb2Ycc(const itpp::Vec< kind >& rgb);
  template < typename kind >
  itpp::Vec< kind > Ycc2Rgb(const itpp::Vec< kind >& ycc);

  template < typename kind >
  std::vector< Vector_2D< kind > > Rgb2Ycc(const std::vector< Vector_2D< kind > >& rgb);
  template < typename kind >
  std::vector< Vector_2D< kind > > Ycc2Rgb(const std::vector< Vector_2D< kind > >& ycc);
  template < typename kind >
  std::vector< itpp::Mat< kind > > Rgb2Ycc(const std::vector< itpp::Mat< kind > >& rgb);
  template < typename kind >
  std::vector< itpp::Mat< kind > > Ycc2Rgb(const std::vector< itpp::Mat< kind > >& ycc);
} // end of mylib

#endif

