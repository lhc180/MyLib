// Tiff writer is not implemented.

#ifndef MYTIFF_H                
#define MYTIFF_H

#include <string>
#include <vector>
#include <cstdlib>
#include <tiffio.h>
#include "mymatrix.h"

class TiffReader
{
protected:
  uint32* pixels_;
  uint32 width_;
  uint32 height_;
  uint16 numComp_;
  bool setDone_;

public:
  TiffReader(): setDone_(false)
  {
  }

  TiffReader(const char* filename)
  {
    Set(filename);
  }

  ~TiffReader()
  {
    Close();
  }
  
  void Set(const char* filename)
  {
    TIFF* tif = TIFFOpen(filename, "r");
    if(tif)
      {
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width_);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height_);
        size_t numPixels = width_*height_;
        pixels_ = (uint32*) _TIFFmalloc(numPixels * sizeof(uint32));
        if(pixels_ != NULL)
          {
            if(!TIFFReadRGBAImageOriented(tif, width_, height_, pixels_, ORIENTATION_TOPLEFT, 0))
              {
                std::cerr << "Error: can not read RGBA data from \""
                          << filename
                          << "\"." << std::endl;
                exit(1);
              }
            else
              {
                TIFFClose(tif);
              }
          }
        else
          {
            std::cerr << "Error: can not implement _TIFFmalloc in cTiffReader." << std::endl;
            exit(1);
          }
      }
    else
      {
        std::cerr << "Error: can not read file \"" << filename <<"\"." << std::endl;
        exit(1);
      }
    setDone_ = true;
  }

  mylib::Vector_2D< u_char > GetGray();
  std::vector< mylib::Vector_2D< u_char > > GetRGB();

  void Close()
  {
    _TIFFfree(pixels_);
  }

  int Width()
  {
    return width_;
  }

  int Height()
  {
    return height_;
  }

};

#endif
