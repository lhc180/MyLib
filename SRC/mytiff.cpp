#include "../include/mytiff.h"

namespace mylib{
  // 8 bit unsigned.
  // return Y of YUV.
  itpp::Mat< u_char > TiffReader::GetGray()
  {
    assert(setDone_);

    itpp::Mat< u_char > output(height_, width_);

    for(int h = 0; h < static_cast<int>(height_); h++)
      {
        for(int w = 0; w < static_cast<int>(width_); w++)
          {
            int i = h*width_ + w;
            int r = TIFFGetR(pixels_[i]);
            int g = TIFFGetG(pixels_[i]);
            int b = TIFFGetG(pixels_[i]);
            output(h, w) = static_cast< u_char >( 0.2990*r + 0.5870*g + 0.1140*b);
          }
      }

    return output;
  }

  std::vector< itpp::Mat< u_char > > TiffReader::GetRGB()
  {
    assert(setDone_);

    std::vector< itpp::Mat< u_char > > output(3);
    for(int i = 0; i < 3; i++)
      {
        output[i].set_size(height_, width_);
      }

    for(int h = 0; h < static_cast<int>(height_); h++)
      {
        for(int w = 0; w < static_cast<int>(width_); w++)
          {
            int i = h*width_ + w;
            for(int comp = 0; comp < 3; comp++)
              {
                output[comp](h, w) = ((pixels_[i] >> (8*comp)) & 0xff);
              }
          }
      }

    return output;
  }
}
