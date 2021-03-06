#include <cassert>
#include "../include/mybmp.h"

namespace mylib{
  // numComp_の値によってYUV変換したほうがいいかも
  template itpp::imat BmpReader::GetGray();
  template itpp::Mat< u_char > BmpReader::GetGray();

  template < typename kind >
  itpp::Mat< kind > BmpReader::GetGray()
  {
    assert(setDone_);
    assert(reqComp_ == 1);

    itpp::Mat< kind > output(height_, width_);
  
    for(int h = 0; h < height_; h++){
      for(int w = 0; w < width_; w++){
        int i        = h*width_ + w;
        output(h, w) = static_cast< kind >(pixels_[i]);
      }
    }

    return output;
  }

  template
  std::vector< itpp::Mat< int > > BmpReader::Get();
  template
  std::vector< itpp::Mat< u_char > > BmpReader::Get();

  template < typename kind >
  std::vector< itpp::Mat< kind > > BmpReader::Get()
  {
    assert(setDone_);
  
    int numComp = reqComp_;

    std::vector< itpp::Mat< kind > > output(numComp);
    for(int i = 0; i < numComp; i++){
      output[i].set_size(height_, width_);
    }

    int i                    = 0;
    for(int h = 0; h < height_; h++){
      for(int w = 0; w < width_; w++){
        for(int comp = 0; comp < numComp; comp++){
          output[comp](h, w) = static_cast< kind >(pixels_[i]);
          i++;
        }                         // for comp
      }                           // for w
    }                             // for h

    return output;
  }


  template
  bool MakeBmp(const char* filename, 
               const itpp::Mat< int > &pixelVec);
  template
  bool MakeBmp(const char* filename, 
               const itpp::Mat< u_char > &pixelVec);
    
  // for grey scale
  template < typename kind >
  bool MakeBmp(const char* filename, 
               const itpp::Mat< kind > &pixelVec)
  {
    int height = pixelVec.rows();
    int width  = pixelVec.cols();

    u_char* pixels = new u_char[height*width];
  
    int i         = 0;
    for(int h = 0; h < height; h++){
      for(int w = 0; w < width; w++){
        pixels[i] = static_cast< u_char >(pixelVec(h,w));
        i++;
      } // for w
    }   // for h

    int ret = stbi_write_bmp(filename, width, height, 1, pixels);
    
    stbi_image_free(pixels);

    if(ret == 0){
      return false;
    }
    else{
      return true;
    }
  }

  template
  bool MakeBmp(const char* filename,
               const std::vector< itpp::Mat< int > > &pixelVec);
  template
  bool MakeBmp(const char* filename,
               const std::vector< itpp::Mat< u_char > > &pixelVec);
  
  // for color scale
  template < typename kind >
  bool MakeBmp(const char* filename,
               const std::vector< itpp::Mat< kind > > &pixelVec)
  {
    int numComp = pixelVec.size();

    int height = pixelVec[0].rows();
    int width = pixelVec[0].cols();

    u_char* pixels = new u_char[height*width*numComp];

    int i = 0;
    for(int h = 0; h < height; h++){
      for(int w = 0; w < width; w++){
        for(int c = 0; c < numComp; c++){
          pixels[i] = static_cast< u_char >(pixelVec[c](h,w));
          i++;
        } // for c
      }   // for w
    }     // for h
  
    int ret = stbi_write_bmp(filename, width, height, numComp, pixels);

    stbi_image_free(pixels);

    if(ret == 0){
      return false;
    }
    else{
      return true;
    }
  }

  template
  bool MakeBmpRGB(const char* filename,
                  const itpp::Mat< int > &pixelVec);
  template
  bool MakeBmpRGB(const char* filename,
                  const itpp::Mat< u_char > &pixelVec);

  template < typename kind >
  bool MakeBmpRGB(const char* filename,
                  const itpp::Mat< kind > &pixelVec)
  {
    int numComp = 3;              // RGB

    int height = pixelVec.rows();
    int width = pixelVec.cols();

    u_char* pixels = new u_char[height * width];

    int i = 0;
    for(int h = 0; h < height; h++)
      {
        for(int w = 0; w < width; w++)
          {
            pixels[i] = static_cast< kind >(pixelVec(h,w));
            i++;
          } // for w
      }     // for h

    int ret = stbi_write_bmp(filename, width/numComp, height, numComp, pixels);

    stbi_image_free(pixels);

    if(ret == 0) {
      return false;
    }
    else {
      return true;
    }
  }
  
  // 入力値を0~255に補正する
  template < typename kind >
  inline kind ColorCompensation(double input)
  {
    if(input >= 255){
      return 255;
    }
    else if(input <= 0){
      return 0;
    }
    else{
      return static_cast< kind >(input);
    }
  }

  template
  std::vector< int > Rgb2Ycc(const std::vector< int >& rgb);
  template
  std::vector< u_char > Rgb2Ycc(const std::vector< u_char >& rgb);
    
  template < typename kind >
  std::vector< kind > Rgb2Ycc(const std::vector< kind >& rgb)
  {
    assert(rgb.size() == 3);
    std::vector< kind > ycc(3);
    double r = rgb[0], g = rgb[1], b = rgb[2];
    
    ycc[0] = ColorCompensation< kind >(0.2990*r + 0.5870*g + 0.1140*b);
    ycc[1] = ColorCompensation< kind >(-0.1687*r - 0.3313*g + 0.5000*b + 128);
    ycc[2] = ColorCompensation< kind >(0.5000*r - 0.4187*g - 0.0813*b + 128);

    return ycc;
  }

  template
  itpp::Vec< int > Rgb2Ycc(const itpp::Vec< int >& rgb);
  template
  itpp::Vec< u_char > Rgb2Ycc(const itpp::Vec< u_char >& rgb);
    
  template < typename kind >
  itpp::Vec< kind > Rgb2Ycc(const itpp::Vec< kind >& rgb)
  {
    assert(rgb.size() == 3);
    itpp::Vec< kind > ycc(3);
    double r = rgb[0], g = rgb[1], b = rgb[2];
    
    ycc[0] = ColorCompensation< kind >(0.2990*r + 0.5870*g + 0.1140*b);
    ycc[1] = ColorCompensation< kind >(-0.1687*r - 0.3313*g + 0.5000*b + 128);
    ycc[2] = ColorCompensation< kind >(0.5000*r - 0.4187*g - 0.0813*b + 128);

    return ycc;
  }

  template
  std::vector< int > Ycc2Rgb(const std::vector< int >& ycc);
  template
  std::vector< u_char > Ycc2Rgb(const std::vector< u_char >& ycc);
    
  template < typename kind >
  std::vector< kind > Ycc2Rgb(const std::vector< kind >& ycc)
  {
    assert(ycc.size() == 3);
    std::vector< kind > rgb(3);
    double y = ycc[0], cb = ycc[1], cr = ycc[2];

    rgb[0] = ColorCompensation< kind >(y + 1.40200*(cr - 128));
    rgb[1] = ColorCompensation< kind >(y - 0.34414*(cb - 128) - 0.71414*(cr - 128));
    rgb[2] = ColorCompensation< kind >(y + 1.77200*(cb - 128));

    return rgb;
  }

  template
  itpp::Vec< int > Ycc2Rgb(const itpp::Vec< int >& ycc);
  template
  itpp::Vec< u_char > Ycc2Rgb(const itpp::Vec< u_char >& ycc);

  template < typename kind >
  itpp::Vec< kind > Ycc2Rgb(const itpp::Vec< kind >& ycc)
  {
    assert(ycc.size() == 3);
    itpp::Vec< kind > rgb(3);
    double y = ycc[0], cb = ycc[1], cr = ycc[2];

    rgb[0] = ColorCompensation< kind >(y + 1.40200*(cr - 128));
    rgb[1] = ColorCompensation< kind >(y - 0.34414*(cb - 128) - 0.71414*(cr - 128));
    rgb[2] = ColorCompensation< kind >(y + 1.77200*(cb - 128));

    return rgb;
  }

  template
  std::vector< Vector_2D< int > > Rgb2Ycc(const std::vector< Vector_2D< int > >& rgb);
  template
  std::vector< Vector_2D< u_char > > Rgb2Ycc(const std::vector< Vector_2D< u_char > >& rgb);
    
  template < typename kind >
  std::vector< Vector_2D< kind > > Rgb2Ycc(const std::vector< Vector_2D< kind > >& rgb)
  {
    assert(rgb.size() == 3);
    // 各要素のサイズが同じかどうかのチェックはめんどくさいのでやらない
    int height = rgb[0].Height();
    int width = rgb[0].Width();

    std::vector< Vector_2D< kind > > ycc(3, Vector_2D< kind >(height, width));
    
    for (int row = 0; row < height; ++row){
      for (int col = 0; col < width; ++col){
        std::vector< kind > components(3);
        for (int i = 0; i < 3; ++i){
          components[i] = rgb[i](row, col);
        } // for i
        std::vector< kind > transformed = Rgb2Ycc(components);
        for (int i = 0; i < 3; ++i){
          ycc[i](row, col) = transformed[i];
        } // for i
      } // for col
    } // for row

    return ycc;
  }


  template
  std::vector< itpp::Mat< int > > Rgb2Ycc(const std::vector< itpp::Mat< int > >& rgb);
  template
  std::vector< itpp::Mat< u_char > > Rgb2Ycc(const std::vector< itpp::Mat< u_char > >& rgb);
    
  template < typename kind >
  std::vector< itpp::Mat< kind > > Rgb2Ycc(const std::vector< itpp::Mat< kind > >& rgb)
  {
    assert(rgb.size() == 3);
    // 各要素のサイズが同じかどうかのチェックはめんどくさいのでやらない
    int height = rgb[0].rows();
    int width = rgb[0].cols();

    std::vector< itpp::Mat< kind > > ycc(3, itpp::Mat< kind >(height, width));
    
    for (int row = 0; row < height; ++row){
      for (int col = 0; col < width; ++col){
        std::vector< kind > components(3);
        for (int i = 0; i < 3; ++i){
          components[i] = rgb[i](row, col);
        } // for i
        std::vector< kind > transformed = Rgb2Ycc(components);
        for (int i = 0; i < 3; ++i){
          ycc[i](row, col) = transformed[i];
        } // for i
      } // for col
    } // for row

    return ycc;
  }

  template
  std::vector< Vector_2D< int > > Ycc2Rgb(const std::vector< Vector_2D< int > >& ycc);
  template
  std::vector< Vector_2D< u_char > > Ycc2Rgb(const std::vector< Vector_2D< u_char > >& ycc);
    
  template < typename kind >
  std::vector< Vector_2D< kind > > Ycc2Rgb(const std::vector< Vector_2D< kind > >& ycc)
  {
    assert(ycc.size() == 3);
    // 各要素のサイズが同じかどうかのチェックはめんどくさいのでやらない
    int height = ycc[0].Height();
    int width = ycc[0].Width();

    std::vector< Vector_2D< kind > > rgb(3, Vector_2D< kind >(height, width));
    
    for (int row = 0; row < height; ++row){
      for (int col = 0; col < width; ++col){
        std::vector< kind > components(3);
        for (int i = 0; i < 3; ++i){
          components[i] = ycc[i](row, col);
        } // for i
        std::vector< kind > transformed = Ycc2Rgb(components);
        for (int i = 0; i < 3; ++i){
          rgb[i](row, col) = transformed[i];
        } // for i
      } // for col
    } // for row

    return rgb;
  }

  template
  std::vector< itpp::Mat< int > > Ycc2Rgb(const std::vector< itpp::Mat< int > >& ycc);
  template
  std::vector< itpp::Mat< u_char > > Ycc2Rgb(const std::vector< itpp::Mat< u_char > >& ycc);
    
  template<typename kind>
  std::vector< itpp::Mat< kind > > Ycc2Rgb(const std::vector< itpp::Mat< kind > >& ycc)
  {
    assert(ycc.size() == 3);
    // 各要素のサイズが同じかどうかのチェックはめんどくさいのでやらない
    int height = ycc[0].rows();
    int width = ycc[0].cols();

    std::vector< itpp::Mat< kind > > rgb(3, itpp::Mat< kind >(height, width));
    
    for (int row = 0; row < height; ++row){
      for (int col = 0; col < width; ++col){
        std::vector< kind > components(3);
        for (int i = 0; i < 3; ++i){
          components[i] = ycc[i](row, col);
        } // for i
        std::vector< kind > transformed = Ycc2Rgb(components);
        for (int i = 0; i < 3; ++i){
          rgb[i](row, col) = transformed[i];
        } // for i
      } // for col
    } // for row

    return rgb;
  }

  
} // end of mylib
