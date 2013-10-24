#include <cassert>
#include "../include/mybmp.h"

namespace mylib{
  // numComp_���ͤˤ�ä�YUV�Ѵ������ۤ�����������
  itpp::Mat< u_char > BmpReader::GetGray()
  {
    assert(setDone_);
    assert(reqComp_ == 1);

    itpp::Mat< u_char > output(height_, width_);
  
    for(int h = 0; h < height_; h++){
      for(int w = 0; w < width_; w++){
        int i        = h*width_ + w;
        output(h, w) = pixels_[i];
      }
    }

    return output;
  }
  
  std::vector< itpp::Mat< u_char > > BmpReader::Get()
  {
    assert(setDone_);
  
    int numComp = reqComp_;

    std::vector< itpp::Mat< u_char > > output(numComp);
    for(int i = 0; i < numComp; i++){
      output[i].set_size(height_, width_);
    }

    int i                    = 0;
    for(int h = 0; h < height_; h++){
      for(int w = 0; w < width_; w++){
        for(int comp = 0; comp < numComp; comp++){
          output[comp](h, w) = pixels_[i];
          i++;
        }                         // for comp
      }                           // for w
    }                             // for h

    return output;
  }

  // for grey scale
  bool MakeBmp(const char* filename, 
               const itpp::Mat< u_char > &pixelVec)
  {
    int height = pixelVec.rows();
    int width  = pixelVec.cols();

    u_char* pixels = new u_char[height*width];
  
    int i         = 0;
    for(int h = 0; h < height; h++){
      for(int w = 0; w < width; w++){
        pixels[i] = pixelVec(h,w);
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
  
  // for color scale
  bool MakeBmp(const char* filename,
               const std::vector< itpp::Mat< u_char > > &pixelVec)
  {
    int numComp = pixelVec.size();

    int height = pixelVec[0].rows();
    int width = pixelVec[0].cols();

    u_char* pixels = new u_char[height*width*numComp];

    int i = 0;
    for(int h = 0; h < height; h++){
      for(int w = 0; w < width; w++){
        for(int c = 0; c < numComp; c++){
          pixels[i] = pixelVec[c](h,w);
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

  bool MakeBmpRGB(const char* filename,
                  const itpp::Mat< u_char > &pixelVec)
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
            pixels[i] = pixelVec(h,w);
            i++;
          } // for w
      }     // for h

    int ret = stbi_write_bmp(filename, width/numComp, height, numComp, pixels);

    stbi_image_free(pixels);

    if(ret == 0)
      {
        return false;
      }
    else
      {
        return true;
      }
  }

  // �����ͤ�0~255����������
  inline u_char ColorCompensation(double input)
  {
    if(input >= 255){
      return 255;
    }
    else if(input <= 0){
      return 0;
    }
    else{
      return static_cast< u_char >(input);
    }
  }
  
  std::vector< u_char > Rgb2Ycc(const std::vector< u_char >& rgb)
  {
    assert(rgb.size() == 3);
    std::vector< u_char > ycc(3);
    double r = rgb[0], g = rgb[1], b = rgb[2];
    
    ycc[0] = ColorCompensation(0.2990*r + 0.5870*g + 0.1140*b);
    ycc[1] = ColorCompensation(-0.1687*r - 0.3313*g + 0.5000*b + 128);
    ycc[2] = ColorCompensation(0.5000*r - 0.4187*g - 0.0813*b + 128);

    return ycc;
  }

  itpp::Vec< u_char > Rgb2Ycc(const itpp::Vec< u_char >& rgb)
  {
    assert(rgb.size() == 3);
    itpp::Vec< u_char > ycc(3);
    double r = rgb[0], g = rgb[1], b = rgb[2];
    
    ycc[0] = ColorCompensation(0.2990*r + 0.5870*g + 0.1140*b);
    ycc[1] = ColorCompensation(-0.1687*r - 0.3313*g + 0.5000*b + 128);
    ycc[2] = ColorCompensation(0.5000*r - 0.4187*g - 0.0813*b + 128);

    return ycc;
  }

  std::vector< u_char > Ycc2Rgb(const std::vector< u_char >& ycc)
  {
    assert(ycc.size() == 3);
    std::vector< u_char > rgb(3);
    double y = ycc[0], cb = ycc[1], cr = ycc[2];

    rgb[0] = ColorCompensation(y + 1.40200*(cr - 128));
    rgb[1] = ColorCompensation(y - 0.34414*(cb - 128) - 0.71414*(cr - 128));
    rgb[2] = ColorCompensation(y + 1.77200*(cb - 128));

    return rgb;
  }

  itpp::Vec< u_char > Ycc2Rgb(const itpp::Vec< u_char >& ycc)
  {
    assert(ycc.size() == 3);
    itpp::Vec< u_char > rgb(3);
    double y = ycc[0], cb = ycc[1], cr = ycc[2];

    rgb[0] = ColorCompensation(y + 1.40200*(cr - 128));
    rgb[1] = ColorCompensation(y - 0.34414*(cb - 128) - 0.71414*(cr - 128));
    rgb[2] = ColorCompensation(y + 1.77200*(cb - 128));

    return rgb;
  }

  
  std::vector< Vector_2D< u_char > > Rgb2Ycc(const std::vector< Vector_2D< u_char > >& rgb)
  {
    assert(rgb.size() == 3);
    // �����ǤΥ�������Ʊ�����ɤ����Υ����å��Ϥ��ɤ������ΤǤ��ʤ�
    int height = rgb[0].Height();
    int width = rgb[0].Width();

    std::vector< Vector_2D< u_char > > ycc(3, Vector_2D< u_char >(height, width));
    
    for (int row = 0; row < height; ++row){
      for (int col = 0; col < width; ++col){
        std::vector< u_char > components(3);
        for (int i = 0; i < 3; ++i){
          components[i] = rgb[i](row, col);
        } // for i
        std::vector< u_char > transformed = Rgb2Ycc(components);
        for (int i = 0; i < 3; ++i){
          ycc[i](row, col) = transformed[i];
        } // for i
      } // for col
    } // for row

    return ycc;
  }

  std::vector< itpp::Mat< u_char > > Rgb2Ycc(const std::vector< itpp::Mat< u_char > >& rgb)
  {
    assert(rgb.size() == 3);
    // �����ǤΥ�������Ʊ�����ɤ����Υ����å��Ϥ��ɤ������ΤǤ��ʤ�
    int height = rgb[0].rows();
    int width = rgb[0].cols();

    std::vector< itpp::Mat< u_char > > ycc(3, itpp::Mat< u_char >(height, width));
    
    for (int row = 0; row < height; ++row){
      for (int col = 0; col < width; ++col){
        std::vector< u_char > components(3);
        for (int i = 0; i < 3; ++i){
          components[i] = rgb[i](row, col);
        } // for i
        std::vector< u_char > transformed = Rgb2Ycc(components);
        for (int i = 0; i < 3; ++i){
          ycc[i](row, col) = transformed[i];
        } // for i
      } // for col
    } // for row

    return ycc;
  }

  
  std::vector< Vector_2D< u_char > > Ycc2Rgb(const std::vector< Vector_2D< u_char > >& ycc)
  {
    assert(ycc.size() == 3);
    // �����ǤΥ�������Ʊ�����ɤ����Υ����å��Ϥ��ɤ������ΤǤ��ʤ�
    int height = ycc[0].Height();
    int width = ycc[0].Width();

    std::vector< Vector_2D< u_char > > rgb(3, Vector_2D< u_char >(height, width));
    
    for (int row = 0; row < height; ++row){
      for (int col = 0; col < width; ++col){
        std::vector< u_char > components(3);
        for (int i = 0; i < 3; ++i){
          components[i] = ycc[i](row, col);
        } // for i
        std::vector< u_char > transformed = Ycc2Rgb(components);
        for (int i = 0; i < 3; ++i){
          rgb[i](row, col) = transformed[i];
        } // for i
      } // for col
    } // for row

    return rgb;
  }

  std::vector< itpp::Mat< u_char > > Ycc2Rgb(const std::vector< itpp::Mat< u_char > >& ycc)
  {
    assert(ycc.size() == 3);
    // �����ǤΥ�������Ʊ�����ɤ����Υ����å��Ϥ��ɤ������ΤǤ��ʤ�
    int height = ycc[0].rows();
    int width = ycc[0].cols();

    std::vector< itpp::Mat< u_char > > rgb(3, itpp::Mat< u_char >(height, width));
    
    for (int row = 0; row < height; ++row){
      for (int col = 0; col < width; ++col){
        std::vector< u_char > components(3);
        for (int i = 0; i < 3; ++i){
          components[i] = ycc[i](row, col);
        } // for i
        std::vector< u_char > transformed = Ycc2Rgb(components);
        for (int i = 0; i < 3; ++i){
          rgb[i](row, col) = transformed[i];
        } // for i
      } // for col
    } // for row

    return rgb;
  }

  
} // end of mylib
