#include <cassert>
#include "../include/mybmp.h"

// numComp_の値によってYUV変換したほうがいいかも
mylib::Vector_2D< u_char > BmpReader::GetGray()
{
  assert(setDone_);
  assert(reqComp_ == 1);

  mylib::Vector_2D< u_char > output(height_, width_);
  
  for(int h = 0; h < height_; h++){
    for(int w = 0; w < width_; w++){
      int i        = h*width_ + w;
      output(h, w) = pixels_[i];
    }
  }

  return output;
}

std::vector< mylib::Vector_2D< u_char > > BmpReader::Get()
{
  assert(setDone_);
  
  int numComp = reqComp_;

  std::vector< mylib::Vector_2D< u_char > > output(numComp);
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
             const mylib::Vector_2D< u_char > &pixelVec)
{
  assert(pixelVec.is_rectangular());
  int height = pixelVec.size_rows();
  int width  = pixelVec.size_cols(0);

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
             const std::vector< mylib::Vector_2D< u_char > > &pixelVec)
{
  int numComp = pixelVec.size();
  for(int c = 0; c < numComp; c++){
    assert(pixelVec[c].is_rectangular());
  }

  int height = pixelVec[0].size_rows();
  int width = pixelVec[0].size_cols(0);

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
                  const mylib::Vector_2D< u_char > &pixelVec)
{
  int numComp = 3;              // RGB
  assert(pixelVec.is_rectangular());

  int height = pixelVec.size_rows();
  int width = pixelVec.size_cols();

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

