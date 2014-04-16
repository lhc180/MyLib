/************************************************************************************
 * spiht.cpp
 *   
 * SPIHTクラスの実装
 *
 *
 * Last Updated: <2014/04/16 17:07:31 from dr-yst-no-pc.local by yoshito>
 ************************************************************************************/


#include "spiht.h"

namespace mylib {

  template < typename kind >
  inline itpp::bin SignBit(kind value)
  {
    if (value >= 0){
      return 0;
    } // if
    else{
      return 1;
    } // else
  }
  
  inline int SignFromBit(const itpp::bin& bit)
  {
    if (bit == 1){
      return -1;
    } // if bit
    else{
      return 1;
    } // else 
  }
  
  void SPIHTencoder::Initialize(const itpp::imat& image, itpp::bvec* bout) {
    lip_.clear();
    lsp_.clear();
    lis_.clear();
    bout->set_size(0);
    int max = 0;
    for (int y = 0; y < image.rows(); y++){
      for (int x = 0; x < image.cols(); x++){
        if (std::abs(image(x,y)) > max) {
          max = std::abs(image(x,y));
        }
      } // for x
    }   // for y
    
    step_ = (int)floor(log2(static_cast< double >(max)));
    itpp::bvec step_bin = itpp::dec2bin(8, step_);
    *bout = itpp::concat(*bout, step_bin);

    for (int y = 0; y < (image.rows()) / (1 << numStages_); y++) {
      for (int x = 0; x < (image.cols()) / (1 << numStages_); x++) {
        lip_.push_back(PixItem(x, y));
        if ((x % 2 != 0) || (y % 2 != 0)) // ## 多分ここ違う
          lis_.push_back(SetItem(x, y, LIS_A));
      } // for x
    }   // for y
    
  }

  void SPIHTencoder::GetSuccessor(const itpp::imat& image, int x, int y, int *sx, int *sy)
  {
    int lx = (image.cols()) / (1 << numStages_);
    int ly = (image.rows()) / (1 << numStages_);
    if (x < lx && y < ly) {
      if (x % 2 == 1){
        *sx = x + lx - 1;
      }
      else {
        *sx = x;
      }
      if (y % 2 == 1) {
        *sy = y + ly - 1;
      }
      else{
        *sy = y;
      }
      
      if (*sx == x && *sy == y) {
        *sx = -1;
        *sy = -1;
      }
    } else {
      *sx = 2 * x;
      *sy = 2 * y;
      if (*sx >= (image.cols()) || *sy >= (image.rows())) {
        *sx = -1;
        *sy = -1;
      }
    }
  }

  bool SPIHTencoder::IsSignificantPixel(const itpp::imat& image, int x, int y)
  {
    return (std::abs((int) image(x, y)) >= (1 << step_));
  }

  bool SPIHTencoder::IsSignificant_SetA(const itpp::imat& image, int x, int y, int count)
  {
    if (count > 1 && IsSignificantPixel(image, x, y)){
      return true;
    }
    
    int sx, sy;
    GetSuccessor(image, x, y, &sx, &sy);
    if (sx == -1 || sy == -1)
      return false;
    if (IsSignificant_SetA(image, sx, sy, count + 1))
      return true;
    else if (IsSignificant_SetA(image, sx + 1, sy, count + 1))
      return true;
    else if (IsSignificant_SetA(image, sx, sy + 1, count + 1))
      return true;
    else if (IsSignificant_SetA(image, sx + 1, sy + 1, count + 1))
      return true;
    return false;
  }

  bool SPIHTencoder::IsSignificant_SetB(const itpp::imat& image, int x, int y, int count)
  {
    if (count > 2 && IsSignificantPixel(image, x, y)){
      return true;
    }
    
    int sx, sy;
    GetSuccessor(image, x, y, &sx, &sy);
    if (sx == -1 || sy == -1)
      return false;
    if (IsSignificant_SetB(image, sx, sy, count + 1))
      return true;
    else if (IsSignificant_SetB(image, sx + 1, sy, count + 1))
      return true;
    else if (IsSignificant_SetB(image, sx, sy + 1, count + 1))
      return true;
    else if (IsSignificant_SetB(image, sx + 1, sy + 1, count + 1))
      return true;
    return false;
  }

  void SPIHTencoder::Encode(const itpp::imat& image, int bits, itpp::bvec *bout)
  {
    if (bits < 0){
      bits = std::numeric_limits<int>::max();
    } // if 
    Initialize(image, bout);
    int bit_cnt = 8;
    if (bits <= 8){
      *bout = bout->left(bits);
      return;
    } // if 
    
    while (step_ >= 0) {
      /* Sorting pass */
      /* first process LIP */
      for (int i = 0; i < static_cast< int >(lip_.size()); i++) {
        bool sig = IsSignificantPixel(image, lip_[i].x, lip_[i].y);
        *bout = itpp::concat(*bout, static_cast< itpp::bin >(sig));
        if (++bit_cnt > bits){
          return;
        } // if
        if (sig) {
          lsp_.push_back(PixItem(lip_[i].x, lip_[i].y));
          itpp::bin s = SignBit(image(lip_[i].x, lip_[i].y)); // sign
          *bout = itpp::concat(*bout, s);
          if (++bit_cnt > bits){
            return; 
          } // if
          lip_.erase(lip_.begin() + i);
          i--;
        } // if sig
      }   // for i
      /* now process LIS */
      for (int i = 0; i < static_cast< int >(lis_.size()); i++) {
        if (lis_[i].type == LIS_A) {
          bool sig = IsSignificant_SetA(image, lis_[i].x, lis_[i].y);
          *bout = itpp::concat(*bout, static_cast< itpp::bin >(sig));
          if (++bit_cnt > bits) return;
          if (sig) {
            int sx, sy;
            GetSuccessor(image, lis_[i].x, lis_[i].y, &sx, &sy);
            /* process the four offsprings */
            sig = IsSignificantPixel(image, sx, sy);
            *bout = itpp::concat(*bout, static_cast< itpp::bin >(sig));

            // for (sx, sy)
            if (++bit_cnt > bits) return;
            if (sig) {
              lsp_.push_back(PixItem(sx, sy));
              itpp::bin s = SignBit(image(sx,sy));
              *bout = itpp::concat(*bout, s);
              if (++bit_cnt > bits) return;
            } else {
              lip_.push_back(PixItem(sx, sy));
            }
            sig = IsSignificantPixel(image, sx + 1, sy);
            *bout = itpp::concat(*bout, static_cast< itpp::bin >(sig));

            // for (sx+1, sy)
            if (++bit_cnt > bits) return;
            if (sig) {
              lsp_.push_back(PixItem(sx + 1, sy));
              itpp::bin s = SignBit(image(sx + 1, sy));
              *bout = itpp::concat(*bout, s);
              if (++bit_cnt > bits) return;
            } else {
              lip_.push_back(PixItem(sx + 1, sy));
            }
            sig = IsSignificantPixel(image, sx, sy + 1);
            *bout = itpp::concat(*bout, static_cast< itpp::bin >(sig));
            
            // for (sx, sy+1)
            if (++bit_cnt > bits) return;
            if (sig) {
              lsp_.push_back(PixItem(sx, sy + 1));
              itpp::bin s = SignBit(image(sx, sy+1));
              *bout = itpp::concat(*bout, s);
              if (++bit_cnt > bits) return;
            } else {
              lip_.push_back(PixItem(sx, sy + 1));
            }
            sig = IsSignificantPixel(image, sx + 1, sy + 1);
            *bout = itpp::concat(*bout, static_cast< itpp::bin >(sig));

            // for (sx+1, sy+1)
            if (++bit_cnt > bits) return;
            if (sig) {
              lsp_.push_back(PixItem(sx + 1, sy + 1));
              itpp::bin s = SignBit(image(sx+1, sy+1));
              *bout = itpp::concat(*bout, s);
              if (++bit_cnt > bits) return;
            } else {
              lip_.push_back(PixItem(sx + 1, sy + 1));
            }
            /* test if L(i, j) != 0 */
            GetSuccessor(image, sx, sy, &sx, &sy);
            if (sx != -1)
              lis_.push_back(SetItem(lis_[i].x, lis_[i].y, LIS_B));
            lis_.erase(lis_.begin() + i);
            i--;
          }
        } else {
          bool sig = IsSignificant_SetB(image, lis_[i].x, lis_[i].y);
          *bout = itpp::concat(*bout, static_cast< itpp::bin >(sig));

          if (++bit_cnt > bits) return;
          if (sig) {
            int sx, sy;
            GetSuccessor(image, lis_[i].x, lis_[i].y, &sx, &sy);
            lis_.push_back(SetItem(sx, sy, LIS_A));
            lis_.push_back(SetItem(sx + 1, sy, LIS_A));
            lis_.push_back(SetItem(sx, sy + 1, LIS_A));
            lis_.push_back(SetItem(sx + 1, sy + 1, LIS_A));
            lis_.erase(lis_.begin() + i);
            i--;
          }
        }
      }
      /* Refinement pass */
      for (int i = 0; i < static_cast< int >(lsp_.size()); i++) {
        if (std::abs(image(lsp_[i].x, lsp_[i].y)) >= (1 << (step_ + 1))) {
          itpp::bin refinement = (std::abs(static_cast< int >(image(lsp_[i].x, lsp_[i].y))) >> step_) & 1;
          *bout = itpp::concat(*bout, refinement);
          if (++bit_cnt > bits) return;
        }
      }
      /* Quantization step update */
      step_--;
    }
  }

/*************************************************************************************************************/

  void SPIHTdecoder::Initialize(const itpp::bvec &bin, itpp::imat *image)
  {
    lip_.clear();
    lsp_.clear();
    lis_.clear();

    image->set_size(rows_, cols_);
    image->zeros();

    if (bin.size() < 8){
      return;
    } // if 

    itpp::bvec stepBin = bin.left(8);
    step_ = itpp::bin2dec(stepBin);
    for (int y = 0; y < rows_ / (1 << numStages_); y++)
      for (int x = 0; x < cols_ / (1 << numStages_); x++) {
        lip_.push_back(PixItem(x, y));
        if ((x % 2 != 0) || (y % 2 != 0))
          lis_.push_back(SetItem(x, y, LIS_A));
      }
  }

  void SPIHTdecoder::GetSuccessor(int x, int y, int* sx, int* sy) {
    int lx = cols_ / (1 << numStages_);
    int ly = rows_ / (1 << numStages_);
    if (x < lx && y < ly) {
      if (x % 2 == 1){
        *sx = x + lx - 1;
      }
      else{
        *sx = x;
      }
      if (y % 2 == 1){
        *sy = y + ly - 1;
      }
      else{
        *sy = y;
      }
      if (*sx == x && *sy == y) {
        *sx = -1;
        *sy = -1;
      }
    } else {
      *sx = 2 * x;
      *sy = 2 * y;
      if (*sx >= cols_ || *sy >= rows_) {
        *sx = -1;
        *sy = -1;
      }
    }
  }

  void SPIHTdecoder::Decode(const itpp::bvec& bin, itpp::imat *imageOut)
  {
    Initialize(bin, imageOut);
    int bits = bin.size();
    int bit_cnt = 8;

    if (bit_cnt > bits){
      return;
    } // if 
    
    while (step_ >= 0) {
      /* Sorting pass */
      /* first process LIP */
      for (int i = 0; i < static_cast< int >(lip_.size()); i++) {
        bool sig = static_cast< bool >(bin[bit_cnt]);
        if (++bit_cnt > bits) return;
        if (sig) {
          lsp_.push_back(PixItem(lip_[i].x, lip_[i].y));
          int s = SignFromBit(bin[bit_cnt]);
          (*imageOut)(lip_[i].x, lip_[i].y) = s * (1 << step_);
          if (++bit_cnt > bits) return;
          lip_.erase(lip_.begin() + i);
          i--;
        }
      }
      /* now process LIS */
      for (int i = 0; i < static_cast< int >(lis_.size()); i++) {
        if (lis_[i].type == LIS_A) {
          bool sig = static_cast< bool >(bin[bit_cnt]);
          if (++bit_cnt > bits) return;
          if (sig) {
            int sx, sy;
            GetSuccessor(lis_[i].x, lis_[i].y, &sx, &sy);
            /* process the four offsprings */
            // for (sx, sy)
            sig = static_cast< bool >(bin[bit_cnt]);
            if (++bit_cnt > bits) return;
            if (sig) {
              lsp_.push_back(PixItem(sx, sy));
              int s = SignFromBit(bin[bit_cnt]);
              (*imageOut)(sx, sy) = s * (1 << step_);
              if (++bit_cnt > bits) return;
            } else {
              lip_.push_back(PixItem(sx, sy));
            }
            
            // for (sx+1, sy)
            sig = static_cast< bool >(bin[bit_cnt]);
            if (++bit_cnt > bits) return;
            if (sig) {
              lsp_.push_back(PixItem(sx + 1, sy));
              int s = SignFromBit(bin[bit_cnt]);
              (*imageOut)(sx + 1, sy) = s * (1 << step_);
              if (++bit_cnt > bits) return;
            } else {
              lip_.push_back(PixItem(sx + 1, sy));
            }

            // for (sx, sy+1)
            sig = static_cast< bool >(bin[bit_cnt]);
            if (++bit_cnt > bits) return;
            if (sig) {
              lsp_.push_back(PixItem(sx, sy + 1));
              int s = SignFromBit(bin[bit_cnt]);
              (*imageOut)(sx, sy + 1) = s * (1 << step_);
              if (++bit_cnt > bits) return;
            } else {
              lip_.push_back(PixItem(sx, sy + 1));
            }

            // for (sx+1, sy+1)
            sig = static_cast< bool >(bin[bit_cnt]);
            if (++bit_cnt > bits) return;
            if (sig) {
              lsp_.push_back(PixItem(sx + 1, sy + 1));
              int s = SignFromBit(bin[bit_cnt]);
              (*imageOut)(sx + 1, sy + 1) = s * (1 << step_);
              if (++bit_cnt > bits) return;
            } else {
              lip_.push_back(PixItem(sx + 1, sy + 1));
            }
            /* test if L(i, j) != 0 */
            GetSuccessor(sx, sy, &sx, &sy);
            if (sx != -1)
              lis_.push_back(SetItem(lis_[i].x, lis_[i].y, LIS_B));
            lis_.erase(lis_.begin() + i);
            i--;
          }
        } else {
          bool sig = static_cast< bool >( bin[bit_cnt]);
          if (++bit_cnt > bits) return;
          if (sig) {
            int sx, sy;
            GetSuccessor(lis_[i].x, lis_[i].y, &sx, &sy);
            lis_.push_back(SetItem(sx, sy, LIS_A));
            lis_.push_back(SetItem(sx + 1, sy, LIS_A));
            lis_.push_back(SetItem(sx, sy + 1, LIS_A));
            lis_.push_back(SetItem(sx + 1, sy + 1, LIS_A));
            lis_.erase(lis_.begin() + i);
            i--;
          }
        }
      }
      /* Refinement pass */
      for (int i = 0; i < static_cast< int >(lsp_.size()); i++) {
        if (std::abs((*imageOut)(lsp_[i].x, lsp_[i].y)) >= (1 << (step_ + 1))) {
          if (static_cast< bool >(bin[bit_cnt])) {
            if ((*imageOut)(lsp_[i].x, lsp_[i].y) >= 0)
              (*imageOut)(lsp_[i].x, lsp_[i].y) = static_cast< int >((*imageOut)(lsp_[i].x, lsp_[i].y)) | (1 << step_);
            else
              (*imageOut)(lsp_[i].x, lsp_[i].y) = static_cast< int >(-(std::abs((*imageOut)(lsp_[i].x, lsp_[i].y)) | (1 << step_)));
          }
          else {
            (*imageOut)(lsp_[i].x, lsp_[i].y) =  static_cast< int >((*imageOut)(lsp_[i].x, lsp_[i].y)) & (~(1 << step_));
          }
          if (++bit_cnt > bits) return;
        } // if abs
      } // for i
      /* Quantization step update */
      step_--;
    }
  }
  
  
}
