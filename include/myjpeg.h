/***************************************
 * This file defines my original JPEG decoder class using libjpeg.
 * 
 * Last Updated: <2014/09/12 15:39:30 from WatanabeYoshito-no-iMac.local by yoshito>
 ***************************************/
// To do list
// -- 行数を指定してデータを読み込んでいく関数を作る.
// jpeg_read_scanlines()を使う
// -- テキストのエンコーダでバイナリデータを作
// らないといけない
// -- JpegFileのclose()は最後かも
// -- Huffman符号データをitpp::binで取り出す
// -- Get Huffman codewords as they are every times appearing RST marker.
// -- raw data (DCT係数)から画像を復元する関数を使う(生成されるのはJPEGファイル)
// -- DCTの式を見直して丸め誤差を見る
// -- もしかしたら、DC係数がDPCMかもしれないのでそれも確認する

#ifndef MYJPEG_H
#define MYJPEG_H

#include <iostream>
#include <jpeglib.h>
#include <csetjmp>
#include <cstdio>
// #include <cassert>
#include "mymatrix.h"


namespace mylib{
  
  static const int IMG_DCT_SIZE      = 8;
  static const int IMG_DCT_SIZE2     = 64;
  static const int ZIGZAG_TABLE_SIZE = 64;

  //  static const int ZRL_INDEX = 151;    // ## ハフマン符号化テーブルによっては違うかも
  static const int ZRL_INDEX = 0xF0;    // ## jpeglibによるとこの値
  
  static const int EOB_INDEX = 0;

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
  std::vector< kind > ForwardZigzag(const std::vector< kind > &input)
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
  itpp::Vec< kind > ForwardZigzag(const itpp::Vec< kind > &input)
  {
    assert(input.size() == 64);
  
    itpp::Vec< kind > output(64);
  
    for(int i = 0; i < 64; i++)
      {
        output[i] = input[ ZIGZAG_TABLE[i] ];
      }
  
    return output;
  }

  template < typename kind >
  itpp::Vec< kind > ForwardZigzagFromMat(const itpp::Mat< kind >& input)
  {
    assert(input.rows() == 8 && input.cols() == 8);

    itpp::Vec< kind > output(64);

    for (int i = 0; i < 64; ++i){
      output[i] = input(ZIGZAG_TABLE[i]/8, ZIGZAG_TABLE[i]%8);
    } // for i

    return output;
  }
  
  template<typename kind>
  std::vector< kind > InverseZigzag(const std::vector< kind > &input)
  {
    assert(input.size() == 64);
  
    std::vector< kind > output(64);

    for(int i = 0; i < 64; i++)
      {
        output[ ZIGZAG_TABLE[i] ] = input[i];
      }
  
    return output;
  }

  template<typename kind>
  itpp::Vec< kind > InverseZigzag(const itpp::Vec< kind > &input)
  {
    assert(input.size() == 64);
  
    itpp::Vec< kind > output(64);

    for(int i = 0; i < 64; i++)
      {
        output[ ZIGZAG_TABLE[i] ] = input[i];
      }
  
    return output;
  }

  template < typename kind >
  itpp::Mat< kind > InverseZigzagToMat(const itpp::Vec< kind >& input)
  {
    assert(input.size() == 64);

    itpp::Mat< kind > output(8,8);

    for (int i = 0; i < 64; ++i){
      output(ZIGZAG_TABLE[i]/8, ZIGZAG_TABLE[i]%8) = input[i];
    } // for i

    return output;
  }

  template < typename kind >
  itpp::Mat< u_char > ToUchar(const itpp::Mat< kind >& input)
  {
    itpp::Mat< u_char > output(input.rows(), input.cols());

    for (int row = 0; row < input.rows(); ++row){
      for (int col = 0; col < input.cols(); ++col){
        kind temp = input(row,col);
        if (temp > 255){
          temp = 255;
        } // if
        if (temp < 0){
          temp = 0;
        } // if
        output(row, col) = static_cast< u_char >(temp);
      } // for col
    } // for row
    
    return output;
  }
  
  enum eDctElement {
    DC_ELEMENT = 0,
    AC_ELEMENT = 1
  };

  // referenced form example.c in jpeg-9
  struct my_jpeg_error_mgr {
    struct jpeg_error_mgr pub;	/* "public" fields */

    jmp_buf setjmp_buffer;	/* for return to caller */
  };

  typedef struct my_jpeg_error_mgr * my_jpeg_error_ptr;

  METHODDEF(void)
  my_error_exit (j_common_ptr cinfo)
  {
    /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
    my_jpeg_error_ptr myerr = (my_jpeg_error_ptr) cinfo->err;

    /* Always display the message. */
    /* We could postpone this until after returning, if we chose. */
    (*cinfo->err->output_message) (cinfo);

    /* Return control to the setjmp point */
    longjmp(myerr->setjmp_buffer, 1);
  }
  
  inline int Components(J_COLOR_SPACE colorSpace)
  {
    if (colorSpace == JCS_GRAYSCALE)
      {
        return 1;
      }
    else
      {
        return 3;
      }
  }

  // ++++++++++ JpegReader ++++++++++
  // 仮想クラス
  class JpegReader
  {
  private:
    int width_, height_;
  protected:
    virtual void SetWidth(int width) { width_ = width; }
    virtual void SetHeight(int height) { height_ = height; }
    
  public:
    JpegReader(): width_(0), height_(0) { }
    virtual ~JpegReader() { }

    int Width() const { return width_; } // public非仮想
    int Height() const { return height_; }
  };
  // ---------- JpegReader ----------

  // JPEGからbitmap形式に変換
  // C++ Coding Standardsのガイドラインによりメンバ変数はprivateにした
  class Jpeg2Bitmap: public JpegReader
  {
  private:
    itpp::Mat< JSAMPLE > pixels_;
    
  public:
    Jpeg2Bitmap(): pixels_(0,0)
    { }

    Jpeg2Bitmap(const char *jpegFileName, J_COLOR_SPACE colorSpace = JCS_RGB)
    {
      Open(jpegFileName, colorSpace);
    }
    
    virtual ~Jpeg2Bitmap()
    {
      Close();
    }

    // 出てくる並びはRGBか?
    itpp::Mat< JSAMPLE > Get()
    {
      return pixels_;
    }
  
    void Close()
    {
      pixels_.clear();
    }

    void Open(const char *jpegFileName, J_COLOR_SPACE colorSpace = JCS_RGB);
    
  };

  // JPEGからDCT係数を取り出す
  class Jpeg2Dct: public JpegReader
  {
  private:
    jpeg_decompress_struct cinfo_;
    FILE *infile_;
    std::vector< int > widthInBlocks_, heightInBlocks_;
    std::vector< itpp::Mat< itpp::Mat< JCOEF > > > dctcoef_; // [componet](h_block, r_block)(h, r)
    std::vector< itpp::Mat< int > > quantizeTable_; // [component](h, r)
    std::vector< JHUFF_TBL > dcHuffmanTable_;    
    std::vector< JHUFF_TBL > acHuffmanTable_;
    int restarInterval_;
    
  public:
    // Jpeg2Dct(): widthInBlocks_(0), heightInBlocks_(0),
    //             dctcoef_(0), quantizeTable_(0), dcHuffmanTable_(0), acHuffmanTable_(0), restarInterval(0)
    // { }

    Jpeg2Dct(const char *jpegFileName, J_COLOR_SPACE colorSpace = JCS_RGB)
    {
      Open(jpegFileName, colorSpace);
    }
    
    virtual ~Jpeg2Dct()
    {
      Close();
    }

    int WidthInBlocks(int component)
    {
      return widthInBlocks_[component];
    }

    int HeightInBlocks(int component)
    {
      return heightInBlocks_[component];
    }
  
    itpp::Mat< JCOEF > Get(int component, int rowBlocks, int colBlocks)
    {
      return dctcoef_[component](rowBlocks, colBlocks);
    }

    itpp::Mat< int > QuantTable(int component)
    {
      return quantizeTable_[component];
    }

    JHUFF_TBL DcHuffmanTable(int component)
    {
      return dcHuffmanTable_[component];
    }

    JHUFF_TBL AcHuffmanTable(int component)
    {
      return acHuffmanTable_[component];
    }

    int RestartInterval()
    {
      return restarInterval_;
    }
    
    void Close()
    {
      (void) jpeg_finish_decompress(&cinfo_);
      jpeg_destroy_decompress(&cinfo_);
      fclose(infile_);
      dctcoef_.clear();
    }

    void Open(const char* jpegFileName, J_COLOR_SPACE colorSpace = JCS_RGB);
    
    j_decompress_ptr Cinfo()
    {
      return &cinfo_;
    }
    
  };

  /************************************************************************************
   * Dct2Jpeg 
   * 
   * jpeglib.hを使ってrawデータからJPEGデータを作る
   ************************************************************************************/
  class Dct2Jpeg
  {
  private:
    jpeg_compress_struct cinfo_;
    my_jpeg_error_mgr      jerr;
    FILE *outfile_;
    jvirt_barray_ptr* coeff_arrays_;
    JBLOCKARRAY coef_buffers_[MAX_COMPONENTS];
    
  public:
    // Dct2Jpeg(const char* jpegFileName, jpeg_compress_struct& cinfo)
    // {
    //   Init(jpegFileName, cinfo);
    // }
    
    Dct2Jpeg(const char* jpegFileName, j_decompress_ptr srcinfo)
    {
      Init(jpegFileName, srcinfo);
    }
    virtual ~Dct2Jpeg()
    {
      Close();
    }

    // void Init(const char* jpegFileName, jpeg_compress_struct& cinfo)
    // {
    //   cinfo_ = cinfo;
    //   Init(jpegFileName);
    // }

    void Init(const char* jpegFileName, j_decompress_ptr srcinfo);
    
    // DCT係数を保存していく
    void Make(std::vector< mylib::Vector_2D< itpp::Mat< JCOEF > > >& dctcoef); // dctcoef[component](blockY, blockX)(y, x)
    
    void Close()
    {
      fclose(outfile_);
      jpeg_destroy_compress(&cinfo_); 
    }
  };
  
}


#endif
