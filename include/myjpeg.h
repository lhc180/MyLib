/***************************************
 * This file defines my original JPEG decoder class using libjpeg.
 * 
 * Last Updated: <2013/05/27 14:00:12 from Yoshitos-iMac.local by yoshito>
 ***************************************/
// To do list
// -- 行数を指定してデータを読み込んでいく関数を作る.
// jpeg_read_scanlines()を使う
// -- テキストのエンコーダでバイナリデータを作
// らないといけない
// -- JpegFileのclose()は最後かも
// -- Huffman符号データをitpp::binで取り出す
// -- Get Huffman codewords as they are every times appearing RST marker.
// -- Huffman Tableの取得はJpeg2Dctのメンバ関数で行う
// -- もしくはJpeg2Dctを継承させて新しいクラスを作るかJpeg2Dctをメンバに持ったクラスを作る
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
#include "myhuffman.h"

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
  inline std::vector< kind > ForwardZigzag(const std::vector< kind > &input)
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
  inline std::vector< kind > InverseZigzag(const std::vector< kind > &input)
  {
    assert(input.size() == 64);
  
    std::vector< kind > output(64);

    for(int i = 0; i < 64; i++)
      {
        output[ ZIGZAG_TABLE[i] ] = input[i];
      }
  
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

  // JPEGからbitmap形式に変換
  // C++ Coding Standardsのガイドラインによりメンバ変数はprivateにした
  class Jpeg2Bitmap: public JpegReader
  {
  private:
    mylib::Vector_2D< JSAMPLE > pixels_;
    
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

    mylib::Vector_2D< JSAMPLE > Get()
    {
      return pixels_;
    }
  
    void Close()
    {
      pixels_.clear();
    }

    void Open(const char *jpegFileName, J_COLOR_SPACE colorSpace = JCS_RGB)
    {
      FILE *infile;
      if ((infile = fopen(jpegFileName, "rb")) == NULL) {
        fprintf(stderr, "can't open %s\n", jpegFileName);
        exit(1);
      }


      jpeg_decompress_struct cinfo;
      my_jpeg_error_mgr      jerr;
  
      /* Step 1: allocate and initialize JPEG decompression object */

      /* We set up the normal JPEG error routines, then override error_exit. */
      cinfo.err           = jpeg_std_error(&jerr.pub);
      jerr.pub.error_exit = my_error_exit;
      /* Establish the setjmp return context for my_error_exit to use. */
      if (setjmp(jerr.setjmp_buffer)) {
        /* If we get here, the JPEG code has signaled an error.
         * We need to clean up the JPEG object, close the input file, and return.
         */
        jpeg_destroy_decompress(&cinfo);
        fclose(infile);
        exit(1);
      }
      /* Now we can initialize the JPEG decompression object. */
      jpeg_create_decompress(&cinfo);

      /* Step 2: specify data source (eg, a file) */

      jpeg_stdio_src(&cinfo, infile);

      /* Step 3: read file parameters with jpeg_read_header() */

      (void) jpeg_read_header(&cinfo, TRUE);
      /* We can ignore the return value from jpeg_read_header since
       *   (a) suspension is not possible with the stdio data source, and
       *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
       * See libjpeg.txt for more info.
       */

      /* Step 4: set parameters for decompression */

      /* In this example, we don't need to change any of the defaults set by
       * jpeg_read_header(), so we do nothing here.
       */

      if (colorSpace != JCS_UNKNOWN)
        {
          cinfo.out_color_space = colorSpace;
        }
    
      /* Step 5: Start decompressor */

      (void) jpeg_start_decompress(&cinfo);
      /* We can ignore the return value since suspension is not possible
       * with the stdio data source.
       */

      /* We may need to do some setup of our own at this point before reading
       * the data.  After jpeg_start_decompress() we have the correct scaled
       * output image dimensions available, as well as the output colormap
       * if we asked for color quantization.
       * In this example, we need to make an output work buffer of the right size.
       */ 
      /* JSAMPLEs per row in output buffer */
      SetWidth(cinfo.output_width);
      SetHeight(cinfo.output_height);
        
      int row_stride   = cinfo.output_width * cinfo.output_components;
      int numComponents = cinfo.output_components;

      JSAMPARRAY buffer = (*cinfo.mem->alloc_sarray)
        ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

      /* Step 6: while (scan lines remain to be read) */
      /*           jpeg_read_scanlines(...); */

      /* Here we use the library's state variable cinfo.output_scanline as the
       * loop counter, so that we don't have to keep track ourselves.
       */

      std::cout << "## height = " << Height() << " width = " << Width() << std::endl;
      std::cout << "## numComponents = " << numComponents << std::endl;
   
    
      pixels_.set_size(cinfo.output_height, row_stride);

    
      std::cout << "## Debug." << std::endl;

      int currrentRow = 0;
      while (cinfo.output_scanline < cinfo.output_height) {
        /* jpeg_read_scanlines expects an array of pointers to scanlines.
         * Here the array is only one element long, but you could ask for
         * more than one scanline at a time if that's more convenient.
         */
        int readRows = jpeg_read_scanlines(&cinfo, buffer, 1);
        /* Assume put_scanline_someplace wants a pointer and sample count. */

        // std::cout << "## readRow = " << readRows << "\tcinfo.output_scanline = "<< cinfo.output_scanline << std::endl;
      
        for (int w = 0; w < static_cast< int >(cinfo.output_width); ++w)
          {
            for (int r = 0; r < readRows; ++r)
              {
                for (int n = 0; n < numComponents; ++n)
                  {
                    pixels_(currrentRow + r, numComponents*w + n) = buffer[r][numComponents*w + n];
                  } // for n
              }     // for w
          }         // for r
        currrentRow += readRows;
      } // while cinfo.output_scanline

      /* Step 7: Finish decompression */

      (void) jpeg_finish_decompress(&cinfo);
      /* We can ignore the return value since suspension is not possible
       * with the stdio data source.
       */

      /* Step 8: Release JPEG decompression object */

      /* This is an important step since it will release a good deal of memory. */
      jpeg_destroy_decompress(&cinfo);

      /* After finish_decompress, we can close the input file.
       * Here we postpone it until after no more JPEG errors are possible,
       * so as to simplify the setjmp error logic above.  (Actually, I don't
       * think that jpeg_destroy can do an error exit, but why assume anything...)
       */
      fclose(infile);

    
      /* At this point you may want to check to see whether any corrupt-data
       * warnings occurred (test whether jerr.pub.num_warnings is nonzero).
       */

      /* And we're done! */

    }
    
  };

  // JPEGからDCT係数を取り出す
  class Jpeg2Dct: public JpegReader
  {
  private:
    std::vector< int > widthInBlocks_, heightInBlocks_;
    std::vector< mylib::Vector_2D< mylib::Vector_2D< JCOEF > > > dctcoef_; // [componet](h_block, r_block)(h, r)
    std::vector< mylib::Vector_2D< u_int > > quantizeTable_; // [component](h, r)
    std::vector< JHUFF_TBL > dcHuffmanTable_;    
    std::vector< JHUFF_TBL > acHuffmanTable_;
    
  public:
    Jpeg2Dct(): widthInBlocks_(0), heightInBlocks_(0),
                dctcoef_(0), quantizeTable_(0), dcHuffmanTable_(0), acHuffmanTable_(0)
    { }

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
  
    mylib::Vector_2D< JCOEF > Get(int component, int rowBlocks, int colBlocks)
    {
      return dctcoef_[component](rowBlocks, colBlocks);
    }

    mylib::Vector_2D< u_int > QuantTable(int component)
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
    
    void Close()
    {
      dctcoef_.clear();
    }

    void Open(const char* jpegFileName, J_COLOR_SPACE colorSpace = JCS_RGB)
    {
      FILE *infile;
      if ((infile = fopen(jpegFileName, "rb")) == NULL) {
        fprintf(stderr, "can't open %s\n", jpegFileName);
        exit(1);
      }

      jpeg_decompress_struct cinfo;
      my_jpeg_error_mgr      jerr;

      cinfo.err           = jpeg_std_error(&jerr.pub);
      jerr.pub.error_exit = my_error_exit;

      if (setjmp(jerr.setjmp_buffer)) {
        jpeg_destroy_decompress(&cinfo);
        fclose(infile);
        exit(1);
      }

      jpeg_create_decompress(&cinfo);
      
      jpeg_stdio_src(&cinfo, infile);
      
      (void) jpeg_read_header(&cinfo, TRUE);

      int numComponents = 3;
      if (colorSpace == JCS_GRAYSCALE)
        {
          cinfo.out_color_space = JCS_GRAYSCALE;
          numComponents          = 1;
        }
      else if (colorSpace != JCS_UNKNOWN)
        {
          cinfo.out_color_space = colorSpace;
        }
    
      jvirt_barray_ptr* coeff_arrays = jpeg_read_coefficients(&cinfo);
    
      SetWidth(cinfo.output_width);      
      SetHeight(cinfo.output_height);
        
      std::cout << "## height = " << Height() << " width = " << Width() << std::endl;
      std::cout << "## numComponents = " << numComponents << std::endl;


      dctcoef_.resize(numComponents);
      heightInBlocks_.resize(numComponents);
      widthInBlocks_.resize(numComponents);

      quantizeTable_.resize(numComponents);
      dcHuffmanTable_.resize(numComponents);
      acHuffmanTable_.resize(numComponents);
      
      jpeg_component_info* compptr;
      JHUFF_TBL* dc_huff_tbl_ptr,* ac_huff_tbl_ptr;
      
      for (int component = 0; component < numComponents; ++component){
        compptr = cinfo.comp_info + component;
        std::cout << "## Debug." << std::endl;

        quantizeTable_[component].set_size(DCTSIZE, DCTSIZE);

        for (int i = 0; i < DCTSIZE2; ++i){
          quantizeTable_[component](i/DCTSIZE, i%DCTSIZE) = compptr->quant_table->quantval[i];
        } // for i

        dc_huff_tbl_ptr = cinfo.dc_huff_tbl_ptrs[component];
        ac_huff_tbl_ptr = cinfo.ac_huff_tbl_ptrs[component];
        assert(dc_huff_tbl_ptr != NULL && ac_huff_tbl_ptr != NULL);
        dcHuffmanTable_[component] = *dc_huff_tbl_ptr;
        acHuffmanTable_[component] = *ac_huff_tbl_ptr;

        heightInBlocks_[component] = compptr->height_in_blocks;
        widthInBlocks_[component]  = compptr->width_in_blocks;

        dctcoef_[component].set_size(heightInBlocks_[component], widthInBlocks_[component]);
                
        std::cout << "## height_in_blocks = " << compptr->height_in_blocks
                  << " width_in_blocks = " << compptr->width_in_blocks << std::endl;

        for (int blockY = 0; blockY < heightInBlocks_[component]; ++blockY){
          JBLOCKARRAY buffer = (cinfo.mem->access_virt_barray)((j_common_ptr)&cinfo, coeff_arrays[component], blockY, JDIMENSION(1), FALSE);

          for (int blockX = 0; blockX < widthInBlocks_[component]; ++blockX){
            dctcoef_[component](blockY, blockX).set_size(DCTSIZE, DCTSIZE);
            JCOEFPTR blockptr = buffer[0][blockX];

            for (int i = 0; i < DCTSIZE2; ++i){
              dctcoef_[component](blockY, blockX)(i/DCTSIZE, i%DCTSIZE) = blockptr[i];
            } // for i
          }   // for blockX
        }     // for blockY
      }       // for component
    
      (void) jpeg_finish_decompress(&cinfo);

      jpeg_destroy_decompress(&cinfo);

      fclose(infile);
    
    }
    
  };

  /************************************************************************************
   * JpegHuffman
   * 
   * jpeglib.hのJHUFF_TBLを使ってハフマン符号器を作る
   ************************************************************************************/
  class JpegHuffman
  {
  private:
    int numElement_;
    std::vector< u_char > sizeTable_;
    std::vector< u_int > codeTable_;
    std::vector< u_char > valueTable_;
    
  public:
    JpegHuffman(): numElement_(0), sizeTable_(0), codeTable_(0), valueTable_(0){ }
    explicit JpegHuffman(const JHUFF_TBL& jhuff_tbl)
    {
      Set(jhuff_tbl);
    }
    virtual ~JpegHuffman(){ }

    void Set(const JHUFF_TBL& jhuff_tbl)
    {
      assert(&jhuff_tbl != NULL);

      int num = 0;
      for (int i = 1; i <= 16; ++i){
        std::cout << "## bits[" << i << "] = " << static_cast< int >(jhuff_tbl.bits[i]) << std::endl;
        num += jhuff_tbl.bits[i];
      } // for i

      numElement_ = num;
      sizeTable_.resize(num);
      codeTable_.resize(num);
      valueTable_.resize(num);

      // サイズテーブルの生成
      for (int i = 1, k = 0; i <= 16; ++i){
        int j = 1;
        while( j <= jhuff_tbl.bits[i]){
          assert(k < num);
          sizeTable_[k] = i;
          ++k;
          ++j;
        } // while j
      } // for i
      
      int k = 0;
      int code = 0;
      int si = sizeTable_[0];
      while (1){
        while (sizeTable_[k] == si){
          codeTable_[k] = code;
          k++;
          code++;
        } // while sizeTable_[k]

        if (k >= num){
          break;
        } // if k

        do {
          code <<= 1;         // 符号語長を1ビット増やす
          si++;               // サイズを1ビット増やす
        } while (sizeTable_[k] != si);
      } // while 1

      for (int i = 0; i < num; ++i){
        valueTable_[i] = jhuff_tbl.huffval[i];
      } // for i 
    }

    // エンコーダ
    itpp::bvec Encode(int input) const
    {
      std::cout << "## numElement_ = " << numElement_ << std::endl;
      assert(input >= 0 && input < numElement_);
      itpp::bvec output = itpp::dec2bin(static_cast< int >(sizeTable_[input]), static_cast< int >(codeTable_[input]));
      return output;
    }

    itpp::bvec Encode(const std::vector< int >& input) const
    {
      itpp::bvec output(0);
      for (std::vector< int >::const_iterator ite = input.begin(); ite != input.end(); ++ite){
        output = itpp::concat(output, Encode(*ite));
      } // for ite
      return output;
    }

    int Decode1Code(const itpp::bvec& input, itpp::bvec* output)
    {
      u_int code = 0;             // ハフマン符号の候補
      u_char length = 0;          // ハフマン符号候補のビット数
      int k = 0;                  // 表のインデックス

      while (k < numElement_ && length <= 16 && length < input.size()){
        length++;
        code <<= 1;
        itpp::bvec temp = input.mid(0, length);
        code = itpp::bin2dec(temp);

        while (sizeTable_[k] == length){
          if (codeTable_[k] == code){
            *output = input.get(length, -1);
            return valueTable_[k];
          } // if code
          k++;
        } // while size
      } // while k
      
      *output = itpp::bvec(0);
      return -1;
    }

    std::vector< int > Decode(const itpp::bvec& input)
    {
      std::vector< int > decoded(0);
      itpp::bvec t_input = input;
      itpp::bvec t_output;
      int res = 0;
      while ((res = Decode1Code(t_input, &t_output)) >= 0){
        decoded.push_back(res);
        t_input = t_output;
      } // while res

      return decoded;
    }

    // friend class JpegEntropy;
  };

  /************************************************************************************
   * JpegEntropy 
   * 
   * ハフマン符号とランレングス符号を組み合わせたJPEGのエントロピー符号化
   * コンポーネントの数だけ用意しなければいけない
   * つまりグレイスケールなら1つ、RGBなどなら3つ
   ************************************************************************************/
  class JpegEntropy
  {
  private:
    JpegHuffman dcHuffman_;
    JpegHuffman acHuffman_;
    
    itpp::bvec DcEncode(int input)
    {
      int absInput = abs(input); // 絶対値

      // まずビット数を求める
      int bits = 0;
      while (absInput > 0){
        absInput >>= 1;
        bits++;
      } // while absInput
      
      itpp::bvec huffman = dcHuffman_.Encode(bits);
      
      itpp::bvec value(0);
      if (bits != 0){
        value = itpp::dec2bin(bits, absInput);
        if (input < 0){         // 負の数の場合ビットを反転させる(テキストとはやり方違う)
          value += 1;
        } // if input
      } // if bits

      itpp::bvec output = itpp::concat(huffman, value);
      
      return output;
    }

    itpp::bvec AcEncode(int input, int* run)
    {
      itpp::bvec output(0);

      int absInput = abs(input);
      while (*run > 15){
        itpp::bvec temp = acHuffman_.Encode(ZRL_INDEX);
        output = itpp::concat(output, temp);
        *run -= 16;
      } // while *run

      // ビット数を求める
      int bits = 0;
      while (absInput > 0){
        absInput >>= 1;
        bits++;
      } // while absInput
      int index = *run * 10 + bits + (*run == 15);
      itpp::bvec huffman = acHuffman_.Encode(index);

      output = itpp::concat(output, huffman);

      itpp::bvec value = itpp::dec2bin(bits, absInput);
      if (input < 0){           // 負の場合はビット反転
        value += itpp::bin(1);
      } // if input
      output = itpp::concat(output, value);

      return output;
    }

    /************************************************************************************
     * DcDecode -- DC成分のエントロピー復号
     * 
     * Arguments:
     *   input -- 入力バイナリデータ
     *   output -- 復号に用いたものを除いた残りのバイナリデータ
     *   result -- 復号結果
     *
     * Return Value:
     *   bool -- 復号成功ならtrue
     ************************************************************************************/    
    bool DcDecode(const itpp::bvec& input, itpp::bvec* output, int* result)
    {
      itpp::bvec t_output(0);
      int category = dcHuffman_.Decode1Code(input, &t_output); // 差分値のビット数
      int diff = 0;                                            // DC成分の差分値
      if (category >= 0){
        itpp::bvec temp = t_output.mid(0, category);
        if (temp[0] == 0){      // 差分が負の場合はビット反転してある
          temp += itpp::bin(1);
          diff -= itpp::bin2dec(temp);
        } // if temp[0]
        else{
          diff = itpp::bin2dec(temp);
        } // else

        *output = t_output.get(category, -1);
        *result = diff;
        return true;
      } // if category
      else{                     // categoryが負の場合
        *output = t_output;
        *result = 0;
        return false;
      }
    }

    bool AcDecode(const itpp::bvec& input, itpp::bvec* output, std::vector< int >* results)
    {
      results->resize(IMG_DCT_SIZE2-1, 0);
      int k = 1;
      itpp::bvec t_input = input;
      itpp::bvec t_output(0);
      while (k < IMG_DCT_SIZE2){
        int category = acHuffman_.Decode1Code(t_input, &t_output);
        if (category == 0){
          while (k < IMG_DCT_SIZE2){
            (*results)[k] = 0;
            k++;
          } // while k
          break;
        } // if category
        else if(category < 0){
          *output = t_output;
          return false;
        }

        int run = category >> 4;
        category &= 0x0f;

        int acValue = 0;
        if (category){
          itpp::bvec temp = t_output.mid(0, category);
          if (temp[0] == 0){
            temp += itpp::bin(1);
            acValue -= itpp::bin2dec(temp);
          } // if temp[0]
          else{
            acValue = itpp::bin2dec(temp);
          } // else
          
          t_output = t_output.get(category, -1);

        } // if category
        else if(run != 15){     // EOBでもZRLでもない
          *output = t_output;
          return false;
        }

        if (run + k > IMG_DCT_SIZE2 - 1){ // 係数が多すぎる
          *output = t_output;
          return false;
        } // if run + k

        while (run > 0){
          (*results)[k] = 0;
          k++;
          run--;
        } // while run

        (*results)[k] = acValue; // ランレングスの後に数値
        
        t_input = t_output;
      } // while k

      *output = t_output;
      return true;
    }
    
  public:
    // コンストラクタは2種類
    JpegEntropy(const JpegHuffman& dcHuffman, const JpegHuffman& acHuffman): dcHuffman_(dcHuffman), acHuffman_(acHuffman)
    { }
    JpegEntropy(const JHUFF_TBL& dcHuffTable, const JHUFF_TBL& acHuffTable):
      dcHuffman_(dcHuffTable), acHuffman_(acHuffTable)
    { }
    
    // 8*8の1ブロックのみの符号化
    itpp::bvec Encode(const std::vector< JCOEF >& input)
    {
      assert(static_cast< int >(input.size()) == IMG_DCT_SIZE2);
      // DC成分
      itpp::bvec dcCode = DcEncode(input[0]);
      
      // AC成分
      itpp::bvec acCode(0);
      int run = 0;
      for (int n = 1; n < IMG_DCT_SIZE2; ++n){
        // 係数が0でなければ
        if (input[n] != 0){
          acCode = itpp::concat(acCode, AcEncode(input[n], &run));
          run = 0;
        } // if input[n]
        else{
          if (n == IMG_DCT_SIZE2-1){
            itpp::bvec temp = acHuffman_.Encode(EOB_INDEX);
            acCode = itpp::concat(acCode, temp);
          } // if n
          else{
            run++;
          } // else n
          
        } // else input[n]
        
      } // for n
      itpp::bvec output = itpp::concat(dcCode, acCode);
      return output;
    }

    /************************************************************************************
     * Decode -- エントロピー復号
     * もし途中で復号エラーが起きたら、それ以降は0で埋めるという処理を行う
     * 
     * Arguments:
     *   input -- 入力バイナリデータ
     *
     * Return Value:
     *   std::vector< int > -- 出力データ系列
     ************************************************************************************/
    std::vector< int > Decode(const itpp::bvec& input)
    {
      itpp::bvec t_output;
      int dcCoeff;
      if (DcDecode(input, &t_output, &dcCoeff) == false){
        return std::vector< int >(IMG_DCT_SIZE2, 0);
      }
      
      itpp::bvec t_input = t_output;
      std::vector< int > acCoeff;
      if (AcDecode(t_input, &t_output, &acCoeff) == false){ // 一応エラーが出たらチェックする
        std::cerr << "Error has occured in AcDecode."  << std::endl;
      }
      std::vector< int > output = mylib::Concat(dcCoeff, acCoeff);
      return output;
    }
    
  };

  
}


#endif
