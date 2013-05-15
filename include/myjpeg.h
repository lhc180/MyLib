/***************************************
 * This file defines my original JPEG decoder class using libjpeg.
 * 
 * Last Updated: <2013/04/19 16:08:31 from Yoshitos-iMac.local by yoshito>
 ***************************************/
// To do list
// -- 行数を指定してデータを読み込んでいく関数を作る.
// jpeg_read_scanlines()を使う
// -- テキストのエンコーダでバイナリデータを作
// らないといけない
// -- JpegFileのclose()は最後かも
// -- Huffman符号データをitpp::binで取り出す
// -- Get Huffman codewords as they are every times appearing RST marker.

#ifndef MYJPEG_H
#define MYJPEG_H

#include <iostream>
#include <jpeglib.h>
#include <csetjmp>
#include <cstdio>
#include <cassert>
#include <mymatrix.h>

namespace mylib{

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

  // JPEGからbitmap形式に変換
  // C++ Coding Standardsのガイドラインによりメンバ変数はprivateにした
  class Jpeg2Bitmap// : public cJpegReader
  {
  private:
    int width_, height_;
    mylib::Vector_2D< JSAMPLE > pixels_;
  
  public:
    Jpeg2Bitmap(): width_(0), height_(0),
                   pixels_(0,0)
    { }

    Jpeg2Bitmap(const char *jpegFileName, J_COLOR_SPACE colorSpace = JCS_RGB)
    {
      Open(jpegFileName, colorSpace);
    }
    
    virtual ~Jpeg2Bitmap()
    {
      Close();
    }

    int Width()
    {
      return width_;
    }

    int Height()
    {
      return height_;
    }

    mylib::Vector_2D< JSAMPLE > Get()
    {
      return pixels_;
    }
  
    void Close()
    {
      pixels_.clear();
    }

    //

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
      width_  = cinfo.output_width;
      height_ = cinfo.output_height;
        
      int row_stride   = cinfo.output_width * cinfo.output_components;
      int numComponets = cinfo.output_components;

      JSAMPARRAY buffer = (*cinfo.mem->alloc_sarray)
        ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

      /* Step 6: while (scan lines remain to be read) */
      /*           jpeg_read_scanlines(...); */

      /* Here we use the library's state variable cinfo.output_scanline as the
       * loop counter, so that we don't have to keep track ourselves.
       */

      std::cout << "## height = " << height_ << " width = " << width_ << std::endl;
      std::cout << "## numComponets = " << numComponets << std::endl;
   
    
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
                for (int n = 0; n < numComponets; ++n)
                  {
                    pixels_(currrentRow + r, numComponets*w + n) = buffer[r][numComponets*w + n];
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
  class Jpeg2Dct// : public cJpegReader
  {
  private:
    int width_, height_;
    std::vector< int > widthInBlocks_, heightInBlocks_;
    std::vector< mylib::Vector_2D< mylib::Vector_2D< JCOEF > > > dctcoef_; // [componet](h_block, r_block)(h, r)
    std::vector< mylib::Vector_2D< u_int > > quantizeTable_; // [component](h, r)
  
  public:
    Jpeg2Dct(): width_(0), height_(0), widthInBlocks_(0), heightInBlocks_(0),
                dctcoef_(0), quantizeTable_(0)
    { }

    Jpeg2Dct(const char *jpegFileName, J_COLOR_SPACE colorSpace = JCS_RGB)
    {
      Open(jpegFileName, colorSpace);
    }
    
    virtual ~Jpeg2Dct()
    {
      Close();
    }

    int Width()
    {
      return width_;
    }

    int Height()
    {
      return height_;
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
    
      width_  = cinfo.output_width;
      height_ = cinfo.output_height;
        
      std::cout << "## height = " << height_ << " width = " << width_ << std::endl;
      std::cout << "## numComponents = " << numComponents << std::endl;


      dctcoef_.resize(numComponents);
      heightInBlocks_.resize(numComponents);
      widthInBlocks_.resize(numComponents);
    
      quantizeTable_.resize(numComponents);
    
    
      for (int component = 0; component < numComponents; ++component)
        {

          jpeg_component_info* compptr = cinfo.comp_info + component;
          std::cout << "## Debug." << std::endl;

          quantizeTable_[component].set_size(DCTSIZE, DCTSIZE);

          for (int i = 0; i < DCTSIZE2; ++i)
            {
              quantizeTable_[component](i/DCTSIZE, i%DCTSIZE) = compptr->quant_table->quantval[i];
            } // for i
        
          heightInBlocks_[component] = compptr->height_in_blocks;
          widthInBlocks_[component]  = compptr->width_in_blocks;

          dctcoef_[component].set_size(heightInBlocks_[component], widthInBlocks_[component]);
                
          std::cout << "## height_in_blocks = " << compptr->height_in_blocks
                    << " width_in_blocks = " << compptr->width_in_blocks << std::endl;

          for (int blockY = 0; blockY < heightInBlocks_[component]; ++blockY)
            {
              JBLOCKARRAY buffer = (cinfo.mem->access_virt_barray)((j_common_ptr)&cinfo, coeff_arrays[component], blockY, JDIMENSION(1), FALSE);

              for (int blockX = 0; blockX < widthInBlocks_[component]; ++blockX)
                {
                  dctcoef_[component](blockY, blockX).set_size(DCTSIZE, DCTSIZE);
                  JCOEFPTR blockptr = buffer[0][blockX];
                  for (int i = 0; i < DCTSIZE2; ++i)
                    {
                      dctcoef_[component](blockY, blockX)(i/DCTSIZE, i%DCTSIZE) = blockptr[i];
                    }             // for i
                }                 // for blockX
            }                     // for blockY
        }                         // for component
    
      (void) jpeg_finish_decompress(&cinfo);

      jpeg_destroy_decompress(&cinfo);

      fclose(infile);
    
    }
  
  
  };
  
}


#endif
