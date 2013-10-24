/************************************************************************************
 * myjpeg.cpp
 *   
 * myjpeg.hの実装部分
 *
 * Contents:
 *   色々
 *
 * Last Updated: <2013/10/11 17:54:05 from Yoshitos-iMac.local by yoshito>
 ************************************************************************************/

#include "../include/myjpeg.h"

namespace mylib{
  
  // ++++++++++ Jpeg2Bitmap ++++++++++
  void Jpeg2Bitmap::Open(const char *jpegFileName, J_COLOR_SPACE colorSpace)
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
  // ---------- Jpeg2Bitmap ----------

  // ++++++++++ Jpeg2Dct ++++++++++
  void Jpeg2Dct::Open(const char* jpegFileName, J_COLOR_SPACE colorSpace)
  {      
    if ((infile_ = fopen(jpegFileName, "rb")) == NULL) {
      fprintf(stderr, "can't open %s\n", jpegFileName);
      exit(1);
    }

    //      jpeg_decompress_struct cinfo;
    my_jpeg_error_mgr      jerr;

    cinfo_.err           = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;

    if (setjmp(jerr.setjmp_buffer)) {
      jpeg_destroy_decompress(&cinfo_);
      fclose(infile_);
      exit(1);
    }

    jpeg_create_decompress(&cinfo_);

    std::cout << "## numComponents = " << cinfo_.num_components << std::endl;
      
    jpeg_stdio_src(&cinfo_, infile_);
      
    (void) jpeg_read_header(&cinfo_, TRUE);

    int numComponents = 3;
    if (colorSpace == JCS_GRAYSCALE)
      {
        cinfo_.out_color_space = JCS_GRAYSCALE;
        numComponents          = 1;
      }
    else if (colorSpace != JCS_UNKNOWN)
      {
        cinfo_.out_color_space = colorSpace;
      }
    
    jvirt_barray_ptr* coeff_arrays = jpeg_read_coefficients(&cinfo_);

    restarInterval_ = cinfo_.restart_interval;
      
    SetWidth(cinfo_.output_width);      
    SetHeight(cinfo_.output_height);
        
    std::cout << "## height = " << Height() << " width = " << Width() << std::endl;
    std::cout << "## input numComponents = " << cinfo_.num_components << std::endl;

    dctcoef_.resize(numComponents);
    heightInBlocks_.resize(numComponents);
    widthInBlocks_.resize(numComponents);

    quantizeTable_.resize(numComponents);
    dcHuffmanTable_.resize(numComponents);
    acHuffmanTable_.resize(numComponents);
      
    jpeg_component_info* compptr;
    JHUFF_TBL* dc_huff_tbl_ptr,* ac_huff_tbl_ptr;
      
    for (int component = 0; component < numComponents; ++component){
      compptr = cinfo_.comp_info + component;

      quantizeTable_[component].set_size(DCTSIZE, DCTSIZE);

      for (int i = 0; i < DCTSIZE2; ++i){
        quantizeTable_[component](i/DCTSIZE, i%DCTSIZE) = compptr->quant_table->quantval[i];
      } // for i

      dc_huff_tbl_ptr = cinfo_.dc_huff_tbl_ptrs[component];
      ac_huff_tbl_ptr = cinfo_.ac_huff_tbl_ptrs[component];
      assert(dc_huff_tbl_ptr != NULL && ac_huff_tbl_ptr != NULL);
      dcHuffmanTable_[component] = *dc_huff_tbl_ptr;
      acHuffmanTable_[component] = *ac_huff_tbl_ptr;

      heightInBlocks_[component] = compptr->height_in_blocks;
      widthInBlocks_[component]  = compptr->width_in_blocks;

      dctcoef_[component].set_size(heightInBlocks_[component], widthInBlocks_[component]);
      
      std::cout << "## height_in_blocks = " << compptr->height_in_blocks
                << " width_in_blocks = " << compptr->width_in_blocks << std::endl;

      for (int blockY = 0; blockY < heightInBlocks_[component]; ++blockY){
        JBLOCKARRAY buffer = (cinfo_.mem->access_virt_barray)((j_common_ptr)&cinfo_, coeff_arrays[component], blockY, JDIMENSION(1), FALSE);
          
        for (int blockX = 0; blockX < widthInBlocks_[component]; ++blockX){
          dctcoef_[component](blockY, blockX).set_size(DCTSIZE, DCTSIZE);
          JCOEFPTR blockptr = buffer[0][blockX];

          for (int i = 0; i < DCTSIZE2; ++i){
            dctcoef_[component](blockY, blockX)(i/DCTSIZE, i%DCTSIZE) = blockptr[i];
          } // for i
        }   // for blockX
      }     // for blockY
    }       // for component      
  }
  // -------------------- Jpeg2Dct --------------------

  // ++++++++++++++++++++ Dct2Jpeg ++++++++++++++++++++
  void Dct2Jpeg::Init(const char* jpegFileName, j_decompress_ptr srcinfo)
  {

    cinfo_.err = jpeg_std_error(&jerr.pub);
    // jerr.pub.error_exit = my_error_exit;

    jpeg_create_compress(&cinfo_);
      
    // DCT係数を蓄える領域を確保
    for (int compnum=0; compnum < srcinfo->num_components; compnum++){
      coef_buffers_[compnum] = ((&cinfo_)->mem->alloc_barray) 
        ((j_common_ptr) &cinfo_, JPOOL_IMAGE,
         srcinfo->comp_info[compnum].width_in_blocks,
         srcinfo->comp_info[compnum].height_in_blocks);
    }
      
    coeff_arrays_ = jpeg_read_coefficients(srcinfo);
    jpeg_copy_critical_parameters(srcinfo, &cinfo_);
      
    if ((outfile_ = fopen(jpegFileName, "wb")) == NULL) {
      fprintf(stderr, "can't open %s\n", jpegFileName);
      exit(1);
    }
  }
  // DCT係数を保存していく
  void Dct2Jpeg::Make(std::vector< mylib::Vector_2D< itpp::Mat< JCOEF > > >& dctcoef) // dctcoef[component](blockY, blockX)(y, x)
  {
      
    int numComponents = cinfo_.num_components;
    std::cout << "## numComponents = " << numComponents << std::endl;
      
    assert(numComponents == static_cast< int >(dctcoef.size()));

      
    jpeg_component_info* compptr;
    JBLOCKARRAY row_ptrs[MAX_COMPONENTS];
      
    for (int component = 0; component < numComponents; ++component){
      compptr = cinfo_.comp_info + component;
        
      int heightInBlocks_ = compptr->height_in_blocks;
      int widthInBlocks_  = compptr->width_in_blocks;
        
      for (int blockY = 0; blockY < heightInBlocks_; ++blockY){
        row_ptrs[component] = (cinfo_.mem->access_virt_barray)((j_common_ptr)&cinfo_, coeff_arrays_[component], blockY, JDIMENSION(1), TRUE);
          
        for (int blockX = 0; blockX < widthInBlocks_; ++blockX){
          JCOEFPTR blockptr = row_ptrs[component][0][blockX];
            
          for (int i = 0; i < DCTSIZE2; ++i){
            blockptr[i] = dctcoef[component](blockY, blockX)(i/DCTSIZE, i%DCTSIZE);
          } // for i
        }   // for blockX
      }     // for blockY
    } // for comp

    jpeg_stdio_dest(&cinfo_, outfile_);

    std::cout << "## OK !!" << std::endl;
    jpeg_write_coefficients(&cinfo_, coeff_arrays_);
      
    jpeg_finish_compress(&cinfo_); // 最後
  }

  // -------------------- Dct2Jpeg --------------------


}

