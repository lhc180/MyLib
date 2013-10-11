/************************************************************************************
 * myentropy.cpp
 *   
 * myentropy.hの実装部分
 *
 * Contents:
 *   EntropyEncoderとEntropyDecoderクラスの実装
 *
 * Last Updated: <2013/10/11 18:23:58 from Yoshitos-iMac.local by yoshito>
 ************************************************************************************/

#include "../include/myentropy.h"

namespace mylib{

      // ++++++++++++++++++++ JpegHuffmanEncoder ++++++++++++++++++++

  void JpegHuffmanEncoder::Set(const JHUFF_TBL& jhuff_tbl)
  {
    assert(&jhuff_tbl != NULL);

    int num = 0;
    for (int i = 1; i <= 16; ++i){
      // std::cout << "## bits[" << i << "] = " << static_cast< int >(jhuff_tbl.bits[i]) << std::endl;
      num += jhuff_tbl.bits[i];
    } // for i

    numElement_ = num;
    // sizeTable_.resize(num);
    // codeTable_.resize(num);
      
    std::vector< u_char > t_size(num, 0);
    std::vector< u_int > t_code(num, 0);
      
    // サイズテーブルの生成
    int p = 0;
    for (int l = 1; l <= 16; ++l){
      int i = jhuff_tbl.bits[l];
      assert(i >= 0 && p + i <= 256);
        
      while(i--){
        t_size[p++] = static_cast< char >(l);
      } // while i
    } // for l

    t_size[p] = 0;
    int lastp = p;
      
    // コードテーブルの生成
      
    int code = 0;
    int si = t_size[0];
    p = 0;
    while (t_size[p]){
      while (t_size[p] == si){
        t_code[p++] = code;
        code++;
      } // while si

        // if (k >= num){
        //   break;
        // } // if k

      assert(code < (1 << si));
        
      code <<= 1;         // 符号語長を1ビット増やす
      si++;               // サイズを1ビット増やす

    } // while 1
      
    for (int i = 0; i < lastp; ++i){
      int v = jhuff_tbl.huffval[i];
      codeTable_[v] = t_code[i];
      sizeTable_[v] = t_size[i];
    } // for i 
  }

  // -------------------- JpegHuffmanEncoder --------------------

  // ++++++++++++++++++++ JpegHuffmanDecoder ++++++++++++++++++++
  void JpegHuffmanDecoder::Set(const JHUFF_TBL& jhuff_tbl)
  {
    assert(&jhuff_tbl != NULL);

    int num = 0;
    for (int i = 1; i <= 16; ++i){
      // std::cout << "## bits[" << i << "] = " << static_cast< int >(jhuff_tbl.bits[i]) << std::endl;
      num += jhuff_tbl.bits[i];
    } // for i

    numElement_ = num;
    sizeTable_.resize(num);
    codeTable_.resize(num);
    valueTable_.resize(num);

    // std::cout << "## numElement_ = " << numElement_ << std::endl;
      
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
        // std::cout << "## k = " << k << " ## code = " << code << " ## size = "<< static_cast< int >(sizeTable_[k]) << std::endl;
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
  // -------------------- JpegHuffmanDecoder --------------------

  // ++++++++++++++++++++ JpegEntropyEncoder ++++++++++++++++++++
  itpp::bvec JpegEntropyEncoder::DcEncode_(int input)
  {
    int absInput = abs(input); // 絶対値

    // まずビット数を求める
    int bits = 0;
    while (absInput > 0){
      absInput >>= 1;
      bits++;
    } // while absInput
      
    itpp::bvec huffman = dcHuffman_(bits);
    
    itpp::bvec value(0);
    if (bits != 0){
      value = itpp::dec2bin(bits, abs(input));
      if (input < 0){         // 負の数の場合ビットを反転させる(テキストとはやり方違う)
        value += 1;
      } // if input
    } // if bits

    itpp::bvec output = itpp::concat(huffman, value);
      
    return output;
  }

  itpp::bvec JpegEntropyEncoder::AcEncode_(int input, int* run)
  {
    itpp::bvec output(0);

    int absInput = abs(input);
    while (*run > 15){
      itpp::bvec temp = acHuffman_(ZRL_INDEX);
      output = itpp::concat(output, temp);
      *run -= 16;
    } // while *run

      // ビット数を求める
    int bits = 0;
    while (absInput > 0){
      absInput >>= 1;
      bits++;
    } // while absInput
      // std::cout << "## bits = " << bits << std::endl;
    int index = (*run << 4) + bits; // ## jpeglib.hに従った。
    itpp::bvec huffman = acHuffman_(index);
      
    output = itpp::concat(output, huffman);

    itpp::bvec value = itpp::dec2bin(bits, abs(input));
    if (input < 0){           // 負の場合はビット反転
      value += itpp::bin(1);
    } // if input
    output = itpp::concat(output, value);

    return output;
  }

  // -------------------- JpegEntropyEncoder --------------------

  // ++++++++++++++++++++ JpegEntropyDecoder ++++++++++++++++++++

  bool JpegEntropyDecoder::DcDecode_(const itpp::bvec& input, itpp::bvec* output, int* result)
  {
    itpp::bvec t_output(0);
    int category = dcHuffman_(input, &t_output); // 差分値のビット数
    int diff = 0;      // DC成分の差分値
    if (category > 0){
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
    else if (category == 0){
      *output = t_output;
      *result = 0;
      return true;
    } // if category
    else{                     // categoryが負の場合
      *output = t_output;
      *result = 0;
      return false;
    }
  }
  
  bool JpegEntropyDecoder::AcDecode_(const itpp::bvec& input, itpp::bvec* output, itpp::ivec* results)
  {
    results->set_size(IMG_DCT_SIZE2-1);
    results->zeros();
    
    int k = 0;
    itpp::bvec t_input = input;
    itpp::bvec t_output(0);
    while (k < IMG_DCT_SIZE2-1){
      int category = acHuffman_(t_input, &t_output);
      //        std::cout << "## category = " << category << std::endl;
      if (category == 0){
        while (k < IMG_DCT_SIZE2-1){
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
        // std::cout << "## category = " << category << std::endl;
        itpp::bvec temp = t_output.mid(0, category); // ## ここでエラー
        // categoryがt_outputの要素数超えていたらダメ
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

      if (run + k > IMG_DCT_SIZE2 - 2){ // 係数が多すぎる
        *output = t_output;
        return false;
      } // if run + k

        // std::cout << "## run = "  << run << std::endl;
        
      while (run > 0){
        (*results)[k] = 0;
        k++;
        run--;
      } // while run

        // std::cout << "## acValue = " << acValue << std::endl;
      (*results)[k] = acValue; // ランレングスの後に数値
      k++;
      t_input = t_output;
    } // while k

    *output = t_output;
    return true;
  }
  // -------------------- JpegEntropyDecoder --------------------

  
  // ++++++++++ EntropyEncoder ++++++++++
  // DC係数を並べていったものに対して符号化する
  itpp::bvec EntropyEncoder::EncodeDc(const itpp::Vec< JCOEF >& input)
  {
    itpp::bvec output(0);
    for (int i = 0; i < input.size(); ++i){
      output = itpp::concat(output, DcEncode_(input[i]));
    } // for ite
    return output;
  }
  
  itpp::bvec EntropyEncoder::EncodeAc(const itpp::Vec< JCOEF >& input)
  {
    itpp::bvec output(0);
    int run = 0;
    for (int n = 0; n < input.size(); ++n){
      // 係数が0でなければ
      if (input[n] != 0){
        itpp::bvec temp = AcEncode_(input[n], &run);
        output = itpp::concat(output, temp);
        run = 0;
      } // if input[n]
      else{
        if (n == input.size() - 1){
          itpp::bvec temp = acHuffman_(EOB_INDEX);
          output = itpp::concat(output, temp);
        } // if n
        else{
          run++;
        } // else n
          
      } // else ite
        
    } // for n
      
    return output;
  }

  itpp::bvec EntropyEncoder::EncodePart(const itpp::Vec< JCOEF >& input, int blocks, int layers)
  {
    assert(layers > 0 && layers <= 64);
    
    itpp::bvec dcEncoded;
    
    dcEncoded = EncodeDc(input);

    if (layers == 1){
      return dcEncoded;
    } // if
    else{
      itpp::bvec acEncoded = EncodeAc(input.get(blocks, -1));
      return itpp::concat(dcEncoded, acEncoded);
    } // else
  }
  
  // ---------- EntropyEncoder ----------
  
  // ++++++++++++++++ EntropyDecoder ++++++++++++++++++++
  int EntropyDecoder::AcDecode_(const itpp::bvec& input,
                                 itpp::bvec* output,
                                 itpp::Vec< int >* results,
                                 int size,
                                 int effectiveSize)
  {
    results->set_size(size);
    results->zeros();
    int inputSize = input.size();
    if (effectiveSize < 0){
      effectiveSize = inputSize;
    } // if effectiveSize

    int k = 0;
    itpp::bvec t_input = input;
    itpp::bvec t_output(0);
    itpp::bvec temp2;

    int countSize = 0;

    while (k < size){
      // std::cout << "## k = " << k << std::endl;
        
      int category = acHuffman_(t_input, &t_output);
      
      if (category == 0){
        // while (k < size){
        //   (*results)[k] = 0;
        //   ++k;
        // } // while k
        countSize = size;
        break;
      } // if category
      else if(category < 0){
        *output = t_output;
        return countSize;
      }
        
      int run = category >> 4;
      category &= 0x0f;

      int acValue = 0;
      if (category > 0){

        if (inputSize - t_output.size() - category > effectiveSize){
          countSize = k;
          break;
        } // if 
        
        if (category > t_output.size()){
          *output = t_output;
          return countSize;
        } // if 
          
        itpp::bvec temp = t_output.mid(0, category); // ## ここでエラー
        // categoryがt_outputの要素数超えていたらダメ
        if (temp[0] == 0){
          temp += itpp::bin(1);
          acValue -= itpp::bin2dec(temp);
        } // if temp[0]
        else{
          acValue = itpp::bin2dec(temp);
        } // else
        
        temp2 = t_output.get(category, -1);
        t_output = temp2;

      } // if category
      else if(run != 15){     // EOBでもZRLでもない
        *output = t_output;
        return countSize;
      }

      if (run + k > size - 1){ // 係数が多すぎる
        *output = t_output;
        return countSize;
      } // if run + k
        
      while (run-- > 0 && k < size){
        (*results)[k] = 0;
        k++;
      } // while run
        
      (*results)[k] = acValue; // ランレングスの後に数値
      k++;
      t_input = t_output;
      countSize = k;
    } // while k

    *output = t_output;
    return countSize;
  }

  itpp::Vec< int > EntropyDecoder::DecodeDc(const itpp::bvec& input, int size)
  {
    itpp::Vec< int > output(size);
    itpp::bvec t_input = input;
    int result;
      
    for (int i = 0; i < size; ++i){
      itpp::bvec t_output;
      if (DcDecode_(t_input, &t_output, &result)){
        output[i] = result;
        t_input = t_output;
      } // if
      else{                   // もしデコードエラー出たらループ抜ける
        break;
      } // else 
    } // for i
    return output;
  }

  itpp::Vec< int > EntropyDecoder::DecodeDc(const itpp::bvec& input, int size, itpp::bvec* rest)
  {
    itpp::Vec< int > output(size);
    itpp::bvec t_input = input;
    int result;
      
    for (int i = 0; i < size; ++i){
      if (DcDecode_(t_input, rest, &result)){
        output[i] = result;
        t_input = *rest;
      } // if
      else{                   // もしデコードエラー出たらループ抜ける
        break;
      } // else 
    } // for i
    return output;
  }

  
  itpp::Vec< int > EntropyDecoder::DecodeAc(const itpp::bvec& input, int size)
  {
    itpp::Vec< int > output;
    itpp::bvec t_output;
    this->AcDecode_(input, &t_output, &output, size); // ## falseなら復号失敗だけどエラー隠蔽している
    return output;
  }

  itpp::Vec< int > EntropyDecoder::DecodeAc(const itpp::bvec& input, int size, itpp::bvec* rest)
  {
    itpp::Vec< int > output;
    this->AcDecode_(input, rest, &output, size); // ## falseなら復号失敗だけどエラー隠蔽している
    return output;
  }

  itpp::Vec< int > EntropyDecoder::DecodeAll(const itpp::bvec& input, int originalSize, int effectiveSize)
  {
    itpp::Vec< int > output;
    DecodeAll(input, &output, originalSize, effectiveSize);
    
    return output;
  }

  int EntropyDecoder::DecodeAll(const itpp::bvec &input, itpp::ivec *output, int originalSize, int effectiveSize)
  {
    assert(originalSize % DCTSIZE2 == 0);
    int inputSize = input.size();
    if (effectiveSize < 0){
      effectiveSize = inputSize;
    } // if effectiveSize
    
    int blocks = originalSize / DCTSIZE2;

    output->set_size(originalSize);
    output->zeros();
    itpp::bvec t_input = input;
    // DCのデコード
    for (int i = 0; i < blocks; ++i){
      int result;
      itpp::bvec t_output;
      if (DcDecode_(t_input, &t_output, &result)){
        if (inputSize - t_output.size() > effectiveSize){ // 有効なバイナリ長を越えていたら
          return i;
        } // if
        else{
          (*output)[i] = result;
          t_input = t_output;
        } // else 
      } // if
      else{
        return i;
      } // else
    } // for i  
      
      // ACのデコード
    itpp::Vec< int > acOut;
    itpp::bvec t_output;
    // std::cout << "## originalSize - blocks = " << originalSize - blocks << std::endl;
      
    int acCount = this->AcDecode_(t_input, &t_output, &acOut, originalSize - blocks,
                    effectiveSize - (inputSize - t_input.size()));

    // std::cout << "## decoded AC" << std::endl;
      
    int b = blocks;
    for (int i = 0; i < acOut.size(); ++b, ++i){
      (*output)[b] = acOut[i];
    }

    return blocks+acCount;
  }

  // --------------- EntropyDecoder -------------------
  
}
