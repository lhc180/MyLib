/************************************************************************************
 * myentropy.cpp
 *   
 * myentropy.hの実装部分
 *
 * Contents:
 *   EntropyEncoderとEntropyDecoderクラスの実装
 *
 * Last Updated: <2013/09/23 01:02:52 from yoshitos-mac-mini.local by yoshito>
 ************************************************************************************/

#include "../include/myentropy.h"

namespace mylib{

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
  
  // ---------- EntropyEncoder ----------
  
  // ++++++++++++++++ EntropyDecoder ++++++++++++++++++++
  bool EntropyDecoder::AcDecode_(const itpp::bvec& input,
                                itpp::bvec* output,
                                itpp::Vec< int >* results,
                                int size)
  {
    results->set_size(size);
    results->zeros();
    int k = 0;
    itpp::bvec t_input = input;
    itpp::bvec t_output(0);
    itpp::bvec temp2;
    while (k < size){
      // std::cout << "## k = " << k << std::endl;
        
      int category = acHuffman_(t_input, &t_output);
      if (category == 0){
        while (k < size){
          (*results)[k] = 0;
          ++k;
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
      if (category > 0){

        if (category > t_output.size()){
          *output = t_output;
          return false;
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
        return false;
      }

      if (run + k > size - 1){ // 係数が多すぎる
        *output = t_output;
        return false;
      } // if run + k
        
      while (run-- > 0 && k < size){
        (*results)[k] = 0;
        k++;
      } // while run
        
      (*results)[k] = acValue; // ランレングスの後に数値
      k++;
      t_input = t_output;
    } // while k

    *output = t_output;
    return true;
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

  itpp::Vec< int > EntropyDecoder::DecodeAll(const itpp::bvec& input, int originalSize)
  {
    assert(originalSize % DCTSIZE2 == 0);
    int blocks = originalSize / DCTSIZE2;

    itpp::Vec< int > output(originalSize);
    output.zeros();
    itpp::bvec t_input = input;
    // DCのデコード
    for (int i = 0; i < blocks; ++i){
      int result;
      itpp::bvec t_output;
      if (DcDecode_(t_input, &t_output, &result)){
        output[i] = result;
        t_input = t_output;
      } // if
      else{
        return output;
      } // else
    } // for i  
      
      // ACのデコード
    itpp::Vec< int > acOut;
    itpp::bvec t_output;
    // std::cout << "## originalSize - blocks = " << originalSize - blocks << std::endl;
      
    this->AcDecode_(t_input, &t_output, &acOut, originalSize - blocks);

    // std::cout << "## decoded AC" << std::endl;
      
    int b = blocks;
    for (int i = 0; i < acOut.size(); ++b, ++i){
      output[b] = acOut[i];
    }
      
    return output;
  }


  // --------------- EntropyDecoder -------------------
  
}
