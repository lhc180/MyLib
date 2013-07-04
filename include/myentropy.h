#ifndef MYENTROPY_H
#define MYENTROPY_H
/************************************************************************************
 * My Entropy Coder
 *
 * JPEGで使われているエントロピー符号化を自分の研究に使うようにしたもの。
 * 汎用ハフマンテーブルを使ったJPEGデータじゃないと上手く動かない。
 * 
 * Contents:
 *   Encode()
 *   Decode()
 *   SetTable()
 * Last Updated: <2013/07/04 15:06:17 from Yoshitos-iMac.local by yoshito>
 ************************************************************************************/
// To Do List
// -- JPEGプログレッシブのやり方を使いたい

#include <algorithm>
#include <cassert>
#include "myjpeg.h"

namespace mylib{
  
  /************************************************************************************
   * SeparateEntropyEncoder 
   * 
   * JpegEntropyEncoderのDC成分とAC成分のエンコードを別々で行えるようにしたもの
   ************************************************************************************/
  class EntropyEncoder: public JpegEntropyEncoder
  {
  public:
    EntropyEncoder(const JpegHuffmanEncoder& dcHuffman, const JpegHuffmanEncoder& acHuffman):
      JpegEntropyEncoder(dcHuffman, acHuffman)
    { }
    EntropyEncoder(const JHUFF_TBL& dcHuffmanTable, const JHUFF_TBL& acHuffmanTable):
      JpegEntropyEncoder(dcHuffmanTable, acHuffmanTable)
    { }
    virtual ~EntropyEncoder(){ }

    // DC係数を並べていったものに対して符号化する
    itpp::bvec EncodeDc(const itpp::Vec< JCOEF >& input)
    {
      itpp::bvec output(0);
      for (int i = 0; i < input.size(); ++i){
        output = itpp::concat(output, DcEncode(input[i]));
      } // for ite
      return output;
    }

    itpp::bvec EncodeAc(const itpp::Vec< JCOEF >& input)
    {
      itpp::bvec output(0);
      int run = 0;
      for (int n = 0; n < input.size(); ++n){
        // 係数が0でなければ
        if (input[n] != 0){
          itpp::bvec temp = AcEncode(input[n], &run);
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

    itpp::bvec EncodeAll(const itpp::Vec< JCOEF >& input, int blocks)
    {      
      itpp::bvec dcEncoded = EncodeDc(input.left(blocks));
      itpp::bvec acEncoded = EncodeAc(input.mid(blocks, -1));

      return itpp::concat(dcEncoded, acEncoded);
    }
  };

  /************************************************************************************
   * SeparateEntropyDecoder 
   * 
   * JpegEntropyDecoderのDC成分とAC成分を別々に行えるようにしたもの
   ************************************************************************************/
  class EntropyDecoder: public JpegEntropyDecoder
  {
  protected:
    virtual bool AcDecode(const itpp::bvec& input, itpp::bvec* output, itpp::Vec< int >* results, u_int size)
    {
      results->set_size(size);
      u_int k = 0;
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

  public:
    EntropyDecoder(const JpegHuffmanDecoder& dcHuffman, const JpegHuffmanDecoder& acHuffman):
      JpegEntropyDecoder(dcHuffman, acHuffman){ }
    EntropyDecoder(const JHUFF_TBL& dcHuffTable, const JHUFF_TBL& acHuffTable):
      JpegEntropyDecoder(dcHuffTable, acHuffTable){ }
    virtual ~EntropyDecoder(){ }

    itpp::Vec< int > DecodeDc(const itpp::bvec& input, int size)
    {
      itpp::Vec< int > output(size);
      itpp::bvec t_input = input;
      int result;
      
      for (int i = 0; i < size; ++i){
        itpp::bvec t_output;
        if (DcDecode(t_input, &t_output, &result)){
          output[i] = result;
          t_input = t_output;
        } // if
        else{                   // もしデコードエラー出たらループ抜ける
          break;
        } // else 
      } // for i
      return output;
    }

    itpp::Vec< int > DecodeAc(const itpp::bvec& input, u_int size)
    {
      itpp::Vec< int > output;
      itpp::bvec t_output;
      this->AcDecode(input, &t_output, &output, size); // ## falseなら復号失敗だけどエラー隠蔽している
      return output;
    }

    itpp::Vec< int > DecodeAll(const itpp::bvec& input, int originalSize)
    {
      assert(originalSize % DCTSIZE2 == 0);
      int blocks = originalSize / DCTSIZE2;

      itpp::Vec< int > output(originalSize);
      itpp::bvec t_input = input;

      // DCのデコード
      for (int i = 0; i < blocks; ++i){
        int result;
        itpp::bvec t_output;
        if (DcDecode(t_input, &t_output, &result)){
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
      
      this->AcDecode(t_input, &t_output, &acOut, originalSize - blocks);

      // std::cout << "## decoded AC" << std::endl;
      
      int b = blocks;
      for (int i = 0; i < acOut.size(); ++b, ++i){
        output[b] = acOut[i];
      }
      
      return output;
    }
  };
  
}

#endif
