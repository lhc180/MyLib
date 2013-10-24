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
 * Last Updated: <2013/10/11 18:27:58 from Yoshitos-iMac.local by yoshito>
 ************************************************************************************/
// To Do List
// -- JPEGプログレッシブのやり方を使いたい

#include <cassert>
#include "myjpeg.h"

namespace mylib{

  
    /************************************************************************************
   * JpegHuffmanEncoder
   * 
   * jpeglib.hのJHUFF_TBLを使ってハフマン符号器を作る
   ************************************************************************************/
  class JpegHuffmanEncoder
  {
  private:
    int numElement_;
    std::vector< u_char > sizeTable_;
    std::vector< u_int > codeTable_;
    
  public:
    //     JpegHuffmanEncoder(): numElement_(0), sizeTable_(256, 0), codeTable_(256, 0){ }
    explicit JpegHuffmanEncoder(const JHUFF_TBL& jhuff_tbl): numElement_(0), sizeTable_(256, 0), codeTable_(256, 0)
    {
      Set(jhuff_tbl);
    }
    virtual ~JpegHuffmanEncoder(){ }

    void Set(const JHUFF_TBL& jhuff_tbl);
    
    // エンコーダ
    itpp::bvec operator ()(int input) const
    {
      //  std::cout << "## numElement_ = " << numElement_ << std::endl;
      assert(input >= 0 && input < 256);
      itpp::bvec output = itpp::dec2bin(static_cast< int >(sizeTable_[input]), static_cast< int >(codeTable_[input]));
      // std::cout << "## codeTable_[input] = " << codeTable_[input] << std::endl;
      return output;
    }

    itpp::bvec operator ()(const itpp::Vec< int >& input) const
    {
      itpp::bvec output(0);
      for (int i = 0; i < input.size(); ++i){
        output = itpp::concat(output, operator()(input[i]));
      } // for ite
      return output;
    }
  };

  /************************************************************************************
   * JpegHuffmanDecoder
   * 
   * jpeglib.hのJHUFF_TBLを使ってハフマン符号器を作る
   ************************************************************************************/
  class JpegHuffmanDecoder
  {
  private:
    int numElement_;
    std::vector< u_char > sizeTable_;
    std::vector< u_int > codeTable_;
    std::vector< u_char > valueTable_;
    
  public:
    // JpegHuffmanDecoder(): numElement_(0), sizeTable_(0), codeTable_(0), valueTable_(0){ }
    explicit JpegHuffmanDecoder(const JHUFF_TBL& jhuff_tbl)
    {
      Set(jhuff_tbl);
    }
    virtual ~JpegHuffmanDecoder(){ }

    void Set(const JHUFF_TBL& jhuff_tbl);
    
    // 1語のみのデコード
    int operator ()(const itpp::bvec& input, itpp::bvec* output)
    {
      u_int code = 0;             // ハフマン符号の候補
      int length = 0;          // ハフマン符号候補のビット数
      int k = 0;                  // 表のインデックス

      // std::cout << "## input.size = " << input.size() << std::endl;
      while ((k < numElement_) && (length <= 16) && (length < input.size())){
        length++;
        // std::cout << "## length = " << static_cast< int >(length) << std::endl;
        //        code <<= 1;
        itpp::bvec temp = input.mid(0, length);
        code = itpp::bin2dec(temp);
        
        while (sizeTable_[k] == length){
          if (codeTable_[k] == code){
            // std::cout << "## code = " << code << std::endl;
            *output = input.get(length, -1);
            return valueTable_[k];
          } // if code
          k++;
        } // while size
      } // while k
      
      *output = itpp::bvec(0);
      return -1;
    }

    itpp::Vec< int > operator ()(const itpp::bvec& input)
    {
      itpp::Vec< int > decoded(0);
      itpp::bvec t_input = input;
      itpp::bvec t_output;
      int res = 0;
      while ((res = operator()(t_input, &t_output)) >= 0){
        decoded = itpp::concat(decoded, res);
        t_input = t_output;
      } // while res

      return decoded;
    }

    void PrintTable()
    {
      for (int i = 0; i < numElement_; ++i){
        std::cout << static_cast< int >(valueTable_[i]) << "\t" <<
          itpp::dec2bin(static_cast< int >(sizeTable_[i]), static_cast< int >(codeTable_[i])) << std::endl;
      } // for i
    }
    
    // friend class JpegEntropy;
  };

  /************************************************************************************
   * JpegEntropyEncoder
   * 
   * ハフマン符号とランレングス符号を組み合わせたJPEGのエントロピー符号化
   * コンポーネントの数だけ用意しなければいけない
   * つまりグレイスケールなら1つ、RGBなどなら3つ
   ************************************************************************************/
  class JpegEntropyEncoder
  {
  protected:
    JpegHuffmanEncoder dcHuffman_;
    JpegHuffmanEncoder acHuffman_;

    virtual itpp::bvec DcEncode_(int input);
    
    virtual itpp::bvec AcEncode_(int input, int* run);
    
  public:
    // コンストラクタは2種類
    JpegEntropyEncoder(const JpegHuffmanEncoder& dcHuffman, const JpegHuffmanEncoder& acHuffman): dcHuffman_(dcHuffman), acHuffman_(acHuffman)
    { }
    JpegEntropyEncoder(const JHUFF_TBL& dcHuffTable, const JHUFF_TBL& acHuffTable):
      dcHuffman_(dcHuffTable), acHuffman_(acHuffTable)
    { }
    
    // 8*8の1ブロックのみの符号化
    itpp::bvec operator ()(const itpp::Vec< JCOEF >& input)
    {
      assert(static_cast< int >(input.size()) == IMG_DCT_SIZE2);
      // DC成分
      itpp::bvec dcCode = DcEncode_(input[0]);
      
      // AC成分
      itpp::bvec acCode(0);
      int run = 0;
      for (int n = 1; n < IMG_DCT_SIZE2; ++n){
        // 係数が0でなければ
        if (input[n] != 0){
          acCode = itpp::concat(acCode, AcEncode_(input[n], &run));
          run = 0;
        } // if input[n]
        else{
          if (n == IMG_DCT_SIZE2-1){
            itpp::bvec temp = acHuffman_(EOB_INDEX);
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

  };

  
  class JpegEntropyDecoder
  {
  protected:
    JpegHuffmanDecoder dcHuffman_;
    JpegHuffmanDecoder acHuffman_;

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
    virtual bool DcDecode_(const itpp::bvec& input, itpp::bvec* output, int* result);
    
    virtual bool AcDecode_(const itpp::bvec& input, itpp::bvec* output, itpp::ivec* results);
    
  public:
    // コンストラクタは2種類
    JpegEntropyDecoder(const JpegHuffmanDecoder& dcHuffman, const JpegHuffmanDecoder& acHuffman): dcHuffman_(dcHuffman), acHuffman_(acHuffman)
    { }
    JpegEntropyDecoder(const JHUFF_TBL& dcHuffTable, const JHUFF_TBL& acHuffTable):
      dcHuffman_(dcHuffTable), acHuffman_(acHuffTable)
    { }
    
    /************************************************************************************
     * Decode -- エントロピー復号
     * もし途中で復号エラーが起きたら、それ以降は0で埋めるという処理を行う
     * 
     * Arguments:
     *   input -- 入力バイナリデータ
     *
     * Return Value:
     *   itpp::Vec< int > -- 出力データ系列
     ************************************************************************************/
    itpp::ivec operator ()(const itpp::bvec& input)
    {
      itpp::bvec t_output;
      int dcCoeff;
      if (DcDecode_(input, &t_output, &dcCoeff) == false){
        return itpp::ivec(IMG_DCT_SIZE2);
      }
      
      itpp::bvec t_input = t_output;
      itpp::ivec acCoeff;
      if (AcDecode_(t_input, &t_output, &acCoeff) == false){ // 一応エラーが出たらチェックする
        std::cerr << "Error has occured in AcDecode."  << std::endl;
      }

      acCoeff.ins(0, dcCoeff);
      return acCoeff;
    }
  };  

  
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
    itpp::bvec EncodeDc(const itpp::Vec< JCOEF >& input);
    
    itpp::bvec EncodeAc(const itpp::Vec< JCOEF >& input);
    
    itpp::bvec EncodeAll(const itpp::Vec< JCOEF >& input, int blocks)
    {      
      itpp::bvec dcEncoded = EncodeDc(input.left(blocks));
      itpp::bvec acEncoded = EncodeAc(input.get(blocks, -1));

      return itpp::concat(dcEncoded, acEncoded);
    }

    // layerのmaxは64で
    itpp::bvec EncodePart(const itpp::Vec< JCOEF >& input, int blocks, int layers);
  };

  /************************************************************************************
   * EntropyDecoder 
   * 
   * JpegEntropyDecoderのDC成分とAC成分を別々に行えるようにしたもの
   ************************************************************************************/
  class EntropyDecoder: public JpegEntropyDecoder
  {
  protected:
    virtual int AcDecode_(const itpp::bvec& input, itpp::bvec* output, itpp::Vec< int >* results, int size,
                           int effectiveSize = -1);
    
  public:
    EntropyDecoder(const JpegHuffmanDecoder& dcHuffman, const JpegHuffmanDecoder& acHuffman):
      JpegEntropyDecoder(dcHuffman, acHuffman){ }
    EntropyDecoder(const JHUFF_TBL& dcHuffTable, const JHUFF_TBL& acHuffTable):
      JpegEntropyDecoder(dcHuffTable, acHuffTable){ }
    virtual ~EntropyDecoder(){ }

    itpp::Vec< int > DecodeDc(const itpp::bvec& input, int size);
    itpp::Vec< int > DecodeDc(const itpp::bvec& input, int size, itpp::bvec* rest);
    
    itpp::Vec< int > DecodeAc(const itpp::bvec& input, int size);
    itpp::Vec< int > DecodeAc(const itpp::bvec& input, int size, itpp::bvec* rest);

    // 復号できた係数の数を返す
    // つまりoutputの中の意味のある数字の数
    int DecodeAll(const itpp::bvec& input, itpp::ivec* output, int originalSize, int effectiveSize = -1);
    
    // effectiveSizeは復号に使えるバイナリのサイズ
    itpp::Vec< int > DecodeAll(const itpp::bvec& input, int originalSize, int effectiveSize = -1);

    
  };
  
}

#endif
