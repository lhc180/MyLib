#ifndef HUFFMAN_H
#define HUFFMAN_H
/************************************************************************************
 * JPEG Entropy Coder
 *   
 * JPEGで使われているエントロピー符号化のエンコード・デコードを行うクラス
 * つまり、ハフマン符号とランレングスの組み合わせ
 *
 * Contents:
 *   Encode()
 *   Decode()
 *   SetTable()
 * Last Updated: <2013/05/23 17:55:52 from Yoshitos-iMac.local by yoshito>
 ************************************************************************************/
#include <cassert>
#include <itpp/itcomm.h>
#include <myutl.h>
#include <jpeglib.h>

namespace mylib{
  
  struct HuffmanTable{
    char sizeMax;
    //     int numElement;
    std::vector< char >  sizeTable; // それぞれの符号が何ビットで構成されているか
    std::vector< int >  codeTable; // 実際の符号
  };


  class Huffman
  {
  private:
    HuffmanTable table_;
        
    // 以下デコード用のメンバ
    struct Convert
    {
      int  code_;
      char size_;
      int  index_;
    };

    //-----------------------------
      // set 用 operator<

    //   bool                    // out: 大きいかどうか
    //                           // in : 比較対象
    //   operator<(const Convert& in) const
    //   {
    //     return code_ < in.code_;
    //   }
    // };
    // std::set< Convert > decodeSet_;
    
  protected:
    // 見つかったらtrueを返す
    bool Search(Convert* candidate)
    {
      while (table_.sizeTable[candidate->index_] == candidate->size_){
        if (table_.codeTable[candidate->index_] == candidate->code_){
          return true;
        } // if code_
        candidate->index_++;
      } // while
      return false; 
    }

   public:
     explicit Huffman(const HuffmanTable& table): table_(table)
     {
       assert(table_.sizeTable.size() == table_.codeTable.size());
       assert(table_.sizeMax != 0);

       // デコーダで使うdecodeSet_を設定する
       // for (int i = 0, max = table_.sizeTable.size(); i < max; ++i){
       //   if (table_.sizeTable[i] == 0){
       //     continue;
       //   } // if table_.sizeTable[i]

       //   Convert cv = {
       //     table_.codeTable[i],
       //     table_.sizeTable[i],
       //     i
       //   };
       //   decodeSet_.insert(cv);
       // } // for i

       // assert(!decodeSet_.empty());

     }
    // virtual ~Huffman() -- デフォルトデストラクタ

    // コピーコンストラクタはデフォルト
    
     // エンコーダ
    itpp::bvec Encode(int input) const
    {
      assert(input >= 0 && input < static_cast< int >(table_.sizeTable.size()));
      itpp::bvec output = itpp::dec2bin(table_.sizeTable[input], table_.codeTable[input]);
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
      Convert candidate = {
        0,                      // code_
        0,                      // size_
        0,                      // index_
      };
      
      while (candidate.index_ < static_cast< int >(table_.sizeTable.size()) &&
             candidate.size_ < table_.sizeMax && 
             candidate.size_ < input.size()){
        candidate.size_++;
        itpp::bvec temp = input.mid(0, candidate.size_);
        candidate.code_ = itpp::bin2dec(temp);

        if(Search(&candidate)){
          *output = input.get(candidate.size_, -1);
          return candidate.index_;
        } // if
      } // while 

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
      
}

#endif
