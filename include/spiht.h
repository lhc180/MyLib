#ifndef SPIHT_H
#define SPIHT_H

#include <vector>
#include <itpp/itbase.h>

namespace mylib {

  struct PixItem 
  {
    int x;
    int y;
    PixItem(int t_x, int t_y): x(t_x), y(t_y)
    { }
  };

  enum LIS_Type {LIS_A, LIS_B};

  struct SetItem
  {
    int x;
    int y;
    LIS_Type type;
    SetItem(int t_x, int t_y, LIS_Type t_type)
      :x(t_x), y(t_y), type(t_type)
    { }
  };

  typedef std::vector<PixItem> LIP;
  typedef std::vector<PixItem> LSP;
  typedef std::vector<SetItem> LIS;

  // +++++++++++++++++ Encoder +++++++++++++++++++
  class SPIHTencoder {
  private:
    int numStages_;
    LIP lip_;
    LSP lsp_;
    LIS lis_;
    int step_;    /* quantization step */
    void GetSuccessor(const itpp::imat& image, int x, int y, int* sx, int* sy);
    bool IsSignificantPixel(const itpp::imat& image, int x, int y);
    bool IsSignificant_SetA(const itpp::imat& image, int x, int y, int count = 1);
    bool IsSignificant_SetB(const itpp::imat& image, int x, int y, int count = 1);
    void Initialize(const itpp::imat& image, itpp::bvec* bout);
  public:
    SPIHTencoder():numStages_(0)
    { }
    SPIHTencoder(int numStages)
    {
      Set(numStages);
    }
    void Set(int numStages)
    {
      numStages_ = numStages;
    }
    // bitsは出力したいビット数。-1は全ビットを表す
    void Encode(const itpp::imat& image, int bits, itpp::bvec* bout);
    itpp::bvec Encode(const itpp::imat& image, int bits){
      itpp::bvec output;
      Encode(image, bits, &output);
      return output;
    }
  };

  // +++++++++++++++++++++ Decoder +++++++++++++++++++
  class SPIHTdecoder {
    int rows_;                  // 画像のサイズ
    int cols_;
    int numStages_;
    LIP lip_;
    LSP lsp_;
    LIS lis_;
    int step_;    /* quantization step */
    void GetSuccessor(int x, int y, int* sx, int* sy);
    void Initialize(const itpp::bvec& bin, itpp::imat* image);
  public:
    SPIHTdecoder(): rows_(0), cols_(0)
    { }
    SPIHTdecoder(int numStages, int rows, int cols)
    {
      Set(numStages, rows, cols);
    }
    void Set(int numStages, int rows, int cols)
    {
      numStages_ = numStages;
      rows_ = rows;
      cols_ = cols;
    }
    void Decode(const itpp::bvec& bin, itpp::imat* imageOut);
    itpp::imat Decode(itpp::bvec& bin)
    {
      itpp::imat output;
      Decode(bin, &output);
      return output;
    }
  };

  
}

#endif
