#ifndef MYMODULATION_H
#define MYMODULATION_H

#include <cassert>
#include <math.h>
#include <itpp/itcomm.h>

namespace mylib {


  /************************************************************************************
   * SetPartitioningPSK
   * 
   * セットパーティショニング(ナチュラル)マッピングのPSK
   ************************************************************************************/
  class SetPartitioningPSK: public itpp::PSK
  {
  protected:
    virtual void Init();
    
  public:
    explicit SetPartitioningPSK(int m): itpp::PSK(m)
    {
      Init();
    }
    virtual ~SetPartitioningPSK()
    { }
  };

  /************************************************************************************
   * ReverseSpPSK 
   * 
   * セットパーティショニングの逆マッピング
   ************************************************************************************/
  class ReverseSpPSK: public itpp::PSK
  {
  protected:
    virtual void Init();
    
  public:
    explicit ReverseSpPSK(int m): itpp::PSK(m)
    {
      Init();
    }
    virtual ~ReverseSpPSK()
    { }
  };
  
  
  // itpp::Modulatorのデフォルトのグレイマッピングがそのままblock partitioning
  class BlockPartitioningPSK: public itpp::PSK
  {
  protected:
    virtual void Init(double angle);
    
  public:
    BlockPartitioningPSK(int m, double angle = M_PI/4.0): itpp::PSK(m)
    {
      Init(angle);
    }
    virtual ~BlockPartitioningPSK()
    { }
  };


  class BlockPartitioningQAM: public itpp::QAM
  {
  protected:
    virtual void Init()
    {
      assert(M == 16);          // まだ16QAMにしか対応していない
      itpp::cvec newSymbols = symbols;
      itpp::ivec newBits2Symbol = "0 1 4 5 2 3 6 7 8 9 12 13 10 11 14 15";

      set(newSymbols, newBits2Symbol);
      
    }
    
  public:
    explicit BlockPartitioningQAM(int m): itpp::QAM(m)
    {
      Init();
    }
    virtual ~BlockPartitioningQAM()
    { }
  };

  // Hybrid Partitioning
  class HybridPartitioningPSK_1: public itpp::PSK
  {
  protected:
    virtual void Init(double angle);
    
  public:
    explicit HybridPartitioningPSK_1(int m, double angle = M_PI/4.0): itpp::PSK(m)
    {
      Init(angle);
    }
    virtual ~HybridPartitioningPSK_1() { }
  };

  class HybridPartitioningPSK_2:public itpp::PSK
  {
  protected:
    virtual void Init(double angle);
    
  public:
    explicit HybridPartitioningPSK_2(int m, double angle = M_PI/4.0): itpp::PSK(m)
    {
      Init(angle);
    }
    virtual ~HybridPartitioningPSK_2() { }
  };

  // class BlockPartitioning16APSK:public itpp::Modulator_2D
  // {
  // protected:
  //   virtual void Init();
    
  // public:
  //   explicit BlockPartitioning16APSK();
  //   virtual ~BlockPartitioning16APSK() { }
  // };

}

#endif
