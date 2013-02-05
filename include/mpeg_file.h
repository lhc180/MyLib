#ifndef MPEG_FILE_H
#define MPEG_FILE_H

class CMpegFile
{
private:
  FILE* m_pFile;

  BOOL DecodeSideInfo( CBitStream* pbs, CMpegFrame* pmf );

public:
  CMpegFile();
  ~CMpegFile();

  BOOL Open( const char* pFilename );
  void Close();
  BOOL IsOpen() const { return m_pFile!=NULL; }
  BOOL SeekMpegFrame( CMpegFrame* pMpegFrame );
};

#endif
