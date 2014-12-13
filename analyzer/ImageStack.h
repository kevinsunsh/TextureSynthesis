#ifndef _IMAGESTACK_H__
#define _IMAGESTACK_H__

#include <OpenImageIO/imagebuf.h>

OIIO_NAMESPACE_USING

class ImageStack
{
public:
  static unsigned int wrapAccess(int c, unsigned int size)
  {
    int ct = c % size;
    return (ct < 0) ? ct + size : ct;
  }
private:
  void bilinear(const ImageBuf* img, const float u, const float v, float* data) const
  {
    int width = img->spec().width;
    int height = img->spec().height;
    // translate into "Array coordinate space"
    float i  = u * width - 0.5f;
    float j  = v * height - 0.5f;

    // get corners of interpolation cell
    int i0   = int(floor(i));
    int i1   = i0 + 1;
    int j0   = int(floor(j));
    int j1   = j0 + 1;

    // apply access policy
    i0 = wrapAccess(i0,width);
    i1 = wrapAccess(i1,width);
    j0 = wrapAccess(j0,height);
    j1 = wrapAccess(j1,height);

    // compute interpolation fractions
    float fi = (i - int(floor(i)));
    float fj = (j - int(floor(j)));

    // interpolate each component and return
    float ij00[3], ij01[3], ij10[3], ij11[3];
    img->getpixel(i0, j0, ij00);
    img->getpixel(i0, j1, ij01);
    img->getpixel(i1, j0, ij10);
    img->getpixel(i1, j1, ij11);
    for (int i = 0; i < 3; ++i)
    {
      data[i] = (1.0f-fi)*((1.0f-fj)*ij00[i]
                + (fj)*           ij01[i])
                + (fi)*((1.0f-fj)*ij10[i]
                + (fj)*           ij11[i]);
    }
  }


protected:

  std::vector<ImageBuf*> m_Levels;

public:

  ImageStack(uint l) { m_Levels.resize(l); }
  ImageStack(std::vector<ImageBuf*> pyr)
  {
    m_Levels.resize( pyr.size() );
    for (int l = 0; l < m_Levels.size(); ++l)
    {
      int width = pyr[0]->spec().width;
      int height = pyr[0]->spec().height;
      int nc = pyr[0]->spec().nchannels;
      ImageSpec specOutput(width, height, nc, TypeDesc::FLOAT);
      m_Levels[l] = new ImageBuf(specOutput);
      for (int i = 0; i < width; ++i)
        for (int j = 0; j < height; ++j)
      {
        float fi = (i+0.5f) / float(width);
        float fj = (j+0.5f) / float(height);
        float clr[3];
        bilinear(pyr[l], fi, fj, clr);
        m_Levels[l]->setpixel(i , j, clr);
      }
    }
  }

  unsigned int numLevels() const {return m_Levels.size();}

  //! Accessors to single Level 
  const ImageBuf* level(const unsigned int l) const { return m_Levels[l];}
  ImageBuf*       level(const unsigned int l)       { return m_Levels[l];}
};

#endif //_IMAGESTACK_H__