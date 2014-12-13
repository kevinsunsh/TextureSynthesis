#ifndef _ANALYZER_H__
#define _ANALYZER_H__


#include <OpenImageIO/imagebuf.h>
#include <OpenEXR/ImathVec.h>
#include <vector>

#include "ImageStack.h"

OIIO_NAMESPACE_USING
#define N 5
#define K 8
#define DIM 3
#define VN_05 25
#define VN 12
#define M_NUM 3
#define D_NUM 4

class Analyzer
{
public:
  class KNearest
  {
  public:
    KNearest()
    {
      for (int i = 0; i < K; ++i)
      {
        coords[i] = Imath::V2s(-1, -1);
      }
    }
    Imath::V2s coords[K];
  };

  class Neighborhood
  {
    static const Imath::V2s Delta[D_NUM];

    static const Imath::V2s MM[M_NUM][2];
    float pixel[DIM * VN];
  public:
    Neighborhood()
    {
      memset(pixel, 0, sizeof(float)* DIM * VN);
    }
    Neighborhood(const Neighborhood& nb)
    {
      *this=nb;
    }

    //static void ForNeighborhood_05(std::function<void(int, int)> func)
    //{
    //  for (int i = 0; i < VN_05; ++i)
    //  {
    //    func(i);
    //  }
    //}

    static void ForNeighborhood(std::function<void(int, int, int)> func)
    {
      for (int i = 0; i < D_NUM; ++i)
      {
        for (int j = 0; j < M_NUM; ++j)
        {
          int x = MM[j][0].dot(Delta[i]);
          int y = MM[j][1].dot(Delta[i]);
          x = Delta[i][0] + x;
          y = Delta[i][1] + y;
          int index = j + M_NUM * i;
          func(x, y, index);
        }
      }
    }

    const Neighborhood& operator=(const Neighborhood &nb)
    {
      memcpy(pixel, nb.pixel, sizeof(float) * VN * DIM);
      return *this;
    }

    Neighborhood
    operator - (const Neighborhood &v) const
    {
      Neighborhood ret(*this);
      for (int i = 0; i < VN; ++i) {
        for (int d = 0; d < DIM; ++d)
          ret.pixel[i*DIM + d] -= v.pixel[i*DIM + d];
      }

      return ret;
    }

    void getPixel(int index, float* clr) const
    {
      clr[0] = pixel[index * DIM + 0];///M_NUM
      clr[1] = pixel[index * DIM + 1];///M_NUM
      clr[2] = pixel[index * DIM + 2];///M_NUM
    }

    void setPixel(int index, const float* clr)
    {
      int d_i = index;///M_NUM
      pixel[d_i*DIM + 0] = clr[0];//+
      pixel[d_i*DIM + 1] = clr[1];//+
      pixel[d_i*DIM + 2] = clr[2];//+
    }

    float sqLength()
    {
      float sum = 0.0f;
      for (int i = 0; i < VN; ++i) {
        for (int d = 0; d < DIM; ++d)
          sum += pow(pixel[i*DIM + d], 2.0f);
      }
      return sum;
    }
  };
  static void analyzeLevel(int level, Analyzer* theAnalyzer);
private:
  //std::string                               m_Name;          // Exemplar name
  ImageBuf*                                              m_Exemplar;      // Exemplar image
  ImageBuf*                                              m_PCAExemplar;
  ImageStack*                                 m_Stack;         // Exemplar stack, computed from the image
  std::vector<std::vector<KNearest> >     m_KNearests;     // k-most similar neighborhoods within same exemplar stack level
  std::vector<std::vector<Neighborhood> > m_Neighborhoods; // All neighborhoods (pre-gathered for efficiency)
  int                                                    m_NumThreads;    // Number of threads to be used

  //! analyzes the exemplar stack, level per level
  void analyzeStack();
  //! analyzes one exemplar stack level
  void analyzeStackLevel  (int l);
  //! gathers all neighborhoods of the exemplar stack level
  void gatherNeighborhoods(int l, std::vector<Neighborhood>& _neighs);
  //! gathers neighborhood at i,j in the stack level l
  Neighborhood gatherNeighborhood (int l,int i,int j) const;
  
  void GenPyramidsEx();

public:
  
  /**
  Constructor - takes exemplar name and image as input
  */
  Analyzer(ImageBuf* ex, ImageBuf* pca);
  ~Analyzer();

  /**
  Runs analysis. Ideally the result would be saved for later reuse. In this 
  simple implementation nothing gets saved.
  */
  void  run();

  /**
  Returns the neighborhood at i,j in the stack level l. This is using pre-gathered neighborhoods.
  It is meant to be called after analysis, during synthesis.
  */
  const Neighborhood& neighborhoodAt(int l, int i, int j) const;

  /**
  Accessors
  */

  const ImageBuf*                            ex()    { return (m_Exemplar);  }
  const ImageStack*         stack() { return (m_Stack);     }
  const std::vector<std::vector<KNearest> >& kNrst() { return (m_KNearests); }
};

#endif // _ANALYZER_H__