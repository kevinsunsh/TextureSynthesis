#pragma warning disable 873
#pragma warning disable 2621
#include <flann/flann.hpp>
#include <tbb/tbb.h>
#include <OpenImageIO/imagebufalgo_util.h>
#include <OpenImageIO/imagebufalgo.h>
#include <math.h>

#include "Analyzer.h"

using namespace std;
using namespace tbb;
// --------------------------------------------------------------
const Imath::V2s Analyzer::Neighborhood::Delta[] = {Imath::V2s(-1,-1), Imath::V2s(1,-1), Imath::V2s(-1,1), Imath::V2s(1,1)};
const Imath::V2s Analyzer::Neighborhood::MM[][2] = {{Imath::V2s(0,0),
                                  Imath::V2s(0,0)},
                                 {Imath::V2s(1,0),
                                  Imath::V2s(0,0)},
                                 {Imath::V2s(0,0),
                                  Imath::V2s(0,1)}};
//! simple function to check that a number is a power of two
static bool isPow2(int v)
{
  int n = 0;
  while (v > 0) {
    if (v & 1) { n ++; }
    v >>= 1;
  }
  return (n == 1);
}

// --------------------------------------------------------------

Analyzer::Analyzer(ImageBuf* ex, ImageBuf* pca)
{
  assert(isPow2(ex->spec().width) || ex->spec().width == ex->spec().height);
  m_Exemplar = ex;
  m_PCAExemplar = pca;
}

// --------------------------------------------------------------
ImageBuf* GeneratePyramid(ImageBuf* pImage)
{
  int size_u = pImage->spec().width/2;
  int size_v = pImage->spec().height/2;
  int nc = pImage->spec().nchannels;
  ImageSpec specOutput(size_u, size_v, nc, TypeDesc::FLOAT);
  ImageBuf* pOutputBuf = new ImageBuf(specOutput);
  ImageBufAlgo::resize(*pOutputBuf, *pImage, "", 0, ROI(0, specOutput.width, 0, specOutput.height));
  return pOutputBuf;
}

// --------------------------------------------------------------

void Analyzer::GenPyramidsEx()
{
  std::vector<ImageBuf*> pyramids;
  int level_count = log2(m_PCAExemplar->spec().width);
  pyramids.push_back(m_PCAExemplar);
  for (int i = 0; i < level_count; ++i)
  {
    pyramids.push_back(GeneratePyramid(pyramids.back()));
  }
  
  m_Stack = new ImageStack(pyramids);
}

// --------------------------------------------------------------

void Analyzer::run()
{
  GenPyramidsEx();
  // analyze stack
  analyzeStack();
}

// --------------------------------------------------------------
void Analyzer::analyzeLevel(int level, Analyzer* theAnalyzer)
{
  ImageBuf* img = theAnalyzer->m_Stack->level(level);
  theAnalyzer->m_KNearests[level].resize( img->spec().width * img->spec().height );
  theAnalyzer->gatherNeighborhoods(level, theAnalyzer->m_Neighborhoods[level]);
  theAnalyzer->analyzeStackLevel(level);
}

// --------------------------------------------------------------
void Analyzer::analyzeStack()
{
  // Analyzes
  int level_count = m_Stack->numLevels();
  m_KNearests    .resize( level_count );
  m_Neighborhoods.resize( level_count );

//for (int i = 0; i < level_count; ++i)
//{
//  Analyzer::analyzeLevel(i, this); 
//}
  parallel_for( blocked_range<size_t>(0,level_count), 
    [&](const blocked_range<size_t>& r) {
      for(size_t i=r.begin(); i!=r.end(); ++i) 
          Analyzer::analyzeLevel(i, this); 
    }
  );
}

void Analyzer::analyzeStackLevel(int l)
{
  int pixel_count = m_Neighborhoods[l].size();
  float* dataset_buf = new float[pixel_count * DIM * VN];
  int* indexs_buf = new int[pixel_count * K];
  float* dists_buf = new float[pixel_count * K];
  int data_count = 0;

  for (int p_index = 0; p_index < pixel_count; ++p_index)
  {
    for(int i = 0; i < VN; ++i)
    {
      float clr[DIM];
      m_Neighborhoods[l][p_index].getPixel(i, clr);
      dataset_buf[data_count++] = clr[0];
      dataset_buf[data_count++] = clr[1];
      dataset_buf[data_count++] = clr[2];
    }
  }
  //int stride = DIM * VN_05;
  int stride = DIM * VN;
  flann::Matrix<float> dataset(dataset_buf, pixel_count, stride);
  flann::Matrix<float> query(dataset_buf, pixel_count, stride);
  flann::Matrix<int> indices(indexs_buf, pixel_count, K);
  flann::Matrix<float> dists(dists_buf, pixel_count, K);

  flann::IndexParams* indexParams = new flann::KDTreeSingleIndexParams();
  if (indexParams == NULL) return;
  
  const float eps = 1.0f; // colors are in [0..255]
  flann::SearchParams* sParams = new flann::SearchParams(128);
  if (sParams == NULL) return;
  sParams->eps = eps;
  sParams->max_neighbors = K;

  flann::Index<flann::L2<float> > index(dataset, *indexParams);
  index.buildIndex();

  // do a knn search, using 128 checks
  int hr = index.knnSearch(query, indices, dists, K, *sParams);

  for (int i = 0; i < pixel_count; ++i)
  {
    for (int j = 0; j < K; ++j)
    {
      int sqrt_count = sqrt(pixel_count);
      if (indexs_buf[i*K + j] == -1)
        break;
      m_KNearests[l][i].coords[j][0] = indexs_buf[i*K + j] % sqrt_count;
      m_KNearests[l][i].coords[j][1] = indexs_buf[i*K + j] / sqrt_count;
    }
  }
  delete indexParams;
  delete sParams;

  delete [] dataset_buf;
  delete [] indexs_buf;
  delete [] dists_buf;
}

// --------------------------------------------------------------

void Analyzer::gatherNeighborhoods(int l, std::vector<Neighborhood>& _neighs)
{
  ImageBuf* img = m_Stack->level(l);
  int width = img->spec().width;
  int height = img->spec().height;
  _neighs.resize(width * height);
  // gather neighborhoods
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      // extract neighborhood
      _neighs[i + j * width] = gatherNeighborhood(l,i,j);
    }
  }
}

// --------------------------------------------------------------

Analyzer::Neighborhood Analyzer::gatherNeighborhood(int l,int i,int j) const
{
  // Gather a neighborhood within the stack. Note that contrary to neighborhoods
  // in a regular image, neighbors are not next to each others in the stack but
  // separated by a level-dependent offset.
  Neighborhood n;
  ImageBuf* img = m_Stack->level(l);
  int width = img->spec().width;
  int height = img->spec().height;
  int          spacing = (1 << l); // level-dependent offset
  Neighborhood::ForNeighborhood([&](int di, int dj, int index)->void {
      int x  = (i + di * spacing);
      int y  = (j + dj * spacing);
      float clr[DIM];
      img->getpixel(x, y, clr);
      n.setPixel(index, clr);
    }
  );
  return (n);
}

// --------------------------------------------------------------

const Analyzer::Neighborhood& Analyzer::neighborhoodAt(int l,int i,int j) const
{
  // Returns the neighborhood at i,j in level l, using pre-gathered neighborhoods (see analyzeStackLevel)
  assert(l >= 0 && l < int(m_Neighborhoods.size()));
  ImageBuf* img = m_Stack->level(l);
  int width = img->spec().width;
  int height = img->spec().height;
  i = ImageStack::wrapAccess(i, width);
  j = ImageStack::wrapAccess(j, height);
  return m_Neighborhoods[l][i + j * width];
}

// --------------------------------------------------------------

Analyzer::~Analyzer()
{
  delete m_Exemplar;
  delete m_PCAExemplar;
}

