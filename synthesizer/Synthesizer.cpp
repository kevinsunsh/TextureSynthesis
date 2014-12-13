#include <float.h>
#include <assert.h>
#include <tbb/tbb.h>

#include "Synthesizer.h"

using namespace std;
using namespace tbb;

// --------------------------------------------------------------

Synthesizer::Synthesizer(Analyzer& a) : m_Analyzer(a)
{

}

// --------------------------------------------------------------

Synthesizer::~Synthesizer()
{

}

// --------------------------------------------------------------

void Synthesizer::init(uint w,uint h,float jitterStrength,float kappa, int subpasslevel)
{
  assert(kappa > 0.0f);
  // initialize coarsest level to obtain desired resolution at finest level
  int nx = int(ceil(w / float(m_Analyzer.ex()->spec().width)));
  int ny = int(ceil(h / float(m_Analyzer.ex()->spec().height)));
  SynthesisData s_data(boost::extents[nx][ny]);
  boost::multi_array_ref<Imath::V2s, 1> s_data_ref(s_data.data(), boost::extents[s_data.num_elements()]);
  std::fill(s_data_ref.begin(), s_data_ref.end(), Imath::V2s(m_Analyzer.ex()->spec().width/2,m_Analyzer.ex()->spec().height/2));
  
  m_Synthesized.push_back(s_data);

  // init parameters
  m_Kappa          = kappa;
  m_JitterStrength = jitterStrength;

  // start level is coarsest
  m_StartLevel = m_Analyzer.stack()->numLevels() - 1;
  m_Subpasslevel = subpasslevel;

  assert(m_StartLevel > 0); // with current algorithm it makes no sense to start at level 0
}

// --------------------------------------------------------------

int Synthesizer::currentExemplarLevel()
{
  // Compute exemplar level from number of synthesis steps (number of calls to synthesizeNextLevel)
  int l = (m_StartLevel - (int(m_Synthesized.size())-1));
  assert(l >= 0);
  return l;
}

// --------------------------------------------------------------

bool Synthesizer::done()
{
  // Done if finest level has been reached.
  // This happens when the number of perfromed synthesis steps 
  // equals the exemplar level at which synthesis started.
  return (m_Synthesized.size() == m_StartLevel+1);
}

// --------------------------------------------------------------

void Synthesizer::synthesizeNextLevel()
{
  // Performs a synthesis step, producing the result at the next level.
  // The new result is added to m_Synthesized
  assert(!done());
  assert(m_Synthesized.size() > 0);
  // add the next level result

  SynthesisData data(boost::extents[m_Synthesized.back().size() * 2][m_Synthesized.back()[0].size() * 2]);
  m_Synthesized.push_back( data );
  /// 1. upsample
  upsample    ( m_Synthesized[m_Synthesized.size()-2] , m_Synthesized.back() );
  /// 2. jitter
  // adapt jitter strength per level - arbitrary, ideally should be per-level 
  // user control. Overall it is often more desirable to add strong jitter at 
  // coarser levels and let synthesis recover at finer resolution levels.
  float strength = m_JitterStrength * (currentExemplarLevel() < 3 ? 0 : currentExemplarLevel()) / (float)m_Analyzer.stack()->numLevels()+1;
  // apply jitter
  jitter      ( strength , m_Synthesized.back() );
  ///// 3. correct
  for (int p = 0; p < 2; ++p) {
    correction( m_Synthesized.back() );
  }
}

// --------------------------------------------------------------

void Synthesizer::upsample(const SynthesisData& parent, SynthesisData& _child)
{

  int l       = currentExemplarLevel();

  int spacing = (1 << l);
  // next level has twice the resolution of the previous one
  int row = parent.size();
  int column = parent[0].size();
  // coordinate inheritence
  for (int j = 0; j < row; ++j) {
    for (int i = 0; i < column; ++i) {
      //printf("i j = (%d,%d)\n", i, j);
      for (int cj = 0; cj < 2; ++cj) {
        for (int ci = 0; ci < 2; ++ci) {
          _child[j*2 + cj][i*2 + ci] = parent[j][i] + Imath::V2s(ci,cj) * spacing;
        }
      }
    }
  }
}

// --------------------------------------------------------------

void Synthesizer::jitter(float strength, SynthesisData& synthesis)
{
  // coordinate inheritence
  boost::multi_array_ref<Imath::V2s, 1> synthesis_ref(synthesis.data(), boost::extents[synthesis.num_elements()]);
  int pixel_count = synthesis_ref.num_elements();
  for (int i = 0; i < pixel_count; ++i) {
    float temp = (rand()%100000) / 100000.0f;
    temp = (temp - 0.5f) * 2.0f;
    synthesis_ref[i] = synthesis_ref[i] + Imath::V2s(short(strength*temp),short(strength*temp));
  }
}

// --------------------------------------------------------------

void Synthesizer::correction(SynthesisData& synthesis)
{
  // Performs one correction pass, made of four sub-passes
  for (int i_row = 0; i_row < m_Subpasslevel; ++i_row) {
    for (int i_column = 0; i_column < m_Subpasslevel; ++i_column) {
      // apply correction sub-pass
      correctionSubpass(Imath::V2s(i_column, i_row), synthesis);
    }
  }
}

// --------------------------------------------------------------

void Synthesizer::correctionSubpass(const Imath::V2s& subpass_index, SynthesisData& synthesis)
{
  // temporary copy of the buffer, all reads occur in _S, all writes in tmp
  SynthesisData tmp = synthesis;
  // _S is copied since only a quarter of the pixels are updated by sub-pass correction.
  // NOTE: Please see original publication for details on how to efficiently implement this 
  //       through pixel re-ordering.
  int level = currentExemplarLevel();
  const std::vector<Analyzer::KNearest>& nrst = m_Analyzer.kNrst()[level];
  int pixel_count = synthesis.num_elements();
  parallel_for( blocked_range<size_t>(0,pixel_count), 
   [&](const blocked_range<size_t>& r) {
    for(size_t i=r.begin(); i!=r.end(); ++i) 
        Synthesizer::correctionSubpassForOne(i, level, subpass_index, this, nrst, synthesis, tmp); 
  }
  );
  // done, store result
  synthesis = tmp;
}

// --------------------------------------------------------------

void Synthesizer::correctionSubpassForOne(int pixel_index, int level,
                                            const Imath::V2s& subpass_index,
                                            const Synthesizer* theSynthesizer,
                                            const std::vector<Analyzer::KNearest>& nrst,
                                            const SynthesisData& synthesis,
                                            SynthesisData& synthesis_out)
{
  int column = synthesis[0].size();
  int row = synthesis.size();
  int i_column = pixel_index % column;
  int i_row = pixel_index / column;

  //sub-pass mechanism
  if ((i_column % theSynthesizer->m_Subpasslevel) != subpass_index[0] ||
      (i_row % theSynthesizer->m_Subpasslevel) != subpass_index[1] )
        return;

  int spacing = (1 << level);
  const int numCand = 9*K+1;
  Imath::V2s kcand[numCand];
  // non-coherent candidates stored in [0 ; (9*K-1)], coherent candidates in [9*K ; 9*(K+1)], self in last
  // it is important to put coherent candidates last so that they are chosen in case of tie
  // ties happen constantly in coherent patches
  /// Gather candidates
  // for each neighbor around the pixel (9 of them, including center)
  for (int nj = -1; nj < 2; nj = nj+1) {
    for (int ni = -1; ni < 2; ni = ni+1) {
      int nid = (ni+1)+(nj+1)*3;
      for (int k = 0; k < K; ++k) {
        int x = ImageStack::wrapAccess(i_column + ni, column);
        int y = ImageStack::wrapAccess(i_row + nj, row);
        Imath::V2s n = synthesis[y][x];
        // n is a coordinate in exemplar stack
        // delta must be multiplied by stack level offset
        int kn_length_2 = nrst.size();
        int kn_length = sqrt(kn_length_2);
        n[0] = ImageStack::wrapAccess(n[0], kn_length);
        n[1] = ImageStack::wrapAccess(n[1], kn_length);
        //if(kNrst[n[0] + n[1] * kn_length].coords[k][0] == -1)
        //  break;
        Imath::V2s c = nrst[n[0] + n[1] * kn_length].coords[k] - Imath::V2s(ni,nj) * spacing;
        if (k > 0) {
          //int index = nid*(K-1) + (k - 1);
          //printf("k = %d index = %d\n", k , index);
          kcand[nid*(K-1) + (k - 1)] = c; // non-coherent candidate
        } else {
          //int index = 9*(K-1) + nid;
          //printf("k = %d index = %d\n", k , index);
          kcand[  9*(K-1) + nid] = c; // coherent candidate - we want them to be treated separately in case of tie
        }
      }
    }
  }
  kcand[numCand-1] = synthesis[i_row][i_column]; // self as last -- VERY IMPORTANT to ensure identity in coherent patches --

  /// Gather current neighborhood in synthesized texture
  Analyzer::Neighborhood syN = theSynthesizer->gatherNeighborhood(int(theSynthesizer->m_Synthesized.size())-1, i_column, i_row);

  /// Find best matching candidate
  float mind = FLT_MAX;
  Imath::V2s best = synthesis[i_row][i_column];

  for (int k = 0; k < numCand; ++k) {
    const Analyzer::Neighborhood& exN = theSynthesizer->m_Analyzer.neighborhoodAt(level, kcand[k][0], kcand[k][1]);
    // compare
    float d = ( exN - syN ).sqLength();
    if (k >= 9*(K-1)) {
      d = d * theSynthesizer->m_Kappa; // favor (or defavor) coherent candidates
    }
    if (d <= mind) {
      mind = d;
      best = kcand[k];
    }
  }
  // replace in output
  synthesis_out[i_row][i_column] = best;
}

// --------------------------------------------------------------

Analyzer::Neighborhood Synthesizer::gatherNeighborhood(int step,int i,int j) const
{
  // Gather a neighborhood in the current synthesis result
  Analyzer::Neighborhood n;
  assert(step>=0 && step<int(m_Synthesized.size()));
  assert(m_Synthesized[step].size()>0);
  int l = m_StartLevel - step;
  assert(l>=0 && l<int(m_Analyzer.stack()->numLevels()));
  const ImageBuf* stackLevel = m_Analyzer.stack()->level(l);
  int column = m_Synthesized[step][0].size();
  int row = m_Synthesized[step].size();
  int width = stackLevel->spec().width;
  Analyzer::Neighborhood::ForNeighborhood([&](int di, int dj, int index)->void {
      int x  = (i + di);
      int y  = (j + dj);
      x = ImageStack::wrapAccess(x, column);
      y = ImageStack::wrapAccess(y, row);
      Imath::V2s s = m_Synthesized[step][y][x];  //   S[p]  (coordinate in exemplar stack)
      float clr[DIM];
      stackLevel->getpixel(s[0], s[1], clr);
      n.setPixel(index, clr);
    }
  );
  return n;
}

// --------------------------------------------------------------

ImageBuf* Synthesizer::colorize(int step)
{
  // Create color version of the synthesis result (which contains coordinates only)
  assert(step < int(m_Synthesized.size()));
  const ImageBuf* src = ((m_StartLevel-step) == 0) ? m_Analyzer.ex() : m_Analyzer.stack()->level(m_StartLevel-step);
  int width = src->spec().width;
  int row = m_Synthesized[step].size();
  int column = m_Synthesized[step][0].size();
  ImageSpec specOutput(column, row, 3, TypeDesc::FLOAT);
  ImageBuf* img = new ImageBuf(specOutput);
  for (int j = 0; j < row; ++j) {
    for (int i = 0; i < column; ++i) {
      Imath::V2s xy = m_Synthesized[step][j][i];
      xy[0] = ImageStack::wrapAccess(xy[0], width);
      xy[1] = ImageStack::wrapAccess(xy[1], width);
      float clr[3];
      src->getpixel( xy[0], xy[1], clr );
      img->setpixel(i, j, clr);
    }
  }
  return img;
}

// --------------------------------------------------------------

ImageBuf* Synthesizer::result()
{
  // Current synthesis result
  return colorize(int(m_Synthesized.size())-1);
}

// --------------------------------------------------------------

ImageBuf* Synthesizer::resultPatches()
{
  // Color-code patches produced through synthesis. For visulization purposes only.
  int spacing = (1 << currentExemplarLevel());
  const ImageBuf* src = m_Analyzer.stack()->level(0);
  int width = src->spec().width;
  int row = m_Synthesized.back().size();
  int column = m_Synthesized.back()[0].size();
  ImageSpec specOutput(column, row, 3, TypeDesc::FLOAT);
  ImageBuf* img = new ImageBuf(specOutput);
  for (int j = 0; j < row; ++j) {
    for (int i = 0; i < column; ++i) {
      Imath::V2s xy = m_Synthesized.back()[j][i];
      img->setpixel(i, j, Imath::V3f((xy[0]%width)/float(width), (xy[1] % width)/float(width), 0.0f).getValue());
    }
  }
  return (img);
}

