#ifndef _SYNTHESIZER_H__
#define _SYNTHESIZER_H__

#include <boost/multi_array.hpp>
#include "../analyzer/Analyzer.h"

class Synthesizer
{
public:
  typedef boost::multi_array<Imath::V2s,2> SynthesisData;

  //! correctionSubpassInRegion processes pixels of a sub-region of the synthesized image
  static void correctionSubpassForOne(int pixel_index, int level,
                                      const Imath::V2s& subpass_index,
                                      const Synthesizer* theSynthesizer,
                                      const std::vector<Analyzer::KNearest>& nrst,
                                      const SynthesisData& synthesis,
                                      SynthesisData& synthesis_out);
private:

  Analyzer&                             m_Analyzer;       // Analyzer holding exemplar data
  std::vector<SynthesisData>            m_Synthesized;    // The number of entries correspond to the number of upsampling steps applied
  int                                   m_StartLevel;     // Exemplar statck level at which synthesis was started
  float                                 m_Kappa;          // Controls whether coherent candidates are favored; 1.0 has no effect, 0.1 has strong effect, 0.0 is invalid.
  float                                 m_JitterStrength; // Controls jitter strength. 
  int                                   m_Subpasslevel;

  /**
  The three main steps of the algorithm
  */
  //! upsamples previous level synthesis result
  void upsample                 (const SynthesisData& parent, SynthesisData& _child);
  //! adds jitter
  void jitter                   (float strength, SynthesisData& synthesis);
  //! correct neighborhoods to ensure result is visually similar to exemplar
  void correction               (SynthesisData& synthesis);

  /**
  Sub-pass mechanism
  */
  //! correctionSubpass processes pixels in an interleaved pattern aligned with ci,cj
  //! m_NumThreads are created, each calling correctionSubpassInRegion
  void correctionSubpass        (const Imath::V2s& index, SynthesisData& synthesis);
 
/**
  Helper methods
  */
  //! gather a neighborhood in the current synthesis result
  Analyzer::Neighborhood gatherNeighborhood(int step,int i,int j) const; 
  //! returns the exemplar level that must be used at the current synthesis step
  int  currentExemplarLevel();
  //! colorizes current synthesis result (synthesis results are made of exemplar pixel coordinates)
  ImageBuf*              colorize(int step);

public:
  /**
  Constructor - an analyzer must be given since it holds pre-computed data on the exemplar
  */
  Synthesizer(Analyzer& a);
  ~Synthesizer();

  /**
  Initializes synthesis - call first!
  
  w,h is the resolution of the top-most level of the pyramid
  It is equivalent to the number of times the exemplar will appear along each axis 
  */
  void         init(unsigned int w = 512, unsigned int h = 512, float jitterStrength = 25.0f, float kappa = 1.0f, int subpasslevel = 2); 

  /**
  Synthesizes the next level of the multi-resolution pyramid.
  - produces an error if done() is true
  */
  void         synthesizeNextLevel();
  
  /**
  Returns true if last level was synthesized.
  */
  bool         done();
  
  //! returns current result
  ImageBuf* result();
  //! returns color-coded patches for the current result
  ImageBuf* resultPatches();
};

#endif // _SYNTHESIZER_H__