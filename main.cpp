#include <iostream>
#include <ctime>
#include <cmath>
#include <OpenImageIO/imagebuf.h>
//#include <pca.h>
// --------------------------------------------------------------
// --------------------------------------------------------------
#include "synthesizer/Synthesizer.h"
#include "analyzer/Analyzer.h"

// --------------------------------------------------------------

OIIO_NAMESPACE_USING
using namespace std;

// --------------------------------------------------------------
//ImageBuf* DoPCA1(ImageBuf* pOutBuf) {
//    int image_size_u = pOutBuf->spec().width;
//    int image_size_v = pOutBuf->spec().height;
//    ImageSpec specOutput(image_size_u, image_size_v, 3, TypeDesc::FLOAT);
//    ImageBuf* pOutBuf_2 = new ImageBuf(specOutput);
//    int n = pOutBuf->spec().width*pOutBuf->spec().height;
//    int dim = 5;
//    int m = 75;
//    ImageBuf::Iterator<float> it_pOutBuf(*pOutBuf);
//    ImageBuf::Iterator<float> it_pOutBuf_2(*pOutBuf_2);
//    float* pBegin_1 = (float*)it_pOutBuf.rawptr();
//    float* pBegin_2 = (float*)it_pOutBuf_2.rawptr();
//    int size = pOutBuf_2->spec().width;
//    stats::pca pca(m);
//    pca.set_do_bootstrap(true, 100);
//    vector<double> record[n];
//    int cnt = 0;
//    for(int b = 0; b < size; b++) {
//        for(int a = 0; a < size; a++) {
//            ++it_pOutBuf;
//            ++it_pOutBuf_2;
//            int u = 0; int v = 0;
//            for(int j = b - dim/2; j <= b + dim/2;j++) {
//                for(int i = a - dim/2; i <= a + dim/2; i++) {
//                    int x = i; int y = j;
//                    if(x < 0) x = abs(i); if(y < 0) y = abs(y);
//                    if(x > size - 1) x = 2*size - x -2 ; if(y > size - 1) y = 2*size - y -2;
//                    record[cnt].push_back(pBegin_1[(x + y*size) * 3 + 0]);
//                    record[cnt].push_back(pBegin_1[(x + y*size) * 3 + 1]);
//                    record[cnt].push_back(pBegin_1[(x + y*size) * 3 + 2]);
//                    u++;
//                }
//                u = 0; v++;
//            }
//            pca.add_record(record[cnt]);
//            cnt++;
//        }
//    }
//
//    pca.solve();
//    std::vector<double> principal[m];
//    for( int i = 0;  i < m; i++) {
//        principal[i] = pca.get_principal(i);
//    }
//    float min = FLT_MAX; float max = -FLT_MAX;
//    for(unsigned int i = 0; i < principal[0].size(); i++) {
//        if(min > principal[0][i]) min = principal[0][i];
//        if(max < principal[0][i]) max = principal[0][i];
//        if(min > principal[1][i]) min = principal[1][i];
//        if(max < principal[1][i]) max = principal[1][i];
//        if(min > principal[2][i]) min = principal[2][i];
//        if(max < principal[2][i]) max = principal[2][i];
//    }
//    float fsize = max - min;
//    for(unsigned int i = 0; i < principal[0].size(); i++) {
//        principal[0][i] = (principal[0][i] - min)/fsize;
//        principal[1][i] = (principal[1][i] - min)/fsize;
//        principal[2][i] = (principal[2][i] - min)/fsize;
//    }
//    for( int j = 0;  j < n; j++) {
//        pBegin_2[j*3 + 0] = (principal[0][j]);
//        pBegin_2[j*3 + 1] = (principal[1][j]);
//        pBegin_2[j*3 + 2] = (principal[2][j]);
//    }
//    return pOutBuf_2;
//}

int main(int argc, char **argv)
{
  int s = 0;

  // load the exemplar

  ImageBuf* ex = new ImageBuf("TestData/stone3_exemplar.png");
  //ImageBuf* pcaBuf = DoPCA1(ex);
  // init the analyzer
  Analyzer    analyzer(ex, ex);
  // init the synthesizer
  Synthesizer synthesizer(analyzer);

  analyzer.run();

  synthesizer.init( 512, 512, 25.0f, 0.2f , 2);
  //    // go down synthesis pyramid until finest level reached
  while (!synthesizer.done()) {
  // synthesize next level
    synthesizer.synthesizeNextLevel();
  }

  synthesizer.result()->save(std::string("testsynth.png"), std::string("png"));
  synthesizer.resultPatches()->save(std::string("testcoords.png"), std::string("png"));

  return (0);
}

/* -------------------------------------------------------- */
