#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include "WienerInterpolator.h"
#include "UniformRNG.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>


using namespace boost::numeric::ublas;
using namespace gplib;

int main()
  {
    std::string filename;
    std::ifstream infile;
    UniformRNG Random;
    /*std::cout << "Infile name: ";
     std::cin >> filename;
     std::vector<double> read;
     infile.open(filename.c_str());
     std::copy(std::istream_iterator<double>(infile),std::istream_iterator<double>(),std::back_inserter(read));*/

    const int datasize = 500;
    const int filterlength = 20;
    const int inputlength = 100;
    gplib::rvec Data(datasize);
    for (int i = 0; i < datasize; ++i)
      Data(i) = sin(static_cast<double> (i) / 0.9) + cos(
          static_cast<double> (i) / 2.7) + (0.5 - Random.GetNumber()) * 0.1;
    //copy(read.begin(),read.end(),Data.begin());
    WienerInterpolator Wiener(filterlength);
    gplib::rvec Input(inputlength);
    gplib::rvec Prediction(datasize);
    for (int i = 0; i < inputlength; ++i)
      {
        Input(i) = Data(i);
        Prediction(i) = Data(i);
      }
    Wiener.AdaptFilter(Input, Input);

    for (int i = inputlength - filterlength; i < datasize - filterlength; ++i)
      {
        vector_range<vector<double> >
            vr(Prediction, range(i, i + filterlength));
        Prediction(i + filterlength) = inner_prod(Wiener.GetWeights(), vr);
      }
    std::ofstream outfile("out");
    for (int i = 0; i < datasize; ++i)
      outfile << i << " " << Data(i) << " " << Prediction(i) << std::endl;
    //copy(Wiener.GetWeights().begin(),Wiener.GetWeights().end(),std::ostream_iterator<double>(std::cout,"\n"));
  }
