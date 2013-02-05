#ifndef C_CPP_MEDIATION_H
#define C_CPP_MEDIATION_H

#include <vector>
#include <itpp/base/vec.h>

template<typename kind>
void arrayToVector(const kind *array, std::vector<kind> &vec, const unsigned data)
{
  vec.resize(data);
  for(unsigned i = 0; i < data; i++){
    vec[i] = array[i];
  }
}

template<typename kind>
void arrayToVector(const kind *array, itpp::Vec<kind> &vec, const unsigned data)
{
  vec.set_size(data);
  for(unsigned i = 0; i < data; i++){
    vec[i] = array[i];
  }
}


template<typename kind>
void vectorToArray(const std::vector<kind> &vec, kind *array)
{
  for(unsigned i = 0; i < vec.size(); i++){
    array[i] = vec[i];
  }
}

template<typename kind>
void vectorToArray(const itpp::Vec<kind> &vec, kind *array)
{
  for(unsigned i = 0; i < vec.size(); i++){
    array[i] = vec[i];
  }
}

#endif
