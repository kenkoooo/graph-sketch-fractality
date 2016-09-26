#pragma once

#include <string>
#include <typeinfo>

namespace agl {
std::string demangle_typename(const char* name);

template <class T>
std::string typename_of(const T& t) {
  return demangle_typename(typeid(t).name());
}
}  // namespace agl
