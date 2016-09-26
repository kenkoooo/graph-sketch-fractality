// http://stackoverflow.com/questions/281818/unmangling-the-result-of-stdtype-infoname

#include "type.h"

#ifdef __GNUG__
#include <cstdlib>
#include <memory>
#include <cxxabi.h>
#endif

namespace agl {
#ifdef __GNUG__
std::string demangle_typename(const char* name) {
  int status = -4;
  std::unique_ptr<char, void (*)(void*)> res {
    abi::__cxa_demangle(name, NULL, NULL, &status), std::free };
  return (status == 0) ? res.get() : name;
}

#else

std::string demangle_typename(const char* name) {
  return name;
}

#endif
}  // namespace agl
