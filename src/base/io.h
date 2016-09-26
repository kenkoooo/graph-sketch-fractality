#pragma once
#include <cstdio>
#include <iostream>
#include <locale>
#include <string>
#include "macros.h"

namespace agl {
inline bool read_text_eof(FILE *fp) {
  for (;;) {
    int c = getc(fp);
    if (c == EOF) {
      return true;
    } else if (isspace(c)) {
      continue;
    } else {
      ungetc(c, fp);
      return false;
    }
  }
  FAIL_MSG("format error, EOF expected");
}

template<typename T = int>
inline long read_text_integer(FILE *fp) {
  int c, s;
  T x;
  while (!isdigit(c = getc(fp)) && c != '-');
  if (c == '-') {
    s = -1;
    x = 0;
  } else {
    s = 1;
    x = c - '0';
  }
  while (isdigit(c = getc(fp))) {
    x = (x * 10 + (c - '0'));
  }
  return s * x;
}

template<typename T>
void read_binary(std::istream &is, T *t) {
  CHECK_PERROR(is.read((char*)t, sizeof(T)));
}

template<typename T>
void write_binary(std::ostream &os, const T &t) {
  CHECK_PERROR(os.write((char*)&t, sizeof(T)));
}
}  // namespace agl
