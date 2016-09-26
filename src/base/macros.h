#ifndef SRC_COMMON_MACROS_H_
#define SRC_COMMON_MACROS_H_

#define CHECK(expr)                                                     \
  if (expr) {                                                           \
  } else {                                                              \
    fprintf(stderr, "CHECK Failed (%s:%d): %s\n",                       \
            __FILE__, __LINE__, #expr);                                 \
    exit(EXIT_FAILURE);                                                 \
  }

#define CHECK_MSG(expr, msg)                                            \
  if (expr) {                                                           \
  } else {                                                              \
    fprintf(stderr, "CHECK Failed (%s:%d): %s (%s)\n",                  \
            __FILE__, __LINE__, msg, #expr);                            \
    exit(EXIT_FAILURE);                                                 \
  }

#define CHECK_PERROR(expr)                                                     \
  if (expr) {                                                           \
  } else {                                                              \
    fprintf(stderr, "CHECK Failed (%s:%d): %s: ",                       \
            __FILE__, __LINE__, #expr);                                 \
    perror(nullptr); \
    exit(EXIT_FAILURE);                                                 \
  }

#define FAIL_PERROR() \
    do {              \
      fprintf(stderr, "Error (%s:%d): ", __FILE__, __LINE__); \
      perror(nullptr); \
      exit(EXIT_FAILURE);                                                 \
    } while (0)

#define FAIL_MSG(msg) \
    do {              \
      fprintf(stderr, "Error (%s:%d): %s\n", __FILE__, __LINE__, std::string(msg).c_str()); \
      exit(EXIT_FAILURE);                                                 \
    } while (0)


#define DISALLOW_COPY_AND_ASSIGN(TypeName)  \
  TypeName(const TypeName&);                \
  void operator=(const TypeName&)

#endif /* SRC_COMMON_MACROS_H_ */
