#include "macros.h"
#include "functions.h"

#include <thread>
#include <future>
#include <chrono>
#include <sys/time.h>
#include <sys/utsname.h>
#include <unistd.h>
using namespace std;

namespace agl {
vector<string> split(const string &str, char splitter) {
  vector<string> res;
  string tmp;
  for (char c : str) {
    if (c == splitter) {
      if (tmp != "") res.emplace_back(move(tmp));
      tmp = "";
    } else {
      tmp += c;
    }
  }
  if (tmp != "") res.emplace_back(move(tmp));
  return res;
}

string strip(const string &s) {
  size_t i = 0;
  while (i < s.length() && isspace(s[i])) ++i;
  size_t j = s.length();
  while (i < j && isspace(s[j - 1])) --j;
  return s.substr(i, j - i);
}

double get_current_time_sec() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

void execute_within_time_limit_or_die(double time_limit_sec,
                                      std::function<void()> task,
                                      std::function<void()> hook_before_death) {
  auto future = std::async(std::launch::async, task);
  double start_time_sec = get_current_time_sec();

  for (;;) {
    auto status = future.wait_for(std::chrono::milliseconds(0));
    if (status == std::future_status::ready) {
      break;
    }

    if (get_current_time_sec() - start_time_sec > time_limit_sec) {
      hook_before_death();
      FAIL_MSG("Time limit exceeded");
    }

    this_thread::sleep_for(chrono::seconds(1));
  }
}
}  // namespace agl
