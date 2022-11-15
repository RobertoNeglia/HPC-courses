#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

template <typename T> bool eq(const T &lhs, const T &rhs) {
  if (lhs.size() != rhs.size()) {
    return false;
  }
  for (auto i = lhs.begin(), j = rhs.begin(); i != lhs.end() && j != rhs.end();
       ++i, ++j) {
    if (*i != *j) {
      return false;
    }
  }
  return true;
}

void print_passed(bool b) {
  std::cout << (b ? "PASSED" : "FAILED") << std::endl;
}

int main() {
  // --------------------------------------------------------------------------
  // ASSIGMENT 1
  // sort and remove the duplicates from the given vector v
  // use just std::set
  const std::vector<int> v = {7, 8, 9, 2, 0, 7, 6, 4, 1, 1, 0, 4, 5, 6, 7, 3};
  std::set<int> v_set(v.begin(), v.end());
  std::vector<int> sorted_and_unique(v_set.begin(), v_set.end());
  print_passed(eq(sorted_and_unique, {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}}));
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  // ASSIGMENT 2
  // given a tokenized sentence in the vector words, compute
  // how many times each word appears, write the result in a std::map
  const std::vector<std::string> words = {
      "and",    "i",      "will", "show",   "you",     "something", "different",
      "from",   "either", "your", "shadow", "at",      "morning",   "striding",
      "behind", "you",    "or",   "your",   "shadow",  "at",        "evening",
      "rising", "to",     "meet", "you",    "i",       "will",      "show",
      "you",    "fear",   "in",   "a",      "handful", "of",        "dust"};
  std::map<std::string, size_t> word_count;
  // WRITE HERE YOUR SOLUTION
  for (const auto s : words) {
    auto i = word_count.find(s);
    if (i == word_count.end())
      // Element is not already present in the map, so i add it
      word_count.insert({s, 1});
    else
      // Element is in the map, i increment the value
      i->second++;
  }

  // PROF SOLUTION
  // The [] operator is used: if a key is not present, the operator inserts it
  // automatically
  //   for(const auto &w : words)
  //     {
  //       word_count[w]++;
  //     }
  print_passed(
      eq(word_count,
         {{"a", 1},         {"and", 1},      {"at", 2},      {"behind", 1},
          {"different", 1}, {"dust", 1},     {"either", 1},  {"evening", 1},
          {"fear", 1},      {"from", 1},     {"handful", 1}, {"i", 2},
          {"in", 1},        {"meet", 1},     {"morning", 1}, {"of", 1},
          {"or", 1},        {"rising", 1},   {"shadow", 2},  {"show", 2},
          {"something", 1}, {"striding", 1}, {"to", 1},      {"will", 2},
          {"you", 4},       {"your", 2}}));
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  // ASSIGMENT 3
  // elvaluate the functions specified in the vector funs at the point x=0.9
  // instert the results fun(0.9) into the vector fx
  // exploit the given unordered map to convert the string into an std::function
  // you will need to add elements to the given fmap
  const std::unordered_map<std::string, std::function<double(double)>> fmap = {
      {"cos", [](auto x) { return std::cos(x); }},
      {"tan", [](auto x) { return std::tan(x); }},
      {"sin", [](auto x) { return std::sin(x); }},
      {"exp", [](auto x) { return std::exp(x); }}};
  const std::vector<std::string> funs = {"tan", "sin", "exp"};
  std::vector<double> fx;
  for (const auto f : funs) {
    fx.emplace_back(fmap.find(f)->second(0.9));
  }
  // fmap.find() returns the pair (key -> value) with the specified key

  // PROF SOLUTION
  // fmap.at() returns only the value corresponding to the specified key
  // for (auto const f : funs) {
  //   fx.emplace_back(fmap.at(f)(0.9));
  // }

  print_passed(eq(fx, {{std::tan(0.9), std::sin(0.9), std::exp(0.9)}}));
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  // ASSIGMENT 4
  // compute the square of each even number between 0 and 11 starting from
  // the vector range10, employ the STL functions:
  // * copy_if
  // * for_each
  std::vector<int> range11(11);
  std::iota(range11.begin(), range11.end(), 0);
  std::vector<int> range11even_squares(6);
  std::copy_if(range11.begin(), range11.end(), range11even_squares.begin(),
               [](const auto x) { return (x % 2 == 0); });
  std::for_each(range11even_squares.begin(), range11even_squares.end(),
                [](auto &x) { x = x * x; });
  print_passed(eq(range11even_squares, {{0, 4, 16, 36, 64, 100}}));
  // --------------------------------------------------------------------------
  return 0;
}