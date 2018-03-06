/***********************************************************************/
/************************* DISCLAIMER **********************************/
/***********************************************************************/
/** This UCR Suite software is copyright protected � 2012 by          **/
/** Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,            **/
/** Gustavo Batista and Eamonn Keogh.                                 **/
/**                                                                   **/
/** Unless stated otherwise, all software is provided free of charge. **/
/** As well, all software is provided on an "as is" basis without     **/
/** warranty of any kind, express or implied. Under no circumstances  **/
/** and under no legal theory, whether in tort, contract,or otherwise,**/
/** shall Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,      **/
/** Gustavo Batista, or Eamonn Keogh be liable to you or to any other **/
/** person for any indirect, special, incidental, or consequential    **/
/** damages of any character including, without limitation, damages   **/
/** for loss of goodwill, work stoppage, computer failure or          **/
/** malfunction, or for any and all other damages or losses.          **/
/**                                                                   **/
/** If you do not agree with these terms, then you you are advised to **/
/** not use this software.                                            **/
/***********************************************************************/
/***********************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

using namespace std;

namespace {

/// If serious error happens, terminate the program.
void error(const string& msg) {
  cout << "ERROR  : " << msg << "!!!"
       << "Usage  : ./uce_ed data_file query_file m" << endl
       << "Example: ./uce_ed data.txt  query.txt  128" << endl;
  exit(1);
}

/// Main function for calculating ED distance between the query, Q, and
/// current data, T. Note that Q is already sorted by absolute
/// z-normalization value, |z_norm(Q[i])|.
double ed_distance(const vector<double>& query,
                   const vector<double>& data,
                   const int j,
                   const size_t length,
                   const double mean,
                   const double std,
                   const vector<unsigned>& order,
                   const double best) {
  double sum = 0;
  for (size_t i = 0; i < length && sum < best; ++i) {
    double x = (data[(order[i] + j)] - mean) / std;
    sum += (x - query[i]) * (x - query[i]);
  }
  return sum;
};

/// Sorts `vec` and returns the new index order.
template <class T>
vector<unsigned> sort_get_order(vector<T>& vec) {
  auto length = vec.size();
  vector<pair<double, unsigned>> tmp(length);
  for (size_t i = 0; i < length; ++i) {
    tmp[i].first = vec[i];
    tmp[i].second = i;
  }
  // The query will be sorted by absolute z-normalization value,
  // |z_norm(Q[i])| from high to low.
  //sort(tmp.begin(), tmp.end(),
       //[](auto& lhs, auto& rhs) { return abs(rhs.first) > abs(lhs.first); });
  // Although this has no effect on the sort time itself, using quicksort
  // halfs the runtime ...
  qsort(tmp.data(), length, sizeof(pair<double, unsigned>),
        [](const void* lhs, const void* rhs) -> int {
    auto x = reinterpret_cast<const pair<double, unsigned>*>(lhs);
    auto y = reinterpret_cast<const pair<double, unsigned>*>(rhs);
    return abs(y->first) - abs(x->first);
  });
  // Ordering of query by |z(q_i)|.
  vector<unsigned> order(length);
  for (size_t i = 0; i < length; ++i) {
    vec[i] = tmp[i].first;
    order[i] = tmp[i].second;
  }
  return order;
}

}  // namespace

int main(int argc, char* argv[]) {
  // best-so-far
  double best = numeric_limits<double>::max();

  // Answer: location of the best-so-far match.
  size_t index = 0;

  double d = 0;
  double ex = 0;
  double ex2 = 0;

  auto start = chrono::system_clock::now();

  if (argc <= 2)
    error("Invalid number of arguments");

  // Stream from data file.
  ifstream fp(argv[1]);
  if (!fp)
    error("Data file not found");

  // Stream from query file.
  ifstream qp(argv[2]);
  if (!qp)
    error("Query file not found");

  // Length of query.
  size_t length = atoll(argv[3]);

  // Array for keeping the query data.
  vector<double> query;
  query.reserve(length);

  // Read the query data from input file and calculate its statistic such
  // as mean, std.
  while (qp >> d && query.size() < length) {
    ex += d;
    ex2 += d * d;
    query.push_back(d);
  }
  qp.close();

  double mean = ex / length;
  double std = sqrt((ex2 / length) - mean * mean);

  // Do z_normalixation on query data.
  //transform(begin(query), end(query), begin(query),
            //[mean, std](auto val) { return (val - mean) / std; });
  for (auto& val : query)
    val = (val - mean) / std;

  // Sort the query data.
  auto order = sort_get_order(query);

  // Array for keeping the current data. Twice the size for removing modulo
  // (circulation) in distance calculation.
  vector<double> data(2 * length);

  double dist = 0;
  size_t j = 0;
  size_t i = 0;
  ex = 0;
  ex2 = 0;

  auto mid = chrono::system_clock::now();

  // Read data file, one value at a time.
  while (fp >> d) {
    ex += d;
    ex2 += d * d;
    data[i % length] = d;
    data[(i % length) + length] = d;

    // If there is enough data in T, the ED distance can be calculated.
    if (i >= length - 1) {
      // The current starting location of T.
      j = (i + 1) % length;

      // Z_norm(T[i]) will be calculated on the fly.
      mean = ex / length;
      std = sqrt((ex2 / length) - mean * mean);

      // Calculate ED distance.
      dist = ed_distance(query, data, j, length, mean, std, order, best);
      if (dist < best) {
        best = dist;
        index = i - length + 1;
      }
      ex -= data[j];
      ex2 -= data[j] * data[j];
    }
    ++i;
  }
  fp.close();
  auto end = chrono::system_clock::now();

  cout << "Location: " << index << endl
       << "Distance: " << sqrt(best) << endl
       << "Data Scanned: " << i << endl
       << "Partial Execution Time: "
       << chrono::duration_cast<chrono::milliseconds>(mid - start).count()
       << " ms" << endl
       << "Total Execution Time: "
       << chrono::duration_cast<chrono::milliseconds>(end - start).count()
       << " ms" << endl;
}
