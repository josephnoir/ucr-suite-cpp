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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <chrono>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

namespace {

// Pseudo infinite number for this code
constexpr double INF = 1e20;

// If serious error happens, terminate the program.
void error(const string& msg) {
  cout << "ERROR  : " << msg << "!!!"
       << "Usage  : ./uce_ed data_file query_file m" << endl
       << "Example: ./uce_ed data.txt  query.txt  128" << endl;
  exit(1);
}
  
} // namespace <anonymous>

int main(int argc, char *argv[]) {
  // Ordering of query by |z(q_i)|.
  vector<int> order;
  // best-so-far
  double best = INF;
  // Answer: location of the best-so-far match.
  size_t index = 0;

  double d = 0;
  double ex = 0;
  double ex2 = 0;

  auto start = chrono::system_clock::now();

  if (argc <= 3)
    error("Invalid number of arguments");

  // Stream from data file.
  fstream fp(argv[1], std::ios_base::in);
  if (!fp)
    error("Data file not found");

  // Stream from query file.
  fstream qp(argv[2], std::ios_base::in);
  if (!qp)
    error("Query file not found");

  // Length of query.
  size_t length = atoll(argv[3]);

  // Array for keeping the query data.
  vector<double> queries(length);

  // Read the query data from input file and calculate its statistic such
  // as mean, std.
  while (qp >> d && queries.size() < length) {
    ex += d;
    ex2 += d * d;
    queries.push_back(d);
  }
  qp.close();
  
  size_t mean = ex / length;
  size_t std = sqrt((ex2 / length) - mean * mean);

  // Do z_normalixation on query data.
  /*
  for (size_t i = 0; i < m; ++i)
    Q[i] = (Q[i] - mean) / std;
  */
  transform(begin(queries), end(queries), begin(queries),
            [mean, std](auto val) { return (val - mean) / std; });

  // Sort the query data.
  vector<pair<double, unsigned>> tmp(length);
  for (size_t i = 0; i < length; ++i) {
    tmp[i].first = queries[i];
    tmp[i].second = i;
  }
  // The query will be sorted by absolute z-normalization value,
  // |z_norm(Q[i])| from high to low.
  sort(tmp.begin(), tmp.end(), [](auto& lhs, auto& rhs) {
    return abs(rhs.first) - abs(lhs.first);
  });
  order.resize(length);
  for (size_t i = 0; i < length; ++i) {
    queries[i] = tmp[i].first;
    order[i] = tmp[i].second;
  }
  tmp.clear();

  // Array for keeping the current data. Twice the size for removing modulo
  // (circulation) in distance calculation.
  vector<double> data(2 * length);

  double dist = 0;
  size_t j = 0;
  size_t i = 0;
  ex = 0;
  ex2 = 0;

  // Main function for calculating ED distance between the query, Q, and
  // current data, T. Note that Q is already sorted by absolute
  // z-normalization value, |z_norm(Q[i])|.
  auto dist_fun = [&queries, &order, length](const vector<double>& T,
                                             const int& j,
                                             const double& mean,
                                             const double& std,
                                             const double& best) {
    double sum = 0;
    for (int i = 0; i < length && sum < best; ++i) {
      double x = (T[(order[i] + j)] - mean) / std;
      sum += (x - queries[i]) * (x - queries[i]);
    }
    return sum;
  };

  // Read data file, one value at a time.
  while (fp >> d) {
    ex += d;
    ex2 += d*d;
    data[i % length] = d;
    data[(i % length) + length] = d;

    // If there is enough data in T, the ED distance can be calculated.
    if (i >= length - 1) {
      // The current starting location of T.
      j = (i + 1) % length;

      // Z_norm(T[i]) will be calculated on the fly.
      mean = ex / length;
      std = ex2 / length;
      std = sqrt(std - mean * mean);

      // Calculate ED distance.
      dist = dist_fun(data, j, mean, std, best);
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

  cout << "Location: " << index << endl;
  cout << "Distance: " << sqrt(best) << endl;
  cout << "Data Scanned: " << i << endl;
  cout << "Total Execution Time: "
       << chrono::duration_cast<chrono::seconds>(end - start).count()
       << " secs" << endl;
}
