#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>
#include <tuple>
#include <vector>

using namespace std;

float NODF(const vector<vector<int>> &A,
           const tuple<vector<int>, const vector<int>> &degrees) {
  int n = A[0].size();
  int m = A.size();
  vector<int> row_degree = get<0>(degrees);
  vector<int> col_degree = get<1>(degrees);

  double N_row = .0;
  for (size_t i = 0; i < m - 1; i++) {
    if (!row_degree[i]) {
      continue;
    }
    for (size_t j = i + 1; j < m; j++) {
      int overlap = 0;

      if (!row_degree[j] || row_degree[i] == row_degree[j]) {
        continue;
      }

      for (size_t k = 0; k < n; k++) {
        if (A[i][k] && A[j][k]) {
          overlap++;
        }
      }
      N_row += (double)overlap / min(row_degree[i], row_degree[j]);
    }
  }

  double N_col = .0;
  for (size_t i = 0; i < n - 1; i++) {
    if (!col_degree[i]) {
      continue;
    }
    for (size_t j = i + 1; j < n; j++) {
      int overlap = 0;

      if (!col_degree[j] || col_degree[i] == col_degree[j]) {
        continue;
      }

      for (size_t k = 0; k < m; k++) {
        if (A[k][i] && A[k][j]) {
          overlap++;
        }
      }
      N_col += (double)overlap / min(col_degree[i], col_degree[j]);
    }
  }

  for (auto &&i : row_degree) {
    if (i == 0) {
      m--;
    }
  }
  for (auto &&i : col_degree) {
    if (i == 0) {
      n--;
    }
  }

  double K = (n * (n - 1)) / 2.0 + (m * (m - 1)) / 2.0;
  return 100.0 * (N_row + N_col) / K;
}

tuple<vector<int>, vector<int>> degree_sequence(const vector<vector<int>> &A) {
  int n = A[0].size();
  int m = A.size();
  vector<int> row(m), col(n);

  for (size_t i = 0; i < m; i++) {
    int degree = 0;
    for (size_t j = 0; j < n; j++) {
      if (A[i][j]) {
        degree++;
      }
    }
    row[i] = degree;
  }
  for (size_t j = 0; j < n; j++) {
    int degree = 0;
    for (size_t i = 0; i < m; i++) {
      if (A[i][j]) {
        degree++;
      }
    }
    col[j] = degree;
  }

  return make_tuple(row, col);
}

double mean(vector<double> N) {
  return accumulate(N.begin(), N.end(), .0) / N.size();
}

double stdev(vector<double> N, double N_mean) {
  vector<double> x(N);
  for (size_t i = 0; i < x.size(); i++) {
    x[i] = pow(x[i] - N_mean, 2.0);
  }

  return sqrt(mean(x));
}

tuple<vector<double>, vector<double>, vector<double>, vector<double>, double> null_model(vector<vector<int>> A,
                                                                                  int samples = 1000,
                                                                                  short type = 0) {
  vector<double> N_sample(samples);
  int n = A[0].size();
  int m = A.size();
  vector<double> c_row(m), c_col(n), node_N_mean(m+n), node_N_std(m+n);
  tuple<vector<int>, vector<int>> degrees = degree_sequence(A);
  vector<int> row_degrees = get<0>(degrees);
  vector<int> col_degrees = get<1>(degrees);
  double N_mean, N_std;
  double N = NODF(A, degrees);
  mt19937_64 generator;
  uniform_real_distribution<double> distribution(0, 1);
  vector<int> unmodified_row(n), unmodified_col(m);
  vector<int> current_column(m); // For null model type 1

  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      unmodified_row[j] = A[i][j];
    }
    for (size_t k = 0; k < samples; k++) {
      if (type == 0) {
        for (size_t j = 0; j < n; j++) {
          double p =
              ((double)row_degrees[i] / n + (double)col_degrees[j] / m) / 2.0;
          if (distribution(generator) < p) {
            A[i][j] = 1;
          } else {
            A[i][j] = 0;
          }
        }
      } else if (type == 1) {
        shuffle(A[i].begin(), A[i].end(), generator);
      }
      N_sample[k] = NODF(A, degree_sequence(A));
    }
    for (size_t j = 0; j < n; j++) {
      A[i][j] = unmodified_row[j];
    }

    N_mean = mean(N_sample);
    N_std = stdev(N_sample, N_mean);
    c_row[i] = (N - N_mean) / N_std;
    node_N_mean[i] = N_mean;
    node_N_std[i] = N_std;
  }

  for (size_t j = 0; j < n; j++) {
    for (size_t i = 0; i < m; i++) {
      unmodified_col[i] = A[i][j];
      current_column[i] = A[i][j];
    }
    for (size_t k = 0; k < samples; k++) {
      if (type == 0) {
        for (size_t i = 0; i < m; i++) {
          double p =
              ((double)row_degrees[i] / n + (double)col_degrees[j] / m) / 2.0;
          if (distribution(generator) < p) {
            A[i][j] = 1;
          } else {
            A[i][j] = 0;
          }
        }
      } else if (type == 1) {
        shuffle(current_column.begin(), current_column.end(), generator);
        
        for (size_t i = 0; i < m; i++) {
          A[i][j] = current_column[i];
        }
      }
      N_sample[k] = NODF(A, degree_sequence(A));
    }
    for (size_t i = 0; i < m; i++) {
      A[i][j] = unmodified_col[i];
    }

    N_mean = mean(N_sample);
    N_std = stdev(N_sample, N_mean);
    c_col[j] = (N - N_mean) / N_std;
    node_N_mean[m + j] = N_mean;
    node_N_std[m + j] = N_std;
  }
  return make_tuple(c_row, c_col, node_N_mean, node_N_std, N);
}

int main(int argc, char const *argv[]) {
  bool return_nodf = false, return_contributions = false;
  short null_model_type = 0;
  string line;
  string i_matrix;
  vector<vector<int>> A;

  for (size_t i = 1; i < argc; i++) {
    string arg = argv[i];
    if (arg == "--nodf") {
      return_nodf = true;
    } else if (arg == "--contributions") {
      return_contributions = true;
    } else if (arg == "-i") {
      i_matrix = argv[i + 1];
      i++;
    } else if (arg == "--null-model") {
      null_model_type = stoi(argv[i + 1]);
    }
  }
  cout << "Run parameters:" << endl;
  cout << "i_matrix: " << i_matrix << endl;
  cout << "return_nodf: " << return_nodf << endl;
  cout << "return_contributions: " << return_contributions << endl;
  cout << "null_model_type: ";
  if (null_model_type == 0)
  {
    cout << "0 (probabilistic)" << endl;
  } else if (null_model_type == 1)
  {
    cout << "1 (shuffle)" << endl;
  }

  ifstream matrix_in("matrix" + i_matrix + ".in.txt");

  while (getline(matrix_in, line)) {
    istringstream lin(line);
    int d;
    while (lin) {
      vector<int> row =
          vector<int>(istream_iterator<int>(lin), istream_iterator<int>());
      A.push_back(row);
    }
  }

  if (return_nodf) {
    ofstream matrix_out("matrix" + i_matrix + ".nodf.txt");
    matrix_out << NODF(A, degree_sequence(A)) << endl;
  }
  if (return_contributions) {
    ofstream matrix_out("matrix" + i_matrix + ".contributions.csv");
    matrix_out.precision(10);
    matrix_out << "type,contribution,N_obs,N_mean,N_std" << endl;
    tuple<vector<double>, vector<double>, vector<double>, vector<double>, double> res = null_model(A, 1000, null_model_type);
    double N_obs = get<4>(res);
    vector<double> N_mean = get<2>(res);
    vector<double> N_std = get<3>(res);
    int count = 0;

    for (auto &&i : get<0>(res)) {
      matrix_out << "row," << i << "," << N_obs << "," << N_mean[count] << "," << N_std[count] << endl;
      ++count;
    }
    for (auto &&i : get<1>(res)) {
      matrix_out << "column," << i << "," << N_obs << "," << N_mean[count] << "," << N_std[count] << endl;
      ++count;
    }
  }

  return 0;
}