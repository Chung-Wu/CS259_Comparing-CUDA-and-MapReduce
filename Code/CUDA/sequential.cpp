#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>

using namespace std;

string sliding_window(const string &s, int start, int size) {
  return s.substr(start, size);
}

int main() {
  ifstream file("../../Data/hg19.fa");
  unordered_map<string, int> map;
  if (!file.is_open()) {
    cerr << "Error: cannot open file.\n";
    return 1;
  }

  string line, data;
  while (getline(file, line)) {
    if (!line.empty() && line[0] != '>')
      data += line; // concatenate sequence lines
  }

  int k = 3; // example window size
  for (int i = 0; i + k <= (unsigned long)data.size(); ++i) {
    string kmer = sliding_window(data, i, k);
    map[kmer] += 1;
  }

  return 0;
}
