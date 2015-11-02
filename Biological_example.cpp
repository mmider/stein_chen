#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
using namespace std;

int find_largest_one_run(vector<int> coincidences, int window){
  int largest = 0;
  int n = coincidences.size();
  
  for (int i = 0; i < window; i++){
    largest += coincidences[i];
  }
  int current = largest;
  for (int i = window; i < n; i++){
    current = current + coincidences[i] - coincidences[i - window];
    if (current > largest){
      largest = current;
    }
  }
  return largest;
}


// [[Rcpp::export]]
int longest_common_run(CharacterVector a, CharacterVector b, int cap){
  int n = a.size();
  int m = b.size();
  if (n != m){
    return 0;
  }
  if (n < cap || m < cap)
    return 0;

  int longest = 0;
  for (int i = 0; i < n-cap+1; i++){
    vector<int> coincidences;
    for (int j = 0; j < n-i; j++){
      int comparison = 0;
      if (a[i + j] == b[j])
        comparison = 1;
      coincidences.push_back(comparison);
    }
    int longest_current = find_largest_one_run(coincidences, cap);
    if (longest_current > longest)
      longest = longest_current;
  }
  return longest;
}