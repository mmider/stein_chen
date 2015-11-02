#include <Rcpp.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;
using namespace std;
/* I apologise for many inconsistencies within this code
 * such as sometimes using vector<int> and sometimes arrays
 * or sometimes using std:: and sometimes accesing objects directly
 * without explicit reference - I am still trying to figure these out myself
 */
 
int no_k_bdays(vector<int> bday_people_count, int k){
  /* test function - checks whether there are no k people
	 * with the same birth day in the room. If there aren't
	 * it returns 1, else returns 0
	 */
	int n = bday_people_count.size();
	for (int i = 0; i < n; i++){
		if (bday_people_count[i]>k-1){
			return 0;
		}
	}
	return 1;
}

double mean(vector<int> vec){
	// calculate mean of elements in a vector
  int n = vec.size();
	int total = 0;
	for (int i = 0; i < n; i++){
		total += vec[i];
	}
	// both total and size are ints, so I need to cast total to double first
	// in order to get the right answer from division
	double average = ((double) total)/n;
	return average;
}

vector<int> sample_indices(int size){
	/* given a vector (0,1,2,...,size-1) it samples from this vector
	 * uniformly with replacement. This function will be used in a bootstrap
	 */
	vector<int> sampled_indices;
	for (int i = 0; i < size; i++){
		// sample uniformly from the vector
		int index = floor(((double) rand()) / RAND_MAX * size);
		// store sampled value
		sampled_indices.push_back(index);
	}
	return sampled_indices;
}

vector<double> bootstrap_mean(vector<int> vec, int iter_num){
	/* This function simulates bootstrap samples for a given vector and
	 * for each of those samples calculates the mean. It then returns
	 * the vector of means
	 */
  int n = vec.size();
	vector<double> boot_sample_mean;
	for (int i = 0; i < iter_num;i++){
		// sample the indices that should be used for bootstrap
		vector<int> indices = sample_indices(n);
		vector<int> sample;
		for (int j = 0; j < n; j++){
			// create a boot sample vactor, which uses elements form vec, whose
			// indices weere drawn at random and were stored in "indices" vector
			sample.push_back(vec[indices[j]]);
		}
		// calculate mean for a given sample
		double sample_mean = mean(sample);
		// and store it
		boot_sample_mean.push_back(sample_mean);
	}
	return boot_sample_mean;
}

vector<double> bootstrap_mean_ci(vector<int> vec, const int iter_num, double certainty_level){
	// calculate the bootstrap confidence intervals for a given vector
	// find the index of the lower percentile from CI
	int lower_index = max(floor((1-certainty_level)/2 * iter_num)-1,0.0);
	// find the index of the upper percentile from CI
	int upper_index = ceil((0.5 + certainty_level/2) * iter_num)-1;
	// find the vector of bootstrapped sample means
	vector<double> boot_sample_mean = bootstrap_mean(vec, iter_num);
	// sort the vector of means in ascending order
	std::sort(boot_sample_mean.begin(), boot_sample_mean.end());
	// copy back to vector<double>
	vector<double> boot_conf_interval;
	// lower confidence bound is just lower_index'th element from the ordered mean vector
	boot_conf_interval.push_back(boot_sample_mean[lower_index]);
	// analogous for the upper confidence bound
	boot_conf_interval.push_back(boot_sample_mean[upper_index]);
	return boot_conf_interval;
}


//[[Rcpp::export]]
NumericVector Bday_MC(NumericVector out, int iterations = 1e4,
int num_days_in_year = 365, int room_size = 187, int iter_num = 100,
int bday_coincidences = 4,double confidence_level = 0.95)
{	
	/* IMPORTANT : you have to pass a vector "out" of length 3. silly, but I am not sure how to do it otherwise. 
   * declare vector of indices, where for each room there will be an indicator
	 * telling whether a particular event tested with "statistics" fuction
	 * has happened in a given room
	 */
	vector<int> statistics;
	for (int i = 0; i< iterations;i++)
	{
		/* declare people in a room and their birthdays. This is slightly different
		 * than in R function, in a sense that we will arange the data into favourable
		 * format from the onset. As we will be letting new people into room we will be checking 
		 * if it has a birthday on the same day as any person that is already in a room
		 * if it has than we will just increment the vector "people" at the position of
		 * corresponding birthday. Otherwise we will add new birthday to our "bday" vector
		 * and add new element equal to 1 to "people" vector, signifying that there is now
		 * 1 person in a room with new birthday. In the end we will  have bday listing all
		 * birthdays for people in a room and vector "people" which will give counts of how
		 * many people in a room share this given birthday. 
		 */
		vector<int> people;
		vector<int> bdays;
		// current size of the vectors above
		int current_size = 0;
		for (int j = 0; j < room_size; j++)
		{
			// simulate bday of the person arriving to a room
			int bday = floor(((double) rand()) / RAND_MAX * num_days_in_year);
			// declare test variable remembering if new person has birthday that
			// someone in a room already has
			bool added = false;
			// for all birthdays of people already in a room
			for (int k = 0; k < current_size; k++){
				// check if new person share the birthday
				if (bdays[k] == bday){
					// if it does than increment number of people that share this birthday
					++people[k];
					// someone has this bday already
					added = true;
					// no need to search any further
					break;
				}
			}
			// if no person in a room had this bday than add new entry
			if (!added == 1){
				// add new bday
				bdays.push_back(bday);
				// add count to number of people sharing this bday (1)
				people.push_back(1);
				// both vectors increased in size by 1
				++current_size;
			}
		}
		// by this point the room is populated, so let us check if the event has happened,
		// store result
		statistics.push_back(no_k_bdays(people, bday_coincidences));
	}
	// calculate the CI
	vector<double> conf_interval = bootstrap_mean_ci(statistics, iter_num, confidence_level);
  
  out[0] = conf_interval[0];
  out[1] = mean(statistics);
  out[2] = conf_interval[1];
	return out;
}
