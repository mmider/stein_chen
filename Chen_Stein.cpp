#include <iostream>
#include <random>
#include <chrono>
#include <cmath>
#include <stdlib.h>
using namespace std;
/* I apologise for many inconsistencies within this code
 * such as sometimes using vector<int> and sometimes arrays
 * or sometimes using std:: and sometimes accesing objects directly
 * without explicit reference - I am still trying to figure these out myself
 */
 
int no_k_bdays(vector<int> bday_people_count, const int k){
	/* test function - checks whether there are no k people
	 * with the same birth day in the room. If there aren't
	 * it returns 1, else returns 0
	 */
	int n = bday_people_count.size();
	for (int i = 0; i < n; i++){
		if (bday_people_count[i]>k-1){
			return(0);
		}
	}
	return 1;
}

int compare(const void * a, const void * b){
	/* helper function for the qsort - quicksort
	 * I will want to sort an array at some point and I need to supply
	 * a function which determines the order for objects to be sorted
	 * in my case I want the array to be in an increasing order
	 */
	if (*(double*)a > *(double*)b) return 1;
	else if (*(double*)a < *(double*)b) return -1;
	else return 0;
}

double mean(int vec[], int size){
	// calculate mean of elements in a vector
	int total = 0;
	for (int i = 0; i < size; i++){
		total += vec[i];
	}
	// both total and size are ints, so I need to cast total to double first
	// in order to get the right answer from division
	double average = ((double) total)/size;
	return average;
}

vector<int> sample_indices(int size){
	/* given a vector (0,1,2,...,size-1) it samples from this vector
	 * uniformly with replacement. This function will be used in a bootstrap
	 */
	 // set a random seed using current time
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// declare a random number engine and set its seed to our random seed from above
	std::default_random_engine generator(seed);
	// declare uniform distribution on (0,1,2,...,size-1)
	std::uniform_int_distribution<int> distribution(0,size-1);
	// declare vector, where we will store our samples
	vector<int> sampled_indices;
	for (int i = 0; i < size; i++){
		// sample uniformly from the vector
		int index = distribution(generator);
		// store sampled value
		sampled_indices.push_back(index);
	}
	return sampled_indices;
}

vector<double> bootstrap_mean(int vec[], const int size, int iter_num){
	/* This function simulates bootstrap samples for a given vector and
	 * for each of those samples calculates the mean. It then returns
	 * the vector of means
	 */
	vector<double> boot_sample_mean;
	for (int i = 0; i < iter_num;i++){
		// sample the indices that should be used for bootstrap
		vector<int> indices = sample_indices(size);
		int sample[size];
		for (int j = 0; j < size; j++){
			// create a boot sample vactor, which uses elements form vec, whose
			// indices weere drawn at random and were stored in "indices" vector
			sample[j] = vec[indices[j]];
		}
		// calculate mean for a given sample
		double sample_mean = mean(sample, size);
		// and store it
		boot_sample_mean.push_back(sample_mean);
	}
	return boot_sample_mean;
}

vector<double> bootstrap_mean_ci(int vec[], const int size, const int iter_num, double certainty_level){
	// calculate the bootstrap confidence intervals for a given vector
	// find the index of the lower percentile from CI
	int lower_index = max(floor((1-certainty_level)/2 * iter_num)-1,0.0);
	// find the index of the upper percentile from CI
	int upper_index = ceil((0.5 + certainty_level/2) * iter_num)-1;
	// find the vector of bootstrapped sample means 	
	vector<double> boot_sample_mean = bootstrap_mean(vec, size, iter_num);
	// quite stupid, but I need to change vector<double> to array of doubles
	// before I can use qsort... that's the price you pay if you use bad programming
	double temp[iter_num];
	for (int i = 0; i < iter_num; i++){
		temp[i] = boot_sample_mean[i];
	}
	// sort the vector of means in ascending order
	qsort(temp, iter_num,sizeof(double), compare);
	// copy back to vector<double>
	for (int i = 0; i < iter_num; i++){
		boot_sample_mean[i] = temp[i];
	}
	vector<double> boot_conf_interval;
	// lower confidence bound is just lower_index'th element from the ordered mean vector
	boot_conf_interval.push_back(boot_sample_mean[lower_index]);
	// analogous for the upper confidence bound
	boot_conf_interval.push_back(boot_sample_mean[upper_index]);
	return boot_conf_interval;
}


int main()
{
	// set the constants for the problem, feel free to play with those
	const int iterations = 100000;
	const int lower_bound = 0;
	const int upper_bound = 364;
	const int room_size = 187;
	const int iter_num = 100;
	const int bday_coincidences = 4;
	const double confidence_level = 0.95;
	
	// set the seed using clock (I know that I have set the seed somewhere else as well,
	// but this function has a different scope, so I need a separate seed here (I think)
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// declare random number engine
	std::default_random_engine generator(seed);
	// declare uniform distribution over days in a year
	std::uniform_int_distribution<int> distribution(lower_bound,upper_bound);
	
	/* declare vector of indices, where for each room there will be an indicator
	 * telling whether a particular event tested with "statistics" fuction
	 * has happened in a given room
	 */
	int statistics[iterations];
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
			int bday = distribution(generator);
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
		statistics[i] = no_k_bdays(people, bday_coincidences);
	}
	// calculate the CI
	vector<double> conf_interval = bootstrap_mean_ci(statistics, iterations, iter_num, confidence_level);
	
	// print results
	std::cout << ((int) (confidence_level*100))<<" \% level confidence interval: (" << conf_interval[0]  << ", " << conf_interval[1] << ")"<< std::endl;
	std::cout << "estimated mean: " << mean(statistics, iterations) << std::endl;
	
	// wait for the user input (used so that the command window won't close immediately
	cin.get();
	
	return 0;
}
