#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <array>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
const int output_length = 25;
const int results_size = 10000;
const int n_x = 200;
const int n_y = 200;
const int num_gen = 20;

//data

#define array_length(array) (sizeof((array))/sizeof((array[0])))
#define random_pick(array) (array[rand() % array_length(array)])
#define random_index(array) (rand() % array_length(array))

const double pi = 3.1415926535897;

inline double rand_doubleRange(double a, double b)
{
    return ((b - a) * ((double)rand() / RAND_MAX)) + a;
}

inline double min(double x,double y,double z){
double small = x;
if (y<small) 
small =  y;
if (z<small) 
small =  z;
return small;}

void compute(double result[output_length],int inputs[num_inputs]) {

//inputs


//equations


bool all_contraints = true;
for( int i = 0; i < array_length(constraints); i = i + 1 ) 
	all_contraints = all_contraints && constraints[i];

//results

for( int i = 0; i < output_length; i = i + 1 )  
	all_contraints = all_contraints && isfinite(result[i]);

result[0] = (double)(int) all_contraints;



return;

}





void model(double results[results_size][output_length]) {

int k = 0;

int inputs[num_inputs] = {0};
int new_inputs[num_inputs] = {0};

double result[output_length] = {0};
double min_x = 1000000;double max_x = -1000000;
double min_y = 1000000;double max_y = -1000000;
double min_z = 1000000;double max_z = -1000000;
double grid[n_x][n_y] = {{0}};
vector<int> selected_indices;


do {

//randoms

compute(results[k],inputs);

if (results[k][x_index]>max_x)
	max_x = results[k][x_index];

if (results[k][x_index]<min_x)
	min_x = results[k][x_index];

if (results[k][y_index]>max_y)
	max_y = results[k][y_index];

if (results[k][y_index]<min_y)
	min_y = results[k][y_index];

if (results[k][z_index]>max_z)
	max_z = results[k][z_index];

if (results[k][z_index]<min_z)
	min_z = results[k][z_index];

if ((int)results[k][0] == 1) {
   	k = k + 1;
}


}  while (k<results_size);

//cout <<"initialization done!" << "\n";



for (int gen = 0; gen < num_gen; ++gen)
{
	
for (int i = 0; i < results_size; ++i)
{

int x_ = floor(n_x*(results[i][x_index]-min_x)/(max_x-min_x));
int y_ = floor(n_y*(results[i][y_index]-min_y)/(max_y-min_y));



if (results[i][z_index]>grid[x_][y_])
{
	grid[x_][y_] = results[i][z_index];
	selected_indices.push_back(i);
}



}

//cout << "gen: " << gen << " , selection phase done!" << "\n";


k = 0;

do {

int i1 = rand() % selected_indices.size();
int i2 = rand() % selected_indices.size();


for (int i = 0; i < num_inputs; ++i)
{
new_inputs[i] = (int)rand_doubleRange(results[selected_indices[i1]][i+num_outputs],results[selected_indices[i1]][i+num_outputs]);
}




compute(results[k],new_inputs);


if (results[k][x_index]>max_x)
	max_x = results[k][x_index];

if (results[k][x_index]<min_x)
	min_x = results[k][x_index];

if (results[k][y_index]>max_y)
	max_y = results[k][y_index];

if (results[k][y_index]<min_y)
	min_y = results[k][y_index];

if (results[k][z_index]>max_z)
	max_z = results[k][z_index];

if (results[k][z_index]<min_z)
	min_z = results[k][z_index];

if ((int)results[k][0] == 1) {
   	k = k + 1;
}


}  while (k<(int)(0.7*results_size));

//cout << "generation: " << gen << "\n";


}




for (int i = 0; i < results_size; ++i)
{
//outputs
}


return;

}

