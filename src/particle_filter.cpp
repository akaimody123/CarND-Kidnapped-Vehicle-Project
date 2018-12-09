/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

#define EPS 0.00001

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	//initialize with 100 particles
	num_particles=100;
	//create normal distribution
	normal_distribution<double> dist_x(x,std[0]);
	normal_distribution<double> dist_y(y,std[1]);
	normal_distribution<double> dist_theta(theta,std[2]);

	for(int i=0;i<num_particles;i++)
	{
		Particle particle;
		particle.id=i;
		particle.x=dist_x(gen);
		particle.y=dist_y(gen);
		particle.theta=dist_theta(gen);
		particle.weight=1;

		particles.push_back(particle);
	}

	is_initialized=true;
	//std::cout<<"ini ok"<<std::endl;



}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/



	for(int i=0;i<num_particles;i++)
	{
		
		if (fabs(yaw_rate)<0.001)
		{
			//std::cout<<"mean"<<particles[i].x+velocity*cos(particles[i].theta)*delta_t<<std::endl;
			normal_distribution<double> pred_x(particles[i].x+velocity*cos(particles[i].theta)*delta_t,std_pos[0]);
			normal_distribution<double> pred_y(particles[i].y+velocity*sin(particles[i].theta)*delta_t,std_pos[1]);
			//std::cout<<"dis ok"<<std::endl;

			particles[i].x=pred_x(gen);
			//std::cout<<"genok"<<std::endl;
			particles[i].y=pred_y(gen);

		}
		else
		{
			//std::cout<<"mean"<<particles[i].x+velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta))<<std::endl;

			normal_distribution<double> pred_x(particles[i].x+velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta)),std_pos[0]);
			normal_distribution<double> pred_y(particles[i].y+velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t)),std_pos[1]);
			//std::cout<<"dis ok"<<std::endl;
			particles[i].x=pred_x(gen);
			//std::cout<<"genok"<<std::endl;
			particles[i].y=pred_y(gen);
		
		}
		normal_distribution<double> pred_theta(particles[i].theta+yaw_rate*delta_t,std_pos[2]);

		particles[i].theta=pred_theta(gen);
		//std::cout<<"x"<<particles[i].x<<" y"<<particles[i].y<<" theta"<<particles[i].theta<<std::endl;


	
	}
	//std::cout<<"prediction ok"<<std::endl;

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	// Initialize min distance as a really big number.
	

	for(int i=0;i<observations.size();i++)
	{
		double min_dist = numeric_limits<double>::max();
		int mapID=-1;

		for(int j=0;j<predicted.size();j++)
		{


			double sum=sqrt((observations[i].x-predicted[j].x)*(observations[i].x-predicted[j].x)+(observations[i].y-predicted[j].y)*(observations[i].y-predicted[j].y));
			if (sum<min_dist)
			{
				min_dist=sum;
				mapID=predicted[j].id;
			}

		}

		observations[i].id=mapID;
		

		
	}
	//std::cout<<"dataAsso ok"<<std::endl;
	

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	

    //std::cout<<"2"<<std::endl;

	//find particles inside the sensor range
	for(int i=0;i<num_particles;i++)
	{
		vector<LandmarkObs> mappedObservations;
		vector<LandmarkObs> mappedLandmarks;

		for(int j=0;j<map_landmarks.landmark_list.size();j++)
		{

			if((particles[i].x-map_landmarks.landmark_list[j].x_f)*(particles[i].x-map_landmarks.landmark_list[j].x_f)+(particles[i].y-map_landmarks.landmark_list[j].y_f)*(particles[i].y-map_landmarks.landmark_list[j].y_f)<sensor_range*sensor_range)
			{
				mappedLandmarks.push_back(LandmarkObs{ map_landmarks.landmark_list[j].id_i,map_landmarks.landmark_list[j].x_f,map_landmarks.landmark_list[j].y_f});
			}
		}
	


		for(int j=0;j<observations.size();j++)
		{


			
			double xx=cos(particles[i].theta)*observations[j].x-sin(particles[i].theta)*observations[j].y+particles[i].x;
			double yy=sin(particles[i].theta)*observations[j].x+cos(particles[i].theta)*observations[j].y+particles[i].y;
			mappedObservations.push_back(LandmarkObs{observations[j].id, xx, yy });
		}

		dataAssociation(mappedLandmarks,mappedObservations);

		// Reseting weight.
	    particles[i].weight = 1.0;
	    // Calculate weights.
	    for(unsigned int j = 0; j < mappedObservations.size(); j++) {
	      double observationX = mappedObservations[j].x;
	      double observationY = mappedObservations[j].y;

	      int landmarkId = mappedObservations[j].id;

	      double landmarkX, landmarkY;
	      unsigned int k = 0;
	      unsigned int nLandmarks = mappedLandmarks.size();
	      bool found = false;

	      while( !found && k < nLandmarks ) {
	        if ( mappedLandmarks[k].id == landmarkId) {
	          found = true;
	          landmarkX = mappedLandmarks[k].x;
	          landmarkY = mappedLandmarks[k].y;
	        }
	        k++;
	      }

	      // Calculating weight.
	     

	      double weight = ( 1/(2*M_PI*std_landmark[0]*std_landmark[1])) * exp( -(( observationX - landmarkX)*( observationX - landmarkX)/(2*std_landmark[0]*std_landmark[0]) + ((observationY - landmarkY)*(observationY - landmarkY)/(2*std_landmark[1]*std_landmark[1])) ) );
	      if (weight == 0) {
	        particles[i].weight *= EPS;
	      } else {
	        particles[i].weight *= weight;
	      }
	    }

	   // std::cout<<"weightUPdated "<<i<<" ok "<<particles[i].weight<<std::endl;
	}



	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  
    //std::cout<<"resample start"<<std::endl;
	std::vector<double> weights;

	for(int i=0;i<num_particles;i++)
	{	
		weights.push_back(particles[i].weight);
	}
    //std::cout<<"weights size"<<weights.size()<<std::endl;
	discrete_distribution<int> reSamp(weights.begin(),weights.end());  //if this function is used, no need to implement the wheel?
    //std::cout<<"resamp dist ok"<<std::endl;
  
	vector<Particle> resampledParticles;

	for(int i=0;i<num_particles;i++)
	{
		int resamp_id=reSamp(gen);
		resampledParticles.push_back(particles[resamp_id]);
	}

	particles=resampledParticles;
    //std::cout<<"resamp ok"<<std::endl;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
