#include<iostream>
#include<cmath>
#include "vector.h"

const double G = 1.0;
const int N = 2;

const double xi = 0.1786178958448091;
const double lambda = -0.2123418310626054;
const double chi = -0.06626458266981849;
const double Um2lambdau2 = (1 - 2 * lambda) / 2;
const double Um2chiplusxi = 1 - 2 * (chi + xi);

class body;
class collider;

//------------------------------------------//

class body {
private:
	vector3D r, V, F;
	double m, R;
public:
	void start(double x0, double y0, double z0, 
		double Vx0, double Vy0, double Vz0, double m, double R0);
	void SumForce(vector3D dF);
	void EraseForce(void) { F.load(0, 0, 0); };
	void Mover(double dt, double coeff);
	void Movev(double dt, double coeff);
	double Getx(void) { return r.x(); }; //inline
	double Gety(void) { return r.y(); }; //inline
	friend class collider;
};

class collider {
private:

public:
	void CalculateEveryForce(body * bodies);
	void CalculateForceBetween(body & planet1, body & planet2);
};
//-----------------------------------------------//
void body::start(double x0, double y0, double z0, double Vx0, 
	double Vy0, double Vz0, double m0, double R0) {
	r.load(x0, y0, z0);
	V.load(Vx0, Vy0, Vz0);
	m = m0;
	R = R0;
}

void body::SumForce(vector3D dF) {
	F += dF;
}

void body::Mover(double dt, double coeff) {
	r += V * coeff * dt ;
}

void body::Movev(double dt, double coeff) {
	V += (F / m) * coeff * dt;
}

void collider::CalculateEveryForce(body * planet){
	int i, j;
	for (i = 0; i < N; i++) {
		planet[i].EraseForce();
	}
	for (i = 0; i < N; i++) {
		for (j = 0; j < i; j++) {
			CalculateForceBetween(planet[i], planet[j]);
		}
	}
}

void collider::CalculateForceBetween(body& planet1, body& planet2) {
	double m1 = planet1.m;
	double m2 = planet2.m;
	vector3D r21 = planet2.r - planet1.r;
	double r2 = r21.norm2();
	double aux = G * m2 * m1 * std::pow(r2, -1.5);
	vector3D F1 = aux * r21;
	planet1.SumForce(F1);
	planet2.SumForce(F1 * (-1));
}

//--------------------------------------------//
int main() {
	//Body 0 = Sun, Body 1 = Jupiter
	double dt = 0.1;                               //Time step  
	double r = 1000.0;	                           //Distance from Sun to Jupiter	
	double m0 = 1047.0, m1 = 1.0;                  //Masses
	double M = m0 + m1;							   //Total mass													
	double mu = (m0 * m1) / M;					   //Reduced mass
	double x0 = -(m1 * r) / M, x1 = (m0 * r) / M;  //Initial positons to place C.M. in (0, 0)
	double omega = sqrt(G * M / (pow(r, 3)));	   //Angular velocity
	double T = 2 * M_PI / omega;				   //Period
	double Ttotal = 20 * T;						   //Time frame			
	double v0 = omega * x0, v1 = omega * x1;	   //Initial velocities
	
	collider newton;
	body bodies[N];

	//x0, y0, z0, Vx0, Vy0, Vz0, m0, R0
	bodies[0].start(x0, 0, 0, 0, v0, 0, m0, 1.0);
	bodies[1].start(x1, 0, 0, 0, v1, 0, m1, 0.1);

	for (double t = 0; t < Ttotal; t += dt) {
		std::cout << bodies[0].Getx() << " " << bodies[0].Gety() << " " << bodies[1].Getx() << " " << bodies[1].Gety() << std::endl;

		for (int i = 0; i < N; i++) {
			bodies[i].Mover(dt, xi);
		}
		newton.CalculateEveryForce(bodies);

		for (int i = 0; i < N; i++) {
			bodies[i].Movev(dt, Um2lambdau2);
		}

		for (int i = 0; i < N; i++) {
			bodies[i].Mover(dt, chi);
		}

		newton.CalculateEveryForce(bodies);

		for (int i = 0; i < N; i++) {
			bodies[i].Movev(dt, lambda);
		}

		for (int i = 0; i < N; i++) {
			bodies[i].Mover(dt, Um2chiplusxi);
		}

		newton.CalculateEveryForce(bodies);

		for (int i = 0; i < N; i++) {
			bodies[i].Movev(dt, lambda);
		}

		for (int i = 0; i < N; i++) {
			bodies[i].Mover(dt, chi);
		}

		newton.CalculateEveryForce(bodies);

		for (int i = 0; i < N; i++) { 
			bodies[i].Movev(dt, Um2lambdau2);
		}

		for (int i = 0; i < N; i++) {
			bodies[i].Mover(dt, xi);
		}
	}

	return 0;
}