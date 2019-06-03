#include <cmath>
#include <string>
#include <iostream>
#include <libnum/libnum.h>
#include <sys/time.h>
using namespace std;
#include "objectsc60.h"

double createPerturbation(double Rmax, double &addx,double &addy, double &addz){
	double t,theta,f,fi,r;
	t=LIBNUM::RngDouble();
	f=LIBNUM::RngDouble();
	r=(LIBNUM::RngDouble())*Rmax;
	theta=t*2.0*2.0*asin(1.0);
	fi=f*2.0*asin(1.0);

	addx=r*cos(theta)*sin(fi);
	addy=r*sin(theta)*sin(fi);
	addz=r*cos(fi);

}


int main(){

	LIBNUM::RngInit();


	string movieFileName = "c60";
	string lastClusterFile = "last_c60";
	const double hartree = 27.2113845;
	int n = 60; // atomos
	int k = 1000*n; // iterações
	double T = 5000000.0;
	double Tmin = 1.0E-4;
	double coolingRate = 0.05; // = (1- alpha)
	double Kb = 8.617343E-5;
	double Rmax = 1.0;
	double RrateOfChange = 0.05; // variação de Rmax
	int frameRate = k/5;
	double E1,E2,iatom,deltaE,oldx,oldy,oldz,cnt,xp,yp,zp,rej;

	struct timeval start, end;
	long mtime, seconds, useconds;

	Cluster cn(n,6.6,6.8); // c60 r = 6.7 , c20 r = 3.85 raios das cages
	cn.initializeSpherical(6.7);

	cn.createMovie(movieFileName);

	E1 = cn.getEnergy();
	while(T>=Tmin){
		cnt = 0;
		for (int i = 0; i < k ; i++){
			cnt++;
			iatom=LIBNUM::RngInt(n);
			oldx=cn.getAtom(iatom).getX();
			oldy=cn.getAtom(iatom).getY();
			oldz=cn.getAtom(iatom).getZ();

			createPerturbation(Rmax,xp,yp,zp);
			cn.getAtom(iatom).setX(oldx+xp);
			cn.getAtom(iatom).setY(oldy+yp);
			cn.getAtom(iatom).setZ(oldz+zp);

			E2 = cn.getEnergy();

			deltaE = E2-E1;

			if (deltaE < 0.0){E1 = E2;} 
			else{
				if( LIBNUM::RngDouble() > exp(-deltaE/(Kb*T)) ){
					cn.getAtom(iatom).setX(oldx);
					cn.getAtom(iatom).setY(oldy);
					cn.getAtom(iatom).setZ(oldz);
					rej++;
				}
				else {E1 = E2;} 
			}

			if(cnt == 100){
				if (rej > 50){ Rmax = Rmax * (1.0 - RrateOfChange);}
				else { Rmax = Rmax * (1.0 + RrateOfChange);}
				rej = 0;
				cnt = 0;
			}
			if((i+1)%frameRate == 0){
				//cn.centerAtoms();
				cn.appendFrame(movieFileName);
				cn.storeInFile(lastClusterFile);
				gettimeofday(&end, NULL);
				seconds  = end.tv_sec  - start.tv_sec;
				useconds = end.tv_usec - start.tv_usec;
				mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;

cout << T << " " << E1*hartree << " " << (0.0415+(E1-n*(-(2.0*0.50097+2.0*0.19930)))/n)*hartree << " " << Rmax << " " << mtime << "\n";
				gettimeofday(&start, NULL);
			}
		}

		T = (1.0 - coolingRate)*T;
	}

	return 0;



}

