#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#define PI 3.14159265

class Vector3D
{
      public:
	double x;
	double y;
	double z;
	double NormSq() // mÃ©todo
	{
		return(x*x+y*y+z*z);
	}

	double Norm() // mÃ©todo
	{
		return(sqrt(NormSq()));
	}
};


class Atom
{
      public:
	string name;
	Vector3D r;
	// constructores nÃ£o sÃ£o obrigatÃ³rios, e podem ser substituÃ­dos por mÃ©todos...
	void SetCoordinates(string name, double x, double y, double z)
	{
		this->name = name;
		r.x = x;
		r.y = y;
		r.z = z;
	}

	void PrintCoordinates() // mÃ©todo
	{
		printf("%s %14lf %14lf %14lf\n",name.c_str(),r.x,r.y,r.z);
	}

	void ReadFromFile(FILE *fin) // mÃ©todo
	{
		char el[3];
		fscanf(fin,"%s %lf %lf %lf\n", &el,&r.x, &r.y, &r.z);
		name = el;
	}
};


class Cluster
{
      public:
	Atom *atom;
	string name;
	int n;

	Cluster(int n) // constructor
	{
		this->n = n;
		atom = new Atom[n];
	}

	~Cluster() // destructor
	{
		delete[] atom;
	}

	void PrintCoordinates() // mÃ©todo
	{
		for(int i=0;i<n;i++)
		{
			atom[i].PrintCoordinates();
		}
	}

	void ReadFromFile(string FileName) // mÃ©todo
	{
		FILE *fin = fopen(FileName.c_str(),"r");
		if(fin==NULL)
		{
			printf("FATAL Cluster::ReadCoordinates, file %s not found\n",FileName.c_str());
			exit(1);
		}

		for(int i=0;i<n;i++)
		{
			atom[i].ReadFromFile(fin); // este ReadFromFile Ã© um mÃ©todo
			// da class Atom !!!
		}

		fclose(fin);
	}

	double GetDistanceSq(int i, int j) // mÃ©todo
	{
		Vector3D rij;
		// Calculate the distance between atom[i] and atom[j]
		if(i==j)
		{
			printf("Warning Cluster::GetDistance, i and j are the same!\n");
			return(0.0);
		}
		else
		{
			rij.x = atom[j].r.x - atom[i].r.x;
			rij.y = atom[j].r.y - atom[i].r.y;
			rij.z = atom[j].r.z - atom[i].r.z;
			return(rij.NormSq()); // usa um mÃ©todo da class Vector3D
		}
	}

	void PrintDistances()
	{
		for(int i=0;i<n;i++)
		{
			for(int j=i+1;j<n;j++)
			{
				printf("|r%2.2d-r%2.2d| = %14lf\n",j,i,0.52917720859*sqrt(GetDistanceSq(i,j)));
			}
		}
	}


	double AddPerturbation(double Rmax, double &addx,double &addy, double &addz)
	{     double t,theta,f,fi,r;
	      t=LIBNUM::RngDouble();
	      f=LIBNUM::RngDouble();
	      r=(LIBNUM::RngDouble())*Rmax;
	      theta=t*2.0*2.0*asin(1.0);
	      fi=f*2.0*asin(1.0);

	      addx=r*cos(theta)*sin(fi);
	      addy=r*sin(theta)*sin(fi);
	      addz=r*cos(fi);

	}

	double Chebyshev (double a,double b, double r,double c[10])
	{

		double y = (r-(b+a)/2.0)/((b-a)/2.0);

		return c[0]/2.0 + c[1]*y + c[2]*(2.0*y*y-1.0) + c[3]*y*(4.0*y*y-3.0)
		+ c[4]*(8.0*y*y*(y*y-1.0)+1.0) + c[5]*(4.0*y*y*(4.0*y*y*y-5.0*y)+5.0*y)
		+ c[6]*(32.0*y*y*y*y*y*y-48.0*y*y*y*y+18.0*y*y-1.0)
		+ c[7]*(64.0*y*y*y*y*y*y*y-112.0*y*y*y*y*y+56.0*y*y*y-7.0*y)
		+ c[8]*(128.0*y*y*y*y*y*y*y*y-256.0*y*y*y*y*y*y+160.0*y*y*y*y-32.0*y*y+1.0)
		+ c[9]*(256.0*y*y*y*y*y*y*y*y*y-576.0*y*y*y*y*y*y*y+432.0*y*y*y*y*y-120.0*y*y*y+9.0*y);

	}

};



class TBCluster: public Cluster // TBCluster â€œherdaâ€ de Cluster
{
public:

	TBCluster(int n):Cluster(n){} // o mesmo â€œnâ€ em LJCluster e Cluster

	void H(gsl_matrix * mat){
		double a1 = 1.0;
		double b1 = 7.0;

		double HCCssSig1_7[10] = {-0.4663805,0.3528951,-0.1402985,0.0050519,0.0269723,-0.015881,0.0036716,0.0010301,-0.0015546,0.0008601};
		double HCCspSig1_7[10] = {0.3395418,-0.2250358,0.0298224,0.0653476,-0.0605786,0.0298962,-0.0099609,0.0020609,0.0001264,-0.0003381};
		double HCCppSig1_7[10] = {0.2422701,-0.1315258,-0.0372696,0.0942352,-0.0673216,0.0316900,-0.0117293,0.0033519,-0.0004838,-0.0000906};
		double HCCppPi1_7[10] = {-0.3793837,0.320447,-0.1956799,0.0883986,-0.0300733,0.0074465,-0.0008563,-0.0004453,0.0003842,-0.0001855};


		for ( int a = 0; a < n ; a++){
			for ( int b = 0; b < n ; b++){
				if (a == b){
					for (int i = 0 ; i < 4 ; i++){
						for (int j = 0 ; j < 4 ; j++){
							if (i != j){gsl_matrix_set(mat,a*4+i,b*4+j,0.0);}
							else if( i == 0 ) {gsl_matrix_set(mat,a*4+i,b*4+j,-0.50097);}
							else {gsl_matrix_set(mat,a*4+i,b*4+j,-0.19930);}
						}
					}
				}
				else{
					double r = sqrt(GetDistanceSq(a,b));

					if (r>b1 || r<a1){
						for (int i = 0 ; i < 4 ; i++){
							for (int j = 0 ; j < 4 ; j++){gsl_matrix_set(mat,a*4+i,b*4+j,0.0);}
						}
					}
					else{
						double c[3] = {(atom[b].r.x-atom[a].r.x)/r,(atom[b].r.y-atom[a].r.y)/r,(atom[b].r.z-atom[a].r.z)/r};

						for (int i = 0 ; i < 4 ; i++){
							for (int j = 0 ; j < 4 ; j++){
								if (i == 0 && j == 0){ gsl_matrix_set(mat,a*4+i,b*4+j,Chebyshev(a1,b1,r,HCCssSig1_7));}
								else if ( i == 0 ) {gsl_matrix_set(mat,a*4+i,b*4+j,c[j-1]*Chebyshev(a1,b1,r,HCCspSig1_7));}
								else if ( j == 0 ) {gsl_matrix_set(mat,a*4+i,b*4+j,-c[i-1]*Chebyshev(a1,b1,r,HCCspSig1_7));}
								else if ( i == j ){ gsl_matrix_set(mat,a*4+i,b*4+j,
										c[i-1]*c[i-1]*Chebyshev(a1,b1,r,HCCppSig1_7)+(1.0-c[i-1]*c[i-1])*Chebyshev(a1,b1,r,HCCppPi1_7));}
								else {gsl_matrix_set(mat,a*4+i,b*4+j,c[i-1]*c[j-1]*(Chebyshev(a1,b1,r,HCCppSig1_7) - Chebyshev(a1,b1,r,HCCppPi1_7)));}

							}
						}
					}
				}
			}
		}
	}

	void S(gsl_matrix * mat){
		double a1 = 1.0;
		double b1 = 7.0;

		double SCCssSig1_7[10] = {0.4728644,-0.3661623,0.1594782,-0.0204934,-0.0170732,0.0096695,-0.0007135,-0.0013826,0.0007849,-0.0002005};
		double SCCspSig1_7[10] = {-0.3662838,0.2490285,-0.0431248,-0.0584391,0.0492775,-0.0150447,-0.0010758,0.0027734,-0.0011214,0.0002303};
		double SCCppPi1_7[10] = {0.3715732,-0.3070867,0.1707304,-0.0581555,0.0061645,0.0051460,-0.0032776,0.0009119,-0.0001265,-0.0000227};
		double SCCppSig1_7[10] = {-0.1359608,0.0226235,0.1406440,-0.1573794,0.0753818,-0.0108677,-0.0075444,0.0051533,-0.0013747,0.0000751};


		for ( int a = 0; a < n ; a++){
			for ( int b = 0; b < n ; b++){
				if (a == b){
					for (int i = 0 ; i < 4 ; i++){
						for (int j = 0 ; j < 4 ; j++){
							if ( i == j){gsl_matrix_set(mat,a*4+i,b*4+j,1.0);}
							else{gsl_matrix_set(mat,a*4+i,b*4+j,0.0);}
						}
					}
				}
				else{
					double r = sqrt(GetDistanceSq(a,b));

					if (r>b1 || r<a1){
						for (int i = 0 ; i < 4 ; i++){
							for (int j = 0 ; j < 4 ; j++){gsl_matrix_set(mat,a*4+i,b*4+j,0.0);}
						}
					}
					else{
						double c[3] = {(atom[b].r.x-atom[a].r.x)/r,(atom[b].r.y-atom[a].r.y)/r,(atom[b].r.z-atom[a].r.z)/r};

						for (int i = 0 ; i < 4 ; i++){
							for (int j = 0 ; j < 4 ; j++){
								if (i == 0 && j == 0){ gsl_matrix_set(mat,a*4+i,b*4+j,Chebyshev(a1,b1,r,SCCssSig1_7));}
								else if ( i == 0 ) {gsl_matrix_set(mat,a*4+i,b*4+j,c[j-1]*Chebyshev(a1,b1,r,SCCspSig1_7));}
								else if ( j == 0 ) {gsl_matrix_set(mat,a*4+i,b*4+j,-c[i-1]*Chebyshev(a1,b1,r,SCCspSig1_7));}
								else if ( i == j ){ gsl_matrix_set(mat,a*4+i,b*4+j,
										c[i-1]*c[i-1]*Chebyshev(a1,b1,r,SCCppSig1_7)+(1.0-c[i-1]*c[i-1])*Chebyshev(a1,b1,r,SCCppPi1_7));}
								else {gsl_matrix_set(mat,a*4+i,b*4+j, c[i-1]*c[j-1]*(Chebyshev(a1,b1,r,SCCppSig1_7) - Chebyshev(a1,b1,r,SCCppPi1_7)));}

							}
						}
					}
				}
			}
		}
	}

	void printmatrix(gsl_matrix * mat, int d1, int d2){
		for (int i = 0 ; i < d1 ; i++){
			for (int j = 0 ; j < d2 ; j++){
				printf("%14f ",	gsl_matrix_get (mat, i, j));	
	
			}
			cout << "\n";
		}
	}

	double getEnergyAtractive(){
		double energy = 0.0;
		gsl_vector_complex * alpha = gsl_vector_complex_calloc(n*4);   //lara
		gsl_vector * beta = gsl_vector_calloc(n*4);   //lara
		gsl_vector * energies = gsl_vector_calloc(n*4);
		gsl_matrix * HH = gsl_matrix_calloc (n*4, n*4);
		gsl_matrix * SS = gsl_matrix_calloc (n*4, n*4);
		gsl_eigen_gen_workspace * workH = gsl_eigen_gen_alloc (n*4);

		H(HH);
		S(SS);

		gsl_eigen_gen (HH,SS, alpha, beta, workH);

		for(int i=0;i<n*4;i++){gsl_vector_set(energies,i,GSL_REAL(gsl_vector_complex_get(alpha,i))/gsl_vector_get(beta,i));}//lara

		gsl_sort_vector(energies);

		for (int i = 0 ; i < n*2 ; i++){
			energy += 2.0*gsl_vector_get(energies,i);
		}
		gsl_matrix_free(HH);
		gsl_matrix_free(SS);
		gsl_vector_complex_free(alpha);   //lara
		gsl_vector_free(beta);   //lara
		gsl_vector_free(energies);
		gsl_eigen_gen_free (workH);

		return energy;
	}

	double getEnergyRepulsive(){
		double r;
		double a1 = 1.0;
		double b1 = 4.1;
		double energy = 0.0;
		double VCCrep[10] = {2.2681036,-1.9157174,1.1677745,-0.5171036,0.1529242,-0.0219294,-0.0000002,-0.0000001,-0.0000005,0.0000009};

		for ( int a = 0; a < n ; a++){
			for ( int b = 0; b < a ; b++){
				double r = sqrt(GetDistanceSq(a,b));

				if (r <= b1) {
					if (r <= a1) {energy += 1000.0;}
					else{energy += Chebyshev(a1,b1,r,VCCrep);}
				}
			}
		}

		return energy;
	}

	double getEnergyBondPerAtom()
	{	double Etot, Eint;

	Etot=getEnergyAtractive()+getEnergyRepulsive();
	Eint=-(2*0.50097+2*0.19930);
	return (Etot-n*Eint)/n;


	}

	double getEnergy()
	{	double Etot, Eint;

	return getEnergyAtractive()+getEnergyRepulsive();
	}

};
