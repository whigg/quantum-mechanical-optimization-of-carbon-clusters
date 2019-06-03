#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>


#define PI 3.14159265

class Vector3D{
	double x,y,z;
public:
	Vector3D(double x, double y, double z){
		this-> x = x;
		this-> y = y;
		this-> z = z;
	}
	Vector3D(){
		this-> x = 0.0;
		this-> y = 0.0;
		this-> z = 0.0;
	}

	double getX(){return x;}

	double getY(){return y;}

	double getZ(){return z;}

	double NormSq(){return(x*x+y*y+z*z);}

	double Norm() {return(sqrt(NormSq()));}

	double dotProduct(Vector3D v){return v.getX()*x+v.getY()*y+v.getZ()*z;}

	void setX(double x){this->x = x;}

	void setY(double y){this->y = y;}

	void setZ(double z){this->z = z;}

	void setCoordinates(double x, double y, double z){
		setX(x);
		setY(y);
		setZ(z);
	}
};


class Atom: public Vector3D{
	string name;
public:
	Atom(string name, double x, double y, double z):Vector3D(x,y,z){this-> name = name;}

	Atom():Vector3D(){this-> name = "NoName";}

	string getName(){return name;}

	double getDistanceTo(Atom at){return  sqrt((at.getX()-getX())*(at.getX()-getX())+
			(at.getY()-getY())*(at.getY()-getY())+
			(at.getZ()-getZ())*(at.getZ()-getZ()));}

	void PrintCoordinates(){printf("%s %14lf %14lf %14lf\n",name.c_str(),getX(),getY(),getZ());}

	void ReadFromFile(FILE *fin){
		char el[10];
		double x,y,z;
		fscanf(fin,"%s %lf %lf %lf\n", &el,&x, &y, &z);
		name = el;
		setCoordinates(x,y,z);
	}
};

class Cluster
{
	int n;
	Atom *atom;
	double minRadius,maxRadius;
	struct timeval start, end;
	long mtime, seconds, useconds;
public:

	Cluster(int n, double minRadius, double maxRadius){
		this->n = n;
		atom = new Atom[n];
		this-> minRadius = minRadius;
		this-> maxRadius = maxRadius;
	}

	Cluster(int n, const char * filename, double minRadius, double maxRadius){
		this->n = n;
		atom = new Atom[n];
		this-> minRadius = minRadius;
		this-> maxRadius = maxRadius;
		FILE * pFile = fopen (filename,"r");
		for(int i=0;i<n;i++){atom[i].ReadFromFile(pFile);}
		fclose (pFile);
	}

	~Cluster(){ delete[] atom; }

	Atom & getAtom(int i){return atom[i];}

	void PrintCoordinates(){
		for(int i=0;i<n;i++){  atom[i].PrintCoordinates(); }
	}

	void ReadFromFile(string FileName){
		FILE *fin = fopen(FileName.c_str(),"r");
		if(fin==NULL){
			printf("FATAL Cluster::ReadFromFile, file %s not found\n",FileName.c_str());
			exit(1);
		}

		for(int i=0;i<n;i++){
			atom[i].ReadFromFile(fin); // este ReadFromFile Ã© um mÃ©todo
			// da class Atom !!!
		}

		fclose(fin);
	}

	double getDistance(int i, int j) // mÃ©todo
	{
		// Calculate the distance between atom[i] and atom[j]
		if(i==j)
		{
			printf("Warning Cluster::getDistance, i and j are the same!\n");
			return(0.0);
		}
		else{ return atom[i].getDistanceTo(atom[j]); }
	}

	void PrintDistances(){
		for(int i=0;i<n;i++){
			for(int j=i+1;j<n;j++){
				printf("|r%2.2d-r%2.2d| = %14lf\n",j,i,getDistance(i,j));
			}
		}
	}

	double Chebyshev (double a,double b, double r,double c[10]){

		double y = (r-(b+a)/2.0)/((b-a)/2.0);

		return c[0]/2.0 + c[1]*y + c[2]*(2.0*y*y-1.0) + c[3]*y*(4.0*y*y-3.0)
		+ c[4]*(8.0*y*y*(y*y-1.0)+1.0) + c[5]*(4.0*y*y*(4.0*y*y*y-5.0*y)+5.0*y)
		+ c[6]*(32.0*y*y*y*y*y*y-48.0*y*y*y*y+18.0*y*y-1.0)
		+ c[7]*(64.0*y*y*y*y*y*y*y-112.0*y*y*y*y*y+56.0*y*y*y-7.0*y)
		+ c[8]*(128.0*y*y*y*y*y*y*y*y-256.0*y*y*y*y*y*y+160.0*y*y*y*y-32.0*y*y+1.0)
		+ c[9]*(256.0*y*y*y*y*y*y*y*y*y-576.0*y*y*y*y*y*y*y+432.0*y*y*y*y*y-120.0*y*y*y+9.0*y);

	}

	void storeInFile(string filename){
		FILE *fout= fopen((filename+".dat").c_str(),"w");
		for ( int i = 0; i < n ; i++){
			fprintf(fout,"%s %14lf %14lf %14lf\n",atom[i].getName().c_str(),atom[i].getX(),atom[i].getY(),atom[i].getZ());
		}
		fclose(fout);
	}

	void createMovie(string filename){
		FILE *fout= fopen((filename+".bs").c_str(),"w");
		FILE *fout1= fopen((filename+".mv").c_str(),"w");

		fprintf(fout,"\n");
		fprintf(fout,"inc 1.0");
		fprintf(fout,"\n\n");
		fprintf(fout,"spec C 0.6 Orange");
		fprintf(fout,"\n\n");
		fprintf(fout,"bonds C C 0.0 3.0 0.1 Blue");
		fprintf(fout,"\n\n");

		fprintf(fout1,"frame \n");

		for ( int i = 0; i < n ; i++){
			fprintf(fout,"atom C %14lf %14lf %14lf\n",atom[i].getX(),atom[i].getY(),atom[i].getZ());
			fprintf(fout1,"%14lf %14lf %14lf\n",atom[i].getX(),atom[i].getY(),atom[i].getZ());
		}

		fclose(fout);
		fclose(fout1);
	}

	void appendFrame(string filename){
		FILE *fout= fopen((filename+".mv").c_str(),"a");
		fprintf(fout,"frame\n");
		for ( int i = 0; i < n ; i++){
			fprintf(fout,"%14lf %14lf %14lf\n",atom[i].getX(),atom[i].getY(),atom[i].getZ());
		}

		fclose(fout);
	}

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
					double r = getDistance(a,b);

					if (r>b1 || r<a1){
						for (int i = 0 ; i < 4 ; i++){
							for (int j = 0 ; j < 4 ; j++){gsl_matrix_set(mat,a*4+i,b*4+j,0.0);}
						}
					}
					else{
						double c[3] = {(atom[b].getX()-atom[a].getX())/r,(atom[b].getY()-atom[a].getY())/r,(atom[b].getZ()-atom[a].getZ())/r};

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
					double r = getDistance(a,b);

					if (r>b1 || r<a1){
						for (int i = 0 ; i < 4 ; i++){
							for (int j = 0 ; j < 4 ; j++){gsl_matrix_set(mat,a*4+i,b*4+j,0.0);}
						}
					}
					else{
						double c[3] = {(atom[b].getX()-atom[a].getX())/r,(atom[b].getY()-atom[a].getY())/r,(atom[b].getZ()-atom[a].getZ())/r};

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
		gsl_vector_complex * alpha = gsl_vector_complex_calloc(n*4);
		gsl_vector * beta = gsl_vector_calloc(n*4);
		gsl_vector * energies = gsl_vector_calloc(n*4);
		gsl_matrix * HH = gsl_matrix_calloc (n*4, n*4);
		gsl_matrix * SS = gsl_matrix_calloc (n*4, n*4);
		gsl_eigen_gen_workspace * workH = gsl_eigen_gen_alloc (n*4);

		H(HH);
		S(SS);
		gsl_eigen_gen (HH,SS, alpha, beta, workH);

		for(int i=0;i<n*4;i++){gsl_vector_set(energies,i,GSL_REAL(gsl_vector_complex_get(alpha,i))/gsl_vector_get(beta,i));}

		gsl_sort_vector(energies);

		for (int i = 0 ; i < n*2 ; i++){
			energy += 2.0*gsl_vector_get(energies,i);
		}
		gsl_matrix_free(HH);
		gsl_matrix_free(SS);
		gsl_vector_complex_free(alpha);
		gsl_vector_free(beta);
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
			double radius = atom[a].Norm();

			if ( radius > maxRadius || radius < minRadius ){energy += 1.0E30;}

			for ( int b = 0; b < a ; b++){
				double r = getDistance(a,b);

				if (r <= b1) {
					if (r < a1) {energy += 1.0E20;}
					else{energy += Chebyshev(a1,b1,r,VCCrep);}
				}
			}
		}
		return energy;
	}

	double centerAtoms(){
		for ( int a = 1; a < n ; a++){
			atom[a].setX(atom[a].getX()-atom[0].getX());
			atom[a].setY(atom[a].getY()-atom[0].getY());
			atom[a].setZ(atom[a].getZ()-atom[0].getZ());
		}
		atom[0].setX(0.0);
		atom[0].setY(0.0);
		atom[0].setZ(0.0);
	}

	double getEnergyBondPerAtom(){	

		double Etot, Eint;
		double Espin = 0.0415; //hartree

		Etot=getEnergyAtractive()+getEnergyRepulsive();
		Eint=-(2.0*0.50097+2.0*0.19930);
		return Espin+(Etot-n*Eint)/n;


	}

	double getEnergy(){
		double Erep = getEnergyRepulsive();
		if (Erep >= 1.0E10){return Erep;}
		else{return Erep+getEnergyAtractive();}

		//return getEnergyAtractive()+getEnergyRepulsive();
	}

	void initializeSpherical(double raio){
		double x,y,z;
		for (int i = 0 ; i < n ; i++){
			x = 2.0*(LIBNUM::RngInt(2)-0.5)*raio*LIBNUM::RngDouble();
			y = 2.0*(LIBNUM::RngInt(2)-0.5)*raio*LIBNUM::RngDouble();
			z = 2.0*(LIBNUM::RngInt(2)-0.5)*sqrt(raio*raio-y*y-x*x);


			while( (x*x+y*y) > raio*raio){
				x = 2.0*(LIBNUM::RngInt(2)-0.5)*raio*LIBNUM::RngDouble();
				y = 2.0*(LIBNUM::RngInt(2)-0.5)*raio*LIBNUM::RngDouble();
				z = 2.0*(LIBNUM::RngInt(2)-0.5)*sqrt(raio*raio-y*y-x*x);
			}
			atom[i].setCoordinates(x, y, z);
		}
	}

	void initializeBox(double side){
		double x,y,z;
		for (int i = 0 ; i < n ; i++){
			x = side*LIBNUM::RngDouble();
			y = side*LIBNUM::RngDouble();
			z = side*LIBNUM::RngDouble();

			atom[i].setCoordinates(x, y, z);
		}
	}

	void initializeChain(double leght){
		double x,y,z;
		for (int i = 0 ; i < n ; i++){

			x = leght*i/n;;
			y = 0.0;
			z = 0.0;

			atom[i].setCoordinates(x, y, z);
		}
	}

	void initializeRing(double raio){
		double x,y,z,tetta;
		for (int i = 0 ; i < n ; i++){
			tetta = i*2.0*PI/n;

			x = raio*sin(tetta);
			y = raio*cos(tetta);
			z = 0.0;

			atom[i].setCoordinates(x, y, z);
		}
	}

	void geometry(int neighbors){

		double atomindex[n][neighbors];

		for (int i = 0; i < n ; i++){
			gsl_vector * vectOrdered = gsl_vector_calloc(n); 
			gsl_vector * vect= gsl_vector_calloc(n);

			for (int j = 0; j < n ; j++){
				double data = round(1000.0*getDistance(i,j))/1000.0;
				gsl_vector_set(vectOrdered,j,data); 
				gsl_vector_set(vect,j,data);
			}

			gsl_sort_vector(vectOrdered);

			int idx = 0;
			for (int j = 1; j <= neighbors ; j++){

				while( gsl_vector_get(vectOrdered,j) != gsl_vector_get(vect,idx)) { idx++; if(idx == n){ idx = 0;}}

				atomindex[i][j-1] = idx;
				cout << "Atom " << i+1 << " Vizinho " << j << " raio " << gsl_vector_get(vect,idx) << " bohrs\n";

				idx++;
				if(idx == n){ idx = 0;}
			}

			gsl_vector_free(vect);
			gsl_vector_free(vectOrdered); 
		}

		for (int i = 0; i < n ; i++){

			for (int l = 0; l < neighbors ; l++){
				for (int m = 0; m < l ; m++){


					Vector3D vect1(getAtom(i).getX()-getAtom(atomindex[i][l]).getX(),getAtom(i).getY()-getAtom(atomindex[i][l]).getY(),getAtom(i).getZ()-getAtom(atomindex[i][l]).getZ());

					Vector3D vect2(getAtom(i).getX()-getAtom(atomindex[i][m]).getX(),getAtom(i).getY()-getAtom(atomindex[i][m]).getY(),getAtom(i).getZ()-getAtom(atomindex[i][m]).getZ());

					cout << "Atom " << i+1 << " AnguloVizinhos(" << l+1 << "," << m+1 << ") " << 360.0*acos(vect1.dotProduct(vect2)/(vect1.Norm()*vect2.Norm()))/(2.0*PI) << "\n";
				}
			}
		}
	}

	void initializeDodecahedron(double r){


		double gamma = (1.0+sqrt(5.0))/2.0;
		double gamma_inv = 2.0/(1.0+sqrt(5.0));

		r = r*gamma/2.0;

		atom[0].setCoordinates(r,r,r);
		atom[1].setCoordinates(r,r,-r);
		atom[2].setCoordinates(r,-r,r);
		atom[3].setCoordinates(r,-r,-r);
		atom[4].setCoordinates(-r,r,r);
		atom[5].setCoordinates(-r,r,-r);
		atom[6].setCoordinates(-r,-r,r);
		atom[7].setCoordinates(-r,-r,-r);

		atom[8].setCoordinates(0.0,r*gamma_inv,r*gamma);
		atom[9].setCoordinates(0.0,r*gamma_inv,-r*gamma);
		atom[10].setCoordinates(0.0,-r*gamma_inv,r*gamma);
		atom[11].setCoordinates(0.0,-r*gamma_inv,-r*gamma);

		atom[12].setCoordinates(r*gamma_inv,r*gamma,0.0);
		atom[13].setCoordinates(r*gamma_inv,-r*gamma,0.0);
		atom[14].setCoordinates(-r*gamma_inv,r*gamma,0.0);
		atom[15].setCoordinates(-r*gamma_inv,-r*gamma,0.0);

		atom[16].setCoordinates(r*gamma,0.0,r*gamma_inv);
		atom[17].setCoordinates(r*gamma,0.0,-r*gamma_inv);
		atom[18].setCoordinates(-r*gamma,0.0,r*gamma_inv);
		atom[19].setCoordinates(-r*gamma,0.0,-r*gamma_inv);

	}

	void intializeTruncatedIcosahedron(double r){

		double gamma = (1.0+sqrt(5.0))/2.0;
		r = r/2.0;

		atom[0].setCoordinates(0.0,r,3.0*r*gamma);
		atom[1].setCoordinates(0.0,r,-3.0*r*gamma);
		atom[2].setCoordinates(0.0,-r,3.0*r*gamma);
		atom[3].setCoordinates(0.0,-r,-3.0*r*gamma);

		atom[4].setCoordinates(r,3.0*r*gamma,0.0);
		atom[5].setCoordinates(r,-3.0*r*gamma,0.0);
		atom[6].setCoordinates(-r,3.0*r*gamma,0.0);
		atom[7].setCoordinates(-r,-3.0*r*gamma,0.0);

		atom[8].setCoordinates(3.0*r*gamma,0.0,r);
		atom[9].setCoordinates(3.0*r*gamma,0.0,-r);
		atom[10].setCoordinates(-3.0*r*gamma,0.0,r);
		atom[11].setCoordinates(-3.0*r*gamma,0.0,-r);

		atom[12].setCoordinates(r*2.0,r*(1.0+2.0*gamma),r*gamma);
		atom[13].setCoordinates(r*2.0,r*(1.0+2.0*gamma),-r*gamma);
		atom[14].setCoordinates(r*2.0,-r*(1.0+2.0*gamma),r*gamma);
		atom[15].setCoordinates(r*2.0,-r*(1.0+2.0*gamma),-r*gamma);
		atom[16].setCoordinates(-r*2.0,r*(1.0+2.0*gamma),r*gamma);
		atom[17].setCoordinates(-r*2.0,r*(1.0+2.0*gamma),-r*gamma);
		atom[18].setCoordinates(-r*2.0,-r*(1.0+2.0*gamma),r*gamma);
		atom[19].setCoordinates(-r*2.0,-r*(1.0+2.0*gamma),-r*gamma);

		atom[20].setCoordinates(r*(1.0+2.0*gamma),r*gamma,r*2.0);
		atom[21].setCoordinates(r*(1.0+2.0*gamma),r*gamma,-r*2.0);
		atom[22].setCoordinates(r*(1.0+2.0*gamma),-r*gamma,r*2.0);
		atom[23].setCoordinates(r*(1.0+2.0*gamma),-r*gamma,-r*2.0);
		atom[24].setCoordinates(-r*(1.0+2.0*gamma),r*gamma,r*2.0);
		atom[25].setCoordinates(-r*(1.0+2.0*gamma),r*gamma,-r*2.0);
		atom[26].setCoordinates(-r*(1.0+2.0*gamma),-r*gamma,r*2.0);
		atom[27].setCoordinates(-r*(1.0+2.0*gamma),-r*gamma,-r*2.0);

		atom[28].setCoordinates(r*gamma,r*2.0,r*(1.0+2.0*gamma));
		atom[29].setCoordinates(r*gamma,r*2.0,-r*(1.0+2.0*gamma));
		atom[30].setCoordinates(r*gamma,-r*2.0,r*(1.0+2.0*gamma));
		atom[31].setCoordinates(r*gamma,-r*2.0,-r*(1.0+2.0*gamma));
		atom[32].setCoordinates(-r*gamma,r*2.0,r*(1.0+2.0*gamma));
		atom[33].setCoordinates(-r*gamma,r*2.0,-r*(1.0+2.0*gamma));
		atom[34].setCoordinates(-r*gamma,-r*2.0,r*(1.0+2.0*gamma));
		atom[35].setCoordinates(-r*gamma,-r*2.0,-r*(1.0+2.0*gamma));

		atom[36].setCoordinates(r*1.0,r*(2.0+gamma),r*2.0*gamma);
		atom[37].setCoordinates(r*1.0,r*(2.0+gamma),-r*2.0*gamma);
		atom[38].setCoordinates(r*1.0,-r*(2.0+gamma),r*2.0*gamma);
		atom[39].setCoordinates(r*1.0,-r*(2.0+gamma),-r*2.0*gamma);
		atom[40].setCoordinates(-r*1.0,r*(2.0+gamma),r*2.0*gamma);
		atom[41].setCoordinates(-r*1.0,r*(2.0+gamma),-r*2.0*gamma);
		atom[42].setCoordinates(-r*1.0,-r*(2.0+gamma),r*2.0*gamma);
		atom[43].setCoordinates(-r*1.0,-r*(2.0+gamma),-r*2.0*gamma);

		atom[44].setCoordinates(r*(2.0+gamma),r*2.0*gamma,r*1.0);
		atom[45].setCoordinates(r*(2.0+gamma),r*2.0*gamma,-r*1.0);
		atom[46].setCoordinates(r*(2.0+gamma),-r*2.0*gamma,r*1.0);
		atom[47].setCoordinates(r*(2.0+gamma),-r*2.0*gamma,-r*1.0);
		atom[48].setCoordinates(-r*(2.0+gamma),r*2.0*gamma,r*1.0);
		atom[49].setCoordinates(-r*(2.0+gamma),r*2.0*gamma,-r*1.0);
		atom[50].setCoordinates(-r*(2.0+gamma),-r*2.0*gamma,r*1.0);
		atom[51].setCoordinates(-r*(2.0+gamma),-r*2.0*gamma,-r*1.0);

		atom[52].setCoordinates(r*2.0*gamma,r*1.0,r*(2.0+gamma));
		atom[53].setCoordinates(r*2.0*gamma,r*1.0,-r*(2.0+gamma));
		atom[54].setCoordinates(r*2.0*gamma,-r*1.0,r*(2.0+gamma));
		atom[55].setCoordinates(r*2.0*gamma,-r*1.0,-r*(2.0+gamma));
		atom[56].setCoordinates(-r*2.0*gamma,r*1.0,r*(2.0+gamma));
		atom[57].setCoordinates(-r*2.0*gamma,r*1.0,-r*(2.0+gamma));
		atom[58].setCoordinates(-r*2.0*gamma,-r*1.0,r*(2.0+gamma));
		atom[59].setCoordinates(-r*2.0*gamma,-r*1.0,-r*(2.0+gamma));

	}

};

