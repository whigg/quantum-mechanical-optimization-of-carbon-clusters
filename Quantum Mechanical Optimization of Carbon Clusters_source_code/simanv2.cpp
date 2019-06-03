#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "libnum/libnum.h"
#include "libnum/displot.h"
#include "atomosv3.h"
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>

#define PI 3.14159265
#define Hartree 27.2113845;

using namespace std;



int main(int argc, char *argv[])
{
  int M,iter,i,t,rej,N,iatom,count,NT,frame;
  double Rmax,Rmax0;
  double oldx,oldy,oldz,xp,yp,zp,xi,yi,zi,xi2,yi2,zi2,A,B,C,q,p;
  double Energia,Energia2,En,Enq,Enmq,Enqm,Cv,alpha;
  double de,r,g,w,Temp,Tmax,dist,e;

  N=0;
  frame=0;
//  printf("here0 %s %s %s\n",argv[0],argv[1],argv[2]);

  Rmax0=atof(argv[1]);  // Diâmetro máximo da distribuição inicial
  iter=atoi(argv[2]);   //nº de iterações do passeio de Metropolis Monte-Carlo
  alpha=atof(argv[3]);  //parâmetro alpha em T=alpha T1 (ele calcula automaticamente o nº de passos nt do simmulated annealing,
                        //para uma temperatura final de 0.00001)
  Tmax=atof(argv[4]);   //temperatura inicial
  M=atoi(argv[5]);      // nº de átomos no cluster

  Rmax=Rmax0;

//Tmax*alpha^e = Tmin, então "e" é o nº de iterações que queremos
  e=log(0.00001/Tmax)/log(alpha);

  int nt=int(e);



  TBCluster Cn(M);
  
  LIBNUM::RngInit();

//ponho o 1º átomo na origem e os outros estão numa rede cúbica à volta deste (é só uma forma de os arranjar inicialmente,
//não é muito importante)
  Cn.atom[0].r.x=0.0;
  Cn.atom[0].r.y=0.0;
  Cn.atom[0].r.z=0.0;

  for(int k=1;k<M;k++)
  {	for(int a=-1; a<=1;a++)
	{  for(int b=-1; b<=1;b++)
	    {  for(int c=-1; c<=1;c++)
		{
		      A=a*1.0*3.5;
		      B=b*1.0*3.5;
		      C=c*1.0*3.5;
		      p=sqrt(A*A+B*B+C*C);
		      count=0;

		  if(p<Rmax0*.5)
		  {  
		      for(int i=0;i<k;i++)
		      {
			          q=(A-Cn.atom[i].r.x)*(A-Cn.atom[i].r.x)+(B-Cn.atom[i].r.y)*(B-Cn.atom[i].r.y)+(C-Cn.atom[i].r.z)*(C-Cn.atom[i].r.z);

			  if(q>2.0)
			  {
			      count+=1;
			  }
			  else
			  {
			  }
			  
		      }	
	      
		      if(count==k)
		      {
			  Cn.atom[k].r.x=A;
			  Cn.atom[k].r.y=B;
			  Cn.atom[k].r.z=C;
    
		      }
		      else
		      {
		      }

		    }

		}
	    }
	 }

  }

  //FILE *fout=fopen("vulcao2.bs","w");

	Cn.PrintDistances();
//frame 0

  NT=0;

    for(int j=0;j<=2;j++)
    {
     Temp=Tmax;
      
//aqui começa o simmulated annealing
      for(t=1;t<=nt;t++)
      {   NT+=1;
          N=0;
          rej=0;
    
          //En=0.0;
    
          //passeio de Metropolis Monte-Carlo, em que em cada iteração escolho um dos átomos e dou novas coordenadas aleatoriamente
          //numa esfera de raio Rmax, e se a energia do novo agregado for mais pequena escolho as novas coordenadas. Se não, há
          //aquela probabilidade de escolher ou não escolher, dada pela função de Boltzman. Após 100 iterações reajusto o raio Rmax.
          for(i=0;i<iter;i++)
          {   
            	N+=1;	
            	Energia=Cn.getEnergyAtractive()+Cn.getEnergyRepulsive();
            	//En+=Energia;
            
            	iatom=LIBNUM::RngInt(M-1)+1;
            
            	oldx=Cn.atom[iatom].r.x;
            	oldy=Cn.atom[iatom].r.y;
            	oldz=Cn.atom[iatom].r.z;
            
            	Cn.AddPerturbation(Rmax,xp,yp,zp);
            
            	Cn.atom[iatom].r.x+=xp;
            	Cn.atom[iatom].r.y+=yp;
            	Cn.atom[iatom].r.z+=zp;
              
            	Energia2=Cn.getEnergyAtractive()+Cn.getEnergyRepulsive();
            
            	de=Energia2-Energia;
            
            	if (de<0.0)
            	{
            	} 
               
            	else
            	{     g=LIBNUM::RngDouble();
                     
            	      w=exp(-de/(8.617343*0.00001*Temp));
            
            	      if(g>w)
            	      {   
            		  Cn.atom[iatom].r.x=oldx;
            		  Cn.atom[iatom].r.y=oldy;
            		  Cn.atom[iatom].r.z=oldz;
            		  rej+=1;
            	      }
            
            	      else
            	      {
            	      }
            	 }
            
            	  if(N==100)
            	  {	if (rej>50)
            		{	Rmax=Rmax*0.95;
            			rej=0;
            			N=0;
            		}
            		else 
            		{       Rmax=Rmax*1.05;
            			rej=0;
            			N=0;
            		}
            	  }
            	  else
            	  {   
            	  }
    
            }
    
    //de 10 em 10 passos do simmulated annealing faz um frame do vídeo
        
          Temp=alpha*Temp;

          }
  
    }

	double Ebond = Cn.getEnergyBondPerAtom()*Hartree;
	double E = Cn.getEnergy()*Hartree;
	double Eatr = Cn.getEnergyAtractive()*Hartree;
	double Erep = Cn.getEnergyRepulsive()*Hartree;

	cout << "E atractiva " << Eatr << "\n";
	cout << "E repulsiva "<< Erep << "\n";
	cout << "E "<< E << "\n";
	cout << "E ligaÃ§Ã£o por atomo "<< Ebond+1.13 << "\n";
	Cn.PrintDistances();

    return 0;

}



