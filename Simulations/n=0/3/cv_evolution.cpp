#include<vector>     /*Library for using vector insted of arrays*/
#include<iostream>   /*Basic functions*/   
#include <cmath>     /*Define some fundamental mathematical operations*/
#include<stdlib.h>   /*Useful for random numbers*/
#include <sstream>
#include <fstream>   /*To read and write into files*/
#include <limits>
#include <stdio.h>   /*Basic functions in particular "sprintf"*/
#include <iomanip>      // std::setprecision

using namespace std;

#define pi 3.141592653589793238462643383279502884197169399375

#define SMALL 1e-20

// function declaration
vector<double> quatproduct(vector<double> P, vector<double> Q);           
vector<double> q_to_ex(vector<double> Q);
vector<double> q_to_ey(vector<double> Q);
vector<double> q_to_ez(vector<double> Q);
double sign(double a);
double length(vector<double> a);
double dotProduct (vector<double> a, vector<double> b);
vector<double> crossProduct(vector<double> a, vector<double> b);
void fix_roundoff(double &a);
vector<vector <double> > matrixOrthogonal(vector<double> P, vector<double> Q, vector<double> R);
vector<vector <double> > matrixTranspose(vector<vector <double> > M);
vector<vector <double> > matrixProduct(vector<vector <double> > P, vector<vector <double> > Q);
double matrixTrace(vector<vector <double> > M);

int main(int argc, char *argv[])
{
    /**********************************************/
    /*Pass the following arguments to the programm*/
    /**********************************************/
    /*To get help screen of the program: ./trajectories help*/
    if(argc<10)
    {
 
        cout << "Eleven arguments are required by the programm:"                      << endl;
        cout << "1.- Directory where the data is"                          << endl;
        cout << "2.- Directory where the output will be stored"            << endl;
        cout << "3.- Are we reading multiple configurations? (1-yes 0-no)" << endl;
        cout << "4.- If yes to (3), which configuration do you want to read? If not to (3) this argument doesn't mean anything" << endl;
        cout << "5.- Starting timestep" << endl; 
        cout << "6.- The dump frequency" <<endl;
        cout << "7.- Last timestep of the simulation" << endl;
        cout << "8.- number of nucleotides" << endl;       
        cout << "9.- oxDNA version" << endl;
        cout << "10.- Visualization (1 for vmd - 2 for ovito):" << endl;
    }


if(argc==10)
{   
    cout << "ok" << endl;
    /*argv[1] Directory where the data is*/
    /*argv[2] Directory where the output will be stored*/

    /*Multiple configurations?*/
    int multipleconf;
    multipleconf = atoi(argv[3]);
    
    /*Which configuration*/
    int config;
    config = atoi(argv[4]);
    
    /*Start the calculation from this timestep*/
    int equiltime;
    equiltime = atoi(argv[5]);
    
    /*Dump frequency*/
    int dumpfreq;
    dumpfreq = atoi(argv[6]);
    
    /*Last timestep*/
    int totrun;
    totrun = atoi(argv[7]);
    
    /*number of nucleotides*/  
    int N;
    N = atoi(argv[8]);

    /*Choose the oxDNA version you are using to compute the hb and bb location*/
    int version;         
    version = atoi(argv[9]); 
                               
                          
                           
    int frames = 1 + (totrun-equiltime)/dumpfreq;
   
    int num1, num2, num3;    

    
    /*the following variables will be useful to read the data file*/
    int id, mol, types, nx, ny, nz;
    double x, y, z, q0, q1, q2, q3;
 
    double Lx, Ly, Lz;
    double Lmaxx,Lminx,
           Lmaxy,Lminy,
           Lmaxz,Lminz;

    vector< vector< vector<double> > > quat(frames, vector<vector<double> >(4, vector<double>(N)));
    vector< vector< vector<double> > > com (frames, vector<vector<double> >(3, vector<double>(N)));
    vector< vector< vector<double> > > hb  (frames, vector<vector<double> >(3, vector<double>(N)));
    vector< vector< vector<double> > > bb  (frames, vector<vector<double> >(3, vector<double>(N)));

    vector<vector <int> > ids(frames, vector<int>(N));
    vector<vector <int> > type(frames, vector<int>(N));

  
    /*Usefull when reading input files*/
    int  timestep;
    char name0[10]    = "conf";
    char name1[10]    = "out";
    char name11[10]    = "_out";
    char name2[10]    = ".data";
    char readFile[250] = {'\0'};
    string dummy; 



    /*Usefull when writing output files*/

    string line;
    /*********************/
    /**Reads input files**/
    /*********************/
    //If the simulation starts from a single initial configuration and you produce out0.dat, out1.data, ...
    if(multipleconf==0)
    {
       sprintf(readFile, "%sout.data",argv[1]);
    }

    //This is useful when you gave sevaral configurations and you want to track only one.
    //For example in the simulations after releasing twist
    if(multipleconf==1)
    {
        sprintf(readFile, "%s%s%d%s%s",argv[1], name0, config, name11, name2);
    }
    
    ifstream read(readFile);
    if(!read){cout << "Error while opening data file!"<<endl;}
    
    //Inside the file, the number of rows assigned per time-frame is: 9, N;
    int nlinespf=9+N;

    if (read.is_open())    
    {  
        for (int t=0; t<frames; t++)
        {
        
            for(int i=0; i<nlinespf; i++)   
            {
                //Read the first 9 lines
                if(i<9)
                {                      
                    //Avoid the first 5 lines containing the timestep, number of particle and boundary conditions.
                    if(i<5){
                        getline(read,line); 
                        //cout << "i=" << i << "    "<< line << endl;
                    }
                
                    //Read The box size
                    if(i==5)
                    {
                        read >> Lminx >> Lmaxx;
                        //cout << "i=" << i << "  Lx " << Lminx << " " << Lmaxx <<endl;
                        Lx = Lmaxx-Lminx;
                    }

                    if(i==6)
                    {
                        read >> Lminy >> Lmaxy;
                        //cout << "i=" << i << "  Ly " << Lminy << " " << Lmaxy <<endl;
                        Ly = Lmaxy-Lminy;
                    }

                    if(i==7)
                    {
                        read >> Lminz >> Lmaxz;
                        //cout << "i=" << i << "  Lz " << Lminz << " " << Lmaxz <<endl;
                        Lz = Lmaxz-Lminz;
                        getline(read,dummy);
                    }
                    
                    if(i>7){
                        getline(read,line); 
                        //cout << "i=" << i << "    "<< line << endl;
                    }
                }
                            
                    
                
                // Continue reading the file
                if(i>=9)
                { 
                    read >> id >> mol >> types >> x >> y >> z >> nx >> ny >> nz >> q0 >> q1 >> q2 >> q3;

                    ids[t][id-1] = id-1;
                    type[t][id-1] = types;
                    com[t][0][id-1] = x+(Lx*nx);
                    com[t][1][id-1] = y+(Ly*ny);
                    com[t][2][id-1] = z+(Lz*nz);

                    quat[t][0][id-1] = q0;
                    quat[t][1][id-1] = q1;
                    quat[t][2][id-1] = q2;
                    quat[t][3][id-1] = q3;
                }
            }
            getline(read,dummy);
        }
    }
    read.close();        

//Check
/*
    for (int t=1; t<frames; t++)
    {
        for (int i=0; i<N; i++)
        {
            int j = N-1-i;
            cout << t << " " << ids[t][i] << endl;
        }
    }
*/


    /**************************************/
    /**Compute HB and phosphate positions**/
    /**************************************/
    double x1, y1, z1, q01, q11, q21, q31;

    for (int t=0; t<frames; t++)
    {
        for(int i=0; i<N; i++ )
        {
            x1  =  com[t][0][i]; y1  =  com[t][1][i]; z1  =  com[t][2][i];
            q01 = quat[t][0][i]; q11 = quat[t][1][i]; q21 = quat[t][2][i]; q31 = quat[t][3][i];

            vector<double> qs1 {q01, q11, q21, q31};
            vector<double> ex1(3), ey1(3), ez1(3);
            ex1 = q_to_ex(qs1); ey1 = q_to_ey(qs1); ez1 = q_to_ez(qs1);


            //oxdna1 backbone
            if(version==1)
            {
                bb[t][0][i] = x1 - 0.40*ex1[0];
                bb[t][1][i] = y1 - 0.40*ex1[1];
                bb[t][2][i] = z1 - 0.40*ex1[2];                
            }

            //oxdna2 backbone
            if(version==2)
            {
                bb[t][0][i] = x1 - 0.34*ex1[0] + 0.3408*ey1[0];
                bb[t][1][i] = y1 - 0.34*ex1[1] + 0.3408*ey1[1];
                bb[t][2][i] = z1 - 0.34*ex1[2] + 0.3408*ey1[2];
            }


            //The hydrogen bond sites
            hb[t][0][i] = x1 + 0.40*ex1[0];
            hb[t][1][i] = y1 + 0.40*ex1[1];
            hb[t][2][i] = z1 + 0.40*ex1[2];            
        }          
    }



    /********************************************************************/
    /*Compute the COM of the three arms near to the core of the molecule*/
    /********************************************************************/
    vector<vector <double> > core(frames, vector<double>(3));
    
    //First bp at the beggining of arm1: 27-79
    int id1a1=27;
    int id2a1=75;

    //First bp at the beggining of arm2: 76-128
    int id1a2=74;
    int id2a2=122;
    
    //First bp at the beggining of arm3: 30-125
    int id1a3=28;
    int id2a3=121;
    
    for (int t=0; t<frames; t++)
    {
        vector <double> c1(3);
        vector <double> c2(3);
        vector <double> c3(3);
        
        //Beggining of the three dsDNA arms (at the core)
        for(int d=0; d<3; d++)
        {            
            c1[d] = (com[t][d][id1a1-1] + com[t][d][id2a1-1])/2.0;
            c2[d] = (com[t][d][id1a2-1] + com[t][d][id2a2-1])/2.0;
            c3[d] = (com[t][d][id1a3-1] + com[t][d][id2a3-1])/2.0;
        }
                
        //COM of the molecule at the core
        for(int d=0; d<3; d++)
        { 
            core[t][d] = (c1[d] + c2[d] + c3[d])/3.0;
        }
    
    }
    
    /**********************************************/
    /*Compute the COM of the three arms at the end*/
    /**********************************************/ 
    //The unitary vector pointing from the core of the molecule til the last base-pair in each arm
    vector<vector <double> > dsdna1(frames, vector<double>(3));
    vector<vector <double> > dsdna2(frames, vector<double>(3));
    vector<vector <double> > dsdna3(frames, vector<double>(3));
      
    //First bp at the end of arm1: 8-98
    int id1a1end=8;
    int id2a1end=94;

    //First bp at the end of arm2: 57-147
    int id1a2end=55;
    int id2a2end=141;
    
    //First bp at the end of arm3: 49-106
    int id1a3end=47;
    int id2a3end=102;

    for (int t=0; t<frames; t++)
    {
        //End of the three dsDNA arms
        vector <double> c1end(3);
        vector <double> c2end(3);
        vector <double> c3end(3);
     
        for(int d=0; d<3; d++)
        {            
            c1end[d] = (com[t][d][id1a1end-1] + com[t][d][id2a1end-1])/2.0;
            c2end[d] = (com[t][d][id1a2end-1] + com[t][d][id2a2end-1])/2.0;
            c3end[d] = (com[t][d][id1a3end-1] + com[t][d][id2a3end-1])/2.0;
        }
                
        //The vector pointing from the core of the molecule to the end of the arms defines the direction of the arms
        vector <double> v1(3);
        vector <double> v2(3);
        vector <double> v3(3);
     
        for(int d=0; d<3; d++)
        {
            v1[d] = c1end[d] - core[t][d];
            v2[d] = c2end[d] - core[t][d];
            v3[d] = c3end[d] - core[t][d];
        }
        
        //The magnitude of the previous vectors
        double dd1=sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
        double dd2=sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
        double dd3=sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
        
        //The unitary vector representing the arms:
        for(int d=0; d<3; d++)
        {
            dsdna1[t][d] = v1[d]/dd1;
            dsdna2[t][d] = v2[d]/dd2;
            dsdna3[t][d] = v3[d]/dd3;
        }
    }

    /******************************/
    /*Volume of the parallelepiped*/
    /******************************/
    //Volume=|(u1 x u2).u3|
    vector<double> volume(frames);
    
    for (int t=0; t<frames; t++)
    {
        double ax=dsdna1[t][0];   double ay=dsdna1[t][1];   double az=dsdna1[t][2];
        double bx=dsdna2[t][0];   double by=dsdna2[t][1];   double bz=dsdna2[t][2];
        double cx=dsdna3[t][0];   double cy=dsdna3[t][1];   double cz=dsdna3[t][2];
        
        vector<double> A1 {ax, ay, az};
        vector<double> A2 {bx, by, bz};
        vector<double> A3 {cx, cy, cz};
        
        vector<double> cp(3);   
        cp=crossProduct(A1,A2);
        
        double dot=dotProduct(cp,A3);
        
        volume[t] = sqrt(dot*dot);
    }
    


    /**************************************************************************/
    /*Compute the distance between nucleotides forming a base-pair at the core*/
    /**************************************************************************/
    vector<double> dist1(frames);
    vector<double> dist2(frames);
    vector<double> dist3(frames);

    for (int t=0; t<frames; t++)
    {
        double l1=0.0;
        double l2=0.0;
        double l3=0.0;

        for(int d=0; d<3; d++)
        {
            l1 += pow(com[t][d][id1a1-1] - com[t][d][id2a1-1],2.0);
            l2 += pow(com[t][d][id1a2-1] - com[t][d][id2a2-1],2.0);
            l3 += pow(com[t][d][id1a3-1] - com[t][d][id2a3-1],2.0);
        }
        
        l1=sqrt(l1);
        l2=sqrt(l2);
        l3=sqrt(l3);
        
        dist1[t]=l1;
        dist2[t]=l2;
        dist3[t]=l3;
    }
    
    /************************************************************/
    /*Compute distance from the plane to the COM of the molecule*/
    /************************************************************/
    vector<double> distance(frames);
    for (int t=0; t<frames; t++)
    {
        double ax=dsdna1[t][0];   double ay=dsdna1[t][1];   double az=dsdna1[t][2];
        double bx=dsdna2[t][0];   double by=dsdna2[t][1];   double bz=dsdna2[t][2];
        double cx=dsdna3[t][0];   double cy=dsdna3[t][1];   double cz=dsdna3[t][2];
        
        //Vector pointing from U1 to U2
        vector<double> A1 {bx-ax, by-ay, bz-az};
        
        //Vector pointing from U1 to U3
        vector<double> A2 {cx-ax, cy-ay, cz-az};
        
        //The cross product (A1 x A2) gives the normal to the plane
        vector<double> cp(3);   
        cp=crossProduct(A1,A2);
        
        //general equation of the plane a*x + b*y + c*z = dp (with a=cp[0], b=cp[1], c=cp[2] and x, y, z one point in the plane (for example, the end of one of the unitary vectors, U2))
        double dp = cp[0]*bx + cp[1]*by + cp[2]*bz;
        
        //The distance can be computed as: (D = ax0 + by0 + cz0 + dp)/(sqrt(a^2 + b^2 + c^2)), where (x0=0, y0=0, z0=0) is the coordinate of the COM.
        double d0 = dp/sqrt(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2]);

        double d0abs=sqrt(d0*d0);
        
        distance[t] = d0abs;
    }
    
    /****************/
    /*Compute angles*/
    /****************/
    vector<double> theta12(frames);
    vector<double> theta13(frames);
    vector<double> theta23(frames);
    for (int t=0; t<frames; t++)
    {
        double ax=dsdna1[t][0];   double ay=dsdna1[t][1];   double az=dsdna1[t][2];
        double bx=dsdna2[t][0];   double by=dsdna2[t][1];   double bz=dsdna2[t][2];
        double cx=dsdna3[t][0];   double cy=dsdna3[t][1];   double cz=dsdna3[t][2];
     
        //unitary dsDNA arm1
        vector<double> A1 {ax, ay, az};
        
        //unitary dsDNA arm2
        vector<double> A2 {bx, by, bz};
        
        //unitary dsDNA arm3
        vector<double> A3 {cx, cy, cz};
        
        //DOt products
        double ctheta12 = dotProduct(A1,A2);
        double ctheta13 = dotProduct(A1,A3);
        double ctheta23 = dotProduct(A2,A3);
        
        theta12[t] = acos(ctheta12)*180./pi;
        theta13[t] = acos(ctheta13)*180./pi;
        theta23[t] = acos(ctheta23)*180./pi;
     
     }
    
    

    //PRINT RESULTS
    char writeFile1 [250] = {'\0'};
    char name3[100] ="colvar_vs_time";
    sprintf(writeFile1, "%sN%d.oxDNA%d.%s", argv[2], N, version, name3);
        
    ofstream write1(writeFile1);
    cout << "writing on .... " << writeFile1 <<endl;
        
    //set precision and the number of decimals to be printed always
    write1.precision(16);
    write1.setf(ios::fixed);
    write1.setf(ios::showpoint);

    for (int t=0; t<frames; t++) 
    { 
        timestep = equiltime+dumpfreq*t;
        write1 << timestep << " " << distance[t] << " " << volume[t] << " " << theta12[t] << " " << theta13[t] << " " << theta23[t] << " " << dist1[t] << " " << dist2[t] << " " << dist3[t] << endl;
    }
    write1.close();
    write1.clear();
    




}








return 0;
}



//Quaternion product
vector<double> quatproduct(vector<double> P, vector<double> Q)
{
    vector<double> R(4);
    R[0] = P[0]*Q[0] - P[1]*Q[1] - P[2]*Q[2] - P[3]*Q[3];
    R[1] = P[0]*Q[1] + Q[0]*P[1] + P[2]*Q[3] - P[3]*Q[2];
    R[2] = P[0]*Q[2] + Q[0]*P[2] - P[1]*Q[3] + P[3]*Q[1];
    R[3] = P[0]*Q[3] + Q[0]*P[3] + P[1]*Q[2] - P[2]*Q[1];

    return R;
}

// converts quaternion DOF into local body reference frame
vector<double> q_to_ex(vector<double> Q)
{
    vector<double> ex(3);

    ex[0]=Q[0]*Q[0]+Q[1]*Q[1]-Q[2]*Q[2]-Q[3]*Q[3];
    ex[1]=2*(Q[1]*Q[2]+Q[0]*Q[3]);
    ex[2]=2*(Q[1]*Q[3]-Q[0]*Q[2]);
    return ex;
}


vector<double> q_to_ey(vector<double> Q)
{
    vector<double> ey(3);

    ey[0]=2*(Q[1]*Q[2]-Q[0]*Q[3]);
    ey[1]=Q[0]*Q[0]-Q[1]*Q[1]+Q[2]*Q[2]-Q[3]*Q[3];
    ey[2]=2*(Q[2]*Q[3]+Q[0]*Q[1]);
    return ey;
}


vector<double> q_to_ez(vector<double> Q)
{
    vector<double> ez(3);

    ez[0]=2*(Q[1]*Q[3]+Q[0]*Q[2]);
    ez[1]=2*(Q[2]*Q[3]-Q[0]*Q[1]);
    ez[2]=Q[0]*Q[0]-Q[1]*Q[1]-Q[2]*Q[2]+Q[3]*Q[3];
    return ez;
}


//Sign function
double sign(double a) { if (a>0) {return 1.0;} else {return -1.0;} } 


//length
double length(vector<double> a)
{
  double r;
  r = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  return r;
}

//Dot product
double dotProduct (vector<double> a, vector<double> b)
{
  double r;
  r = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  return r;
}

//Cross product
vector<double> crossProduct(vector<double> a, vector<double> b)
{
  vector<double> r(3);
  r[0] = a[1]*b[2]-a[2]*b[1];
  r[1] = a[2]*b[0]-a[0]*b[2];
  r[2] = a[0]*b[1]-a[1]*b[0];
  
  return r;
}


// checks for round off in a number in the interval [-1,1] before using asin or acos
void fix_roundoff(double &a)
{
  if (abs(a)>1.001) 
  {
    std::cout<<"Error - number should be in interval [-1,1] but is "<< a << endl;
    exit(0);
  }
  if (a>1){a=1.0;}
  if (a<-1){a=-1.0;}
}


//The orthogonal matrix. This is the rotation matrix transforming the canonical frame into the frame of the respective triad
vector<vector <double> > matrixOrthogonal(vector<double> a, vector<double> b, vector<double> c)
{
  vector<vector <double> > T(3, vector<double>(3));

  T[0][0] = a[0];    T[0][1] = b[0];    T[0][2] = c[0];
  T[1][0] = a[1];    T[1][1] = b[1];    T[1][2] = c[1];
  T[2][0] = a[2];    T[2][1] = b[2];    T[2][2] = c[2];

  return T;
}


vector<vector <double> > matrixTranspose(vector<vector <double> > a)
{
  vector<vector <double> > M(3, vector<double>(3));

  for(int i=0; i<3; i++)
  {
      for(int j=0; j<3; j++)
      {
          M[i][j] = a[j][i];
      }
  }

  return M;
}


vector<vector <double> > matrixProduct(vector<vector <double> > a, vector<vector <double> > b)
{
  vector<vector <double> > M(3, vector<double>(3));

  for(int i=0; i<3; i++)
  {
      for(int j=0; j<3; j++)
      {
          double s=0;
          
          for(int k=0; k<3; k++)
          {
              s = s + a[i][k]*b[k][j];

              M[i][j] = s;
          }
      }
  }

  return M;
}
     

double matrixTrace(vector<vector <double> > a)
{
  double r;
  r = a[0][0] + a[1][1] + a[2][2];
  return r;
}









