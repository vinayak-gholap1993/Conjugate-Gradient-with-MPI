
#include"type.h"
#include"Grid.h"
#include"cg.h"

void write( const Grid& para_storeToFile,const uint& para_nx,const uint& para_ny,const real& para_hx,const real& para_hy)
{
  std::ofstream file;
  file.open("finalU.txt");
  for(uint i = 0; i < para_ny ; ++i)
    {   
	uint term = i * para_nx; 
	for(uint j = 0; j < para_nx ; ++j)
	{
	  file << i*para_hy <<" "<< j*para_hx <<" "<< para_storeToFile.data[term + j] <<"\n";
	 }
	 file<<"\n";
    }
     file.close();
}


int main(int argc, char** argv) 
{
	  
	MPI_Init( &argc , &argv );
	MPI_Comm GRID_COMM_CART;
	
	int numProc , rankProc;
	
	MPI_Comm_size(MPI_COMM_WORLD,&numProc);
	MPI_Comm_rank(MPI_COMM_WORLD,&rankProc);
	

	int proc_dims[ndims] = { 0 , 0 } , periods[ndims] = {0,0} , mycoord[ndims];
	
	MPI_Dims_create(numProc, ndims ,proc_dims); 
	
//	std::cout<<"Number of node in x "<<proc_dims[0]<<"\tNumber of node in y "<<proc_dims[1]<<std::endl;
	
	MPI_Cart_create(MPI_COMM_WORLD,ndims,proc_dims,periods,0,&GRID_COMM_CART);
	
	//grid_... for GRID_COMM_CART
	int grid_proc , grid_rank;
	MPI_Comm_size(GRID_COMM_CART,&grid_proc);
	MPI_Comm_rank(GRID_COMM_CART,&grid_rank);
	
	MPI_Cart_coords(GRID_COMM_CART,grid_rank,ndims,mycoord);

	int numProcX = proc_dims[0] , numProcY = proc_dims[1];
	int cartCordX = mycoord[0] , cartCordY = mycoord[1];
//-----------------------------------------------------------------------------------------------

	///////
	if(argc!=4)	{;}
	
	uint nx = std::atoi(argv[1]);
	uint ny = std::atoi(argv[2]);

	//const uint nx = 10;
	//const uint ny = 10;
	
	const real hx = (x1 - x0) / (nx-1);
	const real hy = (y1 - y0) / (ny-1);

	const real invHX = 1.0 / (hx * hx);
	const real invHY = 1.0 / (hy * hy);
	const real weightC = 2.0 * invHX + 2.0 * invHY + M_PI * M_PI;

	const uint totalGrid = nx * ny;

	//std::cout<<"hx is... "<<hx<<"\t hy is..."<<hy<<std::endl;

	Grid u_exact(nx,ny);
	Grid u(nx,ny);
	Grid rhs(nx,ny);
	Grid res(nx,ny);
	Grid d(nx,ny);
	Grid z(nx,ny);

	cg cgMPI;
	
	int elementX = 0,elementY = 0 , xStart = 0, xEnd = 0 , yStart = 0 , yEnd = 0;
	
	cgMPI.findElementsinXY(elementX,elementY,nx,ny,numProcX,numProcY,cartCordX,cartCordY);
	cgMPI.findGlobalCoord(xStart,xEnd,yStart,yEnd,elementX,elementY,cartCordX,cartCordY,numProcX,numProcY,nx,ny);
	
	//create vector in contigous in memory
	MPI_Datatype rowType,colType;
	MPI_Type_contiguous(elementX,MPI_DOUBLE,&rowType);
	MPI_Type_commit(&rowType);
	MPI_Type_vector(elementY,1,0,rowType,&colType);
	MPI_Type_commit(&colType);
	
	
	if(rankProc == 0)
	{
	  std::cout<<"xstart "<<xStart<<"\t xEnd "<<xEnd<<"\tystart "<<yStart<<"\t yEnd "<<yEnd<<std::endl;
	}

	u_exact.applyBoundaryU(u_exact,nx,ny,hx);
	rhs.applyBoundaryRHS(rhs,nx,ny,hx,hy);
	

	real delta0 = 0.0 , alpha = 0.0 , delta1 = 0.0 , beta = 0.0;

	cgMPI.computeFirstStep(u,rhs,res,d,delta0,nx,ny,invHX,invHY,weightC);

	for(int iter = 1; iter < std::atoi(argv[3]) ; ++iter)
	{
		cgMPI.multiplyVector(d,z,alpha,nx,ny,invHX,invHY,weightC);
		alpha = delta0 / alpha;
		cgMPI.updateU(u,d,alpha,nx,ny);
		cgMPI.updateRes(res,z,alpha,delta1,nx,ny);
		
		std::cout<<"Rank... "<<rankProc<<"\tresidual norm ... "<<sqrt(delta1/totalGrid)<<std::endl;
		beta = delta1 / delta0;
		cgMPI.updateD(d,res,beta,nx,ny);
	
		delta0 = delta1;
		alpha = 0.0 , delta1 = 0.0 , beta = 0.0;
	}

	//write(u,nx,ny,hx,hy);
	
	MPI_Finalize();
	
return 0;
}
