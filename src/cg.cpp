#include "cg.h"

//r   f - A  u
// delta0 =  r * r
// d = r
void cg::computeFirstStep(const Grid& u,const Grid& rhs,Grid& res,Grid& d,real& para_delta0,const uint& para_nx,const uint& para_ny,const real& para_invHx,const real& para_invHy,const real& para_weightC)
{
	for(uint i = 1; i < para_ny-1; ++i)
	{
		for (uint j = 1; j < para_nx-1; j++)
		{
			res(i,j) = rhs(i,j) + para_invHx * (u(i,j-1) + u(i,j+1)) + para_invHy * (u(i-1,j) + u(i+1,j)) - para_weightC*u(i,j);
			d(i,j) = res(i,j);
			para_delta0 += res(i,j) * res(i,j);
		}
	}
}

// z = A * d
// alpha = d * A * d or d * z
void cg::multiplyVector(const Grid& d,Grid& z,real& para_alpha,const uint& para_nx,const uint& para_ny,const real& para_invHx,const real& para_invHy,const real& para_weightC)
{
	for(uint i = 1; i < para_ny-1; ++i)
	{
		for (uint j = 1; j < para_nx-1; j++)
		{
			z(i,j) =  (-1) * para_invHx * (d(i,j-1) + d(i,j+1)) - para_invHy * (d(i-1,j) + d(i+1,j)) + para_weightC*d(i,j);
			para_alpha += d(i,j) * z(i,j);
		}
	}
}

void cg::updateU(Grid& u,const Grid& d,const real& para_alpha,const uint& para_nx,const uint& para_ny)
{
	for(uint i = 1; i < para_ny-1; ++i)
	{
		for (uint j = 1; j < para_nx-1; j++)
		{
			u(i,j) =  u(i,j) + para_alpha * d(i,j);
		}
	}
}

void cg::updateRes(Grid& res,const Grid& z,const real& para_alpha,real& para_delta1,const uint& para_nx,const uint& para_ny)
{
	for(uint i = 1; i < para_ny-1; ++i)
	{
		for (uint j = 1; j < para_nx-1; j++)
		{
			res(i,j) =  res(i,j) - para_alpha * z(i,j);
			para_delta1 += res(i,j) * res(i,j);
		}
	}
}

void cg::updateD(Grid& d,const Grid& res,const real& para_beta,const uint& para_nx,const uint& para_ny)
{
	for(uint i = 1; i < para_ny-1; ++i)
	{
		for (uint j = 1; j < para_nx-1; j++)
		{
			d(i,j) =  res(i,j) + para_beta * d(i,j);
		}
	}
}

//output : number in grid x / processor , same as y
void cg::findElementsinXY(int& elementX, int& elementY, const uint& nx, const uint& ny,const int& numProcX,const int& numProcY,const int& cartCordX,const int& cartCordY)
{
    elementX = nx / numProcX;
    if(cartCordX == numProcX-1)
    {
      elementX += nx % numProcX;
     // std::cout<<"coord in X "<<cartCordX<<"\t coord in Y "<<cartCordY<<"\t elementX "<<elementX<<"\t elementY  "<<elementY<<std::endl;
    }
    
    elementY = ny / numProcY;
    if(cartCordY == numProcY-1)
    {
      elementY += ny % numProcY;
      //std::cout<<"coord in X "<<cartCordX<<"\t coord in Y "<<cartCordY<<"\t elementX "<<elementX<<"\t elementY  "<<elementY<<std::endl;
    }
}

void cg::findGlobalCoord(int& xStart,int& xEnd,int& yStart,int& yEnd,const int& elementX,const int& elementY,const int& mycoordX,const int& mycoordY,const int& numProcX,const int& numProcY,const uint& nx,const uint& ny)
{
    if(mycoordX == numProcX-1)
      {
	xStart = mycoordX * (elementX - (nx % numProcX));
	xEnd = (xStart + elementX)-1;
      }
    else
      {
	xStart = mycoordX * elementX;
	xEnd = xStart + elementX;
      }
  //-------------------------------------------------------------------------------------
    if(mycoordY == numProcY-1)
      {
	yStart = mycoordY * (elementY - (ny % numProcY));
	yEnd = (yStart + elementY)-1;
      }
    
    else
      {
	yStart = mycoordY * elementY;
	yEnd = yStart + elementY; 
      }
//we don't want to iterate at boundary points,so,
    if(mycoordX == 0)
    {
      xStart += 1;
    }
    /*
    if(mycoordX == numProcX-1)
    {
      xEnd -= 1;
    }*/
    
    if(mycoordY == 0)
    {
      yStart += 1;
    }
    /*
    if(mycoordY == numProcY-1)
    {
      yEnd -= 1;
    }*/
      
}

/* 
  dir = 0 , disp = 1 => + x-axis ,
  dir = 0 , disp = -1 => - x-axis,
  dir = 1 , disp = 1 => + y-axis , 
  dir = 1 , disp = -1 => - y-axis 
*/
void cg::sendReceiveMsg(Grid& u,int& elementX,int& elementY,MPI_Comm& GRID_COMM_CART,MPI_Datatype& colType,MPI_Datatype& rowType)
{
   
  int srcRank , dstRank;
  int count = 0 , myrank;
  MPI_Request requestHandle[8];
  MPI_Status statusHandle[8];
  
  //left , right , top , bottom
  int sendBuffer0[] = { 1 , 1  , elementY , 1} , sendBuffer1[] = { elementX , 1 , 1 , 1};
  int recvBuffer0[] = { 1 , 1 , 0 , elementY+1} , recvBuffer1[] = { 0 , elementX+1 , 1 , 1}; 
  
  int tag[] = { 1 , 2 , 3 , 4 }; //for left , right , top , bottom
  
  //we want to transfer/receive data in x and y direction
  for(int dir = 0 ; dir < ndims ; ++dir)
  {
    for(int disp = 1 ; disp > -2 ; disp -= 2)
    {
	MPI_Cart_shift(GRID_COMM_CART,dir,disp,&srcRank,&dstRank);
	
	MPI_Comm_rank(GRID_COMM_CART,&myrank);
	
	
// count < 2 means first, we shift in x-axis.	
//---------------------------------------------------------------------------------------------------------------------------------	
	if(count < 2)	// column shift
	{ 
	  
	  if(srcRank != MPI_PROC_NULL)
	    {
	      
	      // receive message from left boundary
	      MPI_Irecv( &u(recvBuffer0[count],recvBuffer1[count]),1,colType,srcRank,tag[count],GRID_COMM_CART,&requestHandle[count]);
	    }
	  
	  if(dstRank != MPI_PROC_NULL)
	  {
	    //std::cout<<"srcRank "<<srcRank<<"\tdstRank "<<dstRank<<"\t"<<myrank<<"  I am here"<<std::endl;
	    MPI_Isend(&u(sendBuffer0[count],sendBuffer1[count]),1,colType,dstRank,tag[count],GRID_COMM_CART,&requestHandle[count+4]);
	    
	  }
	}
//-------------------------------------------------------------------------------------------------------------------------------------	  
	else	// row shift
	{
	  if(srcRank != MPI_PROC_NULL)
	    {
	      // receive message from left boundary
	      MPI_Irecv( &u(recvBuffer0[count],recvBuffer1[count]),1,rowType,srcRank,tag[count],GRID_COMM_CART,&requestHandle[count]);
	    }
	  
	  if(dstRank != MPI_PROC_NULL)
	  {
	    MPI_Isend(&u(sendBuffer0[count],sendBuffer1[count]),1,rowType,dstRank,tag[count],GRID_COMM_CART,&requestHandle[count+4]);
	  }
	}  
	   
//-------------------------------------------------------------------------------------------------------------------------------------

      if(srcRank != MPI_PROC_NULL)
      {
	MPI_Wait(&requestHandle[count],&statusHandle[count]);
      }
	
      if(dstRank != MPI_PROC_NULL)
      {
	MPI_Wait(&requestHandle[count+4],&statusHandle[count+4]);
      }	
      count++;
    } // disp
    
  } // dir

  
} // function



