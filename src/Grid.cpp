#include "Grid.h"



//constructor
Grid::Grid(const uint& para_totalGrid)
{
	this->data = new real[para_totalGrid];
	
	for(uint i = 0; i < para_totalGrid ; ++i)
	{
		this->data[i] = 0.0;
	}
}

Grid::Grid(const uint& para_nx,const uint& para_ny)
{
	this->nx = para_nx;	this->ny = para_ny;

	const uint totalGrid = para_nx * para_ny;
	
	this->data = new real[totalGrid];
	for(uint i = 0; i < totalGrid ; ++i)
	{
		this->data[i] = 0.0;
	}

}

//apply boundary condition to U
void Grid::applyBoundaryU(Grid& para_u, const uint& para_xStart, const uint& para_xEnd,const uint& para_yStart, const uint& para_yEnd,const real& para_hx,const int& numProcY,const int& mycoordX,const int& mycoordY)
{
  
  // processes that touch bottom boundary 
    if(mycoordY == numProcY - 1)
    {  
	int value = 0;
	
	if(mycoordX == 0) {value = 1;}
      
	real SIN_H = std::sinh(m_pi);
	
	const uint term = para_yEnd - para_yStart;	// i * para_nx
	
	for(uint j = para_xStart, x = value ; j < para_xEnd ; j++, x++)
	{
		//para_u.data[term + j] = std::sin(m_pi * para_hx * j) * SIN_H;
	  para_u(term,x) = std::sin(m_pi * para_hx * j) * SIN_H;
	}

    }
}

//apply boundary condition to RHS
void Grid::applyBoundaryRHS(Grid& para_rhs, const uint& para_xStart, const uint& para_elementX,const uint& para_yStart, const uint& para_elementY,const real& para_hx,const real& para_hy,const int& mycoordX)
{
	const real termC =  m_pi * m_pi;	
	
	int value = para_xStart;
	//in y direction we consider boundary means yEnd = 10.
	if(mycoordX == 0)	{ value -= 1;	}
	for (uint y = 0 ,i = para_yStart; y < para_elementY; i++ , ++y)
	{
		uint term = std::sinh(m_pi * para_hy * i);
		
		for(uint j = value , x = 0 ; x < para_elementX ; j++ ,++x )
		{
        		para_rhs(y,x) = termC * std::sin(m_pi * para_hx * j) * term;
		}
	}
}

//write 
real& Grid::operator()(const uint& x,const uint& y)
{
		return this->data[x * this->nx + y];
}

//read only
real Grid::operator()(const uint& x,const uint& y) const
{
		return this->data[x * this->nx + y];
}

//print matrix just for test
void Grid::print(const Grid& para_print,const uint& para_elementX,const uint& para_elementY)
{
  std::cout<<"I am in print "<<para_elementX<<std::endl;
	for(uint i = 0; i < para_elementY ; ++i)
	{   
		uint term = i * para_elementX; 
		for(uint j = 0; j < para_elementX ; ++j)
		{
			std::cout<<para_print.data[term + j] << "\t";
		}
		std::cout<<"\n";
	}
}

/*
void Grid::writeFile(const std::string& fileName, const Grid& para_storeToFile,const uint& para_nx,const uint& para_ny,const real& para_hx,const real& para_hy)
{
  std::ofstream file(fileName);
  for(uint i = 0; i < para_ny ; ++i)
    {   
	uint term = i * para_nx; 
	for(uint j = 0; j < para_nx ; ++j)
	{
	  file << i*para_hy <<" "<< j*para_hx <<" "<< para_storeToFile.data[term + j] <<"\n";
	 }
	 file<<"\n";
    }
}
*/