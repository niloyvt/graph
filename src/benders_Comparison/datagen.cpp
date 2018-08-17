// datagen.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include "datalap.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
	/*
	The inequalities obtained from the Birkhoff Neumann decomposition can be used to solve the ALSP using Bender's decomposition. However, we found the computational performance of such a method to be vastly inferior to that of the A-formulation. Consider the following formulation:
\begin{eqnarray}
minimize & z\label{eq:LO_obj}
\end{eqnarray}
\begin{eqnarray}
\mbox{subject to:}\eqref{eq:Total_order},\eqref{eq:Transitivity},\eqref{eq:assignment1},\eqref{eq:assign2},\\
\sum_j (j-1)\pi_{ij}=\sum_j y_{ji} &\forall i=1,\cdots,l.
\end{eqnarray}
To this, we can add the cuts \eqref{eq:cuts} with the left side replaced with $z$. While such a formulation converges to the optimal, the solution time required is large. We fixed the upper bound on the number of sublots to 3. For a problem with 10 lots and 10 vendors, the solution time was 6.56s as against 1.34s for the A-formulation. Increasing the number of vendors to 20 results in a CPU-time of 13.47s compared to 1.28s for the A-formulation. The Bender's approach does not obtain an optimal solution within the specified time limit when the number of lots is increased to 20, irrespective of the number of vendors used. This is significantly worse than any of the other formulations discussed.
*/

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
		datalap obj;
		cout<<"Hello Seattle!";
		
	obj.experiment(20,20);

	//	vector<double> difference;
	//	vector<int> badroots;
	//for(int root =14; root< 15; ++root)
	//{
	//	
	//	double nb, b;
	////obj.generatedata(root,4,20,2);
	//	obj.generatedata(root,4,20,2);
	////for obj.generatedata(root,4,20,2); found bad roots 14 (4), 19 (1), and 46 (3) with errors shown in ().
	//
	//
	//
	//
	//
	//

	//nb = obj.otimalLObased();
	//cout<<"optimal sequence LO based \n";
	//obj.printprec();
	//cout<<endl;
	//cout<<"For this sequence givenseq= "<<obj.givensequence()<<endl;
	//
	//cout<<"Tested optimal= "<<nb<<endl;
	////cin>>b;
	//b = obj.benders();
	//cout<<b;

	//cout<<"master - nb= "<<b-nb<<endl;
	//cout<<root<<endl;
	//
	//if(0.1 < abs(b-nb))
	//{
	//	difference.push_back(abs(b-nb));
	//	badroots.push_back(root);
	//}

	////assert(b-nb < 1);
	////assert(nb-b < 1);

	//}

	//for(unsigned short i=0; i<difference.size(); ++i)
	//	cout<<"difference= "<<difference[i]<<" for bad root ="<<badroots[i]<<endl;






	//obj.lpNB();
	
	//obj.otimalLObased();
	////obj.sillygiven();
	//bool foundone=1;
	//int root=0;
	//
	//while(foundone && root<2)
	//{
	//	root++;
	//	obj.generatedata(root,4,2,2);
	//	
	//	double badheuristic=obj.lpformulation(obj.precedence);		
	//	double goodheuristic = obj.otimalLObased();
	//	cout<<"root="<<root<<endl;

	//	
	//	if(badheuristic-goodheuristic>0.01 || badheuristic-goodheuristic<-0.01 )//added tolerance
	//	{
	//		cout<<string(50,'\n');
	//		cout<<endl<<"root="<<root<<endl;
	//		cout<<"correct answer"<<goodheuristic<<endl;
	//		cout<<"wrong answer"<<badheuristic<<endl;
	//		foundone=0;
	//	}
	//}
	//if(!foundone)
	//	cout<<"different answers\n";
	//else
	//	cout<<"both heuristics give the same answer";

	//
	/*obj.experiment(10,2);
	obj.experiment(20,2);
	obj.experiment(30,2);
	obj.experiment(40,2);

	obj.experiment(10,20);
	obj.experiment(20,20);
	obj.experiment(30,20);
	obj.experiment(40,20);

	obj.experiment(10,30);
	obj.experiment(20,30);
	obj.experiment(30,30);
	obj.experiment(40,30);

	obj.experiment(10,40);
	obj.experiment(20,40);
	obj.experiment(30,40);
	obj.experiment(40,40);*/

//	ofstream resultfile;
//	resultfile.open("otto.bat",ios::out);//,std::fstream::app);
//	if(!resultfile)
//		cout<<"Could not open";
//	resultfile<<"Hello World!";
//	resultfile.close();
//	
//	
//datalap obj;
//obj.generatedata(1,12,20,10);
//cout<<"Heuristic value="<<obj.givensequence()<<endl;
//cout<<obj.optimalNetwork();
//obj.optimalNetwork();
//cout<<"Heuristic value="<<obj.heuristicLPbased()<<endl;
//try
//{
////obj.experiment(25,10);//was 21
//	obj.generatedata(1,19,15,10);
//	obj.lpLO();
//	obj.lpNB();
//	obj.optimalNetwork();
////throw 1;
//}
//catch (int e)
//{
//	cout<<"experient failed";
//}
//obj.experiment(10,20);



	return 0;
}

