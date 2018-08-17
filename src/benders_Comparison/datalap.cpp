#include "stdafx.h"
#include "datalap.h"
#include <stdlib.h>
#include <string>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <vector>
#include <assert.h>
#include <algorithm>
#include <cmath>
#include <ilcplex/ilocplex.h>

using namespace std;

typedef IloBoolVarArray BoolVarArray1;
typedef IloArray<BoolVarArray1> BoolVarArray2;
typedef IloArray<BoolVarArray2> BoolVarArray3;
typedef IloArray<BoolVarArray3> BoolVarArray4;

typedef IloNumVarArray FloatVarArray1;
typedef IloArray<FloatVarArray1> FloatVarArray2;
typedef IloArray<FloatVarArray2> FloatVarArray3;
typedef IloArray<FloatVarArray3> FloatVarArray4;
typedef IloArray<FloatVarArray4> FloatVarArray5;
typedef IloArray<FloatVarArray5> FloatVarArray6;
typedef IloArray<FloatVarArray6> FloatVarArray7;

typedef IloArray<IloRangeArray> IloRangeArray2;
typedef IloArray<IloRangeArray2> IloRangeArray3;

typedef IloArray<IloExprArray> IloExprArray2;

ILOUSERCUTCALLBACK4(BendersUserCallback, FloatVarArray2, y, BoolVarArray2, x, IloNumVar, Obj, IloNum, ObjConstant)
{
	
	unsigned short TLots = datalap::master_instance->TLots;
	IloEnv masterEnv = getEnv();


	datalap::master_instance->m_perm.resize(TLots);
	datalap::master_instance->precedence.resize(TLots);
	cout<<"perm received from master:"<<endl;
	for(unsigned short i=0;i<TLots;++i)
	{
		datalap::master_instance->m_perm[i].resize(TLots);
		datalap::master_instance->precedence[i].resize(TLots);

		for(unsigned short j=0;j<TLots;++j)
		{
			

			
			if(i==j)
				datalap::master_instance->precedence[i][j]=0;
			else
				datalap::master_instance->precedence[i][j]= getValue(y[i][j]);
			datalap::master_instance->m_perm[i][j]= getValue(x[i][j]);
			
			cout<<datalap::master_instance->m_perm[i][j];

			

		}
		cout<<endl;

	}	
	cout<<"perm from master:"<<endl;
	datalap::master_instance->printprec();
	datalap::master_instance->clean_perm(datalap::master_instance->m_perm);
	datalap::master_instance->perm_to_prec();
	double LowerBound = datalap::master_instance->givensequence();
	cout<<"Critical Lot="<<datalap::master_instance->critical_lot<<endl;
	cout<<"Critical Sublot="<<datalap::master_instance->critical_sublot<<endl;
	cout<<"Critical Vendor="<<datalap::master_instance->critical_vendor<<endl;
	cout<<LowerBound<<" from given seq.\n";

	IloExpr cutLhs(masterEnv);
	IloNumVar cutRhs = Obj;
	cout<<"obj constant in sub="<<datalap::master_instance->ObjConstant<<endl;
	cout<<"Calling Seperate()                    \n \n \n  WORLD \n\n"<<"\n\n\n Ta Daaaa\n";
	double subobj = 0;

	//subobj = datalap::master_instance->separate( cutLhs,
	//	x, y, LowerBound,
	//datalap::master_instance->critical_lot, 
	//datalap::master_instance->critical_sublot-1,
	//datalap::master_instance->critical_vendor);
	//cout<<"from sub = "<<subobj<<endl;
	//cout<<cutLhs<<endl;
	//assert(abs(subobj  -LowerBound)<0.1);

	vector<unsigned short> crsublotconfigs;
	crsublotconfigs.push_back(datalap::master_instance->submap_inv[datalap::master_instance->critical_lot]
	[datalap::master_instance->critical_vendor]
	[datalap::master_instance->critical_sublot-1]);
	
	subobj = datalap::master_instance->separatesmall( cutLhs,
		x, y, LowerBound, crsublotconfigs);

	cout<<"from sub = "<<subobj<<endl;
	cout<<cutLhs<<endl;
	assert(abs(subobj  -LowerBound)<0.1);

	cout<<"back to master \n";

	add(cutLhs+ datalap::master_instance->ObjConstant <= cutRhs).end();
	cutLhs.end();
	crsublotconfigs.clear();

	


	

	cout<<"HELLO Finished one call to separate()                     \n \n \n                                                         WORLD \n\n"<<"\n\n\n Ta Daaaa";
	
	/*cout<<NumSublots<<endl;
	IloEnv work_env= workerCplex.getEnv();
	IloNumArray init(work_env,l.getSize());
	for(unsigned short i=0;i<l.getSize();i++)
		init[i]=i;
	workerCplex.setVectors(init,0,l,0,0,0);
	workerCplex.setParam(IloCplex::PreInd, false);
	workerCplex.solve();
	cout<<"Dual="<<workerCplex.getDual(pickone)<<endl;
	for(unsigned short i=0;i<l.getSize();i++)
		cout<<workerCplex.getValue(l[i])<<endl;*/
	
	return;
}
ILOLAZYCONSTRAINTCALLBACK4(BendersLazyCallback, BoolVarArray2, y, BoolVarArray2, x, IloNumVar, Obj, IloNum, ObjConstant)
{
	/*
	unsigned short TLots = datalap::master_instance->TLots;
	IloEnv masterEnv = getEnv();
	datalap::master_instance->m_perm.resize(TLots);
	datalap::master_instance->precedence.resize(TLots);
	cout<<"perm received from master:"<<endl;
	for(unsigned short i=0;i<TLots;++i)
	{
		datalap::master_instance->m_perm[i].resize(TLots);
		datalap::master_instance->precedence[i].resize(TLots);

		for(unsigned short j=0;j<TLots;++j)
		{
			

			
			if(i==j)
				datalap::master_instance->precedence[i][j]=0;
			else
				datalap::master_instance->precedence[i][j]= getValue(y[i][j]);
			datalap::master_instance->m_perm[i][j]= getValue(x[i][j]);
			if(datalap::master_instance->m_perm[i][j]<0.1)
				datalap::master_instance->m_perm[i][j]=0;
			else
				datalap::master_instance->m_perm[i][j]=1;
			cout<<datalap::master_instance->m_perm[i][j];

			

		}
		cout<<endl;

	}	
	cout<<"perm from master:"<<endl;
	datalap::master_instance->printprec();
	//datalap::master_instance->clean_perm(datalap::master_instance->m_perm);
	//datalap::master_instance->perm_to_prec();
	double LowerBound = datalap::master_instance->givensequence();
	cout<<"Critical Lot="<<datalap::master_instance->critical_lot<<endl;
	cout<<"Critical Sublot="<<datalap::master_instance->critical_sublot<<endl;
	cout<<"Critical Vendor="<<datalap::master_instance->critical_vendor<<endl;
	cout<<LowerBound<<" from given seq.\n";

	IloExpr cutLhs(masterEnv);
	IloNumVar cutRhs = Obj;
	cout<<"obj constant in sub="<<datalap::master_instance->ObjConstant<<endl;
	cout<<"Calling Seperate()                    \n \n \n  WORLD \n\n"<<"\n\n\n Ta Daaaa\n";
	double subobj = 0;

	//subobj = datalap::master_instance->separate( cutLhs,
	//	x, y, LowerBound,
	//datalap::master_instance->critical_lot, 
	//datalap::master_instance->critical_sublot-1,
	//datalap::master_instance->critical_vendor);
	//cout<<"from sub = "<<subobj<<endl;
	//cout<<cutLhs<<endl;
	//assert(abs(subobj  -LowerBound)<0.1);

	vector<unsigned short> crsublotconfigs;
	crsublotconfigs.push_back(datalap::master_instance->submap_inv[datalap::master_instance->critical_lot]
	[datalap::master_instance->critical_vendor]
	[datalap::master_instance->critical_sublot-1]);
	
	subobj = datalap::master_instance->separatesmallbool( cutLhs,
		x, y, LowerBound, crsublotconfigs);

	cout<<"from sub = "<<subobj<<endl;
	cout<<cutLhs<<endl;
	assert(abs(subobj  -LowerBound)<0.1);

	cout<<"back to master \n";

	add(cutLhs+ datalap::master_instance->ObjConstant <= cutRhs).end();
	cutLhs.end();
	crsublotconfigs.clear();

	


	

	cout<<"HELLO Finished one call to separate()                     \n \n \n                                                         WORLD \n\n"<<"\n\n\n Ta Daaaa";
	
	
	
	return;*/
}

ILOLAZYCONSTRAINTCALLBACK3(BendersLazyCallbacky, BoolVarArray2, y, IloNumVar, Obj, IloNum, ObjConstant)
{
	
	unsigned short TLots = datalap::master_instance->TLots;
	IloEnv masterEnv = getEnv();
	datalap::master_instance->m_perm.resize(TLots);
	datalap::master_instance->precedence.resize(TLots);
	
	for(unsigned short i=0;i<TLots;++i)
	{
		datalap::master_instance->m_perm[i].resize(TLots);
		datalap::master_instance->precedence[i].resize(TLots);

		for(unsigned short j=0;j<TLots;++j)
		{
			

			
			if(i==j)
				datalap::master_instance->precedence[i][j]=0;
			else
				datalap::master_instance->precedence[i][j]= getValue(y[i][j]);
			
			cout<<datalap::master_instance->precedence[i][j]<<" ";
			

		}
		cout<<endl;

	}	

	cout<<endl;
	datalap::master_instance->prec_to_perm();

	for(unsigned short i=0;i<TLots;++i)
	{
		for(unsigned short j=0;j<TLots;++j)
		{
			
			cout<<datalap::master_instance->m_perm[i][j]<< " ";
		} cout<< endl;}cout<<endl;

	cout<<"perm from master:"<<endl;
	datalap::master_instance->printprec();
	//datalap::master_instance->clean_perm(datalap::master_instance->m_perm);
	//datalap::master_instance->perm_to_prec();
	double LowerBound = datalap::master_instance->givensequence();
	cout<<"Critical Lot="<<datalap::master_instance->critical_lot<<endl;
	cout<<"Critical Sublot="<<datalap::master_instance->critical_sublot<<endl;
	cout<<"Critical Vendor="<<datalap::master_instance->critical_vendor<<endl;
	cout<<LowerBound<<" from given seq.\n";

	IloExpr cutLhs(masterEnv);
	IloNumVar cutRhs = Obj;
	cout<<"obj constant in sub="<<datalap::master_instance->ObjConstant<<endl;
	cout<<"Calling Seperate()                    \n \n \n  WORLD \n\n"<<"\n\n\n Ta Daaaa\n";
	double subobj = 0;

	//subobj = datalap::master_instance->separate( cutLhs,
	//	x, y, LowerBound,
	//datalap::master_instance->critical_lot, 
	//datalap::master_instance->critical_sublot-1,
	//datalap::master_instance->critical_vendor);
	//cout<<"from sub = "<<subobj<<endl;
	//cout<<cutLhs<<endl;
	//assert(abs(subobj  -LowerBound)<0.1);

	vector<unsigned short> crsublotconfigs;
	crsublotconfigs.push_back(datalap::master_instance->submap_inv[datalap::master_instance->critical_lot]
	[datalap::master_instance->critical_vendor]
	[datalap::master_instance->critical_sublot-1]);
	
	subobj = datalap::master_instance->separatesmallbooly( cutLhs, y, LowerBound, crsublotconfigs);

	cout<<"from sub = "<<subobj<<endl;
	cout<<cutLhs<<endl;
	assert(abs(subobj  -LowerBound)<0.1);

	cout<<"back to master \n";

	add(cutLhs+ datalap::master_instance->ObjConstant <= cutRhs).end();
	cutLhs.end();
	crsublotconfigs.clear();

	


	

	cout<<"HELLO Finished one call to separate()                     \n \n \n                                                         WORLD \n\n"<<"\n\n\n Ta Daaaa";
	
	/*cout<<NumSublots<<endl;
	IloEnv work_env= workerCplex.getEnv();
	IloNumArray init(work_env,l.getSize());
	for(unsigned short i=0;i<l.getSize();i++)
		init[i]=i;
	workerCplex.setVectors(init,0,l,0,0,0);
	workerCplex.setParam(IloCplex::PreInd, false);
	workerCplex.solve();
	cout<<"Dual="<<workerCplex.getDual(pickone)<<endl;
	for(unsigned short i=0;i<l.getSize();i++)
		cout<<workerCplex.getValue(l[i])<<endl;*/
	
	return;
}

void datalap::prec_to_perm()
{
	vector<unsigned short> position;
	position.resize(TLots);

	for(unsigned short i=0;i<TLots;++i)
	{
		position[i]=0;
		for(unsigned short j=0;j<TLots;++j)
		{
			position[i]+= precedence[i][j];

		}
		//cout<<position[i]<<endl;

	}

	for(unsigned short i=0;i<TLots;++i)
	{
		
		for(unsigned short j=0;j<TLots;++j)
		{
			if(j==position[i])
				m_perm[i][j]=1;
			else
				m_perm[i][j]=0;
			//cout<<precedence[i][j]<<" ";
		}	
		//cout<<endl;

	}
}

double datalap::separate(IloExpr cutLhs, 
	BoolVarArray2 x, FloatVarArray2 y, IloNum GivenObj,
	unsigned short crLot, unsigned short crSublot, unsigned short crVendor)
{
	std::ofstream datafile;
	std::fstream cutfile;
	datafile.open("out.txt", std::fstream::app);
	cutfile.open("cuts.txt", std::fstream::app);

	unsigned short numSublots = submap.size();
	double * ctime = new double[numSublots];
	cout<<endl;
	for(unsigned short i=0;i<TLots;++i)
	{
		
		for(unsigned short j=0;j<TLots;++j)
		{
			
			cout<<precedence[i][j]<<" ";

		}	
		cout<<endl;

	}
	cout<<"critical sublot config="<<submap_inv[crLot][crVendor][crSublot]<<endl;
	cout<<"numsublots="<<numSublots<<endl;
	cout<<"ctime:"<<endl;
	for(unsigned short i=0;i<numSublots ;++i)
	{
		ctime[i] = 0;
		unsigned short crl = submap[i].parentlot;
		unsigned short crv = submap[i].parentvendor;
		unsigned short crs= submap[i].sublotsused;
		
		for(unsigned short kk = 0; kk < TLots; ++kk)
				{
					if(kk != crl)
					ctime[i]+= (lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*precedence[kk][crl];
					
					
				}
		ctime[i]+= lvdata[crl][crv].makespansmp[crs]- lvdata[crl][crv].assemblytime ;
					if(ctime[i]==ctime[submap_inv[crLot][crVendor][crSublot]])
						cout<<i<<endl;
	}
	

	IloEnv workerEnv;
	

	//Declare Variables
	FloatVarArray1 lamda(workerEnv, numSublots,0,1);

	FloatVarArray2 z(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		z[i] = FloatVarArray1(workerEnv, numSublots,0,1);

	FloatVarArray2 ylamda(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		ylamda[i] = FloatVarArray1(workerEnv, TLots,0,1);

	FloatVarArray2 zlamda(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		zlamda[i] = FloatVarArray1(workerEnv, numSublots,0,1);

	FloatVarArray3 zyi(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		zyi[i] = FloatVarArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
	{
		for(unsigned short j=0; j<numSublots; ++j)
		{
			zyi[i][j] = FloatVarArray1(workerEnv, TLots,0,1);
			
		}
	}

	FloatVarArray3 zyj(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		zyj[i] = FloatVarArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
	{
		for(unsigned short j=0; j<numSublots; ++j)
		{
			zyj[i][j] = FloatVarArray1(workerEnv, TLots,0,1);
		}
	}

	// Declare model
	IloModel workerModel(workerEnv);

	//Declare Constraints
	
	IloExpr gammaex(workerEnv);
	IloExpr workerObj(workerEnv);
	
	

	IloRangeArray2 pi1(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
	{
		pi1[i] = IloRangeArray(workerEnv);
		pi1[i].setSize(TLots);
	}

	IloRangeArray2 pi2(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		pi2[i] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 pi3(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		pi3[i] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 w1(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		w1[i] = IloRangeArray(workerEnv, numSublots);

	IloRangeArray2 w2(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		w2[i] = IloRangeArray(workerEnv, numSublots);

	IloRangeArray2 w3(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		w3[i] = IloRangeArray(workerEnv, numSublots);

	IloRangeArray2 chat1(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		chat1[i] = IloRangeArray(workerEnv, numSublots);

	IloRangeArray2 chat2(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		chat2[i] = IloRangeArray(workerEnv, numSublots);

	IloRangeArray3 the1(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		the1[i] = IloRangeArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		for(unsigned short j=0; j<numSublots; ++j)
			the1[i][j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray3 the2(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		the2[i] = IloRangeArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		for(unsigned short j=0; j<numSublots; ++j)
			the2[i][j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray3 the3(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		the3[i] = IloRangeArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		for(unsigned short j=0; j<numSublots; ++j)
			the3[i][j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray3 phi1(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		phi1[i] = IloRangeArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		for(unsigned short j=0; j<numSublots; ++j)
			phi1[i][j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray3 phi2(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		phi2[i] = IloRangeArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		for(unsigned short j=0; j<numSublots; ++j)
			phi2[i][j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray3 phi3(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		phi3[i] = IloRangeArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		for(unsigned short j=0; j<numSublots; ++j)
			phi3[i][j] = IloRangeArray(workerEnv, TLots);
	//Add constraints	
	for(unsigned short i=0;i<numSublots ;++i)
	{
		unsigned short crl = submap[i].parentlot;
		unsigned short crv = submap[i].parentvendor;
		unsigned short crs= submap[i].sublotsused;

		if(crl==crLot && crv==crVendor && crs==crSublot)
			gammaex+= lamda[i];

		//lamda in objective
		workerObj+= (lvdata[crl][crv].makespansmp[crs] 
		-lvdata[crl][crv].assemblytime)*lamda[i];

		//yl in objective
		for(unsigned short lk = 0; lk < TLots; ++ lk)
		{
			if(lk != crl)
			{
				workerObj+= (lvdata[lk][crv].vendortime - lvdata[lk][crv].assemblytime)*ylamda[i][lk];
			}
		}

		//zlamda in objective
		for(unsigned short subk = 0; subk < numSublots; ++ subk)
		{
			unsigned short subk_l = submap[subk].parentlot;
			unsigned short subk_v = submap[subk].parentvendor;
			unsigned short subk_s = submap[subk].sublotsused;

			if(subk != i)
			{
				if(subk_s <= lvdata[subk_l][subk_v].maxsublots - 2)// 2 since -1 is last element index 
					workerObj+= (lvdata[subk_l][subk_v].handlingcost[subk_s +1] -
					lvdata[subk_l][subk_v].handlingcost[subk_s] )*zlamda[i][subk];
				else
					workerObj+= 0*zlamda[i][subk];			

			}
		}
		
	
		//pi
		
	for(unsigned short lk = 0; lk < TLots; ++ lk)
		{			
			if(lk != crl)
			{
				IloExpr temp1(workerEnv);
				IloExpr temp2(workerEnv);
				IloExpr temp3(workerEnv);

				temp1+= lamda[i]-ylamda[i][lk];
				temp2+= -ylamda[i][lk];
				temp3+=-lamda[i]+ylamda[i][lk];
				

				pi1[i][lk] = IloRange(workerEnv,0,temp1,IloInfinity);
				pi2[i][lk] = IloRange(workerEnv,-1*precedence[lk][crl],temp2,IloInfinity);
				pi3[i][lk] = IloRange(workerEnv,-1*precedence[crl][lk],temp3,IloInfinity);

				workerModel.add(pi1[i][lk]);
				workerModel.add(pi2[i][lk]);
				workerModel.add(pi3[i][lk]);

				temp1.end();
				temp2.end();
				temp3.end();

				
			}	

			
		}

	//w
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
		{			
			if(jbu != i)
			{
				

				IloExpr temp1(workerEnv);
				IloExpr temp2(workerEnv);
				IloExpr temp3(workerEnv);

				temp1+= lamda[i]-zlamda[i][jbu];
				temp2+= z[i][jbu]-zlamda[i][jbu];
				if(lvdata[submap[jbu].parentlot][submap[jbu].parentvendor].maxsublots-1 == submap[jbu].sublotsused)
					temp3+=  -lamda[i] -z[i][jbu];
				else
					temp3+= zlamda[i][jbu] -lamda[i] -z[i][jbu];
				

				w1[i][jbu] = IloRange(workerEnv,0,temp1,IloInfinity);
				w2[i][jbu] = IloRange(workerEnv,0,temp2,IloInfinity);
				w3[i][jbu] = IloRange(workerEnv,-1,temp3,IloInfinity);

				workerModel.add(w1[i][jbu]);
				workerModel.add(w2[i][jbu]);
				workerModel.add(w3[i][jbu]);

				temp1.end();
				temp2.end();
				temp3.end();
				
			}	

		}

	

	//chat
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
		{			
			if(jbu != i)
			{
				
				unsigned short jbu_l = submap[jbu].parentlot;
				unsigned short jbu_v = submap[jbu].parentvendor;
				unsigned short jbu_s = submap[jbu].sublotsused;
				

				IloExpr temp1(workerEnv);
				IloNum tempRHS = 0;

				char zname[100];
				sprintf(zname, "z.%d.%d", (int) i, (int) jbu); 
				z[i][jbu].setName(zname);
				
				for(unsigned short kk = 0; kk < TLots; ++kk)
				{	

					if(kk != jbu_l)
					{
						char varNameyj[100];
						sprintf(varNameyj, "zj.%d.%d_y.%d.%d", (int) i, (int) jbu, (int) kk, (int) jbu_l); 
						zyj[i][jbu][kk].setName(varNameyj);

						temp1+= (lvdata[kk][jbu_v].vendortime - lvdata[kk][jbu_v].assemblytime)*zyj[i][jbu][kk];
						tempRHS+= (lvdata[kk][jbu_v].vendortime - lvdata[kk][jbu_v].assemblytime)*precedence[kk][jbu_l];
					}

					if(kk != crl)
					{
						char varNamezyi[100];
						sprintf(varNamezyi, "zi.%d.%d_y.%d.%d", (int) i, (int) jbu, (int) kk, (int) crl); 
						zyi[i][jbu][kk].setName(varNamezyi);

						temp1+= -1*(lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*zyi[i][jbu][kk];
						tempRHS+= -1*(lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*precedence[kk][crl];


					}

					


				}//end kk
				
				
				//cout<<jbu_l<<" "<<jbu_v<<" " <<jbu_s<<endl<<"max subs="<<lvdata[jbu_l][jbu_v].maxsublots<<endl;
				
				temp1+= (lvdata[jbu_l][jbu_v].makespansmp[jbu_s] - lvdata[crl][crv].makespansmp[crs] + 
						-lvdata[jbu_l][jbu_v].assemblytime + lvdata[crl][crv].assemblytime)*z[i][jbu];
				
					tempRHS+= (lvdata[jbu_l][jbu_v].makespansmp[jbu_s] - lvdata[crl][crv].makespansmp[crs] + 
						-lvdata[jbu_l][jbu_v].assemblytime 
						+lvdata[crl][crv].assemblytime);

					

				chat1[i][jbu] = IloRange(workerEnv,tempRHS,temp1,IloInfinity);	
				char varName[100];
				sprintf(varName, "chat1.%d.%d", (int) i, (int) jbu); 
				chat1[i][jbu].setName(varName);

				chat2[i][jbu] = IloRange(workerEnv,0,temp1,IloInfinity);	
				char varName2[100];
				sprintf(varName2, "chat2.%d.%d", (int) i, (int) jbu); 
				chat2[i][jbu].setName(varName2);

				workerModel.add(chat1[i][jbu]);
				workerModel.add(chat2[i][jbu]);

				

				/*for(unsigned short ii=0;ii<TLots;++ii)
				{

					for(unsigned short jjj=0;jjj<TLots;++jjj)
					{

						cout<<precedence[ii][jjj]<<" ";

					}	
					cout<<endl;

				}}*/

				temp1.end();
				
				
			}//endif	

		}//end jbu
	
	//theta
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
	{			
		if(jbu != i)
		{
			for(unsigned short lk = 0; lk < TLots; ++ lk)
			{	
				unsigned short jbu_l = submap[jbu].parentlot;
				if(lk != jbu_l)
				{
					

					IloExpr temp1(workerEnv);
					IloExpr temp2(workerEnv);
					IloExpr temp3(workerEnv);

					temp1+= -zyj[i][jbu][lk];
					temp2+=  -zyj[i][jbu][lk] + z[i][jbu];
					temp3+= zyj[i][jbu][lk] - z[i][jbu];


					the1[i][jbu][lk] = IloRange(workerEnv,-1*precedence[lk][jbu_l],temp1,IloInfinity);
					the2[i][jbu][lk] = IloRange(workerEnv,0,temp2,IloInfinity);
					the3[i][jbu][lk] = IloRange(workerEnv,-1*precedence[jbu_l][lk],temp3,IloInfinity);

					workerModel.add(the1[i][jbu][lk]);
					workerModel.add(the2[i][jbu][lk]);
					workerModel.add(the3[i][jbu][lk]);

					temp1.end();
					temp2.end();
					temp3.end();
				}}

		}}//end theta

	//phi
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
	{			
		if(jbu != i)
		{
			char varName[100];
				sprintf(varName, "z.%d.%d", (int) i, (int) jbu); 
				z[i][jbu].setName(varName);

			for(unsigned short lk = 0; lk < TLots; ++ lk)
			{	
				unsigned short cbu_j = submap[jbu].parentlot;
				if(lk != crl)
				{
					

					IloExpr temp1(workerEnv);
					IloExpr temp2(workerEnv);
					IloExpr temp3(workerEnv);

					temp1+=  -zyi[i][jbu][lk];
					temp2+=  -zyi[i][jbu][lk]  +z[i][jbu];
					temp3+=   zyi[i][jbu][lk] - z[i][jbu];


					phi1[i][jbu][lk] = IloRange(workerEnv,-1*precedence[lk][crl],temp1,IloInfinity);
					phi2[i][jbu][lk] = IloRange(workerEnv,0,temp2,IloInfinity);
					phi3[i][jbu][lk] = IloRange(workerEnv,-1*precedence[crl][lk],temp3,IloInfinity);

					/*if(i == 0 && jbu == 1)
				{
					cout<<"phi 1 "<<phi1[i][jbu][lk].getExpr()<<">="<<phi1[i][jbu][lk].getLB()<<endl;
					cout<<"phi 2 "<<phi2[i][jbu][lk].getExpr()<<">="<<phi2[i][jbu][lk].getLB()<<endl;
					cout<<"phi 3 "<<phi3[i][jbu][lk].getExpr()<<">="<<phi3[i][jbu][lk].getLB()<<endl;
					}*/

					workerModel.add(phi1[i][jbu][lk]);
					workerModel.add(phi2[i][jbu][lk]);
					workerModel.add(phi3[i][jbu][lk]);

					temp1.end();
					temp2.end();
					temp3.end();
				}}

		}}//end phi

	}// end of big subs loop

	


		IloRange gamma(workerEnv,1,gammaex,IloInfinity,"gamma");
		workerModel.add(gamma);
	
		workerModel.add(IloMinimize(workerEnv,workerObj));

	IloCplex cplex(workerModel);
	cplex.setOut(workerEnv.getNullStream());
	cout<<"Attempting worker solution"<<endl;
	if(cplex.solve())
	{
		double y_factor=0;
		cout<<"Solved worker cplex"<<endl;

		cut_yij_coeff.clear();		
		cut_yij_coeff.resize(TLots);		

		for(unsigned short i=0; i<TLots; ++i)
		{
			cut_yij_coeff[i].resize(TLots);
			fill(cut_yij_coeff[i].begin(), cut_yij_coeff[i].end(), 0);
			
		}

		

		unsigned short icr = submap_inv[crLot][crVendor][crSublot];
		
		//add pi
		for(unsigned short k=0; k<TLots; ++k)
		{
			if(k!=crLot)
			{
				double temp= cplex.getDual(pi2[icr][k]);
				cut_yij_coeff[k][crLot]+= -1*temp;

				temp= cplex.getDual(pi3[icr][k]);
				cut_yij_coeff[crLot][k]= -1*temp;
			}
		}

		

		//add chat1
		for(unsigned short jbu=0; jbu<numSublots; ++jbu)
		{
			unsigned short jbu_l = submap[jbu].parentlot;
			unsigned short jbu_v = submap[jbu].parentvendor;

			if(jbu != icr)
			{
				double chat1_jbu_dual = cplex.getDual(chat1[icr][jbu]);
				for(unsigned short k=0; k<TLots; ++k)
				{
					if(jbu_l != k)
					{
						double temp = (lvdata[k][jbu_v].vendortime - lvdata[k][jbu_v].assemblytime);
						cut_yij_coeff[k][jbu_l]+= temp*chat1_jbu_dual;
						if(abs(chat1_jbu_dual)>0)
						{
						cout<<"y co-effs during chat1 first loop \n";
						cout<<"jbu_l="<<jbu_l<<" k="<<k
							<<" temp="<<temp<<" dual="<<chat1_jbu_dual<<endl;
		for(unsigned short i=0; i<TLots; ++i)
		{	
			
			for(unsigned short j=0; j<TLots; ++j)
			cout<<cut_yij_coeff[i][j]<<" ";
			cout<<endl;
		}
		}
					}

					if(k != crLot)
					{
						cut_yij_coeff[k][crLot] += -1*chat1_jbu_dual*(lvdata[k][crVendor].vendortime - lvdata[k][crVendor].assemblytime);
					}
				}
				

			}
		}
		
		
		//add theta and phi
		for(unsigned short jbu=0; jbu<numSublots; ++jbu)
		{
			unsigned short jbu_l = submap[jbu].parentlot;
			unsigned short jbu_v = submap[jbu].parentvendor;

			if(jbu != icr)
			{
				
				for(unsigned short k=0; k<TLots; ++k)
				{
					if(jbu_l != k)
					{
						double temp = -1* cplex.getDual(the1[icr][jbu][k]);
						cut_yij_coeff[k][jbu_l]+= temp;

						temp = -1* cplex.getDual(the3[icr][jbu][k]);
						cut_yij_coeff[jbu_l][k] += temp;
					}

					if(k != crLot)
					{
						double temp = -1* cplex.getDual(phi1[icr][jbu][k]);
						cut_yij_coeff[k][crLot]+= temp;

						temp = -1* cplex.getDual(phi3[icr][jbu][k]);
						cut_yij_coeff[crLot][k] += temp;
					}
				}
				

			}
		}
		
		
		//create cut LHS
		cout<<"y co-effs \n";

		cutfile<< "for sequence: \n";
		for(unsigned short i=0; i<TLots; ++i)
		{
			
			for(unsigned short j=0; j<TLots; ++j)
			{
				
				cutfile<<precedence[i][j]<<" ";
			}
			cutfile<<endl;
		}//add y factors
		cutfile<<"cut added: \n";
		for(unsigned short i=0; i<TLots; ++i)
		{
			
			for(unsigned short j=0; j<TLots; ++j)
			{
				cout<<cut_yij_coeff[i][j]<<" ";
				cutLhs+= cut_yij_coeff[i][j]*y[i][j];
				y_factor+= cut_yij_coeff[i][j]*precedence[i][j];
				if(cut_yij_coeff[i][j] != 0)
					cutfile<<endl<<cut_yij_coeff[i][j]<<"*y["<<i<<"]["<<j<<"] ";
			}
			cout<<endl;
		}//add y factors
		cout<<"y factor="<<y_factor<<endl;
		cout<<"slave objective="<<cplex.getObjValue()<<endl;
		cutLhs+= cplex.getObjValue()-y_factor;
		cutfile<<" + "<<cplex.getObjValue()-y_factor<<" <= z_master \n";
		cout<<cutLhs<<endl;
	}

	//output solution
	cout<<"Objective="<<cplex.getObjValue()<<" + "<<ObjConstant<<endl;
	/*cout<<cplex.getObjValue()+ ObjConstant - m_optimalobjective<<endl;;
	if(abs(cplex.getObjValue()+ ObjConstant - m_optimalobjective) > 0.1)
		objective_check=1;
	else
		objective_check=0;*/

	cout<<"gamma="<<cplex.getDual(gamma)<<endl;
	bool show_out = 1;

	for(unsigned short i = 0; i<numSublots; ++i)
	{
		if(show_out && cplex.getValue(lamda[i])==1)
	{
		cout<<"lamda["<<i<<"]="<<cplex.getValue(lamda[i])<<" ctime= "<<ctime[i]<<" lot="<<submap[i].parentlot<<" "<<"ven="<<submap[i].parentvendor<<" "<<"sub="<<submap		[i].sublotsused<<" "<<endl;
		cout<<"z= \t\t ";
		for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
			if(i!=jbu)
				cout<<cplex.getValue(z[i][jbu])<<" ";
		cout<<endl;

		cout<<"zlam= \t\t ";
		for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
			if(i!=jbu)
				cout<<cplex.getValue(zlamda[i][jbu])<<" ";
				cout<<endl;

		cout<<"ctime diff= \t ";
		for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
			if(i!=jbu)
				cout<<ctime[jbu] - ctime[i]<<" ";
		cout<<endl;

		cout<<"ylamda= ";
		for(unsigned short jbu = 0; jbu < TLots; ++jbu)
			if(submap[i].parentlot != jbu)
				cout<<cplex.getValue(ylamda[i][jbu])<<" ";
		cout<<endl;
	}
	}

	

	
	

	int show_duals=1;
	//cin>>show_duals;
	if(show_duals )
	{
		

		for(unsigned short i = 0; i<numSublots; ++i)
		{
			if( cplex.getValue(lamda[i])>=1)
			{
			datafile<<"lamda["<<i<<"]="<<cplex.getValue(lamda[i])<<" "<<"lot="<<submap[i].parentlot<<" "<<"ven="<<submap[i].parentvendor
				<<"			 "<<"sub="<<submap		[i].sublotsused<<" "<<endl;

			for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
			{			
				if(jbu != i)
				{
					datafile<<"w["<<i<<"]["<<jbu<<"]"<<cplex.getDual(w1[i][jbu])<<" "<<cplex.getDual(w2[i][jbu])<<" "<<cplex.getDual(w3[i][jbu])<<endl;

				}			
			}
			datafile<<endl;
			for(unsigned short lk = 0; lk < TLots; ++lk)
			{			
				if(lk != submap[i].parentlot)
				{
					datafile<<"pi["<<i<<"]["<<lk<<"]"<<cplex.getDual(pi1[i][lk])<<" "<<cplex.getDual(pi2[i][lk])<<" "<<cplex.getDual(pi3[i][lk])<<endl;

				}			
			}

			datafile<<endl;
			for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
			{			
				if(jbu != i)
				{
					datafile<<"chat1["<<i<<"]["<<jbu<<"]"<<cplex.getDual(chat1[i][jbu])<<endl;

				}			
			}
			datafile<<endl;
			datafile<<endl;

			
			for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
			{			
				if(jbu != i)
				{
					datafile<<"chat2["<<i<<"]["<<jbu<<"]"<<cplex.getDual(chat2[i][jbu])<<endl;

				}			
			}
			datafile<<endl;
			datafile<<endl;

			//theta
			for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
			{			
				if(jbu != i)
				{
					for(unsigned short k = 0; k < TLots; ++k)
					{
						if(k!=submap[jbu].parentlot)
							datafile<<"theta["<<i<<"]["<<jbu<<"]"<<cplex.getDual(the1[i][jbu][k])<<" "
							<<cplex.getDual(the2[i][jbu][k])<<" "<<cplex.getDual(the3[i][jbu][k])<<" "<<endl;
					}

				}			
			}
			datafile<<endl;
			datafile<<endl;

			for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
			{			
				if(jbu != i)
				{
					for(unsigned short k = 0; k < TLots; ++k)
					{
						if(k!=submap[i].parentlot)
							datafile<<"phi["<<i<<"]["<<jbu<<"]"<<cplex.getDual(phi1[i][jbu][k])<<" "
							<<cplex.getDual(phi2[i][jbu][k])<<" "<<cplex.getDual(phi3[i][jbu][k])<<" "<<endl;
					}

				}			
			}
			datafile<<endl;
			datafile<<endl;



		}
		}
	}//show duals till this

	cout<<"here\n";
	//workerModel.end();
	double total_obj = cplex.getObjValue()+ObjConstant;
	workerEnv.end();
	cout<<"hear\n";

	datafile.close();
	cutfile.close();
	return total_obj;
}

//solve with all sublots and add duals from all 
//may have some old code
double datalap::separate_allcritical(IloExpr cutLhs, 
	BoolVarArray2 x, BoolVarArray2 y, IloNum GivenObj,
	unsigned short crLot, unsigned short crSublot, unsigned short crVendor)
{
	unsigned short numSublots = submap.size();
	double * ctime = new double[numSublots];
	cout<<endl;
	for(unsigned short i=0;i<TLots;++i)
	{
		
		for(unsigned short j=0;j<TLots;++j)
		{
			
			cout<<precedence[i][j]<<" ";

		}	
		cout<<endl;

	}
	cout<<"critical lot="<<submap_inv[crLot][crVendor][crSublot]<<endl;
	cout<<"numsublots="<<numSublots<<endl;
	cout<<"ctime:"<<endl;
	for(unsigned short i=0;i<numSublots ;++i)
	{
		ctime[i] = 0;
		unsigned short crl = submap[i].parentlot;
		unsigned short crv = submap[i].parentvendor;
		unsigned short crs= submap[i].sublotsused;
		
		for(unsigned short kk = 0; kk < TLots; ++kk)
				{
					if(kk != crl)
					ctime[i]+= (lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*precedence[kk][crl];
					
					
				}
		ctime[i]+= lvdata[crl][crv].makespansmp[crs]- lvdata[crl][crv].assemblytime ;
					if(ctime[i]==ctime[submap_inv[crLot][crVendor][crSublot]])
						cout<<i<<endl;
	}
	

	IloEnv workerEnv;
	

	//Declare Variables
	FloatVarArray1 lamda(workerEnv, numSublots,0,1);

	FloatVarArray2 z(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		z[i] = FloatVarArray1(workerEnv, numSublots,0,1);

	FloatVarArray2 ylamda(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		ylamda[i] = FloatVarArray1(workerEnv, TLots,0,1);

	FloatVarArray2 zlamda(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		zlamda[i] = FloatVarArray1(workerEnv, numSublots,0,1);

	FloatVarArray3 zyi(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		zyi[i] = FloatVarArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
	{
		for(unsigned short j=0; j<numSublots; ++j)
		{
			zyi[i][j] = FloatVarArray1(workerEnv, TLots,0,1);
			
		}
	}

	FloatVarArray3 zyj(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		zyj[i] = FloatVarArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
	{
		for(unsigned short j=0; j<numSublots; ++j)
		{
			zyj[i][j] = FloatVarArray1(workerEnv, TLots,0,1);
		}
	}

	// Declare model
	IloModel workerModel(workerEnv);

	//Declare Constraints
	
	IloExpr gammaex(workerEnv);
	IloExpr workerObj(workerEnv);
	
	

	IloRangeArray2 pi1(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
	{
		pi1[i] = IloRangeArray(workerEnv);
		pi1[i].setSize(TLots);
	}

	IloRangeArray2 pi2(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		pi2[i] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 pi3(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		pi3[i] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 w1(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		w1[i] = IloRangeArray(workerEnv, numSublots);

	IloRangeArray2 w2(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		w2[i] = IloRangeArray(workerEnv, numSublots);

	IloRangeArray2 w3(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		w3[i] = IloRangeArray(workerEnv, numSublots);

	IloRangeArray2 chat1(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		chat1[i] = IloRangeArray(workerEnv, numSublots);

	IloRangeArray2 chat2(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		chat2[i] = IloRangeArray(workerEnv, numSublots);

	IloRangeArray3 the1(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		the1[i] = IloRangeArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		for(unsigned short j=0; j<numSublots; ++j)
			the1[i][j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray3 the2(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		the2[i] = IloRangeArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		for(unsigned short j=0; j<numSublots; ++j)
			the2[i][j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray3 the3(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		the3[i] = IloRangeArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		for(unsigned short j=0; j<numSublots; ++j)
			the3[i][j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray3 phi1(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		phi1[i] = IloRangeArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		for(unsigned short j=0; j<numSublots; ++j)
			phi1[i][j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray3 phi2(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		phi2[i] = IloRangeArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		for(unsigned short j=0; j<numSublots; ++j)
			phi2[i][j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray3 phi3(workerEnv,numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		phi3[i] = IloRangeArray2(workerEnv, numSublots);
	for(unsigned short i=0; i<numSublots; ++i)
		for(unsigned short j=0; j<numSublots; ++j)
			phi3[i][j] = IloRangeArray(workerEnv, TLots);
	//Add constraints	
	for(unsigned short i=0;i<numSublots ;++i)
	{
		unsigned short crl = submap[i].parentlot;
		unsigned short crv = submap[i].parentvendor;
		unsigned short crs= submap[i].sublotsused;

		if(crl==crLot && crv==crVendor && crs==crSublot)
			gammaex+= lamda[i];

		//lamda in objective
		workerObj+= (lvdata[crl][crv].makespansmp[crs] 
		-lvdata[crl][crv].assemblytime)*lamda[i];

		//yl in objective
		for(unsigned short lk = 0; lk < TLots; ++ lk)
		{
			if(lk != crl)
			{
				workerObj+= (lvdata[lk][crv].vendortime - lvdata[lk][crv].assemblytime)*ylamda[i][lk];
			}
		}

		//zlamda in objective
		for(unsigned short subk = 0; subk < numSublots; ++ subk)
		{
			unsigned short subk_l = submap[subk].parentlot;
			unsigned short subk_v = submap[subk].parentvendor;
			unsigned short subk_s = submap[subk].sublotsused;

			if(subk != i)
			{
				if(subk_s <= lvdata[subk_l][subk_v].maxsublots - 2)// 2 since -1 is last element index 
					workerObj+= (lvdata[subk_l][subk_v].handlingcost[subk_s +1] -
					lvdata[subk_l][subk_v].handlingcost[subk_s] )*zlamda[i][subk];
				else
					workerObj+= 0*zlamda[i][subk];			

			}
		}
		
	
		//pi
		
	for(unsigned short lk = 0; lk < TLots; ++ lk)
		{			
			if(lk != crl)
			{
				IloExpr temp1(workerEnv);
				IloExpr temp2(workerEnv);
				IloExpr temp3(workerEnv);

				temp1+= lamda[i]-ylamda[i][lk];
				temp2+= -ylamda[i][lk];
				temp3+=-lamda[i]+ylamda[i][lk];
				

				pi1[i][lk] = IloRange(workerEnv,0,temp1,IloInfinity);
				pi2[i][lk] = IloRange(workerEnv,-1*precedence[lk][crl],temp2,IloInfinity);
				pi3[i][lk] = IloRange(workerEnv,-1*precedence[crl][lk],temp3,IloInfinity);

				workerModel.add(pi1[i][lk]);
				workerModel.add(pi2[i][lk]);
				workerModel.add(pi3[i][lk]);

				temp1.end();
				temp2.end();
				temp3.end();

				
			}	

			
		}

	//w
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
		{			
			if(jbu != i)
			{
				

				IloExpr temp1(workerEnv);
				IloExpr temp2(workerEnv);
				IloExpr temp3(workerEnv);

				temp1+= lamda[i]-zlamda[i][jbu];
				temp2+= z[i][jbu]-zlamda[i][jbu];
				if(lvdata[submap[jbu].parentlot][submap[jbu].parentvendor].maxsublots-1 == submap[jbu].sublotsused)
					temp3+=  -lamda[i] -z[i][jbu];
				else
					temp3+= zlamda[i][jbu] -lamda[i] -z[i][jbu];
				

				w1[i][jbu] = IloRange(workerEnv,0,temp1,IloInfinity);
				w2[i][jbu] = IloRange(workerEnv,0,temp2,IloInfinity);
				w3[i][jbu] = IloRange(workerEnv,-1,temp3,IloInfinity);

				workerModel.add(w1[i][jbu]);
				workerModel.add(w2[i][jbu]);
				workerModel.add(w3[i][jbu]);

				temp1.end();
				temp2.end();
				temp3.end();
				
			}	

		}

	

	//chat
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
		{			
			if(jbu != i)
			{
				
				unsigned short jbu_l = submap[jbu].parentlot;
				unsigned short jbu_v = submap[jbu].parentvendor;
				unsigned short jbu_s = submap[jbu].sublotsused;
				

				IloExpr temp1(workerEnv);
				IloNum tempRHS = 0;

				char zname[100];
				sprintf(zname, "z.%d.%d", (int) i, (int) jbu); 
				z[i][jbu].setName(zname);
				
				for(unsigned short kk = 0; kk < TLots; ++kk)
				{	

					if(kk != jbu_l)
					{
						char varNameyj[100];
						sprintf(varNameyj, "zj.%d.%d_y.%d.%d", (int) i, (int) jbu, (int) kk, (int) jbu_l); 
						zyj[i][jbu][kk].setName(varNameyj);

						temp1+= (lvdata[kk][jbu_v].vendortime - lvdata[kk][jbu_v].assemblytime)*zyj[i][jbu][kk];
						tempRHS+= (lvdata[kk][jbu_v].vendortime - lvdata[kk][jbu_v].assemblytime)*precedence[kk][jbu_l];
					}

					if(kk != crl)
					{
						char varNamezyi[100];
						sprintf(varNamezyi, "zi.%d.%d_y.%d.%d", (int) i, (int) jbu, (int) kk, (int) crl); 
						zyi[i][jbu][kk].setName(varNamezyi);

						temp1+= -1*(lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*zyi[i][jbu][kk];
						tempRHS+= -1*(lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*precedence[kk][crl];


					}

					


				}//end kk
				
				
				//cout<<jbu_l<<" "<<jbu_v<<" " <<jbu_s<<endl<<"max subs="<<lvdata[jbu_l][jbu_v].maxsublots<<endl;
				
				temp1+= (lvdata[jbu_l][jbu_v].makespansmp[jbu_s] - lvdata[crl][crv].makespansmp[crs] + 
						-lvdata[jbu_l][jbu_v].assemblytime + lvdata[crl][crv].assemblytime)*z[i][jbu];
				
					tempRHS+= (lvdata[jbu_l][jbu_v].makespansmp[jbu_s] - lvdata[crl][crv].makespansmp[crs] + 
						-lvdata[jbu_l][jbu_v].assemblytime 
						+lvdata[crl][crv].assemblytime);

					

				chat1[i][jbu] = IloRange(workerEnv,tempRHS,temp1,IloInfinity);	
				char varName[100];
				sprintf(varName, "chat1.%d.%d", (int) i, (int) jbu); 
				chat1[i][jbu].setName(varName);

				chat2[i][jbu] = IloRange(workerEnv,0,temp1,IloInfinity);	
				char varName2[100];
				sprintf(varName2, "chat2.%d.%d", (int) i, (int) jbu); 
				chat2[i][jbu].setName(varName2);

				workerModel.add(chat1[i][jbu]);
				workerModel.add(chat2[i][jbu]);

				

				/*for(unsigned short ii=0;ii<TLots;++ii)
				{

					for(unsigned short jjj=0;jjj<TLots;++jjj)
					{

						cout<<precedence[ii][jjj]<<" ";

					}	
					cout<<endl;

				}}*/

				temp1.end();
				
				
			}//endif	

		}//end jbu
	
	//theta
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
	{			
		if(jbu != i)
		{
			for(unsigned short lk = 0; lk < TLots; ++ lk)
			{	
				unsigned short jbu_l = submap[jbu].parentlot;
				if(lk != jbu_l)
				{
					

					IloExpr temp1(workerEnv);
					IloExpr temp2(workerEnv);
					IloExpr temp3(workerEnv);

					temp1+= -zyj[i][jbu][lk];
					temp2+=  -zyj[i][jbu][lk] + z[i][jbu];
					temp3+= zyj[i][jbu][lk] - z[i][jbu];


					the1[i][jbu][lk] = IloRange(workerEnv,-1*precedence[lk][jbu_l],temp1,IloInfinity);
					the2[i][jbu][lk] = IloRange(workerEnv,0,temp2,IloInfinity);
					the3[i][jbu][lk] = IloRange(workerEnv,-1*precedence[jbu_l][lk],temp3,IloInfinity);

					workerModel.add(the1[i][jbu][lk]);
					workerModel.add(the2[i][jbu][lk]);
					workerModel.add(the3[i][jbu][lk]);

					temp1.end();
					temp2.end();
					temp3.end();
				}}

		}}//end theta

	//phi
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
	{			
		if(jbu != i)
		{
			char varName[100];
				sprintf(varName, "z.%d.%d", (int) i, (int) jbu); 
				z[i][jbu].setName(varName);

			for(unsigned short lk = 0; lk < TLots; ++ lk)
			{	
				unsigned short cbu_j = submap[jbu].parentlot;
				if(lk != crl)
				{
					

					IloExpr temp1(workerEnv);
					IloExpr temp2(workerEnv);
					IloExpr temp3(workerEnv);

					temp1+=  -zyi[i][jbu][lk];
					temp2+=  -zyi[i][jbu][lk]  +z[i][jbu];
					temp3+=   zyi[i][jbu][lk] - z[i][jbu];


					phi1[i][jbu][lk] = IloRange(workerEnv,-1*precedence[lk][crl],temp1,IloInfinity);
					phi2[i][jbu][lk] = IloRange(workerEnv,0,temp2,IloInfinity);
					phi3[i][jbu][lk] = IloRange(workerEnv,-1*precedence[crl][lk],temp3,IloInfinity);

					/*if(i == 0 && jbu == 1)
				{
					cout<<"phi 1 "<<phi1[i][jbu][lk].getExpr()<<">="<<phi1[i][jbu][lk].getLB()<<endl;
					cout<<"phi 2 "<<phi2[i][jbu][lk].getExpr()<<">="<<phi2[i][jbu][lk].getLB()<<endl;
					cout<<"phi 3 "<<phi3[i][jbu][lk].getExpr()<<">="<<phi3[i][jbu][lk].getLB()<<endl;
					}*/

					workerModel.add(phi1[i][jbu][lk]);
					workerModel.add(phi2[i][jbu][lk]);
					workerModel.add(phi3[i][jbu][lk]);

					temp1.end();
					temp2.end();
					temp3.end();
				}}

		}}//end phi

	}// end of big subs loop

	


		IloRange gamma(workerEnv,1,gammaex,IloInfinity,"gamma");
		workerModel.add(gamma);
	
		workerModel.add(IloMinimize(workerEnv,workerObj));

	IloCplex cplex(workerModel);
	cplex.setOut(workerEnv.getNullStream());
	cout<<"Attempting worker solution"<<endl;
	if(cplex.solve())
	{
		double y_factor=0;
		cout<<"Solved worker cplex"<<endl;

		cut_yij_coeff.clear();		
		cut_yij_coeff.resize(TLots);		

		for(unsigned short i=0; i<TLots; ++i)
		{
			cut_yij_coeff[i].resize(TLots);
			fill(cut_yij_coeff[i].begin(), cut_yij_coeff[i].end(), 0);
			
		}

		
		// submap_inv[crLot][crVendor][crSublot] earlier tried to use only this as icr
		for(unsigned short icr = 0; icr < numSublots; ++icr)
		{
			crLot = submap[icr].parentlot;
			crVendor = submap[icr].parentvendor;
			crSublot = submap[icr].sublotsused;
		//add pi
		for(unsigned short k=0; k<TLots; ++k)
		{
			if(k!=crLot)
			{
				double temp= cplex.getDual(pi2[icr][k]);
				cut_yij_coeff[k][crLot]+= -1*temp;

				temp= cplex.getDual(pi3[icr][k]);
				cut_yij_coeff[crLot][k]= -1*temp;
			}
		}

		

		//add chat1
		for(unsigned short jbu=0; jbu<numSublots; ++jbu)
		{
			unsigned short jbu_l = submap[jbu].parentlot;
			unsigned short jbu_v = submap[jbu].parentvendor;

			if(jbu != icr)
			{
				double chat1_jbu_dual = cplex.getDual(chat1[icr][jbu]);
				for(unsigned short k=0; k<TLots; ++k)
				{
					if(jbu_l != k)
					{
						double temp = (lvdata[k][jbu_v].vendortime - lvdata[k][jbu_v].assemblytime);
						cut_yij_coeff[k][jbu_l]+= temp*chat1_jbu_dual;
						
					}

					if(k != crLot)
					{
						cut_yij_coeff[k][crLot] += -1*chat1_jbu_dual*(lvdata[k][crVendor].vendortime - lvdata[k][crVendor].assemblytime);
					}
				}
				

			}
		}
		
		
		//add theta and phi
		for(unsigned short jbu=0; jbu<numSublots; ++jbu)
		{
			unsigned short jbu_l = submap[jbu].parentlot;
			unsigned short jbu_v = submap[jbu].parentvendor;

			if(jbu != icr)
			{
				
				for(unsigned short k=0; k<TLots; ++k)
				{
					if(jbu_l != k)
					{
						double temp = -1* cplex.getDual(the1[icr][jbu][k]);
						cut_yij_coeff[k][jbu_l]+= temp;

						temp = -1* cplex.getDual(the3[icr][jbu][k]);
						cut_yij_coeff[jbu_l][k] += temp;
					}

					if(k != crLot)
					{
						double temp = -1* cplex.getDual(phi1[icr][jbu][k]);
						cut_yij_coeff[k][crLot]+= temp;

						temp = -1* cplex.getDual(phi3[icr][jbu][k]);
						cut_yij_coeff[crLot][k] += temp;
					}
				}
				

			}
		}
		

		}//end of icr
		
		//create cut LHS
		cout<<"y co-effs \n";
		for(unsigned short i=0; i<TLots; ++i)
		{
			
			for(unsigned short j=0; j<TLots; ++j)
			{
				cout<<cut_yij_coeff[i][j]<<" ";
				cutLhs+= cut_yij_coeff[i][j]*y[i][j];
				y_factor+= cut_yij_coeff[i][j]*precedence[i][j];
			}
			cout<<endl;
		}//add y factors
		cout<<"y factor="<<y_factor<<endl;
		cout<<"slave objective="<<cplex.getObjValue()<<endl;
		cutLhs+= cplex.getObjValue()-y_factor;
		cout<<cutLhs<<endl;
	}

	//output solution
	cout<<"Objective="<<cplex.getObjValue()<<" + "<<ObjConstant<<endl;
	/*cout<<cplex.getObjValue()+ ObjConstant - m_optimalobjective<<endl;;
	if(abs(cplex.getObjValue()+ ObjConstant - m_optimalobjective) > 0.1)
		objective_check=1;
	else
		objective_check=0;*/

	cout<<"gamma="<<cplex.getDual(gamma)<<endl;
	bool show_out = 0;

	for(unsigned short i = 0; i<numSublots; ++i)
	{
		if(show_out && cplex.getValue(lamda[i])==1)
	{
		cout<<"lamda["<<i<<"]="<<cplex.getValue(lamda[i])<<" ctime= "<<ctime[i]<<" lot="<<submap[i].parentlot<<" "<<"ven="<<submap[i].parentvendor<<" "<<"sub="<<submap		[i].sublotsused<<" "<<endl;
		cout<<"z= \t\t ";
		for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
			if(i!=jbu)
				cout<<cplex.getValue(z[i][jbu])<<" ";
		cout<<endl;

		cout<<"ctime diff= \t ";
		for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
			if(i!=jbu)
				cout<<ctime[jbu] - ctime[i]<<" ";
		cout<<endl;

		cout<<"ylamda= ";
		for(unsigned short jbu = 0; jbu < TLots; ++jbu)
			if(submap[i].parentlot != jbu)
				cout<<cplex.getValue(ylamda[i][jbu])<<" ";
		cout<<endl;
	}
	}

	

	
	

	int show_duals=0;
	//cin>>show_duals;
	if(show_duals )
	{
		

		for(unsigned short i = 0; i<numSublots; ++i)
		{
			if( cplex.getValue(lamda[i])>=1)
			{
			cout<<"lamda["<<i<<"]="<<cplex.getValue(lamda[i])<<" "<<"lot="<<submap[i].parentlot<<" "<<"ven="<<submap[i].parentvendor
				<<"			 "<<"sub="<<submap		[i].sublotsused<<" "<<endl;

			for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
			{			
				if(jbu != i)
				{
					cout<<"w["<<i<<"]["<<jbu<<"]"<<cplex.getDual(w1[i][jbu])<<" "<<cplex.getDual(w2[i][jbu])<<" "<<cplex.getDual(w3[i][jbu])<<endl;

				}			
			}
			cout<<endl;
			for(unsigned short lk = 0; lk < TLots; ++lk)
			{			
				if(lk != submap[i].parentlot)
				{
					cout<<"pi["<<i<<"]["<<lk<<"]"<<cplex.getDual(pi1[i][lk])<<" "<<cplex.getDual(pi2[i][lk])<<" "<<cplex.getDual(pi3[i][lk])<<endl;

				}			
			}

			cout<<endl;
			for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
			{			
				if(jbu != i)
				{
					cout<<"chat1["<<i<<"]["<<jbu<<"]"<<cplex.getDual(chat1[i][jbu])<<endl;

				}			
			}
			cout<<endl;
			cout<<endl;

			
			for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
			{			
				if(jbu != i)
				{
					cout<<"chat2["<<i<<"]["<<jbu<<"]"<<cplex.getDual(chat2[i][jbu])<<endl;

				}			
			}
			cout<<endl;
			cout<<endl;

			//theta
			for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
			{			
				if(jbu != i)
				{
					for(unsigned short k = 0; k < TLots; ++k)
					{
						if(k!=submap[jbu].parentlot)
							cout<<"theta["<<i<<"]["<<jbu<<"]"<<cplex.getDual(the1[i][jbu][k])<<" "
							<<cplex.getDual(the2[i][jbu][k])<<" "<<cplex.getDual(the3[i][jbu][k])<<" "<<endl;
					}

				}			
			}
			cout<<endl;
			cout<<endl;

			for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
			{			
				if(jbu != i)
				{
					for(unsigned short k = 0; k < TLots; ++k)
					{
						if(k!=submap[i].parentlot)
							cout<<"phi["<<i<<"]["<<jbu<<"]"<<cplex.getDual(phi1[i][jbu][k])<<" "
							<<cplex.getDual(phi2[i][jbu][k])<<" "<<cplex.getDual(phi3[i][jbu][k])<<" "<<endl;
					}

				}			
			}
			cout<<endl;
			cout<<endl;



		}
		}
	}//show duals till this

	double total_obj_i = cplex.getObjValue() +  ObjConstant;
	cout<<"floor "<<floor( cplex.getObjValue())<<endl;
	cout<<"ceil "<<ceil( cplex.getObjValue())<<endl;
	cout<<"default "<< cplex.getObjValue() <<endl;
	
	cout<<"Returning from sub = "<<total_obj_i<<" as"<<cplex.getObjValue()<<" + "<<ObjConstant<<endl;
	
	workerEnv.end();
	return total_obj_i;
}

//solve only critical and add only the duals from the critical lot
double datalap::separatesmall(IloExpr cutLhs, 
	BoolVarArray2 x, FloatVarArray2 y, IloNum GivenObj,
	vector<unsigned short> crsublotconfigs)
{
	std::ofstream datafile;
	std::fstream cutfile;
	datafile.open("out.txt", std::fstream::app);
	cutfile.open("cuts.txt", std::fstream::app);
	cout<<"Entered separatesmall"<<endl;
	unsigned short numSublots = submap.size();
	unsigned short numCrSublots = crsublotconfigs.size();

	unsigned short crLot = submap[crsublotconfigs[0]].parentlot;
	unsigned short crVendor = submap[crsublotconfigs[0]].parentvendor;
	unsigned short crSublot = submap[crsublotconfigs[0]].sublotsused;

	double * ctime = new double[numSublots];
	cout<<endl;
	for(unsigned short i=0;i<TLots;++i)
	{
		
		for(unsigned short j=0;j<TLots;++j)
		{
			
			cout<<precedence[i][j]<<" ";

		}	
		cout<<endl;

	}
	cout<<"critical sublot config="<<submap_inv[crLot][crVendor][crSublot]<<endl;
	cout<<"numsublots="<<numSublots<<endl;
	cout<<"ctime:"<<endl;
	for(unsigned short i=0;i<numSublots ;++i)
	{
		ctime[i] = 0;
		unsigned short crl = submap[i].parentlot;
		unsigned short crv = submap[i].parentvendor;
		unsigned short crs= submap[i].sublotsused;
		
		for(unsigned short kk = 0; kk < TLots; ++kk)
				{
					if(kk != crl)
					ctime[i]+= (lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*precedence[kk][crl];
					
					
				}
		ctime[i]+= lvdata[crl][crv].makespansmp[crs]- lvdata[crl][crv].assemblytime ;
					if(ctime[i]==ctime[submap_inv[crLot][crVendor][crSublot]])
						cout<<i<<endl;
	}
	

	IloEnv workerEnv;
	

	//Declare Variables
	IloNumVar lamda(workerEnv,0,1);
	
	
	FloatVarArray1	z(workerEnv, numSublots,0,1);

	
	FloatVarArray1	ylamda(workerEnv, TLots,0,1);

	FloatVarArray1 zlamda(workerEnv, numSublots,0,1);
	 
	
	FloatVarArray2	zyi(workerEnv, numSublots);	
	for(unsigned short j=0; j<numSublots; ++j)
		{
			
			zyi[j] = FloatVarArray1(workerEnv, TLots,0,1);
			
		}
	

	
	FloatVarArray2	zyj(workerEnv, numSublots);
	
		for(unsigned short j=0; j<numSublots; ++j)
		{
			zyj[j] = FloatVarArray1(workerEnv, TLots,0,1);
		}
	
	cout<<"Declared varaibles"<<endl;
	//create values for setvector
	IloNumVarArray vars(workerEnv);
	IloNumArray vals(workerEnv);
	{
		unsigned short i = crsublotconfigs[0];
		

		vars.add(lamda);
		vals.add(1);
		for(unsigned short j=0; j<numSublots; ++j)
		{
			if(i!=j)
			{
				vars.add(z[j]);
				if(ctime[j] > ctime[i])
					vals.add(1);
				else
					vals.add(0);
					
			}
		}

		for(unsigned short j=0; j<numSublots; ++j)
		{
			if(i!=j)
			{
				vars.add(zlamda[j]);
				if(ctime[j] > ctime[submap_inv[crLot][crVendor][crSublot]])
					vals.add(1);
				else
					vals.add(0);
					
			}
		}

		for(unsigned short j=0; j<TLots; ++j)
		{
			unsigned short crlot_temp = submap[i].parentlot;
			if(crlot_temp !=j)
			{
				vars.add(ylamda[j]);
				vals.add(precedence[j][crlot_temp]);
					
			}
		}

		for(unsigned short j=0; j<numSublots; ++j)
		{
			if(i!=j)
			{
				double tempz;
				if(ctime[j] > ctime[submap_inv[crLot][crVendor][crSublot]])
					tempz =1;
				else
					tempz=0;

				for(unsigned short k=0; k<TLots; ++k)
				{
					unsigned short crlot_temp = submap[i].parentlot;
					if(crlot_temp !=k)
					{
						vars.add(zyi[j][k]);
						vals.add(tempz*precedence[k][crlot_temp]);
					}

					unsigned short crlotj_temp = submap[j].parentlot;
					if(crlot_temp !=k)
					{
						vars.add(zyj[j][k]);
						vals.add(tempz*precedence[k][crlotj_temp]);
					}
				}
			}
		}



	}//end of set vector 

	cout<<"Set initial vector values"<<endl;
	
	// Declare model
	IloModel workerModel(workerEnv);

	//Declare Constraints
	
	IloExpr gammaex(workerEnv);
	IloExpr workerObj(workerEnv);
	
	cout<<"Attempting to create constraints"<<endl;

	IloRangeArray pi1(workerEnv,numSublots);
	

	IloRangeArray pi2(workerEnv,numSublots);
	

	IloRangeArray pi3(workerEnv,numSublots);
	

	IloRangeArray w1(workerEnv,numSublots);
	

	IloRangeArray w2(workerEnv,numSublots);
	

	IloRangeArray w3(workerEnv,numSublots);
	

	IloRangeArray chat1(workerEnv,numSublots);
	

	IloRangeArray chat2(workerEnv,numSublots);
	

	IloRangeArray2 the1(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
		the1[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 the2(workerEnv,numSublots);
		for(unsigned short j=0; j<numSublots; ++j)
			the2[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 the3(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			the3[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 phi1(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			phi1[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 phi2(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			phi2[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 phi3(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			phi3[j] = IloRangeArray(workerEnv, TLots);

	//Add constraints	
	//for(unsigned short i_num=0;i_num<numCrSublots ;++i_num)
	{
		unsigned short i = crsublotconfigs[0];
		unsigned short crl = submap[i].parentlot;
		unsigned short crv = submap[i].parentvendor;
		unsigned short crs= submap[i].sublotsused;

		if(crl==crLot && crv==crVendor && crs==crSublot)
			gammaex+= lamda;

		//lamda in objective
		workerObj+= (lvdata[crl][crv].makespansmp[crs] 
		-lvdata[crl][crv].assemblytime)*lamda;

		//yl in objective
		for(unsigned short lk = 0; lk < TLots; ++ lk)
		{
			if(lk != crl)
			{
				workerObj+= (lvdata[lk][crv].vendortime - lvdata[lk][crv].assemblytime)*ylamda[lk];
			}
		}

		//zlamda in objective
		for(unsigned short subk = 0; subk < numSublots; ++ subk)
		{
			unsigned short subk_l = submap[subk].parentlot;
			unsigned short subk_v = submap[subk].parentvendor;
			unsigned short subk_s = submap[subk].sublotsused;

			if(subk != i)
			{
				if(subk_s <= lvdata[subk_l][subk_v].maxsublots - 2)// 2 since -1 is last element index 
					workerObj+= (lvdata[subk_l][subk_v].handlingcost[subk_s +1] -
					lvdata[subk_l][subk_v].handlingcost[subk_s] )*zlamda[subk];
				else
					workerObj+= 1*zlamda[subk];			

			}
		}
		
	
		//pi
		
	for(unsigned short lk = 0; lk < TLots; ++ lk)
		{			
			if(lk != crl)
			{
				IloExpr temp1(workerEnv);
				IloExpr temp2(workerEnv);
				IloExpr temp3(workerEnv);

				temp1+= lamda-ylamda[lk];
				temp2+= -ylamda[lk];
				temp3+=-lamda+ylamda[lk];
				

				pi1[lk] = IloRange(workerEnv,0,temp1,IloInfinity);
				pi2[lk] = IloRange(workerEnv,-1*precedence[lk][crl],temp2,IloInfinity);
				pi3[lk] = IloRange(workerEnv,-1*precedence[crl][lk],temp3,IloInfinity);

				workerModel.add(pi1[lk]);
				workerModel.add(pi2[lk]);
				workerModel.add(pi3[lk]);

				temp1.end();
				temp2.end();
				temp3.end();

				
			}	

			
		}

	//w
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
		{			
			if(jbu != i)
			{
				

				IloExpr temp1(workerEnv);
				IloExpr temp2(workerEnv);
				IloExpr temp3(workerEnv);

				temp1+= lamda-zlamda[jbu];
				temp2+= z[jbu]-zlamda[jbu];
				if(lvdata[submap[jbu].parentlot][submap[jbu].parentvendor].maxsublots-1 == submap[jbu].sublotsused)
					temp3+=  -lamda -z[jbu];
				else
					temp3+= zlamda[jbu] -lamda -z[jbu];
				

				w1[jbu] = IloRange(workerEnv,0,temp1,IloInfinity);
				w2[jbu] = IloRange(workerEnv,0,temp2,IloInfinity);
				w3[jbu] = IloRange(workerEnv,-1,temp3,IloInfinity);

				workerModel.add(w1[jbu]);
				workerModel.add(w2[jbu]);
				workerModel.add(w3[jbu]);

				temp1.end();
				temp2.end();
				temp3.end();
				
			}	

		}

	

	//chat
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
		{			
			if(jbu != i)
			{
				
				unsigned short jbu_l = submap[jbu].parentlot;
				unsigned short jbu_v = submap[jbu].parentvendor;
				unsigned short jbu_s = submap[jbu].sublotsused;
				

				IloExpr temp1(workerEnv);
				IloNum tempRHS = 0;

				char zname[100];
				sprintf(zname, "z.%d.%d", (int) i, (int) jbu); 
				z[jbu].setName(zname);
				
				for(unsigned short kk = 0; kk < TLots; ++kk)
				{	

					if(kk != jbu_l)
					{
						char varNameyj[100];
						sprintf(varNameyj, "zj.%d.%d_y.%d.%d", (int) i, (int) jbu, (int) kk, (int) jbu_l); 
						zyj[jbu][kk].setName(varNameyj);

						temp1+= (lvdata[kk][jbu_v].vendortime - lvdata[kk][jbu_v].assemblytime)*zyj[jbu][kk];
						tempRHS+= (lvdata[kk][jbu_v].vendortime - lvdata[kk][jbu_v].assemblytime)*precedence[kk][jbu_l];
					}

					if(kk != crl)
					{
						char varNamezyi[100];
						sprintf(varNamezyi, "zi.%d.%d_y.%d.%d", (int) i, (int) jbu, (int) kk, (int) crl); 
						zyi[jbu][kk].setName(varNamezyi);

						temp1+= -1*(lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*zyi[jbu][kk];
						tempRHS+= -1*(lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*precedence[kk][crl];


					}

					


				}//end kk
				
				
				//cout<<jbu_l<<" "<<jbu_v<<" " <<jbu_s<<endl<<"max subs="<<lvdata[jbu_l][jbu_v].maxsublots<<endl;
				
				temp1+= (lvdata[jbu_l][jbu_v].makespansmp[jbu_s] - lvdata[crl][crv].makespansmp[crs] + 
						-lvdata[jbu_l][jbu_v].assemblytime + lvdata[crl][crv].assemblytime)*z[jbu];
				
					tempRHS+= (lvdata[jbu_l][jbu_v].makespansmp[jbu_s] - lvdata[crl][crv].makespansmp[crs] + 
						-lvdata[jbu_l][jbu_v].assemblytime 
						+lvdata[crl][crv].assemblytime);

					

				chat1[jbu] = IloRange(workerEnv,tempRHS,temp1,IloInfinity);	
				char varName[100];
				sprintf(varName, "chat1.%d.%d", (int) i, (int) jbu); 
				chat1[jbu].setName(varName);

				chat2[jbu] = IloRange(workerEnv,0,temp1,IloInfinity);	
				char varName2[100];
				sprintf(varName2, "chat2.%d.%d", (int) i, (int) jbu); 
				chat2[jbu].setName(varName2);

				workerModel.add(chat1[jbu]);
				workerModel.add(chat2[jbu]);

				

				temp1.end();
				
				
			}//endif	

		}//end jbu
	
	//theta
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
	{			
		if(jbu != i)
		{
			for(unsigned short lk = 0; lk < TLots; ++ lk)
			{	
				unsigned short jbu_l = submap[jbu].parentlot;
				if(lk != jbu_l)
				{
					

					IloExpr temp1(workerEnv);
					IloExpr temp2(workerEnv);
					IloExpr temp3(workerEnv);

					temp1+= -zyj[jbu][lk];
					temp2+=  -zyj[jbu][lk] + z[jbu];
					temp3+= zyj[jbu][lk] - z[jbu];


					the1[jbu][lk] = IloRange(workerEnv,-1*precedence[lk][jbu_l],temp1,IloInfinity);
					the2[jbu][lk] = IloRange(workerEnv,0,temp2,IloInfinity);
					the3[jbu][lk] = IloRange(workerEnv,-1*precedence[jbu_l][lk],temp3,IloInfinity);

					workerModel.add(the1[jbu][lk]);
					workerModel.add(the2[jbu][lk]);
					workerModel.add(the3[jbu][lk]);

					temp1.end();
					temp2.end();
					temp3.end();
				}}

		}}//end theta

	//phi
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
	{			
		if(jbu != i)
		{
			char varName[100];
				sprintf(varName, "z.%d.%d", (int) i, (int) jbu); 
				z[jbu].setName(varName);

			for(unsigned short lk = 0; lk < TLots; ++ lk)
			{	
				unsigned short cbu_j = submap[jbu].parentlot;
				if(lk != crl)
				{
					

					IloExpr temp1(workerEnv);
					IloExpr temp2(workerEnv);
					IloExpr temp3(workerEnv);

					temp1+=  -zyi[jbu][lk];
					temp2+=  -zyi[jbu][lk]  +z[jbu];
					temp3+=   zyi[jbu][lk] - z[jbu];


					phi1[jbu][lk] = IloRange(workerEnv,-1*precedence[lk][crl],temp1,IloInfinity);
					phi2[jbu][lk] = IloRange(workerEnv,0,temp2,IloInfinity);
					phi3[jbu][lk] = IloRange(workerEnv,-1*precedence[crl][lk],temp3,IloInfinity);

					/*if(i == 0 && jbu == 1)
				{
					cout<<"phi 1 "<<phi1[i][jbu][lk].getExpr()<<">="<<phi1[i][jbu][lk].getLB()<<endl;
					cout<<"phi 2 "<<phi2[i][jbu][lk].getExpr()<<">="<<phi2[i][jbu][lk].getLB()<<endl;
					cout<<"phi 3 "<<phi3[i][jbu][lk].getExpr()<<">="<<phi3[i][jbu][lk].getLB()<<endl;
					}*/

					workerModel.add(phi1[jbu][lk]);
					workerModel.add(phi2[jbu][lk]);
					workerModel.add(phi3[jbu][lk]);

					temp1.end();
					temp2.end();
					temp3.end();
				}}

		}}//end phi

	}// end of big subs loop

	cout<<"Added all constraints"<<endl;


	IloRange gamma(workerEnv,1,gammaex,IloInfinity,"gamma");
	workerModel.add(gamma);
	
	workerModel.add(IloMinimize(workerEnv,workerObj));

	IloCplex cplex(workerModel);
	cplex.setParam(IloCplex::PreDual,1);
	cplex.setVectors(vals, 0, vars,0,0,0);
	cplex.setParam(IloCplex::PreInd, false);
	cplex.setParam(IloCplex::WorkMem, MaxMemory);
	//cplex.setOut(workerEnv.getNullStream());

	cout<<"Attempting worker solution"<<endl;
	if(cplex.solve())
	{
		double y_factor=0;
		cout<<"Solved worker cplex"<<endl;

		cut_yij_coeff.clear();		
		cut_yij_coeff.resize(TLots);		

		for(unsigned short i=0; i<TLots; ++i)
		{
			cut_yij_coeff[i].resize(TLots);
			fill(cut_yij_coeff[i].begin(), cut_yij_coeff[i].end(), 0);
			
		}

		

		unsigned short icr = submap_inv[crLot][crVendor][crSublot];
		cout<<"trying to add pi duals"<<endl;
		//add pi
		for(unsigned short k=0; k<TLots; ++k)
		{
			if(k!=crLot)
			{
				double temp= cplex.getDual(pi2[k]);
				cut_yij_coeff[k][crLot]+= -1*temp;

				temp= cplex.getDual(pi3[k]);
				cut_yij_coeff[crLot][k]= -1*temp;
			}
		}

		
		cout<<"trying to add chat duals"<<endl;
		//add chat1
		for(unsigned short jbu=0; jbu<numSublots; ++jbu)
		{
			unsigned short jbu_l = submap[jbu].parentlot;
			unsigned short jbu_v = submap[jbu].parentvendor;

			if(jbu != icr)
			{
				double chat1_jbu_dual = cplex.getDual(chat1[jbu]);
				for(unsigned short k=0; k<TLots; ++k)
				{
					if(jbu_l != k)
					{
						double temp = (lvdata[k][jbu_v].vendortime - lvdata[k][jbu_v].assemblytime);
						cut_yij_coeff[k][jbu_l]+= temp*chat1_jbu_dual;
						if(abs(chat1_jbu_dual)>0)
						{
						cout<<"y co-effs during chat1 first loop \n";
						cout<<"jbu_l="<<jbu_l<<" k="<<k
							<<" temp="<<temp<<" dual="<<chat1_jbu_dual<<endl;
		for(unsigned short i=0; i<TLots; ++i)
		{	
			
			for(unsigned short j=0; j<TLots; ++j)
			cout<<cut_yij_coeff[i][j]<<" ";
			cout<<endl;
		}
		}
					}

					if(k != crLot)
					{
						cut_yij_coeff[k][crLot] += -1*chat1_jbu_dual*(lvdata[k][crVendor].vendortime - lvdata[k][crVendor].assemblytime);
					}
				}
				

			}
		}
		
		cout<<"trying to add phi and theta duals"<<endl;
		//add theta and phi
		for(unsigned short jbu=0; jbu<numSublots; ++jbu)
		{
			unsigned short jbu_l = submap[jbu].parentlot;
			unsigned short jbu_v = submap[jbu].parentvendor;

			if(jbu != icr)
			{
				
				for(unsigned short k=0; k<TLots; ++k)
				{
					if(jbu_l != k)
					{
						double temp = -1* cplex.getDual(the1[jbu][k]);
						cut_yij_coeff[k][jbu_l]+= temp;

						temp = -1* cplex.getDual(the3[jbu][k]);
						cut_yij_coeff[jbu_l][k] += temp;
					}

					if(k != crLot)
					{
						double temp = -1* cplex.getDual(phi1[jbu][k]);
						cut_yij_coeff[k][crLot]+= temp;

						temp = -1* cplex.getDual(phi3[jbu][k]);
						cut_yij_coeff[crLot][k] += temp;
					}
				}
				

			}
		}
		
		
		//create cut LHS
		cout<<"y co-effs \n";

		cutfile<< "for sequence: \n";
		for(unsigned short i=0; i<TLots; ++i)
		{
			
			for(unsigned short j=0; j<TLots; ++j)
			{
				
				cutfile<<precedence[i][j]<<" ";
			}
			cutfile<<endl;
		}//add y factors
		cutfile<<"cut added: \n";
		for(unsigned short i=0; i<TLots; ++i)
		{
			
			for(unsigned short j=0; j<TLots; ++j)
			{
				cout<<cut_yij_coeff[i][j]<<" ";
				cutLhs+= cut_yij_coeff[i][j]*y[i][j];
				y_factor+= cut_yij_coeff[i][j]*precedence[i][j];
				if(cut_yij_coeff[i][j] != 0)
					cutfile<<endl<<cut_yij_coeff[i][j]<<"*y["<<i<<"]["<<j<<"] ";
			}
			cout<<endl;
		}//add y factors
		cout<<"y factor="<<y_factor<<endl;
		cout<<"slave objective="<<cplex.getObjValue()<<endl;
		cutLhs+= cplex.getObjValue()-y_factor;
		cutfile<<" + "<<cplex.getObjValue()-y_factor<<" <= z_master \n";
		cout<<cutLhs<<endl;
	}

	

	cout<<"here\n";
	//workerModel.end();
	double total_obj = cplex.getObjValue()+ObjConstant;
	workerEnv.end();
	cout<<"hear\n";

	datafile.close();
	cutfile.close();
	return total_obj;
}

//solve only critical and add only the duals from the critical lot
double datalap::separatesmallbooly(IloExpr cutLhs, 
	 BoolVarArray2 y, IloNum GivenObj,
	vector<unsigned short> crsublotconfigs)
{
	std::ofstream datafile;
	std::fstream cutfile;
	datafile.open("out.txt", std::fstream::app);
	cutfile.open("cuts.txt", std::fstream::app);
	cout<<"Entered separatesmall"<<endl;
	unsigned short numSublots = submap.size();
	unsigned short numCrSublots = crsublotconfigs.size();

	unsigned short crLot = submap[crsublotconfigs[0]].parentlot;
	unsigned short crVendor = submap[crsublotconfigs[0]].parentvendor;
	unsigned short crSublot = submap[crsublotconfigs[0]].sublotsused;

	double * ctime = new double[numSublots];
	cout<<endl;
	for(unsigned short i=0;i<TLots;++i)
	{
		
		for(unsigned short j=0;j<TLots;++j)
		{
			
			cout<<precedence[i][j]<<" ";

		}	
		cout<<endl;

	}
	cout<<"critical sublot config="<<submap_inv[crLot][crVendor][crSublot]<<endl;
	cout<<"numsublots="<<numSublots<<endl;
	cout<<"ctime:"<<endl;
	for(unsigned short i=0;i<numSublots ;++i)
	{
		ctime[i] = 0;
		unsigned short crl = submap[i].parentlot;
		unsigned short crv = submap[i].parentvendor;
		unsigned short crs= submap[i].sublotsused;
		
		for(unsigned short kk = 0; kk < TLots; ++kk)
				{
					if(kk != crl)
					ctime[i]+= (lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*precedence[kk][crl];
					
					
				}
		ctime[i]+= lvdata[crl][crv].makespansmp[crs]- lvdata[crl][crv].assemblytime ;
					if(ctime[i]==ctime[submap_inv[crLot][crVendor][crSublot]])
						cout<<i<<endl;
	}
	

	IloEnv workerEnv;
	

	//Declare Variables
	IloNumVar lamda(workerEnv,0,1);
	
	
	FloatVarArray1	z(workerEnv, numSublots,0,1);

	
	FloatVarArray1	ylamda(workerEnv, TLots,0,1);

	FloatVarArray1 zlamda(workerEnv, numSublots,0,1);
	 
	
	FloatVarArray2	zyi(workerEnv, numSublots);	
	for(unsigned short j=0; j<numSublots; ++j)
		{
			
			zyi[j] = FloatVarArray1(workerEnv, TLots,0,1);
			
		}
	

	
	FloatVarArray2	zyj(workerEnv, numSublots);
	
		for(unsigned short j=0; j<numSublots; ++j)
		{
			zyj[j] = FloatVarArray1(workerEnv, TLots,0,1);
		}
	
	cout<<"Declared varaibles"<<endl;
	//create values for setvector
	IloNumVarArray vars(workerEnv);
	IloNumArray vals(workerEnv);
	{
		unsigned short i = crsublotconfigs[0];
		

		vars.add(lamda);
		vals.add(1);
		for(unsigned short j=0; j<numSublots; ++j)
		{
			if(i!=j)
			{
				vars.add(z[j]);
				if(ctime[j] > ctime[i])
					vals.add(1);
				else
					vals.add(0);
					
			}
		}

		for(unsigned short j=0; j<numSublots; ++j)
		{
			if(i!=j)
			{
				vars.add(zlamda[j]);
				if(ctime[j] > ctime[submap_inv[crLot][crVendor][crSublot]])
					vals.add(1);
				else
					vals.add(0);
					
			}
		}

		for(unsigned short j=0; j<TLots; ++j)
		{
			unsigned short crlot_temp = submap[i].parentlot;
			if(crlot_temp !=j)
			{
				vars.add(ylamda[j]);
				vals.add(precedence[j][crlot_temp]);
					
			}
		}

		for(unsigned short j=0; j<numSublots; ++j)
		{
			if(i!=j)
			{
				double tempz;
				if(ctime[j] > ctime[submap_inv[crLot][crVendor][crSublot]])
					tempz =1;
				else
					tempz=0;

				for(unsigned short k=0; k<TLots; ++k)
				{
					unsigned short crlot_temp = submap[i].parentlot;
					if(crlot_temp !=k)
					{
						vars.add(zyi[j][k]);
						vals.add(tempz*precedence[k][crlot_temp]);
					}

					unsigned short crlotj_temp = submap[j].parentlot;
					if(crlot_temp !=k)
					{
						vars.add(zyj[j][k]);
						vals.add(tempz*precedence[k][crlotj_temp]);
					}
				}
			}
		}



	}//end of set vector 

	cout<<"Set initial vector values"<<endl;
	
	// Declare model
	IloModel workerModel(workerEnv);

	//Declare Constraints
	
	IloExpr gammaex(workerEnv);
	IloExpr workerObj(workerEnv);
	
	cout<<"Attempting to create constraints"<<endl;

	IloRangeArray pi1(workerEnv,numSublots);
	

	IloRangeArray pi2(workerEnv,numSublots);
	

	IloRangeArray pi3(workerEnv,numSublots);
	

	IloRangeArray w1(workerEnv,numSublots);
	

	IloRangeArray w2(workerEnv,numSublots);
	

	IloRangeArray w3(workerEnv,numSublots);
	

	IloRangeArray chat1(workerEnv,numSublots);
	

	IloRangeArray chat2(workerEnv,numSublots);
	

	IloRangeArray2 the1(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
		the1[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 the2(workerEnv,numSublots);
		for(unsigned short j=0; j<numSublots; ++j)
			the2[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 the3(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			the3[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 phi1(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			phi1[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 phi2(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			phi2[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 phi3(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			phi3[j] = IloRangeArray(workerEnv, TLots);

	//Add constraints	
	//for(unsigned short i_num=0;i_num<numCrSublots ;++i_num)
	{
		unsigned short i = crsublotconfigs[0];
		unsigned short crl = submap[i].parentlot;
		unsigned short crv = submap[i].parentvendor;
		unsigned short crs= submap[i].sublotsused;

		if(crl==crLot && crv==crVendor && crs==crSublot)
			gammaex+= lamda;

		//lamda in objective
		workerObj+= (lvdata[crl][crv].makespansmp[crs] 
		-lvdata[crl][crv].assemblytime)*lamda;

		//yl in objective
		for(unsigned short lk = 0; lk < TLots; ++ lk)
		{
			if(lk != crl)
			{
				workerObj+= (lvdata[lk][crv].vendortime - lvdata[lk][crv].assemblytime)*ylamda[lk];
			}
		}

		//zlamda in objective
		for(unsigned short subk = 0; subk < numSublots; ++ subk)
		{
			unsigned short subk_l = submap[subk].parentlot;
			unsigned short subk_v = submap[subk].parentvendor;
			unsigned short subk_s = submap[subk].sublotsused;

			if(subk != i)
			{
				if(subk_s <= lvdata[subk_l][subk_v].maxsublots - 2)// 2 since -1 is last element index 
					workerObj+= (lvdata[subk_l][subk_v].handlingcost[subk_s +1] -
					lvdata[subk_l][subk_v].handlingcost[subk_s] )*zlamda[subk];
				else
					workerObj+= 1*zlamda[subk];			

			}
		}
		
	
		//pi
		
	for(unsigned short lk = 0; lk < TLots; ++ lk)
		{			
			if(lk != crl)
			{
				IloExpr temp1(workerEnv);
				IloExpr temp2(workerEnv);
				IloExpr temp3(workerEnv);

				temp1+= lamda-ylamda[lk];
				temp2+= -ylamda[lk];
				temp3+=-lamda+ylamda[lk];
				

				pi1[lk] = IloRange(workerEnv,0,temp1,IloInfinity);
				pi2[lk] = IloRange(workerEnv,-1*precedence[lk][crl],temp2,IloInfinity);
				pi3[lk] = IloRange(workerEnv,-1*precedence[crl][lk],temp3,IloInfinity);

				workerModel.add(pi1[lk]);
				workerModel.add(pi2[lk]);
				workerModel.add(pi3[lk]);

				temp1.end();
				temp2.end();
				temp3.end();

				
			}	

			
		}

	//w
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
		{			
			if(jbu != i)
			{
				

				IloExpr temp1(workerEnv);
				IloExpr temp2(workerEnv);
				IloExpr temp3(workerEnv);

				temp1+= lamda-zlamda[jbu];
				temp2+= z[jbu]-zlamda[jbu];
				if(lvdata[submap[jbu].parentlot][submap[jbu].parentvendor].maxsublots-1 == submap[jbu].sublotsused)
					temp3+=  -lamda -z[jbu];
				else
					temp3+= zlamda[jbu] -lamda -z[jbu];
				

				w1[jbu] = IloRange(workerEnv,0,temp1,IloInfinity);
				w2[jbu] = IloRange(workerEnv,0,temp2,IloInfinity);
				w3[jbu] = IloRange(workerEnv,-1,temp3,IloInfinity);

				workerModel.add(w1[jbu]);
				workerModel.add(w2[jbu]);
				workerModel.add(w3[jbu]);

				temp1.end();
				temp2.end();
				temp3.end();
				
			}	

		}

	

	//chat
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
		{			
			if(jbu != i)
			{
				
				unsigned short jbu_l = submap[jbu].parentlot;
				unsigned short jbu_v = submap[jbu].parentvendor;
				unsigned short jbu_s = submap[jbu].sublotsused;
				

				IloExpr temp1(workerEnv);
				IloNum tempRHS = 0;

				char zname[100];
				sprintf(zname, "z.%d.%d", (int) i, (int) jbu); 
				z[jbu].setName(zname);
				
				for(unsigned short kk = 0; kk < TLots; ++kk)
				{	

					if(kk != jbu_l)
					{
						char varNameyj[100];
						sprintf(varNameyj, "zj.%d.%d_y.%d.%d", (int) i, (int) jbu, (int) kk, (int) jbu_l); 
						zyj[jbu][kk].setName(varNameyj);

						temp1+= (lvdata[kk][jbu_v].vendortime - lvdata[kk][jbu_v].assemblytime)*zyj[jbu][kk];
						tempRHS+= (lvdata[kk][jbu_v].vendortime - lvdata[kk][jbu_v].assemblytime)*precedence[kk][jbu_l];
					}

					if(kk != crl)
					{
						char varNamezyi[100];
						sprintf(varNamezyi, "zi.%d.%d_y.%d.%d", (int) i, (int) jbu, (int) kk, (int) crl); 
						zyi[jbu][kk].setName(varNamezyi);

						temp1+= -1*(lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*zyi[jbu][kk];
						tempRHS+= -1*(lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*precedence[kk][crl];


					}

					


				}//end kk
				
				
				//cout<<jbu_l<<" "<<jbu_v<<" " <<jbu_s<<endl<<"max subs="<<lvdata[jbu_l][jbu_v].maxsublots<<endl;
				
				temp1+= (lvdata[jbu_l][jbu_v].makespansmp[jbu_s] - lvdata[crl][crv].makespansmp[crs] + 
						-lvdata[jbu_l][jbu_v].assemblytime + lvdata[crl][crv].assemblytime)*z[jbu];
				
					tempRHS+= (lvdata[jbu_l][jbu_v].makespansmp[jbu_s] - lvdata[crl][crv].makespansmp[crs] + 
						-lvdata[jbu_l][jbu_v].assemblytime 
						+lvdata[crl][crv].assemblytime);

					

				chat1[jbu] = IloRange(workerEnv,tempRHS,temp1,IloInfinity);	
				char varName[100];
				sprintf(varName, "chat1.%d.%d", (int) i, (int) jbu); 
				chat1[jbu].setName(varName);

				chat2[jbu] = IloRange(workerEnv,0,temp1,IloInfinity);	
				char varName2[100];
				sprintf(varName2, "chat2.%d.%d", (int) i, (int) jbu); 
				chat2[jbu].setName(varName2);

				workerModel.add(chat1[jbu]);
				workerModel.add(chat2[jbu]);

				

				temp1.end();
				
				
			}//endif	

		}//end jbu
	
	//theta
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
	{			
		if(jbu != i)
		{
			for(unsigned short lk = 0; lk < TLots; ++ lk)
			{	
				unsigned short jbu_l = submap[jbu].parentlot;
				if(lk != jbu_l)
				{
					

					IloExpr temp1(workerEnv);
					IloExpr temp2(workerEnv);
					IloExpr temp3(workerEnv);

					temp1+= -zyj[jbu][lk];
					temp2+=  -zyj[jbu][lk] + z[jbu];
					temp3+= zyj[jbu][lk] - z[jbu];


					the1[jbu][lk] = IloRange(workerEnv,-1*precedence[lk][jbu_l],temp1,IloInfinity);
					the2[jbu][lk] = IloRange(workerEnv,0,temp2,IloInfinity);
					the3[jbu][lk] = IloRange(workerEnv,-1*precedence[jbu_l][lk],temp3,IloInfinity);

					workerModel.add(the1[jbu][lk]);
					workerModel.add(the2[jbu][lk]);
					workerModel.add(the3[jbu][lk]);

					temp1.end();
					temp2.end();
					temp3.end();
				}}

		}}//end theta

	//phi
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
	{			
		if(jbu != i)
		{
			char varName[100];
				sprintf(varName, "z.%d.%d", (int) i, (int) jbu); 
				z[jbu].setName(varName);

			for(unsigned short lk = 0; lk < TLots; ++ lk)
			{	
				unsigned short cbu_j = submap[jbu].parentlot;
				if(lk != crl)
				{
					

					IloExpr temp1(workerEnv);
					IloExpr temp2(workerEnv);
					IloExpr temp3(workerEnv);

					temp1+=  -zyi[jbu][lk];
					temp2+=  -zyi[jbu][lk]  +z[jbu];
					temp3+=   zyi[jbu][lk] - z[jbu];


					phi1[jbu][lk] = IloRange(workerEnv,-1*precedence[lk][crl],temp1,IloInfinity);
					phi2[jbu][lk] = IloRange(workerEnv,0,temp2,IloInfinity);
					phi3[jbu][lk] = IloRange(workerEnv,-1*precedence[crl][lk],temp3,IloInfinity);

					/*if(i == 0 && jbu == 1)
				{
					cout<<"phi 1 "<<phi1[i][jbu][lk].getExpr()<<">="<<phi1[i][jbu][lk].getLB()<<endl;
					cout<<"phi 2 "<<phi2[i][jbu][lk].getExpr()<<">="<<phi2[i][jbu][lk].getLB()<<endl;
					cout<<"phi 3 "<<phi3[i][jbu][lk].getExpr()<<">="<<phi3[i][jbu][lk].getLB()<<endl;
					}*/

					workerModel.add(phi1[jbu][lk]);
					workerModel.add(phi2[jbu][lk]);
					workerModel.add(phi3[jbu][lk]);

					temp1.end();
					temp2.end();
					temp3.end();
				}}

		}}//end phi

	}// end of big subs loop

	cout<<"Added all constraints"<<endl;


	IloRange gamma(workerEnv,1,gammaex,IloInfinity,"gamma");
	workerModel.add(gamma);
	
	workerModel.add(IloMinimize(workerEnv,workerObj));

	IloCplex cplex(workerModel);
	cplex.setParam(IloCplex::PreDual,1);
	cplex.setVectors(vals, 0, vars,0,0,0);
	cplex.setParam(IloCplex::PreInd, false);
	cplex.setParam(IloCplex::WorkMem, MaxMemory);
	//cplex.setOut(workerEnv.getNullStream());

	cout<<"Attempting worker solution"<<endl;
	if(cplex.solve())
	{
		double y_factor=0;
		cout<<"Solved worker cplex"<<endl;

		cut_yij_coeff.clear();		
		cut_yij_coeff.resize(TLots);		

		for(unsigned short i=0; i<TLots; ++i)
		{
			cut_yij_coeff[i].resize(TLots);
			fill(cut_yij_coeff[i].begin(), cut_yij_coeff[i].end(), 0);
			
		}

		

		unsigned short icr = submap_inv[crLot][crVendor][crSublot];
		cout<<"trying to add pi duals"<<endl;
		//add pi
		for(unsigned short k=0; k<TLots; ++k)
		{
			if(k!=crLot)
			{
				double temp= cplex.getDual(pi2[k]);
				cut_yij_coeff[k][crLot]+= -1*temp;

				temp= cplex.getDual(pi3[k]);
				cut_yij_coeff[crLot][k]= -1*temp;
			}
		}

		
		cout<<"trying to add chat duals"<<endl;
		//add chat1
		for(unsigned short jbu=0; jbu<numSublots; ++jbu)
		{
			unsigned short jbu_l = submap[jbu].parentlot;
			unsigned short jbu_v = submap[jbu].parentvendor;

			if(jbu != icr)
			{
				double chat1_jbu_dual = cplex.getDual(chat1[jbu]);
				for(unsigned short k=0; k<TLots; ++k)
				{
					if(jbu_l != k)
					{
						double temp = (lvdata[k][jbu_v].vendortime - lvdata[k][jbu_v].assemblytime);
						cut_yij_coeff[k][jbu_l]+= temp*chat1_jbu_dual;
						if(abs(chat1_jbu_dual)>0)
						{
						cout<<"y co-effs during chat1 first loop \n";
						cout<<"jbu_l="<<jbu_l<<" k="<<k
							<<" temp="<<temp<<" dual="<<chat1_jbu_dual<<endl;
		for(unsigned short i=0; i<TLots; ++i)
		{	
			
			for(unsigned short j=0; j<TLots; ++j)
			cout<<cut_yij_coeff[i][j]<<" ";
			cout<<endl;
		}
		}
					}

					if(k != crLot)
					{
						cut_yij_coeff[k][crLot] += -1*chat1_jbu_dual*(lvdata[k][crVendor].vendortime - lvdata[k][crVendor].assemblytime);
					}
				}
				

			}
		}
		
		cout<<"trying to add phi and theta duals"<<endl;
		//add theta and phi
		for(unsigned short jbu=0; jbu<numSublots; ++jbu)
		{
			unsigned short jbu_l = submap[jbu].parentlot;
			unsigned short jbu_v = submap[jbu].parentvendor;

			if(jbu != icr)
			{
				
				for(unsigned short k=0; k<TLots; ++k)
				{
					if(jbu_l != k)
					{
						double temp = -1* cplex.getDual(the1[jbu][k]);
						cut_yij_coeff[k][jbu_l]+= temp;

						temp = -1* cplex.getDual(the3[jbu][k]);
						cut_yij_coeff[jbu_l][k] += temp;
					}

					if(k != crLot)
					{
						double temp = -1* cplex.getDual(phi1[jbu][k]);
						cut_yij_coeff[k][crLot]+= temp;

						temp = -1* cplex.getDual(phi3[jbu][k]);
						cut_yij_coeff[crLot][k] += temp;
					}
				}
				

			}
		}
		
		
		//create cut LHS
		cout<<"y co-effs \n";

		cutfile<< "for sequence: \n";
		for(unsigned short i=0; i<TLots; ++i)
		{
			
			for(unsigned short j=0; j<TLots; ++j)
			{
				
				cutfile<<precedence[i][j]<<" ";
			}
			cutfile<<endl;
		}//add y factors
		cutfile<<"cut added: \n";
		for(unsigned short i=0; i<TLots; ++i)
		{
			
			for(unsigned short j=0; j<TLots; ++j)
			{
				cout<<cut_yij_coeff[i][j]<<" ";
				cutLhs+= cut_yij_coeff[i][j]*y[i][j];
				y_factor+= cut_yij_coeff[i][j]*precedence[i][j];
				if(cut_yij_coeff[i][j] != 0)
					cutfile<<endl<<cut_yij_coeff[i][j]<<"*y["<<i<<"]["<<j<<"] ";
			}
			cout<<endl;
		}//add y factors
		cout<<"y factor="<<y_factor<<endl;
		cout<<"slave objective="<<cplex.getObjValue()<<endl;
		cutLhs+= cplex.getObjValue()-y_factor;
		cutfile<<" + "<<cplex.getObjValue()-y_factor<<" <= z_master \n";
		cout<<cutLhs<<endl;
	}

	

	cout<<"here\n";
	//workerModel.end();
	double total_obj = cplex.getObjValue()+ObjConstant;
	workerEnv.end();
	cout<<"hear\n";

	datafile.close();
	cutfile.close();
	return total_obj;
}


double datalap::separatepermx(IloExpr cutLhs,  BoolVarArray2 x_master)
	
{
	

	std::ofstream datafile;
	std::fstream cutfile;
	datafile.open("out.txt", std::fstream::app);
	cutfile.open("cuts.txt", std::fstream::app);
	//datalap::master_instance->clean_perm(datalap::master_instance->m_perm);
	precedence.resize(TLots);
	cout<<"perm received from master:"<<endl;
	for(unsigned short i=0;i<TLots;++i)
	{
		
		precedence[i].resize(TLots);

		for(unsigned short j=0;j<TLots;++j)
		{			
			
			if(m_perm[i][j]<0.1)
				m_perm[i][j]=0;
			else
				m_perm[i][j]=1;
			cout<<m_perm[i][j];

		}
		cout<<endl;

	}	
	datalap::master_instance->perm_to_prec();
	cout<<"prec from master:"<<endl;
	datalap::master_instance->printprec();
	
	
	double LowerBound = givensequence();
	cout<<"Critical Lot="<<critical_lot<<endl;
	cout<<"Critical Sublot="<<critical_sublot<<endl;
	cout<<"Critical Vendor="<<critical_vendor<<endl;
	cout<<LowerBound<<" from given seq.\n";
		
	
	
	


	vector<unsigned short> crsublotconfigs;
	crsublotconfigs.push_back(submap_inv[critical_lot][critical_vendor][critical_sublot-1]);
	
	cout<<"Entered solver phase"<<endl;
	unsigned short numSublots = submap.size();
	unsigned short numCrSublots = crsublotconfigs.size();

	unsigned short crLot = submap[crsublotconfigs[0]].parentlot;
	unsigned short crVendor = submap[crsublotconfigs[0]].parentvendor;
	unsigned short crSublot = submap[crsublotconfigs[0]].sublotsused;

	double * ctime = new double[numSublots];
	cout<<endl;
	for(unsigned short i=0;i<TLots;++i)
	{
		
		for(unsigned short j=0;j<TLots;++j)
		{
			
			cout<<precedence[i][j]<<" ";

		}	
		cout<<endl;

	}
	cout<<"critical sublot config="<<submap_inv[crLot][crVendor][crSublot]<<endl;
	cout<<"numsublots="<<numSublots<<endl;
	cout<<"ctime:"<<endl;
	for(unsigned short i=0;i<numSublots ;++i)
	{
		ctime[i] = 0;
		unsigned short crl = submap[i].parentlot;
		unsigned short crv = submap[i].parentvendor;
		unsigned short crs= submap[i].sublotsused;
		
		for(unsigned short kk = 0; kk < TLots; ++kk)
				{
					if(kk != crl)
					ctime[i]+= (lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*precedence[kk][crl];
					
					
				}
		ctime[i]+= lvdata[crl][crv].makespansmp[crs]- lvdata[crl][crv].assemblytime ;
					if(ctime[i]==ctime[submap_inv[crLot][crVendor][crSublot]])
						cout<<i<<endl;
	}
	

	IloEnv workerEnv;
	

	//Declare Variables
	IloNumVar lamda(workerEnv,0,1);
	
	
	FloatVarArray1	z(workerEnv, numSublots,0,1);

	
	FloatVarArray1	ylamda(workerEnv, TLots,0,1);

	FloatVarArray1 zlamda(workerEnv, numSublots,0,1);
	 
	
	FloatVarArray2	zyi(workerEnv, numSublots);	
	for(unsigned short j=0; j<numSublots; ++j)
		{
			
			zyi[j] = FloatVarArray1(workerEnv, TLots,0,1);
			
		}
	

	
	FloatVarArray2	zyj(workerEnv, numSublots);
	
		for(unsigned short j=0; j<numSublots; ++j)
		{
			zyj[j] = FloatVarArray1(workerEnv, TLots,0,1);
		}
	
	cout<<"Declared varaibles"<<endl;
	//create values for setvector
	IloNumVarArray vars(workerEnv);
	IloNumArray vals(workerEnv);
	{
		unsigned short i = crsublotconfigs[0];
		

		vars.add(lamda);
		vals.add(1);
		for(unsigned short j=0; j<numSublots; ++j)
		{
			if(i!=j)
			{
				vars.add(z[j]);
				if(ctime[j] > ctime[i])
					vals.add(1);
				else
					vals.add(0);
					
			}
		}

		for(unsigned short j=0; j<numSublots; ++j)
		{
			if(i!=j)
			{
				vars.add(zlamda[j]);
				if(ctime[j] > ctime[submap_inv[crLot][crVendor][crSublot]])
					vals.add(1);
				else
					vals.add(0);
					
			}
		}

		for(unsigned short j=0; j<TLots; ++j)
		{
			unsigned short crlot_temp = submap[i].parentlot;
			if(crlot_temp !=j)
			{
				vars.add(ylamda[j]);
				vals.add(precedence[j][crlot_temp]);
					
			}
		}

		for(unsigned short j=0; j<numSublots; ++j)
		{
			if(i!=j)
			{
				double tempz;
				if(ctime[j] > ctime[submap_inv[crLot][crVendor][crSublot]])
					tempz =1;
				else
					tempz=0;

				for(unsigned short k=0; k<TLots; ++k)
				{
					unsigned short crlot_temp = submap[i].parentlot;
					if(crlot_temp !=k)
					{
						vars.add(zyi[j][k]);
						vals.add(tempz*precedence[k][crlot_temp]);
					}

					unsigned short crlotj_temp = submap[j].parentlot;
					if(crlot_temp !=k)
					{
						vars.add(zyj[j][k]);
						vals.add(tempz*precedence[k][crlotj_temp]);
					}
				}
			}
		}



	}//end of set vector 

	cout<<"Set initial vector values"<<endl;
	
	// Declare model
	IloModel workerModel(workerEnv);

	//Declare Constraints
	
	IloExpr gammaex(workerEnv);
	IloExpr workerObj(workerEnv);
	
	cout<<"Attempting to create constraints"<<endl;

	IloRangeArray pi1(workerEnv,numSublots);
	

	IloRangeArray pi2(workerEnv,numSublots);
	

	IloRangeArray pi3(workerEnv,numSublots);
	

	IloRangeArray w1(workerEnv,numSublots);
	

	IloRangeArray w2(workerEnv,numSublots);
	

	IloRangeArray w3(workerEnv,numSublots);
	

	IloRangeArray chat1(workerEnv,numSublots);
	

	IloRangeArray chat2(workerEnv,numSublots);
	

	IloRangeArray2 the1(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
		the1[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 the2(workerEnv,numSublots);
		for(unsigned short j=0; j<numSublots; ++j)
			the2[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 the3(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			the3[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 phi1(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			phi1[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 phi2(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			phi2[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 phi3(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			phi3[j] = IloRangeArray(workerEnv, TLots);

	//Add constraints	
	//for(unsigned short i_num=0;i_num<numCrSublots ;++i_num)
	{
		unsigned short i = crsublotconfigs[0];
		unsigned short crl = submap[i].parentlot;
		unsigned short crv = submap[i].parentvendor;
		unsigned short crs= submap[i].sublotsused;

		if(crl==crLot && crv==crVendor && crs==crSublot)
			gammaex+= lamda;

		//lamda in objective
		workerObj+= (lvdata[crl][crv].makespansmp[crs] 
		-lvdata[crl][crv].assemblytime)*lamda;

		//yl in objective
		for(unsigned short lk = 0; lk < TLots; ++ lk)
		{
			if(lk != crl)
			{
				workerObj+= (lvdata[lk][crv].vendortime - lvdata[lk][crv].assemblytime)*ylamda[lk];
			}
		}

		//zlamda in objective
		for(unsigned short subk = 0; subk < numSublots; ++ subk)
		{
			unsigned short subk_l = submap[subk].parentlot;
			unsigned short subk_v = submap[subk].parentvendor;
			unsigned short subk_s = submap[subk].sublotsused;

			if(subk != i)
			{
				if(subk_s <= lvdata[subk_l][subk_v].maxsublots - 2)// 2 since -1 is last element index 
					workerObj+= (lvdata[subk_l][subk_v].handlingcost[subk_s +1] -
					lvdata[subk_l][subk_v].handlingcost[subk_s] )*zlamda[subk];
				else
					workerObj+= 1*zlamda[subk];			

			}
		}
		
	
		//pi
		
	for(unsigned short lk = 0; lk < TLots; ++ lk)
		{			
			if(lk != crl)
			{
				IloExpr temp1(workerEnv);
				IloExpr temp2(workerEnv);
				IloExpr temp3(workerEnv);

				temp1+= lamda-ylamda[lk];
				temp2+= -ylamda[lk];
				temp3+=-lamda+ylamda[lk];
				

				pi1[lk] = IloRange(workerEnv,0,temp1,IloInfinity);
				pi2[lk] = IloRange(workerEnv,-1*precedence[lk][crl],temp2,IloInfinity);
				pi3[lk] = IloRange(workerEnv,-1*precedence[crl][lk],temp3,IloInfinity);

				workerModel.add(pi1[lk]);
				workerModel.add(pi2[lk]);
				workerModel.add(pi3[lk]);

				temp1.end();
				temp2.end();
				temp3.end();

				
			}	

			
		}

	//w
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
		{			
			if(jbu != i)
			{
				

				IloExpr temp1(workerEnv);
				IloExpr temp2(workerEnv);
				IloExpr temp3(workerEnv);

				temp1+= lamda-zlamda[jbu];
				temp2+= z[jbu]-zlamda[jbu];
				if(lvdata[submap[jbu].parentlot][submap[jbu].parentvendor].maxsublots-1 == submap[jbu].sublotsused)
					temp3+=  -lamda -z[jbu];
				else
					temp3+= zlamda[jbu] -lamda -z[jbu];
				

				w1[jbu] = IloRange(workerEnv,0,temp1,IloInfinity);
				w2[jbu] = IloRange(workerEnv,0,temp2,IloInfinity);
				w3[jbu] = IloRange(workerEnv,-1,temp3,IloInfinity);

				workerModel.add(w1[jbu]);
				workerModel.add(w2[jbu]);
				workerModel.add(w3[jbu]);

				temp1.end();
				temp2.end();
				temp3.end();
				
			}	

		}

	

	//chat
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
		{			
			if(jbu != i)
			{
				
				unsigned short jbu_l = submap[jbu].parentlot;
				unsigned short jbu_v = submap[jbu].parentvendor;
				unsigned short jbu_s = submap[jbu].sublotsused;
				

				IloExpr temp1(workerEnv);
				IloNum tempRHS = 0;

				char zname[100];
				sprintf(zname, "z.%d.%d", (int) i, (int) jbu); 
				z[jbu].setName(zname);
				
				for(unsigned short kk = 0; kk < TLots; ++kk)
				{	

					if(kk != jbu_l)
					{
						char varNameyj[100];
						sprintf(varNameyj, "zj.%d.%d_y.%d.%d", (int) i, (int) jbu, (int) kk, (int) jbu_l); 
						zyj[jbu][kk].setName(varNameyj);

						temp1+= (lvdata[kk][jbu_v].vendortime - lvdata[kk][jbu_v].assemblytime)*zyj[jbu][kk];
						tempRHS+= (lvdata[kk][jbu_v].vendortime - lvdata[kk][jbu_v].assemblytime)*precedence[kk][jbu_l];
					}

					if(kk != crl)
					{
						char varNamezyi[100];
						sprintf(varNamezyi, "zi.%d.%d_y.%d.%d", (int) i, (int) jbu, (int) kk, (int) crl); 
						zyi[jbu][kk].setName(varNamezyi);

						temp1+= -1*(lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*zyi[jbu][kk];
						tempRHS+= -1*(lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*precedence[kk][crl];


					}

					


				}//end kk
				
				
				//cout<<jbu_l<<" "<<jbu_v<<" " <<jbu_s<<endl<<"max subs="<<lvdata[jbu_l][jbu_v].maxsublots<<endl;
				
				temp1+= (lvdata[jbu_l][jbu_v].makespansmp[jbu_s] - lvdata[crl][crv].makespansmp[crs] + 
						-lvdata[jbu_l][jbu_v].assemblytime + lvdata[crl][crv].assemblytime)*z[jbu];
				
					tempRHS+= (lvdata[jbu_l][jbu_v].makespansmp[jbu_s] - lvdata[crl][crv].makespansmp[crs] + 
						-lvdata[jbu_l][jbu_v].assemblytime 
						+lvdata[crl][crv].assemblytime);

					

				chat1[jbu] = IloRange(workerEnv,tempRHS,temp1,IloInfinity);	
				char varName[100];
				sprintf(varName, "chat1.%d.%d", (int) i, (int) jbu); 
				chat1[jbu].setName(varName);

				chat2[jbu] = IloRange(workerEnv,0,temp1,IloInfinity);	
				char varName2[100];
				sprintf(varName2, "chat2.%d.%d", (int) i, (int) jbu); 
				chat2[jbu].setName(varName2);

				workerModel.add(chat1[jbu]);
				workerModel.add(chat2[jbu]);

				

				temp1.end();
				
				
			}//endif	

		}//end jbu
	
	//theta
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
	{			
		if(jbu != i)
		{
			for(unsigned short lk = 0; lk < TLots; ++ lk)
			{	
				unsigned short jbu_l = submap[jbu].parentlot;
				if(lk != jbu_l)
				{
					

					IloExpr temp1(workerEnv);
					IloExpr temp2(workerEnv);
					IloExpr temp3(workerEnv);

					temp1+= -zyj[jbu][lk];
					temp2+=  -zyj[jbu][lk] + z[jbu];
					temp3+= zyj[jbu][lk] - z[jbu];


					the1[jbu][lk] = IloRange(workerEnv,-1*precedence[lk][jbu_l],temp1,IloInfinity);
					the2[jbu][lk] = IloRange(workerEnv,0,temp2,IloInfinity);
					the3[jbu][lk] = IloRange(workerEnv,-1*precedence[jbu_l][lk],temp3,IloInfinity);

					workerModel.add(the1[jbu][lk]);
					workerModel.add(the2[jbu][lk]);
					workerModel.add(the3[jbu][lk]);

					temp1.end();
					temp2.end();
					temp3.end();
				}}

		}}//end theta

	//phi
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
	{			
		if(jbu != i)
		{
			char varName[100];
				sprintf(varName, "z.%d.%d", (int) i, (int) jbu); 
				z[jbu].setName(varName);

			for(unsigned short lk = 0; lk < TLots; ++ lk)
			{	
				unsigned short cbu_j = submap[jbu].parentlot;
				if(lk != crl)
				{
					

					IloExpr temp1(workerEnv);
					IloExpr temp2(workerEnv);
					IloExpr temp3(workerEnv);

					temp1+=  -zyi[jbu][lk];
					temp2+=  -zyi[jbu][lk]  +z[jbu];
					temp3+=   zyi[jbu][lk] - z[jbu];


					phi1[jbu][lk] = IloRange(workerEnv,-1*precedence[lk][crl],temp1,IloInfinity);
					phi2[jbu][lk] = IloRange(workerEnv,0,temp2,IloInfinity);
					phi3[jbu][lk] = IloRange(workerEnv,-1*precedence[crl][lk],temp3,IloInfinity);

					/*if(i == 0 && jbu == 1)
				{
					cout<<"phi 1 "<<phi1[i][jbu][lk].getExpr()<<">="<<phi1[i][jbu][lk].getLB()<<endl;
					cout<<"phi 2 "<<phi2[i][jbu][lk].getExpr()<<">="<<phi2[i][jbu][lk].getLB()<<endl;
					cout<<"phi 3 "<<phi3[i][jbu][lk].getExpr()<<">="<<phi3[i][jbu][lk].getLB()<<endl;
					}*/

					workerModel.add(phi1[jbu][lk]);
					workerModel.add(phi2[jbu][lk]);
					workerModel.add(phi3[jbu][lk]);

					temp1.end();
					temp2.end();
					temp3.end();
				}}

		}}//end phi

	}// end of big subs loop

	cout<<"Added all constraints"<<endl;


	IloRange gamma(workerEnv,1,gammaex,IloInfinity,"gamma");
	workerModel.add(gamma);
	
	workerModel.add(IloMinimize(workerEnv,workerObj));

	IloCplex cplex(workerModel);
	cplex.setParam(IloCplex::PreDual,1);
	cplex.setVectors(vals, 0, vars,0,0,0);
	cplex.setParam(IloCplex::PreInd, false);
	cplex.setParam(IloCplex::WorkMem, MaxMemory);
	//cplex.setOut(workerEnv.getNullStream());

	cout<<"Attempting worker solution"<<endl;
	if(cplex.solve())
	{
		double y_factor=0;
		cout<<"Solved worker cplex"<<endl;

		cut_yij_coeff.clear();		
		cut_yij_coeff.resize(TLots);		

		for(unsigned short i=0; i<TLots; ++i)
		{
			cut_yij_coeff[i].resize(TLots);
			fill(cut_yij_coeff[i].begin(), cut_yij_coeff[i].end(), 0);
			
		}

		

		unsigned short icr = submap_inv[crLot][crVendor][crSublot];
		cout<<"trying to add pi duals"<<endl;
		//add pi
		for(unsigned short k=0; k<TLots; ++k)
		{
			if(k!=crLot)
			{
				double temp= cplex.getDual(pi2[k]);
				cut_yij_coeff[k][crLot]+= -1*temp;

				temp= cplex.getDual(pi3[k]);
				cut_yij_coeff[crLot][k]= -1*temp;
			}
		}

		
		cout<<"trying to add chat duals"<<endl;
		//add chat1
		for(unsigned short jbu=0; jbu<numSublots; ++jbu)
		{
			unsigned short jbu_l = submap[jbu].parentlot;
			unsigned short jbu_v = submap[jbu].parentvendor;

			if(jbu != icr)
			{
				double chat1_jbu_dual = cplex.getDual(chat1[jbu]);
				for(unsigned short k=0; k<TLots; ++k)
				{
					if(jbu_l != k)
					{
						double temp = (lvdata[k][jbu_v].vendortime - lvdata[k][jbu_v].assemblytime);
						cut_yij_coeff[k][jbu_l]+= temp*chat1_jbu_dual;
						if(abs(chat1_jbu_dual)>0)
						{
						cout<<"y co-effs during chat1 first loop \n";
						cout<<"jbu_l="<<jbu_l<<" k="<<k
							<<" temp="<<temp<<" dual="<<chat1_jbu_dual<<endl;
		for(unsigned short i=0; i<TLots; ++i)
		{	
			
			for(unsigned short j=0; j<TLots; ++j)
			cout<<cut_yij_coeff[i][j]<<" ";
			cout<<endl;
		}
		}
					}

					if(k != crLot)
					{
						cut_yij_coeff[k][crLot] += -1*chat1_jbu_dual*(lvdata[k][crVendor].vendortime - lvdata[k][crVendor].assemblytime);
					}
				}
				

			}
		}
		
		cout<<"trying to add phi and theta duals"<<endl;
		//add theta and phi
		for(unsigned short jbu=0; jbu<numSublots; ++jbu)
		{
			unsigned short jbu_l = submap[jbu].parentlot;
			unsigned short jbu_v = submap[jbu].parentvendor;

			if(jbu != icr)
			{
				
				for(unsigned short k=0; k<TLots; ++k)
				{
					if(jbu_l != k)
					{
						double temp = -1* cplex.getDual(the1[jbu][k]);
						cut_yij_coeff[k][jbu_l]+= temp;

						temp = -1* cplex.getDual(the3[jbu][k]);
						cut_yij_coeff[jbu_l][k] += temp;
					}

					if(k != crLot)
					{
						double temp = -1* cplex.getDual(phi1[jbu][k]);
						cut_yij_coeff[k][crLot]+= temp;

						temp = -1* cplex.getDual(phi3[jbu][k]);
						cut_yij_coeff[crLot][k] += temp;
					}
				}
				

			}
		}
		
		
		//create cut LHS
		cout<<"y co-effs \n";

		cutfile<< "for sequence: \n";
		for(unsigned short i=0; i<TLots; ++i)
		{
			
			for(unsigned short j=0; j<TLots; ++j)
			{
				
				cutfile<<precedence[i][j]<<" ";
			}
			cutfile<<endl;
		}//add y factors
		cutfile<<"cut added: \n";
		unsigned short imp_locn = 0;
		for(unsigned short i=0; i<TLots; ++i)
		{
			
			for(unsigned short j=0; j<TLots; ++j)
			{
				cout<<cut_yij_coeff[i][j]<<" ";
				if(cut_yij_coeff[i][j] >= 0)
					cutLhs+= cut_yij_coeff[i][j]*x_master[i][imp_locn]+cut_yij_coeff[i][j]*x_master[j][imp_locn+1] - cut_yij_coeff[i][j];
				else
					cutLhs+= -1*cut_yij_coeff[i][j]*x_master[i][imp_locn+1]-cut_yij_coeff[i][j]*x_master[j][imp_locn] + 2*cut_yij_coeff[i][j];
				
				y_factor+= cut_yij_coeff[i][j]*precedence[i][j];
				if(cut_yij_coeff[i][j] != 0)
					cutfile<<endl<<cut_yij_coeff[i][j]<<"*y["<<i<<"]["<<j<<"] ";
			}
			cout<<endl;
		}//add y factors
		cout<<"y factor="<<y_factor<<endl;
		cout<<"sub objective="<<cplex.getObjValue()<<endl;
		cutLhs+= cplex.getObjValue()-y_factor + ObjConstant;
		cutfile<<" + "<<cplex.getObjValue()-y_factor<<" <= z_master \n";
		cout<<cutLhs<<endl;
	}

	

	cout<<"here\n";
	//workerModel.end();
	double total_obj = cplex.getObjValue()+ObjConstant;
	workerEnv.end();
	cout<<"hear\n";

	datafile.close();
	cutfile.close();
	return total_obj;
	}
double datalap::separateperm(IloExpr cutLhs, 
	FloatVarArray2 y)
{
	

	std::ofstream datafile;
	std::fstream cutfile;
	datafile.open("out.txt", std::fstream::app);
	cutfile.open("cuts.txt", std::fstream::app);
	//datalap::master_instance->clean_perm(datalap::master_instance->m_perm);
	precedence.resize(TLots);
	cout<<"perm received from master:"<<endl;
	for(unsigned short i=0;i<TLots;++i)
	{
		
		precedence[i].resize(TLots);

		for(unsigned short j=0;j<TLots;++j)
		{			
			
			if(m_perm[i][j]<0.1)
				m_perm[i][j]=0;
			else
				m_perm[i][j]=1;
			cout<<m_perm[i][j];

		}
		cout<<endl;

	}	
	datalap::master_instance->perm_to_prec();
	cout<<"prec from master:"<<endl;
	datalap::master_instance->printprec();
	
	
	double LowerBound = givensequence();
	cout<<"Critical Lot="<<critical_lot<<endl;
	cout<<"Critical Sublot="<<critical_sublot<<endl;
	cout<<"Critical Vendor="<<critical_vendor<<endl;
	cout<<LowerBound<<" from given seq.\n";
		
	
	
	


	vector<unsigned short> crsublotconfigs;
	crsublotconfigs.push_back(submap_inv[critical_lot][critical_vendor][critical_sublot-1]);
	
	cout<<"Entered solver phase"<<endl;
	unsigned short numSublots = submap.size();
	unsigned short numCrSublots = crsublotconfigs.size();

	unsigned short crLot = submap[crsublotconfigs[0]].parentlot;
	unsigned short crVendor = submap[crsublotconfigs[0]].parentvendor;
	unsigned short crSublot = submap[crsublotconfigs[0]].sublotsused;

	double * ctime = new double[numSublots];
	cout<<endl;
	for(unsigned short i=0;i<TLots;++i)
	{
		
		for(unsigned short j=0;j<TLots;++j)
		{
			
			cout<<precedence[i][j]<<" ";

		}	
		cout<<endl;

	}
	cout<<"critical sublot config="<<submap_inv[crLot][crVendor][crSublot]<<endl;
	cout<<"numsublots="<<numSublots<<endl;
	cout<<"ctime:"<<endl;
	for(unsigned short i=0;i<numSublots ;++i)
	{
		ctime[i] = 0;
		unsigned short crl = submap[i].parentlot;
		unsigned short crv = submap[i].parentvendor;
		unsigned short crs= submap[i].sublotsused;
		
		for(unsigned short kk = 0; kk < TLots; ++kk)
				{
					if(kk != crl)
					ctime[i]+= (lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*precedence[kk][crl];
					
					
				}
		ctime[i]+= lvdata[crl][crv].makespansmp[crs]- lvdata[crl][crv].assemblytime ;
					if(ctime[i]==ctime[submap_inv[crLot][crVendor][crSublot]])
						cout<<i<<endl;
	}
	

	IloEnv workerEnv;
	

	//Declare Variables
	IloNumVar lamda(workerEnv,0,1);
	
	
	FloatVarArray1	z(workerEnv, numSublots,0,1);

	
	FloatVarArray1	ylamda(workerEnv, TLots,0,1);

	FloatVarArray1 zlamda(workerEnv, numSublots,0,1);
	 
	
	FloatVarArray2	zyi(workerEnv, numSublots);	
	for(unsigned short j=0; j<numSublots; ++j)
		{
			
			zyi[j] = FloatVarArray1(workerEnv, TLots,0,1);
			
		}
	

	
	FloatVarArray2	zyj(workerEnv, numSublots);
	
		for(unsigned short j=0; j<numSublots; ++j)
		{
			zyj[j] = FloatVarArray1(workerEnv, TLots,0,1);
		}
	
	cout<<"Declared varaibles"<<endl;
	//create values for setvector
	IloNumVarArray vars(workerEnv);
	IloNumArray vals(workerEnv);
	{
		unsigned short i = crsublotconfigs[0];
		

		vars.add(lamda);
		vals.add(1);
		for(unsigned short j=0; j<numSublots; ++j)
		{
			if(i!=j)
			{
				vars.add(z[j]);
				if(ctime[j] > ctime[i])
					vals.add(1);
				else
					vals.add(0);
					
			}
		}

		for(unsigned short j=0; j<numSublots; ++j)
		{
			if(i!=j)
			{
				vars.add(zlamda[j]);
				if(ctime[j] > ctime[submap_inv[crLot][crVendor][crSublot]])
					vals.add(1);
				else
					vals.add(0);
					
			}
		}

		for(unsigned short j=0; j<TLots; ++j)
		{
			unsigned short crlot_temp = submap[i].parentlot;
			if(crlot_temp !=j)
			{
				vars.add(ylamda[j]);
				vals.add(precedence[j][crlot_temp]);
					
			}
		}

		for(unsigned short j=0; j<numSublots; ++j)
		{
			if(i!=j)
			{
				double tempz;
				if(ctime[j] > ctime[submap_inv[crLot][crVendor][crSublot]])
					tempz =1;
				else
					tempz=0;

				for(unsigned short k=0; k<TLots; ++k)
				{
					unsigned short crlot_temp = submap[i].parentlot;
					if(crlot_temp !=k)
					{
						vars.add(zyi[j][k]);
						vals.add(tempz*precedence[k][crlot_temp]);
					}

					unsigned short crlotj_temp = submap[j].parentlot;
					if(crlot_temp !=k)
					{
						vars.add(zyj[j][k]);
						vals.add(tempz*precedence[k][crlotj_temp]);
					}
				}
			}
		}



	}//end of set vector 

	cout<<"Set initial vector values"<<endl;
	
	// Declare model
	IloModel workerModel(workerEnv);

	//Declare Constraints
	
	IloExpr gammaex(workerEnv);
	IloExpr workerObj(workerEnv);
	
	cout<<"Attempting to create constraints"<<endl;

	IloRangeArray pi1(workerEnv,numSublots);
	

	IloRangeArray pi2(workerEnv,numSublots);
	

	IloRangeArray pi3(workerEnv,numSublots);
	

	IloRangeArray w1(workerEnv,numSublots);
	

	IloRangeArray w2(workerEnv,numSublots);
	

	IloRangeArray w3(workerEnv,numSublots);
	

	IloRangeArray chat1(workerEnv,numSublots);
	

	IloRangeArray chat2(workerEnv,numSublots);
	

	IloRangeArray2 the1(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
		the1[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 the2(workerEnv,numSublots);
		for(unsigned short j=0; j<numSublots; ++j)
			the2[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 the3(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			the3[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 phi1(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			phi1[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 phi2(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			phi2[j] = IloRangeArray(workerEnv, TLots);

	IloRangeArray2 phi3(workerEnv,numSublots);
	for(unsigned short j=0; j<numSublots; ++j)
			phi3[j] = IloRangeArray(workerEnv, TLots);

	//Add constraints	
	//for(unsigned short i_num=0;i_num<numCrSublots ;++i_num)
	{
		unsigned short i = crsublotconfigs[0];
		unsigned short crl = submap[i].parentlot;
		unsigned short crv = submap[i].parentvendor;
		unsigned short crs= submap[i].sublotsused;

		if(crl==crLot && crv==crVendor && crs==crSublot)
			gammaex+= lamda;

		//lamda in objective
		workerObj+= (lvdata[crl][crv].makespansmp[crs] 
		-lvdata[crl][crv].assemblytime)*lamda;

		//yl in objective
		for(unsigned short lk = 0; lk < TLots; ++ lk)
		{
			if(lk != crl)
			{
				workerObj+= (lvdata[lk][crv].vendortime - lvdata[lk][crv].assemblytime)*ylamda[lk];
			}
		}

		//zlamda in objective
		for(unsigned short subk = 0; subk < numSublots; ++ subk)
		{
			unsigned short subk_l = submap[subk].parentlot;
			unsigned short subk_v = submap[subk].parentvendor;
			unsigned short subk_s = submap[subk].sublotsused;

			if(subk != i)
			{
				if(subk_s <= lvdata[subk_l][subk_v].maxsublots - 2)// 2 since -1 is last element index 
					workerObj+= (lvdata[subk_l][subk_v].handlingcost[subk_s +1] -
					lvdata[subk_l][subk_v].handlingcost[subk_s] )*zlamda[subk];
				else
					workerObj+= 1*zlamda[subk];			

			}
		}
		
	
		//pi
		
	for(unsigned short lk = 0; lk < TLots; ++ lk)
		{			
			if(lk != crl)
			{
				IloExpr temp1(workerEnv);
				IloExpr temp2(workerEnv);
				IloExpr temp3(workerEnv);

				temp1+= lamda-ylamda[lk];
				temp2+= -ylamda[lk];
				temp3+=-lamda+ylamda[lk];
				

				pi1[lk] = IloRange(workerEnv,0,temp1,IloInfinity);
				pi2[lk] = IloRange(workerEnv,-1*precedence[lk][crl],temp2,IloInfinity);
				pi3[lk] = IloRange(workerEnv,-1*precedence[crl][lk],temp3,IloInfinity);

				workerModel.add(pi1[lk]);
				workerModel.add(pi2[lk]);
				workerModel.add(pi3[lk]);

				temp1.end();
				temp2.end();
				temp3.end();

				
			}	

			
		}

	//w
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
		{			
			if(jbu != i)
			{
				

				IloExpr temp1(workerEnv);
				IloExpr temp2(workerEnv);
				IloExpr temp3(workerEnv);

				temp1+= lamda-zlamda[jbu];
				temp2+= z[jbu]-zlamda[jbu];
				if(lvdata[submap[jbu].parentlot][submap[jbu].parentvendor].maxsublots-1 == submap[jbu].sublotsused)
					temp3+=  -lamda -z[jbu];
				else
					temp3+= zlamda[jbu] -lamda -z[jbu];
				

				w1[jbu] = IloRange(workerEnv,0,temp1,IloInfinity);
				w2[jbu] = IloRange(workerEnv,0,temp2,IloInfinity);
				w3[jbu] = IloRange(workerEnv,-1,temp3,IloInfinity);

				workerModel.add(w1[jbu]);
				workerModel.add(w2[jbu]);
				workerModel.add(w3[jbu]);

				temp1.end();
				temp2.end();
				temp3.end();
				
			}	

		}

	

	//chat
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
		{			
			if(jbu != i)
			{
				
				unsigned short jbu_l = submap[jbu].parentlot;
				unsigned short jbu_v = submap[jbu].parentvendor;
				unsigned short jbu_s = submap[jbu].sublotsused;
				

				IloExpr temp1(workerEnv);
				IloNum tempRHS = 0;

				char zname[100];
				sprintf(zname, "z.%d.%d", (int) i, (int) jbu); 
				z[jbu].setName(zname);
				
				for(unsigned short kk = 0; kk < TLots; ++kk)
				{	

					if(kk != jbu_l)
					{
						char varNameyj[100];
						sprintf(varNameyj, "zj.%d.%d_y.%d.%d", (int) i, (int) jbu, (int) kk, (int) jbu_l); 
						zyj[jbu][kk].setName(varNameyj);

						temp1+= (lvdata[kk][jbu_v].vendortime - lvdata[kk][jbu_v].assemblytime)*zyj[jbu][kk];
						tempRHS+= (lvdata[kk][jbu_v].vendortime - lvdata[kk][jbu_v].assemblytime)*precedence[kk][jbu_l];
					}

					if(kk != crl)
					{
						char varNamezyi[100];
						sprintf(varNamezyi, "zi.%d.%d_y.%d.%d", (int) i, (int) jbu, (int) kk, (int) crl); 
						zyi[jbu][kk].setName(varNamezyi);

						temp1+= -1*(lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*zyi[jbu][kk];
						tempRHS+= -1*(lvdata[kk][crv].vendortime - lvdata[kk][crv].assemblytime)*precedence[kk][crl];


					}

					


				}//end kk
				
				
				//cout<<jbu_l<<" "<<jbu_v<<" " <<jbu_s<<endl<<"max subs="<<lvdata[jbu_l][jbu_v].maxsublots<<endl;
				
				temp1+= (lvdata[jbu_l][jbu_v].makespansmp[jbu_s] - lvdata[crl][crv].makespansmp[crs] + 
						-lvdata[jbu_l][jbu_v].assemblytime + lvdata[crl][crv].assemblytime)*z[jbu];
				
					tempRHS+= (lvdata[jbu_l][jbu_v].makespansmp[jbu_s] - lvdata[crl][crv].makespansmp[crs] + 
						-lvdata[jbu_l][jbu_v].assemblytime 
						+lvdata[crl][crv].assemblytime);

					

				chat1[jbu] = IloRange(workerEnv,tempRHS,temp1,IloInfinity);	
				char varName[100];
				sprintf(varName, "chat1.%d.%d", (int) i, (int) jbu); 
				chat1[jbu].setName(varName);

				chat2[jbu] = IloRange(workerEnv,0,temp1,IloInfinity);	
				char varName2[100];
				sprintf(varName2, "chat2.%d.%d", (int) i, (int) jbu); 
				chat2[jbu].setName(varName2);

				workerModel.add(chat1[jbu]);
				workerModel.add(chat2[jbu]);

				

				temp1.end();
				
				
			}//endif	

		}//end jbu
	
	//theta
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
	{			
		if(jbu != i)
		{
			for(unsigned short lk = 0; lk < TLots; ++ lk)
			{	
				unsigned short jbu_l = submap[jbu].parentlot;
				if(lk != jbu_l)
				{
					

					IloExpr temp1(workerEnv);
					IloExpr temp2(workerEnv);
					IloExpr temp3(workerEnv);

					temp1+= -zyj[jbu][lk];
					temp2+=  -zyj[jbu][lk] + z[jbu];
					temp3+= zyj[jbu][lk] - z[jbu];


					the1[jbu][lk] = IloRange(workerEnv,-1*precedence[lk][jbu_l],temp1,IloInfinity);
					the2[jbu][lk] = IloRange(workerEnv,0,temp2,IloInfinity);
					the3[jbu][lk] = IloRange(workerEnv,-1*precedence[jbu_l][lk],temp3,IloInfinity);

					workerModel.add(the1[jbu][lk]);
					workerModel.add(the2[jbu][lk]);
					workerModel.add(the3[jbu][lk]);

					temp1.end();
					temp2.end();
					temp3.end();
				}}

		}}//end theta

	//phi
	for(unsigned short jbu = 0; jbu < numSublots; ++jbu)
	{			
		if(jbu != i)
		{
			char varName[100];
				sprintf(varName, "z.%d.%d", (int) i, (int) jbu); 
				z[jbu].setName(varName);

			for(unsigned short lk = 0; lk < TLots; ++ lk)
			{	
				unsigned short cbu_j = submap[jbu].parentlot;
				if(lk != crl)
				{
					

					IloExpr temp1(workerEnv);
					IloExpr temp2(workerEnv);
					IloExpr temp3(workerEnv);

					temp1+=  -zyi[jbu][lk];
					temp2+=  -zyi[jbu][lk]  +z[jbu];
					temp3+=   zyi[jbu][lk] - z[jbu];


					phi1[jbu][lk] = IloRange(workerEnv,-1*precedence[lk][crl],temp1,IloInfinity);
					phi2[jbu][lk] = IloRange(workerEnv,0,temp2,IloInfinity);
					phi3[jbu][lk] = IloRange(workerEnv,-1*precedence[crl][lk],temp3,IloInfinity);

					/*if(i == 0 && jbu == 1)
				{
					cout<<"phi 1 "<<phi1[i][jbu][lk].getExpr()<<">="<<phi1[i][jbu][lk].getLB()<<endl;
					cout<<"phi 2 "<<phi2[i][jbu][lk].getExpr()<<">="<<phi2[i][jbu][lk].getLB()<<endl;
					cout<<"phi 3 "<<phi3[i][jbu][lk].getExpr()<<">="<<phi3[i][jbu][lk].getLB()<<endl;
					}*/

					workerModel.add(phi1[jbu][lk]);
					workerModel.add(phi2[jbu][lk]);
					workerModel.add(phi3[jbu][lk]);

					temp1.end();
					temp2.end();
					temp3.end();
				}}

		}}//end phi

	}// end of big subs loop

	cout<<"Added all constraints"<<endl;


	IloRange gamma(workerEnv,1,gammaex,IloInfinity,"gamma");
	workerModel.add(gamma);
	
	workerModel.add(IloMinimize(workerEnv,workerObj));

	IloCplex cplex(workerModel);
	cplex.setParam(IloCplex::PreDual,1);
	cplex.setVectors(vals, 0, vars,0,0,0);
	cplex.setParam(IloCplex::PreInd, false);
	cplex.setParam(IloCplex::WorkMem, MaxMemory);
	//cplex.setOut(workerEnv.getNullStream());

	cout<<"Attempting worker solution"<<endl;
	if(cplex.solve())
	{
		double y_factor=0;
		cout<<"Solved worker cplex"<<endl;

		cut_yij_coeff.clear();		
		cut_yij_coeff.resize(TLots);		

		for(unsigned short i=0; i<TLots; ++i)
		{
			cut_yij_coeff[i].resize(TLots);
			fill(cut_yij_coeff[i].begin(), cut_yij_coeff[i].end(), 0);
			
		}

		

		unsigned short icr = submap_inv[crLot][crVendor][crSublot];
		cout<<"trying to add pi duals"<<endl;
		//add pi
		for(unsigned short k=0; k<TLots; ++k)
		{
			if(k!=crLot)
			{
				double temp= cplex.getDual(pi2[k]);
				cut_yij_coeff[k][crLot]+= -1*temp;

				temp= cplex.getDual(pi3[k]);
				cut_yij_coeff[crLot][k]= -1*temp;
			}
		}

		
		cout<<"trying to add chat duals"<<endl;
		//add chat1
		for(unsigned short jbu=0; jbu<numSublots; ++jbu)
		{
			unsigned short jbu_l = submap[jbu].parentlot;
			unsigned short jbu_v = submap[jbu].parentvendor;

			if(jbu != icr)
			{
				double chat1_jbu_dual = cplex.getDual(chat1[jbu]);
				for(unsigned short k=0; k<TLots; ++k)
				{
					if(jbu_l != k)
					{
						double temp = (lvdata[k][jbu_v].vendortime - lvdata[k][jbu_v].assemblytime);
						cut_yij_coeff[k][jbu_l]+= temp*chat1_jbu_dual;
						if(abs(chat1_jbu_dual)>0)
						{
						cout<<"y co-effs during chat1 first loop \n";
						cout<<"jbu_l="<<jbu_l<<" k="<<k
							<<" temp="<<temp<<" dual="<<chat1_jbu_dual<<endl;
		for(unsigned short i=0; i<TLots; ++i)
		{	
			
			for(unsigned short j=0; j<TLots; ++j)
			cout<<cut_yij_coeff[i][j]<<" ";
			cout<<endl;
		}
		}
					}

					if(k != crLot)
					{
						cut_yij_coeff[k][crLot] += -1*chat1_jbu_dual*(lvdata[k][crVendor].vendortime - lvdata[k][crVendor].assemblytime);
					}
				}
				

			}
		}
		
		cout<<"trying to add phi and theta duals"<<endl;
		//add theta and phi
		for(unsigned short jbu=0; jbu<numSublots; ++jbu)
		{
			unsigned short jbu_l = submap[jbu].parentlot;
			unsigned short jbu_v = submap[jbu].parentvendor;

			if(jbu != icr)
			{
				
				for(unsigned short k=0; k<TLots; ++k)
				{
					if(jbu_l != k)
					{
						double temp = -1* cplex.getDual(the1[jbu][k]);
						cut_yij_coeff[k][jbu_l]+= temp;

						temp = -1* cplex.getDual(the3[jbu][k]);
						cut_yij_coeff[jbu_l][k] += temp;
					}

					if(k != crLot)
					{
						double temp = -1* cplex.getDual(phi1[jbu][k]);
						cut_yij_coeff[k][crLot]+= temp;

						temp = -1* cplex.getDual(phi3[jbu][k]);
						cut_yij_coeff[crLot][k] += temp;
					}
				}
				

			}
		}
		
		
		//create cut LHS
		cout<<"y co-effs \n";

		cutfile<< "for sequence: \n";
		for(unsigned short i=0; i<TLots; ++i)
		{
			
			for(unsigned short j=0; j<TLots; ++j)
			{
				
				cutfile<<precedence[i][j]<<" ";
			}
			cutfile<<endl;
		}//add y factors
		cutfile<<"cut added: \n";
		for(unsigned short i=0; i<TLots; ++i)
		{
			
			for(unsigned short j=0; j<TLots; ++j)
			{
				cout<<cut_yij_coeff[i][j]<<" ";
				cutLhs+= cut_yij_coeff[i][j]*y[i][j];
				y_factor+= cut_yij_coeff[i][j]*precedence[i][j];
				if(cut_yij_coeff[i][j] != 0)
					cutfile<<endl<<cut_yij_coeff[i][j]<<"*y["<<i<<"]["<<j<<"] ";
			}
			cout<<endl;
		}//add y factors
		cout<<"y factor="<<y_factor<<endl;
		cout<<"sub objective="<<cplex.getObjValue()<<endl;
		cutLhs+= cplex.getObjValue()-y_factor + ObjConstant;
		cutfile<<" + "<<cplex.getObjValue()-y_factor<<" <= z_master \n";
		cout<<cutLhs<<endl;
	}

	

	cout<<"here\n";
	//workerModel.end();
	double total_obj = cplex.getObjValue()+ObjConstant;
	workerEnv.end();
	cout<<"hear\n";

	datafile.close();
	cutfile.close();
	return total_obj;
	}

	
double datalap::optimalNetworkBendersallycuts(void)
{
	datalap::master_instance = this;
	IloEnv env;
	IloModel model(env);
	IloNumVar obj(env,0);
	IloExpr objective(env);	
	


	BoolVarArray2 y(env, TLots);
	for(unsigned short int i=0;i<TLots;i++)
		y[i] =BoolVarArray1(env, TLots);



	IloConstraintArray constraints(env);


	//create transitive constraints
	//y[i][j]+y[j][k]+y[k][i]<=2 i!=j!=k!=i (implemented here as i>j>k)
	for(unsigned int short i=0;i<TLots;i++)
	{
		for(unsigned int short j=0;j<TLots;j++)
		{
			for(unsigned int short k=0;k<TLots;k++)
			{					 
				if(i!=j && j!=k && k!=i)
					constraints.add(y[i][j]+y[j][k]+y[k][i]<=2);
			}}}
	for(unsigned int short i=0;i<TLots;i++)
	{
		for(unsigned int short j=0;j<i;j++)
		{

			constraints.add(y[i][j]+y[j][i]==1);
		}}

	for(unsigned short int i=0;i<TLots;i++)
	{
		constraints.add(y[i][i]==0);
	}







	//objecctive funtion expression
	objective+=obj;

	

	
		/*try
		{
		separatepermx(cutLHS_temp, x);
		}
		catch(...)
		{
			cout<<"could not generate cut";
		}*/
		
	//add constraints to model and solve
	model.add(constraints);
	model.add(IloMinimize(env,obj));
	IloCplex cplex(model);
	//cplex.setOut(env.getNullStream());

	cplex.setParam(IloCplex::WorkMem, MaxMemory);
	cplex.setParam(IloCplex::NodeFileInd, 3);
	cplex.setParam(IloCplex::TreLim,MaxMemory);
	cplex.setParam(IloCplex::TiLim,MaxTime);
	cplex.setParam(IloCplex::PreInd, IloFalse); 
		cplex.setParam(IloCplex::Threads, 1); 
		cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
		cplex.use(BendersLazyCallbacky(env,y,obj,ObjConstant));
		//cplex.use(BendersUserCallback(env,y,x,obj,ObjConstant));
	IloTimer cputime(env);
	cputime.start();
	//cplex.setOut(env.getNullStream());
	//cplex.setWarning(env.getNullStream());

	if(cplex.solve())
	{
		if (cplex.getStatus()==IloAlgorithm::Optimal)
		{

			cout<<"\n solved"<<endl;

			m_cputime=cputime.stop();
			m_optimalobjective=cplex.getObjValue();
			cout<<"optimal objective="<<cplex.getObjValue()<<" CPU time="<<m_cputime<<endl;
		}	
		else
		{
			m_cputime=0;
			m_optimalobjective=0;
		}
	}
	else
	{
		m_cputime=0;
		m_optimalobjective=0;
		cout<<"Something went wrong."<<endl;
		cputime.stop();

	}


	env.end();
	return m_optimalobjective;

}

datalap * datalap::master_instance= NULL;



double datalap::optimalNetworkBenders1(void)
{	
	datalap::master_instance = this;
	IloEnv env;
	IloModel model(env);
	double objectivevaue;


	//Building model by rows//
	//create the 2 dimensional matrix of variables y[i][j]

	IloNumVar obj(env);
	IloNumVar Cmax(env);
	IloExpr objective(env);	 

	FloatVarArray2 y(env,TLots);
	for(unsigned short int i=0;i<TLots;i++)
		y[i] = FloatVarArray1(env, TLots,0,1);

	BoolVarArray2 x(env, TLots);
	for(unsigned short int i=0;i<TLots;i++)
		x[i] =BoolVarArray1(env, TLots);

	FloatVarArray4 zsublots(env,TLots);
	for(unsigned short int i=0;i<TLots;i++)
		zsublots[i]=FloatVarArray3(env,NumVendors);
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			int k=lvdata[i][j].maxsublots;
			zsublots[i][j]=FloatVarArray2(env,k);
		}
	}
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				zsublots[i][j][k]=FloatVarArray1(env,TLots,0,1);
			}
		}
	}

	

	IloConstraintArray constraints(env);


	//create transitive constraints
	//y[i][j]+y[j][k]+y[k][i]<=2 i!=j!=k!=i (implemented here as i>j>k)
	for(unsigned int short i=0;i<TLots;i++)
	{
		IloExpr sumofrow(env),sumofcolumn(env);
		for(unsigned int short j=0;j<TLots;j++)
		{
			sumofrow+=x[i][j];
			sumofcolumn+=x[j][i];

		}
		constraints.add(sumofcolumn==1);
		constraints.add(sumofrow==1);

	}

	IloExpr sumy(env);
	for(unsigned int short i=0;i<TLots;i++)
	{
		
		for(unsigned int short j=0;j<TLots;j++)
		{
			if(i!=j)
			{
				sumy+=y[i][j];
				for(unsigned int short k=0;k<TLots;k++)
				{
					IloExpr sety(env);
					for(unsigned int short kk=k+1;kk<TLots;kk++)
					{
						if(kk<TLots)
							sety+= x[j][kk];

					}

					sety += x[i][k]-1-y[i][j];
					constraints.add(sety <= 0);
					sety.end();

				}
			}

		}
	}
	constraints.add(sumy== 0.5*TLots*(TLots-1));


	//objecctive funtion expression
	objective+=Cmax;

	//create  constraints



	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			IloExpr sublotssumtoone(env);
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				for(unsigned short int l=0;l<TLots;l++)
				{
					sublotssumtoone+=zsublots[i][j][k][l];		
					objective+=lvdata[i][j].handlingcost[k]*zsublots[i][j][k][l];//add to objective
					constraints.add(zsublots[i][j][k][l]<=x[i][l]);//constraints saying sublot assignment and xij are linked

				}				 

			}//made an expression summing subots
			constraints.add(sublotssumtoone==1);
			sublotssumtoone.end();

		}
	}//end i j loop


	for(unsigned short int l=0;l<TLots;l++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			IloExpr assosciatedmakespan(env);

			for(unsigned short int lotto=0;lotto<TLots;lotto++)
				for(unsigned short int i=0;i<l;i++)
					assosciatedmakespan+=lvdata[lotto][j].vendortime*x[lotto][i];

			for(unsigned short int lotto=0;lotto<TLots;lotto++)
				for(unsigned short int i=l+1;i<TLots;i++)
					assosciatedmakespan+=lvdata[lotto][j].assemblytime*x[lotto][i];

			for(unsigned short int i=0;i<TLots;i++)
			{
				for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
				{
					assosciatedmakespan+=lvdata[i][j].makespansmp[k]*zsublots[i][j][k][l];

				}
			}

			constraints.add(assosciatedmakespan<=Cmax);
			assosciatedmakespan.end();

		}//next vendor
	}//next final destination location



	constraints.add(obj >= objective);


	//add constraints to model and solve
	model.add(constraints);
	model.add(IloMinimize(env,obj));
	IloCplex cplex(model);
	//cplex.setOut(env.getNullStream());

	///////////////////////
	//set initial values
	//////////////////////

	lpNB();
	clean_perm(m_perm);
	IloNumVarArray startVar(env);
	IloNumArray startVal(env);
	for (int i = 0; i<TLots;++i){
	for (int j = 0; j<TLots;j++) {
	startVar.add(x[i][j]);
	startVal.add(m_perm[i][j]);
	}
	}
	cplex.addMIPStart(startVar, startVal, IloCplex::MIPStartAuto,  "secondMIPStart");
	startVal.end();
	startVar.end();
	

	////////////////////////
	//finished setting initial values
	////////////////////////
	cplex.setParam(IloCplex::WorkMem, MaxMemory);
	// cplex.setParam(IloCplex::NodeFileInd, 3);
	cplex.setParam(IloCplex::TreLim,MaxMemory);
	cplex.setParam(IloCplex::TiLim,MaxTime);
	cplex.setParam(IloCplex::PreInd, IloFalse); 
		cplex.setParam(IloCplex::Threads, 1); 
		cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
		//cplex.use(BendersLazyCallback(env,y,x,obj,ObjConstant));
		//cplex.use(BendersUserCallback(env,y,x,obj,ObjConstant));
	IloTimer cputime(env);
	cputime.start();
	//cplex.setOut(env.getNullStream());
	//cplex.setWarning(env.getNullStream());

	if(cplex.solve())
	{
		if (cplex.getStatus()==IloAlgorithm::Optimal)
		{

			cout<<"\n solved"<<endl;

			m_cputime=cputime.stop();
			m_optimalobjective=cplex.getObjValue();
			cout<<"optimal objective="<<cplex.getObjValue()<<" CPU time="<<m_cputime<<" Makespan="<<cplex.getValue(Cmax)<<endl;
		}	
		else
		{
			m_cputime=0;
			m_optimalobjective=0;
		}
	}
	else
	{
		m_cputime=0;
		m_optimalobjective=0;
		cout<<"Something went wrong."<<endl;
		cputime.stop();

	}


	env.end();
	return m_optimalobjective;

}

double datalap::optimalNetworkBenders2(void)
{	
	datalap::master_instance = this;
	IloEnv env;
	IloModel model(env);
	double objectivevaue;


	//Building model by rows//
	//create the 2 dimensional matrix of variables y[i][j]

	IloNumVar obj(env);
	IloNumVar Cmax(env);
	IloExpr objective(env);	 

	FloatVarArray2 y(env,TLots);
	for(unsigned short int i=0;i<TLots;i++)
		y[i] = FloatVarArray1(env, TLots,0,1);

	BoolVarArray2 x(env, TLots);
	for(unsigned short int i=0;i<TLots;i++)
		x[i] =BoolVarArray1(env, TLots);

	BoolVarArray4 zsublots(env,TLots);
	for(unsigned short int i=0;i<TLots;i++)
		zsublots[i]=BoolVarArray3(env,NumVendors);
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			int k=lvdata[i][j].maxsublots;
			zsublots[i][j]=BoolVarArray2(env,k);
		}
	}
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				zsublots[i][j][k]=BoolVarArray1(env,TLots);
			}
		}
	}

	

	IloConstraintArray constraints(env);


	//create transitive constraints
	//y[i][j]+y[j][k]+y[k][i]<=2 i!=j!=k!=i (implemented here as i>j>k)
	for(unsigned int short i=0;i<TLots;i++)
	{
		IloExpr sumofrow(env),sumofcolumn(env);
		for(unsigned int short j=0;j<TLots;j++)
		{
			sumofrow+=x[i][j];
			sumofcolumn+=x[j][i];

		}
		constraints.add(sumofcolumn==1);
		constraints.add(sumofrow==1);

	}

	IloExpr sumy(env);
	for(unsigned int short i=0;i<TLots;i++)
	{
		
		for(unsigned int short j=0;j<TLots;j++)
		{
			if(i!=j)
			{
				sumy+=y[i][j];
				for(unsigned int short k=0;k<TLots;k++)
				{
					IloExpr sety(env);
					for(unsigned int short kk=k+1;kk<TLots;kk++)
					{
						if(kk<TLots)
							sety+= x[j][kk];

					}

					sety += x[i][k]-1-y[i][j];
					constraints.add(sety <= 0);
					sety.end();

				}
			}

		}
	}
	constraints.add(sumy== 0.5*TLots*(TLots-1));


	//objecctive funtion expression
	objective+=Cmax;

	//create  constraints



	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			IloExpr sublotssumtoone(env);
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				for(unsigned short int l=0;l<TLots;l++)
				{
					sublotssumtoone+=zsublots[i][j][k][l];		
					objective+=lvdata[i][j].handlingcost[k]*zsublots[i][j][k][l];//add to objective
					constraints.add(zsublots[i][j][k][l]<=x[i][l]);//constraints saying sublot assignment and xij are linked

				}				 

			}//made an expression summing subots
			constraints.add(sublotssumtoone==1);
			sublotssumtoone.end();

		}
	}//end i j loop


	for(unsigned short int l=0;l<TLots;l++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			IloExpr assosciatedmakespan(env);

			for(unsigned short int lotto=0;lotto<TLots;lotto++)
				for(unsigned short int i=0;i<l;i++)
					assosciatedmakespan+=lvdata[lotto][j].vendortime*x[lotto][i];

			for(unsigned short int lotto=0;lotto<TLots;lotto++)
				for(unsigned short int i=l+1;i<TLots;i++)
					assosciatedmakespan+=lvdata[lotto][j].assemblytime*x[lotto][i];

			for(unsigned short int i=0;i<TLots;i++)
			{
				for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
				{
					assosciatedmakespan+=lvdata[i][j].makespansmp[k]*zsublots[i][j][k][l];

				}
			}

			constraints.add(assosciatedmakespan<=Cmax);
			assosciatedmakespan.end();

		}//next vendor
	}//next final destination location



	constraints.add(obj >= objective);

	//START PREPROCESSING
	///////////////////////
	//set initial values
	//////////////////////

	lpNB();
	vector<vector<double>> temp_perm = m_perm;
	//add bender's cuts
	bool addbenders = true;
	while(addbenders)
	{
		clean_perm(temp_perm);
		IloExpr cutLHS_temp(env);
		separateperm(cutLHS_temp, y);
		constraints.add(obj >= cutLHS_temp);
		//find min co-eff of m_perm 
		double min_coeff = 0;
		for(unsigned short int i=0;i<TLots;i++)
		{
			for(unsigned short int j=0;j<TLots;j++)
			{
				if(m_perm[i][j] > 0.99 && min_coeff< temp_perm[i][j])
					min_coeff = temp_perm[i][j];
			}}

		if(min_coeff < 0.05)
			addbenders  = false;
		else
		{
			for(unsigned short int i=0;i<TLots;i++)
			{
				for(unsigned short int j=0;j<TLots;j++)
				{
					if(m_perm[i][j] > 0.99 )
						temp_perm[i][j] = temp_perm[i][j]-min_coeff;
				}}
		}
	}
	//end bender's cuts
	IloNumVarArray startVar(env);
	IloNumArray startVal(env);
	for (int i = 0; i<TLots;++i){
	for (int j = 0; j<TLots;j++) {
	startVar.add(x[i][j]);
	startVal.add(m_perm[i][j]);
	}
	}
	
	//END PREPROCESSING


	//add constraints to model and solve
	model.add(constraints);
	model.add(IloMinimize(env,obj));
	IloCplex cplex(model);
	//cplex.setOut(env.getNullStream());

	
	cplex.addMIPStart(startVar, startVal, IloCplex::MIPStartAuto,  "secondMIPStart");
	startVal.end();
	startVar.end();

	////////////////////////
	//finished setting initial values
	////////////////////////
	cplex.setParam(IloCplex::WorkMem, MaxMemory);
	// cplex.setParam(IloCplex::NodeFileInd, 3);
	cplex.setParam(IloCplex::TreLim,MaxMemory);
	cplex.setParam(IloCplex::TiLim,MaxTime);
	//cplex.setParam(IloCplex::PreInd, IloFalse); 
		//cplex.setParam(IloCplex::Threads, 1); 
		//cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
		//cplex.use(BendersLazyCallback(env,y,x,obj,ObjConstant));
		//cplex.use(BendersUserCallback(env,y,x,obj,ObjConstant));
	IloTimer cputime(env);
	cputime.start();
	//cplex.setOut(env.getNullStream());
	//cplex.setWarning(env.getNullStream());

	if(cplex.solve())
	{
		if (cplex.getStatus()==IloAlgorithm::Optimal)
		{

			cout<<"\n solved"<<endl;

			m_cputime=cputime.stop();
			m_optimalobjective=cplex.getObjValue();
			cout<<"optimal objective="<<cplex.getObjValue()<<" CPU time="<<m_cputime<<" Makespan="<<cplex.getValue(Cmax)<<endl;
		}	
		else
		{
			m_cputime=0;
			m_optimalobjective=0;
		}
	}
	else
	{
		m_cputime=0;
		m_optimalobjective=0;
		cout<<"Something went wrong."<<endl;
		cputime.stop();

	}


	env.end();
	return m_optimalobjective;

}

double datalap::optimalNetworkBenders3(void)
{	
	datalap::master_instance = this;
	IloEnv env;
	IloModel model(env);
	double objectivevaue;


	//Building model by rows//
	//create the 2 dimensional matrix of variables y[i][j]

	IloNumVar obj(env);
	IloNumVar Cmax(env);
	IloExpr objective(env);	 

	/*FloatVarArray2 y(env,TLots);
	for(unsigned short int i=0;i<TLots;i++)
		y[i] = FloatVarArray1(env, TLots,0,1);*/

	BoolVarArray2 x(env, TLots);
	for(unsigned short int i=0;i<TLots;i++)
		x[i] =BoolVarArray1(env, TLots);

	BoolVarArray4 zsublots(env,TLots);
	for(unsigned short int i=0;i<TLots;i++)
		zsublots[i]=BoolVarArray3(env,NumVendors);
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			int k=lvdata[i][j].maxsublots;
			zsublots[i][j]=BoolVarArray2(env,k);
		}
	}
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				zsublots[i][j][k]=BoolVarArray1(env,TLots);
			}
		}
	}

	

	IloConstraintArray constraints(env);


	//create transitive constraints
	//y[i][j]+y[j][k]+y[k][i]<=2 i!=j!=k!=i (implemented here as i>j>k)
	for(unsigned int short i=0;i<TLots;i++)
	{
		IloExpr sumofrow(env),sumofcolumn(env);
		for(unsigned int short j=0;j<TLots;j++)
		{
			sumofrow+=x[i][j];
			sumofcolumn+=x[j][i];

		}
		constraints.add(sumofcolumn==1);
		constraints.add(sumofrow==1);

	}

	//IloExpr sumy(env);
	/*for(unsigned int short i=0;i<TLots;i++)
	{
		
		for(unsigned int short j=0;j<TLots;j++)
		{
			if(i!=j)
			{
				sumy+=y[i][j];
				for(unsigned int short k=0;k<TLots;k++)
				{
					IloExpr sety(env);
					for(unsigned int short kk=k+1;kk<TLots;kk++)
					{
						if(kk<TLots)
							sety+= x[j][kk];

					}

					sety += x[i][k]-1-y[i][j];
					constraints.add(sety <= 0);
					sety.end();

				}
			}

		}
	}
	constraints.add(sumy== 0.5*TLots*(TLots-1));*/


	//objecctive funtion expression
	objective+=Cmax;

	//create  constraints



	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			IloExpr sublotssumtoone(env);
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				for(unsigned short int l=0;l<TLots;l++)
				{
					sublotssumtoone+=zsublots[i][j][k][l];		
					objective+=lvdata[i][j].handlingcost[k]*zsublots[i][j][k][l];//add to objective
					constraints.add(zsublots[i][j][k][l]<=x[i][l]);//constraints saying sublot assignment and xij are linked

				}				 

			}//made an expression summing subots
			constraints.add(sublotssumtoone==1);
			sublotssumtoone.end();

		}
	}//end i j loop


	for(unsigned short int l=0;l<TLots;l++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			IloExpr assosciatedmakespan(env);

			for(unsigned short int lotto=0;lotto<TLots;lotto++)
				for(unsigned short int i=0;i<l;i++)
					assosciatedmakespan+=lvdata[lotto][j].vendortime*x[lotto][i];

			for(unsigned short int lotto=0;lotto<TLots;lotto++)
				for(unsigned short int i=l+1;i<TLots;i++)
					assosciatedmakespan+=lvdata[lotto][j].assemblytime*x[lotto][i];

			for(unsigned short int i=0;i<TLots;i++)
			{
				for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
				{
					assosciatedmakespan+=lvdata[i][j].makespansmp[k]*zsublots[i][j][k][l];

				}
			}

			constraints.add(assosciatedmakespan<=Cmax);
			assosciatedmakespan.end();

		}//next vendor
	}//next final destination location



	constraints.add(obj >= objective);

	//START PREPROCESSING
	///////////////////////
	//set initial values
	//////////////////////

	lpNB();
	vector<vector<double>> temp_perm = m_perm;
	//add bender's cuts
	bool addbenders = true;
	while(addbenders)
	{
		clean_perm(temp_perm);
		IloExpr cutLHS_temp(env);
		try
		{
		separatepermx(cutLHS_temp, x);
		}
		catch(...)
		{
			cout<<"could not generate cut";
		}
		cutLHS_temp.normalize();
		constraints.add(obj >= cutLHS_temp);
		cutLHS_temp.end();
		//find min co-eff of m_perm 
		double min_coeff = 0;
		for(unsigned short int i=0;i<TLots;i++)
		{
			for(unsigned short int j=0;j<TLots;j++)
			{
				if(m_perm[i][j] > 0.99 && min_coeff< temp_perm[i][j])
					min_coeff = temp_perm[i][j];
			}}

		if(min_coeff < 0.05)
			addbenders  = false;
		else
		{
			for(unsigned short int i=0;i<TLots;i++)
			{
				for(unsigned short int j=0;j<TLots;j++)
				{
					if(m_perm[i][j] > 0.99 )
						temp_perm[i][j] = temp_perm[i][j]-min_coeff;
				}}
		}
	}
	//end bender's cuts
	IloNumVarArray startVar(env);
	IloNumArray startVal(env);
	for (int i = 0; i<TLots;++i){
	for (int j = 0; j<TLots;j++) {
	startVar.add(x[i][j]);
	startVal.add(m_perm[i][j]);
	}
	}
	
	//END PREPROCESSING


	//add constraints to model and solve
	model.add(constraints);
	model.add(IloMinimize(env,obj));
	IloCplex cplex(model);
	//cplex.setOut(env.getNullStream());

	
	cplex.addMIPStart(startVar, startVal, IloCplex::MIPStartAuto,  "secondMIPStart");
	startVal.end();
	startVar.end();

	////////////////////////
	//finished setting initial values
	////////////////////////
	cplex.setParam(IloCplex::WorkMem, MaxMemory);
	// cplex.setParam(IloCplex::NodeFileInd, 3);
	cplex.setParam(IloCplex::TreLim,MaxMemory);
	cplex.setParam(IloCplex::TiLim,MaxTime);
	//cplex.setParam(IloCplex::PreInd, IloFalse); 
		//cplex.setParam(IloCplex::Threads, 1); 
		//cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
		//cplex.use(BendersLazyCallback(env,y,x,obj,ObjConstant));
		//cplex.use(BendersUserCallback(env,y,x,obj,ObjConstant));
	IloTimer cputime(env);
	cputime.start();
	//cplex.setOut(env.getNullStream());
	//cplex.setWarning(env.getNullStream());

	if(cplex.solve())
	{
		if (cplex.getStatus()==IloAlgorithm::Optimal)
		{

			cout<<"\n solved"<<endl;

			m_cputime=cputime.stop();
			m_optimalobjective=cplex.getObjValue();
			cout<<"optimal objective="<<cplex.getObjValue()<<" CPU time="<<m_cputime<<" Makespan="<<cplex.getValue(Cmax)<<endl;
		}	
		else
		{
			m_cputime=0;
			m_optimalobjective=0;
		}
	}
	else
	{
		m_cputime=0;
		m_optimalobjective=0;
		cout<<"Something went wrong."<<endl;
		cputime.stop();

	}


	env.end();
	return m_optimalobjective;

}

double datalap::optimalNetworkBendersallcuts(void)
{
	datalap::master_instance = this;
	IloEnv env;
	IloModel model(env);
	IloNumVar obj(env,0);
	IloExpr objective(env);	
	
	BoolVarArray2 x(env, TLots);
	for(unsigned short int i=0;i<TLots;i++)
		x[i] =BoolVarArray1(env, TLots);

	BoolVarArray2 y(env, TLots);
	for(unsigned short int i=0;i<TLots;i++)
		y[i] =BoolVarArray1(env, TLots);



	IloConstraintArray constraints(env);


	//create transitive constraints
	//y[i][j]+y[j][k]+y[k][i]<=2 i!=j!=k!=i (implemented here as i>j>k)
	for(unsigned int short i=0;i<TLots;i++)
	{
		for(unsigned int short j=0;j<TLots;j++)
		{
			for(unsigned int short k=0;k<TLots;k++)
			{					 
				if(i!=j && j!=k && k!=i)
					constraints.add(y[i][j]+y[j][k]+y[k][i]<=2);
			}}}
	for(unsigned int short i=0;i<TLots;i++)
	{
		for(unsigned int short j=0;j<i;j++)
		{

			constraints.add(y[i][j]+y[j][i]==1);
		}}

	for(unsigned short int i=0;i<TLots;i++)
	{
		constraints.add(y[i][i]==0);
	}


		for(unsigned int short i=0;i<TLots;i++)
	{
		IloExpr sumofrow(env),sumofcolumn(env);
		for(unsigned int short j=0;j<TLots;j++)
		{
			sumofrow+=x[i][j];
			sumofcolumn+=x[j][i];

		}
		constraints.add(sumofcolumn==1);
		constraints.add(sumofrow==1);

	}


		//setx
		for(unsigned int short i=0;i<TLots;i++)
	{
		IloExpr setx(env);
		for(unsigned int short j=0;j<TLots;j++)
		{
			setx+=j*x[i][j]-y[j][i];
			

		}
		constraints.add(setx==0);
	}

	//objecctive funtion expression
	objective+=obj;

	

	
		/*try
		{
		separatepermx(cutLHS_temp, x);
		}
		catch(...)
		{
			cout<<"could not generate cut";
		}*/
		
	//add constraints to model and solve
	model.add(constraints);
	model.add(IloMinimize(env,obj));
	IloCplex cplex(model);
	//cplex.setOut(env.getNullStream());

	

	cplex.setParam(IloCplex::WorkMem, MaxMemory);
	 cplex.setParam(IloCplex::NodeFileInd, 3);
	cplex.setParam(IloCplex::TreLim,MaxMemory);
	cplex.setParam(IloCplex::TiLim,MaxTime);
	cplex.setParam(IloCplex::PreInd, IloFalse); 
		cplex.setParam(IloCplex::Threads, 1); 
		cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
		cplex.use(BendersLazyCallback(env,y,x,obj,ObjConstant));
		//cplex.use(BendersUserCallback(env,y,x,obj,ObjConstant));
	IloTimer cputime(env);
	cputime.start();
	//cplex.setOut(env.getNullStream());
	//cplex.setWarning(env.getNullStream());

	if(cplex.solve())
	{
		if (cplex.getStatus()==IloAlgorithm::Optimal)
		{

			cout<<"\n solved"<<endl;

			m_cputime=cputime.stop();
			m_optimalobjective=cplex.getObjValue();
			cout<<"optimal objective="<<cplex.getObjValue()<<" CPU time="<<m_cputime<<endl;
		}	
		else
		{
			m_cputime=0;
			m_optimalobjective=0;
		}
	}
	else
	{
		m_cputime=0;
		m_optimalobjective=0;
		cout<<"Something went wrong."<<endl;
		cputime.stop();

	}


	env.end();
	return m_optimalobjective;

}

void datalap::perm_to_prec(void)
{
	vector<unsigned short> position;
	position.resize(TLots);

	for(unsigned short i=0;i<TLots;++i)
	{
		position[i]=0;
		for(unsigned short j=0;j<TLots;++j)
		{
			position[i]+= m_perm[i][j]*j;

		}
		//cout<<position[i]<<endl;

	}

	for(unsigned short i=0;i<TLots;++i)
	{
		
		for(unsigned short j=0;j<TLots;++j)
		{
			if(i==j)
				precedence[i][j]=0;
			if(position[i]<position[j])
				precedence[i][j]=1;
			else
				precedence[i][j]=0;
			//cout<<precedence[i][j]<<" ";

		}	
		//cout<<endl;

	}


}



void datalap::generatedata(int root, int Lots, int Vendors, int ubsublots)
{
	MaxMemory = 1400;
	MaxTime= 1800;
	
	//First clear any old data
	lvdata.clear();

	//set upper limits on parameter values
	srand(root);
	//ULlots=40;
	ULvendors=10;
	ULmaxsublots=ubsublots;
	ULprocesstime=20; 
	ULtransfercost=ULprocesstime/ULmaxsublots;


	//set parameter values
	TLots=Lots;
	NumVendors=Vendors;

	lvdata.resize(TLots);
	precedence.resize(TLots);
	for(unsigned short int i=0; i<lvdata.size();i++)
	{
		lvdata[i].resize(NumVendors);

	}

	///set default precedence to 1,2,3...
	for(unsigned short int i=0; i<TLots;i++)
	{
		for(unsigned short int j=0; j<TLots;j++)
		{
			if(i<j)
				precedence[i].push_back(1);
			else
				precedence[i].push_back(0);
		}

	}



	//create vectors to store the values created
	vector< vector<int>> maximum_sublots, vendor_processing_time;
	vector<int> assembly_time;

	//output it to a stream that goes to a .dat file
	std::ofstream datafile;
	datafile.open("test.dat");

	datafile<< "\n TLots="<<TLots<<";";
	datafile<<"\n TVendors="<<NumVendors<<";";
	datafile<<"\n TMaxSublots="<<ULmaxsublots<<";";

	//Insert Maxsublots
	datafile<< "\n MaxSublots=[";
	//Each lot has a new row
	for(int i=0;i<TLots;i++)
	{
		vector<int> maxsublotlot;
		maxsublotlot.clear();
		datafile<<"\n [";
		//add a variable for each vendor
		for(int j=0;j<NumVendors;j++)
		{
			int maxsub= rand()%ULmaxsublots+1;
			datafile<<" "<<maxsub;
			maxsublotlot.push_back(maxsub);
			lvdata[i][j].maxsublots=maxsub;
			if(NumVendors-1-j)
				datafile<<",";
		}
		maximum_sublots.push_back(maxsublotlot);
		maxsublotlot.clear();
		datafile<<"],";

	}
	datafile<< "\n ]; \n";
	////////////////finished creating MaxSublots matrix

	//Insert Vendor Processing Times
	vendor_processing_time.clear();
	datafile<< "\n p=[";
	//Each lot has a new row
	for(int i=0;i<TLots;i++)
	{
		datafile<<"\n [";
		vector<int> lotptime;
		//add a variable for each vendor
		lotptime.clear();
		for(int j=0;j<NumVendors;j++)
		{

			int ptemp= rand()%ULprocesstime+1; //upperbound used here!
			lotptime.push_back(ptemp);
			datafile<<" "<<ptemp;
			lvdata[i][j].vendortime=ptemp;
			if(NumVendors-1-j)
				datafile<<",";
		}
		vendor_processing_time.push_back(lotptime);
		lotptime.clear();
		datafile<<"],";

	}
	datafile<< "\n ]; \n";
	////////////////finished creating Processingtimes matrix matrix



	//Insert assembly times
	//Insert Vendor Processing Times
	assembly_time.clear();
	datafile<< "\n pass=[";
	//Each lot has a new row
	for(int i=0;i<TLots;i++)
	{
		int patemp= rand()%ULprocesstime+1; //upperbound used here!
		datafile<<" "<<patemp;
		for(unsigned short int j=0; j<NumVendors;j++)
			lvdata[i][j].assemblytime=patemp;
		assembly_time.push_back(patemp);
		if(TLots-1)
			datafile<<",";

	}
	datafile<< "]; \n";
	////////////////finished creating Processingtimes matrix matrix

	////Writing min makespan
	datafile<<"\n MinMakespan=[";

	//Each lot has a new row
	for(int i=0;i<TLots;i++)
	{
		datafile<<"\n [";
		//add a variable for each vendor

		for(int j=0;j<NumVendors;j++)
		{
			datafile<<"[";
			int sublot=0;
			int tempmake=0;


			while (sublot<ULmaxsublots)
			{
				sublot++;

				if(sublot==1)
				{ 
					tempmake=vendor_processing_time[i][j];
					lvdata[i][j].makespansmp.push_back(tempmake+assembly_time[i]);
				}
				else
					if(sublot<=maximum_sublots[i][j])
					{
						tempmake*=(0.5+((rand()%5)*0.10));
						if(!tempmake)
							tempmake=1;
						lvdata[i][j].makespansmp.push_back(tempmake+assembly_time[i]);
					}
					else
						tempmake=0;


				datafile<<" "<<tempmake+assembly_time[i];
				if(ULmaxsublots>sublot)
					datafile<<",";
			}

			datafile<<"]";
			if(NumVendors-j-1)
				datafile<<",";
		}

		//at the end of each lot close a square bracket
		datafile<<"]";
		if(TLots-i-1)
			datafile<<",";


	}
	datafile<< "\n ]; \n";
	//Finished min-makespan matrix



	//Insert transfer cost
	datafile<<"\n TransferCost=[";

	//Each lot has a new row
	for(int i=0;i<TLots;i++)
	{
		datafile<<"\n [";
		//add a variable for each vendor

		for(int j=0;j<NumVendors;j++)
		{
			datafile<<"[";
			int sublot=0;
			int tempcost=0;
			int singletransfer;
			singletransfer=1+rand()%ULtransfercost;


			while (sublot<ULmaxsublots)
			{
				sublot++;


				if(sublot==1)
				{ 

					tempcost=singletransfer;
					lvdata[i][j].handlingcost.push_back(tempcost);
				}
				else
					if(sublot<=maximum_sublots[i][j])
					{
						int incrementincost= singletransfer*(1+((rand()%10)*0.10));
						if(incrementincost	>1)
							tempcost=tempcost+incrementincost;
						else
							tempcost=tempcost+1;


						lvdata[i][j].handlingcost.push_back(tempcost);
					}
					else
						tempcost=0;


				datafile<<" "<<tempcost;
				if(ULmaxsublots>sublot)
					datafile<<",";
			}

			datafile<<"]";
			if(NumVendors-j-1)
				datafile<<",";
		}

		//at the end of each lot close a square bracket
		datafile<<"]";
		if(TLots-i-1)
			datafile<<",";


	}
	datafile<< "\n ]; \n";
	//Finished min-makespan matrix

	datafile.close();
	submap.clear();
	for(unsigned short int i=0; i<TLots;i++)
	{
		for(unsigned short int j=0; j<NumVendors;j++)
		{
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				lvsubdata temp;
				temp.parentlot=i;
				temp.parentvendor=j;
				temp.sublotsused=k;

				submap.push_back(temp);

				
				unsigned short crl = submap.back().parentlot;
			unsigned short crv = submap.back().parentvendor;
			unsigned short crs= submap.back().sublotsused;
			//cout<<crl<<" " <<crv<<" "<<crs<<" "<< " max sub"<<lvdata[crl][crv].maxsublots<<endl<<" lv max="<<lvdata[i][j].maxsublots<<" for "<<i<<" "<<j<<endl;
				
			}

			
		}

		//create mapping array
		submap_inv.clear();
		submap_inv.resize(TLots);
		for(unsigned short i=0; i<TLots; ++i)
		{
			submap_inv[i].resize(NumVendors);
			for(unsigned short j=0; j<NumVendors; ++j)
			{
				submap_inv[i][j].resize(lvdata[i][j].maxsublots);
			}
		}

		unsigned int numsubs = submap.size();

		for(unsigned short i=0; i<numsubs; ++i)
		{
			unsigned short crl = submap[i].parentlot;
			unsigned short crv = submap[i].parentvendor;
			unsigned short crs = submap[i].sublotsused;

			submap_inv[crl][crv][crs]=i;
			
		}

		

	}

	double P=0;
		for(unsigned short int i=0;i<TLots;i++)
		{
			P+=lvdata[i][0].assemblytime;

		}

		//Min Transfer Cost
		double MinTransfer=0;
		for(unsigned short int i=0;i<TLots;i++)
		{
			for(unsigned short j=0;j<NumVendors;j++)
			{
				MinTransfer+=lvdata[i][j].handlingcost[0];
			}

		}

		ObjConstant= MinTransfer+P;
		//cout<<ObjConstant<<" constant after generate data \n";


	/*for(unsigned short int i=0; i<TLots;i++)
	{
	for(unsigned short int j=0; j<NumVendors;j++)
	{
	for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
	cout<<lvdata[i][j].makespansmp[k]<<" ";

	cout<<"\t";
	}
	cout<<endl;

	}*/

	


}



void datalap::generatewindow(bool prec)//prec or perm
{
	//size window matrix
	lvwindow.resize(TLots);
	for(unsigned int i=0;i<lvwindow.size();i++)
		lvwindow[i].resize(NumVendors);
	/*for(unsigned short int i=0;i<lvdata.size();i++)
	{

	for(unsigned short int k=0;k<lvdata.size();k++)
	{
	cout<<precedence[k][i]<<"\t";
	}
	cout<<endl;
	}*/



	if(prec)
	{


		for(unsigned short int i=0;i<lvdata.size();i++)
		{
			for(unsigned short int j=0;j<lvdata[i].size();j++)
			{
				lvwindow[i][j]=0;
				for(unsigned short int k=0;k<lvdata.size();k++)
				{
					if(k!=i)
					{
						lvwindow[i][j]+=lvdata[k][j].vendortime*precedence[k][i]+lvdata[k][j].assemblytime*precedence[i][k];
						//cout<<"\n considering lot "<<i<<" vendor "<<j<<" relative lot "<<k<<" with makespan "<<lvdata[k][j].vendortime<<" X "<<precedence[k][i]<<" + "<<lvdata[k][j].assemblytime<<" X"<<precedence[i][k]<<endl;
					}
				}




			}
		}//finished generating data

		//to check code
		//for(int i=0;i<precedence.size();i++)
		//{
		//	cout<<"Precedence:";
		//	for(int j=0;j<precedence[i].size();j++)
		//		cout<<precedence[i][j]<<"\t";
		//	cout<<endl;

		//

		//	cout<<endl;
		//}

		//	cout<<"Vendor 1 process time:";
		//	for(int j=0;j<lvdata.size();j++)
		//		cout<<lvdata[j].front().vendortime<<"\t";
		//	cout<<endl;

		//	cout<<"Assembly time:";
		//	for(int j=0;j<lvdata.size();j++)
		//		cout<<lvdata[j].front().assemblytime<<"\t";
		//	cout<<endl;

		//	cout<<"Vendor 1 window time:";
		//	for(int j=0;j<lvdata.size();j++)
		//		cout<<lvwindow[j].front()<<"\t";
		//	cout<<endl;


	}
}

double datalap::givensequence(void)
{


	vector<lvsubdata> criticalsublotcandidates, handlingsum;
	double minfeasiblemakespan=0;
	vector<double> cumulativehandlingcost;
	double mincost=0;
	double minhandlingcost;

	generatewindow(1);//lead in and lead out times of each lv 



	//find min feasible makespan
	for(int i=0;i<lvdata.size();i++)//over lots
	{
		for(int j=0;j<lvdata[i].size();j++)//over vendors
		{
			double lvminmakespan=lvdata[i][j].makespansmp[lvdata[i][j].maxsublots-1]+lvwindow[i][j];
			if(lvminmakespan>minfeasiblemakespan)
				minfeasiblemakespan=lvminmakespan;

		}
	}//found minimum
	//cout<<"min feasible"<<minfeasiblemakespan<<endl;

	//The next few lines evaluate the makespan that would result if a particlar sublot were critical
	//we add critical sublot candidates if the makespan assosciated with them is >= minfeasiblemakespan AND add exactly one sublot with makespan=minfeasiblemakespan
	bool notreachedminimum=true;
	for(int i=0;i<lvdata.size();i++)//over lots
	{//cout<<endl;
		for(int j=0;j<lvdata[i].size();j++)//over vendors
		{//cout<<"  ";
			for(int k=0;k<lvdata[i][j].maxsublots;k++)//over sublots
			{
				lvsubdata newsublot;
				newsublot.assosciated_makespan=lvdata[i][j].makespansmp[k]+lvwindow[i][j];
				//cout<<"\n considering lot "<<i<<" vendor "<<j<<" sublot "<<k<<" with makespan "<<lvdata[i][j].makespansmp[k]<<" + "<<lvwindow[i][j]<<endl;

				if(newsublot.assosciated_makespan>minfeasiblemakespan)
				{
					//cout<<"\n added lot "<<i<<" vendor "<<j<<" sublot "<<k<<endl;
					newsublot.handling_cost=lvdata[i][j].handlingcost[k];
					newsublot.parentlot=i;
					newsublot.parentvendor=j;
					newsublot.parent=&lvdata[i][j];
					newsublot.sublotsused=k+1;//+1 since we are using actual number, not array index
					criticalsublotcandidates.push_back(newsublot);
				}

				if(newsublot.assosciated_makespan==minfeasiblemakespan && notreachedminimum)
				{
					//cout<<"\n added lot "<<i<<" vendor "<<j<<" sublot "<<k<<endl;
					notreachedminimum=false;
					newsublot.handling_cost=lvdata[i][j].handlingcost[k];
					newsublot.parentlot=i;
					newsublot.parentvendor=j;
					newsublot.parent=&lvdata[i][j];
					newsublot.sublotsused=k+1;//+1 since we are using actual number, not array index
					criticalsublotcandidates.push_back(newsublot);
				}

			}
		}
	}//finished creating sublot-based vector



	sort(criticalsublotcandidates.begin(),criticalsublotcandidates.end());//sort the vector based on assosciated makespan




	//Calculate the handling cost at each point

	//initialize handling cost to the cost of having a sigle transfer

	msublotsused.clear();//clear previous configration
	cumulativehandlingcost.resize(criticalsublotcandidates.size());
	minhandlingcost=0;//set handling cost to zero
	for(int i=0;i<TLots;i++)//initialize handling cost to the cost of having a single transfer
	{
		vector<int> tempsub;
		for(int j=0;j<NumVendors;j++)
		{
			minhandlingcost+=lvdata[i][j].handlingcost[0];

			tempsub.push_back(1);//initialize the number of sublots to 1 for all LotXVendors
		}
		msublotsused.push_back(tempsub);//set all sublot numbers to 1
	}


	/*cout<<"numbr of critical sublot candidates"<<criticalsublotcandidates.size()<<endl;
	for(int i=0;i<criticalsublotcandidates.size();i++)
	{
	cout<<endl<<"lot "<<criticalsublotcandidates[i].parentlot<<" "<<" vendor "<<criticalsublotcandidates[i].parentvendor<<" sublots used= "<<criticalsublotcandidates[i].sublotsused<<" assosciated makespan="<<criticalsublotcandidates[i].assosciated_makespan<<" cumulative handling cost"<<cumulativehandlingcost[i];

	if(criticalsublotcandidates[i].sublotsused==criticalsublotcandidates[i].parent->maxsublots)
	cout<<"last subot!"<<endl;
	else
	cout<<endl;
	}*/


	//Now calculate handling cost for every point
	for(int i=criticalsublotcandidates.size()-1;i>=0;i--)
	{
		//cout<<"current i="<<i<<endl;
		if(i==criticalsublotcandidates.size()-1)
			cumulativehandlingcost[i]=minhandlingcost;
		else
		{ 
			//cout<<criticalsublotcandidates[i+1].parent->maxsublots<<" "<<criticalsublotcandidates[i+1].sublotsused<<endl;
			lotxvendordata *prevpointsparent=criticalsublotcandidates[i+1].parent;
			int prevpointsublots=criticalsublotcandidates[i+1].sublotsused;
			assert(prevpointsublots<prevpointsparent->maxsublots);

			cumulativehandlingcost[i]=cumulativehandlingcost[i+1]-prevpointsparent->handlingcost[prevpointsublots-1]+prevpointsparent->handlingcost[prevpointsublots];//offset 
		}

	}//finished assigning cumulative handling cost

	assert(cumulativehandlingcost.size()==criticalsublotcandidates.size());
	mincost=cumulativehandlingcost[0]+criticalsublotcandidates[0].assosciated_makespan;//set one as the minimum;
	//cout<<cumulativehandlingcost[0]<<" + "<<criticalsublotcandidates[0].assosciated_makespan;
	unsigned short int minindex=0;//set 0 as the index of the min
	for(unsigned short int i=0; i<cumulativehandlingcost.size(); i++)
	{
		double costifchosen= cumulativehandlingcost[i]+criticalsublotcandidates[i].assosciated_makespan;
		cout<<costifchosen;
		if(costifchosen<mincost)
		{
			mincost=costifchosen;
			minindex=i;
			criticallv=criticalsublotcandidates[i].parent;
		}
	}

	//Found best sublot, record it
	critical_lot=criticalsublotcandidates[minindex].parentlot;
	critical_vendor=criticalsublotcandidates[minindex].parentvendor;
	critical_sublot=criticalsublotcandidates[minindex].sublotsused;

	emptysolution(1,0);
	//set sublot numbers used correctly
	for(unsigned short int i=minindex+1;i<criticalsublotcandidates.size();i++)
	{
		unsigned short int l,v;
		l=criticalsublotcandidates[i].parentlot;
		v=criticalsublotcandidates[i].parentvendor;
		msublotsused[l][v]++;

	}




	/*cout<<endl;
	cout<<"Given CMax="<<criticalsublotcandidates[minindex].assosciated_makespan<<endl;
	for(int i=0;i<criticalsublotcandidates.size();i++)
	{
	cout<<"Lot="<<criticalsublotcandidates[i].parentlot<<" Vendor="<<criticalsublotcandidates[i].parentvendor<<" Sublots Used="<<criticalsublotcandidates[i].sublotsused<<" Assosciated Makespan"<<criticalsublotcandidates[i].assosciated_makespan<<" Handling Cost="<<cumulativehandlingcost[i];

	if(criticalsublotcandidates[i].sublotsused==criticalsublotcandidates[i].parent->maxsublots)
	cout<<"last subot!"<<endl;
	else
	cout<<endl;
	}*/
	//cout<<"returning "<< mincost<<" for given sequence \n";
	return mincost;

}


void datalap::generatepermutation(void)
{
	for(unsigned short int i=0;i<TLots;i++)
	{

	}
}

double datalap::heuristicVolvo(void)
{   


	return 0;
}
double datalap::heuristicLPbased(void)
{
	lpNB();//generated objective co-eff
	IloEnv env;
	IloModel model(env);



	//Building model by rows//
	//create the 2 dimensional matrix of variables y[i][j]


	IloExpr objective(env);	 

	FloatVarArray2 x(env, TLots);
	for(unsigned short int i=0;i<TLots;i++)
		x[i] =FloatVarArray1(env, TLots,0,1,ILOFLOAT);




	IloConstraintArray constraints(env);



	for(unsigned int short i=0;i<TLots;i++)
	{
		IloExpr sumofrow(env),sumofcolumn(env);
		for(unsigned int short j=0;j<TLots;j++)
		{
			sumofrow+=x[i][j];
			sumofcolumn+=x[j][i];
			objective+=m_perm[i][j]*x[i][j];

		}
		constraints.add(sumofcolumn==1);
		constraints.add(sumofrow==1);

	}

	model.add(constraints);
	model.add(IloMaximize(env,objective));
	IloCplex cplex(model);

	IloTimer cputime(env);
	cputime.start();
	//cplex.setOut(env.getNullStream());
	//cplex.setWarning(env.getNullStream());
	cplex.solve();
	m_cputime=cputime.stop();

	m_perm.clear();
	m_perm.resize(TLots);
	for(unsigned int short i=0;i<TLots;i++)
		m_perm[i].resize(TLots);

	for(unsigned int short i=0;i<TLots;i++)
	{

		for(unsigned int short j=0;j<TLots;j++)
		{
			m_perm[i][j]=cplex.getValue(x[i][j]);


		}	
	}

	vector<int> destination;
	destination.resize(TLots);
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<TLots;j++)
		{
			destination[i]+=j*m_perm[i][j];
		}
	}


	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<TLots;j++)
		{
			precedence[i][j]=(destination[i]<destination[j])? 1:0;
		}
	}



	//finished reading results	


	env.end();
	m_optimalobjective=sillygiven();



	return m_optimalobjective;
}

double datalap::optimalNetwork(void)
{	
	IloEnv env;
	IloModel model(env);
	double objectivevaue;


	//Building model by rows//
	//create the 2 dimensional matrix of variables y[i][j]

	IloNumVar Cmax(env);
	IloExpr objective(env);	 

	BoolVarArray2 x(env, TLots);
	for(unsigned short int i=0;i<TLots;i++)
		x[i] =BoolVarArray1(env, TLots);

	BoolVarArray4 zsublots(env,TLots);
	for(unsigned short int i=0;i<TLots;i++)
		zsublots[i]=BoolVarArray3(env,NumVendors);
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			int k=lvdata[i][j].maxsublots;
			zsublots[i][j]=BoolVarArray2(env,k);
		}
	}
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				zsublots[i][j][k]=BoolVarArray1(env,TLots);
			}
		}
	}

	IloConstraintArray constraints(env);


	//create transitive constraints
	//y[i][j]+y[j][k]+y[k][i]<=2 i!=j!=k!=i (implemented here as i>j>k)
	for(unsigned int short i=0;i<TLots;i++)
	{
		IloExpr sumofrow(env),sumofcolumn(env);
		for(unsigned int short j=0;j<TLots;j++)
		{
			sumofrow+=x[i][j];
			sumofcolumn+=x[j][i];

		}
		constraints.add(sumofcolumn==1);
		constraints.add(sumofrow==1);

	}




	//objecctive funtion expression
	objective+=Cmax;

	//create  constraints



	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			IloExpr sublotssumtoone(env);
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				for(unsigned short int l=0;l<TLots;l++)
				{
					sublotssumtoone+=zsublots[i][j][k][l];		
					objective+=lvdata[i][j].handlingcost[k]*zsublots[i][j][k][l];//add to objective
					constraints.add(zsublots[i][j][k][l]<=x[i][l]);//constraints saying sublot assignment and xij are linked

				}				 

			}//made an expression summing subots
			constraints.add(sublotssumtoone==1);
			sublotssumtoone.end();

		}
	}//end i j loop


	for(unsigned short int l=0;l<TLots;l++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			IloExpr assosciatedmakespan(env);

			for(unsigned short int lotto=0;lotto<TLots;lotto++)
				for(unsigned short int i=0;i<l;i++)
					assosciatedmakespan+=lvdata[lotto][j].vendortime*x[lotto][i];

			for(unsigned short int lotto=0;lotto<TLots;lotto++)
				for(unsigned short int i=l+1;i<TLots;i++)
					assosciatedmakespan+=lvdata[lotto][j].assemblytime*x[lotto][i];

			for(unsigned short int i=0;i<TLots;i++)
			{
				for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
				{
					assosciatedmakespan+=lvdata[i][j].makespansmp[k]*zsublots[i][j][k][l];

				}
			}

			constraints.add(assosciatedmakespan<=Cmax);
			assosciatedmakespan.end();

		}//next vendor
	}//next final destination location






	//add constraints to model and solve
	model.add(constraints);
	model.add(IloMinimize(env,objective));
	IloCplex cplex(model);
	cplex.setOut(env.getNullStream());
	///////////////////////
	//set initial values
	//////////////////////

	/* heuristicLPbased();
	IloNumVarArray startVar(env);
	IloNumArray startVal(env);
	for (int i = 0; i<TLots;++i){
	for (int j = 0; j<TLots;j++) {
	startVar.add(x[i][j]);
	startVal.add(m_perm[i][j]);
	}
	}
	cplex.addMIPStart(startVar, startVal, IloCplex::MIPStartAuto,  "secondMIPStart");
	startVal.end();
	startVar.end();
	*/

	////////////////////////
	//finished setting initial values
	////////////////////////
	cplex.setParam(IloCplex::WorkMem, MaxMemory);
	// cplex.setParam(IloCplex::NodeFileInd, 3);
	cplex.setParam(IloCplex::TreLim,MaxMemory);
	cplex.setParam(IloCplex::TiLim,MaxTime);
	IloTimer cputime(env);
	cputime.start();
	//cplex.setOut(env.getNullStream());
	//cplex.setWarning(env.getNullStream());

	//for benders comparison	
	cplex.setParam(IloCplex::NodeFileInd, 3);	
	cplex.setParam(IloCplex::PreInd, IloFalse); 
	cplex.setParam(IloCplex::Threads, 1); 
	cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
	//finished for benders comparison

	if(cplex.solve())
	{
		if (cplex.getStatus()==IloAlgorithm::Optimal)
		{

			cout<<"\n solved"<<endl;

			m_cputime=cputime.stop();
			m_optimalobjective=cplex.getObjValue();
			cout<<"optimal objective="<<cplex.getObjValue()<<" CPU time="<<m_cputime<<" Makespan="<<cplex.getValue(Cmax)<<endl;


			//transfer results to dvars here
			emptysolution(0,1);//clear previous values
			for(unsigned short int i=0;i<TLots;i++)
			{

				for(unsigned short int j=0;j<NumVendors;j++)
				{
					int ijsublotsused=1;
					for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
					{
						for(unsigned short int l=0;l<TLots;l++)
						{

							ijsublotsused+=k*int(cplex.getValue(zsublots[i][j][k][l]));
						}
						msublotsused[i][j]=ijsublotsused;
					}


				}

			}

			vector<int> destination;
			destination.resize(TLots);
			for(unsigned short int i=0;i<TLots;i++)
			{
				for(unsigned short int j=0;j<TLots;j++)
				{
					destination[i]+=j*cplex.getValue(x[i][j]);
				}
			}


			for(unsigned short int i=0;i<TLots;i++)
			{
				for(unsigned short int j=0;j<TLots;j++)
				{
					precedence[i][j]=(destination[i]<destination[j])? 1:0;
				}
			}

			//finished reading results
			solutionoptimal=true;
		}//Finished todo list of optimal found
		else
		{
			solutionoptimal=false;
		}
	}
	else
	{
		cout<<"Something went wrong."<<endl;
		cputime.stop();

	}


	env.end();
	return m_optimalobjective;

}
double datalap::optimalLiftedNetwork(void)
{	
	IloEnv env;
	IloModel model(env);
	double objectivevaue;


	//Building model by rows//
	//create the 2 dimensional matrix of variables y[i][j]

	IloNumVar Cmax(env);
	IloExpr objective(env);	 

	BoolVarArray2 x(env, TLots);
	for(unsigned short int i=0;i<TLots;i++)
		x[i] =BoolVarArray1(env, TLots);

	BoolVarArray4 zsublots(env,TLots);
	for(unsigned short int i=0;i<TLots;i++)
		zsublots[i]=BoolVarArray3(env,NumVendors);
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			int k=lvdata[i][j].maxsublots;
			zsublots[i][j]=BoolVarArray2(env,k);
		}
	}
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				zsublots[i][j][k]=BoolVarArray1(env,TLots);
			}
		}
	}

	IloConstraintArray constraints(env);


	//create transitive constraints
	//y[i][j]+y[j][k]+y[k][i]<=2 i!=j!=k!=i (implemented here as i>j>k)
	for(unsigned int short i=0;i<TLots;i++)
	{
		IloExpr sumofrow(env),sumofcolumn(env);
		for(unsigned int short j=0;j<TLots;j++)
		{
			sumofrow+=x[i][j];
			sumofcolumn+=x[j][i];

		}
		constraints.add(sumofcolumn==1);
		constraints.add(sumofrow==1);

	}




	//objecctive funtion expression
	objective+=Cmax;

	//create  constraints



	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			IloExpr sublotssumtoone(env);
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				for(unsigned short int l=0;l<TLots;l++)
				{
					sublotssumtoone+=zsublots[i][j][k][l];		
					objective+=lvdata[i][j].handlingcost[k]*zsublots[i][j][k][l];//add to objective
					constraints.add(zsublots[i][j][k][l]<=x[i][l]);//constraints saying sublot assignment and xij are linked

				}				 

			}//made an expression summing subots
			constraints.add(sublotssumtoone==1);
			sublotssumtoone.end();

		}
	}//end i j loop


	for(unsigned short int l=0;l<TLots;l++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			IloExpr assosciatedmakespan(env);

			for(unsigned short int lotto=0;lotto<TLots;lotto++)
				for(unsigned short int i=0;i<l;i++)
					assosciatedmakespan+=lvdata[lotto][j].vendortime*x[lotto][i];

			for(unsigned short int lotto=0;lotto<TLots;lotto++)
				for(unsigned short int i=l+1;i<TLots;i++)
					assosciatedmakespan+=lvdata[lotto][j].assemblytime*x[lotto][i];

			for(unsigned short int i=0;i<TLots;i++)
			{
				for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
				{
					assosciatedmakespan+=lvdata[i][j].makespansmp[k]*zsublots[i][j][k][l];

				}
			}

			constraints.add(assosciatedmakespan<=Cmax);
			assosciatedmakespan.end();

		}//next vendor
	}//next final destination location






	//add constraints to model and solve
	model.add(constraints);
	model.add(IloMinimize(env,objective));
	IloCplex cplex(model);
	///////////////////////
	//set initial values
	//////////////////////

	/* heuristicLPbased();
	IloNumVarArray startVar(env);
	IloNumArray startVal(env);
	for (int i = 0; i<TLots;++i){
	for (int j = 0; j<TLots;j++) {
	startVar.add(x[i][j]);
	startVal.add(m_perm[i][j]);
	}
	}
	cplex.addMIPStart(startVar, startVal, IloCplex::MIPStartAuto,  "secondMIPStart");
	startVal.end();
	startVar.end();
	*/

	////////////////////////
	//finished setting initial values
	////////////////////////
	cplex.setParam(IloCplex::WorkMem, 30000);
	// cplex.setParam(IloCplex::NodeFileInd, 3);
	cplex.setParam(IloCplex::TreLim,24000);
	cplex.setParam(IloCplex::TiLim,900);
	IloTimer cputime(env);
	cputime.start();
	//cplex.setOut(env.getNullStream());
	//cplex.setWarning(env.getNullStream());

	if(cplex.solve())
	{
		if (cplex.getStatus()==IloAlgorithm::Optimal)
		{

			cout<<"\n solved"<<endl;

			m_cputime=cputime.stop();
			m_optimalobjective=cplex.getObjValue();
			cout<<"optimal objective="<<cplex.getObjValue()<<" CPU time="<<m_cputime<<" Makespan="<<cplex.getValue(Cmax)<<endl;


			//transfer results to dvars here
			emptysolution(0,1);//clear previous values
			for(unsigned short int i=0;i<TLots;i++)
			{

				for(unsigned short int j=0;j<NumVendors;j++)
				{
					int ijsublotsused=1;
					for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
					{
						for(unsigned short int l=0;l<TLots;l++)
						{

							ijsublotsused+=k*int(cplex.getValue(zsublots[i][j][k][l]));
						}
						msublotsused[i][j]=ijsublotsused;
					}


				}

			}

			vector<int> destination;
			destination.resize(TLots);
			for(unsigned short int i=0;i<TLots;i++)
			{
				for(unsigned short int j=0;j<TLots;j++)
				{
					destination[i]+=j*cplex.getValue(x[i][j]);
				}
			}


			for(unsigned short int i=0;i<TLots;i++)
			{
				for(unsigned short int j=0;j<TLots;j++)
				{
					precedence[i][j]=(destination[i]<destination[j])? 1:0;
				}
			}

			//finished reading results
			solutionoptimal=true;
		}//Finished todo list of optimal found
		else
		{
			solutionoptimal=false;
		}
	}
	else
	{
		cout<<"Something went wrong."<<endl;
		cputime.stop();

	}


	env.end();
	return m_optimalobjective;

}



double datalap::otimalLObased(void )
{
	IloEnv env;
	IloModel model(env);
	env.setDeleter(IloSafeDeleterMode);


	//Building model by rows//
	//create the 2 dimensional matrix of variables y[i][j]

	IloNumVar Cmax(env);
	IloExpr objective(env);	 

	BoolVarArray2 y(env, TLots);
	for(unsigned short int i=0;i<TLots;i++)
		y[i] =BoolVarArray1(env, TLots);

	BoolVarArray3 sublots(env,TLots);
	for(unsigned short int i=0;i<TLots;i++)
		sublots[i]=BoolVarArray2(env,NumVendors);
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			int k=lvdata[i][j].maxsublots;
			sublots[i][j]=BoolVarArray1(env,k);
		}
	}

	IloConstraintArray constraints(env);


	//create transitive constraints
	//y[i][j]+y[j][k]+y[k][i]<=2 i!=j!=k!=i (implemented here as i>j>k)
	for(unsigned int short i=0;i<TLots;i++)
	{
		for(unsigned int short j=0;j<TLots;j++)
		{
			for(unsigned int short k=0;k<TLots;k++)
			{					 
				if(i!=j && j!=k && k!=i)
					constraints.add(y[i][j]+y[j][k]+y[k][i]<=2);
			}}}
	for(unsigned int short i=0;i<TLots;i++)
	{
		for(unsigned int short j=0;j<i;j++)
		{

			constraints.add(y[i][j]+y[j][i]==1);
		}}



	//objecctive funtion expression
	objective+=Cmax;

	//create  constraints



	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			IloExpr expr(env),expr2(env);
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				expr+=sublots[i][j][k];
				expr2+=lvdata[i][j].makespansmp[k]*sublots[i][j][k];
				objective+=lvdata[i][j].handlingcost[k]*sublots[i][j][k];

			}
			expr2+=-1*Cmax;
			for(unsigned short int lotto=0;lotto<TLots;lotto++)
			{
				if(lotto!=i)
					expr2+=lvdata[lotto][j].vendortime*y[lotto][i]+lvdata[lotto][j].assemblytime*y[i][lotto];
			}
			expr2.normalize();
			constraints.add(expr==1);
			constraints.add(expr2<=0);
			expr.end();
			expr2.end();

		}
	}//end i j looping through L X V

	for(unsigned short int i=0;i<TLots;i++)
	{
		constraints.add(y[i][i]==0);
	}







	//add constraints to model and solve
	model.add(constraints);
	model.add(IloMinimize(env,objective));
	IloCplex cplex(model);
	cplex.setParam(IloCplex::WorkMem, MaxMemory);
	// cplex.setParam(IloCplex::NodeFileInd, 3);
	cplex.setParam(IloCplex::TreLim,MaxMemory);
	cplex.setParam(IloCplex::TiLim,MaxTime);
	cplex.setParam(IloCplex::PreInd, IloFalse); 
	cplex.setParam(IloCplex::Threads, 1); 
	cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);

	IloTimer cputime(env);
	cputime.start();
	//cplex.setOut(env.getNullStream());
	//cplex.setWarning(env.getNullStream());
	if(cplex.solve())
	{
		if (cplex.getStatus()==IloAlgorithm::Optimal)
		{


			m_cputime=cputime.stop();
			m_optimalobjective=cplex.getObjValue();
			cout<<"optimal objective="<<cplex.getObjValue()<<" CPU time="<<m_cputime<<" Makespan="<<cplex.getValue(Cmax)<<endl;


			//transfer results to dvars here
			emptysolution(0,1);//clear previous values
			for(unsigned short int i=0;i<TLots;i++)
			{

				for(unsigned short int j=0;j<NumVendors;j++)
				{
					int ijsublotsused=1;
					for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
					{
						ijsublotsused+=k*int(cplex.getValue(sublots[i][j][k]));
					}
					msublotsused[i][j]=ijsublotsused;


				}

			}

			for(unsigned short int i=0;i<TLots;i++)
			{
				for(unsigned short int j=0;j<TLots;j++)
				{
					precedence[i][j]=cplex.getValue(y[i][j]);
				}
			}

			solutionoptimal=true;
		}//Finished todo list of optimal found
		else
		{
			solutionoptimal=false;
			m_cputime=0;
			m_optimalobjective=cplex.getBestObjValue();
		}
	}
	else
	{
		m_cputime=0;
		m_optimalobjective=0;
		cout<<"Something went wrong."<<endl;
		cputime.stop();

	}

	//finished reading results



	env.end();
	return m_optimalobjective;
}
double datalap::sillygiven(void )
{
	IloEnv env;
	IloModel model(env);


	//Building model by rows//
	//create the 2 dimensional matrix of variables y[i][j]

	IloNumVar Cmax(env);
	IloExpr objective(env);	 

	BoolVarArray2 y(env, TLots);
	for(unsigned short int i=0;i<TLots;i++)
		y[i] =BoolVarArray1(env, TLots);

	BoolVarArray3 sublots(env,TLots);
	for(unsigned short int i=0;i<TLots;i++)
		sublots[i]=BoolVarArray2(env,NumVendors);
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			int k=lvdata[i][j].maxsublots;
			sublots[i][j]=BoolVarArray1(env,k);
		}
	}

	IloConstraintArray constraints(env);


	//create transitive constraints
	//y[i][j]+y[j][k]+y[k][i]<=2 i!=j!=k!=i (implemented here as i>j>k)
	for(unsigned int short i=0;i<TLots;i++)
	{
		for(unsigned int short j=0;j<TLots;j++)
		{
			for(unsigned int short k=0;k<TLots;k++)
			{					 
				if(i!=j && j!=k && k!=i)
					constraints.add(y[i][j]+y[j][k]+y[k][i]<=2);
			}}}
	for(unsigned int short i=0;i<TLots;i++)
	{
		for(unsigned int short j=0;j<i;j++)
		{

			constraints.add(y[i][j]+y[j][i]==1);
		}}



	//objecctive funtion expression
	objective+=Cmax;

	//create  constraints



	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			IloExpr expr(env),expr2(env);
			// cout<<"lot"<<i<<"vendor"<<j<<" vp "<<lvdata[i][j].vendortime<<" ap="<<lvdata[i][j].assemblytime<<endl;
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				expr+=sublots[i][j][k];
				//cout<<"lot"<<i<<"vendor"<<j<<" make "<<lvdata[i][j].makespansmp[k]<<" handling="<<lvdata[i][j].handlingcost[k]<<endl;
				expr2+=lvdata[i][j].makespansmp[k]*sublots[i][j][k];
				objective+=lvdata[i][j].handlingcost[k]*sublots[i][j][k];

			}
			expr2+=-1*Cmax;
			for(unsigned short int lotto=0;lotto<TLots;lotto++)
			{
				if(lotto!=i)
					expr2+=lvdata[lotto][j].vendortime*y[lotto][i]+lvdata[lotto][j].assemblytime*y[i][lotto];
			}
			// expr2.normalize();
			//cout<<expr2<<endl;
			constraints.add(expr==1);
			constraints.add(expr2<=0);
			expr.end();
			expr2.end();

		}
	}//end i j looping through L X V

	for(unsigned short int i=0;i<TLots;i++)
	{
		constraints.add(y[i][i]==0);
	}
	for(unsigned int short i=0;i<TLots;i++)
	{
		for(unsigned int short j=0;j<TLots;j++)
		{
			cout<<precedence[i][j]<<"\t";
			if(i!=j)
				constraints.add(y[i][j]==precedence[i][j]);

		}//cout<<endl;
	}	






	//add constraints to model and solve
	model.add(constraints);
	model.add(IloMinimize(env,objective));
	IloCplex cplex(model);
	cplex.setParam(IloCplex::WorkMem, MaxMemory);
	// cplex.setParam(IloCplex::NodeFileInd, 3);
	cplex.setParam(IloCplex::TreLim,MaxMemory);
	cplex.setParam(IloCplex::TiLim,MaxTime);
	IloTimer cputime(env);
	cputime.start();
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	if(cplex.solve())
	{
		if (cplex.getStatus()==IloAlgorithm::Optimal)
		{



			m_cputime=cputime.stop();
			m_optimalobjective=cplex.getObjValue();
			cout<<"optimal objective="<<cplex.getObjValue()<<" CPU time="<<m_cputime<<" Makespan="<<cplex.getValue(Cmax)<<endl;


			//transfer results to dvars here
			emptysolution(0,1);//clear previous values
			for(unsigned short int i=0;i<TLots;i++)
			{

				for(unsigned short int j=0;j<NumVendors;j++)
				{
					int ijsublotsused=1;
					for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
					{
						ijsublotsused+=k*int(cplex.getValue(sublots[i][j][k]));
					}
					msublotsused[i][j]=ijsublotsused;

					//cout<<"lot "<<i<<"vendor "<< j<<" sulots used="<<ijsublotsused<<endl<<lvdata[i][j].maxsublots<<" "<<lvdata[i][j].assemblytime<<" "<<lvdata[i][j].vendortime<<endl<<cplex.getObjective()<<endl;

					/* for(int k=0;k<lvdata[i][j].maxsublots;k++)
					cout<<lvdata[i][j].handlingcost[k]<<endl;*/



				}

			}

			for(unsigned short int i=0;i<TLots;i++)
			{
				for(unsigned short int j=0;j<TLots;j++)
				{
					precedence[i][j]=cplex.getValue(y[i][j]);
				}
			}

			solutionoptimal=true;
		}//Finished todo list of optimal found
		else
		{
			solutionoptimal=false;
		}
	}
	else
	{
		cout<<"Something went wrong."<<endl;
		cputime.stop();

	}

	//finished reading results



	env.end();
	return m_optimalobjective;
}
double datalap::lpformulation(vector<vector<double>> y)
{
	IloEnv env;
	IloModel model(env);

	IloConstraintArray constraints(env);

	//Parameters by sublots
	vector<lvsubdata> v_sublots;
	int total_sublots=0;
	for(unsigned short int i=0;i<TLots;++i)
	{
		for(unsigned short int j=0;j<NumVendors;++j)
		{
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;++k)
			{
				total_sublots++;
				lvsubdata temp_sublot;
				temp_sublot.parentlot=i;
				temp_sublot.parentvendor=j;
				temp_sublot.sublotsused=k;
				v_sublots.push_back(temp_sublot);
			}}}

	FloatVarArray2 z(env,total_sublots);
	FloatVarArray2 z_lamda(env,total_sublots);
	FloatVarArray1 lamda(env,total_sublots);
	FloatVarArray2 y_lamda(env,total_sublots);
	FloatVarArray3 zy(env,total_sublots);












	/*
	Old way of making model
	//Building model by rows//
	//create the 2 dimensional matrix of variables y[i][j]

	IloNumVar Total_cost(env);
	IloExpr objective(env);	 

	BoolVarArray2 y(env, TLots);

	FloatVarArray3 critical_sublot(env,TLots);
	FloatVarArray3 sublot(env,TLots);//=1 when used

	FloatVarArray4 yc_ikva(env,TLots);//y_ik*c_kva
	FloatVarArray4 yc_kiva(env,TLots);//y_ki*c_kva

	FloatVarArray6 zc_iva__jub(env,TLots);//iva preceeds jub and iva is critical
	FloatVarArray6 z_iva__jub(env,TLots);//iva preceeds jub 

	FloatVarArray7 yz1u(env,TLots);
	FloatVarArray7 yz1d(env,TLots);
	FloatVarArray7 yz2u(env,TLots);//y_ir*z_ujb^via
	FloatVarArray7 yz2d(env,TLots);

	//To map constraints to y[i][j] variables we use the following pointers





	//size y
	for(unsigned short int i=0;i<TLots;i++)
	y[i] =BoolVarArray1(env, TLots);

	//size critical_sublot


	for(unsigned short int i=0;i<TLots;i++)
	{
	critical_sublot[i]=FloatVarArray2(env,NumVendors);
	sublot[i]=FloatVarArray2(env,NumVendors);

	for(unsigned short int j=0;j<NumVendors;j++)
	{
	int k=lvdata[i][j].maxsublots;
	critical_sublot[i][j]=FloatVarArray1(env,k,0,1);
	sublot[i][j]=FloatVarArray1(env,k,0,1);
	}
	}

	//size yz12ud, 	 
	for(unsigned short int r=0;r<TLots;r++)
	{
	yz1u[r]=FloatVarArray6(env,TLots);
	yz1d[r]=FloatVarArray6(env,TLots);
	yz2u[r]=FloatVarArray6(env,TLots);
	yz2d[r]=FloatVarArray6(env,TLots);




	for(unsigned short int i=0;i<TLots;i++)
	{
	yz1u[r][i]=FloatVarArray5(env,NumVendors);
	yz2u[r][i]=FloatVarArray5(env,NumVendors);
	yz1d[r][i]=FloatVarArray5(env,NumVendors);
	yz2d[r][i]=FloatVarArray5(env,NumVendors);

	if(!r)
	zc_iva__jub[i]=FloatVarArray5(env,NumVendors);

	if(!r)
	z_iva__jub[i]=FloatVarArray5(env,NumVendors);


	for(unsigned short int v=0;v<NumVendors;v++)
	{
	yz1u[r][i][v]=FloatVarArray4(env,lvdata[i][v].maxsublots);
	yz1d[r][i][v]=FloatVarArray4(env,lvdata[i][v].maxsublots);
	yz2u[r][i][v]=FloatVarArray4(env,lvdata[i][v].maxsublots);
	yz2d[r][i][v]=FloatVarArray4(env,lvdata[i][v].maxsublots);

	if(!r)
	zc_iva__jub[i][v]=FloatVarArray4(env,lvdata[i][v].maxsublots);

	if(!r)
	z_iva__jub[i][v]=FloatVarArray4(env,lvdata[i][v].maxsublots);


	for(unsigned short int a=0;a<lvdata[i][v].maxsublots;a++)
	{
	yz1u[r][i][v][a]=FloatVarArray3(env,TLots);
	yz1d[r][i][v][a]=FloatVarArray3(env,TLots);
	yz2u[r][i][v][a]=FloatVarArray3(env,TLots);
	yz2d[r][i][v][a]=FloatVarArray3(env,TLots);

	if(!r)
	zc_iva__jub[i][v][a]=FloatVarArray3(env,TLots);

	if(!r)
	z_iva__jub[i][v][a]=FloatVarArray3(env,TLots);

	for(unsigned short int j=0;j<TLots;j++)
	{
	yz1u[r][i][v][a][j]=FloatVarArray2(env,NumVendors);
	yz1d[r][i][v][a][j]=FloatVarArray2(env,NumVendors);
	yz2u[r][i][v][a][j]=FloatVarArray2(env,NumVendors);
	yz2d[r][i][v][a][j]=FloatVarArray2(env,NumVendors);

	if(!r)
	zc_iva__jub[i][v][a][j]=FloatVarArray2(env,NumVendors);

	if(!r)
	z_iva__jub[i][v][a][j]=FloatVarArray2(env,NumVendors);

	for(unsigned short int u=0;u<NumVendors;u++)
	{
	yz1u[r][i][v][a][j][u]=FloatVarArray1(env,lvdata[j][u].maxsublots,0,1);
	yz1d[r][i][v][a][j][u]=FloatVarArray1(env,lvdata[j][u].maxsublots,0,1);
	yz2u[r][i][v][a][j][u]=FloatVarArray1(env,lvdata[j][u].maxsublots,0,1);
	yz2d[r][i][v][a][j][u]=FloatVarArray1(env,lvdata[j][u].maxsublots,0,1);

	if(!r)
	zc_iva__jub[i][v][a][j][u]=FloatVarArray1(env,lvdata[j][u].maxsublots,0,1);

	if(!r)
	z_iva__jub[i][v][a][j][u]=FloatVarArray1(env,lvdata[j][u].maxsublots,0,1);

	}}
	}}}
	}

	//size yc
	for(unsigned short int r=0;r<TLots;r++)
	{
	yc_ikva[r]= FloatVarArray3(env, TLots);
	yc_kiva[r]= FloatVarArray3(env, TLots);

	for(unsigned short int i=0;i<TLots;i++)
	{
	yc_ikva[r][i] =FloatVarArray2(env, NumVendors);
	yc_kiva[r][i] =FloatVarArray2(env, NumVendors);

	for(unsigned short int v=0;v<NumVendors;v++)
	{
	yc_ikva[r][i][v]= FloatVarArray1(env,lvdata[i][v].maxsublots,0,1);
	yc_kiva[r][i][v]= FloatVarArray1(env,lvdata[i][v].maxsublots,0,1);

	}}};






	//create transitive constraints
	//y[i][j]+y[j][k]+y[k][i]<=2 i!=j!=k!=i (implemented here as i>j>k)
	for(unsigned int short i=0;i<TLots;i++)
	{
	for(unsigned int short j=0;j<i;j++)
	{
	for(unsigned int short k=0;k<j;k++)
	{					 

	constraints.add(y[i][j]+y[j][k]+y[k][i]<=2);
	constraints.add(y[i][k]+y[k][j]+y[j][i]<=2);
	}}}
	for(unsigned int short i=0;i<TLots;i++)
	{
	for(unsigned int short j=0;j<i;j++)
	{								 
	constraints.add(y[i][j]+y[j][i]==1);
	}}

	for(unsigned short int i=0;i<TLots;i++)
	{
	constraints.add(y[i][i]==0);
	}



	//objecctive funtion expression
	objective+=Total_cost;

	//create  constraints
	//add constraint only only one sublot is critical
	IloExpr expr(env);
	for(unsigned short int i=0;i<TLots;i++)
	{
	for(unsigned short int v=0;v<NumVendors;v++)
	{
	//expr_onesublotused(env);
	for(unsigned short int a=0;a<lvdata[i][v].maxsublots;a++)
	{
	expr+=critical_sublot[i][v][a];
	}
	}
	}
	constraints.add(expr==1);
	expr.clear();
	assert(expr.getEnv()==env);



	//add constraint only only one sublot is critical
	//IloExpr expr(env);
	for(unsigned short int i=0;i<TLots;i++)
	{
	for(unsigned short int v=0;v<NumVendors;v++)
	{

	for(unsigned short int a=0;a<lvdata[i][v].maxsublots;a++)
	{
	//expr+=critical_sublot[i][v][a];
	}
	}
	}
	constraints.add(expr==1);
	expr.end();









	//add constraints to model and solve
	model.add(constraints);
	model.add(IloMinimize(env,objective));
	IloCplex cplex(model);
	cplex.setParam(IloCplex::WorkMem, 8000);
	IloTimer cputime(env);
	cputime.start();
	//cplex.setOut(env.getNullStream());
	//cplex.setWarning(env.getNullStream());
	cplex.solve();
	//cout<<"optimal objective="<<cplex.getObjValue()<<" CPU time="<<cputime.stop()<<" Makespan="<<cplex.getValue(Cmax)<<endl;

	//transfer results to dvars here
	emptysolution(0,1);//clear previous values


	for(unsigned short int i=0;i<TLots;i++)
	{
	for(unsigned short int j=0;j<TLots;j++)
	{
	precedence[i][j]=cplex.getValue(y[i][j]);
	cout<<precedence[i][j]<<" ";
	}
	cout<<endl;
	}

	//finished reading results



	env.end();*/
	return 0;
}

IloNum datalap::lpLO(void)
{
	IloEnv env;
	IloModel model(env);


	//Building model by rows//
	//create the 2 dimensional matrix of variables y[i][j]

	IloNumVar Cmax(env);
	IloExpr objective(env);	 

	FloatVarArray2 y(env, TLots);
	for(unsigned short int i=0;i<TLots;i++)
		y[i] =FloatVarArray1(env, TLots,0,1,ILOFLOAT);

	FloatVarArray3 sublots(env,TLots);
	for(unsigned short int i=0;i<TLots;i++)
		sublots[i]=FloatVarArray2(env,NumVendors);
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			int k=lvdata[i][j].maxsublots;
			sublots[i][j]=FloatVarArray1(env,k,0,1,ILOFLOAT);
		}
	}

	IloConstraintArray constraints(env);


	//create transitive constraints
	//y[i][j]+y[j][k]+y[k][i]<=2 i!=j!=k!=i (implemented here as i>j>k)
	for(unsigned int short i=0;i<TLots;i++)
	{
		for(unsigned int short j=0;j<TLots;j++)
		{
			for(unsigned int short k=0;k<TLots;k++)
			{					 
				if(i!=j && j!=k && k!=i)
					constraints.add(y[i][j]+y[j][k]+y[k][i]<=2);
			}}}
	for(unsigned int short i=0;i<TLots;i++)
	{
		for(unsigned int short j=0;j<i;j++)
		{

			constraints.add(y[i][j]+y[j][i]==1);
		}}



	//objecctive funtion expression
	objective+=Cmax;

	//create  constraints



	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			IloExpr expr(env),expr2(env);
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				expr+=sublots[i][j][k];
				expr2+=lvdata[i][j].makespansmp[k]*sublots[i][j][k];
				objective+=lvdata[i][j].handlingcost[k]*sublots[i][j][k];

			}
			expr2+=-1*Cmax;
			for(unsigned short int lotto=0;lotto<TLots;lotto++)
			{
				if(lotto!=i)
					expr2+=lvdata[lotto][j].vendortime*y[lotto][i]+lvdata[lotto][j].assemblytime*y[i][lotto];
			}
			expr2.normalize();
			constraints.add(expr==1);
			constraints.add(expr2<=0);
			expr.end();
			expr2.end();

		}
	}//end i j looping through L X V

	for(unsigned short int i=0;i<TLots;i++)
	{
		constraints.add(y[i][i]==0);
	}







	//add constraints to model and solve
	model.add(constraints);
	model.add(IloMinimize(env,objective));
	IloCplex cplex(model);
	//cplex.setParam(IloCplex::WorkMem, 1800);
	//cplex.setParam(IloCplex::NodeFileInd, 3);
	//cplex.setParam(IloCplex::TreLim,1500);
	IloTimer cputime(env);
	cputime.start();
	//cplex.setOut(env.getNullStream());
	//cplex.setWarning(env.getNullStream());
	cplex.solve();


	m_cputime=cputime.stop();
	m_optimalobjective=cplex.getObjValue();
	cout<<"optimal objective="<<cplex.getObjValue()<<" CPU time="<<m_cputime<<" Makespan="<<cplex.getValue(Cmax)<<endl;


	//transfer results to dvars here
	emptysolution(0,1);//clear previous values
	for(unsigned short int i=0;i<TLots;i++)
	{

		for(unsigned short int j=0;j<NumVendors;j++)
		{
			int ijsublotsused=1;
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				ijsublotsused+=k*int(cplex.getValue(sublots[i][j][k]));
			}
			msublotsused[i][j]=ijsublotsused;


		}

	}

	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<TLots;j++)
		{
			precedence[i][j]=cplex.getValue(y[i][j]);
		}
	}

	//finished reading results



	env.end();
	return m_optimalobjective;
}

IloNum datalap::lpNB(void)
{
	IloEnv env;
	IloModel model(env);
	//double objectivevaue;


	//Building model by rows//
	//create the 2 dimensional matrix of variables y[i][j]

	IloNumVar Cmax(env);
	IloExpr objective(env);	 

	FloatVarArray2 x(env, TLots);
	for(unsigned short int i=0;i<TLots;i++)
		x[i] =FloatVarArray1(env, TLots,0,1,ILOFLOAT);

	FloatVarArray4 zsublots(env,TLots);
	for(unsigned short int i=0;i<TLots;i++)
		zsublots[i]=FloatVarArray3(env,NumVendors);
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			int k=lvdata[i][j].maxsublots;
			zsublots[i][j]=FloatVarArray2(env,k);
		}
	}
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				zsublots[i][j][k]=FloatVarArray1(env,TLots,0,1,ILOFLOAT);
			}
		}
	}

	IloConstraintArray constraints(env);


	//create transitive constraints
	//y[i][j]+y[j][k]+y[k][i]<=2 i!=j!=k!=i (implemented here as i>j>k)
	for(unsigned int short i=0;i<TLots;i++)
	{
		IloExpr sumofrow(env),sumofcolumn(env);
		for(unsigned int short j=0;j<TLots;j++)
		{
			sumofrow+=x[i][j];
			sumofcolumn+=x[j][i];

		}
		constraints.add(sumofcolumn==1);
		constraints.add(sumofrow==1);

	}




	//objecctive funtion expression
	objective+=Cmax;

	//create  constraints



	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			IloExpr sublotssumtoone(env);
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				for(unsigned short int l=0;l<TLots;l++)
				{
					sublotssumtoone+=zsublots[i][j][k][l];		
					objective+=lvdata[i][j].handlingcost[k]*zsublots[i][j][k][l];//add to objective
					constraints.add(zsublots[i][j][k][l]<=x[i][l]);//constraints saying sublot assignment and xij are linked

				}				 

			}//made an expression summing subots
			constraints.add(sublotssumtoone==1);
			sublotssumtoone.end();

		}
	}//end i j loop


	for(unsigned short int l=0;l<TLots;l++)
	{
		for(unsigned short int j=0;j<NumVendors;j++)
		{
			IloExpr assosciatedmakespan(env);

			for(unsigned short int lotto=0;lotto<TLots;lotto++)
				for(unsigned short int i=0;i<l;i++)
					assosciatedmakespan+=lvdata[lotto][j].vendortime*x[lotto][i];

			for(unsigned short int lotto=0;lotto<TLots;lotto++)
				for(unsigned short int i=l+1;i<TLots;i++)
					assosciatedmakespan+=lvdata[lotto][j].assemblytime*x[lotto][i];

			for(unsigned short int i=0;i<TLots;i++)
			{
				for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
				{
					assosciatedmakespan+=lvdata[i][j].makespansmp[k]*zsublots[i][j][k][l];

				}
			}

			constraints.add(assosciatedmakespan<=Cmax);
			assosciatedmakespan.end();

		}//next vendor
	}//next final destination location






	//add constraints to model and solve
	model.add(constraints);
	model.add(IloMinimize(env,objective));
	IloCplex cplex(model);
	//cplex.setParam(IloCplex::WorkMem, 1800);
	//cplex.setParam(IloCplex::NodeFileInd, 3);
	//cplex.setParam(IloCplex::TreLim,1500);
	IloTimer cputime(env);
	cputime.start();
	//cplex.setOut(env.getNullStream());
	//cplex.setWarning(env.getNullStream());
	cplex.solve();
	m_cputime=cputime.stop();
	m_optimalobjective=cplex.getObjValue();
	cout<<"optimal objective="<<m_optimalobjective<<" CPU time="<<m_cputime<<" Makespan="<<cplex.getValue(Cmax)<<endl;


	//transfer results to dvars here
	emptysolution(0,1);//clear previous values
	for(unsigned short int i=0;i<TLots;i++)
	{

		for(unsigned short int j=0;j<NumVendors;j++)
		{
			int ijsublotsused=1;
			for(unsigned short int k=0;k<lvdata[i][j].maxsublots;k++)
			{
				for(unsigned short int l=0;l<TLots;l++)
				{

					ijsublotsused+=k*int(cplex.getValue(zsublots[i][j][k][l]));
				}
				msublotsused[i][j]=ijsublotsused;
			}


		}

	}

	vector<int> destination;
	destination.resize(TLots);
	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<TLots;j++)
		{
			destination[i]+=j*cplex.getValue(x[i][j]);
		}
	}


	for(unsigned short int i=0;i<TLots;i++)
	{
		for(unsigned short int j=0;j<TLots;j++)
		{
			precedence[i][j]=(destination[i]<destination[j])? 1:0;
		}
	}

	//read permutation
	m_perm.clear();
	m_perm.resize(TLots);
	for(unsigned int short i=0;i<TLots;i++)
		m_perm[i].resize(TLots);

	for(unsigned int short i=0;i<TLots;i++)
	{

		for(unsigned int short j=0;j<TLots;j++)
		{
			m_perm[i][j]=cplex.getValue(x[i][j]);
		}


	}

	//finished reading results	
	env.end();


	return m_optimalobjective;
}

void datalap::emptysolution(const int defaultsublots=0, bool reset_permutation=0)
{
	msublotsused.clear();


	for(int i=0;i<TLots;i++)
	{
		vector<int> tempsub;
		for(int j=0;j<NumVendors;j++)
		{

			tempsub.push_back(defaultsublots);//initialize the number of sublots to "defaultsublots" for all LotXVendors
		}
		msublotsused.push_back(tempsub);//set all sublot numbers to 0
	}

	if(reset_permutation)
	{
		precedence.clear();
		for(int i=0;i<TLots;i++)
		{
			vector<double> tempsub;
			for(int j=0;j<TLots;j++)
			{
				if(i<j)
					tempsub.push_back(1);
				else
					tempsub.push_back(0);//initialize the number 
			}
			precedence.push_back(tempsub);//set all  numbers to 0
		}
	}//reset permutation to 1,2,3,4..



}


void datalap::experiment(int maxlots, int maxvendors)
{
	cout<<"started experiment\n";
	std::ofstream resultfile;
	
	generatedata(1, maxlots, maxvendors, 4);
	resultfile.open("resuts.txt",std::fstream::app);
	resultfile<<"\n Max time="<<MaxTime<<", Max memory="<<MaxMemory<<"\n";
	resultfile<<"\n root, Tlots, Numvendors, maxmaxsublots, cputime_NBB1, Optimal_NBB1,cputime_NB, Optimal_NB, cputime_LO, Optimal_LO\n";
	
	resultfile.close();
	
	unsigned short root = 10;
	/*for(unsigned short int lots=10;lots<=maxlots;lots+=10)
	{

		for(unsigned short int ven=10;ven<=maxvendors;ven+=10)
		{
			
			for(unsigned short int msubs=3;msubs<=12;msubs+=3)
			{*/

for(unsigned short int lots=10;lots<=41;lots+=10)
	{

		for(unsigned short int ven=10;ven<=41;ven+=10)
		{
			
			for(unsigned short int msubs=3;msubs<=13;msubs+=3)
			{
				
				++root;

				resultfile.open("resuts.txt",std::fstream::app);
				resultfile<<root<<",";
				resultfile<<lots<<",";
				resultfile<<ven<<",";
				resultfile<<msubs<<",";
				resultfile.close();

				
				generatedata(root, lots, ven, msubs);
				
				cout<<"generated data\n";
				
				double opt_temp = 0;
				opt_temp = optimalNetworkBendersallycuts();//try optimalNetworkBenders3() for pre-processing only cuts

				resultfile.open("resuts.txt",std::fstream::app);
				resultfile<<"Benders time,"<<m_cputime<<","<<opt_temp<<",";
				cout<<"FINISHED BENDERS"<<endl;
				resultfile.close();

				resultfile.open("resuts.txt",std::fstream::app);
				opt_temp= optimalNetwork();//change parameters to multi-threaded when not comparing benders
				resultfile<<m_cputime<<","<<opt_temp<<",";
				cout<<"FINISHED NB \n\n";
				resultfile.close();

				/*resultfile.open("resuts.txt",std::fstream::app);
				opt_temp = otimalLObased();
				resultfile<<m_cputime<<","<<opt_temp<<"\n";
				cout<<"FINISHED LO \n\n";
				resultfile.close();*/
			}}}
	
		
	

}
void datalap::printprec(void)
{
	for(unsigned short i=0;i<TLots;++i)
	{
		
		for(unsigned short j=0;j<TLots;++j)
		{
			
			cout<<precedence[i][j]<<" ";

		}	
		cout<<endl;

	}
	cout<<endl;
}
void datalap::changelv(void)
{
	for(unsigned short i=0; i<TLots; ++i)
	{
		for(unsigned short v=0; v<NumVendors; ++v)
		{
			if(lvdata[i][v].assemblytime == lvdata[i][v].vendortime)
			{
				lvdata[i][v].vendortime++;
				lvdata[i][v].makespansmp[0]++;
			}
		}

	}
	return;
}
void datalap::clean_perm(vector<vector<double>> perm)
{
	IloEnv env;
	IloModel model(env);

	
	
	IloConstraintArray constraints(env);
	IloExpr objective(env);	 
	BoolVarArray2 x(env, TLots);
	for(unsigned short int i=0;i<TLots;i++)
		x[i] =BoolVarArray1(env, TLots);

	


	for(unsigned int short i=0;i<TLots;i++)
	{
		IloExpr sumofrow(env),sumofcolumn(env);
		for(unsigned int short j=0;j<TLots;j++)
		{
			sumofrow+=x[i][j];
			sumofcolumn+=x[j][i];
			objective+=perm[i][j]*x[i][j];

		}
		constraints.add(sumofcolumn==1);
		constraints.add(sumofrow==1);

	}

	model.add(IloMinimize(env,objective));
	model.add(constraints);
	IloCplex cplex(model);
	cplex.setOut(env.getNullStream());

	if(cplex.solve())
	{

		m_perm.resize(TLots);
	for(unsigned short int i=0;i<TLots;i++)
		m_perm[i].resize(TLots);

	
		for(unsigned int short i=0;i<TLots;i++)
		{

			for(unsigned int short j=0;j<TLots;j++)
			{
				m_perm[i][j] = cplex.getValue(x[i][j]);

			}
		}
	}




}
datalap::datalap(void)
{
	me=this;
}
datalap::~datalap(void)
{
}
