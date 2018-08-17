#pragma once
#include <vector>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ilcplex/ilocplex.h>

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

//#include <iostream.h>
class datalap
{

private:
	unsigned short int ULlots, ULvendors, ULmaxsublots, ULprocesstime, ULtransfercost;//upper limit on parameter sizes

	
	
public:
	struct lotxvendordata{
		unsigned short int vendortime,assemblytime;
		unsigned short int maxsublots;
		std::vector<int>  makespansmp, handlingcost;//minmakespan for single machine problem with that lotXvendorXassemblymachine
	};

	struct lvsubdata{ //lot X vendor X sublot data
		int parentlot,parentvendor;//0 to L-1 and 0 to M-1
		lotxvendordata *parent;
		int sublotsused;
		double assosciated_makespan; //eventaully calculated overall makespan, if this is critical
		int handling_cost; //cost of using parent lot X vendor with this many sublots

		bool operator<(const lvsubdata &a){
			bool smallerthaninput=false;
		if(assosciated_makespan<a.assosciated_makespan)
			smallerthaninput=true;
		if(parentlot==a.parentlot && parentvendor==a.parentvendor)
		smallerthaninput=sublotsused>a.sublotsused;
		return smallerthaninput;};
		};
	


	//ilocplex related
	int MaxMemory, MaxTime;

	double ObjConstant;
	bool solutionoptimal;
	bool objective_check;
	IloEnv env;
	IloModel LinearOrdering;
	double m_optimalobjective;
	double m_cputime;
	//dvar related storage
	datalap* me;
	
	std::vector<std::vector<double>> m_perm;
	std::vector<std::vector<double>> precedence;//to store sequence
	std::vector<std::vector<std::vector<double>>> m_perm_arr;
	
	std::vector<std::vector<double>> lvwindow;
	std::vector<lvsubdata> submap;
	std::vector<std::vector<std::vector<unsigned short>>> submap_inv;
	std::vector<std::vector<double>> cut_yij_coeff;
	
	std::vector<std::vector<int>> msublotsused;//0 to N-1 //to store handling cost implicitly
	lotxvendordata *criticallv;//don't use; sometimes heuristics stores the makespan implicitly, i.e., given the number of sublots used we can find the solution

	//input data
	 std::vector<std::vector<lotxvendordata>> lvdata;//lot,vendor This is the main data repository, it also represents the index of each lot.
	unsigned short int TLots;
	unsigned short int NumVendors;

	//data for sub problem
	static datalap* master_instance;
	unsigned short critical_lot,critical_sublot, critical_vendor;
	
	//functions
	void generate_birkhoff(void);
	void perm_to_prec(void);
	void prec_to_perm(void);
	void generatedata(int root, int Lots, int Vendors, int ubsublots);//To generate the data
	void generatepermutation(void);
	void emptysolution(const int defaultsublots, const bool reset_permutation);//set all precedence values to zero and msublotsused to 0
	void generatewindow(bool precedence);//precedence==1 then use precedence matrix, permutation o/w
	double givensequence(void);//optimal solution for a given sequence
	void printprec(void);

	double separate( IloExpr cutLhs, 
	BoolVarArray2 x, FloatVarArray2 y, IloNum GivenObj,
	unsigned short crLot, unsigned short crSublot, unsigned short crVendor);

	double separatesmall( IloExpr cutLhs, 
	BoolVarArray2 x, FloatVarArray2 y, IloNum GivenObj,
	std::vector<unsigned short> crsublotconfigs);

	double separatesmallbooly( IloExpr cutLhs, 
	BoolVarArray2 y, IloNum GivenObj,
	std::vector<unsigned short> crsublotconfigs);

	double separatepermx(IloExpr cutLhs,  BoolVarArray2 x);

	double separateperm(IloExpr cutLhs, FloatVarArray2 y);

	double separate_allcritical( IloExpr cutLhs, 
	BoolVarArray2 x, BoolVarArray2 y, IloNum GivenObj,
	unsigned short crLot, unsigned short crSublot, unsigned short crVendor);

	IloBool separate2( IloExpr cutLhs, 
	BoolVarArray1 lamdap, BoolVarArray2 y, IloNum GivenObj,
	unsigned short crLot, unsigned short crSublot, unsigned short crVendor, double lambdap_value);

	//heuristic and optimal solutions to the problems, all return the objective function value.
	double heuristicVolvo(void);
	double heuristicLPbased(void);
	double optimalNetwork(void);
	double optimalLiftedNetwork(void);
	double optimalNetworkBenders1(void);
	double optimalNetworkBenders2(void);
	double optimalNetworkBenders3(void);
	double optimalNetworkBendersallcuts(void);
	double optimalNetworkBendersallycuts(void);

	double otimalLObased(void);
	double lpformulation(std::vector<std::vector<double>>);
	double sillygiven(void);

	static int trial;
	static double test(void);
	
	IloNum lpLO(void);
	IloNum lpNB(void);
	void experiment(int lots, int vendors);
	void changelv(void);
	void datalap::clean_perm(std::vector<std::vector<double>> perm);

	double benders(void);
	double benders2(void);
		
	datalap(void);
	~datalap(void);

	/* Discarded code
	if ( solStatus == IloAlgorithm::Optimal ) {

				cout<<"Duals:\n";
				IloNumArray duals(masterEnv);
					masterCplex.getDuals(duals,ranger);
					for(unsigned int short i=0;i<ranger.getSize();i++)
				{
					//int *index (int *) constraints[i].getObject();
					
					cout<<duals[i];
					
					cout<<"\n";

				}

			//IloNumVarArray lamda(workerEnv,NumSublots,0,1);
		//FloatVarArray2 z(workerEnv, NumSublots);
		//FloatVarArray2 zlamda(workerEnv, NumSublots);
		//FloatVarArray2 lamday(workerEnv, NumSublots);
		//FloatVarArray3 yz(workerEnv, NumSublots);
		//FloatVarArray3 zy(workerEnv, NumSublots);
		//for(unsigned short int i=0;i<NumSublots;i++)
		//{
		//	z[i]= IloNumVarArray(workerEnv,NumSublots,0,1);
		//	zlamda[i]=IloNumVarArray(workerEnv,NumSublots,0,1);
		//	lamday[i]=IloNumVarArray(workerEnv,TLots,0,1);
		//	yz[i]=FloatVarArray2(workerEnv,NumSublots);
		//	zy[i]=FloatVarArray2(workerEnv,NumSublots);


		//	for(unsigned short int j=0;j<NumSublots;j++)
		//	{
		//		yz[i][j] = IloNumVarArray(workerEnv,TLots,0,1);//y_ij*z_i_alpha_v
		//		zy[i][j] = IloNumVarArray(workerEnv,TLots,0,1);//y_ji*z_i_alpha_v

		//	}
		//}
		//
		//Finished declaring variables

		

		//Set Sum of lamda==1
		/*IloExpr lsum(workerEnv);
		IloExpr WorkObjExpr(workerEnv);

		for(unsigned short int i=0;i<NumSublots;i++)
		{
			lsum +=lamda[i];

		}
		IloRange lamda_pick_one_s(workerEnv,1,lsum,1);*/
		//Create Objective funtion
		/*for(unsigned short sub=0; sub<NumSublots;++sub)
		{
			int Lotassosciated=sublot[sub].parentlot;
			int Vendorassosciated=sublot[sub].parentvendor;
			int Sublotassosciated=sublot[sub].sublotsused;
			int SMMassosciated=lvdata[Lotassosciated][Vendorassosciated].makespansmp[Sublotassosciated];
			

			int lamdaM=P+SMMassosciated-lvdata[Lotassosciated][Vendorassosciated].assemblytime+MinTransfer-lvdata[Lotassosciated][Vendorassosciated].handlingcost[0]+lvdata[Lotassosciated][Vendorassosciated].handlingcost[Sublotassosciated];
			WorkObjExpr+=lamdaM*lamda[sub];
		}

		workerMod.add(lamda_pick_one_s);
		workerMod.add(IloMinimize(workerEnv,WorkObjExpr));
		workerCplex.extract(workerMod);

		


		
		//masterCplex.use(BendersLazyCallback(masterEnv, workerCplex,y,NumSublots,lamda_pick_one_s,lamda));

		int CriticalSublot=1;///Choose critical

		FloatVarArray1 pi_c(workerEnv,NumSublots,0,IloInfinity);
		FloatVarArray1 pi_cb(workerEnv,NumSublots,0,IloInfinity);
		FloatVarArray1 piz(workerEnv,NumSublots,0,IloInfinity);
		IloNumVar pi_tc(workerEnv,0,IloInfinity);

		FloatVarArray2 pi1(workerEnv,NumSublots);
		FloatVarArray2 pi2(workerEnv,NumSublots);
		FloatVarArray2 pi3(workerEnv,NumSublots);
		FloatVarArray2 pi4(workerEnv,NumSublots);
		FloatVarArray2 pi5(workerEnv,NumSublots);
		FloatVarArray2 pi6(workerEnv,NumSublots);

		
		FloatVarArray2 piy_ki_z(workerEnv,NumSublots);
		FloatVarArray2 piy_kj_z(workerEnv,NumSublots);

		for(unsigned short i=0; i<NumSublots;++i)
		{
			pi1[i]= FloatVarArray1(workerEnv,TLots,0,IloInfinity);
			pi2[i]= FloatVarArray1(workerEnv,TLots,0,IloInfinity);
			pi3[i]= FloatVarArray1(workerEnv,TLots,0,IloInfinity);
			pi4[i]= FloatVarArray1(workerEnv,TLots,0,IloInfinity);
			pi5[i]= FloatVarArray1(workerEnv,TLots,0,IloInfinity);
			pi6[i]= FloatVarArray1(workerEnv,TLots,0,IloInfinity);

			piy_ki_z[i]=FloatVarArray1(workerEnv,TLots,0,IloInfinity);
			piy_kj_z[i]=FloatVarArray1(workerEnv,TLots,0,IloInfinity);

		}
		//objective
		//First give yy some value
		IloNumArray2 yy(workerEnv,TLots);
			for(unsigned short i=0;i<TLots;i++)
				yy[i]= IloNumArray(workerEnv,TLots);
		IloExpr workerObjective(workerEnv);
		for(unsigned short i=0;i<TLots;i++)
		{
			for(unsigned short j=0;j<TLots;j++)
			{
				if(i<j)
				yy[i][j]=1;
				else
				yy[i][j]=0;
			}
		}//gave yy some value

		for(unsigned short j=0;j<NumSublots;j++)
		{
			if(j!=CriticalSublot)
				workerObjective+=piz[j];
			for(unsigned short k=0;k<TLots;k++)
			{
				if(j!=CriticalSublot && k!=sublot[j].parentlot)
					workerObjective+=piy_kj_z[j][k];

				if(j!=CriticalSublot && k!=sublot[CriticalSublot].parentlot)
					workerObjective+=piy_ki_z[j][k];
			}
		}

		double pi_tc_coeff=0;
		int i_crt=sublot[CriticalSublot].parentlot;
		int v_crt=sublot[CriticalSublot].parentvendor;
		int a_crt=sublot[CriticalSublot].sublotsused;

		pi_tc_coeff+=MinTransfer+ P+lvdata[i_crt][v_crt].makespansmp[a_crt]-lvdata[i_crt][v_crt].assemblytime;
		cout<<"pi_tc_coeff="<<pi_tc_coeff<<endl<<endl;
		for(unsigned short i=0;i<TLots;i++)
		{
			if(i!=i_crt)
				pi_tc_coeff+=yy[i][i_crt]*(lvdata[i][v_crt].vendortime-lvdata[i][v_crt].assemblytime);
		}
		cout<<"pi_tc_coeff="<<pi_tc_coeff<<endl<<endl;
		
		workerObjective+=pi_tc_coeff*pi_tc;

		//pi_c_bar coeff
		for(unsigned short j=0;j<NumSublots;++j)
		{
			double pi_c_bar_coeff=0;
			for(unsigned short k=0;k<TLots;++k)
			{
				if(k!=i_crt && j!=CriticalSublot)
					pi_c_bar_coeff+=yy[k][i_crt]*(lvdata[k][v_crt].vendortime-lvdata[k][v_crt].assemblytime);
			}

			for(unsigned short k=0;k<TLots;++k)
			{
				if(k!=sublot[j].parentlot && j!=CriticalSublot)
					pi_c_bar_coeff+=-1*yy[k][sublot[j].parentlot]*(lvdata[k][sublot[j].parentvendor].vendortime-lvdata[k][sublot[j].parentvendor].assemblytime);
			}

			pi_c_bar_coeff+=lvdata[sublot[j].parentlot][sublot[j].parentvendor].assemblytime-lvdata[i_crt][sublot[j].parentvendor].assemblytime;
			pi_c_bar_coeff+=lvdata[i_crt][v_crt].makespansmp[a_crt]-lvdata[sublot[j].parentlot][sublot[j].parentvendor].makespansmp[sublot[j].sublotsused];

			workerObjective+=pi_c_bar_coeff*pi_cb[j];
			cout<<"pi_c_bar_coeff["<<j<<"]="<<pi_c_bar_coeff<<endl;


		}
		

     	//other obj
		for(unsigned short i=0;i<NumSublots;i++)
		{
			double pi_cb_coeff=0;
			for(unsigned short k=0;k<TLots;k++)
			{
				if(k!=i_crt && i!=CriticalSublot)
				{
					workerObjective+=yy[k][i_crt]*pi1[i][k];
					workerObjective+=yy[i_crt][k]*pi2[i][k];
				}

				int j_this=sublot[i].parentlot;
				if(k!=j_this && i!=CriticalSublot)
				{
					workerObjective+=yy[k][j_this]*pi4[i][k];
					workerObjective+=yy[j_this][k]*pi5[i][k];
				}
			}
		}
		


		workerMod.add(IloMinimize(workerEnv,workerObjective));

		//constraints
		IloRange TC_constraint(workerEnv,1,pi_tc);
		workerMod.add(TC_constraint);

		IloRangeArray  z_iav(workerEnv,NumSublots);
		for(unsigned short i=0; i<NumSublots;++i)
		{
			if(i!=CriticalSublot)
			{
				IloExpr exz_iav(workerEnv);
				exz_iav+=-1*tau_star[i]*pi_tc;
				exz_iav+=piz[i];
				exz_iav+=(sublot[i].parent->assemblytime-sublot[i].parent->makespansmp[sublot[i].sublotsused]-sublot[CriticalSublot].parent->assemblytime+sublot[CriticalSublot].parent->makespansmp[sublot[CriticalSublot].sublotsused])*pi_c[i];

				exz_iav+=(sublot[i].parent->assemblytime-sublot[i].parent->makespansmp[sublot[i].sublotsused]-sublot[CriticalSublot].parent->assemblytime+sublot[CriticalSublot].parent->makespansmp[sublot[CriticalSublot].sublotsused])*pi_cb[i];




				for(unsigned short j=0; j<TLots;++j)
				{
					exz_iav+=pi2[i][j]-pi3[i][j]+pi5[i][j]-pi6[i][j];
				}

				z_iav[i]= IloRange(workerEnv,0,exz_iav,IloInfinity);
				exz_iav.end();
			}

		}
		
		IloRangeArray wk_constraints(workerEnv);
		for(unsigned short i=0; i<NumSublots;++i)
		{
			if(i!=CriticalSublot)
			{
				

				for(unsigned short k=0;k<TLots;k++)
				{
					if(k!=sublot[CriticalSublot].parentlot)
					{
						IloExpr temp(workerEnv);
						unsigned short v_cr=sublot[CriticalSublot].parentvendor;

						temp+=piy_ki_z[i][k]+pi1[i][k]-pi2[i][k]+pi3[i][k]
						+(lvdata[k][v_cr].vendortime
							-lvdata[k][v_cr].assemblytime)*pi_c[i]
						+(lvdata[k][v_cr].vendortime
							-lvdata[k][v_cr].assemblytime)*pi_cb[i];

						wk_constraints.add(IloRange(workerEnv,0,temp,IloInfinity));
						temp.end();
					}

				}
			}
		}

		
		for(unsigned short i=0; i<NumSublots;++i)
		{
			if(i!=CriticalSublot)
			{
				

				for(unsigned short k=0;k<TLots;k++)
				{
					if(k!=sublot[i].parentlot)
					{
						IloExpr temp(workerEnv);
						unsigned short u=sublot[i].parentvendor;

						temp+=piy_kj_z[i][k]+pi4[i][k]-pi5[i][k]+pi6[i][k]-(lvdata[k][u].vendortime-lvdata[k][u].assemblytime)*pi_c[i]-(lvdata[k][u].vendortime-lvdata[k][u].assemblytime)*pi_cb[i];

						wk_constraints.add(IloRange(workerEnv,0,temp,IloInfinity));
						temp.end();
					}

				}
			}
		}


		workerMod.add(wk_constraints);





		//cplex

		workerCplex.extract(workerMod);
		IloTimer worker_time(workerEnv);
		workerCplex.setParam(IloCplex::PreInd, IloFalse);
		worker_time.start();
		if(workerCplex.solve())
			cout<<"Worker solved"<<endl;
		cout<<workerCplex.getStatus()<<endl;
		cout<<"Time taken="<<worker_time.stop()<<endl;

		//print values
		cout<<"pi_tc="<<workerCplex.getValue(pi_tc)<<endl;

		for(unsigned short i=0; i<NumSublots;++i)
		{
			if(i!=CriticalSublot)
				cout<<"pi_c["<<i<<"]="<<workerCplex.getValue(pi_c[i])<<"\t";
		}
		cout<<endl;
		for(unsigned short i=0; i<NumSublots;++i)
		{
			if(i!=CriticalSublot)
				cout<<"pi_cb["<<i<<"]="<<workerCplex.getValue(pi_cb[i])<<"\t";
		}
		for(unsigned short i=0; i<NumSublots;++i)
		{
			cout<<endl;			
			for(unsigned short j=0; j<TLots;++j)
				if(i!=CriticalSublot && j != sublot[CriticalSublot].parentlot)
				cout<<"pi1["<<i<<"]="<<workerCplex.getValue(pi1[i][j])<<"\t";
		}

		//IloFloatVar gamma(workerEnv,0,IloInfinity);

	//FloatVarArray1 pi1(workerEnv,TLots);
	//FloatVarArray1 pi2(workerEnv,TLots);
	//FloatVarArray1 pi3(workerEnv,TLots);

	//FloatVarArray3 omega1(workerEnv,TLots);
	//FloatVarArray3 omega2(workerEnv,TLots);
	//FloatVarArray3 omega3(workerEnv,TLots);

	//FloatVarArray3 Chat(workerEnv,TLots);

	//FloatVarArray4 theta1(workerEnv,TLots);
	//FloatVarArray4 theta2(workerEnv,TLots);
	//FloatVarArray4 theta3(workerEnv,TLots);

	//FloatVarArray4 phi1(workerEnv,TLots);
	//FloatVarArray4 phi2(workerEnv,TLots);
	//FloatVarArray4 phi3(workerEnv,TLots);

	//IloNumVar ulamda(workerEnv,0,IloInfinity);
	//FloatVarArray3 uz(workerEnv, TLots);
	//FloatVarArray3 uzlamda(workerEnv, TLots);
	//FloatVarArray1 uylamda(workerEnv,TLots);
	//FloatVarArray4 uzyi(workerEnv, TLots);
	//FloatVarArray4 uzyj(workerEnv, TLots);

	//
	

	////size and assign range to the decision variables:
	//for(int i=0;i<TLots;i++)//over lots
	//{
	//	omega1[i] = FloatVarArray2(workerEnv,NumVendors);
	//	omega2[i] = FloatVarArray2(workerEnv,NumVendors);
	//	omega3[i] = FloatVarArray2(workerEnv,NumVendors);

	//	Chat[i]  = FloatVarArray2(workerEnv,NumVendors);

	//	theta1[i]   = FloatVarArray3(workerEnv,NumVendors);
	//	theta2[i]   = FloatVarArray3(workerEnv,NumVendors);
	//	theta3[i]   = FloatVarArray3(workerEnv,NumVendors);


	//	phi1[i]   = FloatVarArray3(workerEnv,NumVendors);
	//	phi2[i]  = FloatVarArray3(workerEnv,NumVendors);
	//	phi3[i]   = FloatVarArray3(workerEnv,NumVendors);

	//	uz[i]  = FloatVarArray2(workerEnv,NumVendors);
	//	uzlamda[i]  = FloatVarArray2(workerEnv,NumVendors); 

	//	uzyi[i]  = FloatVarArray3(workerEnv,NumVendors);
	//	uzyj[i] = FloatVarArray3(workerEnv,NumVendors);

	//	for(int j=0;j<NumVendors;j++)//over vendors
	//	{
	//		unsigned short maxsublots_t = lvdata[i][j].maxsublots;
	//		omega1[i][j] = FloatVarArray1(workerEnv,maxsublots_t,0,IloInfinity, ILOFLOAT);
	//		omega2[i][j] = FloatVarArray1(workerEnv,maxsublots_t,0,IloInfinity, ILOFLOAT);
	//		omega3[i][j] = FloatVarArray1(workerEnv,maxsublots_t,0,IloInfinity, ILOFLOAT);

	//		Chat[i][j]  = FloatVarArray1(workerEnv,maxsublots_t,0,IloInfinity, ILOFLOAT);

	//		theta1[i][j]   = FloatVarArray2(workerEnv,maxsublots_t);
	//		theta2[i][j]   = FloatVarArray2(workerEnv,maxsublots_t);
	//		theta3[i][j]   = FloatVarArray2(workerEnv,maxsublots_t);


	//		phi1[i][j]   = FloatVarArray2(workerEnv,maxsublots_t);
	//		phi2[i][j]  = FloatVarArray2(workerEnv,maxsublots_t);
	//		phi3[i][j]   = FloatVarArray2(workerEnv,maxsublots_t);

	//		uz[i][j]  = FloatVarArray1(workerEnv,maxsublots_t,0,IloInfinity, ILOFLOAT);
	//		uzlamda[i][j]  = FloatVarArray1(workerEnv,maxsublots_t,0,IloInfinity, ILOFLOAT); 

	//		uzyi[i][j]  = FloatVarArray2(workerEnv,maxsublots_t);
	//		uzyj[i][j] = FloatVarArray2(workerEnv,maxsublots_t);

	//		for(int k=0;k<maxsublots_t;k++)//over sublots
	//		{
	//			theta1[i][j][k]   = FloatVarArray1(workerEnv,TLots,0,IloInfinity, ILOFLOAT);
	//		theta2[i][j][k]   = FloatVarArray1(workerEnv,TLots,0,IloInfinity, ILOFLOAT);
	//		theta3[i][j][k]   = FloatVarArray1(workerEnv,TLots,0,IloInfinity, ILOFLOAT);


	//		phi1[i][j][k]   = FloatVarArray1(workerEnv,TLots,0,IloInfinity, ILOFLOAT);
	//		phi2[i][j][k]  = FloatVarArray1(workerEnv,TLots,0,IloInfinity, ILOFLOAT);
	//		phi3[i][j][k]   = FloatVarArray1(workerEnv,TLots,0,IloInfinity, ILOFLOAT);

	//		uzyi[i][j][k]  = FloatVarArray1(workerEnv,TLots,0,IloInfinity, ILOFLOAT);
	//		uzyj[i][j][k] = FloatVarArray1(workerEnv,TLots,0,IloInfinity, ILOFLOAT);

	//		}}}
		

	*/
};

//double datalap::benders2(void)
//{
//
//	datalap::master_instance = this;
//	double master_optimum_value;
//	IloEnv masterEnv;
//	unsigned short numsublots = submap.size();
//	try
//	{
//		IloNumVar Obj(masterEnv,0);
//
//		//Y-- rational
//		BoolVarArray2 y(masterEnv, TLots);
//		for(unsigned short int i=0;i<TLots;i++)
//			y[i] =BoolVarArray1(masterEnv, TLots);
//
//		BoolVarArray1 lambdap(masterEnv,numsublots);
//
//		IloObjective masterObj(masterEnv,Obj,IloObjective::Minimize);
//
//		IloNum P=0;
//		for(unsigned short int i=0;i<TLots;i++)
//		{
//			P+=lvdata[i][0].assemblytime;
//
//		}
//
//		//Min Transfer Cost
//		IloNum MinTransfer=0;
//		for(unsigned short int i=0;i<TLots;i++)
//		{
//			for(unsigned short j=0;j<NumVendors;j++)
//			{
//				MinTransfer+=lvdata[i][j].handlingcost[0];
//			}
//
//		}
//
//		ObjConstant= MinTransfer+P;
//		cout<<ObjConstant<<" constant after generate data \n";
//
//		IloConstraintArray constraints(masterEnv);
//		IloRangeArray ranger(masterEnv);
//		//////////////////////////////////////////////////////////////////////////////
//		//Create constraints
//
//		IloExpr temp(masterEnv);
//		for(unsigned int short i=0;i<numsublots;i++)
//			temp+=lambdap[i];
//		constraints.add(temp==1);
//
//		for(unsigned int short i=0;i<TLots;i++)
//		{
//			constraints.add(y[i][i]==0);
//
//			for(unsigned int short j=0;j<i;j++)
//				ranger.add(y[i][j]+y[j][i]==1);				
//		}
//
//
//		for(unsigned int short i=0;i<TLots;i++)
//		{
//
//			for(unsigned int short j=0;j<TLots;j++)
//			{
//				for(unsigned int short k=0;k<TLots;k++)
//				{					 
//					if(i!=j && j!=k && k!=i)
//						constraints.add(y[i][j]+y[j][k]+y[k][i]<=2);
//				}}}
//
//
//
//
//
//
//		//////////////////////////////////////////////////////////////////////
//		//After building model basics
//
//		//create models
//		//MASTER
//		IloModel masterMod(masterEnv, "master");
//		masterMod.add(constraints);
//		masterMod.add(ranger);
//		masterMod.add(masterObj);
//
//		IloCplex masterCplex(masterMod);
//		masterCplex.setParam(IloCplex::PreInd, IloFalse); 
//		masterCplex.setParam(IloCplex::Threads, 1); 
//		masterCplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
//
//		masterCplex.setOut(masterEnv.getNullStream());
//		//masterCplex.use(BendersUserCallback(me));
//		masterCplex.use(BendersLazyCallback2_given(masterEnv,y,lambdap,Obj,ObjConstant));
//		masterCplex.use(BendersLazyCallback2_cri(masterEnv,y,lambdap,Obj,ObjConstant));
//
//
//
//
//
//
//
//
//
//
//		/////////////////////////////////////////
//		//SOLVE AND POST PROCESS
//		IloTimer cputime(masterEnv);
//		cputime.start();
//		if(masterCplex.solve())
//		{
//			master_optimum_value = masterCplex.getObjValue() ;
//			m_cputime=cputime.stop();
//			IloAlgorithm::Status solStatus= masterCplex.getStatus();
//			masterEnv.out() << endl << "Solution status: " << solStatus << endl;
//
//
//
//
//
//
//			cout<<"Transitivity Matrix:\n";
//			for(unsigned int short i=0;i<TLots;i++)
//			{
//
//				for(unsigned int short j=0;j<TLots;j++)
//				{
//					cout<<masterCplex.getValue(y[i][j])<<" ";
//
//				}
//				cout<<"\n";
//
//			}
//
//			cout<<"Objective of Master Problem = "<<masterCplex.getObjValue()<<"\n"<<"Cpu time= "<<m_cputime<<"\n";
//
//
//
//		}
//
//		else {
//			masterEnv.out() << "No solution available" << endl;
//		}
//
//
//
//
//
//
//	}
//
//	catch (const IloException& e) {
//		cerr << "Exception caught: " << e << endl;
//
//
//	}
//	catch (...) {
//		cerr << "Unknown exception caught!" << endl;
//	}
//
//
//
//	//wrap up and free memory
//
//	return master_optimum_value;
//
//
//}


//double datalap::benders(void)
//{
//
//	datalap::master_instance = this;
//	double master_optimum_value;
//	IloEnv masterEnv;
//	IloEnv workerEnv;
//	vector<vector<int>> ConstrIndex(TLots);
//	for(unsigned short int i=0;i<TLots;i++)
//		ConstrIndex[i].resize(TLots);
//
//	try
//	{
//			////Variable definition
//		
//
//		// Objective of master problem--rational
//
//		IloNumVar Obj(masterEnv,0);
//
//		//Y-- rational
//		BoolVarArray2 y(masterEnv, TLots);
//		for(unsigned short int i=0;i<TLots;i++)
//			y[i] =BoolVarArray1(masterEnv, TLots);
//
//		//X-- binary
//		BoolVarArray2 x(masterEnv, TLots);
//		for(unsigned short int i=0;i<TLots;i++)
//			x[i] =BoolVarArray1(masterEnv, TLots);
//
//		IloObjective masterObj(masterEnv,Obj,IloObjective::Minimize);
//
//		IloNum P=0;
//		for(unsigned short int i=0;i<TLots;i++)
//		{
//			P+=lvdata[i][0].assemblytime;
//
//		}
//
//		//Min Transfer Cost
//		IloNum MinTransfer=0;
//		for(unsigned short int i=0;i<TLots;i++)
//		{
//			for(unsigned short j=0;j<NumVendors;j++)
//			{
//				MinTransfer+=lvdata[i][j].handlingcost[0];
//			}
//
//		}
//
//		ObjConstant= MinTransfer+P;
//		
//		
//
//
//		IloConstraintArray constraints(masterEnv);
//		IloRangeArray ranger(masterEnv);
//		//////////////////////////////////////////////////////////////////////////////
//		//Create constraints
//		
//		for(unsigned int short i=0;i<TLots;i++)
//		{
//			constraints.add(y[i][i]==0);
//			
//			for(unsigned int short j=0;j<i;j++)
//			{
//
//				ranger.add(y[i][j]+y[j][i]==1);
//				
//				constraints[constraints.getSize()-1].setObject(&ConstrIndex[i][j]);//////////////
//			}}
//		
//
//	
//
//		bool lift=false;
//		if(lift)
//		{
//			//create permutation
//			for(unsigned int short i=0;i<TLots;i++)
//			{
//				IloExpr sumofrow(masterEnv),sumofcolumn(masterEnv);
//				for(unsigned int short j=0;j<TLots;j++)
//				{
//					sumofrow+=x[i][j];
//					sumofcolumn+=x[j][i];
//
//					char varNameyj[100];
//						sprintf(varNameyj, "y.%d.%d", (int) i, (int) j); 
//						y[i][j].setName(varNameyj);
//
//
//				}
//				constraints.add(sumofcolumn==1);
//				constraints.add(sumofrow==1);
//
//			}
//
//			//match permutation to transitive
//			for(unsigned int short i=0;i<TLots;i++)
//			{
//
//				for(unsigned int short j=0;j<TLots;j++)
//				{
//					if(i!=j)
//					{
//
//						for(unsigned int short k=0;k<TLots;k++)
//						{
//							IloExpr prec_to_perm(masterEnv);
//							prec_to_perm+=x[i][k]-1-y[i][j];
//
//							for(unsigned int short bfr_k=k;bfr_k<TLots;bfr_k++)
//							{
//								prec_to_perm+=x[j][bfr_k];
//
//							}
//
//							constraints.add(prec_to_perm<=0);
//
//						}
//
//
//					}
//
//					
//
//				}
//				
//			}
//
//			
//
//
//
//
//		}
//		else
//		{
//				//create transitive constraints
//
//			for(unsigned int short i=0;i<TLots;i++)
//			{
//
//				for(unsigned int short j=0;j<TLots;j++)
//				{
//					for(unsigned int short k=0;k<TLots;k++)
//					{					 
//						if(i!=j && j!=k && k!=i)
//							constraints.add(y[i][j]+y[j][k]+y[k][i]<=2);
//					}}}
//		}
//
//
//		/*for(unsigned int short i=0;i<TLots;i++)
//			{
//
//				for(unsigned int short j=0;j<TLots;j++)
//				{
//					constraints.add(y[i][j]==precedence[i][j]);
//				}}*/
//
//
//
//
//		//////////////////////////////////////////////////////////////////////
//		//After building model basics
//
//		//create models
//		//MASTER
//		IloModel masterMod(masterEnv, "master");
//		masterMod.add(constraints);
//		masterMod.add(ranger);
//		masterMod.add(masterObj);
//		
//		IloCplex masterCplex(masterMod);
//		masterCplex.setParam(IloCplex::PreInd, IloFalse); 
//		masterCplex.setParam(IloCplex::Threads, 1); 
//		masterCplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
//		
//		masterCplex.setOut(masterEnv.getNullStream());
//		//masterCplex.use(BendersUserCallback(me));
//		masterCplex.use(BendersLazyCallback(masterEnv,y,x,Obj,ObjConstant));
//
//		//WORKER
//
//		//Num Sublots
//		unsigned short int NumSublots=0;
//		vector<lvsubdata> sublot;
//		vector<int> tau_star;
//		
//		for(unsigned short int i=0;i<TLots;i++)
//		{
//			for(unsigned short int j=0;j<NumVendors;j++)
//			{
//				unsigned short maxsublots=lvdata[i][j].maxsublots;
//				NumSublots+=maxsublots;
//				for(unsigned short k=0;k<maxsublots;k++)
//				{
//					lvsubdata temp;
//					temp.parentlot=i;
//					temp.parentvendor=j;
//					temp.sublotsused=k;
//					temp.parent= &lvdata[i][j];
//					sublot.push_back(temp);
//
//					if(k==maxsublots-1)
//						tau_star.push_back(0);
//					else
//						tau_star.push_back((lvdata[i][j].handlingcost[k+1]-lvdata[i][j].handlingcost[k]));
//
//				}
//			}
//		}
//		
//		
//
//
//
//		IloModel workerMod(workerEnv, "worker");
//		IloCplex workerCplex(workerEnv); 
//		
//
//		
//
//
//
//		/////////////////////////////////////////
//		//SOLVE AND POST PROCESS
//		IloTimer cputime(masterEnv);
//		cputime.start();
//		if(masterCplex.solve())
//		{
//			master_optimum_value = masterCplex.getObjValue() ;
//			m_cputime=cputime.stop();
//			IloAlgorithm::Status solStatus= masterCplex.getStatus();
//			masterEnv.out() << endl << "Solution status: " << solStatus << endl;
//
//			masterEnv.out() << "Objective value: "
//				<< masterCplex.getObjValue() << endl;
//
//			
//
//				// Write out the optimal solution
//				/*cout<<"Permutaion Matrix:\n";
//				for(unsigned int short i=0;i<TLots;i++)
//				{
//
//					for(unsigned int short j=0;j<TLots;j++)
//					{
//						cout<<masterCplex.getValue(x[i][j])<<" ";
//
//					}
//					cout<<"\n";
//
//				}*/
//
//				cout<<"Transitivity Matrix:\n";
//				for(unsigned int short i=0;i<TLots;i++)
//				{
//
//					for(unsigned int short j=0;j<TLots;j++)
//					{
//						cout<<masterCplex.getValue(y[i][j])<<" ";
//
//					}
//					cout<<"\n";
//
//				}
//
//				cout<<"Objective of Master Problem = "<<masterCplex.getObjValue()<<"\n"<<"Cpu time= "<<m_cputime<<"\n";
//
//
//
//		}
//			
//		else {
//			masterEnv.out() << "No solution available" << endl;
//		}
//
//
//
//
//
//
//	}
//
//	catch (const IloException& e) {
//		cerr << "Exception caught: " << e << endl;
//
//	
//	}
//	catch (...) {
//		cerr << "Unknown exception caught!" << endl;
//	}
//
//
//
//	//wrap up and free memory
//	masterEnv.end();
//	workerEnv.end();
//	return master_optimum_value;
//
//
//}