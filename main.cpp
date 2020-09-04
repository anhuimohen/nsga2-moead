#include"question.h"
#include"evaluation.h"
#include"NSGA2.h"
#include<ctime>
#include<stdlib.h>
vector< vector<double> > ws;
#define DD 200															//迭代次数

int main()
{
	srand((unsigned int)time(NULL));
	///////////////////////////////////
	//  更改决策空间大小时 先进行初始化写入文件
	/*
	decision_space spa;
	spa.init();															//决策空间的初始化
	spa.print();														//决策空间打印
	spa.write_true_value();												//将决策空间写入txt文件
	spa.write_normalization_value();
	*/

	/*
	///////////////////////////////////
	//  算法IGD值的计算
	vector<individual> select_solution;
	vector<individual> NSGA2_solution;

	vector<individual> pvariation_solution;
	vector<individual> pvariation_select_moead;
	vector<individual> pvariation_select_moead02;
	vector<individual> subpop;
	read_solution(".\\evalution\\aes-5-200_EP_value_.xls",2002,aesmoead_solution);
	read_solution(".\\evalution\\moead-5-200_EP_value_.xls",2603, moead_solution);
	read_solution(".\\evalution\\select-moead--5-200_EP_value_.xls",3348, select_solution);
	read_solution(".\\evalution\\NSGA2-5-200_EP_value_.xls",787, NSGA2_solution);
	read_solution(".\\evalution\\subpop_200_EP_value_.xls",2925, subpop);
	//add_solution(aesmoead_solution,non_dominated_solution);
	add_solution(moead_solution, non_dominated_solution);
	add_solution(subpop, non_dominated_solution);
	//add_solution(select_solution, non_dominated_solution);
	//(NSGA2_solution, non_dominated_solution);
	cout << cal_IGD(moead_solution, non_dominated_solution) << endl;
	cout << cal_IGD(subpop, non_dominated_solution) << endl;
	//cout << cal_IGD(aesmoead_solution, non_dominated_solution) << endl;
	//cout << cal_IGD(select_solution, non_dominated_solution) << endl;
	//cout << cal_IGD(NSGA2_solution, non_dominated_solution) << endl;

	cout << cal_coverage(moead_solution, subpop) << endl;
	cout << cal_coverage(subpop, moead_solution)<<endl;
	cout << cal_coverage(moead_solution, aesmoead_solution) << endl;
	cout << cal_coverage(aesmoead_solution, moead_solution) << endl;
	cout << cal_coverage(aesmoead_solution, select_solution) << endl;
	cout << cal_coverage(select_solution, aesmoead_solution) << endl;
	cout << cal_coverage(pvariation_solution, moead_solution) << endl;
	cout << cal_coverage(moead_solution, pvariation_solution) << endl;
	*/

	//////////////////////////////////
	//  读取决策空间数据
	decision_space space;
	space.read_true_value();
	space.read_normalization_value();
	space.print();


	population pop;
	pop.init(space);													//种群初始化
	cout << pop.mean_vector(objective, 3);
	pop.subpop_assign();
	pop.subpop_nondominant_sort();
	pop.cal_nei_index();
	pop.make_offspring(space);
	pop.merge();
	pop.update_ep_nsga_moead();
	cal_z_min(pop);
	pop.subpop_nondominant_sort();
	pop.eliminate();







	


	//save_pop_data(pop, 1);
	int dd = 1;
	do
	{
		cout << "第" << dd << "次跌代：+++++++++++++++++++";

		pop.make_offspring(space);
		pop.merge();
		pop.update_ep_nsga_moead();
		cal_z_min(pop);
		pop.subpop_nondominant_sort();
		pop.eliminate();
		pop.update_deerta(dd, 10);
		cal_mean_objective();
		
		if (dd % 25 == 0)
		{
			save_pop_data(pop, dd);
			save_ep_data(dd);
		}
		
		dd++;
	} while (dd < DD + 1);
	save_mean_objective();
	save_deerta();

	cout << "";
	pop.make_offspring(space);


	pop.merge();
	cal_z_min(pop);
	pop.subpop_nondominant_sort();
	pop.eliminate();
  }



