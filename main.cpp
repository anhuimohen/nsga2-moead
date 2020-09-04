#include"question.h"
#include"evaluation.h"
#include"NSGA2.h"
#include<ctime>
#include<stdlib.h>
vector< vector<double> > ws;
#define DD 200															//��������

int main()
{
	srand((unsigned int)time(NULL));
	///////////////////////////////////
	//  ���ľ��߿ռ��Сʱ �Ƚ��г�ʼ��д���ļ�
	/*
	decision_space spa;
	spa.init();															//���߿ռ�ĳ�ʼ��
	spa.print();														//���߿ռ��ӡ
	spa.write_true_value();												//�����߿ռ�д��txt�ļ�
	spa.write_normalization_value();
	*/

	/*
	///////////////////////////////////
	//  �㷨IGDֵ�ļ���
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
	//  ��ȡ���߿ռ�����
	decision_space space;
	space.read_true_value();
	space.read_normalization_value();
	space.print();


	population pop;
	pop.init(space);													//��Ⱥ��ʼ��
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
		cout << "��" << dd << "�ε�����+++++++++++++++++++";

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



