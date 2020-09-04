#pragma once
#include"question.h"
#include <algorithm>
#include <fstream>
#include<math.h>
#include<string>
#include<stdlib.h>
#include<iostream>

#define popsize 126                                             //��Ⱥ��ģ
#define num_T  4                                                 //�����ھӵĸ���
#define p_variate 0.2											 //�������
#define SBX_n 2													 //�ֲ�ָ�� nԽ���ʾ�Ӵ��͸������ӽ�
#define apap_x 400												 //����Ӧ�������
#define DE_F 0.5												 //DE���ӵ���������
#define num_vector 35											//�ο������ĸ���

class individual;
class populatioon;
double z_min[objective];										 //Z* ��Сֵ
vector<int> num_update;										 //��Ⱥ�и���ĸ��´���
vector<double> mean_objective[objective];						 //��¼���Ž������������ӷ�֧����ȺĿ��ֵ�ľ�ֵ
vector<individual> EP;											 //�ⲿ��Ⱥ
double p_select_concrete[dimension][concrete_service_num];		 //��¼��ǰEP��ÿάÿ�����������ֵĸ���

vector<individual> pop_temporary;									 //�����Ҫ������ӵ����ѡ����ǲ����
vector<double> deerta_data;										//ÿ���ĵ¶�������

class individual
{
public:
	int select[dimension];										 //ÿһάѡ��ķ�����
	double value[objective];                                     //������Ŀ��ֵ
	void random_select();										 //������ɽ�
	void cal_value(decision_space a);
	void print();

	double crowded_distace[objective + 1];						 //�����ӵ�����룬���һ�����ܺ�
	
	///////////////////////////////////
	// С��Ⱥ
	int source;													 //�����Դ 2-��ʼ���Ľ�  0-�����ھӣ�1-�����ھ� 3-�������		
	int rank;													 //����������Ⱥ�ķ�֧��ȼ�
	double value_g;												 //�������ڵ�ǰ����Ⱥ���б�ѩ����ֵ
};

void individual::random_select()								 //������ɽ�
{
	for (int i = 0; i < dimension; i++)
	{
		select[i] = rand() % concrete_service_num;
	}
}

void individual::cal_value(decision_space a)						     //�����ʼ��
{
	for (int i = 0; i < objective; i++)
	{
		if (ob_cal_mode[i] == 0)								 //0-��� ��ʼΪ0
			value[i] = 0;
		if (ob_cal_mode[i] == 1)                                 //1-��� ��ʼΪ1
			value[i] = 1;
	}
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < objective; j++)
		{
			if (ob_cal_mode[j] == 0)
				value[j] += a.service_unit[i][select[i]].attribute_value[j][1];
			if (ob_cal_mode[j] == 1)
				value[j] = value[j] * a.service_unit[i][select[i]].attribute_value[j][1];
		}
	}
}

void individual::print()                                         //��ӡ����
{
	cout << "����ÿάѡ����ǣ�";
	for (int i = 0; i < dimension; i++)
	{
		cout << setw(5) << select[i];
	}
	cout << "Ŀ��ֵ�ǣ�";
	for (int i = 0; i < objective; i++)
	{
		cout << setw(15) << value[i];
	}
}

class population
{
public:
	individual pop[popsize];
	double value_g[popsize];									 //�����б�ѩ��ʽ�����Ŀ��ֵ g
	vector< vector<double> > ws;								 //Ȩ������
	int neighbor_index[num_vector][num_T];						 //ÿ�������ھ��������
	int mean_vector(int m, double H);						     //���ɾ��ȷֲ�������
	void vector_print();
	void init(decision_space a);
	void print();
	void cal_nei_index();
	void cal_value_g();											 //�б�ѩ��ʽ�����Ŀ��ֵ g
	void update(decision_space space, int dd);					 //����������Ⱥ����
	void select_parents(int k, individual& z, individual& x);

	int nondominant_rank[popsize];								 //ÿ������ķ�֧��ȼ�
	individual pop_Q[popsize];
	individual pop_R[popsize * 2];
	int R_nondominant_rank[popsize * 2];
	void nondominant_sort();									 //��Ⱥ��֧������
	void nondominant_select(decision_space space);				 //��֧��ѡ��
	void set_pop_R();
	void fast_nondominant_sort();
	void make_new_pop();										//NSGAѡ�񸸴�
	void select_individual(individual& a, individual& b);
	void make_new_pop_Q(decision_space space);

	//////////////////////////////////////////////////////////////
	vector<individual> subpop[num_vector];
	vector<individual> offspring_pop[num_vector];
	double deerta;												//���ھ���Ⱥѡ�񸸴�����������Ⱥѡ������Ⱥ��ʼ��Ϊ0.3
	int distance_subpop_index(int x);							//�������������Ⱥi��Զ����Ⱥ������
	//double z_subpop[num_vector][objective];						//ÿ������Ⱥ�Ĳο��㣬����͵�
	void subpop_assign();										//��������������Ⱥ
	void make_offspring(decision_space space);
	void subpop_nondominant_sort();
	void individual_assign_offspring_pop(individual x);				//��������䵽����Ⱥ��
	void merge();												//��subpop��offspring_pop�ϲ���һ��
	void eliminate();											//���ϲ������Ⱥ��С����Ϊpopsize
	void update_ep_nsga_moead();
	void update_deerta(int dd,int L);
};

bool operator<(individual a, individual b);
void individual_cross(individual& a, individual& b);
void variation(individual& a, decision_space space);
void population::nondominant_sort()
{
	int mark[popsize];											 //������飬0-δ���ּ���1-�Ѿ��ּ�
	memset(mark, 0, sizeof(mark));
	memset(nondominant_rank, 0, sizeof(nondominant_rank));
	int rank = 1;
	for (int k = 0; k < popsize; k++)
	{
		for (int i = 0; i < popsize; i++)
		{
			if (nondominant_rank[i] != 0)
				continue;
			else
			{
				int flag = 1;												//1-�ø����Ƿ�֧��ĸ���
				for (int j = 0; j < popsize; j++)
				{
					if (mark[j] == 1 || i == j)
						continue;
					else
					{
						if (pop[j] < pop[i])
						{
							flag = 0;
							break;
						}
					}
				}
				if (flag)
				{
					nondominant_rank[i] = rank;
				}
			}
		}
		for (int i = 0; i < popsize; i++)
		{
			if (nondominant_rank[i] == rank)
				mark[i] = 1;
		}
		rank++;
	}
}

void population::nondominant_select(decision_space space)
{
	individual a, b;
	double fitness[popsize];
	int sum = 0;
	int rank = 0;
	for (int i = 0; i < popsize; i++)
	{
		if (nondominant_rank[i] > rank)
			rank = nondominant_rank[i];
	}
	for (int i = 0; i < popsize; i++)
	{
		fitness[i] = double(rank + 1 - nondominant_rank[i]);
		sum += fitness[i];
	}
	fitness[0] = fitness[0] / sum;
	for (int i = 1; i < popsize; i++)
	{
		fitness[i] = fitness[i] / sum;
		fitness[i] = fitness[i] + fitness[i - 1];
	}
	for (int i = 0; i / 2 < popsize / 2; i = i + 2)
	{
		double shu = rand() % 100 / 100.0;
		for (int j = 0; j < popsize; j++)
		{
			if (shu < fitness[j])
			{
				a = pop[j];
				break;
			}
		}
		double shu2 = rand() % 100 / 100.0;
		for (int j = 0; j < popsize; j++)
		{
			if (shu2 < fitness[j])
			{
				b = pop[j];
				break;
			}
		}
		individual_cross(a, b);
		variation(a, space);
		variation(b, space);
		pop_Q[i] = a;
		pop_Q[i + 1] = b;
	}
	if (popsize % 2 == 1)
	{
		double shu = rand() % 100 / 100.0;
		for (int j = 0; j < popsize; j++)
		{
			if (shu < fitness[j])
			{
				a = pop[j];
				break;
			}
		}
		pop_Q[popsize - 1] = a;
	}
}

void population::init(decision_space a)
{
	deerta = 0.5;
	for (int i = 0; i < objective; i++)							//��ʼz*_min													
	{
		z_min[i] = RAND_MAX;
	}
	for (int i = 0; i < popsize; i++)
	{
		pop[i].source = 2;
		pop[i].random_select();
		pop[i].cal_value(a);
	}

}

void population::print()
{
	cout << "��ӡ��Ⱥ��*******************" << endl;
	for (int i = 0; i < popsize; i++)
	{
		pop[i].print();
		cout << "   value_g  " << value_g[i];
		cout << endl;
	}
}

void cal_z_min(population a)                                    //�������Ŀ�����Сֵ ��z*
{
	for (int i = 0; i < objective; i++)							//��ʼz*_min													
	{
		z_min[i] = RAND_MAX;
	}
	for (int k = 0; k < num_vector; k++)
	{
		for (int i = 0; i < a.subpop[k].size(); i++)
		{
			for (int j = 0; j < objective; j++)
			{
				if (a.subpop[k][i].value[j] < z_min[j])
					z_min[j] = a.subpop[k][i].value[j];
			}
		}
	}
}

void z_min_update(individual a)									//���ø���a����Z*
{
	for (int i = 0; i < objective; i++)
	{
		if (a.value[i] < z_min[i])
			z_min[i] = a.value[i];
	}
}

void z_min_print()
{
	cout << "��ӡ����Ŀ�����СֵZ*��";
	for (int i = 0; i < objective; i++)
	{
		cout << setw(15) << z_min[i];
	}
	cout << endl;
}

int population::mean_vector(int m, double H)//��mάĿ��ռ��ÿά�о��Ȳ���H������
{
	/*
	int m; // the number of objectives
	double stepsize;
	double H; // H = 1 / stepsize
	cout << "Please input the number of objectives (m): \n";
	cin >> m;
	cout << "Please input the stepsize (1/H): \n";
	cin >> stepsize;
	H = 1 / stepsize;
	cout << "H = " << H << endl;
	*/
	vector<int> sequence;
	for (unsigned int i = 0; i < H; i++) // the number of zero is (H)
	{
		sequence.push_back(0);
	}
	for (unsigned int i = 0; i < (m - 1); i++) // the number of 1 is (H + m - 1 - (m - 1))
	{
		sequence.push_back(1);
	}

	do
	{
		int s = -1;
		vector<double> weight;
		for (unsigned int i = 0; i < (H + m - 1); i++)
		{
			if (sequence[i] == 1)
			{
				double w = i - s;
				w = (w - 1) / H;
				s = i;
				weight.push_back(w);
			}
		}
		double w = H + m - 1 - s;
		w = (w - 1) / H;
		weight.push_back(w);
		ws.push_back(weight);
	} while (next_permutation(sequence.begin(), sequence.end()));
	ofstream outfile("weight.txt");
	for (unsigned int i = 0; i < ws.size(); i++)
	{
		for (unsigned int j = 0; j < ws[i].size(); j++)
		{
			outfile << ws[i][j] << " ";
		}
		outfile << "\n";
	}
	return ws.size();
}

void population::vector_print()
{
	cout << "vector:" << endl;
	for (int i = 0; i < ws.size(); i++)
	{
		for (int j = 0; j < ws[0].size(); j++)
		{
			cout << ws[i][j] << '\t';
		}
		cout << endl;
	}
}

void population::cal_nei_index()                                  //�����������ھ�����������
{
	for (int i = 0; i < num_vector; i++)
	{
		int index[num_vector];
		double distance[num_vector];
		for (int j = 0; j < num_vector; j++)
		{
			index[j] = j;
			double dis = 0;
			for (int k = 0; k < objective; k++)
			{
				dis += pow((ws[i][k] - ws[j][k]), 2);
			}
			distance[j] = sqrt(dis);
		}
		for (int j = 0; j < num_vector; j++)                       //�����밴�մ�С��������
		{														//ͬʱ�������Ҳ����λ��
			for (int k = 0; k < num_vector - 1; k++)
			{
				if (distance[k] > distance[k + 1])
				{
					double di = distance[k];
					distance[k] = distance[k + 1];
					distance[k + 1] = di;
					int in = index[k];
					index[k] = index[k + 1];
					index[k + 1] = in;
				}
			}
		}
		for (int j = 0; j < num_T; j++)
		{
			neighbor_index[i][j] = index[j];
		}
	}
	////////////////////////////////////////////
	// ���ÿ���������ھ�����������
	
	for (int i = 0; i < num_vector; i++)
	{
		cout << "��" << i << "���������ھ��ǣ�";
		for (int j = 0; j < num_T; j++)
		{
			cout << neighbor_index[i][j] << '\t';
		}
		cout << endl;
	}
	
}

void population::cal_value_g()										//�б�ѩ��ʽ�����Ŀ��ֵ g
{
	for (int i = 0; i < popsize; i++)
	{
		double temporary_g = 0;
		for (int j = 0; j < objective; j++)
		{
			double tem = (pop[i].value[j] - z_min[j]) * ws[i][j];
			if (tem > temporary_g)
				temporary_g = tem;
		}
		value_g[i] = temporary_g;
	}
}

void DE(individual& a, individual& b, individual& c)											//DE����
{
	int x = rand() % dimension;
	for (int i = x; i < dimension; i++)
	{
		int x1 = round(a.select[i] + DE_F * (b.select[i] - c.select[i]));
		int x2 = round(b.select[i] + DE_F * (a.select[i] - c.select[i]));
		///////////////////////////////////
		//  ������  �������߽�ȡ�ദ��
		if (x1 < 0)
			x1 = 0 + (0 - x1) % (concrete_service_num - 1);
		if (x1 >= concrete_service_num)
			x1 = concrete_service_num - 1 - (x1 - (concrete_service_num - 1)) % (concrete_service_num - 1);
		if (x2 < 0)
			x2 = 0 + (0 - x2) % (concrete_service_num - 1);
		if (x2 >= concrete_service_num)
			x2 = concrete_service_num - 1 - (x2 - (concrete_service_num - 1)) % (concrete_service_num - 1);
		a.select[i] = x2;
		b.select[i] = x1;
	}
}

void SBX(individual& a, individual& b)											//SBX����
{
	int x = rand() % dimension;
	double u = rand() % 100 / 100.0;
	double beita;
	if (u <= 0.5)
		beita = pow(2 * u, 1 / (SBX_n + 1));
	else
		beita = pow(1 / (2 - 2 * u), 1 / SBX_n + 1);
	for (int i = x; i < dimension; i++)
	{
		int c1 = round(0.5 * (double(a.select[i] + b.select[i])) - 0.5 * beita * (double(b.select[i] - a.select[i])));
		int c2 = round(0.5 * (double(a.select[i] + b.select[i])) + 0.5 * beita * (double(b.select[i] - a.select[i])));
		///////////////////////////////////
		//  ������  �������߽�ȡ�ദ��
		if (c1 < 0)
			c1 = 0 + (0 - c1) % (concrete_service_num - 1);
		if (c1 >= concrete_service_num)
			c1 = concrete_service_num - 1 - (c1 - (concrete_service_num - 1)) % (concrete_service_num - 1);
		if (c2 < 0)
			c2 = 0 + (0 - c2) % (concrete_service_num - 1);
		if (c2 >= concrete_service_num)
			c2 = concrete_service_num - 1 - (c2 - (concrete_service_num - 1)) % (concrete_service_num - 1);
		a.select[i] = c2;
		b.select[i] = c1;
	}
}

void adaptive_variation(individual& a, int dd, decision_space space)			//����Ӧ����  dd��������
{
	double p_mu = 0.2 * (1 - exp(-1 * dd / apap_x));
	int flag = 0;															//�Ƿ�����ͻ�����	1-����
	if (dd > 10)
	{
		int count = 0;
		for (int i = dd - 1; i >= dd - 5; i--)
		{
			if (num_update[i] < 5)											//�������5�����´�����<5���� ����0.1�ı������
				count++;
		}
		if (count == 5)
			flag = 1;
	}
	if (flag)
		p_mu += 0.1;
	for (int i = 0; i < dimension; i++)
	{
		if (rand() % 100 / 100.0 < p_mu)
			a.select[i] = rand() % concrete_service_num;
	}
	a.cal_value(space);
}

void adaptive_variation2(individual& a, int dd, decision_space space)			//����Ӧ����  dd��������
{
	double p_mu = 0.2 * (1 - exp(-1 * dd / apap_x));
	int flag = 0;															//�Ƿ�����ͻ�����	1-����
	if (dd > 10)
	{
		int count = 0;
		for (int i = dd - 1; i >= dd - 5; i--)
		{
			if (num_update[i] < 5)											//�������5�����´�����<5���� ����0.1�ı������
				count++;
		}
		if (count == 5)
			flag = 1;
	}
	if (flag)
		p_mu += 0.1;
	for (int i = 0; i < dimension; i++)
	{
		if (rand() % 100 / 100.0 < p_mu)
		{
			double p1 = rand() % 100 / 100.0;
			for (int j = 0; j < concrete_service_num; j++)
			{
				if (p1 < p_select_concrete[i][j])
				{
					a.select[i] = j;
					break;
				}
			}
		}
	}
	a.cal_value(space);
}

void individual_cross(individual& a, individual& b)								//���彻��
{
	int x1 = rand() % dimension;
	int x2;
	do
	{
		x2 = rand() % dimension;
	} while (x1 == x2 || x1 - x2 == dimension - 1 || x2 - x1 == dimension - 1);
	if (x1 > x2)
	{
		int x3 = x1;
		x1 = x2;
		x2 = x3;
	}
	for (int i = x1; i <= x2; i++)
	{
		int x = a.select[i];
		a.select[i] = b.select[i];
		b.select[i] = x;
	}
}

void cal_p_select_concrete()
{
	memset(p_select_concrete, 0, sizeof(p_select_concrete));
	for (int i = 0; i < EP.size(); i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			p_select_concrete[j][EP[i].select[j]]++;
		}
	}
	for (int i = 0; i < dimension; i++)
	{
		int count = 0;
		for (int j = 0; j < concrete_service_num; j++)
		{
			count += p_select_concrete[i][j];
		}
		for (int j = 0; j < concrete_service_num; j++)
		{
			p_select_concrete[i][j] = p_select_concrete[i][j] / (count * 1.0);
		}
		for (int j = 1; j < concrete_service_num; j++)
		{
			p_select_concrete[i][j] += p_select_concrete[i][j - 1];
		}
	}
}

void p_variation(individual& a, decision_space space)						//���ݸ��ʽ������̶ķ�����
{
	for (int i = 0; i < dimension; i++)
	{
		if (rand() % 100 / 100.0 < p_variate)
		{
			double p1 = rand() % 100 / 100.0;
			for (int j = 0; j < concrete_service_num; j++)
			{
				if (p1 < p_select_concrete[i][j])
				{
					a.select[i] = j;
					break;
				}
			}
		}
	}
	a.cal_value(space);
}

void variation(individual& a, decision_space space)							//������� ����Ŀ��ֵ
{
	for (int i = 0; i < dimension; i++)
	{
		if (rand() % 100 / 100.0 < p_variate)
			a.select[i] = rand() % concrete_service_num;
	}
	a.cal_value(space);
}

double cal_y_g(individual a, vector<double> w)								//�������������w�ϵ�Ŀ��ֵvalue_g
{
	double temporary = 0;
	for (int i = 0; i < objective; i++)
	{
		double tem = (a.value[i] - z_min[i]) * w[i];
		if (tem > temporary)
			temporary = tem;
	}
	return temporary;
}

bool operator<(individual a, individual b)						//����a�Ƿ�֧��b
{
	int falg = 1;												//1-ÿһ��Ŀ��ֵa<=b
	int fa = 0;													//1-����һ��Ŀ��ֵa<b
	for (int i = 0; i < objective; i++)
	{
		if (a.value[i] > b.value[i])
		{
			falg = 0;
			break;
		}
	}
	for (int i = 0; i < objective; i++)
	{
		if (a.value[i] < b.value[i])
		{
			fa = 1;
			break;
		}
	}
	if (falg && fa)
		return true;
	return false;
}

void update_ep(individual a)												//�����ⲿ�浵��ȺEP
{
	////////////////////////////////////////////////
	//  ����ɾ��EP�б�a֧��Ľ�
	for (int i = 0; i < EP.size(); i++)
	{
		if (a < EP[i])
		{
			EP.erase(EP.begin() + i);
			i--;
		}
	}
	//  ���û�н�֧��a,�ͽ�a����EP
	int flag = 1;															//1-�޸���֧��a
	for (int i = 0; i < EP.size(); i++)
	{
		int count = 0;														//��ͬ��Ŀ�����
		if (EP[i] < a)
		{
			flag = 0;
			break;
		}
		for (int j = 0; j < objective; j++)
		{
			if (EP[i].value[j] == a.value[j])
				count++;
		}
		if (count == objective)
		{
			flag = 0;
			break;
		}
	}
	if (flag)
		EP.push_back(a);
}

void ep_print()
{
	cout << "��ӡ�浵��ȺEP------------------------";
	for (int i = 0; i < EP.size(); i++)
	{
		cout << "Ŀ��ֵ��:";
		for (int j = 0; j < objective; j++)
		{
			cout << setw(15) << EP[i].value[j];
		}
		cout << endl;
	}
}

void save_num_upate()
{
	string str1 = ".\\data\\num_update\\";
	string str2 = "num_update.xls";
	str2 = str1 + str2;
	ofstream outfile(str2);
	for (int i = 0; i < num_update.size(); i++)
	{
		outfile << num_update[i] << endl;
	}
	outfile.close();
}

void cal_mean_objective()													//����EP��ȺĿ��ֵ�ľ�ֵ
{
	double ob[objective] = { 0 };
	for (int i = 0; i < objective; i++)
	{
		for (int j = 0; j < EP.size(); j++)
			ob[i] += EP[j].value[i];
		ob[i] = ob[i] / EP.size();
		mean_objective[i].push_back(ob[i]);
	}
}

void save_mean_objective()													//���õ���Ŀ���ֵ��������������ݱ�������
{
	string str1 = ".\\data\\num_update\\";
	string str2 = "EP_objective_mean_dd.xls";
	str2 = str1 + str2;
	ofstream out(str2);
	for (int i = 0; i < mean_objective[0].size(); i++)
	{
		for (int j = 0; j < objective; j++)
		{
			out << mean_objective[j][i] << '\t';
		}
		out << endl;
	}
	out.close();
}

void population::select_parents(int k, individual& z, individual& x)
{																			//ѡ�������ѡ����������н��ŵĸ���
	vector<individual> non_dominant;
	for (int i = 0; i < num_T; i++)
	{
		int flag = 1;														//1-��i���ھӸ����Ƿ�֧���
		for (int j = 0; j < num_T; j++)
		{
			if (i != j)
			{
				if (pop[neighbor_index[k][j]] < pop[neighbor_index[k][i]])
				{
					flag = 0;
					break;
				}
			}
		}
		if (flag)
		{
			int fla = 1;
			for (int q = 0; q < non_dominant.size(); q++)
			{
				int count = 0;
				for (int p = 0; p < objective; p++)
				{
					if (pop[neighbor_index[k][i]].value[p] == non_dominant[q].value[p])
						count++;
				}
				if (count == objective)
				{
					fla = 0; break;
				}
			}
			if (fla)
				non_dominant.push_back(pop[neighbor_index[k][i]]);
		}
	}
	if (non_dominant.size() < 2)
	{
		int x1 = rand() % num_T;											//�����������i�������ھӱ��
		int x2;
		do
		{
			x2 = rand() % num_T;
		} while (x1 == x2);
		int w = neighbor_index[k][x1];
		int l = neighbor_index[k][x2];
		z = pop[w];
		x = pop[l];
	}
	else
	{
		int x1 = rand() % non_dominant.size();
		int x2;
		do
		{
			x2 = rand() % non_dominant.size();
		} while (x1 == x2);
		z = non_dominant[x1];
		x = non_dominant[x2];
	}
}

void population::update(decision_space space, int dd)
{
	int number_update = 0;													//ÿ�ָ���ĸ��´���
	for (int i = 0; i < popsize; i++)
	{
		individual a, b, c;
		//////////////////////////////
		//  ���ѡ���ھ��е���������
		/*
		int x1 = rand() % num_T;											//�����������i�������ھӱ��
		int x2,x3;
		do
		{
			x2 = rand() % num_T;
			x3 = rand() % num_T;
		} while (x1==x2||x1==x3||x2==x3);
		int k = neighbor_index[i][x1];
		int l = neighbor_index[i][x2];
		int j = neighbor_index[i][x3];
		a = pop[k];
		b = pop[l];
		c = pop[j];
		*/
		//////////////////////////////////
		//  ѡ���ھ��еķ�֧�������Ϊ����
		select_parents(i, a, b);
		/*
		double update_rate = 0;												//�������Ƶ��
		if (dd > 10)
		{
			for (int j = dd-1; j >= dd - 5; j--)
			{
				update_rate += num_update[j];
			}
			update_rate = update_rate / 5;
			if (update_rate < 5)						//�����10��ƽ������Ƶ��<5�������SBX���ӣ��������DE����
				SBX(a, b);
			else
				DE(a, b,c);
		}
		else
			DE(a, b,c);
		adaptive_variation(a,dd,space);
		adaptive_variation(b, dd, space);
		*/


		individual_cross(a, b);

		/////////////////////////
		//  ���ݸ������̶ı���  --Ч����
		/*
		cal_p_select_concrete();
		p_variation(a, space);
		p_variation(b, space);
		*/
		//////


		variation(a, space);
		variation(b, space);
		z_min_update(a);
		z_min_update(b);
		cal_value_g();
		for (int j = 0; j < num_T; j++)					//���½����ھ�������Ŀ��ֵ��ԭ����Ŀ��Ƚϣ����¸���
		{
			if (cal_y_g(a, ws[neighbor_index[i][j]]) < value_g[neighbor_index[i][j]])
			{
				number_update++;
				pop[neighbor_index[i][j]] = a;
				value_g[neighbor_index[i][j]] = cal_y_g(a, ws[neighbor_index[i][j]]);
			}
			if (cal_y_g(b, ws[neighbor_index[i][j]]) < value_g[neighbor_index[i][j]])
			{
				number_update++;
				pop[neighbor_index[i][j]] = b;
				value_g[neighbor_index[i][j]] = cal_y_g(b, ws[neighbor_index[i][j]]);
			}
		}
		update_ep(a);
		update_ep(b);
	}
	num_update.push_back(number_update);
	cal_mean_objective();
}

void save_pop_data(population a, int dd)
{
	string str1 = ".\\data\\pop\\";
	string str = "_pop_value_.xls";
	string str2 = to_string(dd);
	str = str1 + str2 + str;
	ofstream outfile(str);
	for (int i = 0; i < popsize; i++)
	{
		for (int j = 0; j < objective; j++)
		{
			outfile << a.pop[i].value[j] << '\t';
		}
		outfile << endl;
	}
	outfile.close();
}

void save_ep_data(int dd)
{
	string str1 = ".\\data\\EP\\";
	string str = "_EP_value_.xls";
	string str2 = to_string(dd);
	str = str1 + str2 + str;
	ofstream outfile(str);
	for (int i = 0; i < EP.size(); i++)
	{
		for (int j = 0; j < objective; j++)
		{
			outfile << EP[i].value[j] << '\t';
		}
		outfile << endl;
	}
	outfile.close();
}

void population::set_pop_R()
{
	for (int i = 0; i < popsize; i++)
	{
		pop_R[i] = pop[i];
	}
	for (int i = 0; i < popsize; i++)
	{
		pop_R[popsize + i] = pop_Q[i];
	}
}

void population::fast_nondominant_sort()
{
	int np[popsize * 2];										//֧��ø���ĸ���
	vector<int> sp[popsize * 2];								//������֧��ĸ��弯��
	for (int i = 0; i < popsize * 2; i++)
	{
		int count = 0;
		for (int j = 0; j < popsize * 2; j++)
		{
			if (pop_R[j] < pop_R[i])
			{
				count++;
			}
			if (pop_R[i] < pop_R[j])
			{
				sp[i].push_back(j);
			}
		}
		np[i] = count;
	}
	int rank = 1;
	int mark[popsize * 2] = { 0 };
	memset(R_nondominant_rank, 0, sizeof(R_nondominant_rank));
	int count = 0;
	do
	{
		//cout << "3333";
		for (int i = 0; i < popsize * 2; i++)
		{
			if (np[i] == 0 && mark[i] == 0)
			{
				R_nondominant_rank[i] = rank;
				mark[i] = 1;
				count++;
			}
		}

		cout << "����ȼ�" << endl;
		for (int i = 0; i < popsize * 2; i++)
		{
			cout << R_nondominant_rank[i];
		}
		cout << endl;
		cout << "***count=" << count << "*****";
		for (int i = 0; i < popsize * 2; i++)
		{
			if (R_nondominant_rank[i] == rank)
			{
				for (int j = 0; j < sp[i].size(); j++)
				{
					np[sp[i][j]]--;
				}
			}
		}
		rank++;
		cout << "np" << endl;
		for (int i = 0; i < popsize * 2; i++)
		{
			cout << np[i];
		}
		cout << endl;
	} while (count < popsize * 2);
}

void paixu(vector<individual>& pop, int t)
{
	cout << "paixu";
	individual a;
	for (int i = 0; i < pop.size(); i++)
	{
		for (int j = 0; j < pop.size() - 1; j++)
		{
			if (pop[j].value[t] > pop[j + 1].value[t])
			{
				a = pop[j + 1];
				pop.erase(pop.begin() + j + 1);
				pop.insert(pop.begin() + j, a);
			}
		}
	}
	cout << "pai xu done" << endl;
}

bool paixusum(individual a, individual b)
{
	return a.crowded_distace[objective] > b.crowded_distace[objective];
}

bool paixu0(individual a, individual b)
{
	return a.value[0] < b.value[0];
}

bool paixu1(individual a, individual b)
{
	return a.value[1] < b.value[1];
}

bool paixu2(individual a, individual b)
{
	return a.value[2] < b.value[2];
}

bool paixu3(individual a, individual b)
{
	return a.value[3] < b.value[3];
}

bool paixu4(individual a, individual b)
{
	return a.value[4] < b.value[4];
}

void cal_crowded_distance(vector<individual>& pop)
{
	cout << "cal crowded distance" << endl;
	for (int i = 0; i < objective; i++)
	{
		double min, max;
		//paixu(pop, i);

		if (i == 0)
			sort(pop.begin(), pop.end(), paixu0);
		if (i == 1)
			sort(pop.begin(), pop.end(), paixu1);
		if (i == 2)
			sort(pop.begin(), pop.end(), paixu2);
		if (i == 3)
			sort(pop.begin(), pop.end(), paixu3);
		if (i == 4)
			sort(pop.begin(), pop.end(), paixu4);
		min = pop[0].value[i];
		max = pop[pop.size() - 1].value[i];
		pop[0].crowded_distace[i] = RAND_MAX;
		pop[pop.size() - 1].crowded_distace[i] = RAND_MAX;
		for (int j = 1; j < pop.size() - 1; j++)
		{
			pop[j].crowded_distace[i] = (pop[j + 1].value[i] - pop[j - 1].value[i]) / (max - min);
		}
	}
	cout << "every crowded distance done" << endl;
	for (int i = 0; i < pop.size(); i++)
	{
		double crowded = 0;
		for (int j = 0; j < objective; j++)
		{
			crowded += pop[i].crowded_distace[j];
		}
		pop[i].crowded_distace[objective] = crowded;
	}
	cout << "sum crowded distance done" << endl;
	/*
	individual a;
	for (int i = 0; i < pop.size(); i++)
	{
		for (int j = 0; j < pop.size() - 1; j++)
		{
			if (pop[j].crowded_distace[objective] < pop[j+1].crowded_distace[objective])
			{
				a = pop[j];
				pop[j] = pop[j+1];
				pop[j+1] = a;
			}
		}
	}
	*/
	sort(pop.begin(), pop.end(), paixusum);
	cout << "done";
}

void population::make_new_pop()
{
	cout << endl << "make new pop" << endl;
	int num_rank = 0;
	for (int i = 0; i < popsize * 2; i++)
	{
		if (R_nondominant_rank[i] > num_rank)
			num_rank = R_nondominant_rank[i];
	}
	int getishu[popsize * 2] = { 0 };									//ÿ���ȼ��ĸ�����
	for (int i = 0; i < popsize * 2; i++)
	{
		for (int j = 1; j <= num_rank; j++)
		{
			if (R_nondominant_rank[i] == j)
				getishu[j - 1]++;
		}
	}
	int pop_num = 0;													//�Ѿ���ѡ��Ϊ�����ĸ�����
	int i = 0;															//�ȼ���		
	while (pop_num + getishu[i] <= popsize)
	{
		for (int j = 0; j < popsize * 2; j++)
		{
			if (R_nondominant_rank[j] == i + 1)
			{
				pop[pop_num] = pop_R[j];
				pop_num++;
			}
		}
		i++;
	}
	if (pop_num < popsize)
	{
		for (int j = 0; j < popsize * 2; j++)
		{
			if (R_nondominant_rank[j] == i + 1)
				pop_temporary.push_back(pop_R[j]);
		}
		cal_crowded_distance(pop_temporary);
		int shengxiageshu = popsize - pop_num;
		for (int j = 0; j < shengxiageshu; j++)
		{
			pop[pop_num] = pop_temporary[j];
			pop_num++;
		}
	}
	cout << "done" << endl;
}

void population::select_individual(individual& a, individual& b)
{
	int x1 = rand() % popsize;
	int x2;
	do
	{
		x2 = rand() % popsize;
	} while (x1 == x2);
	a = pop[x1];
	b = pop[x2];
}

void population::make_new_pop_Q(decision_space space)
{
	cout << "make new Q" << endl;
	int count = 0;
	do
	{
		if (popsize % 2 == 1 && count + 1 == popsize)
		{
			int x1 = rand() % popsize;
			individual c = pop[x1];
			variation(c, space);
			pop_Q[count] = c;
			count++;
		}
		else
		{
			individual a, b;
			select_individual(a, b);
			individual_cross(a, b);
			variation(a, space);
			variation(b, space);
			pop_Q[count] = a;
			pop_Q[count + 1] = b;
			count = count + 2;
		}
	} while (count < popsize);
	cout << " Q done" << endl;
}

void NSGA2_update_ep(population p)
{
	for (int i = 0; i < popsize; i++)
	{
		update_ep(p.pop_Q[i]);
	}
}  
         

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//			1.�����ʼ����Ⱥ��
//			2.��������䵽�����Ȩ��������
//			3.�����Ӵ����Ը��ʵ¶���ѡ���Ǵ��ھ�T��Ⱥ������������Ⱥ�������丸����Ȼ�󽻲���칲����N���Ӵ�,���Ӵ�
//			  ���䵽���������Ⱥ��
//			4. ɾ������
//					4.1 �Ƚ�������Ⱥ��֧������
//					N+N=2N���Ӵ�����N�����壬���ȴӶ����ԵĽǶȿ�������С��Ⱥ���������Ӷൽ�ٽ��п��죬���һ��С��Ⱥ
//					�г���1����֧��ȼ�����ô�ʹ����һ����֧��ȼ���ʼɾ�������ø�Ȩ���������б�ѩ�����ɾ�����壬Ȼ
//					�������һ��ɾ��������ת��3.2
//					���ÿ��С��Ⱥȫ�Ƿ�֧��ĸ��壬��ô�ʹӶ�������б�ѩ�������ɾ����ֱ����Ⱥ��ģ�ﵽN
//			5.��һ���Ĵ��������ú���ܱ���������������n1,n2,���¸��ʵ¶���

void population::subpop_assign()
{
	for (int i = 0; i < popsize; i++)
	{
		double cos = -1;										//������ο������ļн�
		int zuijin;												//������ӽ��������ı��
		for (int j = 0; j < num_vector; j++)
		{
			double diancheng = 0;
			double getimo = 0;
			double vectormo = 0;
			for (int k = 0; k < objective; k++)
			{
				diancheng += ws[j][k] * pop[i].value[k];
				getimo += pop[i].value[k] * pop[i].value[k];
				vectormo += ws[j][k] * ws[j][k];
			}
			getimo = sqrt(getimo);
			vectormo = sqrt(vectormo);
			diancheng = diancheng / (getimo * vectormo);
			if (diancheng > cos)
			{
				cos = diancheng;
				zuijin = j;
			}
		}
		subpop[zuijin].push_back(pop[i]);
	}
}

void population::individual_assign_offspring_pop(individual x)
{
	double cos = -1;										//������ο������ļн�
	int zuijin;												//������ӽ��������ı��
	for (int j = 0; j < num_vector; j++)
	{
		double diancheng = 0;
		double getimo = 0;
		double vectormo = 0;
		for (int k = 0; k < objective; k++)
		{
			diancheng += ws[j][k] * x.value[k];
			getimo += x.value[k] * x.value[k];
			vectormo += ws[j][k] * ws[j][k];
		}
		getimo = sqrt(getimo);
		vectormo = sqrt(vectormo);
		diancheng = diancheng / (getimo * vectormo);
		if (diancheng > cos)
		{
			cos = diancheng;
			zuijin = j;
		}
	}
	offspring_pop[zuijin].push_back(x);
}

void population::subpop_nondominant_sort()
{
	for (int k = 0; k < num_vector; k++)
	{
		if (subpop[k].size() > 0)
		{
			int* np = new int[subpop[k].size()];						//֧��ø���ĸ���
			vector<int>* sp = new vector<int>[subpop[k].size()];			//������֧��ĸ��弯��
			for (int i = 0; i < subpop[k].size(); i++)
			{
				int count = 0;
				for (int j = 0; j < subpop[k].size(); j++)
				{
					if (subpop[k][j] < subpop[k][i])
					{
						count++;
					}
					if (subpop[k][i] < subpop[k][j])
					{
						sp[i].push_back(j);
					}
				}
				np[i] = count;
			}
			int rank = 1;
			int* mark = new int[subpop[k].size()];
			for (int i = 0; i < subpop[k].size(); i++)
			{
				mark[i] = 0;
			}
			for (int i = 0; i < subpop[k].size(); i++)
			{
				subpop[k][i].rank = 0;
			}
			int count = 0;
			do
			{
				for (int i = 0; i < subpop[k].size(); i++)
				{
					if (np[i] == 0 && mark[i] == 0)
					{
						subpop[k][i].rank = rank;
						mark[i] = 1;
						count++;
					}
				}
				for (int i = 0; i < subpop[k].size(); i++)
				{
					if (subpop[k][i].rank == rank)
					{
						for (int j = 0; j < sp[i].size(); j++)
						{
							np[sp[i][j]]--;
						}
					}
				}
				rank++;
			} while (count < subpop[k].size());
			delete[]np;
			delete[]sp;
		}
	}
}

int population::distance_subpop_index(int x)
{
	int a;											
	vector<int> index;												//��Ž�Զ������Ⱥ�еĽ�ĸ���>0������Ⱥ���
	for (int i = 0; i < num_vector; i++)
	{
		int flag = 0;
		for (int j = 0; j < num_T; j++)
		{
			if (i == neighbor_index[x][j])
			{
				flag = 1;
				break;
			}
		}
		if (flag)
			continue;
		else
		{
			if (subpop[i].size() > 0)
				index.push_back(i);
		}
	}
	if (index.size() == 0)
		return -2;														//�����Զ��Ⱥ��û�нⷵ��0
	if (index.size() == 1 && subpop[index[0]].size() == 1)				//�����Զ��ֻ��һ����Ⱥ����Ⱥ��ֻ��һ�����򷵻�0
		return -3;
	if (index.size() == 1 && subpop[index[0]].size() > 1)				//��Զֻ��һ����Ⱥ��������Ⱥ���ж����
		return -1;
	if (index.size() > 1)
	{
		a = rand() % index.size();
		return index[a];
	}
}
/*
void population::make_offspring(decision_space space)
{
	for (int i = 0; i < num_vector; i++)
	{
		for (int j = 0; j < subpop[i].size(); j++)
		{
			if (rand() % 100 / 100.0 < deerta)							//���ھ���Ⱥ��ѡ��������
			{
				vector<individual> tempoprary_pop;
				vector<individual> tem_sum_pop;
				for (int k = 0; k < num_T; k++)
				{
					for (int q = 0; q < subpop[neighbor_index[i][k]].size(); q++)
					{
						tem_sum_pop.push_back(subpop[neighbor_index[i][k]][q]);
						if (subpop[neighbor_index[i][k]][q].rank == 1)
							tempoprary_pop.push_back(subpop[neighbor_index[i][k]][q]);
					}
				}
				if (tempoprary_pop.size() >= 2)
				{														//��֧����������ڵ���2��
					int x1 = rand() % tempoprary_pop.size();
					int x2;
					do
					{
						x2 = rand() % tempoprary_pop.size();
						cout << "c";
					} while (x1 == x2);
					individual a = tempoprary_pop[x1];
					individual b = tempoprary_pop[x2];
					cout << "���뽻��ĸ�����" << endl;
					a.print();
					b.print();
					individual_cross(a, b);
					p_variation(a, space);
					p_variation(b, space);
					a.source = 0;
					b.source = 0;
					individual_assign_offspring_pop(a);
					individual_assign_offspring_pop(b);
					cout << "1";
				}
				else                             
				{														//��֧�������С��2���������ڵ���Ⱥ�ϲ������ѡ�񸸴����н���
					if (tem_sum_pop.size() >= 2)
					{													//���������Ⱥ�ĸ�������>=2���Ž��н��䣬���򲻽��н���												
						int x1 = rand() % tem_sum_pop.size();
						int x2;
						do
						{
							x2 = rand() % tem_sum_pop.size();
							cout << "a";
						} while (x1 == x2);
						individual a = tem_sum_pop[x1];
						individual b = tem_sum_pop[x2];
						cout << "���뽻��ĸ�����" << endl;
						a.print();
						b.print();
						individual_cross(a, b);
						p_variation(a, space);
						p_variation(b, space);
						a.source = 0;
						b.source = 0;
						individual_assign_offspring_pop(a);
						individual_assign_offspring_pop(b);
						cout << "2";
					}
				}
				tempoprary_pop.clear();
				tem_sum_pop.clear();
			}
			else
			{																//�ӽ�Զ������Ⱥ��ѡ�����
				if (distance_subpop_index(i) < 0)
				{
					if (distance_subpop_index(i) == -2)						//��Զ��Ⱥ��û�н����ֻ��1���⣬���������һ���½�
					{
						individual a;
						a.source = 2;
						a.random_select();
						a.cal_value(space);
						a.source = 1;
						individual_assign_offspring_pop(a);
						cout << "3";
					}
					if (distance_subpop_index(i) == -1)						//��Զ��ֻ��һ����Ⱥ�н����ж����
					{
						int x;
						for (int p = 0; p < num_vector; p++)
						{
							int flag = 0;
							for (int w = 0; w < num_T; w++)
							{
								if (p == neighbor_index[i][w])
								{
									flag = 1;
									break;
								}
							}
							if (flag)
								continue;
							else
							{
								if (subpop[p].size() > 0)
								{
									x = p;
									break;
								}
							}
						}
						int x1 = rand() % subpop[x].size();
						int x2;
						cout << "===" << subpop[x].size() << "===";
						do
						{
							x2 = rand() % subpop[x].size();
							cout << "q";
						} while (x1 == x2);
						individual a = subpop[x][x1];
						individual b = subpop[x][x2];
						cout << "���뽻��ĸ�����" << endl;
						a.print();
						b.print();
						individual_cross(a, b);
						p_variation(a, space);
						p_variation(b, space);
						a.source = 3;
						b.source = 3;
						individual_assign_offspring_pop(a);
						individual_assign_offspring_pop(b);
						cout << "4";
					}
				}
				else
				{
					int x1 = distance_subpop_index(i);
					int x2;
					do
					{
						x2 = distance_subpop_index(i);
						cout << "b";
					} while (x1 == x2);
					int x3 = rand() % subpop[x1].size();
					int x4 = rand() % subpop[x2].size();
					individual a = subpop[x1][x3];
					individual b = subpop[x2][x4];
					cout << "���뽻��ĸ�����" << endl;
					a.print();
					b.print();
					individual_cross(a, b);
					p_variation(a, space);
					p_variation(b, space);
					a.source = 1;
					b.source = 1;
					individual_assign_offspring_pop(a);
					individual_assign_offspring_pop(b);
					cout << "5";
				}
			}
		}
	}
}
*/

void population::update_ep_nsga_moead()
{
	for (int i = 0; i < num_vector; i++)
	{
		for (int j = 0; j < subpop[i].size(); j++)
			update_ep(subpop[i][j]);
	}
}

void population::make_offspring(decision_space space)
{
	for (int i = 0; i < num_vector; i++)
	{
		for (int j = 0; j < subpop[i].size(); j++)
		{
			if (rand() % 100 / 100.0 < deerta)							//���ھ���Ⱥ��ѡ��������
			{
				vector<individual> tempoprary_pop;
				vector<individual> tem_sum_pop;
				for (int k = 0; k < num_T; k++)
				{
					for (int q = 0; q < subpop[neighbor_index[i][k]].size(); q++)
					{
						tem_sum_pop.push_back(subpop[neighbor_index[i][k]][q]);
						if (subpop[neighbor_index[i][k]][q].rank == 1)
							tempoprary_pop.push_back(subpop[neighbor_index[i][k]][q]);
					}
				}
				if (tempoprary_pop.size() >0)
				{														//��֧���������
					int x1 = rand() % tempoprary_pop.size();
					individual a = subpop[i][j];
					individual b = tempoprary_pop[x1];
					individual_cross(a, b);
					variation(a, space);
					variation(b, space);
					a.source = 0;
					b.source = 0;
					individual_assign_offspring_pop(a);
					individual_assign_offspring_pop(b);
					cout << "1";
				}
				else
				{														//��֧�������С��2���������ڵ���Ⱥ�ϲ������ѡ�񸸴����н���
					if (tem_sum_pop.size() >= 2)
					{													//���������Ⱥ�ĸ�������>=2���Ž��н��䣬���򲻽��н���												
						int x1 = rand() % tem_sum_pop.size();
						individual a = subpop[i][j];
						individual b = tem_sum_pop[x1];
						individual_cross(a, b);
						variation(a, space);
						variation(b, space);
						a.source = 0;
						b.source = 0;
						individual_assign_offspring_pop(a);
						individual_assign_offspring_pop(b);
						cout << "2";
					}
				}
				tempoprary_pop.clear();
				tem_sum_pop.clear();
			}
			else
			{																//�ӽ�Զ������Ⱥ��ѡ�����
				if (distance_subpop_index(i) < 0)
				{
					if (distance_subpop_index(i) == -2)						//��Զ��Ⱥ��û�н⣬���������һ���½�
					{
						individual a;
						a.source = 2;
						a.random_select();
						a.cal_value(space);
						a.source = 1;
						individual_assign_offspring_pop(a);
						cout << "3";
					}
					if (distance_subpop_index(i) == -3)						//��Զ��Ⱥ������ֻ��1���⣬����������⽻��
					{
						int x;
						for (int p = 0; p < num_vector; p++)
						{
							int flag = 0;
							for (int w = 0; w < num_T; w++)
							{
								if (p == neighbor_index[i][w])
								{
									flag = 1;
									break;
								}
							}
							if (flag)
								continue;
							else
							{
								if (subpop[p].size() > 0)
								{
									x = p;
									break;
								}
							}
						}
						individual a = subpop[i][j];
						individual b = subpop[x][0];
						individual_cross(a, b);
						variation(a, space);
						variation(b, space);
						a.source = 1;
						b.source = 1;
						individual_assign_offspring_pop(a);
						individual_assign_offspring_pop(b);
						cout << "4";
					}
					if (distance_subpop_index(i) == -1)						//��Զ��ֻ��һ����Ⱥ�н����ж����
					{
						int x;
						for (int p = 0; p < num_vector; p++)
						{
							int flag = 0;
							for (int w = 0; w < num_T; w++)
							{
								if (p == neighbor_index[i][w])
								{
									flag = 1;
									break;
								}
							}
							if (flag)
								continue;
							else
							{
								if (subpop[p].size() > 0)
								{
									x = p;
									break;
								}
							}
						}
						int x1 = rand() % subpop[x].size();
						individual a = subpop[i][j];
						individual b = subpop[x][x1];
						individual_cross(a, b);
						variation(a, space);
						variation(b, space);
						a.source = 1;
						b.source = 1;
						individual_assign_offspring_pop(a);
						individual_assign_offspring_pop(b);
						cout << "4";
					}
				}
				else
				{
					int x1 = distance_subpop_index(i);
					
					int x3 = rand() % subpop[x1].size();
					individual a = subpop[x1][x3];
					individual b = subpop[i][j];
					individual_cross(a, b);
					variation(a, space);
					variation(b, space);
					a.source = 1;
					b.source = 1;
					individual_assign_offspring_pop(a);
					individual_assign_offspring_pop(b);
					cout << "5";
				}
			}
		}
	}
}

void population::merge()
{
	for (int i = 0; i < num_vector; i++)
	{
		for (int j = 0; j < offspring_pop[i].size(); j++)
		{
			
			int flag = 1;													//1-��ͬ�Ľ⣬���Լ�����Ⱥ  0-������ͬ�Ľ�
			for (int k = 0; k < subpop[i].size(); k++)
			{
				int count = 0;													//����ѡ�����ͬ�ĸ���
				for (int p = 0; p < dimension; p++)
				{
					if (subpop[i][k].select[p] == offspring_pop[i][j].select[p])
						count++;
				}
				if (count == dimension)
				{
					flag = 0;
					break;
				}
			}
			if(flag)
				subpop[i].push_back(offspring_pop[i][j]);
		}
	}
	for (int i = 0; i < num_vector; i++)
	{
		offspring_pop[i].clear();
	}
}

void population::eliminate()
{
	int subpop_num[num_vector];											//��¼ÿ������Ⱥ�Ĵ�С
	int subpop_index[num_vector];										//��¼����Ⱥ�ı��
	int sum_popsize = 0;												//��ǰsubpop���ܸ�����
	int rank[num_vector] = { 0 };										//ÿ����Ⱥ�����ȼ�
	for (int i = 0; i < num_vector; i++)
	{
		sum_popsize += subpop[i].size();
	}
	do
	{
		memset(rank, 0, sizeof(rank));
		for (int i = 0; i < num_vector; i++)							
		{
			subpop_num[i] = subpop[i].size();
			subpop_index[i] = i;
			for (int j = 0; j < subpop[i].size(); j++)					//ÿ����Ⱥ�����ȼ�
			{
				if (subpop[i][j].rank > rank[i])
					rank[i] = subpop[i][j].rank;
			}
		}
		for (int i = 0; i < num_vector; i++)
		{																	//������Ⱥ���մӴ�С����
			for (int j = 0; j < num_vector - 1; j++)
			{
				if (subpop_num[j] < subpop_num[j + 1])
				{
					int x;
					x = subpop_num[j];
					subpop_num[j] = subpop_num[j + 1];
					subpop_num[j + 1] = x;
					x = subpop_index[j];
					subpop_index[j] = subpop_index[j + 1];
					subpop_index[j + 1] = x;
					x = rank[j];
					rank[j] = rank[j + 1];
					rank[j + 1] = x;
				}
			}
		}
		int num_rank_1_pop=0;												//����Ⱥȫ�Ƿ�֧������Ⱥ��
		for (int i = 0; i < num_vector; i++)
		{
			if (rank[i]>1)
			{
				int count = 0;												//���һ���ĸ���
				for (int j = 0; j < subpop[subpop_index[i]].size(); j++)
				{
					if (subpop[subpop_index[i]][j].rank == rank[i])
						count++;
				}
				int* index = new int[count];
				double* value = new double[count];
				int ct = 0;
				for (int j = 0; j < subpop[subpop_index[i]].size(); j++)
				{
					if (subpop[subpop_index[i]][j].rank == rank[i])
					{
						value[ct] = cal_y_g(subpop[subpop_index[i]][j], ws[subpop_index[i]]);
						index[ct] = j;
						ct++;
					}
				}
				int delete_index;
				double va = 0;
				for (int j = 0; j < count; j++)
				{
					if (value[j] > va)
					{
						va = value[j];
						delete_index = index[j];
					}
				}
				subpop[subpop_index[i]].erase(subpop[subpop_index[i]].begin() + delete_index);
				sum_popsize--;
				break;
			}
			if (rank[i] == 1||rank[i]==0)
			{
				num_rank_1_pop++;
			}
		}
		if (num_rank_1_pop == num_vector)
		{
			int count = 0;												//���һ���ĸ���
			for (int j = 0; j < subpop[subpop_index[0]].size(); j++)
			{
				if (subpop[subpop_index[0]][j].rank == 1)
					count++;
			}
			int* index = new int[count];
			double* value = new double[count];
			int ct = 0;
			for (int j = 0; j < subpop[subpop_index[0]].size(); j++)
			{
				if (subpop[subpop_index[0]][j].rank == 1)
				{
					value[ct] = cal_y_g(subpop[subpop_index[0]][j], ws[subpop_index[0]]);
					index[ct] = j;
					ct++;
				}
			}
			int delete_index;
			double va = 0;
			for (int j = 0; j < count; j++)
			{
				if (value[j] > va)
				{
					va = value[j];
					delete_index = index[j];
				}
			}
			subpop[subpop_index[0]].erase(subpop[subpop_index[0]].begin() + delete_index);
			sum_popsize--;
		}
	} while (sum_popsize > popsize);
}

int n0=0;														//�½������Ը���������ĸ���
int n1=0;														//�½������Խ�Զ������ĸ���
void population::update_deerta(int dd,int L)
{
	if (dd % L == 0)
	{
		if (n0 + n1 == 0)
			deerta = 0;
		else
			deerta = n0/((n0+n1)*1.0);
		if (deerta < 0.2)
			deerta = 0.8;
		if (deerta > 0.8)
			deerta = 0.2;
		n0 = 0;
		n1 = 0;
	}
	else
	{
		for (int i = 0; i < num_vector; i++)
		{
			for (int j = 0; j < subpop[i].size(); j++)
			{
				if (subpop[i][j].source == 0)
					n0++;
				if (subpop[i][j].source == 1)
					n1++;
			}
		}
	}
	for (int i = 0; i < num_vector; i++)							//ÿһ������n0,n1�������source
	{
		for (int j = 0; j < subpop[i].size(); j++)
		{
			subpop[i][j].source = 4;
		}
	}
	deerta_data.push_back(deerta);
}

void save_deerta()
{
	string str1 = ".\\data\\num_update\\";
	string str2 = "deerta_data.xls";
	str2 = str1 + str2;
	ofstream out(str2);
	for (int i = 0; i <deerta_data.size(); i++)
	{
		
		out << deerta_data[i] << endl;
	}
	out.close();
}
