#pragma once
#include<fstream>
#include"NSGA2.h"

vector<individual> moead_solution;
vector<individual> aesmoead_solution;
vector<individual> non_dominated_solution;

void read_solution(string str, int hangshu, vector<individual>& solution)	//��ȡstr�����ݵ�solution�����ݵ�������hangshu
{																		//��Щ�ⶼ��EP��ֻȡ��Ŀ��ֵ�Ĵ�С������IGD	
	ifstream infile(str);
	if (!infile)
		cout << "�ļ�û�д�" << endl;
	for (int i = 0; i < hangshu; i++)
	{
		individual a;
		for (int j = 0; j < objective; j++)
		{
			infile >> a.value[j];
		}
		solution.push_back(a);
	}
	infile.close();
}

void add_solution(vector<individual> b, vector<individual>& c)
{													//��a��b�Ľ��һ�������еķ�֧���ŵ�c��

	for (int i = 0; i < b.size(); i++)
	{
		for (int j = 0; j < c.size(); j++)			//��ɾ��c�б�b֧��Ľ�
		{
			if (b[i] < c[j])
			{
				//c[j].print();
				c.erase(c.begin() + j);
				j--;
			}
		}
		//Ȼ���ٿ����Ƿ�b���뵽c��
		int flag = 1;								//1-c��û��֧��b�н�Ľ⣬���Խ�����뵽c��
		for (int j = 0; j < c.size(); j++)
		{
			if (c[j] < b[i])
			{
				flag = 0;
				break;
			}
			int count = 0;							//������Ŀ��ֵ��ȵĸ���
			for (int k = 0; k < objective; k++)
			{
				if (c[j].value[k] == b[i].value[k])
					count++;
			}
			if (count == objective)
				flag = 0;
		}
		if (flag)
			c.push_back(b[i]);
		//else
		//	b[i].print();
	}
}

double cal_IGD(vector<individual> a, vector<individual> b)
{												//�⼯a��b��������ʵPF)�ϵ�IGD
	double igd = 0;
	for (int i = 0; i < b.size(); i++)
	{
		double d_min = RAND_MAX;
		for (int j = 0; j < a.size(); j++)
		{
			double d1[objective] = { 0 };
			for (int k = 0; k < objective; k++)
			{
				d1[k] = pow((b[i].value[k] - a[j].value[k]), 2);
				d1[0] += d1[k];
			}
			if (sqrt(d1[0]) < d_min)
				d_min = d1[0];
		}
		igd += d_min;
	}
	return igd / b.size();
}

double cal_coverage(vector<individual> a, vector<individual> b)
{
	int count = 0;
	for (int i = 0; i < b.size(); i++)
	{
		for (int j = 0; j < a.size(); j++)
		{
			if (a[j] < b[i])
			{
				count++;
				break;
			}
		}
	}
	//cout << count;
	return count / (b.size() * 1.0);
}
