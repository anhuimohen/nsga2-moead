#pragma once
#include<fstream>
#include"NSGA2.h"

vector<individual> moead_solution;
vector<individual> aesmoead_solution;
vector<individual> non_dominated_solution;

void read_solution(string str, int hangshu, vector<individual>& solution)	//读取str的数据到solution，数据的行数是hangshu
{																		//这些解都是EP，只取其目标值的大小，计算IGD	
	ifstream infile(str);
	if (!infile)
		cout << "文件没有打开" << endl;
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
{													//把a和b的解放一起，求其中的非支配解放到c中

	for (int i = 0; i < b.size(); i++)
	{
		for (int j = 0; j < c.size(); j++)			//先删除c中被b支配的解
		{
			if (b[i] < c[j])
			{
				//c[j].print();
				c.erase(c.begin() + j);
				j--;
			}
		}
		//然后再考虑是否将b加入到c中
		int flag = 1;								//1-c中没有支配b中解的解，可以将其加入到c中
		for (int j = 0; j < c.size(); j++)
		{
			if (c[j] < b[i])
			{
				flag = 0;
				break;
			}
			int count = 0;							//两个解目标值相等的个数
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
{												//解集a在b（近似真实PF)上的IGD
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
