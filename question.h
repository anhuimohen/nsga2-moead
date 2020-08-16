#pragma once
#include<vector>
#include<iostream>
#include<fstream>
#include<iomanip>
#define dimension 10                                                //决策空间的维度
#define concrete_service_num 10                                     //具体服务的个数
#define objective 3												   //目标数或者属性数
#define ob1_max 1000                                               //目标的最大值,cost
#define ob2_max 100                                                //time
#define ob3_max 10                                                 //other
using namespace std;
int ob_cal_mode[objective] = { 0,0,1};                            //目标值计算方式 0-相加 1-相乘
class service                                                      //基本服务单元
{
public:
	double attribute_value[objective][2];                          //属性值，0-实际值 1-归一化后的值
	int attribute_flag[objective] = { 0,0,1};                     //属性极性，0-消极属性 1-积极属性
	void init();
	void print();

};

void service::init()                                               //服务初始化
{
	for (int i = 0; i < objective; i++)
	{
		if (i == 0)
			attribute_value[i][0] = rand() % ob1_max;
		if (i == 1)
			attribute_value[i][0] = rand() % ob2_max;
		if (i >= 2)
			attribute_value[i][0] = rand() % ob3_max;
	}
}

void service::print()                    //服务打印
{
	cout << "服务打印：";
	for (int i = 0; i < objective; i++)
	{
		cout << setw(10) << attribute_value[i][0];
	}
	cout << endl;
	if (attribute_value[0][1] <= 1 && attribute_value[0][1] >= 0)
	{
		cout << "归一化后的服务属性值为：";
		for (int i = 0; i < objective; i++)
		{
			cout << setw(10) << attribute_value[i][1];
		}
		cout << endl;
	}
}

class decision_space                                               //决策空间
{
public:
	service service_unit[dimension][concrete_service_num];
	void init();                                                   //决策空间的初始化
	void print();
	void write_true_value();
	void write_normalization_value();
	void read_true_value();
	void read_normalization_value();
};

void decision_space::init()                                        //决策空间的初始化
{
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < concrete_service_num; j++)
		{
			service_unit[i][j].init();
		}
	}
	///////////////////////////////////
	//  决策空间初始化完成后就可以将其归一化到0-1之间
	for (int i = 0; i < dimension; i++)
	{
		double max[objective];                                     //每维服务属性的最大值
		double min[objective];
		memset(max, 0, sizeof(max));
		for (int j = 0; j < objective; j++)
		{
			min[j] = RAND_MAX;
		}
		for (int j = 0; j < concrete_service_num; j++)
		{
			for (int k = 0; k < objective; k++)
			{
				if (service_unit[i][j].attribute_value[k][0] > max[k])
					max[k] = service_unit[i][j].attribute_value[k][0];
				if (service_unit[i][j].attribute_value[k][0] < min[k])
					min[k] = service_unit[i][j].attribute_value[k][0];
			}
		}
		for (int j = 0; j < concrete_service_num; j++)
		{
			for (int k = 0; k < objective; k++)                          //0-消极属性越小越好  1-积极属性越大越好
			{
				if (ob_cal_mode[k] == 0)                                 //目标计算方式是 相加
				{
					if (service_unit[i][j].attribute_flag[k] == 0)       //最终目标 是最小值优化
						service_unit[i][j].attribute_value[k][1] = (service_unit[i][j].attribute_value[k][0] - min[k]) / (max[k] - min[k]);
					if (service_unit[i][j].attribute_flag[k] == 1)
						service_unit[i][j].attribute_value[k][1] = (max[k] - service_unit[i][j].attribute_value[k][0]) / (max[k] - min[k]);
				}
				if (ob_cal_mode[k] == 1)                                 //目标计算方式是 相乘  消除归一化后是0的风险
				{
					if (service_unit[i][j].attribute_flag[k] == 0)       //最终目标 是最小值优化
						service_unit[i][j].attribute_value[k][1] = (service_unit[i][j].attribute_value[k][0] - min[k] + 1) / (max[k] - min[k] + 1);
					if (service_unit[i][j].attribute_flag[k] == 1)
						service_unit[i][j].attribute_value[k][1] = (max[k] - service_unit[i][j].attribute_value[k][0] + 1) / (max[k] - min[k] + 1);
				}
			}
		}
	}
}

void decision_space::print()
{
	cout << "决策空间打印：" << endl;
	for (int i = 0; i < dimension; i++)
	{
		cout << "第" << i << "维空间--------------------------" << endl;
		for (int j = 0; j < concrete_service_num; j++)
		{
			service_unit[i][j].print();
		}
	}
}

void decision_space::write_true_value()
{
	ofstream outfile("decision_space_true_value.txt");
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < concrete_service_num; j++)
		{
			for (int k = 0; k < objective; k++)
			{
				outfile << service_unit[i][j].attribute_value[k][0] << " ";
			}
		}
		outfile << endl;
	}
	outfile.close();
}

void decision_space::write_normalization_value()
{
	ofstream outfile("decision_space_normalization_value.txt");
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < concrete_service_num; j++)
		{
			for (int k = 0; k < objective; k++)
			{
				outfile << service_unit[i][j].attribute_value[k][1] << " ";
			}
		}
		outfile << endl;
	}
	outfile.close();
}

void decision_space::read_true_value()
{
	char a;
	ifstream infile("decision_space_true_value.txt");
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < concrete_service_num; j++)
		{
			for (int k = 0; k < objective; k++)
			{
				infile >> service_unit[i][j].attribute_value[k][0];
			}
		}
	}
	infile.close();
}

void decision_space::read_normalization_value()
{
	ifstream infile("decision_space_normalization_value.txt");
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < concrete_service_num; j++)
		{
			for (int k = 0; k < objective; k++)
			{
				infile >> service_unit[i][j].attribute_value[k][1];
			}
		}
	}
	infile.close();
}

