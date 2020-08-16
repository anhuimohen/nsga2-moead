#pragma once
#include<vector>
#include<iostream>
#include<fstream>
#include<iomanip>
#define dimension 10                                                //���߿ռ��ά��
#define concrete_service_num 10                                     //�������ĸ���
#define objective 3												   //Ŀ��������������
#define ob1_max 1000                                               //Ŀ������ֵ,cost
#define ob2_max 100                                                //time
#define ob3_max 10                                                 //other
using namespace std;
int ob_cal_mode[objective] = { 0,0,1};                            //Ŀ��ֵ���㷽ʽ 0-��� 1-���
class service                                                      //��������Ԫ
{
public:
	double attribute_value[objective][2];                          //����ֵ��0-ʵ��ֵ 1-��һ�����ֵ
	int attribute_flag[objective] = { 0,0,1};                     //���Լ��ԣ�0-�������� 1-��������
	void init();
	void print();

};

void service::init()                                               //�����ʼ��
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

void service::print()                    //�����ӡ
{
	cout << "�����ӡ��";
	for (int i = 0; i < objective; i++)
	{
		cout << setw(10) << attribute_value[i][0];
	}
	cout << endl;
	if (attribute_value[0][1] <= 1 && attribute_value[0][1] >= 0)
	{
		cout << "��һ����ķ�������ֵΪ��";
		for (int i = 0; i < objective; i++)
		{
			cout << setw(10) << attribute_value[i][1];
		}
		cout << endl;
	}
}

class decision_space                                               //���߿ռ�
{
public:
	service service_unit[dimension][concrete_service_num];
	void init();                                                   //���߿ռ�ĳ�ʼ��
	void print();
	void write_true_value();
	void write_normalization_value();
	void read_true_value();
	void read_normalization_value();
};

void decision_space::init()                                        //���߿ռ�ĳ�ʼ��
{
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < concrete_service_num; j++)
		{
			service_unit[i][j].init();
		}
	}
	///////////////////////////////////
	//  ���߿ռ��ʼ����ɺ�Ϳ��Խ����һ����0-1֮��
	for (int i = 0; i < dimension; i++)
	{
		double max[objective];                                     //ÿά�������Ե����ֵ
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
			for (int k = 0; k < objective; k++)                          //0-��������ԽСԽ��  1-��������Խ��Խ��
			{
				if (ob_cal_mode[k] == 0)                                 //Ŀ����㷽ʽ�� ���
				{
					if (service_unit[i][j].attribute_flag[k] == 0)       //����Ŀ�� ����Сֵ�Ż�
						service_unit[i][j].attribute_value[k][1] = (service_unit[i][j].attribute_value[k][0] - min[k]) / (max[k] - min[k]);
					if (service_unit[i][j].attribute_flag[k] == 1)
						service_unit[i][j].attribute_value[k][1] = (max[k] - service_unit[i][j].attribute_value[k][0]) / (max[k] - min[k]);
				}
				if (ob_cal_mode[k] == 1)                                 //Ŀ����㷽ʽ�� ���  ������һ������0�ķ���
				{
					if (service_unit[i][j].attribute_flag[k] == 0)       //����Ŀ�� ����Сֵ�Ż�
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
	cout << "���߿ռ��ӡ��" << endl;
	for (int i = 0; i < dimension; i++)
	{
		cout << "��" << i << "ά�ռ�--------------------------" << endl;
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

