#if 1
#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<string>
using namespace std;
#define MAX(x,y) ((x)>(y)?(x):(y))
template<typename T>
void split_generic(vector<T> &v, const T & str, const T & delimiters) {
	//vector<T> v;
	v.clear();	// clear v to be empty
	typename T::size_type start = 0;
	auto pos = str.find_first_of(delimiters, start);
	while (pos != T::npos) {
		if (pos != start) // ignore empty tokens
			v.emplace_back(str, start, pos - start);
		start = pos + 1;
		pos = str.find_first_of(delimiters, start);
	}
	if (start < str.length()) // ignore trailing delimiter
		v.emplace_back(str, start, str.length() - start); // add what's left of the string
														  //return v;
}
class Instance
{
public:
	class Process {
	public:
		int mach;
		int t;
		Process() {}
		Process(int _t, int _mach) :mach(_mach), t(_t) {}
	};
	class Operation
	{
	public:
		int num_mach;
		vector<Process*> proc_vec;
		Operation() {}
		Operation(int _num_mach, const vector<Process*> &_proc_vec) :num_mach(_num_mach), proc_vec(_proc_vec) {}
	};
	class JobInfo
	{
	public:
		int num_oper;
		vector<Operation*> oper_vec;
		JobInfo() {}
		JobInfo(int _num_oper, const vector<Operation*>&_oper_vec) :num_oper(_num_oper), oper_vec(_oper_vec) {}
	};

	int m, n;	// m and n are the number of machines and jobs
	int avg_num_mach_per_oper;	// average number of machines per operation.
	string file_input;
	vector<JobInfo*> job_vec;
	Instance(string);
	void display_instance();
};
class Solution
{
public:
	class Operation
	{
	public:
		int mach, start_time, end_time;
		Operation *pre_job_oper, *next_job_oper;
		Operation *pre_mach_oper, *next_mach_oper;
		const int job_i, oper_i;	// the job and operation index of this operation
		Operation(int _j, int _o, int _et) :end_time(_et), job_i(_j), oper_i(_o) {}
	};
	class JobInfo
	{
	public:
		vector<Operation*> oper_vec;
	};
	vector<JobInfo*> job_vec;
	vector<Operation*>oper_mach_vec;	// the operations assigned at each machine
	int makespan, makespan_mach;
	Solution(int _m, int _mm) :makespan(_m), makespan_mach(_mm) {}
};
class Solver
{
public:
	int sol_num;
	vector<Solution*> sol_vec;
	const Instance &instance;
	Solver(const Instance &,int);
	void init_solution(int);
	void display_solution(int);
};
Solver::Solver(const Instance &_instance,int _sol_num) :instance(_instance) ,sol_num(_sol_num)
{
	for (int i = 0; i <= sol_num; i++)
	{
		Solution *sol = new Solution(0, 0);
		sol->job_vec.push_back(NULL);	// the first element in job_vec is NULL
		for (int j = 1; j < instance.job_vec.size(); j++)
		{
			Solution::JobInfo *jobinfo = new Solution::JobInfo();
			for (int k = 0; k < instance.job_vec[j]->oper_vec.size(); k++)
			{
				Solution::Operation *oper = new Solution::Operation(j, k, 0);
				jobinfo->oper_vec.push_back(oper);
			}
			sol->job_vec.push_back(jobinfo);
		}
		sol_vec.push_back(sol);
	}
}
void Solver::display_solution(int sol_index)
{
	cout << "*** display solution " << sol_index << " ***, makespan: " << sol_vec[sol_index]->makespan
		<< ", makespan machine: " << sol_vec[sol_index]->makespan_mach << endl;
	for (vector<Solution::JobInfo*>::iterator job_iter = sol_vec[sol_index]->job_vec.begin() + 1;
		job_iter != sol_vec[sol_index]->job_vec.end(); job_iter++)
	{
		cout << (*job_iter)->oper_vec.size() - 1 << " | ";
		for (vector<Solution::Operation*>::iterator oper_iter = (*job_iter)->oper_vec.begin() + 1;
			oper_iter != (*job_iter)->oper_vec.end(); oper_iter++)
		{
			cout << (*oper_iter)->mach << "\t" << (*oper_iter)->start_time << "\t" << (*oper_iter)->end_time << "\t";
		}
		cout << endl;
	}	
}
void Solver::init_solution(int sol_index)
{
	vector<int> mach_ct_vec(instance.job_vec.size(), 0);	// the completion time of each machine
	vector<Solution::Operation*> last_oper_mach(instance.m + 1, NULL);	// the last operation at each machine
	vector<int> oper_to_assign_job(instance.job_vec.size(),1);	// the operation index needed to assign at each job
	bool is_not_finished = true;
	while (is_not_finished)
	{
		is_not_finished = false;
		for (int job_i = 1; job_i < instance.job_vec.size() && oper_to_assign_job[job_i]; job_i++)	// All the operations of this job is completed
		{
			is_not_finished = true;
			Instance::Operation *oper = instance.job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]];
			int min_ct = INT_MAX, min_ct_mach, equ_cnt = 1;
			for (int proc_i = 1; proc_i < oper->proc_vec.size(); proc_i++)	// find the best machine which has the minimum completion time
			{
				int cur_mach_ct = oper->proc_vec[proc_i]->t + MAX(mach_ct_vec[oper->proc_vec[proc_i]->mach],	// the next operation with the same job
					sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i] - 1]->end_time);		// the next operation with the same machine
				if (cur_mach_ct < min_ct)
				{
					min_ct = cur_mach_ct;
					min_ct_mach = oper->proc_vec[proc_i]->mach;
				}
				else if (cur_mach_ct== min_ct)
				{
					equ_cnt += 1;
					if (rand() % equ_cnt)
						min_ct_mach = oper->proc_vec[proc_i]->mach;
				}
			}
			sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]]->mach = min_ct_mach;
			sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]]->start_time = mach_ct_vec[min_ct_mach];
			sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]]->end_time = min_ct;
			sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]]->pre_mach_oper = last_oper_mach[min_ct_mach];
			sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]]->next_mach_oper = NULL;	// the next operation with the same machine is NULL
			mach_ct_vec[min_ct_mach] = min_ct;
			if (last_oper_mach[min_ct_mach] != NULL)	// this is not the first operation of the machine, assgin the next operation of the previous operation
				sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]]->pre_mach_oper->next_mach_oper =
				sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]];
			last_oper_mach[min_ct_mach] = sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]];	// update the last operation at this machine
			oper_to_assign_job[job_i] += 1;	// the next operation is to be assigned
			if (oper_to_assign_job[job_i] >= instance.job_vec[job_i]->oper_vec.size())
				oper_to_assign_job[job_i] = 0;	// All the operations of this job is completed
		}
	}
	for (int i = 1; i <= instance.m; i++)
	{
		if (sol_vec[sol_index]->makespan < mach_ct_vec[i])
		{
			sol_vec[sol_index]->makespan = mach_ct_vec[i];
			sol_vec[sol_index]->makespan_mach = i;
		}
	}
}
Instance::Instance(string _file_input):file_input(_file_input)
{
	ifstream ifs(file_input);
	if (!ifs.is_open())
	{
		cout << file_input << endl;
		perror("file_input"); exit(0);
	}
	string str_line;
	ifs >> n >> m >> avg_num_mach_per_oper;
	cout << n << "\t" << m << "\t" << avg_num_mach_per_oper << endl;
	vector<string>fields_vec;
	getline(ifs, str_line);	// read the '\n'
	job_vec.push_back(NULL);	// the first element in job_vec is NULL
	while (getline(ifs, str_line))
	{
		cout << str_line << endl;
		split_generic<string>(fields_vec, str_line, " \t");
		JobInfo *jobinfo = new JobInfo();
		jobinfo->oper_vec.push_back(NULL);	// the first element in oper_vec is NULL
		int field_cnt = 0;
		jobinfo->num_oper = stoi(fields_vec[field_cnt++]);	// number of operation
		for (int i = 0; i < jobinfo->num_oper; i++)	// for each operation
		{
			Operation *oper = new Operation();
			oper->proc_vec.push_back(NULL);	// the first element in proc_vec is NULL
			oper->num_mach = stoi(fields_vec[field_cnt++]);	// number of process
			for (int j = 0; j < oper->num_mach; j++)
			{
				Process *proc = new Process(stoi(fields_vec[field_cnt++]), stoi(fields_vec[field_cnt++]));
				oper->proc_vec.push_back(proc);
			}
			jobinfo->oper_vec.push_back(oper);
		}
		job_vec.push_back(jobinfo);
	}
	ifs.close();
}
void Instance::display_instance()
{
	cout << "*** display instance: " << file_input << " ***" << endl;
	cout << n << "\t" << m << "\t" << avg_num_mach_per_oper << endl;
	for (vector<JobInfo*>::iterator job_iter = job_vec.begin()+1;
	job_iter != job_vec.end(); job_iter++)
	{
		cout << (*job_iter)->num_oper <<" | ";
		for (vector<Operation*>::iterator oper_iter = (*job_iter)->oper_vec.begin()+1;
		oper_iter != (*job_iter)->oper_vec.end(); oper_iter++)
		{
			cout << (*oper_iter)->num_mach << " : ";
			for (vector<Process*>::iterator proc_iter = (*oper_iter)->proc_vec.begin()+1;
			proc_iter != (*oper_iter)->proc_vec.end(); proc_iter++)
			{
				cout << (*proc_iter)->mach << "\t" << (*proc_iter)->t << "\t";
			}
		}
		cout << endl;
	}
}
int main(int argc, char **argv)
{
	char *argv_win[] = { "",	// 0
		"_ifp", "instances\\Dauzere_Data\\",
		"_ifn", "01a",
		"_sol_num", "5"
	};
	cout << "This is the flexible job shop scheduling problem" << endl;
#ifdef _WIN32
	argc = sizeof(argv_win) / sizeof(argv_win[0]); argv = argv_win;
#endif
	map<string, string> argv_map;
	argv_map["_exe_name"] = argv[0];	// add the exe file name to argv map to append to the output file name
	for (int i = 1; i < argc; i += 2)
		argv_map[string(argv[i])] = string(argv[i + 1]);
	Instance *ins = new Instance(argv_map.at("_ifp") + argv_map.at("_ifn") + ".fjs");
	ins->display_instance();
	Solver *solver = new Solver(*ins, stoi(argv_map.at("_sol_num")));
	solver->init_solution(1);
	solver->display_solution(1);
#ifdef _WIN32
	system("pause");
#endif
	return 0;
}
#endif