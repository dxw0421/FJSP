#if 1
#include<iostream>
#include<fstream>
#include<vector>
#include<stack>
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
	int avg_num_mach_per_oper;	// average number of machines per operation
	int total_num_operation;//
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
		int proc_i, mach, t, start_time, end_time, r, q;	// optimize by removing proc_i or mach and t
		Operation *pre_mach_oper, *next_mach_oper;
		const int job_i, oper_i;	// the job and operation index of this operation
		int in_degree;
		Operation(int _j, int _o, int _et) :end_time(_et), job_i(_j), oper_i(_o), pre_mach_oper(NULL), next_mach_oper(NULL), in_degree(2) {}
	};
	class JobInfo
	{
	public:
		vector<Operation*> oper_vec;
	};
	vector<JobInfo*> job_vec;
	vector<Operation*>head_oper_mach_vec, tail_oper_mach_vec;	// the operations assigned at each machine
	int makespan, makespan_mach;
	Operation *start_dummy_oper, *end_dummy_oper;
	Solution(int, int);
};
Solution::Solution(int _m, int _mm) :makespan(_m), makespan_mach(_mm), start_dummy_oper(NULL), end_dummy_oper(NULL)
{
	start_dummy_oper = new Solution::Operation(0, 0, 0);
	end_dummy_oper = new Solution::Operation(0, 1000, 0);
	start_dummy_oper->r = end_dummy_oper->q = 0;
	start_dummy_oper->in_degree = end_dummy_oper->in_degree = 0;
}
class Solver
{
public:
	int sol_num;
	vector<Solution*> sol_vec;
	const Instance *instance;
	Solver(const Instance &,int);
	void init_solution(int);
	void determine_critical_path(int);
	void display_solution(int);
	void check_solution(int)const;
};
Solver::Solver(const Instance &_instance,int _sol_num) :instance(&_instance) ,sol_num(_sol_num)
{
	for (int i = 0; i <= sol_num; i++)
	{
		Solution *sol = new Solution(0, 0);
		sol->head_oper_mach_vec.resize(instance->m + 1, NULL);	// the initial value is NULL
		sol->tail_oper_mach_vec.resize(instance->m + 1, sol->start_dummy_oper);	// the initial value is start dummy operation
		sol->job_vec.push_back(NULL);	// the first element in job_vec is NULL
		for (int j = 1; j < instance->job_vec.size(); j++)
		{
			Solution::JobInfo *jobinfo = new Solution::JobInfo();
			jobinfo->oper_vec.push_back(sol->start_dummy_oper);	// the first element in oper_vec of this job is dummy operation
			for (int k = 1; k < instance->job_vec[j]->oper_vec.size(); k++)
			{
				Solution::Operation *oper = new Solution::Operation(j, k, 0);
				jobinfo->oper_vec.push_back(oper);
			}
			jobinfo->oper_vec[1]->in_degree -= 1;	// the degree of the first operation of the job is 1
			jobinfo->oper_vec.push_back(sol->end_dummy_oper);	// the last element in oper_vec of this job is dummy operation
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
	// for each machine to display solution
	for (int i = 1; i < sol_vec[sol_index]->head_oper_mach_vec.size(); i++)
	{
		cout << "m " << i << ": ";
		Solution::Operation *oper_same_mach = sol_vec[sol_index]->head_oper_mach_vec[i];
		while (oper_same_mach->job_i != 0)
		{
			cout << oper_same_mach->job_i << ", " << oper_same_mach->oper_i << "\t"
				<< oper_same_mach->start_time << "\t" << oper_same_mach->end_time << "\t";
			oper_same_mach = oper_same_mach->next_mach_oper;
		}
		cout << endl;
	}
	// for each machine to display solution in reverse order
	for (int i = 1; i < sol_vec[sol_index]->tail_oper_mach_vec.size(); i++)
	{
		cout << "m(R) " << i << ": ";
		Solution::Operation *oper_same_mach = sol_vec[sol_index]->tail_oper_mach_vec[i];
		while (oper_same_mach->job_i != 0)	// the last element in the list
		{
			cout << oper_same_mach->job_i << ", " << oper_same_mach->oper_i << "\t"
				<< oper_same_mach->start_time << "\t" << oper_same_mach->end_time << "\t";
			oper_same_mach = oper_same_mach->pre_mach_oper;
		}
		cout << endl;
	}
}
void Solver::check_solution(int sol_index)const
{
	cout << "***check solution " << sol_index << " ***" << endl;
	for (int job_i = 1; job_i < sol_vec[sol_index]->job_vec.size(); job_i++)
	{
		if (sol_vec[sol_index]->job_vec[job_i]->oper_vec[0]->end_time != 0)
		{
			cout << "ERROR: the end time of dummy operation is not 0" << endl;
			system("pause");
		}
		for (int oper_i = 1; oper_i <= instance->job_vec[job_i]->num_oper; oper_i++)
		{
			Solution::Operation *oper = sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_i];
			if (oper->start_time != MAX(sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_i - 1]->end_time, oper->pre_mach_oper->end_time))
			{
				cout << "ERROR: the start time is not max of the end time of previous operations of same job and machine" << endl;
				system("pause");
			}
			if (oper->proc_i == 0)
			{
				cout << "ERROR: the assigned machine is wrong" << endl;
				system("pause");
			}
			if (oper->end_time != oper->start_time + instance->job_vec[job_i]->oper_vec[oper_i]->proc_vec[oper->proc_i]->t)
			{
				cout << "ERROR: the end time is not the start time plus its processing time" << endl;
				system("pause");
			}
		}
	}
	// for each machine to check the solution
	int max_ct = 0, max_ct_mach;
	for (int mach_i = 1; mach_i < sol_vec[sol_index]->head_oper_mach_vec.size(); mach_i++)
	{
		Solution::Operation *oper = sol_vec[sol_index]->head_oper_mach_vec[mach_i];
		while (oper->job_i != 0)
		{
			if (oper->start_time != MAX(sol_vec[sol_index]->job_vec[oper->job_i]->oper_vec[oper->oper_i - 1]->end_time, oper->pre_mach_oper->end_time))
			{
				cout << "ERROR: the start time is not max of the end time of previous operations of same job and machine" << endl;
				system("pause");
			}
			int proc_i_select = 0;
			for (int proc_i = 1; proc_i < instance->job_vec[oper->job_i]->oper_vec[oper->oper_i]->proc_vec.size(); proc_i++)
			{
				if (oper->mach == instance->job_vec[oper->job_i]->oper_vec[oper->oper_i]->proc_vec[proc_i]->mach)
				{
					proc_i_select = proc_i;
					break;
				}
			}
			if (proc_i_select == 0)
			{
				cout << "ERROR: the assigned machine is wrong" << endl;
				system("pause");
			}
			if (oper->end_time != oper->start_time + instance->job_vec[oper->job_i]->oper_vec[oper->oper_i]->proc_vec[proc_i_select]->t)
			{
				cout << "ERROR: the end time is not the start time plus its processing time" << endl;
				system("pause");
			}
			if (max_ct < oper->end_time)
			{
				max_ct = oper->end_time;
				max_ct_mach = oper->mach;
			}
			oper = oper->next_mach_oper;
		}
	}
	if (sol_vec[sol_index]->makespan != max_ct || sol_vec[sol_index]->makespan_mach != max_ct_mach)
	{
		cout << "ERROR: makespan or makespan machine is wrong" << endl;
		system("pause");
	}
}
void Solver::init_solution(int sol_index)
{
	vector<int> mach_ct_vec(instance->job_vec.size(), 0);	// the completion time of each machine
	//vector<Solution::Operation*> last_oper_mach(instance->m + 1,sol_vec[sol_index]->start_dummy_oper);	// the last operation at each machine
	vector<int> oper_to_assign_job(instance->job_vec.size(),1);	// the operation index needed to assign at each job
	bool is_not_finished = true;
	while (is_not_finished)
	{
		is_not_finished = false;
		for (int job_i = 1; job_i < instance->job_vec.size(); job_i++)	// All the operations of this job is completed
		{
			if (oper_to_assign_job[job_i] == 0)	// All the operations of this job is completed
				continue;
			is_not_finished = true;
			Instance::Operation *oper = instance->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]];
			int min_ct = INT_MAX, min_ct_mach, select_proc_i, equ_cnt = 1;
			for (int proc_i = 1; proc_i < oper->proc_vec.size(); proc_i++)	// find the best machine which has the minimum completion time
			{
				int cur_mach_ct = oper->proc_vec[proc_i]->t + MAX(mach_ct_vec[oper->proc_vec[proc_i]->mach],	// the next operation with the same machine
					sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i] - 1]->end_time);		// the next operation with the same job
				if (cur_mach_ct < min_ct)
				{
					min_ct = cur_mach_ct;
					select_proc_i = proc_i;
					min_ct_mach = oper->proc_vec[proc_i]->mach;
				}
				else if (cur_mach_ct== min_ct)
				{
					equ_cnt += 1;
					if (rand() % equ_cnt)
					{
						select_proc_i = proc_i;
						min_ct_mach = oper->proc_vec[proc_i]->mach;
					}
				}
			}
			sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]]->proc_i = select_proc_i;
			sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]]->mach = min_ct_mach;
			sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]]->t = oper->proc_vec[select_proc_i]->t;
			sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]]->start_time = MAX(mach_ct_vec[min_ct_mach],
				sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i] - 1]->end_time); 
			sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]]->end_time = min_ct;
			sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]]->pre_mach_oper = sol_vec[sol_index]->tail_oper_mach_vec[min_ct_mach];
			sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]]->next_mach_oper = sol_vec[sol_index]->end_dummy_oper;	// the next operation with the same machine is dummy end operation
			mach_ct_vec[min_ct_mach] = min_ct;
			//sol_vec[sol_index]->tail_oper_mach_vec[min_ct_mach] = sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]];	// the tail operation assigned to the machine
			if (sol_vec[sol_index]->head_oper_mach_vec[min_ct_mach] == NULL)	// assign the first operation to the machine
			{
				sol_vec[sol_index]->head_oper_mach_vec[min_ct_mach] = sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]];
				sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]]->in_degree -= 1;	// the degree of the first opertion at the machine is 1
			}
			if (sol_vec[sol_index]->tail_oper_mach_vec[min_ct_mach]->job_i != 0)	// this is not the first operation of the machine, assgin the next operation of the previous operation
				sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]]->pre_mach_oper->next_mach_oper =
				sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]];
			sol_vec[sol_index]->tail_oper_mach_vec[min_ct_mach]= sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]];	// update the tail operation at this machine
			oper_to_assign_job[job_i] += 1;	// the next operation is to be assigned
			if (oper_to_assign_job[job_i] >= instance->job_vec[job_i]->oper_vec.size())
				oper_to_assign_job[job_i] = 0;	// All the operations of this job is completed
		}
	}
	for (int i = 1; i <= instance->m; i++)
	{
		if (sol_vec[sol_index]->makespan < mach_ct_vec[i])
		{
			sol_vec[sol_index]->makespan = mach_ct_vec[i];
			sol_vec[sol_index]->makespan_mach = i;
		}
	}
}
void Solver::determine_critical_path(int sol_index)
{
	cout << "***solution " << sol_index << ", critical path***" << endl;
	vector<Solution::Operation *> oper_in_degree_vec;	// store operations with in-degree of 0
	stack<Solution::Operation*> oper_re_top_order_stack;	// store operations in stack to get re-topological order
	stack<Solution::Operation*> critical_oper_stack;	// store critical operations
	for (int mach_i = 1; mach_i < sol_vec[sol_index]->head_oper_mach_vec.size(); mach_i++)
	{
		Solution::Operation *oper = sol_vec[sol_index]->head_oper_mach_vec[mach_i];
		if (oper->in_degree==0)	// or oper->oper_i == 1, which means the first operation at this job, in-degree is 0 
			oper_in_degree_vec.push_back(oper);
	}
	int num_oper_visited = 0;
	while (!oper_in_degree_vec.empty())
	{
		Solution::Operation *oper = oper_in_degree_vec.front();
		cout << oper->job_i << ", " << oper->oper_i << "\t" << oper->mach << "\t"
			<< oper->start_time << "\t" << oper->end_time << "\t";
		oper->r = MAX(oper->pre_mach_oper->r + oper->pre_mach_oper->t,
			sol_vec[sol_index]->job_vec[oper->job_i]->oper_vec[oper->oper_i-1]->r + sol_vec[sol_index]->job_vec[oper->job_i]->oper_vec[oper->oper_i-1]->t);
		if (oper->r != oper->start_time)
			cout << endl;
		num_oper_visited += 1;
		oper->next_mach_oper->in_degree -= 1;	// minus the in-degree of next operation at the same machine by 1
		sol_vec[sol_index]->job_vec[oper->job_i]->oper_vec[oper->oper_i + 1]->in_degree -= 1;	// minus the in-degree of next operation of the same job by 1
		if (oper->next_mach_oper->in_degree == 0)
			oper_in_degree_vec.push_back(oper->next_mach_oper);
		if (sol_vec[sol_index]->job_vec[oper->job_i]->oper_vec[oper->oper_i + 1]->in_degree == 0 &&
			oper->next_mach_oper != sol_vec[sol_index]->job_vec[oper->job_i]->oper_vec[oper->oper_i + 1])	// avoid the next machine operation is the same as the next job operation
			oper_in_degree_vec.push_back(sol_vec[sol_index]->job_vec[oper->job_i]->oper_vec[oper->oper_i + 1]);
		oper_re_top_order_stack.push(oper_in_degree_vec.front());	// push the first element into stack
		oper_in_degree_vec.erase(oper_in_degree_vec.begin());	// delete the first element
	}
	cout << endl;
	while (!oper_re_top_order_stack.empty())
	{
		Solution::Operation *oper = oper_re_top_order_stack.top();
		oper->q = MAX(oper->next_mach_oper->q + oper->next_mach_oper->t,
			sol_vec[sol_index]->job_vec[oper->job_i]->oper_vec[oper->oper_i + 1]->q + sol_vec[sol_index]->job_vec[oper->job_i]->oper_vec[oper->oper_i + 1]->t);
		oper_re_top_order_stack.pop();
		if (oper->r + oper->q + oper->t == sol_vec[sol_index]->makespan)
			critical_oper_stack.push(oper);
	}
	cout << "critical operations: ";
	int pre_mach = !critical_oper_stack.empty() ? critical_oper_stack.top()->mach : 0;
	while (!critical_oper_stack.empty())
	{
		Solution::Operation *oper = critical_oper_stack.top();
		if (oper->mach != pre_mach)
		{
			cout << endl;
			pre_mach = oper->mach;
		}
		cout << oper->job_i << ", " << oper->oper_i << "\t" << oper->mach << "\t"
			<< oper->start_time << "\t" << oper->end_time << "\t";
		critical_oper_stack.pop();
	}
	if (num_oper_visited != instance->total_num_operation)
		cout << "ERROR: not all operations are visited" << endl;
}
Instance::Instance(string _file_input):file_input(_file_input),total_num_operation(0)
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
		total_num_operation += jobinfo->num_oper;
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
	solver->check_solution(1);
	solver->determine_critical_path(1);
#ifdef _WIN32
	system("pause");
#endif
	return 0;
}
#endif