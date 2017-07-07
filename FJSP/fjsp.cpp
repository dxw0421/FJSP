#if 1
#include<iostream>
#include<fstream>
#include<vector>
#include<list>
#include<stack>
#include<map>
#include<string>
#include<algorithm>
#include<time.h>
using namespace std;
#define MAX(x,y) (((x)>(y))?(x):(y))
#define MAXN 52	// max job
#define MAXM 50	// max machine	
#define MAXS 10	// max solution
#define MAXO 60	// max operation
enum MOVE_TYPE { FORWARD_INSERT, BACKWARD_INSERT };
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
}
class Instance
{
public:
	class Process {
	public:
		int mach_i;
		int t;
		Process() {}
		Process(int _t, int _mach) :mach_i(_mach), t(_t) {}
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
	void display_instance()const;
};
class Solution
{
public:
	class Operation
	{
	public:
		const int job_i, oper_job_i;	// the job and operation index of this operation
		int mach_i, oper_mach_i, t, proc_i, start_time, end_time, r, q, apx_r, apx_q, in_degree;	// optimize by removing proc_i or mach and t
		Operation *pre_mach_oper, *next_mach_oper, *pre_job_oper, *next_job_oper;
		Operation(int _j, int _o, int _et) :end_time(_et), job_i(_j), oper_job_i(_o), pre_mach_oper(NULL), next_mach_oper(NULL), in_degree(2) {}
	};
	class JobInfo
	{
	public:
		vector<Operation*> oper_vec;
	};
	class Tabu
	{
	public:
		int tabu_iteration;	// the tabu iteration
		int oper_mach_i;	// operation u at machine position
		MOVE_TYPE move_tpye;	// move type
		vector<Operation *> tabu_oper_vec;	// the tabu block of operations from u to v
	};
	vector<JobInfo*> job_vec;
	vector<Operation*>head_oper_mach_vec, tail_oper_mach_vec;	// the operations assigned at each machine
	list<Tabu*> tabu_list;
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
	vector<Solution::Operation*> critical_block_vec;
	vector<int> critical_block_mach_pos_vec;
	const Instance *instance;

	enum DUMMY_OPER { START, END };
	Solution::Operation *job[MAXS][MAXN][MAXO];
	Solution::Operation *machine[MAXS][MAXM][MAXO];
	Solution::Operation *dummy_oper[MAXS][sizeof(DUMMY_OPER)];
	Solution::Operation *last_oper_makespan_machine[MAXS][MAXM];
	Solution::Operation *crit_path[MAXS][MAXM*MAXO];
	Solution::Operation *queue[MAXN*MAXO];
	int job_oper_num[MAXS][MAXN], machine_oper_num[MAXS][MAXM];
	int makespan[MAXS], last_oper_makespan_machine_num[MAXS], crit_path_num[MAXS],queue_num;
	int critical_flag[MAXN][MAXO];
	int crit_block[MAXN*MAXM][3];

	clock_t start_time, end_time;
	class Tabu
	{
	public:
		int tabu_iteration;	// the tabu iteration
		int tabu_u, tabu_oper_num;	// the number of tabu operation
		Solution::Operation *tabu_oper[MAXO];	// the tabu block of operations from u to v
		Tabu *front, *next;
	};
	Tabu *tabu_list[MAXS][MAXM];

	int tt0, d1, d2;	// for tabu tenure
	int best_known_makespan, iteration, alg_best_makespan, alg_best_iter;
	//int globel_iter;

	Solver(const Instance &, int);
	Solver(const Instance &, int, int);
	void display_machine_operation(int, int)const;
	void read_solution(int, string);
	void init_solution(int);
	void init_solution1(int);
	void determine_critical_path(int);
	void display_solution(int)const;
	void check_solution(int)const;
	void check_solution1(int);
	void backward_insert_move(int);
	void forward_insert_move(int);
	void insert_move(int);
	void insert_move1(int);
	void try_backward_insert_move(int, int&, Solution::Operation*, Solution::Operation*);
	void try_backward_insert_move1(int, int&, int, int, int);
	void try_forward_insert_move(int, int&, Solution::Operation*, Solution::Operation*);
	void try_forward_insert_move1(int, int&, int, int, int);
	void backward_insert(int, Solution::Operation*, Solution::Operation*);
	void forward_insert(int, Solution::Operation*, Solution::Operation*);
	void apply_move(int, int, int, int, MOVE_TYPE);
	void update_tabu(int, int, int, Solution::Operation*, Solution::Operation*, MOVE_TYPE);
	void update_tabu1(int, int, int, int,int);
	bool check_tabu(int, int, int, int, Solution::Operation*, Solution::Operation*, MOVE_TYPE);
	bool check_tabu1(int, int, int, int,int, MOVE_TYPE);
	bool check_cycle(int);
	void calculate_qr_obj(int);
	void calculate_q_crit_block(int);
	void perturb(int, int);
	void perturb1(int, int, int);
};
Instance::Instance(string _file_input) :file_input(_file_input), total_num_operation(0)
{
	ifstream ifs(file_input);
	if (!ifs.is_open())
	{
		cout << file_input << endl;
		perror("file_input"); exit(0);
	}
	string str_line;
	vector<string>fields_vec;
	if (file_input.find(".fjs") != string::npos)	// for flexible job shop instance
	{
		ifs >> n >> m >> avg_num_mach_per_oper;
		cout << n << "\t" << m << "\t" << avg_num_mach_per_oper << endl;
		getline(ifs, str_line);	// read the '\n'
		job_vec.push_back(NULL);	// the first element in job_vec is NULL
		while (getline(ifs, str_line))
		{
			//cout << str_line << endl;
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
	}
	else // for job shop instance
	{
		ifs >> n >> m;
		cout << n << "\t" << m << "\t" << avg_num_mach_per_oper << endl;
		getline(ifs, str_line);	// read the '\n'
		job_vec.push_back(NULL);	// the first element in job_vec is NULL
		while (getline(ifs, str_line))
		{
			cout << str_line << endl;
			split_generic<string>(fields_vec, str_line, " \t");
			JobInfo *jobinfo = new JobInfo();
			jobinfo->oper_vec.push_back(NULL);	// the first element in oper_vec is NULL
			int field_cnt = 0;
			jobinfo->num_oper = m;	// stoi(fields_vec[field_cnt++]);	// number of operation
			total_num_operation += jobinfo->num_oper;
			for (int i = 0; i < jobinfo->num_oper; i++)	// for each operation
			{
				Operation *oper = new Operation();
				oper->proc_vec.push_back(NULL);	// the first element in proc_vec is NULL
				oper->num_mach = 1;	// stoi(fields_vec[field_cnt++]);	// number of process
				for (int j = 0; j < oper->num_mach; j++)
				{
					Process *proc = new Process(stoi(fields_vec[field_cnt++]), stoi(fields_vec[field_cnt++]) + 1);	// machine index +1
					oper->proc_vec.push_back(proc);
				}
				jobinfo->oper_vec.push_back(oper);
			}
			job_vec.push_back(jobinfo);
		}
	}
	ifs.close();
}
void Instance::display_instance()const
{
	cout << "*** display instance: " << file_input << " ***" << endl;
	cout << n << "\t" << m << "\t" << avg_num_mach_per_oper << endl;
	for (vector<JobInfo*>::const_iterator job_iter = job_vec.begin() + 1;
		job_iter != job_vec.end(); job_iter++)
	{
		cout << (*job_iter)->num_oper << " | ";
		for (vector<Operation*>::iterator oper_iter = (*job_iter)->oper_vec.begin() + 1;
			oper_iter != (*job_iter)->oper_vec.end(); oper_iter++)
		{
			cout << (*oper_iter)->num_mach << " : ";
			for (vector<Process*>::iterator proc_iter = (*oper_iter)->proc_vec.begin() + 1;
				proc_iter != (*oper_iter)->proc_vec.end(); proc_iter++)
			{
				cout << (*proc_iter)->mach_i << "\t" << (*proc_iter)->t << "\t";
			}
		}
		cout << endl;
	}
}
Solver::Solver(const Instance &_instance, int _sol_num) :instance(&_instance), sol_num(_sol_num)
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
				oper->pre_job_oper = jobinfo->oper_vec[k - 1];
				if (k != 1)	// except the start dummy operation
					oper->pre_job_oper->next_job_oper = oper;
				jobinfo->oper_vec.push_back(oper);
			}
			jobinfo->oper_vec[1]->in_degree -= 1;	// the degree of the first operation of the job is 1
			jobinfo->oper_vec.back()->next_job_oper = sol->end_dummy_oper;	// assign the next job operation for the last operation of the job 
			jobinfo->oper_vec.push_back(sol->end_dummy_oper);	// the last element in oper_vec of this job is dummy operation
			sol->job_vec.push_back(jobinfo);
		}
		sol_vec.push_back(sol);
	}
}
Solver::Solver(const Instance &_instance, int _sol_num, int _temp) :instance(&_instance), sol_num(_sol_num)
{
	for (int i = 0; i <= sol_num; i++)
	{
		dummy_oper[i][START] = new Solution::Operation(0, 0, 0);
		dummy_oper[i][END] = new Solution::Operation(0, 1000, 0);
		for (int j = 1; j < instance->job_vec.size(); j++)
		{
			job[i][j][0] = dummy_oper[i][START];
			int k;
			for (k = 1; k < instance->job_vec[j]->oper_vec.size(); k++)
			{
				job[i][j][k] = new Solution::Operation(j, k, 0);
				job[i][j][k]->pre_job_oper = job[i][j][k - 1];
				if (k != 1)
					job[i][j][k]->pre_job_oper->next_job_oper = job[i][j][k];
			}
			job[i][j][k] = dummy_oper[i][END];
			job[i][j][k - 1]->next_job_oper = dummy_oper[i][END];
			job_oper_num[i][j] = k - 1;
		}
		for (int j = 1; j <= instance->m; j++)
		{
			machine[i][j][0] = dummy_oper[i][START];
			machine_oper_num[i][j] = 0;
		}
	}
}
void Solver::read_solution(int sol_index, string file_input)
{
	ifstream ifs(file_input);
	if (!ifs.is_open())
	{
		cout << file_input << endl;
		perror("file_input"); exit(0);
	}
	for (int job_i = 1; job_i <= instance->n; job_i++)
	{
		for (int oper_i = 1; oper_i <= instance->job_vec[job_i]->num_oper; oper_i++)
		{
			sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_i]->proc_i = 1;	// only one machine to process the operation
			sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_i]->mach_i = instance->job_vec[job_i]->oper_vec[oper_i]->proc_vec[1]->mach_i;
			sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_i]->t = instance->job_vec[job_i]->oper_vec[oper_i]->proc_vec[1]->t;
		}
	}
	string str_line;
	vector<string>fields_vec;
	getline(ifs, str_line);
	best_known_makespan = stoi(str_line);
	int mach_i = 0, makespan = 0;
	while (getline(ifs, str_line))
	{
		cout << str_line << endl;
		split_generic<string>(fields_vec, str_line, " \t");
		mach_i += 1;
		for (int i = 0; i<fields_vec.size(); i++)
		{
			int job_i = stoi(fields_vec[i]) + 1;
			Solution::Operation *oper = sol_vec[sol_index]->job_vec[job_i]->oper_vec[1];
			while (oper->mach_i != mach_i)
				oper = oper->next_job_oper;
			oper->pre_mach_oper = sol_vec[sol_index]->tail_oper_mach_vec[mach_i];
			if (sol_vec[sol_index]->head_oper_mach_vec[mach_i] == NULL)
				sol_vec[sol_index]->head_oper_mach_vec[mach_i] = oper;
			else
				oper->pre_mach_oper->next_mach_oper = oper;
			sol_vec[sol_index]->tail_oper_mach_vec[mach_i] = oper;
			/*oper->r = MAX(oper->pre_job_oper->r + oper->pre_job_oper->t, oper->pre_mach_oper->r + oper->pre_mach_oper->t);
			oper->start_time = oper->r;
			oper->end_time = oper->start_time + oper->t;
			if (oper->end_time > makespan)
			makespan = oper->end_time;*/
		}
		sol_vec[sol_index]->tail_oper_mach_vec[mach_i]->next_mach_oper = sol_vec[sol_index]->end_dummy_oper;
		/*if (sol_vec[sol_index]->tail_oper_mach_vec[mach_i]->end_time > makespan)
		makespan = sol_vec[sol_index]->tail_oper_mach_vec[mach_i]->end_time;*/
	}
	calculate_qr_obj(sol_index);
	if (best_known_makespan != sol_vec[sol_index]->makespan)
	{
		cout << "ERROR: the solution given by the author is wrong." << endl;
		system("pause");
	}
	ifs.close();
}
void Solver::display_solution(int sol_index)const
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
			cout << (*oper_iter)->mach_i << "\t" << (*oper_iter)->start_time << "\t" << (*oper_iter)->end_time << "\t";
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
			cout << oper_same_mach->job_i << ", " << oper_same_mach->oper_job_i << "\t"
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
			cout << oper_same_mach->job_i << ", " << oper_same_mach->oper_job_i << "\t"
				<< oper_same_mach->start_time << "\t" << oper_same_mach->end_time << "\t";
			oper_same_mach = oper_same_mach->pre_mach_oper;
		}
		cout << endl;
	}
}
void Solver::display_machine_operation(int sol_index, int mach_index)const
{
	cout << "solution " << sol_index << ", machine " << mach_index << ", operations: ";
	for (Solution::Operation *iter_oper = sol_vec[sol_index]->head_oper_mach_vec[mach_index];
		iter_oper != sol_vec[sol_index]->end_dummy_oper; iter_oper = iter_oper->next_mach_oper)
		cout << iter_oper->job_i << ", " << iter_oper->oper_job_i << "\t";
	cout << endl;
}
void Solver::check_solution(int sol_index)const
{
	cout << "***check solution " << sol_index << " *** " << sol_vec[sol_index]->makespan << endl;
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
			if (oper->start_time != MAX(oper->pre_job_oper->end_time, oper->pre_mach_oper->end_time))
			{
				cout << "ERROR: the start time is not max of the end time of previous operations of same job and machine" << endl;
				system("pause");
			}
			int proc_i_select = 0;
			for (int proc_i = 1; proc_i < instance->job_vec[oper->job_i]->oper_vec[oper->oper_job_i]->proc_vec.size(); proc_i++)
			{
				if (oper->mach_i == instance->job_vec[oper->job_i]->oper_vec[oper->oper_job_i]->proc_vec[proc_i]->mach_i)
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
			if (oper->end_time != oper->start_time + instance->job_vec[oper->job_i]->oper_vec[oper->oper_job_i]->proc_vec[proc_i_select]->t)
			{
				cout << "ERROR: the end time is not the start time plus its processing time" << endl;
				system("pause");
			}
			if (max_ct < oper->end_time)
			{
				max_ct = oper->end_time;
				max_ct_mach = oper->mach_i;
			}
			oper = oper->next_mach_oper;
		}
	}
	if (sol_vec[sol_index]->makespan != max_ct /*|| sol_vec[sol_index]->makespan_mach != max_ct_mach*/)
	{
		cout << "ERROR: makespan or makespan machine is wrong" << endl;
		system("pause");
	}
}
void Solver::check_solution1(int sol_index)
{
	cout << "***check solution " << sol_index << " *** " << makespan[sol_index] << endl;
	stack<Solution::Operation*> top_oper_stack;
	memset(critical_flag, 0, sizeof(critical_flag));
	int oper_num = 0;	// optimize by removing queue_num
	int front = 0, rear = 0;
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		if (machine[sol_index][mach_i][1]->oper_job_i==1)	// alternative
		{
			queue[rear++] = machine[sol_index][mach_i][1];
			critical_flag[mach_i][1] = 0;
		}else
			critical_flag[mach_i][1] = 1;
		for (int oper_i = 2; oper_i <= machine_oper_num[sol_index][mach_i]; oper_i++)
		{
			if (machine[sol_index][mach_i][oper_i]->oper_job_i == 1)
				critical_flag[mach_i][oper_i] = 1;
			else
				critical_flag[mach_i][oper_i] = 2;
		}
		oper_num += machine_oper_num[sol_index][mach_i];
	}
	if (oper_num != instance->m*instance->n)
	{
		cout << "ERROR: the number of operations is wrong" << endl;
		system("pause");
	}
	oper_num = 0;
	while (front != rear)
	{
		Solution::Operation *oper = queue[front++];
		int start_tm = MAX(oper->pre_job_oper->end_time, machine[sol_index][oper->mach_i][oper->oper_mach_i - 1]->end_time);
		int end_tm = start_tm + oper->t;
		top_oper_stack.push(oper);
		oper_num += 1;
		if (oper->start_time != start_tm)
		{
			cout << "ERROR: the start time is not max of the end time of previous operations of same job and machine" << endl;
			system("pause");
		}
		if (oper->end_time != end_tm)
		{
			cout << "ERROR: the end time is wrong" << endl;
			system("pause");
		}
		critical_flag[oper->next_job_oper->mach_i][oper->next_job_oper->oper_mach_i] -= 1;
		critical_flag[oper->mach_i][oper->oper_mach_i + 1] -= 1;
		if (critical_flag[oper->next_job_oper->mach_i][oper->next_job_oper->oper_mach_i] == 0)
			queue[rear++] = oper->next_job_oper;
		if (critical_flag[oper->mach_i][oper->oper_mach_i + 1] == 0 &&
			machine[sol_index][oper->mach_i][oper->oper_mach_i + 1] != oper->next_job_oper)
			queue[rear++] = machine[sol_index][oper->mach_i][oper->oper_mach_i + 1];
	}
	if (oper_num != instance->m*instance->n)
	{
		cout << "ERROR: the number of operations is wrong" << endl;
		system("pause");
	}

	while (!top_oper_stack.empty())
	{
		Solution::Operation *oper = top_oper_stack.top();
		top_oper_stack.pop();
		int real_q = MAX(oper->next_job_oper->q + oper->next_job_oper->t, 
			machine[sol_index][oper->mach_i][oper->oper_mach_i + 1]->q + machine[sol_index][oper->mach_i][oper->oper_mach_i + 1]->t);
		if (oper->q != real_q)
		{
			cout << "ERROR: q is wrong" << endl;
			system("pause");
		}
		if (oper->q + oper->end_time == makespan[sol_index])
			critical_flag[oper->mach_i][oper->oper_mach_i] = 1;
	}

	int block_cnt = 0, block_mach_i, first_i, last_i;
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		for (int oper_i = 1; oper_i <= machine_oper_num[sol_index][mach_i]; oper_i++)
		{
			if (critical_flag[mach_i][oper_i] == 1)
			{
				block_cnt += 1;
				block_mach_i = mach_i;
				first_i = oper_i;
				oper_i += 1;
				while (critical_flag[mach_i][oper_i] == 1 && oper_i <= machine_oper_num[sol_index][mach_i])
					oper_i += 1;
				last_i = oper_i - 1;
				if (first_i == last_i)
					block_cnt -= 1;
				else
				{
					if (crit_block[block_cnt][0]!=block_mach_i||crit_block[block_cnt][1] != first_i || crit_block[block_cnt][2] != last_i)
					{
						cout << "ERROR: critical block is wrong" << endl;
						system("pause");
					}
				}
			}
		}
	}
	if (crit_block[0][0] != block_cnt)
	{
		cout << "ERROR: number of critical block is wrong" << endl;
		system("pause");
	}
	int max_ct = 0;
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		if (max_ct < machine[sol_index][mach_i][machine_oper_num[sol_index][mach_i]]->end_time)
			max_ct = machine[sol_index][mach_i][machine_oper_num[sol_index][mach_i]]->end_time;
	}
	if (makespan[sol_index] != max_ct)
	{
		cout << "ERROR: makespan or makespan machine is wrong" << endl;
		system("pause");
	}
}
bool Solver::check_cycle(int sol_index)
{
	cout << "***check solution " << sol_index << " *** " << makespan[sol_index] << endl;
	//memset(critical_flag, 0, sizeof(critical_flag));
	int oper_num = 0;	// optimize by removing queue_num
	int front = 0, rear = 0;
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		if (machine[sol_index][mach_i][1]->pre_job_oper == dummy_oper[sol_index][START])	// alternative
		{
			queue[rear++] = machine[sol_index][mach_i][1];
			critical_flag[mach_i][1] = 0;
		}
		else
			critical_flag[mach_i][1] = 1;
		for (int oper_i = 2; oper_i <= machine_oper_num[sol_index][mach_i]; oper_i++)
		{
			if (machine[sol_index][mach_i][oper_i]->pre_job_oper == dummy_oper[sol_index][START])
				critical_flag[mach_i][oper_i] = 1;
			else
				critical_flag[mach_i][oper_i] = 2;
		}
		oper_num += machine_oper_num[sol_index][mach_i];
	}
	if (oper_num != instance->m*instance->n)
	{
		cout << "ERROR: the number of operations is wrong" << endl;
		return true;
	}
	oper_num = 0;
	while (front != rear)
	{
		Solution::Operation *oper = queue[front++];
		oper_num += 1;
		critical_flag[oper->next_job_oper->mach_i][oper->next_job_oper->oper_mach_i] -= 1;
		critical_flag[oper->mach_i][oper->oper_mach_i + 1] -= 1;
		if (critical_flag[oper->next_job_oper->mach_i][oper->next_job_oper->oper_mach_i] == 0)
			queue[rear++] = oper->next_job_oper;
		if (critical_flag[oper->mach_i][oper->oper_mach_i + 1] == 0 &&
			machine[sol_index][oper->mach_i][oper->oper_mach_i + 1] != oper->next_job_oper)
			queue[rear++] = machine[sol_index][oper->mach_i][oper->oper_mach_i + 1];
	}
	if (oper_num != instance->m*instance->n)
	{
		cout << "ERROR: the number of operations is wrong" << endl;
		return true;
	}
	return false;
}
void Solver::init_solution(int sol_index)
{
	vector<int> mach_ct_vec(instance->job_vec.size(), 0);	// the completion time of each machine
															//vector<Solution::Operation*> last_oper_mach(instance->m + 1,sol_vec[sol_index]->start_dummy_oper);	// the last operation at each machine
	vector<int> oper_to_assign_job(instance->job_vec.size(), 1);	// the operation index needed to assign at each job
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
			Solution::Operation * cur_oper = sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i]];
			int min_ct = INT_MAX, min_ct_mach, select_proc_i, equ_cnt = 1;
			for (int proc_i = 1; proc_i < oper->proc_vec.size(); proc_i++)	// find the best machine which has the minimum completion time
			{
				int cur_mach_ct = oper->proc_vec[proc_i]->t + MAX(mach_ct_vec[oper->proc_vec[proc_i]->mach_i],	// the next operation with the same machine
					sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i] - 1]->end_time);		// the next operation with the same job
				if (cur_mach_ct < min_ct)
				{
					min_ct = cur_mach_ct;
					select_proc_i = proc_i;
					min_ct_mach = oper->proc_vec[proc_i]->mach_i;
				}
				else if (cur_mach_ct == min_ct)
				{
					equ_cnt += 1;
					if (rand() % equ_cnt)
					{
						select_proc_i = proc_i;
						min_ct_mach = oper->proc_vec[proc_i]->mach_i;
					}
				}
			}
			cur_oper->proc_i = select_proc_i;
			cur_oper->mach_i = min_ct_mach;
			cur_oper->t = oper->proc_vec[select_proc_i]->t;
			cur_oper->start_time = MAX(mach_ct_vec[min_ct_mach],
				sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_to_assign_job[job_i] - 1]->end_time);
			cur_oper->end_time = min_ct;
			cur_oper->pre_mach_oper = sol_vec[sol_index]->tail_oper_mach_vec[min_ct_mach];
			cur_oper->next_mach_oper = sol_vec[sol_index]->end_dummy_oper;	// the next operation with the same machine is dummy end operation
			mach_ct_vec[min_ct_mach] = min_ct;
			//sol_vec[sol_index]->tail_oper_mach_vec[min_ct_mach] = cur_oper;	// the tail operation assigned to the machine
			if (sol_vec[sol_index]->head_oper_mach_vec[min_ct_mach] == NULL)	// assign the first operation to the machine
			{
				sol_vec[sol_index]->head_oper_mach_vec[min_ct_mach] = cur_oper;
				cur_oper->in_degree -= 1;	// the degree of the first opertion at the machine is 1
			}
			if (sol_vec[sol_index]->tail_oper_mach_vec[min_ct_mach]->job_i != 0)	// this is not the first operation of the machine, assgin the next operation of the previous operation
				cur_oper->pre_mach_oper->next_mach_oper = cur_oper;
			sol_vec[sol_index]->tail_oper_mach_vec[min_ct_mach] = cur_oper;	// update the tail operation at this machine
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
void Solver::init_solution1(int sol_index)
{
	vector<int> oper_to_assign_job(instance->job_vec.size(), 1);	// the operation index needed to assign at each job
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
			Solution::Operation * cur_oper = job[sol_index][job_i][oper_to_assign_job[job_i]];
			int min_ct = INT_MAX, min_ct_mach, select_proc_i, equ_cnt = 1;
			for (int proc_i = 1; proc_i < oper->proc_vec.size(); proc_i++)	// find the best machine which has the minimum completion time
			{
				int oper_i = machine_oper_num[sol_index][oper->proc_vec[proc_i]->mach_i];
				int cur_mach_ct = oper->proc_vec[proc_i]->t + MAX(machine[sol_index][oper->proc_vec[proc_i]->mach_i][oper_i]->end_time,	// the previous operation with the same machine
					job[sol_index][job_i][oper_to_assign_job[job_i] - 1]->end_time);							// the previous operation with the same job
				if (cur_mach_ct < min_ct)
				{
					min_ct = cur_mach_ct;
					select_proc_i = proc_i;
					min_ct_mach = oper->proc_vec[proc_i]->mach_i;
				}
				else if (cur_mach_ct == min_ct)
				{
					equ_cnt += 1;
					if (rand() % equ_cnt)
					{
						select_proc_i = proc_i;
						min_ct_mach = oper->proc_vec[proc_i]->mach_i;
					}
				}
			}
			cur_oper->proc_i = select_proc_i;
			cur_oper->mach_i = min_ct_mach;
			cur_oper->t = oper->proc_vec[select_proc_i]->t;
			cur_oper->start_time = MAX(machine[sol_index][min_ct_mach][machine_oper_num[sol_index][min_ct_mach]]->end_time,
				job[sol_index][job_i][oper_to_assign_job[job_i] - 1]->end_time);
			cur_oper->end_time = min_ct;

			machine_oper_num[sol_index][min_ct_mach] += 1;
			cur_oper->oper_mach_i = machine_oper_num[sol_index][min_ct_mach];	// the position at the assigned machine of the operation
			machine[sol_index][min_ct_mach][machine_oper_num[sol_index][min_ct_mach]] = cur_oper;

			oper_to_assign_job[job_i] += 1;	// the next operation is to be assigned
			if (oper_to_assign_job[job_i] >= instance->job_vec[job_i]->oper_vec.size())
				oper_to_assign_job[job_i] = 0;	// All the operations of this job is completed
		}
	}
	makespan[sol_index] = 0;
	for (int i = 1; i <= instance->m; i++)
	{
		machine[sol_index][i][machine_oper_num[sol_index][i] + 1] = dummy_oper[sol_index][END];
		if (makespan[sol_index] < machine[sol_index][i][machine_oper_num[sol_index][i]]->end_time)
		{
			makespan[sol_index] = machine[sol_index][i][machine_oper_num[sol_index][i]]->end_time;
			last_oper_makespan_machine_num[sol_index] = 0;	// begin from 1
			last_oper_makespan_machine[sol_index][last_oper_makespan_machine_num[sol_index]++] = machine[sol_index][i][machine_oper_num[sol_index][i]];
		}
		else if (makespan[sol_index] == machine[sol_index][i][machine_oper_num[sol_index][i]]->end_time)
			last_oper_makespan_machine[sol_index][last_oper_makespan_machine_num[sol_index]++] = machine[sol_index][i][machine_oper_num[sol_index][i]];
	}
}
void Solver::determine_critical_path(int sol_index)	// optimize by removing this function
{
	cout << "***solution " << sol_index << ", critical path***" << endl;
	vector<Solution::Operation *> oper_in_degree_vec;	// store operations with in-degree of 0
	stack<Solution::Operation*> oper_re_top_order_stack;	// store operations in stack to get re-topological order
	stack<Solution::Operation*> critical_oper_stack;	// store critical operations
	for (int mach_i = 1; mach_i < sol_vec[sol_index]->head_oper_mach_vec.size(); mach_i++)
	{
		Solution::Operation *oper = sol_vec[sol_index]->head_oper_mach_vec[mach_i];
		if (oper->oper_job_i == 1)	// or oper->oper_job_i == 1, which means the first operation at this job, in-degree is 0, optimize by removing in_degree
		{
			oper->in_degree = 0;
			oper_in_degree_vec.push_back(oper);
		}
		else
			oper->in_degree = 1;
		oper = oper->next_mach_oper;
		while (oper != sol_vec[sol_index]->end_dummy_oper)
		{
			if (oper->oper_job_i == 1)
				oper->in_degree = 1;
			else
				oper->in_degree = 2;
			oper = oper->next_mach_oper;
		}
	}
	int num_oper_visited = 0;
	while (!oper_in_degree_vec.empty())
	{
		Solution::Operation *oper = oper_in_degree_vec.front();
		/*cout << oper->job_i << ", " << oper->oper_job_i << "\t" << oper->mach_i << "\t"
		<< oper->start_time << "\t" << oper->end_time << "\t";*/
		oper->r = MAX(oper->pre_mach_oper->r + oper->pre_mach_oper->t, oper->pre_job_oper->r + oper->pre_job_oper->t);
		if (oper->r != oper->start_time)
		{
			system("pause");
		}
		num_oper_visited += 1;
		oper->next_mach_oper->in_degree -= 1;	// minus the in-degree of next operation at the same machine by 1
		oper->next_job_oper->in_degree -= 1;	// minus the in-degree of next operation of the same job by 1
		if (oper->next_mach_oper->in_degree == 0)
			oper_in_degree_vec.push_back(oper->next_mach_oper);
		if (oper->next_job_oper->in_degree == 0 && oper->next_mach_oper != oper->next_job_oper)	// avoid the next machine operation is the same as the next job operation
			oper_in_degree_vec.push_back(oper->next_job_oper);
		oper_re_top_order_stack.push(oper_in_degree_vec.front());	// push the first element into stack
		oper_in_degree_vec.erase(oper_in_degree_vec.begin());	// delete the first element
	}
	//cout << endl;
	while (!oper_re_top_order_stack.empty())
	{
		Solution::Operation *oper = oper_re_top_order_stack.top();
		oper->q = MAX(oper->next_mach_oper->q + oper->next_mach_oper->t, oper->next_job_oper->q + oper->next_job_oper->t);
		oper_re_top_order_stack.pop();
		if (oper->r + oper->q + oper->t == sol_vec[sol_index]->makespan)
			critical_oper_stack.push(oper);
	}
	cout << "critical blocks in determine: " << endl;
	critical_block_vec.clear();
	Solution::Operation *pre_oper, *cur_oper, *first_oper;
	first_oper = !critical_oper_stack.empty() ? critical_oper_stack.top() : NULL;
	while (!critical_oper_stack.empty())
	{
		cur_oper = critical_oper_stack.top();
		if (cur_oper->mach_i != first_oper->mach_i)
		{
			cout << endl;
			critical_block_vec.push_back(first_oper);
			critical_block_vec.push_back(pre_oper);
			first_oper = cur_oper;
		}
		pre_oper = cur_oper;
		cout << cur_oper->job_i << ", " << cur_oper->oper_job_i << "\t" << cur_oper->mach_i << "\t"
			<< cur_oper->start_time << "\t" << cur_oper->end_time << endl;
		critical_oper_stack.pop();
	}
	cout << endl;
	if (num_oper_visited != instance->total_num_operation)
	{
		cout << "ERROR: not all operations are visited" << endl;
		system("pause");
	}
}
void Solver::try_backward_insert_move(int sol_index, int &makespan, Solution::Operation *oper_u, Solution::Operation *oper_v)	// insert oper_u behind oper_v
{
	oper_u->next_mach_oper->apx_r = MAX(oper_u->next_mach_oper->pre_job_oper->r + oper_u->next_mach_oper->pre_job_oper->t, oper_u->pre_mach_oper->r + oper_u->pre_mach_oper->t);
	for (Solution::Operation *oper = oper_u->next_mach_oper->next_mach_oper; oper != oper_v->next_mach_oper; oper = oper->next_mach_oper)
		oper->apx_r = MAX(oper->pre_job_oper->r + oper->pre_job_oper->t, oper->pre_mach_oper->apx_r + oper->pre_mach_oper->t);
	oper_u->apx_r = MAX(oper_u->pre_job_oper->r + oper_u->pre_job_oper->t, oper_v->apx_r + oper_v->t);
	oper_u->apx_q = MAX(oper_u->next_job_oper->q + oper_u->next_job_oper->t, oper_v->next_mach_oper->q + oper_v->next_mach_oper->t);
	oper_v->apx_q = MAX(oper_v->next_job_oper->q + oper_v->next_job_oper->t, oper_u->apx_q + oper_u->t);
	for (Solution::Operation *oper = oper_v->pre_mach_oper; oper != oper_u; oper = oper->pre_mach_oper)
		oper->apx_q = MAX(oper->next_job_oper->q + oper->next_job_oper->t, oper->next_mach_oper->apx_q + oper->next_mach_oper->t);
	makespan = 0;
	for (Solution::Operation *oper = oper_u; oper != oper_v->next_mach_oper; oper = oper->next_mach_oper)
		if (oper->apx_r + oper->apx_q + oper->t > makespan)
			makespan = oper->apx_r + oper->apx_q + oper->t;
}
void Solver::try_backward_insert_move1(int sol_index, int &makespan, int mach_i, int u, int v)	// insert oper_u behind oper_v
{
	machine[sol_index][mach_i][u + 1]->apx_r = MAX(machine[sol_index][mach_i][u + 1]->pre_job_oper->end_time, machine[sol_index][mach_i][u - 1]->end_time);
	for (int oper_i = u + 2; oper_i <= v; oper_i++)
		machine[sol_index][mach_i][oper_i]->apx_r = MAX(machine[sol_index][mach_i][oper_i]->pre_job_oper->end_time, machine[sol_index][mach_i][oper_i - 1]->apx_r + machine[sol_index][mach_i][oper_i - 1]->t);
	machine[sol_index][mach_i][u]->apx_r = MAX(machine[sol_index][mach_i][u]->pre_job_oper->end_time, machine[sol_index][mach_i][v]->apx_r + machine[sol_index][mach_i][v]->t);
	machine[sol_index][mach_i][u]->apx_q = MAX(machine[sol_index][mach_i][u]->next_job_oper->q + machine[sol_index][mach_i][u]->next_job_oper->t, machine[sol_index][mach_i][v + 1]->q + machine[sol_index][mach_i][v + 1]->t);
	machine[sol_index][mach_i][v]->apx_q = MAX(machine[sol_index][mach_i][v]->next_job_oper->q + machine[sol_index][mach_i][v]->next_job_oper->t,
		machine[sol_index][mach_i][u]->apx_q + machine[sol_index][mach_i][u]->t);
	for (int oper_i = v - 1; oper_i > u; oper_i--)
		machine[sol_index][mach_i][oper_i]->apx_q = MAX(machine[sol_index][mach_i][oper_i]->next_job_oper->q + machine[sol_index][mach_i][oper_i]->next_job_oper->t,
			machine[sol_index][mach_i][oper_i + 1]->apx_q + machine[sol_index][mach_i][oper_i + 1]->t);
	makespan = 0;
	for (int oper_i = u; oper_i <= v; oper_i++)
		if (makespan < machine[sol_index][mach_i][oper_i]->apx_r + machine[sol_index][mach_i][oper_i]->apx_q + machine[sol_index][mach_i][oper_i]->t)
			makespan = machine[sol_index][mach_i][oper_i]->apx_r + machine[sol_index][mach_i][oper_i]->apx_q + machine[sol_index][mach_i][oper_i]->t;
}
void Solver::try_forward_insert_move(int sol_index, int &makespan, Solution::Operation *oper_u, Solution::Operation *oper_v)	// insert oper_v before oper_u
{
	oper_v->apx_r = MAX(oper_v->pre_job_oper->r + oper_v->pre_job_oper->t, oper_u->pre_mach_oper->r + oper_u->pre_mach_oper->t);
	oper_u->apx_r = MAX(oper_u->pre_job_oper->r + oper_v->pre_job_oper->t, oper_v->apx_r + oper_v->t);
	for (Solution::Operation *oper = oper_u->next_mach_oper; oper != oper_v; oper = oper->next_mach_oper)
		oper->apx_r = MAX(oper->pre_job_oper->r + oper->pre_job_oper->t, oper->pre_mach_oper->apx_r + oper->pre_mach_oper->t);
	oper_v->pre_mach_oper->apx_q = MAX(oper_v->pre_mach_oper->next_job_oper->q + oper_v->pre_mach_oper->next_job_oper->t, oper_v->next_mach_oper->q + oper_v->next_mach_oper->t);
	for (Solution::Operation *oper = oper_v->pre_mach_oper->pre_mach_oper; oper != oper_u->pre_mach_oper; oper = oper->pre_mach_oper)
		oper->apx_q = MAX(oper->next_job_oper->q + oper->next_job_oper->t, oper->next_mach_oper->apx_q + oper->next_mach_oper->t);
	oper_v->apx_q = MAX(oper_v->next_job_oper->q + oper_v->next_job_oper->t, oper_u->apx_q + oper_u->t);
	makespan = 0;
	for (Solution::Operation *oper = oper_u; oper != oper_v->next_mach_oper; oper = oper->next_mach_oper)
		if (oper->apx_r + oper->apx_q + oper->t > makespan)
			makespan = oper->apx_r + oper->apx_q + oper->t;
}
void Solver::try_forward_insert_move1(int sol_index, int &makespan, int mach_i, int u, int v)	// insert oper_v before oper_u
{
	machine[sol_index][mach_i][v]->apx_r = MAX(machine[sol_index][mach_i][v]->pre_job_oper->end_time, machine[sol_index][mach_i][u - 1]->end_time);
	machine[sol_index][mach_i][u]->apx_r = MAX(machine[sol_index][mach_i][u]->pre_job_oper->end_time, machine[sol_index][mach_i][v]->apx_r + machine[sol_index][mach_i][v]->t);
	for (int oper_i = u + 1; oper_i < v; oper_i++)
		machine[sol_index][mach_i][oper_i]->apx_r = MAX(machine[sol_index][mach_i][oper_i]->pre_job_oper->end_time, machine[sol_index][mach_i][oper_i - 1]->apx_r + machine[sol_index][mach_i][oper_i - 1]->t);
	machine[sol_index][mach_i][v - 1]->apx_q = MAX(machine[sol_index][mach_i][v - 1]->next_job_oper->q + machine[sol_index][mach_i][v - 1]->next_job_oper->t,
		machine[sol_index][mach_i][v + 1]->q + machine[sol_index][mach_i][v + 1]->t);
	for (int oper_i = v - 2; oper_i >= u; oper_i--)
		machine[sol_index][mach_i][oper_i]->apx_q = MAX(machine[sol_index][mach_i][oper_i]->next_job_oper->q + machine[sol_index][mach_i][oper_i]->next_job_oper->t,
			machine[sol_index][mach_i][oper_i + 1]->apx_q + machine[sol_index][mach_i][oper_i + 1]->t);
	machine[sol_index][mach_i][v]->apx_q = MAX(machine[sol_index][mach_i][v]->next_job_oper->q + machine[sol_index][mach_i][v]->next_job_oper->t,
		machine[sol_index][mach_i][u]->apx_q + machine[sol_index][mach_i][u]->t);
	makespan = 0;
	for (int oper_i = u; oper_i <= v; oper_i++)
		if (makespan < machine[sol_index][mach_i][oper_i]->apx_r + machine[sol_index][mach_i][oper_i]->apx_q + machine[sol_index][mach_i][oper_i]->t)
			makespan = machine[sol_index][mach_i][oper_i]->apx_r + machine[sol_index][mach_i][oper_i]->apx_q + machine[sol_index][mach_i][oper_i]->t;
}
void Solver::backward_insert(int sol_index, Solution::Operation *oper_u, Solution::Operation *oper_v)
{
	//cout << "backward insert 1" << endl;
	if (oper_u->pre_mach_oper == sol_vec[sol_index]->start_dummy_oper)	// if oper_u is the first operation at the machine
		sol_vec[sol_index]->head_oper_mach_vec[oper_u->mach_i] = oper_u->next_mach_oper;
	if (oper_v->next_mach_oper == sol_vec[sol_index]->end_dummy_oper)	// if oper_v is the last operation at the machine
		sol_vec[sol_index]->tail_oper_mach_vec[oper_u->mach_i] = oper_u;
	oper_u->pre_mach_oper->next_mach_oper = oper_u->next_mach_oper;
	oper_u->next_mach_oper->pre_mach_oper = oper_u->pre_mach_oper;	// for removing u
	oper_u->next_mach_oper = oper_v->next_mach_oper;
	oper_v->next_mach_oper->pre_mach_oper = oper_u;
	oper_v->next_mach_oper = oper_u;
	oper_u->pre_mach_oper = oper_v;	// for adding u
}
void Solver::forward_insert(int sol_index, Solution::Operation *oper_u, Solution::Operation *oper_v)
{
	//cout << "forward insert 1" << endl;
	if (oper_u->pre_mach_oper == sol_vec[sol_index]->start_dummy_oper)	// if oper_u is the first operation at the machine
		sol_vec[sol_index]->head_oper_mach_vec[oper_u->mach_i] = oper_v;
	if (oper_v->next_mach_oper == sol_vec[sol_index]->end_dummy_oper)	// if oper_v is the last operation at the machine
		sol_vec[sol_index]->tail_oper_mach_vec[oper_u->mach_i] = oper_v->pre_mach_oper;
	oper_v->pre_mach_oper->next_mach_oper = oper_v->next_mach_oper;
	oper_v->next_mach_oper->pre_mach_oper = oper_v->pre_mach_oper;
	oper_u->pre_mach_oper->next_mach_oper = oper_v;
	oper_v->pre_mach_oper = oper_u->pre_mach_oper;
	oper_v->next_mach_oper = oper_u;
	oper_u->pre_mach_oper = oper_v;
}
// calculate start and end time, makespan, last operation at makespan machine
void Solver::apply_move(int sol_index, int mach_i, int u, int v, MOVE_TYPE move_type)
{
	if (move_type == BACKWARD_INSERT)
	{
		Solution::Operation *oper_u = machine[sol_index][mach_i][u];
		for (int i = u; i < v; i++)
		{
			machine[sol_index][mach_i][i] = machine[sol_index][mach_i][i + 1];
			machine[sol_index][mach_i][i]->oper_mach_i = i;
		}
		machine[sol_index][mach_i][v] = oper_u;
		machine[sol_index][mach_i][v]->oper_mach_i = v;
	}
	else
	{
		Solution::Operation *oper_v = machine[sol_index][mach_i][v];
		for (int i = v; i > u; i--)
		{
			machine[sol_index][mach_i][i] = machine[sol_index][mach_i][i - 1];
			machine[sol_index][mach_i][i]->oper_mach_i = i;
		}
		machine[sol_index][mach_i][u] = oper_v;
		machine[sol_index][mach_i][u]->oper_mach_i = u;
	}
	int start_line[MAXM];
	for (int i = 0; i < MAXM; i++)
		start_line[i] = 100;
	memset(critical_flag, 0, sizeof(critical_flag));

	int front = 0, rear = 0;
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		if (machine[sol_index][mach_i][1]->pre_job_oper == dummy_oper[sol_index][START])	// alternative
		{
			queue[rear++] = machine[sol_index][mach_i][1];
			critical_flag[mach_i][1] = 0;
		}
		else
			critical_flag[mach_i][1] = 1;
		for (int oper_i = 2; oper_i <= machine_oper_num[sol_index][mach_i]; oper_i++)
		{
			if (machine[sol_index][mach_i][oper_i]->pre_job_oper == dummy_oper[sol_index][START])
				critical_flag[mach_i][oper_i] = 1;
			else
				critical_flag[mach_i][oper_i] = 2;
		}
	}
	while (front != rear)
	{
		Solution::Operation *oper = queue[front++];
		oper->start_time = MAX(oper->pre_job_oper->end_time, machine[sol_index][oper->mach_i][oper->oper_mach_i - 1]->end_time);
		oper->end_time = oper->start_time + oper->t;
		critical_flag[oper->next_job_oper->mach_i][oper->next_job_oper->oper_mach_i] -= 1;
		critical_flag[oper->mach_i][oper->oper_mach_i + 1] -= 1;
		if (critical_flag[oper->next_job_oper->mach_i][oper->next_job_oper->oper_mach_i] == 0)
			queue[rear++] = oper->next_job_oper;
		if (critical_flag[oper->mach_i][oper->oper_mach_i + 1] == 0 &&
			machine[sol_index][oper->mach_i][oper->oper_mach_i + 1] != oper->next_job_oper)
			queue[rear++] = machine[sol_index][oper->mach_i][oper->oper_mach_i + 1];
	}
	makespan[sol_index] = 0;
	for (int i = 1; i <= instance->m; i++)
	{
		if (makespan[sol_index] < machine[sol_index][i][machine_oper_num[sol_index][i]]->end_time)
		{
			makespan[sol_index] = machine[sol_index][i][machine_oper_num[sol_index][i]]->end_time;
			last_oper_makespan_machine_num[sol_index] = 0;	
			last_oper_makespan_machine[sol_index][last_oper_makespan_machine_num[sol_index]++] = machine[sol_index][i][machine_oper_num[sol_index][i]];
		}
		else if (makespan[sol_index] == machine[sol_index][i][machine_oper_num[sol_index][i]]->end_time)
			last_oper_makespan_machine[sol_index][last_oper_makespan_machine_num[sol_index]++] = machine[sol_index][i][machine_oper_num[sol_index][i]];
	}
}
void Solver::update_tabu(int sol_index, int cur_iter, int mach_i, Solution::Operation *oper_u, Solution::Operation *oper_v, MOVE_TYPE move_type)
{
	Solution::Tabu *tabu = new Solution::Tabu();
	tabu->move_tpye = move_type;
	tabu->oper_mach_i = mach_i;
	int r1 = rand() / ((sol_vec[sol_index]->makespan - best_known_makespan) / d1);
	tabu->tabu_iteration = cur_iter + MAX(r1, d2);
	for (Solution::Operation *oper = oper_u; oper != oper_v->next_mach_oper; oper = oper->next_mach_oper)
		tabu->tabu_oper_vec.push_back(oper);
	sol_vec[sol_index]->tabu_list.push_back(tabu);
}
void Solver::update_tabu1(int sol_index, int cur_iter, int mach_i, int u,int v)
{
	Tabu *tabu = new Tabu();
	int r1 = MAX((makespan[sol_index] - best_known_makespan) / d1, d2);	// optimize by making it global
	tabu->tabu_iteration = cur_iter + tt0 + rand() / r1;
	tabu->tabu_u = u;
	tabu->tabu_oper_num = v - u + 1;
	for (int i = 1; i <= tabu->tabu_oper_num; i++)
		tabu->tabu_oper[i] = machine[sol_index][mach_i][u + i - 1];
	tabu->front = NULL; 
	tabu->next = tabu_list[sol_index][mach_i];
	tabu_list[sol_index][mach_i] = tabu;
	if (tabu->next != NULL)
		tabu->next->front = tabu;
}
bool Solver::check_tabu(int sol_index, int cur_iter, int oper_u_pos, int oper_v_pos, Solution::Operation* oper_u, Solution::Operation* oper_v, MOVE_TYPE move_type)// if u and v in tabu, return true, else return false
{
	/*if (cur_iter == 151 && oper_u_pos==4&&oper_v_pos==5)
	cout << endl;*/
	for (list<Solution::Tabu*>::iterator tabu_iter = sol_vec[sol_index]->tabu_list.begin();
		tabu_iter != sol_vec[sol_index]->tabu_list.end();)
	{
		if (cur_iter > (*tabu_iter)->tabu_iteration) // not in tabu status
		{
			delete *tabu_iter;
			tabu_iter = sol_vec[sol_index]->tabu_list.erase(tabu_iter);
			continue;
		}
		if ((*tabu_iter)->tabu_oper_vec.front()->mach_i != oper_u->mach_i)	// not the same machine
		{
			tabu_iter++;
			continue;
		}
		// in tabu status
		int f_pos = (*tabu_iter)->oper_mach_i;
		int e_pos = (*tabu_iter)->oper_mach_i + (*tabu_iter)->tabu_oper_vec.size() - 1;
		if (oper_u_pos > e_pos || oper_v_pos < f_pos)
		{
			tabu_iter++;
			continue;
		}
		bool is_same = true;
		if (move_type == BACKWARD_INSERT)
		{
			if (f_pos < oper_u_pos)
			{
				Solution::Operation *cur_oper = oper_u;
				for (int i = oper_u_pos - f_pos - 1; i >= 0; i -= 1)
				{
					cur_oper = cur_oper->pre_mach_oper;
					if ((*tabu_iter)->tabu_oper_vec[i] != cur_oper)
					{
						is_same = false;
						break;
					}
				}
				if (is_same)
				{
					for (int i = oper_u_pos - f_pos; i < (*tabu_iter)->tabu_oper_vec.size(); i++)
					{
						if (i == oper_u_pos - f_pos)
							cur_oper = oper_u->next_mach_oper;
						else if (i > oper_u_pos - f_pos&&i < oper_v_pos - f_pos)
							cur_oper = cur_oper->next_mach_oper;
						else if (i == oper_v_pos - f_pos)
							cur_oper = oper_u;
						else
							cur_oper = cur_oper->next_job_oper;
						if ((*tabu_iter)->tabu_oper_vec[i] != cur_oper)
						{
							is_same = false;
							break;
						}
						if (i == oper_v_pos - f_pos)
							cur_oper = oper_v;
					}
				}
			}
			else if (f_pos == oper_u_pos)
			{
				Solution::Operation *cur_oper = oper_u;
				for (int i = 0; i < (*tabu_iter)->tabu_oper_vec.size(); i++)
				{
					if (f_pos + i < oper_v_pos)
						cur_oper = cur_oper->next_mach_oper;
					else if (f_pos + i == oper_v_pos)
						cur_oper = oper_u;
					else
						cur_oper = cur_oper->next_job_oper;
					if ((*tabu_iter)->tabu_oper_vec[i] != cur_oper)
					{
						is_same = false;
						break;
					}
					if (f_pos + i == oper_v_pos)
						cur_oper = oper_v;
				}
			}
			else //if (f_pos > oper_u_pos)
			{
				Solution::Operation *cur_oper = oper_u;
				for (int i = 0; i < f_pos - oper_u_pos; i += 1)
					cur_oper = cur_oper->next_mach_oper;
				for (int i = 0; i < (*tabu_iter)->tabu_oper_vec.size(); i++)
				{
					if (f_pos + i < oper_v_pos)
						cur_oper = cur_oper->next_mach_oper;
					else if (f_pos + i == oper_v_pos)
						cur_oper = oper_u;
					else
						cur_oper = cur_oper->next_job_oper;
					if ((*tabu_iter)->tabu_oper_vec[i] != cur_oper)
					{
						is_same = false;
						break;
					}
					if (f_pos + i == oper_v_pos)
						cur_oper = oper_v;
				}
			}
		}
		else //	if (move_type == FORWARD_INSERT)
		{
			if (f_pos < oper_u_pos)
			{
				Solution::Operation *cur_oper = oper_u;
				for (int i = oper_u_pos - f_pos - 1; i >= 0; i -= 1)
				{
					cur_oper = cur_oper->pre_mach_oper;
					if ((*tabu_iter)->tabu_oper_vec[i] != cur_oper)
					{
						is_same = false;
						break;
					}
				}
				if (is_same)
				{
					int i = oper_u_pos - f_pos;
					if ((*tabu_iter)->tabu_oper_vec[i] != oper_v)	// v not the smae
						is_same = false;
					if (is_same)
					{
						i += 1;
						if (i < (*tabu_iter)->tabu_oper_vec.size() && (*tabu_iter)->tabu_oper_vec[i] != oper_u)	// u not the same
							is_same = false;
						if (is_same)
						{
							i += 1;
							cur_oper = oper_u->next_mach_oper;
							for (; i < (*tabu_iter)->tabu_oper_vec.size(); i++)
							{
								if ((*tabu_iter)->tabu_oper_vec[i] != cur_oper)
								{
									is_same = false;
									break;
								}
								cur_oper = cur_oper->next_mach_oper;
							}
						}
					}
				}
			}
			else if (f_pos == oper_u_pos)
			{
				int i = 0;
				if ((*tabu_iter)->tabu_oper_vec[i] != oper_v)
					is_same = false;
				if (is_same)
				{
					i += 1;
					if ((*tabu_iter)->tabu_oper_vec[i] != oper_u)
						is_same = false;
					if (is_same)
					{
						i += 1;
						Solution::Operation *cur_oper = oper_u->next_mach_oper;
						for (; i < (*tabu_iter)->tabu_oper_vec.size(); i++)
						{
							if ((*tabu_iter)->tabu_oper_vec[i] != cur_oper)
							{
								is_same = false;
								break;
							}
							cur_oper = cur_oper->next_mach_oper;
						}
					}
				}
			}
			else  //if (f_pos > oper_u_pos)
			{
				Solution::Operation *cur_oper = oper_u;
				for (int i = 0; i < f_pos - oper_u_pos - 1; i += 1)
					cur_oper = cur_oper->next_mach_oper;
				for (int i = 0; i < (*tabu_iter)->tabu_oper_vec.size(); i++)
				{
					if ((*tabu_iter)->tabu_oper_vec[i] != cur_oper)
					{
						is_same = false;
						break;
					}
					cur_oper = cur_oper->next_mach_oper;
				}
			}
		}
		if (is_same)
			return true;
		tabu_iter++;
	}
	return false;
}
bool Solver::check_tabu1(int sol_index, int cur_iter, int mach_i, int u, int v, MOVE_TYPE move_type)// if u and v in tabu, return true, else return false
{
	/*if (cur_iter == 151 && oper_u_pos==4&&oper_v_pos==5)
	cout << endl;*/
	int r1 = MAX((makespan[sol_index] - best_known_makespan) / d1, d2);
	for (Tabu *tabu=tabu_list[sol_index][mach_i];tabu!=NULL;tabu=tabu->next)
	{
		if (tabu->tabu_iteration<=cur_iter-r1) // not in tabu status
		{
			if (tabu == tabu_list[sol_index][mach_i])	// tabu is the first node
				tabu_list[sol_index][mach_i] = NULL;
			else
				tabu->front->next = NULL;
			Tabu *tabu1 = tabu->next;	// optimize here
			while (tabu1 != NULL)
			{
				delete tabu;
				tabu = tabu1;
				tabu1 = tabu1->next;
			}
			if (tabu != NULL)
				delete tabu;
			tabu = NULL;
			break;
		}
		else if (tabu->tabu_iteration < cur_iter)
			continue;
		int tabu_v = tabu->tabu_u + tabu->tabu_oper_num - 1;
		if (u > tabu_v || v < tabu->tabu_u)
			continue;
		bool is_same = true;
		if (move_type == BACKWARD_INSERT)
		{
			for (int i = 1; i <= tabu->tabu_oper_num; i++)
			{
				int tabu_u = tabu->tabu_u + i - 1;
				Solution::Operation *oper;
				if (tabu_u < u)
					oper = machine[sol_index][mach_i][tabu_u];
				else if (tabu_u >= u&&tabu_u < v)
					oper = machine[sol_index][mach_i][tabu_u + 1];
				else if (tabu_u == v)
					oper = machine[sol_index][mach_i][u];
				else  // if (tabu_u > v)
					oper = machine[sol_index][mach_i][tabu_u];
				if (oper != tabu->tabu_oper[i])
				{
					is_same = false;
					break;
				}
			}
		}
		else  // if(move_type==FORWARD_INSERT)
		{
			for (int i = 1; i <= tabu->tabu_oper_num; i++)
			{
				int tabu_u = tabu->tabu_u + i - 1;
				Solution::Operation *oper;
				if (tabu_u < u)
					oper = machine[sol_index][mach_i][tabu_u];
				else if (tabu_u == u)
					oper = machine[sol_index][mach_i][v];
				else if (tabu_u > u&&tabu_u <= v)
					oper = machine[sol_index][mach_i][tabu_u - 1];
				else  // if (tabu_u > v)
					oper = machine[sol_index][mach_i][tabu_u];
				if (oper != tabu->tabu_oper[i])
				{
					is_same = false;
					break;
				}
			}
		}
		if (is_same)
			return is_same;
	}
	return false;
}
void Solver::backward_insert_move(int sol_index)
{
	cout << "backward insert move, solution " << sol_index << ", " << sol_vec[sol_index]->makespan << endl;
	Solution::Operation *min_oper_u = NULL, *min_oper_v = NULL;
	int  min_makespan, equ_cnt;
	alg_best_makespan = sol_vec[sol_index]->makespan;
	for (int cur_iter = 1; cur_iter <= iteration; cur_iter++)
	{
		min_makespan = INT_MAX;
		for (vector<Solution::Operation *>::iterator oper_iter = critical_block_vec.begin(); oper_iter != critical_block_vec.end(); oper_iter += 2)
		{
			for (Solution::Operation *oper_u = *oper_iter; oper_u != *(oper_iter + 1); oper_u = oper_u->next_mach_oper)
			{
				for (Solution::Operation *oper_v = oper_u->next_mach_oper; oper_v != (*(oper_iter + 1))->next_mach_oper; oper_v = oper_v->next_mach_oper)	// move u behind v
				{
					if (oper_v->q >= oper_u->next_job_oper->q &&	// q[v]>=q[JS[u]]
						oper_v->job_i != oper_u->job_i)	// u, v do not belong to the same job
					{
						/*cout << oper_u->job_i << "*, " << oper_u->oper_job_i << "\t"
						<< oper_v->job_i << ", " << oper_v->oper_job_i << "\t";*/
						int makespan;
						try_backward_insert_move(sol_index, makespan, oper_u, oper_v);
						if (makespan <= min_makespan && (makespan < sol_vec[sol_index]->makespan/* ||
																								!check_tabu(sol_index, oper_u, oper_v, cur_iter, BACKWARD_INSERT)*/))
						{
							if (makespan < min_makespan)
							{
								min_makespan = makespan;
								min_oper_u = oper_u;
								min_oper_v = oper_v;
								equ_cnt = 1;
							}
							else if (makespan == min_makespan)
							{
								equ_cnt += 1;
								if (rand() % equ_cnt)
								{
									min_oper_u = oper_u;
									min_oper_v = oper_v;
								}
							}
						}
						//cout << makespan << "\t";
					}
				}
				//cout << endl;
			}
			//cout << endl;
		}
		if (min_makespan == INT_MAX)	// still_improve==false
			continue;
		cout << "The best move: ";
		cout << cur_iter << "\t"
			<< min_oper_u->job_i << ", " << min_oper_u->oper_job_i << "\t"
			<< min_oper_v->job_i << ", " << min_oper_v->oper_job_i << "\t"
			<< min_makespan << "\t";
		/*update_tabu(sol_index, min_oper_u, min_oper_v, cur_iter,BACKWARD_INSERT);*/
		//display_machine_operation(sol_index, min_oper_u->mach_i);
		backward_insert(sol_index, min_oper_u, min_oper_v);
		calculate_qr_obj(sol_index);
		if (alg_best_makespan> sol_vec[sol_index]->makespan)
			alg_best_makespan = sol_vec[sol_index]->makespan;
		//display_machine_operation(sol_index, min_oper_u->mach_i);
		cout << sol_vec[sol_index]->makespan << "\t" << alg_best_makespan << endl;
		check_solution(sol_index);
	}
}
void Solver::forward_insert_move(int sol_index)
{
	cout << "forward insert move, solution " << sol_index << ", " << sol_vec[sol_index]->makespan << endl;
	Solution::Operation *min_oper_u = NULL, *min_oper_v = NULL;
	int  min_makespan, equ_cnt;
	alg_best_makespan = sol_vec[sol_index]->makespan;
	for (int cur_iter = 1; cur_iter <= iteration; cur_iter++)
	{
		min_makespan = INT_MAX;
		for (vector<Solution::Operation *>::iterator oper_iter = critical_block_vec.begin(); oper_iter != critical_block_vec.end(); oper_iter += 2)
		{
			for (Solution::Operation *oper_u = *oper_iter; oper_u != *(oper_iter + 1); oper_u = oper_u->next_mach_oper)
			{
				for (Solution::Operation *oper_v = oper_u->next_mach_oper; oper_v != (*(oper_iter + 1))->next_mach_oper; oper_v = oper_v->next_mach_oper)	// move u behind v
				{
					if (oper_u->r + oper_u->t >= oper_v->pre_job_oper->r + oper_v->pre_job_oper->t &&	// r[u]+t[u]>=r[JP[v]]+t[JP[v]]
						oper_v->job_i != oper_u->job_i)	// u, v do not belong to the same job
					{
						/*cout << oper_u->job_i << "*, " << oper_u->oper_job_i << "\t"
						<< oper_v->job_i << ", " << oper_v->oper_job_i << "\t";*/
						int makespan;
						try_forward_insert_move(sol_index, makespan, oper_u, oper_v);
						if (makespan <= min_makespan && (makespan < sol_vec[sol_index]->makespan
							/*|| !check_tabu(sol_index, oper_u, oper_v, cur_iter, FORWARD_INSERT)*/))
						{
							if (makespan < min_makespan)
							{
								min_makespan = makespan;
								min_oper_u = oper_u;
								min_oper_v = oper_v;
								equ_cnt = 1;
							}
							else if (makespan == min_makespan)
							{
								equ_cnt += 1;
								if (rand() % equ_cnt)
								{
									min_oper_u = oper_u;
									min_oper_v = oper_v;
								}
							}
						}
						//cout << makespan << "\t";
					}
				}
				//cout << endl;
			}
			//cout << endl;
		}
		if (min_makespan == INT_MAX)	// still_improve==false
			continue;
		cout << "The best move: ";
		cout << cur_iter << "\t"
			<< min_oper_u->job_i << ", " << min_oper_u->oper_job_i << "\t"
			<< min_oper_v->job_i << ", " << min_oper_v->oper_job_i << "\t"
			<< min_makespan << "\t";
		/*update_tabu(sol_index, min_oper_u, min_oper_v, cur_iter,FORWARD_INSERT);*/
		//display_machine_operation(sol_index, min_oper_u->mach_i);
		forward_insert(sol_index, min_oper_u, min_oper_v);
		calculate_qr_obj(sol_index);
		if (alg_best_makespan> sol_vec[sol_index]->makespan)
			alg_best_makespan = sol_vec[sol_index]->makespan;
		//display_machine_operation(sol_index, min_oper_u->mach_i);
		cout << sol_vec[sol_index]->makespan << "\t" << alg_best_makespan << endl;
		//check_solution(sol_index);
	}
}
void Solver::insert_move(int sol_index)
{
	cout << "insert move, solution " << sol_index << ", " << sol_vec[sol_index]->makespan << endl;
	Solution::Operation *min_oper_u = NULL, *min_oper_v = NULL;
	int  min_makespan, equ_cnt, min_oper_mach_pos;
	MOVE_TYPE move_type;
	for (int cur_iter = 1; cur_iter <= iteration; cur_iter++)
	{
		min_makespan = INT_MAX;
		int critical_block_pos = 0;
		/*if (cur_iter == 3527)
		{
		cout << endl;
		display_machine_operation(sol_index, 3);
		}*/
		for (vector<Solution::Operation *>::iterator oper_iter = critical_block_vec.begin(); oper_iter != critical_block_vec.end(); oper_iter += 2, critical_block_pos += 2)
		{
			Solution::Operation *oper_u = *oper_iter;	// move u behind v
			int oper_v_pos = critical_block_mach_pos_vec[critical_block_pos] + 1;
			for (Solution::Operation *oper_v = oper_u->next_mach_oper; oper_v != (*(oper_iter + 1))->next_mach_oper; oper_v = oper_v->next_mach_oper, oper_v_pos += 1)
			{
				if (oper_v->q >= oper_u->next_job_oper->q &&	// q[v]>=q[JS[u]]
					oper_v->job_i != oper_u->job_i&&	// u, v do not belong to the same job
					oper_u->pre_job_oper->q + oper_u->pre_job_oper->r + oper_u->pre_job_oper->t == sol_vec[sol_index]->makespan)
				{
					/*cout << oper_u->job_i << "*, " << oper_u->oper_job_i << "\t"
					<< oper_v->job_i << ", " << oper_v->oper_job_i << "\t";*/
					int makespan;
					try_backward_insert_move(sol_index, makespan, oper_u, oper_v);
					if (makespan <= min_makespan && (makespan<alg_best_makespan ||
						!check_tabu(sol_index, cur_iter, critical_block_mach_pos_vec[critical_block_pos], oper_v_pos, oper_u, oper_v, BACKWARD_INSERT)))
					{
						if (makespan < min_makespan)
						{
							min_makespan = makespan;
							min_oper_u = oper_u;
							min_oper_v = oper_v;
							min_oper_mach_pos = critical_block_mach_pos_vec[critical_block_pos];
							move_type = BACKWARD_INSERT;
							equ_cnt = 1;
						}
						else if (makespan == min_makespan)
						{
							equ_cnt += 1;
							if (rand() % equ_cnt)
							{
								min_oper_u = oper_u;
								min_oper_v = oper_v;
								min_oper_mach_pos = critical_block_mach_pos_vec[critical_block_pos];
								move_type = BACKWARD_INSERT;
							}
						}
					}
					//cout << makespan << "\t";
				}
				// move v before u
				if (oper_u->r + oper_u->t >= oper_v->pre_job_oper->r + oper_v->pre_job_oper->t &&	// r[u]+t[u]>=r[JP[v]]+t[JP[v]]
					oper_v->job_i != oper_u->job_i&&	// u, v do not belong to the same job
					oper_u->pre_job_oper->q + oper_u->pre_job_oper->r + oper_u->pre_job_oper->t == sol_vec[sol_index]->makespan)
				{
					/*cout << oper_u->job_i << "*, " << oper_u->oper_job_i << "\t"
					<< oper_v->job_i << ", " << oper_v->oper_job_i << "\t";*/
					int makespan;
					try_forward_insert_move(sol_index, makespan, oper_u, oper_v);
					if (makespan <= min_makespan && (makespan<alg_best_makespan ||
						!check_tabu(sol_index, cur_iter, critical_block_mach_pos_vec[critical_block_pos],
							oper_v_pos, oper_u, oper_v, FORWARD_INSERT)))
					{
						if (makespan < min_makespan)
						{
							min_makespan = makespan;
							min_oper_u = oper_u;
							min_oper_v = oper_v;
							min_oper_mach_pos = critical_block_mach_pos_vec[critical_block_pos];
							move_type = FORWARD_INSERT;
							equ_cnt = 1;
						}
						else if (makespan == min_makespan)
						{
							equ_cnt += 1;
							if (rand() % equ_cnt)
							{
								min_oper_u = oper_u;
								min_oper_v = oper_v;
								min_oper_mach_pos = critical_block_mach_pos_vec[critical_block_pos];
								move_type = FORWARD_INSERT;
							}
						}
					}
					//cout << makespan << "\t";
				}
				//cout << endl;
			}
			//cout << endl;
			if ((*oper_iter)->next_mach_oper == *(oper_iter + 1))
				continue;
			Solution::Operation *oper_v = *(oper_iter + 1);
			int oper_u_pos = critical_block_mach_pos_vec[critical_block_pos] + 1;
			for (Solution::Operation *oper_u = (*oper_iter)->next_mach_oper; oper_u != oper_v; oper_u = oper_u->next_mach_oper, oper_u_pos += 1)
			{
				if (oper_v->q >= oper_u->next_job_oper->q &&	// q[v]>=q[JS[u]]
					oper_v->job_i != oper_u->job_i&&	// u, v do not belong to the same job
					oper_v->next_job_oper->q + oper_v->next_job_oper->r + oper_v->next_job_oper->t == sol_vec[sol_index]->makespan)
				{
					/*cout << oper_u->job_i << "*, " << oper_u->oper_job_i << "\t"
					<< oper_v->job_i << ", " << oper_v->oper_job_i << "\t";*/
					int makespan;
					try_backward_insert_move(sol_index, makespan, oper_u, oper_v);
					if (makespan <= min_makespan && (makespan<alg_best_makespan ||
						!check_tabu(sol_index, cur_iter, oper_u_pos, critical_block_mach_pos_vec[critical_block_pos + 1], oper_u, oper_v, BACKWARD_INSERT)))
					{
						if (makespan < min_makespan)
						{
							min_makespan = makespan;
							min_oper_u = oper_u;
							min_oper_v = oper_v;
							min_oper_mach_pos = oper_u_pos;
							move_type = BACKWARD_INSERT;
							equ_cnt = 1;
						}
						else if (makespan == min_makespan)
						{
							equ_cnt += 1;
							if (rand() % equ_cnt)
							{
								min_oper_u = oper_u;
								min_oper_v = oper_v;
								min_oper_mach_pos = oper_u_pos;
								move_type = BACKWARD_INSERT;
							}
						}
					}
					//cout << makespan << "\t";
				}
				// move v before u
				if (oper_u->r + oper_u->t >= oper_v->pre_job_oper->r + oper_v->pre_job_oper->t &&	// r[u]+t[u]>=r[JP[v]]+t[JP[v]]
					oper_v->job_i != oper_u->job_i&&	// u, v do not belong to the same job
					oper_v->next_job_oper->q + oper_v->next_job_oper->r + oper_v->next_job_oper->t == sol_vec[sol_index]->makespan)
				{
					/*cout << oper_u->job_i << "*, " << oper_u->oper_job_i << "\t"
					<< oper_v->job_i << ", " << oper_v->oper_job_i << "\t";*/
					int makespan;
					try_forward_insert_move(sol_index, makespan, oper_u, oper_v);
					if (makespan <= min_makespan && (makespan<alg_best_makespan ||
						!check_tabu(sol_index, cur_iter, oper_u_pos, critical_block_mach_pos_vec[critical_block_pos + 1], oper_u, oper_v, FORWARD_INSERT)))
					{
						if (makespan < min_makespan)
						{
							min_makespan = makespan;
							min_oper_u = oper_u;
							min_oper_v = oper_v;
							min_oper_mach_pos = oper_u_pos;
							move_type = FORWARD_INSERT;
							equ_cnt = 1;
						}
						else if (makespan == min_makespan)
						{
							equ_cnt += 1;
							if (rand() % equ_cnt)
							{
								min_oper_u = oper_u;
								min_oper_v = oper_v;
								min_oper_mach_pos = oper_u_pos;
								move_type = FORWARD_INSERT;
							}
						}
					}
					//cout << makespan << "\t";
				}
				//cout << endl;
			}
		}
		if (min_makespan == INT_MAX)	// still_improve==false
		{
			perturb(sol_index, 2);
			continue;
		}
		cout << "The best move: ";
		cout << cur_iter << "\t"
			<< min_oper_u->job_i << ", " << min_oper_u->oper_job_i << "\t"
			<< min_oper_v->job_i << ", " << min_oper_v->oper_job_i << "\t" << min_oper_v->mach_i << "\t"
			<< (move_type == BACKWARD_INSERT ? "b" : "f") << "\t"
			<< min_makespan << "\t";
		update_tabu(sol_index, cur_iter, min_oper_mach_pos, min_oper_u, min_oper_v, move_type);
		if (move_type == BACKWARD_INSERT)
			backward_insert(sol_index, min_oper_u, min_oper_v);
		else
			forward_insert(sol_index, min_oper_u, min_oper_v);
		calculate_qr_obj(sol_index);
		/*if (cur_iter >= 3527 && cur_iter <= 3545)
		display_machine_operation(sol_index, 3);*/
		if (alg_best_makespan > sol_vec[sol_index]->makespan)
		{
			//check_solution(sol_index);
			alg_best_makespan = sol_vec[sol_index]->makespan;
			alg_best_iter = cur_iter;
			cout << cur_iter << "\t"
				<< sol_vec[sol_index]->makespan << "\t" << alg_best_makespan << endl;
		}
		cout << sol_vec[sol_index]->makespan << "\t" << alg_best_makespan << endl;
	}
}
void Solver::insert_move1(int sol_index)
{
	cout << "insert move, solution " << sol_index << ", " << makespan[sol_index] << endl;
	int min_u, min_v,min_mach_i, mkspan, min_makespan, equ_cnt;
	MOVE_TYPE move_type, min_move_type;
	alg_best_makespan = makespan[sol_index];
	start_time = clock();
	for (int cur_iter = 1; cur_iter <= iteration; cur_iter++)
	{
		min_makespan = INT_MAX;
		/*if (cur_iter == 3527)
		{
		cout << endl;
		display_machine_operation(sol_index, 3);
		}*/
		for (int block_i = 1; block_i <= crit_block[0][0]; block_i++)
		{
			int mach_i = crit_block[block_i][0];
			int u = crit_block[block_i][1];
			for (int v = u + 1; v <= crit_block[block_i][2]; v++)
			{
				if (machine[sol_index][mach_i][v]->q + machine[sol_index][mach_i][v]->t>= machine[sol_index][mach_i][u]->next_job_oper->q+ machine[sol_index][mach_i][u]->next_job_oper->t&&	// q[v]>=q[JS[u]]
					machine[sol_index][mach_i][u]->pre_job_oper!=dummy_oper[sol_index][START]&&
					machine[sol_index][mach_i][u]->job_i != machine[sol_index][mach_i][v]->job_i)
				{
					// move u behind v
					move_type = BACKWARD_INSERT;
					try_backward_insert_move1(sol_index, mkspan, mach_i, u, v);
					if (mkspan <= min_makespan && (mkspan < alg_best_makespan || !check_tabu1(sol_index, cur_iter, mach_i, u, v, move_type)))
					{
						if (mkspan < min_makespan)
						{
							min_makespan = mkspan;
							min_u = u;
							min_v = v;
							min_mach_i = mach_i;
							min_move_type = move_type;
							equ_cnt = 1;
						}
						else if (mkspan == min_makespan)
						{
							equ_cnt += 1;
							if (rand() % equ_cnt == 0)
							{
								min_u = u;
								min_v = v;
								min_mach_i = mach_i;
								min_move_type = move_type;
							}
						}
					}
				}
				if (v!=u+1&&machine[sol_index][mach_i][u]->end_time >= machine[sol_index][mach_i][v]->pre_job_oper->end_time&&	// r[u] + t[u] >= r[JP[v]] + t[JP[v]]
					machine[sol_index][mach_i][u]->pre_job_oper != dummy_oper[sol_index][START] &&
					machine[sol_index][mach_i][u]->job_i != machine[sol_index][mach_i][v]->job_i)
				{
					// move v before u
					move_type = FORWARD_INSERT;
					try_forward_insert_move1(sol_index, mkspan, mach_i, u, v);
					if (mkspan <= min_makespan && (mkspan < alg_best_makespan || !check_tabu1(sol_index, cur_iter, mach_i, u, v, move_type)))
					{
						if (mkspan < min_makespan)
						{
							min_makespan = mkspan;
							min_u = u;
							min_v = v;
							min_mach_i = mach_i;
							min_move_type = move_type;
							equ_cnt = 1;
						}
						else if (mkspan == min_makespan)
						{
							equ_cnt += 1;
							if (rand() % equ_cnt == 0)
							{
								min_u = u;
								min_v = v;
								min_mach_i = mach_i;
								min_move_type = move_type;
							}
						}
					}
				}
			}

			int v = crit_block[block_i][2];
			for (u = crit_block[block_i][1] + 1; u < v; u++)
			{
				if (machine[sol_index][mach_i][v]->q + machine[sol_index][mach_i][v]->t>= machine[sol_index][mach_i][u]->next_job_oper->q+ machine[sol_index][mach_i][u]->next_job_oper->t&&	// q[v]>=q[JS[u]]
					machine[sol_index][mach_i][v]->next_job_oper!=dummy_oper[sol_index][END]&&
					machine[sol_index][mach_i][u]->job_i != machine[sol_index][mach_i][v]->job_i)
				{
					// move u behind v
					move_type = BACKWARD_INSERT;
					try_backward_insert_move1(sol_index, mkspan, mach_i, u, v);
					if (mkspan <= min_makespan && (mkspan < alg_best_makespan || !check_tabu1(sol_index, cur_iter, mach_i, u, v, move_type)))
					{
						if (mkspan < min_makespan)
						{
							min_makespan = mkspan;
							min_u = u;
							min_v = v;
							min_mach_i = mach_i;
							min_move_type = move_type;
							equ_cnt = 1;
						}
						else if (mkspan == min_makespan)
						{
							equ_cnt += 1;
							if (rand() % equ_cnt == 0)
							{
								min_u = u;
								min_v = v;
								min_mach_i = mach_i;
								min_move_type = move_type;
							}
						}
					}
				}
				if (v!=u+1&&machine[sol_index][mach_i][u]->end_time >= machine[sol_index][mach_i][v]->pre_job_oper->end_time&&	// r[u] + t[u] >= r[JP[v]] + t[JP[v]]
					machine[sol_index][mach_i][u]->job_i != machine[sol_index][mach_i][v]->job_i)
				{
					// move v before u
					move_type = FORWARD_INSERT;
					try_forward_insert_move1(sol_index, mkspan, mach_i, u, v);
					if (mkspan <= min_makespan && (mkspan < alg_best_makespan || !check_tabu1(sol_index, cur_iter, mach_i, u, v, move_type)))
					{
						if (mkspan < min_makespan)
						{
							min_makespan = mkspan;
							min_u = u;
							min_v = v;
							min_mach_i = mach_i;
							min_move_type = move_type;
							equ_cnt = 1;
						}
						else if (mkspan == min_makespan)
						{
							equ_cnt += 1;
							if (rand() % equ_cnt == 0)
							{
								min_u = u;
								min_v = v;
								min_mach_i = mach_i;
								min_move_type = move_type;
							}
						}
					}
				}
			}

		}
		if (min_makespan == INT_MAX)
		{
			perturb1(sol_index,cur_iter,2);
			continue;
		}
		//check_solution1(sol_index);
		//bool cyc = check_cycle(sol_index);
		update_tabu1(sol_index, cur_iter, min_mach_i, min_u, min_v);
		apply_move(sol_index, min_mach_i, min_u, min_v, min_move_type);
		calculate_q_crit_block(sol_index);
		//cyc = check_cycle(sol_index);
		//check_solution1(sol_index);
		/*if (cur_iter >= 3527 && cur_iter <= 3545)
		display_machine_operation(sol_index, 3);*/
		if (alg_best_makespan > makespan[sol_index])
		{
			end_time = clock();
			check_solution1(sol_index);
			alg_best_makespan = makespan[sol_index];
			alg_best_iter = cur_iter;

			cout << cur_iter << "\t" << min_mach_i << "\t"
				<< min_u << "\t" << min_v << "\t"
				<< (min_move_type == BACKWARD_INSERT ? "b" : "f") << "\t"
				<< min_makespan << "\t"
				<< makespan[sol_index] << "\t" << alg_best_makespan <<"\t"
				<<crit_block[0][0]<<"\t"
				<< (end_time - start_time) / CLOCKS_PER_SEC
				<< endl;
		}
		/*cout << cur_iter << "\t" << min_mach_i << "\t"
			<< min_u << "\t" << min_v << "\t"
			<< (min_move_type == BACKWARD_INSERT ? "b" : "f") << "\t"
			<< min_makespan << "\t"
			<< makespan[sol_index] << "\t" << alg_best_makespan << endl;*/
	}
}
void Solver::perturb(int sol_index, int ptr_len)
{
	Solution::Operation *min_oper_u = NULL, *min_oper_v = NULL;
	int  min_makespan = 0, equ_cnt = 0;
	MOVE_TYPE move_type;
	for (int cur_iter = 1; cur_iter <= ptr_len; cur_iter++)
	{
		equ_cnt = 0;
		for (vector<Solution::Operation *>::iterator oper_iter = critical_block_vec.begin(); oper_iter != critical_block_vec.end(); oper_iter += 2)
		{
			if (*oper_iter == *(oper_iter + 1))
				continue;
			Solution::Operation *oper_u = *oper_iter;	// move u behind v
			for (Solution::Operation *oper_v = oper_u->next_mach_oper; oper_v != (*(oper_iter + 1))->next_mach_oper; oper_v = oper_v->next_mach_oper)
			{
				if (oper_v->q >= oper_u->next_job_oper->q &&	// q[v]>=q[JS[u]]
					oper_v->job_i != oper_u->job_i)	// u, v do not belong to the same job
				{
					equ_cnt += 1;
					if (rand() % equ_cnt == 0)
					{
						min_oper_u = oper_u;
						min_oper_v = oper_v;
						move_type = BACKWARD_INSERT;
					}
				}
				// move v before u
				if (oper_u->r + oper_u->t >= oper_v->pre_job_oper->r + oper_v->pre_job_oper->t &&	// r[u]+t[u]>=r[JP[v]]+t[JP[v]]
					oper_v->job_i != oper_u->job_i)	// u, v do not belong to the same job
				{
					equ_cnt += 1;
					if (rand() % equ_cnt == 0)
					{
						min_oper_u = oper_u;
						min_oper_v = oper_v;
						move_type = FORWARD_INSERT;
					}
				}
			}

			Solution::Operation *oper_v = *(oper_iter + 1);
			for (Solution::Operation *oper_u = (*oper_iter)->next_mach_oper; oper_u != oper_v; oper_u = oper_u->next_mach_oper)
			{
				if (oper_v->q >= oper_u->next_job_oper->q &&	// q[v]>=q[JS[u]]
					oper_v->job_i != oper_u->job_i)	// u, v do not belong to the same job
				{
					equ_cnt += 1;
					if (rand() % equ_cnt == 0)
					{
						min_oper_u = oper_u;
						min_oper_v = oper_v;
						move_type = BACKWARD_INSERT;
					}
				}
				// move v before u
				if (oper_u->r + oper_u->t >= oper_v->pre_job_oper->r + oper_v->pre_job_oper->t &&	// r[u]+t[u]>=r[JP[v]]+t[JP[v]]
					oper_v->job_i != oper_u->job_i)	// u, v do not belong to the same job
				{
					equ_cnt += 1;
					if (rand() % equ_cnt == 0)
					{
						min_oper_u = oper_u;
						min_oper_v = oper_v;
						move_type = FORWARD_INSERT;
					}
				}
			}
		}
		/*cout << "The best move*: ";
		cout << cur_iter << "\t"
		<< min_oper_u->job_i << ", " << min_oper_u->oper_job_i << "\t"
		<< min_oper_v->job_i << ", " << min_oper_v->oper_job_i << "\t"
		<< (move_type == BACKWARD_INSERT ? "1" : "0") << "\t"
		<< min_makespan << "\t";*/
		if (move_type == BACKWARD_INSERT)
			backward_insert(sol_index, min_oper_u, min_oper_v);
		else
			forward_insert(sol_index, min_oper_u, min_oper_v);
		calculate_qr_obj(sol_index);
		if (alg_best_makespan> sol_vec[sol_index]->makespan)
			alg_best_makespan = sol_vec[sol_index]->makespan;
		//cout << sol_vec[sol_index]->makespan << "\t" << alg_best_makespan << endl;
		//check_solution(sol_index);
	}
}
void Solver::perturb1(int sol_index, int cur_iter, int ptr_len)
{
	int  equ_cnt, min_u, min_v, min_mach_i;
	MOVE_TYPE min_move_type;
	for (int ptr_iter = 1; ptr_iter <= ptr_len; ptr_iter++)
	{
		equ_cnt = 0;
		for (int block_i = 1; block_i <= crit_block[0][0]; block_i++)
		{
			int mach_i = crit_block[block_i][0];
			int u = crit_block[block_i][1];
			for (int v = u + 1; v <= crit_block[block_i][2]; v++)
			{
				if(machine[sol_index][mach_i][v]->q + machine[sol_index][mach_i][v]->t >= machine[sol_index][mach_i][u]->next_job_oper->q + machine[sol_index][mach_i][u]->next_job_oper->t&&	// q[v]>=q[JS[u]]
					machine[sol_index][mach_i][u]->job_i != machine[sol_index][mach_i][v]->job_i)
				{
					// move u behind v
					equ_cnt += 1;
					if (rand() % equ_cnt == 0)
					{
						min_u = u;
						min_v = v;
						min_mach_i = mach_i;
						min_move_type = BACKWARD_INSERT;
					}
				}
				if (machine[sol_index][mach_i][u]->end_time >= machine[sol_index][mach_i][v]->pre_job_oper->end_time&&	// r[u] + t[u] >= r[JP[v]] + t[JP[v]]
					machine[sol_index][mach_i][u]->job_i != machine[sol_index][mach_i][v]->job_i)
				{
					// move v before u
					equ_cnt += 1;
					if (rand() % equ_cnt == 0)
					{
						min_u = u;
						min_v = v;
						min_mach_i = mach_i;
						min_move_type = FORWARD_INSERT;
					}
				}
			}

			int v = crit_block[block_i][2];
			for (u = crit_block[block_i][1] + 1; u < v; u++)
			{
				if (machine[sol_index][mach_i][v]->q + machine[sol_index][mach_i][v]->t >= machine[sol_index][mach_i][u]->next_job_oper->q + machine[sol_index][mach_i][u]->next_job_oper->t&&	// q[v]>=q[JS[u]]
					machine[sol_index][mach_i][u]->job_i != machine[sol_index][mach_i][v]->job_i)
				{
					// move u behind v
					equ_cnt += 1;
					if (rand() % equ_cnt == 0)
					{
						min_u = u;
						min_v = v;
						min_mach_i = mach_i;
						min_move_type = BACKWARD_INSERT;
					}
				}
				if (machine[sol_index][mach_i][u]->end_time >= machine[sol_index][mach_i][v]->pre_job_oper->end_time&&	// r[u] + t[u] >= r[JP[v]] + t[JP[v]]
					machine[sol_index][mach_i][u]->job_i != machine[sol_index][mach_i][v]->job_i)
				{
					// move v before u
					equ_cnt += 1;
					if (rand() % equ_cnt == 0)
					{
						min_u = u;
						min_v = v;
						min_mach_i = mach_i;
						min_move_type = FORWARD_INSERT;
					}
				}
			}
		}
		apply_move(sol_index, min_mach_i, min_u, min_v, min_move_type);
		calculate_q_crit_block(sol_index);
		if (alg_best_makespan> makespan[sol_index])
			alg_best_makespan = makespan[sol_index];
		//check_solution1(sol_index);
		
		/*cout << cur_iter << "*\t" << min_mach_i << "\t"
			<< min_u << "\t" << min_v << "\t"
			<< (min_move_type == BACKWARD_INSERT ? "b" : "f") << "\t"
			<< makespan[sol_index] << "\t" << alg_best_makespan << endl;*/
	}
}
void Solver::calculate_qr_obj(int sol_index)
{
	//cout << "calculate q, r, and makespan of solution " << sol_index << endl;
	int num_oper_visited = 0;
	vector<Solution::Operation *> oper_in_degree_vec;	// store operations with in-degree of 0
	stack<Solution::Operation*> oper_re_top_order_stack;	// store operations in stack to get re-topological order
	stack<Solution::Operation*> critical_oper_stack;	// store critical operations
	vector<stack<Solution::Operation*>> critical_oper_stack1(instance->m + 1);	// store critical operations
	for (int mach_i = 1; mach_i < sol_vec[sol_index]->head_oper_mach_vec.size(); mach_i++)
	{
		Solution::Operation *oper = sol_vec[sol_index]->head_oper_mach_vec[mach_i];
		if (oper->oper_job_i == 1)	// or oper->oper_job_i == 1, which means the first operation at this job, in-degree is 0, optimize by removing in_degree
		{
			oper->in_degree = 0;
			oper_in_degree_vec.push_back(oper);
		}
		else
			oper->in_degree = 1;
		oper = oper->next_mach_oper;
		while (oper != sol_vec[sol_index]->end_dummy_oper)
		{
			if (oper->oper_job_i == 1)
				oper->in_degree = 1;
			else
				oper->in_degree = 2;
			oper = oper->next_mach_oper;
		}
	}
	while (!oper_in_degree_vec.empty())
	{
		Solution::Operation *oper = oper_in_degree_vec.front();
		oper->r = MAX(oper->pre_mach_oper->r + oper->pre_mach_oper->t, oper->pre_job_oper->r + oper->pre_job_oper->t);
		oper->start_time = oper->r;
		oper->end_time = oper->start_time + oper->t;
		num_oper_visited += 1;
		oper->next_mach_oper->in_degree -= 1;	// minus the in-degree of next operation at the same machine by 1
		oper->next_job_oper->in_degree -= 1;	// minus the in-degree of next operation of the same job by 1
		if (oper->next_mach_oper->in_degree == 0)
			oper_in_degree_vec.push_back(oper->next_mach_oper);
		if (oper->next_job_oper->in_degree == 0 && oper->next_mach_oper != oper->next_job_oper)	// avoid the next machine operation is the same as the next job operation
			oper_in_degree_vec.push_back(oper->next_job_oper);
		oper_re_top_order_stack.push(oper_in_degree_vec.front());	// push the first element into stack
		oper_in_degree_vec.erase(oper_in_degree_vec.begin());	// delete the first element
	}
	sol_vec[sol_index]->makespan = 0;
	for (int mach_i = 1; mach_i < sol_vec[sol_index]->tail_oper_mach_vec.size(); mach_i++)
	{
		if (sol_vec[sol_index]->makespan < sol_vec[sol_index]->tail_oper_mach_vec[mach_i]->end_time)
			sol_vec[sol_index]->makespan = sol_vec[sol_index]->tail_oper_mach_vec[mach_i]->end_time;
	}
	//cout << endl;
	while (!oper_re_top_order_stack.empty())
	{
		Solution::Operation *oper = oper_re_top_order_stack.top();
		oper->q = MAX(oper->next_mach_oper->q + oper->next_mach_oper->t, oper->next_job_oper->q + oper->next_job_oper->t);
		oper_re_top_order_stack.pop();
		if (oper->r + oper->q + oper->t == sol_vec[sol_index]->makespan)
		{
			critical_oper_stack.push(oper);
			critical_oper_stack1[oper->mach_i].push(oper);
		}
	}
	//cout << "critical blocks: " << endl;
	critical_block_vec.clear();
	critical_block_mach_pos_vec.clear();
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		if (critical_oper_stack1[mach_i].empty())
			continue;
		Solution::Operation *pre_oper, *cur_oper, *first_oper;
		pre_oper = first_oper = critical_oper_stack1[mach_i].top();
		critical_oper_stack1[mach_i].pop();
		/*cout << first_oper->job_i << ", " << first_oper->oper_job_i << "\t" << first_oper->mach_i << "\t"
		<< first_oper->start_time << "\t" << first_oper->end_time << endl;*/
		while (!critical_oper_stack1[mach_i].empty())
		{
			cur_oper = critical_oper_stack1[mach_i].top();
			if (cur_oper->start_time != pre_oper->end_time)
			{
				//cout << endl;
				if (first_oper != pre_oper)
				{
					critical_block_vec.push_back(first_oper);	// optimize by ingoring first_oper == pre_oper to speedup
					critical_block_vec.push_back(pre_oper);

					int first_pos = 0, pre_pos = 0;
					for (Solution::Operation *oper = pre_oper; oper != sol_vec[sol_index]->start_dummy_oper;
						oper = oper->pre_mach_oper)
					{
						pre_pos += 1;
						if (oper == first_oper)
							first_pos = pre_pos;
					}
					critical_block_mach_pos_vec.push_back(pre_pos - first_pos + 1);
					critical_block_mach_pos_vec.push_back(pre_pos);
				}
				first_oper = cur_oper;
			}
			pre_oper = cur_oper;
			/*cout << cur_oper->job_i << ", " << cur_oper->oper_job_i << "\t" << cur_oper->mach_i << "\t"
			<< cur_oper->start_time << "\t" << cur_oper->end_time << endl;*/
			critical_oper_stack1[mach_i].pop();
		}
		if (first_oper != pre_oper)
		{
			critical_block_vec.push_back(first_oper);
			critical_block_vec.push_back(pre_oper);

			int first_pos = 0, pre_pos = 0;
			for (Solution::Operation *oper = pre_oper; oper != sol_vec[sol_index]->start_dummy_oper;
				oper = oper->pre_mach_oper)
			{
				pre_pos += 1;
				if (oper == first_oper)
					first_pos = pre_pos;
			}
			critical_block_mach_pos_vec.push_back(pre_pos - first_pos + 1);
			critical_block_mach_pos_vec.push_back(pre_pos);
		}
		//cout << endl;
	}
	/*cout << "blocks:" << endl;
	for (vector<Solution::Operation*>::iterator oper_iter = critical_block_vec.begin();
	oper_iter != critical_block_vec.end(); oper_iter++)
	{
	Solution::Operation *cur_oper = *oper_iter;
	cout << cur_oper->job_i << ", " << cur_oper->oper_job_i << "\t" << cur_oper->mach_i << "\t"
	<< cur_oper->start_time << "\t" << cur_oper->end_time << endl;
	}*/
	if (num_oper_visited != instance->total_num_operation)
	{
		cout << "ERROR: not all operations are visited" << endl;
		system("pause");
	}
}
// determine critical path, and calculate q for all operations
void Solver::calculate_q_crit_block(int sol_index)
{
	//determine critical path
	int front = 0, rear = last_oper_makespan_machine_num[sol_index];
	memset(critical_flag, 0, sizeof(critical_flag));
	while (front != rear)
	{
		Solution::Operation *cur_oper = last_oper_makespan_machine[sol_index][front++];	
		/*cur_oper->q = MAX(cur_oper->next_job_oper->q + cur_oper->next_job_oper->t,
			machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i + 1]->q + machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i + 1]->t);*/
		if (critical_flag[cur_oper->mach_i][cur_oper->oper_mach_i] == 0)
		{
			critical_flag[cur_oper->mach_i][cur_oper->oper_mach_i] = 1;
			if (cur_oper->pre_job_oper != dummy_oper[sol_index][START] &&
				cur_oper->start_time == cur_oper->pre_job_oper->end_time)
				last_oper_makespan_machine[sol_index][rear++] = cur_oper->pre_job_oper;
			if (machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i - 1] != dummy_oper[sol_index][START] &&
				cur_oper->start_time == machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i - 1]->end_time&&
				machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i - 1] != cur_oper->pre_job_oper)
				// avoid the previous machine operation is the same as the previous job operation
				last_oper_makespan_machine[sol_index][rear++] = machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i - 1];
		}

	}
	int block_cnt = 0;
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		for (int oper_i = 1; oper_i <= machine_oper_num[sol_index][mach_i]; oper_i++)	// according to the second index of critical
		{
			if (critical_flag[mach_i][oper_i] == 1)
			{
				crit_block[++block_cnt][0] = mach_i;
				crit_block[block_cnt][1] = oper_i;
				oper_i += 1;
				while (critical_flag[mach_i][oper_i] == 1 && oper_i <= machine_oper_num[sol_index][mach_i])
					oper_i += 1;
				crit_block[block_cnt][2] = oper_i - 1;
				if (crit_block[block_cnt][1] == crit_block[block_cnt][2])
					block_cnt -= 1;	// a block at least holds two operations
			}
		}
	}
	crit_block[0][0] = block_cnt;

	// calculate q for all operations
	front = 0, rear = 0;
	memset(critical_flag, 0, sizeof(critical_flag));
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
		critical_flag[mach_i][machine_oper_num[sol_index][mach_i]] += 1;
	for (int job_i = 1; job_i <= instance->n; job_i++)
	{
		Solution::Operation *oper = job[sol_index][job_i][job_oper_num[sol_index][job_i]];
		critical_flag[oper->mach_i][oper->oper_mach_i] += 1;
		if (critical_flag[oper->mach_i][oper->oper_mach_i] == 2)
			queue[rear++] = oper;
	}
	while (front != rear)
	{
		Solution::Operation *cur_oper = queue[front++];
		cur_oper->q = MAX(cur_oper->next_job_oper->q + cur_oper->next_job_oper->t,
			machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i + 1]->q + machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i + 1]->t);
		if (cur_oper->pre_job_oper != dummy_oper[sol_index][START])
		{
			critical_flag[cur_oper->pre_job_oper->mach_i][cur_oper->pre_job_oper->oper_mach_i] += 1;
			if (critical_flag[cur_oper->pre_job_oper->mach_i][cur_oper->pre_job_oper->oper_mach_i] == 2)
				queue[rear++] = cur_oper->pre_job_oper;
		}
		if (machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i - 1] != dummy_oper[sol_index][START])
		{
			critical_flag[cur_oper->mach_i][cur_oper->oper_mach_i - 1] += 1;
			if (critical_flag[cur_oper->mach_i][cur_oper->oper_mach_i - 1] == 2)
				queue[rear++] = machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i - 1];
		}
	}
}
int main(int argc, char **argv)
{
	int rs = time(NULL);
	//rs = 1499240941;
	srand(rs);
	char *argv_win[] = { "",	// 0
		"_ifp", "instances\\DemirkolBenchmarksJobShop\\",	//"instances\\Dauzere_Data\\",
		"_sfp","solutions\\best_solutions\\",	// solution file path
		"_ifn", "rcmax_30_15_1",	"_suffix",".txt",	// 01a, .fjs cscmax_20_15_1 rcmax_40_15_5 rcmax_20_15_4
		"_sfn","dmu15_rcmax_30_15_1",	// solution file name dmu45_cscmax_20_15_1 dmu21_rcmax_40_15_5 dmu01_rcmax_20_15_4
		"_sol_num", "6",
		"_tt0","2", "_d1","5", "_d2", "12",
		"_itr","40000000","_best_obj","3343"
	};
	cout << "This is the flexible job shop scheduling problem" << endl;
#ifdef _WIN32
	argc = sizeof(argv_win) / sizeof(argv_win[0]); argv = argv_win;
#endif
	map<string, string> argv_map;
	argv_map["_exe_name"] = argv[0];	// add the exe file name to argv map to append to the output file name
	for (int i = 1; i < argc; i += 2)
		argv_map[string(argv[i])] = string(argv[i + 1]);
	Instance *ins = new Instance(argv_map.at("_ifp") + argv_map.at("_ifn") + argv_map.at("_suffix"));
	//ins->display_instance();
	Solver *solver = new Solver(*ins, stoi(argv_map.at("_sol_num")), 0);
	solver->tt0 = stoi(argv_map.at("_tt0"));
	solver->d1 = stoi(argv_map.at("_d1"));
	solver->d2 = stoi(argv_map.at("_d2"));
	solver->iteration = stoi(argv_map.at("_itr"));
	solver->best_known_makespan = stoi(argv_map.at("_best_obj"));	// read solution
	solver->init_solution1(1);
	solver->calculate_q_crit_block(1);
	//solver->read_solution(2, argv_map.at("_sfp") + argv_map.at("_sfn") + argv_map.at("_suffix"));

	//solver->display_solution(1);
	solver->check_solution1(1);
	//solver->determine_critical_path(1);
	//solver->calculate_qr_obj(2);
	//solver->forward_insert_move(1); 
	//solver->backward_insert_move(1);
	solver->insert_move1(1);
#ifdef _WIN32
	system("pause");
#endif
	return 0;
}
#endif