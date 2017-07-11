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
	//vector<Solution*> sol_vec;
	vector<Solution::Operation*> critical_block_vec;
	vector<int> critical_block_mach_pos_vec;
	const Instance *instance;

	enum DUMMY_OPER { START, END };
	Solution::Operation *job[MAXS][MAXN][MAXO];
	Solution::Operation *machine[MAXS][MAXM][MAXO];
	Solution::Operation *dummy_oper[MAXS][sizeof(DUMMY_OPER)];
	Solution::Operation *crit_path[MAXS][MAXM*MAXO];
	Solution::Operation *queue[MAXN*MAXO];
	int job_oper_num[MAXS][MAXN], machine_oper_num[MAXS][MAXM];
	int makespan[MAXS], crit_path_num[MAXS], queue_num;
	int critical_flag[MAXN][MAXO];
	int crit_block[MAXN*MAXM][3],globel_iteration;

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
	int best_known_makespan, iteration, alg_best_iter;
	//int globel_iter;

	Solver(const Instance &, int);
	void display_machine_operation(int, int)const;
	void read_solution(int, string);
	void init_solution1(int);
	void display_solution(int)const;
	void check_solution1(int);
	void insert_move1(int,int);
	void tabu_search(int);
	void try_backward_insert_move1(int, int&, int, int, int);
	void try_forward_insert_move1(int, int&, int, int, int);
	void apply_move(int, int, int, int, MOVE_TYPE);
	void update_tabu1(int, int, int, int, int);
	bool check_tabu1(int, int, int, int, int, MOVE_TYPE);
	bool check_cycle(int);
	void calculate_q_crit_block(int);
	void perturb1(int, int, int);
	void replace_solution(int, int);
	void clear_tabu_list(int);
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
			//sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_i]->proc_i = 1;	// only one machine to process the operation
			//sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_i]->mach_i = instance->job_vec[job_i]->oper_vec[oper_i]->proc_vec[1]->mach_i;
			//sol_vec[sol_index]->job_vec[job_i]->oper_vec[oper_i]->t = instance->job_vec[job_i]->oper_vec[oper_i]->proc_vec[1]->t;
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
		for (int i = 0; i < fields_vec.size(); i++)
		{
			int job_i = stoi(fields_vec[i]) + 1;
			/*Solution::Operation *oper = sol_vec[sol_index]->job_vec[job_i]->oper_vec[1];
			while (oper->mach_i != mach_i)
				oper = oper->next_job_oper;
			oper->pre_mach_oper = sol_vec[sol_index]->tail_oper_mach_vec[mach_i];
			if (sol_vec[sol_index]->head_oper_mach_vec[mach_i] == NULL)
				sol_vec[sol_index]->head_oper_mach_vec[mach_i] = oper;
			else
				oper->pre_mach_oper->next_mach_oper = oper;
			sol_vec[sol_index]->tail_oper_mach_vec[mach_i] = oper;	*/		
		}
		//sol_vec[sol_index]->tail_oper_mach_vec[mach_i]->next_mach_oper = sol_vec[sol_index]->end_dummy_oper;
		/*if (sol_vec[sol_index]->tail_oper_mach_vec[mach_i]->end_time > makespan)
		makespan = sol_vec[sol_index]->tail_oper_mach_vec[mach_i]->end_time;*/
	}
	//calculate_qr_obj(sol_index);
	/*if (best_known_makespan != sol_vec[sol_index]->makespan)
	{
		cout << "ERROR: the solution given by the author is wrong." << endl;
		system("pause");
	}*/
	ifs.close();
}
void Solver::display_solution(int sol_index)const
{
	cout << "*** display solution " << sol_index << " ***, makespan: " << makespan[sol_index] << endl;		
	for (int job_i = 1; job_i <= instance->n; job_i++)
	{
		cout << job_oper_num[sol_index][job_i] << " | ";
		for (int oper_i = 1; oper_i <= job_oper_num[sol_index][job_i]; oper_i++)
			cout << job[sol_index][job_i][oper_i]->mach_i << "\t" << job[sol_index][job_i][oper_i]->oper_mach_i << "\t"
			<< job[sol_index][job_i][oper_i]->start_time << "\t" << job[sol_index][job_i][oper_i]->end_time << "\t";
		cout << endl;
	}

	// for each machine to display solution
	for (int mach_i = 1; mach_i <=instance->m; mach_i++)
	{
		cout << "m " << mach_i << ", " << machine_oper_num[sol_index][mach_i] << ": ";
		for (int oper_i = 1; oper_i <= machine_oper_num[sol_index][mach_i]; oper_i++)
			cout << machine[sol_index][mach_i][oper_i]->job_i << ", " << machine[sol_index][mach_i][oper_i]->oper_job_i << "\t"
			<< machine[sol_index][mach_i][oper_i]->start_time << "\t" << machine[sol_index][mach_i][oper_i]->end_time << "\t";
		cout << endl;
	}
}
void Solver::display_machine_operation(int sol_index, int mach_index)const
{
	cout << "solution " << sol_index << ", machine " << mach_index << ", operations: ";
	for (int oper_i = 1; oper_i <= machine_oper_num[sol_index][mach_index]; oper_i++)
		cout << machine[sol_index][mach_index][oper_i]->job_i << ", " << machine[sol_index][mach_index][oper_i]->oper_job_i << "\t"
		<< machine[sol_index][mach_index][oper_i]->start_time << "\t" << machine[sol_index][mach_index][oper_i]->end_time << "\t";
	cout << endl;
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
		if (machine[sol_index][mach_i][1]->oper_job_i == 1)	// alternative
		{
			queue[rear++] = machine[sol_index][mach_i][1];
			critical_flag[mach_i][1] = 0;
		}
		else
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
	memset(critical_flag, 0, sizeof(critical_flag));
	while (!top_oper_stack.empty())
	{
		Solution::Operation *oper = top_oper_stack.top();
		top_oper_stack.pop();
		int real_q = oper->t + MAX(oper->next_job_oper->q, machine[sol_index][oper->mach_i][oper->oper_mach_i + 1]->q);
		if (oper->q != real_q)
		{
			cout << "ERROR: q is wrong" << endl;
			system("pause");
		}
		if (oper->start_time+oper->q== makespan[sol_index])
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
					if (crit_block[block_cnt][0] != block_mach_i || crit_block[block_cnt][1] != first_i || crit_block[block_cnt][2] != last_i)
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
			makespan[sol_index] = machine[sol_index][i][machine_oper_num[sol_index][i]]->end_time;
	}
}
void Solver::try_backward_insert_move1(int sol_index, int &makespan, int mach_i, int u, int v)	// insert oper_u behind oper_v
{
	if (v==machine_oper_num[sol_index][mach_i] && machine[sol_index][mach_i][v + 1]->q != 0)
		cout << endl;

	machine[sol_index][mach_i][u + 1]->apx_r = MAX(machine[sol_index][mach_i][u + 1]->pre_job_oper->end_time, machine[sol_index][mach_i][u - 1]->end_time);
	for (int oper_i = u + 2; oper_i <= v; oper_i++)
		machine[sol_index][mach_i][oper_i]->apx_r = MAX(machine[sol_index][mach_i][oper_i]->pre_job_oper->end_time, machine[sol_index][mach_i][oper_i - 1]->apx_r + machine[sol_index][mach_i][oper_i - 1]->t);
	machine[sol_index][mach_i][u]->apx_r = MAX(machine[sol_index][mach_i][u]->pre_job_oper->end_time, machine[sol_index][mach_i][v]->apx_r + machine[sol_index][mach_i][v]->t);
	machine[sol_index][mach_i][u]->apx_q = machine[sol_index][mach_i][u]->t + MAX(machine[sol_index][mach_i][u]->next_job_oper->q, machine[sol_index][mach_i][v + 1]->q);
	machine[sol_index][mach_i][v]->apx_q = machine[sol_index][mach_i][v]->t + MAX(machine[sol_index][mach_i][v]->next_job_oper->q, machine[sol_index][mach_i][u]->apx_q);
	for (int oper_i = v - 1; oper_i > u; oper_i--)
		machine[sol_index][mach_i][oper_i]->apx_q = machine[sol_index][mach_i][oper_i]->t + MAX(machine[sol_index][mach_i][oper_i]->next_job_oper->q, machine[sol_index][mach_i][oper_i + 1]->apx_q);
	makespan = 0;
	for (int oper_i = u; oper_i <= v; oper_i++)
		if (makespan < machine[sol_index][mach_i][oper_i]->apx_r + machine[sol_index][mach_i][oper_i]->apx_q)
			makespan = machine[sol_index][mach_i][oper_i]->apx_r + machine[sol_index][mach_i][oper_i]->apx_q;
}
void Solver::try_forward_insert_move1(int sol_index, int &makespan, int mach_i, int u, int v)	// insert oper_v before oper_u
{
	machine[sol_index][mach_i][v]->apx_r = MAX(machine[sol_index][mach_i][v]->pre_job_oper->end_time, machine[sol_index][mach_i][u - 1]->end_time);
	machine[sol_index][mach_i][u]->apx_r = MAX(machine[sol_index][mach_i][u]->pre_job_oper->end_time, machine[sol_index][mach_i][v]->apx_r + machine[sol_index][mach_i][v]->t);
	for (int oper_i = u + 1; oper_i < v; oper_i++)
		machine[sol_index][mach_i][oper_i]->apx_r = MAX(machine[sol_index][mach_i][oper_i]->pre_job_oper->end_time, machine[sol_index][mach_i][oper_i - 1]->apx_r + machine[sol_index][mach_i][oper_i - 1]->t);
	machine[sol_index][mach_i][v - 1]->apx_q = machine[sol_index][mach_i][v - 1]->t + MAX(machine[sol_index][mach_i][v - 1]->next_job_oper->q, machine[sol_index][mach_i][v + 1]->q);
	for (int oper_i = v - 2; oper_i >= u; oper_i--)
		machine[sol_index][mach_i][oper_i]->apx_q = machine[sol_index][mach_i][oper_i]->t + MAX(machine[sol_index][mach_i][oper_i]->next_job_oper->q, machine[sol_index][mach_i][oper_i + 1]->apx_q);
	machine[sol_index][mach_i][v]->apx_q = machine[sol_index][mach_i][v]->t + MAX(machine[sol_index][mach_i][v]->next_job_oper->q, machine[sol_index][mach_i][u]->apx_q);
	makespan = 0;
	for (int oper_i = u; oper_i <= v; oper_i++)
		if (makespan < machine[sol_index][mach_i][oper_i]->apx_r + machine[sol_index][mach_i][oper_i]->apx_q)
			makespan = machine[sol_index][mach_i][oper_i]->apx_r + machine[sol_index][mach_i][oper_i]->apx_q;
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
			makespan[sol_index] = machine[sol_index][i][machine_oper_num[sol_index][i]]->end_time;
	}
}
void Solver::update_tabu1(int sol_index, int cur_iter, int mach_i, int u, int v)
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
bool Solver::check_tabu1(int sol_index, int cur_iter, int mach_i, int u, int v, MOVE_TYPE move_type)// if u and v in tabu, return true, else return false
{
	/*if (cur_iter == 151 && oper_u_pos==4&&oper_v_pos==5)
	cout << endl;*/
	int r1 = MAX((makespan[sol_index] - best_known_makespan) / d1, d2);
	for (Tabu *tabu = tabu_list[sol_index][mach_i]; tabu != NULL; tabu = tabu->next)
	{
		if (tabu->tabu_iteration <= cur_iter - r1) // not in tabu status
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
void Solver::clear_tabu_list(int sol_index)
{
	Tabu *tabu1, *tabu2;
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		if (tabu_list[sol_index][mach_i] != NULL)
		{
			tabu1 = tabu_list[sol_index][mach_i];
			tabu2 = tabu1->next;
			while (tabu2 != NULL)
			{
				delete tabu1;
				tabu1 = tabu2;
				tabu2 = tabu2->next;
			}
			if (tabu1 != NULL)
				delete tabu1;
			tabu1 = NULL;
			tabu_list[sol_index][mach_i] = NULL;
		}
	}
}
void Solver::insert_move1(int sol_index,int best_sol_index)
{
	cout << "insert move, solution " << sol_index << ", " << makespan[sol_index] << endl;
	int min_u, min_v, min_mach_i, mkspan, min_makespan, equ_cnt;
	int min_tb_u, min_tb_v, min_tb_mach_i, equ_tb_cnt, min_tb_makespan;
	MOVE_TYPE move_type, min_move_type, min_tb_move_type;
	int ts_iter = iteration;
	for (int local_iter = 1; local_iter <= iteration; local_iter++)
	{
		globel_iteration += 1;
		min_makespan = min_tb_makespan = INT_MAX;
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
				if (machine[sol_index][mach_i][v]->q >= machine[sol_index][mach_i][u]->next_job_oper->q &&	// q[v]>=q[JS[u]]
					machine[sol_index][mach_i][u]->pre_job_oper != dummy_oper[sol_index][START] )
				{
					// move u behind v
					move_type = BACKWARD_INSERT;
					try_backward_insert_move1(sol_index, mkspan, mach_i, u, v);
					/*if (mkspan < alg_best_makespan&&mkspan == min_makespan)
						cout << endl;*/
					bool is_tabu = check_tabu1(sol_index, globel_iteration, mach_i, u, v, move_type);
					if (mkspan <= min_makespan && !is_tabu)
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
					else if (mkspan < min_makespan&&mkspan < makespan[best_sol_index])
					{
						min_makespan = mkspan;
						min_u = u;
						min_v = v;
						min_mach_i = mach_i;
						min_move_type = move_type;
						equ_cnt = 1;
					}
					if (mkspan < min_tb_makespan&&is_tabu)
					{
						min_tb_makespan = mkspan;
						min_tb_u = u;
						min_tb_v = v;
						min_tb_mach_i = mach_i;
						min_tb_move_type = move_type;
						equ_tb_cnt = 1;
					}
					else if (mkspan == min_tb_makespan&&is_tabu)
					{
						equ_tb_cnt += 1;
						if (rand() % equ_tb_cnt == 0)
						{
							min_tb_u = u;
							min_tb_v = v;
							min_tb_mach_i = mach_i;
							min_tb_move_type = move_type;
						}
					}
				}
				if (/*v != u + 1 && */machine[sol_index][mach_i][u]->end_time >= machine[sol_index][mach_i][v]->pre_job_oper->end_time&&	// r[u] + t[u] >= r[JP[v]] + t[JP[v]]
					machine[sol_index][mach_i][u]->pre_job_oper != dummy_oper[sol_index][START])
				{
					// move v before u
					move_type = FORWARD_INSERT;
					try_forward_insert_move1(sol_index, mkspan, mach_i, u, v);
					bool is_tabu = check_tabu1(sol_index, globel_iteration, mach_i, u, v, move_type);
					if (mkspan <= min_makespan && !is_tabu)
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
					else if (mkspan < min_makespan&&mkspan < makespan[best_sol_index])
					{
						min_makespan = mkspan;
						min_u = u;
						min_v = v;
						min_mach_i = mach_i;
						min_move_type = move_type;
						equ_cnt = 1;
					}
					if (mkspan < min_tb_makespan&&is_tabu)
					{
						min_tb_makespan = mkspan;
						min_tb_u = u;
						min_tb_v = v;
						min_tb_mach_i = mach_i;
						min_tb_move_type = move_type;
						equ_tb_cnt = 1;
					}
					else if (mkspan == min_tb_makespan&&is_tabu)
					{
						equ_tb_cnt += 1;
						if (rand() % equ_tb_cnt == 0)
						{
							min_tb_u = u;
							min_tb_v = v;
							min_tb_mach_i = mach_i;
							min_tb_move_type = move_type;
						}
					}
				}
			}

			int v = crit_block[block_i][2];
			for (u = crit_block[block_i][1]; u < v; u++)
			{
				if (machine[sol_index][mach_i][v]->q >= machine[sol_index][mach_i][u]->next_job_oper->q &&	// q[v]>=q[JS[u]]
					machine[sol_index][mach_i][v]->next_job_oper != dummy_oper[sol_index][END])
				{
					// move u behind v
					move_type = BACKWARD_INSERT;
					try_backward_insert_move1(sol_index, mkspan, mach_i, u, v);
					bool is_tabu = check_tabu1(sol_index, globel_iteration, mach_i, u, v, move_type);
					if (mkspan <= min_makespan && !is_tabu)
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
					else if (mkspan < min_makespan&&mkspan <  makespan[best_sol_index])
					{
						min_makespan = mkspan;
						min_u = u;
						min_v = v;
						min_mach_i = mach_i;
						min_move_type = move_type;
						equ_cnt = 1;
					}
					if (mkspan < min_tb_makespan&&is_tabu)
					{
						min_tb_makespan = mkspan;
						min_tb_u = u;
						min_tb_v = v;
						min_tb_mach_i = mach_i;
						min_tb_move_type = move_type;
						equ_tb_cnt = 1;
					}
					else if (mkspan == min_tb_makespan&&is_tabu)
					{
						equ_tb_cnt += 1;
						if (rand() % equ_tb_cnt == 0)
						{
							min_tb_u = u;
							min_tb_v = v;
							min_tb_mach_i = mach_i;
							min_tb_move_type = move_type;
						}
					}
				}
				if (/*v != u + 1 &&*/ machine[sol_index][mach_i][u]->end_time >= machine[sol_index][mach_i][v]->pre_job_oper->end_time&&	// r[u] + t[u] >= r[JP[v]] + t[JP[v]]
					machine[sol_index][mach_i][v]->next_job_oper != dummy_oper[sol_index][END])
				{
					// move v before u
					move_type = FORWARD_INSERT;
					try_forward_insert_move1(sol_index, mkspan, mach_i, u, v);
					bool is_tabu = check_tabu1(sol_index, globel_iteration, mach_i, u, v, move_type);
					if (mkspan <= min_makespan && !is_tabu)
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
					else if (mkspan < min_makespan&&mkspan <  makespan[best_sol_index])
					{
						min_makespan = mkspan;
						min_u = u;
						min_v = v;
						min_mach_i = mach_i;
						min_move_type = move_type;
						equ_cnt = 1;
					}
					if (mkspan < min_tb_makespan&&is_tabu)
					{
						min_tb_makespan = mkspan;
						min_tb_u = u;
						min_tb_v = v;
						min_tb_mach_i = mach_i;
						min_tb_move_type = move_type;
						equ_tb_cnt = 1;
					}
					else if (mkspan == min_tb_makespan&&is_tabu)
					{
						equ_tb_cnt += 1;
						if (rand() % equ_tb_cnt == 0)
						{
							min_tb_u = u;
							min_tb_v = v;
							min_tb_mach_i = mach_i;
							min_tb_move_type = move_type;
						}
					}
				}
			}

		}
		if (min_makespan == INT_MAX)
		{
			perturb1(sol_index, best_sol_index, 2);
			continue;
		}
		//check_solution1(sol_index); 
		//bool cyc = check_cycle(sol_index);
		update_tabu1(sol_index, globel_iteration, min_mach_i, min_u, min_v);
		apply_move(sol_index, min_mach_i, min_u, min_v, min_move_type);
		calculate_q_crit_block(sol_index);
		//cyc = check_cycle(sol_index);
		//check_solution1(sol_index);
		/*if (cur_iter >= 3527 && cur_iter <= 3545)
		display_machine_operation(sol_index, 3);*/
		if (makespan[best_sol_index] > makespan[sol_index])
		{
			end_time = clock();
			check_solution1(sol_index);
			replace_solution(best_sol_index, sol_index);
			alg_best_iter = globel_iteration;
			local_iter = 1;

			cout << globel_iteration << "\t" << local_iter << "\t" << min_mach_i << "\t"
				<< min_u << "\t" << min_v << "\t"
				<< (min_move_type == BACKWARD_INSERT ? "b" : "f") << "\t"
				<< min_makespan << "\t"
				<< makespan[sol_index] << "\t" << makespan[best_sol_index] << "\t"
				<< crit_block[0][0] << "\t"
				<< (end_time - start_time) / CLOCKS_PER_SEC
				<< endl;
		}
		/*if (makespan[best_sol_index] <= 3367)
			cout << globel_iteration << "\t" << local_iter << "\t" << min_mach_i << "\t"
			<< min_u << "\t" << min_v << "\t"
			<< (min_move_type == BACKWARD_INSERT ? "b" : "f") << "\t"
			<< min_makespan << "\t" << min_tb_makespan << "\t"
			<< makespan[sol_index] << "\t" << makespan[best_sol_index] << "\t"
			<< crit_block[0][0] << "\t"
			<< endl;*/
	}
}
void Solver::replace_solution(int dest, int src)
{
	for (int job_i = 1; job_i <= instance->n; job_i++)
	{
		for (int oper_i = 1; oper_i <= job_oper_num[src][job_i]; oper_i++)
		{
			job[dest][job_i][oper_i]->start_time = job[src][job_i][oper_i]->start_time;
			job[dest][job_i][oper_i]->end_time = job[src][job_i][oper_i]->end_time;
			job[dest][job_i][oper_i]->q = job[src][job_i][oper_i]->q;
			job[dest][job_i][oper_i]->t = job[src][job_i][oper_i]->t;
			job[dest][job_i][oper_i]->mach_i = job[src][job_i][oper_i]->mach_i;
			job[dest][job_i][oper_i]->oper_mach_i = job[src][job_i][oper_i]->oper_mach_i;
		}
	}
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		for (int oper_i = 1; oper_i <= machine_oper_num[src][mach_i]; oper_i++)
			machine[dest][mach_i][oper_i] =job[dest][machine[src][mach_i][oper_i]->job_i][machine[src][mach_i][oper_i]->oper_job_i];
		machine[dest][mach_i][0] = dummy_oper[dest][START];
		machine[dest][mach_i][machine_oper_num[src][mach_i]+1] = dummy_oper[dest][END];
		machine_oper_num[dest][mach_i]=machine_oper_num[src][mach_i];
	}
	makespan[dest] = makespan[src];
}
void Solver::tabu_search(int sol_index)
{
	int best_sol = 0, cur_sol = 1;
	init_solution1(best_sol);
	calculate_q_crit_block(best_sol);
	globel_iteration = 0;
	
	check_solution1(best_sol);

	start_time = clock();
	for (int cur_iter = 1; cur_iter <= 200; cur_iter++)
	{
		replace_solution(cur_sol, best_sol);
		calculate_q_crit_block(cur_sol);
		insert_move1(cur_sol, best_sol);
		clear_tabu_list(cur_sol);
		cout << cur_iter << " **********************************" << endl;
	}
}
void Solver::perturb1(int sol_index, int best_sol_index, int ptr_len)
{
	int  equ_cnt, min_u, min_v, min_mach_i;
	MOVE_TYPE min_move_type;
	for (int local_iter = 1; local_iter <= ptr_len; local_iter++)
	{
		equ_cnt = 0;
		for (int block_i = 1; block_i <= crit_block[0][0]; block_i++)
		{
			int mach_i = crit_block[block_i][0];
			int u = crit_block[block_i][1];
			for (int v = u + 1; v <= crit_block[block_i][2]; v++)
			{
				if (machine[sol_index][mach_i][v]->q >= machine[sol_index][mach_i][u]->next_job_oper->q &&	// q[v]>=q[JS[u]]
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
				if (machine[sol_index][mach_i][v]->q >= machine[sol_index][mach_i][u]->next_job_oper->q &&	// q[v]>=q[JS[u]]
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
		if (makespan[best_sol_index] > makespan[sol_index])
			replace_solution(best_sol_index, sol_index);
		//check_solution1(sol_index);

		/*cout << globel_iteration << "*\t" << local_iter << "\t" << min_mach_i << "\t"
			<< min_u << "\t" << min_v << "\t"
			<< (min_move_type == BACKWARD_INSERT ? "b" : "f") << "\t"
			<< makespan[sol_index] << "\t" << makespan[best_sol_index] << "\t"
			<< crit_block[0][0] << "\t"
			<< endl;*/
	}
}
// determine critical path, and calculate q for all operations
void Solver::calculate_q_crit_block(int sol_index)
{
	//determine critical path
	int front = 0, rear = 0;
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		if (makespan[sol_index] == machine[sol_index][mach_i][machine_oper_num[sol_index][mach_i]]->end_time)
			queue[rear++] = machine[sol_index][mach_i][machine_oper_num[sol_index][mach_i]];
	}
	memset(critical_flag, 0, sizeof(critical_flag));
	while (front != rear)
	{
		Solution::Operation *cur_oper = queue[front++];
		if (critical_flag[cur_oper->mach_i][cur_oper->oper_mach_i] == 0)
		{
			critical_flag[cur_oper->mach_i][cur_oper->oper_mach_i] = 1;
			if (cur_oper->pre_job_oper != dummy_oper[sol_index][START] &&
				cur_oper->start_time == cur_oper->pre_job_oper->end_time)
				queue[rear++] = cur_oper->pre_job_oper;
			if (machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i - 1] != dummy_oper[sol_index][START] &&
				cur_oper->start_time == machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i - 1]->end_time&&
				machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i - 1] != cur_oper->pre_job_oper)
				// avoid the previous machine operation is the same as the previous job operation
				queue[rear++] = machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i - 1];
		}

	}
	int block_cnt = 0;
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		for (int oper_i = 1; oper_i <= machine_oper_num[sol_index][mach_i]; oper_i++)	// according to the second index of critical
		{
			if (critical_flag[mach_i][oper_i] == 1 && oper_i <= machine_oper_num[sol_index][mach_i])
			{
				crit_block[++block_cnt][0] = mach_i;
				crit_block[block_cnt][1] = oper_i;
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
		cur_oper->q = cur_oper->t + MAX(cur_oper->next_job_oper->q, machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i + 1]->q);
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
	//rs = 1499240941;1499650432
	srand(rs);
	char *argv_win[] = { "",	// 0
		"_ifp", "instances\\DemirkolBenchmarksJobShop\\",	//"instances\\Dauzere_Data\\",
		"_sfp","solutions\\best_solutions\\",	// solution file path
		"_ifn", "rcmax_30_15_1",	"_suffix",".txt",	// 01a, .fjs cscmax_20_15_1 rcmax_40_15_5 rcmax_20_15_4
		"_sfn","dmu15_rcmax_30_15_1",	// solution file name dmu45_cscmax_20_15_1 dmu21_rcmax_40_15_5 dmu01_rcmax_20_15_4
		"_sol_num", "6",
		"_tt0","2", "_d1","5", "_d2", "12",
		"_itr","12500","_best_obj","3343"
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
	Solver *solver = new Solver(*ins, stoi(argv_map.at("_sol_num")));
	solver->tt0 = stoi(argv_map.at("_tt0"));
	solver->d1 = stoi(argv_map.at("_d1"));
	solver->d2 = stoi(argv_map.at("_d2"));
	solver->iteration = stoi(argv_map.at("_itr"));
	solver->best_known_makespan = stoi(argv_map.at("_best_obj"));	// read solution
	
	solver->tabu_search(1);
#ifdef _WIN32
	system("pause");
#endif
	return 0;
}
#endif