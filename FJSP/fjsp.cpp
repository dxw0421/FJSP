#if 1
#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<stack>
#include<map>
#include<algorithm>
#include<cstring>
#include <iomanip>
#include <ctime>
#include <climits>

using namespace std;

#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAXN 52	// max job
#define MAXM 50	// max machine	
#define MAXS 12	// max solution
#define MAXO 120	// max operation
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
// replace all occurance of t in s to w  
void replace_all(std::string & s, std::string const & t, std::string const & w)
{
	string::size_type pos = s.find(t), t_size = t.size(), w_size = w.size();
	while (pos != std::string::npos) { // found   
		s.replace(pos, t_size, w);
		pos = s.find(t, pos + w_size);
	}
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

	int m, n;	// m and n are the number of machines and jobs
	int avg_num_mach_per_oper, total_num_operation;	// average number of machines per operation
	string ifn;
	vector<vector<vector<Process*>>> job_vec;
	Instance(string);
	~Instance();
	void display_instance()const;
};

class Solver
{
public:
	class Operation
	{
	public:
		const int job_i, oper_job_i;	// the job and operation index of this operation
		int mach_i, oper_mach_i, t, proc_i, start_time, end_time, q, apx_r, apx_q;	// optimize by removing proc_i or mach and t
		Operation /**pre_mach_oper, *next_mach_oper,*/ *pre_job_oper, *next_job_oper;
		Operation(int _job_i, int _oper_job_i) :job_i(_job_i), oper_job_i(_oper_job_i),
			mach_i(0), oper_mach_i(0), proc_i(0), start_time(0), end_time(0), q(0), apx_r(0), apx_q(0) {}
		friend ostream& operator<<(ostream& out, const Operation& oper)
		{
			out << oper.job_i << ", " << oper.oper_job_i << ", "
				<< oper.mach_i << ", " << oper.oper_mach_i << ", "
				<< oper.start_time << ", " << oper.end_time;
			return out;
		}
	};
	class Tabu
	{
	public:
		int tabu_iteration;	// the tabu iteration
		int tabu_u, tabu_oper_num;	// the number of tabu operation
		Operation *tabu_oper[MAXO];	// the tabu block of operations from u to v
		Tabu *front, *next;
	};
	int sol_num;
	//vector<Solution*> sol_vec;
	vector<Operation*> critical_block_vec;
	vector<int> critical_block_mach_pos_vec;
	const Instance *instance;
	//int m, n, total_num_oper;

	enum DUMMY_OPER { START, END };
	enum MOVE_TYPE { FORWARD_INSERT, BACKWARD_INSERT };
	Operation *job[MAXS][MAXN][MAXO];
	Operation *machine[MAXS][MAXM][MAXO];
	Operation *dummy_oper[MAXS][sizeof(DUMMY_OPER)];
	Operation *queue[MAXN*MAXO];
	int job_oper_num[MAXS][MAXN], machine_oper_num[MAXS][MAXM];
	int makespan[MAXS];
	int critical_flag[MAXN][MAXO], critical_flag_g[MAXN][MAXO];
	int crit_block[MAXS][MAXN*MAXM][3], global_iteration, permutation_iteration, assignment_iteration;

	clock_t start_time, end_time;
	Tabu *permutation_tabu_list[MAXS][MAXM];
	int assignment_tabu_list[MAXS][MAXN][MAXO][MAXM];	// sol_index, job_i, job_oper_i, assigned mach_i

	int common_seq[2][MAXM][2 * MAXO];	// in path relinking

	int tt0, d1, d2, tt1, d3, d4, d5, d6;	// for tabu tenure
	int ts_iteraion, ts_restart;
	int big_ptr, big_ptr_type, small_ptr, small_ptr_type;
	int pr_a, pr_p, iter_cyc, gen, run_cnt, run_i, max_t;
	map<string, string> argv_map;
	vector<string> best_opt_header;
	map<string, map<string, string>> best_opt_map;
	int imp_lb;	// the impossible low bound to reach

	int ns_p_cnt, ns_a_cnt;
	ofstream ofs;
	string ofn, instance_name;
	Solver(const Instance &, int, string, map<string, string> &, string);
	~Solver();
	void display_machine_operation(int, int)const;
	void read_solution(int, string);
	void init_solution(int);
	void init_solution1(int);
	void display_solution(int)const;
	void check_solution(int);
	void change_permutation(int, int, int&, int&, int&, int&, MOVE_TYPE&);
	void tabu_search(int, int);
	void ITS();
	void H2O();
	void try_backward_insert_move(int, int&, int, int, int);
	void try_forward_insert_move(int, int&, int, int, int);
	void apply_permutation_move(int, int, int, int, MOVE_TYPE);
	void update_permutation_tabu(int, int, int, int);
	bool check_permutation_tabu(int, int, int, int, MOVE_TYPE);
	void calculate_q_crit_block(int);
	void perturb(int, int, int, int);
	void replace_solution(int, int);
	void clear_tabu_list(int);
	void calculate_r(int);
	void change_machine(int, int, int&, int&, int&, int&, int&);
	void apply_assign_move(int, int, int, int, int);
	void calculate_distance(int, int, int&, int&);
	void path_relinking(int, int, int, int);
	void find_common_seq(int, int);
	void save_best(int, int, int);
	void save_solution(int);
};
Instance::Instance(string _ifn) :ifn(_ifn), total_num_operation(0)
{
	ifstream ifs(ifn);
	if (!ifs.is_open())
	{
		cout << ifn << endl;
		perror("ifn"); exit(0);
	}
	string str_line;
	vector<string>fields_vec;
	if (ifn.find(".fjs") != string::npos)	// for flexible job shop instance
	{
		ifs >> n >> m >> avg_num_mach_per_oper;
		//cout << n << "\t" << m << "\t" << avg_num_mach_per_oper << endl;
		getline(ifs, str_line);	// read the '\n'
		job_vec.resize(n + 1, vector<vector<Process*>>(NULL));
		int job_i = 0;
		while (getline(ifs, str_line))
		{
			if (str_line == "")
				continue;
			if (job_i == n)
				break;
			split_generic<string>(fields_vec, str_line, " \t");
			int field_cnt = 0;
			int num_oper = stoi(fields_vec[field_cnt++]);	// number of operation
			job_i += 1;
			job_vec[job_i].resize(num_oper + 1, vector<Process*>(NULL));
			total_num_operation += num_oper;
			for (int oper_i = 1; oper_i <= num_oper; oper_i++)	// for each operation
			{
				int num_mach = stoi(fields_vec[field_cnt++]);	// number of process
				job_vec[job_i][oper_i].resize(num_mach + 1, NULL);
				for (int proc_i = 1; proc_i <= num_mach; proc_i++)
				{
					int mach = stoi(fields_vec[field_cnt]);
					field_cnt += 1;
					int t = stoi(fields_vec[field_cnt]);
					field_cnt += 1;
					job_vec[job_i][oper_i][proc_i] = new Process(t, mach);	// machine index +1
				}
			}
		}
	}
	else // for job shop instance
	{
		ifs >> n >> m;
		cout << n << "\t" << m << "\t" << avg_num_mach_per_oper << endl;
		getline(ifs, str_line);	// read the '\n'
		job_vec.resize(n + 1, vector<vector<Process*>>(NULL));
		int job_i = 0;
		while (getline(ifs, str_line))
		{
			cout << str_line << endl;
			split_generic<string>(fields_vec, str_line, " \t");
			int field_cnt = 0;
			int num_oper = m;
			job_i += 1;
			job_vec[job_i].resize(num_oper + 1, vector<Process*>(NULL));
			total_num_operation += num_oper;
			for (int oper_i = 0; oper_i < num_oper; oper_i++)	// for each operation
			{
				int num_mach = 1;	// stoi(fields_vec[field_cnt++]);	// number of process
				job_vec[job_i][oper_i].resize(num_mach + 1, NULL);
				for (int proc_i = 0; proc_i < num_mach; proc_i++)
				{
					int mach = stoi(fields_vec[field_cnt]) + 1;
					field_cnt += 1;
					int t = stoi(fields_vec[field_cnt]);
					field_cnt += 1;
					job_vec[job_i][oper_i][proc_i] = new Process(t, mach);	// machine index +1
				}
			}
		}
	}
	ifs.close();
	//display_instance();
}
Instance::~Instance()
{
	/*for (int job_i = 1; job_i <= n; job_i++)
	{
	for (int oper_i = 1; oper_i <= job_vec[job_i]->num_oper; oper_i++)
	{
	for (int proc_i = 1; proc_i <= job_vec[job_i]->oper_vec[oper_i]->num_mach; proc_i++)
	{
	delete job_vec[job_i]->oper_vec[oper_i]->proc_vec[proc_i];
	}
	delete job_vec[job_i]->oper_vec[oper_i];
	}
	delete job_vec[job_i];
	}*/
	for (int job_i = 1; job_i <= job_vec.size() - 1; job_i++)
		for (int oper_i = 1; oper_i <= job_vec[job_i].size() - 1; oper_i++)
			for (int proc_i = 1; proc_i <= job_vec[job_i][oper_i].size() - 1; proc_i++)
				delete job_vec[job_i][oper_i][proc_i];
}
void Instance::display_instance()const
{
	cout << "*** display instance: " << ifn << " ***" << endl;
	cout << n << "\t" << m << "\t" << avg_num_mach_per_oper << endl;
	for (int job_i = 1; job_i <= job_vec.size() - 1; job_i++)
	{
		cout << job_vec[job_i].size() - 1 << " | ";
		for (int oper_i = 1; oper_i <= job_vec[job_i].size() - 1; oper_i++)
		{
			cout << job_vec[job_i][oper_i].size() - 1 << " : ";
			for (int proc_i = 1; proc_i <= job_vec[job_i][oper_i].size() - 1; proc_i++)
			{
				cout << job_vec[job_i][oper_i][proc_i]->mach_i << "\t"
					<< job_vec[job_i][oper_i][proc_i]->t << "\t";
			}
		}
	}
}
Solver::Solver(const Instance &_instance, int _sol_num, string _ofn, map<string, string> &_argv_map, string _qt_ss)
	:instance(&_instance), sol_num(_sol_num), ofn(_ofn), argv_map(_argv_map)
{
	for (int i = 0; i <= sol_num; i++)
	{
		dummy_oper[i][START] = new Operation(0, 0);
		dummy_oper[i][END] = new Operation(0, 1000);
		for (int j = 1; j < instance->job_vec.size(); j++)
		{
			job[i][j][0] = dummy_oper[i][START];
			int k;
			for (k = 1; k < instance->job_vec[j].size(); k++)
			{
				job[i][j][k] = new Operation(j, k);
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
	ifstream ifs_best(_qt_ss);
	if (!ifs_best.is_open())
	{
		cout << _qt_ss << endl;
		perror("_best_opt"); exit(0);
	}
	string str_line;
	vector<string>fields_vec;
	getline(ifs_best, str_line);
	split_generic<string>(best_opt_header, str_line, "*\t");
	while (getline(ifs_best, str_line))
	{
		split_generic<string>(fields_vec, str_line, "*\t");
		if (str_line == "" || fields_vec.size() <= 1 || fields_vec[0] == "Instance")
			continue;
		map<string, string> sol_map;
		for (int i = 0; i < fields_vec.size(); i++)
			sol_map[best_opt_header[i]] = fields_vec[i];
		best_opt_map[fields_vec[0]] = sol_map;
	}
	ifs_best.close();
	bool ofp_is_exist = true;
	ifs_best.open(ofn);
	if (!ifs_best.is_open())
		ofp_is_exist = false;
	ifs_best.close();
	instance_name = argv_map.at("_ifn");
	ofs.open(ofn, ios::app | ios::out);
	ofs.setf(ios::fixed, ios::floatfield);
	ofs.precision(6);
	if (!ofp_is_exist) {
		for (int i = 0; i < best_opt_header.size(); i++)
			ofs << best_opt_header[i] << "\t";
		for (map<string, string>::iterator iter = argv_map.begin(); iter != argv_map.end(); iter++)
			ofs << iter->first << "\t";
		ofs << "run_start_time" << endl;
		for (int i = 0; i < best_opt_header.size(); i++)
			ofs << best_opt_map.at(instance_name).at(best_opt_header[i]) << "\t";
		for (map<string, string>::iterator iter = argv_map.begin(); iter != argv_map.end(); iter++)
			ofs << iter->second << "\t";
		auto t = std::time(nullptr);
		auto tm = *std::localtime(&t);
		ofs << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << std::endl;
		ofs << "run_i\tgen\tglobal_iteration\t"
			<< "permutation_iteration\tassignment_iteration\t"
			<< "time\t"
			<< " makespan"
			<< endl;
	}
	imp_lb = stoi(best_opt_map.at(instance_name).at("imp_lb"));
}
Solver::~Solver()
{
	for (int i = 0; i <= sol_num; i++)
	{
		delete dummy_oper[i][START];
		delete dummy_oper[i][END];
		for (int j = 1; j < instance->job_vec.size(); j++)
			for (int k = 1; k < instance->job_vec[j].size(); k++)
				delete job[i][j][k];
		clear_tabu_list(i);
	}
	ofs.close();
}
void Solver::read_solution(int sol_index, string ifn)
{
	ifstream ifs(ifn);
	if (!ifs.is_open())
	{
		cout << ifn << endl;
		perror("ifn"); exit(0);
	}
	for (int job_i = 1; job_i < instance->job_vec.size(); job_i++)
	{
		for (int oper_i = 1; oper_i < instance->job_vec[job_i].size(); oper_i++)
		{
			cout << job[sol_index][job_i][oper_i]->job_i << ", "
				<< job[sol_index][job_i][oper_i]->oper_job_i << " ";
			job[sol_index][job_i][oper_i]->mach_i = instance->job_vec[job_i][oper_i][1]->mach_i;
			job[sol_index][job_i][oper_i]->t = instance->job_vec[job_i][oper_i][1]->t;
		}
		cout << endl;
	}
	string str_line;
	vector<string>fields_vec;
	getline(ifs, str_line);
	//best_known_makespan = stoi(str_line);
	int mach_i = 0;
	while (getline(ifs, str_line))
	{
		cout << str_line << endl;
		split_generic<string>(fields_vec, str_line, " \t");
		mach_i += 1;
		for (int i = 0; i < fields_vec.size(); i++)
		{
			int job_i = stoi(fields_vec[i]) + 1;
			Operation *oper = job[sol_index][job_i][1];
			while (oper->mach_i != mach_i)
				oper = oper->next_job_oper;
			machine[sol_index][mach_i][i + 1] = oper;
			oper->oper_mach_i = i + 1;
		}
		machine[sol_index][mach_i][fields_vec.size() + 1] = dummy_oper[sol_index][END];
		machine_oper_num[sol_index][mach_i] = fields_vec.size();
	}
	calculate_r(sol_index);
	if (imp_lb != makespan[sol_index])
	{
		cout << "ERROR: the solution given by the author is wrong." << endl;
		system("pause");
	}
	ifs.close();
}
void Solver::display_solution(int sol_index)const
{
	cout << "*** display solution " << sol_index << " ***, makespan: " << makespan[sol_index] << endl;
	/*for (int job_i = 1; job_i <= instance->n; job_i++)
	{
	cout << job_oper_num[sol_index][job_i] << " | ";
	for (int oper_i = 1; oper_i <= job_oper_num[sol_index][job_i]; oper_i++)
	cout << job[sol_index][job_i][oper_i]->mach_i << "\t" << job[sol_index][job_i][oper_i]->oper_mach_i << "\t"
	<< job[sol_index][job_i][oper_i]->start_time << "\t" << job[sol_index][job_i][oper_i]->end_time << "\t";
	cout << endl;
	}*/

	// for each machine to display solution
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		cout << "m " << mach_i << ", " << machine_oper_num[sol_index][mach_i] << ":\t";
		for (int oper_i = 1; oper_i <= machine_oper_num[sol_index][mach_i]; oper_i++)
			cout << machine[sol_index][mach_i][oper_i]->job_i << "," << machine[sol_index][mach_i][oper_i]->oper_job_i << "\t"
			<< machine[sol_index][mach_i][oper_i]->start_time << "\t" << machine[sol_index][mach_i][oper_i]->end_time << "\t" << machine[sol_index][mach_i][oper_i]->t << "\t";
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
void Solver::check_solution(int sol_index)
{
	cout << "***check solution " << sol_index << " *** " << makespan[sol_index] << endl;
	stack<Operation*> top_oper_stack;
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
	if (oper_num != instance->total_num_operation)
	{
		cout << "ERROR: the number of operations is wrong" << endl;
		system("pause");
	}
	oper_num = 0;
	while (front != rear)
	{
		Operation *oper = queue[front++];
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
	if (oper_num != instance->total_num_operation)
	{
		cout << "ERROR: the number of operations is wrong" << endl;
		system("pause");
	}
	memset(critical_flag, 0, sizeof(critical_flag));
	while (!top_oper_stack.empty())
	{
		Operation *oper = top_oper_stack.top();
		top_oper_stack.pop();
		int real_q = MAX(oper->next_job_oper->q + oper->next_job_oper->t,
			machine[sol_index][oper->mach_i][oper->oper_mach_i + 1]->q + machine[sol_index][oper->mach_i][oper->oper_mach_i + 1]->t);
		if (oper->q != real_q)
		{
			cout << "ERROR: q is wrong" << endl;
			system("pause");
		}
		if (oper->end_time + oper->q == makespan[sol_index])
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
					;// block_cnt -= 1;
				else
				{
					if (crit_block[sol_index][block_cnt][0] != block_mach_i || crit_block[sol_index][block_cnt][1] != first_i || crit_block[sol_index][block_cnt][2] != last_i)
					{
						cout << "ERROR: critical block is wrong" << endl;
						system("pause");
					}
				}
			}
		}
	}
	if (crit_block[sol_index][0][0] != block_cnt)
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
void Solver::init_solution(int sol_index)
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
			Instance::Process *min_proc = NULL;
			Operation * cur_oper = job[sol_index][job_i][oper_to_assign_job[job_i]];
			int min_ct = INT_MAX, min_ct_mach, select_proc_i, equ_cnt = 1;
			for (int proc_i = 1; proc_i < instance->job_vec[job_i][oper_to_assign_job[job_i]].size(); proc_i++)	// find the best machine which has the minimum completion time
			{
				Instance::Process *proc = instance->job_vec[job_i][oper_to_assign_job[job_i]][proc_i];
				int oper_i = machine_oper_num[sol_index][proc->mach_i];
				int cur_mach_ct = proc->t + MAX(machine[sol_index][proc->mach_i][oper_i]->end_time,	// the previous operation with the same machine
					job[sol_index][job_i][oper_to_assign_job[job_i] - 1]->end_time);							// the previous operation with the same job
				if (cur_mach_ct < min_ct)
				{
					min_ct = cur_mach_ct;
					select_proc_i = proc_i;
					min_ct_mach = proc->mach_i;
					min_proc = proc;
				}
				else if (cur_mach_ct == min_ct)
				{
					equ_cnt += 1;
					if (rand() % equ_cnt)
					{
						select_proc_i = proc_i;
						min_ct_mach = proc->mach_i;
						min_proc = proc;
					}
				}
			}
			cur_oper->proc_i = select_proc_i;
			cur_oper->mach_i = min_ct_mach;
			cur_oper->t = min_proc->t;
			cur_oper->start_time = MAX(machine[sol_index][min_ct_mach][machine_oper_num[sol_index][min_ct_mach]]->end_time,
				job[sol_index][job_i][oper_to_assign_job[job_i] - 1]->end_time);
			cur_oper->end_time = min_ct;

			machine_oper_num[sol_index][min_ct_mach] += 1;
			cur_oper->oper_mach_i = machine_oper_num[sol_index][min_ct_mach];	// the position at the assigned machine of the operation
			machine[sol_index][min_ct_mach][machine_oper_num[sol_index][min_ct_mach]] = cur_oper;

			oper_to_assign_job[job_i] += 1;	// the next operation is to be assigned
			if (oper_to_assign_job[job_i] >= instance->job_vec[job_i].size())
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
void Solver::init_solution1(int sol_index)
{
	int job_line[MAXN];
	memset(job_line, 0, sizeof(job_line));
	memset(machine_oper_num[sol_index], 0, sizeof(machine_oper_num[sol_index]));
	for (int arrange_cnt = 1; arrange_cnt <= instance->total_num_operation; arrange_cnt++)
	{
		int rand_job = rand() % instance->n + 1;
		while (job_line[rand_job] + 1 > instance->job_vec[rand_job].size() - 1)
			rand_job = rand() % instance->n + 1;
		job_line[rand_job] += 1;
		int rand_proc = rand() % (instance->job_vec[rand_job][job_line[rand_job]].size() - 1) + 1;
		int rand_machine = instance->job_vec[rand_job][job_line[rand_job]][rand_proc]->mach_i;
		machine_oper_num[sol_index][rand_machine] += 1;
		machine[sol_index][rand_machine][machine_oper_num[sol_index][rand_machine]] = job[sol_index][rand_job][job_line[rand_job]];
		machine[sol_index][rand_machine][machine_oper_num[sol_index][rand_machine]]->mach_i = rand_machine;
		machine[sol_index][rand_machine][machine_oper_num[sol_index][rand_machine]]->oper_mach_i = machine_oper_num[sol_index][rand_machine];
		machine[sol_index][rand_machine][machine_oper_num[sol_index][rand_machine]]->t = instance->job_vec[rand_job][job_line[rand_job]][rand_proc]->t;
	}
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		machine[sol_index][mach_i][0] = dummy_oper[sol_index][START];
		machine[sol_index][mach_i][machine_oper_num[sol_index][mach_i] + 1] = dummy_oper[sol_index][END];

	}
	calculate_r(sol_index);
	calculate_q_crit_block(sol_index);
}
// u is less than v, insert u after v, return makespan
void Solver::try_backward_insert_move(int sol_index, int &makespan, int mach_i, int u, int v)
{
	machine[sol_index][mach_i][u + 1]->apx_r = MAX(machine[sol_index][mach_i][u + 1]->pre_job_oper->end_time, machine[sol_index][mach_i][u - 1]->end_time);
	for (int oper_i = u + 2; oper_i <= v; oper_i++)
		machine[sol_index][mach_i][oper_i]->apx_r = MAX(machine[sol_index][mach_i][oper_i]->pre_job_oper->end_time, machine[sol_index][mach_i][oper_i - 1]->apx_r + machine[sol_index][mach_i][oper_i - 1]->t);
	machine[sol_index][mach_i][u]->apx_r = MAX(machine[sol_index][mach_i][u]->pre_job_oper->end_time, machine[sol_index][mach_i][v]->apx_r + machine[sol_index][mach_i][v]->t);
	machine[sol_index][mach_i][u]->apx_q = MAX(machine[sol_index][mach_i][u]->next_job_oper->q + machine[sol_index][mach_i][u]->next_job_oper->t, machine[sol_index][mach_i][v + 1]->q + machine[sol_index][mach_i][v + 1]->t);
	machine[sol_index][mach_i][v]->apx_q = MAX(machine[sol_index][mach_i][v]->next_job_oper->q + machine[sol_index][mach_i][v]->next_job_oper->t, machine[sol_index][mach_i][u]->apx_q + machine[sol_index][mach_i][u]->t);
	for (int oper_i = v - 1; oper_i > u; oper_i--)
		machine[sol_index][mach_i][oper_i]->apx_q = MAX(machine[sol_index][mach_i][oper_i]->next_job_oper->q + machine[sol_index][mach_i][oper_i]->next_job_oper->t, machine[sol_index][mach_i][oper_i + 1]->apx_q + machine[sol_index][mach_i][oper_i + 1]->t);
	makespan = 0;
	for (int oper_i = u; oper_i <= v; oper_i++)
		if (makespan < machine[sol_index][mach_i][oper_i]->apx_r + machine[sol_index][mach_i][oper_i]->apx_q + machine[sol_index][mach_i][oper_i]->t)
			makespan = machine[sol_index][mach_i][oper_i]->apx_r + machine[sol_index][mach_i][oper_i]->apx_q + machine[sol_index][mach_i][oper_i]->t;
}
// u is less than v, insert v before u, return makespan
void Solver::try_forward_insert_move(int sol_index, int &makespan, int mach_i, int u, int v)
{
	machine[sol_index][mach_i][v]->apx_r = MAX(machine[sol_index][mach_i][v]->pre_job_oper->end_time, machine[sol_index][mach_i][u - 1]->end_time);
	machine[sol_index][mach_i][u]->apx_r = MAX(machine[sol_index][mach_i][u]->pre_job_oper->end_time, machine[sol_index][mach_i][v]->apx_r + machine[sol_index][mach_i][v]->t);
	for (int oper_i = u + 1; oper_i < v; oper_i++)
		machine[sol_index][mach_i][oper_i]->apx_r = MAX(machine[sol_index][mach_i][oper_i]->pre_job_oper->end_time, machine[sol_index][mach_i][oper_i - 1]->apx_r + machine[sol_index][mach_i][oper_i - 1]->t);
	machine[sol_index][mach_i][v - 1]->apx_q = MAX(machine[sol_index][mach_i][v - 1]->next_job_oper->q + machine[sol_index][mach_i][v - 1]->next_job_oper->t, machine[sol_index][mach_i][v + 1]->q + machine[sol_index][mach_i][v + 1]->t);
	for (int oper_i = v - 2; oper_i >= u; oper_i--)
		machine[sol_index][mach_i][oper_i]->apx_q = MAX(machine[sol_index][mach_i][oper_i]->next_job_oper->q + machine[sol_index][mach_i][oper_i]->next_job_oper->t, machine[sol_index][mach_i][oper_i + 1]->apx_q + machine[sol_index][mach_i][oper_i + 1]->t);
	machine[sol_index][mach_i][v]->apx_q = MAX(machine[sol_index][mach_i][v]->next_job_oper->q + machine[sol_index][mach_i][v]->next_job_oper->t, machine[sol_index][mach_i][u]->apx_q + machine[sol_index][mach_i][u]->t);
	makespan = 0;
	for (int oper_i = u; oper_i <= v; oper_i++)
		if (makespan < machine[sol_index][mach_i][oper_i]->apx_r + machine[sol_index][mach_i][oper_i]->apx_q + machine[sol_index][mach_i][oper_i]->t)
			makespan = machine[sol_index][mach_i][oper_i]->apx_r + machine[sol_index][mach_i][oper_i]->apx_q + machine[sol_index][mach_i][oper_i]->t;
}
// u is less than v, move u after v, or move v before u
void Solver::apply_permutation_move(int sol_index, int mach_i, int u, int v, MOVE_TYPE move_type)
{
	if (move_type == BACKWARD_INSERT)
	{
		Operation *oper_u = machine[sol_index][mach_i][u];
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
		Operation *oper_v = machine[sol_index][mach_i][v];
		for (int i = v; i > u; i--)
		{
			machine[sol_index][mach_i][i] = machine[sol_index][mach_i][i - 1];
			machine[sol_index][mach_i][i]->oper_mach_i = i;
		}
		machine[sol_index][mach_i][u] = oper_v;
		machine[sol_index][mach_i][u]->oper_mach_i = u;
	}

	//int front = 0, rear = 0;
	//memset(critical_flag, 0, sizeof(critical_flag));
	//queue[rear++] = machine[sol_index][mach_i][u];
	//critical_flag[queue[0]->mach_i][queue[0]->oper_mach_i] = 0;	// in-degree is 0
	//while (front != rear)
	//{
	//	Operation *oper = queue[front++];
	//	/*if (oper->mach_i == 2 && oper->oper_mach_i == 23)
	//		cout << endl;
	//	if (oper->mach_i == 3 && oper->oper_mach_i == 21)
	//		cout << endl;*/
	//	if (oper->next_job_oper != dummy_oper[sol_index][END])
	//	{
	//		if (critical_flag[oper->next_job_oper->mach_i][oper->next_job_oper->oper_mach_i] == 0)
	//		{
	//			critical_flag[oper->next_job_oper->mach_i][oper->next_job_oper->oper_mach_i] = 1;
	//			queue[rear++] = oper->next_job_oper;
	//		}
	//		else if (critical_flag[oper->next_job_oper->mach_i][oper->next_job_oper->oper_mach_i] == 1)
	//			critical_flag[oper->next_job_oper->mach_i][oper->next_job_oper->oper_mach_i] = 2;
	//	}
	//	if (/*machine[sol_index][oper->mach_i][oper->oper_mach_i + 1] != oper->next_job_oper&&*/
	//		machine[sol_index][oper->mach_i][oper->oper_mach_i + 1] != dummy_oper[sol_index][END])
	//	{
	//		if (critical_flag[oper->mach_i][oper->oper_mach_i + 1] == 0)
	//		{
	//			critical_flag[oper->mach_i][oper->oper_mach_i + 1] = 1;
	//			queue[rear++] = machine[sol_index][oper->mach_i][oper->oper_mach_i + 1];
	//		}
	//		else if (critical_flag[oper->mach_i][oper->oper_mach_i + 1] == 1)
	//			critical_flag[oper->mach_i][oper->oper_mach_i + 1] = 2;
	//	}
	//}
	//front = 0, rear = 0;
	//queue[rear++] = machine[sol_index][mach_i][u];


	//while (front != rear)
	//{
	//	Operation *oper = queue[front++];
	//	/*if (oper->mach_i == 2 && oper->oper_mach_i == 22)
	//	cout << endl;
	//	if (oper->mach_i == 3 && oper->oper_mach_i == 21)
	//	cout << endl;*/
	//	oper->start_time = MAX(oper->pre_job_oper->end_time, machine[sol_index][oper->mach_i][oper->oper_mach_i - 1]->end_time);
	//	oper->end_time = oper->start_time + oper->t;
	//	critical_flag[oper->next_job_oper->mach_i][oper->next_job_oper->oper_mach_i] -= 1;
	//	critical_flag[oper->mach_i][oper->oper_mach_i + 1] -= 1;
	//	if (critical_flag[oper->next_job_oper->mach_i][oper->next_job_oper->oper_mach_i] == 0)
	//		queue[rear++] = oper->next_job_oper;
	//	if (critical_flag[oper->mach_i][oper->oper_mach_i + 1] == 0 &&
	//		machine[sol_index][oper->mach_i][oper->oper_mach_i + 1] != oper->next_job_oper)
	//		queue[rear++] = machine[sol_index][oper->mach_i][oper->oper_mach_i + 1];
	//}
	//makespan[sol_index] = 0;
	//for (int i = 1; i <= instance->m; i++)
	//{
	//	if (makespan[sol_index] < machine[sol_index][i][machine_oper_num[sol_index][i]]->end_time)
	//		makespan[sol_index] = machine[sol_index][i][machine_oper_num[sol_index][i]]->end_time;
	//}
}
void Solver::update_permutation_tabu(int sol_index, int mach_i, int u, int v)
{
	Tabu *tabu = new Tabu();
	int r1 = MAX((makespan[sol_index] - imp_lb) / d1, d2);	// optimize by making it global
	tabu->tabu_iteration = permutation_iteration + tt0 + rand() % r1;
	tabu->tabu_u = u;
	tabu->tabu_oper_num = v - u + 1;
	for (int i = 1; i <= tabu->tabu_oper_num; i++)
		tabu->tabu_oper[i] = machine[sol_index][mach_i][u + i - 1];
	tabu->front = NULL;
	tabu->next = permutation_tabu_list[sol_index][mach_i];
	permutation_tabu_list[sol_index][mach_i] = tabu;
	if (tabu->next != NULL)
		tabu->next->front = tabu;
}
bool Solver::check_permutation_tabu(int sol_index, int mach_i, int u, int v, MOVE_TYPE move_type)// if u and v in tabu, return true, else return false
{
	int r1 = MAX((makespan[sol_index] - imp_lb) / d1, d2);
	for (Tabu *tabu = permutation_tabu_list[sol_index][mach_i]; tabu != NULL; tabu = tabu->next)
	{
		if (tabu->tabu_iteration <= permutation_iteration - r1) // not in tabu status
		{
			if (tabu == permutation_tabu_list[sol_index][mach_i])	// tabu is the first node
				permutation_tabu_list[sol_index][mach_i] = NULL;
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
		else if (tabu->tabu_iteration < permutation_iteration)
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
				Operation *oper;
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
				Operation *oper;
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
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		Tabu *tabu1 = permutation_tabu_list[sol_index][mach_i];
		while (tabu1 != NULL)
		{
			Tabu *tabu2 = tabu1->next;
			delete tabu1;
			tabu1 = tabu2;
		}
		permutation_tabu_list[sol_index][mach_i] = NULL;
	}
}
void Solver::change_permutation(int sol_cur, int sol_best, int &min_u, int &min_v, int &min_mach_i, int &min_makespan, MOVE_TYPE &min_move_type)
{
	int mkspan, equ_cnt, equ_g_cnt;
	MOVE_TYPE move_type;
	ns_p_cnt = 0;
	min_makespan = INT_MAX;
	for (int block_i = 1; block_i <= crit_block[sol_cur][0][0]; block_i++)
	{
		int mach_i = crit_block[sol_cur][block_i][0];
		int u = crit_block[sol_cur][block_i][1];
		for (int v = u + 1; v <= crit_block[sol_cur][block_i][2]; v++)
		{
			if (machine[sol_cur][mach_i][u]->job_i == machine[sol_cur][mach_i][v]->job_i)
				continue;
			if (machine[sol_cur][mach_i][v]->start_time < machine[sol_cur][mach_i][u]->next_job_oper->end_time /*&&*/	// q[v] + t[v] > q[JS[u]] || r[JS[u]]+t[JS[u]]>r[v]
				/*machine[sol_cur][mach_i][u]->pre_job_oper != dummy_oper[sol_cur][START]*/)
			{
				// move u behind v
				ns_p_cnt += 1;
				move_type = BACKWARD_INSERT;
				try_backward_insert_move(sol_cur, mkspan, mach_i, u, v);
				if (mkspan <= min_makespan && (mkspan < makespan[sol_best] || !check_permutation_tabu(sol_cur, mach_i, u, v, move_type)))
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
			if (v != u + 1 && machine[sol_cur][mach_i][u]->end_time >
				machine[sol_cur][mach_i][v]->pre_job_oper->start_time/*&&*/	// r[u] + t[u] > r[JP[v]]
				/*machine[sol_cur][mach_i][u]->pre_job_oper != dummy_oper[sol_cur][START]*/)
			{
				// move v before u
				ns_p_cnt += 1;
				move_type = FORWARD_INSERT;
				try_forward_insert_move(sol_cur, mkspan, mach_i, u, v);
				if (mkspan <= min_makespan && (mkspan < makespan[sol_best] || !check_permutation_tabu(sol_cur, mach_i, u, v, move_type)))
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

		int v = crit_block[sol_cur][block_i][2];
		for (u = crit_block[sol_cur][block_i][1]; u < v; u++)
		{
			if (machine[sol_cur][mach_i][u]->job_i == machine[sol_cur][mach_i][v]->job_i)
				continue;
			if (machine[sol_cur][mach_i][v]->start_time < machine[sol_cur][mach_i][u]->next_job_oper->end_time /*&&*/	// q[v] + t[v] > q[JS[u]]
				/*machine[sol_cur][mach_i][v]->next_job_oper != dummy_oper[sol_cur][END]*/)
			{
				// move u behind v
				ns_p_cnt += 1;
				move_type = BACKWARD_INSERT;
				try_backward_insert_move(sol_cur, mkspan, mach_i, u, v);
				if (mkspan <= min_makespan && (mkspan < makespan[sol_best] || !check_permutation_tabu(sol_cur, mach_i, u, v, move_type)))
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
			if (v != u + 1 && machine[sol_cur][mach_i][u]->end_time >
				machine[sol_cur][mach_i][v]->pre_job_oper->start_time/*&&*/	// r[u] + t[u] > r[JP[v]]
				/*machine[sol_cur][mach_i][v]->next_job_oper != dummy_oper[sol_cur][END]*/)
			{
				// move v before u
				ns_p_cnt += 1;
				move_type = FORWARD_INSERT;
				try_forward_insert_move(sol_cur, mkspan, mach_i, u, v);
				if (mkspan <= min_makespan && (mkspan < makespan[sol_best] || !check_permutation_tabu(sol_cur, mach_i, u, v, move_type)))
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
}
void Solver::change_machine(int sol_cur, int sol_best, int &min_u, int &min_mach_u, int &min_v, int &min_mach_v, int &min_makespan)
{
	int equ_cnt, equ_g_cnt;
	min_makespan = INT_MAX;
	ns_a_cnt = 0;
	for (int block_i = 1; block_i <= crit_block[sol_cur][0][0]; block_i++)	// for each block
	{
		int mach_u = crit_block[sol_cur][block_i][0];
		for (int u = crit_block[sol_cur][block_i][1]; u <= crit_block[sol_cur][block_i][2]; u++)	// for each critical operation
		{
			int job_i = machine[sol_cur][mach_u][u]->job_i;
			int oper_job_i = machine[sol_cur][mach_u][u]->oper_job_i;
			if (instance->job_vec[job_i][oper_job_i].size() == 2)		// only one machine
				continue;
			/*if (global_iteration == 1364828)
			cout << "job_i: " << job_i << ", oper_job_i: " << oper_job_i
			<< " mach_i: " << mach_u << ", oper_mach_i: " << u << " | ";*/

			for (int proc_i = 1; proc_i < instance->job_vec[job_i][oper_job_i].size(); proc_i++)	// for each process
			{
				if (instance->job_vec[job_i][oper_job_i][proc_i]->mach_i == mach_u)	// the same machine
					continue;
				int mach_v = instance->job_vec[job_i][oper_job_i][proc_i]->mach_i;
				int left_i = 0, right_i = machine_oper_num[sol_cur][mach_v] + 1;
				Operation *oper_u = machine[sol_cur][mach_u][u];

				for (int m_i = machine_oper_num[sol_cur][mach_v]; m_i >= 1; m_i--)	//determine right_i
				{
					Operation *oper_v = machine[sol_cur][mach_v][m_i];
					if (oper_v->end_time > oper_u->pre_job_oper->end_time)	// move u before v
																			//if (v->end_time <= oper_u->pre_job_oper->start_time)
					{
						right_i = m_i;
						//break;
					}
					else break;
				}

				for (int m_i = 1; m_i <= machine_oper_num[sol_cur][mach_v]; m_i++)	// determine left_i
				{
					Operation *oper_v = machine[sol_cur][mach_v][m_i];	// move u after v
					if (oper_v->q + oper_v->t > oper_u->next_job_oper->q + oper_u->next_job_oper->t)
						//if (oper_u->next_job_oper->q >= v->q + v->t)
					{
						left_i = m_i;
						//break;
					}
					else break;
				}
				if (left_i > right_i)
				{
					int temp = left_i;
					left_i = right_i;
					right_i = temp;
				}
				else if (left_i == right_i)
				{
					left_i -= 1;
					right_i += 1;
				}
				int v_t = instance->job_vec[job_i][oper_job_i][proc_i]->t;

				for (int v = left_i; v < right_i; v++)	// for each feasible position
				{
					Operation *oper_v = machine[sol_cur][mach_v][v];	// insert oper_u after v
					oper_u->apx_r = MAX(oper_u->pre_job_oper->end_time, oper_v->end_time);
					oper_u->apx_q = MAX(oper_u->next_job_oper->q + oper_u->next_job_oper->t,
						machine[sol_cur][mach_v][v + 1]->q + machine[sol_cur][mach_v][v + 1]->t);	// v->next_mach_oper
					ns_a_cnt += 1;
					int mkspan = oper_u->apx_r + oper_u->apx_q + v_t;
					if (mkspan <= min_makespan && (mkspan < makespan[sol_best] || assignment_iteration >= assignment_tabu_list[sol_cur][oper_u->job_i][oper_u->oper_job_i][mach_v]))
					{
						if (mkspan < min_makespan)
						{
							min_makespan = mkspan;
							min_u = u;
							min_mach_u = mach_u;
							min_v = v;
							min_mach_v = mach_v;
							equ_cnt = 1;
						}
						else if (min_makespan == mkspan)
						{
							equ_cnt += 1;
							if (rand() % equ_cnt == 0)
							{
								min_u = u;
								min_mach_u = mach_u;
								min_v = v;
								min_mach_v = mach_v;
							}
						}
					}
					/*if (global_iteration == 1364828)
					cout << v << ", " << mkspan << "\t";*/
				}
			}
			/*if (global_iteration == 1364828)
			cout << endl;*/
		}
	}
	/*if (makespan[sol_best] <= 2515)
	cout << ns_cnt << endl;*/
}
// move u of mach_u after v of mach_v
void Solver::apply_assign_move(int sol_index, int mach_u, int u, int mach_v, int v)
{
	for (int i = machine_oper_num[sol_index][mach_v] + 1; i > v; i--)
	{
		machine[sol_index][mach_v][i + 1] = machine[sol_index][mach_v][i];
		machine[sol_index][mach_v][i + 1]->oper_mach_i = i + 1;
	}
	machine[sol_index][mach_v][v + 1] = machine[sol_index][mach_u][u];
	machine[sol_index][mach_v][v + 1]->mach_i = mach_v;
	machine[sol_index][mach_v][v + 1]->oper_mach_i = v + 1;
	machine_oper_num[sol_index][mach_v] += 1;

	for (int i = u; i <= machine_oper_num[sol_index][mach_u]; i++)
	{
		machine[sol_index][mach_u][i] = machine[sol_index][mach_u][i + 1];
		machine[sol_index][mach_u][i]->oper_mach_i = i;
	}
	machine_oper_num[sol_index][mach_u] -= 1;
}
void Solver::replace_solution(int dest, int src)
{
	for (int job_i = 1; job_i <= instance->n; job_i++)
	{
		for (int oper_i = 1; oper_i <= job_oper_num[src][job_i]; oper_i++)	// copy construct?
		{
			job[dest][job_i][oper_i]->start_time = job[src][job_i][oper_i]->start_time;	// optimize by class copy construction
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
			machine[dest][mach_i][oper_i] = job[dest][machine[src][mach_i][oper_i]->job_i][machine[src][mach_i][oper_i]->oper_job_i];
		machine[dest][mach_i][0] = dummy_oper[dest][START];
		machine[dest][mach_i][machine_oper_num[src][mach_i] + 1] = dummy_oper[dest][END];
		machine_oper_num[dest][mach_i] = machine_oper_num[src][mach_i];
	}
	for (int block_i = 0; block_i <= crit_block[src][0][0]; block_i++)
		memcpy(crit_block[dest][block_i], crit_block[src][block_i], 3 * sizeof(int));
	makespan[dest] = makespan[src];
}
void Solver::tabu_search(int sol_cur, int sol_best)
{
	//cout << "TS: " << sol_cur << "\t" << makespan[sol_cur] << "\t" << makespan[sol_best] << endl;
	int sol_base = 0;
	replace_solution(sol_base, sol_cur);
	for (int local_iter = 1; local_iter <= ts_iteraion; local_iter++)
	{
		global_iteration += 1;
		int min_u_p, min_v_p, min_mach_i, min_makespan_p;
		MOVE_TYPE min_move_type;
		change_permutation(sol_cur, sol_best, min_u_p, min_v_p, min_mach_i, min_makespan_p, min_move_type);

		int min_u_a, min_mach_u_a, min_v_a, min_mach_v_a, min_makespan_a;
		change_machine(sol_cur, sol_best, min_u_a, min_mach_u_a, min_v_a, min_mach_v_a, min_makespan_a);
		if (min_makespan_p < min_makespan_a)
		{
			permutation_iteration += 1;
			update_permutation_tabu(sol_cur, min_mach_i, min_u_p, min_v_p);
			apply_permutation_move(sol_cur, min_mach_i, min_u_p, min_v_p, min_move_type);
		}
		else if (min_makespan_p >= min_makespan_a&&min_makespan_a != INT_MAX)
		{
			/*if (global_iteration == 450226)
			check_solution(sol_cur);*/
			assignment_iteration += 1;
			int r1 = MAX((makespan[sol_cur] - imp_lb) / d3, d4);	// optimize by making it global
			assignment_tabu_list[sol_cur][machine[sol_cur][min_mach_u_a][min_u_a]->job_i][machine[sol_cur][min_mach_u_a][min_u_a]->oper_job_i][min_mach_v_a] =
				assignment_iteration + tt1 + rand() % r1;
			apply_assign_move(sol_cur, min_mach_u_a, min_u_a, min_mach_v_a, min_v_a);
			//cout << global_iteration << "\t" << assignment_iteration << "\t" << min_mach_u_a << "\t" << min_u_a << "\t"
			//	<< min_mach_v_a << "\t" << min_v_a << "\t"
			//	<< min_makespan_a << "\t"
			//	<< makespan[sol_cur] << "\t" << makespan[sol_best] << "\t"
			//	//<< crit_block[sol_cur][0][0] << "\t"
			//	<< (end_time - start_time) / CLOCKS_PER_SEC
			//	<< "\tA"
			//	<< endl;
		}
		else if (min_makespan_p == min_makespan_a&&min_makespan_a == INT_MAX)
			perturb(sol_cur, sol_best, small_ptr, small_ptr_type);
		calculate_r(sol_cur);
		calculate_q_crit_block(sol_cur);
		if (makespan[sol_base] > makespan[sol_cur])
		{
			end_time = clock();
			replace_solution(sol_base, sol_cur);
			//check_solution(sol_base);
			local_iter = 1;
			/*cout << global_iteration << "\t" << permutation_iteration << "\t" << min_mach_i << "\t"
			<< min_u_p << "\t" << min_v_p << "\t"
			<< (min_move_type == BACKWARD_INSERT ? "b" : "f") << "\t"
			<< min_makespan_p << "\t"
			<< makespan[sol_cur] << "\t" << makespan[sol_base] << "\t" << makespan[sol_best] << "\t"
			<< (end_time - start_time) / CLOCKS_PER_SEC
			<< endl;*/
			/*if (makespan[sol_best] > makespan[sol_cur])
			{
			replace_solution(sol_best, sol_cur);
			cout << global_iteration << "\t" << permutation_iteration << "\t" << min_mach_i << "\t"
			<< min_u_p << "\t" << min_v_p << "\t"
			<< (min_move_type == BACKWARD_INSERT ? "b" : "f") << "\t"
			<< min_makespan_p << "\t"
			<< makespan[sol_cur] << "\t" << makespan[sol_base] << "\t" << makespan[sol_best] << "\t"
			<< (end_time - start_time) / CLOCKS_PER_SEC
			<< endl;
			}*/
		}
		/*if (ts_iter >= 4)
		cout << ns_p_cnt << "\t" << ns_a_cnt << endl;*/
		/*if (ns_p_cnt == 0 && ns_a_cnt == 0)
		cout << endl;*/
		/*if (global_iteration == 2147946)
		cout << endl;*/
	}
	replace_solution(sol_cur, sol_base);
	memset(assignment_tabu_list[sol_cur], 0, sizeof(assignment_tabu_list[sol_cur]));
	clear_tabu_list(sol_cur);
}
void Solver::save_best(int best, int c1, int c2)
{
	if (makespan[best] > makespan[c1])
	{
		replace_solution(best, c1);
	}
	else if (makespan[best] == makespan[c1])
	{
		if (rand() % 2)
			replace_solution(best, c1);
	}
	if (c1 != c2)
	{
		if (makespan[best] > makespan[c2])
		{
			replace_solution(best, c2);
		}
		else if (makespan[best] == makespan[c2])
		{
			if (rand() % 2)
				replace_solution(best, c2);
		}
	}
}
void Solver::save_solution(int sol_index)
{
	end_time = clock();
	ofs << run_i << "\t" << gen << "\t" << global_iteration << "\t"
		<< permutation_iteration << "\t" << assignment_iteration << "\t"
		<< (double)(end_time - start_time) / CLOCKS_PER_SEC << "\t"
		<< makespan[sol_index]
		<< endl;
	ofs.flush();
}
void Solver::H2O()
{
	int sol_p1 = 1, sol_p2 = 2, sol_elite1 = 3, sol_elite2 = 4, sol_pr = 5, sol_best = 6, sol_temp = 7;
	for (run_i = 1; run_i <= run_cnt; run_i++)
	{
		start_time = clock();
		global_iteration = permutation_iteration = assignment_iteration = 0;
		init_solution1(sol_p1);
		init_solution1(sol_p2);
		init_solution1(sol_elite1);
		init_solution1(sol_elite2);
		replace_solution(sol_best, sol_p1);

		tabu_search(sol_p1, sol_best);
		tabu_search(sol_p2, sol_best);
		save_best(sol_best, sol_p1, sol_p2);
		cout << "makespan: "
			<< makespan[sol_p1] << "\t" << makespan[sol_p2] << endl;
		int no_move_cnt = 0, non_consec_imp = 0;
		for (gen = 1; gen <= ts_restart; gen++)
		{
			path_relinking(sol_p1, sol_p2, sol_pr, 0);
			replace_solution(sol_temp, sol_pr);
			//path_relinking(sol_p2, sol_p1, sol_pr, 0);
			replace_solution(sol_p2, sol_pr);
			replace_solution(sol_p1, sol_temp);

			int a_dis, p_dis;
			calculate_distance(sol_p1, sol_p2, a_dis, p_dis);
			//cout << a_dis << "\t" << p_dis << endl;
			/*calculate_distance(sol_p1, sol_pr, a_dis, p_dis);
			cout << a_dis << "\t" << p_dis << endl;
			calculate_distance(sol_pr, sol_p2, a_dis, p_dis);
			cout << a_dis << "\t" << p_dis << endl;*/


			tabu_search(sol_p1, sol_best);
			tabu_search(sol_p2, sol_best);

			save_best(sol_elite1, sol_p1, sol_p2);
			//save_best(sol_best, sol_elite1, sol_elite1);

			if (makespan[sol_best] > makespan[sol_elite1])
			{
				non_consec_imp = 0;
				replace_solution(sol_best, sol_elite1);
				end_time = clock();
				save_solution(sol_best);
			}
			else if (makespan[sol_best] == makespan[sol_elite1])
			{
				if (rand() % 2)
					replace_solution(sol_best, sol_elite1);
			}
			else
				non_consec_imp += 1;
			if ((non_consec_imp + 1) % (2 * iter_cyc) == 0)
				init_solution1(sol_p2);
			if (gen%iter_cyc == 0)
			{
				replace_solution(sol_p1, sol_elite2);
				replace_solution(sol_elite2, sol_elite1);
				init_solution(sol_elite1);
				//cout << "*";
			}

			end_time = clock();
			if ((end_time - start_time) / CLOCKS_PER_SEC > max_t ||
				makespan[sol_best] == imp_lb)
				break;
			/*cout << gen << "\t" << (end_time - start_time) / CLOCKS_PER_SEC << "\t" << no_move_cnt << "\t"
			<< makespan[sol_best] << "\t"
			<< " **********************************" << endl;*/
		}
		/*end_time = clock();
		cout << permutation_iteration << "\t"
		<< makespan[sol_p1] << "\t" << makespan[sol_best] << "\t" << best_known_makespan << "\t"
		<< (end_time - start_time) / CLOCKS_PER_SEC
		<< endl;*/
	}
}
void Solver::ITS()
{
	int sol_p1 = 1, sol_p2 = 2, sol_best = 3, sol_pr = 4;
	for (run_i = 1; run_i <= run_cnt; run_i++)
	{
		start_time = clock();
		global_iteration = permutation_iteration = assignment_iteration = 0;
		init_solution1(sol_p1);
		replace_solution(sol_best, sol_p1);
		tabu_search(sol_p1, sol_best);
		int non_consec_imp = 0;
		for (gen = 1; gen <= ts_restart; gen++)
		{
			tabu_search(sol_p1, sol_best);
			if (makespan[sol_best] > makespan[sol_p1])
			{
				replace_solution(sol_best, sol_p1);
				save_solution(sol_best);
				non_consec_imp = 0;
			}
			non_consec_imp += 1;
			if (non_consec_imp != 0 && non_consec_imp % d6 == 0)
			{
#ifdef _WIN32
				cout << "***init\t" << non_consec_imp << endl;
#endif
				init_solution1(sol_p1);
			}
			else if (non_consec_imp != 0 && non_consec_imp % d5 == 0)
			{
#ifdef _WIN32
				cout << "***replace best\t" << non_consec_imp << endl;
#endif
				replace_solution(sol_p1, sol_best);
				//perturb(sol_p1, sol_best, 10, 0);
			}
			perturb(sol_p1, sol_best, big_ptr, big_ptr_type);
			if (makespan[sol_best] > makespan[sol_p1])
			{
				replace_solution(sol_best, sol_p1);
				save_solution(sol_best);
				non_consec_imp = 0;
			}
			//path_relinking(sol_p1, sol_best, sol_pr, 0);
			//replace_solution(sol_p1, sol_pr);
			if ((end_time - start_time) / CLOCKS_PER_SEC > max_t)
			{
#ifdef _WIN32
				cout << "Time exceeds." << endl;
#endif
				break;
			}
			if (makespan[sol_best] == imp_lb)
			{
#ifdef _WIN32
				cout << "Obtain the lower bound." << endl;
#endif
				break;
			}
#ifdef _WIN32
			end_time = clock();
			cout << run_i << "\t"
				<< gen << "\t"
				<< (end_time - start_time) / CLOCKS_PER_SEC << "\t"
				<< makespan[sol_best] << "\t"
				<< " *****" << endl;
#endif
		}
#ifdef _WIN32
		end_time = clock();
		cout << permutation_iteration << "\t"
			<< makespan[sol_p1] << "\t" << makespan[sol_best] << "\t" << imp_lb << "\t"
			<< (end_time - start_time) / CLOCKS_PER_SEC
			<< endl;
#endif
	}
}
void Solver::calculate_distance(int sol_a, int sol_b, int &a_dis, int &p_dis)
{
	a_dis = p_dis = 0;
	for (int mach_a = 1; mach_a <= instance->m; mach_a++)
	{
		for (int op_a = 1; op_a <= machine_oper_num[sol_a][mach_a]; op_a++)
		{
			Operation *oper_a = machine[sol_a][mach_a][op_a];
			if (job[sol_b][oper_a->job_i][oper_a->oper_job_i]->mach_i != mach_a)
			{
				a_dis += 1;
				//cout << oper_a->job_i << "," << oper_a->oper_job_i << endl;
			}
			else
				p_dis += abs(job[sol_b][oper_a->job_i][oper_a->oper_job_i]->oper_mach_i - op_a);
		}
	}
}
void Solver::find_common_seq(int sol_s, int sol_g)
{

	memset(common_seq, 0, sizeof(common_seq));
	int com_s = 0, com_g = 1;
	// find the common operation sequences
	for (int mach_s = 1; mach_s <= instance->m; mach_s++)
	{
		int consec_cnt = 0, seq_cnt = 0;
		for (int op_s = 1, op_g = 1; op_s <= machine_oper_num[sol_s][mach_s] && op_g <= machine_oper_num[sol_g][mach_s];)
		{
			Operation *oper_s = machine[sol_s][mach_s][op_s];
			Operation *oper_g = job[sol_g][oper_s->job_i][oper_s->oper_job_i];
			if (oper_g->mach_i != mach_s || oper_g->oper_mach_i < op_g)
			{
				consec_cnt = 0;
				op_s += 1;
				continue;
			}
			op_g = oper_g->oper_mach_i;

			common_seq[com_s][mach_s][seq_cnt * 2 + 1] = op_s;
			common_seq[com_g][mach_s][seq_cnt * 2 + 1] = op_g;

			consec_cnt += 1;
			op_s += 1;
			op_g += 1;
			while (machine[sol_g][mach_s][op_g]->job_i == machine[sol_s][mach_s][op_s]->job_i&&
				machine[sol_g][mach_s][op_g]->oper_job_i == machine[sol_s][mach_s][op_s]->oper_job_i&&
				machine[sol_g][mach_s][op_g] != dummy_oper[sol_g][END] && machine[sol_s][mach_s][op_s] != dummy_oper[sol_s][END])
			{
				consec_cnt += 1;
				op_s += 1;
				op_g += 1;
			}
			if (consec_cnt == 1)
			{
				if (seq_cnt > 0)
					op_g = common_seq[com_g][mach_s][(seq_cnt - 1) * 2 + 2] + 1;	// previous seq of g
				else
					op_g = 1;
			}
			else
			{
				common_seq[com_s][mach_s][seq_cnt * 2 + 2] = op_s - 1;
				common_seq[com_g][mach_s][seq_cnt * 2 + 2] = op_g - 1;
				seq_cnt += 1;
				op_g -= 1;
			}
			consec_cnt = 0;
		}
		common_seq[com_s][mach_s][0] = common_seq[com_g][mach_s][0] = seq_cnt;
	}
	/*cout << "common seq:" << endl;
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
	for (int seq_i = 1; seq_i <= common_seq[com_s][mach_i][0]; seq_i++)
	{
	for (int op_i = common_seq[com_s][mach_i][(seq_i - 1) * 2 + 1]; op_i <= common_seq[com_s][mach_i][(seq_i - 1) * 2 + 2]; op_i++)
	{
	Operation *oper = machine[sol_s][mach_i][op_i];
	cout << oper->job_i << "," << oper->oper_job_i << "\t";
	}
	cout << " | ";
	}
	cout << endl << endl;
	}

	cout << "common seq g:" << endl;
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
	for (int seq_i = 1; seq_i <= common_seq[com_g][mach_i][0]; seq_i++)
	{
	for (int op_i = common_seq[com_g][mach_i][(seq_i - 1) * 2 + 1]; op_i <= common_seq[com_g][mach_i][(seq_i - 1) * 2 + 2]; op_i++)
	{
	Operation *oper = machine[sol_g][mach_i][op_i];
	cout << oper->job_i << "," << oper->oper_job_i << "\t";
	}
	cout << " | ";
	}
	cout << endl << endl;
	}*/
}
void Solver::path_relinking(int sol_s, int sol_g, int sol_pr, int temp)
{
	replace_solution(sol_pr, sol_s);
	int pl = 0;
	int mach_p_dis = 1, mach_a_dis = 0;
	while (mach_p_dis!=0)
	{
		pl += 1;
		vector<Operation*> oper_ass_vec;
		int p_dis = 0;
		for (int mach_i = 1; mach_i <= instance->m; mach_i++)
		{
			mach_p_dis = 0, mach_a_dis = 0;
			for (int oper_i = 1; oper_i <= machine_oper_num[sol_pr][mach_i]; oper_i++)
			{
				Operation *oper_u = machine[sol_pr][mach_i][oper_i];
				if (job[sol_g][oper_u->job_i][oper_u->oper_job_i]->mach_i != mach_i)
				{
					oper_ass_vec.push_back(oper_u);
					mach_a_dis += 1;
				}
				else
				{
					if (job[sol_g][oper_u->job_i][oper_u->oper_job_i]->oper_mach_i != oper_i)
					{
						p_dis += 1;
						mach_p_dis += 1;
					}
				}
			}
			cout << mach_i << ": " << mach_p_dis << "\t" << mach_a_dis << endl;
		}
		cout <<pl<< "\ta_dis: " << oper_ass_vec.size() << "\t"
			<< "p_dis: " << p_dis << "\t"
			<< endl;
		int min_ass_dis = INT_MAX, ass_equ_cnt=1, min_ass_mach_u, min_ass_u, min_ass_mach_v, min_ass_v;
		for (int diff_i = 0; diff_i < oper_ass_vec.size(); diff_i++)
		{
			Operation *oper_u = oper_ass_vec[diff_i];
			int mach_v = job[sol_g][oper_u->job_i][oper_u->oper_job_i]->mach_i;
			int left_i = 0, right_i = machine_oper_num[sol_pr][mach_v] + 1;
			for (int m_i = machine_oper_num[sol_pr][mach_v]; m_i >= 1; m_i--)	//determine right_i
			{
				Operation *oper_v = machine[sol_pr][mach_v][m_i];
				if (oper_v->end_time > oper_u->pre_job_oper->end_time)	// move u before v
																		//if (v->end_time <= oper_u->pre_job_oper->start_time)
				{
					right_i = m_i;
					//break;
				}
				else break;
			}
			for (int m_i = 1; m_i <= machine_oper_num[sol_pr][mach_v]; m_i++)	// determine left_i
			{
				Operation *oper_v = machine[sol_pr][mach_v][m_i];	// move u after v
				if (oper_v->q + oper_v->t > oper_u->next_job_oper->q + oper_u->next_job_oper->t)
					//if (oper_u->next_job_oper->q >= v->q + v->t)
				{
					left_i = m_i;
					//break;
				}
				else break;
			}
			if (left_i > right_i)
			{
				int temp = left_i;
				left_i = right_i;
				right_i = temp;
			}
			else if (left_i == right_i)
			{
				left_i -= 1;
				right_i += 1;
			}

			int seq_dis_rm = 0;
			// for u in the first machine
			for (int oper_i = oper_u->oper_mach_i + 1; oper_i <= machine_oper_num[sol_pr][oper_u->mach_i]; oper_i++)
			{
				Operation *oper_w = machine[sol_pr][oper_u->mach_i][oper_i];
				if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->mach_i != oper_u->mach_i)
					continue;
				if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->oper_mach_i == oper_i)
					seq_dis_rm += 1;
				else if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->oper_mach_i == oper_i - 1)
					seq_dis_rm -= 1;
			}
			//cout << left_i << "\t" << right_i << "\tseq_dis_rm: " << seq_dis_rm << endl;
			for (int v = left_i; v <right_i; v++)
			{
				int seq_dis_add = 0;
				// for u moved to pos of the second machine
				if (job[sol_g][oper_u->job_i][oper_u->oper_job_i]->oper_mach_i != v+1)
					seq_dis_add += 1;
				for (int pos_i = v + 1; pos_i <= machine_oper_num[sol_pr][mach_v]; pos_i++)
				{
					Operation *oper_w = machine[sol_pr][mach_v][pos_i];
					if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->mach_i != mach_v)
						continue;
					if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->oper_mach_i == pos_i)
						seq_dis_add += 1;
					else if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->oper_mach_i == pos_i + 1)
						seq_dis_add -= 1;
				}
				if (min_ass_dis > seq_dis_rm + seq_dis_add)
				{
					min_ass_dis = seq_dis_rm + seq_dis_add;
					min_ass_mach_u = oper_ass_vec[diff_i]->mach_i;
					min_ass_u = oper_ass_vec[diff_i]->oper_mach_i;
					min_ass_mach_v = mach_v;
					min_ass_v = v;	// move u after pos
					ass_equ_cnt = 1;
				}
				else if (min_ass_dis == seq_dis_rm + seq_dis_add)
				{
					ass_equ_cnt += 1;
					if (rand() % ass_equ_cnt == 0)
					{
						min_ass_mach_u = oper_ass_vec[diff_i]->mach_i;
						min_ass_u = oper_ass_vec[diff_i]->oper_mach_i;
						min_ass_mach_v = mach_v;
						min_ass_v = v;	// move u after pos
						ass_equ_cnt = 1;
					}
				}
			}
		}


		//********* find the best sequence move
		int min_seq_dis = INT_MAX, equ_cnt = 1, min_mach, min_u, min_v;
		MOVE_TYPE move_type, min_move_type;
		for (int mach_i = 1; mach_i <= instance->m; mach_i++)
		{
			for (int u = 1; u < machine_oper_num[sol_pr][mach_i]; u++)
			{
				for (int v = u + 1; v <= machine_oper_num[sol_pr][mach_i]; v++)
				{
					if (machine[sol_pr][mach_i][u]->job_i == machine[sol_pr][mach_i][v]->job_i)
						continue;
					// move u after v
					if (machine[sol_pr][mach_i][v]->start_time < machine[sol_pr][mach_i][u]->next_job_oper->end_time) 	// q[v] + t[v] > q[JS[u]] || r[JS[u]]+t[JS[u]]>r[v]
					{
						move_type = BACKWARD_INSERT;
						int seq_dis_back_insert = 0;
						Operation *oper_w = machine[sol_pr][mach_i][u];
						if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->mach_i == mach_i)
						{
							if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->oper_mach_i == u)
								seq_dis_back_insert += 1;
							else if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->oper_mach_i == v)
								seq_dis_back_insert -= 1;
						}
						for (int w = u + 1; w <= v; w++)
						{
							oper_w = machine[sol_pr][mach_i][w];
							if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->mach_i != mach_i)
								continue;
							if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->oper_mach_i == w)
								seq_dis_back_insert += 1;
							else if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->oper_mach_i == w - 1)
								seq_dis_back_insert -= 1;
						}
						if (min_seq_dis > seq_dis_back_insert)
						{
							min_seq_dis = seq_dis_back_insert;
							min_move_type = move_type;
							min_mach = mach_i;
							min_u = u;
							min_v = v;
							equ_cnt = 1;
						}
						else if (min_seq_dis == seq_dis_back_insert)
						{
							equ_cnt += 1;
							if (rand() % equ_cnt == 0)
							{
								min_move_type = move_type;
								min_mach = mach_i;
								min_u = u;
								min_v = v;
							}
						}
					}
					if (v == u + 1)
						continue;
					// move v before u
					if (machine[sol_pr][mach_i][u]->end_time > machine[sol_pr][mach_i][v]->pre_job_oper->start_time)
					{
						move_type = FORWARD_INSERT;
						int seq_dis_fore_insert = 0;
						Operation *oper_w = machine[sol_pr][mach_i][v];
						if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->mach_i == mach_i)
						{
							if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->oper_mach_i == v)
								seq_dis_fore_insert += 1;
							else if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->oper_mach_i == u)
								seq_dis_fore_insert -= 1;
						}
						for (int w = u; w < v; w++)
						{
							oper_w = machine[sol_pr][mach_i][w];
							if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->oper_mach_i == w)
								seq_dis_fore_insert += 1;
							else if (job[sol_g][oper_w->job_i][oper_w->oper_job_i]->oper_mach_i == w + 1)
								seq_dis_fore_insert -= 1;
						}
						if (min_seq_dis > seq_dis_fore_insert)
						{
							min_seq_dis = seq_dis_fore_insert;
							min_move_type = move_type;
							min_mach = mach_i;
							min_u = u;
							min_v = v;
							equ_cnt = 1;
						}
						else if (min_seq_dis == seq_dis_fore_insert)
						{
							equ_cnt += 1;
							if (rand() % equ_cnt == 0)
							{
								min_move_type = move_type;
								min_mach = mach_i;
								min_u = u;
								min_v = v;
							}
						}
					}
				}
			}
		}
		if ((min_seq_dis < min_ass_dis&&min_seq_dis<0)||min_ass_dis==INT_MAX)
		{
			apply_permutation_move(sol_pr, min_mach, min_u, min_v, min_move_type);
			cout << min_seq_dis << "\t" << min_ass_dis << "\t"
				<< "PERMU: " << min_mach << "\t" << min_u << "\t" << min_v << "\t" << min_move_type << endl;
		}
		else if (min_seq_dis > min_ass_dis&&min_ass_dis<0)
		{
			apply_assign_move(sol_pr, min_ass_mach_u, min_ass_u, min_ass_mach_v, min_ass_v);
			cout << min_seq_dis << "\t" << min_ass_dis << "\t" 
				<< "ASSI: " << min_ass_mach_u << "\t" << min_ass_u << "\t" << min_ass_mach_v << "\t" << min_ass_v << endl;
		}
		else
		{
			Operation *oper_u = machine[sol_pr][min_ass_mach_u][min_ass_u];	// insert oper_u after v
			Operation *oper_v = machine[sol_pr][min_ass_mach_v][min_ass_v];
			oper_u->apx_r = MAX(oper_u->pre_job_oper->end_time, oper_v->end_time);
			oper_u->apx_q = MAX(oper_u->next_job_oper->q + oper_u->next_job_oper->t,
				machine[sol_pr][min_ass_mach_v][min_ass_v + 1]->q + machine[sol_pr][min_ass_mach_v][min_ass_v + 1]->t);	// v->next_mach_oper
			int u_t = job[sol_g][oper_u->job_i][oper_u->oper_job_i]->t;
			int makespan_ass = oper_u->apx_r + oper_u->apx_q + u_t;
			int makespan_seq;
			if (min_move_type == BACKWARD_INSERT)
				try_backward_insert_move(sol_pr, makespan_seq, min_mach, min_u, min_v);
			else
				try_forward_insert_move(sol_pr, makespan_seq, min_mach, min_u, min_v);
			if (makespan_seq < makespan_ass)
			{
				apply_permutation_move(sol_pr, min_mach, min_u, min_v, min_move_type);
				cout << min_seq_dis << "\t" << min_ass_dis << "\t" 
					<< makespan_seq << "\t" << makespan_ass << "\t"
					<< "PERMU mm: " << min_mach << "\t" << min_u << "\t" << min_v << "\t" << min_move_type << endl;
			}
			else
			{
				apply_assign_move(sol_pr, min_ass_mach_u, min_ass_u, min_ass_mach_v, min_ass_v);
				cout << min_seq_dis << "\t" << min_ass_dis << "\t" 
					<< makespan_seq << "\t" << makespan_ass << "\t"
					<< "ASSI mm: " << min_ass_mach_u << "\t" << min_ass_u << "\t" << min_ass_mach_v << "\t" << min_ass_v << endl;
			}
		}
		calculate_r(sol_pr);
		calculate_q_crit_block(sol_pr);
		//cout << makespan[sol_pr] << endl;
		check_solution(sol_pr);
	}
	/*int a_dis, p_dis;
	calculate_distance(sol_s, sol_g, a_dis, p_dis);
	cout << a_dis << "\t" << p_dis << endl;
	calculate_distance(sol_s, sol_pr, a_dis, p_dis);
	cout << a_dis << "\t" << p_dis << endl;
	calculate_distance(sol_pr, sol_g, a_dis, p_dis);
	cout << a_dis << "\t" << p_dis << endl;*/

}
void Solver::perturb(int sol_ptr, int sol_best, int ptr_len, int ptr_type)
{
	int  equ_cnt, min_u, min_v, min_mach_i;
	int min_u_a, min_mach_u_a, min_v_a, min_mach_v_a;
	MOVE_TYPE min_move_type;
	bool change_permutation;
	for (int local_iter = 1; local_iter <= ptr_len; local_iter++)
	{
		equ_cnt = 0;
		for (int block_i = 1; block_i <= crit_block[sol_ptr][0][0]; block_i++)
		{
			int mach_i = crit_block[sol_ptr][block_i][0];
			int u = crit_block[sol_ptr][block_i][1];
			for (int v = u + 1; v <= crit_block[sol_ptr][block_i][2]; v++)
			{
				if (machine[sol_ptr][mach_i][u]->job_i == machine[sol_ptr][mach_i][v]->job_i)
					continue;
				if (machine[sol_ptr][mach_i][v]->q + machine[sol_ptr][mach_i][v]->t >
					machine[sol_ptr][mach_i][u]->next_job_oper->q)	// q[v] + t[v] > q[JS[u]]
				{
					// move u behind v
					equ_cnt += 1;
					if (rand() % equ_cnt == 0)
					{
						min_u = u;
						min_v = v;
						min_mach_i = mach_i;
						min_move_type = BACKWARD_INSERT;
						change_permutation = true;
					}
				}
				if (machine[sol_ptr][mach_i][u]->end_time > machine[sol_ptr][mach_i][v]->pre_job_oper->start_time)	// r[u] + t[u] > r[JP[v]]
				{
					// move v before u
					equ_cnt += 1;
					if (rand() % equ_cnt == 0)
					{
						min_u = u;
						min_v = v;
						min_mach_i = mach_i;
						min_move_type = FORWARD_INSERT;
						change_permutation = true;
					}
				}
			}

			int v = crit_block[sol_ptr][block_i][2];
			for (u = crit_block[sol_ptr][block_i][1] + 1; u < v; u++)
			{
				if (machine[sol_ptr][mach_i][u]->job_i == machine[sol_ptr][mach_i][v]->job_i)
					continue;
				if (machine[sol_ptr][mach_i][v]->q + machine[sol_ptr][mach_i][v]->t >
					machine[sol_ptr][mach_i][u]->next_job_oper->q)	// q[v] + t[v] > q[JS[u]]
				{
					// move u behind v
					equ_cnt += 1;
					if (rand() % equ_cnt == 0)
					{
						min_u = u;
						min_v = v;
						min_mach_i = mach_i;
						min_move_type = BACKWARD_INSERT;
						change_permutation = true;
					}
				}
				if (machine[sol_ptr][mach_i][u]->end_time > machine[sol_ptr][mach_i][v]->pre_job_oper->start_time)	// r[u] + t[u] > r[JP[v]]
				{
					// move v before u
					equ_cnt += 1;
					if (rand() % equ_cnt == 0)
					{
						min_u = u;
						min_v = v;
						min_mach_i = mach_i;
						min_move_type = FORWARD_INSERT;
						change_permutation = true;
					}
				}
			}
			if (ptr_type == 1)
				continue;
			// change machine 
			for (int u = crit_block[sol_ptr][block_i][1]; u <= crit_block[sol_ptr][block_i][2]; u++)	// for each critical operation
			{
				int job_i = machine[sol_ptr][mach_i][u]->job_i;
				int oper_job_i = machine[sol_ptr][mach_i][u]->oper_job_i;
				if (instance->job_vec[job_i][oper_job_i].size() == 2)
					continue;
				for (int proc_i = 1; proc_i < instance->job_vec[job_i][oper_job_i].size(); proc_i++)	// for each process
				{
					if (instance->job_vec[job_i][oper_job_i][proc_i]->mach_i == mach_i)	// the same machine
						continue;
					int mach_v = instance->job_vec[job_i][oper_job_i][proc_i]->mach_i;
					int left_i = 0, right_i = machine_oper_num[sol_ptr][mach_v] + 1;
					Operation *oper_u = machine[sol_ptr][mach_i][u];

					for (int m_i = machine_oper_num[sol_ptr][mach_v]; m_i >= 1; m_i--)	//determine right_i
					{
						Operation *oper_v = machine[sol_ptr][mach_v][m_i];
						if (oper_v->end_time > oper_u->pre_job_oper->end_time)	// move u before v
																				//if (v->end_time <= oper_u->pre_job_oper->start_time)
							right_i = m_i;
						else break;
					}

					for (int m_i = 1; m_i <= machine_oper_num[sol_ptr][mach_v]; m_i++)	// determine left_i
					{
						Operation *oper_v = machine[sol_ptr][mach_v][m_i];	// move u after v
						if (oper_v->q + oper_v->t > oper_u->next_job_oper->q + oper_u->next_job_oper->t)
							//if (oper_u->next_job_oper->q >= v->q + v->t)
							left_i = m_i;
						else break;
					}
					if (left_i > right_i)
					{
						int temp = left_i;
						left_i = right_i;
						right_i = temp;
						//left_i -= 1;
					}

					for (int v = left_i; v < right_i; v++)
					{
						equ_cnt += 1;
						if (rand() % equ_cnt == 0)
						{
							min_mach_u_a = mach_i;
							min_u_a = u;
							min_mach_v_a = mach_v;
							min_v_a = left_i + rand() % (right_i - left_i);
							change_permutation = false;
						}
					}
				}
			}
		}
		if (equ_cnt == 0)
		{
			//change_machine(sol_ptr, sol_best);
			continue;
		}
		if (change_permutation == true)
		{
			/*permutation_iteration += 1;*/
			//update_permutation_tabu(sol_ptr, min_mach_i, min_u, min_v);
			apply_permutation_move(sol_ptr, min_mach_i, min_u, min_v, min_move_type);
		}
		else
		{
			/*assignment_iteration += 1;*/
			/*assignment_tabu_list[sol_ptr][machine[sol_ptr][min_mach_u_a][min_u_a]->job_i][machine[sol_ptr][min_mach_u_a][min_u_a]->oper_job_i][min_mach_v_a] =
			assignment_iteration + tt1 + rand() % d4;*/
			apply_assign_move(sol_ptr, min_mach_u_a, min_u_a, min_mach_v_a, min_v_a);

		}
		calculate_r(sol_ptr);
		calculate_q_crit_block(sol_ptr);
		//check_solution(sol_ptr);
	}
}// calculate r and makespan
void Solver::calculate_r(int sol_index)
{
	memset(critical_flag, -1, sizeof(critical_flag));
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
	if (oper_num != instance->total_num_operation)
	{
		cout << "ERROR: the number of operations is wrong pr1" << endl;
		system("pause");
	}
	oper_num = 0;
	while (front != rear)
	{
		Operation *oper = queue[front++];
		oper->start_time = MAX(oper->pre_job_oper->end_time, machine[sol_index][oper->mach_i][oper->oper_mach_i - 1]->end_time);
		oper->end_time = oper->start_time + oper->t;
		oper_num += 1;
		//if (oper->next_job_oper != dummy_oper[sol_index][END])
		critical_flag[oper->next_job_oper->mach_i][oper->next_job_oper->oper_mach_i] -= 1;
		//if(machine[sol_index][oper->mach_i][oper->oper_mach_i+1]!=dummy_oper[sol_index][END])
		critical_flag[oper->mach_i][oper->oper_mach_i + 1] -= 1;
		if (critical_flag[oper->next_job_oper->mach_i][oper->next_job_oper->oper_mach_i] == 0)
			queue[rear++] = oper->next_job_oper;
		if (critical_flag[oper->mach_i][oper->oper_mach_i + 1] == 0 &&
			machine[sol_index][oper->mach_i][oper->oper_mach_i + 1] != oper->next_job_oper)
			queue[rear++] = machine[sol_index][oper->mach_i][oper->oper_mach_i + 1];

	}
	/*if (front != instance->total_num_operation)
	cout << endl;*/
	if (oper_num != instance->total_num_operation)
	{
		cout << "ERROR: the number of operations is wrong cr2\t" << argv_map.at("_rand_seed") << endl;
		system("pause");
	}
	makespan[sol_index] = 0;
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		if (makespan[sol_index] < machine[sol_index][mach_i][machine_oper_num[sol_index][mach_i]]->end_time)
			makespan[sol_index] = machine[sol_index][mach_i][machine_oper_num[sol_index][mach_i]]->end_time;
	}
}
// determine critical path, and calculate q for all operations
void Solver::calculate_q_crit_block(int sol_index)
{

	// calculate q for all operations
	/*int front = 0, rear = 0;
	front = 0, rear = 0;
	memset(critical_flag, 0, sizeof(critical_flag));
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	critical_flag[mach_i][machine_oper_num[sol_index][mach_i]] += 1;
	for (int job_i = 1; job_i <= instance->n; job_i++)
	{
	Operation *oper = job[sol_index][job_i][job_oper_num[sol_index][job_i]];
	critical_flag[oper->mach_i][oper->oper_mach_i] += 1;
	if (critical_flag[oper->mach_i][oper->oper_mach_i] == 2)
	queue[rear++] = oper;
	}
	while (front != rear)
	{
	Operation *cur_oper = queue[front++];
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
	}*/

	memset(critical_flag, 0, sizeof(critical_flag));
	for (int oper_i = instance->total_num_operation - 1; oper_i >= 0; oper_i--)
	{
		Operation *cur_oper = queue[oper_i];
		cur_oper->q = MAX(cur_oper->next_job_oper->q + cur_oper->next_job_oper->t,
			machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i + 1]->q + machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i + 1]->t);
		if (cur_oper->end_time + cur_oper->q == makespan[sol_index])
			critical_flag[cur_oper->mach_i][cur_oper->oper_mach_i] = 1;
	}


	//determine critical path
	//int front = 0, rear = 0;
	//for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	//{
	//	if (makespan[sol_index] == machine[sol_index][mach_i][machine_oper_num[sol_index][mach_i]]->end_time)
	//		queue[rear++] = machine[sol_index][mach_i][machine_oper_num[sol_index][mach_i]];
	//}
	//memset(critical_flag, 0, sizeof(critical_flag));
	//while (front != rear)
	//{
	//	Operation *cur_oper = queue[front++];
	//	if (critical_flag[cur_oper->mach_i][cur_oper->oper_mach_i] == 0)
	//	{
	//		critical_flag[cur_oper->mach_i][cur_oper->oper_mach_i] = 1;
	//		if (cur_oper->pre_job_oper != dummy_oper[sol_index][START] &&
	//			cur_oper->start_time == cur_oper->pre_job_oper->end_time)
	//			queue[rear++] = cur_oper->pre_job_oper;
	//		if (machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i - 1] != dummy_oper[sol_index][START] &&
	//			cur_oper->start_time == machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i - 1]->end_time&&
	//			machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i - 1] != cur_oper->pre_job_oper)
	//			// avoid the previous machine operation is the same as the previous job operation
	//			queue[rear++] = machine[sol_index][cur_oper->mach_i][cur_oper->oper_mach_i - 1];
	//	}

	//}
	int block_cnt = 0;
	for (int mach_i = 1; mach_i <= instance->m; mach_i++)
	{
		for (int oper_i = 1; oper_i <= machine_oper_num[sol_index][mach_i]; oper_i++)	// according to the second index of critical
		{
			if (critical_flag[mach_i][oper_i] == 1 && oper_i <= machine_oper_num[sol_index][mach_i])
			{
				crit_block[sol_index][++block_cnt][0] = mach_i;
				crit_block[sol_index][block_cnt][1] = oper_i;
				while (critical_flag[mach_i][oper_i] == 1 && oper_i <= machine_oper_num[sol_index][mach_i])
					oper_i += 1;
				crit_block[sol_index][block_cnt][2] = oper_i - 1;
				//if (crit_block[sol_index][block_cnt][1] == crit_block[sol_index][block_cnt][2])
				//	block_cnt -= 1;	// a block at least holds two operations
			}
		}
	}
	crit_block[sol_index][0][0] = block_cnt;
}
int main(int argc, char **argv)
{
	int rs = time(NULL);
	rs = 1510173144;// 1510173144  1509981230
	srand(rs);
	string argv_win[] = { "",	// 0
		"_px","ITS_",
		"_ifp", "instances\\",
		"_sfp","solutions\\best_solutions\\",	// solution file path
		"_ofp","solutions\\",	// solution output file path
		"_ifn", "01a",
		"_suffix_ins",".fjs", "_suffix_sol",".txt",
		"_qt_ss","quintiq_sspr",
		"_sol_num", "10", "_tt0","10", "_d1","5", "_d2", "12",
		"_tt1","10", "_d3","5", "_d4", "12","_d5","20","_d6","40",
		"_bptr","10","_bptrt","0","_sptr","4","_sptrt","0",
		"_pr_a","50","_pr_p","100","_iter_cyc","10",
		"_itr","10000","_ts_rs","200000000","_run_cnt","20","_max_t","1800"
	};
	std::map<string, string> argv_map;
	argv_map["_exe_name"] = argv[0];
	argv_map["_job_id"] = argv_win[0];
	argv_map["_rand_seed"] = to_string(rs);
	for (int i = 1; i < sizeof(argv_win) / sizeof(argv_win[0]); i += 2)
		argv_map[string(argv_win[i])] = string(argv_win[i + 1]);
	for (int i = 1; i < argc; i += 2)
		argv_map[string(argv[i])] = string(argv[i + 1]);
	string dir[6] = { "\\","Barnes","Brandimarte_Data","Dauzere_Data","Hurink_Data","DemirkolBenchmarksJobShop" };
#ifdef __linux__
	replace_all(argv_map["_ifp"], "\\", "//");
	replace_all(argv_map["_sfp"], "\\", "//");
	replace_all(argv_map["_ofp"], "\\", "//");
	replace_all(dir[0], "\\", "//");
	argv_map["_job_id"] = argv_map.at("_job_id") + ".";
#endif
	string ifn_ins = argv_map.at("_ifp");
	string ifn_qt_ss = argv_map.at("_sfp") + argv_map.at("_qt_ss") + argv_map.at("_suffix_sol");
	if (argv_map.at("_ifn").find("data-") != string::npos)
		ifn_ins += dir[4];
	else if (argv_map.at("_ifn").find("Mk") != string::npos&&argv_map.at("_ifn").length() == 4)
		ifn_ins += dir[2];
	else if (argv_map.at("_ifn").find("mt") != string::npos || argv_map.at("_ifn").find("set") != string::npos)
		ifn_ins += dir[1];
	else if (argv_map.at("_ifn").find("max_") != string::npos)
		ifn_ins += dir[5];
	else  if (argv_map.at("_ifn").find("a") != string::npos&&argv_map.at("_ifn").length() == 3)
		ifn_ins += dir[3];
	ifn_ins += dir[0] + argv_map.at("_ifn");
	if (argv_map.at("_ifn").find("max_") != string::npos)
		ifn_ins += argv_map.at("_suffix_jsp");
	else
		ifn_ins += argv_map.at("_suffix_ins");
	Instance *ins = new Instance(ifn_ins);
	string fnw = argv_map.at("_ofp") + argv_map.at("_job_id") +
		argv_map.at("_px") + argv_map.at("_ifn");
#ifdef _WIN32
	fnw += "_itc" + argv_map.at("_iter_cyc") +
		"_itr" + argv_map.at("_itr") +
		"_tt0" + argv_map.at("_tt0") +
		"_d1" + argv_map.at("_d1") +
		"_d2" + argv_map.at("_d2") +
		"_tt1" + argv_map.at("_tt1") +
		"_d3" + argv_map.at("_d3") +
		"_d4" + argv_map.at("_d4") +
		"_max_t" + argv_map.at("_max_t") +
		"_ts_rs" + argv_map.at("_ts_rs");
#endif
	fnw += argv_map.at("_suffix_sol");
	Solver *solver = new Solver(*ins, stoi(argv_map.at("_sol_num")), fnw, argv_map, ifn_qt_ss);
	solver->tt0 = stoi(argv_map.at("_tt0"));
	solver->d1 = stoi(argv_map.at("_d1"));
	solver->d2 = stoi(argv_map.at("_d2"));
	solver->tt1 = stoi(argv_map.at("_tt1"));
	solver->d3 = stoi(argv_map.at("_d3"));
	solver->d4 = stoi(argv_map.at("_d4"));
	solver->d5 = stoi(argv_map.at("_d5"));
	solver->d6 = stoi(argv_map.at("_d6"));
	solver->big_ptr = stoi(argv_map.at("_bptr"));
	solver->big_ptr_type = stoi(argv_map.at("_bptrt"));
	solver->small_ptr = stoi(argv_map.at("_sptr"));
	solver->small_ptr_type = stoi(argv_map.at("_sptrt"));
	solver->ts_iteraion = stoi(argv_map.at("_itr"));
	solver->ts_restart = stoi(argv_map.at("_ts_rs"));
	solver->pr_a = stoi(argv_map.at("_pr_a"));
	solver->pr_p = stoi(argv_map.at("_pr_p"));
	solver->iter_cyc = stoi(argv_map.at("_iter_cyc"));
	solver->run_cnt = stoi(argv_map.at("_run_cnt"));
	solver->max_t = stoi(argv_map.at("_max_t"));
	//solver->read_solution(1, argv_map.at("_sfp") + argv_map.at("_sfn") + argv_map.at("_suffix"));
	//solver->init_solution1(1);
	//solver->ts(10000);
	//solver->ITS();
	solver->H2O();
	delete solver;
	delete ins;
#ifdef _WIN32
	system("pause");
#endif
	return 0;
}
#endif
