#if 1
#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<string>
using namespace std;
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
private:
	int m, n;	// m and n are the number of machines and jobs
	int avg_num_mach_per_oper;	// average number of machines per operation.
	class JobInfo 
	{
	public:
		int num_oper;
		class Operation 
		{
			int num_mach;
			class Process {
				int mach;
				int t;
				Process(int _mach, int _t) :mach(_mach), t(_t) {}
			};
			vector<Process*> proc_vec;
			Operation(int _num_mach, vector<Process*> &_proc_vec) :num_mach(_num_mach), proc_vec(_proc_vec) {}
		};
		vector<Operation*> oper_vec;
		JobInfo() {}
		JobInfo(int _num_oper,vector<Operation*>_oper_vec):num_oper(_num_oper),oper_vec(_oper_vec) {}
	};
	vector<JobInfo*> job_vec;
public:
	Instance(string);
};
Instance::Instance(string file_input)
{
	ifstream ifs(file_input);
	if (!ifs.is_open())
	{
		cout << file_input << endl;
		perror("file_input"); exit(0);
	}
	string str_line;
	ifs >> n >> m >> avg_num_mach_per_oper;
	//job_vec = new vector<JobInfo*>;
	vector<string>fields_vec;
	while (getline(ifs,str_line))
	{
		cout << str_line << endl;
		split_generic<string>(fields_vec, str_line, "\t");
		JobInfo *jobinfo = new JobInfo();
		//jobinfo->num_oper = stoi(fields_vec[0]);
		//ifs>>op->num_oper;
	}
	ifs.close();
}
int main(int argc, char **argv)
{
	char *argv_win[] = {"",	// 0
		"_ifp", "instances\\Dauzere_Data\\",
		"_ifn", "01a"
	};
	cout << "This is the flexible job shop scheduling problem" << endl;
#ifdef _WIN32
	argc = sizeof(argv_win) / sizeof(argv_win[0]); argv = argv_win;
#endif
	map<string, string> argv_map;
	argv_map["_exe_name"] = argv[0];	// add the exe file name to argv map to append to the output file name
	for (int i = 1; i < argc; i += 2)
		argv_map[string(argv[i])] = string(argv[i + 1]);
	Instance *ins = new Instance(argv_map.at("_ifp")  + argv_map.at("_ifn") + ".fjs");
#ifdef _WIN32
	system("pause");
#endif
	return 0;
}
#endif