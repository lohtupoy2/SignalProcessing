#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <TLatex.h>
#include <TSpectrum.h>
//----------------------//
using namespace std;
// ---------------------//
//som

//командный пункт
//----------------------------------------------------------------------------------------
int ch_fast = 5;					// номер каннала
int ch_slow = 2;

int run_start = 84;
int run_stop = run_start + 34;
int event_start = 20;
int event_stop = event_start + 20;    // выводит определненные события, сигнал
int run_slow_vs_fast = run_start;			// выводит форму сигнала данного run, для быстрых и медленных
int evt_slow_vs_fast_start = 60;
int evt_slow_vs_fast_stop = 80 -1;

int threshold_slow = /*200*/ 40 /*40*/;
int threshold_fast = 20;

int one_alfa_l = 420E3;   // граница cut для одного пика alfa по total_area
int one_alfa_r = 900E3;
int one_peak_l = -10E3; 	// граница cut для одного пика phe для калибровки
int one_peak_r = 10000E3;
int colib_l = 100E3;
int colib_r = 100000E3;
int left_line_integrate = 32000;  /*32000;*/	// границы cut для ворот калибровки
int right_line_integrate = 32700; /*32700;*/
int seperate_coeff = 1;				//коэффициент разделения площадей

int lefts = 32150; /*32110;*/
int rights = 32360; /*32200;*//*32400; 32250*/
int count_peak = 1; 		//количество пиков выше которого - альфа

bool colibration_mode = 0;
bool alfa = 0;

string file_name = "D:\\Data\\new_setup\\250513_part1_caen_archive\\f3";               				// путь к файлу
string data_path = "C:\\Users\Mikheev\\Desktop\\code_root\\241112\\out_runNumb_eventNumb.txt";		// путь к записи данных в файл
FILE* outFile = NULL;

//----------------------//

const int points_per_event = 9999;
const int events_per_file = 1000;
const int sec_per_point = 16;
const int time_to_integrate = 32050;


const int period_stBaseLine = 10000/sec_per_point;

//вспомогательные структуры - тип данных
//----------------------------------------------------------------------------------------

struct Channel
{
	int fast;
	int slow;

};

Channel ch;
Channel TRIG;
TRIG.fast = 100;
TRIG.slow = 25;


struct Point
{
	int x;
	int y;
};

struct Events
{
	int runs;
	int events;

};

struct RunMode
{
	int start;
	int	stop;
};

// структура для хранения данных для гистограмм
//----------------------------------------------------------------------------------------
struct data_bank_HIST
{
	//#vlad на событие
	vector <double> base_line_avr;
	vector <double> base_line_sigma;
	vector <double>	total_area;
	vector <double> integrated_area;//#vlad gated_area area_in_gate
	//#vlad base_line_y_max

	//#vlad на один пик
	vector <double> area;
	vector <double> x_max;
	vector <double> y_max;
};


struct data_bank_EVENT						// структура для хранения данных для вывода событий
{
	vector <vector<Point>> peaks;
	vector <vector<Point>> peaks_max;
	vector <Events> runs_events;

};

struct data_bank_TIME						// струкутра для временных калибровок
{
	vector <double> area_in_rigeon;
	vector <int> t_start;
	vector <int> t_end;
	vector <int> max;

	TH1F* hist_tau_start = new TH1F("hist_tau_start", "hist_tau_start", 16000, 0, 160000);
	TH1F* hist_tau_max = new TH1F("hist_tau_max", "hist_tau_max", 16000, 0, 160000);
	TH1F* hist_start_end = new TH1F("hist_start_end", "hist_start_end", 350, 0, 3500);
	TH1F* hist_max_start = new TH1F("hist_max_start", "hist_max_start", 600, -1000, 5000);
	TH1F* hist_max_end = new TH1F("hist_max_end", "hist_max_end", 500, 0, 5000);
	TH1F* hist_area_in_region = new TH1F("hist_area_in_region", "hist_area_in_region", 150 , -10E3,1500E3);


};

struct data_bank_CORE
{
	int cut_size = 200;
	int no_cut_size = 500;
	int size = cut_size;

	TH2F* cor_start_vs_area = new TH2F("cor_start_vs_area", "cor_start_vs_area", size*50, 20E3, 60E3, size/2, 0, 1000E3);
	TH2F* cor_max_vs_area = new TH2F("cor_max_vs_area", "cor_max_vs_area", size*50, 20E3, 60E3, size/2, 0, 1000E3);
	TH2F* cor_max_vs_start = new TH2F("cor_max_vs_start", "cor_max_vs_start", 1000, 0, 160E3, 1000, 0, 160E3);
	TH2F* cor_max_start_vs_area = new TH2F("cor_max_start_vs_area", "cor_max_start_vs_area", size*10, -1000, 5000, size, 0, 1000E3);
	TH2F* cor_start_end_vs_area = new TH2F("cor_start_end_vs_area", "cor_start_end_vs_area", size, 0, 5000, size, 0, 1000E3);

	void fillHistograms(const data_bank_TIME& TYPE_CUT)
	{
        for (int i = 0; i < TYPE_CUT.area_in_rigeon.size(); i++)
		{
            cor_start_vs_area->Fill(TYPE_CUT.t_start[i], TYPE_CUT.area_in_rigeon[i]);
            cor_max_vs_area->Fill(TYPE_CUT.max[i], TYPE_CUT.area_in_rigeon[i]);
			cor_max_vs_start->Fill(TYPE_CUT.max[i], TYPE_CUT.t_start[i]);
			cor_max_start_vs_area->Fill(abs(TYPE_CUT.max[i] - TYPE_CUT.t_start[i]), TYPE_CUT.area_in_rigeon[i]);
			cor_start_end_vs_area->Fill(abs(TYPE_CUT.t_end[i] - TYPE_CUT.t_start[i]), TYPE_CUT.area_in_rigeon[i]);
		}
    }

	void ClearHistograms()
	{
		cor_start_vs_area->Reset();
		cor_max_vs_area->Reset();
		cor_max_vs_start->Reset();
		cor_max_start_vs_area->Reset();
		cor_start_end_vs_area->Reset();
	}

};

struct data_bank_PEAKS
{
	vector <int> peaks_number;
};

// объекты структур например:  data_bank_HIST COLIB объект COLIB данных гистограмм для калибровок, SHORT_LIVE - объект для данных, которые стираем копируем для калибровок

data_bank_HIST HIST;
data_bank_HIST COLIB;
data_bank_HIST SHORT_LIVE;

data_bank_EVENT FAST_EVENT;
data_bank_EVENT SINGLE_EVENT;
data_bank_EVENT MULTY_EVENT;
data_bank_EVENT PEAK_COUNT_EVENT;
data_bank_EVENT ALPHA_EVENT;
data_bank_EVENT PEAK;

data_bank_TIME CUT;
data_bank_TIME NO_CUT;
data_bank_TIME WITH_PIEDISTAL;
data_bank_TIME DISSPERSION;
data_bank_TIME ALPHA;

data_bank_CORE CORELATION;

data_bank_PEAKS NUMBER_PEAKS;

//  вспомогательные функции для записи данных в структурный банк
//----------------------------------------------------------------------------------------
void get_runs_events(int run, int event, data_bank_EVENT& EVENT_TYPE)
{
	Events s;
	s.runs = run;
	s.events = event;
	EVENT_TYPE.runs_events.push_back(s);
}

//#vlad убрать
void get_base_line_data(double base_line_avr,double base_line_sigma)
{
	HIST.base_line_avr.push_back(base_line_avr);
	HIST.base_line_sigma.push_back(base_line_sigma);

}
void get_peaks(vector <vector<Point>> p,vector <vector<Point>> p_max)
{
	PEAK.peaks = p;
	PEAK.peaks_max = p_max;

}

//#vlad я бы убрал функцию и писал сразу
void get_area_max(double area, double max_x, double max_y)
{
	HIST.area.push_back(area);
	SHORT_LIVE.area.push_back(area);

	HIST.x_max.push_back(max_x);
	SHORT_LIVE.x_max.push_back(max_x);

	HIST.y_max.push_back(max_y);
	SHORT_LIVE.y_max.push_back(max_y);

}

void get_tstart_tend_tmax(int start, int end, int max, double area, data_bank_TIME& TYPE_CUT)
{
	TYPE_CUT.t_start.push_back(start);
	TYPE_CUT.t_end.push_back(end);
	TYPE_CUT.max.push_back(max);
	TYPE_CUT.area_in_rigeon.push_back(area);

}

void get_integrated_area(double area_integrated)
{
	HIST.integrated_area.push_back(area_integrated);

}


void get_total_area(double total_area, data_bank_HIST& HIST_TYPE)
{
	HIST_TYPE.total_area.push_back(total_area);
}


void copy_data_colibr(void)  	// функция копирования SHORT_LIVE в COLIB для калибровки
{
	 COLIB.area.insert(COLIB.area.end(), SHORT_LIVE.area.begin(), SHORT_LIVE.area.end());
	 COLIB.x_max.insert(COLIB.x_max.end(), SHORT_LIVE.x_max.begin(), SHORT_LIVE.x_max.end());
	 COLIB.y_max.insert(COLIB.y_max.end(), SHORT_LIVE.y_max.begin(), SHORT_LIVE.y_max.end());
}
//----------------------------------------------------------------------------------------

// Вспомогательные функции для чтения, расчета, вывода
//----------------------------------------------------------------------------------------
vector <Events> read_data_file(string _data_path)
{
ifstream inputFile(_data_path);
vector <Events> v_adress;
string line;
	while(getline(inputFile, line))
	{
		istringstream iss(line);
		Events adress;
		if(iss >> adress.runs >> adress.events)
		{
			v_adress.push_back(adress);
		}
	}
	inputFile.close();

return v_adress;
}

int columns (int number_events)
{
	int col;
	if (number_events == 1) col = 1;
	else col = 3;
	if (number_events > 9) col = 4;
	if (number_events > 12) col = 5;
	return col;
}

double* extract_xs(vector<Point> event)
{
	double* xs = new double[event.size()];
	for (int i = 0; i < event.size(); i++)
	{
		xs[i] = sec_per_point*event[i].x;
	}
	return xs;
}

double* extract_ys(vector<Point> event)
{
	double* ys = new double[event.size()];
	for (int i = 0; i < event.size(); i++)
	{
		ys[i] = event[i].y;
	}
	return ys;
}

void alfa_mode(string mode)
{
	if (mode == "FAST")
	{
		alfa = 0;
		cout << "mode_cut is = " << alfa << endl;
	}
	else
	{
		cout << "alfa mode is = ? " << endl
		<< " '1' - on " << endl << " '0' - off " << endl;
		cin >> alfa;
		cout << " mode_cut is = " << alfa << endl;
	}

}
void colibration(string mode)
{
	if(mode == "FAST")
	{
		colibration_mode = 0;
		cout << "colibration_mode = " << colibration_mode << endl;
	}
	else
	{
		cout << "colibration_mode = ?" << endl;
		cin >> colibration_mode;
		cout << "colibration_mode is = " << colibration_mode << endl;
	}

}


RunMode conditions_mode(string choose_mode)
{
	RunMode Run;

	if (choose_mode == "HIST")
	{
		Run.start = run_start;
		Run.stop = run_stop;
	}
	if (choose_mode == "EVENT")
	{
		Run.start = run_start;
		Run.stop  = run_start;
	}
	return Run;
}


bool positive(double numb){
	bool answer = 0;
	if (numb > 0)
	answer = 1;
	return answer;
}

bool negative(double numb){
	bool answer = 0;
	if (numb < 0)
	answer = 1;
	return answer;
}

bool Form_of_signal(vector <Point> signal_peak)		// Сборочная функция по отображению сигнала
{
	int N_segments = 6;
	size_t segment_size = signal_peak.size()/N_segments;
	vector <double> derivetive;

	derivetive.clear();
	//cout << signal_peak.size() <<endl;
	for (int i = 0; i < N_segments; i++)
	{
		size_t start = segment_size*i;
		size_t end = segment_size*(1+i);

		if (i == N_segments - 1) end = signal_peak.size();

		vector<Point> segment(signal_peak.begin() + start, signal_peak.begin() + end);

		if (segment.size() < 2) {
		//cerr << "Segment size is less than 2 for segment index " << i << endl;
		return 1;
		break;
		//continue;
		}

		if (segment.back().x == segment.front().x) {
		cerr << "Division by zero for segment index " << i << endl;
		return 1;
		break;
		//continue;
		}

		derivetive.push_back((segment.back().y - segment.front().y)/(segment.back().x - segment.front().x));
	}

	int counter = 0;
	for (int i = 1; i < derivetive.size(); i++)
	{
		if(positive(derivetive[i - 1]) && negative(derivetive[i]) && abs(derivetive[i] / derivetive[i - 1]) > 1.2)
		{
			counter++;
		}
	}

	if (counter > 1)
		return 0;
	else
		return 1;

}

//функции очистки структур
//----------------------------------------------------------------------------------------
void clear_vector()
	{
		HIST.total_area.clear();
		HIST.total_area.shrink_to_fit();
		HIST.x_max.clear();
		HIST.x_max.shrink_to_fit();
		HIST.base_line_sigma.clear();
		HIST.base_line_sigma.shrink_to_fit();
		HIST.base_line_avr.clear();
		HIST.base_line_avr.shrink_to_fit();
		COLIB.total_area.clear();
		COLIB.total_area.shrink_to_fit();
	}

void clear_vector_colib()
{
		SHORT_LIVE.area.clear();
		SHORT_LIVE.area.shrink_to_fit();
		SHORT_LIVE.x_max.clear();
		SHORT_LIVE.x_max.shrink_to_fit();
		SHORT_LIVE.y_max.clear();
		SHORT_LIVE.y_max.shrink_to_fit();
}

void clear_vector_area()
{
	HIST.area.clear();
	HIST.area.shrink_to_fit();
	COLIB.area.clear();
	COLIB.area.shrink_to_fit();
}

void out_data_file_hist(string outFile_name, TH1F* hist){
	string adress_file = "C:\\Users\\Mikheev\\Desktop\\code_root\\";
	string outFile_path = adress_file + outFile_name;
	ofstream outFile(outFile_path);
	if(outFile.is_open()){
		cout << "file " << outFile_name << " is open" << endl;
	}
	for (int bin = 1; bin <= hist->GetNbinsX(); bin++){
		double bin_сontent = hist->GetBinContent(bin);
		double bin_сenter  = hist->GetBinCenter(bin);
		outFile << static_cast<int>(bin_сenter) << "\t" << static_cast<int>(bin_сontent) << endl;
	}
	outFile.close();
	cout << "data save in "<< outFile_name << " complete " << endl;
}

//----------------------------------------------------------------------------------------



//Тело программы: основные функции чтения, расчета, отбора событий, выреза данных из сигнала
//----------------------------------------------------------------------------------------

vector <Point> read_data(int ch, int run, string choose_mode, string preHIST = {}) // функция чтения для медленных сигналов
{
	vector <Point> signal;
	signal.reserve((event_stop-event_start+1)*points_per_event);

	int coeff = 1;
	if(ch == 5||ch == 6) coeff = -1;
	if(ch == 14||ch == 15) coeff = -1;
	if (ch >= 32)  coeff = -1;

	FILE *f = NULL;
	string filepath;

	ostringstream _ch;
	_ch << ch;
	string channel =  "__ch_" + _ch.str() + ".dat";

	ostringstream _run;
	_run << run;
	filepath = file_name + "\\run_" + _run.str() + channel;

	f = fopen(filepath.c_str(), "rb");

	if (!f) {
		cout << "Error opening file: " << filepath << endl;
	}

	int start; int stop;
	if (choose_mode == "EVENT")
	{
		fseek(f, 2 * event_start * points_per_event, SEEK_SET);
		start = event_start;
		stop = event_stop;
	}
	if (choose_mode == "HIST")
	{
		start = 0;
		stop = events_per_file-1;
	}

	for (int event = start; event <= stop; event++)
	{
		vector<int16_t> buffer(points_per_event);
		size_t read_count = fread(buffer.data(), 2, points_per_event, f);

		for (int j = 0; j < read_count; j++)
		{
			Point point;
			point.x = j;
			point.y = coeff*buffer[j];
			signal.push_back(point);
		}
	}

	fclose(f);

	if (preHIST != "preHIST")
	cout << filepath << endl;

	return signal;

}

vector<Point> read_data_fast_tau(int run_vector, int ch, vector <Events> eventss) // функция чтения для быстрых сигналов (изветны события которые нужно прочитать)
{
	vector <Point> signal;

	int coeff = 1;
	if(ch == 5||ch == 6) coeff = -1;
	if(ch == 14||ch == 15) coeff = -1;
	if (ch >= 32)  coeff = -1;

	FILE *f = NULL;
	string filepath;

	ostringstream _ch;
	_ch << ch;
	string channel =  "__ch_" + _ch.str() + ".dat";

	//cout << "i = " << run_vector << endl;
	ostringstream _run;
	_run << eventss[run_vector].runs;
	filepath = file_name + "\\run_" + _run.str() + channel;

	f = fopen(filepath.c_str(), "rb");

	if (!f)
	{
		cout << "Error opening file: " << filepath << endl;
	}

	vector<int16_t> buffer(points_per_event);
	fseek(f, 2 * eventss[run_vector].events * points_per_event, SEEK_SET);
	size_t read_count = fread(buffer.data(), 2, points_per_event, f);

	for (int j = 0; j < read_count; j++)
	{
		Point point;
		point.x = j;
		point.y = coeff*buffer[j];
		signal.push_back(point);
	}

	fclose(f);

	return signal;
}

vector<vector<Point>> partition_data(vector <Point> Signal, string choose_mode)		//функция делит дату на events and runs
{
	vector <vector<Point>> signal_events;

	int size_vector;

	if (choose_mode == "EVENT")
	{
		size_vector = event_stop - event_start + 1;

		for (int event = 0; event < size_vector; event++)
		{
			vector <Point> temp(Signal.begin() + event*points_per_event, Signal.begin() + (event+1)*points_per_event);
			signal_events.push_back(std::move(temp));
		}
	}
	if (choose_mode == "HIST")
	{
		size_vector = events_per_file;

		for (int event = 0; event < size_vector; event++)
		{
			vector <Point> temp(Signal.begin() + event*points_per_event, Signal.begin() + (event+1)*points_per_event);
			signal_events.push_back(std::move(temp));
		}
	}

	if (choose_mode == "SELECTION")
	{
		signal_events.push_back(Signal);
	}

	return signal_events;
}

//#vlad убрать
bool check_floating_baseLine(const vector<Point>& checked_normal_baseline, const int& period)
{

	int trigger_level = TRIG.slow;

	for (int point = 0; point < period; point++)
	{
		if (abs(checked_normal_baseline[point].y) > trigger_level)
		{
			return true;
		}
	}
	return false;
}

vector <vector<Point>> normalize_baseLine(vector <vector<Point>> Signal, string choose_mode, const int& period) 	//функция нормормировки базовой линии
{
	vector <vector<Point>> normalize_signal = Signal;
	double base_line_avr = 0;
	double base_line_sigma = 0;
	int sum = 0;
	int sum_sigma = 0;

	for (int event = 0; event < Signal.size(); event++)
	{
		sum = 0;
		sum_sigma = 0;

		for (int points = 0; points < period; points++)
		{
			sum += normalize_signal[event][points].y;
		}

		base_line_avr = sum/period;

		for (int points = 0; points < period; points++)
		{
			sum_sigma += pow(normalize_signal[event][points].y - base_line_avr, 2);
		}

		base_line_sigma = sqrt(sum_sigma/period);

		for (int points = 0; points < points_per_event; points++)
		{
			normalize_signal[event][points].y -= base_line_avr;
		}

		// if (choose_main_read != NULL)
		// {
			// if (check_floating_baseLine(normalize_signal[event], period)) continue;
		// }

		//#vlad добавить СЮДА функцию вычисления max y

		if (choose_mode == "HIST")
		{
			//#vlad убрать
			get_base_line_data(base_line_avr, base_line_sigma);

			//#vlad добавть ... HIST.base_line_avr.push_back(base_line_avr);
		}
	}

	return normalize_signal;
}


double integrate_area_in_diapozone(int coeff_separation, vector<vector<Point>>& signal, int threshold, const int& event)	// функция интегрирования сигнала для калибровок по методу ворот
{
	int coeff = 1;
	double integrated_area = 0;

	for (int i = left_line_integrate/sec_per_point; i < right_line_integrate/sec_per_point; i++)
	{
		if (signal[event][i-1].y < threshold && signal[event][i].y > threshold)
		{
			coeff = coeff_separation;
			break;
		}
		else
			coeff = 1;
	}

	for (int j = left_line_integrate/sec_per_point; j <= right_line_integrate/sec_per_point; j++)
	{
		integrated_area += sec_per_point*(signal[event][j].y + signal[event][j-1].y)/2;
	}

	return coeff*integrated_area;
}

bool alfa_ID(int run, int event, vector <Events> alfa_events)
{
	bool status_pass_event = 1;

	for (int i = 0; i < alfa_events.size(); i++)
	{
		if (alfa_events[i].runs == run && alfa_events[i].events == event)
			status_pass_event = 0;
	}
	return  status_pass_event;
}

void Find_peaks(vector <vector<Point>> norm_signal, string choose_mode, int run, string come_again, int threshold, bool alfa_out ={})   // функция нахождения пиков, в этой функции основные действия по сбору и вырезу данных - нужно разделить на несколько
{
	double area = 0;
	double total_area = 0;
	double integrated_area = 0;
	double prevmin = 0;
	double prevmax = 0;

	vector <Point> peak;
	vector <Point> peak_short;
	vector <Point> peak_max;
	vector <vector<Point>> peakls;
	vector <vector<Point>> peakls_max;

	int counter_peak = 0;
	int counter_alpha = 0;

	if (come_again == "SLOW")
	 outFile = fopen(data_path.c_str(), "a+");

	for (int event = 0; event < norm_signal.size(); event++) 	// цикл событий
	{
		bool condition_is_done = 0;

		if (come_again == "SLOW")
			integrated_area = integrate_area_in_diapozone(seperate_coeff, norm_signal, threshold, event);

		int p_min, p_max;
		for (int p = 1; p < norm_signal[0].size(); p++)		// цикл точек в событии
		{
			Point max = {0,0};

		  //	peak_fineder(norm_signal, choose_mode, threshold, peak, peak_max, area, total_area, max, event, p);

			if (norm_signal[event][p-1].y < threshold && norm_signal[event][p].y > threshold)			//часть поиска пиков, записи данных в структуры
			{
				p_min = p;

				while ( norm_signal[event][p_min].y > threshold/6 && p_min > 0) p_min--;

				p_max = p;

				while ( norm_signal[event][p_max].y > threshold/5 && p_max < points_per_event) p_max++;

				if (p_min != prevmin && p_max != prevmax)
				{
					for (int i = p_min; i <= p_max; i++)
					{
						area += sec_per_point*(norm_signal[event][i].y + norm_signal[event][i-1].y)/2;

						if (max.y < norm_signal[event][i].y)
							max = norm_signal[event][i];

						if (choose_mode == "EVENT")
							peak.push_back(norm_signal[event][i]);

						//#vlad зачем это?
						peak_short.push_back(norm_signal[event][i]);
					}
				}
				else
					area = 0;

				if (choose_mode == "EVENT")
					peak_max.push_back(max);

				if (choose_mode == "HIST" && area > 0)
					total_area += area;

				if (choose_mode == "HIST" && area > 0)
				{
					//#vlad убрать cut в другую функцию
					if (/*(max.y > threshold+1) &&  */((area > one_peak_l && area < one_peak_r) || !alfa))
						get_area_max(area, sec_per_point * max.x, max.y);
				}


				prevmin = p_min;
				prevmax = p_max;
			}

			//#vlad убрать в другую функцию
			if (choose_mode == "preHIST")								//для сбора данных событий с альфа частицей
			{
				if(area > 100E3 && max.x*sec_per_point >= lefts && max.x*sec_per_point <= rights)
				{
					get_runs_events(run, event, PEAK_COUNT_EVENT);
					get_tstart_tend_tmax(p_min*sec_per_point, p_max*sec_per_point, max.x*sec_per_point, area, ALPHA);
				}
			}

			if (choose_mode == "HIST" && come_again == "SLOW")			// отбор событий калиброки по времени и по площади для наложения сигнала
			{
				if (area > 0)
				{
					if (total_area > one_alfa_l && total_area < one_alfa_r || !alfa)
						get_tstart_tend_tmax(p_min*sec_per_point, p_max*sec_per_point, max.x*sec_per_point, area, NO_CUT);

				}

				if( area > 0 /*&& area < 380E3*/ && max.x*sec_per_point >= lefts && max.x*sec_per_point <= rights)
				{
					condition_is_done = 1;

					get_tstart_tend_tmax(p_min*sec_per_point, p_max*sec_per_point, max.x*sec_per_point, area, WITH_PIEDISTAL);




						if(alfa_ID(run, event, ALPHA_EVENT.runs_events))
						{
							get_tstart_tend_tmax(p_min, p_max, max.x, area, DISSPERSION);
							get_tstart_tend_tmax(p_min*sec_per_point, p_max*sec_per_point, max.x*sec_per_point, area, CUT);
							get_runs_events(run, event, MULTY_EVENT);

							if( area > 100E3)
							get_runs_events(run, event, SINGLE_EVENT);
						}
						else{counter_alpha ++;}

				}



			}

				// if ((max.x > (left_line_integrate)/sec_per_point && max.x < (right_line_integrate)/sec_per_point))
				// {
					// if (!colibration_mode || (total_area > 0 && (total_area > colib_r || total_area < colib_l))){

						// get_runs_events(run, event, MULTY_EVENT);

					// }

				// }


			peak_short.clear();

			area = 0;

		} //цикл точек

		if (choose_mode == "HIST" && come_again == "SLOW")
		{
			if (!condition_is_done)
				get_tstart_tend_tmax(20, 20, 20, 0.00, WITH_PIEDISTAL);
		}


		if (choose_mode == "EVENT")		// очистка данных для событий
		{
			peakls.push_back(peak);
			peakls_max.push_back(peak_max);

			get_peaks(peakls, peakls_max);

			peak.clear();
			peak_max.clear();
			peak.shrink_to_fit();
			peak_max.shrink_to_fit();

		}

		if (choose_mode == "HIST" && total_area > 0)//&& total_area > 0
		{
			if (total_area > one_alfa_l && total_area < one_alfa_r || !alfa) 	// Условие отбора событий для БЫСТРЫХ сигналов
			{
				get_total_area(total_area, HIST);

				if (come_again == "SLOW")
				{
					ostringstream output_data;

					output_data << run << "\t" << event << "\n";
					fwrite(output_data.str().c_str(), 1, output_data.str().length(), outFile);

					get_runs_events(run, event, FAST_EVENT);
				}
			}
		}

		if (!colibration_mode || (total_area > 0 && (total_area > colib_r || total_area < colib_l)))  // условие отбора событий для калибровки по площади
		{
			get_total_area(total_area, COLIB);

			if (come_again == "SLOW")
			{
				copy_data_colibr();

				get_integrated_area(integrated_area);

			}

		}

		clear_vector_colib(); 		// очистка калибровочных векторов SHORT_LIVE
		total_area = 0;


	}//цикл событий

	if (come_again == "SLOW" && !(alfa_out == NULL))
	{
		fclose(outFile);

	cout << " counter_alpha = " << counter_alpha << endl;

	}

}

//----------------------------------------------------------------------------------------

// функции отрисовки гисторамм, событий, функции отрисовки для калибровок
//----------------------------------------------------------------------------------------

void Plot_EVENT(vector<vector<Point>> divided_data, int columns, int threshold, TCanvas* canvas, const vector <Events>& run_event = {})  	// функция для отображения выбранных событий
{
	int rows = ceil((float)divided_data.size() / columns);
	canvas->Divide(columns,rows,0.01,0.01); //columns and rows

	int size = divided_data.size();
	if (divided_data.size() >= 10)
		int size = 25;

	for(int i = 0; i < divided_data.size(); i++) //size to divided.data.size()
	{
		canvas->cd(i+1);
		vector<Point> event = divided_data[i];
		double* x_signal = extract_xs(event);
		double* y_signal = extract_ys(event);
		if (event.size() > 0)
		{
			TGraph *signal = new TGraph(event.size(), x_signal, y_signal);
			signal->SetMarkerStyle(10);
			signal->Draw();

			ostringstream title;
			signal->GetXaxis()->SetTitle("Time");
			signal->GetYaxis()->SetTitle("signal level");


			if (run_event.empty())
			{
				title << "Event: " << event_start + i;
				signal->SetTitle(title.str().c_str());

			}
			else
			{
				title << "Run: " << PEAK_COUNT_EVENT.runs_events[i].runs
				<< "  Event: " << PEAK_COUNT_EVENT.runs_events[i].events;
				signal->SetTitle(title.str().c_str());
			}

			TF1 *line = new TF1("line","[0] + [1]*x", x_signal[0], x_signal[event.size()-1]);
			line->SetParameter(0, threshold);
			line->SetParameter(1, 0);
			line->SetNpx(10000);
			line->Draw("L same");
		}

		vector <Point> event_peaks = PEAK.peaks[i];
		if (PEAK.peaks.size() > 0)
		{
			TGraph*signal_peaks = new TGraph(PEAK.peaks[i].size(), extract_xs(event_peaks),  extract_ys(event_peaks));
			signal_peaks->Draw("P same");
			signal_peaks->SetMarkerColor(kOrange+1);
			signal_peaks->SetMarkerStyle(20);
		}

		vector <Point> event_peaks_max = PEAK.peaks_max[i];
		if (event_peaks_max.size() > 0)
		{
			TGraph*signal_peaks_max = new TGraph(PEAK.peaks_max[i].size(),
			extract_xs(event_peaks_max),  extract_ys(event_peaks_max));

			signal_peaks_max->Draw("P same");
			signal_peaks_max->SetMarkerColor(kBlue+2);
			signal_peaks_max->SetMarkerStyle(20);
		}

	}

}

void Plot_HIST(const char* canvas_name)				//функция для отображения основных гистограмм (плошадь, tau, baseLine)
{
	int rangeSignal = 1200E3;
	if (!(strcasecmp(canvas_name, "FAST_signal")))
	{
		rangeSignal = 100E3;
	}
	TH1* h1_base_line_sigma = new TH1F("h1_base_line_sigma", "h1_base_line_sigma", 1000, -50, 200);
	TH1* h1_base_line_avr = new TH1F("h1_base_line_avr", "h1_base_line_avr", 1000, -5000, 5000);
	TH1F* h1_peak_time_weight = new TH1F("h1","h1_avr_weight",1600,0,160000);
	TH1F* h1_peak_area = new TH1F("h1", "h1_peak_area", 550*5, -100000, 1.5E6);
	TH1F* h1_total_area = new TH1F("h1", "h1_total_area", 2000*5, -200000, 10E6);
	TH2F* hist_peak_area_ev_vs_evnum = new TH2F("h2", "h2", 110, 0, HIST.total_area.size(), 100, 0, rangeSignal);

	TCanvas *c1 = new TCanvas(canvas_name,canvas_name,1000,1000);
	c1->Divide(3,2,0.01,0.01);

	string str(canvas_name);



	for (int run = 0; run < HIST.total_area.size(); run++)
	{
		hist_peak_area_ev_vs_evnum->Fill(run, HIST.total_area[run]);
	}

	c1->cd(1);
	c1->cd(1)->SetTickx();
	c1->cd(1)->SetTicky();
	hist_peak_area_ev_vs_evnum->GetXaxis()->SetTitle("event");
	hist_peak_area_ev_vs_evnum->GetYaxis()->SetTitle("peak_area_ev [ADC_code*ns]");
	hist_peak_area_ev_vs_evnum->Draw("colz");

	TProfile* prof_hist_peak_area_ev_vs_evnum = hist_peak_area_ev_vs_evnum->ProfileX();
	//prof_hist_peak_area_ev_vs_evnum->Draw("same");
	prof_hist_peak_area_ev_vs_evnum->SetMarkerStyle(20);
	prof_hist_peak_area_ev_vs_evnum->SetMarkerColor(kRed);

	c1->cd(2);
	c1->cd(2)->SetTickx();
	c1->cd(2)->SetTicky();

	if (!colibration_mode)
	{
		for (int i = 0; i < HIST.total_area.size(); i++)
		{
			h1_total_area->Fill(HIST.total_area[i]);
		}

	}
	else
	{
		for (int i = 0; i < COLIB.total_area.size(); i++)
		{
			h1_total_area->Fill(COLIB.total_area[i]);
		}

	}

	h1_total_area->GetXaxis()->SetTitle("Area");
	h1_total_area->GetYaxis()->SetTitle("Number");
	h1_total_area->SetFillColor(kBlue-1.5);
	h1_total_area->Draw();

	c1->cd(3);
	c1->cd(3)->SetTickx();
	c1->cd(3)->SetTicky();

	if (!colibration_mode)
	{
		for (int i = 0; i < HIST.area.size(); i++)
		{
			h1_peak_area->Fill(HIST.area[i]);
		}
	}
	else
	{
		for (int i = 0; i < COLIB.area.size(); i++)
		{
			h1_peak_area->Fill(COLIB.area[i]);
		}
	}

	h1_peak_area->GetXaxis()->SetTitle("Area");
	h1_peak_area->GetYaxis()->SetTitle("Number of peaks");
	h1_peak_area->SetFillColor(kBlue-1.5);
	h1_peak_area->Draw();



	// if (!(strcasecmp(canvas_name, "SLOW_signal")))
	// {
		// ofstream outFile("C:\\Users\\Mikheev\\Desktop\\code_root\\phe_Hist.txt");
		// if (outFile.is_open())
		// {
			// cout << "file phe_Hist.txt is open" << endl;

			// for (int bin = 1; bin <= h1_peak_area->GetNbinsX(); bin++)
			// {
				// double binContent = h1_peak_area->GetBinContent(bin);
				// int intBinContent = static_cast<int>(binContent);
				// double binCenter = h1_peak_area->GetBinCenter(bin);
				// outFile << static_cast<int>(binCenter) << "\t" << intBinContent << endl;

			// }
		// }
		// outFile.close();
		// cout << "data save complete " << endl;

	// }


	c1->cd(4);
	c1->cd(4)->SetTickx();
	c1->cd(4)->SetTicky();

	if (!colibration_mode)
	{
		for (int i = 0; i < HIST.area.size(); i++)
		{
			h1_peak_time_weight->Fill(HIST.x_max[i], HIST.area[i]);
		}
	}
	else
	{
		for (int i = 0; i < COLIB.area.size(); i++) // как сделать ?
		{
			h1_peak_time_weight->Fill(COLIB.x_max[i], COLIB.area[i]);
		}

	}
	
	
	h1_peak_time_weight->SetFillStyle(0);
	h1_peak_time_weight->Draw("HIST");
	gPad->SetLogy();
	if (!(strcasecmp(canvas_name, "FAST_signal"))){
		string tau_file = "tau_file.txt";
		out_data_file_hist(tau_file, h1_peak_time_weight);
	}
	
	c1->cd(5);
	c1->cd(5)->SetTickx();
	c1->cd(5)->SetTicky();
	for (int i = 0; i < HIST.base_line_sigma.size(); i++)
	{
		h1_base_line_sigma->Fill(HIST.base_line_sigma[i]);
	}
	h1_base_line_sigma->SetFillColor(kBlue-1.5);
	h1_base_line_sigma->Draw();
	gPad->SetLogy();

	c1->cd(6);
	c1->cd(6)->SetTickx();
	c1->cd(6)->SetTicky();
	for (int i = 0; i < HIST.base_line_avr.size(); i++)
	{
		h1_base_line_avr->Fill(HIST.base_line_avr[i]);
	}
	h1_base_line_avr->Draw();
	gPad->SetLogy();

}

void Plot_multy_signals(vector <vector<Point>> normalize_signal_multy, TCanvas* canvas, int left_line_integrate, int left_right_integrate)		// функция для отображения наложения сигнала
{
	canvas->cd(3);
	canvas->cd(3)->SetTickx();
	canvas->cd(3)->SetTicky();

	int roof_signal = 2000;
	vector <Point> events;

	int start_events_runs = 0;
	int end_events_runs = 0;
	int number_of_signals = 20;


	if (normalize_signal_multy.size() >= number_of_signals)
	end_events_runs = number_of_signals;
	else
	end_events_runs = normalize_signal_multy.size();


	for (int i = start_events_runs; i <= end_events_runs ; i++) /*_normalize_signal_multy.size()/points_per_event*/
	{
		events = normalize_signal_multy[i];
		TGraph *signal = new TGraph(events.size(), extract_xs(events), extract_ys(events));

		// string number_events = to_string(end_events_runs-start_events_runs + 1);
		// string title = "number_of_Signals: " + number_events;
		// signal->SetTitle(title.c_str());
		// signal->SetMarkerStyle(15);


		int color = TColor::GetColor(50 + (i * 50) % 200, 50 + (i * 100) % 200, 50 + (i * 150) % 200); // Генерация цвета

		signal->SetLineColor(color);
		signal->SetMarkerColor(color);

		if (i == start_events_runs) signal->Draw();
		if (i > start_events_runs) signal->Draw("SAME");

		signal->GetXaxis()->SetTitle("Time");
		signal->GetYaxis()->SetTitle("signal level");
		signal->GetYaxis()->SetRangeUser(0, roof_signal); // Устанавливаем границы по оси Y
		signal->GetXaxis()->SetRangeUser(31E3, 34E3); // Устанавливаем границы по оси Y

	}

	TLine *line_lefty = new TLine(left_line_integrate, 0, left_line_integrate, roof_signal);
	line_lefty->SetLineColor(kRed);
	line_lefty->SetLineStyle(2);
	line_lefty->SetLineWidth(2);

	TLine *line_rigty = new TLine(right_line_integrate, 0, right_line_integrate, roof_signal);
	line_rigty->SetLineColor(kRed);
	line_rigty->SetLineStyle(2);
	line_rigty->SetLineWidth(2);


}


void Multy_Signals(TCanvas* canvas, int ch)		// Сборочная функция по отображению сигнала
{
	vector <Point> multy_signal;
	vector <vector<Point>> partition_signal;
	vector <vector<Point>> normalize_signal;
	vector <vector<Point>> colect_multy_signal;

	//cout << "HIST.runs_events_multy = " << HIST.runs_events_multy.size() << endl;

	for (int i = 0; i < MULTY_EVENT.runs_events.size(); i++)
	{
		multy_signal = read_data_fast_tau(i, ch, MULTY_EVENT.runs_events);
		partition_signal = partition_data(multy_signal, "SELECTION");
		normalize_signal = normalize_baseLine(partition_signal, "NO_GIST", period_stBaseLine);
		colect_multy_signal.push_back(normalize_signal[0]);
	}
	Plot_multy_signals(colect_multy_signal, canvas, left_line_integrate, right_line_integrate);

}

void Plot_area_Amp(TCanvas* canvas)				// функция отрисовки гистограммы зависимости амплитуды сигнала от площади пика
{

	if (!colibration_mode)
	{

		TH2F* AMP_AREA = new TH2F("h2", "AMP_AREA", 100, 0, 150000, 100, 0, 700);

		canvas->cd(1);
		canvas->cd(1)->SetTickx();
		canvas->cd(1)->SetTicky();

		for (int i = 0; i < HIST.area.size(); i++)
		{
			AMP_AREA->Fill( HIST.area[i], HIST.y_max[i]);
		}
		AMP_AREA->GetXaxis()->SetTitle("area_of_phe");
		AMP_AREA->GetYaxis()->SetTitle("Amplitude_of_peak");
		AMP_AREA->Draw("COLZ");

	}
	else
	{
		TH2F* AMP_AREA = new TH2F("h2", "AMP_AREA", 100, 0, 150000 , 100, 0, 700);

		canvas->cd(1);
		canvas->cd(1)->SetTickx();
		canvas->cd(1)->SetTicky();

		for (int i = 0; i < COLIB.area.size(); i++)
		{
			AMP_AREA->Fill(COLIB.area[i], COLIB.y_max[i]);
		}

		AMP_AREA->GetXaxis()->SetTitle("area_of_phe");
		AMP_AREA->GetYaxis()->SetTitle("Amplitude_of_peak");
		AMP_AREA->Draw("COLZ");

	}

}


void Plot_intagrate_Area(TCanvas* canvas, int seperate_coeff)		// функция для отрисоввки гистораммы площади пика калибровка по воротам с методом пьедистала
{
	TH1F* h1_hist_inegrate_area = new TH1F("h3", "h1_hist_inegrate_area", 4088, -1000000, 7E6);

	canvas->cd(2);
	for (int i = 0; i < HIST.integrated_area.size(); i++)
	{
		h1_hist_inegrate_area->Fill(HIST.integrated_area[i]);
	}

	int seperate_line = 22000;																		//
	int lowbound = 0;
	int upbound = h1_hist_inegrate_area->FindBin(seperate_line);

	for (int i = 1; i <= h1_hist_inegrate_area->GetNbinsX(); ++i) {
        if (h1_hist_inegrate_area->GetBinContent(i) > 0) {
            lowbound = i;
            break;
        }
    }

	double Integral_all = h1_hist_inegrate_area->Integral();
	double Integral_piedistal = h1_hist_inegrate_area->Integral(lowbound, upbound);
	double Probability_of_Piedistal = Integral_piedistal/Integral_all;
	double mean_S = 0;
	double S_piedistal = 0;
	double S_signal = 0;
	double seperate_S = 0;
	double mean_phe = 0;
	double G = 0;
	double seperete_G = 0;
	int total_Events = 0;
	int total_Events_signal = 0;

	for (int i = h1_hist_inegrate_area->FindBin(-20000); i < h1_hist_inegrate_area->FindBin(seperate_line); i++)
	{
		S_piedistal += h1_hist_inegrate_area->GetBinContent(i) * h1_hist_inegrate_area->GetBinCenter(i);
		total_Events += h1_hist_inegrate_area->GetBinContent(i);
	}

	for (int i = h1_hist_inegrate_area->FindBin(seperate_line); i < h1_hist_inegrate_area->GetNbinsX(); i++)
	{
		S_signal += h1_hist_inegrate_area->GetBinContent(i) * h1_hist_inegrate_area->GetBinCenter(i);
		total_Events_signal += h1_hist_inegrate_area->GetBinContent(i);
	}

	S_piedistal = S_piedistal/total_Events;
	S_signal = S_signal/total_Events_signal;
	S_signal = S_signal/seperate_coeff;

	seperate_S = (total_Events*S_piedistal + total_Events_signal*S_signal)/(total_Events + total_Events_signal);
	mean_S = h1_hist_inegrate_area-> GetMean();
	mean_phe = -log(Probability_of_Piedistal);

	G = mean_S/mean_phe;
	seperete_G = seperate_S/mean_phe;

	// cout << "Integral_all = " << Integral_all << endl;
	// cout << "Integral_piedistal = " << Integral_piedistal << endl;
	// cout << "Probability_of_Piedistal = " << Probability_of_Piedistal << endl;
	// cout << "mean_S = " << mean_S <<endl;
	// cout << "piedistal_S = " << S_piedistal <<endl;
	// cout << "signal_S = " << S_signal <<endl;
	 cout << "seperate_S = " << seperate_S <<endl;
	 cout << "seperete_G = " << seperete_G <<endl;
	 //cout << "G = " << G <<endl;
	 cout << "mean_phe = " << mean_phe <<endl;


	TLatex latex;
	latex.SetTextSize(1);
	latex.DrawLatex(1000, 200, to_string(static_cast<double>(G)).c_str());

	h1_hist_inegrate_area->GetXaxis()->SetTitle("area_of_phe");
	h1_hist_inegrate_area->GetYaxis()->SetTitle("event");
	h1_hist_inegrate_area->Draw();

	TLine *line_seperate = new TLine(seperate_line, 0 , seperate_line, h1_hist_inegrate_area->GetMaximum());
	line_seperate->Draw();
	line_seperate->SetLineColor(kRed);


}

vector<Point> cut_normalize_signal(vector <vector<Point>> normalize_signal, int left_border, int right_border)
{
	vector<Point> cut_norm_signal;

	for (int i = 0; i < normalize_signal[0].size(); i++)
	{
		if (i >= int(left_border/sec_per_point) && i <= int(right_border/sec_per_point))
		{
			cut_norm_signal.push_back(normalize_signal[0][i]);
		}
	}
	return cut_norm_signal;
}

void PLOT_Core_hist(data_bank_CORE& CORELATION, data_bank_TIME& TYPE_CUT, TCanvas* canvas)
{
	canvas-> Divide(3,2,0.01,0.01);

	CORELATION.fillHistograms(TYPE_CUT);

	canvas->cd(1);
	CORELATION.cor_start_vs_area->Draw("COLZ");
	canvas->cd(2);
	CORELATION.cor_max_vs_area->Draw("COLZ");
	canvas->cd(3);
	CORELATION.cor_max_vs_start->Draw("COLZ");
	canvas->cd(4);
	CORELATION.cor_start_end_vs_area->Draw("COLZ");
	canvas->cd(5);
	CORELATION.cor_max_start_vs_area->Draw("COLZ");

	CORELATION.cor_start_vs_area->GetXaxis()->SetRangeUser(31000,33000);
	CORELATION.cor_max_vs_area->GetXaxis()->SetRangeUser(30000,35000);

}

void Find_G_for_histogram(vector<TH1F*> cascade_hist) // функция вычисления G для калибровки по методу cut сигнала по времени
{
	double entries;
	double mean;
	double S;
	double G;
	double piedistal_probability;
	int Events;

   	for (int i = 0; i < cascade_hist[5]->GetNbinsX(); i++)  //cascade_hist[5] - гистграмма area пиков медленных сигналов в характерном времени
	{
		S += cascade_hist[5]->GetBinContent(i) * cascade_hist[5]->GetBinCenter(i);
	}

	entries = cascade_hist[5]->GetEntries();				// число зарегестрированных событий, то есть события вне пьедестала
	Events = (run_stop - run_start + 1)*events_per_file;	// число всех событий
	S = S/Events;
	piedistal_probability = 1 - (entries/Events);			// вероятность пьедестала
	mean = -log(piedistal_probability);
	G = S/mean;

	cout << "Events = " << Events << "  piedistal_probability = " << piedistal_probability << endl << "mean = " << mean << endl << "S = " << S << endl << "G = " << G << endl;

}

vector<TH1F*> Fill_histograms(data_bank_TIME& TYPE_CUT)  		// функция записи данных для калибровок по времени
{
	vector <TH1F*> cascade_hist;

	for (int i = 0; i < TYPE_CUT.t_start.size(); i++)
	{
		TYPE_CUT.hist_tau_start->Fill(TYPE_CUT.t_start[i]);
		TYPE_CUT.hist_tau_max->Fill(TYPE_CUT.max[i]);
		TYPE_CUT.hist_start_end->Fill(TYPE_CUT.t_end[i] - TYPE_CUT.t_start[i]);
		TYPE_CUT.hist_max_start->Fill(TYPE_CUT.max[i] - TYPE_CUT.t_start[i]);
		TYPE_CUT.hist_max_end->Fill(TYPE_CUT.t_end[i] - TYPE_CUT.max[i]);
		TYPE_CUT.hist_area_in_region->Fill(TYPE_CUT.area_in_rigeon[i]);
	}

	cascade_hist.push_back(TYPE_CUT.hist_tau_start);
	cascade_hist.push_back(TYPE_CUT.hist_tau_max);
	cascade_hist.push_back(TYPE_CUT.hist_start_end);
	cascade_hist.push_back(TYPE_CUT.hist_max_start);
	cascade_hist.push_back(TYPE_CUT.hist_max_end);
	cascade_hist.push_back(TYPE_CUT.hist_area_in_region);


	return cascade_hist;
}



void Set_Histograms_single(TCanvas* canvas, data_bank_TIME& TYPE_CUT, vector<TH1F*> cascade_hist) 	// функция установки пвраметров рисовки гистограмм
{
	canvas-> Divide(3,2,0.01,0.01);

	string canvasName = canvas->GetName();


	if (canvasName == "ccut")
	{
		Find_G_for_histogram(cascade_hist);
	}


    for (int i = 0; i < cascade_hist.size(); i++)
	{
		canvas->cd(i+1);
		cascade_hist[i]->Draw();

		if (canvasName == "ccut")
		{
			cascade_hist[i]->SetLineColor(kPink - 0.1);
			cascade_hist[i]->SetLineWidth(2);
		}

		if(canvasName == "cnocut")
		{
			cascade_hist[i]->SetLineColor(kGreen - 0.2);
			cascade_hist[i]->SetLineWidth(2);
		}

		if(canvasName == "piedistal")
		{
			cascade_hist[i]->SetLineColor(kBlack + 0.1);
			cascade_hist[i]->SetLineWidth(2);
		}
	}


	TYPE_CUT.hist_max_start->GetXaxis()->SetRangeUser(0,500);
	TYPE_CUT.hist_tau_start->GetXaxis()->SetRangeUser(31000,33000);
	TYPE_CUT.hist_tau_max->GetXaxis()->SetRangeUser(31000,36000);
}



void Plot_time_signals()		//Сброрчная функция для отрисовки временных гистограмм сигнала
{
	TCanvas* ccut = new TCanvas("ccut", "ccut", 1000, 10000);
	TCanvas* cnocut = new TCanvas("cnocut", "cnocut", 1000, 10000);
	TCanvas* piedistal = new TCanvas("piedistal", "piedistal", 1000, 10000);

	Set_Histograms_single(ccut, CUT, Fill_histograms(CUT));
	Set_Histograms_single(cnocut, NO_CUT, Fill_histograms(NO_CUT));
	Set_Histograms_single(piedistal, WITH_PIEDISTAL, Fill_histograms(WITH_PIEDISTAL));

}

double min(vector <double> data)
	{
		double min = data[0];
		for (int i = 0; i < data.size(); i++)
		{
			if(data[i] < min)
				min = data[i];
		}
		//cout << "min = " <<  min << endl;
		return min;
	}

	double max(vector <double> data)
	{
		double max = data[0];
		for (int i = 0; i < data.size(); i++)
		{
			if(data[i] > max)
				max = data[i];
		}

		//cout << "max = " <<  max << endl;
		return max;
	}

void Plot_3D_corelation(vector <double> disspersion, vector <double> time, vector <double> area, TCanvas* canvas)
{
	canvas->Divide(2, 1);


	TH3F* sigma_vs_time_vs_area = new TH3F("sigma_vs_time_vs_area", "sigma_vs_time_vs_area",area.size()/100, min(area), max(area),
	time.size()/100, min(time), max(time), disspersion.size()/100, min(disspersion), max(disspersion));
	TH2F* sigma_vs_area = new TH2F("sigma_vs_area", "sigma_vs_area", 20, min(area), max(area),20, min(disspersion), max(disspersion));


	for (int i = 0; i < disspersion.size(); i++)
	{
		sigma_vs_area->Fill(area[i], disspersion[i]);
		sigma_vs_time_vs_area->Fill(area[i], time[i], disspersion[i]);
	}

	canvas->cd(1);
	sigma_vs_time_vs_area->Draw("COLZ");

	canvas->cd(2);
	sigma_vs_area->Draw("COLZ");

}

double get_disspersion(vector<Point> signal_peak, double area, int time_start, int time_end)
{
	int N_points = 0;
	double average_y = 0;
	double disspersion = 0;
	double cv;

	for(int i = signal_peak[time_start].x; i <= signal_peak[time_end].x; i++)
	{
		average_y += signal_peak[i].y;
		N_points++;
	}

	for (int i = signal_peak[time_start].x; i <= signal_peak[time_end].x; i++)
	{
		disspersion += (average_y - signal_peak[i].y)*(average_y - signal_peak[i].y);
	}

	disspersion = sqrt(disspersion/N_points);

	cv = disspersion/average_y;

	return cv;
}


//----------------------------------------------------------------------------------------

// Основные сборочные функции, функции которые содержат подфункции
//----------------------------------------------------------------------------------------


void Plot_slow_vs_fast_signals(vector<vector<vector<Point>>> signal_multy, TCanvas* canvas)
{
	int column = columns(evt_slow_vs_fast_stop-evt_slow_vs_fast_start+1);
	int rows = ceil((float)(evt_slow_vs_fast_stop-evt_slow_vs_fast_start+1)/column);

	canvas->Divide(column,rows,0.01,0.01); //columns and rows
	cout << "yes2" << endl;

	for (int i = 0; i < column*rows/*signal_multy.size()*/; i++){

		canvas->cd(i+1);

		TGraph *signal = new TGraph(signal_multy[i][0].size(), extract_xs(signal_multy[i][0]), extract_ys(signal_multy[i][0]));

		ostringstream title;

		title << "event: " << evt_slow_vs_fast_start + i;
		signal->SetTitle(title.str().c_str());

		signal->SetLineColor(kBlue);
		signal->SetMarkerColor(kBlue);
		//signal->SetMarkerStyle(20);
		//signal->SetMarkerSize(0.4);
		signal->Draw();

		TGraph *signal1 = new TGraph(signal_multy[i][1].size(), extract_xs(signal_multy[i][1]), extract_ys(signal_multy[i][1]));

		signal1->SetLineColor(kRed);
		signal1->SetMarkerColor(kRed);

		int coefficient = 10.0;

		for (int i = 0; i < signal1->GetN(); ++i)
		{
			double xValue, yValue;
			signal1->GetPoint(i, xValue, yValue);
			yValue *= coefficient;
			signal1->SetPoint(i, xValue, yValue);
		}

		signal1->Draw("SAME");

	}


}

void Bounded_slow_fast_signals(vector<vector<Events>> runs_events, int ch_slow, int ch_fast)
{
	vector<vector<Point>> colect_slow_fast_signal;
	vector<vector<vector<Point>>> collect;

	TCanvas* slow_fast_sig_canvas = new TCanvas("slow_fast_sig_canvas", "slow_fast_sig_canvas", 1000, 1000);

	vector <Point> Signal;
	vector<vector<Point>> partition_signal;
	vector<vector<Point>> normalize_signal;

	for (int i = 0; i < runs_events.size(); i++){

		Signal = read_data_fast_tau(0, ch_slow, runs_events[i]);
		partition_signal = partition_data(Signal, "SELECTION");
		normalize_signal = normalize_baseLine(partition_signal, "EVENT", period_stBaseLine);
		colect_slow_fast_signal.push_back(normalize_signal[0]);

		Signal = read_data_fast_tau(0, ch_fast, runs_events[i]);
		partition_signal = partition_data(Signal, "SELECTION");
		normalize_signal = normalize_baseLine(partition_signal, "EVENT", period_stBaseLine);
		colect_slow_fast_signal.push_back(normalize_signal[0]);
		collect.push_back(colect_slow_fast_signal);
		colect_slow_fast_signal.clear();

	}

		Plot_slow_vs_fast_signals(collect, slow_fast_sig_canvas);

}

void Other_PARAMETERS(int ch_slow, int ch_fast)		// Сборочная функция вызова функций отрисовки гистограмм калировок, сигналов наложения, временных гистограмм
{
	//TCanvas *c4 = new TCanvas("OTHER PARAMETERS","OTHER PARAMETERS",1000,1000);
	//c4 ->Divide(3,2,0.01,0.01);

	//Plot_area_Amp(c4);
	//Plot_intagrate_Area(c4, seperate_coeff);
	//Multy_Signals(c4, ch);
	Plot_time_signals();

	vector<vector<Events>> slow_vs_fast_sig;
	Events run_event;

	// for (int i = evt_slow_vs_fast_start; i <= evt_slow_vs_fast_stop; i++)
	// {
		// run_event.runs = run_slow_vs_fast;
		// run_event.events = i;

		// vector <Events> r_ev;
		// r_ev.push_back(run_event);
		// slow_vs_fast_sig.push_back(r_ev);
	// }

	for (int i = 0; i < SINGLE_EVENT.runs_events.size(); i++)
	{
		run_event.runs = SINGLE_EVENT.runs_events[i].runs;
		run_event.events = SINGLE_EVENT.runs_events[i].events;

		vector <Events> r_ev;
		r_ev.push_back(run_event);
		slow_vs_fast_sig.push_back(r_ev);
	}

	Bounded_slow_fast_signals(slow_vs_fast_sig, ch_slow, ch_fast);
	TCanvas* core_time_area = new TCanvas("core_time_area", "core_time_area", 1000, 10000);
	//PLOT_Core_hist(CORELATION, CUT, core_time_area);
	PLOT_Core_hist(CORELATION, NO_CUT, core_time_area);

}



void pre_finder_peaks(string HIST, int ch)
{
	vector <Point> Signal;
	vector<vector<Point>> partition_signal;
	vector<vector<Point>> normalize_signal;

	RunMode Run = conditions_mode(HIST);
	for (int run = Run.start; run <= Run.stop; run++)
	{
		Signal = read_data(ch, run, "HIST", "preHIST");
		partition_signal = partition_data(Signal, "HIST");
		normalize_signal = normalize_baseLine(partition_signal, "noHIST", period_stBaseLine);
		Find_peaks(normalize_signal, "preHIST", run, "noSLOW", threshold_slow);
	}
}

void func_HIST(string HIST, int ch, bool alfaout) 		// Сборочая функция для основных гистограмм МЕДЛЕННОГО сигнала (tau , total area, peak area .. )
{
	colibration("SLOW");
	alfa_mode("SLOW");
	vector <Point> Signal;
	vector<vector<Point>> partition_signal;
	vector<vector<Point>> normalize_signal;
	vector<vector<Point>> cutty_norm_signal;

	RunMode Run = conditions_mode(HIST);
	for (int run = Run.start; run <= Run.stop; run++)
	{
		Signal = read_data(ch, run, HIST);
		partition_signal = partition_data(Signal, HIST);
		normalize_signal = normalize_baseLine(partition_signal, "HIST", period_stBaseLine);
		cutty_norm_signal.push_back(cut_normalize_signal(normalize_signal, left_line_integrate, right_line_integrate));
		Find_peaks(normalize_signal, HIST, run, "SLOW", threshold_slow, alfaout);
		//Find_peaks(cutty_norm_signal, HIST, run, "0", threshold_slow);
	}

	Plot_HIST("SLOW_signal");

	clear_vector();
}

void func_EVENT(string EVENT, int ch)			// Сборочная функция для отображения сигнала по событиям
{
	TCanvas *multy_event = new TCanvas("event","event",1000,1000);

	int threshold_temp = threshold_slow;
	if (ch == 5 || ch == 6) threshold_temp = threshold_fast;
	if (ch == 2 || ch == 3 || ch == 1) threshold_temp = threshold_slow;
	if (ch == 14 || ch == 15) threshold_temp = threshold_fast;

	colibration("FAST");
	alfa_mode("FAST");
	vector <Point> Signal;
	vector<vector<Point>> partition_signal;
	vector<vector<Point>> normalize_signal;

	RunMode Run = conditions_mode(EVENT);
	for (int run = Run.start; run <= Run.stop; run++)
	{
		Signal = read_data(ch, run, EVENT);
		partition_signal = partition_data(Signal, EVENT);
		normalize_signal = normalize_baseLine(partition_signal, "EVENT", period_stBaseLine);
		Find_peaks(normalize_signal, EVENT, run, "0", threshold_temp);
		Plot_EVENT(normalize_signal, columns(event_stop-event_start+1), threshold_temp, multy_event);
	}
}


void watch_single_EVENT(string EVENT, int ch, vector <Events> events, TCanvas* canvas)			// Сборочная функция для просмотра сигнала по событиям с условием отбора пиков (по времени )
{
	int threshold_temp = threshold_slow;
	if (ch == 5 || ch == 6) threshold_temp = threshold_fast;
	if (ch == 2 || ch == 3 || ch == 1) threshold_temp = threshold_slow;
	if (ch == 14 || ch == 15) threshold_temp = threshold_fast;

	colibration("FAST");
	alfa_mode("FAST");

	vector <Point> Signal;
	vector<vector<Point>> partition_signal;
	vector<vector<Point>> normalize_signal;
	vector <vector<Point>> collect_signal;

	int max_event_per_list = 20;
	int temp_size_single_event;

	if (events.size() >= max_event_per_list)
		temp_size_single_event = max_event_per_list;
	else
		temp_size_single_event = events.size();

	for (int i = 0; i < temp_size_single_event; i++)
	{
		Signal = read_data_fast_tau(i, ch, events);
		partition_signal = partition_data(Signal, "SELECTION");
		normalize_signal = normalize_baseLine(partition_signal, "EVENT", period_stBaseLine);
		collect_signal.push_back(normalize_signal[0]);
	}

		Find_peaks(collect_signal, EVENT, 0, "0", threshold_temp);
		Plot_EVENT(collect_signal, columns(events.size()), threshold_temp, canvas, events);
}



void watch_disperssion_peaks(string EVENT, int ch)			// Сборочная функция для просмотра диссперсии пика во временных воротах
{
	TCanvas *disspersion_vs_time_peaks_area_peaks = new TCanvas("disspersion_vs_time_peaks_area_peaks","disspersion_vs_time_peaks_area_peaks",1000,1000);

	int threshold_temp = threshold_slow;
	if (ch == 5 || ch == 6) threshold_temp = threshold_fast;
	if (ch == 2 || ch == 3 || ch == 1) threshold_temp = threshold_slow;
	if (ch == 14 || ch == 15) threshold_temp = threshold_fast;

	colibration("FAST");
	alfa_mode("FAST");

	vector <Point> Signal;
	vector <double> disspersion_peak;
	vector <double> max_peak;
	vector <double> area_peak;

	vector<vector<Point>> partition_signal;
	vector<vector<Point>> normalize_signal;

	for (int i = 0; i < SINGLE_EVENT.runs_events.size(); i++)
	{
		Signal = read_data_fast_tau(i, ch, SINGLE_EVENT.runs_events);
		partition_signal = partition_data(Signal, "SELECTION");
		normalize_signal = normalize_baseLine(partition_signal, "EVENT", period_stBaseLine);
		disspersion_peak.push_back(get_disspersion(normalize_signal[0], DISSPERSION.area_in_rigeon[i], DISSPERSION.t_start[i], DISSPERSION.t_end[i]));
		max_peak.push_back(DISSPERSION.max[i]*sec_per_point);
		area_peak.push_back(DISSPERSION.area_in_rigeon[i]);

	}

	Plot_3D_corelation(disspersion_peak, max_peak, area_peak, disspersion_vs_time_peaks_area_peaks);
}


Events finder_alpha(vector<Point> normal_signal, int time_start, int time_end, Events Run_event)
{
	Events alfa_runs_events;
	int peak_count = 0;
	int shift_time = 300;
	for (int i = (time_start - shift_time)/sec_per_point ; i <= (time_end + shift_time)/sec_per_point; i++)
	{
		if (normal_signal[i-1].y < 10 && normal_signal[i].y > 10)
		//if (abs(double(normal_signal[i].y / normal_signal[i-1].y)) >= threshold_fast)
			peak_count ++;
	}
	NUMBER_PEAKS.peaks_number.push_back(peak_count);

	if (peak_count > count_peak)
	alfa_runs_events = Run_event;
	else
	alfa_runs_events = {-1, -1};


	return alfa_runs_events;
}

void number_peak_spectrum()
{
	TCanvas* hist_peak_number = new TCanvas("hist_peak_number", "hist_peak_number", 1000, 1000);
	TH1F* peak_spectrum = new TH1F("peak_spectrum", "peak_spectrum", 20, 0, 20);

	for (int i = 0; i < NUMBER_PEAKS.peaks_number.size(); i ++)
	{
		peak_spectrum->Fill(NUMBER_PEAKS.peaks_number[i]);
	}

	for (int bin = 1; bin < peak_spectrum->GetNbinsX(); bin++)
	{
		double variation = peak_spectrum->GetBinContent(bin)/NUMBER_PEAKS.peaks_number.size();
		double binContent = peak_spectrum->GetBinContent(bin);
		double binCenter = peak_spectrum->GetBinCenter(bin);
		variation = variation *100;
		string var = to_string(variation) + " % ";
		cout << bin << " varioation = " << var << endl;
		TLatex latex;
		latex.SetTextSize(0.03);
		latex.DrawLatex(binCenter, binContent, (var).c_str());
	}

	hist_peak_number->cd(1);
	peak_spectrum->Draw("COLZ");

}

void search_alfa_fast_tau(vector <Events> peaks_in_time, data_bank_TIME TYPE, int ch, int ch_watch)
{
	TCanvas* canvas_alfa = new TCanvas("canvas_alfa", "canvas_alfa", 1000, 1000);

	vector <Point> signal;
	vector <vector<Point>> partition_signal;
	vector <vector<Point>> normalize_signal;
	vector <Events> alfa_events;

	for (int i = 0; i < peaks_in_time.size(); i++)
	{
		signal = read_data_fast_tau(i, ch, peaks_in_time);
		partition_signal = partition_data(signal, "SELECTION");
		normalize_signal = normalize_baseLine(partition_signal, " noHIST", period_stBaseLine);
		alfa_events.push_back(finder_alpha(normalize_signal[0], TYPE.t_start[i], TYPE.t_end[i], peaks_in_time[i]));
	}

	vector <Events> nonZero_event;
	for (int i = 0; i < alfa_events.size(); i++)
	{
		if (!(alfa_events[i].runs == -1 && alfa_events[i].events == -1))
			nonZero_event.push_back(alfa_events[i]);
	}
	alfa_events = nonZero_event;


	watch_single_EVENT("EVENT", ch_watch, alfa_events, canvas_alfa);

	ALPHA_EVENT.runs_events = alfa_events;
}


void func_HIST_fast_tau(int ch)			// Сборочная функция для гистограмм БЫСТРОГО сигнала
{
	colibration("FAST");
	alfa_mode("FAST");
	vector<Point> signal;
	vector <vector<Point>> partition_signal;
	vector <vector<Point>> normalize_signal;


	for (int i = 0; i < FAST_EVENT.runs_events.size(); i++)
	{
		signal = read_data_fast_tau(i, ch, FAST_EVENT.runs_events);
		partition_signal = partition_data(signal, "SELECTION");
		normalize_signal = normalize_baseLine(partition_signal, "HIST", period_stBaseLine);
		Find_peaks(normalize_signal, "HIST", 0, "0", threshold_fast);
	}
	Plot_HIST("FAST_signal");

	clear_vector();
}


void Calibr_func(int slow, int fast, bool alfa_out)
{
	pre_finder_peaks("HIST", slow);
	search_alfa_fast_tau(PEAK_COUNT_EVENT.runs_events, ALPHA, fast, slow);
	number_peak_spectrum();
	func_HIST("HIST", slow, alfa_out);
	Other_PARAMETERS(slow, fast);
	clear_vector_area();
	func_EVENT("EVENT", fast);
}

void Tau_func(int slow, int fast, bool alfa_out)
{
	func_HIST("HIST", slow, alfa_out);
	clear_vector_area();
	func_HIST_fast_tau(fast);
}


void calibration_or_tau(int slow, int fast)
{
	string temp_parametr;
	bool alfa_out = NULL;

	cout << endl << "calib / tau ?" << endl << "enter 'c' or 't' " << endl;
	cin >> temp_parametr;
	cout << "mode = " << temp_parametr << endl;

	if (temp_parametr == "c") alfa_out = 1;

	if (temp_parametr == "c") Calibr_func(slow, fast, alfa_out);
	if (temp_parametr == "t") Tau_func(slow, fast, alfa_out);
}

//----------------------------------------------------------------------------------------

int SignalProcessing()
{
	ch.slow = ch_slow;
	ch.fast = ch_fast;

	outFile = fopen(data_path.c_str(), "w");
	fclose(outFile);

	calibration_or_tau(ch.slow, ch.fast);


	//TCanvas* single_event = new TCanvas("single_event","single_event",1000,1000); 	//  функция для вызова просмотра событий отобранных, второй параметр - канал

	//watch_single_EVENT("EVENT", ch.slow , SINGLE_EVENT.runs_events, single_event);

	//watch_single_EVENT("EVENT", 2 , SINGLE_EVENT.runs_events, single_event);		//  функция для вызова просмотра событий отобранных, второй параметр - канал

	//watch_disperssion_peaks("EVENT", 2);

	//func_EVENT("EVENT", ch.slow); 			//  функция для вызова просмотра событий , второй параметр - канал

	 return 0;
}
