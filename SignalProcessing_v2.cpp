 #include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <malloc.h>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <math.h>

#include <TFile.h>
#include <TTree.h>


using namespace std;
namespace fs = std::filesystem; // Just an alias

const size_t DATA_OBJECT_SIZE = 2; // number of bytes per data object
typedef int16_t data_t; // sizeof(DATA_OBJECT_TYPE) == DATA_OBJECT_SIZE

const int RUN_START = 1;
const int RUN_STOP = 1;

const fs::path DATA_DIR = "C:\\code_root\\data_in";
const int EVENTS_PER_RUN = 1000;
const int POINTS_PER_EVENT = 9999;
const int SEC_PER_POINT = 16;
const int PERIOD = 10000;

typedef enum {
    CHANNEL_SLOW,
    CHANNEL_FAST,
} CHANNEL_TYPE;

typedef enum {
	THRESH_SLOW = 40,
	THRESH_FAST = 20,
}THRESHOLD_TYPE;

typedef struct {
    CHANNEL_TYPE t;
    int n;
	double threshhold;
} Channel;

typedef struct {
    int x;
    int y;
} Point;

typedef vector<Point> Event;

struct Events {
    int event;
    int run;

    bool operator<(const Events& other) const {
        return tie(run, event) < tie(other.run, other.event);
    }

};

typedef struct {
	int left;
	int right;
} Borders;


// two trees for each channel which fill with struct Event_data 
// Event_data include struct Peak_data as vector of peaks
struct Baseline_data{
double_t base_avr = 0.0,
	   base_sigma = 0.0, 
	   base_fluct = 0.0;
};
struct Peak_data{
	double_t area;
	Point max;
	Borders bord;
};
struct Event_data{
	///pushing fields by func find peaks
	double_t total_area = 0.0;
	int run = 1;
	vector <Peak_data> peaks;
	///lately pushing baselineinfo
	Baseline_data baseline;
	
	
};

void setup_branches(TTree* tree, Event_data& data){
	tree->Branch("total_area", &data.total_area);
	tree->Branch("run", &data.run);
    tree->Branch("baseline", &data.baseline);
	tree->Branch("peaks", &data.peaks); // its complicated branch for vector data
}
    

/*void creat_tree(){
	TFile *file = new TFile("event.root", "RECREATE");
	TTree *Chan_slow = new TTree("Chan_slow", "CHANNEL SLOW data");
	TTree *Chan_fast = new TTree("Chan_fast", "CHANNEL FAST data");
	// two trees of fastslow signals
	Event_data CHANNEL_SLOW, CHANNEL_FAST;
	setup_branches(Chan_slow, CHANNEL_SLOW);
	setup_branches(Chan_fast, CHANNEL_FAST);
	
	file->Write(); 
    file->Close(); 
}*/


vector<data_t> read_data(fs::path filepath, size_t n_obj) {
    vector<data_t> data(n_obj);

    FILE* f = fopen(filepath.string().c_str(), "rb");
    if (!f) {
        cout << "Error opening file: " << filepath << endl;
        exit(EXIT_FAILURE);
    }
    fread(data.data(), DATA_OBJECT_SIZE, n_obj, f);
    fclose(f);

    return data;
}



// Read data from channel `ch` and run `run`
vector<data_t> read_data_by_run(Channel chan, int run) {
    char filename[32];
    sprintf(filename, "run_%d__ch_%d.dat", run, chan.n);
    fs::path filepath = DATA_DIR / filename;

    return read_data(filepath, POINTS_PER_EVENT * EVENTS_PER_RUN);
}


vector<Event> split_data_to_events(Channel chan, vector<data_t> data) {
    int coeff;
    if (chan.t == CHANNEL_SLOW) coeff = 1;
    else if (chan.t == CHANNEL_FAST) coeff = -1;

    vector<Event> events(EVENTS_PER_RUN);
    for (size_t i = 0; i < EVENTS_PER_RUN; i++) {
        Event event(POINTS_PER_EVENT);
        for (size_t j = 0; j < POINTS_PER_EVENT; j++) {
            event[j] = Point {
                .x = (int)j,
                .y = coeff*data[i * POINTS_PER_EVENT + j],
            };
        }
        events[i] = event;
    }

    return events;
}


Event_data find_peaks(Channel chan, Event data, const int run, Borders left_right = {}){
	Event event = data;
	Borders t_bords;
	double_t area = 0, total_area = 0;
	int temp_left, temp_right;
	int depth = 5;	
	vector <Peak_data> peaks_info;
	if (left_right.left != NULL ) {t_bords.left = 1; t_bords.right = data.size();}
	else {t_bords.left = left_right.left; t_bords.right = left_right.right;}  
	if (t_bords.left < 0 || t_bords.right >event.size()) {
		cout << "out of limits" << endl;
		t_bords.left = 1; t_bords.right = event.size();
	} // conditions for diff borders
	for (size_t i = t_bords.left + 1; i < t_bords.right; i++){
		Point max = {0, 0}; Borders p_bords = {0, 0}; area = 0;
		//conditions of finding peak
		if(event[i-1].y < chan.threshhold && event[i].y > chan.threshhold){
			p_bords.left = i;
			while (event[i].y > chan.threshhold/depth && p_bords.left > 0) p_bords.left--;
			p_bords.right = i;
			while (event[i].y > chan.threshhold/depth && p_bords.right < POINTS_PER_EVENT) p_bords.right++;
			if (p_bords.left != temp_left && p_bords.right != temp_right){
				for (size_t i = p_bords.left; i <= p_bords.right; i++){
					area += SEC_PER_POINT*(event[i-1].y + event[i].y)/2;
					if (max.y < event[i].y)
						max = event[i];
				} //loop of one peak
				peaks_info.push_back(Peak_data{
					.area = area,
					.max = max,
					.bord = p_bords
				});	// get fill vector of peak data	
			}
			else area = 0;
			total_area += area;
			temp_left  = p_bords.left;
			temp_right = p_bords.right;
		}
	}//loop of one event/segment
	Event_data Event_info = (Event_data{
		.total_area = total_area,
		.run = run,
		.peaks = peaks_info 
	}); // get fill of event info data
return Event_info;
}


vector<Long64_t> find_events_by_condition(TTree* tree, const string& cut) {
    vector<Long64_t> selected_events;
    TTreeFormula cutFormula("cut", cut.c_str(), tree);
    
    for (Long64_t i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        if (cutFormula.EvalInstance()) {
            selected_events.push_back(i);
        }
    }
    
    return selected_events;
}// get a formalure and make cut on events this numbers of events will be needed

vector<Event> normalize_baseline(vector<Event> data){
	double_t baseline_avr;
	int start_line = PERIOD/SEC_PER_POINT;
	int  sum; 
	vector <Event> normal_signal = data;
	for (size_t i = 0; i < normal_signal.size(); i++){
		sum = 0;
		for (size_t j = 0; j < start_line; j++){
			sum += normal_signal[i][j].y;
		}
		baseline_avr = sum/start_line;
		for (size_t j = 0; j < POINTS_PER_EVENT; j++){
			normal_signal[i][j].y -= baseline_avr;
		}	
	}
	return normal_signal;
} 

Baseline_data get_baseline_info(Event data){
	double_t baseline_avr, baseline_sigma, maxy_start;
	int start_line = PERIOD/SEC_PER_POINT;
	Event normal_signal = data;
	int sum = 0; int sum_sigma = 0;
	for (size_t i = 0; i < start_line; i++){
		sum += normal_signal[i].y;
	}
	baseline_avr = sum/start_line;	
	for (size_t i = 0; i < POINTS_PER_EVENT; i++){
		normal_signal[i].y -= baseline_avr;
	}
	for (size_t i = 0; i < start_line; i++){
		sum_sigma += pow((normal_signal[i].y - baseline_avr), 2);
		if ( maxy_start > normal_signal[i].y)
		maxy_start = normal_signal[i].y; 	
	}
	baseline_sigma = sqrt(sum_sigma/start_line);

	Baseline_data base = (Baseline_data{
		.base_avr = baseline_avr,
		.base_sigma = baseline_sigma,
		.base_fluct = maxy_start
	});
	return base;
}



//.L Signal_processing.C
 int Signal_processing_v2() {
    Channel chan = {
        .t = CHANNEL_SLOW,
        .n = 2,
		.threshhold = THRESH_SLOW, 
		//for fast channel for ex it would be .n = 5 .treashold = 20; 
    };

	Event_data SLOW;
	TFile *file = new TFile("event.root", "RECREATE");
	TTree *Chan_slow = new TTree("Chan_slow", "CHANNEL SLOW data");

	for (int run = RUN_START; run <=RUN_STOP; run++){
		vector<data_t> data = read_data_by_run(chan, run);
    	vector<Event> events = split_data_to_events(chan, data);
		vector<Event> nr_events = normalize_baseline(events);
		for (int event = 0; event < nr_events.size(); event++){
			SLOW = find_peaks(chan, nr_events[event], run);
			SLOW.baseline = get_baseline_info(events[event]);
			setup_branches(Chan_slow, SLOW);
			Chan_slow->Fill();
		}
		
	}
    file->Write(); 
    file->Close(); 
    // TODO: Discuss this thingy
 
    return 0;
}