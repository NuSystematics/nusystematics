#include "RwFramework/GReWeight.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepRecord.h"
#include <vector>

using namespace std;
using namespace genie;
using namespace genie::rew;

std::vector<double> Flatness(genie::EventRecord const &GenieGHep, systtools::paramId_t pid, std::unique_ptr<genie::rew::GReWeight> & grw, std::vector<double> node, std::vector<double> weights){//const genie::EventRecord & GenieGHep, genie::rew::GReWeight & rw, genie::rew::GSystSet & syst, double non_flat_point, double non_flat_weight, double first_flat_point, double first_flat_weight){
    bool is_flat=false;
    double genie_non_flat_point;
    double genie_non_flat_weight;
    double genie_first_flat_point;
    double genie_first_flat_weight;
    //Point instance of first flat weight
    for(int i=1; i<weights.size(); i++){
        if(weights[i-1]==weights[i] && is_flat==false){
            genie_non_flat_point=node[i-2];
            genie_non_flat_weight=weights[i-2];
            genie_first_flat_point=node[i-1];
            genie_first_flat_weight=weights[i-1];
            is_flat=true;
        }
    }
    //if point is on left side i.e. less than 0, then recalc
    if(is_flat==true && genie_first_flat_point<0){
        bool no_longer_flat=false;
        for(int i=1; i<weights.size(); i++){
            if(weights[i-1]!=weights[i] && no_longer_flat==false){
                genie_non_flat_point=node[i];
                genie_non_flat_weight=weights[i];
                genie_first_flat_point=node[i-1];
                genie_first_flat_weight=weights[i-1];
                no_longer_flat=true;
            }
        }
    }

    if(is_flat==true){
    
    // double genie_non_flat_point = non_flat_point;
    // double genie_first_flat_point = first_flat_point;
    // double genie_non_flat_weight = non_flat_weight;
    // double genie_first_flat_weight = first_flat_weight;
        vector<double> test_points; //plot line
        double spacing=0.5*abs(genie_first_flat_point-genie_non_flat_point);
        while(spacing>=0.1){
            spacing = 0.5*abs(genie_first_flat_point-genie_non_flat_point);
            double point_in_middle;
            if(genie_non_flat_point>0){ //Flatness on rightside
                point_in_middle = genie_non_flat_point+spacing;
            }
            else{ //Flatness on leftside
                point_in_middle = genie_non_flat_point-spacing;
            }
            grw->Systematics().Set(grw->Systematics().AllIncluded()[0], point_in_middle); //Generate GENIE RW
            grw->Reconfigure();
            double test_middle_weight=grw->CalcWeight(GenieGHep);//Calculate Point
            if(test_middle_weight==genie_first_flat_weight){
            genie_first_flat_weight=test_middle_weight;
            genie_first_flat_point=point_in_middle;
            }
            else{
                genie_non_flat_weight=test_middle_weight;
                genie_non_flat_point=point_in_middle;
            }
            test_points.push_back(point_in_middle); //plot line
            test_points.push_back(test_middle_weight); //plot line
        }
        vector<double> return_vec = {genie_non_flat_point, genie_non_flat_weight};
        return return_vec;
        }
    else{
        // will return a vector of size 1 if never flat
        vector<double> np_flat_vec={1};
        return np_flat_vec;
    }
}
