#include<vector>
#include<TSpline.h>

std::array<double, 3> RoundRobin(const std::vector<double>& new_points, const std::vector<double>& genie_weights){
    double passed=0.0;
    auto const nweights = genie_weights.size();

    std::vector<double> test_points; //holds all test points
    std::vector<double> new_test_deviation;

    for(std::size_t skip = 1; skip < nweights - 1; ++skip){ //Create all N-1 node points (N-2 Total Splines)
        std::vector<double> weights;
        std::vector<double> points;
        for(std::size_t i = 0; i < nweights; ++i){
            if(i!=skip){
                weights.push_back(genie_weights[i]);
                points.push_back(new_points[i]);
            }
            else{
                test_points.push_back(new_points[i]);
            }
        }
        const char* opt=0;
        const TSpline3 spline("spline", &points[0], &weights[0], points.size(), opt); //take gaps and eval at that gap
        // OR
        // const TSpline3 spline("spline", points.data(), weights.data(), points.size(), opt); //take gaps and eval at that gap

		const double test_val = spline.Eval(test_points[skip-1]);
        auto const next_weight = genie_weights[skip+1];
        if( next_weight==0 )
        {
            new_test_deviation.push_back(std::abs(test_val));
            continue;
        }
        double ratio = std::abs(test_val-next_weight)/std::abs(next_weight);
        new_test_deviation.push_back(ratio);
    }

    //checkif max dev is above 0.01
    // And the position of the max_dev
    double max_dev_new_rr = new_test_deviation[0];
    double max_dev_pos_new_rr = test_points[0];
    for(std::size_t i = 1; i < test_points.size(); i++){
        double current_dev_new_rr = new_test_deviation[i];
        if(current_dev_new_rr > max_dev_new_rr){
            max_dev_new_rr = current_dev_new_rr;
            max_dev_pos_new_rr = test_points[i];
        }
    }
    if(max_dev_new_rr < 0.01)
    {
        passed=1.0;
    }

    std::array<double, 3> return_vector = {passed, max_dev_pos_new_rr, max_dev_new_rr}; //return if passed, give the param (x-value) of highest dev
    return return_vector;
}

// std::vector<double> RoundRobin(std::vector<double> new_points, std::vector<double> genie_weights){    
//     double passed=0.0;
//     std::vector<double> test_points; //holds all test points
//     std::vector<std::vector<double>> skippable_genie_weights; //holds genie weight with gaps
//     std::vector<std::vector<double>> skippable_genie_param; //holds genie params with gaps
//     for(int skip=1; skip<genie_weights.size()-1; skip++){ //Create all N-1 node points (N-2 Total Splines)
// 	    std::vector<double> temp_genie;
// 	    std::vector<double> temp_point;
//         for(int i=0; i<genie_weights.size(); i++){
//             if(i!=skip){
//                 temp_genie.push_back(genie_weights[i]);
//                 temp_point.push_back(new_points[i]);
//             }
//             else{
//                 test_points.push_back(new_points[i]);
//                 //cout << new_points[i] << endl;
//             }
//         }
//         skippable_genie_weights.push_back(temp_genie);
//         skippable_genie_param.push_back(temp_point);
//     }
//     double skip_genie_weight_arr[skippable_genie_weights.size()][skippable_genie_weights[0].size()];
//     double skip_genie_param_arr[skippable_genie_weights.size()][skippable_genie_weights[0].size()];
//     for(int i=0; i<skippable_genie_weights.size(); i++){
//         for(int j=0; j<skippable_genie_weights[0].size(); j++){
//             skip_genie_weight_arr[i][j]=skippable_genie_weights[i][j];
//             skip_genie_param_arr[i][j]=skippable_genie_param[i][j];
//         }
//     }

//     //How are the vectors organized as of now?
//     //genie_weight--> all genie weights for all points -3 to 3
//     //skippable_genie_weight--> first vec is first (excluding -3) point pissing, second is the second, and so on till +3(excluding)
//     //skippable_genie_points--> first vec is first point missing(excluding 3), till 3(excluding)
//     std::vector<double> new_test_deviation;
//     std::vector<double> new_test_points;
//     for(int i=0; i<skippable_genie_weights.size(); i++){
//         const char* opt=0;
//         TSpline3* spline = new TSpline3("spline", skip_genie_param_arr[i], skip_genie_weight_arr[i], skippable_genie_param[i].size(), opt); //take gaps and eval at that gap
//         double test_val = spline->Eval(test_points[i]);
//         double ratio = abs(test_val-genie_weights[i+1])/abs(genie_weights[i+1]);
//         if(genie_weights[i+1]==0){
//             ratio=abs(test_val);
//         }
//         new_test_deviation.push_back(ratio);
//         delete spline;
//     }
//     //checkif max dev is above 0.01
//     // And the position of the max_dev
//     double max_dev_new_rr=new_test_deviation[0];
//     double max_dev_pos_new_rr = test_points[0];
//     for(int i=1; i<test_points.size(); i++){
//         double current_dev_new_rr=new_test_deviation[i];
//         if(current_dev_new_rr>max_dev_new_rr){
//             max_dev_new_rr=current_dev_new_rr;
//             max_dev_pos_new_rr=test_points[i];
//         }
//     }
//     if(max_dev_new_rr<0.01){
//         passed=1.0;
//     }
//     std::vector<double> return_vector = {passed, max_dev_pos_new_rr};//,max_dev_new_rr}; //return if passed, give the param (x-value) of highest dev
//     return return_vector;
// }
