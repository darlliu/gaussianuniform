#include"gucompute.h"

void gurunner::load (const std::string& s) {

    std::istringstream is(s);
    char ss[256];

    if (!is.getline( ss,256, '\t')) {
        std ::cerr << "Error reading line "<< s<<std::endl;
        return;
    }

    std::istringstream iss (ss);

    iss >> gene;

    data.clear(); // clear the data vector
    double d;
    while (is.getline(ss, 256,'\t')){
        std::istringstream iss1 (ss);
        iss1 >> d;
        data.push_back(d);
    }
    sz_data = data.size();
    //stub
    Pg.resize(sz_data);
    Pu.resize(sz_data);
    W.resize(sz_data);
    current_run_idx = 0; // reset index for new run
#if LOG2
    std::cout << "Loaded "<<gene <<std::endl;
#endif
    return;

};

void gurunner::get_member_likelihood (){

#if LOG
    std::cout << "Now assigning weights...";
#endif

    double c = 1/(sigma*pi2), c2 = 2*(sigma*sigma);

    //auto normpdf = [] (auto& d) {return c*exp(-1*(d-u)*(d-u)/c2); } ;
    // named lambdas with non concrete type, C++14
    auto normpdf = [&] (const double & d) {return c*exp(-1*(d-mu)*(d-mu)/c2); } ;

    auto uniformpdf = [&] (const double & d){
        if (d>=a && d<=b){
            return 1/(b-a);
        } else {
            return 0.0;
        }
    };
    for (int i=0; i<sz_data; ++i){
        Pg[i] = normpdf(data[i]);
        Pu[i] = uniformpdf(data[i]);
        W[i] = Pg[i]*w/ (Pg[i]*w + (1-w)*Pu[i]);
#if LOG
        std::cout << W[i] << ",";
#endif
    }
#if LOG
    std::cout << std::endl;
#endif

    return;
};


void gurunner::get_w(){
    w=0;
    for (auto &d: W){
        w+=d;
    }

    w /= sz_data;

#if LOG
    std::cout<< "Weight for Gaussian is "<< w << std::endl;

#endif
    return;
};

void gurunner::get_uni_params(){
    if (uniform_fixed) {
        return;
    }
    a=*std::min_element(data.begin(),data.end());
    b=*std::max_element(data.begin(),data.end());
    double mm=b, mx=a, d ;
    for (int i=0; i<sz_data; ++i){
        if (Pu[i]*(1-w) >= Pg[i]* w ){
            d= data[i];
            if (d<mm) {
                mm=d;
            } else if (d>mx){
                mx=d;
            }
        }
    }
    if (mm < mx){
        a=mm;
        b=mx;
    }
#if LOG
    std::cout<< "New a, b is "<< a  << " , "<< b << std::endl;
#endif
    return;
};


void gurunner::get_normal_params(){
    mu = 0;
    sigma = 0;
    for (int i =0; i < sz_data; ++i){
        mu+=W[i]*data[i];
        sigma+=W[i]*(data[i]-mu)*(data[i]-mu);
    }
    mu= mu/w/sz_data;
    sigma = sigma/w/sz_data;
    sigma = sqrt(sigma);
#if LOG
    std::cout << " New mu, sigma is " << mu <<" , " <<sigma <<std::endl;
#endif
    return;
};


void gurunner::get_total_likelihood(){
    res = 0;
    for (int i=0; i<sz_data; ++i){
        res += log(Pg[i]*w + Pu[i]*(1-w));
    }
#if LOG
    std::cout << " Log-Likelihood is now "<< res << std::endl;
#endif
    return;
};

void gurunner::run(double tol =1e-5, int mx = 500){

    double r=10;
    int i = 1;
    while (fabs(res-r)>tol && i < mx){
        r=res;
        get_member_likelihood();
        get_w();
        get_uni_params();
        get_normal_params();
        get_total_likelihood(); //res gets updated here
        ++i;
    }
#if LOG
    std::cout << "After " << i << " runs, likelihood at " << res << std::endl;
    std::cout << " New mu, sigma is " << mu <<" , " <<sigma <<std::endl;
    std::cout<< "New a, b is "<< a  << " , "<< b << std::endl;
    std::cout<< "Weight for Gaussian is "<< w << std::endl;
#endif
    return;
};


void gurunner::init(){
// get w
    std::normal_distribution<double> g1(0.5, 0.2);
    w = g1 (rndgen);
    //initialize around 0.5
// get a, b if not bounded
    a=*std::min_element(data.begin(),data.end());
    b=*std::max_element(data.begin(),data.end());
    if (!uniform_fixed) {
        std::uniform_real_distribution <double> g2(a, a+ (b-a)/2);
        std::uniform_real_distribution <double> g3(b-(b-a)/2, b);
        a = g2(rndgen);
        b = g3(rndgen);
        //tentative
    }
    mu = 0;
    for (auto & d: data){
        mu+=d;
    }
    mu/=sz_data;
    sigma=0;
    for (auto & d: data){
        sigma +=(d-mu)*(d-mu);
    }
    sigma/=sz_data;
    sigma = sqrt(sigma);
// get mu from conjugate prior
    std::normal_distribution <double> g4(mu, sigma);
    mu = g4(rndgen);
#if LOG
    std::cout << "Assigned random params, w, a, b, mu, sigma: " << w << "," <<\
        a<<","<<b<<","<<mu<<","<<sigma<<std::endl;
#endif
    return;
};

void gurunner::record(){
    if (current_run_idx>=restart_num) return;
    as[current_run_idx] = a;
    bs[current_run_idx] = b;
    ws[current_run_idx] = w;
    mus[current_run_idx] = mu;
    sigmas[current_run_idx] = sigma;
    ress[current_run_idx] = res;
    ++current_run_idx;
    return;
};

void gurunner::writeout(const unsigned& idx ){
#if LOG
    std::cout << "min element is "<<idx <<std::endl;
#endif
    std::ofstream f;
    f.open(prefix+"_EMResults.txt",std::ofstream::app);
    if (f.good()){
        f << prefix<<","<<gene << ","<<ws[idx]<<","<<as[idx]<<","<<bs[idx]\
            <<","<<mus[idx]<<","<<sigmas[idx]<<","<<ress[idx]<<","<<idx<<std::endl;
        f.flush();
        f.close();
    } else {
        std::cerr << " Failed to open file : "<<prefix <<"_EMResults.txt" <<std::endl;
    }
    return;
};

void gurunner::train(){
    for (int i = 0; i <restart_num; ++i){
        init();
        run();
        record();
    }
    auto mi = std::max_element(ress.begin(),ress.end());
    writeout((unsigned) (mi- ress.begin()));
    return;
};
