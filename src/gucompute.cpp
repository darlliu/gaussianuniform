#include"gucompute.h"

void gurunner::load (const char* fname) {
    sz_data=5;
    //stub
    data.resize(sz_data);
    Pg.resize(sz_data);
    Pu.resize(sz_data);
    W.resize(sz_data);

    data[0]=1;
    data[1]=2;
    data[3]=1;
    data[2]=50;
    data[4]=10;
    a=*std::min_element(data.begin(),data.end());
    b=*std::max_element(data.begin(),data.end());

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
        if (Pu[i]*(1-w) > Pg[i]* w ){
            d= data[i];
            if (d<mm) {
                mm=d;
            } else if (d>mx){
                mx=d;
            }
        }
    }
#if LOG
    std::cout<< "New a, b is "<<mm << " , "<<mx << std::endl;
#endif
    a=mm;
    b=mx;
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
#endif
    return;
};
