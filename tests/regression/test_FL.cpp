// Regression test: the strict conformal estimator must reproduce the frozen
// reference values at the two benchmark points. Regression, not independent
// physics validation (see tests/physics and docs/validation.md).
#include <SMVacuumDecay/FubiniLipatov.hpp>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

int main() {
    // Standard Model Point
    double SM_Mh = 125.1;
    double SM_Mt = 173.1;
    double SM_ref_conf = 2120.340020693041;
    int SM_ref_Status = 2;

    // Low-Mass Point (deeply metastable toy region of the scan range)
    double Low_Mh = 5.0;
    double Low_Mt = 105.0;
    double Low_ref_conf = 11571.930525635938;
    int Low_ref_Status = 2;

    auto sm_res = SMVacuumDecay::classify_conformal(SM_Mh, SM_Mt);
    int sm_status = std::get<0>(sm_res);
    double sm_action = std::get<1>(sm_res);

    auto low_res = SMVacuumDecay::classify_conformal(Low_Mh, Low_Mt);
    int low_status = std::get<0>(low_res);
    double low_action = std::get<1>(low_res);

    double err_sm = std::abs((sm_action - SM_ref_conf)/SM_ref_conf);
    double err_low = std::abs((low_action - Low_ref_conf)/Low_ref_conf);

    std::cout << std::scientific << std::setprecision(12);
    std::cout << "SM Error: " << err_sm << " Status: " << sm_status << " (Ref: " << SM_ref_Status << ")\n";
    std::cout << "Low Error: " << err_low << " Status: " << low_status << " (Ref: " << Low_ref_Status << ")\n";

    double max_err = std::max(err_sm, err_low);
    std::cout << "Conformal Regression Max Error: " << max_err << std::endl;

    if (err_sm < 1e-12 && err_low < 1e-12 && sm_status == SM_ref_Status && low_status == Low_ref_Status) {
        std::cout << "PASS" << std::endl;
        return 0;
    } else {
        std::cout << "FAIL" << std::endl;
        return 1;
    }
}
