//brew install nlohmann-json
//# or
//apt install nlohmann-json3-dev
#include <TFile.h>
#include <TFitResult.h>
#include <TMatrixTSym.h>

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

using json = nlohmann::json;

int main()
{
    TFile file(
        " /Users/zhaokangkang/Downloads/Bayesian0nbbAnalysis-2Tyr_brad/CUORE-0nbb-Analysis/ScalingFit/residual_and_width_vs_energy_ds3021.root",
        "READ"
    );

    // ---------- Q-value fit ----------
    TFitResult* fitQ =
        (TFitResult*)file.Get("fit_result_Q");
    TMatrixTSym<double>* covQ =
        (TMatrixTSym<double>*)file.Get("covariance_matrix_Q");

    json j;

    // Fit metadata
    j["Q"]["chi2"] = fitQ->Chi2();
    j["Q"]["ndf"]  = fitQ->Ndf();

    // Parameters
    int npar = fitQ->NPar();
    for (int i = 0; i < npar; i++)
    {
        j["Q"]["params"].push_back(fitQ->Parameter(i));
        j["Q"]["errors"].push_back(fitQ->ParError(i));
    }

    // Covariance matrix
    for (int i = 0; i < npar; i++)
    {
        std::vector<double> row;
        for (int jcol = 0; jcol < npar; jcol++)
            row.push_back((*covQ)(i, jcol));
        j["Q"]["cov"].push_back(row);
    }


    // ---------- Resolution sigma fit ----------
    TFitResult* fitSigma =
        (TFitResult*)file.Get("fit_result_sigma");
    TMatrixTSym<double>* covSigma =
        (TMatrixTSym<double>*)file.Get("covariance_matrix_sigma");

    // json j;

    // Fit metadata
    j["Sigma"]["chi2"] = fitSigma->Chi2();
    j["Sigma"]["ndf"]  = fitSigma->Ndf();

    // Parameters
    int nparSigma = fitSigma->NPar();
    for (int i = 0; i < nparSigma; i++)
    {
        j["Sigma"]["params"].push_back(fitSigma->Parameter(i));
        j["Sigma"]["errors"].push_back(fitSigma->ParError(i));
    }

    // Covariance matrix
    for (int i = 0; i < nparSigma; i++)
    {
        std::vector<double> rowSigma;
        for (int jcol = 0; jcol < nparSigma; jcol++)
            rowSigma.push_back((*covSigma)(i, jcol));
        j["Sigma"]["cov"].push_back(rowSigma);
    }


    // ---------- write JSON ----------
    std::ofstream out("lineshape_scaling_ds3021.json");
    out << j.dump(4);   // pretty print
    out.close();

    std::cout << "JSON file written successfully\n";
    return 0;
}

int converter(){
    main();
    return 1;
}
