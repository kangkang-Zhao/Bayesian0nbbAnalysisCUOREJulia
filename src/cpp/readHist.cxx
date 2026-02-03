// readHist.cpp
#include <TFile.h>
#include <TH1.h>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: ./read_hist <file.root> <histname>\n";
        return 1;
    }

    const char* filename = argv[1];
    const char* histname = argv[2];

    TFile file(filename, "READ");
    if (file.IsZombie()) {
        std::cerr << "Error opening file " << filename << "\n";
        return 1;
    }

    TH1* h = dynamic_cast<TH1*>(file.Get(histname));
    if (!h) {
        std::cerr << "Histogram " << histname << " not found\n";
        return 1;
    }

    int nbins = h->GetNbinsX();
    std::cout << "# bin_low, bin_high, content\n";
    for (int i = 1; i <= nbins; i++) {
        double low = h->GetBinLowEdge(i);
        double high = h->GetBinLowEdge(i+1);
        double val = h->GetBinContent(i);
        std::cout << low << "," << high << "," << val << "\n";
    }

    return 0;
}


