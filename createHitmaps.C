#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TH2D.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <algorithm>

/**
 * @brief Compute row sums per layer (Y axis) for a 2D histogram.
 * 
 * @param h Pointer to a TH2D histogram.
 * @return std::vector<double> Vector of sums per row/layer.
 */
static std::vector<double> rowSums(TH2D* h) {
    int nx = h->GetNbinsX(), ny = h->GetNbinsY();
    std::vector<double> s(ny, 0.0);
    for (int ly = 1; ly <= ny; ++ly)
        for (int fx = 1; fx <= nx; ++fx)
            s[ly-1] += h->GetBinContent(fx, ly);
    return s;
}

/**
 * @brief Apply per-layer multiplicative weights to a 2D histogram.
 * 
 * Optionally applies a special boost for the first layer.
 * 
 * @param h Pointer to a TH2D histogram.
 * @param w Vector of weights per layer.
 */
static void applyRowWeights(TH2D* h, const std::vector<double>& w) {
    int nx = h->GetNbinsX(), ny = h->GetNbinsY();
    for (int ly = 1; ly <= ny; ++ly) {
        double scale = w[ly-1];
        if (scale == 1.0) continue;
        for (int fx = 1; fx <= nx; ++fx) {
            if (ly == 1) {
                h->SetBinContent(fx, ly, h->GetBinContent(fx, ly) * scale * 1.0658);
            } else {
                h->SetBinContent(fx, ly, h->GetBinContent(fx, ly) * scale);
            }
        }
    }
}

/**
 * @brief Generate physics-motivated exponential weights per layer.
 * 
 * @param nLayers Number of layers in the histogram.
 * @param beta Attenuation coefficient (default 0.3).
 * @return std::vector<double> Vector of exponential weights.
 */
static std::vector<double> weightsExponential(int nLayers, double beta = 0.3) {
    std::vector<double> w(nLayers, 1.0);
    for (int l = 0; l < nLayers; ++l) {
        w[l] = std::exp(-(beta * l * 2)*(beta * l * 2));
    }
    return w;
}

/**
 * @brief Generate monotonic weights to enforce decreasing counts per layer.
 * 
 * Useful after boosting first layer to avoid unphysical increase in later layers.
 * 
 * @param sums Vector of row sums per layer.
 * @param eps Small number to prevent division by zero.
 * @return std::vector<double> Monotonically decreasing weights.
 */
static std::vector<double> weightsMonotonic(const std::vector<double>& sums, double eps=1e-12) {
    std::vector<double> w(sums.size(), 1.0);
    if (sums.empty()) return w;

    double last = sums[0];
    for (size_t l = 1; l < sums.size(); ++l) {
        if (sums[l] > last) {
            w[l] = last / (sums[l] + eps);
        }
        last = sums[l] * w[l];
    }
    return w;
}

/**
 * @brief Normalize a histogram so that its total integral is 1.
 * 
 * @param h Pointer to TH2D histogram.
 */
static void normalizeTotal(TH2D* h) {
    double total = h->Integral();
    if (total > 0) h->Scale(1.0 / total);
}

/**
 * @brief Normalize a histogram such that the first layer sum is 1.
 * 
 * @param h Pointer to TH2D histogram.
 */
static void normalizeFirstLayer(TH2D* h) {
    double sumL0 = 0;
    int nx = h->GetNbinsX();
    for (int fx = 1; fx <= nx; ++fx) sumL0 += h->GetBinContent(fx, 1);
    if (sumL0 > 0) h->Scale(1.0 / sumL0);
}

// Function declarations
TString processFile(const TString& filepath, const TString& outputFilePath, bool energyDepo);
void processAllRootFilesInDir(const TString& inputDir, const TString& outputDir, bool energyDepo);

/**
 * @brief Process all ROOT files in a directory and create hitmap histograms.
 * 
 * @param inputDir Input directory containing ROOT files.
 * @param outputDir Output directory to save hitmap ROOT files.
 * @param energyDepo If true, fill energy deposition instead of simple counts.
 */
void processAllRootFilesInDir(const TString& inputDir, const TString& outputDir, bool energyDepo) {
    // Ensure output directory exists
    if (gSystem->AccessPathName(outputDir, kWritePermission)) {
        gSystem->mkdir(outputDir, true);
    }

    void* dirp = gSystem->OpenDirectory(inputDir);
    if (!dirp) {
        std::cerr << "❌ Cannot open directory: " << inputDir << std::endl;
        return;
    }

    const char* fileName = nullptr;
    while ((fileName = gSystem->GetDirEntry(dirp)) != nullptr) {
        TString fname(fileName);
        if (!fname.EndsWith(".root")) continue;

        TString inputFilePath = inputDir + "/" + fname;
        TString baseName = fname;
        baseName.ReplaceAll(".root", "");
        TString outputFilePath = outputDir + "/" + baseName + "_hitmap.root";

        std::cout << "➡️ Processing: " << inputFilePath << std::endl;
        processFile(inputFilePath, outputFilePath, energyDepo);
    }
    gSystem->FreeDirectory(dirp);
}

/**
 * @brief Create hitmaps from a single ROOT file and save to output file.
 * 
 * @param filepath Path to input ROOT file.
 * @param outputFilePath Path to save output hitmap ROOT file.
 * @param energyDepo Fill energy deposition (true) or counts (false).
 * @return TString Path to output file.
 */
TString processFile(const TString& filepath, const TString& outputFilePath, bool energyDepo = false) {
    TFile* file = TFile::Open(filepath);
    if (!file || file->IsZombie()) {
        std::cerr << "❌ Cannot open file: " << filepath << std::endl;
        return "";
    }

    TTree* tree = (TTree*)file->Get("deposits");
    if (!tree) {
        std::cerr << "⚠️ No tree named 'deposits' in: " << filepath << std::endl;
        file->Close();
        return "";
    }

    // Branch variables
    Int_t fiberID = 0;
    TVector3* position = nullptr;
    double_t energy = 0.0;

    tree->SetBranchAddress("position", &position);
    tree->SetBranchAddress("fiberID", &fiberID);
    tree->SetBranchAddress("energy", &energy);

    // Detector geometry
    const int nLayers = 7;
    const int nFibers = 55;
    double layerThickness = 2;
    double zMin = 226.0;
    double layerZ[nLayers];
    for (int l = 0; l < nLayers; ++l) {
        layerZ[l] = zMin + l * layerThickness + layerThickness / 2.0;
    }
    double halfThickness = layerThickness / 2.0;

    // Energy thresholds (MeV)
    std::vector<double> thresholds = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5};

    // Create histograms per threshold
    std::vector<TH2D*> hHitmaps;
    for (double thr : thresholds) {
        int thrInt = (int)(thr * 1000);
        TString thrStr = Form("%04d", thrInt);
        TString name = energyDepo ? Form("hEnergyDepoMap_t%s", thrStr.Data())
                                  : Form("hFiberHitMap_t%s", thrStr.Data());
        TString title = Form("%s (E > %.1f);FiberID;Layer",
                             energyDepo ? "Energy deposition map" : "Fiber hitmap", thr);
        hHitmaps.push_back(new TH2D(name, title, nFibers, 0, nFibers, nLayers, 0, nLayers));
    }

    // Fill histograms
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (!position) continue;

        for (size_t t = 0; t < thresholds.size(); ++t) {
            if (energy > thresholds[t]) {
                for (int l = 0; l < nLayers; ++l) {
                    if (std::abs(position->Z() - layerZ[l]) <= halfThickness) {
                        if (energyDepo) hHitmaps[t]->Fill(fiberID, l, energy);
                        else hHitmaps[t]->Fill(fiberID, l);
                        break;
                    }
                }
            }
        }
    }

    // Save histograms
    TFile* fout = new TFile(outputFilePath, "RECREATE");
    for (auto* h : hHitmaps) {
        h->Write();
        auto sums = rowSums(h);

        // Variant A: exponential attenuation
        TH2D* hExp = (TH2D*)h->Clone(TString(h->GetName()) + "_expAtt");
        applyRowWeights(hExp, weightsExponential(nLayers, 0.05));
        hExp->Write();
    }

    fout->Close();
    file->Close();

    std::cout << "✅ Saved: " << outputFilePath << std::endl;
    return outputFilePath;
}
